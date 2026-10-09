#pragma once

/** @file compressed.h
 * @brief Compressed rank/select support with metadata-validated source
 * prefetch.
 */

#include <pixie/rank_select.h>
#include <pixie/rank_select/detail/compressed/packed_prefix_search.h>
#include <pixie/rank_select/detail/compressed/position_samples.h>
#include <pixie/rank_select/detail/compressed/samples.h>
#include <pixie/rank_select/detail/compressed/select_loaded_block.h>
#include <pixie/rank_select/detail/compressed/select_metadata_search.h>
#include <pixie/storage/aligned.h>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie {
/**
 * @brief Non-owning rank/select index with interleaved compressed metadata.
 * @details Source words must remain alive and immutable for the index lifetime.
 * Copying duplicates metadata and continues to borrow the same source words.
 * Moving transfers metadata and leaves the source empty with select disabled.
 * Bit numbering, half-open rank prefixes, one-based select, and sentinels
 * follow RankSelectBase. Rank boundaries beyond size() clamp to the total
 * count.
 *
 * Each 16,384-bit superblock has one 64-byte record: a global one rank
 * followed by 32 monotone prefixes packed uniformly in fourteen bits. A zero
 * terminal encoding represents either zero or 16,384 ones: the first prefix
 * is then exactly zero or 512, so multiplying it by 32 reconstructs the total.
 * All 448 local bits fit beside the global rank. The default AlignedStorage
 * aligns every rank record to 64 bytes; another backend need only align its
 * allocation to uint64_t. On little-endian hosts each
 * nonterminal prefix decodes from one unaligned 32-bit read within the record;
 * other hosts use word extraction. Full-superblock rank overhead is 3.125%.
 *
 * Select samples store superblock IDs every SampleFrequency occurrences, using
 * block Elias--Fano. During construction, neighboring samples delimit rank
 * searching. A span below SearchThreshold uses a short forward rank scan.
 * Other intervals have
 * a separate sparse directory, compressed occupied-superblock IDs and bounded
 * relative rank endpoints; select searches those endpoints without probing
 * interleaved rank records, then decodes just the chosen superblock ID.
 * If every wide interval has contiguous occupied superblocks, no occupied-ID
 * array is allocated; all target IDs follow arithmetically from their samples.
 * Mixed directions store a 32-bit wide-interval bitmap and prefix count per
 * 32 samples. A bit test and popcount locate the exact sparse-directory entry;
 * queries never binary-search sample IDs. Endpoint selection uses lower_bound.
 * Ordinary samples carry no child pointers. A
 * terminal sample follows the last occupied superblock, excluding trailing
 * empty blocks. Construction classifies each direction as homogeneous, narrow,
 * wide-prefix, or mixed. If all superblocks are homogeneous and sampling
 * divides their full population, the starting sample gives an arithmetic
 * answer. Such indexes allocate neither occupied IDs nor endpoint exceptions.
 * The retained arithmetic shortcut obtains the within-superblock offset from
 * (rank - 1) modulo 16,384, avoiding a rank-record fetch on that path.
 * Proven shapes need no bitmap allocation. Wide endpoint search probes
 * existing group ends, then compares up to 64 bytes together; no group index
 * is stored. Singleton intervals need no endpoint load. Endpoint values are
 * zero-based, allowing a sampling frequency of 65,536 to use sixteen bits.
 * Indexes of at most 1,024 bits allocate no select samples. An outlined path
 * chooses between at most two basic blocks using their local prefix; a fully
 * matching index returns rank - 1 directly.
 *
 * Dense, nonhomogeneous sources replace coarse primary samples with exact
 * selected bit positions in fixed-stride Elias--Fano groups of 32. Each group
 * interleaves its base, high bitmap, optional wide-interval map and low fields;
 * a shared maximum low width pads groups uniformly without an offset directory.
 * There is no duplicate primary sample array. Every such query decodes both
 * neighboring positions, predicts a 512-bit source block and prefetches its
 * words while validating against rank metadata. Exact validation may correct
 * the estimate by one adjacent basic block; otherwise the original narrow
 * scan or wide-directory lookup finds the answer. Prediction is skipped when
 * the sampled bit span exceeds 65,536. Sources with either bit population below
 * floor(size()/16), and homogeneous directions, retain coarse samples and their
 * arithmetic shortcuts. The terminal position is one past the last selected
 * bit, excluding trailing empty bits. Speculation never changes correctness.
 *
 * On AVX-512 VBMI targets, local search reconstructs all packed prefixes with
 * two byte permutations, then compares them together. Other targets use
 * interpolation with exact correction bounds. Homogeneous blocks return an
 * arithmetic position. The 512-bit select primitive
 * keeps source words in registers through word localization; other targets
 * use the established SIMD/scalar fallback. Construction builds rank records
 * from source bits, then samples and exceptions from those records. Exact
 * position samples additionally resolve selected source blocks. The index owns
 * metadata, not source bits. All retained metadata arrays use MetadataStorage;
 * construction may use temporary vectors. Serialization is not provided.
 * @tparam MetadataStorage Default-empty, contiguous storage for index metadata.
 * Construction requires resize(bit_count) and writable_bytes(); read access
 * uses as_bytes(), and memory accounting uses allocated_bytes(). Allocations
 * must expose exactly the requested bytes with uint64_t alignment. The storage
 * must own its bytes, preserve contents when copied, and transfer ownership
 * when moved. Its alignment does not constrain the borrowed source words.
 * @tparam SampleFrequency Positive power-of-two occurrence sampling interval.
 * @tparam SearchThreshold Positive superblock-span threshold for extra
 * indexing.
 */
template <StorageImplementation MetadataStorage = AlignedStorage,
          std::size_t SampleFrequency = 4096,
          std::size_t SearchThreshold = 16>
class RankSelectSupportCompressed
    : public RankSelectBase<RankSelectSupportCompressed<MetadataStorage,
                                                        SampleFrequency,
                                                        SearchThreshold>> {
  static_assert(std::has_single_bit(SampleFrequency));
  static_assert(SearchThreshold > 0);
  static constexpr bool kNothrowMove =
      std::is_nothrow_default_constructible_v<MetadataStorage> &&
      std::is_nothrow_move_constructible_v<MetadataStorage> &&
      std::is_nothrow_move_assignable_v<MetadataStorage>;

 public:
  /** @brief Optional select directions; rank is always indexed. */
  enum class SelectSupport : std::uint8_t {
    kNone = 0,
    kSelect1 = 1,
    kSelect0 = 2,
    kBoth = 3,
  };

  /** @brief Construct an empty index with no select directions enabled. */
  RankSelectSupportCompressed() = default;

  /** @brief Duplicate metadata while borrowing the same immutable source. */
  RankSelectSupportCompressed(const RankSelectSupportCompressed&) = default;

  /** @brief Transfer metadata, leaving the source empty with select disabled.
   */
  RankSelectSupportCompressed(RankSelectSupportCompressed&& other) noexcept(
      kNothrowMove) {
    swap_state(other);
  }

  /**
   * @brief Duplicate metadata while borrowing the assigned source.
   * @details Allocation failure while copying metadata leaves this index
   * unchanged. The strong exception guarantee also covers assignment when
   * MetadataStorage has nonthrowing move construction and move assignment.
   */
  RankSelectSupportCompressed& operator=(
      const RankSelectSupportCompressed& other) {
    if (this != &other) {
      RankSelectSupportCompressed copy(other);
      swap_state(copy);
    }
    return *this;
  }

  /**
   * @brief Transfer metadata, leaving the source empty with select disabled.
   * @details Self-move assignment preserves the index.
   */
  RankSelectSupportCompressed& operator=(
      RankSelectSupportCompressed&& other) noexcept(kNothrowMove) {
    if (this != &other) {
      RankSelectSupportCompressed moved(std::move(other));
      swap_state(moved);
    }
    return *this;
  }

  /**
   * @brief Index caller-owned packed words, ignoring padding after valid bits.
   * @param source_words External immutable complete 64-bit words.
   * @param num_bits Valid bit count, clamped to the available source capacity.
   * @param select_support Directions to build; disabled select returns size().
   * @throws std::invalid_argument If select_support is not a defined mode.
   * Also thrown if the metadata backend exposes an invalid size or alignment.
   * @throws std::length_error If an auxiliary allocation size is
   * unrepresentable.
   * @details Allocation failure propagates. Select rank zero returns zero even
   * when its direction is disabled. The source lifetime is caller-managed.
   */
  explicit RankSelectSupportCompressed(
      std::span<const std::uint64_t> source_words,
      std::size_t num_bits,
      SelectSupport select_support = SelectSupport::kBoth)
    requires detail::compressed::MutableMetadataStorage<MetadataStorage>
      : bits_(source_words),
        num_bits_(num_bits),
        select_support_(select_support) {
    if (static_cast<std::uint8_t>(select_support) > 3) {
      throw std::invalid_argument("Invalid compressed select support mode");
    }
    if (source_words.size() <= std::numeric_limits<std::size_t>::max() / 64) {
      num_bits_ = std::min(num_bits_, source_words.size() * 64);
    }
    build_rank();
    if (supports_select1_impl()) {
      select1_ = build_samples<false>();
    }
    if (supports_select0_impl()) {
      select0_ = build_samples<true>();
    }
  }

  /** @brief Return the logical source bit count. */
  std::size_t size_impl() const { return num_bits_; }

  /** @brief Read a zero-based bit; position must be below size(). */
  int bit_impl(std::size_t position) const {
    return (bits_[position / 64] >> (position % 64)) & 1;
  }

  /** @brief Count ones in [0, end), clamping end at size(). */
  // Keep the shared prefix calculation inside the directional query.
  [[gnu::always_inline]] inline std::uint64_t rank_impl(std::size_t end) const {
    if (end >= num_bits_) {
      return one_count_;
    }
    const std::size_t super = end / kSuperBits;
    const std::size_t basic = (end % kSuperBits) / kBasicBits;
    const auto* record = rank_record(super);
    std::uint64_t result = record[0] + local_rank(record, basic);
    const std::size_t first_word = (end / kBasicBits) * 8;
    const std::size_t offset = end % kBasicBits;
    if (offset == 0) {
      return result;
    }
    if (bits_.size() - first_word >= 8) {
      return result + ::rank_512(bits_.data() + first_word, offset);
    }
    const std::size_t last_word = end / 64;
    for (std::size_t word = first_word; word < last_word; ++word) {
      result += std::popcount(bits_[word]);
    }
    if (end % 64 != 0) {
      result += std::popcount(bits_[last_word] & ::first_bits_mask(end % 64));
    }
    return result;
  }

  /** @brief Count zeros in [0, end), clamping to the logical bit length. */
  std::uint64_t rank0_impl(std::size_t end) const {
    if (end >= num_bits_) {
      return num_bits_ - one_count_;
    }
    return end - rank_impl(end);
  }

  /**
   * @brief Locate a one-based one occurrence as a zero-based bit position.
   * @return Zero for rank zero; size() if absent or one-select is disabled.
   */
  std::uint64_t select_impl(std::size_t rank) const {
    return select_direction<false>(rank, select1_, supports_select1_impl());
  }

  /**
   * @brief Locate a one-based zero occurrence as a zero-based bit position.
   * @return Zero for rank zero; size() if absent or zero-select is disabled.
   */
  std::uint64_t select0_impl(std::size_t rank) const {
    return select_direction<true>(rank, select0_, supports_select0_impl());
  }

  /** @brief Report whether one-select metadata was requested. */
  bool supports_select1_impl() const {
    return (static_cast<std::uint8_t>(select_support_) & 1) != 0;
  }

  /** @brief Report whether zero-select metadata was requested. */
  bool supports_select0_impl() const {
    return (static_cast<std::uint8_t>(select_support_) & 2) != 0;
  }

  /** @brief Object plus allocated auxiliary bytes, excluding borrowed words. */
  std::size_t memory_usage_bytes_impl() const
    requires requires(const MetadataStorage& storage) {
      storage.allocated_bytes();
    }
  {
    return sizeof(*this) + ranks_.allocated_bytes() +
           select1_.allocated_bytes() + select0_.allocated_bytes();
  }

 private:
  static constexpr std::size_t kSuperBits = 16384;
  static constexpr std::size_t kBasicBits = 512;
  static constexpr std::size_t kRecordBytes = 64;
  // An interval contains at most ceil(SIZE_MAX / 2^14) superblocks, leaving
  // the high size_t bit available for this contiguity marker.
  static constexpr std::size_t kContiguous =
      std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);

  struct WideEntry {
    std::size_t begin;
    std::size_t count;
  };
  using PrimarySamples = detail::compressed::BitmapMonotoneSamples<
      detail::compressed::BlockedMonotoneSamples<MetadataStorage>>;

  using PositionSamples = detail::compressed::PositionSamples<MetadataStorage>;
  template <class T>
  using MetadataArray = detail::compressed::MetadataArray<T, MetadataStorage>;

  struct SelectIndex {
    enum class Shape : std::uint8_t {
      kMixed,
      kNarrow,
      kWidePrefix,
      kHomogeneous,
    };
    using Endpoint = std::conditional_t<
        (SampleFrequency - 1 <= std::numeric_limits<std::uint16_t>::max()),
        std::uint16_t,
        std::conditional_t<(SampleFrequency - 1 <=
                            std::numeric_limits<std::uint32_t>::max()),
                           std::uint32_t,
                           std::uint64_t>>;
    PrimarySamples samples;
    PositionSamples positions;
    bool predicted = false;
    detail::compressed::BlockedMonotoneSamples<MetadataStorage> occupied;
    MetadataArray<WideEntry> wide;
    MetadataArray<Endpoint> endpoints;
    Shape shape{};

    std::size_t allocated_bytes() const {
      return samples.allocated_bytes() + positions.allocated_bytes() +
             occupied.allocated_bytes() + wide.allocated_bytes() +
             endpoints.allocated_bytes();
    }
  };

  std::span<const std::uint64_t> bits_;
  std::size_t num_bits_ = 0;
  std::uint64_t one_count_ = 0;
  MetadataArray<std::uint64_t> ranks_;
  SelectIndex select1_;
  SelectIndex select0_;
  SelectSupport select_support_ = SelectSupport::kNone;

  void swap_state(RankSelectSupportCompressed& other) noexcept(kNothrowMove) {
    using std::swap;
    swap(bits_, other.bits_);
    swap(num_bits_, other.num_bits_);
    swap(one_count_, other.one_count_);
    swap(ranks_, other.ranks_);
    swap(select1_, other.select1_);
    swap(select0_, other.select0_);
    swap(select_support_, other.select_support_);
  }

  std::size_t super_count() const { return ranks_.size() / (kRecordBytes / 8); }

  const std::uint64_t* rank_record(std::size_t super) const {
    // Complete word-aligned records are guaranteed by build_rank(). Avoid
    // repeating storage alignment/size validation in query paths.
    return ranks_.data() + super * (kRecordBytes / 8);
  }

  static std::uint64_t local_rank(const std::uint64_t* record,
                                  std::size_t blocks) {
    if (blocks == 0) {
      return 0;
    }

    if (blocks == 32) {
      const auto terminal = record[7] >> 50;
      // With a zero terminal encoding, the first prefix is exactly 0 or
      // 512; shifting it reconstructs 0 or 16384 without a second test.
      return terminal != 0 ? terminal : (record[1] & 0x3fff) << 5;
    }
    const std::size_t bit = (blocks - 1) * 14;
    constexpr auto mask = 0x3fffu;
    if constexpr (std::endian::native == std::endian::little) {
      std::uint32_t packed;
      // The largest byte offset is 52. A four-byte load stays inside the seven
      // packed words even at the final nonterminal prefix (block 31).
      std::memcpy(&packed,
                  reinterpret_cast<const std::byte*>(record + 1) + bit / 8,
                  sizeof(packed));
      return (packed >> (bit % 8)) & mask;
    } else {
      auto packed = record[1 + bit / 64] >> (bit % 64);
      constexpr unsigned width = 14;
      if (bit % 64 + width > 64) {
        packed |= record[2 + bit / 64] << (64 - bit % 64);
      }
      return packed & mask;
    }
  }

  static void write_local_rank(std::uint64_t* record,
                               std::size_t blocks,
                               std::uint64_t value) {
    const std::size_t bit = (blocks - 1) * 14;
    constexpr unsigned width = 14;
    value &= 0x3fff;
    record[1 + bit / 64] |= value << (bit % 64);
    if (bit % 64 + width > 64) {
      record[2 + bit / 64] |= value >> (64 - bit % 64);
    }
  }

  std::size_t super_size(std::size_t super) const {
    return std::min(kSuperBits, num_bits_ - super * kSuperBits);
  }

  void build_rank() {
    const std::size_t count =
        num_bits_ == 0 ? 0 : 1 + (num_bits_ - 1) / kSuperBits;
    ranks_ = MetadataArray<std::uint64_t>(count * (kRecordBytes / 8));
    auto* words = ranks_.data();
    for (std::size_t super = 0; super < count; ++super) {
      auto* record = words + super * (kRecordBytes / 8);
      record[0] = one_count_;
      const std::size_t base = super * kSuperBits;
      const std::size_t valid = super_size(super);
      std::uint64_t local = 0;

      for (std::size_t basic = 0; basic < 32; ++basic) {
        for (std::size_t lane = 0; lane < 8; ++lane) {
          const std::size_t offset = basic * kBasicBits + lane * 64;
          if (offset >= valid) {
            break;
          }
          const auto word = bits_[base / 64 + offset / 64];
          local += std::popcount(word & ::first_bits_mask(std::min<std::size_t>(
                                            64, valid - offset)));
        }

        write_local_rank(record, basic + 1, local);
      }

      one_count_ += local;
    }
  }

  template <bool Zero>
  std::uint64_t prefix_before(std::size_t super) const {
    const std::uint64_t ones = rank_record(super)[0];
    if constexpr (Zero) {
      return super * kSuperBits - ones;
    }
    return ones;
  }

  template <bool Zero>
  std::uint64_t prefix_after(std::size_t super) const {
    const auto* record = rank_record(super);
    const std::uint64_t ones = record[0] + local_rank(record, 32);
    if constexpr (Zero) {
      return super * kSuperBits + super_size(super) - ones;
    }
    return ones;
  }

  template <bool Zero>
  SelectIndex build_samples() const {
    const std::uint64_t total = Zero ? num_bits_ - one_count_ : one_count_;
    if (total == 0) {
      return {};
    }

    if (num_bits_ <= 1024) {
      // All queries fit one superblock and at most two basic blocks.
      // The first local prefix supplies the complete localization boundary.
      return {};
    }

    const std::size_t sample_count = 1 + (total - 1) / SampleFrequency;
    std::vector<std::uint64_t> samples;
    if (sample_count == std::numeric_limits<std::size_t>::max()) {
      throw std::length_error("Too many select samples");
    }
    samples.reserve(sample_count + 1);
    std::size_t last_occupied = 0;
    [[maybe_unused]] bool homogeneous = true;
    for (std::size_t super = 0; super < super_count(); ++super) {
      const auto end = prefix_after<Zero>(super);
      const auto before = prefix_before<Zero>(super);
      if constexpr (SampleFrequency <= kSuperBits) {
        homogeneous &= end == before || end - before == super_size(super);
      }
      if (end != before) {
        last_occupied = super;
      }
      while (samples.size() < sample_count &&
             samples.size() * SampleFrequency < end) {
        samples.push_back(super);
      }
    }
    samples.push_back(last_occupied + 1);

    SelectIndex result;
    if constexpr (SampleFrequency <= kSuperBits) {
      if (homogeneous) {
        result.samples = PrimarySamples(samples);
        result.shape = SelectIndex::Shape::kHomogeneous;
        return result;
      }
    }
    std::vector<std::uint64_t> occupied;
    std::vector<std::uint64_t> wide_samples;
    std::vector<WideEntry> wide;
    std::vector<typename SelectIndex::Endpoint> endpoints;
    bool needs_occupied_ids = false;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
      if (samples[sample + 1] - samples[sample] < SearchThreshold) {
        continue;
      }
      const std::uint64_t begin_rank = sample * SampleFrequency;
      const std::uint64_t end_rank =
          begin_rank +
          std::min<std::uint64_t>(SampleFrequency, total - begin_rank);
      const std::size_t begin = occupied.size();
      const std::size_t last =
          std::min<std::size_t>(samples[sample + 1], super_count() - 1);
      for (std::size_t super = samples[sample]; super <= last; ++super) {
        const auto before = prefix_before<Zero>(super);
        const auto after = prefix_after<Zero>(super);
        if (after > before && after > begin_rank && before < end_rank) {
          occupied.push_back(super);
          endpoints.push_back(static_cast<typename SelectIndex::Endpoint>(
              std::min<std::uint64_t>(after - begin_rank, SampleFrequency) -
              1));
        }
      }
      const std::size_t count = occupied.size() - begin;
      const bool contiguous = occupied.back() - occupied[begin] + 1 == count;
      needs_occupied_ids |= !contiguous;

      wide.push_back({begin, count | (contiguous ? kContiguous : 0)});
      wide_samples.push_back(sample);
    }
    if (needs_occupied_ids) {
      result.occupied =
          detail::compressed::BlockedMonotoneSamples<MetadataStorage>(occupied);
    }
    result.wide = MetadataArray<WideEntry>(wide);
    result.endpoints = MetadataArray<typename SelectIndex::Endpoint>(endpoints);

    const auto last_wide = wide_samples.empty() ? 0 : wide_samples.back();
    result.shape = result.wide.empty() ? SelectIndex::Shape::kNarrow
                   : last_wide == result.wide.size() - 1
                       ? SelectIndex::Shape::kWidePrefix
                       : SelectIndex::Shape::kMixed;

    if (result.shape != SelectIndex::Shape::kMixed) {
      wide_samples.clear();
    }

    if (one_count_ >= num_bits_ / 16 &&
        num_bits_ - one_count_ >= num_bits_ / 16) {
      std::vector<std::uint64_t> positions;
      positions.reserve(sample_count + 1);
      for (std::size_t sample = 0; sample < sample_count; ++sample) {
        positions.push_back(select_in_super<Zero>(sample * SampleFrequency + 1,
                                                  samples[sample]));
      }
      // Use one past the last selected position, excluding trailing empty bits.
      positions.push_back(select_in_super<Zero>(total, last_occupied) + 1);
      result.positions = PositionSamples(positions, wide_samples);
      result.predicted = true;
    } else {
      result.samples = PrimarySamples(samples, wide_samples);
    }

    return result;
  }

  template <bool Zero>
  std::uint64_t select_direction(std::size_t rank,
                                 const SelectIndex& index,
                                 bool enabled) const {
    if (rank == 0) {
      return 0;
    }
    const std::uint64_t total = Zero ? num_bits_ - one_count_ : one_count_;
    if (!enabled || rank > total) {
      return num_bits_;
    }

    if (num_bits_ <= 1024) {
      return total == num_bits_ ? rank - 1
                                : select_small_index<Zero>(rank, total);
    }

    const std::size_t sample = (rank - 1) / SampleFrequency;
    if constexpr (SampleFrequency <= kSuperBits) {
      if (index.shape == SelectIndex::Shape::kHomogeneous) {
        // The sample interval divides a full superblock's population. Every
        // earlier occupied block contributes exactly kSuperBits occurrences;
        // only the final block can be partial. The sample's first superblock
        // therefore contains the complete query interval, despite empty gaps.
        const auto super = index.samples[sample];

        // Every earlier occupied block contributes exactly kSuperBits, so
        // the local offset needs neither a global nor a local rank load.
        return super * kSuperBits + ((rank - 1) & (kSuperBits - 1));
      }
    }
    std::uint64_t predicted_first = 0;
    if (index.predicted) {
      const auto [first, next] = index.positions.bounds(sample);
      predicted_first = first / kSuperBits;
      if (const auto answer = predict<Zero>(rank, first, next)) {
        return *answer;
      }
    }
    const auto first_super = [&]() {
      return index.predicted ? predicted_first : index.samples[sample];
    };
    std::uint64_t first;
    bool narrow;
    std::optional<std::size_t> mapped_wide;

    if (index.shape == SelectIndex::Shape::kMixed) {
      mapped_wide = index.predicted ? index.positions.wide_entry(sample)
                                    : index.samples.wide_entry(sample);
    }

    if (index.shape != SelectIndex::Shape::kMixed &&
        (index.shape == SelectIndex::Shape::kNarrow ||
         sample < index.wide.size())) {
      narrow = index.shape == SelectIndex::Shape::kNarrow;
      first = narrow || (index.wide[sample].count & kContiguous) != 0
                  ? first_super()
                  : 0;
    } else {
      narrow = !mapped_wide;
      first = narrow || (index.wide[*mapped_wide].count & kContiguous) != 0
                  ? first_super()
                  : 0;
    }
    std::size_t super = first;
    if (narrow) {
      while (prefix_after<Zero>(super) < rank) {
        ++super;
      }
    } else {
      // Wide-prefix shapes index directly; mixed shapes use the exact bitmap
      // map.
      const auto* entry = &index.wide[mapped_wide ? *mapped_wide : sample];

      std::size_t low = entry->begin;
      const std::size_t count = entry->count & ~kContiguous;
      const auto relative_rank = (rank - 1) % SampleFrequency;

      // A single occupied superblock already proves the answer; its stored
      // endpoint contributes no localization information to this query.
      if (count != 1) {
        low += detail::compressed::grouped_lower_bound(
            index.endpoints.data() + low, count, relative_rank);
      }

      assert(low < entry->begin + count);
      super = (entry->count & kContiguous) != 0 ? super + low - entry->begin
                                                : index.occupied[low];
    }

    return select_in_super<Zero>(rank, super);
  }

  template <bool Zero>
  std::uint64_t select_in_super(std::size_t rank, std::size_t super) const {
    const auto* record = rank_record(super);
    const std::uint64_t target = rank - prefix_before<Zero>(super);
    const std::size_t valid = super_size(super);
    const auto ones_local = local_rank(record, 32);
    const std::uint64_t total_local = Zero ? valid - ones_local : ones_local;
    if (total_local == valid) {
      return super * kSuperBits + target - 1;
    }
    const std::size_t blocks = 1 + (valid - 1) / kBasicBits;
    const auto packed_location =
        detail::compressed::packed_prefix_search<Zero>(record, target, valid);

    std::size_t basic;
    std::uint64_t before;
    std::uint64_t after;
    if (packed_location) {
      basic = packed_location->basic;
      before = packed_location->before;
      after = packed_location->after;
    } else {
      basic = (target - 1) * blocks / total_local;
      const auto count_at = [&](std::size_t boundary) {
        const std::uint64_t ones = local_rank(record, boundary);
        return Zero ? std::min(boundary * kBasicBits, valid) - ones : ones;
      };
      // The estimate is only a starting point. Exact corrections also handle
      // arbitrary clustering and repeated local counts within this cache line.
      before = count_at(basic);
      after = before;
      if (before >= target) {
        do {
          after = before;
          before = count_at(--basic);
        } while (before >= target);
      } else {
        after = count_at(basic + 1);
        while (after < target) {
          before = after;
          after = count_at(++basic + 1);
        }
      }
    }
    const std::size_t base = super * kSuperBits + basic * kBasicBits;
    return select_basic<Zero>(base, target - before, after - before);
  }

  template <bool Zero>
  std::uint64_t basic_prefix(std::size_t end) const {
    if constexpr (Zero) {
      return rank0_impl(end);
    } else {
      return rank_impl(end);
    }
  }

  template <bool Zero>
  std::optional<std::uint64_t> predict(std::size_t rank,
                                       std::uint64_t first,
                                       std::uint64_t next) const {
    const auto span = next - first;
    if (span > 65536) {
      return std::nullopt;
    }
    const auto guess =
        first + static_cast<std::uint64_t>(static_cast<unsigned __int128>(
                                               (rank - 1) % SampleFrequency) *
                                           span / SampleFrequency);
    const std::size_t base = (guess / kBasicBits) * kBasicBits;
    if (num_bits_ - base < kBasicBits) {
      return std::nullopt;
    }
    __builtin_prefetch(bits_.data() + base / 64, 0, 3);
    const auto before = basic_prefix<Zero>(base);
    const auto after = basic_prefix<Zero>(base + kBasicBits);
    if (before < rank && rank <= after) {
      return select_basic<Zero>(base, rank - before, after - before);
    }
    if (rank <= before && base >= kBasicBits) {
      const auto lower = basic_prefix<Zero>(base - kBasicBits);
      if (lower < rank) {
        return select_basic<Zero>(base - kBasicBits, rank - lower,
                                  before - lower);
      }
    } else if (rank > after && num_bits_ - base >= 2 * kBasicBits) {
      const auto upper = basic_prefix<Zero>(base + 2 * kBasicBits);
      if (rank <= upper) {
        return select_basic<Zero>(base + kBasicBits, rank - after,
                                  upper - after);
      }
    }
    return std::nullopt;
  }

  template <bool Zero>
  [[gnu::noinline]] std::uint64_t select_small_index(
      std::size_t rank,
      std::uint64_t total) const {
    const std::uint64_t first_ones = local_rank(rank_record(0), 1);
    const std::uint64_t first_end =
        Zero ? std::min(kBasicBits, num_bits_) - first_ones : first_ones;
    const std::size_t basic = rank > first_end;
    const std::uint64_t before = basic ? first_end : 0;
    const std::uint64_t after = basic ? total : first_end;
    return select_basic<Zero>(basic * kBasicBits, rank - before,
                              after - before);
  }

  template <bool Zero>
  [[gnu::always_inline]] inline std::uint64_t select_basic(
      std::size_t base,
      std::uint64_t remaining,
      std::uint64_t population) const {
    if (population == kBasicBits) {
      return base + remaining - 1;
    }
    if (num_bits_ - base >= kBasicBits) {
      return base + detail::compressed::select_loaded_block<Zero>(
                        bits_.data() + base / 64, remaining - 1);
    }
    for (std::size_t lane = 0; lane < 8; ++lane) {
      const std::size_t position = base + lane * 64;
      if (position >= num_bits_) {
        break;
      }
      std::uint64_t word = bits_[position / 64];
      if constexpr (Zero) {
        word = ~word;
      }
      word &=
          ::first_bits_mask(std::min<std::size_t>(64, num_bits_ - position));
      const unsigned count = std::popcount(word);
      if (remaining <= count) {
        return position + ::select_64(word, remaining - 1);
      }
      remaining -= count;
    }
    assert(false && "Select localization must contain the target");
    return num_bits_;
  }
};

}  // namespace pixie
