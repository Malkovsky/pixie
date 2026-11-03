#pragma once

/** @file compressed.h
 * @brief Compressed rank/select support with metadata-validated source
 * prefetch.
 */

#include <pixie/detail/serialization_frame.h>
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
 * @brief Non-owning rank/select index with compressed interleaved metadata.
 * @details Source words must remain alive and immutable. Copies own independent
 * metadata and borrow the same source; moves leave an empty, disabled index.
 * Half-open rank, one-based select and sentinels follow RankSelectBase; rank
 * boundaries beyond size() clamp to the total. Serialization retains metadata
 * only; the source encoding is supplied independently when loading.
 *
 * @code{.text}
 * Rank: one 64-byte record per 16,384 source bits (3.125% overhead)
 *   +------------------+-------------------------------------------+
 *   | global ones: 64b | 32 local prefix ends: 14b each (448b)       |
 *   +------------------+-------------------------------------------+
 *   local prefix ends correspond to 512-bit source-block boundaries
 *
 * Select: independent samples every SampleFrequency ones / zeros
 *   coarse EF groups: descriptor -> [base | unary highs | low bits]
 *   dense EF groups:  [base | unary highs | wide map* | 32 low fields]
 *                      <---------- fixed stride ------------------>
 *   wide map*: [directory prefix | bitmap], present for mixed shapes
 *   coarse groups keep the wide map separately
 *
 *   wide directory: [begin, count] -> rank endpoints + occupied IDs*
 *   occupied IDs* are omitted when all wide intervals are contiguous
 * @endcode
 *
 * Short intervals scan global ranks; wide intervals search stored endpoints.
 * Dense samples store exact positions: decode both neighbors, prefetch an
 * estimated source block, then validate against rank. A miss corrects by one
 * block or falls back to the directory. Homogeneous regions use arithmetic.
 * AlignedStorage aligns rank records to 64 bytes; custom backends need only
 * uint64_t alignment. All retained arrays use MetadataStorage.
 *
 * @tparam MetadataStorage Default-empty, owning contiguous storage with deep
 * copies and nonthrowing move construction/assignment. Construction requires
 * resize(bit_count), writable_bytes() and as_bytes(), exposing exactly the
 * requested bytes with uint64_t alignment. Memory reporting requires
 * allocated_bytes(). The borrowed source needs no cache-line alignment.
 * @tparam SampleFrequency Positive power-of-two occurrence sampling interval.
 * @tparam SearchThreshold Positive superblock-span threshold for wide indexing.
 */
template <StorageImplementation MetadataStorage = AlignedStorage,
          std::size_t SampleFrequency = 4096,
          std::size_t SearchThreshold = 16>
class RankSelectSupportCompressed
    : public RankSelectBase<RankSelectSupportCompressed<MetadataStorage,
                                                        SampleFrequency,
                                                        SearchThreshold>>,
      public SerializationBase<RankSelectSupportCompressed<MetadataStorage,
                                                           SampleFrequency,
                                                           SearchThreshold>> {
  static_assert(std::has_single_bit(SampleFrequency));
  static_assert(SearchThreshold > 0);
  static_assert(std::is_nothrow_move_constructible_v<MetadataStorage> &&
                    std::is_nothrow_move_assignable_v<MetadataStorage>,
                "Compressed metadata storage must have nonthrowing moves");
  static constexpr bool kNothrowDefaultConstruction =
      std::is_nothrow_default_constructible_v<MetadataStorage>;

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

  /** @brief Write versioned, little-endian compressed metadata, excluding bits.
   * @details Template sampling parameters are recorded and must match on load.
   */
  void serialize_impl(BinaryWriter& writer) const {
    const auto begin = detail::begin_frame(writer, kSerializationTag);
    writer.write_size(SampleFrequency);
    writer.write_size(SearchThreshold);
    writer.write_size(num_bits_);
    writer.write_u64(one_count_);
    writer.write_u64(static_cast<uint8_t>(select_support_));
    ranks_.serialize(writer);
    select1_.serialize(writer);
    select0_.serialize(writer);
    detail::end_frame(writer, begin);
  }

  /** @brief Restore owning metadata over an independently supplied encoding.
   * @param reader Cursor advanced by one artifact only on success.
   * @param words Immutable source words, which must outlive the index.
   * @param validation Quick checks metadata without scanning the source; full
   * validation rebuilds and compares metadata against the supplied source.
   * @details The caller must pair quick-loaded metadata with its original
   * encoding. No hash is computed. Metadata is copied, not mapped.
   * @throws std::invalid_argument for incompatible, truncated or inconsistent
   * metadata or a short source. Allocation and size exceptions propagate.
   */
  static RankSelectSupportCompressed deserialize_impl(
      BinaryReader& reader,
      std::span<const uint64_t> words,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires detail::compressed::MutableMetadataStorage<MetadataStorage>
  {
    auto candidate = reader;
    auto payload = detail::read_frame(candidate, kSerializationTag);
    if (payload.read_size() != SampleFrequency ||
        payload.read_size() != SearchThreshold) {
      throw std::invalid_argument(
          "Incompatible compressed rank/select configuration");
    }
    RankSelectSupportCompressed result;
    result.bits_ = words;
    result.num_bits_ = payload.read_size();
    result.one_count_ = payload.read_u64();
    const auto support = payload.read_u64();
    if (support > 3 || result.one_count_ > result.num_bits_ ||
        result.num_bits_ / 64 + (result.num_bits_ % 64 != 0) > words.size()) {
      throw std::invalid_argument("Invalid compressed rank/select dimensions");
    }
    result.select_support_ = static_cast<SelectSupport>(support);
    result.ranks_ = decltype(ranks_)::deserialize(payload);
    result.select1_ = SelectIndex::deserialize(payload);
    result.select0_ = SelectIndex::deserialize(payload);
    payload.require_zero_padding(7);
    result.validate_metadata();
    if (validation == DeserializationValidation::kFull) {
      RankSelectSupportCompressed expected(words, result.num_bits_,
                                           result.select_support_);
      if (result.one_count_ != expected.one_count_ ||
          result.ranks_ != expected.ranks_ ||
          result.select1_ != expected.select1_ ||
          result.select0_ != expected.select0_) {
        throw std::invalid_argument(
            "Compressed rank/select metadata does not match source");
      }
    }
    reader = candidate;
    return result;
  }

  /** @brief Duplicate metadata while borrowing the same immutable source. */
  RankSelectSupportCompressed(const RankSelectSupportCompressed&) = default;

  /** @brief Transfer metadata, leaving the source empty with select disabled.
   */
  RankSelectSupportCompressed(RankSelectSupportCompressed&& other) noexcept(
      kNothrowDefaultConstruction) {
    swap_state(other);
  }

  /**
   * @brief Duplicate metadata while borrowing the assigned source.
   * @details Copy failure leaves this index unchanged.
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
  RankSelectSupportCompressed&
  operator=(RankSelectSupportCompressed&& other) noexcept(
      kNothrowDefaultConstruction) {
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
  // Inlining exposes block-aligned prefix calls in select's prediction path.
#if __has_cpp_attribute(gnu::always_inline)
  [[gnu::always_inline]]
#endif
  inline std::uint64_t rank_impl(std::size_t end) const {
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
  static constexpr uint64_t kSerializationTag = 0x3153435258495050ULL;
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
    bool operator==(const WideEntry&) const = default;
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

    bool operator==(const SelectIndex&) const = default;

    void serialize(BinaryWriter& writer) const {
      writer.write_u64(predicted);
      writer.write_u64(static_cast<uint8_t>(shape));
      samples.serialize(writer);
      positions.serialize(writer);
      occupied.serialize(writer);
      wide.serialize(writer, [](BinaryWriter& out, WideEntry entry) {
        out.write_size(entry.begin);
        out.write_size(entry.count);
      });
      endpoints.serialize(writer);
    }

    static SelectIndex deserialize(BinaryReader& reader) {
      SelectIndex result;
      const auto predicted = reader.read_u64();
      const auto shape = reader.read_u64();
      if (predicted > 1 || shape > 3) {
        throw std::invalid_argument("Invalid compressed select shape");
      }
      result.predicted = predicted;
      result.shape = static_cast<Shape>(shape);
      result.samples = PrimarySamples::deserialize(reader);
      result.positions = PositionSamples::deserialize(reader);
      result.occupied = decltype(occupied)::deserialize(reader);
      result.wide =
          decltype(wide)::deserialize(reader, 16, [](BinaryReader& in) {
            return WideEntry{in.read_size(), in.read_size()};
          });
      result.endpoints = decltype(endpoints)::deserialize(reader);
      return result;
    }

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

  void swap_state(RankSelectSupportCompressed& other) noexcept {
    using std::swap;
    swap(bits_, other.bits_);
    swap(num_bits_, other.num_bits_);
    swap(one_count_, other.one_count_);
    swap(ranks_, other.ranks_);
    swap(select1_, other.select1_);
    swap(select0_, other.select0_);
    swap(select_support_, other.select_support_);
  }

  void validate_metadata() const {
    const auto supers = num_bits_ / kSuperBits + (num_bits_ % kSuperBits != 0);
    if (ranks_.size() != supers * 8) {
      throw std::invalid_argument("Invalid rank record count");
    }
    uint64_t global = 0;
    for (size_t s = 0; s < supers; ++s) {
      const auto* record = rank_record(s);
      if (record[0] != global) {
        throw std::invalid_argument("Invalid global rank prefix");
      }
      uint64_t previous = 0;
      size_t previous_bits = 0;
      for (size_t b = 1; b <= 32; ++b) {
        const auto end = std::min(b * kBasicBits, super_size(s));
        const auto rank = local_rank(record, b);
        if (rank < previous || rank - previous > end - previous_bits) {
          throw std::invalid_argument("Invalid local rank prefix");
        }
        previous = rank;
        previous_bits = end;
      }
      global += previous;
    }
    if (global != one_count_) {
      throw std::invalid_argument("Invalid total rank");
    }
    validate_select<false>(select1_, supports_select1_impl());
    validate_select<true>(select0_, supports_select0_impl());
  }

  template <bool Zero>
  void validate_select(const SelectIndex& index, bool enabled) const {
    const uint64_t total = Zero ? num_bits_ - one_count_ : one_count_;
    if (!enabled || !total || num_bits_ <= 1024) {
      if (index != SelectIndex{}) {
        throw std::invalid_argument("Unexpected select metadata");
      }
      return;
    }
    const size_t count = 1 + (total - 1) / SampleFrequency;
    if (count == SIZE_MAX) {
      throw std::length_error("Too many select samples");
    }
    const bool mixed = index.shape == SelectIndex::Shape::kMixed;
    const auto mapped = mixed ? index.wide.size() : 0;
    index.samples.validate(index.predicted ? 0 : count + 1, super_count(),
                           index.predicted ? 0 : mapped);
    index.positions.validate(index.predicted ? count + 1 : 0, num_bits_,
                             index.predicted ? mapped : 0);
    if (index.wide.size() > count ||
        ((index.shape == SelectIndex::Shape::kNarrow ||
          index.shape == SelectIndex::Shape::kHomogeneous) &&
         !index.wide.empty())) {
      throw std::invalid_argument("Invalid wide select directory");
    }
    size_t endpoint_count = 0;
    bool needs_occupied = false;
    for (size_t i = 0; i < index.wide.size(); ++i) {
      const auto entry = index.wide[i];
      const auto n = entry.count & ~kContiguous;
      if (!n || entry.begin != endpoint_count ||
          n > index.endpoints.size() - endpoint_count) {
        throw std::invalid_argument("Invalid wide select interval");
      }
      endpoint_count += n;
      needs_occupied |= !(entry.count & kContiguous);
    }
    if (endpoint_count != index.endpoints.size()) {
      throw std::invalid_argument("Invalid select endpoint count");
    }
    index.occupied.validate(needs_occupied ? endpoint_count : 0,
                            super_count() - 1);
    for (size_t i = 0; i < count; ++i) {
      const uint64_t begin_rank = i * SampleFrequency;
      const uint64_t end_rank =
          begin_rank + std::min<uint64_t>(SampleFrequency, total - begin_rank);
      const auto first =
          index.predicted ? index.positions[i] / kSuperBits : index.samples[i];
      if (first >= super_count() || prefix_before<Zero>(first) > begin_rank ||
          prefix_after<Zero>(first) <= begin_rank) {
        throw std::invalid_argument("Invalid select sample localization");
      }
      if (index.shape == SelectIndex::Shape::kHomogeneous) {
        if (index.predicted || SampleFrequency > kSuperBits ||
            prefix_after<Zero>(first) - prefix_before<Zero>(first) !=
                super_size(first)) {
          throw std::invalid_argument("Invalid homogeneous select interval");
        }
      }
      std::optional<size_t> w;
      if (mixed) {
        w = index.predicted ? index.positions.wide_entry(i)
                            : index.samples.wide_entry(i);
      } else if (index.shape == SelectIndex::Shape::kWidePrefix &&
                 i < index.wide.size()) {
        w = i;
      }
      if (!w) {
        continue;
      }
      if (*w >= index.wide.size()) {
        throw std::invalid_argument("Invalid wide select reference");
      }
      const auto entry = index.wide[*w];
      const auto n = entry.count & ~kContiguous;
      uint64_t previous_rank = begin_rank;
      size_t previous_super = first;
      for (size_t j = 0; j < n; ++j) {
        const auto s = (entry.count & kContiguous)
                           ? first + j
                           : index.occupied[entry.begin + j];
        if (s >= super_count() || (j && s <= previous_super) ||
            prefix_before<Zero>(s) != previous_rank ||
            prefix_after<Zero>(s) <= previous_rank) {
          // The first superblock may begin before this sampled interval.
          if (j || s != first || prefix_before<Zero>(s) > begin_rank ||
              prefix_after<Zero>(s) <= begin_rank) {
            throw std::invalid_argument("Invalid occupied select interval");
          }
        }
        previous_rank = std::min(end_rank, prefix_after<Zero>(s));
        if (index.endpoints[entry.begin + j] !=
            previous_rank - begin_rank - 1) {
          throw std::invalid_argument("Invalid select endpoint");
        }
        previous_super = s;
      }
      if (previous_rank != end_rank) {
        throw std::invalid_argument("Incomplete select interval");
      }
    }
  }

  std::size_t super_count() const { return ranks_.size() / (kRecordBytes / 8); }

  const std::uint64_t* rank_record(std::size_t super) const {
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
        // Sampling divides each full occupied superblock's population, so
        // the local offset follows arithmetically without reading rank
        // metadata.
        const auto super = index.samples[sample];
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
    const WideEntry* entry = nullptr;
    switch (index.shape) {
      case SelectIndex::Shape::kNarrow:
        break;
      case SelectIndex::Shape::kWidePrefix:
        if (sample < index.wide.size()) {
          entry = &index.wide[sample];
        }
        break;
      case SelectIndex::Shape::kMixed:
        if (const auto mapped = index.predicted
                                    ? index.positions.wide_entry(sample)
                                    : index.samples.wide_entry(sample)) {
          entry = &index.wide[*mapped];
        }
        break;
      case SelectIndex::Shape::kHomogeneous:
        // The arithmetic path above handles all homogeneous indexes.
        assert(false);
        break;
    }

    std::size_t super;
    if (entry == nullptr) {
      super = index.predicted ? predicted_first : index.samples[sample];
      while (prefix_after<Zero>(super) < rank) {
        ++super;
      }
    } else {
      const std::size_t count = entry->count & ~kContiguous;
      const auto relative_rank = (rank - 1) % SampleFrequency;
      // A singleton interval identifies the superblock without an endpoint
      // load.
      const std::size_t offset =
          count == 1 ? 0
                     : detail::compressed::grouped_lower_bound(
                           index.endpoints.data() + entry->begin, count,
                           relative_rank);
      assert(offset < count);
      if ((entry->count & kContiguous) != 0) {
        super = (index.predicted ? predicted_first : index.samples[sample]) +
                offset;
      } else {
        super = index.occupied[entry->begin + offset];
      }
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
#if __has_cpp_attribute(gnu::noinline)
  [[gnu::noinline]]
#endif
  std::uint64_t select_small_index(std::size_t rank,
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
#if __has_cpp_attribute(gnu::always_inline)
  [[gnu::always_inline]]
#endif
  inline std::uint64_t select_basic(std::size_t base,
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
