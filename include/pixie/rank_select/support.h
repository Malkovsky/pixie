#pragma once

#include <pixie/bit_stream.h>
#include <pixie/bits.h>
#include <pixie/rank_select.h>
#include <pixie/storage/aligned.h>
#include <pixie/storage/read_only_view.h>

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>

#ifdef PIXIE_DIAGNOSTICS
#include <spdlog/spdlog.h>
#endif

namespace pixie {

/**
 * @brief Rank/select support over an external packed bit sequence.
 * @details This object does not own its source bits. It owns rank/select
 * metadata through @p MetadataStorage and keeps a read-only storage view plus
 * a validated 64-bit-word view of the caller-owned source.
 *
 * Structure overview:
 * - Super blocks of 2^16 bits with 64-bit ranks (0.098% overhead).
 * - Basic blocks of 512 bits with 16-bit ranks (3.125% overhead).
 * - Optional 64-bit select samples every 2^14 occurrences. One enabled
 * direction uses at most 0.391% overhead; when both directions are enabled,
 * their combined asymptotic overhead is also 0.391%.
 *
 * The rank metadata plus both select directions therefore use 3.613%
 * asymptotic metadata overhead. These ratios exclude the object itself,
 * sentinel entries, and 64-byte allocation rounding, which are included by
 * memory_usage_bytes_impl().
 *
 * Rank: 2 table lookups plus SIMD popcount in the 512-bit block.
 *
 * Select:
 * - Start from a sampled super block.
 * - SIMD linear scan to find the super block.
 * - SIMD linear scan to find the basic block.
 *
 * This variant does not interleave data and index, favoring simpler scans.
 * Rank metadata and enabled select samples are built in one scan over the
 * source words.
 */
template <StorageImplementation MetadataStorage = AlignedStorage>
class RankSelectSupport
    : public RankSelectBase<RankSelectSupport<MetadataStorage>> {
 public:
  /**
   * @brief Select directions to index during construction.
   */
  enum class SelectSupport : uint8_t {
    kNone = 0,
    kSelect1 = 1,
    kSelect0 = 2,
    kBoth = 3,
  };

 private:
  constexpr static size_t kWordSize = 64;
  constexpr static size_t kSuperBlockRankIntSize = 64;
  constexpr static size_t kBasicBlockRankIntSize = 16;
  constexpr static size_t kBasicBlockSize = 512;
  constexpr static size_t kWordsPerBlock = 8;
  constexpr static size_t kSuperBlockSize = 65536;
  constexpr static size_t kBlocksPerSuperBlock = 128;
  constexpr static size_t kSelectSampleFrequency = 16384;

  alignas(64) uint64_t delta_super[8]{};
  alignas(64) uint16_t delta_basic[32]{};

  MetadataStorage super_block_rank_;  // 64-bit global prefix sums
  MetadataStorage basic_block_rank_;  // 16-bit local prefix sums
  MetadataStorage select_samples_;    // 64-bit global positions
  ReadOnlyStorageView source_storage_;
  std::span<const uint64_t> bits_;
  size_t num_bits_{};
  size_t padded_size_{};
  size_t max_rank_{};
  size_t select1_sample_begin_{};
  size_t select1_sample_count_{};
  size_t select0_sample_begin_{};
  size_t select0_sample_count_{};
  SelectSupport select_support_ = SelectSupport::kNone;
  bool select0_samples_reversed_ = false;

  static bool builds_select1(SelectSupport support) {
    return (static_cast<uint8_t>(support) &
            static_cast<uint8_t>(SelectSupport::kSelect1)) != 0;
  }

  static bool builds_select0(SelectSupport support) {
    return (static_cast<uint8_t>(support) &
            static_cast<uint8_t>(SelectSupport::kSelect0)) != 0;
  }

  size_t logical_word_count() const {
    return (num_bits_ + kWordSize - 1) / kWordSize;
  }

  size_t logical_word_bits(size_t word_index) const {
    const size_t begin = word_index * kWordSize;
    if (begin >= num_bits_) {
      return 0;
    }
    return std::min(kWordSize, num_bits_ - begin);
  }

  uint64_t logical_word(size_t word_index) const {
    if (word_index >= bits_.size()) {
      return 0;
    }
    const size_t bits = logical_word_bits(word_index);
    if (bits == 0) {
      return 0;
    }
    if (bits == kWordSize) {
      return bits_[word_index];
    }
    return bits_[word_index] & first_bits_mask(bits);
  }

  uint64_t rank_in_basic_block(size_t basic_block, size_t offset) const {
    if (offset == 0) {
      return 0;
    }
    const size_t first_word = basic_block * kWordsPerBlock;
    if (first_word + kWordsPerBlock <= bits_.size()) {
      return rank_512(&bits_[first_word], offset);
    }

    uint64_t result = 0;
    size_t word_index = first_word;
    while (offset >= kWordSize) {
      result += std::popcount(logical_word(word_index));
      offset -= kWordSize;
      ++word_index;
    }
    if (offset != 0) {
      result +=
          std::popcount(logical_word(word_index) & first_bits_mask(offset));
    }
    return result;
  }

  uint64_t select_in_words(size_t first_word, size_t rank, bool value) const {
    const size_t first_bit = first_word * kWordSize;
    if (first_bit + kBasicBlockSize <= num_bits_ &&
        first_word + kWordsPerBlock <= bits_.size()) {
      return value ? first_bit + select_512(&bits_[first_word], rank - 1)
                   : first_bit + select0_512(&bits_[first_word], rank - 1);
    }

    for (size_t word_index = first_word; word_index < logical_word_count();
         ++word_index) {
      const uint64_t word = logical_word(word_index);
      const uint64_t candidates =
          value ? word
                : (~word & first_bits_mask(logical_word_bits(word_index)));
      const size_t count = std::popcount(candidates);
      if (rank > count) {
        rank -= count;
        continue;
      }
      return word_index * kWordSize + select_64(candidates, rank - 1);
    }
    return num_bits_;
  }

  static size_t select_sample_count_for_rank(size_t rank_count) {
    return 1 + rank_count / kSelectSampleFrequency;
  }

  static size_t select_sample_upper_bound(size_t bit_count) {
    return select_sample_count_for_rank(bit_count);
  }

  struct SelectSampleWriter {
    std::span<uint64_t> words;
    size_t next = 0;
    size_t count = 0;
    size_t capacity = 0;
    bool enabled = false;
    bool reversed = false;

    SelectSampleWriter() = default;

    SelectSampleWriter(std::span<uint64_t> words,
                       size_t begin,
                       size_t capacity,
                       bool enabled,
                       bool reversed)
        : words(words),
          next(begin),
          capacity(capacity),
          enabled(enabled),
          reversed(reversed) {}

    void append(uint64_t sample) {
      if (!enabled) {
        return;
      }
      if (count >= capacity) [[unlikely]] {
        throw std::invalid_argument(
            "RankSelectSupport one_count hint is inconsistent with input bits");
      }
      words[next] = sample;
      ++count;
      if (reversed) {
        if (next != 0) {
          --next;
        }
      } else {
        ++next;
      }
    }
  };

  struct SelectSampleWriters {
    SelectSampleWriter ones;
    SelectSampleWriter zeros;
    bool shrink_after_build = false;
  };

  SelectSampleWriters initialize_select_sample_writers(
      bool need_select1,
      bool need_select0,
      std::optional<size_t> one_count) {
    select1_sample_begin_ = 0;
    select1_sample_count_ = 0;
    select0_sample_begin_ = 0;
    select0_sample_count_ = 0;
    select0_samples_reversed_ = false;
    select_samples_.resize(0);

    SelectSampleWriters writers;
    if (!need_select1 && !need_select0) {
      return writers;
    }

    const std::optional<size_t> zero_count =
        one_count ? std::optional<size_t>(num_bits_ - *one_count)
                  : std::nullopt;
    if (need_select1 && need_select0) {
      const size_t one_sample_capacity =
          one_count ? select_sample_count_for_rank(*one_count)
                    : 2 + num_bits_ / kSelectSampleFrequency;
      const size_t zero_sample_capacity =
          zero_count ? select_sample_count_for_rank(*zero_count)
                     : 2 + num_bits_ / kSelectSampleFrequency;
      const size_t total_samples =
          one_count ? one_sample_capacity + zero_sample_capacity
                    : 2 + num_bits_ / kSelectSampleFrequency;
      select_samples_.resize(total_samples * kWordSize);
      auto samples = select_samples_.writable_words64();
      select1_sample_begin_ = 0;
      select0_samples_reversed_ = true;
      writers.ones =
          SelectSampleWriter(samples, 0, one_sample_capacity, true, false);
      writers.zeros = SelectSampleWriter(samples, total_samples - 1,
                                         zero_sample_capacity, true, true);
      writers.ones.append(0);
      writers.zeros.append(0);
      return writers;
    }

    const size_t sample_capacity =
        need_select1 ? (one_count ? select_sample_count_for_rank(*one_count)
                                  : select_sample_upper_bound(num_bits_))
                     : (zero_count ? select_sample_count_for_rank(*zero_count)
                                   : select_sample_upper_bound(num_bits_));
    select_samples_.resize(sample_capacity * kWordSize);
    auto samples = select_samples_.writable_words64();
    writers.shrink_after_build = !one_count;
    if (need_select1) {
      select1_sample_begin_ = 0;
      writers.ones =
          SelectSampleWriter(samples, 0, sample_capacity, true, false);
      writers.ones.append(0);
    } else {
      select0_sample_begin_ = 0;
      writers.zeros =
          SelectSampleWriter(samples, 0, sample_capacity, true, false);
      writers.zeros.append(0);
    }
    return writers;
  }

  void finalize_select_sample_writers(SelectSampleWriters writers) {
    select1_sample_count_ = writers.ones.count;
    select0_sample_count_ = writers.zeros.count;
    if (writers.zeros.reversed) {
      select0_sample_begin_ = writers.zeros.next + 1;
      auto zero_samples = writers.zeros.words.subspan(select0_sample_begin_,
                                                      select0_sample_count_);
      std::reverse(zero_samples.begin(), zero_samples.end());
      select0_samples_reversed_ = false;
    }
    if (writers.shrink_after_build) {
      const size_t sample_count = select1_sample_count_ != 0
                                      ? select1_sample_count_
                                      : select0_sample_count_;
      select_samples_.resize(sample_count * kWordSize);
      select_samples_.shrink_to_fit();
    }
  }

  uint64_t select1_sample(size_t sample_index) const {
    auto samples = select_samples_.as_words64();
    return samples[select1_sample_begin_ + sample_index];
  }

  uint64_t select0_sample(size_t sample_index) const {
    auto samples = select_samples_.as_words64();
    if (select0_samples_reversed_) {
      return samples[select0_sample_begin_ + select0_sample_count_ - 1 -
                     sample_index];
    }
    return samples[select0_sample_begin_ + sample_index];
  }

  /**
   * @brief Precompute rank and requested select samples in one word scan.
   */
  void build_rank_select(SelectSupport support,
                         std::optional<size_t> one_count) {
    select_support_ = support;
    size_t num_superblocks =
        8 + (padded_size_ == 0 ? 0 : (padded_size_ - 1) / kSuperBlockSize);
    // Add more blocks to ease SIMD processing
    // num_basicblocks to fully cover superblock, i.e. 128
    // This reduces branching in select
    num_superblocks = ((num_superblocks + 7) / 8) * 8;
    size_t num_basicblocks = num_superblocks * kBlocksPerSuperBlock;
    super_block_rank_.resize(num_superblocks * 64);
    basic_block_rank_.resize(num_basicblocks * 16);

    auto super_block_rank = super_block_rank_.writable_words64();
    auto basic_block_rank = basic_block_rank_.writable_words16();

    const bool need_select1 = builds_select1(support);
    const bool need_select0 = builds_select0(support);
    if (one_count && *one_count > num_bits_) {
      throw std::invalid_argument(
          "RankSelectSupport one_count hint cannot exceed num_bits");
    }
    auto select_writers =
        initialize_select_sample_writers(need_select1, need_select0, one_count);

    uint64_t super_block_sum = 0;
    uint64_t basic_block_sum = 0;
    uint64_t milestone = kSelectSampleFrequency;
    uint64_t milestone0 = kSelectSampleFrequency;
    uint64_t rank = 0;
    uint64_t rank0 = 0;

    for (size_t i = 0; i / kBasicBlockSize < basic_block_rank.size();
         i += kWordSize) {
      if (i % kSuperBlockSize == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank[i / kSuperBlockSize] = super_block_sum;
        basic_block_sum = 0;
      }
      if (i % kBasicBlockSize == 0) {
        basic_block_rank[i / kBasicBlockSize] =
            static_cast<uint16_t>(basic_block_sum);
      }
      if (i / kWordSize < logical_word_count()) {
        const size_t word_index = i / kWordSize;
        const uint64_t word = logical_word(word_index);
        const size_t word_bits = logical_word_bits(word_index);
        const uint64_t ones = std::popcount(word);
        const uint64_t zeros = word_bits - ones;
        if (need_select1 && rank + ones >= milestone) {
          const auto pos = select_64(word, milestone - rank - 1);
          // TODO: try including global rank into select samples to save
          //       a cache miss on global rank scan
          select_writers.ones.append((64 * word_index + pos) / kSuperBlockSize);
          milestone += kSelectSampleFrequency;
        }
        if (need_select0 && rank0 + zeros >= milestone0) {
          const uint64_t zero_word = ~word & first_bits_mask(word_bits);
          const auto pos = select_64(zero_word, milestone0 - rank0 - 1);
          select_writers.zeros.append((64 * word_index + pos) /
                                      kSuperBlockSize);
          milestone0 += kSelectSampleFrequency;
        }
        basic_block_sum += ones;
        rank += ones;
        rank0 += zeros;
      }
    }
    max_rank_ = super_block_sum + basic_block_sum;
    finalize_select_sample_writers(select_writers);

    for (size_t i = 0; i < 8; ++i) {
      delta_super[i] = i * kSuperBlockSize;
    }
    for (size_t i = 0; i < 32; ++i) {
      delta_basic[i] = i * kBasicBlockSize;
    }
  }

  /**
   * @brief First step of the select operation.
   * @param rank 1-based
   * rank of the 1-bit to locate.
   */
  uint64_t find_superblock(uint64_t rank) const {
    auto super_block_rank = super_block_rank_.as_words64();

    uint64_t left = select1_sample(rank / kSelectSampleFrequency);

    while (left + 7 < super_block_rank.size()) {
      auto len = lower_bound_8x64(&super_block_rank[left], rank);
      if (len < 8) {
        return left + len - 1;
      }
      left += 8;
    }
    if (left + 3 < super_block_rank.size()) {
      auto len = lower_bound_4x64(&super_block_rank[left], rank);
      if (len < 4) {
        return left + len - 1;
      }
      left += 4;
    }
    while (left < super_block_rank.size() && super_block_rank[left] < rank) {
      left++;
    }
    return left - 1;
  }

  /**
   * @brief First step of the select0 operation.
   * @param rank0 1-based
   * rank of the 0-bit to locate.
   */
  uint64_t find_superblock_zeros(uint64_t rank0) const {
    auto super_block_rank = super_block_rank_.as_words64();

    uint64_t left = select0_sample(rank0 / kSelectSampleFrequency);

    while (left + 7 < super_block_rank.size()) {
      auto len = lower_bound_delta_8x64(&super_block_rank[left], rank0,
                                        delta_super, kSuperBlockSize * left);
      if (len < 8) {
        return left + len - 1;
      }
      left += 8;
    }
    if (left + 3 < super_block_rank.size()) {
      auto len = lower_bound_delta_4x64(&super_block_rank[left], rank0,
                                        delta_super, kSuperBlockSize * left);
      if (len < 4) {
        return left + len - 1;
      }
      left += 4;
    }
    while (left < super_block_rank.size() &&
           kSuperBlockSize * left - super_block_rank[left] < rank0) {
      left++;
    }
    return left - 1;
  }

  /**
   * @brief SIMD-optimized linear scan.
   * @param local_rank Rank within
   * the super block.
   * @param s_block Super block index.
   * @details
   *
   * Processes 32 16-bit entries at once (full cache line), so there is at most

   * * 4 iterations.
   */
  uint64_t find_basicblock(uint16_t local_rank, uint64_t s_block) const {
    auto basic_block_rank = basic_block_rank_.as_words16();

    for (size_t pos = 0; pos < kBlocksPerSuperBlock; pos += 32) {
      auto count = lower_bound_32x16(
          &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank);
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
    }
    return kBlocksPerSuperBlock * s_block + kBlocksPerSuperBlock - 1;
  }

  /**
   * @brief SIMD-optimized linear scan.
   * @param local_rank0 Rank of
   * zeros within the super block.
   * @param s_block Super block index.
   *
   * @details
   * Processes 32 16-bit entries at once (full cache line), so
   * there is at most
   * 4 iterations.
   */
  uint64_t find_basicblock_zeros(uint16_t local_rank0, uint64_t s_block) const {
    auto basic_block_rank = basic_block_rank_.as_words16();
    for (size_t pos = 0; pos < kBlocksPerSuperBlock; pos += 32) {
      auto count = lower_bound_delta_32x16(
          &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank0,
          delta_basic, kBasicBlockSize * pos);
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
    }
    return kBlocksPerSuperBlock * s_block + kBlocksPerSuperBlock - 1;
  }

  /**
   * @brief Interpolation search with SIMD optimization.
   * @param
   * local_rank Rank within the super block.
   * @param s_block Super block
   * index.
   * @details
   * Similar to find_basicblock but initial guess is
   * based on linear
   * interpolation, for random data it should make initial
   * guess correct
   * most of the times, we start from the 32 wide block with
   * interpolation
   * guess at the center, if we see that select result lie in
   * lower blocks
   * we backoff to find_basicblock
   */
  uint64_t find_basicblock_is(uint16_t local_rank, uint64_t s_block) const {
    auto super_block_rank = super_block_rank_.as_words64();
    auto basic_block_rank = basic_block_rank_.as_words16();

    auto lower = super_block_rank[s_block];
    auto upper = super_block_rank[s_block + 1];

    uint64_t pos = kBlocksPerSuperBlock * local_rank / (upper - lower);
    pos = pos + 16 < 32 ? 0 : (pos - 16);
    pos = pos > 96 ? 96 : pos;
    while (pos < 96) {
      auto count = lower_bound_32x16(
          &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank);
      if (count == 0) {
        return find_basicblock(local_rank, s_block);
      }
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
      pos += 32;
    }
    pos = 96;
    auto count = lower_bound_32x16(
        &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank);
    if (count == 0) {
      return find_basicblock(local_rank, s_block);
    }
    return kBlocksPerSuperBlock * s_block + pos + count - 1;
  }

  /**
   * @brief Locate a zero target basic block with a validated interpolation
   * estimate.
   * @param local_rank0 1-based zero rank within the super block.
   * @param s_block Super-block index.
   * @details The estimate is checked using the existing one-prefix metadata.
   * Non-uniform inputs fall back to the delta-SIMD search, so no zero-prefix
   * metadata is stored.
   */
  uint64_t find_basicblock_is_zeros(uint16_t local_rank0,
                                    uint64_t s_block) const {
    auto super_block_rank = super_block_rank_.as_words64();
    auto basic_block_rank = basic_block_rank_.as_words16();

    auto lower = kSuperBlockSize * s_block - super_block_rank[s_block];
    auto upper =
        kSuperBlockSize * (s_block + 1) - super_block_rank[s_block + 1];

    uint64_t interpolation =
        kBlocksPerSuperBlock * local_rank0 / (upper - lower);
    // Random data usually places the interpolation estimate in the target
    // block. Validate it from existing one-prefix metadata before the SIMD
    // derived-zero scan.
    const uint64_t block_offset =
        std::min(interpolation, kBlocksPerSuperBlock - 1);
    const uint64_t block = kBlocksPerSuperBlock * s_block + block_offset;
    const uint64_t zero_before =
        kBasicBlockSize * block_offset - basic_block_rank[block];
    const uint64_t zero_after = block_offset + 1 == kBlocksPerSuperBlock
                                    ? upper - lower
                                    : kBasicBlockSize * (block_offset + 1) -
                                          basic_block_rank[block + 1];
    if (zero_before < local_rank0 && local_rank0 <= zero_after) {
      return block;
    }

    uint64_t pos = interpolation;
    pos = pos + 16 < 32 ? 0 : (pos - 16);
    pos = pos > 96 ? 96 : pos;
    while (pos < 96) {
      auto count = lower_bound_delta_32x16(
          &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank0,
          delta_basic, kBasicBlockSize * pos);
      if (count == 0) {
        return find_basicblock_zeros(local_rank0, s_block);
      }
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
      pos += 32;
    }
    pos = 96;
    auto count = lower_bound_delta_32x16(
        &basic_block_rank[kBlocksPerSuperBlock * s_block + pos], local_rank0,
        delta_basic, kBasicBlockSize * pos);
    if (count == 0) {
      return find_basicblock_zeros(local_rank0, s_block);
    }
    return kBlocksPerSuperBlock * s_block + pos + count - 1;
  }

 public:
  RankSelectSupport() = default;
  RankSelectSupport(const RankSelectSupport&) = default;
  RankSelectSupport(RankSelectSupport&&) noexcept = default;
  RankSelectSupport& operator=(const RankSelectSupport&) = default;
  RankSelectSupport& operator=(RankSelectSupport&&) noexcept = default;

#ifdef PIXIE_DIAGNOSTICS
  struct DiagnosticsBytes {
    size_t source_bit_sequence_bytes = 0;
    size_t super_block_rank_bytes = 0;
    size_t basic_block_rank_bytes = 0;
    size_t select1_samples_bytes = 0;
    size_t select0_samples_bytes = 0;
    size_t total_bytes = 0;
  };

  /**
   * @brief Returns the number of bytes used by each internal component.
   */
  DiagnosticsBytes diagnostics_bytes() const {
    DiagnosticsBytes result;
    result.source_bit_sequence_bytes = (num_bits_ + 7) / 8;
    result.super_block_rank_bytes = super_block_rank_.as_bytes().size();
    result.basic_block_rank_bytes = basic_block_rank_.as_bytes().size();
    result.select1_samples_bytes = select1_sample_count_ * sizeof(uint64_t);
    result.select0_samples_bytes = select0_sample_count_ * sizeof(uint64_t);
    result.total_bytes = result.super_block_rank_bytes +
                         result.basic_block_rank_bytes +
                         select_samples_.as_bytes().size();
    return result;
  }

  /**
   * @brief Log memory usage of internal components.
   */
  void memory_report() const {
    const auto diagnostics = diagnostics_bytes();
    const double source_bytes =
        static_cast<double>(diagnostics.source_bit_sequence_bytes);
    const auto log_bytes = [&](std::string_view label, size_t bytes) {
      const double percentage =
          source_bytes > 0.0 ? 100.0 * static_cast<double>(bytes) / source_bytes
                             : 0.0;
      spdlog::info("RankSelectSupport {}: {} bytes ({:.2f}% of source)", label,
                   bytes, percentage);
    };
    log_bytes("source_bit_sequence", diagnostics.source_bit_sequence_bytes);
    log_bytes("super_block_rank", diagnostics.super_block_rank_bytes);
    log_bytes("basic_block_rank", diagnostics.basic_block_rank_bytes);
    log_bytes("select1_samples", diagnostics.select1_samples_bytes);
    log_bytes("select0_samples", diagnostics.select0_samples_bytes);
    log_bytes("total", diagnostics.total_bytes);
  }
#endif
  /**
   * @brief Construct support over a read-only storage view.
   * @param source_storage Packed source bytes, not owned.
   * @param num_bits Number of valid source bits.
   * @param select_support Which select sample tables to build. Rank and rank0
   * remain available in all modes.
   * @param one_count Optional exact number of 1-bits. When supplied,
   * construction can allocate select sample storage exactly.
   */
  explicit RankSelectSupport(
      ReadOnlyStorageView source_storage,
      size_t num_bits,
      SelectSupport select_support = SelectSupport::kBoth,
      std::optional<size_t> one_count = std::nullopt)
      : source_storage_(source_storage),
        bits_(source_storage_.as_words64()),
        num_bits_(std::min(num_bits, bits_.size() * kWordSize)),
        padded_size_(((num_bits_ + kWordSize - 1) / kWordSize) * kWordSize) {
    build_rank_select(select_support, one_count);
  }

  /**
   * @brief Construct support over caller-owned packed 64-bit words.
   * @param source_words Packed source words, not owned.
   * @param num_bits Number of valid source bits.
   * @param select_support Select sample tables to build.
   * @param one_count Optional exact number of one bits.
   */
  explicit RankSelectSupport(
      std::span<const uint64_t> source_words,
      size_t num_bits,
      SelectSupport select_support = SelectSupport::kBoth,
      std::optional<size_t> one_count = std::nullopt)
      : RankSelectSupport(ReadOnlyStorageView(std::as_bytes(source_words)),
                          num_bits,
                          select_support,
                          one_count) {}

  /**
   * @brief Construct support over a Pixie storage implementation.
   * @details The support retains a non-owning view. The source storage must
   * remain alive and must not be resized while this object is used.
   * @param source_storage Packed source storage, not owned.
   * @param num_bits Number of valid source bits.
   * @param select_support Select sample tables to build.
   * @param one_count Optional exact number of one bits.
   */
  template <StorageImplementation SourceStorage>
  explicit RankSelectSupport(
      const SourceStorage& source_storage,
      size_t num_bits,
      SelectSupport select_support = SelectSupport::kBoth,
      std::optional<size_t> one_count = std::nullopt)
      : RankSelectSupport(source_storage.view(),
                          num_bits,
                          select_support,
                          one_count) {}

  /**
   * @brief Returns the number of valid bits.
   */
  size_t size_impl() const { return num_bits_; }

  /**
   * @brief Whether this index stores samples for select1 queries.
   */
  bool supports_select1_impl() const { return builds_select1(select_support_); }

  /**
   * @brief Whether this index stores samples for select0 queries.
   */
  bool supports_select0_impl() const { return builds_select0(select_support_); }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this rank/select object and its owned aligned metadata
   * buffers. The external source bit words are not owned and are excluded.
   */
  size_t memory_usage_bytes_impl() const
    requires requires(const MetadataStorage& storage) {
      storage.allocated_bytes();
    }
  {
    return sizeof(*this) + super_block_rank_.allocated_bytes() +
           basic_block_rank_.allocated_bytes() +
           select_samples_.allocated_bytes();
  }

  /**
   * @brief Returns the bit at the given position.
   * @param pos Bit
   * index in [0, size()).
   */
  int bit_impl(size_t pos) const {
    size_t word_idx = pos / kWordSize;
    size_t bit_off = pos % kWordSize;

    return (bits_[word_idx] >> bit_off) & 1;
  }

  /**
   * @brief Rank of 1s up to position pos (exclusive).
   * @param pos Bit
   * index in [0, size()].
   * @return Number of 1s in [0, pos).
   */
  uint64_t rank_impl(size_t pos) const {
    if (pos >= num_bits_) [[unlikely]] {
      return max_rank_;
    }

    auto super_block_rank = super_block_rank_.as_words64();
    auto basic_block_rank = basic_block_rank_.as_words16();

    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = pos / kSuperBlockSize;
    // Precomputed rank
    uint64_t result = super_block_rank[s_block] + basic_block_rank[b_block];
    // Basic block tail
    result += rank_in_basic_block(b_block, pos - (b_block * kBasicBlockSize));
    return result;
  }

  /**
   * @brief Rank of 0s up to position pos (exclusive).
   * @param pos Bit
   * index in [0, size()].
   * @return Number of 0s in [0, pos).
   */
  uint64_t select_impl(size_t rank) const {
    if (rank == 0) [[unlikely]] {
      return 0;
    }
    if (!supports_select1_impl()) [[unlikely]] {
      return num_bits_;
    }
    if (rank > max_rank_) [[unlikely]] {
      return num_bits_;
    }
    auto super_block_rank = super_block_rank_.as_words64();
    auto basic_block_rank = basic_block_rank_.as_words16();

    uint64_t s_block = find_superblock(rank);
    rank -= super_block_rank[s_block];
    auto pos = find_basicblock_is(rank, s_block);
    rank -= basic_block_rank[pos];
    return select_in_words(pos * kWordsPerBlock, rank, true);
  }

  /**
   * @brief Select the position of the rank0-th 0-bit (1-indexed).
   *
   * @param rank0 1-based rank of the 0-bit to select.
   * @return Bit index,
   * or size() if rank0 is out of range.
   */
  uint64_t select0_impl(size_t rank0) const {
    if (rank0 == 0) [[unlikely]] {
      return 0;
    }
    if (!supports_select0_impl()) [[unlikely]] {
      return num_bits_;
    }
    if (rank0 > num_bits_ - max_rank_) [[unlikely]] {
      return num_bits_;
    }
    auto super_block_rank = super_block_rank_.as_words64();
    auto basic_block_rank = basic_block_rank_.as_words16();

    uint64_t s_block = find_superblock_zeros(rank0);
    rank0 -= kSuperBlockSize * s_block - super_block_rank[s_block];
    auto pos = find_basicblock_is_zeros(rank0, s_block);
    auto pos_in_super_block = pos & (kBlocksPerSuperBlock - 1);
    rank0 -= kBasicBlockSize * pos_in_super_block - basic_block_rank[pos];
    return select_in_words(pos * kWordsPerBlock, rank0, false);
  }

  void serialize(pixie::OutputBitStream& bs) const {
    bs << num_bits_ << padded_size_ << max_rank_ << select1_sample_begin_
       << select1_sample_count_ << select0_sample_begin_
       << select0_sample_count_ << static_cast<uint32_t>(select_support_)
       << static_cast<uint32_t>(select0_samples_reversed_);
    for (const uint64_t delta : delta_super) {
      bs << delta;
    }
    for (const uint16_t delta : delta_basic) {
      bs << delta;
    }
    super_block_rank_.serialize(bs);
    basic_block_rank_.serialize(bs);
    select_samples_.serialize(bs);
  }

  static RankSelectSupport deserialize(std::span<const uint64_t> source_bits,
                                       std::span<const std::byte>& data)
    requires std::same_as<MetadataStorage, ReadOnlyStorageView>
  {
    RankSelectSupport result;
    result.source_storage_ = ReadOnlyStorageView(std::as_bytes(source_bits));
    result.bits_ = result.source_storage_.as_words64();
    auto read = [&data](auto& value) {
      constexpr size_t length = sizeof(value);
      std::memcpy(&value, data.data(), length);
      data = data.subspan(length);
    };
    read(result.num_bits_);
    read(result.padded_size_);
    read(result.max_rank_);
    read(result.select1_sample_begin_);
    read(result.select1_sample_count_);
    read(result.select0_sample_begin_);
    read(result.select0_sample_count_);
    uint32_t buf;
    read(buf);
    result.select_support_ = static_cast<SelectSupport>(buf);
    read(buf);
    result.select0_samples_reversed_ = static_cast<bool>(buf);
    for (uint64_t& delta : result.delta_super) {
      read(delta);
    }
    for (uint16_t& delta : result.delta_basic) {
      read(delta);
    }
    result.super_block_rank_ = ReadOnlyStorageView::deserialize(data);
    result.basic_block_rank_ = ReadOnlyStorageView::deserialize(data);
    result.select_samples_ = ReadOnlyStorageView::deserialize(data);
    return result;
  }
};

}  // namespace pixie
