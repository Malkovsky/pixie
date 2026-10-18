#pragma once

#ifndef PIXIE_3STAR_SUPPORT
#error "ThreeStarRankSelectSupport requires PIXIE_THIRD_PARTY_BACKENDS"
#endif

#include <pixie/rank_select.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <span>

// The upstream header emits unconditional construction diagnostics. Suppress
// only those diagnostics so benchmark construction measures the index rather
// than terminal I/O.
#define printf(...) static_cast<void>(0)
#include <m3.hpp>
#undef printf

static void print_overhead_line(const std::string&, int64_t) {}

static void print_separator() {}

namespace pixie {

/**
 * @brief Optional rank/select implementation over the 3-star backend.
 *
 * @details This adapter exposes only the upstream operations that are complete
 * in the pinned upstream source: one-rank and one-select. Its backing input and
 * metadata are owned by the upstream process-global implementation, so only
 * one instance may be live at a time. The comparison target supplies its
 * upstream global state. Upstream construction diagnostics are suppressed to
 * keep terminal I/O out of the timed build operation. The source copy includes
 * zero padding for complete superblocks and sentinels.
 */
class ThreeStarRankSelectSupport
    : public RankSelectBase<ThreeStarRankSelectSupport> {
 public:
  /**
   * @brief Build the upstream uncompressed 3-star configuration.
   * @param source_words Packed least-significant-bit-first source words.
   * @param num_bits Logical source length, clamped to @p source_words.
   */
  explicit ThreeStarRankSelectSupport(
      std::span<const std::uint64_t> source_words,
      std::size_t num_bits)
      : num_bits_(std::min(num_bits, source_words.size() * kWordBits)) {
    const std::size_t source_word_count =
        (num_bits_ + kWordBits - 1) / kWordBits;
    const std::size_t padded_word_count =
        ((source_word_count + kWordsPerSuperblock - 1) / kWordsPerSuperblock) *
        kWordsPerSuperblock;
    const std::size_t stored_words =
        padded_word_count + kSentinelSuperblocks * kWordsPerSuperblock;
    N_bits = num_bits_;
    // Upstream builds complete superblocks and copies sentinel blocks without
    // checking the input vector's end. Supply initialized padding for both.
    N_words = padded_word_count;
    bits.assign(stored_words, 0);
    one_count_ = 0;
    for (std::size_t word_index = 0; word_index < source_word_count;
         ++word_index) {
      const std::uint64_t word =
          logical_word(source_words, word_index, num_bits_);
      bits[word_index] = word;
      one_count_ += std::popcount(word);
    }
    N_ones = one_count_;
    N_zeros = num_bits_ - one_count_;

    // Match the paper's uncompressed robust configuration: L0=2048, a*,
    // alpha=16, and its two-level L0/L1 summary tree.
    ALPHA = 16;
    SUMMARY_LEVELS = 1;
    TREE3_STRATEGY = GET_LIFTED_CUBIC_THEORY_PARAMS;
    BV_COMPRESSION = 0;

    support_.build_auxiliaries();

    source_copy_bytes_ = stored_words * sizeof(std::uint64_t);
  }

  /** @brief Return the logical source length. */
  std::size_t size_impl() const { return num_bits_; }

  /** @brief Read a source bit; position must be below size(). */
  int bit_impl(std::size_t position) const {
    return (bits[position / kWordBits] >> (position % kWordBits)) & 1;
  }

  /** @brief Return the one rank in `[0, end_position)`. */
  std::uint64_t rank_impl(std::size_t end_position) const {
    return end_position >= num_bits_ ? one_count_
                                     : support_.rank_1(end_position);
  }

  /** @brief Return the zero-based position of a one-based one rank. */
  std::uint64_t select_impl(std::size_t rank) const {
    if (rank == 0) {
      return 0;
    }
    return rank > one_count_ ? num_bits_ : support_.select_1(rank);
  }

  /** @brief 3-star implements one-select, but not zero-select. */
  bool supports_select1_impl() const { return true; }

  /** @brief 3-star's upstream zero-select implementation is incomplete. */
  bool supports_select0_impl() const { return false; }

  /** @brief Return zero for occurrence zero, otherwise size(): disabled. */
  std::uint64_t select0_impl(std::size_t occurrence) const {
    return occurrence == 0 ? 0 : size();
  }

  /** @brief Return upstream-reported logical bytes for source and metadata. */
  std::size_t memory_usage_bytes_impl() const {
    return static_cast<std::size_t>(support_.space_in_bits() / 8);
  }

  /** @brief Return logical bytes in 3-star's uncompressed source copy. */
  std::size_t source_copy_bytes() const { return source_copy_bytes_; }

  /** @brief Return upstream-reported logical metadata bytes. */
  std::size_t index_logical_bytes() const {
    return memory_usage_bytes() - source_copy_bytes_;
  }

 private:
  static constexpr std::size_t kWordBits = 64;
  static constexpr std::size_t kWordsPerSuperblock = 1024;
  static constexpr std::size_t kSentinelSuperblocks = 5;

  static std::uint64_t logical_word(std::span<const std::uint64_t> source,
                                    std::size_t word_index,
                                    std::size_t num_bits) {
    std::uint64_t word = source[word_index];
    const std::size_t remaining = num_bits - word_index * kWordBits;
    if (remaining < kWordBits) {
      word &= (std::uint64_t{1} << remaining) - 1;
    }
    return word;
  }

  m3 support_;
  std::size_t num_bits_ = 0;
  std::size_t one_count_ = 0;
  std::size_t source_copy_bytes_ = 0;
};

}  // namespace pixie
