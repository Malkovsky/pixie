#pragma once

#ifndef PIXIE_3STAR_SUPPORT
#error "ThreeStarRankSelectSupport requires PIXIE_3STAR_SOURCE_DIR"
#endif

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

namespace pixie::benchmarks {

/**
 * @brief Benchmark-only adapter for a permitted local 3-star checkout.
 *
 * @details This adapter exposes only the upstream operations that are complete
 * in the pinned local source: one-rank and one-select. Its backing input and
 * metadata are owned by the upstream process-global implementation, so only
 * one instance may be live at a time. It is intentionally not a public Pixie
 * rank/select implementation. Upstream construction diagnostics are suppressed
 * to keep terminal I/O out of the timed build operation.
 */
class ThreeStarRankSelectSupport {
 public:
  /**
   * @brief Build the upstream uncompressed 3-star configuration.
   * @param source_words Packed least-significant-bit-first source words.
   * @param num_bits Logical source length, clamped to @p source_words.
   */
  explicit ThreeStarRankSelectSupport(
      std::span<const std::uint64_t> source_words,
      std::size_t num_bits);

  /** @brief Return the logical source length. */
  std::size_t size() const { return num_bits_; }

  /** @brief Return the one rank in `[0, end_position)`. */
  std::uint64_t rank(std::size_t end_position) const {
    return end_position >= num_bits_ ? one_count_
                                     : support_.rank_1(end_position);
  }

  /** @brief Return the derived zero rank in `[0, end_position)`. */
  std::uint64_t rank0(std::size_t end_position) const {
    return end_position >= num_bits_ ? num_bits_ - one_count_
                                     : end_position - rank(end_position);
  }

  /** @brief Return the zero-based position of a one-based one rank. */
  std::uint64_t select(std::size_t rank) const {
    if (rank == 0) {
      return 0;
    }
    return rank > one_count_ ? num_bits_ : support_.select_1(rank);
  }

  /** @brief 3-star implements one-select, but not zero-select. */
  bool supports_select1() const { return true; }

  /** @brief 3-star's upstream zero-select implementation is incomplete. */
  bool supports_select0() const { return false; }

  /** @brief Return upstream-reported logical bytes for source and metadata. */
  std::size_t memory_usage_bytes() const {
    return static_cast<std::size_t>(support_.space_in_bits() / 8);
  }

  /** @brief Return logical bytes in 3-star's uncompressed source copy. */
  std::size_t source_copy_bytes() const { return source_copy_bytes_; }

  /** @brief Return upstream-reported logical metadata bytes. */
  std::size_t index_logical_bytes() const {
    return memory_usage_bytes() - source_copy_bytes_;
  }

 private:
  m3 support_;
  std::size_t num_bits_ = 0;
  std::size_t one_count_ = 0;
  std::size_t source_copy_bytes_ = 0;
};

}  // namespace pixie::benchmarks
