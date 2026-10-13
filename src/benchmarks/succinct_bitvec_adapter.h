#pragma once

#ifndef PIXIE_SUCCINCT_BITVEC_SUPPORT
#error "SuccinctFastRankSelectSupport requires PIXIE_THIRD_PARTY_BACKENDS"
#endif

#include <algorithm>
#include <bitvec.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>

namespace pixie::benchmarks {

/**
 * @brief Benchmark-only adapter for succinct-bitvec-cpp's FastBitVec.
 * @details Owns a masked source copy and the upstream index. Upstream ranks
 * include their bit position and selects have zero-based occurrence numbers;
 * this adapter uses half-open prefixes and one-based occurrences. Rank zero
 * selects return zero; absent one occurrences return size(). Zero rank is
 * derived from one rank; zero select is unavailable. Upstream's uint32_t sample
 * positions limit the logical source to 2^32 bits. Copy and move are disabled
 * because FastBitVec retains a reference to its RawBits object.
 */
class SuccinctFastRankSelectSupport {
 public:
  /** @brief Largest supported source; its last bit fits a uint32_t sample. */
  static constexpr std::uint64_t kMaxBits = std::uint64_t{1} << 32;

  /**
   * @brief Build FastBitVec over an owned least-significant-bit-first source.
   * @param source_words Packed complete 64-bit words, including optional
   * padding.
   * @param num_bits Logical bit count, clamped to the available words.
   * @throws std::length_error If the clamped source exceeds 2^32 bits.
   */
  explicit SuccinctFastRankSelectSupport(
      std::span<const std::uint64_t> source_words,
      std::size_t num_bits)
      : bits_(checked_size(source_words, num_bits)),
        support_(copy_source(bits_, source_words)) {}

  SuccinctFastRankSelectSupport(const SuccinctFastRankSelectSupport&) = delete;
  SuccinctFastRankSelectSupport& operator=(
      const SuccinctFastRankSelectSupport&) = delete;
  SuccinctFastRankSelectSupport(SuccinctFastRankSelectSupport&&) = delete;
  SuccinctFastRankSelectSupport& operator=(SuccinctFastRankSelectSupport&&) =
      delete;

  /** @brief Return the logical source length. */
  std::size_t size() const { return bits_.len; }

  /** @brief Return a source bit; position must be below size(). */
  int operator[](std::size_t position) const { return bits_.get(position); }

  /** @brief Count ones in [0, end), clamping end to size(). */
  std::uint64_t rank(std::size_t end) const {
    if (end >= size()) {
      return support_.get_total_ones();
    }
    return end == 0 ? 0 : support_.rank1(end - 1);
  }

  /** @brief Count zeros in [0, end), clamping end to size(). */
  std::uint64_t rank0(std::size_t end) const {
    return std::min(end, size()) - rank(end);
  }

  /** @brief Return a one-based one's position, zero for zero, size() if absent.
   */
  std::uint64_t select(std::size_t occurrence) const {
    if (occurrence == 0) {
      return 0;
    }
    if (occurrence > support_.get_total_ones()) {
      return size();
    }
    return static_cast<std::uint64_t>(support_.select1(occurrence - 1));
  }

  /** @brief FastBitVec always builds one-select metadata. */
  bool supports_select1() const { return true; }
  /** @brief The pinned upstream implementation has no zero-select API. */
  bool supports_select0() const { return false; }

  /** @brief Return allocated bytes in the normalized source copy. */
  std::size_t source_copy_bytes() const {
    return bits_.words.capacity() * sizeof(std::uint64_t);
  }

  /**
   * @brief Return object plus upstream-reported logical index bytes.
   * @details FastBitVec exposes directory sizes, not allocator capacities.
   */
  std::size_t index_logical_bytes() const {
    return sizeof(*this) + support_.overhead_bytes();
  }

  /** @brief Return reported index bytes plus the allocated source copy. */
  std::size_t memory_usage_bytes() const {
    return index_logical_bytes() + source_copy_bytes();
  }

 private:
  static std::size_t checked_size(std::span<const std::uint64_t> source,
                                  std::size_t num_bits) {
    if (source.size() <= std::numeric_limits<std::size_t>::max() / 64) {
      num_bits = std::min(num_bits, source.size() * 64);
    }
    if (num_bits > kMaxBits) {
      throw std::length_error(
          "FastBitVec select samples support at most 2^32 bits");
    }
    return num_bits;
  }

  static const succinct::RawBits& copy_source(
      succinct::RawBits& bits,
      std::span<const std::uint64_t> source) {
    const std::size_t word_count = std::min(bits.num_words, source.size());
    std::copy_n(source.begin(), word_count, bits.words.begin());
    if (bits.len % 64 != 0) {
      bits.words.back() &= (std::uint64_t{1} << (bits.len % 64)) - 1;
    }
    return bits;
  }

  succinct::RawBits bits_;
  succinct::FastBitVec support_;
};

}  // namespace pixie::benchmarks
