#pragma once

// Current experiment snapshot, 2026-09-16. Ryzen 7 8845HS, Linux/WSL,
// GCC 13.3, -O3 -DNDEBUG -march=native, Google Benchmark 1.9.4, CPU 0.
// Command: cpp-bench native bit_sequence_benchmarks <filter> 5 0.15s
// Filter: '^Block(Rotate|Read)<(DirectBlock|MappedBlock).*'
// Median CPU ns/op; randomized repetition interleaving. Hot single block,
// 256 precomputed rotations (seed 128), no construction/reset in timing.
// Modes 0/1: chunk-aligned/unaligned proper subranges; 2: alternating 50/50;
// 3: whole block by 65; 4: unaligned, N=2045, initial h=2038 (both variants
// retain their circular origin). Modes 1/4 retain identity maps; mode 2 keeps
// changing nonidentity maps across partial rotations without normalization.
// Read setup: two aligned subrange rotations, then whole-block rotation by 65;
// indexed reads use 1024 precomputed positions (seed 123). Flatten includes
// producing all 32 output words with zero padding, without mutating the block.
//
// | Workload         | N bits | Direct ns | Mapped ns |
// | ---------------- | -----: | --------: | --------: |
// | Rotate aligned   |   2048 |     16.85 |      3.82 |
// | Rotate unaligned |   2048 |     28.15 |     29.70 |
// | Rotate mixed     |   2048 |     22.03 |     48.95 |
// | Rotate whole     |   2048 |      3.29 |      3.29 |
// | Rotate wrapped   |   2045 |     29.65 |     31.11 |
// | Indexed access   |   2048 |      1.50 |      1.58 |
// | Indexed access   |   2045 |      1.50 |      1.56 |
// | Flatten          |   2048 |     10.57 |     99.69 |
// | Flatten          |   2045 |     11.36 |     99.55 |
//
// Direct/Mapped use the same 256-byte payload and packed copy primitive.
// Object bytes: Direct 272, Mapped 280. Map adds 8 bytes (3.125% of payload).
// Retain only as a bounded aligned-heavy experiment: mixed rotations and
// flattening remain substantially slower than direct. These are block-only
// results, not tree measurements. Chunk-wise gather/scatter avoids
// normalization, not map cost. Independent hot indexed reads measure
// throughput, not dependent latency; WSL timings do not establish cold-cache or
// large-sequence performance.

#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace pixie::experimental {

/**
 * @brief Bounded owning experiment, independent of sequence split/merge policy.
 *
 * @details Holds 0..2048 bits in 16 physical 128-bit chunks. With
 * Permuted=true, nibble j of a uint64_t maps canonical chunk j to its physical
 * chunk. Circular origin h is applied modulo the valid length BEFORE this
 * mapping. Whole-block rotation changes only h. At h=0, chunk-aligned subranges
 * rotate map nibbles without touching payload (including in partially filled
 * blocks). Other rotations gather only the affected interval in rotated order,
 * then scatter through the same map, preserving both the permutation and
 * origin. Copies translate contiguous chunk pieces, never individual bits.
 *
 * Permuted=false is the matched direct-movement control: the same packed copy
 * primitive and circular-origin fast path, without a map or descriptor lookup.
 * Value copies own independent payloads. No allocation or external lifetime.
 * @tparam Permuted Whether to map canonical chunks to physical chunks; false
 * selects the direct-movement control.
 */
template <bool Permuted = true>
class PermutedBitBlock {
 public:
  /** @brief Logical element type returned by value. */
  using value_type = bool;
  /** @brief Maximum logical element count. */
  static constexpr std::size_t capacity = 2048;
  /** @brief Fixed owning payload, in least-significant-bit-first word order. */
  using Payload = std::array<std::uint64_t, 32>;

  /** @brief Construct empty; all storage is initialized. */
  PermutedBitBlock() noexcept = default;

  /**
   * @brief Copy n bits, ignoring extra words and high padding bits.
   * @details The source need not outlive the block. Bits are copied in
   * least-significant-bit-first word order, with initial circular origin h = 0.
   * @param words Source words containing the bits to copy.
   * @param n Valid logical bit count, from 0 through 2048 inclusive.
   * @throws std::invalid_argument If n > 2048 or words does not cover n bits.
   */
  PermutedBitBlock(std::span<const std::uint64_t> words, std::size_t n)
      : n_(n) {
    if (n > 2048 || (n + 63) / 64 > words.size()) {
      throw std::invalid_argument("PermutedBitBlock: input size");
    }
    std::copy_n(words.begin(), (n + 63) / 64, bits_.begin());
  }

  /**
   * @brief Query the logical length.
   * @return Valid logical bit count, excluding capacity and padding.
   */
  std::size_t size() const noexcept { return n_; }
  /**
   * @brief Query whether the block is empty.
   * @return True if there are no valid bits.
   */
  bool empty() const noexcept { return n_ == 0; }

  /**
   * @brief Read a zero-based logical bit by value.
   * @param i Bit position in [0, size()).
   * @return The logical bit value.
   * @throws std::out_of_range If i >= size().
   */
  bool operator[](std::size_t i) const {
    if (i >= n_) {
      throw std::out_of_range("PermutedBitBlock: index");
    }
    const auto p = physical((h_ + i) % n_);
    return (bits_[p / 64] >> (p % 64)) & 1;
  }

  /**
   * @brief Copy logical bits to linear words without mutating the block.
   * @details Copies at most 256 bytes and does not allocate.
   * @return Owning LSB-first payload with all bits beyond size() set to zero.
   */
  Payload flatten() const noexcept {
    Payload result{};
    copy(0, result, 0, n_);
    return result;
  }

  /**
   * @brief Rotate [left, right) left by distance modulo its length.
   * @details Empty ranges and zero effective distances do nothing. Bits outside
   * the range remain unchanged.
   * @param left Inclusive zero-based start of the range.
   * @param right Exclusive end of the range.
   * @param distance Left rotation distance in bits, reduced modulo range
   * length.
   * @pre left <= right <= size(); invalid primitive arguments are asserted.
   */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= n_);
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return;
    }
    if (length == n_) {
      h_ = (h_ + distance) % n_;
      return;
    }
    if constexpr (Permuted) {
      if (h_ == 0 && ((left | right | distance) % 128 == 0)) {
        const auto shift = left / 128 * 4;
        const auto width = length / 128 * 4;
        const auto cut = distance / 128 * 4;
        // Proper subrange: width < 64; distance is strictly inside the range.
        const auto mask = (std::uint64_t{1} << width) - 1;
        const auto slice = (map_ >> shift) & mask;
        const auto rotated = (slice >> cut) | (slice << (width - cut));
        map_ = (map_ & ~(mask << shift)) | ((rotated & mask) << shift);
        return;
      }
    }
    Payload scratch{};
    copy(left + distance, scratch, 0, length - distance);
    copy(left, scratch, length - distance, distance);
    const auto position = (h_ + left) % n_;
    if constexpr (Permuted) {
      if (map_ != identity) {
        auto p = position;
        std::size_t source = 0;
        while (source != length) {
          const auto run = std::min({length - source, n_ - p, 128 - p % 128});
          copy_packed_bits(scratch.data(), source, bits_.data(), physical(p),
                           run);
          source += run;
          p = (p + run == n_) ? 0 : p + run;
        }
        return;
      }
    }
    const auto prefix = std::min(length, n_ - position);
    copy_packed_bits(scratch.data(), 0, bits_.data(), position, prefix);
    copy_packed_bits(scratch.data(), prefix, bits_.data(), 0, length - prefix);
  }

  /**
   * @brief Query fixed payload capacity.
   * @return Payload bytes, including unused capacity.
   */
  std::size_t payload_capacity_bytes() const noexcept { return sizeof(bits_); }
  /**
   * @brief Query object metadata size; there is no dynamic allocation.
   * @return Object bytes excluding the payload, including any padding.
   */
  std::size_t metadata_bytes() const noexcept {
    return sizeof(*this) - sizeof(bits_);
  }

  /**
   * @brief Repartition two distinct blocks without allocation or exceptions.
   * @details Preserves concatenated logical contents; length changes normalize
   * the maps and origins. Ordinary rotations retain their existing kernels.
   * Requires left_size <= capacity and left_size <= size()+rhs.size() <=
   * left_size+capacity. Invalid primitive arguments are asserted.
   * @param rhs Right block, receiving the suffix.
   * @param left_size Desired valid count in this block.
   */
  void redistribute(PermutedBitBlock& rhs, std::size_t left_size) noexcept {
    assert(this != &rhs);
    const auto total = n_ + rhs.n_;
    assert(left_size <= capacity && left_size <= total &&
           total - left_size <= capacity);
    Payload a{}, b{};
    auto gather = [&](std::size_t start, Payload& out, std::size_t count) {
      const auto first = start < n_ ? std::min(count, n_ - start) : 0;
      copy(start, out, 0, first);
      rhs.copy(start + first - std::min(start + first, n_), out, first,
               count - first);
    };
    gather(0, a, left_size);
    gather(left_size, b, total - left_size);
    bits_ = a;
    rhs.bits_ = b;
    n_ = left_size;
    rhs.n_ = total - left_size;
    h_ = rhs.h_ = 0;
    if constexpr (Permuted) {
      map_ = rhs.map_ = identity;
    }
  }

 private:
  static constexpr std::uint64_t identity = 0xfedcba9876543210ULL;
  struct NoMap {};

  std::size_t physical(std::size_t p) const noexcept {
    if constexpr (Permuted) {
      return ((map_ >> (p / 128 * 4)) & 15) * 128 + p % 128;
    }
    return p;
  }

  void copy(std::size_t offset,
            Payload& out,
            std::size_t dest,
            std::size_t count) const noexcept {
    if (count == 0) {
      return;
    }
    auto p = (h_ + offset) % n_;
    while (count != 0) {
      auto run = std::min(count, n_ - p);
      if constexpr (Permuted) {
        if (map_ != identity) {
          run = std::min(run, 128 - p % 128);
        }
      }
      copy_packed_bits(bits_.data(), physical(p), out.data(), dest, run);
      count -= run;
      dest += run;
      p = (p + run == n_) ? 0 : p + run;
    }
  }

  Payload bits_{};
  std::size_t n_ = 0;
  std::size_t h_ = 0;
  [[no_unique_address]] std::conditional_t<Permuted, std::uint64_t, NoMap>
      map_ = []() noexcept {
        if constexpr (Permuted) {
          return identity;
        } else {
          return NoMap{};
        }
      }();
};

}  // namespace pixie::experimental
