#pragma once

/**
 * @file slot_order.h
 * @brief Experimental split-bit and byte orders for 32/64 physical slots.
 * @details Split order encoding occupies 20/48 bytes; complete objects occupy
 * 24/56 bytes, versus 33/65 for byte orders. Physical occupancy is dense, so
 * count represents it without a stored bitmap. Split-bit rotation is scalar;
 * byte rotation dispatches locally to supported SIMD or a scalar fallback.
 * Primitive probes exclude tree traversal, subtree measures, and allocation.
 * Current primitive and integrated measurements are retained in
 * benchmarks/sequence_snapshot.md.
 */

#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace pixie::experimental {
namespace slot_order_detail {

template <class Word>
constexpr Word mask(std::size_t bits) noexcept {
  return bits >= sizeof(Word) * 8 ? ~Word{0} : (Word{1} << bits) - 1;
}

// Rotate a bounded field stream. Unchanged bits, including padding, survive.
template <class Word, std::size_t Words>
void rotate_bits(std::array<Word, Words>& words,
                 std::size_t left,
                 std::size_t length,
                 std::size_t distance) noexcept {
  constexpr auto width = sizeof(Word) * 8;
  if constexpr (Words == 1) {
    const auto region_mask = mask<Word>(length);
    const auto region = (words[0] >> left) & region_mask;
    const auto rotated = (region >> distance) | (region << (length - distance));
    words[0] =
        (words[0] & ~(region_mask << left)) | ((rotated & region_mask) << left);
  }
#ifdef __SIZEOF_INT128__
  else if constexpr (sizeof(Word) == 8 && Words == 2) {
    using Wide = __uint128_t;
    const auto original = Wide(words[0]) | (Wide(words[1]) << 64);
    const auto region_mask = mask<Wide>(length);
    const auto region = (original >> left) & region_mask;
    const auto rotated = (region >> distance) | (region << (length - distance));
    const auto result =
        (original & ~(region_mask << left)) | ((rotated & region_mask) << left);
    words[0] = Word(result);
    words[1] = Word(result >> 64);
  }
#endif
  else {
    const auto original = words;
    auto copy = [&](std::size_t source, std::size_t target, std::size_t count) {
      while (count != 0) {
        const auto source_bit = source % width;
        const auto target_bit = target % width;
        const auto chunk =
            std::min({count, width - source_bit, width - target_bit});
        const auto bits = mask<Word>(chunk);
        const auto value = (original[source / width] >> source_bit) & bits;
        auto& destination = words[target / width];
        destination =
            (destination & ~(bits << target_bit)) | (value << target_bit);
        source += chunk;
        target += chunk;
        count -= chunk;
      }
    };
    copy(left + distance, left, length - distance);
    copy(left, left + length - distance, distance);
  }
}

template <std::size_t N>
inline constexpr auto identity = [] {
  std::array<std::uint8_t, N> result{};
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = i;
  }
  return result;
}();

// SIMD indices select source bytes for the same half-open rotation as
// std::rotate.
template <std::size_t N>
void rotate_bytes(std::array<std::uint8_t, N>& bytes,
                  std::size_t left,
                  std::size_t right,
                  std::size_t distance) noexcept {
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI__)
  if constexpr (N == 64) {
    const auto positions = _mm512_loadu_si512(identity<N>.data());
    auto source = _mm512_add_epi8(positions, _mm512_set1_epi8(distance));
    source = _mm512_mask_sub_epi8(
        source, _mm512_cmpge_epi8_mask(source, _mm512_set1_epi8(right)), source,
        _mm512_set1_epi8(right - left));
    const auto active =
        _mm512_cmpge_epi8_mask(positions, _mm512_set1_epi8(left)) &
        _mm512_cmplt_epi8_mask(positions, _mm512_set1_epi8(right));
    source = _mm512_mask_mov_epi8(positions, active, source);
    const auto input = _mm512_loadu_si512(bytes.data());
    _mm512_storeu_si512(bytes.data(), _mm512_permutexvar_epi8(source, input));
    return;
  }
#endif
#ifdef PIXIE_AVX2_SUPPORT
  const auto input0 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bytes.data()));
  const auto bank0 = _mm256_permute2x128_si256(input0, input0, 0x00);
  const auto bank1 = _mm256_permute2x128_si256(input0, input0, 0x11);
  __m256i bank2{}, bank3{};
  if constexpr (N == 64) {
    const auto input1 =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(bytes.data() + 32));
    bank2 = _mm256_permute2x128_si256(input1, input1, 0x00);
    bank3 = _mm256_permute2x128_si256(input1, input1, 0x11);
  }
  for (std::size_t offset = 0; offset < N; offset += 32) {
    const auto positions = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(identity<N>.data() + offset));
    auto source = _mm256_add_epi8(positions, _mm256_set1_epi8(distance));
    source = _mm256_sub_epi8(
        source,
        _mm256_and_si256(_mm256_cmpgt_epi8(source, _mm256_set1_epi8(right - 1)),
                         _mm256_set1_epi8(right - left)));
    const auto active = _mm256_and_si256(
        _mm256_cmpgt_epi8(positions,
                          _mm256_set1_epi8(static_cast<int>(left) - 1)),
        _mm256_cmpgt_epi8(_mm256_set1_epi8(right), positions));
    source = _mm256_blendv_epi8(positions, source, active);
    const auto upper = _mm256_cmpgt_epi8(
        _mm256_and_si256(source, _mm256_set1_epi8(16)), _mm256_setzero_si256());
    auto result = _mm256_blendv_epi8(_mm256_shuffle_epi8(bank0, source),
                                     _mm256_shuffle_epi8(bank1, source), upper);
    if constexpr (N == 64) {
      const auto second =
          _mm256_blendv_epi8(_mm256_shuffle_epi8(bank2, source),
                             _mm256_shuffle_epi8(bank3, source), upper);
      result = _mm256_blendv_epi8(
          result, second, _mm256_cmpgt_epi8(source, _mm256_set1_epi8(31)));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(bytes.data() + offset),
                        result);
  }
#elif defined(PIXIE_SSE41_SUPPORT)
  __m128i banks[N / 16];
  for (std::size_t i = 0; i < N / 16; ++i) {
    banks[i] = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(bytes.data() + 16 * i));
  }
  for (std::size_t offset = 0; offset < N; offset += 16) {
    const auto positions = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(identity<N>.data() + offset));
    auto source = _mm_add_epi8(positions, _mm_set1_epi8(distance));
    source = _mm_sub_epi8(
        source, _mm_and_si128(_mm_cmpgt_epi8(source, _mm_set1_epi8(right - 1)),
                              _mm_set1_epi8(right - left)));
    const auto active = _mm_and_si128(
        _mm_cmpgt_epi8(positions, _mm_set1_epi8(static_cast<int>(left) - 1)),
        _mm_cmpgt_epi8(_mm_set1_epi8(right), positions));
    source = _mm_blendv_epi8(positions, source, active);
    auto result = _mm_shuffle_epi8(banks[0], source);
    for (std::size_t i = 1; i < N / 16; ++i) {
      result =
          _mm_blendv_epi8(result, _mm_shuffle_epi8(banks[i], source),
                          _mm_cmpgt_epi8(source, _mm_set1_epi8(16 * i - 1)));
    }
    _mm_storeu_si128(reinterpret_cast<__m128i*>(bytes.data() + offset), result);
  }
#else
  std::rotate(bytes.begin() + left, bytes.begin() + left + distance,
              bytes.begin() + right);
#endif
}
}  // namespace slot_order_detail

/**
 * @brief Experimental 32/64-slot order with separate low nibbles and high bits.
 * @details Owns only ordering metadata. Physical slots [0,count) stay occupied;
 * rotations never move payload. No holes are introduced. Invalid indices or
 * ranges violate assertion-checked preconditions. Every operation is noexcept
 * and allocation-free. This initial split-bit implementation is scalar.
 * @tparam N Capacity, either 32 (5 bits/slot) or 64 (6 bits/slot).
 */
template <std::size_t N>
class SplitSlotOrder {
  static_assert(N == 32 || N == 64);
  static constexpr std::size_t high_bits = N == 32 ? 1 : 2;
  using HighWord = std::conditional_t<N == 32, std::uint32_t, std::uint64_t>;
  std::array<std::uint64_t, N / 16> low_{};
  std::array<HighWord, N * high_bits / (sizeof(HighWord) * 8)> high_{};
  std::uint8_t count_;

 public:
  /** @brief Maximum number of occupied slots. */
  static constexpr std::size_t capacity = N;
  /** @brief Packed order bytes, excluding count and C++ alignment padding. */
  static constexpr std::size_t order_bytes = N * (4 + high_bits) / 8;
  /** @brief Identity order for physical [0,count); requires count <= N. */
  explicit SplitSlotOrder(std::size_t count = 0) noexcept : count_(count) {
    assert(count <= N);
    for (std::size_t i = 0; i < count; ++i) {
      low_[i / 16] |= std::uint64_t(i & 15) << (4 * (i % 16));
      constexpr auto width = sizeof(HighWord) * 8;
      high_[i * high_bits / width] |= HighWord(i >> 4)
                                      << (i * high_bits % width);
    }
  }
  /** @brief Occupied count, in [0,N]. */
  std::size_t size() const noexcept { return count_; }
  /** @brief Physical slot for a logical index strictly below size(). */
  unsigned operator[](std::size_t i) const noexcept {
    assert(i < size());
    constexpr auto width = sizeof(HighWord) * 8;
    const auto high =
        (high_[i * high_bits / width] >> (i * high_bits % width)) &
        ((1u << high_bits) - 1);
    return unsigned((low_[i / 16] >> (4 * (i % 16))) & 15) |
           (unsigned(high) << 4);
  }
  /** @brief Whether physical slot i is occupied; requires i < N. */
  bool occupied(std::size_t i) const noexcept {
    assert(i < N);
    return i < size();
  }
  /** @brief Check mapping uniqueness and occupancy. */
  bool valid() const noexcept {
    std::uint64_t seen = 0;
    for (std::size_t i = 0; i < size(); ++i) {
      const auto value = (*this)[i];
      if (value >= size() || (seen & (std::uint64_t{1} << value))) {
        return false;
      }
      seen |= std::uint64_t{1} << value;
    }
    return size() <= N;
  }
  /** @brief Rotate [left,right) by distance modulo its length; empty is a
   * no-op.
   * @pre left <= right <= size(). Physical slots remain unchanged. */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return;
    }
    slot_order_detail::rotate_bits(low_, 4 * left, 4 * length, 4 * distance);
    slot_order_detail::rotate_bits(high_, high_bits * left, high_bits * length,
                                   high_bits * distance);
  }
  /** @brief Apply rotate_left's contract without touching caller-owned
   * pointers.
   * @return Zero physical pointer writes. */
  std::size_t rotate_children(std::array<void*, N>&,
                              std::size_t left,
                              std::size_t right,
                              std::size_t distance) noexcept {
    rotate_left(left, right, distance);
    return 0;
  }
};

/**
 * @brief Byte-index control for the split-bit slot-order experiment.
 * @details Shares SplitSlotOrder's ranges, occupancy, ownership and exception
 * contracts. SIMD dispatch stays in this primitive; false selects the scalar
 * comparison implementation on every host.
 * @tparam N Capacity, 32 or 64. @tparam Simd Permit available SIMD shuffles.
 */
template <std::size_t N, bool Simd = true>
class ByteSlotOrder {
  static_assert(N == 32 || N == 64);
  std::array<std::uint8_t, N> order_ = slot_order_detail::identity<N>;
  std::uint8_t count_;

 public:
  /** @brief Maximum occupied count. */
  static constexpr std::size_t capacity = N;
  /** @brief Order bytes, excluding count and alignment padding. */
  static constexpr std::size_t order_bytes = N;
  /** @brief Construct identity order; requires count <= N. */
  explicit ByteSlotOrder(std::size_t count = 0) noexcept : count_(count) {
    assert(count <= N);
  }
  /** @brief Occupied count in [0,N]. */
  std::size_t size() const noexcept { return count_; }
  /** @brief Physical index; requires i < size(). */
  unsigned operator[](std::size_t i) const noexcept {
    assert(i < size());
    return order_[i];
  }
  /** @brief Whether physical slot i is occupied; requires i < N. */
  bool occupied(std::size_t i) const noexcept {
    assert(i < N);
    return i < size();
  }
  /** @brief Check mapping uniqueness and occupancy. */
  bool valid() const noexcept {
    std::uint64_t seen = 0;
    for (std::size_t i = 0; i < size(); ++i) {
      const auto value = (*this)[i];
      if (value >= size() || (seen & (std::uint64_t{1} << value))) {
        return false;
      }
      seen |= std::uint64_t{1} << value;
    }
    return size() <= N;
  }
  /** @brief Rotate [left,right); reduce distance modulo length; empty is a
   * no-op.
   * @pre left <= right <= size(). */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return;
    }
    if constexpr (Simd) {
      slot_order_detail::rotate_bytes(order_, left, right, distance);
    } else {
      std::rotate(order_.begin() + left, order_.begin() + left + distance,
                  order_.begin() + right);
    }
  }
  /** @brief Rotate metadata using rotate_left's contract; pointers are
   * untouched.
   * @return Zero physical pointer writes. */
  std::size_t rotate_children(std::array<void*, N>&,
                              std::size_t left,
                              std::size_t right,
                              std::size_t distance) noexcept {
    rotate_left(left, right, distance);
    return 0;
  }
};
}  // namespace pixie::experimental
