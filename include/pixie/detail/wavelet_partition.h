#pragma once

/**
 * @file wavelet_partition.h
 * @brief Bulk direction-mask and stable-partition primitives.
 */

#include <pixie/bits.h>
#include <pixie/packed_bit_builder.h>

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>

#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
#include <immintrin.h>
#endif

namespace pixie::detail {

template <std::unsigned_integral Symbol>
std::uint64_t wavelet_direction_mask(std::span<const Symbol> ranks,
                                     Symbol middle) {
  std::uint64_t mask = 0;
  std::size_t offset = 0;

  if constexpr (std::same_as<Symbol, std::uint8_t>) {
#if defined(PIXIE_AVX512_SUPPORT)
    if (ranks.size() == 64) {
      const __m512i values = _mm512_loadu_si512(ranks.data());
      const __m512i threshold = _mm512_set1_epi8(static_cast<char>(middle));
      return _mm512_cmpge_epu8_mask(values, threshold);
    }
#endif
#if defined(PIXIE_AVX2_SUPPORT)
    const __m256i threshold = _mm256_set1_epi8(static_cast<char>(middle));
    for (; offset + 32 <= ranks.size(); offset += 32) {
      const __m256i values = _mm256_loadu_si256(
          reinterpret_cast<const __m256i*>(ranks.data() + offset));
      const __m256i greater_or_equal =
          _mm256_cmpeq_epi8(_mm256_max_epu8(values, threshold), values);
      mask |= static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                  _mm256_movemask_epi8(greater_or_equal)))
              << offset;
    }
#endif
  }

  for (; offset < ranks.size(); ++offset) {
    mask |= static_cast<std::uint64_t>(ranks[offset] >= middle) << offset;
  }
  return mask;
}

/**
 * @brief Build packed node directions and stably partition ranks.
 * @param input Ranks in the node's original subsequence order.
 * @param middle First rank belonging to the right child.
 * @param output Destination split into left then right subsequences.
 * @param expected_left Number of ranks expected in the left child.
 * @param write_left Whether the left subsequence is needed by a child node.
 * @param write_right Whether the right subsequence is needed by a child node.
 * @param directions Packed destination for one direction bit per input rank.
 */
template <std::unsigned_integral Symbol>
void partition_wavelet_ranks(std::span<const Symbol> input,
                             Symbol middle,
                             std::span<Symbol> output,
                             std::size_t expected_left,
                             bool write_left,
                             bool write_right,
                             PackedBitBuilder& directions) {
  const bool output_required = write_left || write_right;
  if ((output_required && input.size() != output.size()) ||
      expected_left > input.size()) {
    throw std::invalid_argument("Invalid wavelet partition buffers");
  }

  std::size_t left = 0;
  std::size_t right = expected_left;
  for (std::size_t offset = 0; offset < input.size(); offset += 64) {
    const std::size_t width = std::min<std::size_t>(64, input.size() - offset);
    const std::span<const Symbol> block = input.subspan(offset, width);
    const std::uint64_t right_mask = wavelet_direction_mask(block, middle);
    directions.write_bits(right_mask, width);

    const std::uint64_t valid_mask =
        width == 64 ? std::numeric_limits<std::uint64_t>::max()
                    : (std::uint64_t{1} << width) - 1;
    std::uint64_t left_mask = (~right_mask) & valid_mask;
    if (write_left) {
      while (left_mask != 0) {
        const unsigned index = std::countr_zero(left_mask);
        output[left++] = block[index];
        left_mask &= left_mask - 1;
      }
    } else {
      left += std::popcount(left_mask);
    }

    std::uint64_t remaining_right = right_mask & valid_mask;
    if (write_right) {
      while (remaining_right != 0) {
        const unsigned index = std::countr_zero(remaining_right);
        output[right++] = block[index];
        remaining_right &= remaining_right - 1;
      }
    } else {
      right += std::popcount(remaining_right);
    }
  }

  if (left != expected_left || right != input.size()) {
    throw std::invalid_argument(
        "Wavelet partition does not match the supplied symbol counts");
  }
}

}  // namespace pixie::detail
