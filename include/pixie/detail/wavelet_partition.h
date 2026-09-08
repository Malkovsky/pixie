#pragma once

/**
 * @file wavelet_partition.h
 * @brief Bulk direction-mask and stable-partition primitives.
 */

#include <pixie/bits.h>
#include <pixie/packed_bit_builder.h>

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>

#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
#include <immintrin.h>
#endif

namespace pixie::detail {

#if defined(PIXIE_AVX2_SUPPORT)

consteval auto make_wavelet_compaction_table() {
  std::array<std::array<std::uint8_t, 16>, 256> table{};
  for (std::size_t mask = 0; mask < table.size(); ++mask) {
    auto& shuffle = table[mask];
    shuffle.fill(0x80);
    std::size_t left = 0;
    std::size_t right = 8;
    for (std::size_t index = 0; index < 8; ++index) {
      if ((mask & (std::size_t{1} << index)) == 0) {
        shuffle[left++] = static_cast<std::uint8_t>(index);
      } else {
        shuffle[right++] = static_cast<std::uint8_t>(index);
      }
    }
  }
  return table;
}

alignas(64) inline constexpr auto kWaveletCompactionTable =
    make_wavelet_compaction_table();

inline void store_compacted_bytes(std::uint8_t* destination,
                                  std::uint64_t packed,
                                  unsigned count) {
  switch (count) {
    case 8:
      std::memcpy(destination, &packed, 8);
      return;
    case 7: {
      const std::uint32_t suffix = static_cast<std::uint32_t>(packed >> 24);
      std::memcpy(destination, &packed, 4);
      std::memcpy(destination + 3, &suffix, 4);
      return;
    }
    case 6: {
      const std::uint16_t suffix = static_cast<std::uint16_t>(packed >> 32);
      std::memcpy(destination, &packed, 4);
      std::memcpy(destination + 4, &suffix, 2);
      return;
    }
    case 5:
      std::memcpy(destination, &packed, 4);
      destination[4] = static_cast<std::uint8_t>(packed >> 32);
      return;
    case 4:
      std::memcpy(destination, &packed, 4);
      return;
    case 3: {
      const std::uint16_t prefix = static_cast<std::uint16_t>(packed);
      std::memcpy(destination, &prefix, 2);
      destination[2] = static_cast<std::uint8_t>(packed >> 16);
      return;
    }
    case 2: {
      const std::uint16_t prefix = static_cast<std::uint16_t>(packed);
      std::memcpy(destination, &prefix, 2);
      return;
    }
    case 1:
      destination[0] = static_cast<std::uint8_t>(packed);
      return;
    default:
      return;
  }
}

// Full-width stores avoid a variable-size store in the hot loop. The logical
// cursor advances by only the valid byte count, so the next store replaces the
// preceding store's unused suffix. The allocation boundary still uses an exact
// store.
class OverlappingByteWriter {
 public:
  OverlappingByteWriter(std::uint8_t* destination, std::uint8_t* end)
      : destination_(destination), end_(end) {}

  void append(std::uint64_t packed, unsigned count) {
    if (count == 0) {
      return;
    }
    if (end_ - destination_ >= 8) {
      std::memcpy(destination_, &packed, 8);
    } else {
      store_compacted_bytes(destination_, packed, count);
    }
    destination_ += count;
  }

  void append(std::uint8_t value) { append(value, 1); }

 private:
  std::uint8_t* destination_;
  std::uint8_t* end_;
};

// A left-partition store can overlap the start of its adjacent right partition
// by at most seven bytes. Retain enough of the right stream to restore it once
// partitioning is complete.
class RightPartitionPrefix {
 public:
  void append(std::uint64_t packed, unsigned count) {
    if (size_ == bytes_.size() || count == 0) {
      return;
    }
    const std::size_t copied =
        std::min<std::size_t>(count, bytes_.size() - size_);
    std::memcpy(bytes_.data() + size_, &packed, copied);
    size_ += copied;
  }

  void append(std::uint8_t value) { append(value, 1); }

  void restore(std::uint8_t* destination) const {
    std::memcpy(destination, bytes_.data(), size_);
  }

 private:
  std::array<std::uint8_t, 8> bytes_{};
  std::size_t size_ = 0;
};

template <bool WriteLeft, bool WriteRight>
std::size_t compact_wavelet_block_avx2(std::span<const std::uint8_t> input,
                                       std::uint64_t right_mask,
                                       OverlappingByteWriter& left_output,
                                       OverlappingByteWriter& right_output,
                                       RightPartitionPrefix& right_prefix,
                                       std::size_t& left,
                                       std::size_t& right) {
  std::size_t offset = 0;
  for (; offset + 8 <= input.size(); offset += 8) {
    const auto mask = static_cast<std::uint8_t>(right_mask >> offset);
    const unsigned right_count = std::popcount(mask);
    const unsigned left_count = 8 - right_count;

    if constexpr (WriteLeft || WriteRight) {
      const __m128i values = _mm_loadl_epi64(
          reinterpret_cast<const __m128i*>(input.data() + offset));
      const __m128i shuffle = _mm_load_si128(reinterpret_cast<const __m128i*>(
          kWaveletCompactionTable[mask].data()));
      const __m128i compacted = _mm_shuffle_epi8(values, shuffle);
      if constexpr (WriteLeft) {
        const std::uint64_t packed =
            static_cast<std::uint64_t>(_mm_cvtsi128_si64(compacted));
        left_output.append(packed, left_count);
      }
      if constexpr (WriteRight) {
        const std::uint64_t packed = static_cast<std::uint64_t>(
            _mm_cvtsi128_si64(_mm_srli_si128(compacted, 8)));
        if constexpr (WriteLeft) {
          right_prefix.append(packed, right_count);
        }
        right_output.append(packed, right_count);
      }
    }

    left += left_count;
    right += right_count;
  }
  return offset;
}

#endif

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

template <bool WriteLeft, bool WriteRight, std::unsigned_integral Symbol>
void partition_wavelet_ranks_scalar(std::span<const Symbol> input,
                                    Symbol middle,
                                    std::span<Symbol> output,
                                    PackedBitBuilder& directions,
                                    std::size_t& left,
                                    std::size_t& right) {
  for (std::size_t offset = 0; offset < input.size(); offset += 64) {
    const std::size_t width = std::min<std::size_t>(64, input.size() - offset);
    const std::span<const Symbol> block = input.subspan(offset, width);
    const std::uint64_t right_mask = wavelet_direction_mask(block, middle);
    directions.write_bits(right_mask, width);

    const std::uint64_t valid_mask =
        width == 64 ? std::numeric_limits<std::uint64_t>::max()
                    : (std::uint64_t{1} << width) - 1;
    std::uint64_t left_mask = (~right_mask) & valid_mask;
    if constexpr (WriteLeft) {
      while (left_mask != 0) {
        const unsigned index = std::countr_zero(left_mask);
        output[left++] = block[index];
        left_mask &= left_mask - 1;
      }
    } else {
      left += std::popcount(left_mask);
    }

    std::uint64_t remaining_right = right_mask & valid_mask;
    if constexpr (WriteRight) {
      while (remaining_right != 0) {
        const unsigned index = std::countr_zero(remaining_right);
        output[right++] = block[index];
        remaining_right &= remaining_right - 1;
      }
    } else {
      right += std::popcount(remaining_right);
    }
  }
}

#if defined(PIXIE_AVX2_SUPPORT)

template <bool WriteLeft, bool WriteRight>
void partition_wavelet_bytes_avx2(std::span<const std::uint8_t> input,
                                  std::uint8_t middle,
                                  std::span<std::uint8_t> output,
                                  PackedBitBuilder& directions,
                                  std::size_t& left,
                                  std::size_t& right) {
  std::uint8_t* const output_end =
      output.empty() ? nullptr : output.data() + output.size();
  OverlappingByteWriter left_output(
      WriteLeft && !output.empty() ? output.data() + left : nullptr,
      output_end);
  OverlappingByteWriter right_output(
      WriteRight && !output.empty() ? output.data() + right : nullptr,
      output_end);
  RightPartitionPrefix right_prefix;
  const std::size_t right_begin = right;

  for (std::size_t offset = 0; offset < input.size(); offset += 64) {
    const std::size_t width = std::min<std::size_t>(64, input.size() - offset);
    const std::span<const std::uint8_t> block = input.subspan(offset, width);
    const std::uint64_t right_mask = wavelet_direction_mask(block, middle);
    directions.write_bits(right_mask, width);

    const std::size_t compacted =
        compact_wavelet_block_avx2<WriteLeft, WriteRight>(
            block, right_mask, left_output, right_output, right_prefix, left,
            right);
    const std::uint64_t valid_mask =
        width == 64 ? std::numeric_limits<std::uint64_t>::max()
                    : (std::uint64_t{1} << width) - 1;
    const std::uint64_t remaining_mask =
        compacted == 64
            ? 0
            : valid_mask &
                  (std::numeric_limits<std::uint64_t>::max() << compacted);
    std::uint64_t left_mask = (~right_mask) & remaining_mask;
    if constexpr (WriteLeft) {
      while (left_mask != 0) {
        const unsigned index = std::countr_zero(left_mask);
        left_output.append(block[index]);
        ++left;
        left_mask &= left_mask - 1;
      }
    } else {
      left += std::popcount(left_mask);
    }

    std::uint64_t remaining_right = right_mask & remaining_mask;
    if constexpr (WriteRight) {
      while (remaining_right != 0) {
        const unsigned index = std::countr_zero(remaining_right);
        if constexpr (WriteLeft) {
          right_prefix.append(block[index]);
        }
        right_output.append(block[index]);
        ++right;
        remaining_right &= remaining_right - 1;
      }
    } else {
      right += std::popcount(remaining_right);
    }
  }

  if constexpr (WriteLeft && WriteRight) {
    if (!output.empty()) {
      right_prefix.restore(output.data() + right_begin);
    }
  }
}

#endif

template <bool WriteLeft, bool WriteRight, std::unsigned_integral Symbol>
void partition_wavelet_ranks_impl(std::span<const Symbol> input,
                                  Symbol middle,
                                  std::span<Symbol> output,
                                  PackedBitBuilder& directions,
                                  std::size_t& left,
                                  std::size_t& right) {
#if defined(PIXIE_AVX2_SUPPORT)
  if constexpr (std::same_as<Symbol, std::uint8_t>) {
    partition_wavelet_bytes_avx2<WriteLeft, WriteRight>(
        input, middle, output, directions, left, right);
    return;
  }
#endif
  partition_wavelet_ranks_scalar<WriteLeft, WriteRight>(
      input, middle, output, directions, left, right);
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
  if (write_left) {
    if (write_right) {
      partition_wavelet_ranks_impl<true, true>(input, middle, output,
                                               directions, left, right);
    } else {
      partition_wavelet_ranks_impl<true, false>(input, middle, output,
                                                directions, left, right);
    }
  } else if (write_right) {
    partition_wavelet_ranks_impl<false, true>(input, middle, output, directions,
                                              left, right);
  } else {
    partition_wavelet_ranks_impl<false, false>(input, middle, output,
                                               directions, left, right);
  }

  if (left != expected_left || right != input.size()) {
    throw std::invalid_argument(
        "Wavelet partition does not match the supplied symbol counts");
  }
}

}  // namespace pixie::detail
