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

#if defined(__aarch64__) && defined(__ARM_NEON)
#define PIXIE_WAVELET_NEON_SUPPORT
#include <arm_neon.h>
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

struct WaveletCompaction16Tables {
  std::array<std::array<std::uint8_t, 16>, 256> low{};
  std::array<std::array<std::uint8_t, 32>, 256> high{};
};

// Adapted from PivCo's p16rev partition tables. One shuffle places the left
// ranks forward at the front and the right ranks in reverse at the back. A
// fixed reverse shuffle then recovers the stable right partition.
consteval WaveletCompaction16Tables make_wavelet_compaction_16_tables() {
  WaveletCompaction16Tables tables;
  for (std::size_t mask = 0; mask < 256; ++mask) {
    std::size_t left = 0;
    std::size_t right = 15;
    for (std::size_t index = 0; index < 8; ++index) {
      if ((mask & (std::size_t{1} << index)) == 0) {
        tables.low[mask][left++] = static_cast<std::uint8_t>(index);
      } else {
        tables.low[mask][right--] = static_cast<std::uint8_t>(index);
      }
    }

    left = 8;
    right = 15;
    for (std::size_t index = 0; index < 8; ++index) {
      if ((mask & (std::size_t{1} << index)) == 0) {
        tables.high[mask][left++] = static_cast<std::uint8_t>(index + 8);
      } else {
        tables.high[mask][right--] = static_cast<std::uint8_t>(index + 8);
      }
    }
  }
  return tables;
}

alignas(64) inline constexpr auto kWaveletCompaction16Tables =
    make_wavelet_compaction_16_tables();

inline std::uint16_t wavelet_direction_mask_16(__m128i ranks,
                                               __m128i first_right) {
  const __m128i right =
      _mm_cmpeq_epi8(_mm_min_epu8(ranks, first_right), first_right);
  return static_cast<std::uint16_t>(_mm_movemask_epi8(right));
}

inline __m128i compact_wavelet_16(__m128i ranks, std::uint16_t right_mask) {
  const std::uint8_t low_mask = static_cast<std::uint8_t>(right_mask);
  const std::uint8_t high_mask = static_cast<std::uint8_t>(right_mask >> 8);
  const unsigned low_right = std::popcount(low_mask);
  const __m128i low = _mm_load_si128(reinterpret_cast<const __m128i*>(
      kWaveletCompaction16Tables.low[low_mask].data()));
  const __m128i high = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
      kWaveletCompaction16Tables.high[high_mask].data() + low_right));
  return _mm_shuffle_epi8(ranks, _mm_or_si128(low, high));
}

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

#if defined(PIXIE_WAVELET_NEON_SUPPORT)

struct NeonWaveletCompactionTables {
  std::array<std::array<std::uint8_t, 16>, 256> low{};
  std::array<std::array<std::uint8_t, 32>, 256> high{};
};

consteval NeonWaveletCompactionTables make_neon_wavelet_compaction_tables() {
  NeonWaveletCompactionTables tables;
  for (std::size_t mask = 0; mask < 256; ++mask) {
    std::size_t left = 0;
    std::size_t right = 15;
    for (std::size_t index = 0; index < 8; ++index) {
      if ((mask & (std::size_t{1} << index)) == 0) {
        tables.low[mask][left++] = static_cast<std::uint8_t>(index);
      } else {
        tables.low[mask][right--] = static_cast<std::uint8_t>(index);
      }
    }
    left = 8;
    right = 15;
    for (std::size_t index = 0; index < 8; ++index) {
      if ((mask & (std::size_t{1} << index)) == 0) {
        tables.high[mask][left++] = static_cast<std::uint8_t>(index + 8);
      } else {
        tables.high[mask][right--] = static_cast<std::uint8_t>(index + 8);
      }
    }
  }
  return tables;
}

alignas(64) inline constexpr auto kNeonWaveletCompactionTables =
    make_neon_wavelet_compaction_tables();

inline std::uint16_t wavelet_direction_mask_16_neon(uint8x16_t ranks,
                                                    uint8x16_t first_right) {
  static constexpr std::array<std::uint8_t, 16> kBitWeights = {
      1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128};
  const uint8x16_t bits =
      vandq_u8(vcgeq_u8(ranks, first_right), vld1q_u8(kBitWeights.data()));
  const std::uint8_t low = vaddv_u8(vget_low_u8(bits));
  const std::uint8_t high = vaddv_u8(vget_high_u8(bits));
  return static_cast<std::uint16_t>(low | (std::uint16_t{high} << 8));
}

inline uint8x16_t compact_wavelet_16_neon(uint8x16_t ranks,
                                          std::uint16_t right_mask) {
  const std::uint8_t low_mask = static_cast<std::uint8_t>(right_mask);
  const std::uint8_t high_mask = static_cast<std::uint8_t>(right_mask >> 8);
  const unsigned low_right = std::popcount(low_mask);
  const uint8x16_t low =
      vld1q_u8(kNeonWaveletCompactionTables.low[low_mask].data());
  const uint8x16_t high =
      vld1q_u8(kNeonWaveletCompactionTables.high[high_mask].data() + low_right);
  return vqtbl1q_u8(ranks, vorrq_u8(low, high));
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

/** Convert byte symbols to in-order leaf ranks. */
inline void map_wavelet_byte_ranks(
    std::span<const std::uint8_t> symbols,
    std::span<std::uint8_t> ranks,
    const std::array<std::uint8_t, 256>& symbol_to_rank,
    const std::array<std::uint16_t, 256>& symbol_to_high_rank) {
  if (symbols.size() != ranks.size()) {
    throw std::invalid_argument("Invalid wavelet rank mapping buffers");
  }
  std::size_t offset = 0;
#if defined(PIXIE_WAVELET_NEON_SUPPORT)
  if (symbols.size() >= 20) {
    uint8x16x4_t table0;
    uint8x16x4_t table1;
    uint8x16x4_t table2;
    uint8x16x4_t table3;
    for (std::size_t lane = 0; lane < 4; ++lane) {
      table0.val[lane] = vld1q_u8(symbol_to_rank.data() + lane * 16);
      table1.val[lane] = vld1q_u8(symbol_to_rank.data() + 64 + lane * 16);
      table2.val[lane] = vld1q_u8(symbol_to_rank.data() + 128 + lane * 16);
      table3.val[lane] = vld1q_u8(symbol_to_rank.data() + 192 + lane * 16);
    }
    const uint8x16_t offset64 = vdupq_n_u8(64);
    const uint8x16_t offset128 = vdupq_n_u8(128);
    const uint8x16_t offset192 = vdupq_n_u8(192);
    for (; offset + 20 <= symbols.size(); offset += 20) {
      const uint8x16_t input = vld1q_u8(symbols.data() + offset);
      std::uint32_t tail;
      std::memcpy(&tail, symbols.data() + offset + 16, sizeof(tail));
      uint8x16_t mapped = vqtbl4q_u8(table0, input);
      const unsigned rank0 = symbol_to_rank[static_cast<std::uint8_t>(tail)];
      mapped = vqtbx4q_u8(mapped, table1, vsubq_u8(input, offset64));
      const unsigned rank1 =
          symbol_to_rank[static_cast<std::uint8_t>(tail >> 8)];
      mapped = vqtbx4q_u8(mapped, table2, vsubq_u8(input, offset128));
      const unsigned rank2 =
          symbol_to_rank[static_cast<std::uint8_t>(tail >> 16)];
      mapped = vqtbx4q_u8(mapped, table3, vsubq_u8(input, offset192));
      const unsigned rank3 =
          symbol_to_rank[static_cast<std::uint8_t>(tail >> 24)];
      vst1q_u8(ranks.data() + offset, mapped);
      const std::uint32_t mapped_tail =
          rank0 | (rank1 << 8) | (rank2 << 16) | (rank3 << 24);
      std::memcpy(ranks.data() + offset + 16, &mapped_tail,
                  sizeof(mapped_tail));
    }
  }
  (void)symbol_to_high_rank;
#else
  if constexpr (std::endian::native == std::endian::little) {
    for (; offset + 16 <= symbols.size(); offset += 16) {
      std::uint64_t low_symbols;
      std::uint64_t high_symbols;
      std::memcpy(&low_symbols, symbols.data() + offset, sizeof(low_symbols));
      std::memcpy(&high_symbols, symbols.data() + offset + 8,
                  sizeof(high_symbols));

      std::array<std::uint16_t, 8> pairs;
      for (std::size_t pair = 0; pair < 4; ++pair) {
        pairs[pair] =
            symbol_to_rank[static_cast<std::uint8_t>(low_symbols)] +
            symbol_to_high_rank[static_cast<std::uint8_t>(low_symbols >> 8)];
        low_symbols >>= 16;
        pairs[pair + 4] =
            symbol_to_rank[static_cast<std::uint8_t>(high_symbols)] +
            symbol_to_high_rank[static_cast<std::uint8_t>(high_symbols >> 8)];
        high_symbols >>= 16;
      }
      std::memcpy(ranks.data() + offset, pairs.data(), 16);
    }
  }
#endif
  for (; offset < symbols.size(); ++offset) {
    ranks[offset] = symbol_to_rank[symbols[offset]];
  }
}

/** Convert one byte block to in-order leaf ranks in place. */
inline void map_wavelet_byte_ranks(
    std::span<std::uint8_t> symbols,
    const std::array<std::uint8_t, 256>& symbol_to_rank,
    const std::array<std::uint16_t, 256>& symbol_to_high_rank) {
  map_wavelet_byte_ranks(std::span<const std::uint8_t>(symbols), symbols,
                         symbol_to_rank, symbol_to_high_rank);
}

template <bool WriteLeft, bool WriteRight>
std::size_t partition_wavelet_byte_block_scalar(
    std::span<std::uint8_t> ranks,
    std::uint8_t first_right,
    std::span<std::uint8_t> right_output,
    PackedBitBuilder& directions) {
  std::size_t left = 0;
  std::size_t right = 0;
  for (std::size_t offset = 0; offset < ranks.size(); offset += 64) {
    const std::size_t width = std::min<std::size_t>(64, ranks.size() - offset);
    std::uint64_t right_mask = 0;
    for (std::size_t index = 0; index < width; ++index) {
      const std::uint8_t rank = ranks[offset + index];
      const bool goes_right = rank >= first_right;
      right_mask |= static_cast<std::uint64_t>(goes_right) << index;
      if (goes_right) {
        if constexpr (WriteRight) {
          right_output[right] = rank;
        }
        ++right;
      } else {
        if constexpr (WriteLeft) {
          ranks[left] = rank;
        }
        ++left;
      }
    }
    directions.write_bits(right_mask, width);
  }
  return left;
}

#if defined(PIXIE_WAVELET_NEON_SUPPORT)

// PivCo's p16rev idea expressed with AArch64 TBL. The two halves are loaded
// before either in-place store, and the template removes dead leaf scatters.
template <bool WriteLeft, bool WriteRight>
std::size_t partition_wavelet_byte_block_neon(
    std::span<std::uint8_t> ranks,
    std::uint8_t first_right,
    std::span<std::uint8_t> right_output,
    PackedBitBuilder& directions) {
  static constexpr std::array<std::uint8_t, 16> kReverse = {
      15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  const uint8x16_t reverse = vld1q_u8(kReverse.data());
  const uint8x16_t threshold = vdupq_n_u8(first_right);
  std::size_t left = 0;
  std::size_t right = 0;
  std::size_t offset = 0;

  for (; offset + 32 <= ranks.size(); offset += 32) {
    const uint8x16_t low = vld1q_u8(ranks.data() + offset);
    const uint8x16_t high = vld1q_u8(ranks.data() + offset + 16);
    const std::uint16_t low_mask =
        wavelet_direction_mask_16_neon(low, threshold);
    const std::uint16_t high_mask =
        wavelet_direction_mask_16_neon(high, threshold);
    const std::uint32_t right_mask =
        low_mask | (static_cast<std::uint32_t>(high_mask) << 16);
    const unsigned low_right = std::popcount(low_mask);
    const unsigned low_left = 16 - low_right;
    const unsigned right_count = std::popcount(right_mask);
    directions.write_bits(right_mask, 32);

    if constexpr (WriteLeft || WriteRight) {
      const uint8x16_t compacted_low = compact_wavelet_16_neon(low, low_mask);
      const uint8x16_t compacted_high =
          compact_wavelet_16_neon(high, high_mask);
      if constexpr (WriteLeft) {
        vst1q_u8(ranks.data() + left, compacted_low);
        vst1q_u8(ranks.data() + left + low_left, compacted_high);
      }
      if constexpr (WriteRight) {
        vst1q_u8(right_output.data() + right,
                 vqtbl1q_u8(compacted_low, reverse));
        vst1q_u8(right_output.data() + right + low_right,
                 vqtbl1q_u8(compacted_high, reverse));
      }
    }
    left += 32 - right_count;
    right += right_count;
  }
  for (; offset + 16 <= ranks.size(); offset += 16) {
    const uint8x16_t values = vld1q_u8(ranks.data() + offset);
    const std::uint16_t right_mask =
        wavelet_direction_mask_16_neon(values, threshold);
    const unsigned right_count = std::popcount(right_mask);
    directions.write_bits(right_mask, 16);
    if constexpr (WriteLeft || WriteRight) {
      const uint8x16_t compacted = compact_wavelet_16_neon(values, right_mask);
      if constexpr (WriteLeft) {
        vst1q_u8(ranks.data() + left, compacted);
      }
      if constexpr (WriteRight) {
        vst1q_u8(right_output.data() + right, vqtbl1q_u8(compacted, reverse));
      }
    }
    left += 16 - right_count;
    right += right_count;
  }
  for (; offset < ranks.size(); ++offset) {
    const std::uint8_t rank = ranks[offset];
    const bool goes_right = rank >= first_right;
    directions.write_bit(goes_right);
    if (goes_right) {
      if constexpr (WriteRight) {
        right_output[right] = rank;
      }
      ++right;
    } else {
      if constexpr (WriteLeft) {
        ranks[left] = rank;
      }
      ++left;
    }
  }
  return left;
}

#endif

#if defined(PIXIE_AVX2_SUPPORT)

// PivCo's non-flat p16rev kernel adapted to append direction masks directly to
// a Pixie node stream. This deliberately keeps the full/right/none compile-time
// specializations so leaf children do not cause dead scatter traffic.
template <bool WriteLeft, bool WriteRight>
std::size_t partition_wavelet_byte_block_x86(
    std::span<std::uint8_t> ranks,
    std::uint8_t first_right,
    std::span<std::uint8_t> right_output,
    PackedBitBuilder& directions) {
  static constexpr std::array<std::uint8_t, 16> kReverse = {
      15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  const __m128i reverse =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(kReverse.data()));
  const __m128i threshold = _mm_set1_epi8(static_cast<char>(first_right));
  std::size_t left = 0;
  std::size_t right = 0;
  std::size_t offset = 0;

#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
  const __m512i threshold_512 =
      _mm512_set1_epi8(static_cast<char>(first_right - 1));
  for (; offset + 64 <= ranks.size(); offset += 64) {
    const __m512i values = _mm512_loadu_si512(ranks.data() + offset);
    const __mmask64 right_mask = _mm512_cmpgt_epu8_mask(values, threshold_512);
    const unsigned right_count = std::popcount(right_mask);
    directions.write_bits(right_mask, 64);
    if constexpr (WriteLeft) {
      _mm512_mask_compressstoreu_epi8(ranks.data() + left, ~right_mask, values);
    }
    if constexpr (WriteRight) {
      _mm512_mask_compressstoreu_epi8(right_output.data() + right, right_mask,
                                      values);
    }
    left += 64 - right_count;
    right += right_count;
  }
#endif

  for (; offset + 32 <= ranks.size(); offset += 32) {
    const __m128i low = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(ranks.data() + offset));
    const __m128i high = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(ranks.data() + offset + 16));
    const std::uint16_t low_mask = wavelet_direction_mask_16(low, threshold);
    const std::uint16_t high_mask = wavelet_direction_mask_16(high, threshold);
    const std::uint32_t right_mask =
        low_mask | (static_cast<std::uint32_t>(high_mask) << 16);
    const unsigned low_right = std::popcount(low_mask);
    const unsigned low_left = 16 - low_right;
    const unsigned right_count = std::popcount(right_mask);
    directions.write_bits(right_mask, 32);

    if constexpr (WriteLeft || WriteRight) {
      const __m128i compacted_low = compact_wavelet_16(low, low_mask);
      const __m128i compacted_high = compact_wavelet_16(high, high_mask);
      if constexpr (WriteLeft) {
        _mm_storeu_si128(reinterpret_cast<__m128i*>(ranks.data() + left),
                         compacted_low);
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(ranks.data() + left + low_left),
            compacted_high);
      }
      if constexpr (WriteRight) {
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(right_output.data() + right),
            _mm_shuffle_epi8(compacted_low, reverse));
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(right_output.data() + right + low_right),
            _mm_shuffle_epi8(compacted_high, reverse));
      }
    }
    left += 32 - right_count;
    right += right_count;
  }

  for (; offset + 16 <= ranks.size(); offset += 16) {
    const __m128i values = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(ranks.data() + offset));
    const std::uint16_t right_mask =
        wavelet_direction_mask_16(values, threshold);
    const unsigned right_count = std::popcount(right_mask);
    directions.write_bits(right_mask, 16);
    if constexpr (WriteLeft || WriteRight) {
      const __m128i compacted = compact_wavelet_16(values, right_mask);
      if constexpr (WriteLeft) {
        _mm_storeu_si128(reinterpret_cast<__m128i*>(ranks.data() + left),
                         compacted);
      }
      if constexpr (WriteRight) {
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(right_output.data() + right),
            _mm_shuffle_epi8(compacted, reverse));
      }
    }
    left += 16 - right_count;
    right += right_count;
  }

  for (; offset < ranks.size(); ++offset) {
    const std::uint8_t rank = ranks[offset];
    const bool goes_right = rank >= first_right;
    directions.write_bit(goes_right);
    if (goes_right) {
      if constexpr (WriteRight) {
        right_output[right] = rank;
      }
      ++right;
    } else {
      if constexpr (WriteLeft) {
        ranks[left] = rank;
      }
      ++left;
    }
  }
  return left;
}

#endif

/**
 * @brief Emit a cache-sized byte block's directions and live child ranks.
 * @details The left child is compacted in place and the right child is written
 * to disjoint scratch. Node streams are appended bit-contiguously across
 * blocks; PivCo's wire framing and per-block bitmap padding are not retained.
 */
inline std::size_t partition_wavelet_byte_block(
    std::span<std::uint8_t> ranks,
    std::uint8_t first_right,
    std::span<std::uint8_t> right_output,
    bool write_left,
    bool write_right,
    PackedBitBuilder& directions) {
  if (write_right && right_output.size() < ranks.size()) {
    throw std::invalid_argument("Invalid wavelet block partition buffers");
  }

  const auto partition = [&]<bool WriteLeft, bool WriteRight>() {
#if defined(PIXIE_AVX2_SUPPORT)
    return partition_wavelet_byte_block_x86<WriteLeft, WriteRight>(
        ranks, first_right, right_output, directions);
#elif defined(PIXIE_WAVELET_NEON_SUPPORT)
    return partition_wavelet_byte_block_neon<WriteLeft, WriteRight>(
        ranks, first_right, right_output, directions);
#else
    return partition_wavelet_byte_block_scalar<WriteLeft, WriteRight>(
        ranks, first_right, right_output, directions);
#endif
  };
  if (write_left) {
    return write_right ? partition.template operator()<true, true>()
                       : partition.template operator()<true, false>();
  }
  return write_right ? partition.template operator()<false, true>()
                     : partition.template operator()<false, false>();
}

}  // namespace pixie::detail

#if defined(PIXIE_WAVELET_NEON_SUPPORT)
#undef PIXIE_WAVELET_NEON_SUPPORT
#endif
