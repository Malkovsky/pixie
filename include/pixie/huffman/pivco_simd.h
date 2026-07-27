#pragma once

/**
 * @file pivco_simd.h
 * @brief SIMD kernels for PivCo-Huffman encode/decode primitives.
 *
 * Implements the two core PivCo-Huffman primitives from the paper with
 * architecture-specific kernels and a portable scalar fallback:
 *
 * - `pivco_merge_decode` (bottom-up decode, `merge_vec_vec`): interleaves two
 *   dense child symbol streams under a 1-bit-per-symbol routing bitmap into a
 *   single dense output. A set bitmap bit selects the right child, a clear bit
 *   the left child.
 * - `pivco_partition_encode` (top-down encode partition): the inverse — given
 *   per-symbol direction bits, splits a dense input into left/right child
 *   streams and writes the routing bitmap.
 * - `pivco_flat_decode` (bottom-up `merge_flat_D`): unpacks fixed-width leaf
 *   indices and translates them through a flat-subtree symbol table.
 *
 * Bitmap convention (matches `pivco_huffman.h`): bits are packed LSB-first
 * into 64-bit words; bit 0 of word 0 is the first symbol.
 *
 * Backends, selected by compiler feature macros (`-march`):
 * - AVX-512 VBMI + VBMI2: byte-granular `vpermb`, `vpmultishiftqb`,
 *   `vpexpandb`, and `vpcompressb` operating on 64 symbols per iteration.
 *   These extensions implement vector direction lookup, partition, and the
 *   full `merge_flat_D` family without changing the wire format.
 * - AVX2 + SSE4.1: 8-byte-at-a-time `_mm_shuffle_epi8` with 256-entry lookup
 *   tables, mirroring the paper's NEON `vqtbl1q_u8` technique.
 * - Scalar: the original bit-by-bit / symbol-by-symbol loops.
 *
 * Every SIMD backend produces byte-identical output to the scalar fallback
 * (the operation is a deterministic permutation). The existing round-trip tests
 * under release and AddressSanitizer serve as the differential oracle for every
 * backend the test matrix builds and runs; the AVX-512 backend is verified only
 * when built on VBMI/VBMI2-capable hardware, since there is no runtime
 * dispatch.
 */

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(__AVX2__) || defined(__AVX512VBMI2__)
#include <immintrin.h>
#endif

namespace pixie::pivco_simd_detail {

/// @brief Bitmap bit-ordering convention used throughout this header.
/// @details Every routing bitmap packs one bit per symbol, LSB-first within
/// each
///          64-bit word: bit `k` is the direction bit of symbol `k`, where a
///          set bit selects the right child and a clear bit the left child. All
///          mask builders (`build_dir_mask8`, the AVX-512 scalar loop), LUT
///          tables, merge kernels, and partition tails MUST agree on this
///          convention — flipping it in one backend only would make the
///          serialized bitmap backend-dependent and silently break
///          cross-`-march` encode/decode.

// ---------------------------------------------------------------------------
// Scalar helpers (always compiled; used as tail + fallback).
// ---------------------------------------------------------------------------

/// @brief Build an 8-bit routing mask for symbols @p s8 (one byte per symbol)
///        under the LSB-first convention (bit k = dir[s8[k]]).
/// @details Shared source of truth for the 8-symbol direction mask; the AVX2
///          partition kernel and (transitively) the AVX-512 tail both build
///          their per-group mask here.
inline std::uint8_t build_dir_mask8(const std::uint8_t* s8,
                                    const std::uint8_t* dir) noexcept {
  return static_cast<std::uint8_t>((dir[s8[0]] << 0) | (dir[s8[1]] << 1) |
                                   (dir[s8[2]] << 2) | (dir[s8[3]] << 3) |
                                   (dir[s8[4]] << 4) | (dir[s8[5]] << 5) |
                                   (dir[s8[6]] << 6) | (dir[s8[7]] << 7));
}

/// @brief Scalar merge tail: process the remaining `< 8` symbols bit-by-bit.
template <bool kLeftConstant, bool kRightConstant>
inline void merge_scalar_tail(std::uint8_t* dst,
                              std::size_t& out_i,
                              const std::uint8_t* left,
                              std::uint8_t left_symbol,
                              std::size_t& c_left,
                              const std::uint8_t* right,
                              std::uint8_t right_symbol,
                              std::size_t& c_right,
                              const std::uint64_t* bits,
                              std::size_t i,
                              std::size_t count) noexcept {
  for (; i < count; ++i) {
    const std::uint8_t bit =
        static_cast<std::uint8_t>((bits[i / 64] >> (i % 64)) & 1u);
    if (bit) {
      if constexpr (kRightConstant) {
        dst[out_i++] = right_symbol;
      } else {
        dst[out_i++] = right[c_right];
      }
      ++c_right;
    } else {
      if constexpr (kLeftConstant) {
        dst[out_i++] = left_symbol;
      } else {
        dst[out_i++] = left[c_left];
      }
      ++c_left;
    }
  }
}

/// @brief Write only the routing bitmap for a scalar partition tail.
inline void partition_bitmap_tail(std::uint64_t* bits,
                                  const std::uint8_t* src,
                                  const std::uint8_t* dir,
                                  std::size_t start_i,
                                  std::size_t weight) noexcept {
  std::size_t word_idx = start_i / 64;
  std::size_t bit_pos = start_i % 64;
  for (std::size_t i = start_i; i < weight; ++i) {
    if (dir[src[i]]) {
      bits[word_idx] |= (1ull << bit_pos);
    }
    if (++bit_pos == 64) {
      bit_pos = 0;
      ++word_idx;
    }
  }
}

/// @brief Scalar partition tail: route the remaining symbols (fewer than the
///        vector width) one at a time, writing the bitmap bit and appending to
///        the left/right child stream.
/// @details Shared by the AVX2 kernel tail, the AVX-512 kernel tail, and the
///          scalar fallback, so the LSB-first `|=`-into-pre-zeroed-arena
///          convention lives in exactly one place. @p bits is the bitmap base
///          (word array); @p start_bit is the absolute bit offset where the
///          tail begins (must be a multiple of the full-group width so it is
///          already byte-aligned within the word packing).
template <bool kWriteLeft, bool kWriteRight>
inline void partition_scalar_tail(std::uint8_t* left_dst,
                                  std::uint8_t* right_dst,
                                  std::uint64_t* bits,
                                  const std::uint8_t* src,
                                  const std::uint8_t* dir,
                                  std::size_t& c_left,
                                  std::size_t& c_right,
                                  std::size_t start_i,
                                  std::size_t weight) noexcept {
  std::size_t i = start_i;
  std::size_t word_idx = i / 64;
  std::size_t bit_pos = i % 64;
  for (; i < weight; ++i) {
    const std::uint8_t s = src[i];
    if (dir[s]) {
      bits[word_idx] |= (1ull << bit_pos);
      if constexpr (kWriteRight) {
        right_dst[c_right] = s;
      }
      ++c_right;
    } else {
      if constexpr (kWriteLeft) {
        left_dst[c_left] = s;
      }
      ++c_left;
    }
    if (++bit_pos == 64) {
      bit_pos = 0;
      ++word_idx;
    }
  }
}

// ===========================================================================
// AVX2 backend (SSE4.1 byte-shuffle lookup tables; AVX2 implies SSE4.1).
// ===========================================================================
#if defined(__AVX2__) && defined(__SSE4_1__)

/// @brief 256-entry × 16-byte shuffle-control table, 16-byte aligned per entry.
/// @details Stored as plain bytes (not `__m128i`) to avoid dropping the
///          `__may_alias__` attribute that `std::array<__m128i>` would warn
///          about; entries are loaded as aligned vectors at use time.
struct ShufTable {
  alignas(16) std::uint8_t entry[256][16];
};

/// @brief Load shuffle entry @p m as an aligned 128-bit vector.
inline __m128i load_shuf(const ShufTable& t, std::uint8_t m) noexcept {
  return _mm_load_si128(
      reinterpret_cast<const __m128i*>(&t.entry[m]));  // NOLINT
}

/// @brief 256-entry shuffle table for the 8-byte merge.
/// @details For an 8-bit mask `m` (bit i = 1 -> right), entry `m` is a
///          16-byte shuffle control that, applied to a 128-bit vector holding
///          `left[0..7]` in its low qword and `right[0..7]` in its high qword,
///          produces the merged 8-byte output. Mirrors the paper's NEON
///          `expand_tab`.
inline const ShufTable& merge_lut() {
  static const ShufTable lut = [] {
    ShufTable t{};
    for (int m = 0; m < 256; ++m) {
      int nl = 0;
      int nr = 0;
      for (int i = 0; i < 8; ++i) {
        if ((m >> i) & 1) {
          t.entry[m][i] = static_cast<std::uint8_t>(8 + nr);
          ++nr;
        } else {
          t.entry[m][i] = static_cast<std::uint8_t>(nl);
          ++nl;
        }
      }
      for (int i = 8; i < 16; ++i) {
        t.entry[m][i] = 0x80;  // zeroed lanes (unused, stored as full 16 bytes)
      }
    }
    return t;
  }();
  return lut;
}

/// @brief 256-entry shuffle tables for the 8-byte compress (partition).
/// @details `compress_right[m]` gathers the elements at the 1-bit positions of
///          mask `m` (in order) to the front of the 128-bit result;
///          `compress_left[m]` does the same for the 0-bit positions. Trailing
///          lanes are zeroed. Mirrors the paper's NEON `compress_tab`.
inline const ShufTable& compress_right_lut() {
  static const ShufTable lut = [] {
    ShufTable t{};
    for (int m = 0; m < 256; ++m) {
      int p = 0;
      for (int i = 0; i < 8; ++i) {
        if ((m >> i) & 1) {
          t.entry[m][p++] = static_cast<std::uint8_t>(i);
        }
      }
      for (int i = p; i < 16; ++i) {
        t.entry[m][i] = 0x80;
      }
    }
    return t;
  }();
  return lut;
}

inline const ShufTable& compress_left_lut() {
  static const ShufTable lut = [] {
    ShufTable t{};
    for (int m = 0; m < 256; ++m) {
      int p = 0;
      for (int i = 0; i < 8; ++i) {
        if (((m >> i) & 1) == 0) {
          t.entry[m][p++] = static_cast<std::uint8_t>(i);
        }
      }
      for (int i = p; i < 16; ++i) {
        t.entry[m][i] = 0x80;
      }
    }
    return t;
  }();
  return lut;
}

/// @brief Merge two constant leaves 32 symbols at a time with AVX2.
/// @details Four routing bytes are widened to four qwords, each byte is
///          broadcast over its eight output lanes, and a fixed bit-position
///          vector expands the bitmap without a shuffle-table lookup.
inline void merge_decode_cst_cst_avx2(std::uint8_t* dst,
                                      std::uint8_t left_symbol,
                                      std::uint8_t right_symbol,
                                      const std::uint64_t* bits,
                                      std::size_t count) noexcept {
  const auto* bit_bytes = reinterpret_cast<const std::uint8_t*>(bits);
  const __m256i broadcast_control =
      _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0,
                       0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8);
  const __m256i bit_positions = _mm256_setr_epi8(
      1, 2, 4, 8, 16, 32, 64, static_cast<char>(0x80), 1, 2, 4, 8, 16, 32, 64,
      static_cast<char>(0x80), 1, 2, 4, 8, 16, 32, 64, static_cast<char>(0x80),
      1, 2, 4, 8, 16, 32, 64, static_cast<char>(0x80));
  const __m256i zero = _mm256_setzero_si256();
  const __m256i left = _mm256_set1_epi8(static_cast<char>(left_symbol));
  const __m256i delta =
      _mm256_set1_epi8(static_cast<char>(left_symbol ^ right_symbol));

  std::size_t i = 0;
  for (; i + 32 <= count; i += 32) {
    std::uint32_t routing = 0;
    std::memcpy(&routing, bit_bytes + i / 8, sizeof(routing));
    const __m128i routing_bytes = _mm_cvtsi32_si128(static_cast<int>(routing));
    const __m256i routing_qwords = _mm256_cvtepu8_epi64(routing_bytes);
    const __m256i broadcast =
        _mm256_shuffle_epi8(routing_qwords, broadcast_control);
    const __m256i clear_bits =
        _mm256_cmpeq_epi8(_mm256_and_si256(broadcast, bit_positions), zero);
    const __m256i symbols =
        _mm256_xor_si256(left, _mm256_andnot_si256(clear_bits, delta));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + i), symbols);
  }
  for (; i < count; ++i) {
    dst[i] =
        ((bits[i / 64] >> (i % 64)) & 1u) != 0 ? right_symbol : left_symbol;
  }
}

template <bool kLeftConstant, bool kRightConstant>
inline void merge_decode_avx2(std::uint8_t* dst,
                              const std::uint8_t* left,
                              std::uint8_t left_symbol,
                              const std::uint8_t* right,
                              std::uint8_t right_symbol,
                              const std::uint64_t* bits,
                              std::size_t count) noexcept {
  if constexpr (kLeftConstant && kRightConstant) {
    merge_decode_cst_cst_avx2(dst, left_symbol, right_symbol, bits, count);
    return;
  }

  const auto& mlut = merge_lut();
  const auto* bit_bytes = reinterpret_cast<const std::uint8_t*>(bits);
  std::size_t out_i = 0;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  std::size_t i = 0;
  for (; i + 8 <= count; i += 8) {
    // AVX2 is little-endian on every supported x86 target. Eight-symbol
    // groups are byte-aligned in the LSB-first bitmap, so loading the routing
    // byte directly avoids a word-index/shift sequence in every merge group.
    const std::uint8_t m8 = bit_bytes[i / 8];
    __m128i left_values;
    if constexpr (kLeftConstant) {
      left_values = _mm_set1_epi8(static_cast<char>(left_symbol));
    } else {
      left_values =
          _mm_loadl_epi64(reinterpret_cast<const __m128i*>(left + c_left));
    }
    __m128i right_values;
    if constexpr (kRightConstant) {
      right_values = _mm_set1_epi8(static_cast<char>(right_symbol));
    } else {
      right_values =
          _mm_loadl_epi64(reinterpret_cast<const __m128i*>(right + c_right));
    }
    const __m128i combined = _mm_unpacklo_epi64(left_values, right_values);
    const __m128i out = _mm_shuffle_epi8(combined, load_shuf(mlut, m8));
    _mm_storel_epi64(reinterpret_cast<__m128i*>(dst + out_i), out);
    out_i += 8;
    const int nr = std::popcount(static_cast<unsigned>(m8));
    c_right += static_cast<std::size_t>(nr);
    c_left += static_cast<std::size_t>(8 - nr);
  }
  merge_scalar_tail<kLeftConstant, kRightConstant>(
      dst, out_i, left, left_symbol, c_left, right, right_symbol, c_right, bits,
      i, count);
}

/// @brief Store @p v (low @p n bytes valid) at @p dst without writing past
///        @p dst + @p n.
/// @details The codec partitions into a ping-pong workspace where a node's left
///          and right child regions are contiguous, and a node's right region
///          is immediately followed by its sibling's region (already written by
///          the parent, consumed later in pre-order). A full-width vector store
///          would clobber the start of the adjacent region, so when the valid
///          byte count is less than the vector width the store is masked to
///          exactly
///          @p n bytes.
inline void store_masked_epi64(std::uint8_t* dst, __m128i v, int n) noexcept {
  if (n >= 8) {
    _mm_storel_epi64(reinterpret_cast<__m128i*>(dst), v);
  } else if (n > 0) {
    alignas(8) std::uint8_t tmp[8];
    _mm_storel_epi64(reinterpret_cast<__m128i*>(tmp), v);
    std::memcpy(dst, tmp, static_cast<std::size_t>(n));
  }
}

template <bool kWriteLeft, bool kWriteRight>
inline void partition_encode_avx2(std::uint8_t* left_dst,
                                  std::uint8_t* right_dst,
                                  std::uint64_t* bits,
                                  const std::uint8_t* src,
                                  const std::uint8_t* dir,
                                  std::size_t weight,
                                  std::size_t left_weight,
                                  std::size_t right_weight) noexcept {
  std::uint8_t* bits_bytes = reinterpret_cast<std::uint8_t*>(bits);

  const std::size_t full = weight / 8;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  std::size_t g = 0;
  for (; g < full; ++g) {
    const std::uint8_t* s8 = src + g * 8;
    const std::uint8_t m8 = build_dir_mask8(s8, dir);
    // Byte-aligned write: full groups always start on a byte boundary in the
    // LSB-first bitmap, so byte g holds exactly symbols [8g, 8g+8).
    bits_bytes[g] = m8;
    if constexpr (kWriteLeft || kWriteRight) {
      const int nr = std::popcount(static_cast<unsigned>(m8));
      const int nl = 8 - nr;
      const __m128i data =
          _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s8));
      if constexpr (kWriteRight) {
        const __m128i right_values =
            _mm_shuffle_epi8(data, load_shuf(compress_right_lut(), m8));
        store_masked_epi64(right_dst + c_right, right_values,
                           static_cast<int>(right_weight - c_right));
      }
      if constexpr (kWriteLeft) {
        const __m128i left_values =
            _mm_shuffle_epi8(data, load_shuf(compress_left_lut(), m8));
        store_masked_epi64(left_dst + c_left, left_values,
                           static_cast<int>(left_weight - c_left));
      }
      c_right += static_cast<std::size_t>(nr);
      c_left += static_cast<std::size_t>(nl);
    }
  }

  // Scalar tail for the trailing < 8 symbols (byte-unaligned bitmap bits).
  if constexpr (kWriteLeft || kWriteRight) {
    partition_scalar_tail<kWriteLeft, kWriteRight>(
        left_dst, right_dst, bits, src, dir, c_left, c_right, full * 8, weight);
  } else {
    partition_bitmap_tail(bits, src, dir, full * 8, weight);
  }
}

/// @brief AVX2 partition when one child is a constant leaf.
/// @tparam kLeftConstant True when the left child is the constant; the dense
///                       right stream is written. False writes the dense left
///                       stream and treats the right child as constant.
template <bool kLeftConstant>
inline void partition_encode_one_constant_avx2(
    std::uint8_t* output,
    std::uint64_t* bits,
    const std::uint8_t* src,
    std::uint8_t constant_symbol,
    std::size_t weight,
    std::size_t output_weight) noexcept {
  auto* bit_bytes = reinterpret_cast<std::uint8_t*>(bits);
  const __m128i constant = _mm_set1_epi8(static_cast<char>(constant_symbol));
  std::size_t output_count = 0;
  const std::size_t full = weight / 8;
  std::size_t group = 0;
  for (; group < full; ++group) {
    const std::uint8_t* s8 = src + group * 8;
    const __m128i data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s8));
    const std::uint8_t equal_mask = static_cast<std::uint8_t>(
        _mm_movemask_epi8(_mm_cmpeq_epi8(data, constant)));
    const std::uint8_t routing =
        kLeftConstant ? static_cast<std::uint8_t>(~equal_mask) : equal_mask;
    bit_bytes[group] = routing;

    const auto& lut =
        kLeftConstant ? compress_right_lut() : compress_left_lut();
    const __m128i packed = _mm_shuffle_epi8(data, load_shuf(lut, routing));
    store_masked_epi64(output + output_count, packed,
                       static_cast<int>(output_weight - output_count));
    const std::size_t right_count =
        std::popcount(static_cast<unsigned>(routing));
    output_count += kLeftConstant ? right_count : 8 - right_count;
  }

  std::size_t word_index = (group * 8) / 64;
  std::size_t bit_position = (group * 8) % 64;
  for (std::size_t i = group * 8; i < weight; ++i) {
    const bool is_constant = src[i] == constant_symbol;
    const bool goes_right = kLeftConstant ? !is_constant : is_constant;
    if (goes_right) {
      bits[word_index] |= std::uint64_t{1} << bit_position;
    }
    if (kLeftConstant ? goes_right : !goes_right) {
      output[output_count++] = src[i];
    }
    if (++bit_position == 64) {
      bit_position = 0;
      ++word_index;
    }
  }
}

/// @brief Decode a 2-bit flat subtree, 32 symbols per iteration.
inline void flat_decode_d2_avx2(std::uint8_t* dst,
                                const std::uint8_t* packed,
                                const std::uint8_t* table,
                                std::size_t count) noexcept {
  alignas(16) std::uint8_t table_lanes[16]{};
  std::memcpy(table_lanes, table, 4);
  const __m128i lookup =
      _mm_load_si128(reinterpret_cast<const __m128i*>(table_lanes));
  const __m128i mask = _mm_set1_epi8(3);

  std::size_t i = 0;
  for (; i + 32 <= count; i += 32) {
    const __m128i data =
        _mm_loadl_epi64(reinterpret_cast<const __m128i*>(packed + i / 4));
    const __m128i code0 = _mm_and_si128(data, mask);
    const __m128i code1 = _mm_and_si128(_mm_srli_epi16(data, 2), mask);
    const __m128i code2 = _mm_and_si128(_mm_srli_epi16(data, 4), mask);
    const __m128i code3 = _mm_and_si128(_mm_srli_epi16(data, 6), mask);
    const __m128i code01 = _mm_unpacklo_epi8(code0, code1);
    const __m128i code23 = _mm_unpacklo_epi8(code2, code3);
    const __m128i indices0 = _mm_unpacklo_epi16(code01, code23);
    const __m128i indices1 = _mm_unpackhi_epi16(code01, code23);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i),
                     _mm_shuffle_epi8(lookup, indices0));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i + 16),
                     _mm_shuffle_epi8(lookup, indices1));
  }
  for (; i < count; ++i) {
    dst[i] = table[(packed[i / 4] >> (2 * (i % 4))) & 3u];
  }
}

/// @brief Decode a 4-bit flat subtree, 32 symbols per iteration.
inline void flat_decode_d4_avx2(std::uint8_t* dst,
                                const std::uint8_t* packed,
                                const std::uint8_t* table,
                                std::size_t count) noexcept {
  const __m128i lookup =
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(table));
  const __m128i mask = _mm_set1_epi8(0x0f);

  std::size_t i = 0;
  for (; i + 32 <= count; i += 32) {
    const __m128i data =
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(packed + i / 2));
    const __m128i low = _mm_and_si128(data, mask);
    const __m128i high = _mm_and_si128(_mm_srli_epi16(data, 4), mask);
    const __m128i low_symbols = _mm_shuffle_epi8(lookup, low);
    const __m128i high_symbols = _mm_shuffle_epi8(lookup, high);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i),
                     _mm_unpacklo_epi8(low_symbols, high_symbols));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i + 16),
                     _mm_unpackhi_epi8(low_symbols, high_symbols));
  }
  for (; i < count; ++i) {
    dst[i] = table[(packed[i / 2] >> (4 * (i % 2))) & 0x0fu];
  }
}

/// @brief Reverse the low @p width bits of one byte.
inline std::uint8_t reverse_low_bits_avx2(std::uint8_t value,
                                          std::uint8_t width) noexcept {
  value = static_cast<std::uint8_t>(((value & 0x55u) << 1) |
                                    ((value & 0xaau) >> 1));
  value = static_cast<std::uint8_t>(((value & 0x33u) << 2) |
                                    ((value & 0xccu) >> 2));
  value = static_cast<std::uint8_t>((value << 4) | (value >> 4));
  return value >> (8 - width);
}

/// @brief Decode a known full-alphabet depth-8 bit-reversal table.
/// @details A balanced tree built from symbols 0..255 stores
///          `table[reverse8(code)] = code`. AVX2 can apply that permutation
///          with two nibble `pshufb` operations, avoiding the arbitrary byte
///          lookup that AVX2 cannot otherwise vectorize efficiently.
inline void flat_decode_d8_bitreverse_unchecked_avx2(
    std::uint8_t* dst,
    const std::uint8_t* packed,
    const std::uint8_t* table,
    std::size_t count) noexcept {
  alignas(16) static constexpr std::uint8_t kReverseNibble[16] = {
      0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
  const __m128i reverse128 =
      _mm_load_si128(reinterpret_cast<const __m128i*>(kReverseNibble));
  const __m256i reverse = _mm256_broadcastsi128_si256(reverse128);
  const __m256i nibble_mask = _mm256_set1_epi8(0x0f);
  std::size_t i = 0;
  for (; i + 32 <= count; i += 32) {
    const __m256i codes =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(packed + i));
    const __m256i low = _mm256_and_si256(codes, nibble_mask);
    const __m256i high =
        _mm256_and_si256(_mm256_srli_epi16(codes, 4), nibble_mask);
    const __m256i reversed_low = _mm256_shuffle_epi8(reverse, low);
    const __m256i reversed_high = _mm256_shuffle_epi8(reverse, high);
    const __m256i symbols =
        _mm256_or_si256(_mm256_slli_epi16(reversed_low, 4), reversed_high);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + i), symbols);
  }
  for (; i < count; ++i) {
    dst[i] = table[packed[i]];
  }
}

/// @brief Decode the common full-alphabet depth-8 bit-reversal table.
/// @returns True when @p table has the reshaped tree's bit-reversal layout and
///          the output was produced; false for an arbitrary 256-byte table.
inline bool flat_decode_d8_bitreverse_avx2(std::uint8_t* dst,
                                           const std::uint8_t* packed,
                                           const std::uint8_t* table,
                                           std::size_t count) noexcept {
  for (std::size_t code = 0; code < 256; ++code) {
    if (table[code] !=
        reverse_low_bits_avx2(static_cast<std::uint8_t>(code), 8)) {
      return false;
    }
  }
  flat_decode_d8_bitreverse_unchecked_avx2(dst, packed, table, count);
  return true;
}

/// @brief Expand eight tightly packed D-bit codes into eight byte lanes.
/// @details Eight codes occupy exactly D bytes, so every eight-code group is
///          byte-aligned. BMI2 `pdep` inserts the padding bits between codes in
///          one instruction. The constexpr fallback is fully unrolled for
///          AVX2 CPUs that do not also expose BMI2.
template <std::uint8_t kDepth>
inline std::uint64_t expand_codes8_avx2(const std::uint8_t* packed) noexcept {
  static_assert(kDepth >= 3 && kDepth <= 8);
  std::uint64_t dense = 0;
  std::memcpy(&dense, packed, kDepth);
  if constexpr (kDepth == 8) {
    return dense;
  }

  constexpr std::uint64_t kCodeMask = (std::uint64_t{1} << kDepth) - 1;
#if defined(__BMI2__)
  constexpr std::uint64_t kByteLaneMask =
      kCodeMask * UINT64_C(0x0101010101010101);
  return _pdep_u64(dense, kByteLaneMask);
#else
  std::uint64_t expanded = 0;
  for (std::size_t lane = 0; lane < 8; ++lane) {
    expanded |= ((dense >> (lane * kDepth)) & kCodeMask) << (lane * 8);
  }
  return expanded;
#endif
}

/// @brief Reverse the low D bits independently in 16 byte lanes.
template <std::uint8_t kDepth>
inline __m128i reverse_low_bits16_avx2(__m128i values) noexcept {
  static_assert(kDepth >= 2 && kDepth <= 7);
  alignas(16) static constexpr std::uint8_t kReverseNibble[16] = {
      0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
  const __m128i reverse =
      _mm_load_si128(reinterpret_cast<const __m128i*>(kReverseNibble));
  const __m128i nibble_mask = _mm_set1_epi8(0x0f);
  const __m128i low = _mm_and_si128(values, nibble_mask);
  const __m128i high = _mm_and_si128(_mm_srli_epi16(values, 4), nibble_mask);
  const __m128i reversed =
      _mm_or_si128(_mm_slli_epi16(_mm_shuffle_epi8(reverse, low), 4),
                   _mm_shuffle_epi8(reverse, high));
  constexpr std::uint8_t kCodeMask = (std::uint8_t{1} << kDepth) - 1;
  return _mm_and_si128(_mm_srli_epi16(reversed, 8 - kDepth),
                       _mm_set1_epi8(static_cast<char>(kCodeMask)));
}

/// @brief Decode a known balanced low-alphabet table as bit reversal.
template <std::uint8_t kDepth>
inline void flat_decode_bitreverse_unchecked_avx2(std::uint8_t* dst,
                                                  const std::uint8_t* packed,
                                                  const std::uint8_t* table,
                                                  std::size_t count) noexcept {
  static_assert(kDepth >= 6 && kDepth <= 7);
  std::size_t i = 0;
  for (; i + 16 <= count; i += 16) {
    const std::size_t byte_offset = (i * kDepth) / 8;
    const std::uint64_t codes0 =
        expand_codes8_avx2<kDepth>(packed + byte_offset);
    const std::uint64_t codes1 =
        expand_codes8_avx2<kDepth>(packed + byte_offset + kDepth);
    const __m128i codes = _mm_set_epi64x(static_cast<long long>(codes1),
                                         static_cast<long long>(codes0));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i),
                     reverse_low_bits16_avx2<kDepth>(codes));
  }
  constexpr std::uint64_t kCodeMask = (std::uint64_t{1} << kDepth) - 1;
  const auto* words = reinterpret_cast<const std::uint64_t*>(packed);
  for (; i < count; ++i) {
    const std::size_t bit_position = i * kDepth;
    const std::size_t word_index = bit_position / 64;
    const std::size_t bit_offset = bit_position % 64;
    std::uint64_t code = words[word_index] >> bit_offset;
    if (bit_offset + kDepth > 64) {
      code |= words[word_index + 1] << (64 - bit_offset);
    }
    dst[i] = table[code & kCodeMask];
  }
}

/// @brief Decode a balanced low-alphabet flat table as bit reversal.
/// @returns True when the table has the expected permutation and was decoded.
template <std::uint8_t kDepth>
inline bool flat_decode_bitreverse_avx2(std::uint8_t* dst,
                                        const std::uint8_t* packed,
                                        const std::uint8_t* table,
                                        std::size_t count) noexcept {
  static_assert(kDepth >= 6 && kDepth <= 7);
  constexpr std::size_t kTableSize = std::size_t{1} << kDepth;
  for (std::size_t code = 0; code < kTableSize; ++code) {
    if (table[code] !=
        reverse_low_bits_avx2(static_cast<std::uint8_t>(code), kDepth)) {
      return false;
    }
  }
  flat_decode_bitreverse_unchecked_avx2<kDepth>(dst, packed, table, count);
  return true;
}

/// @brief Lookup byte indices in a power-of-two collection of 16-byte banks.
/// @details AVX2 has no arbitrary byte permutation across 256 entries. Each
///          leaf performs one `pshufb` using the low nibble; balanced
///          `blendv` levels select the bank with the remaining high bits.
///          The recursion is compile-time and emits straight-line SIMD code.
template <std::size_t kFirstBank, std::size_t kBankCount>
inline __m128i lookup_banked_bytes_avx2(const std::uint8_t* table,
                                        __m128i low_nibbles,
                                        __m128i codes) noexcept {
  static_assert(kBankCount > 0 && std::has_single_bit(kBankCount));
  if constexpr (kBankCount == 1) {
    const __m128i bank = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table + kFirstBank * 16));
    return _mm_shuffle_epi8(bank, low_nibbles);
  } else {
    constexpr std::size_t kHalf = kBankCount / 2;
    constexpr int kSelectBit = 3 + std::countr_zero(kBankCount);
    const __m128i low =
        lookup_banked_bytes_avx2<kFirstBank, kHalf>(table, low_nibbles, codes);
    const __m128i high = lookup_banked_bytes_avx2<kFirstBank + kHalf, kHalf>(
        table, low_nibbles, codes);
    const __m128i select = _mm_slli_epi16(codes, 7 - kSelectBit);
    return _mm_blendv_epi8(low, high, select);
  }
}

/// @brief Decode flat depths 3 and 5..7 in 16-symbol AVX2 batches.
/// @details The existing depth-2/depth-4 kernels unpack 32 symbols per batch
///          and remain faster for those widths. Other depths use BMI2 to
///          expand two groups of eight packed codes, followed by a banked
///          `pshufb` lookup. This avoids the former per-symbol shift, branch,
///          and dependent table load loop.
template <std::uint8_t kDepth>
inline void flat_decode_banked_avx2(std::uint8_t* dst,
                                    const std::uint8_t* packed,
                                    const std::uint8_t* table,
                                    std::size_t count) noexcept {
  static_assert(kDepth == 3 || (kDepth >= 5 && kDepth <= 7));
  constexpr std::size_t kTableSize = std::size_t{1} << kDepth;
  constexpr std::size_t kBankCount = kTableSize <= 16 ? 1 : kTableSize / 16;
  const __m128i nibble_mask = _mm_set1_epi8(0x0f);

  // A depth-3 table has only eight valid bytes, while the banked lookup loads
  // a complete vector. Pad that one small table locally; depths >=5 already
  // contain one or more complete 16-byte banks.
  alignas(16) std::uint8_t padded_depth3[16]{};
  const std::uint8_t* lookup_table = table;
  if constexpr (kDepth == 3) {
    std::memcpy(padded_depth3, table, kTableSize);
    lookup_table = padded_depth3;
  }

  std::size_t i = 0;
  for (; i + 16 <= count; i += 16) {
    const std::size_t byte_offset = (i * kDepth) / 8;
    const std::uint64_t codes0 =
        expand_codes8_avx2<kDepth>(packed + byte_offset);
    const std::uint64_t codes1 =
        expand_codes8_avx2<kDepth>(packed + byte_offset + kDepth);
    const __m128i codes = _mm_set_epi64x(static_cast<long long>(codes1),
                                         static_cast<long long>(codes0));
    const __m128i low_nibbles = _mm_and_si128(codes, nibble_mask);
    const __m128i symbols = lookup_banked_bytes_avx2<0, kBankCount>(
        lookup_table, low_nibbles, codes);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i), symbols);
  }

  const auto* words = reinterpret_cast<const std::uint64_t*>(packed);
  constexpr std::uint64_t kCodeMask = (std::uint64_t{1} << kDepth) - 1;
  for (; i < count; ++i) {
    const std::size_t bit_position = i * kDepth;
    const std::size_t word_index = bit_position / 64;
    const std::size_t bit_offset = bit_position % 64;
    std::uint64_t code = words[word_index] >> bit_offset;
    if (bit_offset + kDepth > 64) {
      code |= words[word_index + 1] << (64 - bit_offset);
    }
    dst[i] = table[code & kCodeMask];
  }
}

#if defined(__BMI2__)
/// @brief Encode one depth-3 flat subtree with AVX2 gather and BMI2 packing.
/// @details AVX2 has no byte gather, so eight input symbols are widened to
///          dword indices and looked up together. The gathered codes are
///          narrowed back to bytes, then `pext` packs their low three bits.
///          Wider flat depths were benchmarked and deliberately retain their
///          faster scalar encoders. The caller provides three readable pad
///          bytes after the 256-entry code table because dword gathers overlap
///          adjacent byte entries.
inline void flat_encode_d3_avx2(std::uint64_t* bits,
                                const std::uint8_t* src,
                                const std::uint8_t* wire_code,
                                std::size_t count) noexcept {
  constexpr std::uint8_t kDepth = 3;
  constexpr std::uint64_t kByteLaneMask = UINT64_C(0x0707070707070707);
  auto* packed = reinterpret_cast<std::uint8_t*>(bits);
  const __m256i byte_mask = _mm256_set1_epi32(0xff);
  const __m128i zero = _mm_setzero_si128();

  std::size_t i = 0;
  for (; i + 8 <= count; i += 8) {
    const __m128i symbols =
        _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src + i));
    const __m256i indices = _mm256_cvtepu8_epi32(symbols);
    const __m256i gathered = _mm256_and_si256(
        _mm256_i32gather_epi32(reinterpret_cast<const int*>(wire_code), indices,
                               1),
        byte_mask);
    const __m128i gathered_low = _mm256_castsi256_si128(gathered);
    const __m128i gathered_high = _mm256_extracti128_si256(gathered, 1);
    const __m128i codes16 = _mm_packus_epi32(gathered_low, gathered_high);
    const __m128i codes8 = _mm_packus_epi16(codes16, zero);
    const std::uint64_t codes =
        static_cast<std::uint64_t>(_mm_cvtsi128_si64(codes8));
    const std::uint64_t dense = _pext_u64(codes, kByteLaneMask);
    std::memcpy(packed + (i * kDepth) / 8, &dense, kDepth);
  }

  std::size_t word_index = (i * kDepth) / 64;
  std::size_t bit_offset = (i * kDepth) % 64;
  for (; i < count; ++i) {
    const std::uint64_t code = wire_code[src[i]];
    bits[word_index] |= code << bit_offset;
    if (bit_offset + kDepth > 64) {
      bits[word_index + 1] |= code >> (64 - bit_offset);
    }
    bit_offset += kDepth;
    if (bit_offset >= 64) {
      bit_offset -= 64;
      ++word_index;
    }
  }
}

/// @brief Encode a known balanced low-alphabet table as bit reversal.
template <std::uint8_t kDepth>
inline void flat_encode_bitreverse_unchecked_avx2(std::uint64_t* bits,
                                                  const std::uint8_t* src,
                                                  const std::uint8_t* wire_code,
                                                  std::size_t count) noexcept {
  static_assert(kDepth >= 2 && kDepth <= 7);
  constexpr std::uint64_t kCodeMask = (std::uint64_t{1} << kDepth) - 1;
  constexpr std::uint64_t kByteLaneMask =
      kCodeMask * UINT64_C(0x0101010101010101);
  auto* packed = reinterpret_cast<std::uint8_t*>(bits);
  std::size_t i = 0;
  for (; i + 8 <= count; i += 8) {
    const __m128i symbols =
        _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src + i));
    const __m128i codes = reverse_low_bits16_avx2<kDepth>(symbols);
    const std::uint64_t byte_codes =
        static_cast<std::uint64_t>(_mm_cvtsi128_si64(codes));
    const std::uint64_t dense = _pext_u64(byte_codes, kByteLaneMask);
    std::memcpy(packed + (i * kDepth) / 8, &dense, kDepth);
  }

  std::size_t word_index = (i * kDepth) / 64;
  std::size_t bit_offset = (i * kDepth) % 64;
  for (; i < count; ++i) {
    const std::uint64_t code = wire_code[src[i]];
    bits[word_index] |= code << bit_offset;
    if (bit_offset + kDepth > 64) {
      bits[word_index + 1] |= code >> (64 - bit_offset);
    }
    bit_offset += kDepth;
    if (bit_offset >= 64) {
      bit_offset -= 64;
      ++word_index;
    }
  }
}

/// @brief Encode a balanced low-alphabet table as vectorized bit reversal.
/// @returns True when @p wire_code has the expected mapping and was encoded.
template <std::uint8_t kDepth>
inline bool flat_encode_bitreverse_avx2(std::uint64_t* bits,
                                        const std::uint8_t* src,
                                        const std::uint8_t* wire_code,
                                        std::size_t count) noexcept {
  static_assert(kDepth >= 2 && kDepth <= 7);
  constexpr std::size_t kTableSize = std::size_t{1} << kDepth;
  for (std::size_t symbol = 0; symbol < kTableSize; ++symbol) {
    if (wire_code[symbol] !=
        reverse_low_bits_avx2(static_cast<std::uint8_t>(symbol), kDepth)) {
      return false;
    }
  }
  flat_encode_bitreverse_unchecked_avx2<kDepth>(bits, src, wire_code, count);
  return true;
}

/// @brief Select checked or preclassified balanced-table encoding.
template <std::uint8_t kDepth>
inline bool flat_encode_bitreverse_maybe_avx2(std::uint64_t* bits,
                                              const std::uint8_t* src,
                                              const std::uint8_t* wire_code,
                                              std::size_t count,
                                              bool known_bitreverse) noexcept {
  if (known_bitreverse) {
    flat_encode_bitreverse_unchecked_avx2<kDepth>(bits, src, wire_code, count);
    return true;
  }
  return flat_encode_bitreverse_avx2<kDepth>(bits, src, wire_code, count);
}
#endif

#endif  // AVX2 + SSE4.1

// ===========================================================================
// AVX-512 VBMI2 backend (byte-granular expand/compress, 64 bytes/iteration).
// ===========================================================================
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)

#if defined(__AVX512VBMI__)

/// @brief Four-bank byte lookup used for 64-, 128-, and 256-entry tables.
/// @details `vpermb` addresses one 64-byte bank. The high two index bits select
///          the bank after all required bank permutations run in parallel.
template <std::uint8_t kDepth>
struct ByteLookupAvx512 {
  __m512i bank0;
  __m512i bank1;
  __m512i bank2;
  __m512i bank3;
};

/// @brief Load a `2^kDepth` byte table without reading beyond its allocation.
template <std::uint8_t kDepth>
inline ByteLookupAvx512<kDepth> load_byte_lookup_avx512(
    const std::uint8_t* table) noexcept {
  static_assert(kDepth >= 1 && kDepth <= 8);
  constexpr std::size_t kTableSize = std::size_t{1} << kDepth;
  constexpr __mmask64 kFirstBankMask =
      kTableSize >= 64 ? ~__mmask64{0}
                       : (__mmask64{1} << kTableSize) - __mmask64{1};
  const __m512i zero = _mm512_setzero_si512();
  ByteLookupAvx512<kDepth> lookup{
      _mm512_maskz_loadu_epi8(kFirstBankMask, table), zero, zero, zero};
  if constexpr (kDepth >= 7) {
    lookup.bank1 = _mm512_loadu_si512(table + 64);
  }
  if constexpr (kDepth == 8) {
    lookup.bank2 = _mm512_loadu_si512(table + 128);
    lookup.bank3 = _mm512_loadu_si512(table + 192);
  }
  return lookup;
}

/// @brief Translate 64 byte indices through a preloaded `2^kDepth` table.
template <std::uint8_t kDepth>
inline __m512i lookup_bytes_avx512(
    __m512i indices,
    const ByteLookupAvx512<kDepth>& lookup) noexcept {
  const __m512i low_indices = _mm512_and_si512(indices, _mm512_set1_epi8(0x3f));
  __m512i result = _mm512_permutexvar_epi8(low_indices, lookup.bank0);
  if constexpr (kDepth >= 7) {
    const __m512i banks =
        _mm512_and_si512(_mm512_srli_epi16(indices, 6), _mm512_set1_epi8(3));
    const __mmask64 bank1 = _mm512_cmpeq_epi8_mask(banks, _mm512_set1_epi8(1));
    result = _mm512_mask_mov_epi8(
        result, bank1, _mm512_permutexvar_epi8(low_indices, lookup.bank1));
    if constexpr (kDepth == 8) {
      const __mmask64 bank2 =
          _mm512_cmpeq_epi8_mask(banks, _mm512_set1_epi8(2));
      const __mmask64 bank3 =
          _mm512_cmpeq_epi8_mask(banks, _mm512_set1_epi8(3));
      result = _mm512_mask_mov_epi8(
          result, bank2, _mm512_permutexvar_epi8(low_indices, lookup.bank2));
      result = _mm512_mask_mov_epi8(
          result, bank3, _mm512_permutexvar_epi8(low_indices, lookup.bank3));
    }
  }
  return result;
}

/// @brief Decode one scalar tail after full 64-symbol flat-D groups.
template <std::uint8_t kDepth>
inline void flat_decode_tail_avx512(std::uint8_t* dst,
                                    const std::uint64_t* bits,
                                    const std::uint8_t* table,
                                    std::size_t start,
                                    std::size_t count) noexcept {
  constexpr std::uint64_t kMask = (std::uint64_t{1} << kDepth) - 1;
  for (std::size_t i = start; i < count; ++i) {
    const std::size_t bit_position = i * kDepth;
    const std::size_t word_index = bit_position / 64;
    const std::size_t bit_offset = bit_position % 64;
    std::uint64_t packed = bits[word_index] >> bit_offset;
    if (bit_offset + kDepth > 64) {
      packed |= bits[word_index + 1] << (64 - bit_offset);
    }
    dst[i] = table[packed & kMask];
  }
}

/// @brief AVX-512 implementation of the paper's `merge_flat_D` primitive.
/// @details Each 64-symbol group occupies exactly @p kDepth 64-bit words. The
///          words are permuted into eight 8-code windows, aligned with vector
///          shifts, unpacked by `vpmultishiftqb`, and translated by `vpermb`.
template <std::uint8_t kDepth>
inline void flat_decode_avx512(std::uint8_t* dst,
                               const std::uint64_t* bits,
                               const std::uint8_t* table,
                               std::size_t count) noexcept {
  static_assert(kDepth >= 2 && kDepth <= 8);
  const ByteLookupAvx512<kDepth> lookup =
      load_byte_lookup_avx512<kDepth>(table);

  alignas(64) static constexpr auto kWordIndices = [] {
    std::array<std::uint64_t, 8> values{};
    for (std::size_t lane = 0; lane < values.size(); ++lane) {
      values[lane] = (lane * 8 * kDepth) / 64;
    }
    return values;
  }();
  alignas(64) static constexpr auto kNextWordIndices = [] {
    std::array<std::uint64_t, 8> values{};
    for (std::size_t lane = 0; lane < values.size(); ++lane) {
      values[lane] = ((lane * 8 * kDepth) / 64 + 1) % 8;
    }
    return values;
  }();
  alignas(64) static constexpr auto kRightShifts = [] {
    std::array<std::uint64_t, 8> values{};
    for (std::size_t lane = 0; lane < values.size(); ++lane) {
      values[lane] = (lane * 8 * kDepth) % 64;
    }
    return values;
  }();
  alignas(64) static constexpr auto kLeftShifts = [] {
    std::array<std::uint64_t, 8> values{};
    for (std::size_t lane = 0; lane < values.size(); ++lane) {
      values[lane] = 64 - ((lane * 8 * kDepth) % 64);
    }
    return values;
  }();
  alignas(64) static constexpr auto kMultishiftControl = [] {
    std::array<std::uint8_t, 64> values{};
    for (std::size_t lane = 0; lane < 8; ++lane) {
      for (std::size_t code = 0; code < 8; ++code) {
        values[lane * 8 + code] = static_cast<std::uint8_t>(code * kDepth);
      }
    }
    return values;
  }();

  const __m512i word_indices = _mm512_load_si512(kWordIndices.data());
  const __m512i next_word_indices = _mm512_load_si512(kNextWordIndices.data());
  const __m512i right_shifts = _mm512_load_si512(kRightShifts.data());
  const __m512i left_shifts = _mm512_load_si512(kLeftShifts.data());
  const __m512i multishift_control =
      _mm512_load_si512(kMultishiftControl.data());
  const __m512i code_mask =
      _mm512_set1_epi8(static_cast<char>((1u << kDepth) - 1u));
  constexpr __mmask8 kPackedWordMask =
      static_cast<__mmask8>((1u << kDepth) - 1u);

  std::size_t i = 0;
  for (; i + 64 <= count; i += 64) {
    __m512i indices;
    if constexpr (kDepth == 8) {
      // Eight-bit codes are already byte indices; only lookup is required.
      indices =
          _mm512_loadu_si512(reinterpret_cast<const std::uint8_t*>(bits) + i);
    } else {
      const __m512i packed =
          _mm512_maskz_loadu_epi64(kPackedWordMask, bits + (i / 64) * kDepth);
      const __m512i current_words =
          _mm512_permutexvar_epi64(word_indices, packed);
      const __m512i next_words =
          _mm512_permutexvar_epi64(next_word_indices, packed);
      const __m512i windows =
          _mm512_or_si512(_mm512_srlv_epi64(current_words, right_shifts),
                          _mm512_sllv_epi64(next_words, left_shifts));
      indices = _mm512_and_si512(
          _mm512_multishift_epi64_epi8(multishift_control, windows), code_mask);
    }
    _mm512_storeu_si512(dst + i, lookup_bytes_avx512(indices, lookup));
  }
  flat_decode_tail_avx512<kDepth>(dst, bits, table, i, count);
}

#endif  // AVX-512 VBMI

template <bool kLeftConstant, bool kRightConstant>
inline void merge_decode_avx512(std::uint8_t* dst,
                                const std::uint8_t* left,
                                std::uint8_t left_symbol,
                                const std::uint8_t* right,
                                std::uint8_t right_symbol,
                                const std::uint64_t* bits,
                                std::size_t count) noexcept {
  std::size_t out_i = 0;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  const std::size_t full = count / 64;
  for (std::size_t w = 0; w < full; ++w) {
    const __mmask64 km = bits[w];
    __m512i left_values;
    if constexpr (kLeftConstant) {
      left_values = _mm512_set1_epi8(static_cast<char>(left_symbol));
    } else {
      left_values = _mm512_loadu_si512(left + c_left);
    }
    __m512i right_values;
    if constexpr (kRightConstant) {
      right_values = _mm512_set1_epi8(static_cast<char>(right_symbol));
    } else {
      right_values = _mm512_loadu_si512(right + c_right);
    }
    // Fill clear-bit positions with left symbols, then set-bit positions with
    // right symbols. maskz zeroes the right positions first; the second expand
    // keeps the left values there.
    __m512i d;
    if constexpr (kLeftConstant) {
      d = left_values;
    } else {
      d = _mm512_maskz_expand_epi8(~km, left_values);
    }
    if constexpr (kRightConstant) {
      d = _mm512_mask_mov_epi8(d, km, right_values);
    } else {
      d = _mm512_mask_expand_epi8(d, km, right_values);
    }
    _mm512_storeu_si512(dst + out_i, d);
    const std::size_t nr = std::popcount(bits[w]);
    out_i += 64;
    c_right += nr;
    c_left += 64 - nr;
  }
  // Tail handled by the AVX2 path if available, otherwise scalar.
#if defined(__AVX2__) && defined(__SSE4_1__)
  const std::size_t done = full * 64;
  const std::uint8_t* tail_left = kLeftConstant ? left : left + c_left;
  const std::uint8_t* tail_right = kRightConstant ? right : right + c_right;
  merge_decode_avx2<kLeftConstant, kRightConstant>(
      dst + done, tail_left, left_symbol, tail_right, right_symbol, bits + full,
      count - done);
  (void)out_i;
#else
  merge_scalar_tail<kLeftConstant, kRightConstant>(
      dst, out_i, left, left_symbol, c_left, right, right_symbol, c_right, bits,
      full * 64, count);
#endif
}

template <bool kWriteLeft, bool kWriteRight>
inline void partition_encode_avx512(std::uint8_t* left_dst,
                                    std::uint8_t* right_dst,
                                    std::uint64_t* bits,
                                    const std::uint8_t* src,
                                    const std::uint8_t* dir,
                                    std::size_t weight,
                                    std::size_t left_weight,
                                    std::size_t right_weight) noexcept {
  // Build the routing mask 64 symbols at a time, store the word, then compress
  // the 64 input symbols into left/right streams with vpcompressb.
  const std::size_t full = weight / 64;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
#if defined(__AVX512VBMI__)
  const ByteLookupAvx512<8> direction_lookup = load_byte_lookup_avx512<8>(dir);
#endif
  std::size_t g = 0;
  for (; g < full; ++g) {
    const std::uint8_t* s64 = src + g * 64;
    __m512i data;
    __mmask64 km = 0;
#if defined(__AVX512VBMI__)
    data = _mm512_loadu_si512(s64);
    const __m512i directions = lookup_bytes_avx512(data, direction_lookup);
    km = _mm512_cmpneq_epi8_mask(directions, _mm512_setzero_si512());
#else
    for (int k = 0; k < 64; ++k) {
      km |= static_cast<__mmask64>(dir[s64[k]]) << k;
    }
#endif
    bits[g] = km;
    if constexpr (kWriteLeft || kWriteRight) {
#if !defined(__AVX512VBMI__)
      data = _mm512_loadu_si512(s64);
#endif
      const std::size_t nr = std::popcount(km);
      const std::size_t nl = 64 - nr;
      if constexpr (kWriteRight) {
        const __m512i right_packed = _mm512_maskz_compress_epi8(km, data);
        const std::size_t remaining = right_weight - c_right;
        if (remaining >= 64) {
          _mm512_storeu_si512(right_dst + c_right, right_packed);
        } else {
          const __mmask64 store_mask = (1ull << remaining) - 1;
          _mm512_mask_storeu_epi8(right_dst + c_right, store_mask,
                                  right_packed);
        }
      }
      if constexpr (kWriteLeft) {
        const __m512i left_packed = _mm512_maskz_compress_epi8(~km, data);
        const std::size_t remaining = left_weight - c_left;
        if (remaining >= 64) {
          _mm512_storeu_si512(left_dst + c_left, left_packed);
        } else {
          const __mmask64 store_mask = (1ull << remaining) - 1;
          _mm512_mask_storeu_epi8(left_dst + c_left, store_mask, left_packed);
        }
      }
      c_right += nr;
      c_left += nl;
    }
  }
  // Tail for the trailing < 64 symbols. `__AVX512VBMI2__` implies `__AVX2__`
  // on every real target, so delegate to the AVX2 kernel (which itself defers
  // the final < 8 symbols to the shared `partition_scalar_tail`). There is no
  // reachable `#else` branch, so no scalar tail is duplicated here.
  const std::size_t done_words = g;
  const std::size_t done_syms = g * 64;
  std::uint8_t* tail_left = kWriteLeft ? left_dst + c_left : left_dst;
  std::uint8_t* tail_right = kWriteRight ? right_dst + c_right : right_dst;
  partition_encode_avx2<kWriteLeft, kWriteRight>(
      tail_left, tail_right, bits + done_words, src + done_syms, dir,
      weight - done_syms, left_weight - c_left, right_weight - c_right);
}

#endif  // AVX-512 VBMI2

// ===========================================================================
// Scalar fallback (no SIMD).
// ===========================================================================

template <bool kLeftConstant, bool kRightConstant>
inline void merge_decode_scalar(std::uint8_t* dst,
                                const std::uint8_t* left,
                                std::uint8_t left_symbol,
                                const std::uint8_t* right,
                                std::uint8_t right_symbol,
                                const std::uint64_t* bits,
                                std::size_t count) noexcept {
  std::size_t out_i = 0;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  const std::size_t words = (count + 63) / 64;
  for (std::size_t w = 0; w < words; ++w) {
    std::uint64_t word = bits[w];
    std::size_t limit = 64;
    if (w + 1 == words) {
      limit = count - w * 64;
    }
    for (std::size_t i = 0; i < limit; ++i) {
      if (word & 1ull) {
        if constexpr (kRightConstant) {
          dst[out_i++] = right_symbol;
        } else {
          dst[out_i++] = right[c_right];
        }
        ++c_right;
      } else {
        if constexpr (kLeftConstant) {
          dst[out_i++] = left_symbol;
        } else {
          dst[out_i++] = left[c_left];
        }
        ++c_left;
      }
      word >>= 1;
    }
  }
}

template <bool kWriteLeft, bool kWriteRight>
inline void partition_encode_scalar(std::uint8_t* left_dst,
                                    std::uint8_t* right_dst,
                                    std::uint64_t* bits,
                                    const std::uint8_t* src,
                                    const std::uint8_t* dir,
                                    std::size_t weight) noexcept {
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  partition_scalar_tail<kWriteLeft, kWriteRight>(
      left_dst, right_dst, bits, src, dir, c_left, c_right, 0, weight);
}

/// @brief Portable fixed-width unpack and lookup for flat subtrees.
inline void flat_decode_scalar(std::uint8_t* dst,
                               const std::uint64_t* bits,
                               const std::uint8_t* table,
                               std::size_t count,
                               std::uint8_t depth) noexcept {
  if (depth == 8) {
    const auto* bytes = reinterpret_cast<const std::uint8_t*>(bits);
    for (std::size_t i = 0; i < count; ++i) {
      dst[i] = table[bytes[i]];
    }
    return;
  }

  const std::uint64_t mask = (std::uint64_t{1} << depth) - 1;
  std::size_t word_index = 0;
  std::size_t bit_offset = 0;
  for (std::size_t i = 0; i < count; ++i) {
    std::uint64_t packed = bits[word_index] >> bit_offset;
    if (bit_offset + depth > 64) {
      packed |= bits[word_index + 1] << (64 - bit_offset);
    }
    dst[i] = table[packed & mask];
    bit_offset += depth;
    if (bit_offset >= 64) {
      bit_offset -= 64;
      ++word_index;
    }
  }
}

template <bool kLeftConstant, bool kRightConstant>
inline void merge_decode_dispatch(std::uint8_t* dst,
                                  const std::uint8_t* left,
                                  std::uint8_t left_symbol,
                                  const std::uint8_t* right,
                                  std::uint8_t right_symbol,
                                  const std::uint64_t* bits,
                                  std::size_t count) noexcept {
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  merge_decode_avx512<kLeftConstant, kRightConstant>(
      dst, left, left_symbol, right, right_symbol, bits, count);
#elif defined(__AVX2__) && defined(__SSE4_1__)
  merge_decode_avx2<kLeftConstant, kRightConstant>(
      dst, left, left_symbol, right, right_symbol, bits, count);
#else
  merge_decode_scalar<kLeftConstant, kRightConstant>(
      dst, left, left_symbol, right, right_symbol, bits, count);
#endif
}

template <bool kWriteLeft, bool kWriteRight>
inline void partition_encode_dispatch(std::uint8_t* left_dst,
                                      std::uint8_t* right_dst,
                                      std::uint64_t* bits,
                                      const std::uint8_t* src,
                                      const std::uint8_t* dir,
                                      std::size_t weight,
                                      std::size_t left_weight,
                                      std::size_t right_weight) noexcept {
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  partition_encode_avx512<kWriteLeft, kWriteRight>(
      left_dst, right_dst, bits, src, dir, weight, left_weight, right_weight);
#elif defined(__AVX2__) && defined(__SSE4_1__)
  partition_encode_avx2<kWriteLeft, kWriteRight>(
      left_dst, right_dst, bits, src, dir, weight, left_weight, right_weight);
#else
  partition_encode_scalar<kWriteLeft, kWriteRight>(left_dst, right_dst, bits,
                                                   src, dir, weight);
  (void)left_weight;
  (void)right_weight;
#endif
}

/// @brief Dispatch the paper's merge-flat-D primitive family.
inline void flat_decode_dispatch(std::uint8_t* dst,
                                 const std::uint64_t* bits,
                                 const std::uint8_t* table,
                                 std::size_t count,
                                 std::uint8_t depth,
                                 bool known_bitreverse) noexcept {
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  if (depth == 1) {
    merge_decode_avx512<true, true>(dst, nullptr, table[0], nullptr, table[1],
                                    bits, count);
    return;
  }
#elif defined(__AVX2__) && defined(__SSE4_1__)
  if (depth == 1) {
    merge_decode_avx2<true, true>(dst, nullptr, table[0], nullptr, table[1],
                                  bits, count);
    return;
  }
#endif

#if defined(__AVX512VBMI2__) && defined(__AVX512BW__) && defined(__AVX512VBMI__)
  switch (depth) {
    case 2:
      flat_decode_avx512<2>(dst, bits, table, count);
      return;
    case 3:
      flat_decode_avx512<3>(dst, bits, table, count);
      return;
    case 4:
      flat_decode_avx512<4>(dst, bits, table, count);
      return;
    case 5:
      flat_decode_avx512<5>(dst, bits, table, count);
      return;
    case 6:
      flat_decode_avx512<6>(dst, bits, table, count);
      return;
    case 7:
      flat_decode_avx512<7>(dst, bits, table, count);
      return;
    case 8:
      flat_decode_avx512<8>(dst, bits, table, count);
      return;
    default:
      break;
  }
#endif

#if defined(__AVX2__) && defined(__SSE4_1__)
  const auto* packed = reinterpret_cast<const std::uint8_t*>(bits);
  switch (depth) {
    case 2:
      flat_decode_d2_avx2(dst, packed, table, count);
      return;
    case 3:
      flat_decode_banked_avx2<3>(dst, packed, table, count);
      return;
    case 4:
      flat_decode_d4_avx2(dst, packed, table, count);
      return;
    case 5:
      flat_decode_banked_avx2<5>(dst, packed, table, count);
      return;
    case 6:
      if (known_bitreverse) {
        flat_decode_bitreverse_unchecked_avx2<6>(dst, packed, table, count);
        return;
      }
      if (flat_decode_bitreverse_avx2<6>(dst, packed, table, count)) {
        return;
      }
      flat_decode_banked_avx2<6>(dst, packed, table, count);
      return;
    case 7:
      if (known_bitreverse) {
        flat_decode_bitreverse_unchecked_avx2<7>(dst, packed, table, count);
        return;
      }
      if (flat_decode_bitreverse_avx2<7>(dst, packed, table, count)) {
        return;
      }
      flat_decode_banked_avx2<7>(dst, packed, table, count);
      return;
    case 8:
      if (known_bitreverse) {
        flat_decode_d8_bitreverse_unchecked_avx2(dst, packed, table, count);
        return;
      }
      if (flat_decode_d8_bitreverse_avx2(dst, packed, table, count)) {
        return;
      }
      break;
    default:
      break;
  }
#endif
  flat_decode_scalar(dst, bits, table, count, depth);
}

}  // namespace pixie::pivco_simd_detail

// ---------------------------------------------------------------------------
// Public dispatch entry points.
// ---------------------------------------------------------------------------

namespace pixie {

/// @brief Bottom-up decode merge: interleave @p left / @p right child streams
///        under @p bits into @p dst, producing @p count symbols.
/// @param dst    Output buffer of (at least) @p count bytes.
/// @param left   Left child output stream (consumed left_weight symbols).
/// @param right  Right child output stream (consumed right_weight symbols).
/// @param bits   Packed 1-bit-per-symbol routing bitmap, LSB-first per word.
///               Bit value 1 selects @p right, 0 selects @p left.
/// @param count  Number of symbols (== bitmap bit count).
/// @details The caller must guarantee a few bytes of slack past the end of
///          @p dst, @p left, and @p right for vectorized over-reads/over-writes
///          (the codec's workspace is padded accordingly).
inline void pivco_merge_decode(std::uint8_t* dst,
                               const std::uint8_t* left,
                               const std::uint8_t* right,
                               const std::uint64_t* bits,
                               std::size_t count) noexcept {
  pivco_simd_detail::merge_decode_dispatch<false, false>(dst, left, 0, right, 0,
                                                         bits, count);
}

/// @brief Merge two constant-symbol leaves under a routing bitmap.
inline void pivco_merge_decode_cst_cst(std::uint8_t* dst,
                                       std::uint8_t left_symbol,
                                       std::uint8_t right_symbol,
                                       const std::uint64_t* bits,
                                       std::size_t count) noexcept {
  pivco_simd_detail::merge_decode_dispatch<true, true>(
      dst, nullptr, left_symbol, nullptr, right_symbol, bits, count);
}

/// @brief Merge a constant left leaf with a dense right child stream.
inline void pivco_merge_decode_cst_vec(std::uint8_t* dst,
                                       std::uint8_t left_symbol,
                                       const std::uint8_t* right,
                                       const std::uint64_t* bits,
                                       std::size_t count) noexcept {
  pivco_simd_detail::merge_decode_dispatch<true, false>(
      dst, nullptr, left_symbol, right, 0, bits, count);
}

/// @brief Merge a dense left child stream with a constant right leaf.
inline void pivco_merge_decode_vec_cst(std::uint8_t* dst,
                                       const std::uint8_t* left,
                                       std::uint8_t right_symbol,
                                       const std::uint64_t* bits,
                                       std::size_t count) noexcept {
  pivco_simd_detail::merge_decode_dispatch<false, true>(
      dst, left, 0, nullptr, right_symbol, bits, count);
}

/// @brief Top-down encode partition: split @p src into left/right child streams
///        and write the routing bitmap, using per-symbol direction bits.
/// @param left_dst     Output for symbols whose direction bit is 0 (exactly
///                     @p left_weight symbols).
/// @param right_dst    Output for symbols whose direction bit is 1. Must equal
///                     @p left_dst + @p left_weight: the two streams are stored
///                     contiguously, and each is immediately followed by a
///                     later-consumed sibling region, so both stores are masked
///                     to never spill past their own region end.
/// @param bits         Destination bitmap (pre-zeroed), LSB-first per word.
/// @param src          Input symbols to partition.
/// @param dir          Per-symbol direction bit table (0 or 1), indexed by
///                     symbol.
/// @param weight       Number of input symbols (== bitmap bit count).
/// @param left_weight  Number of symbols routed left (== bitmap zero-count).
/// @param right_weight Number of symbols routed right (== bitmap one-count).
inline void pivco_partition_encode(std::uint8_t* left_dst,
                                   std::uint8_t* right_dst,
                                   std::uint64_t* bits,
                                   const std::uint8_t* src,
                                   const std::uint8_t* dir,
                                   std::size_t weight,
                                   std::size_t left_weight,
                                   std::size_t right_weight) noexcept {
  pivco_simd_detail::partition_encode_dispatch<true, true>(
      left_dst, right_dst, bits, src, dir, weight, left_weight, right_weight);
}

/// @brief Build a routing bitmap and materialize only the left child stream.
inline void pivco_partition_encode_left(std::uint8_t* left_dst,
                                        std::uint64_t* bits,
                                        const std::uint8_t* src,
                                        const std::uint8_t* dir,
                                        std::size_t weight,
                                        std::size_t left_weight) noexcept {
  pivco_simd_detail::partition_encode_dispatch<true, false>(
      left_dst, nullptr, bits, src, dir, weight, left_weight, 0);
}

/// @brief Build a routing bitmap and materialize only the right child stream.
inline void pivco_partition_encode_right(std::uint8_t* right_dst,
                                         std::uint64_t* bits,
                                         const std::uint8_t* src,
                                         const std::uint8_t* dir,
                                         std::size_t weight,
                                         std::size_t right_weight) noexcept {
  pivco_simd_detail::partition_encode_dispatch<false, true>(
      nullptr, right_dst, bits, src, dir, weight, 0, right_weight);
}

/// @brief Build only the routing bitmap when both children are leaves.
inline void pivco_partition_encode_none(std::uint64_t* bits,
                                        const std::uint8_t* src,
                                        const std::uint8_t* dir,
                                        std::size_t weight) noexcept {
  pivco_simd_detail::partition_encode_dispatch<false, false>(
      nullptr, nullptr, bits, src, dir, weight, 0, 0);
}

/// @brief Build a bitmap for a node whose two children are constant leaves.
/// @details A set bit means @p src equals @p right_symbol. Unlike the generic
///          direction-table path, AVX2 can compare 32 symbols directly and
///          obtain the complete routing word fragment with `vpmovmskb`.
inline void pivco_partition_encode_cst_cst(std::uint64_t* bits,
                                           const std::uint8_t* src,
                                           std::uint8_t right_symbol,
                                           std::size_t weight) noexcept {
  std::size_t i = 0;
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  const __m512i right = _mm512_set1_epi8(static_cast<char>(right_symbol));
  for (; i + 64 <= weight; i += 64) {
    const __m512i symbols = _mm512_loadu_si512(src + i);
    bits[i / 64] = _mm512_cmpeq_epi8_mask(symbols, right);
  }
#elif defined(__AVX2__)
  const __m256i right = _mm256_set1_epi8(static_cast<char>(right_symbol));
  auto* bit_bytes = reinterpret_cast<std::uint8_t*>(bits);
  for (; i + 32 <= weight; i += 32) {
    const __m256i symbols =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + i));
    const std::uint32_t routing = static_cast<std::uint32_t>(
        _mm256_movemask_epi8(_mm256_cmpeq_epi8(symbols, right)));
    std::memcpy(bit_bytes + i / 8, &routing, sizeof(routing));
  }
#endif
  for (; i < weight; ++i) {
    if (src[i] == right_symbol) {
      bits[i / 64] |= std::uint64_t{1} << (i % 64);
    }
  }
}

#if defined(__AVX2__) && defined(__SSE4_1__)
/// @brief Partition a constant left leaf and materialize the right stream.
inline void pivco_partition_encode_cst_vec(std::uint8_t* right_dst,
                                           std::uint64_t* bits,
                                           const std::uint8_t* src,
                                           std::uint8_t left_symbol,
                                           std::size_t weight,
                                           std::size_t right_weight) noexcept {
  pivco_simd_detail::partition_encode_one_constant_avx2<true>(
      right_dst, bits, src, left_symbol, weight, right_weight);
}

/// @brief Partition a dense left stream and a constant right leaf.
inline void pivco_partition_encode_vec_cst(std::uint8_t* left_dst,
                                           std::uint64_t* bits,
                                           const std::uint8_t* src,
                                           std::uint8_t right_symbol,
                                           std::size_t weight,
                                           std::size_t left_weight) noexcept {
  pivco_simd_detail::partition_encode_one_constant_avx2<false>(
      left_dst, bits, src, right_symbol, weight, left_weight);
}
#endif

#if defined(__AVX2__) && defined(__SSE4_1__) && defined(__BMI2__)
/// @brief Encode a balanced depth-2..7 table as AVX2 bit reversal.
/// @param known_bitreverse Skip the table scan when construction already
///        classified the mapping.
inline bool pivco_flat_encode_bitreverse(
    std::uint64_t* bits,
    const std::uint8_t* src,
    const std::uint8_t* wire_code,
    std::size_t count,
    std::uint8_t depth,
    bool known_bitreverse = false) noexcept {
  switch (depth) {
    case 2:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<2>(
          bits, src, wire_code, count, known_bitreverse);
    case 3:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<3>(
          bits, src, wire_code, count, known_bitreverse);
    case 4:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<4>(
          bits, src, wire_code, count, known_bitreverse);
    case 5:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<5>(
          bits, src, wire_code, count, known_bitreverse);
    case 6:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<6>(
          bits, src, wire_code, count, known_bitreverse);
    case 7:
      return pivco_simd_detail::flat_encode_bitreverse_maybe_avx2<7>(
          bits, src, wire_code, count, known_bitreverse);
    default:
      return false;
  }
}

/// @brief Pack a depth-3 flat subtree through the AVX2/BMI2 kernel.
/// @pre @p wire_code has 259 readable bytes.
inline void pivco_flat_encode_d3(std::uint64_t* bits,
                                 const std::uint8_t* src,
                                 const std::uint8_t* wire_code,
                                 std::size_t count) noexcept {
  pivco_simd_detail::flat_encode_d3_avx2(bits, src, wire_code, count);
}
#endif

#if defined(__AVX2__) && defined(__SSE4_1__)
/// @brief Encode the common full-alphabet bit-reversal table with AVX2.
/// @returns True when @p wire_code is the expected bit-reversal permutation.
/// @param known_bitreverse Skip the table scan when construction already
///        classified the mapping.
inline bool pivco_flat_encode_d8_bitreverse(
    std::uint64_t* bits,
    const std::uint8_t* src,
    const std::uint8_t* wire_code,
    std::size_t count,
    bool known_bitreverse = false) noexcept {
  if (known_bitreverse) {
    pivco_simd_detail::flat_decode_d8_bitreverse_unchecked_avx2(
        reinterpret_cast<std::uint8_t*>(bits), src, wire_code, count);
    return true;
  }
  return pivco_simd_detail::flat_decode_d8_bitreverse_avx2(
      reinterpret_cast<std::uint8_t*>(bits), src, wire_code, count);
}
#endif

/// @brief Unpack and translate a flat subtree's fixed-width wire codes.
/// @param known_bitreverse Skip table classification when the caller already
///        established the balanced mapping during construction.
inline void pivco_flat_decode(std::uint8_t* dst,
                              const std::uint64_t* bits,
                              const std::uint8_t* table,
                              std::size_t count,
                              std::uint8_t depth,
                              bool known_bitreverse = false) noexcept {
  pivco_simd_detail::flat_decode_dispatch(dst, bits, table, count, depth,
                                          known_bitreverse);
}

}  // namespace pixie
