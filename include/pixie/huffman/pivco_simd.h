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
 *
 * Bitmap convention (matches `pivco_huffman.h`): bits are packed LSB-first
 * into 64-bit words; bit 0 of word 0 is the first symbol.
 *
 * Backends, selected by compiler feature macros (`-march`):
 * - AVX-512 VBMI2 (`__AVX512VBMI2__`): byte-granular `vpexpandb`/`vpcompressb`
 *   operating on 64 bytes per iteration. This is the instruction the paper
 *   (Sect. 6.6) identifies as mapping directly to the partition kernel.
 * - AVX2 + SSE4.1: 8-byte-at-a-time `_mm_shuffle_epi8` with 256-entry lookup
 *   tables, mirroring the paper's NEON `vqtbl1q_u8` technique.
 * - Scalar: the original bit-by-bit / symbol-by-symbol loops.
 *
 * Every SIMD backend produces byte-identical output to the scalar fallback
 * (the operation is a deterministic permutation). The existing round-trip tests
 * under release and AddressSanitizer serve as the differential oracle for every
 * backend the test matrix builds and runs; the AVX-512 VBMI2 backend is
 * verified only when built on VBMI2-capable hardware, since there is no runtime
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

/// @brief Read 8 bits starting at bit position @p pos from a packed word array.
/// @pre @p pos + 8 is within the bitmap's bit range; callers stop the SIMD
///      loop once fewer than 8 symbols remain, so the span never exceeds the
///      allocated words.
inline std::uint8_t extract_byte(const std::uint64_t* bits,
                                 std::size_t pos) noexcept {
  const std::size_t w = pos / 64;
  const std::size_t b = pos % 64;
  std::uint64_t v = bits[w] >> b;
  if (b + 8 > 64) {
    v |= bits[w + 1] << (64 - b);
  }
  return static_cast<std::uint8_t>(v);
}

/// @brief Scalar merge tail: process the remaining `< 8` symbols bit-by-bit.
inline void merge_scalar_tail(std::uint8_t* dst,
                              std::size_t& out_i,
                              const std::uint8_t* left,
                              std::size_t& c_left,
                              const std::uint8_t* right,
                              std::size_t& c_right,
                              const std::uint64_t* bits,
                              std::size_t i,
                              std::size_t count) noexcept {
  for (; i < count; ++i) {
    const std::uint8_t bit =
        static_cast<std::uint8_t>((bits[i / 64] >> (i % 64)) & 1u);
    dst[out_i++] = bit ? right[c_right++] : left[c_left++];
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
      right_dst[c_right++] = s;
    } else {
      left_dst[c_left++] = s;
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

inline void merge_decode_avx2(std::uint8_t* dst,
                              const std::uint8_t* left,
                              const std::uint8_t* right,
                              const std::uint64_t* bits,
                              std::size_t count) noexcept {
  const auto& mlut = merge_lut();
  std::size_t out_i = 0;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  std::size_t i = 0;
  for (; i + 8 <= count; i += 8) {
    const std::uint8_t m8 = extract_byte(bits, i);
    std::uint64_t lv = 0;
    std::uint64_t rv = 0;
    std::memcpy(&lv, left + c_left, 8);
    std::memcpy(&rv, right + c_right, 8);
    const __m128i combined =
        _mm_set_epi64x(static_cast<long long>(rv), static_cast<long long>(lv));
    const __m128i out = _mm_shuffle_epi8(combined, load_shuf(mlut, m8));
    _mm_storel_epi64(reinterpret_cast<__m128i*>(dst + out_i), out);
    out_i += 8;
    const int nr = std::popcount(static_cast<unsigned>(m8));
    c_right += static_cast<std::size_t>(nr);
    c_left += static_cast<std::size_t>(8 - nr);
  }
  merge_scalar_tail(dst, out_i, left, c_left, right, c_right, bits, i, count);
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

inline void partition_encode_avx2(std::uint8_t* left_dst,
                                  std::uint8_t* right_dst,
                                  std::uint64_t* bits,
                                  const std::uint8_t* src,
                                  const std::uint8_t* dir,
                                  std::size_t weight,
                                  std::size_t left_weight,
                                  std::size_t right_weight) noexcept {
  const auto& cr = compress_right_lut();
  const auto& cl = compress_left_lut();
  std::uint8_t* bits_bytes = reinterpret_cast<std::uint8_t*>(bits);

  const std::size_t full = weight / 8;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  std::size_t g = 0;
  for (; g < full; ++g) {
    const std::uint8_t* s8 = src + g * 8;
    const __m128i data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s8));
    const std::uint8_t m8 = build_dir_mask8(s8, dir);
    // Byte-aligned write: full groups always start on a byte boundary in the
    // LSB-first bitmap, so byte g holds exactly symbols [8g, 8g+8).
    bits_bytes[g] = m8;
    const __m128i rv = _mm_shuffle_epi8(data, load_shuf(cr, m8));
    const __m128i lv = _mm_shuffle_epi8(data, load_shuf(cl, m8));
    const int nr = std::popcount(static_cast<unsigned>(m8));
    const int nl = 8 - nr;
    // Both child regions are contiguous with a later-consumed sibling region,
    // so neither store may overflow past its own region end. Masked stores
    // keep the trailing invalid lanes inside the region when they would
    // otherwise spill into the neighbor.
    store_masked_epi64(right_dst + c_right, rv,
                       static_cast<int>(right_weight - c_right));
    store_masked_epi64(left_dst + c_left, lv,
                       static_cast<int>(left_weight - c_left));
    c_right += static_cast<std::size_t>(nr);
    c_left += static_cast<std::size_t>(nl);
  }

  // Scalar tail for the trailing < 8 symbols (byte-unaligned bitmap bits).
  partition_scalar_tail(left_dst, right_dst, bits, src, dir, c_left, c_right,
                        full * 8, weight);
}

#endif  // AVX2 + SSE4.1

// ===========================================================================
// AVX-512 VBMI2 backend (byte-granular expand/compress, 64 bytes/iteration).
// ===========================================================================
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)

inline void merge_decode_avx512(std::uint8_t* dst,
                                const std::uint8_t* left,
                                const std::uint8_t* right,
                                const std::uint64_t* bits,
                                std::size_t count) noexcept {
  std::size_t out_i = 0;
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  const std::size_t full = count / 64;
  for (std::size_t w = 0; w < full; ++w) {
    const __mmask64 km = bits[w];
    const __m512i lv = _mm512_loadu_si512(left + c_left);
    const __m512i rv = _mm512_loadu_si512(right + c_right);
    // Fill clear-bit positions with left symbols, then set-bit positions with
    // right symbols. maskz zeroes the right positions first; the second expand
    // keeps the left values there.
    __m512i d = _mm512_maskz_expand_epi8(~km, lv);
    d = _mm512_mask_expand_epi8(d, km, rv);
    _mm512_storeu_si512(dst + out_i, d);
    const std::size_t nr = std::popcount(bits[w]);
    out_i += 64;
    c_right += nr;
    c_left += 64 - nr;
  }
  // Tail handled by the AVX2 path if available, otherwise scalar.
#if defined(__AVX2__) && defined(__SSE4_1__)
  const std::size_t done = full * 64;
  merge_decode_avx2(dst + done, left + c_left, right + c_right, bits + full,
                    count - done);
  (void)out_i;
#else
  merge_scalar_tail(dst, out_i, left, c_left, right, c_right, bits, full * 64,
                    count);
#endif
}

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
  std::size_t g = 0;
  for (; g < full; ++g) {
    const std::uint8_t* s64 = src + g * 64;
    const __m512i data = _mm512_loadu_si512(s64);
    // Direction mask: bit k = dir[symbol k], LSB-first. x86 has no
    // byte-granular gather, so the 64-bit mask is built with a scalar OR loop
    // over the dir table. (A qword gather + vpmovb2m would be faster in
    // principle but reads all 64 byte sign-bits, corrupting the mask with
    // neighboring dir entries.)
    __mmask64 km = 0;
    for (int k = 0; k < 64; ++k) {
      km |= static_cast<__mmask64>(dir[s64[k]]) << k;
    }
    bits[g] = km;
    // Compress: pack right-bound (set-bit) elements to the front of one stream
    // and left-bound (clear-bit) elements to the front of another.
    const __m512i right_packed = _mm512_maskz_compress_epi8(km, data);
    const __m512i left_packed = _mm512_maskz_compress_epi8(~km, data);
    const std::size_t nr = std::popcount(km);
    const std::size_t nl = 64 - nr;
    // Both child regions are contiguous with a later-consumed sibling region,
    // so neither store may overflow its own region. Use masked stores that
    // write only the valid lanes when the full vector would spill into the
    // neighbor.
    {
      const std::size_t rem = right_weight - c_right;
      if (rem >= 64) {
        _mm512_storeu_si512(right_dst + c_right, right_packed);
      } else {
        const __mmask64 rmask = (1ull << rem) - 1;
        _mm512_mask_storeu_epi8(right_dst + c_right, rmask, right_packed);
      }
    }
    {
      const std::size_t rem = left_weight - c_left;
      if (rem >= 64) {
        _mm512_storeu_si512(left_dst + c_left, left_packed);
      } else {
        const __mmask64 lmask = (1ull << rem) - 1;
        _mm512_mask_storeu_epi8(left_dst + c_left, lmask, left_packed);
      }
    }
    c_right += nr;
    c_left += nl;
  }
  // Tail for the trailing < 64 symbols. `__AVX512VBMI2__` implies `__AVX2__`
  // on every real target, so delegate to the AVX2 kernel (which itself defers
  // the final < 8 symbols to the shared `partition_scalar_tail`). There is no
  // reachable `#else` branch, so no scalar tail is duplicated here.
  const std::size_t done_words = g;
  const std::size_t done_syms = g * 64;
  partition_encode_avx2(left_dst + c_left, right_dst + c_right,
                        bits + done_words, src + done_syms, dir,
                        weight - done_syms, left_weight - c_left,
                        right_weight - c_right);
}

#endif  // AVX-512 VBMI2

// ===========================================================================
// Scalar fallback (no SIMD).
// ===========================================================================

inline void merge_decode_scalar(std::uint8_t* dst,
                                const std::uint8_t* left,
                                const std::uint8_t* right,
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
      dst[out_i++] = (word & 1ull) ? right[c_right++] : left[c_left++];
      word >>= 1;
    }
  }
}

inline void partition_encode_scalar(std::uint8_t* left_dst,
                                    std::uint8_t* right_dst,
                                    std::uint64_t* bits,
                                    const std::uint8_t* src,
                                    const std::uint8_t* dir,
                                    std::size_t weight,
                                    std::size_t /*left_weight*/,
                                    std::size_t /*right_weight*/) noexcept {
  std::size_t c_left = 0;
  std::size_t c_right = 0;
  partition_scalar_tail(left_dst, right_dst, bits, src, dir, c_left, c_right, 0,
                        weight);
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
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  pivco_simd_detail::merge_decode_avx512(dst, left, right, bits, count);
#elif defined(__AVX2__) && defined(__SSE4_1__)
  pivco_simd_detail::merge_decode_avx2(dst, left, right, bits, count);
#else
  pivco_simd_detail::merge_decode_scalar(dst, left, right, bits, count);
#endif
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
#if defined(__AVX512VBMI2__) && defined(__AVX512BW__)
  pivco_simd_detail::partition_encode_avx512(
      left_dst, right_dst, bits, src, dir, weight, left_weight, right_weight);
#elif defined(__AVX2__) && defined(__SSE4_1__)
  pivco_simd_detail::partition_encode_avx2(left_dst, right_dst, bits, src, dir,
                                           weight, left_weight, right_weight);
#else
  pivco_simd_detail::partition_encode_scalar(
      left_dst, right_dst, bits, src, dir, weight, left_weight, right_weight);
#endif
}

}  // namespace pixie
