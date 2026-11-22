#pragma once

/**
 * @file forest256_simd.h
 * @brief Internal 512-bit level-order forest encodings and SIMD unpackers.
 *
 * @details Provides two 64-byte (512-bit) representations for ordered forests
 * of `1 <= count <= 256` nodes numbered in breadth-first (level) order
 * `0 .. count - 1`, together with AVX-512, AVX2/SSE2, and portable C++20
 * kernels that bulk-unpack either encoding into the same sibling-contiguous
 * byte-width CSR layout (`LevelCsr256`):
 *
 * 1. `LoudsForest256` (truncated 512-bit LOUDS): concatenates `1^r 0` for the
 *    `r` roots followed by `1^deg(u) 0` for each node `u` in level order,
 *    omitting bit 512 for `count == 256` (or tail-padding bits
 *    `2 * count + 1 .. 511` with `1010...1` for `count < 256` so every block
 *    has 256 ones and 256 zeros).
 * 2. `TwinForest256` (dual-bitplane level-order `(H, L)`): stores two 256-bit
 *    masks `words[0..3] = H` (`H[u] = 1` iff `deg(u) > 0`) and
 *    `words[4..7] = L` (`L[v] = 1` iff `v` is the last node in its sibling
 *    group, including `r - 1` for the root group).
 *
 * Dense streams omit padding and use exactly two bits per node. See
 * `tree/experimental/agentic_notes/forest256-handoff.md` for algorithm
 * derivations and Zen 4 benchmark measurements.
 */

#include <pixie/bits.h>

#if defined(__SSE2__) && !defined(PIXIE_AVX2_SUPPORT) &&      \
    !defined(PIXIE_SSE41_SUPPORT) && !defined(__AVX512F__) && \
    !defined(PIXIE_BMI2_SUPPORT)
#include <emmintrin.h>
#endif

#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
#define PIXIE_FOREST256_AVX512_VBMI2
#endif

#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI__) && \
    defined(__AVX512VBMI2__) && defined(__AVX512BITALG__)
#define PIXIE_FOREST256_AVX512_FULL
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

namespace pixie::tree::detail::forest256 {

/// Root sentinel in `LevelCsr256::parent`: valid BFS parent IDs lie in
/// `0..254`.
inline constexpr uint8_t kNoParent = 0xFF;

/**
 * @brief Unpacked sibling-contiguous CSR + parent table for up to 256 BFS
 * nodes.
 *
 * @details For a forest of `1 <= count <= 256` nodes with `r` roots `0 .. r-1`:
 * - `parent[v]` (`0 <= v < count`) is `kNoParent` (`0xFF`) for roots `0 .. r-1`
 *   and the BFS node ID (`0 .. v-1`) of `v`'s parent otherwise.
 * - `offset_m1[u]` stores `first_child_offset(u) - 1`; `offset_m1[0] = r - 1`.
 * - Short units have `offset_m1[count] = count - 1` and a `parent[count]`
 *   distinct from the preceding parent. For a full unit, node 255 is always
 *   a leaf and last in its sibling group; neither array has a 257th entry.
 */
struct alignas(64) LevelCsr256 {
  uint8_t parent[256];
  uint8_t offset_m1[256];

  /** @brief Number of roots `r` (`1 <= r <= count`) in the unpacked forest. */
  constexpr size_t root_count() const noexcept {
    return static_cast<size_t>(offset_m1[0]) + 1;
  }

  /** @brief True iff node @p v (`0 <= v < count`) is a root of the forest. */
  constexpr bool is_root(uint8_t v) const noexcept {
    return parent[v] == kNoParent;
  }

  /** @brief Parent of @p v (`0 <= v < count`), or `kNoParent` if @p v is a
   * root. */
  constexpr uint8_t parent_of(uint8_t v) const noexcept { return parent[v]; }

  /** @brief Number of children of node @p u (`0 <= u < count`). */
  constexpr size_t degree(uint8_t u) const noexcept {
    if (u == 255) {
      return 0;
    }
    return static_cast<size_t>(offset_m1[static_cast<size_t>(u) + 1]) -
           static_cast<size_t>(offset_m1[u]);
  }

  /** @brief True iff node @p u (`0 <= u < count`) has zero children. */
  constexpr bool is_leaf(uint8_t u) const noexcept {
    return u == 255 || offset_m1[static_cast<size_t>(u) + 1] == offset_m1[u];
  }

  /**
   * @brief First child of non-leaf node @p u (`0 <= u < count`, `!is_leaf(u)`).
   */
  constexpr uint8_t first_child(uint8_t u) const noexcept {
    return static_cast<uint8_t>(offset_m1[u] + 1);
  }

  /**
   * @brief Zero-based @p k-th child (`0 <= k < degree(u)`) of node @p u.
   */
  constexpr uint8_t child(uint8_t u, size_t k) const noexcept {
    return static_cast<uint8_t>(offset_m1[u] + 1 + k);
  }

  /** @brief True iff @p v (`0 <= v < count`) is the last node in its sibling
   * group. */
  constexpr bool is_last_child(uint8_t v) const noexcept {
    return v == 255 || parent[v] != parent[static_cast<size_t>(v) + 1];
  }

  /**
   * @brief Next sibling of @p v (`0 <= v < count`, `!is_last_child(v)`).
   */
  constexpr uint8_t next_sibling(uint8_t v) const noexcept {
    return static_cast<uint8_t>(v + 1);
  }
};
static_assert(sizeof(LevelCsr256) == 512);

/**
 * @brief 64-byte (512-bit) truncated LOUDS representation of a `<= 256`-node
 * forest.
 */
struct alignas(64) LoudsForest256 {
  uint64_t words[8]{};
};

/**
 * @brief 64-byte (512-bit) Twin-256 `(H, L)` representation of a `<= 256`-node
 * forest.
 * @details `words[0..3]` store `H` (`deg(u) > 0`); `words[4..7]` store `L`
 * (`is_last_sibling(v)`).
 */
struct alignas(64) TwinForest256 {
  uint64_t words[8]{};
};

namespace detail {

constexpr void set_bit(uint64_t* words, size_t pos) noexcept {
  words[pos >> 6] |= (uint64_t{1} << (pos & 63));
}

inline void validate_bfs_degrees(size_t root_count,
                                 std::span<const uint16_t> bfs_degrees) {
  const size_t count = bfs_degrees.size();
  if (count == 0 || count > 256) {
    throw std::invalid_argument("Forest256 requires 1 <= count <= 256");
  }
  if (root_count == 0 || root_count > count) {
    throw std::invalid_argument("Forest256 requires 1 <= root_count <= count");
  }
  size_t frontier_end = root_count;
  for (size_t u = 0; u < count; ++u) {
    if (u >= frontier_end) {
      throw std::invalid_argument("Disconnected BFS degree sequence");
    }
    const size_t deg = bfs_degrees[u];
    if (deg > count - frontier_end) {
      throw std::invalid_argument("BFS degree sum exceeds node count");
    }
    frontier_end += deg;
  }
  if (frontier_end != count) {
    throw std::invalid_argument("BFS degree sum does not match count - roots");
  }
}

struct LoudsZeroBeforeOneTable {
  std::array<uint64_t, 256> packed{};
  constexpr LoudsZeroBeforeOneTable() noexcept {
    for (unsigned m = 0; m < 256; ++m) {
      uint64_t val = 0;
      unsigned ones = 0;
      for (unsigned b = 0; b < 8; ++b) {
        if ((m >> b) & 1u) {
          val |= (static_cast<uint64_t>(b - ones) << (8 * ones));
          ++ones;
        }
      }
      packed[m] = val;
    }
  }
};

inline constexpr LoudsZeroBeforeOneTable kLoudsDiffLut{};

struct ByteCompressTable {
  std::array<uint64_t, 256> packed{};
  constexpr ByteCompressTable() noexcept {
    for (unsigned m = 0; m < 256; ++m) {
      uint64_t val = 0;
      unsigned k = 0;
      for (unsigned b = 0; b < 8; ++b) {
        if ((m >> b) & 1u) {
          val |= (static_cast<uint64_t>(b) << (8 * k));
          ++k;
        }
      }
      packed[m] = val;
    }
  }
};

inline constexpr ByteCompressTable kByteCompressLut{};

[[gnu::always_inline]] inline uint64_t mask_low_bytes_u64(
    uint64_t val,
    unsigned cnt) noexcept {
#ifdef PIXIE_BMI2_SUPPORT
  return _bzhi_u64(val, cnt * 8);
#else
  if (cnt == 0) {
    return 0;
  }
  return val & (~uint64_t{0} >> (64 - cnt * 8));
#endif
}

[[gnu::always_inline]] inline uint64_t mask_bits_from_u64(unsigned k) noexcept {
#ifdef PIXIE_BMI2_SUPPORT
  return ~_bzhi_u64(~uint64_t{0}, k);
#else
  return (k >= 64) ? 0ULL : (~uint64_t{0} << k);
#endif
}

#ifdef PIXIE_FOREST256_AVX512_FULL
struct ZmmTable256 {
  __m512i t0;
  __m512i t1;
  __m512i t2;
  __m512i t3;
};

[[gnu::always_inline]] inline __m512i rank1_excl_64_avx512(
    uint64_t w,
    uint8_t wbase) noexcept {
  const __m512i kBytePrefixMask = _mm512_setr_epi64(
      0x0000000000000000ULL, 0x00000000000000FFULL, 0x000000000000FFFFULL,
      0x0000000000FFFFFFULL, 0x00000000FFFFFFFFULL, 0x000000FFFFFFFFFFULL,
      0x0000FFFFFFFFFFFFULL, 0x00FFFFFFFFFFFFFFULL);
  const __m512i kBroadcastLowBytePerQword = _mm512_setr_epi64(
      0x0000000000000000ULL, 0x0808080808080808ULL, 0x0000000000000000ULL,
      0x0808080808080808ULL, 0x0000000000000000ULL, 0x0808080808080808ULL,
      0x0000000000000000ULL, 0x0808080808080808ULL);
  const __m512i kSelectEachByteOfQword = _mm512_setr_epi64(
      0x0000000000000000ULL, 0x0909090909090909ULL, 0x0202020202020202ULL,
      0x0B0B0B0B0B0B0B0BULL, 0x0404040404040404ULL, 0x0D0D0D0D0D0D0D0DULL,
      0x0606060606060606ULL, 0x0F0F0F0F0F0F0F0FULL);
  const __m512i kBitPrefixMask = _mm512_set1_epi64(0x7F3F1F0F07030100ULL);

  const __m512i w64 = _mm512_set1_epi64(static_cast<long long>(w));
  const __m512i byte_pop64 =
      _mm512_popcnt_epi64(_mm512_and_si512(w64, kBytePrefixMask));
  const __m512i byte_base =
      _mm512_shuffle_epi8(byte_pop64, kBroadcastLowBytePerQword);
  const __m512i curr_byte = _mm512_shuffle_epi8(w64, kSelectEachByteOfQword);
  const __m512i bit_pop =
      _mm512_popcnt_epi8(_mm512_and_si512(curr_byte, kBitPrefixMask));
  return _mm512_add_epi8(_mm512_add_epi8(byte_base, bit_pop),
                         _mm512_set1_epi8(static_cast<char>(wbase)));
}

[[gnu::always_inline]] inline ZmmTable256 compress_256_branchless_avx512(
    const uint64_t m[4]) noexcept {
  const __m512i iota0 = _mm512_setr_epi64(
      0x0706050403020100ULL, 0x0F0E0D0C0B0A0908ULL, 0x1716151413121110ULL,
      0x1F1E1D1C1B1A1918ULL, 0x2726252423222120ULL, 0x2F2E2D2C2B2A2928ULL,
      0x3736353433323130ULL, 0x3F3E3D3C3B3A3938ULL);
  const __m512i iota64 = _mm512_add_epi8(iota0, _mm512_set1_epi8(64));
  const __m512i iota128 = _mm512_add_epi8(iota0, _mm512_set1_epi8(-128));
  const __m512i iota192 = _mm512_add_epi8(iota0, _mm512_set1_epi8(-64));

  const __m512i z0 = _mm512_maskz_compress_epi8(m[0], iota0);
  const __m512i z1 = _mm512_maskz_compress_epi8(m[1], iota64);
  const __m512i z2 = _mm512_maskz_compress_epi8(m[2], iota128);
  const __m512i z3 = _mm512_maskz_compress_epi8(m[3], iota192);

  const unsigned c0 = static_cast<unsigned>(std::popcount(m[0]));
  const unsigned c1 = static_cast<unsigned>(std::popcount(m[1]));
  const unsigned c2 = static_cast<unsigned>(std::popcount(m[2]));

  const __m512i idx_p =
      _mm512_sub_epi8(iota0, _mm512_set1_epi8(static_cast<char>(c0)));
  const __mmask64 mask_p = mask_bits_from_u64(c0);
  const __m512i p0 = _mm512_mask_permutexvar_epi8(z0, mask_p, idx_p, z1);
  const __m512i p1 = _mm512_permutexvar_epi8(idx_p, z1);

  const __m512i idx_q =
      _mm512_sub_epi8(iota0, _mm512_set1_epi8(static_cast<char>(c2)));
  const __mmask64 mask_q = mask_bits_from_u64(c2);
  const __m512i q0 = _mm512_mask_permutexvar_epi8(z2, mask_q, idx_q, z3);
  const __m512i q1 = _mm512_permutexvar_epi8(idx_q, z3);

  const unsigned c01 = c0 + c1;
  const __m512i vc01 = _mm512_set1_epi8(static_cast<char>(c01));
  const __mmask64 mq0 = mask_bits_from_u64(c01);
  const unsigned rem1 = (c01 > 64) ? (c01 - 64) : 0u;
  const __mmask64 mq1 = mask_bits_from_u64(rem1);

  ZmmTable256 res;
  res.t0 =
      _mm512_mask_permutexvar_epi8(p0, mq0, _mm512_sub_epi8(iota0, vc01), q0);
  res.t1 = _mm512_mask_mov_epi8(
      p1, mq1, _mm512_permutex2var_epi8(q0, _mm512_sub_epi8(iota64, vc01), q1));
  res.t2 = _mm512_permutex2var_epi8(q0, _mm512_sub_epi8(iota0, vc01), q1);
  res.t3 = _mm512_permutex2var_epi8(q0, _mm512_sub_epi8(iota64, vc01), q1);
  return res;
}

[[gnu::always_inline]] inline __m512i lookup_zmm_branchless_avx512(
    const ZmmTable256& tbl,
    __m512i idx) noexcept {
  const __m512i lo = _mm512_permutex2var_epi8(tbl.t0, idx, tbl.t1);
  const __m512i hi = _mm512_permutex2var_epi8(tbl.t2, idx, tbl.t3);
  return _mm512_mask_mov_epi8(lo, _mm512_movepi8_mask(idx), hi);
}
#endif  // PIXIE_FOREST256_AVX512_FULL

}  // namespace detail

// ============================================================================
// Reference 512-bit Encoders (1 <= count <= 256, 1 <= root_count <= count)
// ============================================================================

/**
 * @brief Encode a level-order forest into 512-bit truncated LOUDS.
 *
 * @param root_count Number of forest roots (`1 <= root_count <= count`).
 * @param bfs_degrees Degree of each BFS node `0 .. count - 1` (`1 <= count <=
 * 256`).
 * @return 64-byte `LoudsForest256` with 256 ones and 256 zeros.
 */
inline LoudsForest256 encode_louds_forest256(
    size_t root_count,
    std::span<const uint16_t> bfs_degrees) {
  detail::validate_bfs_degrees(root_count, bfs_degrees);
  const size_t count = bfs_degrees.size();
  LoudsForest256 out{};
  size_t bit = 0;
  for (size_t i = 0; i < root_count; ++i) {
    detail::set_bit(out.words, bit++);
  }
  ++bit;  // super-root 0 delimiter
  for (size_t u = 0; u < count; ++u) {
    const size_t deg = bfs_degrees[u];
    for (size_t j = 0; j < deg; ++j) {
      detail::set_bit(out.words, bit++);
    }
    if (bit < 512) {
      ++bit;  // node 0 delimiter (omitted at bit 512 when count == 256)
    }
  }
  // When count < 256, bit == 2 * count + 1: pad remaining odd bit range
  // 2 * count + 1 .. 511 with alternating 1, 0, ..., 1 so total ones == 256,
  // total zeros == 256, and offset_m1[count] == count - 1. Unpackers set
  // a distinct parent[count] for these short units.
  for (size_t b = bit; b < 512; b += 2) {
    detail::set_bit(out.words, b);
  }
  return out;
}

/**
 * @brief Encode a level-order forest into 512-bit Twin-256 `(H, L)`.
 *
 * @param root_count Number of forest roots (`1 <= root_count <= count`).
 * @param bfs_degrees Degree of each BFS node `0 .. count - 1` (`1 <= count <=
 * 256`).
 * @return 64-byte `TwinForest256` (`words[0..3] = H`, `words[4..7] = L`).
 */
inline TwinForest256 encode_twin_forest256(
    size_t root_count,
    std::span<const uint16_t> bfs_degrees) {
  detail::validate_bfs_degrees(root_count, bfs_degrees);
  const size_t count = bfs_degrees.size();
  TwinForest256 out{};
  uint64_t* H = out.words;
  uint64_t* L = out.words + 4;

  size_t end = root_count;
  detail::set_bit(L, end - 1);
  size_t internal_count = 0;
  for (size_t u = 0; u < count; ++u) {
    const size_t deg = bfs_degrees[u];
    if (deg != 0) {
      detail::set_bit(H, u);
      ++internal_count;
      end += deg;
      detail::set_bit(L, end - 1);
    }
  }
  if (count < 255) {
    detail::set_bit(H, count);
    detail::set_bit(L, 255);
  } else if (count == 255) {
    if (internal_count != 0) {
      detail::set_bit(H, 255);
    }
    detail::set_bit(L, 255);
  }
  return out;
}

// ============================================================================
// Truncated LOUDS Unpacking Kernels
// ============================================================================

/**
 * @brief Portable scalar reference unpacker for `LoudsForest256`.
 */
inline void unpack_louds_scalar(const uint64_t words[8],
                                LevelCsr256& out,
                                size_t count) noexcept {
  size_t ones = 0;
  size_t zeros = 0;
  for (size_t w = 0; w < 8; ++w) {
    const uint64_t bits = words[w];
    for (unsigned b = 0; b < 64; ++b) {
      if ((bits >> b) & 1ULL) {
        out.parent[ones] =
            (zeros == 0) ? kNoParent : static_cast<uint8_t>(zeros - 1);
        ++ones;
      } else {
        out.offset_m1[zeros] = static_cast<uint8_t>(ones - 1);
        ++zeros;
      }
    }
  }
  if (count < 256) {
    out.parent[count] = out.parent[count - 1] ^ 1u;
  }
}

#ifdef __SSE2__
/**
 * @brief SSE2/AVX2 2 KiB LUT unpacker for `LoudsForest256`.
 *
 * @details Uses `_mm_add_epi8` to avoid byte-carry wraparound when `n0 == 0`
 * (`n0 - 1 == 0xFF`) and writes 8-byte unaligned chunks directly into
 * `out.parent` and `out.offset_m1`.
 */
inline void unpack_louds_swar_lut(const uint64_t words[8],
                                  LevelCsr256& out,
                                  size_t count) noexcept {
  size_t n1 = 0;
  size_t n0 = 0;
  for (size_t w = 0; w < 8; ++w) {
    uint64_t x = words[w];
    for (unsigned b = 0; b < 8; ++b) {
      const uint8_t m1 = static_cast<uint8_t>(x & 0xFFu);
      const uint8_t m0 = static_cast<uint8_t>(~m1);
      const __m128i base1 = _mm_set1_epi8(static_cast<char>(n0 - 1));
      const __m128i base0 = _mm_set1_epi8(static_cast<char>(n1 - 1));
      const __m128i v1 = _mm_add_epi8(_mm_cvtsi64_si128(static_cast<long long>(
                                          detail::kLoudsDiffLut.packed[m1])),
                                      base1);
      const __m128i v0 = _mm_add_epi8(_mm_cvtsi64_si128(static_cast<long long>(
                                          detail::kLoudsDiffLut.packed[m0])),
                                      base0);
      auto store = [](uint8_t* dst, size_t n, __m128i v) {
        if (n <= 248) {
          _mm_storeu_si64(dst + n, v);
        } else {
          uint64_t bytes = _mm_cvtsi128_si64(v);
          for (; n < 256; ++n, bytes >>= 8) {
            dst[n] = static_cast<uint8_t>(bytes);
          }
        }
      };
      store(out.parent, n1, v1);
      store(out.offset_m1, n0, v0);
      const size_t c1 = static_cast<size_t>(std::popcount(m1));
      n1 += c1;
      n0 += 8 - c1;
      x >>= 8;
    }
  }
  if (count < 256) {
    out.parent[count] = out.parent[count - 1] ^ 1u;
  }
}
#endif  // __SSE2__

#ifdef PIXIE_FOREST256_AVX512_VBMI2
namespace detail {
// The two output streams each contain 256 bytes. Overlapping stores are fine
// within a stream: the next store replaces every unused compressed lane.
inline void store_louds_block(uint8_t (&dst)[256], size_t start, __m512i v) {
  const auto remaining =
      static_cast<unsigned>(std::min(size_t{64}, 256 - start));
  _mm512_mask_storeu_epi8(dst + start, ~mask_bits_from_u64(remaining), v);
}
}  // namespace detail
/**
 * @brief Two-stage AVX-512 VBMI2 unpacker for `LoudsForest256`.
 *
 * @details Compresses 1-bit and 0-bit positions modulo 256 into the output
 * arrays via `_mm512_maskz_compress_epi8`, then computes `parent[v] = (p1[v] -
 * (v + 1)) mod 256` and `offset_m1[u] = (p0[u] - (u + 1)) mod 256`.
 */
inline void unpack_louds_avx512(const uint64_t words[8],
                                LevelCsr256& out,
                                size_t count) noexcept {
  const __m512i iota0 = _mm512_setr_epi64(
      0x0706050403020100ULL, 0x0F0E0D0C0B0A0908ULL, 0x1716151413121110ULL,
      0x1F1E1D1C1B1A1918ULL, 0x2726252423222120ULL, 0x2F2E2D2C2B2A2928ULL,
      0x3736353433323130ULL, 0x3F3E3D3C3B3A3938ULL);
  const __m512i step64 = _mm512_set1_epi8(64);
  __m512i idx = iota0;
  size_t n1 = 0;
  size_t n0 = 0;
  for (size_t w = 0; w < 8; ++w) {
    const uint64_t m1 = words[w];
    const uint64_t m0 = ~m1;
    detail::store_louds_block(out.parent, n1,
                              _mm512_maskz_compress_epi8(m1, idx));
    detail::store_louds_block(out.offset_m1, n0,
                              _mm512_maskz_compress_epi8(m0, idx));
    const size_t c1 = static_cast<size_t>(std::popcount(m1));
    n1 += c1;
    n0 += 64 - c1;
    idx = _mm512_add_epi8(idx, step64);
  }
  __m512i sub = _mm512_add_epi8(iota0, _mm512_set1_epi8(1));
  for (size_t c = 0; c < 4; ++c) {
    const __m512i v1 = _mm512_load_si512(out.parent + 64 * c);
    const __m512i v0 = _mm512_load_si512(out.offset_m1 + 64 * c);
    _mm512_store_si512(out.parent + 64 * c, _mm512_sub_epi8(v1, sub));
    _mm512_store_si512(out.offset_m1 + 64 * c, _mm512_sub_epi8(v0, sub));
    sub = _mm512_add_epi8(sub, step64);
  }
  if (count < 256) {
    out.parent[count] = out.parent[count - 1] ^ 1u;
  }
}

/**
 * @brief Single-pass zero-scratch AVX-512 VBMI2 unpacker for `LoudsForest256`.
 *
 * @details Fuses `vpcompressb` and modular byte subtraction directly into
 * unaligned stores to `out.parent + n1` and `out.offset_m1 + n0`, avoiding all
 * temporary buffers and 64-byte reload forwarding stalls.
 */
inline void unpack_louds_fused_avx512(const uint64_t words[8],
                                      LevelCsr256& out,
                                      size_t count) noexcept {
  const __m512i iota0 = _mm512_setr_epi64(
      0x0706050403020100ULL, 0x0F0E0D0C0B0A0908ULL, 0x1716151413121110ULL,
      0x1F1E1D1C1B1A1918ULL, 0x2726252423222120ULL, 0x2F2E2D2C2B2A2928ULL,
      0x3736353433323130ULL, 0x3F3E3D3C3B3A3938ULL);
  size_t n1 = 0;
  size_t n0 = 0;
  for (size_t w = 0; w < 8; ++w) {
    const uint64_t m1 = words[w];
    const uint64_t m0 = ~m1;
    const __m512i z1 = _mm512_maskz_compress_epi8(m1, iota0);
    const __m512i z0 = _mm512_maskz_compress_epi8(m0, iota0);
    const __m512i sub1 =
        _mm512_add_epi8(iota0, _mm512_set1_epi8(static_cast<char>(1 - n0)));
    const __m512i sub0 =
        _mm512_add_epi8(iota0, _mm512_set1_epi8(static_cast<char>(1 - n1)));
    detail::store_louds_block(out.parent, n1, _mm512_sub_epi8(z1, sub1));
    detail::store_louds_block(out.offset_m1, n0, _mm512_sub_epi8(z0, sub0));
    const size_t c1 = static_cast<size_t>(std::popcount(m1));
    n1 += c1;
    n0 += 64 - c1;
  }
  if (count < 256) {
    out.parent[count] = out.parent[count - 1] ^ 1u;
  }
}
#endif  // PIXIE_FOREST256_AVX512_VBMI2

/**
 * @brief Unconditional dispatcher for `LoudsForest256` bulk unpacking.
 */
inline void unpack_louds(const LoudsForest256& encoded,
                         LevelCsr256& out,
                         size_t count) noexcept {
#if defined(PIXIE_FOREST256_AVX512_VBMI2)
  unpack_louds_avx512(encoded.words, out, count);
#elif defined(__SSE2__)
  unpack_louds_swar_lut(encoded.words, out, count);
#else
  unpack_louds_scalar(encoded.words, out, count);
#endif
}

/**
 * @brief Unconditional dispatcher preferring the zero-scratch fused LOUDS path.
 */
inline void unpack_louds_fused(const LoudsForest256& encoded,
                               LevelCsr256& out,
                               size_t count) noexcept {
#if defined(PIXIE_FOREST256_AVX512_VBMI2)
  unpack_louds_fused_avx512(encoded.words, out, count);
#elif defined(__SSE2__)
  unpack_louds_swar_lut(encoded.words, out, count);
#else
  unpack_louds_scalar(encoded.words, out, count);
#endif
}

// ============================================================================
// Twin-256 (H, L) Unpacking Kernels
// ============================================================================

/**
 * @brief Portable scalar reference unpacker for `TwinForest256`.
 */
inline void unpack_twin256_scalar(const uint64_t words[8],
                                  LevelCsr256& out) noexcept {
  const uint64_t* H = words;
  const uint64_t* L = words + 4;
  uint8_t u_ext[256]{};
  uint8_t e_pos[256]{};
  u_ext[0] = kNoParent;
  size_t u_len = 1;
  size_t e_len = 0;
  for (size_t w = 0; w < 4; ++w) {
    for (uint64_t x = H[w]; x != 0; x &= x - 1) {
      if (u_len < 256) {
        u_ext[u_len++] = static_cast<uint8_t>(64 * w + std::countr_zero(x));
      }
    }
    for (uint64_t x = L[w]; x != 0; x &= x - 1) {
      e_pos[e_len++] = static_cast<uint8_t>(64 * w + std::countr_zero(x));
    }
  }
  size_t g = 0;
  size_t r = 0;
  for (size_t w = 0; w < 4; ++w) {
    const uint64_t hw = H[w];
    const uint64_t lw = L[w];
    for (unsigned b = 0; b < 64; ++b) {
      const size_t idx = 64 * w + b;
      out.parent[idx] = u_ext[g];
      out.offset_m1[idx] = e_pos[r];
      g += static_cast<size_t>((lw >> b) & 1ULL);
      r += static_cast<size_t>((hw >> b) & 1ULL);
    }
  }
}

#ifdef PIXIE_AVX2_SUPPORT
/**
 * @brief AVX2 unpacker for `TwinForest256`.
 *
 * @details Streams compressed set-bit tables `U_ext` and `E` through a 128-bit
 * integer accumulator into 16-byte-aligned XMM stores, then exploits the
 * monotone 16-bit window rank bound (`0 <= rank1_excl_local <= 15`) so each
 * 16-node half-lane reads at most two aligned 16-byte table vectors resolved by
 * two `_mm256_shuffle_epi8` instructions with zero store-forwarding stalls.
 */
inline void unpack_twin256_avx2(const uint64_t words[8],
                                LevelCsr256& out) noexcept {
  const uint64_t* H = words;
  const uint64_t* L = words + 4;

  alignas(32) __m128i U_vec[16];
  alignas(32) __m128i E_vec[16];

  auto compress_256_aligned16 = [](const uint64_t m[4], __m128i* dst16,
                                   uint64_t init_val, unsigned init_fill) {
    __uint128_t buf = init_val;
    unsigned fill = init_fill;
    uint64_t lo_half = 0;
    unsigned phase = 0;
    size_t out_idx = 0;

    for (size_t w = 0; w < 4; ++w) {
      uint64_t x = m[w];
      if (x == 0) {
        continue;
      }
      uint64_t base_vec = 0x0101010101010101ULL * (64 * w);
      for (unsigned b = 0; b < 8; ++b) {
        const unsigned byte = static_cast<unsigned>(x & 0xFFu);
        const unsigned cnt = static_cast<unsigned>(std::popcount(byte));
        const uint64_t p = detail::mask_low_bytes_u64(
            detail::kByteCompressLut.packed[byte] + base_vec, cnt);
        buf |= (static_cast<__uint128_t>(p) << (fill * 8));
        fill += cnt;
        if (fill >= 8) {
          const uint64_t chunk = static_cast<uint64_t>(buf);
          buf >>= 64;
          fill -= 8;
          if (phase == 0) {
            lo_half = chunk;
            phase = 1;
          } else {
            _mm_store_si128(dst16 + out_idx++,
                            _mm_set_epi64x(static_cast<long long>(chunk),
                                           static_cast<long long>(lo_half)));
            phase = 0;
          }
        }
        x >>= 8;
        base_vec += 0x0808080808080808ULL;
      }
    }
    const uint64_t tail0 = static_cast<uint64_t>(buf);
    const uint64_t tail1 = static_cast<uint64_t>(buf >> 64);
    if (out_idx == 16) {
      return;
    }
    if (phase == 0) {
      _mm_store_si128(dst16 + out_idx++,
                      _mm_set_epi64x(static_cast<long long>(tail1),
                                     static_cast<long long>(tail0)));
      if (out_idx < 16) {
        _mm_store_si128(dst16 + out_idx, _mm_setzero_si128());
      }
    } else {
      _mm_store_si128(dst16 + out_idx++,
                      _mm_set_epi64x(static_cast<long long>(tail0),
                                     static_cast<long long>(lo_half)));
      if (out_idx < 16) {
        _mm_store_si128(dst16 + out_idx,
                        _mm_set_epi64x(0, static_cast<long long>(tail1)));
      }
    }
  };

  compress_256_aligned16(H, U_vec, 0xFFULL, 1);
  compress_256_aligned16(L, E_vec, 0ULL, 0);

  const __m256i kRepBytes32 =
      _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                       2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
  const __m256i kBitPref = _mm256_set1_epi64x(0x7F3F1F0F07030100ULL);
  const __m256i mask0f = _mm256_set1_epi8(0x0F);
  const __m256i x70 = _mm256_set1_epi8(0x70);
  const __m256i x16 = _mm256_set1_epi8(16);

  auto expand_plane_aligned_avx2 = [&](const uint64_t bits[4],
                                       const __m128i* tbl16, uint8_t* dst) {
    size_t r = 0;
    for (size_t w = 0; w < 4; ++w) {
      const uint64_t word = bits[w];
      for (unsigned half = 0; half < 2; ++half) {
        const uint32_t h32 = static_cast<uint32_t>(word >> (32 * half));
        const uint16_t m0 = static_cast<uint16_t>(h32 & 0xFFFFu);
        const uint16_t m1 = static_cast<uint16_t>(h32 >> 16);
        const size_t r0 = r;
        const size_t r1 = r0 + static_cast<size_t>(std::popcount(m0));
        r = r1 + static_cast<size_t>(std::popcount(m1));

        const uint8_t off0 = static_cast<uint8_t>(r0 & 15u);
        const uint8_t off1 = static_cast<uint8_t>(r1 & 15u);
        const uint8_t c0 = static_cast<uint8_t>(
            off0 + std::popcount(static_cast<uint8_t>(m0)));
        const uint8_t c1 = static_cast<uint8_t>(
            off1 + std::popcount(static_cast<uint8_t>(m1)));
        const __m256i base_off = _mm256_setr_epi64x(
            0x0101010101010101LL * off0, 0x0101010101010101LL * c0,
            0x0101010101010101LL * off1, 0x0101010101010101LL * c1);

        const __m256i raw = _mm256_set1_epi32(static_cast<int>(h32));
        const __m256i curr =
            _mm256_and_si256(_mm256_shuffle_epi8(raw, kRepBytes32), kBitPref);
        const __m256i lo_nib = _mm256_and_si256(curr, mask0f);
        const __m256i hi_nib =
            _mm256_and_si256(_mm256_srli_epi16(curr, 4), mask0f);
        const __m256i idx = _mm256_add_epi8(
            base_off,
            _mm256_add_epi8(_mm256_shuffle_epi8(lookup_popcount_4, lo_nib),
                            _mm256_shuffle_epi8(lookup_popcount_4, hi_nib)));

        const size_t q0 = r0 >> 4;
        const size_t q1 = r1 >> 4;
        const __m256i lo256 = _mm256_set_m128i(_mm_load_si128(tbl16 + q1),
                                               _mm_load_si128(tbl16 + q0));
        // A rank fits in a byte. At block 15 no live index selects the
        // following block; wrap those unused high-lane loads to block zero.
        const __m256i hi256 =
            _mm256_set_m128i(_mm_load_si128(tbl16 + ((q1 + 1) & 15)),
                             _mm_load_si128(tbl16 + ((q0 + 1) & 15)));
        const __m256i res = _mm256_or_si256(
            _mm256_shuffle_epi8(lo256, _mm256_add_epi8(idx, x70)),
            _mm256_shuffle_epi8(hi256, _mm256_sub_epi8(idx, x16)));
        _mm256_store_si256(reinterpret_cast<__m256i*>(dst + 64 * w + 32 * half),
                           res);
      }
    }
  };

  expand_plane_aligned_avx2(L, U_vec, out.parent);
  expand_plane_aligned_avx2(H, E_vec, out.offset_m1);
}
#endif  // PIXIE_AVX2_SUPPORT

#ifdef PIXIE_FOREST256_AVX512_FULL
/**
 * @brief Branchless in-register AVX-512 unpacker for `TwinForest256`.
 *
 * @details Compresses `H` and `L` into ZMM register tables (`U`, `E`) with zero
 * stack stores, places `kNoParent` (`0xFF`) at `U[255]`, evaluates 64-lane
 * exclusive bit-ranks via `vpopcntq` + `vpshufb` + `vpopcntb`, and gathers
 * `parent[v] = U[(rank1_excl(L, v) - 1) mod 256]` and
 * `offset_m1[u] = E[rank1_excl(H, u)]` via `vpermt2b`.
 */
inline void unpack_twin256_avx512(const uint64_t words[8],
                                  LevelCsr256& out) noexcept {
  const uint64_t* H = words;
  const uint64_t* L = words + 4;

  detail::ZmmTable256 U = detail::compress_256_branchless_avx512(H);
  U.t3 = _mm512_mask_mov_epi8(U.t3, uint64_t{1} << 63, _mm512_set1_epi8(-1));
  const detail::ZmmTable256 E = detail::compress_256_branchless_avx512(L);

  const uint8_t l0 = 0;
  const uint8_t l1 = static_cast<uint8_t>(std::popcount(L[0]));
  const uint8_t l2 = static_cast<uint8_t>(l1 + std::popcount(L[1]));
  const uint8_t l3 = static_cast<uint8_t>(l2 + std::popcount(L[2]));

  _mm512_store_si512(
      out.parent + 0,
      detail::lookup_zmm_branchless_avx512(
          U, detail::rank1_excl_64_avx512(L[0], static_cast<uint8_t>(l0 - 1))));
  _mm512_store_si512(
      out.parent + 64,
      detail::lookup_zmm_branchless_avx512(
          U, detail::rank1_excl_64_avx512(L[1], static_cast<uint8_t>(l1 - 1))));
  _mm512_store_si512(
      out.parent + 128,
      detail::lookup_zmm_branchless_avx512(
          U, detail::rank1_excl_64_avx512(L[2], static_cast<uint8_t>(l2 - 1))));
  _mm512_store_si512(
      out.parent + 192,
      detail::lookup_zmm_branchless_avx512(
          U, detail::rank1_excl_64_avx512(L[3], static_cast<uint8_t>(l3 - 1))));

  const uint8_t h0 = 0;
  const uint8_t h1 = static_cast<uint8_t>(std::popcount(H[0]));
  const uint8_t h2 = static_cast<uint8_t>(h1 + std::popcount(H[1]));
  const uint8_t h3 = static_cast<uint8_t>(h2 + std::popcount(H[2]));

  _mm512_store_si512(out.offset_m1 + 0,
                     detail::lookup_zmm_branchless_avx512(
                         E, detail::rank1_excl_64_avx512(H[0], h0)));
  _mm512_store_si512(out.offset_m1 + 64,
                     detail::lookup_zmm_branchless_avx512(
                         E, detail::rank1_excl_64_avx512(H[1], h1)));
  _mm512_store_si512(out.offset_m1 + 128,
                     detail::lookup_zmm_branchless_avx512(
                         E, detail::rank1_excl_64_avx512(H[2], h2)));
  _mm512_store_si512(out.offset_m1 + 192,
                     detail::lookup_zmm_branchless_avx512(
                         E, detail::rank1_excl_64_avx512(H[3], h3)));
}
#endif  // PIXIE_FOREST256_AVX512_FULL

/**
 * @brief Unconditional dispatcher for `TwinForest256` bulk unpacking.
 */
inline void unpack_twin256(const TwinForest256& encoded,
                           LevelCsr256& out) noexcept {
#if defined(PIXIE_FOREST256_AVX512_FULL)
  unpack_twin256_avx512(encoded.words, out);
#elif defined(PIXIE_AVX2_SUPPORT)
  unpack_twin256_avx2(encoded.words, out);
#else
  unpack_twin256_scalar(encoded.words, out);
#endif
}

enum class Codec { louds, twin };

namespace detail {
// Caller checked the full slice. Read a second word only if it contributes
// actual bits, including when the final source word ends at a page boundary.
inline void extract_plane(std::span<const uint64_t> src,
                          unsigned shift,
                          size_t bits,
                          uint64_t* dst,
                          uint64_t padding) {
  while (bits != 0) {
    const size_t take = std::min(bits, size_t{64});
    uint64_t value = src[0] >> shift;
    if (take > 64 - shift) {
      value |= src[1] << (64 - shift);
    }
    const uint64_t mask = ~uint64_t{0} >> (64 - take);
    *dst++ = (value & mask) | (padding & ~mask);
    src = src.subspan(1);
    bits -= take;
  }
}

template <bool Scalar>
inline void decode_dense(Codec codec,
                         std::span<const uint64_t> src,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  if (count == 0 || count > 256) {
    throw std::invalid_argument("Forest256 count");
  }
  const size_t start = offset / 64;
  const unsigned shift = offset % 64;
  const size_t needed = (shift + 2 * count + 63) / 64;
  if (start > src.size() || needed > src.size() - start) {
    throw std::out_of_range("Forest256 packed span");
  }
  src = src.subspan(start, needed);
  if (codec == Codec::louds) {
    // The omitted last delimiter is at 2m (even). All subsequent odd bits
    // are ones, exactly as in the reference padded encoder.
    constexpr uint64_t pad = 0xaaaaaaaaaaaaaaaaULL;
    LoudsForest256 encoded;
    std::fill_n(encoded.words, 8, pad);
    extract_plane(src, shift, 2 * count, encoded.words, pad);
    // Both output streams must have exactly 256 entries.
    unsigned ones = 0;
    for (auto word : encoded.words) {
      ones += std::popcount(word);
    }
    if (ones != 256) {
      throw std::invalid_argument("Forest256 LOUDS population");
    }
    if constexpr (Scalar) {
      unpack_louds_scalar(encoded.words, out, count);
    } else {
      unpack_louds_fused(encoded, out, count);
    }
  } else {
    TwinForest256 encoded{};
    extract_plane(src, shift, count, encoded.words, 0);
    const size_t l = shift + count;
    extract_plane(src.subspan(l / 64), l % 64, count, encoded.words + 4, 0);
    if (count < 256) {
      // For 255 isolated roots, synthetic parent 255 would alias kNoParent.
      if (count != 255 || (encoded.words[0] | encoded.words[1] |
                           encoded.words[2] | encoded.words[3])) {
        set_bit(encoded.words, count);
      }
      set_bit(encoded.words + 4, 255);
    }
    if constexpr (Scalar) {
      unpack_twin256_scalar(encoded.words, out);
    } else {
      unpack_twin256(encoded, out);
    }
  }
}
}  // namespace detail

/// Decode a valid dense BFS forest of 1..256 nodes at an LSB-first bit offset.
/// LOUDS stores its first 2m bits; Twin stores m H bits followed by m L bits.
/// Checks count and packed extent, not semantic validity of arbitrary topology.
inline void unpack_dense(Codec codec,
                         std::span<const uint64_t> stream,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  detail::decode_dense<false>(codec, stream, offset, count, out);
}

/// Same dense input and extent checks, using only the scalar reference kernel.
inline void unpack_dense_scalar(Codec codec,
                                std::span<const uint64_t> stream,
                                size_t offset,
                                size_t count,
                                LevelCsr256& out) {
  detail::decode_dense<true>(codec, stream, offset, count, out);
}
}  // namespace pixie::tree::detail::forest256
