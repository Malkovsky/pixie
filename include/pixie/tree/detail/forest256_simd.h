#pragma once

/**
 * @file forest256_simd.h
 * @brief Decode LOUDS forest units into a sibling-contiguous CSR and parent
 * table.
 * @details A dense unit uses two bits per node for up to 256 breadth-first
 * nodes. The decoder pads the unit to 512 bits before running a bulk kernel.
 */

#include <pixie/bits.h>

#ifdef __SSE2__
#include <emmintrin.h>
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

[[gnu::always_inline]] inline uint64_t mask_bits_from_u64(unsigned k) noexcept {
#ifdef __BMI2__
  return ~_bzhi_u64(~uint64_t{0}, k);
#else
  return (k >= 64) ? 0ULL : (~uint64_t{0} << k);
#endif
}

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

#if defined(__AVX512BW__) && defined(__AVX512VBMI2__)
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
#endif

/**
 * @brief Unconditional dispatcher for `LoudsForest256` bulk unpacking.
 */
inline void unpack_louds(const LoudsForest256& encoded,
                         LevelCsr256& out,
                         size_t count) noexcept {
#if defined(__AVX512BW__) && defined(__AVX512VBMI2__)
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
#if defined(__AVX512BW__) && defined(__AVX512VBMI2__)
  unpack_louds_fused_avx512(encoded.words, out, count);
#elif defined(__SSE2__)
  unpack_louds_swar_lut(encoded.words, out, count);
#else
  unpack_louds_scalar(encoded.words, out, count);
#endif
}

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

inline std::span<const uint64_t> checked_dense_words(
    std::span<const uint64_t> src,
    size_t offset,
    size_t count) {
  if (count == 0 || count > 256) {
    throw std::invalid_argument("Forest256 count");
  }
  const size_t start = offset / 64;
  const unsigned shift = offset % 64;
  const size_t needed = (shift + 2 * count + 63) / 64;
  if (start > src.size() || needed > src.size() - start) {
    throw std::out_of_range("Forest256 packed span");
  }
  return src.subspan(start, needed);
}

template <bool Scalar>
inline void decode_dense(std::span<const uint64_t> src,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  src = checked_dense_words(src, offset, count);
  const unsigned shift = offset % 64;
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
}
}  // namespace detail

/// Decode a valid dense BFS forest of 1..256 nodes at an LSB-first bit offset.
/// Checks count and packed extent, not semantic validity of arbitrary topology.
inline void unpack_dense(std::span<const uint64_t> stream,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  detail::decode_dense<false>(stream, offset, count, out);
}

/// Same dense input and extent checks, using only the scalar reference kernel.
inline void unpack_dense_scalar(std::span<const uint64_t> stream,
                                size_t offset,
                                size_t count,
                                LevelCsr256& out) {
  detail::decode_dense<true>(stream, offset, count, out);
}
}  // namespace pixie::tree::detail::forest256
