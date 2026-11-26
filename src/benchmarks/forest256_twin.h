#pragma once

#include <pixie/tree/detail/forest256_cache.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

#if defined(__AVX512BW__) && defined(__AVX512VBMI__) &&      \
    defined(__AVX512VBMI2__) && defined(__AVX512BITALG__) && \
    defined(__AVX512VPOPCNTDQ__)
#define PIXIE_FOREST256_TWIN_AVX512
#endif

namespace pixie::benchmarks::forest256_twin {

using tree::detail::forest256::kNoParent;
using tree::detail::forest256::LevelCsr256;

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

using tree::detail::forest256::detail::checked_dense_words;
using tree::detail::forest256::detail::extract_plane;
using tree::detail::forest256::detail::mask_bits_from_u64;
using tree::detail::forest256::detail::set_bit;
using tree::detail::forest256::detail::validate_bfs_degrees;

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
#ifdef __BMI2__
  return _bzhi_u64(val, cnt * 8);
#else
  if (cnt == 0) {
    return 0;
  }
  return val & (~uint64_t{0} >> (64 - cnt * 8));
#endif
}

#ifdef PIXIE_FOREST256_TWIN_AVX512
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
#endif  // PIXIE_FOREST256_TWIN_AVX512

}  // namespace detail

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

#ifdef __AVX2__
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
#endif  // __AVX2__

#ifdef PIXIE_FOREST256_TWIN_AVX512
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
#endif  // PIXIE_FOREST256_TWIN_AVX512

/**
 * @brief Unconditional dispatcher for `TwinForest256` bulk unpacking.
 */
inline void unpack_twin256(const TwinForest256& encoded,
                           LevelCsr256& out) noexcept {
#if defined(PIXIE_FOREST256_TWIN_AVX512)
  unpack_twin256_avx512(encoded.words, out);
#elif defined(__AVX2__)
  unpack_twin256_avx2(encoded.words, out);
#else
  unpack_twin256_scalar(encoded.words, out);
#endif
}

namespace detail {

template <bool Scalar>
inline void decode_dense(std::span<const uint64_t> src,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  src = checked_dense_words(src, offset, count);
  const unsigned shift = offset % 64;
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
}  // namespace detail

/// Decode a dense Twin unit with count H bits followed by count L bits.
inline void unpack_dense(std::span<const uint64_t> stream,
                         size_t offset,
                         size_t count,
                         LevelCsr256& out) {
  detail::decode_dense<false>(stream, offset, count, out);
}

inline void unpack_dense_scalar(std::span<const uint64_t> stream,
                                size_t offset,
                                size_t count,
                                LevelCsr256& out) {
  detail::decode_dense<true>(stream, offset, count, out);
}

using Cache = tree::detail::forest256::Forest256CacheFor<unpack_dense>;

}  // namespace pixie::benchmarks::forest256_twin

#undef PIXIE_FOREST256_TWIN_AVX512
