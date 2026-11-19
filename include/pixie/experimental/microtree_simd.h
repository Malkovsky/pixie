#pragma once

#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

namespace pixie::tree::experimental::microtree_simd {

/**
 * @brief 128-bit bitmask returned by 128-bit SIMD microtree kernels.
 */
struct Mask128 {
  uint64_t lo = 0;
  uint64_t hi = 0;
};

/**
 * @brief Read up to 64 bits starting at bit offset @p pos, zero-filling past
 * @p words.size().
 */
[[gnu::always_inline]] inline uint64_t load_bits_64(
    std::span<const uint64_t> words,
    size_t pos) noexcept {
  const size_t idx = pos >> 6;
  if (idx >= words.size()) {
    return 0;
  }
  const unsigned shift = static_cast<unsigned>(pos & 63);
  uint64_t lo = words[idx] >> shift;
  if (shift != 0 && idx + 1 < words.size()) {
    lo |= words[idx + 1] << (64 - shift);
  }
  return lo;
}

/**
 * @brief Read up to 128 bits starting at bit offset @p pos into @p out[0..1],
 * zero-filling past @p words.size().
 */
[[gnu::always_inline]] inline void load_bits_128(
    std::span<const uint64_t> words,
    size_t pos,
    uint64_t out[2]) noexcept {
  const size_t idx = pos >> 6;
  const size_t n = words.size();
  if (idx >= n) {
    out[0] = 0;
    out[1] = 0;
    return;
  }
  const unsigned shift = static_cast<unsigned>(pos & 63);
  const uint64_t w0 = words[idx];
  const uint64_t w1 = (idx + 1 < n) ? words[idx + 1] : 0;
  if (shift == 0) {
    out[0] = w0;
    out[1] = w1;
    return;
  }
  const uint64_t w2 = (idx + 2 < n) ? words[idx + 2] : 0;
  out[0] = (w0 >> shift) | (w1 << (64 - shift));
  out[1] = (w1 >> shift) | (w2 << (64 - shift));
}

/**
 * @brief Count consecutive BP leaf encodings ('10', little-endian 0b01) from
 * bit 0 of @p bits, up to 32 leaves.
 */
[[gnu::always_inline]] inline unsigned bp_leaf_run_64(uint64_t bits) noexcept {
  constexpr uint64_t kAlternatingLeaves = 0x5555555555555555ULL;
  return static_cast<unsigned>(std::countr_zero(bits ^ kAlternatingLeaves)) >>
         1;
}

namespace detail {

constexpr uint8_t nibble_record_lows_mask(unsigned nib, int gap) noexcept {
  int cur = 0;
  int best = -gap;
  uint8_t mask = 0;
  for (int b = 0; b < 4; ++b) {
    cur += ((nib >> b) & 1u) ? 1 : -1;
    if (cur < best) {
      best = cur;
      mask |= static_cast<uint8_t>(1u << b);
    }
  }
  return mask;
}

template <int Gap>
constexpr std::array<uint8_t, 16> make_nibble_rl_table() noexcept {
  std::array<uint8_t, 16> t{};
  for (unsigned nib = 0; nib < 16; ++nib) {
    t[nib] = nibble_record_lows_mask(nib, Gap);
  }
  return t;
}

inline constexpr auto kNibbleRlGap0 = make_nibble_rl_table<0>();
inline constexpr auto kNibbleRlGap1 = make_nibble_rl_table<1>();
inline constexpr auto kNibbleRlGap2 = make_nibble_rl_table<2>();
inline constexpr auto kNibbleRlGap3 = make_nibble_rl_table<3>();

inline constexpr std::array<int8_t, 16> kNibbleDelta = [] {
  std::array<int8_t, 16> t{};
  for (unsigned nib = 0; nib < 16; ++nib) {
    int cur = 0;
    for (int b = 0; b < 4; ++b) {
      cur += ((nib >> b) & 1u) ? 1 : -1;
    }
    t[nib] = static_cast<int8_t>(cur);
  }
  return t;
}();

inline constexpr std::array<int8_t, 16> kNibbleMin = [] {
  std::array<int8_t, 16> t{};
  for (unsigned nib = 0; nib < 16; ++nib) {
    int cur = 0;
    int best = 4;
    for (int b = 0; b < 4; ++b) {
      cur += ((nib >> b) & 1u) ? 1 : -1;
      best = std::min(best, cur);
    }
    t[nib] = static_cast<int8_t>(best);
  }
  return t;
}();

}  // namespace detail

/**
 * @brief Compute 32-bit bitmask of 2-bit pair boundaries where BP half-excess
 * equals zero across 64 bits (32 pairs).
 *
 * @details Partition @p bits into 32 consecutive 2-bit pairs P_0..P_31 where
 * P_j = (bit 2j, bit 2j+1) contributes half-delta (bit 2j + bit 2j+1 - 1) in
 * {-1, 0, +1}. Bit j of the return value is 1 iff sum_{t=0..j} half_delta(P_t)
 * == 0. When @p bits starts at `u.pos + 1` of a non-leaf BP node whose subtree
 * has `J + 1 <= 33` nodes (`J` interior pairs), the set bits of
 * `bp_pair_zero_mask_64(bits) & ((1ULL << J) - 1)` are the exact end-pair
 * indices of all immediate children of `u`.
 */
[[gnu::always_inline]] inline uint32_t bp_pair_zero_mask_64(
    uint64_t bits) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  const __m128i raw = _mm_cvtsi64_si128(static_cast<long long>(bits));
  const __m256i z = _mm256_cvtepu8_epi32(raw);
  const __m256i a = _mm256_or_si256(z, _mm256_slli_epi32(z, 12));
  const __m256i b = _mm256_slli_epi32(a, 6);
  const __m256i c =
      _mm256_or_si256(_mm256_and_si256(a, _mm256_set1_epi32(0x00030003)),
                      _mm256_and_si256(b, _mm256_set1_epi32(0x03000300)));
  const __m256i lut_hd =
      _mm256_setr_epi8(-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
                       0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  __m256i ps = _mm256_shuffle_epi8(lut_hd, c);
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 1));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 2));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 4));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 8));
  const __m256i carry_src = _mm256_permute2x128_si256(ps, ps, 0x08);
  const __m256i carry = _mm256_shuffle_epi8(
      carry_src, _mm256_setr_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, 15, 15, 15, 15, 15, 15,
                                  15, 15, 15, 15, 15, 15, 15, 15, 15, 15));
  ps = _mm256_add_epi8(ps, carry);
  return static_cast<uint32_t>(
      _mm256_movemask_epi8(_mm256_cmpeq_epi8(ps, _mm256_setzero_si256())));
#else
  uint32_t mask = 0;
  int h = 0;
  for (unsigned j = 0; j < 32; ++j) {
    const unsigned pair = static_cast<unsigned>((bits >> (2 * j)) & 3u);
    h += static_cast<int>(pair & 1u) + static_cast<int>((pair >> 1) & 1u) - 1;
    if (h == 0) {
      mask |= (1u << j);
    }
  }
  return mask;
#endif
}

/**
 * @brief Compute 64-bit bitmask of 2-bit pair boundaries where BP half-excess
 * equals zero across 128 bits (64 pairs), starting from @p initial_half_excess.
 *
 * @param s 2 little-endian uint64_t words (128 bits = 64 2-bit pairs).
 * @param initial_half_excess Starting half-excess prior to pair 0.
 * @param next_half_excess Optional output receiving half-excess after pair 63.
 * @return 64-bit mask where bit j is 1 iff H_j == 0.
 */
[[gnu::always_inline]] inline uint64_t bp_pair_zero_mask_128(
    const uint64_t s[2],
    int initial_half_excess = 0,
    int* next_half_excess = nullptr) noexcept {
  const int block_delta = std::popcount(s[0]) + std::popcount(s[1]) - 64;
  if (next_half_excess != nullptr) {
    *next_half_excess = initial_half_excess + block_delta;
  }
  if (initial_half_excess < -64 || initial_half_excess > 64) {
    return 0;
  }

#ifdef PIXIE_AVX512_SUPPORT
  const __m128i raw = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s));
  const __m512i z = _mm512_cvtepu8_epi32(raw);
  const __m512i a = _mm512_or_si512(z, _mm512_slli_epi32(z, 12));
  const __m512i b = _mm512_slli_epi32(a, 6);
  const __m512i pairs =
      _mm512_or_si512(_mm512_and_si512(a, _mm512_set1_epi32(0x00030003)),
                      _mm512_and_si512(b, _mm512_set1_epi32(0x03000300)));
  const __m128i lut16 =
      _mm_setr_epi8(-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  __m512i ps = _mm512_shuffle_epi8(_mm512_broadcast_i32x4(lut16), pairs);
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 1));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 2));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 4));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 8));

  const __m512i zero = _mm512_setzero_si512();
  const __m512i lane_tot = _mm512_shuffle_epi8(ps, _mm512_set1_epi8(15));
  __m512i t = _mm512_add_epi8(lane_tot, _mm512_alignr_epi64(lane_tot, zero, 6));
  t = _mm512_add_epi8(t, _mm512_alignr_epi64(t, zero, 4));
  const __m512i carry = _mm512_alignr_epi64(t, zero, 6);
  ps = _mm512_add_epi8(ps, carry);

  const __m512i target =
      _mm512_set1_epi8(static_cast<int8_t>(-initial_half_excess));
  return static_cast<uint64_t>(_mm512_cmpeq_epi8_mask(ps, target));
#elif defined(PIXIE_AVX2_SUPPORT)
  const __m256i lut_hd =
      _mm256_setr_epi8(-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
                       0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  const __m256i mask_lo_carry = _mm256_setr_epi8(
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 15,
      15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15);
  auto half64 = [&](uint64_t w, int init_h) -> uint32_t {
    if (init_h < -32 || init_h > 32) {
      return 0;
    }
    const __m128i raw = _mm_cvtsi64_si128(static_cast<long long>(w));
    const __m256i z = _mm256_cvtepu8_epi32(raw);
    const __m256i a = _mm256_or_si256(z, _mm256_slli_epi32(z, 12));
    const __m256i b = _mm256_slli_epi32(a, 6);
    const __m256i pairs =
        _mm256_or_si256(_mm256_and_si256(a, _mm256_set1_epi32(0x00030003)),
                        _mm256_and_si256(b, _mm256_set1_epi32(0x03000300)));
    __m256i ps = _mm256_shuffle_epi8(lut_hd, pairs);
    ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 1));
    ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 2));
    ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 4));
    ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 8));
    const __m256i carry_src = _mm256_permute2x128_si256(ps, ps, 0x08);
    ps = _mm256_add_epi8(ps, _mm256_shuffle_epi8(carry_src, mask_lo_carry));
    const __m256i target = _mm256_set1_epi8(static_cast<int8_t>(-init_h));
    return static_cast<uint32_t>(
        _mm256_movemask_epi8(_mm256_cmpeq_epi8(ps, target)));
  };
  const int lo_delta = std::popcount(s[0]) - 32;
  const uint64_t lo_mask = half64(s[0], initial_half_excess);
  const uint64_t hi_mask = half64(s[1], initial_half_excess + lo_delta);
  return lo_mask | (hi_mask << 32);
#else
  uint64_t mask = 0;
  int h = initial_half_excess;
  for (unsigned w = 0; w < 2; ++w) {
    const uint64_t bits = s[w];
    for (unsigned j = 0; j < 32; ++j) {
      const unsigned pair = static_cast<unsigned>((bits >> (2 * j)) & 3u);
      h += static_cast<int>(pair & 1u) + static_cast<int>((pair >> 1) & 1u) - 1;
      if (h == 0) {
        mask |= (1ULL << (32 * w + j));
      }
    }
  }
  return mask;
#endif
}

/**
 * @brief Single-pass BP sibling & parent-close scan over 128 bits (64 pairs)
 * when the parent's closing parenthesis is not known in advance.
 *
 * @details Starting at `u.pos + 1`, `u` closes at `u.pos + 1 + 2 * close_pair`
 * at the first pair J where `H_{J-1} == 0 && s[2*J] == 0`.
 *
 * @return true if `u` closes within these 64 pairs (`0 <= close_pair < 64`),
 * with `child_ends_mask` masked to `< close_pair`; false otherwise
 * (`child_ends_mask` holds all 64 pairs' `H_j == 0` flags and
 * `next_half_excess` is updated).
 */
[[gnu::always_inline]] inline bool bp_children_128(
    const uint64_t s[2],
    int initial_half_excess,
    uint64_t& child_ends_mask,
    size_t& close_pair,
    int& next_half_excess) noexcept {
  const int block_delta = std::popcount(s[0]) + std::popcount(s[1]) - 64;
  next_half_excess = initial_half_excess + block_delta;
  if (initial_half_excess < -64 || initial_half_excess > 64) {
    child_ends_mask = 0;
    return false;
  }

#ifdef PIXIE_AVX512_SUPPORT
  const __m128i raw = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s));
  const __m512i z = _mm512_cvtepu8_epi32(raw);
  const __m512i a = _mm512_or_si512(z, _mm512_slli_epi32(z, 12));
  const __m512i b = _mm512_slli_epi32(a, 6);
  const __m512i pairs =
      _mm512_or_si512(_mm512_and_si512(a, _mm512_set1_epi32(0x00030003)),
                      _mm512_and_si512(b, _mm512_set1_epi32(0x03000300)));
  const uint64_t s_even =
      static_cast<uint64_t>(_mm512_test_epi8_mask(pairs, _mm512_set1_epi8(1)));

  const __m128i lut16 =
      _mm_setr_epi8(-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  __m512i ps = _mm512_shuffle_epi8(_mm512_broadcast_i32x4(lut16), pairs);
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 1));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 2));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 4));
  ps = _mm512_add_epi8(ps, _mm512_bslli_epi128(ps, 8));

  const __m512i zero = _mm512_setzero_si512();
  const __m512i lane_tot = _mm512_shuffle_epi8(ps, _mm512_set1_epi8(15));
  __m512i t = _mm512_add_epi8(lane_tot, _mm512_alignr_epi64(lane_tot, zero, 6));
  t = _mm512_add_epi8(t, _mm512_alignr_epi64(t, zero, 4));
  ps = _mm512_add_epi8(ps, _mm512_alignr_epi64(t, zero, 6));

  const uint64_t eq0 = static_cast<uint64_t>(_mm512_cmpeq_epi8_mask(
      ps, _mm512_set1_epi8(static_cast<int8_t>(-initial_half_excess))));
  const uint64_t prev_eq0 =
      (eq0 << 1) | static_cast<uint64_t>(initial_half_excess == 0);
  const uint64_t close_cand = prev_eq0 & ~s_even;
  if (close_cand != 0) {
    close_pair = static_cast<size_t>(std::countr_zero(close_cand));
    child_ends_mask =
        (close_pair == 0) ? 0ULL : (eq0 & ((1ULL << close_pair) - 1));
    return true;
  }
  child_ends_mask = eq0;
  return false;
#elif defined(PIXIE_AVX2_SUPPORT)
  const uint64_t eq0 = bp_pair_zero_mask_128(s, initial_half_excess, nullptr);
  uint64_t s_even = 0;
#ifdef PIXIE_BMI2_SUPPORT
  s_even = _pext_u64(s[0], 0x5555555555555555ULL) |
           (_pext_u64(s[1], 0x5555555555555555ULL) << 32);
#else
  auto compact_even_32 = [](uint64_t x) -> uint64_t {
    x &= 0x5555555555555555ULL;
    x = (x | (x >> 1)) & 0x3333333333333333ULL;
    x = (x | (x >> 2)) & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x >> 4)) & 0x00FF00FF00FF00FFULL;
    x = (x | (x >> 8)) & 0x0000FFFF0000FFFFULL;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;
    return x;
  };
  s_even = compact_even_32(s[0]) | (compact_even_32(s[1]) << 32);
#endif
  const uint64_t prev_eq0 =
      (eq0 << 1) | static_cast<uint64_t>(initial_half_excess == 0);
  const uint64_t close_cand = prev_eq0 & ~s_even;
  if (close_cand != 0) {
    close_pair = static_cast<size_t>(std::countr_zero(close_cand));
    child_ends_mask =
        (close_pair == 0) ? 0ULL : (eq0 & ((1ULL << close_pair) - 1));
    return true;
  }
  child_ends_mask = eq0;
  return false;
#else
  uint64_t eq0 = 0;
  int h = initial_half_excess;
  for (unsigned w = 0; w < 2; ++w) {
    const uint64_t bits = s[w];
    for (unsigned j = 0; j < 32; ++j) {
      const size_t pair_idx = 32 * w + j;
      const unsigned pair = static_cast<unsigned>((bits >> (2 * j)) & 3u);
      if (h == 0 && (pair & 1u) == 0) {
        close_pair = pair_idx;
        child_ends_mask = eq0;
        return true;
      }
      h += static_cast<int>(pair & 1u) + static_cast<int>((pair >> 1) & 1u) - 1;
      if (h == 0) {
        eq0 |= (1ULL << pair_idx);
      }
    }
  }
  child_ends_mask = eq0;
  return false;
#endif
}

/**
 * @brief Compute 64-bit strict record-low bitmask for 64 bits (starting at
 * excess 0, running minimum 0).
 *
 * @details In DFUDS, starting at `delimiter + 1` of a degree-`d` node whose
 * children span `<= 64` bits, the `d` set bits are the exact inclusive subtree
 * end offsets of children `0 .. d - 1`.
 */
[[gnu::always_inline]] inline uint64_t dfuds_record_lows_64(
    uint64_t bits) noexcept {
#ifdef PIXIE_SSE41_SUPPORT
  const __m128i raw = _mm_cvtsi64_si128(static_cast<long long>(bits));
  const __m128i mask0f = _mm_set1_epi8(0x0F);
  const __m128i lo_nib = _mm_and_si128(raw, mask0f);
  const __m128i hi_nib = _mm_and_si128(_mm_srli_epi16(raw, 4), mask0f);
  const __m128i nibs = _mm_unpacklo_epi8(lo_nib, hi_nib);

  const __m128i vdelta = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleDelta.data()));
  const __m128i vmin = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleMin.data()));

  const __m128i tot = _mm_shuffle_epi8(vdelta, nibs);
  __m128i ps = tot;
  ps = _mm_add_epi8(ps, _mm_bslli_si128(ps, 1));
  ps = _mm_add_epi8(ps, _mm_bslli_si128(ps, 2));
  ps = _mm_add_epi8(ps, _mm_bslli_si128(ps, 4));
  ps = _mm_add_epi8(ps, _mm_bslli_si128(ps, 8));
  const __m128i excl_ps = _mm_sub_epi8(ps, tot);

  const __m128i cand = _mm_add_epi8(excl_ps, _mm_shuffle_epi8(vmin, nibs));
  const __m128i inf = _mm_set1_epi8(127);
  __m128i inc_min = cand;
  inc_min = _mm_min_epi8(inc_min, _mm_alignr_epi8(inc_min, inf, 15));
  inc_min = _mm_min_epi8(inc_min, _mm_alignr_epi8(inc_min, inf, 14));
  inc_min = _mm_min_epi8(inc_min, _mm_alignr_epi8(inc_min, inf, 12));
  inc_min = _mm_min_epi8(inc_min, _mm_alignr_epi8(inc_min, inf, 8));
  const __m128i min_before =
      _mm_min_epi8(_mm_alignr_epi8(inc_min, inf, 15), _mm_setzero_si128());

  const __m128i gap = _mm_sub_epi8(excl_ps, min_before);
  const __m128i lut0 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap0.data()));
  const __m128i lut1 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap1.data()));
  const __m128i lut2 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap2.data()));
  const __m128i lut3 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap3.data()));

  __m128i m = _mm_and_si128(_mm_cmpeq_epi8(gap, _mm_setzero_si128()),
                            _mm_shuffle_epi8(lut0, nibs));
  m = _mm_or_si128(m, _mm_and_si128(_mm_cmpeq_epi8(gap, _mm_set1_epi8(1)),
                                    _mm_shuffle_epi8(lut1, nibs)));
  m = _mm_or_si128(m, _mm_and_si128(_mm_cmpeq_epi8(gap, _mm_set1_epi8(2)),
                                    _mm_shuffle_epi8(lut2, nibs)));
  m = _mm_or_si128(m, _mm_and_si128(_mm_cmpeq_epi8(gap, _mm_set1_epi8(3)),
                                    _mm_shuffle_epi8(lut3, nibs)));

  const __m128i paired = _mm_maddubs_epi16(m, _mm_set1_epi16(0x1001));
  const __m128i packed = _mm_packus_epi16(paired, paired);
  return static_cast<uint64_t>(_mm_cvtsi128_si64(packed));
#else
  uint64_t mask = 0;
  int cur = 0;
  int best = 0;
  for (unsigned b = 0; b < 64; ++b) {
    cur += ((bits >> b) & 1ULL) ? 1 : -1;
    if (cur < best) {
      best = cur;
      mask |= (1ULL << b);
    }
  }
  return mask;
#endif
}

/**
 * @brief Compute 128-bit strict record-low bitmask for 128 bits, carrying the
 * non-negative `gap = current_excess - running_min` across chained blocks.
 */
[[gnu::always_inline]] inline Mask128 dfuds_record_lows_128(
    const uint64_t s[2],
    int initial_gap = 0,
    int* next_gap = nullptr) noexcept {
  const int block_delta = 2 * (std::popcount(s[0]) + std::popcount(s[1])) - 128;
  if (initial_gap >= 128) {
    if (next_gap != nullptr) {
      *next_gap = initial_gap + block_delta;
    }
    return {0, 0};
  }

#ifdef PIXIE_AVX2_SUPPORT
  const __m128i raw = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s));
  const __m256i mask0f = _mm256_set1_epi8(0x0F);
  const __m128i lo_raw = _mm_unpacklo_epi64(raw, raw);
  const __m128i hi_raw = _mm_unpackhi_epi64(raw, raw);
  const __m256i wide =
      _mm256_inserti128_si256(_mm256_castsi128_si256(lo_raw), hi_raw, 1);
  // Unpack each 128-bit lane's low 8 bytes into 16 nibbles (32 nibbles total).
  const __m256i lo_nib = _mm256_and_si256(wide, mask0f);
  const __m256i hi_nib = _mm256_and_si256(_mm256_srli_epi16(wide, 4), mask0f);
  const __m256i nibs = _mm256_unpacklo_epi8(lo_nib, hi_nib);

  const __m256i vdelta = _mm256_broadcastsi128_si256(_mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleDelta.data())));
  const __m256i vmin = _mm256_broadcastsi128_si256(_mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleMin.data())));

  const __m256i tot = _mm256_shuffle_epi8(vdelta, nibs);
  __m256i ps = tot;
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 1));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 2));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 4));
  ps = _mm256_add_epi8(ps, _mm256_bslli_epi128(ps, 8));
  const __m256i carry_mask = _mm256_setr_epi8(
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 15,
      15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15);
  const __m256i ps_carry =
      _mm256_shuffle_epi8(_mm256_permute2x128_si256(ps, ps, 0x08), carry_mask);
  ps = _mm256_add_epi8(ps, ps_carry);
  const __m256i excl_ps = _mm256_sub_epi8(ps, tot);

  const __m256i cand =
      _mm256_add_epi8(excl_ps, _mm256_shuffle_epi8(vmin, nibs));
  const __m256i inf = _mm256_set1_epi8(127);
  __m256i inc_min = cand;
  inc_min = _mm256_min_epi8(inc_min, _mm256_alignr_epi8(inc_min, inf, 15));
  inc_min = _mm256_min_epi8(inc_min, _mm256_alignr_epi8(inc_min, inf, 14));
  inc_min = _mm256_min_epi8(inc_min, _mm256_alignr_epi8(inc_min, inf, 12));
  inc_min = _mm256_min_epi8(inc_min, _mm256_alignr_epi8(inc_min, inf, 8));

  const __m128i lo_min128 = _mm256_castsi256_si128(inc_min);
  const int lo_block_min = static_cast<int8_t>(_mm_extract_epi8(lo_min128, 15));
  const __m256i hi_carry_min = _mm256_inserti128_si256(
      inf, _mm_set1_epi8(static_cast<int8_t>(lo_block_min)), 1);
  inc_min = _mm256_min_epi8(inc_min, hi_carry_min);

  if (next_gap != nullptr) {
    const int hi_block_min = static_cast<int8_t>(
        _mm_extract_epi8(_mm256_extracti128_si256(inc_min, 1), 15));
    *next_gap = block_delta + std::max(initial_gap, -hi_block_min);
  }

  // Shift inc_min right by 1 nibble across the full 32-byte vector, inserting
  // +127 at nibble 0 and lo_block_min at nibble 16.
  const __m256i shift_fill = _mm256_inserti128_si256(
      inf, _mm_set1_epi8(static_cast<int8_t>(lo_block_min)), 1);
  const __m256i min_before =
      _mm256_min_epi8(_mm256_alignr_epi8(inc_min, shift_fill, 15),
                      _mm256_set1_epi8(static_cast<int8_t>(-initial_gap)));

  const __m256i gap = _mm256_sub_epi8(excl_ps, min_before);
  const __m128i l0 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap0.data()));
  const __m128i l1 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap1.data()));
  const __m128i l2 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap2.data()));
  const __m128i l3 = _mm_loadu_si128(
      reinterpret_cast<const __m128i*>(detail::kNibbleRlGap3.data()));
  const __m256i lut0 = _mm256_broadcastsi128_si256(l0);
  const __m256i lut1 = _mm256_broadcastsi128_si256(l1);
  const __m256i lut2 = _mm256_broadcastsi128_si256(l2);
  const __m256i lut3 = _mm256_broadcastsi128_si256(l3);

  __m256i m = _mm256_and_si256(_mm256_cmpeq_epi8(gap, _mm256_setzero_si256()),
                               _mm256_shuffle_epi8(lut0, nibs));
  m = _mm256_or_si256(
      m, _mm256_and_si256(_mm256_cmpeq_epi8(gap, _mm256_set1_epi8(1)),
                          _mm256_shuffle_epi8(lut1, nibs)));
  m = _mm256_or_si256(
      m, _mm256_and_si256(_mm256_cmpeq_epi8(gap, _mm256_set1_epi8(2)),
                          _mm256_shuffle_epi8(lut2, nibs)));
  m = _mm256_or_si256(
      m, _mm256_and_si256(_mm256_cmpeq_epi8(gap, _mm256_set1_epi8(3)),
                          _mm256_shuffle_epi8(lut3, nibs)));

  const __m256i paired = _mm256_maddubs_epi16(m, _mm256_set1_epi16(0x1001));
  const __m128i bytes_out = _mm_packus_epi16(
      _mm256_castsi256_si128(paired), _mm256_extracti128_si256(paired, 1));
  Mask128 res;
  _mm_storeu_si128(reinterpret_cast<__m128i*>(&res), bytes_out);
  return res;
#else
  Mask128 res{0, 0};
  int cur = 0;
  int best = -initial_gap;
  for (unsigned b = 0; b < 64; ++b) {
    cur += ((s[0] >> b) & 1ULL) ? 1 : -1;
    if (cur < best) {
      best = cur;
      res.lo |= (1ULL << b);
    }
  }
  for (unsigned b = 0; b < 64; ++b) {
    cur += ((s[1] >> b) & 1ULL) ? 1 : -1;
    if (cur < best) {
      best = cur;
      res.hi |= (1ULL << b);
    }
  }
  if (next_gap != nullptr) {
    *next_gap = cur - best;
  }
  return res;
#endif
}

}  // namespace pixie::tree::experimental::microtree_simd
