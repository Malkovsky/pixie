#pragma once

#include <pixie/bits.h>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace pixie::experimental {

#ifdef PIXIE_AVX2_SUPPORT
// clang-format off
static inline const __m256i excess_branch_lut_em4 = _mm256_setr_epi8(
    0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em3 = _mm256_setr_epi8(
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em2 = _mm256_setr_epi8(
    0x02, 0x08, 0x08, 0x00, 0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00,
    0x02, 0x08, 0x08, 0x00, 0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em1 = _mm256_setr_epi8(
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00);

static inline const __m256i excess_branch_lut_e0 = _mm256_setr_epi8(
    0x00, 0x02, 0x02, 0x08, 0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00, 0x08, 0x02, 0x02, 0x00,
    0x00, 0x02, 0x02, 0x08, 0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00, 0x08, 0x02, 0x02, 0x00);

static inline const __m256i excess_branch_lut_e1 = _mm256_setr_epi8(
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01);

static inline const __m256i excess_branch_lut_e2 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A, 0x00, 0x08, 0x08, 0x02,
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A, 0x00, 0x08, 0x08, 0x02);

static inline const __m256i excess_branch_lut_e3 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04);

static inline const __m256i excess_branch_lut_e4 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08);
// clang-format on

static inline __m256i excess_bit_masks_16x() noexcept {
  return _mm256_setr_epi16(0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020,
                           0x0040, 0x0080, 0x0100, 0x0200, 0x0400, 0x0800,
                           0x1000, 0x2000, 0x4000, (int16_t)0x8000);
}

static inline __m256i excess_prefix_sum_16x_i16(__m256i v) noexcept {
  __m256i x = v;
  __m256i t = _mm256_slli_si256(x, 2);
  x = _mm256_add_epi16(x, t);
  t = _mm256_slli_si256(x, 4);
  x = _mm256_add_epi16(x, t);
  t = _mm256_slli_si256(x, 8);
  x = _mm256_add_epi16(x, t);

  __m128i lo = _mm256_extracti128_si256(x, 0);
  __m128i hi = _mm256_extracti128_si256(x, 1);
  const int16_t carry = (int16_t)_mm_extract_epi16(lo, 7);
  hi = _mm_add_epi16(hi, _mm_set1_epi16(carry));

  __m256i out = _mm256_castsi128_si256(lo);
  out = _mm256_inserti128_si256(out, hi, 1);
  return out;
}

static inline int16_t excess_last_prefix_16x_i16(__m256i pref) noexcept {
  __m128i hi = _mm256_extracti128_si256(pref, 1);
  return (int16_t)_mm_extract_epi16(hi, 7);
}

static inline __m256i excess_bit_masks_32x8() noexcept {
  return _mm256_setr_epi8(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80);
}

static inline __m256i excess_byte_selectors_32x8() noexcept {
  return _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                          2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
}

static inline __m256i excess_prefix_sum_32x_i8(__m256i v) noexcept {
  __m256i x = v;
  __m256i t = _mm256_slli_si256(x, 1);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 2);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 4);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 8);
  x = _mm256_add_epi8(x, t);

  __m128i lo = _mm256_extracti128_si256(x, 0);
  __m128i hi = _mm256_extracti128_si256(x, 1);
  const int8_t carry = (int8_t)_mm_extract_epi8(lo, 15);
  hi = _mm_add_epi8(hi, _mm_set1_epi8(carry));

  __m256i out = _mm256_castsi128_si256(lo);
  out = _mm256_inserti128_si256(out, hi, 1);
  return out;
}

static inline int8_t excess_last_prefix_32x_i8(__m256i pref) noexcept {
  __m128i hi = _mm256_extracti128_si256(pref, 1);
  return (int8_t)_mm_extract_epi8(hi, 15);
}

static inline void excess_positions_512_branching_lut(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  int cur = 0;
  const __m256i vdelta =
      _mm256_setr_epi8(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -4,
                       -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4);
  const __m256i vmult = _mm256_set1_epi16(0x1001);
  const __m128i vnibble_mask = _mm_set1_epi8(0x0F);

  for (int k = 0; k < 4; ++k) {
    __m128i word_vec = _mm_loadu_si128((const __m128i*)&s[2 * k]);
    __m128i lo_nibbles = _mm_and_si128(word_vec, vnibble_mask);
    __m128i hi_nibbles =
        _mm_and_si128(_mm_srli_epi16(word_vec, 4), vnibble_mask);

    __m128i unpack_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
    __m128i unpack_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);
    __m256i nibbles = _mm256_inserti128_si256(_mm256_castsi128_si256(unpack_lo),
                                              unpack_hi, 1);

    __m256i ps = _mm256_shuffle_epi8(vdelta, nibbles);
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 1));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 2));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 4));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 8));

    __m128i ps_lo = _mm256_castsi256_si128(ps);
    __m128i ps_hi = _mm256_extracti128_si256(ps, 1);
    __m128i carry = _mm_set1_epi8((int8_t)_mm_extract_epi8(ps_lo, 15));
    ps_hi = _mm_add_epi8(ps_hi, carry);
    ps = _mm256_inserti128_si256(_mm256_castsi128_si256(ps_lo), ps_hi, 1);

    __m256i b = _mm256_permute2x128_si256(ps, ps, 0x08);
    __m256i excl_ps = _mm256_alignr_epi8(ps, b, 15);

    int target_rel = target_x - cur;
    int block_delta =
        2 * (std::popcount(s[2 * k]) + std::popcount(s[2 * k + 1])) - 128;

    const int d = 2 * target_rel - block_delta;
    if (d < -128 || d > 128) {
      cur += block_delta;
      continue;
    }

    if (target_rel == 128 || target_rel == -128) {
      out[2 * k + 1] |= (uint64_t{1} << 63);
      cur += block_delta;
      continue;
    }

    __m256i t = _mm256_sub_epi8(_mm256_set1_epi8((int8_t)target_rel), excl_ps);
    __m256i total_match = _mm256_setzero_si256();
    __m256i t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-4));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em4, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-3));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em3, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-2));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em2, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-1));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em1, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(0));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e0, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(1));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e1, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(2));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e2, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(3));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e3, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(4));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e4, nibbles)));

    __m256i res = _mm256_maddubs_epi16(total_match, vmult);
    __m128i packed = _mm_packus_epi16(_mm256_castsi256_si128(res),
                                      _mm256_extracti128_si256(res, 1));
    _mm_storeu_si128((__m128i*)&out[2 * k], packed);

    cur += block_delta;
  }
}
#else
static inline void excess_positions_512_branching_lut(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  excess_positions_512(s, target_x, out);
}
#endif

#ifdef PIXIE_AVX512_SUPPORT
static inline __m512i excess_lut_delta_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4));
}

static inline __m512i excess_lut_pos0_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1));
}

static inline __m512i excess_lut_pos1_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-2, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 2));
}

static inline __m512i excess_lut_pos2_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-3, -1, -1, 1, -1, 1, 1, 3, -3, -1, -1, 1, -1, 1, 1, 3));
}

static inline __m512i excess_bit_masks_64x8() noexcept {
  alignas(64) static constexpr int8_t masks[64] = {
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80};
  return _mm512_load_si512((const void*)masks);
}

static inline __m512i excess_byte_selectors_64x8() noexcept {
  alignas(64) static constexpr int8_t selectors[64] = {
      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
      2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5,
      5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7};
  return _mm512_load_si512((const void*)selectors);
}

static inline __m512i excess_prefix_sum_64x_i8(__m512i v) noexcept {
  __m512i x = v;
  __m512i t = _mm512_bslli_epi128(x, 1);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 2);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 4);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 8);
  x = _mm512_add_epi8(x, t);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(x, last_byte);
  const __m512i shift1_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11);
  const __m512i shift2_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7);
  const __m512i shift3_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3);

  __m512i lane_base =
      _mm512_maskz_permutexvar_epi32(0xFFF0, shift1_idx, lane_carry);
  lane_base = _mm512_add_epi8(lane_base, _mm512_maskz_permutexvar_epi32(
                                             0xFF00, shift2_idx, lane_carry));
  lane_base = _mm512_add_epi8(lane_base, _mm512_maskz_permutexvar_epi32(
                                             0xF000, shift3_idx, lane_carry));
  return _mm512_add_epi8(x, lane_base);
}

static inline __m512i excess_prefix_sum_2x32_i8(__m512i v) noexcept {
  __m512i x = v;
  __m512i t = _mm512_bslli_epi128(x, 1);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 2);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 4);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 8);
  x = _mm512_add_epi8(x, t);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(x, last_byte);
  const __m512i prev_lane_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8);
  const __mmask16 carry_to_second_lane_of_each_half = 0xF0F0;
  const __m512i lane_base = _mm512_maskz_permutexvar_epi32(
      carry_to_second_lane_of_each_half, prev_lane_idx, lane_carry);
  return _mm512_add_epi8(x, lane_base);
}

static inline __m512i excess_nibbles_64x_from_256(__m256i words) noexcept {
  const __m512i bytes16 = _mm512_cvtepu8_epi16(words);
  const __m512i low = _mm512_and_si512(bytes16, _mm512_set1_epi16(0x000F));
  const __m512i high = _mm512_and_si512(_mm512_srli_epi16(bytes16, 4),
                                        _mm512_set1_epi16(0x000F));
  return _mm512_or_si512(low, _mm512_slli_epi16(high, 8));
}

static inline __m512i excess_exclusive_prefix_2x32_i8(__m512i pref) noexcept {
  const __m512i zero = _mm512_setzero_si512();
  __m512i out = _mm512_alignr_epi8(pref, zero, 15);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(pref, last_byte);
  const __m512i prev_lane_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8);
  const __m512i carry_dwords =
      _mm512_permutexvar_epi32(prev_lane_idx, lane_carry);
  const __mmask64 first_byte_of_second_lane_in_each_half =
      (uint64_t{1} << 16) | (uint64_t{1} << 48);
  return _mm512_or_si512(
      out, _mm512_maskz_mov_epi8(first_byte_of_second_lane_in_each_half,
                                 carry_dwords));
}

static inline uint64_t excess_repeat_byte(int value) noexcept {
  return uint64_t{0x0101010101010101} *
         static_cast<uint8_t>(static_cast<int8_t>(value));
}

static inline void excess_positions_512_lut_avx512(const uint64_t* s,
                                                   int target_x,
                                                   uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  static const __m512i vdelta = excess_lut_delta_64x();
  static const __m512i vpos0 = excess_lut_pos0_64x();
  static const __m512i vpos1 = excess_lut_pos1_64x();
  static const __m512i vpos2 = excess_lut_pos2_64x();
  static const __m512i vbit0 = _mm512_set1_epi8(1);
  static const __m512i vbit1 = _mm512_set1_epi8(2);
  static const __m512i vbit2 = _mm512_set1_epi8(4);
  static const __m512i vbit3 = _mm512_set1_epi8(8);
  static const __m512i vmult = _mm512_set1_epi16(0x1001);

  for (int k = 0; k < 2; ++k) {
    const int base_word = 4 * k;
    const int delta0 =
        2 * (std::popcount(s[base_word]) + std::popcount(s[base_word + 1])) -
        128;
    const int delta1 = 2 * (std::popcount(s[base_word + 2]) +
                            std::popcount(s[base_word + 3])) -
                       128;
    const int target0 = target_x;
    const int target1 = target_x - delta0;
    const bool reachable0 = [&] {
      const int d = 2 * target0 - delta0;
      return -128 <= d && d <= 128;
    }();
    const bool reachable1 = [&] {
      const int d = 2 * target1 - delta1;
      return -128 <= d && d <= 128;
    }();

    if (!reachable0 && !reachable1) {
      target_x -= delta0 + delta1;
      continue;
    }

    const __m256i words =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&s[base_word]));
    const __m512i nibbles = excess_nibbles_64x_from_256(words);
    const __m512i ps =
        excess_prefix_sum_2x32_i8(_mm512_shuffle_epi8(vdelta, nibbles));
    const __m512i excl_ps = excess_exclusive_prefix_2x32_i8(ps);
    const uint64_t repeated0 = excess_repeat_byte(target0);
    const uint64_t repeated1 = excess_repeat_byte(target1);
    const __m512i vtgt =
        _mm512_setr_epi64(repeated0, repeated0, repeated0, repeated0, repeated1,
                          repeated1, repeated1, repeated1);
    const __m512i t = _mm512_sub_epi8(vtgt, excl_ps);

    const __mmask64 cmp0 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos0, nibbles), t);
    const __mmask64 cmp1 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos1, nibbles), t);
    const __mmask64 cmp2 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos2, nibbles), t);
    const __mmask64 cmp3 = _mm512_cmpeq_epi8_mask(ps, vtgt);
    __m512i total_match = _mm512_maskz_mov_epi8(cmp0, vbit0);
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp1, vbit1));
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp2, vbit2));
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp3, vbit3));

    const __mmask64 active =
        (reachable0 ? __mmask64{0x00000000FFFFFFFFull} : __mmask64{0}) |
        (reachable1 ? __mmask64{0xFFFFFFFF00000000ull} : __mmask64{0});
    total_match = _mm512_maskz_mov_epi8(active, total_match);

    const __m512i res = _mm512_maddubs_epi16(total_match, vmult);
    const __m256i packed = _mm512_cvtepi16_epi8(res);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&out[base_word]), packed);

    target_x -= delta0 + delta1;
  }
}
#else
static inline void excess_positions_512_lut_avx512(const uint64_t* s,
                                                   int target_x,
                                                   uint64_t* out) noexcept {
  excess_positions_512(s, target_x, out);
}
#endif

static inline void excess_positions_512_expand(const uint64_t* s,
                                               int target_x,
                                               uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX2_SUPPORT
  static const __m256i masks = excess_bit_masks_16x();
  static const __m256i vzero = _mm256_setzero_si256();
  static const __m256i vallones = _mm256_cmpeq_epi16(vzero, vzero);
  static const __m256i vminus1 = _mm256_set1_epi16(-1);
  static const __m256i vtwo = _mm256_set1_epi16(2);
  const __m256i vtarget = _mm256_set1_epi16((int16_t)target_x);

  int cur = 0;
  for (int block = 0; block < 4; ++block) {
    const int target_rel = target_x - cur;
    if (target_rel <= -64 || target_rel >= 64) {
      const int block_delta =
          2 * (std::popcount(s[2 * block]) + std::popcount(s[2 * block + 1])) -
          128;
      const int reachability = 2 * target_rel - block_delta;
      if (reachability < -128 || reachability > 128) {
        cur += block_delta;
        continue;
      }
    }

    for (int j = 0; j < 8; ++j) {
      const int k = 8 * block + j;
      const size_t word_idx = size_t(k) >> 2;
      const size_t shift = size_t(k & 3) * 16;
      const uint16_t bits16 =
          static_cast<uint16_t>((s[word_idx] >> shift) & 0xFFFFull);

      const __m256i vb = _mm256_set1_epi16((int16_t)bits16);
      const __m256i m = _mm256_and_si256(vb, masks);
      const __m256i is_zero = _mm256_cmpeq_epi16(m, vzero);
      const __m256i is_set = _mm256_andnot_si256(is_zero, vallones);
      const __m256i steps =
          _mm256_add_epi16(vminus1, _mm256_and_si256(is_set, vtwo));

      const __m256i pref_rel = excess_prefix_sum_16x_i16(steps);
      const __m256i base = _mm256_set1_epi16((int16_t)cur);
      const __m256i pref_abs = _mm256_add_epi16(pref_rel, base);
      const __m256i cmp = _mm256_cmpeq_epi16(pref_abs, vtarget);

      const uint32_t m32 = (uint32_t)_mm256_movemask_epi8(cmp);
      const uint16_t m16 = (uint16_t)_pext_u32(m32, 0xAAAAAAAAu);

      out[word_idx] |= uint64_t(m16) << shift;
      cur += (int)excess_last_prefix_16x_i16(pref_rel);
    }
  }
#else
  int cur = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
#endif
}

static inline void excess_positions_512_expand8(const uint64_t* s,
                                                int target_x,
                                                uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX2_SUPPORT
  static const __m256i byte_selectors = excess_byte_selectors_32x8();
  static const __m256i masks = excess_bit_masks_32x8();
  static const __m256i vzero = _mm256_setzero_si256();
  static const __m256i vallones = _mm256_cmpeq_epi8(vzero, vzero);
  static const __m256i vminus1 = _mm256_set1_epi8(-1);
  static const __m256i vtwo = _mm256_set1_epi8(2);

  int cur = 0;
  for (int k = 0; k < 16; ++k) {
    const size_t word_idx = size_t(k) >> 1;
    const size_t shift = size_t(k & 1) * 32;
    const uint32_t bits32 =
        static_cast<uint32_t>((s[word_idx] >> shift) & 0xFFFFFFFFull);

    const int target_rel = target_x - cur;
    if (target_rel < -32 || target_rel > 32) {
      cur += 2 * static_cast<int>(std::popcount(bits32)) - 32;
      continue;
    }

    const __m256i src = _mm256_set1_epi32((int)bits32);
    const __m256i bytes = _mm256_shuffle_epi8(src, byte_selectors);
    const __m256i m = _mm256_and_si256(bytes, masks);
    const __m256i is_zero = _mm256_cmpeq_epi8(m, vzero);
    const __m256i is_set = _mm256_andnot_si256(is_zero, vallones);
    const __m256i steps =
        _mm256_add_epi8(vminus1, _mm256_and_si256(is_set, vtwo));

    const __m256i pref_rel = excess_prefix_sum_32x_i8(steps);
    const __m256i vtarget = _mm256_set1_epi8((int8_t)target_rel);
    const __m256i cmp = _mm256_cmpeq_epi8(pref_rel, vtarget);
    const uint32_t mask = static_cast<uint32_t>(_mm256_movemask_epi8(cmp));

    out[word_idx] |= uint64_t(mask) << shift;
    cur += static_cast<int>(excess_last_prefix_32x_i8(pref_rel));
  }
#else
  int cur = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
#endif
}

static inline void excess_positions_512_expand_avx512(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX512_SUPPORT
  static const __m512i byte_selectors = excess_byte_selectors_64x8();
  static const __m512i masks = excess_bit_masks_64x8();
  static const __m512i vzero = _mm512_setzero_si512();
  static const __m512i vallones = _mm512_set1_epi8(-1);
  static const __m512i vminus1 = _mm512_set1_epi8(-1);
  static const __m512i vtwo = _mm512_set1_epi8(2);

  int cur = 0;
  for (int k = 0; k < 8; ++k) {
    const uint64_t bits64 = s[k];
    const int target_rel = target_x - cur;
    if (target_rel < -64 || target_rel > 64) {
      cur += 2 * static_cast<int>(std::popcount(bits64)) - 64;
      continue;
    }

    const __m512i src = _mm512_set1_epi64(static_cast<int64_t>(bits64));
    const __m512i bytes = _mm512_shuffle_epi8(src, byte_selectors);
    const __m512i m = _mm512_and_si512(bytes, masks);
    const __mmask64 is_zero = _mm512_cmpeq_epi8_mask(m, vzero);
    const __m512i is_set = _mm512_maskz_mov_epi8(~is_zero, vallones);
    const __m512i steps =
        _mm512_add_epi8(vminus1, _mm512_and_si512(is_set, vtwo));

    const __m512i pref_rel = excess_prefix_sum_64x_i8(steps);
    const __mmask64 match =
        _mm512_cmpeq_epi8_mask(pref_rel, _mm512_set1_epi8((int8_t)target_rel));
    out[k] = static_cast<uint64_t>(match);
    cur += 2 * static_cast<int>(std::popcount(bits64)) - 64;
  }
#else
  excess_positions_512_expand8(s, target_x, out);
#endif
}

}  // namespace pixie::experimental
