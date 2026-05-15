#pragma once

#include <pixie/bits.h>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace pixie::experimental {

#ifdef PIXIE_AVX2_SUPPORT
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
#endif

#ifdef PIXIE_AVX512_SUPPORT
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
