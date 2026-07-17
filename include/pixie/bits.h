#pragma once

#include <immintrin.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>

#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__) && \
    defined(__AVX512BW__)
#define PIXIE_AVX512_SUPPORT
#endif

#if defined(__BMI__) && defined(__BMI2__)
#define PIXIE_BMI_SUPPORT
#endif

#ifdef __AVX2__
#define PIXIE_AVX2_SUPPORT
// Lookup table for 4-bit popcount
// This table maps each 4-bit value (0-15) to its population count
// clang-format off
static inline const __m256i lookup_popcount_4 = _mm256_setr_epi8(
    0, 1, 1, 2,  // 0000, 0001, 0010, 0011
    1, 2, 2, 3,  // 0100, 0101, 0110, 0111
    1, 2, 2, 3,  // 1000, 1001, 1010, 1011
    2, 3, 3, 4,  // 1100, 1101, 1110, 1111
    
    // Same table repeated for high 128 bits
    0, 1, 1, 2,  // 0000, 0001, 0010, 0011
    1, 2, 2, 3,  // 0100, 0101, 0110, 0111
    1, 2, 2, 3,  // 1000, 1001, 1010, 1011
    2, 3, 3, 4   // 1100, 1101, 1110, 1111
);

static inline const __m256i mask_first_half = _mm256_setr_epi8(
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
);

// clang-format on
#endif

static inline constexpr int8_t excess_nibble_min_offset[16] = {
    4, 4, 4, 4, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1};

#if defined(__SSSE3__) && defined(__SSE4_1__)
#define PIXIE_SSE41_SUPPORT
// clang-format off
static inline const __m128i excess_lut_delta_sse = _mm_setr_epi8(
    -4, -2, -2,  0,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
     0,  2,  2,  4);
static inline const __m128i excess_lut_pos0_sse = _mm_setr_epi8(
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1);
static inline const __m128i excess_lut_pos1_sse = _mm_setr_epi8(
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2);
static inline const __m128i excess_lut_pos2_sse = _mm_setr_epi8(
    -3, -1, -1,  1,
    -1,  1,  1,  3,
    -3, -1, -1,  1,
    -1,  1,  1,  3);
static inline const __m128i excess_lut_min_sse = _mm_setr_epi8(
    -4, -2, -2,  0,
    -2,  0, -1,  1,
    -3, -1, -1,  1,
    -2,  0, -1,  1);
static inline const __m128i excess_lut_nibble_index_sse = _mm_setr_epi8(
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12, 13, 14, 15);
static inline const __m128i excess_lut_low_nibble_index_sse = _mm_setr_epi8(
     0,  2,  4,  6,
     8, 10, 12, 14,
    16, 18, 20, 22,
    24, 26, 28, 30);
static inline const __m128i excess_lut_high_nibble_index_sse = _mm_setr_epi8(
     1,  3,  5,  7,
     9, 11, 13, 15,
    17, 19, 21, 23,
    25, 27, 29, 31);
static inline const __m128i excess_lut_nibble_mask_sse = _mm_set1_epi8(0x0F);
// clang-format on
#endif

#ifdef PIXIE_SSE41_SUPPORT
static inline __m128i excess_nibbles_64_sse(const uint64_t* s) noexcept {
  const __m128i word_vec = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(s));
  const __m128i lo_nibbles =
      _mm_and_si128(word_vec, excess_lut_nibble_mask_sse);
  const __m128i hi_nibbles =
      _mm_and_si128(_mm_srli_epi16(word_vec, 4), excess_lut_nibble_mask_sse);
  return _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
}

static inline __m128i excess_prefix_sum_16x_i8(__m128i v) noexcept {
  __m128i x = v;
  __m128i t = _mm_slli_si128(x, 1);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 2);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 4);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 8);
  return _mm_add_epi8(x, t);
}

static inline int excess_horizontal_min_i8(__m128i v) noexcept {
  v = _mm_min_epi8(v, _mm_alignr_epi8(v, v, 8));
  v = _mm_min_epi8(v, _mm_alignr_epi8(v, v, 4));
  v = _mm_min_epi8(v, _mm_alignr_epi8(v, v, 2));
  v = _mm_min_epi8(v, _mm_alignr_epi8(v, v, 1));
  return static_cast<int>(static_cast<int8_t>(_mm_extract_epi8(v, 0)));
}
#endif

/**
 * @brief Test 16 int16 RmM btree child ranges for a node-local target.
 * @details Each lane represents one child summary. The function checks whether
 * @p target lies in `[prefix_before[i] + min_excess[i],
 * prefix_before[i] + max_excess[i]]`. When @p include_zero_boundary is true,
 * it also accepts `target == prefix_before[i]`, which represents the left
 * boundary match used by backward search.
 * @param prefix_before Exclusive prefix excess before each child.
 * @param min_excess Per-child minimum excess relative to the child start.
 * @param max_excess Per-child maximum excess relative to the child start.
 * @param target Target excess relative to the start of the parent node.
 * @param include_zero_boundary Whether to accept child left-boundary matches.
 * @return Bit mask with bit `i` set when lane `i` can contain the target.
 */
static inline uint32_t rmm_btree_match_mask_i16x16(const int16_t* prefix_before,
                                                   const int16_t* min_excess,
                                                   const int16_t* max_excess,
                                                   int16_t target,
                                                   bool include_zero_boundary) {
#ifdef PIXIE_AVX2_SUPPORT
  const __m256i vtarget = _mm256_set1_epi16(target);
  const __m256i vprefix =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(prefix_before));
  const __m256i vmin =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(min_excess));
  const __m256i vmax =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(max_excess));

  const __m256i lower = _mm256_adds_epi16(vprefix, vmin);
  const __m256i upper = _mm256_adds_epi16(vprefix, vmax);
  const __m256i ge_lower = _mm256_or_si256(_mm256_cmpgt_epi16(vtarget, lower),
                                           _mm256_cmpeq_epi16(vtarget, lower));
  const __m256i le_upper = _mm256_or_si256(_mm256_cmpgt_epi16(upper, vtarget),
                                           _mm256_cmpeq_epi16(upper, vtarget));
  __m256i matched = _mm256_and_si256(ge_lower, le_upper);
  if (include_zero_boundary) {
    matched = _mm256_or_si256(matched, _mm256_cmpeq_epi16(vtarget, vprefix));
  }

  const uint32_t byte_mask =
      static_cast<uint32_t>(_mm256_movemask_epi8(matched));
  uint32_t result = 0;
  for (size_t lane = 0; lane < 16; ++lane) {
    const uint32_t lane_mask = 0x3u << (lane * 2);
    if ((byte_mask & lane_mask) == lane_mask) {
      result |= uint32_t{1} << lane;
    }
  }
  return result;
#else
  uint32_t result = 0;
  for (size_t lane = 0; lane < 16; ++lane) {
    const int lower = prefix_before[lane] + min_excess[lane];
    const int upper = prefix_before[lane] + max_excess[lane];
    const bool found = (lower <= target && target <= upper) ||
                       (include_zero_boundary && target == prefix_before[lane]);
    if (found) {
      result |= uint32_t{1} << lane;
    }
  }
  return result;
#endif
}

/**
 * @brief Test 4 int64 RmM btree child ranges for a node-local target.
 * @details Each lane represents one child summary. The function subtracts the
 * child prefix from @p target to form a child-relative target, then checks it
 * against the child's `[min_excess, max_excess]` range. When
 * @p include_zero_boundary is true, it also accepts a zero relative target for
 * the left-boundary match used by backward search.
 * @param prefix_before Exclusive prefix excess before each child.
 * @param min_excess Per-child minimum excess relative to the child start.
 * @param max_excess Per-child maximum excess relative to the child start.
 * @param target Target excess relative to the start of the parent node.
 * @param include_zero_boundary Whether to accept child left-boundary matches.
 * @return Bit mask with bit `i` set when lane `i` can contain the target.
 */
static inline uint32_t rmm_btree_match_mask_i64x4(const int64_t* prefix_before,
                                                  const int64_t* min_excess,
                                                  const int64_t* max_excess,
                                                  int64_t target,
                                                  bool include_zero_boundary) {
#ifdef PIXIE_AVX2_SUPPORT
  const __m256i vtarget = _mm256_set1_epi64x(target);
  const __m256i vprefix =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(prefix_before));
  const __m256i vmin =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(min_excess));
  const __m256i vmax =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(max_excess));

  const __m256i relative = _mm256_sub_epi64(vtarget, vprefix);
  const __m256i ge_min = _mm256_or_si256(_mm256_cmpgt_epi64(relative, vmin),
                                         _mm256_cmpeq_epi64(relative, vmin));
  const __m256i le_max = _mm256_or_si256(_mm256_cmpgt_epi64(vmax, relative),
                                         _mm256_cmpeq_epi64(vmax, relative));
  __m256i matched = _mm256_and_si256(ge_min, le_max);
  if (include_zero_boundary) {
    matched = _mm256_or_si256(matched, _mm256_cmpeq_epi64(vtarget, vprefix));
  }

  const uint32_t byte_mask =
      static_cast<uint32_t>(_mm256_movemask_epi8(matched));
  uint32_t result = 0;
  for (size_t lane = 0; lane < 4; ++lane) {
    const uint32_t lane_mask = 0xffu << (lane * 8);
    if ((byte_mask & lane_mask) == lane_mask) {
      result |= uint32_t{1} << lane;
    }
  }
  return result;
#else
  uint32_t result = 0;
  for (size_t lane = 0; lane < 4; ++lane) {
    const int64_t relative = target - prefix_before[lane];
    const bool found =
        (min_excess[lane] <= relative && relative <= max_excess[lane]) ||
        (include_zero_boundary && relative == 0);
    if (found) {
      result |= uint32_t{1} << lane;
    }
  }
  return result;
#endif
}

/**
 * @brief Return a mask with the lowest @p num bits set.
 * @details Values greater than or equal to 64 produce an all-ones mask, which
 * avoids undefined behavior from shifting by the word width.
 * @param num Number of low bits to set.
 * @return A 64-bit mask containing ones in positions `[0, num)`.
 */
static inline uint64_t first_bits_mask(size_t num) {
  return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
}

/**
 * @brief Number of 1 bits in positions 0 .. count - 1
 * @details Assumes count
 * < 512.
 * @details Surprisingly one cannot just do (1 << L) - 1 for L > 64 to
 * produce mask of ones of length L. The best we can do with a single
 * instruction is (1 << L) - 1 for L=8k using maskz_set1.
 *
 * To make mask for arbitrary L we use shldv instuction, it doesn't
 * really matter what epi is used the recipe is the following:
 * - produce mask with x * (count / x) ones and store in a
 * - produce mask with x * ((count / x) + 1) ones and store in b
 * - use shldv(a, b, count % x) to blend a and b producing required mask
 *
 * Essentially shldv_epiX(a, b, k) takes k high bits of b and (X - k) bits of a
 * and concats them (bits of b become lower bits of the result).
 *
 * The rest is standard, i.e. popcount_epi64 to perform popcount on
 * 64 bits and then reduce_add to sum the result.
 */
static inline uint64_t rank_512(const uint64_t* x, uint64_t count) {
#ifdef PIXIE_AVX512_SUPPORT

  __m512i a = _mm512_maskz_set1_epi64((1ull << ((count >> 6))) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i b = _mm512_maskz_set1_epi64((1ull << ((count >> 6) + 1)) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i mask = _mm512_shldv_epi64(a, b, _mm512_set1_epi64(count % 64));

  __m512i res = _mm512_loadu_epi64(x);
  res = _mm512_and_epi64(res, mask);
  __m512i cnt = _mm512_popcnt_epi64(res);
  return _mm512_reduce_add_epi64(cnt);

#else

  uint64_t last_uint = count < 512 ? count >> 6 : 8;

  uint64_t pop_val = 0;

  for (int i = 0; i < last_uint; i++) {
    pop_val += std::popcount(x[i]);
  }

  pop_val += count < 512
                 ? std::popcount(x[last_uint] & first_bits_mask(count & 63))
                 : 0;
  return pop_val;

#endif
}

/**
 * @brief Return position of @p rank 1 bit in @p x
 * @details Uses BMI/BMI2 when enabled and a portable fallback otherwise.
 */
static inline uint64_t select_64(uint64_t x, uint64_t rank) {
#ifdef PIXIE_BMI_SUPPORT
  return _tzcnt_u64(_pdep_u64(1ull << rank, x));
#else
  while (rank != 0) {
    x &= x - 1;
    --rank;
  }
  return std::countr_zero(x);
#endif
}

/**
 * @brief Return position of @p rank 1 bit in @p x
 * @details Selecting within 64-bit word can be done
 * using combination of _tzcnt_u64(_pdep_u64(1 << rank, x))
 * See Pandey P., Bender M. A., Johnson R. A fast x86 implementation of select
 * https://arxiv.org/abs/1706.00990
 *
 * To find a 64-bit word inside a 512-bit region we use
 * We first popcounts of all 8 64-bit words with _mm512_popcnt_epi64
 * and then perform a linear scan.
 *
 * Notably a SWAR algorithm for parallel binary search
 * http://www-graphics.stanford.edu/~seander/bithacks.html#SelectPosFromMSBRank
 * might be used as a backoff algorithm for selecting in a 64-bit word.
 * It can also be used as an alternative for linear search but i don't
 * see a proper SIMD algorithm to make it faster.
 */
static inline uint64_t select_512(const uint64_t* x, uint64_t rank) {
#ifdef PIXIE_AVX512_SUPPORT

  __m512i res = _mm512_loadu_epi64(x);
  __m512i counts = _mm512_popcnt_epi64(res);
  __m512i prefix = counts;

  const __m512i idx_shift1 = _mm512_set_epi64(6, 5, 4, 3, 2, 1, 0, 0);
  const __m512i idx_shift2 = _mm512_set_epi64(5, 4, 3, 2, 1, 0, 0, 0);
  const __m512i idx_shift4 = _mm512_set_epi64(3, 2, 1, 0, 0, 0, 0, 0);

  __m512i tmp = _mm512_maskz_permutexvar_epi64(0xFE, idx_shift1, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xFC, idx_shift2, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xF0, idx_shift4, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);

  __mmask8 mask = _mm512_cmpgt_epu64_mask(prefix, _mm512_set1_epi64(rank));
  uint32_t i = _tzcnt_u32(static_cast<uint32_t>(mask));
  uint64_t prev = 0;
  if (i != 0) {
    __m512i idx_prev = _mm512_set1_epi64(static_cast<int64_t>(i - 1));
    __m512i prev_vec = _mm512_permutexvar_epi64(idx_prev, prefix);
    prev = static_cast<uint64_t>(
        _mm_cvtsi128_si64(_mm512_castsi512_si128(prev_vec)));
  }
  return i * 64 + select_64(x[i], rank - prev);

#else

  size_t i = 0;
  int popcount = std::popcount(x[0]);
  while (i < 7 && popcount <= rank) {
    rank -= popcount;
    popcount = std::popcount(x[++i]);
  }
  return i * 64 + select_64(x[i], rank);

#endif
}

#ifdef PIXIE_AVX512_SUPPORT
static inline uint64_t select0_512_from_inverted_words(const uint64_t* x,
                                                       uint64_t rank0,
                                                       __m512i res) {
  __m512i counts = _mm512_popcnt_epi64(res);
  __m512i prefix = counts;

  const __m512i idx_shift1 = _mm512_set_epi64(6, 5, 4, 3, 2, 1, 0, 0);
  const __m512i idx_shift2 = _mm512_set_epi64(5, 4, 3, 2, 1, 0, 0, 0);
  const __m512i idx_shift4 = _mm512_set_epi64(3, 2, 1, 0, 0, 0, 0, 0);

  __m512i tmp = _mm512_maskz_permutexvar_epi64(0xFE, idx_shift1, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xFC, idx_shift2, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xF0, idx_shift4, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);

  __mmask8 mask = _mm512_cmpgt_epu64_mask(prefix, _mm512_set1_epi64(rank0));
  uint32_t i = _tzcnt_u32(static_cast<uint32_t>(mask));
  uint64_t prev = 0;
  if (i != 0) {
    __m512i idx_prev = _mm512_set1_epi64(static_cast<int64_t>(i - 1));
    __m512i prev_vec = _mm512_permutexvar_epi64(idx_prev, prefix);
    prev = static_cast<uint64_t>(
        _mm_cvtsi128_si64(_mm512_castsi512_si128(prev_vec)));
  }
  return i * 64 + select_64(~x[i], rank0 - prev);
}
#endif

/**
 * @brief Return position of @p rank0 0 bit in @p x.
 * @details select_512 with bit inversion.
 */
static inline uint64_t select0_512(const uint64_t* x, uint64_t rank0) {
#ifdef PIXIE_AVX512_SUPPORT

  __m512i res = _mm512_loadu_epi64(x);
  res = _mm512_ternarylogic_epi64(res, res, res, 0x55);
  return select0_512_from_inverted_words(x, rank0, res);

#else

  size_t i = 0;
  int popcount = std::popcount(~x[0]);
  while (i < 7 && popcount <= rank0) {
    rank0 -= popcount;
    popcount = std::popcount(~x[++i]);
  }
  return i * 64 + select_64(~x[i], rank0);

#endif
}

/**
 * @brief Compare 4 64-bit numbers of @p x with @p y and
 * return the length of the prefix where @p y is less then @p x
 */
static inline uint16_t lower_bound_4x64(const uint64_t* x, uint64_t y) {
#ifdef PIXIE_AVX512_SUPPORT

  auto y_4 = _mm256_set1_epi64x(y);
  auto reg_256 = _mm256_loadu_epi64(x);
  auto cmp = _mm256_cmpge_epu64_mask(reg_256, y_4);

  return _tzcnt_u16(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  auto y_4 = _mm256_set1_epi64x(y);
  __m256i reg_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));

  const __m256i offset = _mm256_set1_epi64x(0x8000000000000000ULL);
  __m256i x_offset = _mm256_xor_si256(reg_256, offset);
  __m256i y_offset = _mm256_xor_si256(y_4, offset);
  auto mask = _mm256_movemask_epi8(_mm256_cmpgt_epi64(
      x_offset, _mm256_sub_epi64(y_offset, _mm256_set1_epi64x(1))));

  return _tzcnt_u32(mask) >> 3;

#else

  for (uint16_t i = 0; i < 4; ++i) {
    if (x[i] >= y) {
      return i;
    }
  }
  return 4;

#endif
#endif
}

/**
 * @brief Compare 4 64-bit numbers of ( @p delta_array + @p delta_scalar - @p x

 * * ) with @p y and return the length of the prefix
 * where @p y is less then
 * ( @p delta_array + @p delta_scalar - @p x )
 * @param x Base input array.
 *
 * @param y Threshold value for comparison.
 * @param delta_array Per-lane delta
 * offsets.
 * @param delta_scalar Shared delta offset.
 */
static inline uint16_t lower_bound_delta_4x64(const uint64_t* x,
                                              uint64_t y,
                                              const uint64_t* delta_array,
                                              uint64_t delta_scalar) {
#ifdef PIXIE_AVX512_SUPPORT

  const __m256i dlt_256 = _mm256_loadu_epi64(delta_array);
  auto x_256 = _mm256_loadu_epi64(x);
  auto dlt_4 = _mm256_set1_epi64x(delta_scalar);
  auto y_4 = _mm256_set1_epi64x(y);

  auto tmp = _mm256_add_epi64(dlt_4, dlt_256);
  auto reg_256 = _mm256_sub_epi64(tmp, x_256);
  auto cmp = _mm256_cmpge_epu64_mask(reg_256, y_4);

  return _tzcnt_u16(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  const __m256i dlt_256 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(delta_array));
  auto x_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));
  auto dlt_4 = _mm256_set1_epi64x(delta_scalar);
  auto y_4 = _mm256_set1_epi64x(y);

  auto tmp = _mm256_add_epi64(dlt_4, dlt_256);
  auto reg_256 = _mm256_sub_epi64(tmp, x_256);

  const __m256i offset = _mm256_set1_epi64x(0x8000000000000000ULL);
  __m256i x_offset = _mm256_xor_si256(reg_256, offset);
  __m256i y_offset = _mm256_xor_si256(y_4, offset);
  auto mask = _mm256_movemask_epi8(_mm256_cmpgt_epi64(
      x_offset, _mm256_sub_epi64(y_offset, _mm256_set1_epi64x(1))));

  return _tzcnt_u32(mask) >> 3;

#else

  for (uint16_t i = 0; i < 4; ++i) {
    if (delta_array[i] + delta_scalar - x[i] >= y) {
      return i;
    }
  }
  return 4;

#endif
#endif
}

/**
 * @brief Compare 8 64-bit numbers of @p x with @p y and
 * return the length of the prefix where @p y is less then @p x
 */
static inline uint16_t lower_bound_8x64(const uint64_t* x, uint64_t y) {
#ifdef PIXIE_AVX512_SUPPORT

  auto y_8 = _mm512_set1_epi64(y);
  auto reg_512 = _mm512_loadu_epi64(x);
  auto cmp = _mm512_cmpge_epu64_mask(reg_512, y_8);

  return _tzcnt_u16(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  uint16_t len = lower_bound_4x64(x, y);

  if (len < 4) {
    return len;
  }

  return len + lower_bound_4x64(x + 4, y);

#else

  for (uint16_t i = 0; i < 8; ++i) {
    if (x[i] >= y) {
      return i;
    }
  }
  return 8;

#endif
#endif
}

/**
 * @brief Compare 8 64-bit numbers of ( @p delta_array + @p delta_scalar - @p x

 * * ) with @p y and return the length of the prefix
 * where @p y is less then
 * ( @p delta_array + @p delta_scalar - @p x )
 * @param x Base input array.
 *
 * @param y Threshold value for comparison.
 * @param delta_array Per-lane delta
 * offsets.
 * @param delta_scalar Shared delta offset.
 */
static inline uint16_t lower_bound_delta_8x64(const uint64_t* x,
                                              uint64_t y,
                                              const uint64_t* delta_array,
                                              uint64_t delta_scalar) {
#ifdef PIXIE_AVX512_SUPPORT

  const __m512i dlt_512 = _mm512_loadu_epi64(delta_array);
  auto x_512 = _mm512_loadu_epi64(x);
  auto dlt_8 = _mm512_set1_epi64(delta_scalar);
  auto y_8 = _mm512_set1_epi64(y);

  auto tmp = _mm512_add_epi64(dlt_8, dlt_512);
  auto reg_512 = _mm512_sub_epi64(tmp, x_512);
  auto cmp = _mm512_cmpge_epu64_mask(reg_512, y_8);

  return _tzcnt_u16(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  uint16_t len = lower_bound_delta_4x64(x, y, delta_array, delta_scalar);

  if (len < 4) {
    return len;
  }

  return len + lower_bound_delta_4x64(x + 4, y, delta_array + 4, delta_scalar);

#else

  for (uint16_t i = 0; i < 8; ++i) {
    if (delta_array[i] + delta_scalar - x[i] >= y) {
      return i;
    }
  }
  return 8;

#endif
#endif
}

/**
 * @brief Compare 32 16-bit numbers of @p x with @p y and
 * return the count of numbers where @p x is less then @p y
 */
static inline uint16_t lower_bound_32x16(const uint16_t* x, uint16_t y) {
#ifdef PIXIE_AVX512_SUPPORT

  auto y_32 = _mm512_set1_epi16(y);
  auto reg_512 = _mm512_loadu_epi16(x);
  auto cmp = _mm512_cmplt_epu16_mask(reg_512, y_32);
  return std::popcount(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  auto y_16 = _mm256_set1_epi16(y);
  __m256i reg_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));

  const __m256i offset = _mm256_set1_epi16(0x8000);
  __m256i x_offset = _mm256_xor_si256(reg_256, offset);
  __m256i y_offset = _mm256_xor_si256(y_16, offset);
  uint32_t mask = _mm256_movemask_epi8(_mm256_cmpgt_epi16(y_offset, x_offset));

  uint16_t count = std::popcount(mask) >> 1;

  reg_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + 16));

  x_offset = _mm256_xor_si256(reg_256, offset);
  mask = _mm256_movemask_epi8(_mm256_cmpgt_epi16(y_offset, x_offset));

  return count + (std::popcount(mask) >> 1);

#else

  uint16_t cnt = 0;
  for (uint16_t i = 0; i < 32; ++i) {
    if (x[i] < y) {
      cnt++;
    }
  }
  return cnt;

#endif
#endif
}

/**
 * @brief Compare 32 16-bit numbers of ( @p delta_array + @p delta_scalar - @p
 * x
 * ) with @p y and return the count of numbers where
 * ( @p delta_array +
 * @p delta_scalar - @p x ) is less then @p y
 * @param x Base input array.
 *
 * @param y Threshold value for comparison.
 * @param delta_array Per-lane delta
 * offsets.
 * @param delta_scalar Shared delta offset.
 */
static inline uint16_t lower_bound_delta_32x16(const uint16_t* x,
                                               uint16_t y,
                                               const uint16_t* delta_array,
                                               uint16_t delta_scalar) {
#ifdef PIXIE_AVX512_SUPPORT

  const __m512i dlt_512 = _mm512_loadu_epi64(delta_array);
  auto x_512 = _mm512_loadu_epi64(x);
  auto dlt_32 = _mm512_set1_epi16(delta_scalar);
  auto y_32 = _mm512_set1_epi16(y);

  auto tmp = _mm512_add_epi16(dlt_32, dlt_512);
  auto reg_512 = _mm512_sub_epi16(tmp, x_512);
  auto cmp = _mm512_cmplt_epu16_mask(reg_512, y_32);
  return std::popcount(cmp);

#else
#ifdef PIXIE_AVX2_SUPPORT

  auto dlt_256 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(delta_array));
  auto x_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));
  auto dlt_16 = _mm256_set1_epi16(delta_scalar);
  auto y_16 = _mm256_set1_epi16(y);

  auto tmp = _mm256_add_epi16(dlt_16, dlt_256);
  auto reg_256 = _mm256_sub_epi16(tmp, x_256);

  const __m256i offset = _mm256_set1_epi16(0x8000);
  __m256i x_offset = _mm256_xor_si256(reg_256, offset);
  __m256i y_offset = _mm256_xor_si256(y_16, offset);
  uint32_t mask = _mm256_movemask_epi8(_mm256_cmpgt_epi16(y_offset, x_offset));

  uint16_t count = std::popcount(mask) >> 1;

  dlt_256 =
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(delta_array + 16));
  x_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + 16));

  tmp = _mm256_add_epi16(dlt_16, dlt_256);
  reg_256 = _mm256_sub_epi16(tmp, x_256);

  x_offset = _mm256_xor_si256(reg_256, offset);
  mask = _mm256_movemask_epi8(_mm256_cmpgt_epi16(y_offset, x_offset));

  return count + (std::popcount(mask) >> 1);

#else

  uint16_t cnt = 0;
  for (uint16_t i = 0; i < 32; ++i) {
    if (delta_array[i] + delta_scalar - x[i] < y) {
      cnt++;
    }
  }
  return cnt;

#endif
#endif
}

/**
 * @brief Calculates 64 popcounts of 4-bits integers and stores as 64 4-bits
 * integers (packed into 32 bytes)
 *
 * This function counts the number of set bits in each 4-bit value
 * using an efficient shuffle-based approach with lookup tables.
 *
 * @param x Pointer to 32 bytes containing 64 4-bit integers
 * @param result Pointer to store the 64 resulting 4-bit popcount values (packed
 * in 32 bytes)
 */
static inline void popcount_64x4(const uint8_t* x, uint8_t* result) {
#ifdef PIXIE_AVX512_SUPPORT
  __m256i data = _mm256_loadu_si256((__m256i const*)x);

  // Masks for extracting the lower and upper nibbles
  const __m256i low_bits_mask = _mm256_set1_epi8(0x0F);

  // Count bits in the lower half
  __m256i low_bits = _mm256_and_si256(data, low_bits_mask);
  __m256i low_count = _mm256_shuffle_epi8(lookup_popcount_4, low_bits);

  // Count bits in the upper half
  __m256i high_bits = _mm256_srli_epi16(data, 4);
  high_bits = _mm256_and_si256(high_bits, low_bits_mask);
  __m256i high_count = _mm256_shuffle_epi8(lookup_popcount_4, high_bits);

  // Pack the results into a single output vector
  __m256i result_vec =
      _mm256_or_si256(low_count, _mm256_slli_epi16(high_count, 4));
  _mm256_storeu_epi8(result, result_vec);
#else
  // Fallback implementation for non-AVX2 platforms
  for (size_t i = 0; i < 32; i++) {
    // Count bits in the lower half
    uint8_t a = x[i] & 0x0F;
    uint8_t low_count = std::popcount(a);
    // Count bits in the upper half
    a = (x[i] >> 4) & 0x0F;
    uint8_t high_count = std::popcount(a);

    // Pack the counts into the output byte
    result[i] = low_count | (high_count << 4);
  }
#endif
}

/**
 * @brief Calculates 64 popcounts of 4-bits integers and stores as 64 4-bits
 * integers (packed into 32 bytes)
 *
 * This function counts the number of set bits in each 4-bit value
 * using an efficient shuffle-based approach with lookup tables.
 *
 * @param x Pointer to 32 bytes containing 64 4-bit integers
 * @param result Pointer to store the 64 resulting 4-bit popcount values
 * (packed in 32 bytes)
 */
static inline void popcount_32x8(const uint8_t* x, uint8_t* result) {
#ifdef PIXIE_AVX512_SUPPORT
  // Load 64 4-bit integers (256 bits total)
  __m256i data = _mm256_loadu_si256((__m256i const*)x);
  auto popcount_8 = _mm256_popcnt_epi8(data);
  _mm256_storeu_si256((__m256i*)result, popcount_8);
#else
#ifdef PIXIE_AVX2_SUPPORT
  // Load 64 4-bit integers (256 bits total)
  __m256i data = _mm256_loadu_si256((__m256i const*)x);

  // Masks for extracting the lower and upper nibbles
  const __m256i low_bits_mask = _mm256_set1_epi8(0x0F);

  // Count bits in lower half
  __m256i low_bits = _mm256_and_si256(data, low_bits_mask);
  __m256i low_count = _mm256_shuffle_epi8(lookup_popcount_4, low_bits);

  // Count bits upper half
  __m256i high_bits = _mm256_srli_epi16(data, 4);
  high_bits = _mm256_and_si256(high_bits, low_bits_mask);
  __m256i high_count = _mm256_shuffle_epi8(lookup_popcount_4, high_bits);

  __m256i result_vec = _mm256_add_epi8(low_count, high_count);
  _mm256_storeu_si256((__m256i*)result, result_vec);
#else
  // Fallback implementation for non-AVX2 platforms
  for (size_t i = 0; i < 32; i++) {
    result[i] = std::popcount(x[i]);
  }
#endif
#endif
}

#ifdef PIXIE_AVX2_SUPPORT
// clang-format off
// LUT for total excess change across a 4-bit nibble
static inline const __m256i excess_lut_delta = _mm256_setr_epi8(
    -4, -2, -2,  0,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
     0,  2,  2,  4,
    -4, -2, -2,  0,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
     0,  2,  2,  4);

// LUTs for target relative excess positions
static inline const __m256i excess_lut_pos0 = _mm256_setr_epi8(
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1,
    -1,  1, -1,  1);

static inline const __m256i excess_lut_pos1 = _mm256_setr_epi8(
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
    -2,  0,  0,  2);

static inline const __m256i excess_lut_pos2 = _mm256_setr_epi8(
    -3, -1, -1,  1,
    -1,  1,  1,  3,
    -3, -1, -1,  1,
    -1,  1,  1,  3,
    -3, -1, -1,  1,
    -1,  1,  1,  3,
    -3, -1, -1,  1,
    -1,  1,  1,  3);
static inline const __m256i excess_lut_min = _mm256_setr_epi8(
    -4, -2, -2,  0,
    -2,  0, -1,  1,
    -3, -1, -1,  1,
    -2,  0, -1,  1,
    -4, -2, -2,  0,
    -2,  0, -1,  1,
    -3, -1, -1,  1,
    -2,  0, -1,  1);
static inline constexpr int8_t excess_lut_min_offset[16] = {
    4, 4, 4, 4, 2, 2, 1, 1, 3, 3, 1, 1, 2, 2, 1, 1};
static inline const __m256i excess_lut_pack_multiplier =
    _mm256_set1_epi16(0x1001);
static inline const __m256i excess_lut_bit0 = _mm256_set1_epi8(1);
static inline const __m256i excess_lut_bit1 = _mm256_set1_epi8(2);
static inline const __m256i excess_lut_bit2 = _mm256_set1_epi8(4);
static inline const __m256i excess_lut_bit3 = _mm256_set1_epi8(8);
static inline const __m256i excess_lut_nibble_index = _mm256_setr_epi8(
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12, 13, 14, 15,
    16, 17, 18, 19,
    20, 21, 22, 23,
    24, 25, 26, 27,
    28, 29, 30, 31);
static inline const __m128i excess_lut_nibble_mask = _mm_set1_epi8(0x0F);
// clang-format on

static inline __m256i excess_nibbles_128_avx2(const uint64_t* s) noexcept {
  __m128i word_vec = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s));
  __m128i lo_nibbles = _mm_and_si128(word_vec, excess_lut_nibble_mask);
  __m128i hi_nibbles =
      _mm_and_si128(_mm_srli_epi16(word_vec, 4), excess_lut_nibble_mask);

  __m128i unpack_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
  __m128i unpack_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);

  return _mm256_inserti128_si256(_mm256_castsi128_si256(unpack_lo), unpack_hi,
                                 1);
}

static inline __m256i excess_bit_masks_16x_i16() noexcept {
  return _mm256_setr_epi16(0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020,
                           0x0040, 0x0080, 0x0100, 0x0200, 0x0400, 0x0800,
                           0x1000, 0x2000, 0x4000,
                           static_cast<int16_t>(0x8000));
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
  const int16_t carry = static_cast<int16_t>(_mm_extract_epi16(lo, 7));
  hi = _mm_add_epi16(hi, _mm_set1_epi16(carry));

  __m256i out = _mm256_castsi128_si256(lo);
  return _mm256_inserti128_si256(out, hi, 1);
}
#endif

/**
 * @brief Minimum prefix excess in a 128-bit bitstring range.
 * @details Prefix positions are offsets in `[0, 128]`; position 0 is the
 * empty prefix and position `k` is the excess after consuming the first `k`
 * bits. The query range `[left, right]` is inclusive. Ties return the first
 * offset attaining the minimum. Invalid ranges return offset 128 as a
 * sentinel.
 */
struct ExcessResult {
  int min_excess = 0;
  size_t offset = 128;
};

/**
 * @brief Pair of boundary minimum results for adjacent BP query blocks.
 *
 * @details `suffix` is local to the left/suffix block and `prefix` is local to
 * the right/prefix block.
 */
struct ExcessBoundaryPairResult {
  ExcessResult suffix;
  ExcessResult prefix;
};

constexpr int8_t excess_byte_delta_value(uint8_t x) {
  return static_cast<int8_t>(2 * std::popcount(x) - 8);
}

constexpr int8_t excess_byte_min_prefix_value(uint8_t x) {
  int cur = 0;
  int best = 0;
  for (int bit = 0; bit < 8; ++bit) {
    cur += ((x >> bit) & 1u) != 0 ? 1 : -1;
    if (bit == 0 || cur < best) {
      best = cur;
    }
  }
  return static_cast<int8_t>(best);
}

constexpr int8_t excess_byte_min_prefix_offset_value(uint8_t x) {
  int cur = 0;
  int best = 0;
  int best_offset = 1;
  for (int bit = 0; bit < 8; ++bit) {
    cur += ((x >> bit) & 1u) != 0 ? 1 : -1;
    if (bit == 0 || cur < best) {
      best = cur;
      best_offset = bit + 1;
    }
  }
  return static_cast<int8_t>(best_offset);
}

constexpr int8_t excess_nibble_min_prefix_offset_value(uint8_t x, int bits) {
  int cur = 0;
  int best = 0;
  int best_offset = 1;
  for (int bit = 0; bit < bits; ++bit) {
    cur += ((x >> bit) & 1u) != 0 ? 1 : -1;
    if (bit == 0 || cur < best) {
      best = cur;
      best_offset = bit + 1;
    }
  }
  return static_cast<int8_t>(best_offset);
}

template <typename Fn>
constexpr std::array<int8_t, 256> excess_make_byte_lut(Fn fn) {
  std::array<int8_t, 256> out{};
  for (size_t i = 0; i < out.size(); ++i) {
    out[i] = fn(static_cast<uint8_t>(i));
  }
  return out;
}

static inline constexpr std::array<int8_t, 256> excess_byte_delta_lut =
    excess_make_byte_lut([](uint8_t x) { return excess_byte_delta_value(x); });
static inline constexpr std::array<int8_t, 256> excess_byte_min_lut =
    excess_make_byte_lut(
        [](uint8_t x) { return excess_byte_min_prefix_value(x); });
static inline constexpr std::array<int8_t, 256> excess_byte_min_offset_lut =
    excess_make_byte_lut(
        [](uint8_t x) { return excess_byte_min_prefix_offset_value(x); });
static inline constexpr std::array<std::array<int8_t, 16>, 4>
    excess_partial_nibble_min_offset_lut = [] {
      std::array<std::array<int8_t, 16>, 4> out{};
      for (size_t width = 1; width < out.size(); ++width) {
        for (size_t nibble = 0; nibble < out[width].size(); ++nibble) {
          out[width][nibble] = excess_nibble_min_prefix_offset_value(
              static_cast<uint8_t>(nibble), static_cast<int>(width));
        }
      }
      return out;
    }();

/**
 * @brief Compute record-low mask for a single byte relative to a threshold.
 *
 * @details For each bit position in the byte (0..7), computes the local
 * excess starting from 0. Sets the corresponding bit in the result mask
 * iff the local excess is strictly less than @p threshold.
 */
constexpr uint8_t excess_byte_record_lows_mask(uint8_t byte, int threshold) {
  int cur = 0;
  uint8_t mask = 0;
  for (int bit = 0; bit < 8; ++bit) {
    cur += ((byte >> bit) & 1u) ? 1 : -1;
    if (cur < threshold) {
      mask |= static_cast<uint8_t>(1u << bit);
    }
  }
  return mask;
}

/**
 * @brief LUT for record-low masks within a single byte.
 *
 * @details For each byte value (256) and each possible gap g in [0, 7],
 * stores a bitmask of positions whose local excess is strictly less than
 * -g.  Gap g = start_excess - best_excess; when g >= 8 no new record low
 * is possible because the byte's minimum local excess is -8.
 */
static inline constexpr std::array<std::array<uint8_t, 8>, 256>
    excess_byte_record_lows_lut = [] {
      std::array<std::array<uint8_t, 8>, 256> out{};
      for (size_t byte = 0; byte < 256; ++byte) {
        for (int g = 0; g < 8; ++g) {
          out[byte][g] =
              excess_byte_record_lows_mask(static_cast<uint8_t>(byte), -g);
        }
      }
      return out;
    }();

/**
 * @brief Find every prefix whose excess equals target_x in a 128-bit bitstring.
 *
 * Excess(i) = 2*popcount(bits[0..i-1]) - i   for i in [0..128].
 * Bit (w*64 + b) of out[w] is set iff excess(w*64 + b + 1) == target_x.
 * I.e. out bit index b corresponds to prefix length (b+1).
 *
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param target_x Target excess value in [-128, 128]; outside this range out is
 * zeroed.
 * @param out 2 uint64_t words receiving the result bitmask.
 * @return Total excess change across the 128-bit bitstring.
 */
static inline int excess_positions_128(const uint64_t* s,
                                       int target_x,
                                       uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  const int block_delta = 2 * (std::popcount(s[0]) + std::popcount(s[1])) - 128;

  if (target_x < -128 || target_x > 128) {
    return block_delta;
  }

#ifdef PIXIE_AVX2_SUPPORT
  const __m256i vdelta = excess_lut_delta;
  const __m256i vpos0 = excess_lut_pos0;
  const __m256i vpos1 = excess_lut_pos1;
  const __m256i vpos2 = excess_lut_pos2;
  const __m256i vmult = excess_lut_pack_multiplier;
  const __m256i vbit0 = excess_lut_bit0;
  const __m256i vbit1 = excess_lut_bit1;
  const __m256i vbit2 = excess_lut_bit2;
  const __m256i vbit3 = excess_lut_bit3;

  const int d = 2 * target_x - block_delta;
  if (d < -128 || d > 128) {
    return block_delta;
  }

  __m256i nibbles = excess_nibbles_128_avx2(s);

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

  __m256i vtgt = _mm256_set1_epi8((int8_t)target_x);
  __m256i t = _mm256_sub_epi8(vtgt, excl_ps);

  __m256i cmp0 = _mm256_cmpeq_epi8(_mm256_shuffle_epi8(vpos0, nibbles), t);
  __m256i cmp1 = _mm256_cmpeq_epi8(_mm256_shuffle_epi8(vpos1, nibbles), t);
  __m256i cmp2 = _mm256_cmpeq_epi8(_mm256_shuffle_epi8(vpos2, nibbles), t);
  __m256i cmp3 = _mm256_cmpeq_epi8(ps, vtgt);

  __m256i bit0 = _mm256_and_si256(cmp0, vbit0);
  __m256i bit1 = _mm256_and_si256(cmp1, vbit1);
  __m256i bit2 = _mm256_and_si256(cmp2, vbit2);
  __m256i bit3 = _mm256_and_si256(cmp3, vbit3);

  __m256i total_match =
      _mm256_or_si256(_mm256_or_si256(bit0, bit1), _mm256_or_si256(bit2, bit3));

  __m256i res = _mm256_maddubs_epi16(total_match, vmult);
  __m128i res_lo = _mm256_castsi256_si128(res);
  __m128i res_hi = _mm256_extracti128_si256(res, 1);
  __m128i packed = _mm_packus_epi16(res_lo, res_hi);

  _mm_storeu_si128((__m128i*)out, packed);
#else
  int cur = 0;
  for (size_t i = 0; i < 128; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
#endif
  return block_delta;
}

/**
 * @brief Prefix excess in a 128-bit bitstring.
 *
 * Excess(i) = 2*popcount(bits[0..i-1]) - i for i in [0, 128].
 *
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param end_offset Exclusive prefix boundary, clamped to [0, 128].
 * @return Prefix excess on [0, end_offset).
 */
static inline int prefix_excess_128(const uint64_t* s,
                                    size_t end_offset) noexcept {
  end_offset = end_offset > 128 ? 128 : end_offset;
  if (end_offset == 0) {
    return 0;
  }
  if (end_offset <= 64) {
    const int ones = static_cast<int>(std::popcount(
        s[0] & first_bits_mask(static_cast<uint32_t>(end_offset))));
    return 2 * ones - static_cast<int>(end_offset);
  }
  const int ones = static_cast<int>(
      std::popcount(s[0]) +
      std::popcount(s[1] &
                    first_bits_mask(static_cast<uint32_t>(end_offset - 64))));
  return 2 * ones - static_cast<int>(end_offset);
}

/**
 * @brief Prefix excess in a 64-bit bitstring.
 *
 * @details `excess(i) = 2 * popcount(bits[0..i)) - i` for `i` in `[0, 64]`.
 *
 * @param s One little-endian 64-bit word (bit 0 of `s[0]` is first).
 * @param end_offset Exclusive prefix boundary, clamped to `[0, 64]`.
 * @return Prefix excess on `[0, end_offset)`.
 */
static inline int prefix_excess_64(const uint64_t* s,
                                   size_t end_offset) noexcept {
  end_offset = end_offset > 64 ? 64 : end_offset;
  if (end_offset == 0) {
    return 0;
  }
  const int ones = static_cast<int>(
      std::popcount(s[0] & first_bits_mask(static_cast<uint32_t>(end_offset))));
  return 2 * ones - static_cast<int>(end_offset);
}

static inline ExcessResult excess_min_128_byte_lut_short(
    const uint64_t* s,
    size_t left,
    size_t right) noexcept {
  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 7u) != 0; ++bit) {
    current += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }

  for (; bit + 8 <= right; bit += 8) {
    const uint8_t byte =
        static_cast<uint8_t>((s[bit >> 6] >> (bit & 63)) & 0xFFu);
    const int candidate = current + excess_byte_min_lut[byte];
    if (candidate < best) {
      best = candidate;
      best_offset = bit + static_cast<size_t>(excess_byte_min_offset_lut[byte]);
    }
    current += excess_byte_delta_lut[byte];
  }

  for (; bit < right; ++bit) {
    current += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }

  return {best, best_offset};
}

/**
 * @brief Return the minimum prefix excess and first attaining offset.
 *
 * @details Prefix positions are offsets in `[0, 64]`; position 0 is the empty
 * prefix and position `k` is the excess after consuming the first `k` bits.
 * The query range `[left, right]` is inclusive. Ties return the first offset
 * attaining the minimum. Invalid ranges return the default `ExcessResult`
 * sentinel.
 *
 * @param s One little-endian 64-bit word (bit 0 of `s[0]` is first).
 * @param left First prefix position to consider, inclusive.
 * @param right Last prefix position to consider, inclusive.
 * @return Minimum excess and first local offset attaining it.
 */
static inline ExcessResult excess_min_64(const uint64_t* s,
                                         size_t left,
                                         size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 64);
  right = std::min<size_t>(right, 64);

  int best = prefix_excess_64(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

#ifdef PIXIE_SSE41_SUPPORT
  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 3u) != 0; ++bit) {
    current += ((s[0] >> bit) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }

  const size_t first_full_nibble = bit >> 2;
  const size_t last_full_nibble = right >> 2;
  const size_t right_partial_width = bit < right ? (right & 3u) : 0;
  const size_t end_nibble =
      last_full_nibble + (right_partial_width == 0 ? 0 : 1);
  if (first_full_nibble < end_nibble) {
    const __m128i nibbles = excess_nibbles_64_sse(s);

    __m128i ps = _mm_shuffle_epi8(excess_lut_delta_sse, nibbles);
    ps = _mm_add_epi8(ps, _mm_slli_si128(ps, 1));
    ps = _mm_add_epi8(ps, _mm_slli_si128(ps, 2));
    ps = _mm_add_epi8(ps, _mm_slli_si128(ps, 4));
    ps = _mm_add_epi8(ps, _mm_slli_si128(ps, 8));

    const __m128i excl_ps = _mm_slli_si128(ps, 1);
    __m128i local_min = _mm_shuffle_epi8(excess_lut_min_sse, nibbles);
    if (right_partial_width != 0) {
      __m128i partial_min = _mm_shuffle_epi8(excess_lut_pos0_sse, nibbles);
      if (right_partial_width >= 2) {
        partial_min = _mm_min_epi8(
            partial_min, _mm_shuffle_epi8(excess_lut_pos1_sse, nibbles));
      }
      if (right_partial_width >= 3) {
        partial_min = _mm_min_epi8(
            partial_min, _mm_shuffle_epi8(excess_lut_pos2_sse, nibbles));
      }
      local_min = _mm_blendv_epi8(
          local_min, partial_min,
          _mm_cmpeq_epi8(excess_lut_nibble_index_sse,
                         _mm_set1_epi8(static_cast<int8_t>(last_full_nibble))));
    }
    const __m128i partial_candidates = _mm_add_epi8(excl_ps, local_min);

    const __m128i idx = excess_lut_nibble_index_sse;
    const int first_minus_one_value = static_cast<int>(first_full_nibble) - 1;
    const __m128i first_minus_one =
        _mm_set1_epi8(static_cast<int8_t>(first_minus_one_value));
    const __m128i last = _mm_set1_epi8(static_cast<int8_t>(end_nibble));
    const __m128i active = _mm_and_si128(_mm_cmpgt_epi8(idx, first_minus_one),
                                         _mm_cmpgt_epi8(last, idx));
    const __m128i masked_candidates =
        _mm_blendv_epi8(_mm_set1_epi8(127), partial_candidates, active);

    __m128i min128 = masked_candidates;
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 8));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 4));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 2));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 1));

    const int candidate_min =
        static_cast<int>(static_cast<int8_t>(_mm_extract_epi8(min128, 0)));
    if (candidate_min < best) {
      const __m128i equal_min = _mm_cmpeq_epi8(
          masked_candidates, _mm_set1_epi8(static_cast<int8_t>(candidate_min)));
      const uint32_t equal_mask =
          static_cast<uint32_t>(_mm_movemask_epi8(equal_min));
      const uint32_t nibble_index = std::countr_zero(equal_mask);
      const uint8_t nibble =
          static_cast<uint8_t>((s[0] >> (nibble_index * 4u)) & 0xFu);
      best = candidate_min;
      if (right_partial_width != 0 && nibble_index == last_full_nibble) {
        int local = 0;
        int local_best = 0;
        size_t local_offset = 1;
        for (size_t i = 0; i < right_partial_width; ++i) {
          local += ((nibble >> i) & 1u) != 0 ? 1 : -1;
          if (i == 0 || local < local_best) {
            local_best = local;
            local_offset = i + 1;
          }
        }
        best_offset = static_cast<size_t>(nibble_index) * 4u + local_offset;
      } else {
        best_offset = static_cast<size_t>(nibble_index) * 4u +
                      static_cast<size_t>(excess_nibble_min_offset[nibble]);
      }
    }

    bit = end_nibble * 4;
  }

  for (; bit < right; ++bit) {
    current += ((s[0] >> bit) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }
#else
  int current = best;
  for (size_t bit = left; bit < right; ++bit) {
    current += ((s[0] >> bit) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }
#endif

  return {best, best_offset};
}

/**
 * @brief Return the minimum prefix excess and first attaining offset.
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param left First prefix position to consider, inclusive.
 * @param right Last prefix position to consider, inclusive.
 */
static inline ExcessResult excess_min_128(const uint64_t* s,
                                          size_t left,
                                          size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  if (right - left <= 32 && (left & 7u) == 0 && (right & 7u) == 0)
      [[unlikely]] {
    return excess_min_128_byte_lut_short(s, left, right);
  }

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

#ifdef PIXIE_AVX2_SUPPORT
  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 3u) != 0; ++bit) {
    current += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (current < best) {
      best = current;
      best_offset = offset;
    }
  }

  const size_t first_nibble = bit >> 2;
  const size_t last_full_nibble = right >> 2;
  const size_t right_partial_width = bit < right ? (right & 3u) : 0;
  const size_t end_nibble =
      last_full_nibble + (right_partial_width == 0 ? 0 : 1);
  if (first_nibble < end_nibble) {
    const __m128i bytes = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s));
    const __m128i lo_nibbles = _mm_and_si128(bytes, excess_lut_nibble_mask_sse);
    const __m128i hi_nibbles =
        _mm_and_si128(_mm_srli_epi16(bytes, 4), excess_lut_nibble_mask_sse);
    const __m128i lo_delta = _mm_shuffle_epi8(excess_lut_delta_sse, lo_nibbles);
    const __m128i hi_delta = _mm_shuffle_epi8(excess_lut_delta_sse, hi_nibbles);
    const __m128i byte_delta = _mm_add_epi8(lo_delta, hi_delta);
    const __m128i byte_prefix = excess_prefix_sum_16x_i8(byte_delta);
    const __m128i byte_prefix_before = _mm_slli_si128(byte_prefix, 1);

    __m128i lo_local_min = _mm_shuffle_epi8(excess_lut_min_sse, lo_nibbles);
    __m128i hi_local_min = _mm_shuffle_epi8(excess_lut_min_sse, hi_nibbles);

    const __m128i byte_index = excess_lut_nibble_index_sse;
    if (right_partial_width != 0) {
      const bool partial_is_high = (last_full_nibble & 1u) != 0;
      const size_t partial_byte = last_full_nibble >> 1;
      const __m128i partial_source = partial_is_high ? hi_nibbles : lo_nibbles;
      __m128i partial_min =
          _mm_shuffle_epi8(excess_lut_pos0_sse, partial_source);
      if (right_partial_width >= 2) {
        partial_min = _mm_min_epi8(
            partial_min, _mm_shuffle_epi8(excess_lut_pos1_sse, partial_source));
      }
      if (right_partial_width >= 3) {
        partial_min = _mm_min_epi8(
            partial_min, _mm_shuffle_epi8(excess_lut_pos2_sse, partial_source));
      }
      const __m128i partial_lane = _mm_cmpeq_epi8(
          byte_index, _mm_set1_epi8(static_cast<int8_t>(partial_byte)));
      if (partial_is_high) {
        hi_local_min = _mm_blendv_epi8(hi_local_min, partial_min, partial_lane);
      } else {
        lo_local_min = _mm_blendv_epi8(lo_local_min, partial_min, partial_lane);
      }
    }

    const __m128i lo_candidates =
        _mm_add_epi8(byte_prefix_before, lo_local_min);
    const __m128i hi_candidates =
        _mm_add_epi8(_mm_add_epi8(byte_prefix_before, lo_delta), hi_local_min);

    __m128i masked_lo = lo_candidates;
    __m128i masked_hi = hi_candidates;
    if (first_nibble != 0 || end_nibble != 32) {
      const __m128i first_minus_one = _mm_set1_epi8(
          static_cast<int8_t>(static_cast<int>(first_nibble) - 1));
      const __m128i last = _mm_set1_epi8(static_cast<int8_t>(end_nibble));
      const __m128i lo_active = _mm_and_si128(
          _mm_cmpgt_epi8(excess_lut_low_nibble_index_sse, first_minus_one),
          _mm_cmpgt_epi8(last, excess_lut_low_nibble_index_sse));
      const __m128i hi_active = _mm_and_si128(
          _mm_cmpgt_epi8(excess_lut_high_nibble_index_sse, first_minus_one),
          _mm_cmpgt_epi8(last, excess_lut_high_nibble_index_sse));
      masked_lo = _mm_blendv_epi8(_mm_set1_epi8(127), lo_candidates, lo_active);
      masked_hi = _mm_blendv_epi8(_mm_set1_epi8(127), hi_candidates, hi_active);
    }

    const int candidate_min =
        excess_horizontal_min_i8(_mm_min_epi8(masked_lo, masked_hi));
    if (candidate_min < best) {
      const __m128i min_vec = _mm_set1_epi8(static_cast<int8_t>(candidate_min));
      const uint32_t lo_equal_mask = static_cast<uint32_t>(
          _mm_movemask_epi8(_mm_cmpeq_epi8(masked_lo, min_vec)));
      const uint32_t hi_equal_mask = static_cast<uint32_t>(
          _mm_movemask_epi8(_mm_cmpeq_epi8(masked_hi, min_vec)));
      const uint32_t lo_nibble_index =
          lo_equal_mask == 0
              ? 32u
              : static_cast<uint32_t>(std::countr_zero(lo_equal_mask)) * 2u;
      const uint32_t hi_nibble_index =
          hi_equal_mask == 0
              ? 32u
              : static_cast<uint32_t>(std::countr_zero(hi_equal_mask)) * 2u +
                    1u;
      const uint32_t nibble_index = std::min(lo_nibble_index, hi_nibble_index);
      const uint32_t byte_offset = nibble_index >> 1u;
      const uint64_t byte_word = s[byte_offset >> 3u];
      const uint8_t byte = static_cast<uint8_t>(
          (byte_word >> ((byte_offset & 7u) * 8u)) & 0xFFu);
      const uint8_t nibble = (nibble_index & 1u) == 0
                                 ? static_cast<uint8_t>(byte & 0xFu)
                                 : static_cast<uint8_t>((byte >> 4u) & 0xFu);
      const size_t local_offset =
          right_partial_width != 0 && nibble_index == last_full_nibble
              ? static_cast<size_t>(
                    excess_partial_nibble_min_offset_lut[right_partial_width]
                                                        [nibble])
              : static_cast<size_t>(excess_lut_min_offset[nibble]);
      best = candidate_min;
      best_offset = static_cast<size_t>(nibble_index) * 4u + local_offset;
    }
  }
#else
  int current = 0;
  for (size_t bit = 0; bit < right; ++bit) {
    current += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (offset >= left && current < best) {
      best = current;
      best_offset = offset;
    }
  }
#endif

  return {best, best_offset};
}

/**
 * @brief Compute disjoint suffix/prefix boundary minima for two 64-bit blocks.
 *
 * @details Computes `excess_min_64(suffix_s, suffix_left, 63)` and
 * `excess_min_64(prefix_s, 0, prefix_right)`. The individual calls use the
 * SSE nibble-LUT path when available.
 *
 * @param suffix_s Left boundary block.
 * @param suffix_left First local prefix offset in the suffix range.
 * @param prefix_s Right boundary block.
 * @param prefix_right Last local prefix offset in the prefix range.
 * @return Pair of local minimum results.
 */
static inline ExcessBoundaryPairResult excess_min_64_disjoint_suffix_prefix(
    const uint64_t* suffix_s,
    size_t suffix_left,
    const uint64_t* prefix_s,
    size_t prefix_right) noexcept {
  suffix_left = std::min<size_t>(suffix_left, 64);
  prefix_right = std::min<size_t>(prefix_right, 64);
  return {excess_min_64(suffix_s, suffix_left, 63),
          excess_min_64(prefix_s, 0, prefix_right)};
}

/**
 * @brief Compute disjoint suffix/prefix boundary minima for two 128-bit blocks.
 *
 * @details Computes `excess_min_128(suffix_s, suffix_left, 127)` and
 * `excess_min_128(prefix_s, 0, prefix_right)`. When `suffix_left >
 * prefix_right`, the AVX2 path blends active whole nibbles from both blocks
 * into one synthetic lane vector, fills the gap with neutral zero-delta
 * nibbles, and shares the prefix-sum/reduction work. Boundary fragments that
 * do not cover a whole nibble are handled by scalar or partial-nibble logic.
 * Non-disjoint or out-of-shape inputs fall back to the two independent
 * production calls.
 *
 * @param suffix_s Left boundary block.
 * @param suffix_left First local prefix offset in the suffix range.
 * @param prefix_s Right boundary block.
 * @param prefix_right Last local prefix offset in the prefix range.
 * @return Pair of local minimum results.
 */
static inline ExcessBoundaryPairResult excess_min_128_disjoint_suffix_prefix(
    const uint64_t* suffix_s,
    size_t suffix_left,
    const uint64_t* prefix_s,
    size_t prefix_right) noexcept {
  suffix_left = std::min<size_t>(suffix_left, 128);
  prefix_right = std::min<size_t>(prefix_right, 128);
  if (suffix_left <= prefix_right || suffix_left > 127 || prefix_right > 127) {
    return {excess_min_128(suffix_s, suffix_left, 127),
            excess_min_128(prefix_s, 0, prefix_right)};
  }

#ifdef PIXIE_AVX2_SUPPORT
  ExcessResult prefix{0, 0};

  int suffix_best = prefix_excess_128(suffix_s, suffix_left);
  ExcessResult suffix{suffix_best, suffix_left};
  size_t suffix_bit = suffix_left;
  int suffix_current = suffix_best;
  for (; suffix_bit < 127 && (suffix_bit & 3u) != 0; ++suffix_bit) {
    suffix_current +=
        ((suffix_s[suffix_bit >> 6] >> (suffix_bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = suffix_bit + 1;
    if (suffix_current < suffix.min_excess) {
      suffix = {suffix_current, offset};
    }
  }

  const size_t prefix_last_nibble = prefix_right >> 2;
  const size_t prefix_partial_width = prefix_right & 3u;
  const size_t prefix_end_nibble =
      prefix_last_nibble + (prefix_partial_width == 0 ? 0 : 1);
  const size_t suffix_first_nibble = suffix_bit < 127 ? suffix_bit >> 2 : 32;
  const int prefix_artificial_delta =
      prefix_excess_128(prefix_s, prefix_end_nibble * 4);

  const __m256i idx = excess_lut_nibble_index;
  const __m256i prefix_active = _mm256_cmpgt_epi8(
      _mm256_set1_epi8(static_cast<int8_t>(prefix_end_nibble)), idx);
  const __m256i suffix_active = _mm256_cmpgt_epi8(
      idx, _mm256_set1_epi8(static_cast<int8_t>(suffix_first_nibble) - 1));

  const __m256i prefix_nibbles = excess_nibbles_128_avx2(prefix_s);
  const __m256i suffix_nibbles = excess_nibbles_128_avx2(suffix_s);
  __m256i nibbles =
      _mm256_blendv_epi8(_mm256_set1_epi8(3), prefix_nibbles, prefix_active);
  nibbles = _mm256_blendv_epi8(nibbles, suffix_nibbles, suffix_active);

  __m256i ps = _mm256_shuffle_epi8(excess_lut_delta, nibbles);
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 1));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 2));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 4));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 8));

  __m128i ps_lo = _mm256_castsi256_si128(ps);
  __m128i ps_hi = _mm256_extracti128_si256(ps, 1);
  __m128i carry =
      _mm_set1_epi8(static_cast<int8_t>(_mm_extract_epi8(ps_lo, 15)));
  ps_hi = _mm_add_epi8(ps_hi, carry);
  ps = _mm256_inserti128_si256(_mm256_castsi128_si256(ps_lo), ps_hi, 1);

  __m256i b = _mm256_permute2x128_si256(ps, ps, 0x08);
  const __m256i excl_ps = _mm256_alignr_epi8(ps, b, 15);

  __m256i local_min = _mm256_shuffle_epi8(excess_lut_min, nibbles);
  if (prefix_partial_width != 0) {
    __m256i partial_min = _mm256_shuffle_epi8(excess_lut_pos0, nibbles);
    if (prefix_partial_width >= 2) {
      partial_min = _mm256_min_epi8(
          partial_min, _mm256_shuffle_epi8(excess_lut_pos1, nibbles));
    }
    if (prefix_partial_width >= 3) {
      partial_min = _mm256_min_epi8(
          partial_min, _mm256_shuffle_epi8(excess_lut_pos2, nibbles));
    }
    local_min = _mm256_blendv_epi8(
        local_min, partial_min,
        _mm256_cmpeq_epi8(
            idx, _mm256_set1_epi8(static_cast<int8_t>(prefix_last_nibble))));
  }
  __m256i suffix_partial_min = _mm256_min_epi8(
      _mm256_shuffle_epi8(excess_lut_pos0, nibbles),
      _mm256_min_epi8(_mm256_shuffle_epi8(excess_lut_pos1, nibbles),
                      _mm256_shuffle_epi8(excess_lut_pos2, nibbles)));
  const __m256i suffix_tail_active = _mm256_and_si256(
      suffix_active, _mm256_cmpeq_epi8(idx, _mm256_set1_epi8(31)));
  local_min =
      _mm256_blendv_epi8(local_min, suffix_partial_min, suffix_tail_active);

  const __m256i base_candidates = _mm256_add_epi8(excl_ps, local_min);
  const __m256i sentinel = _mm256_set1_epi8(127);

  const __m256i prefix_candidates =
      _mm256_blendv_epi8(sentinel, base_candidates, prefix_active);
  const __m256i suffix_candidates = _mm256_blendv_epi8(
      sentinel,
      _mm256_add_epi8(base_candidates,
                      _mm256_set1_epi8(static_cast<int8_t>(
                          suffix_current - prefix_artificial_delta))),
      suffix_active);

  auto reduce_min = [](__m256i values) {
    __m128i min128 = _mm_min_epi8(_mm256_castsi256_si128(values),
                                  _mm256_extracti128_si256(values, 1));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 8));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 4));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 2));
    min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 1));
    return static_cast<int>(static_cast<int8_t>(_mm_extract_epi8(min128, 0)));
  };

  auto local_offset = [](uint8_t nibble, size_t width_value) {
    if (width_value == 0 || width_value == 4) {
      return static_cast<size_t>(excess_lut_min_offset[nibble]);
    }
    int current = 0;
    int best = 0;
    size_t best_offset = 1;
    for (size_t i = 0; i < width_value; ++i) {
      current += ((nibble >> i) & 1u) != 0 ? 1 : -1;
      if (i == 0 || current < best) {
        best = current;
        best_offset = i + 1;
      }
    }
    return best_offset;
  };

  const int prefix_min = reduce_min(prefix_candidates);
  if (prefix_min < prefix.min_excess) {
    const uint32_t mask = static_cast<uint32_t>(_mm256_movemask_epi8(
        _mm256_cmpeq_epi8(prefix_candidates,
                          _mm256_set1_epi8(static_cast<int8_t>(prefix_min)))));
    const uint32_t prefix_lane = std::countr_zero(mask);
    const uint64_t word = prefix_s[prefix_lane >> 4];
    const uint8_t nibble =
        static_cast<uint8_t>((word >> ((prefix_lane & 15u) * 4u)) & 0xFu);
    prefix.min_excess = prefix_min;
    const size_t width =
        prefix_partial_width != 0 && prefix_lane == prefix_last_nibble
            ? prefix_partial_width
            : 4;
    prefix.offset =
        static_cast<size_t>(prefix_lane) * 4u + local_offset(nibble, width);
  }

  const int suffix_min = reduce_min(suffix_candidates);
  if (suffix_min < suffix.min_excess) {
    const uint32_t mask = static_cast<uint32_t>(_mm256_movemask_epi8(
        _mm256_cmpeq_epi8(suffix_candidates,
                          _mm256_set1_epi8(static_cast<int8_t>(suffix_min)))));
    const uint32_t suffix_lane = std::countr_zero(mask);
    const uint64_t word = suffix_s[suffix_lane >> 4];
    const uint8_t nibble =
        static_cast<uint8_t>((word >> ((suffix_lane & 15u) * 4u)) & 0xFu);
    suffix.min_excess = suffix_min;
    const size_t width = suffix_lane == 31 ? 3 : 4;
    suffix.offset =
        static_cast<size_t>(suffix_lane) * 4u + local_offset(nibble, width);
  }

  return {suffix, prefix};
#else
  return {excess_min_128(suffix_s, suffix_left, 127),
          excess_min_128(prefix_s, 0, prefix_right)};
#endif
}

/**
 * @brief Find the first prefix reaching target_x in a 128-bit bitstring.
 *
 * Searches the prefix excess values represented by excess_positions_128 and
 * ignores matches before start_offset. The returned offset is the bit position
 * whose inclusive prefix reaches target_x.
 *
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param target_x Target excess value relative to the beginning of this
 * 128-bit bitstring.
 * @param start_offset First bit position eligible for a match, in [0, 128].
 * @param block_excess Optional output for the total excess change across the
 * 128-bit bitstring.
 * @return Matching bit offset in [0, 127], or 128 if no match exists.
 */
static inline size_t forward_search_128(const uint64_t* s,
                                        int target_x,
                                        size_t start_offset,
                                        int* block_excess = nullptr) noexcept {
  uint64_t out[2];
  const int delta = excess_positions_128(s, target_x, out);
  if (block_excess != nullptr) {
    *block_excess = delta;
  }
  if (start_offset >= 128) {
    return 128;
  }

  const size_t first_word = start_offset >> 6;
  const size_t first_bit = start_offset & 63;
  for (size_t word = first_word; word < 2; ++word) {
    uint64_t mask = out[word];
    if (word == first_word && first_bit != 0) {
      mask &= ~first_bits_mask(first_bit);
    }
    if (mask != 0) {
      return word * 64 + std::countr_zero(mask);
    }
  }
  return 128;
}

/**
 * @brief Find the last prefix before end_offset reaching target_x in a 128-bit
 * bitstring.
 *
 * Searches prefix boundary positions strictly before end_offset, matching the
 * RmM backward-search convention. A return value of 0 is a valid match for the
 * chunk-start boundary when target_x is zero.
 *
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param target_x Target excess value relative to the beginning of this
 * 128-bit bitstring.
 * @param end_offset Exclusive right boundary for the search, in [0, 128].
 * @param block_excess Optional output for the total excess change across the
 * 128-bit bitstring.
 * @return Matching prefix boundary offset in [0, 127], or 128 if no match
 * exists.
 */
static inline size_t backward_search_128(const uint64_t* s,
                                         int target_x,
                                         size_t end_offset,
                                         int* block_excess = nullptr) noexcept {
  uint64_t out[2];
  const int delta = excess_positions_128(s, target_x, out);
  if (block_excess != nullptr) {
    *block_excess = delta;
  }
  if (end_offset == 0) {
    return 128;
  }

  const size_t max_prefix_length = end_offset - 1;
  if (max_prefix_length > 0) {
    const size_t last_bit_index = max_prefix_length - 1;
    size_t word = last_bit_index >> 6;
    const size_t bit_in_word = last_bit_index & 63;
    uint64_t mask = out[word] & first_bits_mask(bit_in_word + 1);
    while (true) {
      if (mask != 0) {
        return word * 64 + (63 - std::countl_zero(mask)) + 1;
      }
      if (word == 0) {
        break;
      }
      --word;
      mask = out[word];
    }
  }
  return target_x == 0 ? 0 : 128;
}

/**
 * @brief Find every prefix whose excess equals target_x in a 512-bit bitstring.
 *
 * Excess(i) = 2*popcount(bits[0..i-1]) - i   for i in [0..512].
 * Bit (w*64 + b) of out[w] is set iff excess(w*64 + b + 1) == target_x.
 * I.e. out bit index b corresponds to prefix length (b+1).
 *
 * @param s 8 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param target_x Target excess value in [-512, 512]; outside this range out is
 * zeroed.
 * @param out 8 uint64_t words receiving the result bitmask.
 */
static inline void excess_positions_512(const uint64_t* s,
                                        int target_x,
                                        uint64_t* out) noexcept {
  if (target_x < -512 || target_x > 512) {
    out[0] = out[1] = out[2] = out[3] = 0;
    out[4] = out[5] = out[6] = out[7] = 0;
    return;
  }

  for (int k = 0; k < 4; ++k) {
    target_x -= excess_positions_128(s + 2 * k, target_x, out + 2 * k);
  }
}

/**
 * @brief Find every strict record-low prefix excess position in a 128-bit
 * bitstring.
 *
 * @details A position i (1 <= i <= 128) is a strict record low if its prefix
 * excess is strictly smaller than every prefix excess at positions j < i.
 * Position 0 (empty prefix, excess 0) is always a record low vacuously but is
 * not represented in the output bitmask. Bit b of out[w] is set iff position
 * (w*64 + b + 1) is a strict record low.
 *
 * @param s 2 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param out 2 uint64_t words receiving the result bitmask.
 */
static inline void excess_record_lows_128(const uint64_t* s,
                                          uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;
  for (size_t i = 0; i < 128; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = static_cast<int>((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur < best) {
      best = cur;
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
}

/**
 * @brief Byte-LUT accelerated record-low scan for 128-bit blocks.
 *
 * @details Uses precomputed per-byte minimum and delta LUTs to skip whole
 * bytes that cannot contain a record low. When a byte's local minimum is
 * below the running global best, the byte is scanned bit-by-bit; otherwise
 * the byte is skipped and only the running excess is updated.
 */
static inline void excess_record_lows_128_byte_lut(const uint64_t* s,
                                                   uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;
  for (size_t byte_idx = 0; byte_idx < 16; ++byte_idx) {
    const size_t bit_base = byte_idx * 8;
    const uint8_t byte =
        static_cast<uint8_t>((s[bit_base >> 6] >> (bit_base & 63)) & 0xFFu);
    const int byte_min = excess_byte_min_lut[byte];
    if (cur + byte_min < best) {
      for (size_t i = 0; i < 8; ++i) {
        const int bit = static_cast<int>((byte >> i) & 1u);
        cur += bit ? +1 : -1;
        if (cur < best) {
          best = cur;
          const size_t pos = bit_base + i;
          out[pos >> 6] |= (uint64_t{1} << (pos & 63));
        }
      }
    } else {
      cur += excess_byte_delta_lut[byte];
    }
  }
}

/**
 * @brief LUT-based record-low scan for 128-bit blocks (branchless per byte).
 *
 * @details Precomputed per-byte LUT: for each byte value and each possible
 * threshold (cur - best, clamped to 0..7), returns a bitmask of positions
 * within the byte whose excess is strictly below the threshold. The LUT
 * eliminates all bit-by-bit scanning and branching inside the hot loop.
 */
static inline void excess_record_lows_128_lut(const uint64_t* s,
                                              uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;
  for (size_t byte_idx = 0; byte_idx < 16; ++byte_idx) {
    const size_t bit_base = byte_idx * 8;
    const uint8_t byte =
        static_cast<uint8_t>((s[bit_base >> 6] >> (bit_base & 63)) & 0xFFu);
    const int gap = cur - best;
    const int idx = gap > 7 ? 7 : (gap < 0 ? 0 : gap);
    const uint8_t mask =
        excess_byte_record_lows_lut[byte][static_cast<size_t>(idx)];
    if (mask != 0) {
      // Recompute absolute excesses for masked positions to update cur/best.
      int local = 0;
      int local_best = 0;
      uint8_t local_mask = 0;
      for (int bit = 0; bit < 8; ++bit) {
        local += ((byte >> bit) & 1u) ? 1 : -1;
        if (local < local_best) {
          local_best = local;
          local_mask |= static_cast<uint8_t>(1u << bit);
        }
      }
      // Only output positions whose absolute excess is < best.
      uint8_t out_mask = 0;
      local = 0;
      for (int bit = 0; bit < 8; ++bit) {
        local += ((byte >> bit) & 1u) ? 1 : -1;
        if (cur + local < best) {
          out_mask |= static_cast<uint8_t>(1u << bit);
          best = cur + local;
        }
      }
      if (out_mask != 0) {
        const uint64_t word = static_cast<uint64_t>(out_mask)
                              << (bit_base & 63);
        out[bit_base >> 6] |= word;
      }
    }
    cur += excess_byte_delta_lut[byte];
  }
}

#ifdef PIXIE_AVX2_SUPPORT
/**
 * @brief AVX2 record-low scan for 128-bit blocks.
 *
 * @details Expands the 128-bit block into eight 16-bit chunks, computes
 * vector prefix sums, then extracts record lows with a running scalar
 * minimum. This avoids the unpredictable branch in the scalar loop.
 */
static inline void excess_record_lows_128_avx2(const uint64_t* s,
                                               uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;

  const __m256i masks = excess_bit_masks_16x_i16();
  const __m256i zero = _mm256_setzero_si256();
  const __m256i pos = _mm256_set1_epi16(1);
  const __m256i neg = _mm256_set1_epi16(-1);

  for (size_t chunk = 0; chunk < 8; ++chunk) {
    const size_t chunk_bit = chunk * 16;
    const uint16_t bits =
        chunk < 4
            ? static_cast<uint16_t>((s[0] >> (chunk * 16)) & 0xFFFFu)
            : static_cast<uint16_t>((s[1] >> ((chunk - 4) * 16)) & 0xFFFFu);

    const __m256i vb = _mm256_set1_epi16(static_cast<int16_t>(bits));
    const __m256i m = _mm256_and_si256(vb, masks);
    const __m256i is_zero = _mm256_cmpeq_epi16(m, zero);
    const __m256i steps = _mm256_blendv_epi8(pos, neg, is_zero);
    const __m256i pref_rel = excess_prefix_sum_16x_i16(steps);
    const __m256i pref_abs = _mm256_add_epi16(
        pref_rel, _mm256_set1_epi16(static_cast<int16_t>(cur)));

    alignas(32) int16_t vals[16];
    _mm256_store_si256(reinterpret_cast<__m256i*>(vals), pref_abs);

    for (size_t lane = 0; lane < 16; ++lane) {
      const int val = vals[lane];
      if (val < best) {
        best = val;
        const size_t pos_idx = chunk_bit + lane;
        out[pos_idx >> 6] |= (uint64_t{1} << (pos_idx & 63));
      }
    }

    cur += 2 * static_cast<int>(std::popcount(bits)) - 16;
  }
}

/**
 * @brief AVX2 nibble-LUT record-low scan for 128-bit blocks.
 *
 * @details Loads 128 bits as 32 nibbles in one AVX2 register.
 * Uses precomputed LUTs to get per-nibble excess values, prefix-sums
 * them with 32-lane i8 operations, then compares against the running
 * scalar minimum. Record-low positions are written into the output
 * bitmask bit-by-bit.
 */
static inline void excess_record_lows_128_nibble_lut(const uint64_t* s,
                                                     uint64_t* out) noexcept {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;

  const __m256i vdelta = excess_lut_delta;
  const __m256i vmin = excess_lut_min;

  __m256i nibbles = excess_nibbles_128_avx2(s);

  // Compute inclusive prefix sums of per-nibble total excess changes.
  __m256i ps = _mm256_shuffle_epi8(vdelta, nibbles);
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 1));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 2));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 4));
  ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 8));

  __m128i ps_lo = _mm256_castsi256_si128(ps);
  __m128i ps_hi = _mm256_extracti128_si256(ps, 1);
  __m128i carry =
      _mm_set1_epi8(static_cast<int8_t>(_mm_extract_epi8(ps_lo, 15)));
  ps_hi = _mm_add_epi8(ps_hi, carry);
  ps = _mm256_inserti128_si256(_mm256_castsi128_si256(ps_lo), ps_hi, 1);

  // Exclusive prefix: shift in zero at start.
  __m256i b = _mm256_permute2x128_si256(ps, ps, 0x08);
  __m256i excl_ps = _mm256_alignr_epi8(ps, b, 15);

  // Local minima relative to each nibble start.
  __m256i local_min = _mm256_shuffle_epi8(vmin, nibbles);

  alignas(32) int8_t excl[32];
  _mm256_store_si256(reinterpret_cast<__m256i*>(excl), excl_ps);

  alignas(32) int8_t nibble_min[32];
  _mm256_store_si256(reinterpret_cast<__m256i*>(nibble_min), local_min);

  alignas(32) int8_t nibble_vals[32];
  _mm256_store_si256(reinterpret_cast<__m256i*>(nibble_vals), nibbles);

  for (int n = 0; n < 32; ++n) {
    const int nibble_base = cur + excl[n];
    const int nibble_best = nibble_base + nibble_min[n];
    if (nibble_best < best) {
      // This nibble contains at least one record low — scan bit-by-bit.
      const uint8_t nibble = static_cast<uint8_t>(nibble_vals[n]);
      int local = 0;
      for (int bit = 0; bit < 4; ++bit) {
        local += ((nibble >> bit) & 1u) ? 1 : -1;
        const int val = nibble_base + local;
        if (val < best) {
          best = val;
          const size_t pos = static_cast<size_t>(n) * 4 + bit;
          out[pos >> 6] |= (uint64_t{1} << (pos & 63));
        }
      }
    }
  }
}
#endif

/**
 * @brief Find every strict record-low prefix excess position in a 512-bit
 * bitstring.
 *
 * @details A position i (1 <= i <= 512) is a strict record low if its prefix
 * excess is strictly smaller than every prefix excess at positions j < i.
 * Position 0 is not represented in the output bitmask. Bit b of out[w] is set
 * iff position (w*64 + b + 1) is a strict record low.
 *
 * @param s 8 little-endian uint64_t words (bit 0 of s[0] is the first bit).
 * @param out 8 uint64_t words receiving the result bitmask.
 */
static inline void excess_record_lows_512(const uint64_t* s,
                                          uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;
  int best = 0;
  int global = 0;
  for (int k = 0; k < 4; ++k) {
    const uint64_t* block = s + 2 * k;
    uint64_t* block_out = out + 2 * k;
    for (size_t i = 0; i < 128; ++i) {
      const uint64_t w = block[i >> 6];
      const int bit = static_cast<int>((w >> (i & 63)) & 1ull);
      global += bit ? +1 : -1;
      if (global < best) {
        best = global;
        block_out[i >> 6] |= (uint64_t{1} << (i & 63));
      }
    }
  }
}

/**
 * @brief Calculates 32 bit ranks of every 8th bit, result is stored as 32

 * * 8-bit integers.
 * @details Prefix sums are computed modulo 256 (uint8_t
 * wraparound).
 *
 * @param x Pointer to 32 input 8-bit integers
 * @param
 * result Pointer to store the resulting 32 8-bit integers
 */
static inline void rank_32x8(const uint8_t* x, uint8_t* result) {
#ifdef PIXIE_AVX512_SUPPORT
  // Step 1: Calculate popcount of each byte
  popcount_32x8(x, result);
  __m256i prefix_sums = _mm256_loadu_si256((__m256i const*)result);
  const __m256i zero = _mm256_setzero_si256();

  prefix_sums = _mm256_add_epi8(prefix_sums,
                                _mm256_alignr_epi8(prefix_sums, zero, 16 - 1));
  prefix_sums = _mm256_add_epi8(prefix_sums,
                                _mm256_alignr_epi8(prefix_sums, zero, 16 - 2));
  prefix_sums = _mm256_add_epi8(prefix_sums,
                                _mm256_alignr_epi8(prefix_sums, zero, 16 - 4));
  prefix_sums = _mm256_add_epi8(prefix_sums,
                                _mm256_alignr_epi8(prefix_sums, zero, 16 - 8));

  // At this point we have prefix sums for two halfs, the last step is to
  // extract 16-th value and add it to the whole second half
  __m128i low_lane = _mm256_extracti128_si256(prefix_sums, 0);
  __m128i high_lane = _mm256_extracti128_si256(prefix_sums, 1);
  auto last_val_low = _mm_extract_epi8(low_lane, 15);
  __m128i add_to_high = _mm_set1_epi8(last_val_low);
  high_lane = _mm_add_epi8(high_lane, add_to_high);
  prefix_sums = _mm256_set_m128i(high_lane, low_lane);
  _mm256_storeu_epi8(result, prefix_sums);
#else
  // Scalar fallback implementation
  result[0] = std::popcount(x[0]);
  for (size_t i = 1; i < 32; ++i) {
    result[i] = std::popcount(x[i]) + result[i - 1];
  }
#endif
}
