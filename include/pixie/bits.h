#pragma once

#include <immintrin.h>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>

#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__) && \
    defined(__AVX512BW__)
#define PIXIE_AVX512_SUPPORT
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

#if defined(__BMI2__) && !defined(PIXIE_DISABLE_BMI2)
#define PIXIE_BMI2_SUPPORT
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

#ifndef PIXIE_BMI2_SUPPORT
struct PixieSelectByteLut {
  uint8_t popcounts[256];
  uint8_t select[256][8];

  constexpr PixieSelectByteLut() : popcounts{}, select{} {
    for (int byte = 0; byte < 256; ++byte) {
      for (int rank = 0; rank < 8; ++rank) {
        select[byte][rank] = 8;
      }

      int count = 0;
      for (int bit = 0; bit < 8; ++bit) {
        if (((byte >> bit) & 1) != 0) {
          select[byte][count++] = static_cast<uint8_t>(bit);
        }
      }
      popcounts[byte] = static_cast<uint8_t>(count);
    }
  }
};

static inline constexpr PixieSelectByteLut pixie_select_byte_lut;

static inline uint64_t select_64_no_bmi2(uint64_t x, uint64_t rank) {
  uint64_t offset = 0;

  uint64_t count = std::popcount(static_cast<uint32_t>(x));
  if (rank >= count) {
    rank -= count;
    x >>= 32;
    offset += 32;
  }

  count = std::popcount(static_cast<uint16_t>(x));
  if (rank >= count) {
    rank -= count;
    x >>= 16;
    offset += 16;
  }

  const auto low_byte = static_cast<uint8_t>(x);
  count = pixie_select_byte_lut.popcounts[low_byte];
  if (rank >= count) {
    rank -= count;
    x >>= 8;
    offset += 8;
  }

  return offset + pixie_select_byte_lut.select[static_cast<uint8_t>(x)][rank];
}
#endif

/**
 * @brief Return position of @p rank 1 bit in @p x
 */
static inline uint64_t select_64(uint64_t x, uint64_t rank) {
#ifdef PIXIE_BMI2_SUPPORT
  return std::countr_zero(_pdep_u64(1ull << rank, x));
#else
  return select_64_no_bmi2(x, rank);
#endif
}

template <bool Invert>
static inline uint64_t select_512_word_count(uint64_t word) {
  if constexpr (Invert) {
    return std::popcount(~word);
  } else {
    return std::popcount(word);
  }
}

template <bool Invert>
static inline uint64_t select_512_selected_word(uint64_t word) {
  if constexpr (Invert) {
    return ~word;
  } else {
    return word;
  }
}

template <bool Invert>
static inline uint64_t select_512_scalar_impl(const uint64_t* x,
                                              uint64_t rank) {
  for (size_t i = 0; i < 8; ++i) {
    const uint64_t count = select_512_word_count<Invert>(x[i]);
    if (rank < count) {
      return i * 64 + select_64(select_512_selected_word<Invert>(x[i]), rank);
    }
    rank -= count;
  }
  return 512;
}

#ifdef PIXIE_AVX2_SUPPORT
template <bool Invert>
static inline void select_512_avx2_counts(const uint64_t* x, uint64_t* counts) {
  const __m256i low_mask = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();
  const __m256i sixty_four = _mm256_set1_epi64x(64);

  for (int half = 0; half < 2; ++half) {
    const __m256i words =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + 4 * half));

    const __m256i low_nibbles = _mm256_and_si256(words, low_mask);
    const __m256i high_nibbles =
        _mm256_and_si256(_mm256_srli_epi16(words, 4), low_mask);
    const __m256i byte_counts =
        _mm256_add_epi8(_mm256_shuffle_epi8(lookup_popcount_4, low_nibbles),
                        _mm256_shuffle_epi8(lookup_popcount_4, high_nibbles));
    __m256i word_counts = _mm256_sad_epu8(byte_counts, zero);
    if constexpr (Invert) {
      word_counts = _mm256_sub_epi64(sixty_four, word_counts);
    }
    _mm256_store_si256(reinterpret_cast<__m256i*>(counts + 4 * half),
                       word_counts);
  }
}

template <bool Invert>
static inline uint64_t select_512_avx2_impl(const uint64_t* x, uint64_t rank) {
  alignas(32) uint64_t counts[8];
  select_512_avx2_counts<Invert>(x, counts);

  for (size_t i = 0; i < 8; ++i) {
    if (rank < counts[i]) {
      return i * 64 + select_64(select_512_selected_word<Invert>(x[i]), rank);
    }
    rank -= counts[i];
  }
  return 512;
}
#endif

#ifdef PIXIE_AVX512_SUPPORT
static inline __m512i select_512_avx512_prefix_sum_u64(__m512i prefix) {
  const __m512i idx_shift1 = _mm512_set_epi64(6, 5, 4, 3, 2, 1, 0, 0);
  const __m512i idx_shift2 = _mm512_set_epi64(5, 4, 3, 2, 1, 0, 0, 0);
  const __m512i idx_shift4 = _mm512_set_epi64(3, 2, 1, 0, 0, 0, 0, 0);

  __m512i tmp = _mm512_maskz_permutexvar_epi64(0xFE, idx_shift1, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xFC, idx_shift2, prefix);
  prefix = _mm512_add_epi64(prefix, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xF0, idx_shift4, prefix);
  return _mm512_add_epi64(prefix, tmp);
}

static inline uint64_t select_512_avx512_previous_prefix(__m512i prefix,
                                                         uint32_t lane) {
  if (lane == 0) {
    return 0;
  }
  const __m512i idx_previous =
      _mm512_set1_epi64(static_cast<int64_t>(lane - 1));
  const __m512i previous_vec = _mm512_permutexvar_epi64(idx_previous, prefix);
  return static_cast<uint64_t>(
      _mm_cvtsi128_si64(_mm512_castsi512_si128(previous_vec)));
}

template <bool Invert>
static inline uint64_t select_512_avx512_impl(const uint64_t* x,
                                              uint64_t rank) {
  const __m512i words = _mm512_loadu_epi64(x);
  __m512i prefix = _mm512_popcnt_epi64(words);
  if constexpr (Invert) {
    prefix = _mm512_sub_epi64(_mm512_set1_epi64(64), prefix);
  }
  prefix = select_512_avx512_prefix_sum_u64(prefix);

  const __mmask8 mask = _mm512_cmpgt_epu64_mask(
      prefix, _mm512_set1_epi64(static_cast<int64_t>(rank)));
  const uint32_t lane = std::countr_zero(static_cast<uint32_t>(mask));
  const uint64_t previous = select_512_avx512_previous_prefix(prefix, lane);
  return lane * 64 +
         select_64(select_512_selected_word<Invert>(x[lane]), rank - previous);
}
#endif

/**
 * @brief Return position of @p rank 1 bit in @p x.
 * @details Uses AVX-512, then AVX2, then scalar fallback. The 64-bit in-word
 * select uses BMI2 PDEP unless PIXIE_DISABLE_BMI2 is defined or BMI2 is not
 * available.
 */
static inline uint64_t select_512(const uint64_t* x, uint64_t rank) {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<false>(x, rank);
#elif defined(PIXIE_AVX2_SUPPORT)
  return select_512_avx2_impl<false>(x, rank);
#else
  return select_512_scalar_impl<false>(x, rank);
#endif
}

/**
 * @brief Return position of @p rank0 0 bit in @p x.
 * @details select_512 with bit inversion.
 */
static inline uint64_t select0_512(const uint64_t* x, uint64_t rank0) {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<true>(x, rank0);
#elif defined(PIXIE_AVX2_SUPPORT)
  return select_512_avx2_impl<true>(x, rank0);
#else
  return select_512_scalar_impl<true>(x, rank0);
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
static inline const __m256i excess_lut_pack_multiplier =
    _mm256_set1_epi16(0x1001);
static inline const __m256i excess_lut_bit0 = _mm256_set1_epi8(1);
static inline const __m256i excess_lut_bit1 = _mm256_set1_epi8(2);
static inline const __m256i excess_lut_bit2 = _mm256_set1_epi8(4);
static inline const __m256i excess_lut_bit3 = _mm256_set1_epi8(8);
static inline const __m128i excess_lut_nibble_mask = _mm_set1_epi8(0x0F);
// clang-format on
#endif

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
  const __m128i vnibble_mask = excess_lut_nibble_mask;

  const int d = 2 * target_x - block_delta;
  if (d < -128 || d > 128) {
    return block_delta;
  }

  __m128i word_vec = _mm_loadu_si128((const __m128i*)s);
  __m128i lo_nibbles = _mm_and_si128(word_vec, vnibble_mask);
  __m128i hi_nibbles = _mm_and_si128(_mm_srli_epi16(word_vec, 4), vnibble_mask);

  __m128i unpack_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
  __m128i unpack_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);

  __m256i nibbles =
      _mm256_inserti128_si256(_mm256_castsi128_si256(unpack_lo), unpack_hi, 1);

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
