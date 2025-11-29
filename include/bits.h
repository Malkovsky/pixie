#pragma once

#include <immintrin.h>

#include <array>
#include <bit>
#include <cstdint>
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
const __m256i lookup_popcount_4 = _mm256_setr_epi8(
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

const __m256i mask_first_half = _mm256_setr_epi8(
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

inline uint64_t first_bits_mask(size_t num) {
  return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
}

/**
 * @brief Number of 1 bits in positions 0 .. count - 1
 * @details
 * Surprisingly one cannot just do (1 << L) - 1 for
 * L > 64 to produce mask of ones of length L.
 * The best we can do with a single instruction is (1 << L) - 1 for L=8k
 * using maskz_set1.
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
uint64_t rank_512(const uint64_t* x, uint64_t count) {
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
 */
uint64_t select_64(uint64_t x, uint64_t rank) {
  return _tzcnt_u64(_pdep_u64(1ull << rank, x));
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
uint64_t select_512(const uint64_t* x, uint64_t rank) {
#ifdef PIXIE_AVX512_SUPPORT

  __m512i res = _mm512_loadu_epi64(x);
  std::array<uint64_t, 8> counts;
  _mm512_storeu_epi64(counts.data(), _mm512_popcnt_epi64(res));

  size_t i = 0;
  while (i < 8 && counts[i] <= rank) {
    rank -= counts[i++];
  }
  return i * 64 + select_64(x[i], rank);

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

/**
 * @brief Compare 4 64-bit numbers of @p x with @p y and
 * return the length of the prefix where @p y is less then @p x
 */
uint16_t lower_bound_4x64(const uint64_t* x, uint64_t y) {
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
 * @brief Compare 8 64-bit numbers of @p x with @p y and
 * return the length of the prefix where @p y is less then @p x
 */
uint16_t lower_bound_8x64(const uint64_t* x, uint64_t y) {
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
 * @brief Compare 32 16-bit numbers of @p x with @p y and
 * return the count of numbers where @p x is less then @p y
 */
uint16_t lower_bound_32x16(const uint16_t* x, uint64_t y) {
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
void popcount_64x4(const uint8_t* x, uint8_t* result) {
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
void popcount_32x8(const uint8_t* x, uint8_t* result) {
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

/**
 * @brief Calculates 32 bit ranks of every 8th bit, result is stored as 32
 * 8-bit integers.
 *
 * @param x Pointer to 32 input 8-bit integers
 * @param result Pointer to store the resulting 32 8-bit integers
 */
void rank_32x8(const uint8_t* x, uint8_t* result) {
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

/**
 * @brief Efficiently searches for the first occurrence of a 16-bit value in
 * the range [@p begin, @p end_excl) using AVX2 when available.
 * @details Loads 16 consecutive int16_t elements (256 bits) per iteration.
 * Compares them against the @p target value using vectorized equality.
 * If any match is found, extracts the index of the first matching lane from
 * the comparison mask. Falls back to a scalar tail loop for leftover
 * elements, or to a fully scalar search if AVX2 is not supported.
 * @returns The index of the first match, or @p npos if the value is not found.
 */
static inline size_t find_forward_equal_i16_avx2(const int16_t* arr,
                                                 const size_t& begin,
                                                 const size_t& end_excl,
                                                 const int16_t& target,
                                                 const size_t& npos) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  static constexpr size_t STEP = 16;
  __m256i vtarget = _mm256_set1_epi16(target);
  size_t i = begin;
  size_t n = end_excl;
  for (; i + STEP <= n; i += STEP) {
    unsigned mask = _mm256_movemask_epi8(_mm256_cmpeq_epi16(
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + i)),
        vtarget));
    if (mask) {
      return i + (std::countr_zero(mask) >> 1);
    }
  }
  for (; i < n; ++i) {
    if (arr[i] == target) {
      return i;
    }
  }
#else
  for (size_t i = begin; i < end_excl; ++i) {
    if (arr[i] == target) {
      return i;
    }
  }
#endif
  return npos;
}

/**
 * @brief Performs a backward search for a 16-bit value in a given range.
 * @details Scans the array segment [@p begin .. @p end_incl] from right to
 * left.
 * If AVX2 is available, processes data in 256-bit blocks (16 × int16_t) using
 * vectorized equality comparison for higher throughput. Falls back to a
 * scalar backward scan when AVX2 is not supported. Returns the index of the
 * rightmost occurrence of @p target, or @p npos if no match is found.
 */
static inline size_t find_backward_equal_i16_avx2(const int16_t* arr,
                                                  const size_t& begin,
                                                  const size_t& end_incl,
                                                  const int16_t& target,
                                                  const size_t& npos) noexcept {
  if (begin > end_incl) {
    return npos;
  }
#ifdef PIXIE_AVX2_SUPPORT
  static constexpr size_t STEP = 16;
  size_t len = end_incl + 1 - begin;
  size_t nblocks = len / STEP;
  __m256i vtarget = _mm256_set1_epi16(target);
  if (nblocks > 0) {
    size_t first_block = begin + (len % STEP);
    for (size_t p = first_block + (nblocks - 1) * STEP;;) {
      unsigned mask = _mm256_movemask_epi8(_mm256_cmpeq_epi16(
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(arr + p)),
          vtarget));
      if (mask) {
        return p + ((31u - std::countl_zero(mask)) >> 1);
      }
      if (p == first_block) {
        break;
      }
      p -= STEP;
    }

    for (size_t i = first_block; i > begin;) {
      --i;
      if (arr[i] == target) {
        return i;
      }
    }
  } else {
    for (size_t i = end_incl + 1; i > begin;) {
      --i;
      if (arr[i] == target) {
        return i;
      }
    }
  }
#else
  for (size_t i = end_incl + 1; i > begin;) {
    --i;
    if (arr[i] == target) {
      return i;
    }
  }
#endif
  return npos;
}
