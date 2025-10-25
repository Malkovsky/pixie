#pragma once

#include <array>
#include <cstdint>
#include <immintrin.h>
#include <numeric>


#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__) && defined(__AVX512BW__)
#define AVX512SUPP
#endif


inline uint64_t first_bits_mask(size_t num) {
	return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
}


/**
 * TODO: compile/runtime checks fror SIMD
 */

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
#ifdef AVX512SUPP

uint64_t rank_512(const uint64_t *x, uint64_t count) {
  __m512i a = _mm512_maskz_set1_epi64((1ull << ((count >> 6))) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i b = _mm512_maskz_set1_epi64((1ull << ((count >> 6) + 1)) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i mask = _mm512_shldv_epi64(a, b, _mm512_set1_epi64(count % 64));

  __m512i res = _mm512_loadu_epi64(x);
  res = _mm512_and_epi64(res, mask);
  __m512i cnt = _mm512_popcnt_epi64(res);
  return _mm512_reduce_add_epi64(cnt);
}

#else

uint64_t rank_512(const uint64_t *x, uint64_t count) {

	uint64_t last_uint = count >> 6;

	uint64_t pop_val = 0;

	for (int i = 0; i < last_uint; i++) {
		pop_val += std::popcount(x[i]);
	}

	uint64_t final = x[last_uint] & first_bits_mask(count & 63);

	pop_val += std::popcount(final);
	return pop_val;
}

#endif

/**
 * @brief Return position of @p rank 1 bit in @p x
 */
uint64_t select_64(uint64_t x, uint64_t rank) {
  auto a = _pdep_u64(1ull << rank, x);
  return _tzcnt_u64(a);
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
#ifdef AVX512SUPP

uint64_t select_512(const uint64_t *x, uint64_t rank) {
    __m512i res = _mm512_loadu_epi64(x);
  std::array<uint64_t, 8> counts;
  _mm512_storeu_epi64(counts.data(), _mm512_popcnt_epi64(res));

  size_t i = 0;
  while (i < 8 && counts[i] <= rank) {
    rank -= counts[i++];
  }
  return i * 64 + select_64(x[i], rank);
}

#else

uint64_t select_512(const uint64_t* bits, uint64_t rank) {
  size_t i = 0;
	int popcount = std::popcount(bits[0]);
	while (i < 7 && popcount <= rank) {
		rank -= popcount;
		popcount = std::popcount(bits[++i]);
	}
	return i * 64 + select_64(bits[i], rank);
}

#endif


/**
 * @brief Compare 4 64-bit numbers of @p x with @p y and 
 * return the length of the prefix where @p y is less then
 */
#ifdef AVX512SUPP

uint16_t cmpl_pref_len_256(const uint64_t* x, uint64_t y) {

  auto y_4 = _mm256_set1_epi64x(y);
  auto reg_256 = _mm256_loadu_epi64(x);
  auto cmp = _mm256_cmpge_epu64_mask(reg_256, y_4);

  return _tzcnt_u16(cmp);
}

#else

uint16_t cmpl_pref_len_256(const uint64_t* x, uint64_t y) {

  auto y_4 = _mm256_set1_epi64x(y);
  __m256i reg_256 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));

  const __m256i offset = _mm256_set1_epi64x(0x8000000000000000ULL);
  __m256i x_offset = _mm256_xor_si256(reg_256, offset);
  __m256i y_offset = _mm256_xor_si256(y_4, offset);
  auto mask = _mm256_movemask_epi8(_mm256_cmpgt_epi64(x_offset, _mm256_sub_epi64(y_offset, _mm256_set1_epi64x(1))));

  return _tzcnt_u32(mask) >> 3;
}

#endif


/**
 * @brief Compare 8 64-bit numbers of @p x with @p y and 
 * return the length of the prefix where @p y is less then @p x
 */
#ifdef AVX512SUPP

uint16_t cmpl_pref_len_512(const uint64_t* x, uint64_t y) {

  auto y_8 = _mm512_set1_epi64(y);
  auto reg_512 = _mm512_loadu_epi64(x);
  auto cmp = _mm512_cmpge_epu64_mask(reg_512, y_8);

  return _tzcnt_u16(cmp);
}

#else

uint16_t cmpl_pref_len_512(const uint64_t* x, uint64_t y) {

  uint16_t len = cmpl_pref_len_256(x, y);

  if (len < 4)
    return len;

  return len + cmpl_pref_len_256(x + 4, y);
}

#endif


/**
 * @brief Compare 32 16-bit numbers of @p x with @p y and 
 * return the count of numbers where @p x is less then @p y
 */
#ifdef AVX512SUPP

uint16_t cmpl_count_512(const uint16_t* x, uint64_t y) {

  auto y_32 = _mm512_set1_epi16(y);
  auto reg_512 = _mm512_loadu_epi16(x);
  auto cmp = _mm512_cmplt_epu16_mask(reg_512, y_32);
  return std::popcount(cmp);
}

#else

uint16_t cmpl_count_512(const uint16_t* x, uint64_t y) {

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
}

#endif