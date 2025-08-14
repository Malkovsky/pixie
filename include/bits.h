#pragma once

#include <cstdint>
#include <immintrin.h>
#include <numeric>

/**
 * TODO: compile/runtime checks fror SIMD
 */

/**
 * @brief Number of 1 bits in positions 0 .. count - 1
 * @details
 * Surprisingly one cannot just do (1 << L) - 1 for
 * L > 64 to produce mask of ones of length count.
 * The best we can do with a single instruction is (1 << L) - 1 for l=8k
 * using maskz_set1.
 *
 * To make mask for arbitrary L we use shldv instuction, it doesn't
 * really matter what epi is used the recipe is the following:
 * - produce mask with x * (count / x) ones and store in a
 * - produce mask with x * ((count / x) + 1) ones and store in b
 * - use shldv(a, b, count % x) to blend a and b producing requred mask
 *
 * Essentially shldv_epix(a, b, k) takes k high bits of b and (x - k) bits of a
 * and concats them (bits of b become lower bits of the result).
 *
 * The rest is standard, i.e. popcount_epi64 to perform popcount on
 * 64 bits and ten reduce_add to sum the result.
 */
uint64_t popcount_512(const uint64_t *x, uint64_t count) {
  __m512i a = _mm512_maskz_set1_epi64((1ull << ((count >> 6))) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i b = _mm512_maskz_set1_epi64((1ull << ((count >> 6) + 1)) - 1,
                                      std::numeric_limits<uint64_t>::max());
  __m512i mask = _mm512_shldv_epi64(a, b, _mm512_set1_epi64(count % 64));

  __m512i res = _mm512_load_epi64(x);
  res = _mm512_and_epi64(res, mask);
  __m512i cnt = _mm512_popcnt_epi64(res);
  return _mm512_reduce_add_epi64(cnt);
}

uint64_t select_512(uint64_t *x, uint64_t rank) { return 0; }