#pragma once

#include <array>
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
uint64_t rank_512(const uint64_t *x, uint64_t count) {
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
uint64_t select_512(const uint64_t *x, uint64_t rank) {
  __m512i res = _mm512_load_epi64(x);
  std::array<uint64_t, 8> counts;
  _mm512_store_epi64(counts.data(), _mm512_popcnt_epi64(res));

  size_t i = 0;
  while (i < 8 && counts[i] <= rank) {
    rank -= counts[i++];
  }
  return i * 64 + select_64(x[i], rank);
}