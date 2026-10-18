#pragma once

#include <pixie/bits.h>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace pixie::detail::compressed {

// Select a zero-based occurrence in eight readable words. The caller guarantees
// that the occurrence exists. Keep the loaded block in registers through word
// localization; extracting the winner must not reload source memory.
template <bool Zero>
inline std::size_t select_loaded_block(const std::uint64_t* words,
                                       std::size_t rank) {
#ifdef PIXIE_AVX512_SUPPORT
  auto bits = _mm512_loadu_si512(words);
  if constexpr (Zero) {
    bits = _mm512_xor_si512(bits, _mm512_set1_epi64(-1));
  }
  const auto zero = _mm512_setzero_si512();
  auto prefixes = _mm512_popcnt_epi64(bits);
  prefixes = _mm512_add_epi64(prefixes, _mm512_alignr_epi64(prefixes, zero, 7));
  prefixes = _mm512_add_epi64(prefixes, _mm512_alignr_epi64(prefixes, zero, 6));
  prefixes = _mm512_add_epi64(prefixes, _mm512_alignr_epi64(prefixes, zero, 4));
  const auto candidates =
      _mm512_cmp_epu64_mask(prefixes, _mm512_set1_epi64(rank), _MM_CMPINT_GT);
  const auto lane = std::countr_zero(static_cast<unsigned>(candidates));
  const auto selector = _mm512_set1_epi64(lane);
  const auto word = static_cast<std::uint64_t>(_mm_cvtsi128_si64(
      _mm512_castsi512_si128(_mm512_permutexvar_epi64(selector, bits))));
  const auto end = static_cast<std::uint64_t>(_mm_cvtsi128_si64(
      _mm512_castsi512_si128(_mm512_permutexvar_epi64(selector, prefixes))));
  const auto before = end - std::popcount(word);
  return lane * 64 + ::select_64(word, rank - before);
#else
  // The established primitive supplies AVX2 and scalar dispatch on other
  // targets. Callers need no SIMD feature guards.
  if constexpr (Zero) {
    return ::select0_512(words, rank);
  } else {
    return ::select_512(words, rank);
  }
#endif
}

}  // namespace pixie::detail::compressed
