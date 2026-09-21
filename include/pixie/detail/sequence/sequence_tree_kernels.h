#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

#if SIZE_MAX == UINT64_MAX && (defined(__AVX512F__) || defined(__AVX2__))
#include <immintrin.h>
#endif

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

/**
 * @brief Scalar child selection from unsigned cumulative exclusive ends.
 * @details Linear reference and fallback for node_select; retains no storage.
 * @param ends Nondecreasing stored ends, excluding the implicit last child end.
 * The span must remain valid for this call; no alignment or padding is
 * required.
 * @param index Zero-based position to locate, with the full size_t range
 * accepted.
 * @return First position whose end is strictly greater than index, or
 * ends.size() for the implicit last child. Equivalently, the number of ends <=
 * index. An empty span, including null backing storage, returns zero without
 * access.
 * @pre ends is nondecreasing; violating this precondition is unsupported.
 */
inline std::size_t node_select_scalar(std::span<const std::size_t> ends,
                                      std::size_t index) noexcept {
  std::size_t i = 0;
  while (i < ends.size() && ends[i] <= index) {
    ++i;
  }
  return i;
}

/**
 * @brief Compile-time SIMD dispatch for node_select_scalar's exact contract.
 * @param ends Nondecreasing exclusive ends, omitting the implicit last child.
 * Only [0, ends.size()) is read; the span needs no alignment or extra padding.
 * @param index Unsigned zero-based position, including values with the top bit
 * set.
 * @return First end > index, or ends.size() if none; zero without memory access
 * for an empty span. Equality selects the following child.
 * @pre ends remains valid during the call and is nondecreasing.
 * @details Uses only AVX-512F, otherwise AVX2, when size_t is 64 bits; all
 * other configurations use the scalar control. No backing storage is retained.
 */
inline std::size_t node_select(std::span<const std::size_t> ends,
                               std::size_t index) noexcept {
  if (ends.empty()) {
    return 0;
  }
#if SIZE_MAX == UINT64_MAX && defined(__AVX512F__)
  const auto query =
      _mm512_set1_epi64(std::bit_cast<long long>(std::uint64_t{index}));
  for (std::size_t i = 0; i < ends.size();) {
    const auto remaining = ends.size() - i;
    const auto count = remaining < 8 ? remaining : 8;
    const auto active = static_cast<__mmask8>((1u << count) - 1);
    const auto values = _mm512_maskz_loadu_epi64(active, ends.data() + i);
    const auto greater = static_cast<unsigned>(
        _mm512_mask_cmp_epu64_mask(active, values, query, _MM_CMPINT_GT));
    if (greater != 0) {
      return i + std::countr_zero(greater);
    }
    i += count;
  }
  return ends.size();
#elif SIZE_MAX == UINT64_MAX && defined(__AVX2__)
  const auto sign = _mm256_set1_epi64x(std::numeric_limits<long long>::min());
  const auto query = _mm256_xor_si256(
      _mm256_set1_epi64x(std::bit_cast<long long>(std::uint64_t{index})), sign);
  std::size_t i = 0;
  for (; ends.size() - i >= 4; i += 4) {
    const auto values =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ends.data() + i));
    const auto greater =
        static_cast<unsigned>(_mm256_movemask_pd(_mm256_castsi256_pd(
            _mm256_cmpgt_epi64(_mm256_xor_si256(values, sign), query))));
    if (greater != 0) {
      return i + std::countr_zero(greater);
    }
  }
  return i + node_select_scalar(ends.subspan(i), index);
#else
  return node_select_scalar(ends, index);
#endif
}

}  // namespace pixie::detail::sequence
/// @endcond
