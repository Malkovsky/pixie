#pragma once

#include <pixie/bits.h>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#if defined(PIXIE_AVX512_SUPPORT) || defined(PIXIE_AVX2_SUPPORT)
#include <immintrin.h>
#endif

namespace pixie::detail::compressed {

// Return the first endpoint >= key, or count. Unaligned and partial arrays are
// supported; masked loads never read beyond their logical storage. count is at
// most 64 / sizeof(Endpoint), and key fits the unsigned endpoint type.
template <class Endpoint>
inline std::size_t first_metadata_ge(const Endpoint* endpoints,
                                     std::size_t count,
                                     std::uint64_t key) {
  static_assert(std::is_unsigned_v<Endpoint> &&
                (sizeof(Endpoint) == 2 || sizeof(Endpoint) == 4 ||
                 sizeof(Endpoint) == 8));
  assert(count <= 64 / sizeof(Endpoint));
  assert(key <= std::numeric_limits<Endpoint>::max());
#ifdef PIXIE_AVX512_SUPPORT
  if constexpr (sizeof(Endpoint) == 2) {
    const auto valid = static_cast<__mmask32>((std::uint64_t{1} << count) - 1);
    const auto values = _mm512_maskz_loadu_epi16(valid, endpoints);
    const std::uint32_t matches =
        valid & _mm512_cmpge_epu16_mask(values, _mm512_set1_epi16(key));
    return matches == 0 ? count : std::countr_zero(matches);
  } else if constexpr (sizeof(Endpoint) == 4) {
    const auto valid = static_cast<__mmask16>((std::uint64_t{1} << count) - 1);
    const auto values = _mm512_maskz_loadu_epi32(valid, endpoints);
    const std::uint32_t matches =
        valid & _mm512_cmpge_epu32_mask(values, _mm512_set1_epi32(key));
    return matches == 0 ? count : std::countr_zero(matches);
  } else {
    const auto valid = static_cast<__mmask8>((std::uint64_t{1} << count) - 1);
    const auto values = _mm512_maskz_loadu_epi64(valid, endpoints);
    const std::uint32_t matches =
        valid & _mm512_cmpge_epu64_mask(values, _mm512_set1_epi64(key));
    return matches == 0 ? count : std::countr_zero(matches);
  }
#elif defined(PIXIE_AVX2_SUPPORT)
  constexpr std::size_t lanes = 32 / sizeof(Endpoint);
  std::size_t begin = 0;
  while (count - begin >= lanes) {
    const auto values =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(endpoints + begin));
    std::uint32_t matches;
    if constexpr (sizeof(Endpoint) == 2) {
      const auto sign = _mm256_set1_epi16(-32768);
      const auto less =
          _mm256_cmpgt_epi16(_mm256_xor_si256(_mm256_set1_epi16(key), sign),
                             _mm256_xor_si256(values, sign));
      matches = ~static_cast<std::uint32_t>(_mm256_movemask_epi8(less));
      if (matches != 0) {
        return begin + std::countr_zero(matches) / 2;
      }
    } else if constexpr (sizeof(Endpoint) == 4) {
      const auto sign =
          _mm256_set1_epi32(std::numeric_limits<std::int32_t>::min());
      const auto less =
          _mm256_cmpgt_epi32(_mm256_xor_si256(_mm256_set1_epi32(key), sign),
                             _mm256_xor_si256(values, sign));
      matches = ~static_cast<std::uint32_t>(
                    _mm256_movemask_ps(_mm256_castsi256_ps(less))) &
                0xff;
      if (matches != 0) {
        return begin + std::countr_zero(matches);
      }
    } else {
      const auto sign =
          _mm256_set1_epi64x(std::numeric_limits<std::int64_t>::min());
      const auto less =
          _mm256_cmpgt_epi64(_mm256_xor_si256(_mm256_set1_epi64x(key), sign),
                             _mm256_xor_si256(values, sign));
      matches = ~static_cast<std::uint32_t>(
                    _mm256_movemask_pd(_mm256_castsi256_pd(less))) &
                0xf;
      if (matches != 0) {
        return begin + std::countr_zero(matches);
      }
    }
    begin += lanes;
  }
  for (; begin < count; ++begin) {
    if (endpoints[begin] >= key) {
      return begin;
    }
  }
  return count;
#else
  for (std::size_t i = 0; i < count; ++i) {
    if (endpoints[i] >= key) {
      return i;
    }
  }
  return count;
#endif
}

template <class Endpoint>
inline std::size_t grouped_lower_bound(const Endpoint* endpoints,
                                       std::size_t count,
                                       std::uint64_t key) {
  constexpr std::size_t lanes = 64 / sizeof(Endpoint);
  std::size_t low = 0;
  std::size_t high = (count + lanes - 1) / lanes;
  while (low < high) {
    const auto mid = low + (high - low) / 2;
    const auto end = std::min(count, (mid + 1) * lanes);
    if (endpoints[end - 1] < key) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  const auto begin = low * lanes;
  return begin >= count
             ? count
             : begin + first_metadata_ge(endpoints + begin,
                                         std::min(lanes, count - begin), key);
}

}  // namespace pixie::detail::compressed
