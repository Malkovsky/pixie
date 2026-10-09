#pragma once

#include <pixie/rank_select/detail/compressed/packed_prefix_search.h>

namespace pixie::experimental {

// Four 14-bit values occupy seven bytes. Gather each group into a 64-bit
// lane; the eighth byte is ignored after the final 14-bit mask.
alignas(64) inline constexpr auto kPrefixGroupBytes = [] {
  std::array<std::uint8_t, 64> indexes{};
  for (std::size_t i = 0; i < 64; ++i) {
    indexes[i] = (8 + (i / 8) * 7 + i % 8) & 63;
  }
  return indexes;
}();

/**
 * @brief Experimental AVX-512 VBMI multishift search of packed rank prefixes.
 * @details Uses the production record layout and preconditions: a
 * nonhomogeneous superblock, a valid one-based local target, and its valid bit
 * count. Returns no result on unsupported targets. This kernel is not selected
 * by RankSelectSupportCompressed.
 */
template <bool Zero>
[[gnu::always_inline]] inline std::optional<
    detail::compressed::PackedPrefixLocation>
packed_prefix_search_multishift(const std::uint64_t* record,
                                std::uint64_t target,
                                std::size_t valid_bits) {
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI__)
  const auto packed = _mm512_loadu_si512(record);
  const auto groups = _mm512_permutexvar_epi8(
      _mm512_load_si512(kPrefixGroupBytes.data()), packed);
  auto ends =
      _mm512_and_si512(_mm512_multishift_epi64_epi8(
                           _mm512_set1_epi64(0x322a241c160e0800ULL), groups),
                       _mm512_set1_epi16(0x3fff));
  if constexpr (Zero) {
    const auto positions = _mm512_set_epi16(
        16384, 15872, 15360, 14848, 14336, 13824, 13312, 12800, 12288, 11776,
        11264, 10752, 10240, 9728, 9216, 8704, 8192, 7680, 7168, 6656, 6144,
        5632, 5120, 4608, 4096, 3584, 3072, 2560, 2048, 1536, 1024, 512);
    ends = _mm512_sub_epi16(
        _mm512_min_epu16(positions, _mm512_set1_epi16(valid_bits)), ends);
  }
  const std::uint32_t matches =
      _mm512_cmpge_epu16_mask(ends, _mm512_set1_epi16(target));
  assert(matches != 0);
  const std::uint32_t basic = std::countr_zero(matches);
  const auto indexes =
      _mm512_set1_epi32((basic << 16) | ((basic - 1) & 0xffff));
  const std::uint32_t bounds = _mm_cvtsi128_si32(
      _mm512_castsi512_si128(_mm512_permutexvar_epi16(indexes, ends)));
  return detail::compressed::PackedPrefixLocation{
      basic, basic == 0 ? 0 : bounds & 0xffff, bounds >> 16};
#else
  (void)record;
  (void)target;
  (void)valid_bits;
  return std::nullopt;
#endif
}

}  // namespace pixie::experimental
