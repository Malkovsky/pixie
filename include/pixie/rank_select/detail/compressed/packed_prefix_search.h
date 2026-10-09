#pragma once

#include <pixie/bits.h>

#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>

#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI__)
#include <immintrin.h>
#endif

namespace pixie::detail::compressed {

struct PackedPrefixLocation {
  std::uint32_t basic;
  std::uint32_t before;
  std::uint32_t after;
};

// Two byte permutations gather every fourteen-bit field into a word lane.
// Offsets repeat 0, 6, 4, 2. Only those with offsets 6 and 4 need a third byte.
alignas(64) inline constexpr auto kPrefixLowBytes = [] {
  std::array<std::uint8_t, 64> indexes{};
  for (std::size_t i = 0; i < 32; ++i) {
    indexes[2 * i] = 8 + i * 14 / 8;
    indexes[2 * i + 1] = indexes[2 * i] + 1;
  }
  return indexes;
}();

alignas(64) inline constexpr auto kPrefixHighBytes = [] {
  std::array<std::uint8_t, 64> indexes{};
  for (std::size_t i = 0; i < 32; ++i) {
    indexes[2 * i] = (10 + i * 14 / 8) & 63;
  }
  return indexes;
}();

// The selected superblock must be nonhomogeneous, with uniformly packed local
// ends. This excludes the zero terminal alias. Unsupported targets return no
// result so callers retain the established scalar interpolation/correction.
template <bool Zero>
[[gnu::always_inline]] inline std::optional<PackedPrefixLocation>
packed_prefix_search(const std::uint64_t* record,
                     std::uint64_t target,
                     std::size_t valid_bits) {
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI__)
  const auto packed = _mm512_loadu_si512(record);
  const auto low = _mm512_permutexvar_epi8(
      _mm512_load_si512(kPrefixLowBytes.data()), packed);
  const auto high = _mm512_permutexvar_epi8(
      _mm512_load_si512(kPrefixHighBytes.data()), packed);
  const auto offsets = _mm512_set1_epi64(0x0002000400060000);
  auto ends = _mm512_and_si512(
      _mm512_or_si512(
          _mm512_srlv_epi16(low, offsets),
          _mm512_sllv_epi16(high,
                            _mm512_sub_epi16(_mm512_set1_epi16(16), offsets))),
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
  return PackedPrefixLocation{basic, basic == 0 ? 0 : bounds & 0xffff,
                              bounds >> 16};
#else
  (void)record;
  (void)target;
  (void)valid_bits;
  return std::nullopt;
#endif
}

}  // namespace pixie::detail::compressed
