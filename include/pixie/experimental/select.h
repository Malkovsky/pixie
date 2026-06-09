#pragma once

#include <pixie/bits.h>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace pixie::experimental {

struct SelectByteLut {
  uint8_t popcounts[256];
  uint8_t select[256][8];

  constexpr SelectByteLut() : popcounts{}, select{} {
    for (int b = 0; b < 256; ++b) {
      for (int r = 0; r < 8; ++r) {
        select[b][r] = 8;
      }

      int rank = 0;
      for (int i = 0; i < 8; ++i) {
        if (((b >> i) & 1) != 0) {
          select[b][rank++] = static_cast<uint8_t>(i);
        }
      }
      popcounts[b] = static_cast<uint8_t>(rank);
    }
  }
};

inline constexpr SelectByteLut kSelectByteLut;

static inline uint64_t select_64_pdep(uint64_t x, uint64_t rank) noexcept {
  return ::select_64(x, rank);
}

static inline uint64_t select_64_vigna_broadword(uint64_t x,
                                                 uint64_t rank) noexcept {
  constexpr uint64_t kOnesStep4 = 0x1111111111111111ull;
  constexpr uint64_t kOnesStep8 = 0x0101010101010101ull;
  constexpr uint64_t kHighBitsStep8 = 0x8080808080808080ull;

  uint64_t sums = x;
  sums = sums - ((sums & (0xAull * kOnesStep4)) >> 1);
  sums = (sums & (0x3ull * kOnesStep4)) + ((sums >> 2) & (0x3ull * kOnesStep4));
  sums = (sums + (sums >> 4)) & (0xFull * kOnesStep8);

  const uint64_t byte_sums = sums * kOnesStep8;
  const uint64_t rank_bytes = rank * kOnesStep8;
  const uint64_t ge_rank =
      ((rank_bytes | kHighBitsStep8) - byte_sums) & kHighBitsStep8;
  const uint64_t byte_index = std::popcount(ge_rank);
  const uint64_t byte_rank =
      rank - (((byte_sums << 8) >> (8 * byte_index)) & 0xFFull);
  const auto byte = static_cast<uint8_t>(x >> (8 * byte_index));
  return 8 * byte_index + kSelectByteLut.select[byte][byte_rank];
}

static inline uint64_t select_64_byte_lut(uint64_t x, uint64_t rank) noexcept {
  for (uint64_t byte_index = 0; byte_index < 8; ++byte_index) {
    const auto byte = static_cast<uint8_t>(x >> (8 * byte_index));
    const uint8_t count = kSelectByteLut.popcounts[byte];
    if (rank < count) {
      return 8 * byte_index + kSelectByteLut.select[byte][rank];
    }
    rank -= count;
  }
  return 64;
}

static inline uint64_t select_64_broadword_lut(uint64_t x,
                                               uint64_t rank) noexcept {
  constexpr uint64_t kOnes = 0x0101010101010101ull;
  constexpr uint64_t kHighBits = 0x8080808080808080ull;

  uint64_t byte_counts = x - ((x >> 1) & 0x5555555555555555ull);
  byte_counts = (byte_counts & 0x3333333333333333ull) +
                ((byte_counts >> 2) & 0x3333333333333333ull);
  byte_counts = (byte_counts + (byte_counts >> 4)) & 0x0F0F0F0F0F0F0F0Full;

  const uint64_t prefix = byte_counts * kOnes;
  const uint64_t threshold = rank + 1;
  const uint64_t prefix_less_than_threshold =
      (kOnes * (127 + threshold) - (prefix & 0x7F7F7F7F7F7F7F7Full)) & ~prefix &
      kHighBits;
  const uint64_t byte_index = std::popcount(prefix_less_than_threshold);
  const uint64_t previous =
      byte_index == 0 ? 0 : ((prefix >> (8 * (byte_index - 1))) & 0xFFull);
  const auto byte = static_cast<uint8_t>(x >> (8 * byte_index));
  return 8 * byte_index + kSelectByteLut.select[byte][rank - previous];
}

static inline uint64_t select_64_binary_lut(uint64_t x,
                                            uint64_t rank) noexcept {
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
  count = kSelectByteLut.popcounts[low_byte];
  if (rank >= count) {
    rank -= count;
    x >>= 8;
    offset += 8;
  }

  return offset + kSelectByteLut.select[static_cast<uint8_t>(x)][rank];
}

static inline uint64_t select_64_no_pdep(uint64_t x, uint64_t rank) noexcept {
  return select_64_binary_lut(x, rank);
}

static inline uint64_t select_64_experimental_default(uint64_t x,
                                                      uint64_t rank) noexcept {
#ifdef PIXIE_EXPERIMENTAL_SELECT_NO_PDEP
  return select_64_no_pdep(x, rank);
#else
  return select_64_pdep(x, rank);
#endif
}

template <bool Invert>
static inline uint64_t select_512_word_count(uint64_t word) noexcept {
  if constexpr (Invert) {
    return std::popcount(~word);
  } else {
    return std::popcount(word);
  }
}

template <bool Invert>
static inline uint64_t select_512_select_word(uint64_t word) noexcept {
  if constexpr (Invert) {
    return ~word;
  } else {
    return word;
  }
}

template <uint64_t (*SelectWord)(uint64_t, uint64_t), bool Invert>
static inline uint64_t select_512_scalar_impl(const uint64_t* x,
                                              uint64_t rank) noexcept {
  for (size_t i = 0; i < 8; ++i) {
    const uint64_t count = select_512_word_count<Invert>(x[i]);
    if (rank < count) {
      const uint64_t word = select_512_select_word<Invert>(x[i]);
      return i * 64 + SelectWord(word, rank);
    }
    rank -= count;
  }
  return 512;
}

static inline uint64_t select_512_scalar_pdep(const uint64_t* x,
                                              uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_pdep, false>(x, rank);
}

static inline uint64_t select0_512_scalar_pdep(const uint64_t* x,
                                               uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_pdep, true>(x, rank);
}

static inline uint64_t select_512_scalar_byte_lut(const uint64_t* x,
                                                  uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_byte_lut, false>(x, rank);
}

static inline uint64_t select0_512_scalar_byte_lut(const uint64_t* x,
                                                   uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_byte_lut, true>(x, rank);
}

static inline uint64_t select_512_scalar_broadword_lut(const uint64_t* x,
                                                       uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_broadword_lut, false>(x, rank);
}

static inline uint64_t select0_512_scalar_broadword_lut(
    const uint64_t* x,
    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_broadword_lut, true>(x, rank);
}

static inline uint64_t select_512_scalar_binary_lut(const uint64_t* x,
                                                    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_binary_lut, false>(x, rank);
}

static inline uint64_t select0_512_scalar_binary_lut(const uint64_t* x,
                                                     uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_binary_lut, true>(x, rank);
}

static inline uint64_t select_512_scalar_vigna_broadword(
    const uint64_t* x,
    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_vigna_broadword, false>(x, rank);
}

static inline uint64_t select0_512_scalar_vigna_broadword(
    const uint64_t* x,
    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_vigna_broadword, true>(x, rank);
}

static inline uint64_t select_512_scalar_no_pdep(const uint64_t* x,
                                                 uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_no_pdep, false>(x, rank);
}

static inline uint64_t select0_512_scalar_no_pdep(const uint64_t* x,
                                                  uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_no_pdep, true>(x, rank);
}

static inline uint64_t select_512_scalar_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_experimental_default, false>(x, rank);
}

static inline uint64_t select0_512_scalar_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
  return select_512_scalar_impl<select_64_experimental_default, true>(x, rank);
}

#ifdef PIXIE_AVX2_SUPPORT
template <bool Invert>
static inline void select_512_avx2_counts(const uint64_t* x,
                                          uint64_t* counts) noexcept {
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

template <uint64_t (*SelectWord)(uint64_t, uint64_t), bool Invert>
static inline uint64_t select_512_avx2_impl(const uint64_t* x,
                                            uint64_t rank) noexcept {
  alignas(32) uint64_t counts[8];
  select_512_avx2_counts<Invert>(x, counts);

  for (size_t i = 0; i < 8; ++i) {
    if (rank < counts[i]) {
      const uint64_t word = select_512_select_word<Invert>(x[i]);
      return i * 64 + SelectWord(word, rank);
    }
    rank -= counts[i];
  }
  return 512;
}
#endif

static inline uint64_t select_512_avx2_pdep(const uint64_t* x,
                                            uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_pdep, false>(x, rank);
#else
  return select_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_pdep(const uint64_t* x,
                                             uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_pdep, true>(x, rank);
#else
  return select0_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx2_byte_lut(const uint64_t* x,
                                                uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_byte_lut, false>(x, rank);
#else
  return select_512_scalar_byte_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_byte_lut(const uint64_t* x,
                                                 uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_byte_lut, true>(x, rank);
#else
  return select0_512_scalar_byte_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx2_broadword_lut(const uint64_t* x,
                                                     uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_broadword_lut, false>(x, rank);
#else
  return select_512_scalar_broadword_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_broadword_lut(const uint64_t* x,
                                                      uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_broadword_lut, true>(x, rank);
#else
  return select0_512_scalar_broadword_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx2_binary_lut(const uint64_t* x,
                                                  uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_binary_lut, false>(x, rank);
#else
  return select_512_scalar_binary_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_binary_lut(const uint64_t* x,
                                                   uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_binary_lut, true>(x, rank);
#else
  return select0_512_scalar_binary_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx2_vigna_broadword(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_vigna_broadword, false>(x, rank);
#else
  return select_512_scalar_vigna_broadword(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_vigna_broadword(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_vigna_broadword, true>(x, rank);
#else
  return select0_512_scalar_vigna_broadword(x, rank);
#endif
}

static inline uint64_t select_512_avx2_no_pdep(const uint64_t* x,
                                               uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_no_pdep, false>(x, rank);
#else
  return select_512_scalar_no_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_no_pdep(const uint64_t* x,
                                                uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_no_pdep, true>(x, rank);
#else
  return select0_512_scalar_no_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx2_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_experimental_default, false>(x, rank);
#else
  return select_512_scalar_experimental_default(x, rank);
#endif
}

static inline uint64_t select0_512_avx2_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX2_SUPPORT
  return select_512_avx2_impl<select_64_experimental_default, true>(x, rank);
#else
  return select0_512_scalar_experimental_default(x, rank);
#endif
}

#ifdef PIXIE_AVX512_SUPPORT
static inline __m512i select_512_avx512_prefix_sum_u64(
    __m512i counts) noexcept {
  const __m512i idx_shift1 = _mm512_set_epi64(6, 5, 4, 3, 2, 1, 0, 0);
  const __m512i idx_shift2 = _mm512_set_epi64(5, 4, 3, 2, 1, 0, 0, 0);
  const __m512i idx_shift4 = _mm512_set_epi64(3, 2, 1, 0, 0, 0, 0, 0);

  __m512i tmp = _mm512_maskz_permutexvar_epi64(0xFE, idx_shift1, counts);
  counts = _mm512_add_epi64(counts, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xFC, idx_shift2, counts);
  counts = _mm512_add_epi64(counts, tmp);
  tmp = _mm512_maskz_permutexvar_epi64(0xF0, idx_shift4, counts);
  return _mm512_add_epi64(counts, tmp);
}

static inline uint64_t select_512_avx512_previous_prefix(
    __m512i prefix,
    uint32_t lane) noexcept {
  if (lane == 0) {
    return 0;
  }
  const __m512i idx_previous =
      _mm512_set1_epi64(static_cast<int64_t>(lane - 1));
  const __m512i previous_vec = _mm512_permutexvar_epi64(idx_previous, prefix);
  return static_cast<uint64_t>(
      _mm_cvtsi128_si64(_mm512_castsi512_si128(previous_vec)));
}

template <uint64_t (*SelectWord)(uint64_t, uint64_t), bool Invert>
static inline uint64_t select_512_avx512_impl(const uint64_t* x,
                                              uint64_t rank) noexcept {
  __m512i words = _mm512_loadu_epi64(x);
  __m512i prefix = _mm512_popcnt_epi64(words);
  if constexpr (Invert) {
    prefix = _mm512_sub_epi64(_mm512_set1_epi64(64), prefix);
  }

  prefix = select_512_avx512_prefix_sum_u64(prefix);

  const __mmask8 mask = _mm512_cmpgt_epu64_mask(
      prefix, _mm512_set1_epi64(static_cast<int64_t>(rank)));
  const uint32_t lane = _tzcnt_u32(static_cast<uint32_t>(mask));

  const uint64_t previous = select_512_avx512_previous_prefix(prefix, lane);
  const uint64_t word = select_512_select_word<Invert>(x[lane]);
  return lane * 64 + SelectWord(word, rank - previous);
}

template <uint64_t (*SelectWord)(uint64_t, uint64_t), bool Invert>
static inline uint64_t select_512_avx512_tail_impl(
    const uint64_t* x,
    uint64_t rank,
    uint32_t base_word,
    uint32_t word_count) noexcept {
  const __mmask8 active_mask =
      static_cast<__mmask8>((uint32_t{1} << word_count) - 1);
  const __m512i words = _mm512_maskz_loadu_epi64(active_mask, x + base_word);
  __m512i counts = _mm512_popcnt_epi64(words);
  if constexpr (Invert) {
    counts = _mm512_maskz_sub_epi64(active_mask, _mm512_set1_epi64(64), counts);
  }

  const __m512i prefix = select_512_avx512_prefix_sum_u64(counts);
  const __mmask8 mask = _mm512_mask_cmpgt_epu64_mask(
      active_mask, prefix, _mm512_set1_epi64(static_cast<int64_t>(rank)));
  if (mask == 0) {
    return 512;
  }

  const uint32_t lane = _tzcnt_u32(static_cast<uint32_t>(mask));
  const uint64_t previous = select_512_avx512_previous_prefix(prefix, lane);
  const uint64_t word = select_512_select_word<Invert>(x[base_word + lane]);
  return (base_word + lane) * 64 + SelectWord(word, rank - previous);
}

template <uint64_t (*SelectWord)(uint64_t, uint64_t),
          bool Invert,
          uint32_t ProbeWords>
static inline uint64_t select_512_avx512_hybrid_impl(const uint64_t* x,
                                                     uint64_t rank) noexcept {
  for (uint32_t i = 0; i < ProbeWords; ++i) {
    const uint64_t count = select_512_word_count<Invert>(x[i]);
    if (rank < count) {
      const uint64_t word = select_512_select_word<Invert>(x[i]);
      return i * 64 + SelectWord(word, rank);
    }
    rank -= count;
  }
  return select_512_avx512_tail_impl<SelectWord, Invert>(x, rank, ProbeWords,
                                                         8 - ProbeWords);
}
#endif

static inline uint64_t select_512_avx512_pdep(const uint64_t* x,
                                              uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_pdep, false>(x, rank);
#else
  return select_512_avx2_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_pdep(const uint64_t* x,
                                               uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_pdep, true>(x, rank);
#else
  return select0_512_avx2_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_byte_lut(const uint64_t* x,
                                                  uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_byte_lut, false>(x, rank);
#else
  return select_512_avx2_byte_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_byte_lut(const uint64_t* x,
                                                   uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_byte_lut, true>(x, rank);
#else
  return select0_512_avx2_byte_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx512_broadword_lut(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_broadword_lut, false>(x, rank);
#else
  return select_512_avx2_broadword_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_broadword_lut(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_broadword_lut, true>(x, rank);
#else
  return select0_512_avx2_broadword_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx512_binary_lut(const uint64_t* x,
                                                    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_binary_lut, false>(x, rank);
#else
  return select_512_avx2_binary_lut(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_binary_lut(const uint64_t* x,
                                                     uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_binary_lut, true>(x, rank);
#else
  return select0_512_avx2_binary_lut(x, rank);
#endif
}

static inline uint64_t select_512_avx512_vigna_broadword(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_vigna_broadword, false>(x, rank);
#else
  return select_512_avx2_vigna_broadword(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_vigna_broadword(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_vigna_broadword, true>(x, rank);
#else
  return select0_512_avx2_vigna_broadword(x, rank);
#endif
}

static inline uint64_t select_512_avx512_no_pdep(const uint64_t* x,
                                                 uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_no_pdep, false>(x, rank);
#else
  return select_512_avx2_no_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_no_pdep(const uint64_t* x,
                                                  uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_no_pdep, true>(x, rank);
#else
  return select0_512_avx2_no_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_hybrid1_pdep(const uint64_t* x,
                                                      uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, false, 1>(x, rank);
#else
  return select_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_hybrid1_pdep(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, true, 1>(x, rank);
#else
  return select0_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_hybrid2_pdep(const uint64_t* x,
                                                      uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, false, 2>(x, rank);
#else
  return select_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_hybrid2_pdep(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, true, 2>(x, rank);
#else
  return select0_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_hybrid4_pdep(const uint64_t* x,
                                                      uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, false, 4>(x, rank);
#else
  return select_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_hybrid4_pdep(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_pdep, true, 4>(x, rank);
#else
  return select0_512_scalar_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_hybrid2_no_pdep(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_no_pdep, false, 2>(x, rank);
#else
  return select_512_scalar_no_pdep(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_hybrid2_no_pdep(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_hybrid_impl<select_64_no_pdep, true, 2>(x, rank);
#else
  return select0_512_scalar_no_pdep(x, rank);
#endif
}

static inline uint64_t select_512_avx512_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_experimental_default, false>(x, rank);
#else
  return select_512_avx2_experimental_default(x, rank);
#endif
}

static inline uint64_t select0_512_avx512_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_AVX512_SUPPORT
  return select_512_avx512_impl<select_64_experimental_default, true>(x, rank);
#else
  return select0_512_avx2_experimental_default(x, rank);
#endif
}

static inline uint64_t select_512_experimental_default(const uint64_t* x,
                                                       uint64_t rank) noexcept {
#ifdef PIXIE_EXPERIMENTAL_SELECT_NO_PDEP
  return select_512_scalar_no_pdep(x, rank);
#elif defined(PIXIE_AVX512_SUPPORT)
  return select_512_avx512_experimental_default(x, rank);
#else
  return select_512_scalar_experimental_default(x, rank);
#endif
}

static inline uint64_t select0_512_experimental_default(
    const uint64_t* x,
    uint64_t rank) noexcept {
#ifdef PIXIE_EXPERIMENTAL_SELECT_NO_PDEP
  return select0_512_scalar_no_pdep(x, rank);
#elif defined(PIXIE_AVX512_SUPPORT)
  return select0_512_avx512_experimental_default(x, rank);
#else
  return select0_512_scalar_experimental_default(x, rank);
#endif
}

}  // namespace pixie::experimental
