#pragma once

#include <pixie/bits.h>

#if defined(PIXIE_AVX512_SUPPORT) || defined(PIXIE_AVX2_SUPPORT)
#include <immintrin.h>
#endif

#include <array>
#include <cstdint>
#include <cstring>
#include <span>

namespace pixie::detail::tree_unpack {
// Read at most 64 bits, treating bits beyond the source as zero. No padding
// words are required. Callers separately bound consumption by the logical size.
inline uint64_t word_at(std::span<const uint64_t> words, size_t p) {
  const size_t index = p / 64;
  if (index >= words.size()) {
    return 0;
  }
  const unsigned shift = p % 64;
  uint64_t result = words[index] >> shift;
  if (shift && index + 1 < words.size()) {
    result |= words[index + 1] << (64 - shift);
  }
  return result;
}
// Fill pending-child parent slots. Small variable degrees are common; a
// masked store avoids scalar loop exits and vector-tail dispatch per node.
// The caller provides exactly count writable uint16_t elements (zero is valid).
inline void fill_parent_slots(uint16_t* destination,
                              unsigned count,
                              uint16_t parent) {
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512BW__)
  const __m512i values = _mm512_set1_epi16(parent);
  while (count >= 32) {
    _mm512_storeu_si512(destination, values);
    destination += 32;
    count -= 32;
  }
  const __mmask32 mask = (uint32_t{1} << count) - 1;
  _mm512_mask_storeu_epi16(destination, mask, values);
#elif defined(PIXIE_AVX2_SUPPORT)
  const auto values = _mm256_set1_epi16(parent);
  while (count >= 16) {
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), values);
    destination += 16;
    count -= 16;
  }
  // AVX2 has no masked 16-bit store. Overlapping head/tail stores stay
  // within the requested range, including for unaligned destinations.
  const auto half = _mm256_castsi256_si128(values);
  if (count >= 8) {
    _mm_storeu_si128(reinterpret_cast<__m128i*>(destination), half);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(destination + count - 8), half);
  } else if (count >= 4) {
    _mm_storel_epi64(reinterpret_cast<__m128i*>(destination), half);
    _mm_storel_epi64(reinterpret_cast<__m128i*>(destination + count - 4), half);
  } else if (count >= 2) {
    const uint32_t pair = uint32_t(parent) * 0x00010001u;
    std::memcpy(destination, &pair, sizeof(pair));
    std::memcpy(destination + count - 2, &pair, sizeof(pair));
  } else if (count) {
    *destination = parent;
  }
#else
  for (unsigned i = 0; i < count; ++i) {
    destination[i] = parent;
  }
#endif
}

// Uniform lookup over every possible byte and its following bit.
// Entries describe events and surviving local opens; no tree-dependent table.
struct alignas(32) ByteProgram {
  uint8_t event[8], rank[8], update[8];
  uint8_t count, opens, kept;
  int8_t delta;
};
inline constexpr auto byte_programs = [] {
  std::array<ByteProgram, 512> result{};
  for (unsigned pattern = 0; pattern < 512; ++pattern) {
    auto& p = result[pattern];
    int local[17];
    for (int& v : local) {
      v = -1;
    }
    int depth = 0;
    unsigned rank = 0;
    for (int i = 0; i < 8; ++i) {
      p.rank[i] = rank;
      if ((pattern >> i) & 1) {
        local[depth + 8] = i;
        ++depth;
        p.event[p.count++] = 8 + i;
        ++rank;
      } else {
        --depth;
        local[depth + 8] = -1;
        if (!((pattern >> (i + 1)) & 1)) {
          int slot = depth - 1;
          int opening = slot >= -8 ? local[slot + 8] : -1;
          // -9 occurs only for a zero byte; that byte takes the scalar path.
          p.event[p.count++] = opening >= 0 ? 8 + opening : slot + 8;
        }
      }
    }
    p.opens = rank;
    p.delta = depth;
    for (int slot = -8; slot < 8; ++slot) {
      if (local[slot + 8] >= 0) {
        p.update[p.kept++] = local[slot + 8];
      }
    }
  }
  return result;
}();

inline void decode_byte(const ByteProgram& p,
                        size_t base,
                        size_t& number,
                        size_t& depth,
                        std::array<uint64_t, 136>& numbers,
                        std::array<uint64_t, 136>& positions,
                        uint64_t* event_numbers,
                        uint64_t* event_positions) {
#if defined(PIXIE_AVX512_SUPPORT)
  auto local_num = _mm512_add_epi64(
      _mm512_set1_epi64(number),
      _mm512_cvtepu8_epi64(
          _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.rank))));
  auto local_pos = _mm512_add_epi64(_mm512_set1_epi64(base),
                                    _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7));
  auto indices = _mm512_cvtepu8_epi64(
      _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.event)));
  auto num = _mm512_permutex2var_epi64(
      _mm512_loadu_si512(numbers.data() + depth), indices, local_num);
  auto pos = _mm512_permutex2var_epi64(
      _mm512_loadu_si512(positions.data() + depth), indices, local_pos);
  _mm512_store_si512(event_numbers, num);
  _mm512_store_si512(event_positions, pos);
  depth += p.delta;
  number += p.opens;
  auto keep = _mm512_cvtepu8_epi64(
      _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.update)));
  auto mask = __mmask8((1u << p.kept) - 1);
  _mm512_mask_storeu_epi64(numbers.data() + 8 + depth - p.kept, mask,
                           _mm512_permutexvar_epi64(keep, local_num));
  _mm512_mask_storeu_epi64(positions.data() + 8 + depth - p.kept, mask,
                           _mm512_permutexvar_epi64(keep, local_pos));
#elif defined(PIXIE_AVX2_SUPPORT)
  // AVX2 permutes 32-bit lanes. Duplicate each event index to move both
  // halves of a 64-bit coordinate together, then select local or saved nodes.
  const auto ranks = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.rank));
  const auto nbase = _mm256_set1_epi64x(number);
  const auto pbase = _mm256_set1_epi64x(base);
  const auto events =
      _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.event));
  const auto num_low = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(numbers.data() + depth));
  const auto num_high = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(numbers.data() + depth + 4));
  const auto pos_low = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(positions.data() + depth));
  const auto pos_high = _mm256_loadu_si256(
      reinterpret_cast<const __m256i*>(positions.data() + depth + 4));
  auto emit = [&](const __m128i codes, unsigned offset) {
    const auto local = _mm_and_si128(codes, _mm_set1_epi8(7));
    const auto indices = _mm_cvtepu8_epi32(codes);
    const auto duplicated =
        _mm256_set_m128i(_mm_unpackhi_epi32(indices, indices),
                         _mm_unpacklo_epi32(indices, indices));
    const auto lanes =
        _mm256_or_si256(_mm256_slli_epi32(duplicated, 1),
                        _mm256_setr_epi32(0, 1, 0, 1, 0, 1, 0, 1));
    const auto upper = _mm256_srai_epi32(_mm256_slli_epi32(duplicated, 29), 31);
    const auto ancestor = _mm256_cmpgt_epi32(_mm256_set1_epi32(8), duplicated);
    const auto local_num = _mm256_add_epi64(
        nbase, _mm256_cvtepu8_epi64(_mm_shuffle_epi8(ranks, local)));
    const auto local_pos = _mm256_add_epi64(pbase, _mm256_cvtepu8_epi64(local));
    const auto previous_num =
        _mm256_blendv_epi8(_mm256_permutevar8x32_epi32(num_low, lanes),
                           _mm256_permutevar8x32_epi32(num_high, lanes), upper);
    const auto previous_pos =
        _mm256_blendv_epi8(_mm256_permutevar8x32_epi32(pos_low, lanes),
                           _mm256_permutevar8x32_epi32(pos_high, lanes), upper);
    const auto num = _mm256_blendv_epi8(local_num, previous_num, ancestor);
    const auto pos = _mm256_blendv_epi8(local_pos, previous_pos, ancestor);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(event_numbers + offset),
                        num);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(event_positions + offset),
                        pos);
  };
  emit(events, 0);
  emit(_mm_srli_si128(events, 4), 4);
  depth += p.delta;
  const auto updates =
      _mm_loadl_epi64(reinterpret_cast<const __m128i*>(p.update));
  auto keep = [&](const __m128i codes, unsigned offset) {
    const auto mask =
        _mm256_cmpgt_epi64(_mm256_set1_epi64x(p.kept),
                           _mm256_add_epi64(_mm256_set1_epi64x(offset),
                                            _mm256_setr_epi64x(0, 1, 2, 3)));
    const auto num = _mm256_add_epi64(
        nbase, _mm256_cvtepu8_epi64(_mm_shuffle_epi8(ranks, codes)));
    const auto pos = _mm256_add_epi64(pbase, _mm256_cvtepu8_epi64(codes));
    const auto start = 8 + depth - p.kept + offset;
    _mm256_maskstore_epi64(reinterpret_cast<long long*>(numbers.data() + start),
                           mask, num);
    _mm256_maskstore_epi64(
        reinterpret_cast<long long*>(positions.data() + start), mask, pos);
  };
  if (p.kept) {
    keep(updates, 0);
    if (p.kept > 4) {
      keep(_mm_srli_si128(updates, 4), 4);
    }
  }
  number += p.opens;
#else
  for (unsigned i = 0; i < p.count; ++i) {
    unsigned s = p.event[i];
    event_numbers[i] = s < 8 ? numbers[depth + s] : number + p.rank[s - 8];
    event_positions[i] = s < 8 ? positions[depth + s] : base + s - 8;
  }
  depth += p.delta;
  for (unsigned i = 0; i < p.kept; ++i) {
    unsigned b = p.update[i];
    numbers[depth + 8 - p.kept + i] = number + p.rank[b];
    positions[depth + 8 - p.kept + i] = base + b;
  }
  number += p.opens;
#endif
}
// Expand relative event codes into interleaved 64-bit node-number/position
// pairs. Codes encode a 3-bit node offset and a 4-bit position at bits 3..6;
// bit 7 belongs to the caller. Input and output are padded to eight events.
// Output may be unaligned. Only the first count pairs are logically valid.
inline void decode_relative_events(const uint8_t* codes,
                                   size_t count,
                                   uint64_t number,
                                   uint64_t position,
                                   uint64_t* output) {
#if defined(PIXIE_AVX512_SUPPORT)
  const auto nbase = _mm512_set1_epi64(number),
             pbase = _mm512_set1_epi64(position);
  const auto low = _mm512_setr_epi64(0, 8, 1, 9, 2, 10, 3, 11);
  const auto high = _mm512_setr_epi64(4, 12, 5, 13, 6, 14, 7, 15);
  for (size_t i = 0; i < count; i += 8) {
    const auto packed = _mm512_cvtepu8_epi64(
        _mm_loadl_epi64(reinterpret_cast<const __m128i*>(codes + i)));
    const auto nums =
        _mm512_add_epi64(nbase, _mm512_and_si512(packed, _mm512_set1_epi64(7)));
    const auto positions = _mm512_add_epi64(
        pbase,
        _mm512_and_si512(_mm512_srli_epi64(packed, 3), _mm512_set1_epi64(15)));
    _mm512_storeu_si512(output + 2 * i,
                        _mm512_permutex2var_epi64(nums, low, positions));
    _mm512_storeu_si512(output + 2 * i + 8,
                        _mm512_permutex2var_epi64(nums, high, positions));
  }
#elif defined(PIXIE_AVX2_SUPPORT)
  const auto nbase = _mm256_set1_epi64x(number),
             pbase = _mm256_set1_epi64x(position);
  for (size_t i = 0; i < count; i += 4) {
    uint32_t bytes;
    std::memcpy(&bytes, codes + i, 4);
    const auto packed = _mm256_cvtepu8_epi64(_mm_cvtsi32_si128(bytes));
    const auto nums = _mm256_add_epi64(
        nbase, _mm256_and_si256(packed, _mm256_set1_epi64x(7)));
    const auto positions = _mm256_add_epi64(
        pbase,
        _mm256_and_si256(_mm256_srli_epi64(packed, 3), _mm256_set1_epi64x(15)));
    const auto low = _mm256_unpacklo_epi64(nums, positions),
               high = _mm256_unpackhi_epi64(nums, positions);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + 2 * i),
                        _mm256_permute2x128_si256(low, high, 0x20));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(output + 2 * i + 4),
                        _mm256_permute2x128_si256(low, high, 0x31));
  }
#else
  for (size_t i = 0; i < count; ++i) {
    output[2 * i] = number + (codes[i] & 7);
    output[2 * i + 1] = position + ((codes[i] >> 3) & 15);
  }
#endif
}

#if defined(PIXIE_AVX2_SUPPORT) && \
    !(defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__))
// Keep the packing loop out of recursive traversal instantiations to limit
// code growth and preserve the compiler's recursive DFS unrolling budget.
[[gnu::noinline]] inline unsigned zero_offsets_avx2(uint64_t mask,
                                                    unsigned valid,
                                                    uint8_t* out) {
  // Each table word contains the selected bit offsets of one byte. Four
  // independent byte lookups share one gather and one offset addition.
  static constexpr auto offsets = [] {
    std::array<uint64_t, 256> result{};
    for (unsigned bits = 0; bits < result.size(); ++bits) {
      unsigned count = 0;
      for (unsigned bit = 0; bit < 8; ++bit) {
        if ((bits >> bit) & 1) {
          result[bits] |= uint64_t(bit) << (8 * count++);
        }
      }
    }
    return result;
  }();
  unsigned count = 0;
  for (unsigned base = 0; base < valid; base += 32) {
    const uint32_t bits = uint32_t(mask >> base);
    const auto indices = _mm_cvtepu8_epi32(_mm_cvtsi32_si128(bits));
    const auto local = _mm256_i32gather_epi64(
        reinterpret_cast<const long long*>(offsets.data()), indices, 8);
    const auto packed = _mm256_add_epi8(
        local, _mm256_add_epi8(_mm256_set1_epi8(base),
                               _mm256_setr_epi64x(0, 0x0808080808080808LL,
                                                  0x1010101010101010LL,
                                                  0x1818181818181818LL)));
    auto append = [&]<unsigned I>() {
      const uint64_t positions = _mm256_extract_epi64(packed, I);
      std::memcpy(out + count, &positions, sizeof(positions));
      count += std::popcount((bits >> (8 * I)) & 255u);
    };
    append.template operator()<0>();
    append.template operator()<1>();
    append.template operator()<2>();
    append.template operator()<3>();
  }
  return count;
}
#endif

// Extract zero positions below valid (0..64). The output has 64 writable
// bytes; only the returned number of entries is meaningful.
inline unsigned zero_offsets(uint64_t word, unsigned valid, uint8_t* out) {
  uint64_t mask = ~word;
  if (valid < 64) {
    mask &= (uint64_t(1) << valid) - 1;
  }
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
  alignas(64) static constexpr auto indices = [] {
    std::array<uint8_t, 64> x{};
    for (unsigned i = 0; i < 64; ++i) {
      x[i] = i;
    }
    return x;
  }();
  auto packed =
      _mm512_maskz_compress_epi8(mask, _mm512_load_si512(indices.data()));
  _mm512_storeu_si512(out, packed);
  const unsigned count = std::popcount(mask);
  return count;
#elif defined(PIXIE_AVX2_SUPPORT)
  return zero_offsets_avx2(mask, valid, out);
#else
  unsigned count = 0;
  while (mask) {
    out[count++] = std::countr_zero(mask);
    mask &= mask - 1;
  }
  return count;
#endif
}

}  // namespace pixie::detail::tree_unpack
