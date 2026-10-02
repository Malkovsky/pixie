#ifndef PIXIE_BENCHMARKS_PACKED_COPY_KERNELS_H_
#define PIXIE_BENCHMARKS_PACKED_COPY_KERNELS_H_

#include <pixie/bits.h>

// Unrolling experiment, 2026-09-16: Ryzen 7 8845HS, Linux/WSL, GCC 13.3,
// -O3 -march=native, CPU 0. Median CPU ns, five randomized 0.2s repetitions;
// confirmation run after the full aligned/shifted/mixed sweep. Warm buffers,
// runtime indirect calls, setup excluded. Layout 1: source bit 17, destination
// bit 0; layout 2: mixed offsets and lengths. Registered in
// bit_sequence_benchmarks: cpp-bench native bit_sequence_benchmarks
//   '^PackedCopy/(shift1|funnel[124])/bytes:[0-9]+/layout:[12]$' 5 0.2s
//
// | bytes | layout | shift1 ns | funnel1 ns | funnel2 ns | funnel4 ns |
// | ----: | -----: | --------: | ---------: | ---------: | ---------: |
// |   128 |      1 |       4.8 |        4.7 |        5.5 |        6.1 |
// |   256 |      1 |       6.8 |        7.9 |        8.8 |        7.9 |
// |   256 |      2 |      11.0 |       10.4 |       10.6 |       10.6 |
// |  4096 |      2 |      68.0 |       63.3 |       66.7 |       64.5 |
//
// No consistent unrolling gain across both runs; production keeps rolled
// funnel. Controls retained for Intel testing, not as public library
// implementations. nm/objdump: funnel 1x/2x/4x = 1011/1085/1282 bytes, 1/2/4
// VPSHRDVQ per main iteration (7/13/21 total instructions). GCC also unrolls
// the bounded remainder; the rolled main loop is NOT automatically unrolled. No
// hardware counters taken.

namespace packed_copy_benchmark {

// Benchmark-only controls for Intel/AMD comparisons. Same contract and boundary
// work as copy_packed_bits; only the shifted AVX-512 loop varies. Keep out of
// the public API. Noinline gives every variant the same runtime-call boundary
// and leaves inspectable symbols for checking actual unrolling and code size.
#if defined(__AVX512F__)
template <bool Funnel, size_t Unroll>
__attribute__((noinline)) void Copy(const uint64_t* source,
                                    size_t source_bit,
                                    uint64_t* destination,
                                    size_t destination_bit,
                                    size_t count) {
  static_assert(Unroll == 1 || Unroll == 2 || Unroll == 4);
  const auto boundary = [&](size_t width) {
    const auto shift = source_bit % 64;
    uint64_t value = source[source_bit / 64] >> shift;
    if (width > 64 - shift) {
      value |= source[source_bit / 64 + 1] << (64 - shift);
    }
    const auto offset = destination_bit % 64;
    const auto mask = first_bits_mask(width) << offset;
    auto& word = destination[destination_bit / 64];
    word = (word & ~mask) | ((value << offset) & mask);
    source_bit += width;
    destination_bit += width;
    count -= width;
  };
  if (count == 0) {
    return;
  }
  if (destination_bit % 64 != 0) {
    boundary(std::min(count, 64 - destination_bit % 64));
  }
  const auto shift = source_bit % 64;
  const auto* input = source + source_bit / 64;
  auto* output = destination + destination_bit / 64;
  auto words = count / 64;
  if (shift == 0) {
    std::copy_n(input, words, output);
  } else {
    const auto low_shift = _mm_cvtsi64_si128(shift);
    const auto high_shift = _mm_cvtsi64_si128(64 - shift);
#if defined(__AVX512VBMI2__)
    const auto shifts = _mm512_set1_epi64(shift);
#else
    static_assert(!Funnel, "Funnel control requires AVX512VBMI2");
#endif
    const auto vector = [&](size_t offset) {
      const auto low = _mm512_loadu_si512(input + offset);
      const auto high = _mm512_loadu_si512(input + offset + 1);
#if defined(__AVX512VBMI2__)
      if constexpr (Funnel) {
        _mm512_storeu_si512(output + offset,
                            _mm512_shrdv_epi64(low, high, shifts));
      } else
#endif
      {
        _mm512_storeu_si512(
            output + offset,
            _mm512_or_si512(_mm512_srl_epi64(low, low_shift),
                            _mm512_sll_epi64(high, high_shift)));
      }
    };
    for (; words >= 8 * Unroll;
         words -= 8 * Unroll, input += 8 * Unroll, output += 8 * Unroll) {
      vector(0);
      if constexpr (Unroll >= 2) {
        vector(8);
      }
      if constexpr (Unroll >= 4) {
        vector(16);
        vector(24);
      }
    }
    if constexpr (Unroll != 1) {
      for (; words >= 8; words -= 8, input += 8, output += 8) {
        vector(0);
      }
    }
    for (size_t i = 0; i < words; ++i) {
      output[i] = (input[i] >> shift) | (input[i + 1] << (64 - shift));
    }
  }
  source_bit += count / 64 * 64;
  destination_bit += count / 64 * 64;
  count %= 64;
  if (count != 0) {
    boundary(count);
  }
}
#endif

__attribute__((noinline)) inline void Production(const uint64_t* source,
                                                 size_t source_bit,
                                                 uint64_t* destination,
                                                 size_t destination_bit,
                                                 size_t count) {
  pixie::copy_packed_bits(source, source_bit, destination, destination_bit,
                          count);
}

}  // namespace packed_copy_benchmark

#endif  // PIXIE_BENCHMARKS_PACKED_COPY_KERNELS_H_
