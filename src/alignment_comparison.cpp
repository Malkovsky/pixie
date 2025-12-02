#include <benchmark/benchmark.h>
#include <immintrin.h>
#include <vector>

#include "bits.h"

#ifdef PIXIE_AVX512_SUPPORT

static void BM_Loadu512(benchmark::State& state) {
    alignas(64) uint8_t data[128];

    const __m512i* ptr = reinterpret_cast<const __m512i*>(data);

    for (auto _ : state) {
        benchmark::DoNotOptimize(_mm256_loadu_si512(ptr));
    }
}


static void BM_Load512(benchmark::State& state) {
    alignas(64) uint8_t data[128];

    const __m512i* ptr = reinterpret_cast<const __m512i*>(data);

    for (auto _ : state) {
        benchmark::DoNotOptimize(_mm512_load_si512(ptr));
    }
}


static void BM_Storeu512(benchmark::State& state) {
    alignas(64) uint8_t data[128];
    __m512i value = _mm512_setzero_si512();
    __m512i* ptr = reinterpret_cast<__m512i*>(data);

    for (auto _ : state) {
        _mm512_storeu_si512(ptr, value);
        benchmark::DoNotOptimize(*ptr);
    }
}


static void BM_Store512(benchmark::State& state) {
    alignas(64) uint8_t data[128];
    __m512i value = _mm512_setzero_si512();
    __m512i* ptr = reinterpret_cast<__m512i*>(data);

    for (auto _ : state) {
        _mm512_store_si512(ptr, value);
        benchmark::DoNotOptimize(*ptr);
    }
}


BENCHMARK(BM_Loadu512);
BENCHMARK(BM_Load512);
BENCHMARK(BM_Storeu512);
BENCHMARK(BM_Store512);

#else

static void BM_Loadu256(benchmark::State& state) {
    alignas(32) uint8_t data[64];

    const __m256i* ptr = reinterpret_cast<const __m256i*>(data);

    for (auto _ : state) {
        benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
    }
}


static void BM_Load256(benchmark::State& state) {
    alignas(32) uint8_t data[64];

    const __m256i* ptr = reinterpret_cast<const __m256i*>(data);

    for (auto _ : state) {
        benchmark::DoNotOptimize(_mm256_load_si256(ptr));
    }
}


static void BM_Storeu256(benchmark::State& state) {
    alignas(32) uint8_t data[64];
    __m256i value = _mm256_setzero_si256();
    __m256i* ptr = reinterpret_cast<__m256i*>(data);

    for (auto _ : state) {
        _mm256_storeu_si256(ptr, value);
        benchmark::DoNotOptimize(*ptr);
    }
}


static void BM_Store256(benchmark::State& state) {
    alignas(32) uint8_t data[64];
    __m256i value = _mm256_setzero_si256();
    __m256i* ptr = reinterpret_cast<__m256i*>(data);

    for (auto _ : state) {
        _mm256_store_si256(ptr, value);
        benchmark::DoNotOptimize(*ptr);
    }
}


BENCHMARK(BM_Loadu256);
BENCHMARK(BM_Load256);
BENCHMARK(BM_Storeu256);
BENCHMARK(BM_Store256);

#endif