#include <benchmark/benchmark.h>
#include <immintrin.h>

#include <chrono>
#include <random>
#include <vector>

#include "bits.h"

alignas(64) uint8_t data[(1 << 29) + 1];

#ifdef PIXIE_AVX512_SUPPORT

static void BM_Loadu512_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    const __m512i* ptr = reinterpret_cast<const __m512i*>(data + idx);

    benchmark::DoNotOptimize(_mm512_loadu_si512(ptr));
  }
}

static void BM_Loadu512_unaligned_crossing_64byte_border(
    benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 48;
    const __m512i* ptr = reinterpret_cast<const __m512i*>(data + idx);

    benchmark::DoNotOptimize(_mm512_loadu_si512(ptr));
  }
}

static void BM_Load512_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    const __m512i* ptr = reinterpret_cast<const __m512i*>(data + idx);

    benchmark::DoNotOptimize(_mm512_load_si512(ptr));
  }
}

static void BM_Storeu512_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    __m512i* ptr = reinterpret_cast<__m512i*>(data + idx);
    __m512i value = _mm512_setzero_si512();

    _mm512_storeu_si512(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

static void BM_Storeu512_unaligned_crossing_64byte_border(
    benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 48;
    __m512i* ptr = reinterpret_cast<__m512i*>(data + idx);
    __m512i value = _mm512_setzero_si512();

    _mm512_storeu_si512(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

static void BM_Store512_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    __m512i* ptr = reinterpret_cast<__m512i*>(data + idx);
    __m512i value = _mm512_setzero_si512();

    _mm512_store_si512(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

BENCHMARK(BM_Loadu512_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Loadu512_unaligned_crossing_64byte_border)
    ->ArgNames({"k"})
    ->DenseRange(1, 23, 2);

BENCHMARK(BM_Load512_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Storeu512_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Storeu512_unaligned_crossing_64byte_border)
    ->ArgNames({"k"})
    ->DenseRange(1, 23, 2);

BENCHMARK(BM_Store512_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

#else

static void BM_Loadu256_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    const __m256i* ptr = reinterpret_cast<const __m256i*>(data + idx);

    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Loadu256_unaligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 16;
    const __m256i* ptr = reinterpret_cast<const __m256i*>(data + idx);

    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Loadu256_unaligned_crossing_64byte_border(
    benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 48;
    const __m256i* ptr = reinterpret_cast<const __m256i*>(data + idx);

    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Load256_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    const __m256i* ptr = reinterpret_cast<const __m256i*>(data + idx);

    benchmark::DoNotOptimize(_mm256_load_si256(ptr));
  }
}

static void BM_Storeu256_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    __m256i* ptr = reinterpret_cast<__m256i*>(data + idx);
    __m256i value = _mm256_setzero_si256();

    _mm256_storeu_si256(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

static void BM_Storeu256_unaligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 16;
    __m256i* ptr = reinterpret_cast<__m256i*>(data + idx);
    __m256i value = _mm256_setzero_si256();

    _mm256_storeu_si256(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

static void BM_Storeu256_unaligned_crossing_64byte_border(
    benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = ((rng() & ((1 << k) - 1)) << 6) + 48;
    __m256i* ptr = reinterpret_cast<__m256i*>(data + idx);
    __m256i value = _mm256_setzero_si256();

    _mm256_storeu_si256(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

static void BM_Store256_aligned(benchmark::State& state) {
  size_t k = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    size_t idx = (rng() & ((1 << k) - 1)) << 6;
    __m256i* ptr = reinterpret_cast<__m256i*>(data + idx);
    __m256i value = _mm256_setzero_si256();

    _mm256_store_si256(ptr, value);

    benchmark::DoNotOptimize(ptr);
  }
}

BENCHMARK(BM_Loadu256_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Loadu256_unaligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Loadu256_unaligned_crossing_64byte_border)
    ->ArgNames({"k"})
    ->DenseRange(1, 23, 2);

BENCHMARK(BM_Load256_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Storeu256_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Storeu256_unaligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

BENCHMARK(BM_Storeu256_unaligned_crossing_64byte_border)
    ->ArgNames({"k"})
    ->DenseRange(1, 23, 2);

BENCHMARK(BM_Store256_aligned)->ArgNames({"k"})->DenseRange(1, 23, 2);

#endif
