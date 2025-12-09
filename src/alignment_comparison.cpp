#include <benchmark/benchmark.h>
#include <immintrin.h>

#include <vector>

#include "bits.h"

#ifdef PIXIE_AVX512_SUPPORT

alignas(64) uint8_t data[128];

static void BM_Loadu512_shift63(benchmark::State& state) {
  const __m512i* ptr = reinterpret_cast<const __m512i*>(data + 63);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm512_loadu_si512(ptr));
  }
}

static void BM_Loadu512_shift31(benchmark::State& state) {
  const __m512i* ptr = reinterpret_cast<const __m512i*>(data + 31);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm512_loadu_si512(ptr));
  }
}

static void BM_Loadu512_shift0(benchmark::State& state) {
  const __m512i* ptr = reinterpret_cast<const __m512i*>(data);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm512_loadu_si512(ptr));
  }
}

static void BM_Load512_shift0(benchmark::State& state) {
  const __m512i* ptr = reinterpret_cast<const __m512i*>(data);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm512_load_si512(ptr));
  }
}

static void BM_Storeu512_shift63(benchmark::State& state) {
  __m512i value = _mm512_setzero_si512();
  __m512i* ptr = reinterpret_cast<__m512i*>(data + 63);

  for (auto _ : state) {
    _mm512_storeu_si512(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Storeu512_shift31(benchmark::State& state) {
  __m512i value = _mm512_setzero_si512();
  __m512i* ptr = reinterpret_cast<__m512i*>(data + 31);

  for (auto _ : state) {
    _mm512_storeu_si512(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Storeu512_shift0(benchmark::State& state) {
  __m512i value = _mm512_setzero_si512();
  __m512i* ptr = reinterpret_cast<__m512i*>(data);

  for (auto _ : state) {
    _mm512_storeu_si512(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Store512_shift0(benchmark::State& state) {
  __m512i value = _mm512_setzero_si512();
  __m512i* ptr = reinterpret_cast<__m512i*>(data);

  for (auto _ : state) {
    _mm512_store_si512(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

BENCHMARK(BM_Loadu512_shift63);
BENCHMARK(BM_Loadu512_shift31);
BENCHMARK(BM_Loadu512_shift0);
BENCHMARK(BM_Load512_shift0);
BENCHMARK(BM_Storeu512_shift63);
BENCHMARK(BM_Storeu512_shift31);
BENCHMARK(BM_Storeu512_shift0);
BENCHMARK(BM_Store512_shift0);

#else

alignas(64) uint8_t data[128];

static void BM_Loadu256_shift63(benchmark::State& state) {
  const __m256i* ptr = reinterpret_cast<const __m256i*>(data + 63);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Loadu256_shift31(benchmark::State& state) {
  const __m256i* ptr = reinterpret_cast<const __m256i*>(data + 31);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Loadu256_shift0(benchmark::State& state) {
  const __m256i* ptr = reinterpret_cast<const __m256i*>(data);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm256_loadu_si256(ptr));
  }
}

static void BM_Load256_shift0(benchmark::State& state) {
  const __m256i* ptr = reinterpret_cast<const __m256i*>(data);

  for (auto _ : state) {
    benchmark::DoNotOptimize(_mm256_load_si256(ptr));
  }
}

static void BM_Storeu256_shift63(benchmark::State& state) {
  __m256i value = _mm256_setzero_si256();
  __m256i* ptr = reinterpret_cast<__m256i*>(data + 63);

  for (auto _ : state) {
    _mm256_storeu_si256(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Storeu256_shift31(benchmark::State& state) {
  __m256i value = _mm256_setzero_si256();
  __m256i* ptr = reinterpret_cast<__m256i*>(data + 31);

  for (auto _ : state) {
    _mm256_storeu_si256(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Storeu256_shift0(benchmark::State& state) {
  __m256i value = _mm256_setzero_si256();
  __m256i* ptr = reinterpret_cast<__m256i*>(data);

  for (auto _ : state) {
    _mm256_storeu_si256(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

static void BM_Store256_shift0(benchmark::State& state) {
  __m256i value = _mm256_setzero_si256();
  __m256i* ptr = reinterpret_cast<__m256i*>(data);

  for (auto _ : state) {
    _mm256_store_si256(ptr, value);
    benchmark::DoNotOptimize(*ptr);
  }
}

BENCHMARK(BM_Loadu256_shift63);
BENCHMARK(BM_Loadu256_shift31);
BENCHMARK(BM_Loadu256_shift0);
BENCHMARK(BM_Load256_shift0);
BENCHMARK(BM_Storeu256_shift63);
BENCHMARK(BM_Storeu256_shift31);
BENCHMARK(BM_Storeu256_shift0);
BENCHMARK(BM_Store256_shift0);

#endif
