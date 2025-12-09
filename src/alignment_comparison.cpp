#include <benchmark/benchmark.h>
#include <immintrin.h>

#include <chrono>
#include <random>
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

static void BM_Loadu512_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 1 + 64 * (rng() % (n - 1));
    const __m512i* ptr = reinterpret_cast<const __m512i*>(
        reinterpret_cast<uint8_t*>(data) + idx);

    auto start = std::chrono::high_resolution_clock::now();

    _mm512_loadu_si512(ptr);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    benchmark::DoNotOptimize(ptr);

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Load512_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 64 * (rng() % (n - 1));
    const __m512i* ptr = reinterpret_cast<const __m512i*>(
        reinterpret_cast<uint8_t*>(data) + idx);

    auto start = std::chrono::high_resolution_clock::now();

    _mm512_load_si512(ptr);

    auto end = std::chrono::high_resolution_clock::now();

    benchmark::DoNotOptimize(ptr);

    std::chrono::duration<double> duration = end - start;

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Storeu512_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 1 + 64 * (rng() % (n - 1));
    __m512i* ptr =
        reinterpret_cast<__m512i*>(reinterpret_cast<uint8_t*>(data) + idx);
    __m512i value = _mm512_setzero_si512();

    auto start = std::chrono::high_resolution_clock::now();

    _mm512_storeu_si512(ptr, value);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    benchmark::DoNotOptimize(ptr);

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Store512_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 64 * (rng() % (n - 1));
    __m512i* ptr =
        reinterpret_cast<__m512i*>(reinterpret_cast<uint8_t*>(data) + idx);
    __m512i value = _mm512_setzero_si512();

    auto start = std::chrono::high_resolution_clock::now();

    _mm512_store_si512(ptr, value);

    auto end = std::chrono::high_resolution_clock::now();

    benchmark::DoNotOptimize(ptr);

    std::chrono::duration<double> duration = end - start;

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

BENCHMARK(BM_Loadu512_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Load512_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Storeu512_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Store512_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

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

static void BM_Loadu256_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 1 + 32 * (rng() % (n - 1));
    const __m256i* ptr = reinterpret_cast<const __m256i*>(
        reinterpret_cast<uint8_t*>(data) + idx);

    auto start = std::chrono::high_resolution_clock::now();

    _mm256_loadu_si256(ptr);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    benchmark::DoNotOptimize(ptr);

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Load256_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 32 * (rng() % (n - 1));
    const __m256i* ptr = reinterpret_cast<const __m256i*>(
        reinterpret_cast<uint8_t*>(data) + idx);

    auto start = std::chrono::high_resolution_clock::now();

    _mm256_load_si256(ptr);

    auto end = std::chrono::high_resolution_clock::now();

    benchmark::DoNotOptimize(ptr);

    std::chrono::duration<double> duration = end - start;

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Storeu256_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 1 + 32 * (rng() % (n - 1));
    __m256i* ptr =
        reinterpret_cast<__m256i*>(reinterpret_cast<uint8_t*>(data) + idx);
    __m256i value = _mm256_setzero_si256();

    auto start = std::chrono::high_resolution_clock::now();

    _mm256_storeu_si256(ptr, value);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> duration = end - start;

    benchmark::DoNotOptimize(ptr);

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

static void BM_Store256_Random(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  size_t alignment = 64;
  size_t size = 64 * n;

  void* data = std::aligned_alloc(alignment, size);

  for (auto _ : state) {
    size_t idx = 32 * (rng() % (n - 1));
    __m256i* ptr =
        reinterpret_cast<__m256i*>(reinterpret_cast<uint8_t*>(data) + idx);
    __m256i value = _mm256_setzero_si256();

    auto start = std::chrono::high_resolution_clock::now();

    _mm256_store_si256(ptr, value);

    auto end = std::chrono::high_resolution_clock::now();

    benchmark::DoNotOptimize(ptr);

    std::chrono::duration<double> duration = end - start;

    state.SetIterationTime(duration.count());
  }

  std::free(data);
}

BENCHMARK(BM_Loadu256_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Load256_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Storeu256_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

BENCHMARK(BM_Store256_Random)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(2, 1 << 20)
    ->UseManualTime();

#endif
