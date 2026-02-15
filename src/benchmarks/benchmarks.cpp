#include <benchmark/benchmark.h>
#include <pixie/bitvector.h>
#include <pixie/cache_line.h>

#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
// TODO: change the pasta/bit_vector usage of std::aligned_alloc
#include <pasta/bit_vector/bit_vector.hpp>
#endif
#include <bit>
#include <cstdlib>
#include <random>
#include <vector>

constexpr size_t kBenchmarkRandomCopies = 8;

#ifdef _WIN32
#include <windows.h>
void platform_setup() {
  // Pin to core 3
  SetThreadAffinityMask(GetCurrentThread(), DWORD_PTR(1) << 3);
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
  SetProcessPriorityBoost(GetCurrentProcess(), TRUE);
}
#else
#include <sched.h>
void platform_setup() {
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(3, &cpuset);
  sched_setaffinity(0, sizeof(cpuset), &cpuset);
}
#endif

static void PrepareRandomBits50pFill(std::span<uint64_t> bits) {
  auto randomly_filled_length = bits.size() / kBenchmarkRandomCopies;
  std::mt19937_64 rng(42);

  for (int i = 0; i < randomly_filled_length; ++i) {
    bits[i] = rng();
  }
  for (int i = 1; i < kBenchmarkRandomCopies; ++i) {
    std::copy_n(bits.begin(), randomly_filled_length,
                bits.begin() + i * randomly_filled_length);
  }
  for (int i = kBenchmarkRandomCopies * randomly_filled_length; i < bits.size();
       ++i) {
    bits[i] = rng();
  }
}

static void PrepareRandomBits12p5Fill(std::span<uint64_t> bits) {
  auto randomly_filled_length = bits.size() / kBenchmarkRandomCopies;
  std::mt19937_64 rng(42);

  for (int i = 0; i < randomly_filled_length; ++i) {
    bits[i] = rng() & rng() & rng();
  }
  for (int i = 1; i < kBenchmarkRandomCopies; ++i) {
    std::copy_n(bits.begin(), randomly_filled_length,
                bits.begin() + i * randomly_filled_length);
  }
  for (int i = kBenchmarkRandomCopies * randomly_filled_length; i < bits.size();
       ++i) {
    bits[i] = rng() & rng() & rng();
  }
}

static void PrepareRandomBits87p5Fill(std::span<uint64_t> bits) {
  auto randomly_filled_length = bits.size() / kBenchmarkRandomCopies;
  std::mt19937_64 rng(42);

  for (int i = 0; i < randomly_filled_length; ++i) {
    bits[i] = rng() | rng() | rng();
  }
  for (int i = 1; i < kBenchmarkRandomCopies; ++i) {
    std::copy_n(bits.begin(), randomly_filled_length,
                bits.begin() + i * randomly_filled_length);
  }
  for (int i = kBenchmarkRandomCopies * randomly_filled_length; i < bits.size();
       ++i) {
    bits[i] = rng() | rng() | rng();
  }
}

static void BM_RankNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  AlignedStorage bits(n);
  auto bits_as_words = bits.As64BitInts();
  PrepareRandomBits50pFill(bits_as_words);
  pixie::BitVector bv(bits_as_words, n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits50pFill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_RankInterleaved(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits50pFill(bits.As64BitInts());
  pixie::BitVectorInterleaved bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_SelectNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits50pFill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits50pFill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved12p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits12p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  if (std::popcount(n) == 1) {
    for (auto _ : state) {
      uint64_t pos = rng() & (n - 1);
      benchmark::DoNotOptimize(bv.rank(pos));
    }
  } else {
    for (auto _ : state) {
      uint64_t pos = rng() % n;
      benchmark::DoNotOptimize(bv.rank(pos));
    }
  }
}

static void BM_RankZeroNonInterleaved12p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits12p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved12p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits12p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved12p5PercentFill(
    benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits12p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits87p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits87p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits87p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved87p5PercentFill(
    benchmark::State& state) {
  size_t n = state.range(0);

  AlignedStorage bits(n);
  PrepareRandomBits87p5Fill(bits.As64BitInts());
  pixie::BitVector bv(bits.As64BitInts(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  std::mt19937_64 rng(42);
  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select(rank0));
  }
}

BENCHMARK(BM_RankInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankZeroNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectZeroNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankZeroNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectZeroNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_RankZeroNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);

BENCHMARK(BM_SelectZeroNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000)
    ->Repetitions(100);
