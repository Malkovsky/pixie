#include <benchmark/benchmark.h>
#include <pixie/bitvector.h>
#include <pixie/cache_line.h>

#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
// TODO: change the pasta/bit_vector usage of std::aligned_alloc
#include <pasta/bit_vector/bit_vector.hpp>
#endif
#include <cstdlib>
#include <random>
#include <vector>

static void BM_RankNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  AlignedStorage bits(n);

  std::mt19937_64 rng(42);
  for (auto& x : bits.AsWords()) {
    x = rng();
  }
  pixie::BitVector bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng();
  }
  pixie::BitVector bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_RankInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng();
  }
  pixie::BitVectorInterleaved bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_SelectNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng();
  }
  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng();
  }
  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved12p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() & rng() & rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

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
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() & rng() & rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved12p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() & rng() & rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved12p5PercentFill(
    benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() & rng() & rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() | rng() | rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() | rng() | rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved87p5PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() | rng() | rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved87p5PercentFill(
    benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  AlignedStorage bits(n);
  for (auto& x : bits.AsWords()) {
    x = rng() | rng() | rng();
  }

  pixie::BitVector bv(bits.AsWords(), n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select(rank0));
  }
}

BENCHMARK(BM_RankInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankZeroNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectZeroNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankZeroNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectZeroNonInterleaved12p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankZeroNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectZeroNonInterleaved87p5PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);
