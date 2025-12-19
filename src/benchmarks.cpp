#include <benchmark/benchmark.h>
#include <bitvector.h>
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
// TODO: change the pasta/bit_vector usage of std::aligned_alloc
#include <pasta/bit_vector/bit_vector.hpp>
#endif
#include <cstdlib>
#include <random>
#include <vector>

static void BM_RankNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  for (auto& x : bits) {
    x = rng();
  }
  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  for (auto& x : bits) {
    x = rng();
  }
  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_RankInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  for (auto& x : bits) {
    x = rng();
  }
  pixie::BitVectorInterleaved bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_SelectNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  for (auto& x : bits) {
    x = rng();
  }
  pixie::BitVector bv(bits, n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  for (auto& x : bits) {
    x = rng();
  }
  pixie::BitVector bv(bits, n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved10PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.1;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved10PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.1;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved10PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.1;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved10PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.1;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  auto max_rank0 = bv.rank0(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank0 = rng() % max_rank0;
    benchmark::DoNotOptimize(bv.select0(rank0));
  }
}

static void BM_RankNonInterleaved90PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.9;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankZeroNonInterleaved90PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.9;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank0(pos));
  }
}

static void BM_SelectNonInterleaved90PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.9;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

  auto max_rank = bv.rank(bv.size()) + 1;

  for (auto _ : state) {
    uint64_t rank = rng() % max_rank;
    benchmark::DoNotOptimize(bv.select(rank));
  }
}

static void BM_SelectZeroNonInterleaved90PercentFill(benchmark::State& state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
  size_t num_ones = n * 0.9;
  for (int i = 0; i < num_ones; i++) {
    uint64_t pos = rng() % n;
    bits[pos / 64] |= (1ULL << pos % 64);
  }

  pixie::BitVector bv(bits, n);

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

BENCHMARK(BM_RankNonInterleaved10PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankZeroNonInterleaved10PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectNonInterleaved10PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectZeroNonInterleaved10PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankNonInterleaved90PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_RankZeroNonInterleaved90PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectNonInterleaved90PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);

BENCHMARK(BM_SelectZeroNonInterleaved90PercentFill)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34)
    ->Iterations(100000);