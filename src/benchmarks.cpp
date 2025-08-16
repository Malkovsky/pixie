#include <benchmark/benchmark.h>

#include <bitvector.h>
#include <random>
#include <vector>


static void BM_RankNonInterleaved(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(1 + n / 64);
  for (auto &x : bits) {
    x = rng();
  }
  pixie::BitVector bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

static void BM_RankInterleaved(benchmark::State &state) {
  size_t n = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<uint64_t> bits(1 + n / 64);
  for (auto &x : bits) {
    x = rng();
  }
  pixie::BitVectorInterleaved bv(bits, n);

  for (auto _ : state) {
    uint64_t pos = rng() % n;
    benchmark::DoNotOptimize(bv.rank(pos));
  }
}

BENCHMARK(BM_RankNonInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34);

BENCHMARK(BM_RankInterleaved)
    ->ArgNames({"n"})
    ->RangeMultiplier(4)
    ->Range(8, 1ull << 34);
