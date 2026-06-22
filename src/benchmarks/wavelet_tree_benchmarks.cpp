#include <benchmark/benchmark.h>
#include <pixie/utils.h>
#include <pixie/wavelet_tree.h>

#include <random>

using pixie::WaveletTree;

static void BM_WaveletTreeSelect(benchmark::State& state) {
  size_t data_size = state.range(0), alphabet_size = 1024, query = data_size;
  std::mt19937_64 rng(239);

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);
    std::vector<uint64_t> query_symbol =
        generate_random_data(query, alphabet_size, rng);
    std::vector<size_t> count(alphabet_size), query_pos(query);
    for (auto symb : data) {
      count[symb]++;
    }
    for (size_t i = 0; i < query; i++) {
      query_pos[i] = 1 + std::uniform_int_distribution<size_t>(
                             0, count[query_symbol[i]])(rng);
    }

    state.ResumeTiming();

    WaveletTree wavelet_tree(alphabet_size, data);
    benchmark::DoNotOptimize(wavelet_tree);

    for (size_t i = 0; i < query; i++) {
      size_t select = wavelet_tree.select(query_symbol[i], query_pos[i]);
      benchmark::DoNotOptimize(select);
    }
  }
}

static void BM_WaveletTreeRank(benchmark::State& state) {
  size_t data_size = state.range(0), alphabet_size = 1024, query = data_size;
  std::mt19937_64 rng(239);

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);
    std::vector<uint64_t> query_symbol =
                              generate_random_data(query, alphabet_size, rng),
                          query_pos =
                              generate_random_data(query, data_size + 1, rng);

    state.ResumeTiming();

    WaveletTree wavelet_tree(alphabet_size, data);
    benchmark::DoNotOptimize(wavelet_tree);

    for (size_t i = 0; i < query; i++) {
      size_t rank = wavelet_tree.rank(query_symbol[i], query_pos[i]);
      benchmark::DoNotOptimize(rank);
    }
  }
}

BENCHMARK(BM_WaveletTreeSelect)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_WaveletTreeSelect)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);

BENCHMARK(BM_WaveletTreeRank)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_WaveletTreeRank)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);
