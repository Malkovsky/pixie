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

static void BM_WaveletTreeMmapSelect(benchmark::State& state) {
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

    WaveletTree orig_tree(alphabet_size, data);
    pixie::OutputBitStream bs;
    orig_tree.serialize(bs);
    std::vector<uint64_t> serialized_data = bs.extract();
    std::span<const std::byte> byte_span(
        reinterpret_cast<const std::byte*>(serialized_data.data()),
        serialized_data.size() * sizeof(uint64_t));

    state.ResumeTiming();

    auto mmap_tree =
        pixie::WaveletTreeBase<pixie::MmapViewStorage>::deserialize(byte_span);
    benchmark::DoNotOptimize(mmap_tree);

    for (size_t i = 0; i < query; i++) {
      size_t select = mmap_tree.select(query_symbol[i], query_pos[i]);
      benchmark::DoNotOptimize(select);
    }
  }
}

static void BM_WaveletTreeMmapRank(benchmark::State& state) {
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

    WaveletTree orig_tree(alphabet_size, data);
    pixie::OutputBitStream bs;
    orig_tree.serialize(bs);
    std::vector<uint64_t> serialized_data = bs.extract();
    std::span<const std::byte> byte_span(
        reinterpret_cast<const std::byte*>(serialized_data.data()),
        serialized_data.size() * sizeof(uint64_t));

    state.ResumeTiming();

    auto mmap_tree =
        pixie::WaveletTreeBase<pixie::MmapViewStorage>::deserialize(byte_span);
    benchmark::DoNotOptimize(mmap_tree);

    for (size_t i = 0; i < query; i++) {
      size_t rank = mmap_tree.rank(query_symbol[i], query_pos[i]);
      benchmark::DoNotOptimize(rank);
    }
  }
}

static void BM_WaveletTreeMmapSegment(benchmark::State& state) {
  size_t data_size = state.range(0), alphabet_size = 1024, query = data_size/128, length = 128;
  std::mt19937_64 rng(239);

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);
    std::vector<uint64_t> query_pos =
                              generate_random_data(query, data_size + 1 - length, rng);

    WaveletTree orig_tree(alphabet_size, data);
    pixie::OutputBitStream bs;
    orig_tree.serialize(bs);
    std::vector<uint64_t> serialized_data = bs.extract();
    std::span<const std::byte> byte_span(
        reinterpret_cast<const std::byte*>(serialized_data.data()),
        serialized_data.size() * sizeof(uint64_t));

    state.ResumeTiming();

    auto mmap_tree =
        pixie::WaveletTreeBase<pixie::MmapViewStorage>::deserialize(byte_span);
    benchmark::DoNotOptimize(mmap_tree);

    for (size_t i = 0; i < query; i++) {
      auto segment = mmap_tree.get_segment(query_pos[i], query_pos[i] + length);
      benchmark::DoNotOptimize(segment);
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

BENCHMARK(BM_WaveletTreeMmapSelect)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_WaveletTreeMmapSelect)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);

BENCHMARK(BM_WaveletTreeMmapRank)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_WaveletTreeMmapRank)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);

BENCHMARK(BM_WaveletTreeMmapSegment)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_WaveletTreeMmapSegment)
    ->ArgNames({"data_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);
