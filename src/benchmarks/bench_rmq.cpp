#include <benchmark/benchmark.h>
#include <pixie/rmq.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kQueryCount = 32768;
using Index = std::size_t;

struct Dataset {
  std::size_t size = 0;
  std::size_t max_width = 0;
  std::vector<std::int64_t> values;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
};

struct DepthDataset {
  std::size_t size = 0;
  std::size_t max_width = 0;
  std::vector<std::int64_t> depths;
  std::vector<std::uint64_t> bits;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
};

Dataset make_dataset(std::size_t size, std::size_t max_width) {
  Dataset dataset;
  dataset.size = size;
  dataset.max_width = max_width;
  dataset.values.resize(size);
  dataset.ranges.resize(kQueryCount);

  std::mt19937_64 rng(kSeed ^ (size * 0x9E3779B185EBCA87ull) ^
                      (max_width * 0xBF58476D1CE4E5B9ull));
  std::uniform_int_distribution<std::int64_t> value_dist(-1'000'000, 1'000'000);
  std::generate(dataset.values.begin(), dataset.values.end(),
                [&] { return value_dist(rng); });

  std::uniform_int_distribution<std::size_t> left_dist(0, size - 1);
  for (auto& [left, right] : dataset.ranges) {
    left = left_dist(rng);
    const std::size_t available = size - left;
    const std::size_t width_limit = std::min(max_width, available);
    std::uniform_int_distribution<std::size_t> width_dist(1, width_limit);
    right = left + width_dist(rng) - 1;
  }
  return dataset;
}

DepthDataset make_depth_dataset(std::size_t size, std::size_t max_width) {
  DepthDataset dataset;
  dataset.size = size;
  dataset.max_width = max_width;
  dataset.depths.resize(size);
  dataset.ranges.resize(kQueryCount);

  std::mt19937_64 rng(kSeed ^ (size * 0xD6E8FEB86659FD93ull) ^
                      (max_width * 0xA5A3564E27F88695ull));
  for (std::size_t i = 1; i < dataset.depths.size(); ++i) {
    dataset.depths[i] = dataset.depths[i - 1] + ((rng() & 1u) ? 1 : -1);
  }
  dataset.bits.assign((size - 1 + 63) / 64, 0);
  for (std::size_t i = 1; i < dataset.depths.size(); ++i) {
    if (dataset.depths[i] - dataset.depths[i - 1] == 1) {
      dataset.bits[(i - 1) >> 6] |= std::uint64_t{1} << ((i - 1) & 63);
    }
  }

  std::uniform_int_distribution<std::size_t> left_dist(0, size - 1);
  for (auto& [left, right] : dataset.ranges) {
    left = left_dist(rng);
    const std::size_t available = size - left;
    const std::size_t width_limit = std::min(max_width, available);
    std::uniform_int_distribution<std::size_t> width_dist(1, width_limit);
    right = left + width_dist(rng) - 1;
  }
  return dataset;
}

template <class Rmq>
void run_queries(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const Dataset dataset = make_dataset(size, max_width);
  const Rmq rmq(std::span<const std::int64_t>(dataset.values));

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] =
        dataset.ranges[query_index++ % dataset.ranges.size()];
    std::size_t result = rmq.arg_min(left, right);
    benchmark::DoNotOptimize(result);
  }

  state.counters["N"] = static_cast<double>(size);
  state.counters["max_width"] = static_cast<double>(max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
}

template <class Rmq>
void run_depth_queries(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const Rmq rmq(std::span<const std::uint64_t>(dataset.bits),
                dataset.depths.size());

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] =
        dataset.ranges[query_index++ % dataset.ranges.size()];
    std::size_t result = rmq.arg_min(left, right);
    benchmark::DoNotOptimize(result);
  }

  state.counters["N"] = static_cast<double>(size);
  state.counters["max_width"] = static_cast<double>(max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
}

void register_benchmarks() {
  const std::vector<std::size_t> sizes = {1ull << 10, 1ull << 14, 1ull << 18,
                                          1ull << 22, 1ull << 26};
  const std::vector<std::size_t> widths = {64, 4096, 1ull << 18, 1ull << 22,
                                           1ull << 26};

  for (const std::size_t size : sizes) {
    std::vector<std::size_t> effective_widths;
    for (const std::size_t width : widths) {
      if (width > size) {
        continue;
      }
      effective_widths.push_back(width);
    }
    effective_widths.push_back(size);
    std::sort(effective_widths.begin(), effective_widths.end());
    effective_widths.erase(
        std::unique(effective_widths.begin(), effective_widths.end()),
        effective_widths.end());

    for (const std::size_t width : effective_widths) {
      benchmark::RegisterBenchmark(
          "rmq_sparse_table",
          &run_queries<pixie::rmq::SparseTable<std::int64_t,
                                               std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_segment_tree",
          &run_queries<pixie::rmq::SegmentTree<std::int64_t,
                                               std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_cartesian_tree",
          &run_queries<pixie::rmq::CartesianTreeRmq<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_plus_minus_one",
          &run_depth_queries<pixie::rmq::BpPlusMinusOneRmq<Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  benchmark::MaybeReenterWithoutASLR(argc, argv);
  benchmark::Initialize(&argc, argv);
  register_benchmarks();
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
