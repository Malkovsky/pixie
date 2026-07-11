#include <benchmark/benchmark.h>
#include <pixie/rmq.h>
#include <pixie/rmq/cartesian_hybrid_btree.h>
#include <pixie/rmq/cartesian_rmm.h>

#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
#include <pixie/rmq/sdsl_sct.h>
#endif

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <functional>
#include <limits>
#include <mutex>
#include <random>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kQueryPoolBytes = 512 * 1024;
constexpr std::size_t kQueryCount =
    kQueryPoolBytes / sizeof(std::pair<std::size_t, std::size_t>);
constexpr std::int64_t kValueDatasetGlobalMinimum =
    std::numeric_limits<std::int64_t>::min() / 4;
constexpr double kBenchmarkWarmupSeconds = 0.2;
constexpr double kBenchmarkMinSeconds = 2.0;
using Index = std::size_t;

static_assert(kQueryPoolBytes % sizeof(std::pair<std::size_t, std::size_t>) ==
              0);

template <class Rmq>
concept HasRmqMemoryUsage = requires(const Rmq& rmq) {
  { rmq.memory_usage_bytes() } -> std::convertible_to<std::size_t>;
};

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9E3779B97F4A7C15ull;
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ull;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBull;
  return value ^ (value >> 31);
}

std::uint64_t mix_seed(std::uint64_t seed, std::uint64_t value) {
  return splitmix64(seed ^ splitmix64(value));
}

class BenchmarkRepetitionSeeds {
 public:
  std::uint64_t repetition_index(const benchmark::State& state,
                                 std::size_t size,
                                 std::size_t max_width) {
    std::lock_guard<std::mutex> lock(mutex_);
    Entry& entry = entries_[key_for(state, size, max_width)];
    const benchmark::IterationCount max_iterations = state.max_iterations;
    if (!entry.seen) {
      entry.seen = true;
      entry.last_max_iterations = max_iterations;
      return entry.repetition_index;
    }

    if (max_iterations == entry.last_max_iterations) {
      if (!entry.saw_iteration_count_change &&
          !entry.skipped_initial_warmup_equal) {
        entry.skipped_initial_warmup_equal = true;
      } else {
        ++entry.repetition_index;
      }
    } else {
      entry.saw_iteration_count_change = true;
    }
    entry.last_max_iterations = max_iterations;
    return entry.repetition_index;
  }

 private:
  struct Entry {
    benchmark::IterationCount last_max_iterations = 0;
    std::uint64_t repetition_index = 0;
    bool seen = false;
    bool saw_iteration_count_change = false;
    bool skipped_initial_warmup_equal = kBenchmarkWarmupSeconds == 0.0;
  };

  std::string key_for(const benchmark::State& state,
                      std::size_t size,
                      std::size_t max_width) const {
    return state.name() + "/" + std::to_string(size) + "/" +
           std::to_string(max_width);
  }

  std::mutex mutex_;
  std::unordered_map<std::string, Entry> entries_;
};

BenchmarkRepetitionSeeds& repetition_seeds() {
  static BenchmarkRepetitionSeeds seeds;
  return seeds;
}

struct SeedContext {
  std::uint64_t repetition_index = 0;
  std::uint64_t value_seed = 0;
  std::uint64_t query_seed = 0;
};

SeedContext make_seed_context(const benchmark::State& state,
                              std::size_t size,
                              std::size_t max_width = 0) {
  const std::uint64_t repetition_index =
      repetition_seeds().repetition_index(state, size, max_width);
  std::uint64_t seed = mix_seed(kSeed, static_cast<std::uint64_t>(size));
  seed = mix_seed(seed, repetition_index);

  SeedContext context;
  context.repetition_index = repetition_index;
  context.value_seed = mix_seed(seed, 0x4F1BBCDCBFA54005ull);
  context.query_seed =
      mix_seed(mix_seed(seed, static_cast<std::uint64_t>(max_width)),
               0x9D2C5680D3F6C58Bull);
  return context;
}

class RandomQueryGenerator {
 public:
  RandomQueryGenerator(std::size_t size,
                       std::size_t max_width,
                       std::uint64_t seed)
      : size_(size), rng_state_(seed), pool_(kQueryCount) {
    std::mt19937_64 rng(seed);
    const std::size_t width_limit = std::min(max_width, size_);
    std::uniform_int_distribution<std::size_t> width_dist(1, width_limit);
    for (auto& [left, width] : pool_) {
      width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size_ - width);
      left = left_dist(rng);
    }
    shift_ = next_shift();
  }

  std::pair<std::size_t, std::size_t> get_random_query() {
    if (index_ == pool_.size()) {
      index_ = 0;
      shift_ = next_shift();
    }
    const std::size_t query_number = query_number_++;
    const auto [left, width] = pool_[index_++];
    const std::size_t shifted_left =
        shifted_left_position(left, width, query_number);
    return {shifted_left, shifted_left + width};
  }

 private:
  std::uint64_t splitmix64() {
    std::uint64_t z = (rng_state_ += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  std::size_t next_shift() { return static_cast<std::size_t>(splitmix64()); }

  std::size_t shifted_left_position(std::size_t left,
                                    std::size_t width,
                                    std::size_t query_number) const {
    const std::size_t valid_left_count = size_ - width + 1;
    const std::size_t offset =
        (shift_ % valid_left_count + query_number % valid_left_count) %
        valid_left_count;
    return (left + offset) % valid_left_count;
  }

  std::size_t size_ = 0;
  std::uint64_t rng_state_ = 0;
  std::vector<std::pair<std::size_t, std::size_t>> pool_;
  std::size_t index_ = 0;
  std::size_t query_number_ = 0;
  std::size_t shift_ = 0;
};

/**
 * @brief Single generated value dataset for value-RMQ benchmarks.
 */
class ValueDataset {
 public:
  explicit ValueDataset(std::size_t size, std::uint64_t seed) : values_(size) {
    size_ = size;

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> value_dist(-1'000'000,
                                                           1'000'000);
    std::generate(values_.begin(), values_.end(),
                  [&] { return value_dist(rng); });

    if (!values_.empty()) {
      values_[values_.size() / 2] = kValueDatasetGlobalMinimum;
    }
  }

  std::size_t size() const { return size_; }

  std::span<const std::int64_t> values() const { return values_; }

 private:
  std::size_t size_ = 0;
  std::vector<std::int64_t> values_;
};

void set_aux_memory_counters(benchmark::State& state,
                             std::size_t size,
                             double aux_bytes);

template <class Rmq>
void set_query_aux_memory_counters(benchmark::State& state,
                                   std::size_t size,
                                   const Rmq& rmq);

/**
 * @brief Run value-RMQ query benchmarks over one live RMQ instance.
 */
template <class Rmq>
void run_queries(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const SeedContext seeds = make_seed_context(state, size, max_width);
  ValueDataset dataset(size, seeds.value_seed);
  Rmq rmq(dataset.values());
  RandomQueryGenerator queries(size, max_width, seeds.query_seed);

  for (auto _ : state) {
    const auto [left, right] = queries.get_random_query();
    std::size_t result = rmq.arg_min(left, right);
    benchmark::DoNotOptimize(result);
  }

  state.counters["N"] = static_cast<double>(size);
  state.counters["max_width"] = static_cast<double>(max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
  state.counters["seed_repetition"] =
      static_cast<double>(seeds.repetition_index);
  state.counters["value_arrays"] = 1.0;
  set_query_aux_memory_counters(state, size, rmq);
}

void set_build_counters(benchmark::State& state,
                        std::size_t size,
                        std::size_t input_bytes) {
  state.counters["N"] = static_cast<double>(size);
  state.counters["input_bytes"] = static_cast<double>(input_bytes);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

void set_aux_memory_counters(benchmark::State& state,
                             std::size_t size,
                             double aux_bytes) {
  state.counters["aux_bytes"] = aux_bytes;
  state.counters["aux_mib"] = aux_bytes / (1024.0 * 1024.0);
  state.counters["aux_bits_per_value"] =
      size == 0 ? 0.0 : (8.0 * aux_bytes) / static_cast<double>(size);
}

template <class Rmq>
void set_query_aux_memory_counters(benchmark::State& state,
                                   std::size_t size,
                                   const Rmq& rmq) {
  if constexpr (HasRmqMemoryUsage<Rmq>) {
    set_aux_memory_counters(state, size,
                            static_cast<double>(rmq.memory_usage_bytes()));
  }
}

/**
 * @brief Run value-RMQ construction benchmarks with one value array.
 */
template <class Rmq>
void run_value_rmq_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size);
  ValueDataset dataset(size, seeds.value_seed);

  for (auto _ : state) {
    Rmq rmq(dataset.values());
    std::size_t built_size = rmq.size();
    benchmark::DoNotOptimize(built_size);
    benchmark::ClobberMemory();
  }

  set_build_counters(state, size, dataset.size() * sizeof(std::int64_t));
  state.counters["seed_repetition"] =
      static_cast<double>(seeds.repetition_index);
  state.counters["value_arrays"] = 1.0;
  if constexpr (HasRmqMemoryUsage<Rmq>) {
    const Rmq rmq(dataset.values());
    const std::size_t aux_bytes = rmq.memory_usage_bytes();
    set_aux_memory_counters(state, size, static_cast<double>(aux_bytes));
  }
}

void register_benchmarks() {
  using CartesianRmM =
      pixie::rmq::CartesianRmM<std::int64_t, std::less<std::int64_t>, Index>;
  using CartesianHybrid =
      pixie::rmq::CartesianHybridBTree<std::int64_t, std::less<std::int64_t>,
                                       Index>;
  using CartesianBTree =
      pixie::rmq::CartesianBTree<std::int64_t, std::less<std::int64_t>, Index>;

  const std::vector<std::size_t> sizes = {1ull << 10, 1ull << 14, 1ull << 18,
                                          1ull << 22, 1ull << 24, 1ull << 26};
  const std::vector<std::size_t> build_sizes = {
      1ull << 10, 1ull << 14, 1ull << 18, 1ull << 22, 1ull << 24, 1ull << 26};
  const std::vector<std::size_t> widths = {64, 4096, 1ull << 18, 1ull << 22,
                                           1ull << 26};
  // Pure sparse-table rows materialize O(n log n) indexes; above 2^22 they
  // dominate memory/runtime and make large benchmark passes noisy.
  constexpr std::size_t kSparseTableBenchmarkMaxSize = 1ull << 22;
  constexpr std::size_t kFocusedCartesianRmMLargeSize = 1ull << 28;
  constexpr std::size_t kFocusedCartesianHybridLargeSize = 1ull << 28;
  constexpr std::size_t kFocusedCartesianBTreeLargeSize = 1ull << 28;
  constexpr std::size_t kFocusedHybridBTreeLargeSize = 1ull << 28;
  // Focused 2^28/2^30 rows are intentionally limited to the RMQ variants whose
  // construction and query memory stay bounded at those sizes. CartesianRmM is
  // included here because its Cartesian BP construction now uses a succinct
  // monotone bit-stack instead of an n-entry index stack.
  constexpr std::size_t kVeryLargeHybridSize = 1ull << 30;

  auto effective_widths_for = [&](std::size_t size) {
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
    return effective_widths;
  };

  for (const std::size_t size : build_sizes) {
    if (size <= kSparseTableBenchmarkMaxSize) {
      benchmark::RegisterBenchmark(
          "rmq_build_sparse_table",
          &run_value_rmq_build<pixie::rmq::SparseTable<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Arg(static_cast<std::int64_t>(size))
          ->Unit(benchmark::kMillisecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
    }
    benchmark::RegisterBenchmark(
        "rmq_build_segment_tree",
        &run_value_rmq_build<pixie::rmq::SegmentTree<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
    benchmark::RegisterBenchmark("rmq_build_cartesian_rmm",
                                 &run_value_rmq_build<CartesianRmM>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
    benchmark::RegisterBenchmark("rmq_build_cartesian_hybrid_btree",
                                 &run_value_rmq_build<CartesianHybrid>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
    benchmark::RegisterBenchmark("rmq_build_cartesian_btree",
                                 &run_value_rmq_build<CartesianBTree>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
    benchmark::RegisterBenchmark(
        "rmq_build_sdsl_sct",
        &run_value_rmq_build<
            pixie::rmq::SdslSct<std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
#endif
    benchmark::RegisterBenchmark(
        "rmq_build_hybrid_btree",
        &run_value_rmq_build<pixie::rmq::HybridBTree<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }

  benchmark::RegisterBenchmark("rmq_build_cartesian_rmm",
                               &run_value_rmq_build<CartesianRmM>)
      ->Arg(static_cast<std::int64_t>(kFocusedCartesianRmMLargeSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark("rmq_build_cartesian_rmm",
                               &run_value_rmq_build<CartesianRmM>)
      ->Arg(static_cast<std::int64_t>(kVeryLargeHybridSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark("rmq_build_cartesian_hybrid_btree",
                               &run_value_rmq_build<CartesianHybrid>)
      ->Arg(static_cast<std::int64_t>(kFocusedCartesianHybridLargeSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark("rmq_build_cartesian_hybrid_btree",
                               &run_value_rmq_build<CartesianHybrid>)
      ->Arg(static_cast<std::int64_t>(kVeryLargeHybridSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark("rmq_build_cartesian_btree",
                               &run_value_rmq_build<CartesianBTree>)
      ->Arg(static_cast<std::int64_t>(kFocusedCartesianBTreeLargeSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark("rmq_build_cartesian_btree",
                               &run_value_rmq_build<CartesianBTree>)
      ->Arg(static_cast<std::int64_t>(kVeryLargeHybridSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
  benchmark::RegisterBenchmark(
      "rmq_build_sdsl_sct",
      &run_value_rmq_build<
          pixie::rmq::SdslSct<std::int64_t, std::less<std::int64_t>, Index>>)
      ->Arg(static_cast<std::int64_t>(kFocusedCartesianRmMLargeSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
#endif
  benchmark::RegisterBenchmark(
      "rmq_build_hybrid_btree",
      &run_value_rmq_build<pixie::rmq::HybridBTree<
          std::int64_t, std::less<std::int64_t>, Index>>)
      ->Arg(static_cast<std::int64_t>(kFocusedHybridBTreeLargeSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
  benchmark::RegisterBenchmark(
      "rmq_build_hybrid_btree",
      &run_value_rmq_build<pixie::rmq::HybridBTree<
          std::int64_t, std::less<std::int64_t>, Index>>)
      ->Arg(static_cast<std::int64_t>(kVeryLargeHybridSize))
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);

  for (const std::size_t size : sizes) {
    const std::vector<std::size_t> effective_widths =
        effective_widths_for(size);

    for (const std::size_t width : effective_widths) {
      if (size <= kSparseTableBenchmarkMaxSize) {
        benchmark::RegisterBenchmark(
            "rmq_sparse_table",
            &run_queries<pixie::rmq::SparseTable<
                std::int64_t, std::less<std::int64_t>, Index>>)
            ->Args({static_cast<std::int64_t>(size),
                    static_cast<std::int64_t>(width)})
            ->Unit(benchmark::kNanosecond)
            ->MinWarmUpTime(kBenchmarkWarmupSeconds)
            ->MinTime(kBenchmarkMinSeconds);
      }
      benchmark::RegisterBenchmark(
          "rmq_segment_tree",
          &run_queries<pixie::rmq::SegmentTree<std::int64_t,
                                               std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
      benchmark::RegisterBenchmark("rmq_cartesian_rmm",
                                   &run_queries<CartesianRmM>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
      benchmark::RegisterBenchmark("rmq_cartesian_hybrid_btree",
                                   &run_queries<CartesianHybrid>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
      benchmark::RegisterBenchmark("rmq_cartesian_btree",
                                   &run_queries<CartesianBTree>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
      benchmark::RegisterBenchmark(
          "rmq_sdsl_sct",
          &run_queries<pixie::rmq::SdslSct<std::int64_t,
                                           std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
#endif
      benchmark::RegisterBenchmark(
          "rmq_hybrid_btree",
          &run_queries<pixie::rmq::HybridBTree<std::int64_t,
                                               std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond)
          ->MinWarmUpTime(kBenchmarkWarmupSeconds)
          ->MinTime(kBenchmarkMinSeconds);
    }
  }

  for (const std::size_t width :
       effective_widths_for(kFocusedCartesianRmMLargeSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_rmm",
                                 &run_queries<CartesianRmM>)
        ->Args({static_cast<std::int64_t>(kFocusedCartesianRmMLargeSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width : effective_widths_for(kVeryLargeHybridSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_rmm",
                                 &run_queries<CartesianRmM>)
        ->Args({static_cast<std::int64_t>(kVeryLargeHybridSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width :
       effective_widths_for(kFocusedCartesianHybridLargeSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_hybrid_btree",
                                 &run_queries<CartesianHybrid>)
        ->Args({static_cast<std::int64_t>(kFocusedCartesianHybridLargeSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width : effective_widths_for(kVeryLargeHybridSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_hybrid_btree",
                                 &run_queries<CartesianHybrid>)
        ->Args({static_cast<std::int64_t>(kVeryLargeHybridSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width :
       effective_widths_for(kFocusedCartesianBTreeLargeSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_btree",
                                 &run_queries<CartesianBTree>)
        ->Args({static_cast<std::int64_t>(kFocusedCartesianBTreeLargeSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width : effective_widths_for(kVeryLargeHybridSize)) {
    benchmark::RegisterBenchmark("rmq_cartesian_btree",
                                 &run_queries<CartesianBTree>)
        ->Args({static_cast<std::int64_t>(kVeryLargeHybridSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width :
       effective_widths_for(kFocusedHybridBTreeLargeSize)) {
    benchmark::RegisterBenchmark(
        "rmq_hybrid_btree",
        &run_queries<pixie::rmq::HybridBTree<std::int64_t,
                                             std::less<std::int64_t>, Index>>)
        ->Args({static_cast<std::int64_t>(kFocusedHybridBTreeLargeSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
  }
  for (const std::size_t width : effective_widths_for(kVeryLargeHybridSize)) {
    benchmark::RegisterBenchmark(
        "rmq_hybrid_btree",
        &run_queries<pixie::rmq::HybridBTree<std::int64_t,
                                             std::less<std::int64_t>, Index>>)
        ->Args({static_cast<std::int64_t>(kVeryLargeHybridSize),
                static_cast<std::int64_t>(width)})
        ->Unit(benchmark::kNanosecond)
        ->MinWarmUpTime(kBenchmarkWarmupSeconds)
        ->MinTime(kBenchmarkMinSeconds);
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
