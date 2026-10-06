#include <benchmark/benchmark.h>
#include <pixie/rank_select/implementations.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>
#include <pixie/storage/implementations.h>

#include "perf_query_window.h"
#include "shifted_query_pool.h"

#ifdef PIXIE_3STAR_SUPPORT
#include "three_star_adapter.h"
#endif

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr double kBenchmarkWarmupSeconds = 0.2;
constexpr double kBenchmarkMinSeconds = 1.0;

// Fix repetition zero for reproducible cross-backend or independent-process
// comparisons; ordinary repeated runs continue varying the generated input.
bool fixed_repetition_seed() {
  return std::getenv("PIXIE_RANK_SELECT_FIXED_SEED") != nullptr;
}

double query_warmup_seconds() {
  return pixie::benchmarks::PerfQueryWindow::requested()
             ? 0.0
             : kBenchmarkWarmupSeconds;
}
constexpr std::array<std::size_t, 7> kSizes = {
    1ull << 10, 1ull << 14, 1ull << 18, 1ull << 22,
    1ull << 26, 1ull << 30, 1ull << 34};
// All source families share the same scaling points, including the final
// large-source FNBP comparison.
constexpr auto kFnbpSizes = kSizes;

enum class Fill {
  k12p5,
  k50,
  k87p5,
  kSparse,
  kSparseZeros,
  kRuns,
  kMixed,
  kMixedZeros,
};

enum class QueryOperation {
  kRank1,
  kRank0,
  kSelect1,
  kSelect0,
};

struct FillSpec {
  std::string_view name;
  double expected_one_fill_percent;
  std::uint64_t seed_tag;
};

constexpr FillSpec fill_spec(Fill fill) {
  switch (fill) {
    case Fill::k12p5:
      return {"12p5", 12.5, 0xC6A4A7935BD1E995ull};
    case Fill::k50:
      return {"50", 50.0, 0x9E3779B97F4A7C15ull};
    case Fill::k87p5:
      return {"87p5", 87.5, 0xD6E8FEB86659FD93ull};
    case Fill::kSparse:
      return {"sparse", 100.0 / 1024, 0x8CB92BA72F3D8DD7ull};
    case Fill::kSparseZeros:
      return {"sparse0", 100.0 - 100.0 / 1024, 0xDB4F0B9175AE2165ull};
    case Fill::kRuns:
      return {"runs", 50.0, 0xBBE0563303A4615Full};
    case Fill::kMixed:
      return {"mixed", (50.0 + 64.0 * 100.0 / 1024) / 65,
              0xD42E3F8A111B1201ull};
    case Fill::kMixedZeros:
      return {"mixed0", 100.0 - (50.0 + 64.0 * 100.0 / 1024) / 65,
              0xD42E3F8A111B1201ull};
  }
  return {"unknown", 0.0, 0};
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9E3779B97F4A7C15ull;
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ull;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBull;
  return value ^ (value >> 31);
}

std::uint64_t mix_seed(std::uint64_t seed, std::uint64_t value) {
  return splitmix64(seed ^ splitmix64(value));
}

std::uint64_t fnbp_priority(std::uint64_t seed, std::size_t value_index) {
  return splitmix64(seed + static_cast<std::uint64_t>(value_index));
}

void fill_fnbp_words(std::span<std::uint64_t> words,
                     std::size_t bit_count,
                     std::uint64_t seed) {
  std::fill(words.begin(), words.end(), 0);

  // This is the same backward monotone-stack construction used for the
  // Ferrada-Navarro BP encoding in CartesianHybridBTree.
  const std::size_t value_count = bit_count / 2;
  pixie::rmq::utils::SuccinctIncreasingStack monotone_stack(value_count);
  std::size_t write_position = bit_count;

  const auto prepend_open = [&] {
    --write_position;
    words[write_position >> 6] |= std::uint64_t{1} << (write_position & 63);
  };

  for (std::size_t value_index = value_count; value_index > 0; --value_index) {
    const std::size_t current_index = value_index - 1;
    const std::uint64_t value = fnbp_priority(seed, current_index);
    while (!monotone_stack.empty() &&
           fnbp_priority(seed, value_count - 1 - monotone_stack.top()) >=
               value) {
      monotone_stack.pop();
      prepend_open();
    }
    monotone_stack.push(value_count - 1 - current_index);
    --write_position;  // Closing parenthesis.
  }

  while (write_position != 0) {
    prepend_open();
  }
}

class BenchmarkRepetitionSeeds {
 public:
  std::uint64_t repetition_index(const benchmark::State& state,
                                 std::size_t size,
                                 Fill fill) {
    std::lock_guard<std::mutex> lock(mutex_);
    Entry& entry = entries_[key_for(state, size, fill)];
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
                      Fill fill) const {
    return state.name() + "/" + std::to_string(size) + "/" +
           std::string(fill_spec(fill).name);
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
  std::uint64_t source_seed = 0;
  std::uint64_t query_seed = 0;
};

SeedContext make_seed_context(const benchmark::State& state,
                              std::size_t size,
                              Fill fill) {
  const std::uint64_t repetition_index =
      fixed_repetition_seed()
          ? 0
          : repetition_seeds().repetition_index(state, size, fill);
  std::uint64_t seed = mix_seed(kSeed, static_cast<std::uint64_t>(size));
  seed = mix_seed(seed, fill_spec(fill).seed_tag);
  seed = mix_seed(seed, repetition_index);

  return {
      .repetition_index = repetition_index,
      .source_seed = mix_seed(seed, 0x4F1BBCDCBFA54005ull),
      .query_seed = mix_seed(seed, 0x9D2C5680D3F6C58Bull),
  };
}

void fill_words(std::span<std::uint64_t> words, Fill fill, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::size_t word_index = 0;
  const std::size_t run_words =
      std::max<std::size_t>(1, std::min<std::size_t>(16384, words.size() / 2));
  for (std::uint64_t& word : words) {
    switch (fill) {
      case Fill::k12p5:
        word = rng() & rng() & rng();
        break;
      case Fill::k50:
        word = rng();
        break;
      case Fill::k87p5:
        word = rng() | rng() | rng();
        break;
      case Fill::kSparse:
      case Fill::kSparseZeros:
        word = rng();
        for (unsigned i = 1; i < 10; ++i) {
          word &= rng();
        }
        if (fill == Fill::kSparseZeros) {
          word = ~word;
        }
        break;
      case Fill::kRuns:
        word = ((word_index++ / run_words) & 1) == 0 ? ~std::uint64_t{0} : 0;
        break;
      case Fill::kMixed:
      case Fill::kMixedZeros:
        // One half-full superblock, then 64 sparse superblocks. Both narrow
        // and wide intervals recur; complementing exercises zero-select too.
        word = rng();
        if ((word_index++ / 256) % 65 != 0) {
          for (unsigned i = 1; i < 10; ++i) {
            word &= rng();
          }
        }
        if (fill == Fill::kMixedZeros) {
          word = ~word;
        }
        break;
    }
  }
}

class BitDataset {
 public:
  BitDataset(std::size_t size, Fill fill, std::uint64_t seed)
      : source_(size), size_(size) {
    fill_words(source_.writable_words64(), fill, seed);
  }

  const pixie::AlignedStorage& source() const { return source_; }
  std::size_t size() const { return size_; }

 private:
  pixie::AlignedStorage source_;
  std::size_t size_ = 0;
};

class FnbpDataset {
 public:
  FnbpDataset(std::size_t size, std::uint64_t seed)
      : source_(size), size_(size) {
    fill_fnbp_words(source_.writable_words64(), size_, seed);
  }

  const pixie::AlignedStorage& source() const { return source_; }
  std::size_t size() const { return size_; }

 private:
  pixie::AlignedStorage source_;
  std::size_t size_ = 0;
};

template <class Support, class Dataset>
Support make_support(const Dataset& dataset) {
  return Support(dataset.source().padded_view().as_words64(), dataset.size());
}

template <class Support>
constexpr bool backend_owns_source() {
  return requires(const Support& support) { support.source_copy_bytes(); };
}

template <class Support>
std::size_t index_logical_bytes(const Support& support) {
  if constexpr (requires { support.index_logical_bytes(); }) {
    return support.index_logical_bytes();
  }
  return support.memory_usage_bytes();
}

template <class Support>
std::size_t source_copy_bytes(const Support& support) {
  if constexpr (requires { support.source_copy_bytes(); }) {
    return support.source_copy_bytes();
  }
  return 0;
}

template <class Support>
bool supports_select1(const Support& support) {
  return support.supports_select1();
}

template <class Support>
bool supports_select0(const Support& support) {
  return support.supports_select0();
}

void set_common_counters(benchmark::State& state,
                         std::size_t size,
                         Fill fill,
                         std::size_t index_bytes,
                         std::size_t owned_bytes,
                         std::size_t source_copy_bytes,
                         bool owns_source,
                         bool select1_enabled,
                         bool select0_enabled,
                         std::uint64_t repetition_index) {
  const double input_bytes = static_cast<double>((size + 7) / 8);
  const double index_bytes_as_double = static_cast<double>(index_bytes);
  state.counters["N"] = static_cast<double>(size);
  state.counters["one_fill_percent"] =
      fill_spec(fill).expected_one_fill_percent;
  state.counters["input_bytes"] = input_bytes;
  state.counters["aux_bytes"] = index_bytes_as_double;
  state.counters["aux_mib"] = index_bytes_as_double / (1024.0 * 1024.0);
  state.counters["aux_bits_per_input_bit"] =
      size == 0 ? 0.0 : 8.0 * index_bytes_as_double / size;
  state.counters["backend_owned_bytes"] = static_cast<double>(owned_bytes);
  state.counters["source_copy_bytes"] = static_cast<double>(source_copy_bytes);
  state.counters["select1_enabled"] = select1_enabled ? 1.0 : 0.0;
  state.counters["select0_enabled"] = select0_enabled ? 1.0 : 0.0;
  state.counters["seed_repetition"] = static_cast<double>(repetition_index);
  state.counters["backend_owns_source"] = owns_source ? 1.0 : 0.0;
}

void set_fnbp_counters(benchmark::State& state,
                       std::size_t size,
                       std::size_t index_bytes,
                       std::size_t owned_bytes,
                       std::size_t source_copy_bytes,
                       bool owns_source,
                       bool select1_enabled,
                       bool select0_enabled,
                       std::uint64_t repetition_index) {
  set_common_counters(state, size, Fill::k50, index_bytes, owned_bytes,
                      source_copy_bytes, owns_source, select1_enabled,
                      select0_enabled, repetition_index);
  state.counters["fnbp_source"] = 1.0;
}

template <class Support, Fill fill>
void run_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, fill);
  const BitDataset dataset(size, fill, seeds.source_seed);
  std::size_t index_bytes = 0;
  std::size_t owned_bytes = 0;
  std::size_t copied_source_bytes = 0;
  bool select1_enabled = false;
  bool select0_enabled = false;

  for (auto _ : state) {
    Support support = make_support<Support>(dataset);
    index_bytes = index_logical_bytes(support);
    owned_bytes = support.memory_usage_bytes();
    copied_source_bytes = source_copy_bytes(support);
    select1_enabled = supports_select1(support);
    select0_enabled = supports_select0(support);
    benchmark::DoNotOptimize(owned_bytes);
    benchmark::ClobberMemory();
  }

  set_common_counters(state, size, fill, index_bytes, owned_bytes,
                      copied_source_bytes, backend_owns_source<Support>(),
                      select1_enabled, select0_enabled, seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

template <class Support, Fill fill, QueryOperation operation>
void run_query(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, fill);
  const BitDataset dataset(size, fill, seeds.source_seed);
  const Support support = make_support<Support>(dataset);
  const std::size_t one_count = support.rank(support.size());
  const std::size_t zero_count = support.rank0(support.size());

  // Rank queries include both prefix boundaries: [0, size].
  std::size_t first_query = 0;
  std::size_t last_query = size;
  if constexpr (operation == QueryOperation::kSelect1) {
    if (one_count == 0) {
      state.SkipWithError("input has no one bits");
      return;
    }
    first_query = 1;
    last_query = one_count;
  } else if constexpr (operation == QueryOperation::kSelect0) {
    if (zero_count == 0) {
      state.SkipWithError("input has no zero bits");
      return;
    }
    first_query = 1;
    last_query = zero_count;
  }

  pixie::benchmarks::ShiftedQueryPool queries(first_query, last_query,
                                              seeds.query_seed);
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  std::size_t query_index = 0;
  perf.start();
  for (auto _ : state) {
    const std::size_t query = queries.next();
    ++query_index;
    if constexpr (operation == QueryOperation::kRank1) {
      benchmark::DoNotOptimize(support.rank(query));
    } else if constexpr (operation == QueryOperation::kRank0) {
      benchmark::DoNotOptimize(support.rank0(query));
    } else if constexpr (operation == QueryOperation::kSelect1) {
      benchmark::DoNotOptimize(support.select(query));
    } else {
      benchmark::DoNotOptimize(support.select0(query));
    }
  }

  perf.stop(query_index);
  if (perf.active()) {
    state.counters["perf_queries"] = static_cast<double>(perf.queries());
    state.counters["perf_phases"] = static_cast<double>(perf.phases());
  }
  set_common_counters(state, size, fill, index_logical_bytes(support),
                      support.memory_usage_bytes(), source_copy_bytes(support),
                      backend_owns_source<Support>(), supports_select1(support),
                      supports_select0(support), seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
}

template <class Support>
void run_fnbp_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, Fill::k50);
  const FnbpDataset dataset(size, seeds.source_seed);
  std::size_t index_bytes = 0;
  std::size_t owned_bytes = 0;
  std::size_t copied_source_bytes = 0;
  bool select1_enabled = false;
  bool select0_enabled = false;

  for (auto _ : state) {
    Support support = make_support<Support>(dataset);
    index_bytes = index_logical_bytes(support);
    owned_bytes = support.memory_usage_bytes();
    copied_source_bytes = source_copy_bytes(support);
    select1_enabled = supports_select1(support);
    select0_enabled = supports_select0(support);
    benchmark::DoNotOptimize(owned_bytes);
    benchmark::ClobberMemory();
  }

  set_fnbp_counters(state, size, index_bytes, owned_bytes, copied_source_bytes,
                    backend_owns_source<Support>(), select1_enabled,
                    select0_enabled, seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

template <class Support, QueryOperation operation>
void run_fnbp_query(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, Fill::k50);
  const FnbpDataset dataset(size, seeds.source_seed);
  const Support support = make_support<Support>(dataset);
  const std::size_t one_count = support.rank(support.size());
  const std::size_t zero_count = support.rank0(support.size());
  if (one_count != size / 2 || zero_count != size / 2) {
    state.SkipWithError("FNBP source is not balanced");
    return;
  }

  // Rank queries include both prefix boundaries: [0, size].
  std::size_t first_query = 0;
  std::size_t last_query = size;
  if constexpr (operation == QueryOperation::kSelect1) {
    first_query = 1;
    last_query = one_count;
  } else if constexpr (operation == QueryOperation::kSelect0) {
    first_query = 1;
    last_query = zero_count;
  }

  pixie::benchmarks::ShiftedQueryPool queries(first_query, last_query,
                                              seeds.query_seed);
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  std::size_t query_index = 0;
  perf.start();
  for (auto _ : state) {
    const std::size_t query = queries.next();
    ++query_index;
    if constexpr (operation == QueryOperation::kRank1) {
      benchmark::DoNotOptimize(support.rank(query));
    } else if constexpr (operation == QueryOperation::kRank0) {
      benchmark::DoNotOptimize(support.rank0(query));
    } else if constexpr (operation == QueryOperation::kSelect1) {
      benchmark::DoNotOptimize(support.select(query));
    } else {
      benchmark::DoNotOptimize(support.select0(query));
    }
  }

  perf.stop(query_index);
  if (perf.active()) {
    state.counters["perf_queries"] = static_cast<double>(perf.queries());
    state.counters["perf_phases"] = static_cast<double>(perf.phases());
  }
  set_fnbp_counters(state, size, index_logical_bytes(support),
                    support.memory_usage_bytes(), source_copy_bytes(support),
                    backend_owns_source<Support>(), supports_select1(support),
                    supports_select0(support), seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
}

template <class Support, Fill fill>
void register_build_row(std::string_view backend_name) {
  const std::string name = std::string(backend_name) + "_build_both_" +
                           std::string(fill_spec(fill).name);
  auto* row =
      benchmark::RegisterBenchmark(name.c_str(), &run_build<Support, fill>);
  for (const std::size_t size : kSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

template <class Support, Fill fill, QueryOperation operation>
void register_query_row(std::string_view backend_name,
                        std::string_view operation_name) {
  const std::string name = std::string(backend_name) + "_" +
                           std::string(operation_name) + "_" +
                           std::string(fill_spec(fill).name);
  auto* row = benchmark::RegisterBenchmark(
      name.c_str(), &run_query<Support, fill, operation>);
  for (const std::size_t size : kSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kNanosecond)
      ->MinWarmUpTime(query_warmup_seconds())
      ->MinTime(kBenchmarkMinSeconds);
}

template <class Support, Fill fill>
void register_fill_rows(std::string_view backend_name) {
  constexpr auto kRank1 = QueryOperation::kRank1;
  constexpr auto kRank0 = QueryOperation::kRank0;
  constexpr auto kSelect1 = QueryOperation::kSelect1;
  constexpr auto kSelect0 = QueryOperation::kSelect0;
  register_build_row<Support, fill>(backend_name);
  register_query_row<Support, fill, kRank1>(backend_name, "rank1");
  register_query_row<Support, fill, kRank0>(backend_name, "rank0");
  register_query_row<Support, fill, kSelect1>(backend_name, "select1");
  register_query_row<Support, fill, kSelect0>(backend_name, "select0");
}

template <class Support>
void register_fnbp_build_row(std::string_view backend_name) {
  const std::string name = std::string(backend_name) + "_fnbp_build_both";
  auto* row =
      benchmark::RegisterBenchmark(name.c_str(), &run_fnbp_build<Support>);
  for (const std::size_t size : kFnbpSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

template <class Support, QueryOperation operation>
void register_fnbp_query_row(std::string_view backend_name,
                             std::string_view operation_name) {
  const std::string name =
      std::string(backend_name) + "_fnbp_" + std::string(operation_name);
  const auto callback = &run_fnbp_query<Support, operation>;
  auto* row = benchmark::RegisterBenchmark(name.c_str(), callback);
  for (const std::size_t size : kFnbpSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kNanosecond)
      ->MinWarmUpTime(query_warmup_seconds())
      ->MinTime(kBenchmarkMinSeconds);
}

template <class Support, bool select0_enabled = true>
void register_fnbp_rows(std::string_view backend_name) {
  constexpr auto kRank1 = QueryOperation::kRank1;
  constexpr auto kRank0 = QueryOperation::kRank0;
  constexpr auto kSelect1 = QueryOperation::kSelect1;
  constexpr auto kSelect0 = QueryOperation::kSelect0;
  register_fnbp_build_row<Support>(backend_name);
  register_fnbp_query_row<Support, kRank1>(backend_name, "rank1");
  register_fnbp_query_row<Support, kRank0>(backend_name, "rank0");
  register_fnbp_query_row<Support, kSelect1>(backend_name, "select1");
  if constexpr (select0_enabled) {
    register_fnbp_query_row<Support, kSelect0>(backend_name, "select0");
  }
}

template <class Support, Fill fill>
void register_select1_only_fill_rows(std::string_view backend_name) {
  constexpr auto kRank1 = QueryOperation::kRank1;
  constexpr auto kRank0 = QueryOperation::kRank0;
  constexpr auto kSelect1 = QueryOperation::kSelect1;
  register_build_row<Support, fill>(backend_name);
  register_query_row<Support, fill, kRank1>(backend_name, "rank1");
  register_query_row<Support, fill, kRank0>(backend_name, "rank0");
  register_query_row<Support, fill, kSelect1>(backend_name, "select1");
}

template <class Support>
void register_backend(std::string_view name) {
  register_fill_rows<Support, Fill::k12p5>(name);
  register_fill_rows<Support, Fill::k50>(name);
  register_fill_rows<Support, Fill::k87p5>(name);
  register_fill_rows<Support, Fill::kSparse>(name);
  register_fill_rows<Support, Fill::kSparseZeros>(name);
  register_fill_rows<Support, Fill::kRuns>(name);
  register_fill_rows<Support, Fill::kMixed>(name);
  register_fill_rows<Support, Fill::kMixedZeros>(name);
  register_fnbp_rows<Support>(name);
}

void register_benchmarks() {
  register_backend<pixie::RankSelectSupport<>>("rank_select");
  register_backend<pixie::RankSelectSupportCompressed<>>(
      "rank_select_compressed");
#ifdef PIXIE_PASTA_SUPPORT
  register_backend<pixie::PastaRankSelectSupport>("rank_select_pasta");
#endif
#ifdef PIXIE_3STAR_SUPPORT
  using Support = pixie::benchmarks::ThreeStarRankSelectSupport;
  register_select1_only_fill_rows<Support, Fill::k12p5>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::k50>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::k87p5>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::kSparse>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::kSparseZeros>(
      "rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::kRuns>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::kMixed>("rank_select_3star");
  register_select1_only_fill_rows<Support, Fill::kMixedZeros>(
      "rank_select_3star");
  register_fnbp_rows<Support, false>("rank_select_3star");
#endif
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
