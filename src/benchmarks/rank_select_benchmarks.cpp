#include <benchmark/benchmark.h>
#include <pixie/rank_select/implementations.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>
#include <pixie/storage/implementations.h>

#ifdef PIXIE_3STAR_SUPPORT
#include "three_star_adapter.h"
#endif

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kQueryPoolBytes = 512 * 1024;
constexpr std::size_t kQueryCount = kQueryPoolBytes / sizeof(std::size_t);
constexpr double kBenchmarkWarmupSeconds = 0.2;
constexpr double kBenchmarkMinSeconds = 1.0;
constexpr std::array<std::size_t, 7> kSizes = {
    1ull << 10, 1ull << 14, 1ull << 18, 1ull << 22,
    1ull << 26, 1ull << 30, 1ull << 34};
// FNBP construction is linear in encoded bits. Include a 2^30-bit case to
// expose large-source behavior while retaining the existing scaling points.
constexpr std::array<std::size_t, 6> kFnbpSizes = {
    1ull << 10, 1ull << 14, 1ull << 18, 1ull << 22, 1ull << 26, 1ull << 30};
using PixieRankSelect = pixie::RankSelectSupport<>;

static_assert(kQueryCount > 0 && (kQueryCount & (kQueryCount - 1)) == 0);

enum class Fill {
  k12p5,
  k50,
  k87p5,
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
      repetition_seeds().repetition_index(state, size, fill);
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

std::vector<std::size_t> make_query_pool(std::size_t first,
                                         std::size_t last,
                                         std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<std::size_t> distribution(first, last);
  std::vector<std::size_t> queries(kQueryCount);
  for (std::size_t& query : queries) {
    query = distribution(rng);
  }
  return queries;
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

  std::vector<std::size_t> queries;
  if constexpr (operation == QueryOperation::kRank1 ||
                operation == QueryOperation::kRank0) {
    queries = make_query_pool(0, size, seeds.query_seed);
  } else if constexpr (operation == QueryOperation::kSelect1) {
    if (one_count == 0) {
      state.SkipWithError("input has no one bits");
      return;
    }
    queries = make_query_pool(1, one_count, seeds.query_seed);
  } else {
    if (zero_count == 0) {
      state.SkipWithError("input has no zero bits");
      return;
    }
    queries = make_query_pool(1, zero_count, seeds.query_seed);
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const std::size_t query = queries[query_index++ & (kQueryCount - 1)];
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

  std::vector<std::size_t> queries;
  if constexpr (operation == QueryOperation::kRank1 ||
                operation == QueryOperation::kRank0) {
    queries = make_query_pool(0, size, seeds.query_seed);
  } else if constexpr (operation == QueryOperation::kSelect1) {
    queries = make_query_pool(1, one_count, seeds.query_seed);
  } else {
    queries = make_query_pool(1, zero_count, seeds.query_seed);
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const std::size_t query = queries[query_index++ & (kQueryCount - 1)];
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
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
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
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
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

void register_benchmarks() {
  register_fill_rows<PixieRankSelect, Fill::k12p5>("rank_select");
  register_fill_rows<PixieRankSelect, Fill::k50>("rank_select");
  register_fill_rows<PixieRankSelect, Fill::k87p5>("rank_select");
  register_fnbp_rows<PixieRankSelect>("rank_select");
#ifdef PIXIE_PASTA_SUPPORT
  register_fill_rows<pixie::PastaRankSelectSupport, Fill::k12p5>(
      "rank_select_pasta");
  register_fill_rows<pixie::PastaRankSelectSupport, Fill::k50>(
      "rank_select_pasta");
  register_fill_rows<pixie::PastaRankSelectSupport, Fill::k87p5>(
      "rank_select_pasta");
  register_fnbp_rows<pixie::PastaRankSelectSupport>("rank_select_pasta");
#endif
#ifdef PIXIE_3STAR_SUPPORT
  using ThreeStarRankSelect = pixie::benchmarks::ThreeStarRankSelectSupport;
  register_select1_only_fill_rows<ThreeStarRankSelect, Fill::k12p5>(
      "rank_select_3star");
  register_select1_only_fill_rows<ThreeStarRankSelect, Fill::k50>(
      "rank_select_3star");
  register_select1_only_fill_rows<ThreeStarRankSelect, Fill::k87p5>(
      "rank_select_3star");
  register_fnbp_rows<ThreeStarRankSelect, false>("rank_select_3star");
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
