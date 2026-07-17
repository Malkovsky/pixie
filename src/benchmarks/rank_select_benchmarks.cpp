#include <benchmark/benchmark.h>
#include <pixie/rank_select/implementations.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>
#include <pixie/storage/implementations.h>

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
using RankSelect = pixie::RankSelectSupport<>;
constexpr RankSelect::SelectSupport kMeasuredSelectSupport =
    RankSelect::SelectSupport::kBoth;

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
                         std::size_t auxiliary_bytes,
                         std::uint64_t repetition_index) {
  const double input_bytes = static_cast<double>((size + 7) / 8);
  const double auxiliary_bytes_as_double = static_cast<double>(auxiliary_bytes);
  state.counters["N"] = static_cast<double>(size);
  state.counters["one_fill_percent"] =
      fill_spec(fill).expected_one_fill_percent;
  state.counters["input_bytes"] = input_bytes;
  state.counters["aux_bytes"] = auxiliary_bytes_as_double;
  state.counters["aux_mib"] = auxiliary_bytes_as_double / (1024.0 * 1024.0);
  state.counters["aux_bits_per_input_bit"] =
      size == 0 ? 0.0 : 8.0 * auxiliary_bytes_as_double / size;
  state.counters["select1_enabled"] = 1.0;
  state.counters["select0_enabled"] = 1.0;
  state.counters["seed_repetition"] = static_cast<double>(repetition_index);
}

void set_fnbp_counters(benchmark::State& state,
                       std::size_t size,
                       std::size_t auxiliary_bytes,
                       std::uint64_t repetition_index) {
  set_common_counters(state, size, Fill::k50, auxiliary_bytes,
                      repetition_index);
  state.counters["fnbp_source"] = 1.0;
}

template <Fill fill>
void run_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, fill);
  const BitDataset dataset(size, fill, seeds.source_seed);
  std::size_t auxiliary_bytes = 0;

  for (auto _ : state) {
    RankSelect support(dataset.source(), dataset.size(),
                       kMeasuredSelectSupport);
    auxiliary_bytes = support.memory_usage_bytes();
    benchmark::DoNotOptimize(auxiliary_bytes);
    benchmark::ClobberMemory();
  }

  set_common_counters(state, size, fill, auxiliary_bytes,
                      seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

template <Fill fill, QueryOperation operation>
void run_query(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, fill);
  const BitDataset dataset(size, fill, seeds.source_seed);
  const RankSelect support(dataset.source(), dataset.size(),
                           kMeasuredSelectSupport);
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

  set_common_counters(state, size, fill, support.memory_usage_bytes(),
                      seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
}

void run_fnbp_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, Fill::k50);
  const FnbpDataset dataset(size, seeds.source_seed);
  std::size_t auxiliary_bytes = 0;

  for (auto _ : state) {
    RankSelect support(dataset.source(), dataset.size(),
                       kMeasuredSelectSupport);
    auxiliary_bytes = support.memory_usage_bytes();
    benchmark::DoNotOptimize(auxiliary_bytes);
    benchmark::ClobberMemory();
  }

  set_fnbp_counters(state, size, auxiliary_bytes, seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

template <QueryOperation operation>
void run_fnbp_query(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const SeedContext seeds = make_seed_context(state, size, Fill::k50);
  const FnbpDataset dataset(size, seeds.source_seed);
  const RankSelect support(dataset.source(), dataset.size(),
                           kMeasuredSelectSupport);
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

  set_fnbp_counters(state, size, support.memory_usage_bytes(),
                    seeds.repetition_index);
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()));
}

template <Fill fill>
void register_build_row() {
  const std::string name =
      "rank_select_build_both_" + std::string(fill_spec(fill).name);
  auto* row = benchmark::RegisterBenchmark(name.c_str(), &run_build<fill>);
  for (const std::size_t size : kSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

template <Fill fill, QueryOperation operation>
void register_query_row(std::string_view operation_name) {
  const std::string name = "rank_select_" + std::string(operation_name) + "_" +
                           std::string(fill_spec(fill).name);
  auto* row =
      benchmark::RegisterBenchmark(name.c_str(), &run_query<fill, operation>);
  for (const std::size_t size : kSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kNanosecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

template <Fill fill>
void register_fill_rows() {
  register_build_row<fill>();
  register_query_row<fill, QueryOperation::kRank1>("rank1");
  register_query_row<fill, QueryOperation::kRank0>("rank0");
  register_query_row<fill, QueryOperation::kSelect1>("select1");
  register_query_row<fill, QueryOperation::kSelect0>("select0");
}

void register_fnbp_build_row() {
  auto* row = benchmark::RegisterBenchmark("rank_select_fnbp_build_both",
                                           &run_fnbp_build);
  for (const std::size_t size : kFnbpSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kMillisecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

template <QueryOperation operation>
void register_fnbp_query_row(std::string_view operation_name) {
  const std::string name = "rank_select_fnbp_" + std::string(operation_name);
  auto* row =
      benchmark::RegisterBenchmark(name.c_str(), &run_fnbp_query<operation>);
  for (const std::size_t size : kFnbpSizes) {
    row->Arg(static_cast<std::int64_t>(size));
  }
  row->ArgNames({"N"})
      ->Unit(benchmark::kNanosecond)
      ->MinWarmUpTime(kBenchmarkWarmupSeconds)
      ->MinTime(kBenchmarkMinSeconds);
}

void register_fnbp_rows() {
  register_fnbp_build_row();
  register_fnbp_query_row<QueryOperation::kRank1>("rank1");
  register_fnbp_query_row<QueryOperation::kRank0>("rank0");
  register_fnbp_query_row<QueryOperation::kSelect1>("select1");
  register_fnbp_query_row<QueryOperation::kSelect0>("select0");
}

void register_benchmarks() {
  register_fill_rows<Fill::k12p5>();
  register_fill_rows<Fill::k50>();
  register_fill_rows<Fill::k87p5>();
  register_fnbp_rows();
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
