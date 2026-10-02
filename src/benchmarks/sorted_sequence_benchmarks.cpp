#include <benchmark/benchmark.h>
#include <pixie/permutations/sorted_sequence.h>
#include <pixie/permutations/sorted_vector.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kLookups = 128;
constexpr std::size_t kMiB = 1024 * 1024;

std::uint64_t Mix(std::uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

std::size_t EnvironmentMiB(const char* name, std::size_t fallback) {
  const char* text = std::getenv(name);
  if (text == nullptr || *text < '0' || *text > '9') {
    return fallback * kMiB;
  }
  errno = 0;
  char* end = nullptr;
  const auto value = std::strtoull(text, &end, 10);
  if (errno != 0 || *end != '\0' || value == 0 ||
      value > std::numeric_limits<std::size_t>::max() / kMiB) {
    return fallback * kMiB;
  }
  return static_cast<std::size_t>(value) * kMiB;
}

std::size_t AvailableBytes() {
  std::ifstream input("/proc/meminfo");
  std::string key, rest;
  std::size_t kib = 0;
  while (input >> key >> kib) {
    if (key == "MemAvailable:") {
      return kib * 1024;
    }
    std::getline(input, rest);
  }
  return 0;
}

template <class Model>
bool Preflight(benchmark::State& state, std::size_t n) {
  const auto cap = EnvironmentMiB("PIXIE_SORTED_SEQUENCE_RAM_MIB", 4096);
  const auto available = std::min(AvailableBytes(), cap);
  // Conservative for every representation: std::set normally dominates with
  // one separately allocated node per key; Pixie variants need substantially
  // less. The benchmark owns one fixture at a time.
  const long double estimate = 16.0L * kMiB + static_cast<long double>(n) *
                                                  Model::planning_bytes_per_key;
  state.counters["available_budget_bytes"] = available;
  state.counters["planning_bytes"] = static_cast<double>(estimate);
  if (available == 0 || estimate > available / 2) {
    state.SkipWithError(
        "RAM preflight refused row; review PIXIE_SORTED_SEQUENCE_RAM_MIB");
    return false;
  }
  return true;
}

using Tree = pixie::SortedPermutableSequence<std::uint64_t>;
using Vector = pixie::SortedPermutationVector<std::uint64_t>;

struct TreeModel {
  using storage_type = Tree;
  static constexpr long double planning_bytes_per_key = 64;
  static void insert(storage_type& values, std::uint64_t key) {
    values.insert(key);
  }
  static void lower_bound(const storage_type& values, std::uint64_t key) {
    benchmark::DoNotOptimize(values.lower_bound_index(key));
  }
  static std::size_t size(const storage_type& values) { return values.size(); }
  static std::size_t memory_usage(const storage_type& values) {
    return values.memory_usage_bytes();
  }
};

struct VectorModel {
  using storage_type = Vector;
  static constexpr long double planning_bytes_per_key = 64;
  static void insert(storage_type& values, std::uint64_t key) {
    values.insert(key);
  }
  static void lower_bound(const storage_type& values, std::uint64_t key) {
    benchmark::DoNotOptimize(values.lower_bound_index(key));
  }
  static std::size_t size(const storage_type& values) { return values.size(); }
  static std::size_t memory_usage(const storage_type& values) {
    return values.memory_usage_bytes();
  }
};

struct SetModel {
  using storage_type = std::set<std::uint64_t>;
  static constexpr long double planning_bytes_per_key = 96;
  static void insert(storage_type& values, std::uint64_t key) {
    values.insert(key);
  }
  static void lower_bound(const storage_type& values, std::uint64_t key) {
    benchmark::DoNotOptimize(values.lower_bound(key));
  }
  static std::size_t size(const storage_type& values) { return values.size(); }
};

template <class Model>
typename Model::storage_type Build(std::size_t n) {
  typename Model::storage_type result;
  for (std::size_t i = 0; i < n; ++i) {
    Model::insert(result, Mix(i));
  }
  return result;
}

std::array<std::uint64_t, kLookups> LookupKeys(std::size_t n) {
  std::array<std::uint64_t, kLookups> result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    // Alternate exact hits and absent, uniformly distributed probes.
    result[i] = i % 2 == 0 ? Mix(Mix(i + 41) % n) : Mix(n + Mix(i + 91) % n);
  }
  return result;
}

template <class Model>
void Report(benchmark::State& state,
            const typename Model::storage_type& values) {
  state.counters["elements"] = Model::size(values);
  if constexpr (requires { Model::memory_usage(values); }) {
    state.counters["owned_bytes"] = Model::memory_usage(values);
  }
}

template <class Model>
void RandomInsert(benchmark::State& state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  if (!Preflight<Model>(state, n)) {
    return;
  }
  typename Model::storage_type values;
  state.SetLabel("N unique random inserts; lower bound included per insert");
  for (auto _ : state) {
    state.PauseTiming();
    values = typename Model::storage_type{};
    state.ResumeTiming();
    for (std::size_t i = 0; i < n; ++i) {
      Model::insert(values, Mix(i));
    }
    benchmark::DoNotOptimize(values);
    benchmark::ClobberMemory();
  }
  if (Model::size(values) != n) {
    state.SkipWithError("Random insertion lost a unique key");
    return;
  }
  Report<Model>(state, values);
  state.counters["operations_per_iteration"] = n;
  state.SetItemsProcessed(state.iterations() * n);
}

template <class Model>
void LowerBound(benchmark::State& state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  if (!Preflight<Model>(state, n)) {
    return;
  }
  const auto values = Build<Model>(n);
  const auto keys = LookupKeys(n);
  state.SetLabel("128 lower-bound calls; 50% hits; fixture build excluded");
  for (auto _ : state) {
    for (const auto key : keys) {
      Model::lower_bound(values, key);
    }
  }
  Report<Model>(state, values);
  state.counters["operations_per_iteration"] = kLookups;
  state.counters["lookup_hits_percent"] = 50;
  state.SetItemsProcessed(state.iterations() * kLookups);
}

template <class Model>
void MixedInsertLowerBound(benchmark::State& state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  if (!Preflight<Model>(state, n)) {
    return;
  }
  std::vector<std::uint64_t> probes(n);
  for (std::size_t i = 0; i < n; ++i) {
    // The query follows its insertion. Half query an already inserted key and
    // half use a disjoint key, while all insertions stay std::set-compatible.
    probes[i] =
        i % 2 == 0 ? Mix(Mix(i + 101) % (i + 1)) : Mix(n + Mix(i + 151) % n);
  }
  typename Model::storage_type values;
  state.SetLabel("N insert+lower-bound pairs; 50% hits after each insert");
  for (auto _ : state) {
    state.PauseTiming();
    values = typename Model::storage_type{};
    state.ResumeTiming();
    for (std::size_t i = 0; i < n; ++i) {
      Model::insert(values, Mix(i));
      Model::lower_bound(values, probes[i]);
    }
    benchmark::DoNotOptimize(values);
    benchmark::ClobberMemory();
  }
  if (Model::size(values) != n) {
    state.SkipWithError("Mixed insertion lost a unique key");
    return;
  }
  Report<Model>(state, values);
  state.counters["operations_per_iteration"] = n;
  state.counters["lookup_hits_percent"] = 50;
  state.SetItemsProcessed(state.iterations() * n);
}

template <class Model>
void Register(const char* name) {
  const auto add = [=](const char* operation, auto function) {
    const auto row = std::string(operation) + "/" + name;
    benchmark::RegisterBenchmark(row.c_str(), function)
        ->ArgName("N")
        ->Arg(4096)
        ->Arg(65536)
        ->Arg(1048576)
        ->Arg(4194304);
  };
  add("SortedRandomInsert", RandomInsert<Model>);
  add("SortedLowerBound", LowerBound<Model>);
  add("SortedMixedInsertLowerBound", MixedInsertLowerBound<Model>);
}

const bool registered = [] {
  Register<TreeModel>("PermutableSequence64_B256_F8_Prefix");
  Register<VectorModel>("VectorPermutation64_B256_F8_Prefix");
  Register<SetModel>("StdSet64");
  return true;
}();

}  // namespace
