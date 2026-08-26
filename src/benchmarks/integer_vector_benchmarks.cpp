#include <benchmark/benchmark.h>
#include <pixie/integer_vector/implementations.h>
#include <pixie/serialization.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kQueryCount = 1 << 16;
constexpr double kWarmupSeconds = 0.1;
constexpr double kMinimumSeconds = 0.5;

enum class Dataset : std::int64_t {
  kAllZero,
  kFixedWidth,
  kSkewed,
  kDenseMonotone,
  kDuplicateMonotone,
  kSparseMonotone,
};

constexpr std::array kDatasets = {
    Dataset::kAllZero,
    Dataset::kFixedWidth,
    Dataset::kSkewed,
    Dataset::kDenseMonotone,
    Dataset::kDuplicateMonotone,
    Dataset::kSparseMonotone,
};

constexpr std::array kMonotoneDatasets = {
    Dataset::kDenseMonotone,
    Dataset::kDuplicateMonotone,
    Dataset::kSparseMonotone,
};

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

const char* dataset_name(Dataset dataset) {
  switch (dataset) {
    case Dataset::kAllZero:
      return "all_zero";
    case Dataset::kFixedWidth:
      return "fixed_width";
    case Dataset::kSkewed:
      return "skewed";
    case Dataset::kDenseMonotone:
      return "dense_monotone";
    case Dataset::kDuplicateMonotone:
      return "duplicate_monotone";
    case Dataset::kSparseMonotone:
      return "sparse_monotone";
  }
  return "unknown";
}

std::vector<std::uint64_t> make_values(std::size_t size, Dataset dataset) {
  std::vector<std::uint64_t> values(size);
  std::uint64_t sparse_value = 0;
  for (std::size_t i = 0; i < size; ++i) {
    const std::uint64_t random = splitmix64(42 + i);
    switch (dataset) {
      case Dataset::kAllZero:
        values[i] = 0;
        break;
      case Dataset::kFixedWidth:
        values[i] = random & ((std::uint64_t{1} << 20) - 1);
        break;
      case Dataset::kSkewed:
        values[i] =
            i % 32 == 0 ? random & ((std::uint64_t{1} << 20) - 1) : random & 15;
        break;
      case Dataset::kDenseMonotone:
        values[i] = i;
        break;
      case Dataset::kDuplicateMonotone:
        values[i] = i / 16;
        break;
      case Dataset::kSparseMonotone:
        sparse_value += 1 + (random & 1023);
        values[i] = sparse_value;
        break;
    }
  }
  return values;
}

template <class Serializable>
std::vector<std::byte> serialize_to_vector(const Serializable& value) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  value.serialize(writer);
  writer.finish();
  return sink.take();
}

class AlignedArtifact {
 public:
  explicit AlignedArtifact(std::span<const std::byte> bytes)
      : words_((bytes.size() + 7) / 8), size_(bytes.size()) {
    std::ranges::copy(bytes, writable_bytes().begin());
  }

  std::span<const std::byte> bytes() const {
    return std::as_bytes(std::span<const std::uint64_t>(words_)).first(size_);
  }

 private:
  std::span<std::byte> writable_bytes() {
    return std::as_writable_bytes(std::span<std::uint64_t>(words_))
        .first(size_);
  }

  std::vector<std::uint64_t> words_;
  std::size_t size_;
};

Dataset selected_dataset(const benchmark::State& state) {
  return static_cast<Dataset>(state.range(1));
}

std::vector<std::size_t> make_positions(std::size_t size) {
  std::vector<std::size_t> positions(kQueryCount);
  for (std::size_t i = 0; i < positions.size(); ++i) {
    positions[i] = static_cast<std::size_t>(splitmix64(117 + i) % size);
  }
  return positions;
}

void set_counters(benchmark::State& state,
                  std::size_t size,
                  std::size_t width,
                  std::size_t artifact_bytes,
                  std::size_t allocated_bytes) {
  state.counters["N"] = static_cast<double>(size);
  state.counters["width"] = static_cast<double>(width);
  state.counters["artifact_bytes"] = static_cast<double>(artifact_bytes);
  state.counters["allocated_bytes"] = static_cast<double>(allocated_bytes);
  state.counters["allocated_bits_per_value"] =
      size == 0 ? 0.0 : static_cast<double>(allocated_bytes) * 8.0 / size;
}

void BM_PackedConstruction(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const pixie::PackedIntegerVector<> sample(values);
  const std::size_t artifact_bytes = serialize_to_vector(sample).size();
  for (auto _ : state) {
    pixie::PackedIntegerVector<> vector(values);
    benchmark::DoNotOptimize(vector);
    benchmark::ClobberMemory();
  }
  set_counters(state, size, sample.width(), artifact_bytes,
               sample.memory_usage_bytes() - sizeof(sample));
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(size));
}

void BM_PackedViewConstruction(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const pixie::PackedIntegerVector<> owner(values);
  const auto serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(serialized);
  for (auto _ : state) {
    pixie::BinaryReader reader(artifact.bytes());
    auto view = pixie::PackedIntegerVectorView<>::deserialize(reader);
    benchmark::DoNotOptimize(view);
  }
  set_counters(state, size, owner.width(), artifact.bytes().size(), 0);
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(size));
}

template <bool View>
void BM_PackedAccess(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const pixie::PackedIntegerVector<> owner(values);
  const auto serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(serialized);
  pixie::BinaryReader reader(artifact.bytes());
  const auto view = pixie::PackedIntegerVectorView<>::deserialize(reader);
  const auto positions = make_positions(size);
  std::size_t query = 0;
  for (auto _ : state) {
    const std::size_t position = positions[query++ & (kQueryCount - 1)];
    if constexpr (View) {
      benchmark::DoNotOptimize(view[position]);
    } else {
      benchmark::DoNotOptimize(owner[position]);
    }
  }
  set_counters(state, size, owner.width(), artifact.bytes().size(),
               View ? 0 : owner.memory_usage_bytes() - sizeof(owner));
  state.SetItemsProcessed(state.iterations());
}

void BM_RawSpanAccess(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const std::span<const std::uint64_t> span(values);
  const auto positions = make_positions(size);
  std::size_t query = 0;
  for (auto _ : state) {
    std::uint64_t value = span[positions[query++ & (kQueryCount - 1)]];
    benchmark::DoNotOptimize(value);
  }
  const std::size_t width =
      values.empty() ? 0 : std::bit_width(*std::ranges::max_element(values));
  set_counters(state, size, width, values.size() * sizeof(std::uint64_t),
               values.capacity() * sizeof(std::uint64_t));
  state.SetItemsProcessed(state.iterations());
}

template <bool View>
void BM_PackedCopy(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const pixie::PackedIntegerVector<> owner(values);
  const auto serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(serialized);
  pixie::BinaryReader reader(artifact.bytes());
  const auto view = pixie::PackedIntegerVectorView<>::deserialize(reader);
  constexpr std::size_t kCopyCount = 256;
  std::vector<std::uint64_t> output(kCopyCount);
  const auto positions = make_positions(size - kCopyCount + 1);
  std::size_t query = 0;
  for (auto _ : state) {
    const std::size_t begin = positions[query++ & (kQueryCount - 1)];
    if constexpr (View) {
      view.copy_to(begin, output);
    } else {
      owner.copy_to(begin, output);
    }
    benchmark::ClobberMemory();
  }
  set_counters(state, size, owner.width(), artifact.bytes().size(),
               View ? 0 : owner.memory_usage_bytes() - sizeof(owner));
  state.SetItemsProcessed(state.iterations() * kCopyCount);
}

template <bool View, bool Upper>
void BM_MonotoneBound(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const auto values = make_values(size, selected_dataset(state));
  const pixie::PackedMonotoneIntegerVector<> owner{
      pixie::PackedIntegerVector<>(values)};
  const auto serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(serialized);
  pixie::BinaryReader reader(artifact.bytes());
  const auto view =
      pixie::PackedMonotoneIntegerVectorView<>::deserialize(reader);
  std::vector<std::uint64_t> queries(kQueryCount);
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const std::uint64_t value = values[splitmix64(i + 991) % size];
    queries[i] =
        i % 2 == 0 || value == std::numeric_limits<std::uint64_t>::max()
            ? value
            : value + 1;
  }
  std::size_t query = 0;
  for (auto _ : state) {
    const std::uint64_t value = queries[query++ & (kQueryCount - 1)];
    if constexpr (View && Upper) {
      benchmark::DoNotOptimize(view.upper_bound_index(value));
    } else if constexpr (View) {
      benchmark::DoNotOptimize(view.lower_bound_index(value));
    } else if constexpr (Upper) {
      benchmark::DoNotOptimize(owner.upper_bound_index(value));
    } else {
      benchmark::DoNotOptimize(owner.lower_bound_index(value));
    }
  }
  const std::size_t width = values.empty() ? 0 : std::bit_width(values.back());
  set_counters(state, size, width, artifact.bytes().size(),
               View ? 0 : owner.memory_usage_bytes() - sizeof(owner));
  state.SetItemsProcessed(state.iterations());
}

void configure(benchmark::internal::Benchmark* row,
               std::size_t size,
               Dataset dataset,
               benchmark::TimeUnit unit) {
  row->Args(
         {static_cast<std::int64_t>(size), static_cast<std::int64_t>(dataset)})
      ->ArgNames({"N", "dataset"})
      ->Unit(unit)
      ->MinWarmUpTime(kWarmupSeconds)
      ->MinTime(kMinimumSeconds);
}

void register_benchmarks() {
  constexpr std::array<std::size_t, 3> general_sizes = {1 << 10, 1 << 16,
                                                        1 << 20};
  for (const std::size_t size : general_sizes) {
    for (const Dataset dataset : kDatasets) {
      const std::string suffix = std::string("/") + dataset_name(dataset);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/construct_owner" + suffix).c_str(),
                    &BM_PackedConstruction),
                size, dataset, benchmark::kMillisecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/construct_view" + suffix).c_str(),
                    &BM_PackedViewConstruction),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/access_owner" + suffix).c_str(),
                    &BM_PackedAccess<false>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/access_view" + suffix).c_str(),
                    &BM_PackedAccess<true>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/access_raw_span" + suffix).c_str(),
                    &BM_RawSpanAccess),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/copy_owner" + suffix).c_str(),
                    &BM_PackedCopy<false>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/copy_view" + suffix).c_str(),
                    &BM_PackedCopy<true>),
                size, dataset, benchmark::kNanosecond);
    }
  }

  constexpr std::array<std::size_t, 5> bound_sizes = {
      std::size_t{1} << 10, std::size_t{1} << 14, std::size_t{1} << 18,
      std::size_t{1} << 22, std::size_t{1} << 26};
  for (const std::size_t size : bound_sizes) {
    for (const Dataset dataset : kMonotoneDatasets) {
      const std::string suffix = std::string("/") + dataset_name(dataset);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/lower_bound_owner" + suffix).c_str(),
                    &BM_MonotoneBound<false, false>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/lower_bound_view" + suffix).c_str(),
                    &BM_MonotoneBound<true, false>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/upper_bound_owner" + suffix).c_str(),
                    &BM_MonotoneBound<false, true>),
                size, dataset, benchmark::kNanosecond);
      configure(benchmark::RegisterBenchmark(
                    ("integer_vector/upper_bound_view" + suffix).c_str(),
                    &BM_MonotoneBound<true, true>),
                size, dataset, benchmark::kNanosecond);
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
