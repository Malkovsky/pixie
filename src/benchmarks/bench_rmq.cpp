#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/rmq.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kQueryCount = 32768;
constexpr std::size_t kBpBlockSize = 128;
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

struct BpQueryShape {
  std::size_t same_block = 0;
  std::size_t cross_block = 0;
  std::size_t middle_block = 0;
  std::size_t disjoint_boundary = 0;
  std::size_t fused_boundary = 0;
  std::size_t excess_calls = 0;
  std::size_t left_boundary_width = 0;
  std::size_t right_boundary_width = 0;
  std::size_t same_block_width = 0;
};

struct BlockRange {
  std::size_t block = 0;
  std::size_t left = 0;
  std::size_t right = 0;
};

struct MacroRange {
  std::size_t left = 0;
  std::size_t right = 0;
};

struct alignas(16) BenchBlockSummary {
  static constexpr std::uint64_t kSigned48Mask = (std::uint64_t{1} << 48) - 1;
  static constexpr std::uint64_t kSigned48SignBit = std::uint64_t{1} << 47;

  std::uint64_t word0 = 0;
  std::uint64_t word1 = 0;

  static BenchBlockSummary make(std::int64_t base_depth,
                                std::int64_t min_value,
                                std::size_t min_offset) {
    const std::uint64_t packed_base = pack_signed48(base_depth);
    const std::uint64_t ordered_min = pack_ordered48(min_value);
    return {ordered_min | (static_cast<std::uint64_t>(min_offset) << 48),
            packed_base};
  }

  std::int64_t min_value() const {
    return unpack_ordered48(ordered_min_value());
  }

  std::uint64_t ordered_min_value() const { return word0 & kSigned48Mask; }

 private:
  static std::uint64_t pack_signed48(std::int64_t value) {
    return static_cast<std::uint64_t>(value) & kSigned48Mask;
  }

  static std::uint64_t pack_ordered48(std::int64_t value) {
    return pack_signed48(value) ^ kSigned48SignBit;
  }

  static std::int64_t unpack_signed48(std::uint64_t value) {
    if ((value & kSigned48SignBit) != 0) {
      value |= ~kSigned48Mask;
    }
    return static_cast<std::int64_t>(value);
  }

  static std::int64_t unpack_ordered48(std::uint64_t value) {
    return unpack_signed48(value ^ kSigned48SignBit);
  }
};

static_assert(sizeof(BenchBlockSummary) == 16);
static_assert(alignof(BenchBlockSummary) == 16);

struct BenchBlockSummaryMinLess {
  bool operator()(const BenchBlockSummary& left,
                  const BenchBlockSummary& right) const {
    return left.ordered_min_value() < right.ordered_min_value();
  }
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

std::size_t bp_block_size(const DepthDataset& dataset, std::size_t block) {
  const std::size_t begin = block * kBpBlockSize;
  return std::min(kBpBlockSize, dataset.depths.size() - begin);
}

std::array<std::uint64_t, 2> bp_block_bits_stack(const DepthDataset& dataset,
                                                 std::size_t block) {
  const std::size_t first_word = block * (kBpBlockSize / 64);
  const std::uint64_t lo =
      first_word < dataset.bits.size() ? dataset.bits[first_word] : 0;
  const std::uint64_t hi =
      first_word + 1 < dataset.bits.size() ? dataset.bits[first_word + 1] : 0;
  return {lo, hi};
}

const std::uint64_t* bp_block_bits_direct(const DepthDataset& dataset,
                                          std::size_t block) {
  return dataset.bits.data() + block * (kBpBlockSize / 64);
}

BpQueryShape compute_bp_query_shape(const DepthDataset& dataset) {
  BpQueryShape shape;
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block == right_block) {
      ++shape.same_block;
      ++shape.excess_calls;
      shape.same_block_width += right - left + 1;
      continue;
    }

    ++shape.cross_block;
    const std::size_t left_offset = left % kBpBlockSize;
    const std::size_t right_offset = right % kBpBlockSize;
    const std::size_t left_width = kBpBlockSize - left_offset;
    const std::size_t right_width = right_offset + 1;
    if (left_offset > right_offset) {
      ++shape.disjoint_boundary;
      if ((kBpBlockSize - 1 - left_offset) + right_offset >= 32) {
        ++shape.fused_boundary;
      }
    }
    shape.excess_calls += 2;
    shape.left_boundary_width += left_width;
    shape.right_boundary_width += right_width;
    if (left_block + 1 < right_block) {
      ++shape.middle_block;
    }
  }
  return shape;
}

void set_depth_counters(benchmark::State& state,
                        const DepthDataset& dataset,
                        bool include_shape) {
  state.counters["N"] = static_cast<double>(dataset.size);
  state.counters["max_width"] = static_cast<double>(dataset.max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));

  if (!include_shape) {
    return;
  }

  const BpQueryShape shape = compute_bp_query_shape(dataset);
  const double query_count = static_cast<double>(dataset.ranges.size());
  const double cross_count = static_cast<double>(shape.cross_block);
  const double same_count = static_cast<double>(shape.same_block);
  state.counters["same_block_ratio"] =
      static_cast<double>(shape.same_block) / query_count;
  state.counters["cross_block_ratio"] =
      static_cast<double>(shape.cross_block) / query_count;
  state.counters["middle_block_ratio"] =
      static_cast<double>(shape.middle_block) / query_count;
  state.counters["disjoint_boundary_ratio"] =
      static_cast<double>(shape.disjoint_boundary) / query_count;
  state.counters["fused_boundary_eligible_ratio"] =
      static_cast<double>(shape.fused_boundary) / query_count;
  state.counters["excess_calls_per_query"] =
      static_cast<double>(shape.excess_calls) / query_count;
  state.counters["avg_same_width"] =
      same_count == 0
          ? 0.0
          : static_cast<double>(shape.same_block_width) / same_count;
  state.counters["avg_left_boundary_width"] =
      cross_count == 0
          ? 0.0
          : static_cast<double>(shape.left_boundary_width) / cross_count;
  state.counters["avg_right_boundary_width"] =
      cross_count == 0
          ? 0.0
          : static_cast<double>(shape.right_boundary_width) / cross_count;
}

std::vector<BlockRange> make_same_block_ranges(const DepthDataset& dataset) {
  std::vector<BlockRange> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block == right_block) {
      out.push_back({left_block, left % kBpBlockSize, right % kBpBlockSize});
    }
  }
  return out;
}

std::vector<BlockRange> make_left_boundary_ranges(const DepthDataset& dataset) {
  std::vector<BlockRange> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block != right_block) {
      out.push_back({left_block, left % kBpBlockSize,
                     bp_block_size(dataset, left_block) - 1});
    }
  }
  return out;
}

std::vector<BlockRange> make_right_boundary_ranges(
    const DepthDataset& dataset) {
  std::vector<BlockRange> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block != right_block) {
      out.push_back({right_block, 0, right % kBpBlockSize});
    }
  }
  return out;
}

std::vector<std::pair<BlockRange, BlockRange>> make_boundary_pairs(
    const DepthDataset& dataset) {
  std::vector<std::pair<BlockRange, BlockRange>> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block != right_block) {
      out.push_back({{left_block, left % kBpBlockSize,
                      bp_block_size(dataset, left_block) - 1},
                     {right_block, 0, right % kBpBlockSize}});
    }
  }
  return out;
}

std::vector<std::pair<BlockRange, BlockRange>> make_disjoint_boundary_pairs(
    const DepthDataset& dataset) {
  std::vector<std::pair<BlockRange, BlockRange>> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    const std::size_t left_offset = left % kBpBlockSize;
    const std::size_t right_offset = right % kBpBlockSize;
    if (left_block != right_block && left_offset > right_offset) {
      out.push_back(
          {{left_block, left_offset, bp_block_size(dataset, left_block) - 1},
           {right_block, 0, right_offset}});
    }
  }
  return out;
}

std::vector<MacroRange> make_macro_ranges(const DepthDataset& dataset) {
  std::vector<MacroRange> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = right / kBpBlockSize;
    if (left_block + 1 < right_block) {
      out.push_back({left_block + 1, right_block - 1});
    }
  }
  return out;
}

std::vector<BenchBlockSummary> make_block_summaries(
    const DepthDataset& dataset) {
  const std::size_t block_count =
      (dataset.depths.size() + kBpBlockSize - 1) / kBpBlockSize;
  std::vector<BenchBlockSummary> summaries;
  summaries.reserve(block_count);
  for (std::size_t block = 0; block < block_count; ++block) {
    const std::size_t begin = block * kBpBlockSize;
    const std::size_t end =
        std::min(begin + kBpBlockSize, dataset.depths.size());
    auto min_it = std::min_element(dataset.depths.begin() + begin,
                                   dataset.depths.begin() + end);
    summaries.push_back(BenchBlockSummary::make(
        dataset.depths[begin], *min_it,
        static_cast<std::size_t>(min_it - dataset.depths.begin() - begin)));
  }
  return summaries;
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

void run_bp_boundary_stack(
    benchmark::State& state,
    std::vector<BlockRange> (*make_ranges)(const DepthDataset&)) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const std::vector<BlockRange> ranges = make_ranges(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no diagnostic ranges for this size/width");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const BlockRange range = ranges[query_index++ % ranges.size()];
    const auto bits = bp_block_bits_stack(dataset, range.block);
    ExcessResult result = excess_min_128(bits.data(), range.left, range.right);
    benchmark::DoNotOptimize(result.min_excess);
    benchmark::DoNotOptimize(result.offset);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
}

void run_bp_boundary_pair_stack(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const auto ranges = make_boundary_pairs(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no cross-block boundary pairs for this size/width");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] = ranges[query_index++ % ranges.size()];
    const auto left_bits = bp_block_bits_stack(dataset, left.block);
    const auto right_bits = bp_block_bits_stack(dataset, right.block);
    ExcessResult left_result =
        excess_min_128(left_bits.data(), left.left, left.right);
    ExcessResult right_result =
        excess_min_128(right_bits.data(), right.left, right.right);
    benchmark::DoNotOptimize(left_result.min_excess);
    benchmark::DoNotOptimize(left_result.offset);
    benchmark::DoNotOptimize(right_result.min_excess);
    benchmark::DoNotOptimize(right_result.offset);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
}

void run_bp_boundary_pair_direct(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const auto ranges = make_boundary_pairs(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no cross-block boundary pairs for this size/width");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] = ranges[query_index++ % ranges.size()];
    ExcessResult left_result = excess_min_128(
        bp_block_bits_direct(dataset, left.block), left.left, left.right);
    ExcessResult right_result = excess_min_128(
        bp_block_bits_direct(dataset, right.block), right.left, right.right);
    benchmark::DoNotOptimize(left_result.min_excess);
    benchmark::DoNotOptimize(left_result.offset);
    benchmark::DoNotOptimize(right_result.min_excess);
    benchmark::DoNotOptimize(right_result.offset);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
}

void run_bp_boundary_pair_disjoint_direct(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const auto ranges = make_disjoint_boundary_pairs(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no disjoint cross-block boundary pairs");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] = ranges[query_index++ % ranges.size()];
    ExcessResult left_result = excess_min_128(
        bp_block_bits_direct(dataset, left.block), left.left, left.right);
    ExcessResult right_result = excess_min_128(
        bp_block_bits_direct(dataset, right.block), right.left, right.right);
    benchmark::DoNotOptimize(left_result.min_excess);
    benchmark::DoNotOptimize(left_result.offset);
    benchmark::DoNotOptimize(right_result.min_excess);
    benchmark::DoNotOptimize(right_result.offset);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
}

void run_bp_boundary_pair_fused(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const auto ranges = make_disjoint_boundary_pairs(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no disjoint cross-block boundary pairs");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const auto [left, right] = ranges[query_index++ % ranges.size()];
    ExcessBoundaryPairResult result = excess_min_128_disjoint_suffix_prefix(
        bp_block_bits_direct(dataset, left.block), left.left,
        bp_block_bits_direct(dataset, right.block), right.right);
    benchmark::DoNotOptimize(result.suffix.min_excess);
    benchmark::DoNotOptimize(result.suffix.offset);
    benchmark::DoNotOptimize(result.prefix.min_excess);
    benchmark::DoNotOptimize(result.prefix.offset);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
}

void run_bp_macro_only(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const std::vector<BenchBlockSummary> block_summaries =
      make_block_summaries(dataset);
  const pixie::rmq::SparseTable<BenchBlockSummary, BenchBlockSummaryMinLess,
                                Index>
      macro_rmq{std::span<const BenchBlockSummary>(block_summaries)};
  const std::vector<MacroRange> ranges = make_macro_ranges(dataset);
  if (ranges.empty()) {
    state.SkipWithError("no middle-block macro ranges for this size/width");
    return;
  }

  std::size_t query_index = 0;
  for (auto _ : state) {
    const MacroRange range = ranges[query_index++ % ranges.size()];
    std::size_t result = macro_rmq.arg_min(range.left, range.right);
    benchmark::DoNotOptimize(result);
  }

  set_depth_counters(state, dataset, false);
  state.counters["diagnostic_ranges"] = static_cast<double>(ranges.size());
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

  set_depth_counters(state, dataset, true);
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
      benchmark::RegisterBenchmark("rmq_bp_diag_same_block_boundary",
                                   [](benchmark::State& state) {
                                     run_bp_boundary_stack(
                                         state, make_same_block_ranges);
                                   })
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_left_boundary",
                                   [](benchmark::State& state) {
                                     run_bp_boundary_stack(
                                         state, make_left_boundary_ranges);
                                   })
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_right_boundary",
                                   [](benchmark::State& state) {
                                     run_bp_boundary_stack(
                                         state, make_right_boundary_ranges);
                                   })
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_boundary_pair_stack",
                                   run_bp_boundary_pair_stack)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_boundary_pair_direct",
                                   run_bp_boundary_pair_direct)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_boundary_pair_disjoint_direct",
                                   run_bp_boundary_pair_disjoint_direct)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_boundary_pair_fused",
                                   run_bp_boundary_pair_fused)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark("rmq_bp_diag_macro_only", run_bp_macro_only)
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
