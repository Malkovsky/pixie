#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/rmq.h>
#include <pixie/rmq/experimental/cartesian_tree_rmm_btree_rmq.h>
#include <pixie/rmq/experimental/node_euler_btree_rmq.h>
#include <pixie/rmq/node_euler_btree_rmq.h>

#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
#include <pixie/rmq/sdsl_sct_rmq.h>
#endif

#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr std::size_t kQueryPoolBytes = 512 * 1024;
constexpr std::size_t kQueryCount =
    kQueryPoolBytes / sizeof(std::pair<std::size_t, std::size_t>);
constexpr std::size_t kBpBlockSize = 128;
constexpr std::size_t kCacheLineBytes = 64;
constexpr std::size_t kSparseFootprintCycles = 16;
// Value-RMQ rows rotate between several arrays so query/build timings do not
// depend on one accidental global-minimum position. Large rows use fewer
// variants to keep multi-index query benchmarks within practical memory limits.
constexpr std::size_t kDefaultValueDatasetVariants = 4;
constexpr std::size_t kLargeValueDatasetVariants = 2;
constexpr std::size_t kValueDatasetVariantLargeSize = 1ull << 18;
constexpr std::size_t kMaxValueDatasetVariants = 16;
constexpr std::int64_t kValueDatasetGlobalMinimum =
    std::numeric_limits<std::int64_t>::min() / 4;
constexpr const char* kValueDatasetVariantsEnv = "PIXIE_BENCH_VALUE_VARIANTS";
using Index = std::size_t;

static_assert(kQueryPoolBytes % sizeof(std::pair<std::size_t, std::size_t>) ==
              0);

struct Dataset {
  std::size_t size = 0;
  std::size_t max_width = 0;
  std::vector<std::int64_t> values;
};

struct DepthDataset {
  std::size_t size = 0;
  std::size_t max_width = 0;
  std::vector<std::int64_t> depths;
  std::vector<std::uint64_t> bits;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
};

struct DeltaBitDataset {
  std::size_t size = 0;
  std::vector<std::uint64_t> bits;
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

class PerfCounterGroup {
 public:
  struct Values {
    std::uint64_t cycles = 0;
    std::uint64_t instructions = 0;
    std::uint64_t cache_misses = 0;
  };

  PerfCounterGroup() {
#ifdef __linux__
    if (!enabled()) {
      return;
    }

    leader_fd_ = open_hardware_event(PERF_COUNT_HW_CPU_CYCLES, -1, true);
    if (leader_fd_ < 0) {
      return;
    }
    instructions_fd_ =
        open_hardware_event(PERF_COUNT_HW_INSTRUCTIONS, leader_fd_, false);
    cache_misses_fd_ =
        open_hardware_event(PERF_COUNT_HW_CACHE_MISSES, leader_fd_, false);
    if (instructions_fd_ < 0 || cache_misses_fd_ < 0) {
      close_all();
      return;
    }

    ioctl(leader_fd_, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
    ioctl(leader_fd_, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    active_ = true;
#endif
  }

  PerfCounterGroup(const PerfCounterGroup&) = delete;
  PerfCounterGroup& operator=(const PerfCounterGroup&) = delete;

  ~PerfCounterGroup() { close_all(); }

  bool stop(Values* values) {
#ifdef __linux__
    if (!active_) {
      return false;
    }
    ioctl(leader_fd_, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
    active_ = false;

    struct ReadGroup {
      std::uint64_t count = 0;
      std::uint64_t values[3] = {};
    };
    ReadGroup group;
    const ssize_t bytes = read(leader_fd_, &group, sizeof(group));
    if (bytes != static_cast<ssize_t>(sizeof(group)) || group.count != 3) {
      close_all();
      return false;
    }

    values->cycles = group.values[0];
    values->instructions = group.values[1];
    values->cache_misses = group.values[2];
    close_all();
    return true;
#else
    (void)values;
    return false;
#endif
  }

 private:
  static bool enabled() {
    const char* value = std::getenv("PIXIE_BENCH_PERF_COUNTERS");
    return value != nullptr && value[0] != '\0' &&
           !(value[0] == '0' && value[1] == '\0');
  }

#ifdef __linux__
  static int open_hardware_event(std::uint64_t config,
                                 int group_fd,
                                 bool disabled) {
    perf_event_attr attr = {};
    attr.type = PERF_TYPE_HARDWARE;
    attr.size = sizeof(attr);
    attr.config = config;
    attr.disabled = disabled ? 1 : 0;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.read_format = PERF_FORMAT_GROUP;
    return static_cast<int>(
        syscall(__NR_perf_event_open, &attr, 0, -1, group_fd, 0));
  }

  void close_fd(int* fd) {
    if (*fd >= 0) {
      close(*fd);
      *fd = -1;
    }
  }

  void close_all() {
    if (active_ && leader_fd_ >= 0) {
      ioctl(leader_fd_, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
      active_ = false;
    }
    close_fd(&cache_misses_fd_);
    close_fd(&instructions_fd_);
    close_fd(&leader_fd_);
  }

  int leader_fd_ = -1;
  int instructions_fd_ = -1;
  int cache_misses_fd_ = -1;
  bool active_ = false;
#else
  void close_all() {}
#endif
};

void set_perf_counters(benchmark::State& state,
                       const PerfCounterGroup::Values& values) {
  const double iterations = static_cast<double>(state.iterations());
  if (iterations == 0.0) {
    return;
  }
  state.counters["cycles_per_query"] =
      static_cast<double>(values.cycles) / iterations;
  state.counters["instructions_per_query"] =
      static_cast<double>(values.instructions) / iterations;
  state.counters["ipc"] = values.cycles == 0
                              ? 0.0
                              : static_cast<double>(values.instructions) /
                                    static_cast<double>(values.cycles);
  state.counters["cache_misses_per_query"] =
      static_cast<double>(values.cache_misses) / iterations;
}

Dataset make_dataset(std::size_t size, std::size_t max_width) {
  Dataset dataset;
  dataset.size = size;
  dataset.max_width = max_width;
  dataset.values.resize(size);

  std::mt19937_64 rng(kSeed ^ (size * 0x9E3779B185EBCA87ull) ^
                      (max_width * 0xBF58476D1CE4E5B9ull));
  std::uniform_int_distribution<std::int64_t> value_dist(-1'000'000, 1'000'000);
  std::generate(dataset.values.begin(), dataset.values.end(),
                [&] { return value_dist(rng); });
  return dataset;
}

std::size_t value_dataset_variant_count(std::size_t size) {
  if (const char* value = std::getenv(kValueDatasetVariantsEnv);
      value != nullptr && value[0] != '\0') {
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value, &end, 10);
    if (end != value && *end == '\0' && parsed != 0) {
      return std::min<std::size_t>(static_cast<std::size_t>(parsed),
                                   kMaxValueDatasetVariants);
    }
  }
  return size > kValueDatasetVariantLargeSize ? kLargeValueDatasetVariants
                                              : kDefaultValueDatasetVariants;
}

std::size_t value_dataset_global_min_position(std::size_t size,
                                              std::size_t variant,
                                              std::size_t variant_count) {
  return ((2 * variant + 1) * size) / (2 * variant_count);
}

std::vector<Dataset> make_value_dataset_variants(std::size_t size,
                                                 std::size_t max_width,
                                                 std::size_t variant_count) {
  variant_count = std::max<std::size_t>(
      1, std::min<std::size_t>(variant_count, std::max<std::size_t>(size, 1)));

  std::vector<std::int64_t> base_values(size);
  std::mt19937_64 rng(kSeed ^ (size * 0x9E3779B185EBCA87ull));
  std::uniform_int_distribution<std::int64_t> value_dist(-1'000'000, 1'000'000);
  std::generate(base_values.begin(), base_values.end(),
                [&] { return value_dist(rng); });

  std::vector<Dataset> variants;
  variants.reserve(variant_count);
  for (std::size_t variant = 0; variant < variant_count; ++variant) {
    Dataset dataset;
    dataset.size = size;
    dataset.max_width = max_width;
    dataset.values = base_values;
    if (!dataset.values.empty()) {
      dataset.values[value_dataset_global_min_position(
          size, variant, variant_count)] = kValueDatasetGlobalMinimum;
    }
    variants.push_back(std::move(dataset));
  }
  return variants;
}

DeltaBitDataset make_delta_bit_dataset(std::size_t size) {
  DeltaBitDataset dataset;
  dataset.size = size;
  if (size == 0) {
    return dataset;
  }
  dataset.bits.assign((size - 1 + 63) / 64, 0);

  std::mt19937_64 rng(kSeed ^ (size * 0xD6E8FEB86659FD93ull));
  for (std::size_t i = 0; i + 1 < size; ++i) {
    if ((rng() & 1u) != 0) {
      dataset.bits[i >> 6] |= std::uint64_t{1} << (i & 63);
    }
  }
  return dataset;
}

struct SparseTableFootprintStats {
  std::size_t query_count = 0;
  std::size_t max_observed_width = 0;
  std::size_t over_max_width = 0;
  std::uint64_t total_width = 0;
  std::uint64_t total_level = 0;
  std::uint64_t same_table_line = 0;
  std::uint64_t same_value_line = 0;
  std::array<std::uint64_t, 64> level_counts = {};
  std::vector<std::uint64_t> table_lines;
  std::vector<std::uint64_t> value_lines;

  void finalize() {
    sort_unique(table_lines);
    sort_unique(value_lines);
  }

  double avg_width() const {
    return query_count == 0 ? 0.0
                            : static_cast<double>(total_width) / query_count;
  }

  double avg_level() const {
    return query_count == 0 ? 0.0
                            : static_cast<double>(total_level) / query_count;
  }

  double ratio(std::uint64_t count) const {
    return query_count == 0 ? 0.0 : static_cast<double>(count) / query_count;
  }

  std::size_t percentile_level(std::uint64_t numerator,
                               std::uint64_t denominator) const {
    if (query_count == 0) {
      return 0;
    }
    const std::uint64_t target =
        (static_cast<std::uint64_t>(query_count) * numerator + denominator -
         1) /
        denominator;
    std::uint64_t cumulative = 0;
    for (std::size_t level = 0; level < level_counts.size(); ++level) {
      cumulative += level_counts[level];
      if (cumulative >= target) {
        return level;
      }
    }
    return level_counts.size() - 1;
  }

 private:
  static void sort_unique(std::vector<std::uint64_t>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
  }
};

SparseTableFootprintStats compute_sparse_table_footprint(
    const pixie::rmq::SparseTable<std::int64_t, std::less<std::int64_t>, Index>&
        rmq,
    std::size_t size,
    std::size_t max_width) {
  constexpr std::size_t kIndexEntriesPerLine = kCacheLineBytes / sizeof(Index);
  constexpr std::size_t kValueEntriesPerLine =
      kCacheLineBytes / sizeof(std::int64_t);

  SparseTableFootprintStats stats;
  stats.query_count = kSparseFootprintCycles * kQueryCount;
  stats.table_lines.reserve(stats.query_count * 2);
  stats.value_lines.reserve(stats.query_count * 2);

  RandomQueryGenerator queries(size, max_width,
                               kSeed ^ (size * 0xC2B2AE3D27D4EB4Full) ^
                                   (max_width * 0x165667B19E3779F9ull));

  for (std::size_t i = 0; i < stats.query_count; ++i) {
    const auto [left, right] = queries.get_random_query();
    const std::size_t width = right - left;
    const std::size_t level = std::bit_width(width) - 1;
    const std::size_t span = std::size_t{1} << level;
    const std::size_t first_table_position = left;
    const std::size_t second_table_position = right - span;
    const std::uint64_t first_table_line =
        (static_cast<std::uint64_t>(level) << 56) |
        (first_table_position / kIndexEntriesPerLine);
    const std::uint64_t second_table_line =
        (static_cast<std::uint64_t>(level) << 56) |
        (second_table_position / kIndexEntriesPerLine);

    const std::size_t first_value_position =
        rmq.arg_min(first_table_position, first_table_position + span);
    const std::size_t second_value_position =
        rmq.arg_min(second_table_position, second_table_position + span);
    const std::uint64_t first_value_line =
        first_value_position / kValueEntriesPerLine;
    const std::uint64_t second_value_line =
        second_value_position / kValueEntriesPerLine;

    stats.total_width += width;
    stats.max_observed_width = std::max(stats.max_observed_width, width);
    stats.over_max_width += width > max_width ? 1 : 0;
    stats.total_level += level;
    ++stats.level_counts[level];
    stats.same_table_line += first_table_line == second_table_line ? 1 : 0;
    stats.same_value_line += first_value_line == second_value_line ? 1 : 0;
    stats.table_lines.push_back(first_table_line);
    stats.table_lines.push_back(second_table_line);
    stats.value_lines.push_back(first_value_line);
    stats.value_lines.push_back(second_value_line);
  }

  stats.finalize();
  return stats;
}

void set_sparse_table_footprint_counters(
    benchmark::State& state,
    const SparseTableFootprintStats& stats) {
  const double query_count = static_cast<double>(stats.query_count);
  state.counters["sample_queries"] = query_count;
  state.counters["sample_cycles"] = static_cast<double>(kSparseFootprintCycles);
  state.counters["avg_width"] = stats.avg_width();
  state.counters["max_observed_width"] =
      static_cast<double>(stats.max_observed_width);
  state.counters["over_max_width_ratio"] = stats.ratio(stats.over_max_width);
  state.counters["avg_level"] = stats.avg_level();
  state.counters["level_p50"] =
      static_cast<double>(stats.percentile_level(50, 100));
  state.counters["level_p90"] =
      static_cast<double>(stats.percentile_level(90, 100));
  state.counters["level_p99"] =
      static_cast<double>(stats.percentile_level(99, 100));
  state.counters["table_distinct_lines_per_query"] =
      2.0 - stats.ratio(stats.same_table_line);
  state.counters["value_distinct_lines_per_query"] =
      2.0 - stats.ratio(stats.same_value_line);
  state.counters["table_working_set_lines_per_query"] =
      static_cast<double>(stats.table_lines.size()) / query_count;
  state.counters["value_working_set_lines_per_query"] =
      static_cast<double>(stats.value_lines.size()) / query_count;
  state.counters["total_working_set_lines_per_query"] =
      static_cast<double>(stats.table_lines.size() + stats.value_lines.size()) /
      query_count;

  for (std::size_t level = 0; level < stats.level_counts.size(); ++level) {
    if (stats.level_counts[level] == 0) {
      continue;
    }
    state.counters["level_" + std::to_string(level) + "_ratio"] =
        stats.ratio(stats.level_counts[level]);
  }
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

  RandomQueryGenerator queries(size, max_width,
                               kSeed ^ (size * 0xA24BAED4963EE407ull) ^
                                   (max_width * 0x9FB21C651E98DF25ull));
  for (auto& range : dataset.ranges) {
    range = queries.get_random_query();
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

BpQueryShape compute_bp_query_shape(const DepthDataset& dataset,
                                    std::size_t block_size) {
  BpQueryShape shape;
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t last = right - 1;
    const std::size_t left_block = left / block_size;
    const std::size_t right_block = last / block_size;
    if (left_block == right_block) {
      ++shape.same_block;
      ++shape.excess_calls;
      shape.same_block_width += right - left;
      continue;
    }

    ++shape.cross_block;
    const std::size_t left_offset = left % block_size;
    const std::size_t right_offset = last % block_size;
    const std::size_t left_width = block_size - left_offset;
    const std::size_t right_width = right_offset + 1;
    if (left_offset > right_offset) {
      ++shape.disjoint_boundary;
      if (left_offset < block_size - 1 && right_offset > 0) {
        ++shape.fused_boundary;
        ++shape.excess_calls;
      } else {
        shape.excess_calls += 2;
      }
    } else {
      shape.excess_calls += 2;
    }
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
                        bool include_shape,
                        std::size_t block_size = kBpBlockSize) {
  state.counters["N"] = static_cast<double>(dataset.size);
  state.counters["max_width"] = static_cast<double>(dataset.max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));

  if (!include_shape) {
    return;
  }

  const BpQueryShape shape = compute_bp_query_shape(dataset, block_size);
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
  state.counters["fused_boundary_ratio"] =
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
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
    if (left_block == right_block) {
      out.push_back({left_block, left % kBpBlockSize, last % kBpBlockSize});
    }
  }
  return out;
}

std::vector<BlockRange> make_left_boundary_ranges(const DepthDataset& dataset) {
  std::vector<BlockRange> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
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
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
    if (left_block != right_block) {
      out.push_back({right_block, 0, last % kBpBlockSize});
    }
  }
  return out;
}

std::vector<std::pair<BlockRange, BlockRange>> make_boundary_pairs(
    const DepthDataset& dataset) {
  std::vector<std::pair<BlockRange, BlockRange>> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
    if (left_block != right_block) {
      out.push_back({{left_block, left % kBpBlockSize,
                      bp_block_size(dataset, left_block) - 1},
                     {right_block, 0, last % kBpBlockSize}});
    }
  }
  return out;
}

std::vector<std::pair<BlockRange, BlockRange>> make_disjoint_boundary_pairs(
    const DepthDataset& dataset) {
  std::vector<std::pair<BlockRange, BlockRange>> out;
  out.reserve(dataset.ranges.size());
  for (const auto [left, right] : dataset.ranges) {
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
    const std::size_t left_offset = left % kBpBlockSize;
    const std::size_t right_offset = last % kBpBlockSize;
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
    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBpBlockSize;
    const std::size_t right_block = last / kBpBlockSize;
    if (left_block + 1 < right_block) {
      out.push_back({left_block + 1, right_block});
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
  const std::vector<Dataset> datasets = make_value_dataset_variants(
      size, max_width, value_dataset_variant_count(size));
  std::vector<Rmq> rmqs;
  rmqs.reserve(datasets.size());
  for (const Dataset& dataset : datasets) {
    rmqs.emplace_back(std::span<const std::int64_t>(dataset.values));
  }

  const std::uint64_t query_seed = kSeed ^ (size * 0xC2B2AE3D27D4EB4Full) ^
                                   (max_width * 0x165667B19E3779F9ull);
  std::vector<RandomQueryGenerator> queries;
  queries.reserve(datasets.size());
  for (std::size_t i = 0; i < datasets.size(); ++i) {
    queries.emplace_back(size, max_width, query_seed);
  }

  std::size_t variant = 0;
  PerfCounterGroup perf_counters;
  for (auto _ : state) {
    const auto [left, right] = queries[variant].get_random_query();
    std::size_t result = rmqs[variant].arg_min(left, right);
    benchmark::DoNotOptimize(result);
    ++variant;
    if (variant == rmqs.size()) {
      variant = 0;
    }
  }
  PerfCounterGroup::Values perf_values;
  if (perf_counters.stop(&perf_values)) {
    set_perf_counters(state, perf_values);
  }

  state.counters["N"] = static_cast<double>(size);
  state.counters["max_width"] = static_cast<double>(max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
  state.counters["value_arrays"] = static_cast<double>(datasets.size());
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

template <class Rmq>
void run_value_rmq_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::vector<Dataset> datasets = make_value_dataset_variants(
      size, size, value_dataset_variant_count(size));

  std::size_t variant = 0;
  for (auto _ : state) {
    const Dataset& dataset = datasets[variant];
    Rmq rmq(std::span<const std::int64_t>(dataset.values));
    std::size_t built_size = rmq.size();
    benchmark::DoNotOptimize(built_size);
    benchmark::ClobberMemory();
    ++variant;
    if (variant == datasets.size()) {
      variant = 0;
    }
  }

  set_build_counters(state, size,
                     datasets.front().values.size() * sizeof(std::int64_t));
  state.counters["value_arrays"] = static_cast<double>(datasets.size());
}

template <class Rmq>
void run_bp_rmq_build(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const DeltaBitDataset dataset = make_delta_bit_dataset(size);

  for (auto _ : state) {
    Rmq rmq(std::span<const std::uint64_t>(dataset.bits), dataset.size);
    std::size_t built_size = rmq.size();
    benchmark::DoNotOptimize(built_size);
    benchmark::ClobberMemory();
  }

  set_build_counters(state, size, dataset.bits.size() * sizeof(std::uint64_t));
}

void run_sparse_table_footprint(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const Dataset dataset = make_dataset(size, max_width);
  const pixie::rmq::SparseTable<std::int64_t, std::less<std::int64_t>, Index>
      rmq(std::span<const std::int64_t>(dataset.values));
  const SparseTableFootprintStats stats =
      compute_sparse_table_footprint(rmq, size, max_width);

  for (auto _ : state) {
    std::size_t query_count = stats.query_count;
    benchmark::DoNotOptimize(query_count);
  }

  state.counters["N"] = static_cast<double>(size);
  state.counters["max_width"] = static_cast<double>(max_width);
  state.counters["index_bytes"] = static_cast<double>(sizeof(Index));
  set_sparse_table_footprint_counters(state, stats);
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

template <class Rmq, std::size_t BlockSize = kBpBlockSize>
void run_depth_queries(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const std::size_t max_width = static_cast<std::size_t>(state.range(1));
  const DepthDataset dataset = make_depth_dataset(size, max_width);
  const Rmq rmq(std::span<const std::uint64_t>(dataset.bits),
                dataset.depths.size());
  RandomQueryGenerator queries(size, max_width,
                               kSeed ^ (size * 0xD1B54A32D192ED03ull) ^
                                   (max_width * 0x94D049BB133111EBull));

  PerfCounterGroup perf_counters;
  for (auto _ : state) {
    const auto [left, right] = queries.get_random_query();
    std::size_t result = rmq.arg_min(left, right);
    benchmark::DoNotOptimize(result);
  }
  PerfCounterGroup::Values perf_values;
  if (perf_counters.stop(&perf_values)) {
    set_perf_counters(state, perf_values);
  }

  set_depth_counters(state, dataset, true, BlockSize);
}

void register_benchmarks() {
  const std::vector<std::size_t> sizes = {1ull << 10, 1ull << 14, 1ull << 18,
                                          1ull << 22, 1ull << 24, 1ull << 26};
  const std::vector<std::size_t> build_sizes = {
      1ull << 10, 1ull << 14, 1ull << 18, 1ull << 22, 1ull << 24};
  const std::vector<std::size_t> widths = {64, 4096, 1ull << 18, 1ull << 22,
                                           1ull << 26};
  // Pure sparse-table rows materialize O(n log n) indexes; above 2^22 they
  // dominate memory/runtime and make large benchmark passes noisy.
  constexpr std::size_t kSparseTableBenchmarkMaxSize = 1ull << 22;

  for (const std::size_t size : build_sizes) {
    if (size <= kSparseTableBenchmarkMaxSize) {
      benchmark::RegisterBenchmark(
          "rmq_build_sparse_table",
          &run_value_rmq_build<pixie::rmq::SparseTable<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Arg(static_cast<std::int64_t>(size))
          ->Unit(benchmark::kMillisecond);
    }
    benchmark::RegisterBenchmark(
        "rmq_build_segment_tree",
        &run_value_rmq_build<pixie::rmq::SegmentTree<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_cartesian_tree",
        &run_value_rmq_build<pixie::rmq::CartesianTreeRmq<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_cartesian_tree_rmm_btree",
        &run_value_rmq_build<pixie::rmq::experimental::CartesianTreeRmMBTreeRmq<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
    benchmark::RegisterBenchmark(
        "rmq_build_sdsl_sct",
        &run_value_rmq_build<pixie::rmq::SdslSctRmq<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
#endif
    benchmark::RegisterBenchmark(
        "rmq_build_node_euler_btree",
        &run_value_rmq_build<pixie::rmq::NodeEulerBTreeRmq<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_node_euler_btree_exp",
        &run_value_rmq_build<pixie::rmq::experimental::NodeEulerBTreeRmq<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_node_euler_btree_exp_mask_leaf",
        &run_value_rmq_build<pixie::rmq::experimental::NodeEulerBTreeRmq<
            std::int64_t, std::less<std::int64_t>, Index, 248, 256,
            pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_segment_btree_xl",
        &run_value_rmq_build<pixie::rmq::SegmentBTreeXL<
            std::int64_t, std::less<std::int64_t>, Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_bp_plus_minus_one",
        &run_bp_rmq_build<pixie::rmq::BpPlusMinusOneRmq<Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_bp_rmm_btree",
        &run_bp_rmq_build<
            pixie::rmq::experimental::RmMBTreePlusMinusOneRmq<Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
    benchmark::RegisterBenchmark(
        "rmq_build_bp_one_interval_btree",
        &run_bp_rmq_build<pixie::rmq::OneIntervalBTreeRmq<Index>>)
        ->Arg(static_cast<std::int64_t>(size))
        ->Unit(benchmark::kMillisecond);
  }

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
      if (size <= kSparseTableBenchmarkMaxSize) {
        if (std::getenv("PIXIE_BENCH_SPARSE_FOOTPRINT") != nullptr) {
          benchmark::RegisterBenchmark("rmq_sparse_table_footprint",
                                       run_sparse_table_footprint)
              ->Args({static_cast<std::int64_t>(size),
                      static_cast<std::int64_t>(width)})
              ->Iterations(1)
              ->Unit(benchmark::kNanosecond);
        }
        benchmark::RegisterBenchmark(
            "rmq_sparse_table",
            &run_queries<pixie::rmq::SparseTable<
                std::int64_t, std::less<std::int64_t>, Index>>)
            ->Args({static_cast<std::int64_t>(size),
                    static_cast<std::int64_t>(width)})
            ->Unit(benchmark::kNanosecond);
        benchmark::RegisterBenchmark(
            "rmq_sparse_table_unaligned",
            &run_queries<pixie::rmq::SparseTable<
                std::int64_t, std::less<std::int64_t>, Index, 0>>)
            ->Args({static_cast<std::int64_t>(size),
                    static_cast<std::int64_t>(width)})
            ->Unit(benchmark::kNanosecond);
      }
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
          "rmq_cartesian_tree_rmm_btree",
          &run_queries<pixie::rmq::experimental::CartesianTreeRmMBTreeRmq<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
#ifdef PIXIE_THIRD_PARTY_BENCHMARKS
      benchmark::RegisterBenchmark(
          "rmq_sdsl_sct",
          &run_queries<pixie::rmq::SdslSctRmq<std::int64_t,
                                              std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
#endif
      benchmark::RegisterBenchmark(
          "rmq_node_euler_btree",
          &run_queries<pixie::rmq::NodeEulerBTreeRmq<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_node_euler_btree_exp",
          &run_queries<pixie::rmq::experimental::NodeEulerBTreeRmq<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_node_euler_btree_exp_mask_leaf",
          &run_queries<pixie::rmq::experimental::NodeEulerBTreeRmq<
              std::int64_t, std::less<std::int64_t>, Index, 248, 256,
              pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_segment_btree_xl",
          &run_queries<pixie::rmq::SegmentBTreeXL<
              std::int64_t, std::less<std::int64_t>, Index>>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_plus_minus_one",
          &run_depth_queries<pixie::rmq::BpPlusMinusOneRmq<Index>, 128>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_rmm_btree",
          &run_depth_queries<
              pixie::rmq::experimental::RmMBTreePlusMinusOneRmq<Index>,
              pixie::rmq::experimental::RmMBTreePlusMinusOneRmq<
                  Index>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmq<Index>,
              pixie::rmq::OneIntervalBTreeRmq<Index>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_h1",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmq<Index, 1>,
              pixie::rmq::OneIntervalBTreeRmq<Index, 1>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_h2",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmq<Index, 2>,
              pixie::rmq::OneIntervalBTreeRmq<Index, 2>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_h4",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmq<Index, 4>,
              pixie::rmq::OneIntervalBTreeRmq<Index, 4>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_records",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index>,
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<
                  Index>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_records_h1",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index, 1>,
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index,
                                                             1>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_records_h2",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index, 2>,
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index,
                                                             2>::kBlockSize>)
          ->Args({static_cast<std::int64_t>(size),
                  static_cast<std::int64_t>(width)})
          ->Unit(benchmark::kNanosecond);
      benchmark::RegisterBenchmark(
          "rmq_bp_one_interval_btree_records_h4",
          &run_depth_queries<
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index, 4>,
              pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<Index,
                                                             4>::kBlockSize>)
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
