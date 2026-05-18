#pragma once

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <regex>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pixie_bench {

struct RmMBenchmarkArgs {
  int min_exp = 14;
  int max_exp = 22;
  std::size_t Q = 200000;
  double p1 = 0.5;
  std::uint64_t seed = 42;
  std::size_t block_bits = 64;
  int per_octave = 0;
  std::vector<std::size_t> explicit_sizes;
  std::set<std::string> ops;
  std::string benchmark_filter;
};

struct RmMBenchmarkPools {
  std::vector<std::size_t> inds_any;
  std::vector<std::size_t> rank10_end_positions;
  std::vector<std::size_t> open_positions_zero_based;
  std::vector<std::size_t> open_positions_one_based;
  std::vector<std::size_t> close_positions_one_based;
  std::vector<std::size_t> inds;
  std::vector<std::size_t> inds_1N;
  std::vector<int> deltas;
  std::vector<std::pair<std::size_t, std::size_t>> segs;

  std::vector<std::size_t> ks1;
  std::vector<std::size_t> ks0;
  std::vector<std::size_t> ks10;
  std::vector<std::size_t> minselect_q;
};

inline const std::vector<std::string>& AllRmMOps() {
  static const std::vector<std::string> ops = {
      "rank1",
      "rank0",
      "select1",
      "select0",
      "rank10",
      "select10",
      "excess",
      "fwdsearch",
      "bwdsearch",
      "range_min_query_pos",
      "range_min_query_val",
      "range_max_query_pos",
      "range_max_query_val",
      "mincount",
      "minselect",
      "close",
      "open",
      "enclose",
  };
  return ops;
}

inline std::vector<std::string> SplitCsv(const std::string& text) {
  std::vector<std::string> out;
  std::size_t position = 0;
  while (position < text.size()) {
    while (position < text.size() &&
           (text[position] == ',' ||
            std::isspace(static_cast<unsigned char>(text[position])))) {
      ++position;
    }
    const std::size_t start = position;
    while (position < text.size() && text[position] != ',') {
      ++position;
    }
    std::string value = text.substr(start, position - start);
    while (!value.empty() &&
           std::isspace(static_cast<unsigned char>(value.back()))) {
      value.pop_back();
    }
    if (!value.empty()) {
      out.push_back(value);
    }
  }
  return out;
}

inline std::string RandomDyckBits(std::size_t N,
                                  double p1,
                                  std::mt19937_64& rng) {
  const std::size_t pairs = N >> 1;
  const std::size_t L = pairs << 1;
  std::string s;
  s.resize(L);
  if (L == 0) {
    return s;
  }
  p1 = std::clamp(p1, 0., 1.);
  std::size_t opens_left = pairs;
  std::size_t closes_left = pairs;
  int balance = 0;
  std::uniform_real_distribution<double> U(0., 1.);
  for (std::size_t i = 0; i < L; ++i) {
    if (opens_left == 0) {
      s[i] = '0';
      --closes_left;
      --balance;
      continue;
    }

    if (closes_left == 0) {
      s[i] = '1';
      --opens_left;
      ++balance;
      continue;
    }

    if (balance == 0 || U(rng) < p1) {
      s[i] = '1';
      --opens_left;
      ++balance;
    } else {
      s[i] = '0';
      --closes_left;
      --balance;
    }
  }
  return s;
}

inline std::size_t Count10(const std::string& s) {
  std::size_t count = 0;
  if (s.size() < 2) {
    return 0;
  }
  for (std::size_t i = 0; i + 1 < s.size(); ++i) {
    if (s[i] == '1' && s[i + 1] == '0') {
      ++count;
    }
  }
  return count;
}

inline std::vector<std::uint64_t> PackWordsLsbFirst(const std::string& bits) {
  std::vector<std::uint64_t> words((bits.size() + 63) / 64, 0);
  for (std::size_t i = 0; i < bits.size(); ++i) {
    if (bits[i] == '1') {
      words[i >> 6] |= (std::uint64_t(1) << (i & 63));
    }
  }
  return words;
}

template <class Tree>
struct RmMBenchmarkTraits {
  static bool SupportsOp(std::string_view) { return true; }
};

template <class Tree>
struct RmMBenchmarkDataset {
  std::size_t N{};
  std::string bits;
  std::vector<std::uint64_t> words;
  Tree tree;

  std::size_t cnt1{};
  std::size_t cnt0{};
  std::size_t cnt10{};
  RmMBenchmarkPools pool;
};

template <class Tree>
class RmMBenchmark {
 public:
  int Run(int argc, char** argv) {
    benchmark::MaybeReenterWithoutASLR(argc, argv);
    ParseArgsAndStrip(argc, argv);
    AddDefaultBenchmarkArgs(argc, argv);

    benchmark::Initialize(&argc, argv);
    RegisterAll();
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
  }

  const RmMBenchmarkArgs& args() const { return args_; }

 protected:
  bool OpEnabled(std::string_view op) const {
    const std::string op_string(op);
    if (!args_.ops.empty()) {
      return args_.ops.contains(op_string);
    }
    if (!args_.benchmark_filter.empty() && args_.benchmark_filter != "all" &&
        args_.benchmark_filter != ".*") {
      try {
        return std::regex_search(op_string, std::regex(args_.benchmark_filter));
      } catch (const std::regex_error&) {
        return true;
      }
    }
    return true;
  }

 private:
  using Dataset = RmMBenchmarkDataset<Tree>;

  bool SupportsOp(std::string_view op) const {
    return RmMBenchmarkTraits<Tree>::SupportsOp(op);
  }

  bool ActiveOp(std::string_view op) const {
    return OpEnabled(op) && SupportsOp(op);
  }

  void ParseArgsAndStrip(int& argc, char**& argv) {
    auto is_value_token = [&](int candidate_index) -> bool {
      if (candidate_index >= argc) {
        return false;
      }
      std::string candidate = argv[candidate_index];
      return !(candidate.size() >= 2 && candidate[0] == '-' &&
               candidate[1] == '-');
    };

    auto get_value = [&](const std::string& key) -> std::optional<std::string> {
      const std::string prefix = "--" + key;
      const std::string prefix_with_equals = prefix + "=";
      for (int argument_index = 1; argument_index < argc; ++argument_index) {
        std::string argument = argv[argument_index];
        if (argument.rfind(prefix_with_equals, 0) == 0) {
          return argument.substr(prefix_with_equals.size());
        }
        if (argument == prefix && is_value_token(argument_index + 1)) {
          return std::string(argv[argument_index + 1]);
        }
      }
      return std::nullopt;
    };

    auto strip = [&](const std::string& key) {
      const std::string prefix = "--" + key;
      const std::string prefix_with_equals = prefix + "=";
      int write_index = 0;
      for (int argument_index = 0; argument_index < argc; ++argument_index) {
        std::string argument = argv[argument_index];
        if (argument_index > 0 && argument.rfind(prefix_with_equals, 0) == 0) {
          continue;
        }
        if (argument_index > 0 && argument == prefix &&
            is_value_token(argument_index + 1)) {
          ++argument_index;
          continue;
        }
        argv[write_index++] = argv[argument_index];
      }
      argc = write_index;
    };

    if (auto argument_value = get_value("min_exp")) {
      args_.min_exp = std::stoi(*argument_value);
      strip("min_exp");
    }
    if (auto argument_value = get_value("max_exp")) {
      args_.max_exp = std::stoi(*argument_value);
      strip("max_exp");
    }
    if (auto argument_value = get_value("Q")) {
      args_.Q = std::stoull(*argument_value);
      strip("Q");
    }
    if (auto argument_value = get_value("p1")) {
      args_.p1 = std::stod(*argument_value);
      strip("p1");
    }
    if (auto argument_value = get_value("seed")) {
      args_.seed = std::stoull(*argument_value);
      strip("seed");
    }
    if (auto argument_value = get_value("block_bits")) {
      args_.block_bits = std::stoull(*argument_value);
      strip("block_bits");
    }
    if (auto argument_value = get_value("per_octave")) {
      args_.per_octave = std::stoi(*argument_value);
      strip("per_octave");
    }
    if (auto argument_value = get_value("ops")) {
      for (const auto& op : SplitCsv(*argument_value)) {
        args_.ops.insert(op);
      }
      strip("ops");
    }
    if (auto argument_value = get_value("benchmark_filter")) {
      args_.benchmark_filter = *argument_value;
    }

    if (auto argument_value = get_value("explicit_sizes")) {
      std::string sizes_text = *argument_value;
      strip("explicit_sizes");
      std::size_t position = 0;
      while (position < sizes_text.size()) {
        while (
            position < sizes_text.size() &&
            (sizes_text[position] == ',' ||
             std::isspace(static_cast<unsigned char>(sizes_text[position])))) {
          ++position;
        }
        const std::size_t start_index = position;
        while (position < sizes_text.size() &&
               std::isdigit(static_cast<unsigned char>(sizes_text[position]))) {
          ++position;
        }
        if (start_index < position) {
          args_.explicit_sizes.push_back(std::stoull(
              sizes_text.substr(start_index, position - start_index)));
        }
      }
      std::sort(args_.explicit_sizes.begin(), args_.explicit_sizes.end());
      args_.explicit_sizes.erase(
          std::unique(args_.explicit_sizes.begin(), args_.explicit_sizes.end()),
          args_.explicit_sizes.end());
    }
  }

  void AddDefaultBenchmarkArgs(int& argc, char**& argv) {
    auto has = [&](const char* key) {
      const std::string prefix = std::string(key) + "=";
      for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == key || argument.rfind(prefix, 0) == 0) {
          return true;
        }
      }
      return false;
    };

    if (!has("--benchmark_out_format")) {
      extra_args_.emplace_back("--benchmark_out_format=json");
    }
    if (!has("--benchmark_counters_tabular")) {
      extra_args_.emplace_back("--benchmark_counters_tabular=true");
    }
    if (!has("--benchmark_time_unit")) {
      extra_args_.emplace_back("--benchmark_time_unit=ns");
    }

    argv_vec_.assign(argv, argv + argc);
    for (auto& argument : extra_args_) {
      argv_vec_.push_back(argument.data());
    }
    argv_vec_.push_back(nullptr);
    argc = static_cast<int>(argv_vec_.size()) - 1;
    argv = argv_vec_.data();
  }

  std::vector<std::size_t> BuildSizeGrid() const {
    if (!args_.explicit_sizes.empty()) {
      std::vector<std::size_t> sizes = args_.explicit_sizes;
      sizes.erase(std::remove(sizes.begin(), sizes.end(), 0), sizes.end());
      return sizes;
    }

    const int low_exp = args_.min_exp;
    const int high_exp = args_.max_exp;
    if (args_.per_octave <= 0) {
      std::vector<std::size_t> sizes;
      for (int exp = low_exp; exp <= high_exp; ++exp) {
        sizes.push_back(std::size_t(1) << exp);
      }
      return sizes;
    }

    std::set<std::size_t> candidates;
    const std::size_t min_size = std::size_t(1) << low_exp;
    const std::size_t max_size = std::size_t(1) << high_exp;

    for (int exp = low_exp; exp <= high_exp; ++exp) {
      for (int step = 0; step <= args_.per_octave; ++step) {
        const long double x =
            exp + static_cast<long double>(step) / args_.per_octave;
        std::size_t size =
            static_cast<std::size_t>(std::llround(std::pow(2.0L, x)));
        size = std::clamp(size, min_size, max_size);
        if (size != 0) {
          candidates.insert(size);
        }
      }
    }

    std::vector<std::size_t> sizes(candidates.begin(), candidates.end());
    if (sizes.front() != min_size) {
      sizes.insert(sizes.begin(), min_size);
    }
    if (sizes.back() != max_size) {
      sizes.push_back(max_size);
    }
    return sizes;
  }

  std::shared_ptr<Dataset> BuildDataset(std::size_t N) {
    std::mt19937_64 rng(args_.seed ^
                        static_cast<std::uint64_t>(N) * 0x9E3779B185EBCA87ull);

    auto data = std::make_shared<Dataset>();
    data->bits = RandomDyckBits(N, args_.p1, rng);
    data->N = data->bits.size();
    N = data->N;

    data->words = PackWordsLsbFirst(data->bits);
    data->tree =
        Tree(std::span<const std::uint64_t>(data->words), N, args_.block_bits);
    FillPools(*data, rng);
    return data;
  }

  void FillPools(Dataset& data, std::mt19937_64& rng) {
    const std::size_t N = data.N;
    const bool need_rank10 = ActiveOp("rank10");
    const bool need_select1 = ActiveOp("select1");
    const bool need_select0 = ActiveOp("select0");
    const bool need_select10 = ActiveOp("select10");
    const bool need_fwdsearch = ActiveOp("fwdsearch");
    const bool need_bwdsearch = ActiveOp("bwdsearch");
    const bool need_range_queries =
        ActiveOp("range_min_query_pos") || ActiveOp("range_min_query_val") ||
        ActiveOp("range_max_query_pos") || ActiveOp("range_max_query_val") ||
        ActiveOp("mincount") || ActiveOp("minselect");
    const bool need_minselect = ActiveOp("minselect");
    const bool need_open_positions = ActiveOp("close") || ActiveOp("enclose");
    const bool need_close_positions = ActiveOp("open");

    if (need_select1 || need_select0 || ActiveOp("rank0")) {
      data.cnt1 = std::count(data.bits.begin(), data.bits.end(), '1');
      data.cnt0 = N - data.cnt1;
    }
    if (need_select10) {
      data.cnt10 = Count10(data.bits);
    }

    std::vector<std::size_t> open_positions_zero_based;
    std::vector<std::size_t> open_positions_one_based;
    std::vector<std::size_t> close_positions_one_based;
    if (need_open_positions || need_close_positions) {
      open_positions_zero_based.reserve(N >> 1);
      open_positions_one_based.reserve(N >> 1);
      close_positions_one_based.reserve(N >> 1);
      for (std::size_t i = 0; i < N; ++i) {
        if (data.bits[i] == '1') {
          open_positions_zero_based.push_back(i);
          open_positions_one_based.push_back(i + 1);
        } else {
          close_positions_one_based.push_back(i + 1);
        }
      }
    }

    const std::size_t sample_count =
        std::max<std::size_t>(1, std::min<std::size_t>(args_.Q, 32768));

    std::uniform_int_distribution<std::size_t> nonedge_index_distribution(
        N > 1 ? 1 : 0, N > 1 ? (N - 1) : 0);
    std::uniform_int_distribution<std::size_t> ind_dist(0, N ? (N - 1) : 0);
    std::uniform_int_distribution<std::size_t> ind_dist_1N(1, N ? N : 1);
    std::uniform_int_distribution<std::size_t> rank10_end_position_distribution(
        N >= 2 ? 2 : 0, N >= 2 ? N : 0);
    std::uniform_int_distribution<int> delta_distribution(-8, +8);

    if (ActiveOp("rank1") || ActiveOp("rank0") || ActiveOp("excess")) {
      data.pool.inds_any.resize(sample_count);
    }
    if (need_rank10) {
      data.pool.rank10_end_positions.resize(sample_count);
    }
    if (need_fwdsearch) {
      data.pool.inds.resize(sample_count);
    }
    if (need_bwdsearch) {
      data.pool.inds_1N.resize(sample_count);
    }
    if (need_fwdsearch || need_bwdsearch) {
      data.pool.deltas.resize(sample_count);
    }
    if (need_range_queries) {
      data.pool.segs.resize(sample_count);
    }

    auto rand_ij = [&]() -> std::pair<std::size_t, std::size_t> {
      if (N == 0) {
        return {0, 0};
      }
      std::size_t left = ind_dist(rng);
      std::size_t right = ind_dist(rng);
      if (left > right) {
        std::swap(left, right);
      }
      return {left, right};
    };

    for (std::size_t i = 0; i < sample_count; ++i) {
      if (!data.pool.inds_any.empty()) {
        data.pool.inds_any[i] = nonedge_index_distribution(rng);
      }
      if (!data.pool.rank10_end_positions.empty()) {
        data.pool.rank10_end_positions[i] =
            rank10_end_position_distribution(rng);
      }
      if (!data.pool.inds.empty()) {
        data.pool.inds[i] = (N ? ind_dist(rng) : 0);
      }
      if (!data.pool.inds_1N.empty()) {
        data.pool.inds_1N[i] = (N ? ind_dist_1N(rng) : 0);
      }
      if (!data.pool.deltas.empty()) {
        data.pool.deltas[i] = delta_distribution(rng);
      }
      if (!data.pool.segs.empty()) {
        data.pool.segs[i] = rand_ij();
      }
    }

    auto fill_from_candidates = [&](const std::vector<std::size_t>& candidates,
                                    std::vector<std::size_t>& destination,
                                    std::size_t fallback_value) {
      destination.resize(sample_count);
      if (candidates.empty()) {
        std::fill(destination.begin(), destination.end(), fallback_value);
        return;
      }
      std::uniform_int_distribution<std::size_t> index_distribution(
          0, candidates.size() - 1);
      for (std::size_t sample_index = 0; sample_index < sample_count;
           ++sample_index) {
        destination[sample_index] = candidates[index_distribution(rng)];
      }
    };

    const std::size_t one_based_fallback = (N > 0 ? 1 : 0);
    if (ActiveOp("close")) {
      fill_from_candidates(open_positions_zero_based,
                           data.pool.open_positions_zero_based, 0);
    }
    if (ActiveOp("enclose")) {
      fill_from_candidates(open_positions_one_based,
                           data.pool.open_positions_one_based,
                           one_based_fallback);
    }
    if (ActiveOp("open")) {
      fill_from_candidates(close_positions_one_based,
                           data.pool.close_positions_one_based,
                           one_based_fallback);
    }

    auto fill_ks = [&](std::size_t total, std::vector<std::size_t>& out) {
      out.resize(sample_count);
      if (total == 0) {
        std::fill(out.begin(), out.end(), 1);
        return;
      }
      std::uniform_int_distribution<std::size_t> dist(1, total);
      for (std::size_t i = 0; i < sample_count; ++i) {
        out[i] = dist(rng);
      }
    };

    if (need_select1) {
      fill_ks(data.cnt1, data.pool.ks1);
    }
    if (need_select0) {
      fill_ks(data.cnt0, data.pool.ks0);
    }
    if (need_select10) {
      fill_ks(data.cnt10, data.pool.ks10);
    }

    if (need_minselect) {
      data.pool.minselect_q.resize(sample_count);
      for (std::size_t i = 0; i < sample_count; ++i) {
        auto [left, right] = data.pool.segs[i];
        std::size_t count = data.tree.mincount(left, right);
        if (count == 0) {
          count = 1;
        }
        std::uniform_int_distribution<std::size_t> query_distribution(1, count);
        data.pool.minselect_q[i] = query_distribution(rng);
      }
    }
  }

  void RegisterAll() {
    const auto sizes = BuildSizeGrid();
    bool any_op = false;
    for (const auto& op : AllRmMOps()) {
      any_op = any_op || ActiveOp(op);
    }
    if (!any_op) {
      std::cerr << "No RmM benchmark operations match the requested filter.\n";
      return;
    }

    for (std::size_t size : sizes) {
      auto data = BuildDataset(size);
      keepalive_.push_back(data);

      RegisterOp("rank1", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.rank1(
                       pool.inds_any[sample_index % pool.inds_any.size()]);
                 });
      RegisterOp("rank0", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.rank0(
                       pool.inds_any[sample_index % pool.inds_any.size()]);
                 });
      RegisterOp("select1", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.select1(
                       pool.ks1[sample_index % pool.ks1.size()]);
                 });
      RegisterOp("select0", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.select0(
                       pool.ks0[sample_index % pool.ks0.size()]);
                 });
      RegisterOp(
          "rank10", data, [](const Dataset& dataset, std::size_t sample_index) {
            const auto& pool = dataset.pool;
            return dataset.tree.rank10(
                pool.rank10_end_positions[sample_index %
                                          pool.rank10_end_positions.size()]);
          });
      RegisterOp("select10", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.select10(
                       pool.ks10[sample_index % pool.ks10.size()]);
                 });
      RegisterOp("excess", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.excess(
                       pool.inds_any[sample_index % pool.inds_any.size()]);
                 });
      RegisterOp("fwdsearch", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.fwdsearch(
                       pool.inds[sample_index % pool.inds.size()],
                       pool.deltas[sample_index % pool.deltas.size()]);
                 });
      RegisterOp("bwdsearch", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   return dataset.tree.bwdsearch(
                       pool.inds_1N[sample_index % pool.inds_1N.size()],
                       pool.deltas[sample_index % pool.deltas.size()]);
                 });
      RegisterOp("range_min_query_pos", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const auto [left, right] =
                       pool.segs[sample_index % pool.segs.size()];
                   return dataset.tree.range_min_query_pos(left, right);
                 });
      RegisterOp("range_min_query_val", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const auto [left, right] =
                       pool.segs[sample_index % pool.segs.size()];
                   return dataset.tree.range_min_query_val(left, right);
                 });
      RegisterOp("range_max_query_pos", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const auto [left, right] =
                       pool.segs[sample_index % pool.segs.size()];
                   return dataset.tree.range_max_query_pos(left, right);
                 });
      RegisterOp("range_max_query_val", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const auto [left, right] =
                       pool.segs[sample_index % pool.segs.size()];
                   return dataset.tree.range_max_query_val(left, right);
                 });
      RegisterOp("mincount", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const auto [left, right] =
                       pool.segs[sample_index % pool.segs.size()];
                   return dataset.tree.mincount(left, right);
                 });
      RegisterOp("minselect", data,
                 [](const Dataset& dataset, std::size_t sample_index) {
                   const auto& pool = dataset.pool;
                   const std::size_t index = sample_index % pool.segs.size();
                   const auto [left, right] = pool.segs[index];
                   return dataset.tree.minselect(left, right,
                                                 pool.minselect_q[index]);
                 });
      RegisterOp(
          "close", data, [](const Dataset& dataset, std::size_t sample_index) {
            if (dataset.N == 0) {
              return std::size_t{0};
            }
            const auto& pool = dataset.pool;
            const std::size_t open_position =
                pool.open_positions_zero_based
                    [sample_index % pool.open_positions_zero_based.size()];
            return dataset.tree.close(open_position);
          });
      RegisterOp(
          "open", data, [](const Dataset& dataset, std::size_t sample_index) {
            const auto& pool = dataset.pool;
            const std::size_t close_position_one_based =
                pool.close_positions_one_based
                    [sample_index % pool.close_positions_one_based.size()];
            return dataset.tree.open(close_position_one_based);
          });
      RegisterOp(
          "enclose", data,
          [](const Dataset& dataset, std::size_t sample_index) {
            const auto& pool = dataset.pool;
            const std::size_t open_position_one_based =
                pool.open_positions_one_based
                    [sample_index % pool.open_positions_one_based.size()];
            return dataset.tree.enclose(open_position_one_based);
          });
    }
  }

  template <class Fn>
  void RegisterOp(const std::string& op,
                  std::shared_ptr<Dataset> data,
                  Fn&& body) {
    if (!ActiveOp(op)) {
      return;
    }

    auto idx_ptr = std::make_shared<std::size_t>(0);

    auto* benchmark = benchmark::RegisterBenchmark(
        op.c_str(), [this, data, idx_ptr,
                     body = std::forward<Fn>(body)](benchmark::State& state) {
          const Dataset& dataset = *data;
          for (auto _ : state) {
            const std::size_t sample_index = (*idx_ptr)++;
            auto result = body(dataset, sample_index);
            benchmark::DoNotOptimize(result);
          }
          state.counters["N"] = static_cast<double>(dataset.N);
          state.counters["seed"] = static_cast<double>(args_.seed);
          state.counters["block_bits"] = static_cast<double>(args_.block_bits);
        });

    benchmark->Unit(benchmark::kNanosecond);
  }

  RmMBenchmarkArgs args_;
  std::vector<std::shared_ptr<Dataset>> keepalive_;
  std::vector<std::string> extra_args_;
  std::vector<char*> argv_vec_;
};

}  // namespace pixie_bench
