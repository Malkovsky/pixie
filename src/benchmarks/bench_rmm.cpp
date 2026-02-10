#include <benchmark/benchmark.h>
#include <pixie/rmm_tree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

using pixie::RmMTree;

struct Args {
  int min_exp = 14;
  int max_exp = 22;
  std::size_t Q = 200000;
  double p1 = 0.5;
  std::uint64_t seed = 42;
  std::size_t block_bits = 64;
  int per_octave = 10;
  std::vector<std::size_t> explicit_sizes;
};

static Args args;

static void parse_args_and_strip(int& argc, char**& argv) {
  auto is_value_token = [&](int candidate_index) -> bool {
    if (candidate_index >= argc) {
      return false;
    }
    std::string candidate = argv[candidate_index];
    return !(candidate.size() >= 2 && candidate[0] == '-' &&
             candidate[1] == '-');
  };

  auto get_value = [&](const std::string& key) -> std::optional<std::string> {
    std::string prefix = "--" + key;
    std::string prefix_with_equals = prefix + "=";
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
    std::string prefix = "--" + key;
    std::string prefix_with_equals = prefix + "=";
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
    args.min_exp = std::stoi(*argument_value);
    strip("min_exp");
  }
  if (auto argument_value = get_value("max_exp")) {
    args.max_exp = std::stoi(*argument_value);
    strip("max_exp");
  }
  if (auto argument_value = get_value("Q")) {
    args.Q = std::stoull(*argument_value);
    strip("Q");
  }
  if (auto argument_value = get_value("p1")) {
    args.p1 = std::stod(*argument_value);
    strip("p1");
  }
  if (auto argument_value = get_value("seed")) {
    args.seed = std::stoull(*argument_value);
    strip("seed");
  }
  if (auto argument_value = get_value("block_bits")) {
    args.block_bits = std::stoull(*argument_value);
    strip("block_bits");
  }
  if (auto argument_value = get_value("per_octave")) {
    args.per_octave = std::stoi(*argument_value);
    strip("per_octave");
  }

  if (auto argument_value = get_value("explicit_sizes")) {
    std::string sizes_text = *argument_value;
    strip("explicit_sizes");
    std::size_t position = 0;
    while (position < sizes_text.size()) {
      while (position < sizes_text.size() &&
             (sizes_text[position] == ',' ||
              std::isspace(static_cast<unsigned char>(sizes_text[position])))) {
        ++position;
      }
      std::size_t start_index = position;
      while (position < sizes_text.size() &&
             std::isdigit(static_cast<unsigned char>(sizes_text[position]))) {
        ++position;
      }
      if (start_index < position) {
        args.explicit_sizes.push_back(std::stoull(
            sizes_text.substr(start_index, position - start_index)));
      }
    }
    std::sort(args.explicit_sizes.begin(), args.explicit_sizes.end());
    args.explicit_sizes.erase(
        std::unique(args.explicit_sizes.begin(), args.explicit_sizes.end()),
        args.explicit_sizes.end());
  }
}

static std::string random_dyck_bits(std::size_t N,
                                    double p1,
                                    std::mt19937_64& rng) {
  std::size_t pairs = N >> 1;
  std::size_t L = pairs << 1;
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

    if (balance == 0) {
      s[i] = '1';
      --opens_left;
      ++balance;
      continue;
    }
    if (U(rng) < p1) {
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

static std::vector<std::size_t> build_size_grid() {
  if (!args.explicit_sizes.empty()) {
    std::vector<std::size_t> v = args.explicit_sizes;
    v.erase(std::remove(v.begin(), v.end(), 0), v.end());
    return v;
  }

  const int lo = args.min_exp, hi = args.max_exp;
  if (args.per_octave <= 0) {
    std::vector<std::size_t> v;
    for (int e = lo; e <= hi; ++e) {
      v.push_back(std::size_t(1) << e);
    }
    return v;
  }

  int per_octave_steps = args.per_octave;
  std::set<std::size_t> N_cands;
  std::size_t min_N = std::size_t(1) << lo, max_N = std::size_t(1) << hi;

  for (int e = lo; e <= hi; ++e) {
    for (int t = 0; t <= per_octave_steps; ++t) {
      long double x = e + (long double)t / (long double)per_octave_steps;
      std::size_t N = (std::size_t)std::llround(std::pow(2.0L, x));
      if (N < min_N) {
        N = min_N;
      }
      if (N > max_N) {
        N = max_N;
      }
      if (N) {
        N_cands.insert(N);
      }
    }
  }

  std::vector<std::size_t> v(N_cands.begin(), N_cands.end());
  if (v.front() != min_N) {
    v.insert(v.begin(), min_N);
  }
  if (v.back() != max_N) {
    v.push_back(max_N);
  }
  return v;
}

struct Pools {
  std::vector<std::size_t> inds_any;
  std::vector<std::size_t> rank10_end_positions;
  std::vector<std::size_t> open_positions_zero_based;
  std::vector<std::size_t> open_positions_one_based;
  std::vector<std::size_t> close_positions_one_based;
  std::vector<std::size_t> inds;
  std::vector<std::size_t> inds_1N;
  std::vector<int> deltas;
  std::vector<std::pair<std::size_t, std::size_t>> segs;

  std::vector<std::size_t> ks1, ks0, ks10;
  std::vector<std::size_t> minselect_q;
};

struct Dataset {
  std::size_t N{};
  std::string bits;
  RmMTree t;

  std::size_t cnt1{}, cnt0{}, cnt10{};
  Pools pool;
};

static std::vector<std::shared_ptr<Dataset>> keepalive;

static std::size_t count10(const std::string& s) {
  std::size_t c = 0;
  if (s.size() < 2) {
    return 0;
  }
  for (std::size_t i = 0; i + 1 < s.size(); ++i) {
    if (s[i] == '1' && s[i + 1] == '0') {
      ++c;
    }
  }
  return c;
}

static Dataset build_dataset(std::size_t N) {
  std::mt19937_64 rng(args.seed ^
                      static_cast<std::uint64_t>(N) * 0x9E3779B185EBCA87ull);

  Dataset d;
  d.bits = random_dyck_bits(N, args.p1, rng);
  d.N = d.bits.size();
  N = d.N;

  if (args.block_bits == 0) {
    d.t = RmMTree(d.bits, 0, -1.0f);
  } else {
    d.t = RmMTree(d.bits, args.block_bits);
  }

  d.cnt1 = std::count(d.bits.begin(), d.bits.end(), '1');
  d.cnt0 = N - d.cnt1;
  d.cnt10 = count10(d.bits);

  std::vector<std::size_t> open_positions_zero_based;
  std::vector<std::size_t> open_positions_one_based;
  std::vector<std::size_t> close_positions_one_based;
  open_positions_zero_based.reserve(d.N >> 1);
  open_positions_one_based.reserve(d.N >> 1);
  close_positions_one_based.reserve(d.N >> 1);
  for (std::size_t i = 0; i < N; ++i) {
    if (d.bits[i] == '1') {
      open_positions_zero_based.push_back(i);
      open_positions_one_based.push_back(i + 1);
    } else {
      close_positions_one_based.push_back(i + 1);
    }
  }

  std::size_t sample_count = std::max<std::size_t>(
      1, std::min<std::size_t>(args.Q, std::size_t(32768)));

  std::uniform_int_distribution<std::size_t> nonedge_index_distribution(
      N > 1 ? 1 : 0, N > 1 ? (N - 1) : 0);
  std::uniform_int_distribution<std::size_t> ind_dist(0, N ? (N - 1) : 0);
  std::uniform_int_distribution<std::size_t> ind_dist_1N(1, N ? N : 1);
  std::uniform_int_distribution<std::size_t> rank10_end_position_distribution(
      N >= 2 ? 2 : 0, N >= 2 ? N : 0);
  std::uniform_int_distribution<int> delta_distribution(-8, +8);

  d.pool.inds_any.resize(sample_count);
  d.pool.rank10_end_positions.resize(sample_count);
  d.pool.inds.resize(sample_count);
  d.pool.inds_1N.resize(sample_count);
  d.pool.deltas.resize(sample_count);
  d.pool.segs.resize(sample_count);

  auto rand_ij = [&]() -> std::pair<std::size_t, std::size_t> {
    if (N == 0) {
      return {0, 0};
    }
    std::size_t a = ind_dist(rng), b = ind_dist(rng);
    if (a > b) {
      std::swap(a, b);
    }
    return {a, b};
  };

  for (std::size_t i = 0; i < sample_count; ++i) {
    d.pool.inds_any[i] = nonedge_index_distribution(rng);
    d.pool.rank10_end_positions[i] = rank10_end_position_distribution(rng);
    d.pool.inds[i] = (N ? ind_dist(rng) : 0);
    d.pool.inds_1N[i] = (N ? ind_dist_1N(rng) : 0);
    d.pool.deltas[i] = delta_distribution(rng);
    d.pool.segs[i] = rand_ij();
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
  fill_from_candidates(open_positions_zero_based,
                       d.pool.open_positions_zero_based, 0);
  fill_from_candidates(open_positions_one_based,
                       d.pool.open_positions_one_based, one_based_fallback);
  fill_from_candidates(close_positions_one_based,
                       d.pool.close_positions_one_based, one_based_fallback);

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
  fill_ks(d.cnt1, d.pool.ks1);
  fill_ks(d.cnt0, d.pool.ks0);
  fill_ks(d.cnt10, d.pool.ks10);

  d.pool.minselect_q.resize(sample_count);
  for (std::size_t i = 0; i < sample_count; ++i) {
    auto [l, r] = d.pool.segs[i];
    std::size_t c = d.t.mincount(l, r);
    if (c == 0) {
      c = 1;
    }
    std::uniform_int_distribution<std::size_t> uq(1, c);
    d.pool.minselect_q[i] = uq(rng);
  }

  return d;
}

template <class Fn>
static void register_op(const std::string& op,
                        std::shared_ptr<Dataset> data,
                        Fn&& body) {
  auto idx_ptr = std::make_shared<std::size_t>(0);

  auto* b = benchmark::RegisterBenchmark(
      op.c_str(), [data, idx_ptr, body](benchmark::State& state) {
        const Dataset& D = *data;
        for (auto _ : state) {
          std::size_t i = (*idx_ptr)++;
          body(state, D, i);
        }
        state.counters["N"] = static_cast<double>(D.N);
        state.counters["seed"] = static_cast<double>(args.seed);
        state.counters["block_bits"] = static_cast<double>(args.block_bits);
      });

  b->Unit(benchmark::kNanosecond);
}

static void register_all() {
  auto Ns = build_size_grid();
  for (std::size_t N : Ns) {
    auto data = std::make_shared<Dataset>(build_dataset(N));
    keepalive.push_back(data);

    const auto& P = data->pool;

    register_op("rank1", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.rank1(P.inds_any[k % P.inds_any.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("rank0", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.rank0(P.inds_any[k % P.inds_any.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("select1", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.select1(P.ks1[k % P.ks1.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("select0", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.select0(P.ks0[k % P.ks0.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op(
        "rank10", data,
        [&](benchmark::State&, const Dataset& dataset,
            std::size_t sample_index) {
          benchmark::DoNotOptimize(dataset.t.rank10(
              dataset.pool.rank10_end_positions
                  [sample_index % dataset.pool.rank10_end_positions.size()]));
        });

    register_op("select10", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.select10(P.ks10[k % P.ks10.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("excess", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.excess(P.inds_any[k % P.inds_any.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("fwdsearch", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.fwdsearch(P.inds[k % P.inds.size()],
                                         P.deltas[k % P.deltas.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("bwdsearch", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.bwdsearch(P.inds_1N[k % P.inds_1N.size()],
                                         P.deltas[k % P.deltas.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("range_min_query_pos", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  const std::size_t segment_index = k % P.segs.size();
                  auto [range_begin, range_end] = P.segs[segment_index];
                  auto min_position =
                      D.t.range_min_query_pos(range_begin, range_end);
                  benchmark::DoNotOptimize(min_position);
                });

    register_op("range_min_query_val", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  const std::size_t segment_index = k % P.segs.size();
                  auto [range_begin, range_end] = P.segs[segment_index];
                  auto min_value =
                      D.t.range_min_query_val(range_begin, range_end);
                  benchmark::DoNotOptimize(min_value);
                });

    register_op("range_max_query_pos", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto [i, j] = P.segs[k % P.segs.size()];
                  auto r = D.t.range_max_query_pos(i, j);
                  benchmark::DoNotOptimize(r);
                });

    register_op("range_max_query_val", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto [i, j] = P.segs[k % P.segs.size()];
                  auto r = D.t.range_max_query_val(i, j);
                  benchmark::DoNotOptimize(r);
                });

    register_op("mincount", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto [i, j] = P.segs[k % P.segs.size()];
                  auto r = D.t.mincount(i, j);
                  benchmark::DoNotOptimize(r);
                });

    register_op("minselect", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  const std::size_t idx = k % P.segs.size();
                  auto [i, j] = P.segs[idx];
                  auto q = P.minselect_q[idx];
                  auto r = D.t.minselect(i, j, q);
                  benchmark::DoNotOptimize(r);
                });

    register_op("close", data,
                [&](benchmark::State&, const Dataset& dataset,
                    std::size_t sample_index) {
                  if (dataset.N == 0) {
                    benchmark::DoNotOptimize(0);
                    return;
                  }
                  const std::size_t open_position =
                      dataset.pool.open_positions_zero_based
                          [sample_index %
                           dataset.pool.open_positions_zero_based.size()];
                  benchmark::DoNotOptimize(dataset.t.close(open_position));
                });

    register_op(
        "open", data,
        [&](benchmark::State&, const Dataset& dataset,
            std::size_t sample_index) {
          const std::size_t close_position_one_based =
              dataset.pool.close_positions_one_based
                  [sample_index %
                   dataset.pool.close_positions_one_based.size()];
          benchmark::DoNotOptimize(dataset.t.open(close_position_one_based));
        });

    register_op(
        "enclose", data,
        [&](benchmark::State&, const Dataset& dataset,
            std::size_t sample_index) {
          const std::size_t open_position_one_based =
              dataset.pool.open_positions_one_based
                  [sample_index % dataset.pool.open_positions_one_based.size()];
          benchmark::DoNotOptimize(dataset.t.enclose(open_position_one_based));
        });
  }
}

int main(int argc, char** argv) {
  benchmark::MaybeReenterWithoutASLR(argc, argv);
  parse_args_and_strip(argc, argv);

  auto has = [&](const char* key) {
    std::string p1 = std::string(key) + "=";
    for (int i = 1; i < argc; ++i) {
      std::string s = argv[i];
      if (s == key || s.rfind(p1, 0) == 0) {
        return true;
      }
    }
    return false;
  };

  static std::vector<std::string> extra;
  if (!has("--benchmark_out_format")) {
    extra.emplace_back("--benchmark_out_format=json");
  }
  if (!has("--benchmark_counters_tabular")) {
    extra.emplace_back("--benchmark_counters_tabular=true");
  }
  if (!has("--benchmark_time_unit")) {
    extra.emplace_back("--benchmark_time_unit=ns");
  }

  static std::vector<char*> argv_vec;
  argv_vec.assign(argv, argv + argc);
  for (auto& s : extra) {
    argv_vec.push_back(s.data());
  }
  argv_vec.push_back(nullptr);
  argc = (int)argv_vec.size() - 1;
  argv = argv_vec.data();

  benchmark::Initialize(&argc, argv);
  register_all();
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
