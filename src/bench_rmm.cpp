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
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "rmm_tree.h"

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
  auto getv = [&](const std::string& key) -> std::optional<std::string> {
    std::string pref = "--" + key + "=";
    for (int i = 1; i < argc; ++i) {
      std::string s = argv[i];
      if (s.rfind(pref, 0) == 0) {
        return s.substr(pref.size());
      }
    }
    return std::nullopt;
  };
  auto strip = [&](const std::string& key) {
    std::string pref = "--" + key + "=";
    int w = 0;
    for (int i = 0; i < argc; ++i) {
      if (i > 0 && std::string(argv[i]).rfind(pref, 0) == 0) {
        continue;
      }
      argv[w++] = argv[i];
    }
    argc = w;
  };

  if (auto v = getv("min-exp")) {
    args.min_exp = std::stoi(*v), strip("min-exp");
  }
  if (auto v = getv("max-exp")) {
    args.max_exp = std::stoi(*v), strip("max-exp");
  }
  if (auto v = getv("q")) {
    args.Q = std::stoull(*v), strip("q");
  }
  if (auto v = getv("p")) {
    args.p1 = std::stod(*v), strip("p");
  }
  if (auto v = getv("seed")) {
    args.seed = std::stoull(*v), strip("seed");
  }
  if (auto v = getv("block")) {
    args.block_bits = std::stoull(*v), strip("block");
  }
  if (auto v = getv("per-octave")) {
    args.per_octave = std::stoi(*v), strip("per-octave");
  }

  if (auto v = getv("sizes")) {
    std::string s = *v;
    strip("sizes");
    std::size_t pos = 0;
    while (pos < s.size()) {
      while (pos < s.size() &&
             (s[pos] == ',' || std::isspace((unsigned char)s[pos]))) {
        ++pos;
      }
      std::size_t start = pos;
      while (pos < s.size() && std::isdigit((unsigned char)s[pos])) {
        ++pos;
      }
      if (start < pos) {
        args.explicit_sizes.push_back(
            std::stoull(s.substr(start, pos - start)));
      }
    }
    std::sort(args.explicit_sizes.begin(), args.explicit_sizes.end());
    args.explicit_sizes.erase(
        std::unique(args.explicit_sizes.begin(), args.explicit_sizes.end()),
        args.explicit_sizes.end());
  }
}

static std::string make_random_bits(std::size_t N,
                                    double p1,
                                    std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::string s;
  s.resize(N);
  for (std::size_t i = 0; i < N; ++i) {
    s[i] = (U(rng) < p1 ? '1' : '0');
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
  d.N = N;
  d.bits = make_random_bits(N, args.p1, rng);

  if (args.block_bits == 0) {
    d.t = RmMTree(d.bits, 0, -1.0f);
  } else {
    d.t = RmMTree(d.bits, args.block_bits);
  }

  d.cnt1 = d.t.rank1(N);
  d.cnt0 = N - d.cnt1;
  d.cnt10 = count10(d.bits);

  const std::size_t L = std::min<std::size_t>(args.Q, 32768);

  std::uniform_int_distribution<std::size_t> ind_dist_incl(0, N ? N : 0);
  std::uniform_int_distribution<std::size_t> ind_dist(0, N ? (N - 1) : 0);
  std::uniform_int_distribution<std::size_t> ind_dist_1N(1, N ? N : 1);
  std::uniform_int_distribution<int> d_dist(-8, +8);

  d.pool.inds_any.resize(L);
  d.pool.inds.resize(L);
  d.pool.inds_1N.resize(L);
  d.pool.deltas.resize(L);
  d.pool.segs.resize(L);

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

  for (std::size_t i = 0; i < L; ++i) {
    d.pool.inds_any[i] = ind_dist_incl(rng);
    d.pool.inds[i] = (N ? ind_dist(rng) : 0);
    d.pool.inds_1N[i] = (N ? ind_dist_1N(rng) : 0);
    d.pool.deltas[i] = d_dist(rng);
    d.pool.segs[i] = rand_ij();
  }

  auto fill_ks = [&](std::size_t total, std::vector<std::size_t>& out) {
    out.resize(L);
    if (total == 0) {
      std::fill(out.begin(), out.end(), 1);
      return;
    }
    std::uniform_int_distribution<std::size_t> dist(1, total);
    for (std::size_t i = 0; i < L; ++i) {
      out[i] = dist(rng);
    }
  };
  fill_ks(d.cnt1, d.pool.ks1);
  fill_ks(d.cnt0, d.pool.ks0);
  fill_ks(d.cnt10, d.pool.ks10);

  d.pool.minselect_q.resize(L);
  for (std::size_t i = 0; i < L; ++i) {
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

    register_op("rank10", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.rank10(P.inds_any[k % P.inds_any.size()]);
                  benchmark::DoNotOptimize(r);
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
                  auto [i, j] = P.segs[k % P.segs.size()];
                  auto r = D.t.range_min_query_pos(i, j);
                  benchmark::DoNotOptimize(r);
                });

    register_op("range_min_query_val", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto [i, j] = P.segs[k % P.segs.size()];
                  auto r = D.t.range_min_query_val(i, j);
                  benchmark::DoNotOptimize(r);
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
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  if (D.N == 0) {
                    benchmark::DoNotOptimize(0);
                    return;
                  }
                  auto i = P.inds[k % P.inds.size()];
                  auto r = D.t.close(i);
                  benchmark::DoNotOptimize(r);
                });

    register_op("open", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.open(P.inds_1N[k % P.inds_1N.size()]);
                  benchmark::DoNotOptimize(r);
                });

    register_op("enclose", data,
                [&](benchmark::State&, const Dataset& D, std::size_t k) {
                  auto r = D.t.enclose(P.inds_1N[k % P.inds_1N.size()]);
                  benchmark::DoNotOptimize(r);
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
