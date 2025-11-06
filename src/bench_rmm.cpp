#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cctype>
#include <iostream>
#include <optional>
#include <memory>
#include "rmm_tree.h"

using namespace std;
using pixie::RmMTree;

struct Args
{
    int min_exp = 14;
    int max_exp = 22;
    size_t Q = 200000;
    double p1 = 0.5;
    uint64_t seed = 42;
    size_t block_bits = 64;
    int per_octave = 10;
    vector<size_t> explicit_sizes;
};

static Args args;

static void parse_args_and_strip(int &argc, char **&argv)
{
    auto getv = [&](const string &key) -> optional<string>
    {
        string pref = "--" + key + "=";
        for (int i = 1; i < argc; ++i)
        {
            string s = argv[i];
            if (s.rfind(pref, 0) == 0)
                return s.substr(pref.size());
        }
        return nullopt;
    };
    auto strip = [&](const string &key)
    {
        string pref = "--" + key + "=";
        int w = 0;
        for (int i = 0; i < argc; ++i)
        {
            if (i > 0 && string(argv[i]).rfind(pref, 0) == 0)
                continue;
            argv[w++] = argv[i];
        }
        argc = w;
    };

    if (auto v = getv("min-exp"))
        args.min_exp = stoi(*v), strip("min-exp");
    if (auto v = getv("max-exp"))
        args.max_exp = stoi(*v), strip("max-exp");
    if (auto v = getv("q"))
        args.Q = stoull(*v), strip("q");
    if (auto v = getv("p"))
        args.p1 = stod(*v), strip("p");
    if (auto v = getv("seed"))
        args.seed = stoull(*v), strip("seed");
    if (auto v = getv("block"))
        args.block_bits = stoull(*v), strip("block");
    if (auto v = getv("per-octave"))
        args.per_octave = stoi(*v), strip("per-octave");

    if (auto v = getv("sizes"))
    {
        string s = *v;
        strip("sizes");
        size_t pos = 0;
        while (pos < s.size())
        {
            while (pos < s.size() && (s[pos] == ',' || isspace((unsigned char)s[pos])))
                ++pos;
            size_t start = pos;
            while (pos < s.size() && isdigit((unsigned char)s[pos]))
                ++pos;
            if (start < pos)
                args.explicit_sizes.push_back(stoull(s.substr(start, pos - start)));
        }
        sort(args.explicit_sizes.begin(), args.explicit_sizes.end());
        args.explicit_sizes.erase(unique(args.explicit_sizes.begin(), args.explicit_sizes.end()), args.explicit_sizes.end());
    }
}

static string make_random_bits(size_t N, double p1, mt19937_64 &rng)
{
    uniform_real_distribution<double> U(0.0, 1.0);
    string s;
    s.resize(N);
    for (size_t i = 0; i < N; ++i)
        s[i] = (U(rng) < p1 ? '1' : '0');
    return s;
}

static vector<size_t> build_size_grid()
{
    if (!args.explicit_sizes.empty())
    {
        vector<size_t> v = args.explicit_sizes;
        v.erase(remove(v.begin(), v.end(), 0), v.end());
        return v;
    }

    const int lo = args.min_exp, hi = args.max_exp;
    if (args.per_octave <= 0)
    {
        vector<size_t> v;
        for (int e = lo; e <= hi; ++e)
            v.push_back(size_t(1) << e);
        return v;
    }

    int per_octave_steps = args.per_octave;
    set<size_t> N_cands;
    size_t min_N = size_t(1) << lo, max_N = size_t(1) << hi;

    for (int e = lo; e <= hi; ++e)
    {
        for (int t = 0; t <= per_octave_steps; ++t)
        {
            long double x = e + (long double)t / (long double)per_octave_steps;
            size_t N = (size_t)llround(powl(2.0L, x));
            if (N < min_N)
                N = min_N;
            if (N > max_N)
                N = max_N;
            if (N)
                N_cands.insert(N);
        }
    }

    vector<size_t> v(N_cands.begin(), N_cands.end());
    if (v.front() != min_N)
        v.insert(v.begin(), min_N);
    if (v.back() != max_N)
        v.push_back(max_N);
    return v;
}

struct Pools
{
    vector<size_t> inds_any;
    vector<size_t> inds;
    vector<size_t> inds_1N;
    vector<int> deltas;
    vector<pair<size_t, size_t>> segs;

    vector<size_t> ks1, ks0, ks10;
    vector<size_t> minselect_q;
};

struct Dataset
{
    size_t N{};
    string bits;
    RmMTree t;

    size_t cnt1{}, cnt0{}, cnt10{};
    Pools pool;
};

static vector<shared_ptr<Dataset>> keepalive;

static size_t count10(const string &s)
{
    size_t c = 0;
    if (s.size() < 2)
        return 0;
    for (size_t i = 0; i + 1 < s.size(); ++i)
        if (s[i] == '1' && s[i + 1] == '0')
            ++c;
    return c;
}

static Dataset build_dataset(size_t N)
{
    mt19937_64 rng(args.seed ^ (uint64_t)N * 0x9E3779B185EBCA87ull);

    Dataset d;
    d.N = N;
    d.bits = make_random_bits(N, args.p1, rng);

    if (args.block_bits == 0)
        d.t = RmMTree(d.bits, 0, -1.0f);
    else
        d.t = RmMTree(d.bits, args.block_bits);

    d.cnt1 = d.t.rank1(N);
    d.cnt0 = N - d.cnt1;
    d.cnt10 = count10(d.bits);

    const size_t L = min<size_t>(args.Q, 32768);

    uniform_int_distribution<size_t> ind_dist_incl(0, N ? N : 0);
    uniform_int_distribution<size_t> ind_dist(0, N ? (N - 1) : 0);
    uniform_int_distribution<size_t> ind_dist_1N(1, N ? N : 1);
    uniform_int_distribution<int> d_dist(-8, +8);

    d.pool.inds_any.resize(L);
    d.pool.inds.resize(L);
    d.pool.inds_1N.resize(L);
    d.pool.deltas.resize(L);
    d.pool.segs.resize(L);

    auto rand_ij = [&]() -> pair<size_t, size_t>
    {
        if (N == 0)
            return {0, 0};
        size_t a = ind_dist(rng), b = ind_dist(rng);
        if (a > b)
            swap(a, b);
        return {a, b};
    };

    for (size_t i = 0; i < L; ++i)
    {
        d.pool.inds_any[i] = ind_dist_incl(rng);
        d.pool.inds[i] = (N ? ind_dist(rng) : 0);
        d.pool.inds_1N[i] = (N ? ind_dist_1N(rng) : 0);
        d.pool.deltas[i] = d_dist(rng);
        d.pool.segs[i] = rand_ij();
    }

    auto fill_ks = [&](size_t total, vector<size_t> &out)
    {
        out.resize(L);
        if (total == 0)
        {
            fill(out.begin(), out.end(), 1);
            return;
        }
        uniform_int_distribution<size_t> dist(1, total);
        for (size_t i = 0; i < L; ++i)
            out[i] = dist(rng);
    };
    fill_ks(d.cnt1, d.pool.ks1);
    fill_ks(d.cnt0, d.pool.ks0);
    fill_ks(d.cnt10, d.pool.ks10);

    d.pool.minselect_q.resize(L);
    for (size_t i = 0; i < L; ++i)
    {
        auto [l, r] = d.pool.segs[i];
        size_t c = d.t.mincount(l, r);
        if (c == 0)
            c = 1;
        uniform_int_distribution<size_t> uq(1, c);
        d.pool.minselect_q[i] = uq(rng);
    }

    return d;
}

template <class Fn>
static void register_op(const string &op, shared_ptr<Dataset> data, Fn &&body)
{
    auto idx_ptr = make_shared<size_t>(0);

    auto *b = benchmark::RegisterBenchmark(
        op.c_str(),
        [data, idx_ptr, body](benchmark::State &state)
        {
            const Dataset &D = *data;
            for (auto _ : state)
            {
                size_t i = (*idx_ptr)++;
                body(state, D, i);
            }
            state.counters["N"] = static_cast<double>(D.N);
            state.counters["seed"] = static_cast<double>(args.seed);
            state.counters["block_bits"] = static_cast<double>(args.block_bits);
        });

    b->Unit(benchmark::kNanosecond);
}

static void register_all()
{
    auto Ns = build_size_grid();
    for (size_t N : Ns)
    {
        auto data = make_shared<Dataset>(build_dataset(N));
        keepalive.push_back(data);

        const auto &P = data->pool;

        register_op("rank1", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.rank1(P.inds_any[k % P.inds_any.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("rank0", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.rank0(P.inds_any[k % P.inds_any.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("select1", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.select1(P.ks1[k % P.ks1.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("select0", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.select0(P.ks0[k % P.ks0.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("rank10", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.rank10(P.inds_any[k % P.inds_any.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("select10", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.select10(P.ks10[k % P.ks10.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("excess", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.excess(P.inds_any[k % P.inds_any.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("fwdsearch", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.fwdsearch(P.inds[k % P.inds.size()], P.deltas[k % P.deltas.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("bwdsearch", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.bwdsearch(P.inds_1N[k % P.inds_1N.size()], P.deltas[k % P.deltas.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("range_min_query_pos", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto [i, j] = P.segs[k % P.segs.size()];
            auto r = D.t.range_min_query_pos(i, j);
            benchmark::DoNotOptimize(r); });

        register_op("range_min_query_val", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto [i, j] = P.segs[k % P.segs.size()];
            auto r = D.t.range_min_query_val(i, j);
            benchmark::DoNotOptimize(r); });

        register_op("range_max_query_pos", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto [i, j] = P.segs[k % P.segs.size()];
            auto r = D.t.range_max_query_pos(i, j);
            benchmark::DoNotOptimize(r); });

        register_op("range_max_query_val", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto [i, j] = P.segs[k % P.segs.size()];
            auto r = D.t.range_max_query_val(i, j);
            benchmark::DoNotOptimize(r); });

        register_op("mincount", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto [i, j] = P.segs[k % P.segs.size()];
            auto r = D.t.mincount(i, j);
            benchmark::DoNotOptimize(r); });

        register_op("minselect", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            const size_t idx = k % P.segs.size();
            auto [i, j] = P.segs[idx];
            auto q = P.minselect_q[idx];
            auto r = D.t.minselect(i, j, q);
            benchmark::DoNotOptimize(r); });

        register_op("close", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            if (D.N == 0) { benchmark::DoNotOptimize(0); return; }
            auto i = P.inds[k % P.inds.size()];
            auto r = D.t.close(i);
            benchmark::DoNotOptimize(r); });

        register_op("open", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.open(P.inds_1N[k % P.inds_1N.size()]);
            benchmark::DoNotOptimize(r); });

        register_op("enclose", data, [&](benchmark::State &, const Dataset &D, size_t k)
                    {
            auto r = D.t.enclose(P.inds_1N[k % P.inds_1N.size()]);
            benchmark::DoNotOptimize(r); });
    }
}

int main(int argc, char **argv)
{
    parse_args_and_strip(argc, argv);

    auto has = [&](const char *key)
    {
        string p1 = string(key) + "=";
        for (int i = 1; i < argc; ++i)
        {
            string s = argv[i];
            if (s == key || s.rfind(p1, 0) == 0)
                return true;
        }
        return false;
    };

    static vector<string> extra;
    if (!has("--benchmark_format"))
        extra.emplace_back("--benchmark_format=csv");
    if (!has("--benchmark_counters_tabular"))
        extra.emplace_back("--benchmark_counters_tabular=true");
    if (!has("--benchmark_time_unit"))
        extra.emplace_back("--benchmark_time_unit=ns");

    static vector<char *> argv_vec;
    argv_vec.assign(argv, argv + argc);
    for (auto &s : extra)
        argv_vec.push_back(s.data());
    argv_vec.push_back(nullptr);
    argc = (int)argv_vec.size() - 1;
    argv = argv_vec.data();

    benchmark::Initialize(&argc, argv);
    register_all();
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
