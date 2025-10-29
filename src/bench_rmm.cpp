#include <bits/stdc++.h>
#include "rmm_tree.h"

using namespace std;

static volatile size_t BH_sz = 0;
static volatile int BH_i = 0;

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

static Args parse_args(int argc, char **argv)
{
    Args a;
    auto getv = [&](string key) -> optional<string>
    {
        string pref = "--" + key + "=";
        for (int i = 1; i < argc; i++)
        {
            string s = argv[i];
            if (s.rfind(pref, 0) == 0)
                return s.substr(pref.size());
        }
        return nullopt;
    };
    if (auto v = getv("min-exp"))
        a.min_exp = stoi(*v);
    if (auto v = getv("max-exp"))
        a.max_exp = stoi(*v);
    if (auto v = getv("q"))
        a.Q = stoull(*v);
    if (auto v = getv("p"))
        a.p1 = stod(*v);
    if (auto v = getv("seed"))
        a.seed = stoull(*v);
    if (auto v = getv("block"))
        a.block_bits = stoull(*v);
    if (auto v = getv("per-octave"))
        a.per_octave = stoi(*v);
    if (auto v = getv("sizes"))
    {
        string s = *v;
        size_t pos = 0;
        while (pos < s.size())
        {
            while (pos < s.size() && (s[pos] == ',' || isspace((unsigned char)s[pos])))
                ++pos;
            if (pos >= s.size())
                break;
            size_t start = pos;
            while (pos < s.size() && isdigit((unsigned char)s[pos]))
                ++pos;
            if (start < pos)
                a.explicit_sizes.push_back(stoull(s.substr(start, pos - start)));
        }
    }
    return a;
}

static string make_random_bits(size_t N, double p1, mt19937_64 &rng)
{
    uniform_real_distribution<double> U(0.0, 1.0);
    string s;
    s.resize(N);
    for (size_t i = 0; i < N; i++)
        s[i] = (U(rng) < p1 ? '1' : '0');
    return s;
}

static vector<size_t> build_size_grid(const Args &args)
{
    if (!args.explicit_sizes.empty())
    {
        vector<size_t> v = args.explicit_sizes;
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());
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
    for (int e = lo; e <= hi; ++e)
    {
        for (int t = 0; t <= per_octave_steps; ++t)
        {
            double x = e + double(t) / double(per_octave_steps);
            long double size_float = powl(2.0L, x);
            size_t N = (size_t)llround(size_float);
            if (N == 0)
                continue;
            size_t min_N = size_t(1) << lo;
            size_t max_N = size_t(1) << hi;
            if (N < min_N)
                N = min_N;
            if (N > max_N)
                N = max_N;
            N_cands.insert(N);
        }
    }
    vector<size_t> v(N_cands.begin(), N_cands.end());

    size_t min_N = size_t(1) << lo, max_N = size_t(1) << hi;
    if (v.front() != min_N)
        v.insert(v.begin(), min_N);
    if (v.back() != max_N)
        v.push_back(max_N);
    return v;
}

template <class F>
static double time_ns_per_call(F &&f, size_t q)
{
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    for (size_t i = 0; i < q; i++)
        f(i);
    auto t1 = clk::now();
    double ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    return ns / (double)q;
}

int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Args args = parse_args(argc, argv);
    mt19937_64 rng(args.seed);

    vector<size_t> Ns = build_size_grid(args);

    cout << "N,op,ns_per_op,queries,hits,misses,seed,block_bits\n";

    for (size_t N : Ns)
    {

        string bits = make_random_bits(N, args.p1, rng);
        pixie::RmMTree t(bits, args.block_bits);
        // pixie::RmMTree t(bits, 0, 0.15f);

        uniform_int_distribution<size_t> ind_dist_incl(0, N ? N : 0);
        uniform_int_distribution<size_t> ind_dist(0, N ? (N - 1) : 0);
        auto rand_i_any = [&]()
        { return (N ? ind_dist_incl(rng) : 0); };
        auto rand_i = [&]()
        { return (N ? ind_dist(rng) : 0); };
        auto rand_ij = [&]()
        {
            if (N == 0)
                return pair<size_t, size_t>(0, 0);
            size_t a = ind_dist(rng), b = ind_dist(rng);
            if (a > b)
                swap(a, b);
            return pair<size_t, size_t>(a, b);
        };
        uniform_int_distribution<int> d_dist(-8, +8);

        size_t cnt1 = t.rank1(N);
        size_t cnt0 = N - cnt1;

        size_t cnt10 = 0;
        for (size_t p = 0; p + 1 < N; ++p)
            if (bits[p] == '1' && bits[p + 1] == '0')
                ++cnt10;

        // rank1
        {
            vector<size_t> inds(args.Q);
            for (size_t k = 0; k < args.Q; k++)
                inds[k] = rand_i_any();
            auto fn = [&](size_t it)
            { BH_sz += t.rank1(inds[it]); };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rank1," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // rank0
        {
            vector<size_t> inds(args.Q);
            for (size_t k = 0; k < args.Q; k++)
                inds[k] = rand_i_any();
            auto fn = [&](size_t it)
            { BH_sz += t.rank0(inds[it]); };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rank0," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // select1
        {
            size_t q = min(args.Q, max<size_t>(cnt1, (size_t)1));
            vector<size_t> ks;
            ks.reserve(q);
            if (cnt1 > 0)
            {
                uniform_int_distribution<size_t> select_dist(1, cnt1);
                for (size_t i = 0; i < q; i++)
                    ks.push_back(select_dist(rng));
            }
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                size_t k = (cnt1 ? ks[it] : (size_t)1);
                auto r = t.select1(k);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",select1," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // select0
        {
            size_t q = min(args.Q, max<size_t>(cnt0, (size_t)1));
            vector<size_t> ks;
            ks.reserve(q);
            if (cnt0 > 0)
            {
                uniform_int_distribution<size_t> select_dist(1, cnt0);
                for (size_t i = 0; i < q; i++)
                    ks.push_back(select_dist(rng));
            }
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                size_t k = (cnt0 ? ks[it] : (size_t)1);
                auto r = t.select0(k);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",select0," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // rank10
        {
            vector<size_t> inds(args.Q);
            for (size_t k = 0; k < args.Q; k++)
                inds[k] = rand_i_any();
            auto fn = [&](size_t it)
            { BH_sz += t.rank10(inds[it]); };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rank10," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // select10
        {
            size_t q = min(args.Q, max<size_t>(cnt10, (size_t)1));
            vector<size_t> ks;
            ks.reserve(q);
            if (cnt10 > 0)
            {
                uniform_int_distribution<size_t> select_dist(1, cnt10);
                for (size_t i = 0; i < q; i++)
                    ks.push_back(select_dist(rng));
            }
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                size_t k = (cnt10 ? ks[it] : (size_t)1);
                auto r = t.select10(k);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",select10," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // excess
        {
            vector<size_t> inds(args.Q);
            for (size_t k = 0; k < args.Q; k++)
                inds[k] = rand_i_any();
            auto fn = [&](size_t it)
            { BH_i += t.excess(inds[it]); };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",excess," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // fwdsearch
        {
            size_t q = args.Q;
            vector<size_t> inds(q);
            vector<int> ds(q);
            for (size_t k = 0; k < q; k++)
            {
                inds[k] = (N ? ind_dist(rng) : 0);
                ds[k] = d_dist(rng);
            }
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                auto r = t.fwdsearch(inds[it], ds[it]);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",fwdsearch," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // bwdsearch
        {
            size_t q = args.Q;
            vector<size_t> inds(q);
            vector<int> ds(q);
            for (size_t k = 0; k < q; k++)
            {
                inds[k] = (N ? ind_dist(rng) : 0);
                ds[k] = d_dist(rng);
            }
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                auto r = t.bwdsearch(inds[it], ds[it]);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",bwdsearch," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // rmq_pos / rmq_val / mincount / minselect
        vector<pair<size_t, size_t>> segs(args.Q);
        for (size_t k = 0; k < args.Q; k++)
            segs[k] = rand_ij();

        // rmq_pos
        {
            auto fn = [&](size_t it)
            {
                auto [i, j] = segs[it];
                BH_sz += t.rmq_pos(i, j);
            };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rmq_pos," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // rmq_val
        {
            auto fn = [&](size_t it)
            {
                auto [i, j] = segs[it];
                BH_i += t.rmq_val(i, j);
            };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rmq_val," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // rMq_pos
        {
            auto fn = [&](size_t it)
            {
                auto [i, j] = segs[it];
                BH_sz += t.rMq_pos(i, j);
            };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rMq_pos," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // rMq_val
        {
            auto fn = [&](size_t it)
            {
                auto [i, j] = segs[it];
                BH_i += t.rMq_val(i, j);
            };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",rMq_val," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // mincount
        vector<size_t> seg_counts(args.Q, 0);
        {
            auto fn = [&](size_t it)
            {
                auto [i, j] = segs[it];
                auto c = t.mincount(i, j);
                seg_counts[it] = c;
                BH_sz += c;
            };
            double ns = time_ns_per_call(fn, args.Q);
            cout << N << ",mincount," << ns << "," << args.Q << ",0,0," << args.seed << "," << args.block_bits << "\n";
        }

        // minselect
        {
            vector<size_t> good_inds;
            good_inds.reserve(args.Q);
            for (size_t it = 0; it < args.Q; ++it)
                if (seg_counts[it] > 0)
                    good_inds.push_back(it);

            size_t q = good_inds.size();
            vector<size_t> qs;
            qs.reserve(q);
            for (size_t k = 0; k < q; k++)
            {
                size_t it = good_inds[k];
                uniform_int_distribution<size_t> query_dist(1, seg_counts[it]);
                qs.push_back(query_dist(rng));
            }

            size_t hits = 0;
            auto fn = [&](size_t r)
            {
                size_t it = good_inds[r];
                auto [i, j] = segs[it];
                auto pos = t.minselect(i, j, qs[r]);
                hits += (pos != pixie::RmMTree::npos);
                BH_sz += pos;
            };
            double ns = (q ? time_ns_per_call(fn, q) : numeric_limits<double>::quiet_NaN());
            cout << N << ",minselect," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // close
        {
            if (N > 0)
            {
                size_t q = args.Q;
                vector<size_t> inds(q);
                for (size_t k = 0; k < q; k++)
                    inds[k] = (N ? ind_dist(rng) : 0); // i in [0..N-1]
                size_t hits = 0;
                auto fn = [&](size_t it)
                {
                    auto r = t.close(inds[it]);
                    hits += (r != pixie::RmMTree::npos);
                    BH_sz += r;
                };
                double ns = time_ns_per_call(fn, q);
                cout << N << ",close," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
            }
            else
            {
                cout << N << ",close," << std::numeric_limits<double>::quiet_NaN()
                     << ",0,0,0," << args.seed << "," << args.block_bits << "\n";
            }
        }

        // open
        {
            size_t q = args.Q;
            vector<size_t> inds(q);
            for (size_t k = 0; k < q; k++)
                inds[k] = (N ? ind_dist_incl(rng) : 0); // i in [0..N]
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                auto r = t.open(inds[it]);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",open," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }

        // enclose
        {
            size_t q = args.Q;
            vector<size_t> inds(q);
            for (size_t k = 0; k < q; k++)
                inds[k] = (N ? ind_dist_incl(rng) : 0); // i in [0..N]
            size_t hits = 0;
            auto fn = [&](size_t it)
            {
                auto r = t.enclose(inds[it]);
                hits += (r != pixie::RmMTree::npos);
                BH_sz += r;
            };
            double ns = time_ns_per_call(fn, q);
            cout << N << ",enclose," << ns << "," << q << "," << hits << "," << (q - hits) << "," << args.seed << "," << args.block_bits << "\n";
        }
    }

    cerr << "blackholes: " << BH_sz << " / " << BH_i << "\n";
    return 0;
}
