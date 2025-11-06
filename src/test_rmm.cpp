#include <bits/stdc++.h>
#include "rmm_tree.h"
#include "naive_rmm_tree.h"

using std::size_t;
static constexpr size_t NPOS = std::numeric_limits<size_t>::max();

static std::string bits_to_parens(const std::string &bits)
{
    std::string s;
    s.reserve(bits.size());
    for (char c : bits)
        s.push_back(c == '1' ? '(' : ')');
    return s;
}
static std::string vecbits_to_string(const std::vector<uint8_t> &v)
{
    std::string s;
    s.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        s[i] = v[i] ? '1' : '0';
    return s;
}

static std::string random_bits(std::mt19937_64 &rng, size_t n)
{
    std::uniform_int_distribution<int> b01(0, 1);
    std::string s;
    s.resize(n);
    for (size_t i = 0; i < n; i++)
        s[i] = char('0' + b01(rng));
    return s;
}

static std::string random_dyck_bits(std::mt19937_64 &rng, size_t m)
{
    std::string s;
    s.resize(2 * m);
    std::bernoulli_distribution coin(0.5);
    int h = 0;
    for (size_t pos = 0; pos < 2 * m; ++pos)
    {
        size_t rem = 2 * m - pos;
        if (h == 0)
        {
            s[pos] = '1';
            ++h;
            continue;
        }
        if ((size_t)h == rem - 1)
        {
            s[pos] = '0';
            --h;
            continue;
        }
        s[pos] = coin(rng) ? '1' : '0';
        if (s[pos] == '1')
            ++h;
        else
            --h;
        if (h < 0)
        {
            s[pos] = '1';
            h += 2;
        }
    }

    if (h > 0)
    {
        for (size_t i = s.size(); i > 0 && h > 0; --i)
        {
            if (s[i - 1] == '1')
            {
                s[i - 1] = '0';
                --h;
            }
        }
    }
    return s;
}

template <class T>
static std::string show_sz_or_npos(T x)
{
    if constexpr (std::is_same_v<T, size_t>)
        return (x == NPOS ? std::string("npos") : std::to_string(x));
    else
        return std::to_string(x);
}

[[noreturn]] static void die_mismatch(
    const std::string &seed_info,
    const std::string &bits,
    const char *op_name,
    const std::vector<std::string> &args,
    const std::string &want,
    const std::string &got)
{
    std::cerr << "\n=== MISMATCH DETECTED ===\n";
    std::cerr << seed_info << "\n";
    std::cerr << "bits:   " << bits << "\n";
    std::cerr << "parens: " << bits_to_parens(bits) << "\n";
    std::cerr << "op: " << op_name << "(";
    for (size_t i = 0; i < args.size(); ++i)
    {
        if (i)
            std::cerr << ", ";
        std::cerr << args[i];
    }
    std::cerr << ")\n";
    std::cerr << "expected: " << want << "\n";
    std::cerr << "got:      " << got << "\n";
    std::cerr << "=========================\n";
    std::exit(1);
}

int main()
{
    std::random_device rd;
    uint64_t seed = ((uint64_t)rd() << 32) ^ (uint64_t)rd() ^ (uint64_t)std::chrono::high_resolution_clock::now().time_since_epoch().count();
    // seed = 4436852161472814998ull;
    std::mt19937_64 rng(seed);
    std::ostringstream oss;
    oss << "seed = " << seed;
    std::string seed_info = oss.str();
    std::cerr << seed_info << "\n";

    const size_t LOG_EVERY = 50;
    const size_t OPS_PER_CASE = 2000;
    const size_t MAX_N = 40960;
    const size_t B = 64;

    size_t iter = 0;
    size_t total_ops = 0;

    for (;;)
    {
        ++iter;

        std::uniform_int_distribution<int> coin(0, 1);
        std::uniform_int_distribution<int> len_u(1, (int)MAX_N);
        std::uniform_int_distribution<int> len_even(0, (int)(MAX_N / 2));

        std::string bits;
        if (coin(rng) == 0)
        {

            bits = random_bits(rng, len_u(rng));
        }
        else
        {

            size_t m = 1 + (size_t)len_even(rng);
            bits = random_dyck_bits(rng, m);
            if (bits.size() == 0)
                bits = "10";
        }

        pixie::RmMTree rm(bits, B);
        // std::cerr << "N=" << rm.N << " B=" << rm.B << " overhead=" << rm.built_overhead << "\n";
        NaiveRmM nv(bits);

        const size_t N = bits.size();
        const size_t ones = nv.rank1(N);
        const size_t zeros = N - ones;
        const size_t pairs10 = (N >= 2 ? nv.rank10(N) : 0);

        std::uniform_int_distribution<size_t> pos_i(0, N);
        std::uniform_int_distribution<size_t> pos_i_nz(0, N ? N - 1 : 0);
        std::uniform_int_distribution<int> d_dist(-(int)std::min<size_t>(N, 200), (int)std::min<size_t>(N, 200));

        for (size_t q = 0; q < OPS_PER_CASE; ++q)
        {
            ++total_ops;
            int which = std::uniform_int_distribution<int>(0, 17)(rng);

            size_t i = 0, j = 0;
            if (N > 0)
            {
                i = pos_i_nz(rng);
                j = pos_i_nz(rng);
                if (i > j)
                    std::swap(i, j);
            }

            switch (which)
            {
            case 0:
            { // rank1
                size_t x = pos_i(rng);
                auto a = nv.rank1(x);
                auto b = rm.rank1(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "rank1", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 1:
            { // rank0
                size_t x = pos_i(rng);
                auto a = nv.rank0(x);
                auto b = rm.rank0(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "rank0", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 2:
            { // select1
                size_t k = std::uniform_int_distribution<size_t>(0, ones + 3)(rng);
                auto a = nv.select1(k);
                auto b = rm.select1(k);
                if (a != b)
                    die_mismatch(seed_info, bits, "select1", {std::to_string(k)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 3:
            { // select0
                size_t k = std::uniform_int_distribution<size_t>(0, zeros + 3)(rng);
                auto a = nv.select0(k);
                auto b = rm.select0(k);
                if (a != b)
                    die_mismatch(seed_info, bits, "select0", {std::to_string(k)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 4:
            { // rank10
                size_t x = (N >= 2 ? std::uniform_int_distribution<size_t>(0, N)(rng) : 0);
                auto a = nv.rank10(x);
                auto b = (N >= 2 ? rm.rank10(x) : 0);
                if (a != b)
                    die_mismatch(seed_info, bits, "rank10", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 5:
            { // select10
                size_t k = std::uniform_int_distribution<size_t>(0, pairs10 + 3)(rng);
                auto a = nv.select10(k);
                auto b = rm.select10(k);
                if (a != b)
                    die_mismatch(seed_info, bits, "select10", {std::to_string(k)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 6:
            { // excess
                size_t x = pos_i(rng);
                auto a = nv.excess(x);
                auto b = rm.excess(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "excess", {std::to_string(x)}, std::to_string(a), std::to_string(b));
            }
            break;
            case 7:
            { // fwdsearch
                if (N == 0)
                    break;
                size_t start = pos_i_nz(rng);
                int d = d_dist(rng);
                auto a = nv.fwdsearch(start, d);
                auto b = rm.fwdsearch(start, d);
                if (a != b)
                    die_mismatch(seed_info, bits, "fwdsearch", {std::to_string(start), std::to_string(d)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 8:
            { // bwdsearch
                if (N == 0)
                    break;
                size_t start = pos_i_nz(rng);
                int d = d_dist(rng);
                auto a = nv.bwdsearch(start, d);
                auto b = rm.bwdsearch(start, d);
                if (a != b)
                    die_mismatch(seed_info, bits, "bwdsearch", {std::to_string(start), std::to_string(d)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 9:
            { // range_min_query_pos
                if (N == 0)
                    break;
                auto a = nv.range_min_query_pos(i, j);
                auto b = rm.range_min_query_pos(i, j);
                if (a != b)
                    die_mismatch(seed_info, bits, "range_min_query_pos", {std::to_string(i), std::to_string(j)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 10:
            { // range_min_query_val
                if (N == 0)
                    break;
                auto a = nv.range_min_query_val(i, j);
                auto b = rm.range_min_query_val(i, j);
                if (a != b)
                    die_mismatch(seed_info, bits, "range_min_query_val", {std::to_string(i), std::to_string(j)}, std::to_string(a), std::to_string(b));
            }
            break;
            case 11:
            { // mincount
                if (N == 0)
                    break;
                auto a = nv.mincount(i, j);
                auto b = rm.mincount(i, j);
                if (a != b)
                    die_mismatch(seed_info, bits, "mincount", {std::to_string(i), std::to_string(j)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 12:
            { // minselect
                if (N == 0)
                    break;
                size_t cnt = nv.mincount(i, j);
                size_t k = cnt == 0 ? 1 : std::uniform_int_distribution<size_t>(1, cnt + 1)(rng);
                auto a = nv.minselect(i, j, k);
                auto b = rm.minselect(i, j, k);
                if (a != b)
                    die_mismatch(seed_info, bits, "minselect", {std::to_string(i), std::to_string(j), std::to_string(k)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 13:
            { // range_max_query_pos
                if (N == 0)
                    break;
                auto a = nv.range_max_query_pos(i, j);
                auto b = rm.range_max_query_pos(i, j);
                if (a != b)
                    die_mismatch(seed_info, bits, "range_max_query_pos",
                                 {std::to_string(i), std::to_string(j)},
                                 show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 14:
            { // range_max_query_val
                if (N == 0)
                    break;
                auto a = nv.range_max_query_val(i, j);
                auto b = rm.range_max_query_val(i, j);
                if (a != b)
                    die_mismatch(seed_info, bits, "range_max_query_val",
                                 {std::to_string(i), std::to_string(j)},
                                 std::to_string(a), std::to_string(b));
            }
            break;
            case 15:
            { // close
                if (N == 0)
                    break;
                size_t x = pos_i_nz(rng); // i in [0..N-1]
                auto a = nv.close(x);
                auto b = rm.close(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "close", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 16:
            { // open
                if (N == 0)
                    break;
                size_t x = pos_i(rng);
                auto a = nv.open(x);
                auto b = rm.open(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "open", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            case 17:
            { // enclose
                if (N == 0)
                    break;
                size_t x = pos_i(rng);
                auto a = nv.enclose(x);
                auto b = rm.enclose(x);
                if (a != b)
                    die_mismatch(seed_info, bits, "enclose", {std::to_string(x)}, show_sz_or_npos(a), show_sz_or_npos(b));
            }
            break;
            }
        }

        if (iter % LOG_EVERY == 0)
        {
            std::cerr << "iter=" << iter
                      << " total_ops=" << total_ops
                      << " last_N=" << bits.size()
                      << " ok\n";
        }
    }
}
