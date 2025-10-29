#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <limits>
#include <algorithm>

struct NaiveRmM
{
    using u64 = std::uint64_t;
    static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

    std::vector<uint8_t> b;
    std::size_t N = 0;

    NaiveRmM() = default;
    explicit NaiveRmM(const std::string &s) { build_from_string(s); }
    NaiveRmM(const std::vector<u64> &words, std::size_t nbits) { build_from_words(words, nbits); }

    void build_from_string(const std::string &s)
    {
        N = s.size();
        b.resize(N);
        for (std::size_t i = 0; i < N; i++)
            b[i] = (s[i] == '1');
    }
    void build_from_words(const std::vector<u64> &words, std::size_t nbits)
    {
        N = nbits;
        b.assign(N, 0);
        for (std::size_t i = 0; i < N; i++)
            b[i] = ((words[i >> 6] >> (i & 63)) & 1u);
    }

    inline int bit(std::size_t i) const { return b[i]; }
    inline std::size_t size() const { return N; }

    std::size_t rank1(std::size_t i) const
    {
        if (i > N)
            i = N;
        std::size_t c = 0;
        for (std::size_t p = 0; p < i; p++)
            c += (b[p] != 0);
        return c;
    }
    std::size_t rank0(std::size_t i) const { return i - rank1(i); }

    // 1-based
    std::size_t select1(std::size_t k) const
    {
        if (k == 0)
            return npos;
        for (std::size_t p = 0; p < N; p++)
            if (b[p])
            {
                if (--k == 0)
                    return p;
            }
        return npos;
    }

    // 1-based
    std::size_t select0(std::size_t k) const
    {
        if (k == 0)
            return npos;
        for (std::size_t p = 0; p < N; p++)
            if (!b[p])
            {
                if (--k == 0)
                    return p;
            }
        return npos;
    }

    std::size_t rank10(std::size_t i) const
    {
        if (i <= 1)
            return 0;
        std::size_t c = 0;
        for (std::size_t p = 0; p + 1 < i; ++p)
            if (b[p] == 1 && b[p + 1] == 0)
                ++c;
        return c;
    }

    // 1-based
    std::size_t select10(std::size_t k) const
    {
        if (k == 0)
            return npos;
        for (std::size_t p = 0; p + 1 < N; ++p)
            if (b[p] == 1 && b[p + 1] == 0)
            {
                if (--k == 0)
                    return p;
            }
        return npos;
    }

    int excess(std::size_t i) const
    {

        return int(2 * static_cast<long long>(rank1(i)) - static_cast<long long>(i));
    }

    std::size_t fwdsearch(std::size_t i, int d) const
    {
        if (i >= N)
            return npos;
        int target = excess(i) + d;
        int cur = excess(i);
        for (std::size_t p = i; p < N; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur == target)
                return p;
        }
        return npos;
    }

    std::size_t bwdsearch(std::size_t i, int d) const
    {

        if (i > N)
            return npos;
        if (i == 0)
            return npos;
        int target = excess(i) + d;

        int cur = excess(i);
        for (std::size_t p = i; p > 0;)
        {
            --p;
            cur += b[p] ? -1 : +1;
            if (cur == target)
                return p;
        }
        return npos;
    }

    std::size_t rmq_pos(std::size_t i, std::size_t j) const
    {
        if (i > j || j >= N)
            return npos;
        int cur = 0, mn = std::numeric_limits<int>::max();
        std::size_t pos = npos;
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur < mn)
            {
                mn = cur;
                pos = p;
            }
        }
        return pos;
    }

    int rmq_val(std::size_t i, std::size_t j) const
    {
        if (i > j || j >= N)
            return 0;
        int cur = 0, mn = std::numeric_limits<int>::max();
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur < mn)
                mn = cur;
        }
        return mn;
    }

    std::size_t rMq_pos(std::size_t i, std::size_t j) const
    {
        if (i > j || j >= N)
            return npos;
        int cur = 0, mx = std::numeric_limits<int>::min();
        std::size_t pos = npos;
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur > mx)
            {
                mx = cur;
                pos = p;
            }
        }
        return pos;
    }

    int rMq_val(std::size_t i, std::size_t j) const
    {
        if (i > j || j >= N)
            return 0;
        int cur = 0, mx = std::numeric_limits<int>::min();
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur > mx)
                mx = cur;
        }
        return mx;
    }

    std::size_t mincount(std::size_t i, std::size_t j) const
    {
        if (i > j || j >= N)
            return 0;
        int cur = 0, mn = std::numeric_limits<int>::max();
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur < mn)
                mn = cur;
        }
        std::size_t cnt = 0;
        cur = 0;
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur == mn)
                ++cnt;
        }
        return cnt;
    }

    // (1-based q)
    std::size_t minselect(std::size_t i, std::size_t j, std::size_t q) const
    {
        if (i > j || j >= N || q == 0)
            return npos;
        int cur = 0, mn = std::numeric_limits<int>::max();
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur < mn)
                mn = cur;
        }
        cur = 0;
        for (std::size_t p = i; p <= j; ++p)
        {
            cur += b[p] ? +1 : -1;
            if (cur == mn)
            {
                if (--q == 0)
                    return p;
            }
        }
        return npos;
    }

    std::size_t close(std::size_t i) const
    {
        if (i >= N)
            return npos;
        return fwdsearch(i, -1);
    }
    std::size_t open(std::size_t i) const
    {
        if (i == 0 || i > N)
            return npos;
        auto r = bwdsearch(i, 0);
        return (r == npos ? npos : r + 1);
    }
    std::size_t enclose(std::size_t i) const
    {
        if (i == 0 || i > N)
            return npos;
        auto r = bwdsearch(i, -2);
        return (r == npos ? npos : r + 1);
    }
};
