#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <climits>
#include <limits>
#include <string>
#include <vector>

namespace pixie
{

#ifdef __GNUC__
#define POPCNT __builtin_popcountll
#else
#include <intrin.h>
    static inline int POPCNT(unsigned long long x) { return __popcnt64(x); }
#endif

    /**
     * @brief Range min–max tree over a bitvector (LSB-first) tailored for balanced-parentheses (BP).
     * @details
     *  Implements the classic RMQ/rMq family and BP navigation in O(log n)
     *  using a perfectly balanced binary tree of blocks. Supported queries:
     *  - rank1 / rank0
     *  - select1 / select0
     *  - rank10 / select10 (starts of "10")
     *  - excess (prefix +1 for '1', −1 for '0')
     *  - fwdsearch / bwdsearch (prefix–sum search)
     *  - rmq_pos / rmq_val (first minimum within a range)
     *  - rMq_pos / rMq_val (first maximum within a range)
     *  - mincount / minselect (count/selection of minima)
     *  - close / open / enclose (BP navigation)
     *
     *  The bitvector is LSB-first inside each 64-bit word.
     */
    class RmMTree
    {
        // ------------ bitvector ------------
        std::vector<std::uint64_t> bits; // LSB-first
        size_t num_bits = 0;             // number of bits

        // ------------ blocking ------------
        size_t block_bits = 64; // block size (bits), leaf covers <= block_bits bits
        size_t leaf_count = 0;  // #leaves = ceil(num_bits/block_bits)

        // ------------ tree arrays (heap order: 1 is root) ------------
        // size of segment (in bits) covered by node
        // needed for: rank1/rank0, select1/select0, select10,
        //             excess, fwdsearch/bwdsearch/close/open/enclose,
        //             rmq/rMq, minselect.
        std::vector<uint32_t> segment_size_bits;

        // node_total_excess = total excess (+1 for '1', -1 for '0') on the node
        // needed for: rank1/rank0, select1/select0, excess,
        //             fwdsearch/bwdsearch/close/open/enclose,
        //             rmq/rMq, mincount/minselect.
        std::vector<int32_t> node_total_excess;

        // node_min_prefix_excess = minimum pref-excess on the node (from 0)
        // needed for: fwdsearch/bwdsearch/close/open/enclose, rmq, mincount/minselect.
        std::vector<int32_t> node_min_prefix_excess;

        // node_max_prefix_excess = maximum pref-excess on the node (from 0)
        // needed for: fwdsearch/bwdsearch/close/open/enclose, rMq.
        std::vector<int32_t> node_max_prefix_excess;

        // node_min_count = number of positions where the minimum is attained
        // needed for: mincount/minselect.
        std::vector<uint32_t> node_min_count;

        // node_pattern10_count = # of "10" pattern occurrences inside the node
        // needed for: rank10, select10.
        std::vector<uint32_t> node_pattern10_count;

        // node_first_bit = first bit (0/1), node_last_bit = last bit (0/1) of the segment (to handle "10" crossing)
        // both needed for: rank10, select10.
        std::vector<uint8_t> node_first_bit, node_last_bit;

    public:
        /**
         * @brief Sentinel for "not found".
         */
        static constexpr size_t npos = std::numeric_limits<size_t>::max();

#ifdef DEBUG
        float built_overhead = 0.0;
#endif

        // --------- construction ----------

        /**
         * @brief Construct empty structure.
         */
        RmMTree() = default;

        /**
         * @brief Build from a '0'/'1' string.
         * @param bp Bitstring (characters '0' and '1').
         * @param leaf_block_bits Desired leaf size (power of two, 0 = auto).
         * @param max_overhead Max allowed overhead fraction (<0 to disable constraint).
         * @details Block size priority: (1) respect @p max_overhead, (2) explicit @p leaf_block_bits,
         *          (3) set to ceil_pow2(log2(num_bits)).
         */
        explicit RmMTree(const std::string &bp,
                         const size_t &leaf_block_bits /*0=auto*/ = 0,
                         const float &max_overhead /*<0=off*/ = -1.0)
        {
            build_from_string(bp, leaf_block_bits, max_overhead);
        }

        /**
         * @brief Build from 64-bit words (LSB-first).
         * @param words Array of words holding bits LSB-first.
         * @param Nbits Number of valid bits.
         * @param leaf_block_bits Desired leaf size (power of two, 0 = auto).
         * @param max_overhead Max allowed overhead fraction (<0 to disable constraint).
         * @details Block size priority: (1) respect @p max_overhead, (2) explicit @p leaf_block_bits,
         *          (3) set to ceil_pow2(log2(num_bits)).
         */
        explicit RmMTree(const std::vector<std::uint64_t> &words, size_t Nbits,
                         const size_t &leaf_block_bits /*0=auto*/ = 0,
                         const float &max_overhead /*<0=off*/ = -1.0)
        {
            build_from_words(words, Nbits, leaf_block_bits, max_overhead);
        }

        // --------- queries: rank/select/excess ----------

        /**
         * @brief Number of ones in prefix [0, i).
         * @details Returns 0 for i==0. Complexity: O(log n) with small constants.
         */
        size_t rank1(const size_t &i) const
        {
            if (i == 0)
                return 0;
            const size_t blk = block_of(i - 1);
            size_t ans = 0;
            if (blk > 0)
            {
                const auto nodes = cover_blocks(0, blk - 1);
                for (const size_t &v : nodes)
                    ans += ones_in_node(v);
            }
            const size_t Lb = blk * block_bits;
            const size_t Rb = std::min(num_bits, Lb + block_bits);
            ans += rank1_in_block(Lb, std::min(i, Rb));
            return ans;
        }

        /**
         * @brief Number of zeros in prefix [0, i).
         * @details Computed as i - rank1(i).
         */
        size_t rank0(const size_t &i) const
        {
            return i - rank1(i);
        }

        /**
         * @brief 1-based select of the k-th one.
         * @return Position of k-th '1' or npos if not found.
         */
        size_t select1(size_t k) const
        {
            if (k == 0)
                return npos;
            size_t v = 1;
            if (ones_in_node(v) < k)
                return npos;
            size_t base = 0;
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1, Rc = Lc | 1;
                const uint32_t o_l = ones_in_node(Lc);
                if (o_l >= k)
                {
                    v = Lc;
                }
                else
                {
                    k -= o_l;
                    base += segment_size_bits[Lc];
                    v = Rc;
                }
            }
            return select1_in_block(base, std::min(base + segment_size_bits[v], num_bits), k);
        }

        /**
         * @brief 1-based select of the k-th zero.
         * @return Position of k-th '0' or npos if not found.
         */
        size_t select0(size_t k) const
        {
            if (k == 0)
                return npos;
            size_t v = 1;
            const auto zeros = [&](const size_t &x) noexcept
            {
                return segment_size_bits[x] - ones_in_node(x);
            };
            if (zeros(v) < k)
                return npos;
            size_t base = 0;
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1, Rc = Lc | 1;
                const size_t z_l = zeros(Lc);
                if (z_l >= k)
                {
                    v = Lc;
                }
                else
                {
                    k -= z_l;
                    base += segment_size_bits[Lc];
                    v = Rc;
                }
            }
            return select0_in_block(base, std::min(base + segment_size_bits[v], num_bits), k);
        }

        /**
         * @brief Rank of the pattern "10" (starts) within [0, i).
         * @details Counts p where bit[p]==1 and bit[p+1]==0 with p+1<i.
         */
        size_t rank10(const size_t &i) const
        {
            if (i <= 1)
                return 0;
            const size_t blk = block_of(i - 1);
            size_t ans = 0;
            int prev_last = -1;

            if (blk > 0)
            {
                const auto nodes = cover_blocks(0, blk - 1);
                for (const size_t &v : nodes)
                {
                    ans += node_pattern10_count[v];
                    if (prev_last != -1 && prev_last == 1 && node_first_bit[v] == 0)
                        ++ans;
                    prev_last = node_last_bit[v];
                }
            }
            const size_t Lb = blk * block_bits;
            ans += rr_in_block(Lb, i);
            // boundary between the last full node and the leaf tail
            if (blk > 0 && i > Lb && prev_last == 1 && bit(Lb) == 0)
                ++ans;
            return ans;
        }

        /**
         * @brief 1-based select of the k-th "10" start.
         * @return Position p such that bits[p..p+1]=="10", or npos if not found.
         */
        size_t select10(size_t k) const
        {
            if (k == 0)
                return npos;
            size_t ind = 1;
            if (node_pattern10_count[ind] < k)
                return npos;
            size_t base = 0;
            const size_t leaf0 = first_leaf_index();
            while (ind < leaf0)
            {
                const size_t Lc = ind << 1, Rc = Lc | 1;
                const size_t cross = (node_last_bit[Lc] == 1 && node_first_bit[Rc] == 0) ? 1u : 0u;
                if (node_pattern10_count[Lc] >= k)
                {
                    ind = Lc;
                    continue;
                }
                size_t rem = k - node_pattern10_count[Lc];
                if (cross)
                {
                    if (rem == 1)
                        return base + segment_size_bits[Lc] - 1;
                    --rem;
                }
                base += segment_size_bits[Lc];
                ind = Rc;
                k = rem;
            }
            return select10_in_block(base, std::min(base + segment_size_bits[ind], num_bits), k);
        }

        /**
         * @brief Prefix excess on [0, i): +1 for '1', −1 for '0'.
         */
        inline int excess(const size_t &i) const
        {
            return int64_t(rank1(i)) * 2 - int64_t(i);
        }

        /**
         * @brief Forward search: first position p ≥ i where excess(p) = excess(i) + d.
         * @details Scans remainder of current leaf, then descends using precomputed bounds.
         *          Returns npos if no such position exists.
         */
        size_t fwdsearch(const size_t &i, const int &d) const
        {
            if (i >= num_bits)
                return npos;

            const int start_excess = excess(i);
            const int target = start_excess + d;

            // 1) scan the remainder of the current leaf
            const size_t blk = block_of(i);
            const size_t Lb = blk * block_bits;
            const size_t Rb = std::min(num_bits, Lb + block_bits);
            int cur = start_excess;
            for (size_t p = i; p < Rb; ++p)
            {
                cur += bit(p) ? +1 : -1;
                if (cur == target)
                    return p;
            }

            // 2) suffix after the leaf: cover full blocks [blk+1 .. leaf_count-1]
            const int excess_at_Rb = excess(Rb);
            int need = target - excess_at_Rb; // target expressed in coordinates of the current node start
            if (blk + 1 <= (leaf_count ? leaf_count - 1 : 0))
            {
                const auto nodes = cover_blocks(blk + 1, leaf_count - 1);
                size_t base = (blk + 1) * block_bits;
                for (const size_t &v : nodes)
                {
                    if (need == 0)
                        return base;
                    if (node_min_prefix_excess[v] <= need && need <= node_max_prefix_excess[v])
                        return descend_fwd(v, need, base);
                    need -= node_total_excess[v];
                    base += segment_size_bits[v];
                }
            }
            return npos;
        }

        /**
         * @brief Backward search: last position p ≤ i where excess(p) = excess(i) + d.
         * @details Scans inside the leaf to the left, then climbs to examine left siblings.
         *          Returns npos if no such position exists.
         */
        size_t bwdsearch(const size_t &i, const int &d) const
        {
            if (i > num_bits || i == 0)
                return npos;
            const int start_excess = excess(i);
            const int target = start_excess + d;

            // 1) scan inside the block
            const size_t blk = block_of(i - 1);
            const size_t Lb = blk * block_bits;
            int cur = start_excess;
            for (size_t p = i; p > Lb;)
            {
                --p;
                cur += bit(p) ? -1 : +1;
                if (cur == target)
                    return p;
                if (p == 0)
                    break;
            }
            const int excess_at_Lb = excess(Lb);
            if (Lb < i && excess_at_Lb == target)
                return Lb;

            // 2) climb up
            int need = target - excess_at_Lb;
            size_t v = leaf_index_of(Lb);
            size_t base = Lb;
            while (v > 1)
            {
                if (v & 1)
                {                                                        // v is the right child
                    const size_t sib = v ^ 1;                            // left sibling
                    const size_t border = base;                          // right border of the sibling (== start(v))
                    const int need_node = need + node_total_excess[sib]; // target in coordinates relative to the start of sib
                    const bool allow_rb = (border != i);                 // j must be < i

                    // try inside the sibling, but return only if a position is found
                    if (need_node == 0 || (node_min_prefix_excess[sib] <= need_node && need_node <= node_max_prefix_excess[sib]))
                    {
                        const size_t ans = descend_bwd(sib, border - segment_size_bits[sib], need_node, border, allow_rb);
                        if (ans != npos)
                            return ans;
                    }
                    // junction between children is a separate branch (allowed only if < i)
                    if (need_node == node_total_excess[sib] && border < i)
                        return border;

                    // stepped over the sibling, shifted the zero point of the coordinates
                    need += node_total_excess[sib];
                    base -= segment_size_bits[sib];
                }
                v >>= 1;
            }
            return npos;
        }

        /**
         * @brief Position of the first minimum of excess on [i, j] (inclusive).
         * @return Position of first occurrence of minimum, or npos on invalid range.
         */
        size_t rmq_pos(const size_t &i, const size_t &j) const
        {
            if (i > j || j >= num_bits)
                return npos;

            const size_t blk_i = block_of(i);
            const size_t Lbi = blk_i * block_bits;
            const size_t Rbi = std::min(num_bits, Lbi + block_bits);
            const size_t blk_j = block_of(j);
            const size_t Lbj = blk_j * block_bits;

            int best_val = INT_MAX;
            size_t best_pos = npos;
            size_t chosen_node = 0;
            int pref = 0, pref_at_choice = 0;

            // prefix
            int mn_first = INT_MAX;
            size_t first_pos = npos;
            const size_t end_first = std::min(j, (size_t)(Rbi ? Rbi - 1 : 0));
            if (i <= end_first)
            {
                first_min_value_pos8(i, end_first, mn_first, first_pos);
                pref = (int64_t)rank1_in_block(i, end_first + 1) * 2 - int64_t(end_first + 1 - i);
                best_val = mn_first;
                best_pos = first_pos;
                chosen_node = 0;
            }

            // middle
            const size_t leaf0 = first_leaf_index();
            if (blk_i + 1 <= blk_j - 1)
            {
                size_t l = leaf0 + (blk_i + 1);
                size_t r = leaf0 + (blk_j - 1);
                size_t Rnodes[64];
                int rn = 0;

                while (l <= r)
                {
                    if (l & 1)
                    {
                        const size_t v = l++;
                        const int cand = pref + node_min_prefix_excess[v];
                        if (cand < best_val)
                        {
                            best_val = cand;
                            best_pos = npos;
                            chosen_node = v;
                            pref_at_choice = pref;
                        }
                        pref += node_total_excess[v];
                    }
                    if ((r & 1) == 0)
                        Rnodes[rn++] = r--;
                    l >>= 1;
                    r >>= 1;
                }
                while (rn--)
                {
                    const size_t v = Rnodes[rn];
                    const int cand = pref + node_min_prefix_excess[v];
                    if (cand < best_val)
                    {
                        best_val = cand;
                        best_pos = npos;
                        chosen_node = v;
                        pref_at_choice = pref;
                    }
                    pref += node_total_excess[v];
                }
            }

            // tail
            if (blk_j != blk_i)
            {
                int mn_last;
                size_t last_pos;
                first_min_value_pos8(Lbj, j, mn_last, last_pos);
                const int cand = pref + mn_last;
                if (cand < best_val)
                {
                    best_val = cand;
                    best_pos = last_pos;
                    chosen_node = 0;
                }
            }

            if (best_pos != npos)
                return best_pos;

            return descend_first_min(chosen_node, best_val - pref_at_choice, node_base(chosen_node));
        }

        /**
         * @brief Value of the minimum prefix excess on [i, j] relative to i.
         * @details Equivalent to min_{t in [i..j]} (excess(t+1) - excess(i)).
         */
        int rmq_val(const size_t &i, const size_t &j) const
        {
            if (i > j || j >= num_bits)
                return 0;
            size_t p = rmq_pos(i, j);
            if (p == npos)
                return 0;
            return excess(p + 1) - excess(i);
        }

        /**
         * @brief Position of the first maximum of excess on [i, j] (inclusive).
         * @return Position of first occurrence of maximum, or npos on invalid range.
         */
        size_t rMq_pos(const size_t &i, const size_t &j) const
        {
            if (i > j || j >= num_bits)
                return npos;

            const size_t blk_i = block_of(i);
            const size_t Lbi = blk_i * block_bits;
            const size_t Rbi = std::min(num_bits, Lbi + block_bits);
            const size_t blk_j = block_of(j);
            const size_t Lbj = blk_j * block_bits;

            int best_val = INT_MIN;
            size_t best_pos = npos;
            size_t chosen_node = 0;
            int pref = 0, pref_at_choice = 0;

            // prefix
            int mx_first = INT_MIN;
            size_t first_pos = npos;
            const size_t end_first = std::min(j, (size_t)(Rbi ? Rbi - 1 : 0));
            if (i <= end_first)
            {
                first_max_value_pos8(i, end_first, mx_first, first_pos);
                pref = (int64_t)rank1_in_block(i, end_first + 1) * 2 - int64_t(end_first + 1 - i);
                best_val = mx_first;
                best_pos = first_pos;
                chosen_node = 0;
            }

            // middle
            const size_t leaf0 = first_leaf_index();
            if (blk_i + 1 <= blk_j - 1)
            {
                size_t l = leaf0 + (blk_i + 1);
                size_t r = leaf0 + (blk_j - 1);
                size_t Rnodes[64];
                int rn = 0;

                while (l <= r)
                {
                    if (l & 1)
                    {
                        const size_t v = l++;
                        const int cand = pref + node_max_prefix_excess[v];
                        if (cand > best_val)
                        {
                            best_val = cand;
                            best_pos = npos;
                            chosen_node = v;
                            pref_at_choice = pref;
                        }
                        pref += node_total_excess[v];
                    }
                    if ((r & 1) == 0)
                        Rnodes[rn++] = r--;
                    l >>= 1;
                    r >>= 1;
                }
                while (rn--)
                {
                    const size_t v = Rnodes[rn];
                    const int cand = pref + node_max_prefix_excess[v];
                    if (cand > best_val)
                    {
                        best_val = cand;
                        best_pos = npos;
                        chosen_node = v;
                        pref_at_choice = pref;
                    }
                    pref += node_total_excess[v];
                }
            }

            // tail
            if (blk_j != blk_i)
            {
                int mx_last;
                size_t last_pos;
                first_max_value_pos8(Lbj, j, mx_last, last_pos);
                const int cand = pref + mx_last;
                if (cand > best_val)
                {
                    best_val = cand;
                    best_pos = last_pos;
                    chosen_node = 0;
                }
            }

            if (best_pos != npos)
                return best_pos;

            return descend_first_max(chosen_node, best_val - pref_at_choice, node_base(chosen_node));
        }

        /**
         * @brief Value of the maximum prefix excess on [i, j] relative to i.
         */
        int rMq_val(const size_t &i, const size_t &j) const
        {
            if (i > j || j >= num_bits)
                return 0;
            size_t p = rMq_pos(i, j);
            if (p == npos)
                return 0;
            return excess(p + 1) - excess(i);
        }

        /**
         * @brief How many times the minimum prefix excess occurs on [i, j].
         */
        size_t mincount(const size_t &i, const size_t &j) const
        {
            if (i > j || j >= num_bits)
                return 0;

            const size_t blk_i = block_of(i);
            const size_t Lbi = blk_i * block_bits;
            const size_t Rbi = std::min(num_bits, Lbi + block_bits);
            const size_t blk_j = block_of(j);
            const size_t Lbj = blk_j * block_bits;

            int best_val = INT_MAX;
            size_t cnt = 0;
            int pref = 0;

            // first chunk
            {
                int cur = 0, mn = INT_MAX, c = 0;
                const size_t end = std::min(j, Rbi - 1);
                for (size_t p = i; p <= end; ++p)
                {
                    cur += bit(p) ? +1 : -1;
                    if (cur < mn)
                    {
                        mn = cur;
                        c = 1;
                    }
                    else if (cur == mn)
                    {
                        ++c;
                    }
                }
                best_val = mn;
                cnt = c;
                pref = cur; // offset toward the middle
            }

            // middle
            if (blk_i + 1 <= blk_j - 1)
            {
                const auto mids = cover_blocks(blk_i + 1, blk_j - 1);
                for (const size_t &v : mids)
                {
                    const int cand = pref + node_min_prefix_excess[v];
                    if (cand < best_val)
                    {
                        best_val = cand;
                        cnt = node_min_count[v];
                    }
                    else if (cand == best_val)
                    {
                        cnt += node_min_count[v];
                    }
                    pref += node_total_excess[v];
                }
            }

            // last chunk
            if (blk_j != blk_i)
            {
                int cur = 0, mn = INT_MAX, c = 0;
                for (size_t p = Lbj; p <= j; ++p)
                {
                    cur += bit(p) ? +1 : -1;
                    if (cur < mn)
                    {
                        mn = cur;
                        c = 1;
                    }
                    else if (cur == mn)
                    {
                        ++c;
                    }
                }
                const int cand = pref + mn;
                if (cand < best_val)
                {
                    best_val = cand;
                    cnt = c;
                }
                else if (cand == best_val)
                {
                    cnt += c;
                }
            }
            return cnt;
        }

        /**
         * @brief Position of the q-th (1-based) occurrence of the minimum on [i, j].
         * @return Position or npos if q exceeds the number of minima.
         */
        size_t minselect(const size_t &i, const size_t &j, size_t q) const
        {
            if (i > j || j >= num_bits || q == 0)
                return npos;

            const size_t blk_i = block_of(i);
            const size_t Lbi = blk_i * block_bits;
            const size_t Rbi = std::min(num_bits, Lbi + block_bits);
            const size_t blk_j = block_of(j);
            const size_t Lbj = blk_j * block_bits;

            // prefix
            const size_t end_first = std::min(j, Rbi - 1);
            int cur_first = 0, mn_first = 0;
            uint32_t c_first = 0;

            if (i <= end_first)
            {
                scan_range_min_count8(i, end_first, cur_first, mn_first, c_first);
            }
            else
            {
                cur_first = 0;
                mn_first = INT_MAX;
                c_first = 0;
            }

            int best_val = (mn_first == INT_MAX ? INT_MAX : mn_first);
            size_t total_cnt = (mn_first == INT_MAX ? 0u : (size_t)c_first);
            int pref = cur_first; // offset for middle

            const size_t leaf0 = first_leaf_index();
            size_t l = leaf0 + blk_i + 1;
            size_t r = leaf0 + blk_j - 1;
            size_t Rnodes[64];
            int rn = 0;

            // middle
            if (blk_i + 1 <= blk_j - 1)
            {
                while (l <= r)
                {
                    if (l & 1)
                    {
                        const int cand = pref + node_min_prefix_excess[l];
                        if (cand < best_val)
                        {
                            best_val = cand;
                            total_cnt = node_min_count[l];
                        }
                        else if (cand == best_val)
                        {
                            total_cnt += node_min_count[l];
                        }
                        pref += node_total_excess[l++];
                    }
                    if ((r & 1) == 0)
                        Rnodes[rn++] = r--;
                    l >>= 1;
                    r >>= 1;
                }
                while (rn--)
                {
                    const size_t v = Rnodes[rn];
                    const int cand = pref + node_min_prefix_excess[v];
                    if (cand < best_val)
                    {
                        best_val = cand;
                        total_cnt = node_min_count[v];
                    }
                    else if (cand == best_val)
                    {
                        total_cnt += node_min_count[v];
                    }
                    pref += node_total_excess[v];
                }
            }

            // tail
            int cur_last = 0, mn_last = INT_MAX;
            uint32_t c_last = 0;
            if (blk_j != blk_i)
            {
                scan_range_min_count8(Lbj, j, cur_last, mn_last, c_last);
                const int cand = pref + mn_last;
                if (cand < best_val)
                {
                    best_val = cand;
                    total_cnt = c_last;
                }
                else if (cand == best_val)
                {
                    total_cnt += c_last;
                }
            }

            if (q > total_cnt)
                return npos;

            // prefix
            if (mn_first == best_val && c_first)
            {
                if (q <= c_first)
                    return qth_min_in_block(i, end_first, q);
                q -= c_first;
            }

            // middle
            pref = cur_first;
            if (blk_i + 1 <= blk_j - 1)
            {
                l = leaf0 + (blk_i + 1);
                r = leaf0 + (blk_j - 1);
                rn = 0;
                while (l <= r)
                {
                    if (l & 1)
                    {
                        const size_t v = l++;
                        const int cand = pref + node_min_prefix_excess[v];
                        if (cand == best_val)
                        {
                            if (q <= node_min_count[v])
                            {
                                return descend_qth_min(v, best_val - pref, q, node_base(v));
                            }
                            q -= node_min_count[v];
                        }
                        pref += node_total_excess[v];
                    }
                    if (!(r & 1))
                        Rnodes[rn++] = r--;
                    l >>= 1;
                    r >>= 1;
                }
                while (rn--)
                {
                    const size_t v = Rnodes[rn];
                    const int cand = pref + node_min_prefix_excess[v];
                    if (cand == best_val)
                    {
                        if (q <= node_min_count[v])
                        {
                            return descend_qth_min(v, best_val - pref, q, node_base(v));
                        }
                        q -= node_min_count[v];
                    }
                    pref += node_total_excess[v];
                }
            }

            // tail
            if (blk_j != blk_i && (pref + mn_last) == best_val)
            {
                return qth_min_in_block(Lbj, j, q);
            }

            return npos;
        }

        // ----- parentheses navigation (BP) -----

        /**
         * @brief close(i): matching ')' for '(' at i.
         * @return Position of matching ')', or npos.
         */
        inline size_t close(const size_t &i) const
        {
            if (i >= num_bits)
                return npos;
            return fwdsearch(i, -1);
        }

        /**
         * @brief open(i): matching '(' for ')' at i.
         * @return Position of matching '(', or npos.
         */
        inline size_t open(const size_t &i) const
        {
            // bwdsearch allows i in [1..num_bits]
            if (i == 0 || i > num_bits)
                return npos;
            const size_t r = bwdsearch(i, 0);
            return (r == npos ? npos : r + 1);
        }

        /**
         * @brief enclose(i): opening '(' that strictly encloses position i.
         * @return Position of enclosing '(', or npos.
         */
        inline size_t enclose(const size_t &i) const
        {
            if (i == 0 || i > num_bits)
                return npos;
            const size_t r = bwdsearch(i, -2);
            return (r == npos ? npos : r + 1);
        }

    private:
        /**
         * @brief Count "10" occurrences inside a 64-bit slice of given logical length.
         * @details Only positions fully inside the slice are counted.
         */
        static inline size_t pop10_in_slice64(const std::uint64_t &slice, const int &len) noexcept
        {
            if (len <= 1)
                return 0;
            std::uint64_t P = slice & ~(slice >> 1); // candidates for "10"
            if (len < 64)
                P &= ((std::uint64_t(1) << (len - 1)) - 1);
            else
                P &= 0x7FFFFFFFFFFFFFFFull;
            return (size_t)POPCNT(P);
        }

        /**
         * @brief Rank of ones within [Lb, Rb).
         * @details Works on word boundaries; Rb may equal Lb.
         */
        size_t rank1_in_block(const size_t &Lb, const size_t &Rb) const noexcept
        {
            if (Rb <= Lb)
                return 0;
            size_t w_l = Lb >> 6;
            const size_t w_r = Rb >> 6;
            size_t off_l = Lb & 63;
            const size_t off_r = Rb & 63;
            size_t cnt = 0;
            if (w_l == w_r)
            {
                const std::uint64_t mask = ((off_r == 0) ? 0 : ((std::uint64_t(1) << off_r) - 1)) & (~std::uint64_t(0) << off_l);
                return (size_t)POPCNT(bits[w_l] & mask);
            }
            if (off_l)
            {
                cnt += (size_t)POPCNT(bits[w_l] & (~std::uint64_t(0) << off_l));
                ++w_l;
            }
            while (w_l < w_r)
            {
                cnt += (size_t)POPCNT(bits[w_l]);
                ++w_l;
            }
            if (off_r)
                cnt += (size_t)POPCNT(bits[w_r] & ((std::uint64_t(1) << off_r) - 1));
            return cnt;
        }

        /**
         * @brief Count "10" starts within [Lb, Rb).
         * @details Accounts for cross-word boundaries.
         */
        size_t rr_in_block(const size_t &Lb, const size_t &Rb) const noexcept
        {
            if (Rb <= Lb + 1)
                return 0;
            size_t w_l = Lb >> 6;
            const size_t w_r = (Rb - 1) >> 6;
            const int off_l = Lb & 63;
            const int off_r = (Rb - 1) & 63;
            size_t cnt = 0;

            if (w_l == w_r)
            {
                const int len = off_r - off_l + 1;
                const std::uint64_t slice = bits[w_l] >> off_l;
                return pop10_in_slice64(slice, len);
            }

            // prefix word
            {
                const int len = 64 - off_l;
                const std::uint64_t slice = bits[w_l] >> off_l;
                cnt += pop10_in_slice64(slice, len);
            }
            // full interior words
            for (size_t w = w_l + 1; w < w_r; ++w)
            {
                const std::uint64_t x = bits[w];
                cnt += pop10_in_slice64(x, 64);
            }
            // suffix word
            {
                const int len = off_r + 1;
                const std::uint64_t mask = (len == 64) ? ~std::uint64_t(0) : ((std::uint64_t(1) << len) - 1);
                const std::uint64_t slice = bits[w_r] & mask;
                cnt += pop10_in_slice64(slice, len);
            }
            // cross-word boundaries (bit 63 of w and bit 0 of w+1)
            for (size_t w = w_l; w < w_r; ++w)
            {
                if (((bits[w] >> 63) & 1u) && ((bits[w + 1] & 1u) == 0))
                    ++cnt;
            }
            return cnt;
        }

        /**
         * @brief 1-based select of k-th "10" within [Lb, Rb).
         * @return Position or npos.
         */
        size_t select10_in_block(const size_t &Lb, const size_t &Rb, size_t k) const noexcept
        {
            if (Rb <= Lb + 1)
                return npos;
            size_t w_l = Lb >> 6;
            const size_t w_r = (Rb - 1) >> 6;
            const int off_l = Lb & 63;
            const int off_r = (Rb - 1) & 63;

            const auto select_in_masked_slice = [&](const std::uint64_t &slice, const int &len, const size_t &kk) noexcept -> int
            {
                if (len <= 1)
                    return -1;
                std::uint64_t P = slice & ~(slice >> 1);
                if (len < 64)
                    P &= ((std::uint64_t(1) << (len - 1)) - 1);
                else
                    P &= 0x7FFFFFFFFFFFFFFFull;
                return select_in_word(P, kk);
            };

            if (w_l == w_r)
            {
                const int len = off_r - off_l + 1;
                const std::uint64_t slice = bits[w_l] >> off_l;
                const int off = select_in_masked_slice(slice, len, k);
                return off >= 0 ? (Lb + (size_t)off) : npos;
            }

            // prefix word
            {
                const int len = 64 - off_l;
                const std::uint64_t slice = bits[w_l] >> off_l;
                std::uint64_t P = slice & ~(slice >> 1);
                P &= ((std::uint64_t(1) << (len - 1)) - 1);
                const int c = POPCNT(P);
                if (k <= (size_t)c)
                {
                    const int off = select_in_masked_slice(slice, len, k);
                    return Lb + (size_t)off;
                }
                k -= c;
            }

            // walk interior boundaries and words
            for (size_t w = w_l; w + 1 < w_r; ++w)
            {
                // boundary between w and w+1
                if (((bits[w] >> 63) & 1u) && ((bits[w + 1] & 1u) == 0))
                {
                    if (--k == 0)
                        return (w << 6) + 63;
                }
                // full word w+1 (positions 0..62)
                const std::uint64_t x = bits[w + 1];
                const std::uint64_t P = (x & ~(x >> 1)) & 0x7FFFFFFFFFFFFFFFull;
                const int c = POPCNT(P);
                if (k <= (size_t)c)
                {
                    const int off = select_in_word(P, k);
                    if (off == -1)
                        return npos;
                    return ((w + 1) << 6) + (size_t)off;
                }
                k -= c;
            }

            // boundary (w_r-1, w_r)
            if (((bits[w_r - 1] >> 63) & 1u) && ((bits[w_r] & 1u) == 0))
            {
                if (--k == 0)
                    return ((w_r - 1) << 6) + 63;
            }

            // suffix word w_r: [0..off_r]
            {
                const int len = off_r + 1;
                const std::uint64_t mask = (len == 64) ? ~std::uint64_t(0) : ((std::uint64_t(1) << len) - 1);
                const std::uint64_t slice = bits[w_r] & mask;
                const int off = select_in_masked_slice(slice, len, k);
                if (off >= 0)
                    return (w_r << 6) + (size_t)off;
            }
            return npos;
        }

        struct ByteAgg
        {
            int8_t excess_total;     // total excess for the byte
            int8_t min_prefix;       // minimum prefix within the byte (from 0)
            int8_t max_prefix;       // maximum prefix within the byte (from 0)
            uint8_t min_count;       // number of positions attaining the minimum in the byte
            uint8_t pattern10_count; // number of "10" patterns inside the byte
            uint8_t first_bit;       // first bit (LSB)
            uint8_t last_bit;        // last bit (MSB)
            uint8_t pos_first_min;   // pos of first minimum in this byte
            uint8_t pos_first_max;   // pos of first maximum in this byte
        };

        /**
         * @brief Returns the static lookup table of byte aggregates.
         */
        static inline const std::array<ByteAgg, 256> &LUT8() noexcept
        {
            static const std::array<ByteAgg, 256> T = []
            {
                std::array<ByteAgg, 256> t{};
                for (int b = 0; b < 256; ++b)
                {
                    int cur = 0, mn = INT_MAX, mx = INT_MIN, cnt = 0, rrc = 0;
                    int pm = 0, pM = 0;
                    const auto get = [&](const int &k)
                    { return (b >> k) & 1; }; // LSB-first
                    for (int k = 0; k < 8; ++k)
                    {
                        int bit = get(k);
                        if (k + 1 < 8 && bit && get(k + 1) == 0)
                            ++rrc;
                        cur += bit ? +1 : -1;
                        if (cur < mn)
                        {
                            mn = cur;
                            cnt = 1;
                            pm = k;
                        }
                        else if (cur == mn)
                        {
                            ++cnt;
                        }
                        if (cur > mx)
                        {
                            mx = cur;
                            pM = k;
                        }
                    }
                    ByteAgg a{};
                    a.excess_total = cur;
                    a.min_prefix = (mn == INT_MAX ? 0 : mn);
                    a.max_prefix = (mx == INT_MIN ? 0 : mx);
                    a.min_count = cnt;
                    a.pattern10_count = rrc;
                    a.first_bit = get(0);
                    a.last_bit = get(7);
                    a.pos_first_min = pm;
                    a.pos_first_max = pM;
                    t[b] = a;
                }
                return t;
            }();
            return T;
        }

        /**
         * @brief Extract 8 bits starting at position pos (LSB-first across words).
         */
        inline uint8_t get_byte(const size_t &pos) const noexcept
        {
            const size_t w = pos >> 6;
            const size_t off = pos & 63;
            const std::uint64_t lo = bits[w] >> off;
            if (off <= 56)
                return uint8_t(lo & 0xFFu);
            const std::uint64_t hi = (w + 1 < bits.size()) ? bits[w + 1] : 0;
            const std::uint64_t x = (lo | (hi << (64 - off))) & 0xFFu;
            return uint8_t(x);
        }

        /**
         * @brief Descend to the first (leftmost) maximum with node-relative prefix equal to d.
         * @param v Node index.
         * @param d Target prefix within node coordinates.
         * @param base Starting global position of node.
         * @return Position or npos.
         */
        size_t descend_first_max(size_t v, int d, size_t base) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1, Rc = Lc | 1;
                const int leftX = node_max_prefix_excess[Lc];
                const int rightX = node_total_excess[Lc] + node_max_prefix_excess[Rc];
                if (leftX >= rightX && leftX == d)
                {
                    v = Lc;
                }
                else if (rightX == d)
                {
                    base += segment_size_bits[Lc];
                    d -= node_total_excess[Lc];
                    v = Rc;
                }
                else
                {
                    return npos;
                }
            }

            const size_t Lb = base;
            const size_t Rb = std::min(base + segment_size_bits[v], num_bits);
            int mx;
            size_t pos;

            first_max_value_pos8(Lb, Rb ? (Rb - 1) : Lb, mx, pos);
            return (mx == d ? pos : npos);
        }

        /**
         * @brief Heap index of the first leaf node.
         */
        inline size_t first_leaf_index() const noexcept
        {
            return std::bit_ceil(std::max<size_t>(1, leaf_count));
        }

        /**
         * @brief Block index containing position i.
         */
        size_t block_of(const size_t &i) const noexcept { return i / block_bits; }

        /**
         * @brief Heap index of the leaf whose segment starts at block_start.
         */
        size_t leaf_index_of(const size_t &block_start) const noexcept
        {
            return first_leaf_index() + block_of(block_start);
        }

        /**
         * @brief Starting global position for node v (0-indexed).
         * @details Walk up to compute base for internal nodes; direct for leaves.
         */
        size_t node_base(size_t v) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            if (v >= leaf0)
                return (v - leaf0) * block_bits;

            size_t base = 0;
            for (; v > 1; v >>= 1)
                if (v & 1)
                    base += segment_size_bits[v - 1];
            return base;
        }

        /**
         * @brief Cover a range of whole blocks [a..b] (inclusive) with O(log n) maximal nodes.
         * @details Returns node indices in left-to-right order.
         */
        std::vector<size_t> cover_blocks(const size_t &a, const size_t &b) const
        {
            const size_t leaf0 = first_leaf_index();
            size_t l = leaf0 + a;
            size_t r = leaf0 + b;
            std::vector<size_t> Lnodes, Rnodes;
            while (l <= r)
            {
                if ((l & 1) == 1)
                    Lnodes.push_back(l++);
                if ((r & 1) == 0)
                    Rnodes.push_back(r--);
                l >>= 1;
                r >>= 1;
            }
            std::reverse(Rnodes.begin(), Rnodes.end());
            Lnodes.insert(Lnodes.end(), Rnodes.begin(), Rnodes.end());
            return Lnodes;
        }

        /**
         * @brief Descend for fwdsearch to find first position where relative prefix equals 'need'.
         */
        size_t descend_fwd(size_t v, int need, size_t base) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1;
                const size_t Rc = Lc | 1;
                if (node_min_prefix_excess[Lc] <= need && need <= node_max_prefix_excess[Lc])
                    v = Lc;
                else
                {
                    need -= node_total_excess[Lc];
                    base += segment_size_bits[Lc];
                    v = Rc;
                }
            }
            const size_t Rb = std::min(base + segment_size_bits[v], num_bits);
            int cur = 0;
            for (size_t p = base; p < Rb; ++p)
            {
                cur += bit(p) ? +1 : -1;
                if (cur == need)
                    return p;
            }
            return npos;
        }

        /**
         * @brief Descend for bwdsearch to return the rightmost solution.
         * @param v Current node.
         * @param base Left border of node v.
         * @param need Target relative prefix inside v.
         * @param right_border Global right limit (exclusive if !allow_rb).
         * @param allow_rb Whether right border is allowed to match.
         */
        size_t descend_bwd(size_t v, const size_t &base, const int &need, const size_t &right_border, const bool &allow_rb) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1;
                const size_t Rc = Lc | 1;
                const int need_r = need - node_total_excess[Lc];

                // 1) try the right child first (to capture the rightmost j)
                if (node_min_prefix_excess[Rc] <= need_r && need_r <= node_max_prefix_excess[Rc])
                {
                    const size_t ans = descend_bwd(Rc, base + segment_size_bits[Lc], need_r, right_border, allow_rb);
                    if (ans != npos)
                        return ans;
                }

                // 2) junction between children (end of the left child)
                const size_t j_border = base + segment_size_bits[Lc];
                if (need == node_total_excess[Lc] && (j_border < right_border || allow_rb))
                {
                    return j_border;
                }

                // 3) can we move left within the range?
                if (node_min_prefix_excess[Lc] <= need && need <= node_max_prefix_excess[Lc])
                {
                    v = Lc;
                    continue;
                }

                // None of (1)-(3) worked. The only possible point is the left border of the node.
                if (need == 0 && (base < right_border || allow_rb))
                {
                    return base;
                }

                return npos;
            }

            const size_t Lb = base;
            const size_t Rb = std::min(base + segment_size_bits[v], num_bits);
            const size_t RB = std::min(right_border, Rb);

            int cur = 0;
            for (size_t p = Lb; p < RB; ++p)
                cur += bit(p) ? +1 : -1;

            if (allow_rb && cur == need)
                return RB;

            for (size_t p = RB; p > Lb;)
            {
                --p;
                cur += bit(p) ? -1 : +1;
                if (cur == need)
                    return p;
                if (p == 0)
                    break;
            }

            if ((Lb < right_border || allow_rb) && cur == need)
                return Lb;

            return npos;
        }

        /**
         * @brief Descend to find first position where node-relative prefix equals d (minimum).
         */
        size_t descend_first_min(size_t v, int d, size_t base) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1, Rc = Lc | 1;
                const int leftm = node_min_prefix_excess[Lc];
                const int rightm = node_total_excess[Lc] + node_min_prefix_excess[Rc];
                if (leftm <= rightm && leftm == d)
                {
                    v = Lc;
                }
                else if (rightm == d)
                {
                    base += segment_size_bits[Lc];
                    d -= node_total_excess[Lc];
                    v = Rc;
                }
                else
                {
                    return npos;
                }
            }

            const size_t Lb = base;
            const size_t Rb = std::min(base + segment_size_bits[v], num_bits);
            int mn;
            size_t pos;

            first_min_value_pos8(Lb, Rb ? (Rb - 1) : Lb, mn, pos);
            return (mn == d ? pos : npos);
        }

        /**
         * @brief Descend to find the q-th minimum (1-based) where node-relative prefix equals d.
         */
        size_t descend_qth_min(size_t v, int d, size_t q, size_t base) const noexcept
        {
            const size_t leaf0 = first_leaf_index();
            while (v < leaf0)
            {
                const size_t Lc = v << 1;
                const size_t Rc = Lc | 1;
                const int leftm = node_min_prefix_excess[Lc];
                const int rightm = node_total_excess[Lc] + node_min_prefix_excess[Rc];
                if (leftm == d)
                {
                    if (node_min_count[Lc] >= q)
                    {
                        v = Lc;
                        continue;
                    }
                    q -= node_min_count[Lc];
                }
                if (rightm == d)
                {
                    base += segment_size_bits[Lc];
                    d -= node_total_excess[Lc];
                    v = Rc;
                    continue;
                }
                return npos;
            }
            return qth_min_in_block(base, std::min(base + segment_size_bits[v], num_bits) - 1, q);
        }

        /**
         * @brief 1-based select of k-th '1' within [Lb, Rb).
         */
        size_t select1_in_block(const size_t &Lb, const size_t &Rb, size_t k) const noexcept
        {
            size_t w_l = Lb >> 6;
            const size_t w_r = (Rb >> 6);
            const size_t off_l = Lb & 63;
            const std::uint64_t mask_l = (off_l ? (~std::uint64_t(0) << off_l) : ~std::uint64_t(0));
            if (w_l == w_r)
            {
                const std::uint64_t w = bits[w_l] & mask_l & ((Rb & 63) ? ((std::uint64_t(1) << (Rb & 63)) - 1) : ~std::uint64_t(0));
                return Lb + select_in_word(w, k);
            }
            // prefix
            if (off_l)
            {
                const std::uint64_t w = bits[w_l] & mask_l;
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                    return Lb + select_in_word(w, k);
                k -= c;
                w_l++;
            }
            // full words
            while (w_l < w_r)
            {
                const std::uint64_t w = bits[w_l];
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                    return (w_l << 6) + select_in_word(w, k);
                k -= c;
                ++w_l;
            }
            // tail
            const size_t off_r = Rb & 63;
            if (off_r)
            {
                const std::uint64_t w = bits[w_l] & ((std::uint64_t(1) << off_r) - 1);
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                    return (w_l << 6) + select_in_word(w, k);
            }
            return npos;
        }

        /**
         * @brief 1-based select of k-th '0' within [Lb, Rb).
         */
        size_t select0_in_block(const size_t &Lb, const size_t &Rb, size_t k) const noexcept
        {
            if (Rb <= Lb)
                return npos;

            size_t w_l = Lb >> 6;
            const size_t w_r = Rb >> 6;
            const size_t off_l = Lb & 63;

            if (w_l == w_r)
            {
                const std::uint64_t mask_l = (off_l ? (~std::uint64_t(0) << off_l) : ~std::uint64_t(0));
                const std::uint64_t mask_r = ((Rb & 63) ? ((std::uint64_t(1) << (Rb & 63)) - 1) : ~std::uint64_t(0));
                const std::uint64_t w = (~bits[w_l]) & mask_l & mask_r;
                const int off = select_in_word(w, k);
                return (off >= 0) ? (Lb + (size_t)off) : npos;
            }

            // prefix
            if (off_l)
            {
                const std::uint64_t w = (~bits[w_l]) & (~std::uint64_t(0) << off_l);
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                {
                    const int off = select_in_word(w, k);
                    return (off >= 0) ? (Lb + (size_t)off) : npos;
                }
                k -= c;
                ++w_l;
            }

            // full words
            while (w_l < w_r)
            {
                const std::uint64_t w = ~bits[w_l];
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                {
                    const int off = select_in_word(w, k);
                    return (off >= 0) ? ((w_l << 6) + (size_t)off) : npos;
                }
                k -= c;
                ++w_l;
            }

            // tail
            const size_t off_r = Rb & 63;
            if (off_r)
            {
                const std::uint64_t w = (~bits[w_l]) & ((std::uint64_t(1) << off_r) - 1);
                const int c = POPCNT(w);
                if (k <= (size_t)c)
                {
                    const int off = select_in_word(w, k);
                    return (off >= 0) ? ((w_l << 6) + (size_t)off) : npos;
                }
            }
            return npos;
        }

        /**
         * @brief 1-based select of k-th set bit inside a 64-bit word.
         * @return Bit index [0..63] or −1 if not found.
         */
        static inline int select_in_word(std::uint64_t w, size_t k) noexcept
        {
#ifdef __GNUC__
            while (w)
            {
                if (--k == 0)
                    return __builtin_ctzll(w);
                w &= (w - 1);
            }
            return -1;
#else
            for (int i = 0; i < 64; ++i)
                if ((w >> i) & 1u)
                    if (--k == 0)
                        return i;
            return -1;
#endif
        }

        /**
         * @brief Ceil division of positive integers.
         */
        static inline size_t ceil_div(const size_t &a, const size_t &b) noexcept
        {
            return (a + b - 1) / b;
        }

        /**
         * @brief Number of node slots needed for @p Nbits with leaf size @p Bpow2.
         * @details Includes internal node space (heap layout) plus leaves.
         */
        static inline size_t nodeslots_for(const size_t &Nbits, const size_t &Bpow2) noexcept
        {
            if (Nbits == 0)
                return 0;
            size_t leaf_count = ceil_div(Nbits, Bpow2);
            return std::bit_ceil(std::max<size_t>(1, leaf_count)) + leaf_count;
        }

        /**
         * @brief Auxiliary overhead in bytes/bitvector bytes for given parameters.
         */
        static inline float overhead_for(const size_t &Nbits, const size_t &Bpow2) noexcept
        {
            static constexpr size_t AUX_SLOT_BYTES = sizeof(uint32_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint8_t) + sizeof(uint8_t);

            size_t bb = ceil_div(Nbits, 64) * 8;
            if (bb == 0)
                return 0;
            size_t slots = nodeslots_for(Nbits, Bpow2);
            size_t aux = slots * AUX_SLOT_BYTES;
            return ((float)aux) / ((float)bb);
        }

        /**
         * @brief Choose minimal block size (power of two) keeping overhead ≤ cap.
         * @details Returns 64 if cap < 0 (no constraint). Clamped to ≤16384 or Nbits.
         */
        static inline size_t choose_block_bits_for_overhead(const size_t &Nbits, const float &cap) noexcept
        {
            if (cap < 0.f)
                return 64;

            const size_t Bmax = std::min<size_t>(Nbits, 16384);
            size_t block_bits = 64;
            while (block_bits < Bmax)
            {
                if (overhead_for(Nbits, block_bits) <= cap)
                    break;
                block_bits <<= 1;
            }
            return block_bits;
        }

        /**
         * @brief Build internal structures from a 0/1 string.
         * @param leaf_block_bits Desired leaf size (0 = auto).
         * @param max_overhead Overhead cap (<0 = no cap).
         */
        void build_from_string(const std::string &s, const size_t &leaf_block_bits = 0, const float &max_overhead = -1.0)
        {
            num_bits = s.size();
            bits.assign((num_bits + 63) / 64, 0);
            for (size_t i = 0; i < num_bits; ++i)
                if (s[i] == '1')
                    set1(i);
            build(leaf_block_bits, max_overhead);
        }

        /**
         * @brief Build internal structures from 64-bit words.
         * @param words Words with LSB-first bits.
         * @param Nbits Number of valid bits.
         * @param leaf_block_bits Desired leaf size (0 = auto).
         * @param max_overhead Overhead cap (<0 = no cap).
         */
        void build_from_words(const std::vector<std::uint64_t> &words, const size_t &Nbits, const size_t &leaf_block_bits = 0, const float &max_overhead = -1.0)
        {
            bits = words;
            num_bits = Nbits;
            if (bits.size() * 64 < num_bits)
                bits.resize((num_bits + 63) / 64);
            build(leaf_block_bits, max_overhead);
        }

        /**
         * @brief Read bit at position i (LSB-first across words).
         */
        inline int bit(const size_t &i) const noexcept
        {
            return (bits[i >> 6] >> (i & 63)) & 1u;
        }

        /**
         * @brief Set bit at position i to 1.
         */
        inline void set1(const size_t &i) noexcept
        {
            bits[i >> 6] |= (std::uint64_t(1) << (i & 63));
        }

        /**
         * @brief Number of ones in node v computed from size and total excess.
         */
        inline uint32_t ones_in_node(const size_t &v) const noexcept
        {
            return ((int64_t)segment_size_bits[v] + (int64_t)node_total_excess[v]) >> 1;
        }

        /**
         * @brief Scan [l..r] inclusive, computing minimum value and how many times it is attained.
         * @details Uses 8-bit LUT for speed; outputs cur at end, mn, and count.
         */
        inline void scan_range_min_count8(size_t l, const size_t &r, int &cur, int &mn, uint32_t &cnt) const noexcept
        {
            cur = 0;
            mn = INT_MAX;
            cnt = 0;
            if (r < l)
            {
                mn = 0;
                return;
            }
            // to byte alignment
            while (l <= r && (l & 7))
            {
                cur += bit(l) ? +1 : -1;
                if (cur < mn)
                {
                    mn = cur;
                    cnt = 1;
                }
                else if (cur == mn)
                {
                    ++cnt;
                }
                ++l;
            }
            // full bytes
            const auto &T = LUT8();
            while (l + 7 <= r)
            {
                const auto &a = T[get_byte(l)];
                const int cand = cur + a.min_prefix;
                if (cand < mn)
                {
                    mn = cand;
                    cnt = a.min_count;
                }
                else if (cand == mn)
                {
                    cnt += a.min_count;
                }
                cur += a.excess_total;
                l += 8;
            }
            // tail
            while (l <= r)
            {
                cur += bit(l) ? +1 : -1;
                if (cur < mn)
                {
                    mn = cur;
                    cnt = 1;
                }
                else if (cur == mn)
                {
                    ++cnt;
                }
                ++l;
            }
            if (mn == INT_MAX)
            {
                mn = cnt = 0;
            }
        }

        /**
         * @brief Select q-th minimum (1-based) inside [l..r] inclusive in two 8-bit passes.
         * @details First pass finds global minimum, second selects the q-th position.
         */
        inline size_t qth_min_in_block(const size_t &l, const size_t &r, size_t q) const noexcept
        {
            if (r < l || q == 0)
                return npos;

            const auto &T = LUT8();

            int cur = 0, mn = INT_MAX;
            size_t p = l;

            while (p <= r && (p & 7))
            {
                cur += bit(p) ? +1 : -1;
                if (cur < mn)
                    mn = cur;
                ++p;
            }
            while (p + 7 <= r)
            {
                const auto &a = T[get_byte(p)];
                mn = std::min(mn, cur + a.min_prefix);
                cur += a.excess_total;
                p += 8;
            }
            while (p <= r)
            {
                cur += bit(p) ? +1 : -1;
                if (cur < mn)
                    mn = cur;
                ++p;
            }

            cur = 0;
            p = l;

            // to byte alignment
            while (p <= r && (p & 7))
            {
                cur += bit(p) ? +1 : -1;
                if (cur == mn)
                {
                    if (--q == 0)
                        return p;
                }
                ++p;
            }

            // full bytes
            while (p + 7 <= r)
            {
                const uint8_t b = get_byte(p);
                const auto &a = T[b];
                const int cand = cur + a.min_prefix;
                if (cand == mn)
                {
                    int s = 0;
                    for (int k = 0; k < 8; ++k)
                    {
                        s += ((b >> k) & 1u) ? +1 : -1;
                        if (s == a.min_prefix)
                        {
                            if (--q == 0)
                                return p + k;
                        }
                    }
                }
                cur += a.excess_total;
                p += 8;
            }

            // tail
            while (p <= r)
            {
                cur += bit(p) ? +1 : -1;
                if (cur == mn)
                {
                    if (--q == 0)
                        return p;
                }
                ++p;
            }

            return npos;
        }

        /**
         * @brief Find first minimum value and its first position in [l..r] inclusive using 8-bit LUT.
         * @param mn_out Output minimum value (0 if empty).
         * @param first_pos Output position of first minimum (npos if none).
         */
        inline void first_min_value_pos8(size_t l, const size_t &r, int &mn_out, size_t &first_pos) const noexcept
        {
            const auto &T = LUT8();
            int cur = 0;
            int mn = INT_MAX;
            first_pos = npos;

            // to byte allignment
            while (l <= r && (l & 7))
            {
                cur += bit(l) ? +1 : -1;
                if (cur < mn)
                {
                    mn = cur;
                    first_pos = l;
                }
                ++l;
            }

            // full bytes
            while (l + 7 <= r)
            {
                const auto &a = T[get_byte(l)];
                const int cand = cur + a.min_prefix;
                if (cand < mn)
                {
                    mn = cand;
                    first_pos = l + a.pos_first_min;
                }
                cur += a.excess_total;
                l += 8;
            }

            // tail
            while (l <= r)
            {
                cur += bit(l) ? +1 : -1;
                if (cur < mn)
                {
                    mn = cur;
                    first_pos = l;
                }
                ++l;
            }

            mn_out = (mn == INT_MAX ? 0 : mn);
        }

        /**
         * @brief Find first maximum value and its first position in [l..r] inclusive using 8-bit LUT.
         * @param mx_out Output maximum value (0 if empty).
         * @param first_pos Output position of first maximum (npos if none).
         */
        inline void first_max_value_pos8(size_t l, const size_t &r, int &mx_out, size_t &first_pos) const noexcept
        {
            const auto &T = LUT8();
            int cur = 0;
            int mx = INT_MIN;
            first_pos = npos;

            while (l <= r && (l & 7))
            {
                cur += bit(l) ? +1 : -1;
                if (cur > mx)
                {
                    mx = cur;
                    first_pos = l;
                }
                ++l;
            }

            while (l + 7 <= r)
            {
                const auto &a = T[get_byte(l)];
                const int cand = cur + a.max_prefix;
                if (cand > mx)
                {
                    mx = cand;
                    first_pos = l + a.pos_first_max;
                }
                cur += a.excess_total;
                l += 8;
            }

            while (l <= r)
            {
                cur += bit(l) ? +1 : -1;
                if (cur > mx)
                {
                    mx = cur;
                    first_pos = l;
                }
                ++l;
            }

            mx_out = (mx == INT_MIN ? 0 : mx);
        }

        /**
         * @brief Build the tree arrays and per-node aggregates.
         * @details Chooses block_bits honoring @p max_overhead or explicit @p leaf_block_bits,
         *          allocates arrays, fills leaves via LUT, and builds internal nodes bottom-up.
         */
        void build(const size_t &leaf_block_bits, const float &max_overhead)
        {
            // the lower clamp depends on the desired overhead fraction; otherwise use 64
            const size_t clamp_by_overhead = (max_overhead >= 0.0 ? choose_block_bits_for_overhead(num_bits, max_overhead) : size_t(64));

            // chosen block_bits: honor an explicit request, but not below clamp_by_overhead
            if (leaf_block_bits == 0)
                block_bits = std::max(clamp_by_overhead, std::bit_ceil<size_t>((num_bits <= 1) ? 1 : std::bit_width(num_bits - 1)));
            else
                block_bits = std::max(clamp_by_overhead, std::bit_ceil(std::max<size_t>(1, leaf_block_bits)));

#ifdef DEBUG
            // finalizes the achieved overhead percentage
            built_overhead = overhead_for(num_bits, block_bits);
#endif

            leaf_count = ceil_div(num_bits, block_bits);
            const size_t leaf0 = first_leaf_index();
            const size_t tree_size = leaf0 + leaf_count - 1;
            segment_size_bits.assign(tree_size + 1, 0);
            node_total_excess.assign(tree_size + 1, 0);
            node_min_prefix_excess.assign(tree_size + 1, 0);
            node_max_prefix_excess.assign(tree_size + 1, 0);
            node_min_count.assign(tree_size + 1, 0);
            node_pattern10_count.assign(tree_size + 1, 0);
            node_first_bit.assign(tree_size + 1, 0);
            node_last_bit.assign(tree_size + 1, 0);

            // leaves
            for (size_t k = 0; k < leaf_count; ++k)
            {
                const size_t v = leaf0 + k;
                const size_t Lb = k * block_bits;
                const size_t Rb = std::min(num_bits, Lb + block_bits);
                segment_size_bits[v] = Rb - Lb;

                if (Lb < Rb)
                {
                    node_first_bit[v] = bit(Lb);
                }

                const auto &T = LUT8();

                int cur = 0, mn = INT_MAX, mx = INT_MIN;
                uint32_t mn_cnt = 0;
                uint32_t rrc = 0;

                uint8_t prev_bit = 0;

                size_t p = Lb;

                // Full bytes
                while (p + 8 <= Rb)
                {
                    const uint8_t b = get_byte(p);
                    const auto &a = T[b];

                    // internal "10" inside the byte
                    rrc += a.pattern10_count;
                    // stitching across the boundary between the previous and current byte (within the segment)
                    if (prev_bit == 1 && a.first_bit == 0)
                        rrc++;

                    // prefix min/max accounting for the current offset
                    const int cand_min = cur + a.min_prefix;
                    if (cand_min < mn)
                    {
                        mn = cand_min;
                        mn_cnt = a.min_count;
                    }
                    else if (cand_min == mn)
                    {
                        mn_cnt += a.min_count;
                    }

                    mx = std::max(mx, cur + a.max_prefix);
                    cur += a.excess_total;
                    prev_bit = a.last_bit;
                    p += 8;
                }

                // Tail < 8 bits
                while (p < Rb)
                {
                    const uint8_t b = bit(p);
                    if (prev_bit == 1 && b == 0)
                        rrc++;
                    const int step = b ? +1 : -1;
                    cur += step;
                    if (cur < mn)
                    {
                        mn = cur;
                        mn_cnt = 1;
                    }
                    else if (cur == mn)
                    {
                        ++mn_cnt;
                    }
                    if (cur > mx)
                        mx = cur;

                    prev_bit = b;
                    ++p;
                }

                if (Lb < Rb)
                    node_last_bit[v] = prev_bit;

                node_total_excess[v] = cur;
                node_min_prefix_excess[v] = (segment_size_bits[v] == 0 ? 0 : mn);
                node_max_prefix_excess[v] = (segment_size_bits[v] == 0 ? 0 : mx);
                node_min_count[v] = mn_cnt;
                node_pattern10_count[v] = (uint32_t)rrc;
            }
            // internal nodes
            for (size_t v = leaf0 - 1; v >= 1; --v)
            {
                const size_t Lc = v << 1;
                const size_t Rc = Lc | 1;
                const bool has_l = (Lc <= tree_size) && segment_size_bits[Lc];
                const bool has_r = (Rc <= tree_size) && segment_size_bits[Rc];
                if (!has_l && !has_r)
                {
                    segment_size_bits[v] = 0;
                    continue;
                }
                if (has_l && !has_r)
                {
                    segment_size_bits[v] = segment_size_bits[Lc];
                    node_total_excess[v] = node_total_excess[Lc];
                    node_min_prefix_excess[v] = node_min_prefix_excess[Lc];
                    node_max_prefix_excess[v] = node_max_prefix_excess[Lc];
                    node_min_count[v] = node_min_count[Lc];
                    node_pattern10_count[v] = node_pattern10_count[Lc];
                    node_first_bit[v] = node_first_bit[Lc];
                    node_last_bit[v] = node_last_bit[Lc];
                }
                else if (!has_l && has_r)
                {
                    segment_size_bits[v] = segment_size_bits[Rc];
                    node_total_excess[v] = node_total_excess[Rc];
                    node_min_prefix_excess[v] = node_min_prefix_excess[Rc];
                    node_max_prefix_excess[v] = node_max_prefix_excess[Rc];
                    node_min_count[v] = node_min_count[Rc];
                    node_pattern10_count[v] = node_pattern10_count[Rc];
                    node_first_bit[v] = node_first_bit[Rc];
                    node_last_bit[v] = node_last_bit[Rc];
                }
                else
                {
                    segment_size_bits[v] = segment_size_bits[Lc] + segment_size_bits[Rc];
                    node_total_excess[v] = node_total_excess[Lc] + node_total_excess[Rc];
                    const int m_r = node_total_excess[Lc] + node_min_prefix_excess[Rc];
                    const int M_R = node_total_excess[Lc] + node_max_prefix_excess[Rc];
                    node_min_prefix_excess[v] = std::min(node_min_prefix_excess[Lc], m_r);
                    node_max_prefix_excess[v] = std::max(node_max_prefix_excess[Lc], M_R);
                    node_min_count[v] = (node_min_prefix_excess[Lc] == node_min_prefix_excess[v] ? node_min_count[Lc] : 0) + (m_r == node_min_prefix_excess[v] ? node_min_count[Rc] : 0);
                    node_pattern10_count[v] = node_pattern10_count[Lc] + node_pattern10_count[Rc] + ((node_last_bit[Lc] == 1 && node_first_bit[Rc] == 0) ? 1u : 0u);
                    node_first_bit[v] = node_first_bit[Lc];
                    node_last_bit[v] = node_last_bit[Rc];
                }
                if (v == 1)
                    break;
            }
        }
    };

} // namespace pixie
