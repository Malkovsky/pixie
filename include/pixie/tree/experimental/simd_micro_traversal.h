#pragma once

#include <pixie/experimental/microtree_simd.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/unpack.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

namespace pixie::tree::experimental {

/**
 * @brief Single-pass SIMD microtree inspected-DFS traversal for `BPTree`.
 *
 * @details Replaces the two-pass Four-Russians lookup decoder + 197-way
 * indirect handler table with on-the-fly packed 2-bit-pair SIMD zero-excess
 * masks (`bp_pair_zero_mask_64` / `bp_pair_zero_mask_128`). Subtrees up to
 * `MaxSimdBits` bits (`MaxSimdBits / 2 + 1` nodes) are decoded directly in
 * registers / small bitmask words without allocating per-node records.
 */
template <size_t MaxSimdBits, class Tree, class Visit, class Inspect>
class BpSimdMicroTraversal {
  static_assert(MaxSimdBits >= 128 && (MaxSimdBits % 128) == 0,
                "MaxSimdBits must be a positive multiple of 128.");
  static constexpr size_t kMaxPairs = MaxSimdBits / 2;
  static constexpr size_t kMaskWords = MaxSimdBits / 128;

  const Tree& tree_;
  std::span<const uint64_t> words_;
  size_t nbits_;
  Visit& visit_;
  Inspect& inspect_;
  size_t max_depth_;

  [[gnu::always_inline]] bool bit_at(size_t pos) const noexcept {
    return ((words_[pos >> 6] >> (pos & 63)) & 1ULL) != 0;
  }

 public:
  static constexpr size_t child_buffer_bytes = 2 * 32 * sizeof(size_t);
  size_t nodes = 0;
  size_t children = 0;

  static constexpr size_t scratch_bytes(size_t max_depth) noexcept {
    return max_depth * (child_buffer_bytes + kMaskWords * sizeof(uint64_t));
  }

  BpSimdMicroTraversal(const Tree& tree,
                       std::span<const uint64_t> words,
                       Visit& visit,
                       Inspect& inspect,
                       size_t max_depth) noexcept
      : tree_(tree),
        words_(words),
        nbits_(2 * tree.size()),
        visit_(visit),
        inspect_(inspect),
        max_depth_(max_depth) {}

  size_t subtree_nodes(TreeNode root) const {
    if (root.number == 0) {
      return tree_.size();
    }
    TreeNode cur = root;
    while (!tree_.is_root(cur) && tree_.is_last_child(cur)) {
      cur = tree_.parent(cur);
    }
    const size_t end_num =
        tree_.is_root(cur) ? tree_.size() : tree_.next_sibling(cur).number;
    return end_num - root.number;
  }

  void walk(TreeNode root) {
    if (tree_.empty()) {
      return;
    }
    const size_t sub_nodes = subtree_nodes(root);
    nodes += sub_nodes;
    children += sub_nodes - 1;
    visit_subtree(root, sub_nodes, 0);
  }

 private:
  [[gnu::noinline]] void visit_subtree(TreeNode u,
                                       size_t sub_nodes,
                                       size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    using namespace pixie::tree::experimental::microtree_simd;

    if (sub_nodes == 1) {
      emit(visit_, u, true);
      return;
    }
    if (sub_nodes == 2) {
      const TreeNode c0(u.number + 1, u.pos + 1);
      emit(visit_, u, true);
      inspect_(c0);
      emit(visit_, c0, true);
      emit(visit_, u, false);
      return;
    }
    if (sub_nodes == 3) {
      const TreeNode c0(u.number + 1, u.pos + 1);
      emit(visit_, u, true);
      if (!bit_at(u.pos + 2)) {
        // Star: 1 10 10 0
        const TreeNode c1(u.number + 2, u.pos + 3);
        inspect_(c0);
        inspect_(c1);
        emit(visit_, c0, true);
        emit(visit_, c1, true);
      } else {
        // Chain: 1 1 10 0 0
        const TreeNode gc(u.number + 2, u.pos + 2);
        inspect_(c0);
        emit(visit_, c0, true);
        inspect_(gc);
        emit(visit_, gc, true);
        emit(visit_, c0, false);
      }
      emit(visit_, u, false);
      return;
    }

    if (depth >= max_depth_) {
      walk_uncached(u);
      return;
    }

    const size_t pairs =
        sub_nodes - 1;  // Interior 2-bit pairs in [u.pos + 1, close(u))
    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + 1;

    if (pairs <= 32) {
      const uint64_t w = load_bits_64(words_, base_pos);
      // Fast check: all `pairs` children are 1-node leaves iff interior is
      // (10)^pairs.
      if (pairs <= 31) {
        const uint64_t interior_mask = (1ULL << (2 * pairs)) - 1ULL;
        if (((w ^ 0x5555555555555555ULL) & interior_mask) == 0) {
          emit(visit_, u, true);
          for (size_t j = 0; j < pairs; ++j) {
            inspect_(TreeNode(base_num + j, base_pos + 2 * j));
          }
          for (size_t j = 0; j < pairs; ++j) {
            emit(visit_, TreeNode(base_num + j, base_pos + 2 * j), true);
          }
          emit(visit_, u, false);
          return;
        }
      }
      uint32_t eq0 = bp_pair_zero_mask_64(w);
      if (pairs < 32) {
        eq0 &= (1u << pairs) - 1u;
      }
      replay_single_word(u, static_cast<uint64_t>(eq0), pairs, depth);
      return;
    }

    if (pairs <= 64) {
      uint64_t s[2];
      load_bits_128(words_, base_pos, s);
      uint64_t eq0 = bp_pair_zero_mask_128(s, 0, nullptr);
      if (pairs < 64) {
        eq0 &= (1ULL << pairs) - 1ULL;
      }
      replay_single_word(u, eq0, pairs, depth);
      return;
    }

    if (pairs <= kMaxPairs) {
      replay_multi_word(u, pairs, depth);
      return;
    }

    visit_large(u, sub_nodes, depth);
  }

  void replay_single_word(TreeNode u,
                          uint64_t eq0,
                          size_t pairs,
                          size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + 1;

    emit(visit_, u, true);
    inspect_(TreeNode(base_num, base_pos));
    uint64_t m_ins = eq0 & ~(1ULL << (pairs - 1));
    while (m_ins) {
      const size_t sp = static_cast<size_t>(std::countr_zero(m_ins)) + 1;
      m_ins &= m_ins - 1;
      inspect_(TreeNode(base_num + sp, base_pos + 2 * sp));
    }

    // Separate leaf children from non-leaf children so consecutive leaf runs
    // stream in a branch-predictable loop without coin-flip branches.
    const uint64_t leaf_ends = eq0 & ((eq0 << 1) | 1ULL);
    uint64_t nonleaf_ends = eq0 & ~leaf_ends;
    size_t cur_pair = 0;
    while (nonleaf_ends) {
      const unsigned end_nl =
          static_cast<unsigned>(std::countr_zero(nonleaf_ends));
      nonleaf_ends &= nonleaf_ends - 1;
      const uint64_t below = eq0 & ((1ULL << end_nl) - 1ULL);
      const size_t start_nl =
          below ? static_cast<size_t>(64 - std::countl_zero(below)) : 0;
      for (size_t p = cur_pair; p < start_nl; ++p) {
        emit(visit_, TreeNode(base_num + p, base_pos + 2 * p), true);
      }
      visit_subtree(TreeNode(base_num + start_nl, base_pos + 2 * start_nl),
                    end_nl - start_nl + 1, depth + 1);
      cur_pair = static_cast<size_t>(end_nl) + 1;
    }
    for (size_t p = cur_pair; p < pairs; ++p) {
      emit(visit_, TreeNode(base_num + p, base_pos + 2 * p), true);
    }
    emit(visit_, u, false);
  }

  void replay_multi_word(TreeNode u, size_t pairs, size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    using namespace pixie::tree::experimental::microtree_simd;

    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + 1;
    const size_t nw = (pairs + 63) / 64;
    std::array<uint64_t, kMaskWords> masks;

    int h = 0;
    for (size_t w = 0; w < nw; ++w) {
      uint64_t s[2];
      load_bits_128(words_, base_pos + w * 128, s);
      masks[w] = bp_pair_zero_mask_128(s, h, &h);
    }
    const unsigned rem = static_cast<unsigned>(pairs & 63);
    if (rem != 0) {
      masks[nw - 1] &= (1ULL << rem) - 1ULL;
    }

    emit(visit_, u, true);
    inspect_(TreeNode(base_num, base_pos));
    for (size_t w = 0; w < nw; ++w) {
      uint64_t m = masks[w];
      if (w + 1 == nw) {
        m &= ~(1ULL << ((pairs - 1) & 63));
      }
      const size_t pair_base = w * 64;
      while (m) {
        const size_t sp =
            pair_base + static_cast<size_t>(std::countr_zero(m)) + 1;
        m &= m - 1;
        inspect_(TreeNode(base_num + sp, base_pos + 2 * sp));
      }
    }

    size_t start_pair = 0;
    for (size_t w = 0; w < nw; ++w) {
      uint64_t m = masks[w];
      const size_t pair_base = w * 64;
      while (m) {
        const size_t end_pair =
            pair_base + static_cast<size_t>(std::countr_zero(m));
        m &= m - 1;
        visit_subtree(
            TreeNode(base_num + start_pair, base_pos + 2 * start_pair),
            end_pair - start_pair + 1, depth + 1);
        start_pair = end_pair + 1;
      }
    }
    emit(visit_, u, false);
  }

  [[gnu::noinline]] void visit_large(TreeNode u,
                                     size_t sub_nodes,
                                     size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;

    emit(visit_, u, true);
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    bool overflow = false;
    tree_.for_each_child(u, [&](TreeNode child) {
      inspect_(child);
      if (count < 32) {
        numbers[count] = child.number;
        positions[count] = child.pos;
        ++count;
      } else {
        overflow = true;
      }
    });

    const size_t subtree_end_num = u.number + sub_nodes;
    if (!overflow) {
      for (size_t j = 0; j < count; ++j) {
        const size_t next_num =
            (j + 1 < count) ? numbers[j + 1] : subtree_end_num;
        visit_subtree(TreeNode(numbers[j], positions[j]), next_num - numbers[j],
                      depth + 1);
      }
    } else {
      for (size_t j = 0; j + 1 < count; ++j) {
        visit_subtree(TreeNode(numbers[j], positions[j]),
                      numbers[j + 1] - numbers[j], depth + 1);
      }
      TreeNode pending(numbers[count - 1], positions[count - 1]);
      auto view = tree_.navigate(pending);
      while (!view.is_last_child()) {
        view.reset(view.next_sibling());
        const TreeNode cur = view.node();
        visit_subtree(pending, cur.number - pending.number, depth + 1);
        pending = cur;
      }
      visit_subtree(pending, subtree_end_num - pending.number, depth + 1);
    }
    emit(visit_, u, false);
  }

  [[gnu::noinline]] void walk_uncached(TreeNode root) {
    using pixie::tree::detail::tree_unpack::emit;

    auto view = tree_.navigate(root);
    bool descending = true;
    while (true) {
      const TreeNode cur = view.node();
      emit(visit_, cur, descending);
      if (descending) {
        tree_.for_each_child(cur, inspect_);
        if (!view.is_leaf()) {
          view.reset(view.first_child());
          continue;
        }
      }
      descending = false;
      if (cur.number == root.number) {
        return;
      }
      if (view.is_last_child()) {
        view.reset(view.parent());
      } else {
        view.reset(view.next_sibling());
        descending = true;
      }
    }
  }
};

/**
 * @brief Single-pass SIMD microtree inspected-DFS traversal for `DFUDSTree`.
 *
 * @details Uses unary degree scanning + `dfuds_record_lows_64` /
 * `dfuds_record_lows_128` SIMD strict-record-low masks to decode immediate
 * child boundaries for subtrees up to `MaxSimdBits` bits without Four-Russians
 * lookup tables or per-node record buffers.
 */
template <size_t MaxSimdBits, class Tree, class Visit, class Inspect>
class DfudsSimdMicroTraversal {
  static_assert(MaxSimdBits >= 128 && (MaxSimdBits % 128) == 0,
                "MaxSimdBits must be a positive multiple of 128.");
  static constexpr size_t kMaskWords = MaxSimdBits / 64;

  const Tree& tree_;
  std::span<const uint64_t> words_;
  size_t nbits_;
  Visit& visit_;
  Inspect& inspect_;
  size_t max_depth_;

  [[gnu::always_inline]] bool bit_at(size_t pos) const noexcept {
    return ((words_[pos >> 6] >> (pos & 63)) & 1ULL) != 0;
  }

  [[gnu::always_inline]] size_t read_degree(size_t pos) const noexcept {
    using namespace pixie::tree::experimental::microtree_simd;
    const uint64_t w0 = load_bits_64(words_, pos);
    const unsigned d0 = static_cast<unsigned>(std::countr_one(w0));
    if (d0 < 64) {
      return d0;
    }
    size_t d = 64;
    size_t cur = pos + 64;
    while (cur < nbits_) {
      const uint64_t w = load_bits_64(words_, cur);
      const unsigned run = static_cast<unsigned>(std::countr_one(w));
      d += run;
      if (run < 64) {
        break;
      }
      cur += 64;
    }
    return d;
  }

 public:
  static constexpr size_t child_buffer_bytes = 2 * 32 * sizeof(size_t);
  size_t nodes = 0;
  size_t children = 0;

  static constexpr size_t scratch_bytes(size_t max_depth) noexcept {
    return max_depth * (child_buffer_bytes + kMaskWords * sizeof(uint64_t));
  }

  DfudsSimdMicroTraversal(const Tree& tree,
                          std::span<const uint64_t> words,
                          Visit& visit,
                          Inspect& inspect,
                          size_t max_depth) noexcept
      : tree_(tree),
        words_(words),
        nbits_(tree.empty() ? 0 : 2 * tree.size() - 1),
        visit_(visit),
        inspect_(inspect),
        max_depth_(max_depth) {}

  size_t subtree_nodes(TreeNode root) const {
    if (root.number == 0) {
      return tree_.size();
    }
    TreeNode cur = root;
    while (!tree_.is_root(cur) && tree_.is_last_child(cur)) {
      cur = tree_.parent(cur);
    }
    const size_t end_num =
        tree_.is_root(cur) ? tree_.size() : tree_.next_sibling(cur).number;
    return end_num - root.number;
  }

  void walk(TreeNode root) {
    if (tree_.empty()) {
      return;
    }
    const size_t sub_nodes = subtree_nodes(root);
    nodes += sub_nodes;
    children += sub_nodes - 1;
    visit_subtree(root, sub_nodes, 0);
  }

 private:
  [[gnu::noinline]] void visit_subtree(TreeNode u,
                                       size_t sub_nodes,
                                       size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    using namespace pixie::tree::experimental::microtree_simd;

    if (sub_nodes == 1) {
      emit(visit_, u, true);
      return;
    }
    if (sub_nodes == 2) {
      const TreeNode c0(u.number + 1, u.pos + 2);
      emit(visit_, u, true);
      inspect_(c0);
      emit(visit_, c0, true);
      emit(visit_, u, false);
      return;
    }
    if (sub_nodes == 3) {
      emit(visit_, u, true);
      if (bit_at(u.pos + 1)) {
        // Star: 110 0 0
        const TreeNode c0(u.number + 1, u.pos + 3);
        const TreeNode c1(u.number + 2, u.pos + 4);
        inspect_(c0);
        inspect_(c1);
        emit(visit_, c0, true);
        emit(visit_, c1, true);
      } else {
        // Chain: 10 10 0
        const TreeNode c0(u.number + 1, u.pos + 2);
        const TreeNode gc(u.number + 2, u.pos + 4);
        inspect_(c0);
        emit(visit_, c0, true);
        inspect_(gc);
        emit(visit_, gc, true);
        emit(visit_, c0, false);
      }
      emit(visit_, u, false);
      return;
    }

    if (depth >= max_depth_) {
      walk_uncached(u);
      return;
    }

    const size_t d = read_degree(u.pos);
    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + d + 1;

    if (d == 1) {
      // Single child inherits the entire remaining subtree of `nodes - 1`
      // nodes.
      const TreeNode c0(base_num, base_pos);
      emit(visit_, u, true);
      inspect_(c0);
      visit_subtree(c0, sub_nodes - 1, depth + 1);
      emit(visit_, u, false);
      return;
    }

    if (sub_nodes == d + 1) {
      // Every child has size 1 (all `d` children are leaves), regardless of
      // `d`.
      emit(visit_, u, true);
      for (size_t j = 0; j < d; ++j) {
        inspect_(TreeNode(base_num + j, base_pos + j));
      }
      for (size_t j = 0; j < d; ++j) {
        emit(visit_, TreeNode(base_num + j, base_pos + j), true);
      }
      emit(visit_, u, false);
      return;
    }

    const size_t child_bits = 2 * sub_nodes - d - 2;

    if (child_bits <= 64) {
      uint64_t rl = dfuds_record_lows_64(load_bits_64(words_, base_pos));
      if (child_bits < 64) {
        rl &= (1ULL << child_bits) - 1ULL;
      }
      replay_single_word(u, rl, child_bits, depth);
      return;
    }

    if (child_bits <= 128) {
      uint64_t s[2];
      load_bits_128(words_, base_pos, s);
      Mask128 rl = dfuds_record_lows_128(s, 0, nullptr);
      if (child_bits < 128) {
        rl.hi &= (1ULL << (child_bits - 64)) - 1ULL;
      }
      std::array<uint64_t, 2> masks{rl.lo, rl.hi};
      replay_words(u, masks.data(), 2, child_bits, depth);
      return;
    }

    if (child_bits <= MaxSimdBits) {
      const size_t blocks128 = (child_bits + 127) / 128;
      const size_t nw = (child_bits + 63) / 64;
      std::array<uint64_t, kMaskWords> masks;
      int gap = 0;
      for (size_t b = 0; b < blocks128; ++b) {
        uint64_t s[2];
        load_bits_128(words_, base_pos + b * 128, s);
        const Mask128 rl = dfuds_record_lows_128(s, gap, &gap);
        masks[2 * b] = rl.lo;
        if (2 * b + 1 < nw) {
          masks[2 * b + 1] = rl.hi;
        }
      }
      const unsigned rem = static_cast<unsigned>(child_bits & 63);
      if (rem != 0) {
        masks[nw - 1] &= (1ULL << rem) - 1ULL;
      }
      replay_words(u, masks.data(), nw, child_bits, depth);
      return;
    }

    visit_large(u, sub_nodes, depth);
  }

  void replay_single_word(TreeNode u,
                          uint64_t rl,
                          size_t child_bits,
                          size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    // Note: in DFUDS `u.pos + d + 1` is `u.pos + popcount(rl) + 1`!
    const size_t d = static_cast<size_t>(std::popcount(rl));
    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + d + 1;

    emit(visit_, u, true);
    inspect_(TreeNode(base_num, base_pos));
    uint64_t m_ins = rl & ~(1ULL << (child_bits - 1));
    size_t child_idx = 1;
    while (m_ins) {
      const size_t sp = static_cast<size_t>(std::countr_zero(m_ins)) + 1;
      m_ins &= m_ins - 1;
      inspect_(TreeNode(base_num + (sp + child_idx) / 2, base_pos + sp));
      ++child_idx;
    }

    // Leaf children in DFUDS have bit length 1, so child r is a leaf iff
    // end_bit == start_bit (adjacent set bits in `rl`, or bit 0 for child 0).
    const uint64_t leaf_ends = rl & ((rl << 1) | 1ULL);
    uint64_t nonleaf_ends = rl & ~leaf_ends;
    size_t cur_bit = 0;
    size_t cur_num = base_num;
    while (nonleaf_ends) {
      const unsigned end_nl =
          static_cast<unsigned>(std::countr_zero(nonleaf_ends));
      nonleaf_ends &= nonleaf_ends - 1;
      const uint64_t below = rl & ((1ULL << end_nl) - 1ULL);
      const size_t start_nl =
          below ? static_cast<size_t>(64 - std::countl_zero(below)) : 0;
      const size_t leaves = start_nl - cur_bit;
      for (size_t j = 0; j < leaves; ++j) {
        emit(visit_, TreeNode(cur_num + j, base_pos + cur_bit + j), true);
      }
      cur_num += leaves;
      const size_t nl_nodes = (static_cast<size_t>(end_nl) - start_nl + 2) / 2;
      visit_subtree(TreeNode(cur_num, base_pos + start_nl), nl_nodes,
                    depth + 1);
      cur_num += nl_nodes;
      cur_bit = static_cast<size_t>(end_nl) + 1;
    }
    const size_t tail_leaves = child_bits - cur_bit;
    for (size_t j = 0; j < tail_leaves; ++j) {
      emit(visit_, TreeNode(cur_num + j, base_pos + cur_bit + j), true);
    }
    emit(visit_, u, false);
  }

  void replay_words(TreeNode u,
                    const uint64_t* masks,
                    size_t nw,
                    size_t child_bits,
                    size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;
    size_t d = 0;
    for (size_t w = 0; w < nw; ++w) {
      d += static_cast<size_t>(std::popcount(masks[w]));
    }
    const size_t base_num = u.number + 1;
    const size_t base_pos = u.pos + d + 1;

    emit(visit_, u, true);
    inspect_(TreeNode(base_num, base_pos));
    size_t child_idx = 1;
    for (size_t w = 0; w < nw; ++w) {
      uint64_t m = masks[w];
      if (w + 1 == nw) {
        m &= ~(1ULL << ((child_bits - 1) & 63));
      }
      const size_t bit_base = w * 64;
      while (m) {
        const size_t sp =
            bit_base + static_cast<size_t>(std::countr_zero(m)) + 1;
        m &= m - 1;
        inspect_(TreeNode(base_num + (sp + child_idx) / 2, base_pos + sp));
        ++child_idx;
      }
    }

    size_t start_bit = 0;
    size_t cur_num = base_num;
    for (size_t w = 0; w < nw; ++w) {
      uint64_t m = masks[w];
      const size_t bit_base = w * 64;
      while (m) {
        const size_t end_bit =
            bit_base + static_cast<size_t>(std::countr_zero(m));
        m &= m - 1;
        const size_t child_nodes = (end_bit - start_bit + 2) / 2;
        visit_subtree(TreeNode(cur_num, base_pos + start_bit), child_nodes,
                      depth + 1);
        cur_num += child_nodes;
        start_bit = end_bit + 1;
      }
    }
    emit(visit_, u, false);
  }

  [[gnu::noinline]] void visit_large(TreeNode u,
                                     size_t sub_nodes,
                                     size_t depth) {
    using pixie::tree::detail::tree_unpack::emit;

    emit(visit_, u, true);
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    bool overflow = false;
    tree_.for_each_child(u, [&](TreeNode child) {
      inspect_(child);
      if (count < 32) {
        numbers[count] = child.number;
        positions[count] = child.pos;
        ++count;
      } else {
        overflow = true;
      }
    });

    const size_t subtree_end_num = u.number + sub_nodes;
    if (!overflow) {
      for (size_t j = 0; j < count; ++j) {
        const size_t next_num =
            (j + 1 < count) ? numbers[j + 1] : subtree_end_num;
        visit_subtree(TreeNode(numbers[j], positions[j]), next_num - numbers[j],
                      depth + 1);
      }
    } else {
      for (size_t j = 0; j + 1 < count; ++j) {
        visit_subtree(TreeNode(numbers[j], positions[j]),
                      numbers[j + 1] - numbers[j], depth + 1);
      }
      TreeNode pending(numbers[count - 1], positions[count - 1]);
      auto view = tree_.navigate(pending);
      while (!view.is_last_child()) {
        view.reset(view.next_sibling());
        const TreeNode cur = view.node();
        visit_subtree(pending, cur.number - pending.number, depth + 1);
        pending = cur;
      }
      visit_subtree(pending, subtree_end_num - pending.number, depth + 1);
    }
    emit(visit_, u, false);
  }

  [[gnu::noinline]] void walk_uncached(TreeNode root) {
    using pixie::tree::detail::tree_unpack::emit;

    auto view = tree_.navigate(root);
    bool descending = true;
    while (true) {
      const TreeNode cur = view.node();
      emit(visit_, cur, descending);
      if (descending) {
        tree_.for_each_child(cur, inspect_);
        if (!view.is_leaf()) {
          view.reset(view.first_child());
          continue;
        }
      }
      descending = false;
      if (cur.number == root.number) {
        return;
      }
      if (view.is_last_child()) {
        view.reset(view.parent());
      } else {
        view.reset(view.next_sibling());
        descending = true;
      }
    }
  }
};

}  // namespace pixie::tree::experimental
