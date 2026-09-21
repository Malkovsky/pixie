#pragma once

#include <pixie/detail/sequence/sequence_tree_kernels.h>
#include <pixie/sequence_options.h>
#include <pixie/storage/aligned.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>

#ifdef PIXIE_SEQUENCE_TREE_TESTING
#include <unordered_set>
#include <vector>
#endif

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie {

template <class Index,
          std::size_t StorageBits,
          std::size_t Fanout,
          LengthLayout Layout>
class Permutation;
}  // namespace pixie

// Unsupported implementation shared by the two public sequence families.
namespace pixie::detail::sequence {

/**
 * @brief Nonallocating local contract for SequenceTree leaves.
 * @details Blocks default-construct empty and own their elements. Capacity is
 * positive and at most SIZE_MAX/2. Reads return value_type by value.
 * Moves, destruction, size(), and valid local mutations cannot throw. Rotation
 * accepts [left,right) within size() and reduces distance modulo its length.
 * redistribute requires distinct blocks and a feasible left count; it preserves
 * concatenated order, leaves that count on the left, and the rest on the right.
 * Neither mutation allocates. These semantic requirements supplement the
 * mechanically checked signatures. Throwing element mutations are unsupported.
 */
template <class B>
concept SequenceBlock =
    std::is_nothrow_default_constructible_v<B> &&
    std::is_nothrow_move_constructible_v<B> &&
    std::is_nothrow_move_assignable_v<B> && std::is_nothrow_destructible_v<B> &&
    requires(B& a, B& b, const B& c, std::size_t n) {
      typename B::value_type;
      requires B::capacity >= 1 &&
                       B::capacity <=
                           std::numeric_limits<std::size_t>::max() / 2;
      { c.size() } noexcept -> std::same_as<std::size_t>;
      { c[n] } -> std::same_as<typename B::value_type>;
      { a.rotate_left(n, n, n) } noexcept -> std::same_as<void>;
      { a.redistribute(b, n) } noexcept -> std::same_as<void>;
    };

/**
 * @brief Exclusively owned positional B+ split/join tree over bounded blocks.
 * @details All leaves have equal depth. Nonroot nodes have F/2..F children;
 * roots have at least two and collapse otherwise. Only the two exterior leaves
 * may be below half capacity. Counts are unsigned size_t, up to SIZE_MAX.
 * Access is O(F*h); split, join and rotation touch O(F*h) structural entries
 * and a constant number of bounded leaf payloads. No normal operation
 * enumerates leaves. Allocation preflight is O(h), bounded by the temporary
 * node deficit with dismantled nodes recycled; commit is nonallocating under
 * SequenceBlock's contract. No sharing, stable references, or concurrent
 * mutation is supported.
 *
 * Internal allocations contain F typed pointer slots and F-1 contiguous
 * measures plus a count, aligned to CacheLine. A level-aware move-only owner
 * controls roots and detached subtrees; occupied node slots own their children.
 * Leaves are separate aligned allocations containing a Block. Optional index
 * bias headers exist only when IndexBias is true. A decoded value equals its
 * stored field plus the nonnegative biases on its root-to-leaf path. Reads
 * accumulate without mutation; dismantling nodes pushes their bias to children,
 * and only redistributed leaves materialize fields. The permutation facade
 * checks the combined domain before attaching a bias, so every partial sum is
 * bounded by the final representable value. No merge-history chain is stored.
 * @tparam Block Nonthrowing bounded local block implementation.
 * @tparam Fanout Even fanout, at least four (4/8/16 are baseline candidates).
 * @tparam Layout Cumulative ends or individual lengths; last length is derived.
 * @tparam IndexBias Enable narrowly scoped lazy unsigned index rebasing. Block
 * must represent its value_type's full unsigned range and provide nonallocating
 * noexcept add_bias(uint64_t), requiring representable resulting values.
 */
template <SequenceBlock Block,
          std::size_t Fanout = 8,
          LengthLayout Layout = LengthLayout::cumulative,
          bool IndexBias = false>
class SequenceTree {
  static_assert(Fanout >= 4 && Fanout % 2 == 0);
  static_assert(!IndexBias ||
                (std::is_unsigned_v<typename Block::value_type> &&
                 !std::is_same_v<typename Block::value_type, bool>));
  static_assert(!IndexBias || requires(Block& block, std::uint64_t bias) {
    { block.add_bias(bias) } noexcept -> std::same_as<void>;
  });
  struct NoBias {};
  struct PendingBias {
    std::uint64_t value = 0;
  };
  using Bias = std::conditional_t<IndexBias, PendingBias, NoBias>;
  struct alignas(std::max(alignof(CacheLine), alignof(Block))) Leaf {
    Block block;
    [[no_unique_address]] Bias bias;
  };
  struct alignas(CacheLine) Node {
    std::array<void*, Fanout> children{};
    std::array<std::size_t, Fanout - 1> measures{};
    std::size_t count = 0;
    [[no_unique_address]] Bias bias;
  };
  template <class, std::size_t, std::size_t, LengthLayout>
  friend class ::pixie::Permutation;
  static constexpr std::size_t max_height =
      std::numeric_limits<std::size_t>::digits;
  static constexpr std::size_t minimum_leaf_size =
      (static_cast<std::size_t>(Block::capacity) + 1) / 2;

 public:
  /** @brief Block type accepted by the consuming factory. */
  using block_type = Block;
  /** @brief Immutable indexed result type. */
  using value_type = typename Block::value_type;
  /** @brief Maximum leaf element count. */
  static constexpr std::size_t block_capacity = Block::capacity;
  /** @brief Actual aligned leaf allocation bytes, including block padding. */
  static constexpr std::size_t block_storage_bytes = sizeof(Leaf);
  /** @brief Actual aligned internal allocation bytes. */
  static constexpr std::size_t node_storage_bytes = sizeof(Node);
  /** @brief Required leaf allocation alignment (payload may differ). */
  static constexpr std::size_t block_alignment = alignof(Leaf);
  /** @brief Required internal allocation alignment. */
  static constexpr std::size_t node_alignment = alignof(Node);

#ifdef PIXIE_SEQUENCE_TREE_TESTING
  /**
   * @brief Test-only transaction counters; absent without the testing macro.
   * @details All translation units using an instantiation must agree on
   * PIXIE_SEQUENCE_TREE_TESTING. Counts are thread-local per instantiation.
   * Payload mutations count calls, not elements. Allocation counts include
   * unused preflight spares; live counts include temporary/spare allocations.
   */
  struct TestCounters {
    std::size_t allocations = 0;
    std::size_t node_visits = 0;
    std::size_t child_transfers = 0;
    std::size_t payload_mutations = 0;
    std::size_t live_allocations = 0;
    std::size_t peak_allocations = 0;
    std::size_t live_bytes = 0;
    std::size_t peak_bytes = 0;
  };
  /** @brief Per-thread test counters, never present in normal builds. */
  inline static thread_local TestCounters test_counters{};
  /**
   * @brief Fail allocation after n successful attempts; -1 disables injection.
   * @param n Number of successful subsequent allocations before bad_alloc.
   */
  static void test_fail_after(std::ptrdiff_t n) noexcept { failure_ = n; }
  /** @brief Reset operation counters, retaining the current live allocation
   * count. */
  static void test_reset_counters() noexcept {
    const auto live = test_counters.live_allocations;
    const auto bytes = test_counters.live_bytes;
    test_counters = {};
    test_counters.live_allocations = test_counters.peak_allocations = live;
    test_counters.live_bytes = test_counters.peak_bytes = bytes;
  }
#endif

 private:
  static void allocation() {
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    if (failure_ == 0) {
      throw std::bad_alloc();
    }
    if (failure_ > 0) {
      --failure_;
    }
#endif
  }
  template <class T>
  static T* allocate() {
    allocation();
    auto* p = new T;
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    ++test_counters.allocations;
    ++test_counters.live_allocations;
    test_counters.peak_allocations = std::max(test_counters.peak_allocations,
                                              test_counters.live_allocations);
    test_counters.live_bytes += sizeof(T);
    test_counters.peak_bytes =
        std::max(test_counters.peak_bytes, test_counters.live_bytes);
#endif
    return p;
  }
  template <class T>
  static void dispose(T* p) noexcept {
    delete p;
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    --test_counters.live_allocations;
    test_counters.live_bytes -= sizeof(T);
#endif
  }
  static void visit() noexcept {
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    ++test_counters.node_visits;
#endif
  }
  static void transfer() noexcept {
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    ++test_counters.child_transfers;
#endif
  }
  static void mutate() noexcept {
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    ++test_counters.payload_mutations;
#endif
  }
  static void destroy(void* p, std::size_t height) noexcept {
    if (!p) {
      return;
    }
    if (height == 0) {
      dispose(static_cast<Leaf*>(p));
    } else {
      auto* node = static_cast<Node*>(p);
      for (std::size_t i = 0; i < node->count; ++i) {
        destroy(node->children[i], height - 1);
      }
      dispose(node);
    }
  }
  struct Owner {
    void* p = nullptr;
    std::size_t total = 0;
    std::size_t height = 0;
    Owner() noexcept = default;
    Owner(void* pointer, std::size_t n, std::size_t h) noexcept
        : p(pointer), total(n), height(h) {}
    Owner(const Owner&) = delete;
    Owner& operator=(const Owner&) = delete;
    Owner(Owner&& other) noexcept
        : p(std::exchange(other.p, nullptr)),
          total(std::exchange(other.total, 0)),
          height(std::exchange(other.height, 0)) {}
    Owner& operator=(Owner&& other) noexcept {
      if (this != &other) {
        destroy(p, height);
        p = std::exchange(other.p, nullptr);
        total = std::exchange(other.total, 0);
        height = std::exchange(other.height, 0);
      }
      return *this;
    }
    ~Owner() { destroy(p, height); }
    void* release() noexcept { return std::exchange(p, nullptr); }
  };

  // Spares use their unoccupied first slot as a free-list link. No separate
  // vector allocations, no live ownership objects copied as numeric arrays.
  // Only this operation owns the pool; dismantled nodes join it after their
  // children have transferred to Owners. Stale slots beyond count own nothing.
  struct Pool {
    Node* nodes = nullptr;
    std::array<Owner, 3> leaves;
    std::size_t leaf_count = 0;
    Pool() = default;
    Pool(const Pool&) = delete;
    Pool& operator=(const Pool&) = delete;
    ~Pool() {
      while (nodes) {
        auto* next = static_cast<Node*>(nodes->children[0]);
        dispose(nodes);
        nodes = next;
      }
    }
    void reserve(std::size_t node_count, std::size_t leaf_spares) {
      for (std::size_t i = 0; i < node_count; ++i) {
        recycle(allocate<Node>());
      }
      for (; leaf_count < leaf_spares; ++leaf_count) {
        leaves[leaf_count] = Owner(allocate<Leaf>(), 0, 0);
      }
    }
    void recycle(Node* node) noexcept {
      node->count = 0;
      if constexpr (IndexBias) {
        node->bias.value = 0;
      }
      node->children[0] = nodes;
      nodes = node;
    }
    Node* node() noexcept {
      assert(nodes);
      assert(nodes->count == 0);
      auto* result = nodes;
      nodes = static_cast<Node*>(nodes->children[0]);
      result->children[0] = nullptr;
      return result;
    }
    Owner leaf() noexcept {
      assert(leaf_count);
      return std::move(leaves[--leaf_count]);
    }
  };
  using Entries = std::array<Owner, 2 * Fanout>;
  struct Pair {
    Owner left, right;
  };

  static void add_bias(void* p, std::size_t height, std::uint64_t bias) noexcept
    requires(IndexBias)
  {
    auto& pending = height == 0 ? static_cast<Leaf*>(p)->bias.value
                                : static_cast<Node*>(p)->bias.value;
    assert(bias <= std::numeric_limits<value_type>::max() - pending);
    pending += bias;
  }
  static void normalize(Leaf& leaf) noexcept {
    if constexpr (IndexBias) {
      if (leaf.bias.value != 0) {
        leaf.block.add_bias(leaf.bias.value);
        leaf.bias.value = 0;
        mutate();
      }
    }
  }
  static void redistribute(Leaf& a, Leaf& b, std::size_t left) noexcept {
    normalize(a);
    normalize(b);
    a.block.redistribute(b.block, left);
    mutate();
  }

  static std::size_t length(const Node& node,
                            std::size_t i,
                            std::size_t total,
                            std::size_t prefix) noexcept {
    if (i + 1 == node.count) {
      return total - prefix;
    }
    if constexpr (Layout == LengthLayout::cumulative) {
      return node.measures[i] - prefix;
    } else {
      return node.measures[i];
    }
  }
  static std::size_t unpack(Owner tree,
                            Entries& entries,
                            Pool& pool,
                            std::size_t offset = 0) noexcept {
    visit();
    auto* node = static_cast<Node*>(tree.p);
    const auto count = node->count;
    std::size_t prefix = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const auto n = length(*node, i, tree.total, prefix);
      if constexpr (IndexBias) {
        add_bias(node->children[i], tree.height - 1, node->bias.value);
      }
      entries[offset + i] = Owner(node->children[i], n, tree.height - 1);
      prefix += n;
      transfer();
    }
    tree.release();
    pool.recycle(node);
    return count;
  }
  static Owner pack(Node* node,
                    Entries& entries,
                    std::size_t begin,
                    std::size_t count) noexcept {
    assert(count >= 2 && count <= Fanout);
    node->count = count;
    const auto height = entries[begin].height + 1;
    std::size_t total = 0;
    for (std::size_t i = 0; i < count; ++i) {
      auto& entry = entries[begin + i];
      assert(entry.p && entry.height + 1 == height);
      total += entry.total;
      node->children[i] = entry.release();
      if (i + 1 < count) {
        node->measures[i] =
            Layout == LengthLayout::cumulative ? total : entry.total;
      }
      transfer();
    }
    return Owner(node, total, height);
  }
  static Owner group(Entries& entries,
                     std::size_t begin,
                     std::size_t count,
                     Pool& pool) noexcept {
    if (count == 0) {
      return {};
    }
    if (count == 1) {
      return std::move(entries[begin]);
    }
    return pack(pool.node(), entries, begin, count);
  }
  static Pair repack(Entries& entries, std::size_t count, Pool& pool) noexcept {
    if (count <= Fanout) {
      return {group(entries, 0, count, pool), {}};
    }
    const auto half = count / 2;
    return {group(entries, 0, half, pool),
            group(entries, half, count - half, pool)};
  }
  // Join compatible roots, carrying at most two nodes at the taller height.
  // Underfull roots are merged/redistributed before becoming nonroot children.
  // Pool deficit: equal-height unpack returns two nodes before repack uses at
  // most two. Each of gap unequal levels returns one before using at most two;
  // join may add one root. Thus gap+1 spares cover every prefix, not just the
  // final node count. Without root growth the bound is gap. A taller root with
  // a spare child slot cannot split, tightening this to gap-1 when gap > 0.
  static Pair join_level(Owner a, Owner b, Pool& pool) noexcept {
    if (a.height == b.height) {
      if (a.height == 0) {
        return {std::move(a), std::move(b)};
      }
      Entries entries;
      const auto n = unpack(std::move(a), entries, pool);
      const auto m = unpack(std::move(b), entries, pool, n);
      return repack(entries, n + m, pool);
    }
    Entries entries;
    if (a.height > b.height) {
      const auto count = unpack(std::move(a), entries, pool);
      auto seam = join_level(std::move(entries[count - 1]), std::move(b), pool);
      entries[count - 1] = std::move(seam.left);
      const bool extra = seam.right.p != nullptr;
      entries[count] = std::move(seam.right);
      return repack(entries, count + extra, pool);
    }
    const auto count = unpack(std::move(b), entries, pool);
    auto seam = join_level(std::move(a), std::move(entries[0]), pool);
    const bool extra = seam.right.p != nullptr;
    if (extra) {
      for (std::size_t i = count; i > 1; --i) {
        entries[i] = std::move(entries[i - 1]);
      }
      entries[1] = std::move(seam.right);
    }
    entries[0] = std::move(seam.left);
    return repack(entries, count + extra, pool);
  }
  static Owner join(Owner a, Owner b, Pool& pool) noexcept {
    if (!a.p) {
      return b;
    }
    if (!b.p) {
      return a;
    }
    auto pair = join_level(std::move(a), std::move(b), pool);
    if (!pair.right.p) {
      return std::move(pair.left);
    }
    Entries entries;
    entries[0] = std::move(pair.left);
    entries[1] = std::move(pair.right);
    return group(entries, 0, 2, pool);
  }

  // One cut descent. On ascent, sibling groups are already balanced subtrees.
  // Each boundary accumulator rises monotonically: spine gaps traversed by
  // join telescope across levels (plus one boundary node per level), O(F*h),
  // not a fresh full-height join at every ancestor. Neither result exceeds the
  // input height: sibling groups have at most F-1 children, so their root has
  // room for a carry; a single sibling joined to a cut grows by at most one.
  //
  /**
   * @brief Split with at most 3h fresh internal nodes, including prefix
   * deficits.
   * @details All d<=h cut ancestors are recycled before ascent. At a level
   * with child height t, a nonempty boundary accumulator has height q<=t.
   * One sibling costs at most t-q+1 join nodes. Two or more siblings cost one
   * group node plus at most (t+1-q)-1 join nodes: the group has <=F-1 children,
   * so its taller root cannot split. Their respective output heights are at
   * least t and t+1. Thus charge each side's height increase plus one only
   * for a singleton sibling. An empty accumulator costs at most one group,
   * covered by its new height; no siblings cost nothing. Heights telescope
   * to <=h per side, and there are <=2d singleton groups. Local prefix-deficit
   * charges sum to <=2h+2d, less the d recycled ancestors: <=3h. Every prefix
   * is covered, including both groups before either join and retained roots;
   * grouping is prepaid by the same height-increase charge. It is not a bound
   * inferred from final growth. Leaf cuts separately need one leaf spare.
   */
  static Pair split(Owner tree, std::size_t p, Pool& pool) noexcept {
    if (p == 0) {
      return {{}, std::move(tree)};
    }
    if (p == tree.total) {
      return {std::move(tree), {}};
    }
    if (tree.height == 0) {
      auto right = pool.leaf();
      redistribute(*static_cast<Leaf*>(tree.p), *static_cast<Leaf*>(right.p),
                   p);
      right.total = tree.total - p;
      tree.total = p;
      return {std::move(tree), std::move(right)};
    }
    Entries entries;
    const auto count = unpack(std::move(tree), entries, pool);
    std::size_t i = 0;
    while (p > entries[i].total) {
      p -= entries[i++].total;
    }
    auto cut = split(std::move(entries[i]), p, pool);
    auto left = group(entries, 0, i, pool);
    auto right = group(entries, i + 1, count - i - 1, pool);
    return {join(std::move(left), std::move(cut.left), pool),
            join(std::move(cut.right), std::move(right), pool)};
  }

  struct Location {
    Leaf* leaf;
    std::size_t offset;
    [[no_unique_address]] Bias bias{};
  };
  /**
   * @brief Rotate within a leaf or reorder complete children at a covering
   * node.
   * @details Search is read-only until all three boundaries match child slots.
   * Only global exterior leaves can be underfull. Check those included in a
   * matched child range before committing; strict interior ranges need no leaf
   * lookups. A descent reaching a leaf uses its nonthrowing local kernel.
   * Typed pointer rotation transfers occupied-slot ownership; numeric lengths
   * rotate separately and rebuild at most F-1 measures. Parent uniform bias
   * stays in place and each child's own bias travels with its allocation.
   * Ancestor totals, heights, and all allocation identities remain unchanged.
   * @return True on nonallocating commit, false without any mutation otherwise.
   */
  bool rotate_local(std::size_t left,
                    std::size_t right,
                    std::size_t distance) noexcept {
    const bool includes_first = left == 0;
    const bool includes_last = right == size();
    void* p = root_.p;
    auto total = root_.total;
    for (auto h = root_.height; h != 0; --h) {
      visit();
      auto& node = *static_cast<Node*>(p);
      std::array<std::size_t, Fanout> lengths;
      std::size_t prefix = 0;
      std::size_t begin = Fanout, middle = Fanout, end = Fanout + 1;
      bool descend = false;
      for (std::size_t i = 0; i < node.count; ++i) {
        const auto n = length(node, i, total, prefix);
        if (left >= prefix && right <= prefix + n) {
          left -= prefix;
          right -= prefix;
          total = n;
          p = node.children[i];
          descend = true;
          break;
        }
        lengths[i] = n;
        if (prefix == left) {
          begin = i;
        }
        if (prefix == left + distance) {
          middle = i;
        }
        prefix += n;
        if (prefix == right) {
          end = i + 1;
        }
      }
      if (descend) {
        continue;
      }
      if (begin == Fanout || middle == Fanout || end == Fanout + 1) {
        return false;
      }
      // The matched slots already identify the endpoint subtrees. Follow only
      // their outer spines, without rescanning ancestors or selecting by rank.
      const auto underfull = [h](void* child, bool back) noexcept {
        for (auto level = h - 1; level != 0; --level) {
          visit();
          const auto& edge = *static_cast<Node*>(child);
          child = edge.children[back ? edge.count - 1 : 0];
        }
        return static_cast<Leaf*>(child)->block.size() < minimum_leaf_size;
      };
      if ((includes_first && underfull(node.children[begin], false)) ||
          (includes_last && underfull(node.children[end - 1], true))) {
        return false;
      }
      std::rotate(node.children.begin() + begin, node.children.begin() + middle,
                  node.children.begin() + end);
      std::rotate(lengths.begin() + begin, lengths.begin() + middle,
                  lengths.begin() + end);
      prefix = 0;
      for (std::size_t i = 0; i + 1 < node.count; ++i) {
        prefix += lengths[i];
        node.measures[i] =
            Layout == LengthLayout::cumulative ? prefix : lengths[i];
      }
      for (auto i = begin; i < end; ++i) {
        transfer();
      }
      return true;
    }
    static_cast<Leaf*>(p)->block.rotate_left(left, right, distance);
    mutate();
    return true;
  }
  template <bool Accumulate = false>
  static Location locate(const Owner& root, std::size_t index) noexcept {
    void* p = root.p;
    auto total = root.total;
    [[maybe_unused]] Bias bias{};
    for (auto h = root.height; h != 0; --h) {
      visit();
      const auto& node = *static_cast<Node*>(p);
      if constexpr (Accumulate && IndexBias) {
        bias.value += node.bias.value;
      }
      if constexpr (Layout == LengthLayout::cumulative) {
        const auto i =
            node_select({node.measures.data(), node.count - 1}, index);
        const auto prefix = i == 0 ? 0 : node.measures[i - 1];
        const auto end = i + 1 == node.count ? total : node.measures[i];
        index -= prefix;
        total = end - prefix;
        p = node.children[i];
      } else {
        std::size_t prefix = 0;
        for (std::size_t i = 0; i < node.count; ++i) {
          const auto n = length(node, i, total, prefix);
          if (index - prefix < n) {
            index -= prefix;
            total = n;
            p = node.children[i];
            break;
          }
          prefix += n;
        }
      }
    }
    if constexpr (Accumulate && IndexBias) {
      bias.value += static_cast<Leaf*>(p)->bias.value;
    }
    return {static_cast<Leaf*>(p), index, bias};
  }
  static Owner pop_edge(Owner& tree, bool back, Pool& pool) noexcept {
    // Removing an entire exterior leaf needs no spares. Inductively, reducing
    // a child from height k-1 to q returns at least k-1-q nodes. Recycling the
    // parent adds one. With >=2 siblings, grouping uses one, but its root has
    // <=F-1 children: joining costs at most k-q-1. With one sibling, no group
    // node is needed and join costs <=k-q (one less without root growth).
    // With none, all credits remain. Every prefix is funded and a height drop
    // from k to r leaves at least k-r credits, completing the induction.
    const auto n = locate(tree, back ? tree.total - 1 : 0).leaf->block.size();
    const auto position = back ? tree.total - n : n;
    auto cut = split(std::move(tree), position, pool);
    if (back) {
      tree = std::move(cut.left);
      return std::move(cut.right);
    }
    tree = std::move(cut.right);
    return std::move(cut.left);
  }
  // Payload seam repair is performed once, not at each structural ancestor.
  // At most two leaves from each side suffice: if anything remains outside,
  // an extracted neighbor was interior and hence at least half full. Compact
  // the four blocks, then balance the final pair. A lone underfull result can
  // therefore only occur when it is an exterior leaf of the whole result.
  // Four pop_edge calls never increase the node deficit. Group the <=4 seam
  // leaves (4<=F) using <=1 node, then join twice, rather than once per leaf.
  // For h>=1, the group has height <=1<=h. The joins cost <=h+1 and <=h+2
  // spares and produce height <=h+2. Thus 1+(h+1)+(h+2)=2h+4 covers every
  // prefix, including retained results. For h=0 both inputs are extracted in
  // full, so only the group can cost a node. No leaf spares are used. An
  // underfull group root is repaired by join_level before becoming nonroot.
  static Owner concatenate(Owner a, Owner b, Pool& pool) noexcept {
    if (!a.p || !b.p) {
      return join(std::move(a), std::move(b), pool);
    }
    if (locate(a, a.total - 1).leaf->block.size() >= minimum_leaf_size &&
        locate(b, 0).leaf->block.size() >= minimum_leaf_size) {
      return join(std::move(a), std::move(b), pool);
    }
    Entries seam;
    seam[1] = pop_edge(a, true, pool);
    if (a.p) {
      seam[0] = pop_edge(a, true, pool);
    }
    seam[2] = pop_edge(b, false, pool);
    if (b.p) {
      seam[3] = pop_edge(b, false, pool);
    }
    std::size_t count = 0;
    for (std::size_t i = 0; i < 4; ++i) {
      if (seam[i].p) {
        if (i != count) {
          seam[count] = std::move(seam[i]);
        }
        ++count;
      }
    }
    for (std::size_t i = 0; i + 1 < count;) {
      auto& x = seam[i];
      auto& y = seam[i + 1];
      const auto total = x.total + y.total;
      const auto first = std::min(block_capacity, total);
      redistribute(*static_cast<Leaf*>(x.p), *static_cast<Leaf*>(y.p), first);
      x.total = first;
      y.total = total - first;
      if (y.total == 0) {
        for (std::size_t j = i + 1; j + 1 < count; ++j) {
          seam[j] = std::move(seam[j + 1]);
        }
        seam[--count] = {};
      } else {
        ++i;
      }
    }
    if (count >= 2 && seam[count - 1].total < minimum_leaf_size) {
      auto& x = seam[count - 2];
      auto& y = seam[count - 1];
      const auto total = x.total + y.total;
      redistribute(*static_cast<Leaf*>(x.p), *static_cast<Leaf*>(y.p),
                   total / 2);
      x.total = total / 2;
      y.total = total - x.total;
    }
    a = join(std::move(a), group(seam, 0, count, pool), pool);
    return join(std::move(a), std::move(b), pool);
  }

 public:
  /** @brief Construct canonical empty, with no allocation. */
  SequenceTree() noexcept = default;
  /** @brief Exclusive ownership disallows copying. */
  SequenceTree(const SequenceTree&) = delete;
  /** @brief Exclusive ownership disallows copy assignment. */
  SequenceTree& operator=(const SequenceTree&) = delete;
  /** @brief Transfer ownership; source becomes canonical empty. */
  SequenceTree(SequenceTree&&) noexcept = default;
  /** @brief Transfer ownership; source becomes empty; self-move is a no-op. */
  SequenceTree& operator=(SequenceTree&&) noexcept = default;
  /** @brief Return logical element count. @return Count in [0, SIZE_MAX]. */
  std::size_t size() const noexcept { return root_.total; }
  /** @brief Test emptiness. @return Whether no elements or allocations exist.
   */
  bool empty() const noexcept { return !root_.p; }
  /** @brief Return internal levels. @return Zero for an empty tree or leaf
   * root. */
  std::size_t height() const noexcept { return root_.height; }
  /**
   * @brief Read a zero-based element by value in O(F*h).
   * @param i Index in [0, size()).
   * @return Immutable logical value.
   * @throws std::out_of_range If i >= size().
   */
  value_type operator[](std::size_t i) const {
    if (i >= size()) {
      throw std::out_of_range("SequenceTree: index");
    }
    const auto position = locate<IndexBias>(root_, i);
    const auto value = std::as_const(position.leaf->block)[position.offset];
    if constexpr (IndexBias) {
      return static_cast<value_type>(value + position.bias.value);
    } else {
      return value;
    }
  }
  /**
   * @brief Consume a single-pass range of blocks in streaming bottom-up order.
   * @details Moves from each input block; empty blocks are ignored. Packs short
   * inputs into full leaves, retaining O(F*h) pending owners, not a directory
   * or second complete payload buffer. Completed levels are built bottom-up;
   * finishing the bounded fringe groups and joins once per level. On iterator,
   * allocation or length failure, already visited input may be consumed; all
   * acquired ownership is reclaimed. This is not a transactional range API.
   * @param blocks Mutable input range whose references can move into Block.
   * @return Exclusively owned tree preserving concatenated element order.
   * @throws std::length_error If total exceeds SIZE_MAX or a block exceeds
   * capacity.
   * @throws std::bad_alloc If construction allocation fails.
   */
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<Block, std::ranges::range_rvalue_reference_t<Range>>
      static SequenceTree from_blocks(Range&& blocks) {
    std::array<Entries, max_height> levels;
    std::array<std::size_t, max_height> counts{};
    std::size_t total = 0;
    Block pending;
    auto emit = [&](Block&& block) {
      Owner carry(allocate<Leaf>(), block.size(), 0);
      static_cast<Leaf*>(carry.p)->block = std::move(block);
      std::size_t level = 0;
      while (true) {
        levels[level][counts[level]++] = std::move(carry);
        if (counts[level] != Fanout) {
          break;
        }
        auto* node = allocate<Node>();
        carry = pack(node, levels[level], 0, Fanout);
        counts[level++] = 0;
        assert(level < max_height);
      }
    };
    auto it = std::ranges::begin(blocks);
    const auto end = std::ranges::end(blocks);
    for (; it != end; ++it) {
      Block incoming(std::ranges::iter_move(it));
      const auto n = incoming.size();
      if (n > Block::capacity ||
          n > std::numeric_limits<std::size_t>::max() - total) {
        throw std::length_error("SequenceTree: construction size");
      }
      total += n;
      if (n == 0) {
        continue;
      }
      const auto combined = pending.size() + n;
      pending.redistribute(incoming, std::min(block_capacity, combined));
      if (pending.size() == Block::capacity) {
        emit(std::move(pending));
        pending = std::move(incoming);
      }
    }
    if (pending.size()) {
      emit(std::move(pending));
    }
    SequenceTree result;
    if (total <= Block::capacity) {
      result.root_ = std::move(levels[0][0]);
      return result;
    }
    // Assemble from low to high: each level's older sibling group precedes the
    // existing fringe. Height gaps telescope just as in split's boundary
    // assembly, rather than repeatedly descending from the tallest root.
    Pool pool;
    std::size_t used_height = 0;
    for (std::size_t i = 0; i < max_height; ++i) {
      if (counts[i]) {
        used_height = i;
      }
    }
    // <=h+1 groups and join roots, plus telescoping gaps <=h+1. This covers
    // every prefix while the unassembled higher-level Owners remain live.
    pool.reserve(3 * (used_height + 1), 0);
    for (std::size_t i = 0; i <= used_height; ++i) {
      auto prefix = group(levels[i], 0, counts[i], pool);
      result.root_ = join(std::move(prefix), std::move(result.root_), pool);
    }
    return result;
  }
  /**
   * @brief Keep [0,p), returning [p,size()) by exclusive ownership transfer.
   * @details Endpoints are accepted. O(F*h) structure and at most one block
   * redistribution. Allocation failure leaves this tree unchanged.
   * @param p Zero-based split boundary in [0,size()].
   * @return Consumed suffix; neither result shares storage.
   * @throws std::out_of_range If p > size().
   * @throws std::bad_alloc On preflight allocation failure.
   */
  SequenceTree split_off(std::size_t p) {
    if (p > size()) {
      throw std::out_of_range("SequenceTree: split position");
    }
    SequenceTree result;
    if (p == size()) {
      return result;
    }
    if (p == 0) {
      return std::move(*this);
    }
    Pool pool;
    pool.reserve(3 * root_.height, locate(root_, p).offset != 0);
    auto pair = split(std::move(root_), p, pool);
    root_ = std::move(pair.left);
    result.root_ = std::move(pair.right);
    return result;
  }
  /**
   * @brief Append and consume a distinct donor in O(F*h).
   * @details Self-merge and empty donor are no-ops. On success donor is
   * canonical empty. Repairs only a bounded leaf seam; off-path ownership stays
   * intact. Both trees remain unchanged on allocation or total-size failure.
   * @param donor Tree whose elements follow this tree's contents.
   * @throws std::length_error If concatenation exceeds SIZE_MAX.
   * @throws std::bad_alloc On preflight allocation failure.
   */
  void merge(SequenceTree& donor) { merge_impl<false>(donor); }

 private:
  // Only Permutation may request rebasing. All allocating preflight paths
  // complete before tagging; from that point every ownership transfer and
  // bounded normalization is noexcept. The structural pool proofs are
  // unchanged.
  void merge_rebased(SequenceTree& donor)
    requires(IndexBias)
  {
    merge_impl<true>(donor);
  }
  template <bool Rebase>
  void merge_impl(SequenceTree& donor) {
    if (this == &donor || donor.empty()) {
      return;
    }
    if (donor.size() > std::numeric_limits<std::size_t>::max() - size()) {
      throw std::length_error("SequenceTree: concatenation size");
    }
    if (empty()) {
      *this = std::move(donor);
      return;
    }
    Pool pool;
    auto tag = [&]() noexcept {
      if constexpr (Rebase) {
        add_bias(donor.root_.p, donor.height(), size());
      }
    };
    if (height() == 0 && donor.height() == 0) {
      // Two leaf roots have no spine to rebuild. A fitting total combines in
      // place; otherwise seam extraction splits only at endpoints and joining
      // the two surviving leaves requires exactly one parent allocation.
      const auto total = size() + donor.size();
      if (total <= block_capacity) {
        tag();
        redistribute(*static_cast<Leaf*>(root_.p),
                     *static_cast<Leaf*>(donor.root_.p), total);
        root_.total = total;
        donor.root_ = {};
      } else {
        pool.reserve(1, 0);
        tag();
        root_ = concatenate(std::move(root_), std::move(donor.root_), pool);
      }
      return;
    }
    const auto h = std::max(height(), donor.height());
    if (std::as_const(locate(root_, size() - 1).leaf->block).size() >=
            minimum_leaf_size &&
        std::as_const(locate(donor.root_, 0).leaf->block).size() >=
            minimum_leaf_size) {
      // Recycled spine nodes cover replacements; only gap+1 additional nodes
      // can be simultaneously needed. No payload repair or leaf spares.
      const auto gap = h - std::min(height(), donor.height());
      pool.reserve(gap + 1, 0);
      tag();
      root_ = join(std::move(root_), std::move(donor.root_), pool);
      return;
    }
    pool.reserve(2 * h + 4, 0);
    tag();
    root_ = concatenate(std::move(root_), std::move(donor.root_), pool);
  }

 public:
  /**
   * @brief Rotate [left,right) left, preserving all outside elements.
   * @details A one-leaf range uses its nonthrowing local kernel. Complete-child
   * ranges at a common covering node reorder in place when both endpoint leaves
   * are at least half full. Otherwise one preflight covers cuts and joins, so
   * no intermediate edit can escape on allocation failure. O(F*h) structure,
   * bounded leaf copying; no persistent spare-node storage.
   * @param left Inclusive start.
   * @param right Exclusive end.
   * @param distance Left distance reduced modulo a nonempty range length.
   * @throws std::out_of_range If left > right or right > size(), even at
   * distance zero.
   * @throws std::bad_alloc On preflight failure, leaving contents unchanged.
   */
  void rotate_left(std::size_t left, std::size_t right, std::size_t distance) {
    if (left > right || right > size()) {
      throw std::out_of_range("SequenceTree: rotation range");
    }
    const auto n = right - left;
    if (n == 0 || (distance %= n) == 0) {
      return;
    }
    if (rotate_local(left, right, distance)) {
      return;
    }
    Pool pool;
    if (left == 0 && right == size()) {
      /** @brief One split (3h) and one seam (2h+4), with cut heights <=h. */
      pool.reserve(5 * height() + 4, locate(root_, distance).offset != 0);
      auto cut = split(std::move(root_), distance, pool);
      root_ = concatenate(std::move(cut.right), std::move(cut.left), pool);
      return;
    }
    /**
     * @brief Reserve the sum of cut and seam temporary-prefix deficits.
     * @details Three splits cost <=9h, retaining all four results. Three seams
     * cost <=(2h+4)+(2(h+2)+4)+(2(h+4)+4), each growing height by at most two.
     * The sum is 15h+24. Only nonaligned cuts need leaf spares; earlier cuts
     * preserve later offsets. This is not a net-growth-only reservation.
     */
    const auto leaf_spares =
        std::size_t(left != 0 && locate(root_, left).offset != 0) +
        (locate(root_, left + distance).offset != 0) +
        (right != size() && locate(root_, right).offset != 0);
    pool.reserve(15 * height() + 24, leaf_spares);
    auto suffix = split(std::move(root_), right, pool);
    auto b = split(std::move(suffix.left), left + distance, pool);
    auto a = split(std::move(b.left), left, pool);
    root_ = concatenate(std::move(a.left), std::move(b.right), pool);
    root_ = concatenate(std::move(root_), std::move(a.right), pool);
    root_ = concatenate(std::move(root_), std::move(suffix.right), pool);
  }

  /** @brief Explicit traversal result for memory accounting, excluding
   * allocator overhead. */
  struct MemoryUsage {
    std::size_t blocks = 0;
    std::size_t nodes = 0;
    std::size_t block_bytes = 0;
    std::size_t node_bytes = 0;
    std::size_t total_bytes = sizeof(SequenceTree);
  };
  /**
   * @brief Enumerate allocations for diagnostics, never called by mutation.
   * @details O(number of allocations). Counts actual typed object sizes,
   * including alignment padding, but excludes allocator overhead and preflight
   * peak. Byte arithmetic saturates at SIZE_MAX if an exotic block represents
   * too many tiny elements for the diagnostic footprint to be representable.
   * @return Live block/node counts and requested bytes including this object.
   */
  MemoryUsage memory_usage() const noexcept {
    MemoryUsage result;
    auto walk = [&](auto&& self, void* p, std::size_t h) -> void {
      if (!p) {
        return;
      }
      if (h == 0) {
        ++result.blocks;
      } else {
        ++result.nodes;
        const auto& node = *static_cast<Node*>(p);
        for (std::size_t i = 0; i < node.count; ++i) {
          self(self, node.children[i], h - 1);
        }
      }
    };
    walk(walk, root_.p, height());
    result.block_bytes = saturated_multiply(result.blocks, sizeof(Leaf));
    result.node_bytes = saturated_multiply(result.nodes, sizeof(Node));
    result.total_bytes = saturated_add(
        sizeof(*this), saturated_add(result.block_bytes, result.node_bytes));
    return result;
  }
  /** @brief Explicit O(allocations) footprint query. @return Requested live
   * bytes. */
  std::size_t memory_usage_bytes() const noexcept {
    return memory_usage().total_bytes;
  }
  /** @brief Explicit O(allocations) leaf count. @return Number of live blocks.
   */
  std::size_t block_count() const noexcept { return memory_usage().blocks; }
  /**
   * @brief Explicit O(allocations) internal-node count.
   * @return Number of live internal nodes, excluding leaves and the root owner.
   */
  std::size_t node_count() const noexcept { return memory_usage().nodes; }
  /** @brief Explicit O(allocations) leaf bytes. @return Includes leaf
   * metadata/padding. */
  std::size_t block_allocation_bytes() const noexcept {
    return memory_usage().block_bytes;
  }
  /** @brief Explicit O(allocations) node bytes. @return Excludes the root
   * owner. */
  std::size_t internal_memory_bytes() const noexcept {
    return memory_usage().node_bytes;
  }
  /**
   * @brief Explicit O(allocations) payload accounting for blocks reporting it.
   * @return Payload capacity bytes including unused elements, excluding
   * metadata.
   */
  std::size_t payload_capacity_bytes() const noexcept
    requires requires(const Block& b) {
      { b.payload_capacity_bytes() } noexcept -> std::same_as<std::size_t>;
    }
  {
    return saturated_multiply(block_count(), Block{}.payload_capacity_bytes());
  }
  /**
   * @brief Explicit O(allocations) metadata accounting for payload-reporting
   * blocks.
   * @return Root owner, internal nodes, block metadata and alignment padding.
   */
  std::size_t metadata_bytes() const noexcept
    requires requires(const Block& b) {
      { b.payload_capacity_bytes() } noexcept -> std::same_as<std::size_t>;
    }
  {
    return memory_usage_bytes() - payload_capacity_bytes();
  }
  /**
   * @brief Visit blocks in logical order without materializing another buffer.
   * @details Explicit O(allocations) traversal; callback receives const Block&.
   * Mutation during traversal is unsupported; callback exceptions propagate.
   * Unavailable on tagged trees: raw stored fields omit pending index biases.
   * @param callback Callable invoked once per nonempty block in sequence order.
   */
  template <class Callback>
  void for_each_block(Callback&& callback) const
    requires(!IndexBias)
  {
    auto walk = [&](auto&& self, void* p, std::size_t h) -> void {
      if (!p) {
        return;
      }
      if (h == 0) {
        callback(std::as_const(static_cast<Leaf*>(p)->block));
      } else {
        const auto& node = *static_cast<Node*>(p);
        for (std::size_t i = 0; i < node.count; ++i) {
          self(self, node.children[i], h - 1);
        }
      }
    };
    walk(walk, root_.p, height());
  }

#ifdef PIXIE_SEQUENCE_TREE_TESTING
  /**
   * @brief Enumerate absolute child boundaries for shape-forcing rotation
   * tests.
   * @details Root-first preorder; each vector starts at the node's first
   * element and ends at its exclusive end. Read-only O(allocations) traversal,
   * absent from normal builds, with no effect on transaction counters or lazy
   * biases.
   * @return One allocating boundary vector per internal node.
   */
  std::vector<std::vector<std::size_t>> test_child_boundaries() const {
    std::vector<std::vector<std::size_t>> result;
    auto walk = [&](auto&& self, const void* p, std::size_t total,
                    std::size_t h, std::size_t start) -> void {
      if (!p || h == 0) {
        return;
      }
      const auto& node = *static_cast<const Node*>(p);
      std::vector<std::size_t> boundaries{start};
      std::size_t prefix = 0;
      for (std::size_t i = 0; i < node.count; ++i) {
        prefix += length(node, i, total, prefix);
        boundaries.push_back(start + prefix);
      }
      result.push_back(boundaries);
      for (std::size_t i = 0; i < node.count; ++i) {
        self(self, node.children[i], boundaries[i + 1] - boundaries[i], h - 1,
             boundaries[i]);
      }
    };
    walk(walk, root_.p, root_.total, root_.height, 0);
    return result;
  }

  /**
   * @brief Test-only recursive validation; never invoked by normal operations.
   * @details Checks lengths, depth, node/leaf occupancy, unique ownership and
   * allocation alignment. O(allocations) space/time, may allocate and throw.
   * @return Whether every representation invariant holds.
   */
  bool test_validate() const {
    if (!root_.p) {
      return size() == 0 && height() == 0;
    }
    std::unordered_set<const void*> seen;
    auto walk = [&](auto&& self, void* p, std::size_t total, std::size_t h,
                    bool first, bool last, bool root,
                    std::uint64_t inherited) -> bool {
      if (!p || !total || !seen.insert(p).second) {
        return false;
      }
      if (h == 0) {
        const auto& b = static_cast<Leaf*>(p)->block;
        if constexpr (IndexBias) {
          const auto pending = static_cast<Leaf*>(p)->bias.value;
          const auto maximum = std::numeric_limits<value_type>::max();
          if (inherited > maximum || pending > maximum - inherited) {
            return false;
          }
          inherited += pending;
          for (std::size_t i = 0; i < b.size(); ++i) {
            if (b[i] > maximum - inherited) {
              return false;
            }
          }
        }
        return reinterpret_cast<std::uintptr_t>(p) % alignof(Leaf) == 0 &&
               b.size() == total && total <= Block::capacity &&
               (first || last || total >= minimum_leaf_size);
      }
      const auto& node = *static_cast<Node*>(p);
      if constexpr (IndexBias) {
        const auto maximum = std::numeric_limits<value_type>::max();
        if (inherited > maximum || node.bias.value > maximum - inherited) {
          return false;
        }
        inherited += node.bias.value;
      }
      if (reinterpret_cast<std::uintptr_t>(p) % alignof(Node) != 0 ||
          node.count < (root ? 2 : Fanout / 2) || node.count > Fanout) {
        return false;
      }
      std::size_t prefix = 0;
      for (std::size_t i = 0; i < node.count; ++i) {
        if (prefix >= total) {
          return false;
        }
        const auto n = length(node, i, total, prefix);
        if (n > total - prefix ||
            !self(self, node.children[i], n, h - 1, first && i == 0,
                  last && i + 1 == node.count, false, inherited)) {
          return false;
        }
        prefix += n;
      }
      return prefix == total;
    };
    return walk(walk, root_.p, size(), height(), true, true, true, 0);
  }
  /** @brief Enumerate stable leaf identities for locality tests. @return
   * Ordered addresses. */
  std::vector<const void*> test_leaf_identities() const {
    std::vector<const void*> result;
    auto walk = [&](auto&& self, const void* p, std::size_t h) -> void {
      if (!p) {
        return;
      }
      if (h == 0) {
        result.push_back(&static_cast<const Leaf*>(p)->block);
      } else {
        const auto& node = *static_cast<const Node*>(p);
        for (std::size_t i = 0; i < node.count; ++i) {
          self(self, node.children[i], h - 1);
        }
      }
    };
    walk(walk, root_.p, height());
    return result;
  }
  /**
   * @brief Snapshot raw pending biases and their allocation identities in
   * tests.
   * @details Tagged trees only; root-first preorder includes every node and
   * leaf, including zero biases. This const traversal neither pushes tags nor
   * changes operation counters. O(allocations) time and returned vector space;
   * vector allocation may throw. Absent from production builds.
   * @return Allocation address and exact unaccumulated pending bias pairs.
   */
  std::vector<std::pair<const void*, std::uint64_t>> test_bias_snapshot() const
    requires(IndexBias)
  {
    std::vector<std::pair<const void*, std::uint64_t>> result;
    auto walk = [&](auto&& self, const void* p, std::size_t h) -> void {
      if (!p) {
        return;
      }
      if (h == 0) {
        result.emplace_back(p, static_cast<const Leaf*>(p)->bias.value);
      } else {
        const auto& node = *static_cast<const Node*>(p);
        result.emplace_back(p, node.bias.value);
        for (std::size_t i = 0; i < node.count; ++i) {
          self(self, node.children[i], h - 1);
        }
      }
    };
    walk(walk, root_.p, height());
    return result;
  }
  /**
   * @brief Attach a synthetic test-only bias without huge identity allocation.
   * @details Requires a nonempty tree and every resulting value representable.
   * Not present in production and never an unchecked public container action.
   */
  void test_add_bias(std::uint64_t bias) noexcept
    requires(IndexBias)
  {
    assert(!empty());
    add_bias(root_.p, height(), bias);
  }
  /**
   * @brief Enumerate internal-node identities for test-only locality checks.
   * @details Explicit O(allocations) traversal, never called by mutation.
   * Addresses follow root-first preorder with children in logical order.
   * The returned vector allocates; callback-free enumeration does not update
   * operation counters. Empty trees and leaf roots return an empty vector.
   * @return Internal allocation addresses, excluding leaves and the root owner.
   */
  std::vector<const void*> test_internal_node_identities() const {
    std::vector<const void*> result;
    auto walk = [&](auto&& self, const void* p, std::size_t h) -> void {
      if (!p || h == 0) {
        return;
      }
      result.push_back(p);
      const auto& node = *static_cast<const Node*>(p);
      for (std::size_t i = 0; i < node.count; ++i) {
        self(self, node.children[i], h - 1);
      }
    };
    walk(walk, root_.p, height());
    return result;
  }
#endif

 private:
  static std::size_t saturated_add(std::size_t a, std::size_t b) noexcept {
    return b > std::numeric_limits<std::size_t>::max() - a
               ? std::numeric_limits<std::size_t>::max()
               : a + b;
  }
  static std::size_t saturated_multiply(std::size_t a, std::size_t b) noexcept {
    return b != 0 && a > std::numeric_limits<std::size_t>::max() / b
               ? std::numeric_limits<std::size_t>::max()
               : a * b;
  }
#ifdef PIXIE_SEQUENCE_TREE_TESTING
  inline static thread_local std::ptrdiff_t failure_ = -1;
#endif
  Owner root_;
};

}  // namespace pixie::detail::sequence
/// @endcond
