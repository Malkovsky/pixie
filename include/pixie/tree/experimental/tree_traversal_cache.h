#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/experimental/tree_dfuds_micro.h>
#include <pixie/tree/experimental/tree_louds_iterative.h>
#include <pixie/tree/experimental/tree_micro_compact.h>
#include <pixie/tree/experimental/tree_navigation_index.h>
#include <pixie/tree/louds.h>

namespace pixie::experimental {
namespace traversal_detail {
template <class Tree>
struct Encoding;
template <class B, class C>
struct Encoding<BPTree<B, C>> {
  static constexpr int value = 0;
};
template <class B, class C>
struct Encoding<DFUDSTree<B, C>> {
  static constexpr int value = 1;
};
template <class B, class C>
struct Encoding<LoudsTree<B, C>> {
  static constexpr int value = 2;
};
struct SizingVisitor {
  void operator()(TreeNode) const {}
};
}  // namespace traversal_detail

/**
 * @brief Runtime-budgeted root index and decoding configuration.
 * @details Borrows a nonempty immutable tree and its matching encoded words;
 * both must outlive this object and its cursors. Requires local navigation so
 * another shared cache cannot consume the same allowance. The tree constructor
 * budget covers this context, the root index and peak explicit traversal
 * buffers. Up to half the available bytes are assigned to decoding, capped at
 * 192 KiB for BP/DFUDS and 80 KiB for LOUDS. Plain streaming needs only 4 KiB
 * and remains available when inspected decoding cannot fit. Remaining bytes
 * hold the root index. The constructor selects compiled chunk/level sizes;
 * BP/DFUDS cached recursion depth is additionally computed at runtime.
 *
 * Zero/tiny budgets use direct navigation. Partial child lists and deeper
 * subtrees retain exact navigation fallbacks. Constexpr tables, source/index
 * storage, allocator metadata, compiler call-frame overhead, user callback
 * storage and construction's temporary priority queue are excluded. The limit
 * covers one active traversal; nested/concurrent walks require independent
 * budgets. Decoded data does not persist between walks.
 */
template <class Tree>
class TreeTraversalCache {
  static_assert(
      std::is_same_v<typename Tree::Cache,
                     typename LocalNavigationCache::template Cache<Tree>>);
  static constexpr int encoding = traversal_detail::Encoding<Tree>::value;
  static constexpr size_t kStreamBytes = 4096;
  struct Plan {
    size_t scratch = 0;
    size_t bits = 0;
    size_t depth = 0;
  };
  template <size_t Bits, class V, class I>
  using ParenthesisWalk = std::conditional_t<
      encoding == 0,
      compact_bp::BpMicroTraversal<Bits, Tree, V, I, 8, true, 1>,
      dfuds_detail::DfudsMicroTraversal<Bits, Tree, V, I, 10>>;
  using SizingVisitor = traversal_detail::SizingVisitor;
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Plan plan_;
  TreeNavigationIndex<Tree> index_;

  static constexpr size_t overhead() {
    return sizeof(TreeTraversalCache) - sizeof(TreeNavigationIndex<Tree>);
  }
  template <size_t Bits>
  static Plan parentheses_plan(size_t allowance) {
    using Walk = ParenthesisWalk<Bits, SizingVisitor, SizingVisitor>;
    if (Walk::scratch_bytes(8) > allowance) {
      return {};
    }
    const size_t depth =
        std::min(size_t{128},
                 (allowance - Walk::scratch_bytes(0)) / (64 * sizeof(size_t)));
    return {std::max(kStreamBytes, Walk::scratch_bytes(depth)), Bits, depth};
  }
  template <size_t Depth, size_t Children>
  static Plan louds_plan(size_t allowance) {
    using Walk = LoudsIterativeTraversal<Tree, SizingVisitor, SizingVisitor,
                                         true, Depth, Children>;
    constexpr size_t bytes = std::max(kStreamBytes, Walk::scratch_bytes());
    return bytes <= allowance ? Plan{bytes, 0, Depth} : Plan{};
  }
  static Plan choose_plan(size_t budget) {
    if (budget < sizeof(TreeTraversalCache) + kStreamBytes) {
      return {};
    }
    const size_t available = budget - sizeof(TreeTraversalCache);
    const size_t allowance = std::min(size_t{encoding == 2 ? 80 : 192} * 1024,
                                      std::max(kStreamBytes, available / 2));
    if constexpr (encoding == 2) {
      for (auto p :
           {louds_plan<128, 32>(allowance), louds_plan<64, 16>(allowance),
            louds_plan<32, 16>(allowance), louds_plan<16, 8>(allowance),
            louds_plan<8, 8>(allowance)}) {
        if (p.depth) {
          return p;
        }
      }
    } else {
      for (auto p :
           {parentheses_plan<8192>(allowance),
            parentheses_plan<2048>(allowance), parentheses_plan<512>(allowance),
            parentheses_plan<128>(allowance)}) {
        if (p.depth) {
          return p;
        }
      }
    }
    return {kStreamBytes, 0, 0};
  }

 public:
  /** @brief Select buffer sizes within tree.cache_limit_bytes(); words match
   * tree. */
  TreeTraversalCache(const Tree& tree, std::span<const uint64_t> words)
      : tree_(tree),
        words_(words),
        plan_(choose_plan(tree.cache_limit_bytes())),
        index_(tree,
               tree.cache_limit_bytes() >= sizeof(*this) + plan_.scratch
                   ? tree.cache_limit_bytes() - overhead() - plan_.scratch
                   : 0) {}
  /** @brief Peak explicit traversal-buffer allowance, including decoder loads.
   */
  size_t scratch_bytes() const { return plan_.scratch; }
  /** @brief Whether plain DFS uses its streaming decoder. */
  bool streaming_enabled() const { return plan_.scratch >= kStreamBytes; }
  /** @brief Whether inspected DFS can use decoded subtrees/level cursors. */
  bool optimized_inspection_enabled() const { return plan_.depth != 0; }
  /** @brief BP/DFUDS chunk size in bits; zero for LOUDS or direct inspection.
   */
  size_t decoded_chunk_bits() const { return plan_.bits; }
  /** @brief Cached traversal depth; BP/DFUDS keep decoding below this depth.
   */
  size_t cached_depth() const { return plan_.depth; }
  /** @brief Context and root-index bytes plus reserved traversal buffers. */
  size_t memory_usage_bytes() const {
    const size_t indexed = index_.memory_usage_bytes();
    return indexed || plan_.scratch
               ? overhead() + std::max(sizeof(index_), indexed) + plan_.scratch
               : 0;
  }
  /** @brief Independent root cursor borrowing this index and its tree. */
  auto cursor() const { return index_.cursor(); }

  /** @brief Emit DFS entry and final-child return events to a borrowed visitor.
   */
  template <class Visitor>
  void traverse(Visitor& visit) const {
    if (streaming_enabled()) {
      traverse_decoded(visit);
    } else {
      pixie::detail::tree_unpack::subtree(tree_, tree_.root(), visit);
    }
  }
  /** @brief DFS events plus immediate children in order before descent.
   * @details Callbacks are borrowed; exceptions propagate. Returns
   * (entered nodes, inspected children). Return events do not inspect children.
   */
  template <class Visitor, class Inspector>
  std::pair<size_t, size_t> inspect(Visitor& visit, Inspector& child) const {
    if constexpr (encoding == 2) {
      switch (plan_.depth) {
        case 128:
          return inspect_louds<128, 32>(visit, child);
        case 64:
          return inspect_louds<64, 16>(visit, child);
        case 32:
          return inspect_louds<32, 16>(visit, child);
        case 16:
          return inspect_louds<16, 8>(visit, child);
        case 8:
          return inspect_louds<8, 8>(visit, child);
      }
    } else {
      switch (plan_.bits) {
        case 8192:
          return inspect_parentheses<8192>(visit, child);
        case 2048:
          return inspect_parentheses<2048>(visit, child);
        case 512:
          return inspect_parentheses<512>(visit, child);
        case 128:
          return inspect_parentheses<128>(visit, child);
      }
    }
    size_t nodes = 0, children = 0;
    auto events = [&](TreeNode n, bool entering) {
      pixie::detail::tree_unpack::emit(visit, n, entering);
      if (entering) {
        ++nodes;
        tree_.for_each_child(n, [&](TreeNode c) {
          ++children;
          child(c);
        });
      }
    };
    pixie::detail::tree_unpack::subtree(tree_, tree_.root(), events);
    return {nodes, children};
  }

 private:
  // Separate frames prevent large-profile buffers appearing on smaller paths.
  template <class Visitor>
  [[gnu::noinline]] void traverse_decoded(Visitor& visit) const {
    if constexpr (encoding == 0) {
      pixie::detail::tree_unpack::bp_byte(tree_, words_, 2 * tree_.size(),
                                          visit);
    } else if constexpr (encoding == 1) {
      pixie::detail::tree_unpack::dfuds_excess_stream(
          tree_, words_, 2 * tree_.size() - 1, visit);
    } else {
      pixie::detail::tree_unpack::louds(tree_, words_, 2 * tree_.size() - 1,
                                        visit);
    }
  }
  template <size_t Bits, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_parentheses(
      Visitor& visit,
      Inspector& child) const {
    ParenthesisWalk<Bits, Visitor, Inspector> walk(tree_, words_, visit, child,
                                                   plan_.depth);
    walk.walk(tree_.root());
    return {walk.nodes, walk.children};
  }
  template <size_t Depth, size_t Children, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_louds(
      Visitor& visit,
      Inspector& child) const {
    LoudsIterativeTraversal<Tree, Visitor, Inspector, true, Depth, Children>
        walk(tree_, words_, visit, child);
    walk.walk(tree_.root(), 0);
    return {walk.node_count, walk.child_count};
  }
};
}  // namespace pixie::experimental
