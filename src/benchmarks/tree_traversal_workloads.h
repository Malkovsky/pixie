#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/detail/unpack_louds.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/path_view.h>
#include <pixie/tree/root_navigation_index.h>
#include <pixie/utils.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

#include "tree_dfs/dfs.h"
#include "tree_dfs/inspect.h"

namespace pixie::benchmarks {

// Dispatch to the concrete implementations' internal unpacking kernels.
// Source words belong to the fixture; they
// must match the tree and need no padding past the final logical word.
template <class Backend, class Storage, class Visitor>
void UnpackedDepthFirst(const BPTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words,
                        Visitor& visitor) {
  detail::tree_unpack::bp_byte(tree, words, 2 * tree.size(), visitor);
}
template <class Backend, class Storage, class Visitor>
void UnpackedDepthFirst(const DFUDSTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words,
                        Visitor& visitor) {
  detail::tree_unpack::dfuds_excess_stream(tree, words, 2 * tree.size() - 1,
                                           visitor);
}
template <class Backend, class Storage, class Visitor>
void UnpackedDepthFirst(const LoudsTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words,
                        Visitor& visitor) {
  detail::tree_unpack::louds(tree, words, 2 * tree.size() - 1, visitor);
}

// Emit entry events and one return event after each internal node's last child.
// Visitors are borrowed; all traversal state is local, including nested calls.
template <bool Unpack = false, class Tree, class Visitor>
void DepthFirst(const Tree& tree,
                std::span<const uint64_t> words,
                Visitor&& visitor) {
  if (tree.empty()) {
    return;
  }
  if constexpr (Unpack) {
    UnpackedDepthFirst(tree, words, visitor);
  } else {
    tree_dfs::cached(tree, visitor);
  }
}

// Use the reference's stored degree; count succinct children during
// enumeration.
template <class Tree, class Node, class Visitor>
std::size_t InspectChildren(const Tree& tree, Node node, Visitor&& visit) {
  if constexpr (std::is_same_v<Tree, AdjListTree>) {
    const auto count = tree.degree(node);
    tree.for_each_child(node, visit);
    return count;
  } else {
    std::size_t count = 0;
    auto inspect = [&](auto child) {
      visit(child);
      ++count;
    };
    if constexpr (requires { tree.navigate(node); }) {
      tree.navigate(node).for_each_child(inspect);
    } else {
      tree.for_each_child(node, inspect);
    }
    return count;
  }
}

// Give the reference the same cursor operations without modifying its API.
template <class Tree>
class ReferencePath {
 public:
  explicit ReferencePath(const Tree& tree) : tree_(tree), node_(tree.root()) {}
  auto node() const { return node_; }
  bool is_leaf() const { return tree_.is_leaf(node_); }
  bool is_root() const { return tree_.is_root(node_); }
  bool is_last_child() const { return tree_.is_last_child(node_); }
  void move_to_first_child() { node_ = tree_.first_child(node_); }
  void move_to_next_sibling() { node_ = tree_.next_sibling(node_); }
  void move_to_parent() { node_ = tree_.parent(node_); }

 private:
  const Tree& tree_;
  decltype(std::declval<const Tree&>().root()) node_;
};

template <class Tree>
auto MakePath(const Tree& tree) {
  if constexpr (requires { tree.navigate(tree.root()); }) {
    return TreePathView<Tree>(tree);
  } else {
    return ReferencePath(tree);
  }
}

struct TraversalWork {
  std::size_t nodes = 0;
  std::size_t children = 0;
};

// Inspect children exactly once on node entry, never on a return event.
// Requires a nonempty tree. One complete traversal is one benchmark iteration.
template <class Tree, class Visit, class Inspect>
TraversalWork DepthFirstInspect(const Tree& tree,
                                Visit&& visit,
                                Inspect&& inspect) {
  auto path = MakePath(tree);
  TraversalWork work;
  bool descending = true;
  while (true) {
    visit(path.node());
    if (descending) {
      ++work.nodes;
      work.children += InspectChildren(tree, path.node(), inspect);
      if (!path.is_leaf()) {
        path.move_to_first_child();
        continue;
      }
    }
    descending = false;
    if (path.is_root()) {
      return work;
    }
    if (path.is_last_child()) {
      path.move_to_parent();
    } else {
      path.move_to_next_sibling();
      descending = true;
    }
  }
}

// One pooled random value per non-leaf, reduced to a child index. This samples
// child choices, not leaves uniformly. Both variants consume the same stream.
// Inspection retains the selected handle instead of fetching that child twice.
template <bool Inspect, class Tree, class Queries, class Visit, class Observe>
TraversalWork RandomRootToLeaf(const Tree& tree,
                               Queries& queries,
                               Visit&& visit,
                               Observe&& observe) {
  auto node = tree.root();
  TraversalWork work;
  while (true) {
    visit(node);
    ++work.nodes;
    if (tree.is_leaf(node)) {
      return work;
    }
    const auto choose = [&](const auto& view) {
      const auto index = queries.next() % view.degree();
      if constexpr (Inspect) {
        auto chosen = node;
        std::size_t i = 0;
        auto inspect = [&](auto child) {
          observe(child);
          if (i++ == index) {
            chosen = child;
          }
        };
        if constexpr (requires { view.for_each_child(inspect); }) {
          view.for_each_child(inspect);
        } else {
          tree.for_each_child(node, inspect);
        }
        work.children += i;
        return chosen;
      } else {
        return view.child(index);
      }
    };
    if constexpr (requires { tree.navigate(node); }) {
      node = choose(tree.navigate(node));
    } else {
      // The reference has no navigation view.
      struct View {
        const Tree& tree;
        decltype(node) current;
        auto degree() const { return tree.degree(current); }
        auto child(std::size_t i) const { return tree.child(current, i); }
      };
      node = choose(View{tree, node});
    }
  }
}
// The auxiliary index is built before timing; each walk gets a fresh cursor.
// Identical random choices and observations to the unindexed workload.
template <bool Inspect, class Index, class Queries, class Visit, class Observe>
TraversalWork IndexedRandomRootToLeaf(const Index& index,
                                      Queries& queries,
                                      Visit&& visit,
                                      Observe&& observe) {
  auto cursor = index.cursor();
  TraversalWork work;
  while (true) {
    visit(cursor.node());
    ++work.nodes;
    const size_t degree = cursor.degree();
    if (!degree) {
      return work;
    }
    const size_t choice = queries.next() % degree;
    if constexpr (Inspect) {
      cursor.inspect_children_and_move_to_child(choice, observe);
      work.children += degree;
    } else {
      cursor.move_to_child(choice);
    }
  }
}

// Per-invocation caches: never retain decoded data between timed walks.
template <bool BP, class Tree>
class ParenthesisUnpackedAccess {
 public:
  ParenthesisUnpackedAccess(const Tree& tree, std::span<const uint64_t> words)
      : tree_(tree), words_(words), view_(tree.navigate(tree.root())) {}
  void reset(TreeNode node) {
    view_.reset(node);
    leaf_ = view_.is_leaf();
    cached_ =
        cache_.contains(node) ||
        (!leaf_ && cache_.template load<BP>(
                       tree_, words_, 2 * tree_.size() - (BP ? 0 : 1), node));
  }
  size_t degree() const {
    return leaf_ ? 0 : cached_ ? cache_.degree(view_.node()) : view_.degree();
  }
  TreeNode child(size_t i) const {
    return cached_ ? cache_.child(view_.node(), i) : view_.child(i);
  }
  template <class Visitor>
  size_t children(Visitor&& visit) const {
    if (leaf_) {
      return 0;
    }
    if (cached_) {
      return cache_.children(view_.node(), visit);
    }
    return InspectChildren(tree_, view_.node(), visit);
  }

 private:
  const Tree& tree_;
  std::span<const uint64_t> words_;
  typename Tree::NodeView view_;
  detail::tree_unpack::SubtreeCache cache_;
  bool cached_ = false;
  bool leaf_ = false;
};
template <class Backend, class Storage>
auto MakeUnpackedAccess(const BPTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words) {
  return ParenthesisUnpackedAccess<true, BPTree<Backend, Storage>>(tree, words);
}
template <class Backend, class Storage>
auto MakeUnpackedAccess(const DFUDSTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words) {
  return ParenthesisUnpackedAccess<false, DFUDSTree<Backend, Storage>>(tree,
                                                                       words);
}
template <class Tree>
class LoudsUnpackedAccess {
 public:
  LoudsUnpackedAccess(const Tree& tree, std::span<const uint64_t> words)
      : tree_(tree), words_(words), view_(tree.navigate(tree.root())) {}
  void reset(TreeNode node) { view_.reset(node); }
  size_t degree() const { return view_.degree(); }
  TreeNode child(size_t i) const {
    // Bound extra child decoding to a group of 64, even for a huge star.
    TreeNode selected = view_.node();
    const size_t first = i / 64 * 64;
    detail::tree_unpack::louds_children(
        tree_, words_, 2 * tree_.size() - 1, view_.node(), first, i - first + 1,
        [&](TreeNode node) { selected = node; });
    return selected;
  }
  template <class Visitor>
  size_t children(Visitor&& visit) const {
    return InspectChildren(tree_, view_.node(), visit);
  }

 private:
  const Tree& tree_;
  std::span<const uint64_t> words_;
  typename Tree::NodeView view_;
};
template <class Backend, class Storage>
auto MakeUnpackedAccess(const LoudsTree<Backend, Storage>& tree,
                        std::span<const uint64_t> words) {
  return LoudsUnpackedAccess<LoudsTree<Backend, Storage>>(tree, words);
}

template <class Backend, class Storage, class Visit, class Inspect>
TraversalWork UnpackedDepthFirstInspect(const LoudsTree<Backend, Storage>& tree,
                                        std::span<const uint64_t> words,
                                        Visit&& visit,
                                        Inspect&& inspect) {
  using Visitor = std::remove_reference_t<Visit>;
  using Inspector = std::remove_reference_t<Inspect>;
  detail::tree_unpack::LoudsTraversal<LoudsTree<Backend, Storage>, Visitor,
                                      Inspector>
      traversal{tree, words, 2 * tree.size() - 1, visit, {}, &inspect};
  traversal.walk(tree.root(), 0);
  return {traversal.node_count, traversal.child_count};
}

template <bool BP, class Tree, class Visit, class Inspect>
TraversalWork ParenthesisDepthFirstInspect(const Tree& tree,
                                           std::span<const uint64_t> words,
                                           Visit& visit,
                                           Inspect& inspect) {
  InspectTraversal<BP, Tree, Visit, Inspect> traversal(tree, words, visit,
                                                       inspect);
  traversal.walk(tree.root());
  return {traversal.node_count, traversal.child_count};
}
template <class Backend, class Storage, class Visit, class Inspect>
TraversalWork UnpackedDepthFirstInspect(const BPTree<Backend, Storage>& tree,
                                        std::span<const uint64_t> words,
                                        Visit&& visit,
                                        Inspect&& inspect) {
  return ParenthesisDepthFirstInspect<true>(tree, words, visit, inspect);
}
template <class Backend, class Storage, class Visit, class Inspect>
TraversalWork UnpackedDepthFirstInspect(const DFUDSTree<Backend, Storage>& tree,
                                        std::span<const uint64_t> words,
                                        Visit&& visit,
                                        Inspect&& inspect) {
  return ParenthesisDepthFirstInspect<false>(tree, words, visit, inspect);
}

template <bool Inspect, class Tree, class Queries, class Visit, class Observe>
TraversalWork UnpackedRandomRootToLeaf(const Tree& tree,
                                       std::span<const uint64_t> words,
                                       Queries& queries,
                                       Visit&& visit,
                                       Observe&& observe) {
  auto access = MakeUnpackedAccess(tree, words);
  auto node = tree.root();
  TraversalWork work;
  while (true) {
    visit(node);
    ++work.nodes;
    access.reset(node);
    const size_t degree = access.degree();
    if (!degree) {
      return work;
    }
    const size_t index = queries.next() % degree;
    if constexpr (Inspect) {
      size_t i = 0;
      work.children += access.children([&](TreeNode child) {
        observe(child);
        if (i++ == index) {
          node = child;
        }
      });
    } else {
      node = access.child(index);
    }
  }
}
}  // namespace pixie::benchmarks
