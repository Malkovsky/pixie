#pragma once

#include <pixie/tree/detail/forest256_cache.h>

namespace pixie::tree::detail {

/// Coordinates belong to their originating tree. Foreign handles unsupported.
struct ForestCoverNode {
  uint32_t unit = 0;   ///< Zero-based forest number.
  uint32_t local = 0;  ///< Zero-based BFS ordinal in that forest.
  bool operator==(const ForestCoverNode&) const = default;
};

/// One decoded row by value. Reset or eviction in a callback cannot
/// invalidate a child run; the original tree must remain alive throughout.
template <class Tree>
class ForestCoverNodeView {
  friend Tree;
  using Node = ForestCoverNode;
  using U = uint32_t;
  const Tree* tree_;
  Node node_;
  forest256::NodeInfo row_;
  U portal() const {
    return tree_->portal(tree_->sizes_.at(node_.unit) + node_.local);
  }
  ForestCoverNodeView(const Tree& tree, Node n, forest256::NodeInfo row)
      : tree_(&tree), node_(n), row_(row) {}

 public:
  Node node() const { return node_; }
  size_t number() const { return tree_->number(node_); }
  bool is_root() const { return node_ == Node{}; }
  void reset(Node n) { *this = tree_->navigate_impl(n); }
  size_t degree() const {
    if (row_.degree) {
      return row_.degree;
    }
    const U p = portal();
    return p == Tree::missing ? 0
                              : tree_->roots_.at(tree_->end(p)) -
                                    tree_->roots_.at(tree_->begin(p));
  }
  bool is_leaf() const { return degree() == 0; }
  Node child(size_t k) const {
    if (k < row_.degree) {
      return {node_.unit, row_.first + U(k)};
    }
    if (!row_.degree) {
      const U p = portal();
      if (p != Tree::missing) {
        const U x = tree_->begin(p), y = tree_->end(p),
                base = tree_->roots_.at(x);
        if (k < tree_->roots_.at(y) - base) {
          const U j = base + k, u = x + 1 == y ? x : tree_->roots_.find(j);
          return {u, j - tree_->roots_.at(u)};
        }
      }
    }
    throw std::out_of_range("child index");
  }
  Node first_child() const { return child(0); }
  Node parent() const {
    if (row_.parent != forest256::kNoParent) {
      return {node_.unit, row_.parent};
    }
    return node_.unit
               ? tree_->find_by_number(tree_->key(tree_->group(node_.unit)))
               : Node{};
  }
  bool is_last_child() const {
    return row_.last &&
           (row_.parent != forest256::kNoParent || !node_.unit ||
            ((tree_->ends_[(node_.unit - 1) / 64] >> ((node_.unit - 1) % 64)) &
             1));
  }
  Node next_sibling() const {
    if (!row_.last) {
      return {node_.unit, node_.local + 1};
    }
    if (!is_last_child()) {
      return {node_.unit + 1, 0};
    }
    throw std::out_of_range("last child");
  }
  template <class F>
  void for_each_child(F&& visit) const {
    if (row_.degree) {
      const auto row = row_;
      const U unit = node_.unit;
      for (U i = 0; i < row.degree; ++i) {
        visit(Node{unit, row.first + i});
      }
      return;
    }
    const auto* tree = tree_;
    const U p = portal();
    if (p != Tree::missing) {
      const U y = tree->end(p);
      for (U u = tree->begin(p); u < y; ++u) {
        for (U i = 0, z = 1u + tree->roots_.delta[u]; i < z; ++i) {
          visit(Node{u, i});
        }
      }
    }
  }
};
}  // namespace pixie::tree::detail
