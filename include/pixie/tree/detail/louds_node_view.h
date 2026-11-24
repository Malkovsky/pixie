#pragma once

#include <pixie/tree.h>

#include <cstddef>
#include <limits>
#include <utility>

namespace pixie::tree::detail {

/**
 * @brief Non-owning navigation view with lazy per-node navigation caches.
 * @details The tree and its source must outlive the view and remain stable.
 * Copies have independent caches. Operations on the same view require
 * external synchronization; independent views do not modify the tree.
 */
template <class Tree>
class LoudsNodeView {
 public:
  /** @brief Return the current valid node handle by value. */
  TreeNode node() const { return node_; }
  /** @brief Return the current node's zero-based level-order number. */
  size_t number() const { return node_.number; }

  /**
   * @brief Move to a valid node of the same tree, invalidating the cache.
   * @param node Valid handle from this view's tree.
   */
  void reset(TreeNode node) {
    node_ = node;
    incoming_ = kUncached;
    next_ = 0;
    first_.number = kUncached;
    parent_.number = kUncached;
  }

  /** @brief Test whether the current node is the root. */
  bool is_root() const { return tree_->is_root(node_); }
  /** @brief Test whether the current node is a leaf. */
  bool is_leaf() const { return tree_->is_leaf(node_); }
  /** @brief Return the current node's child count. */
  size_t degree() const {
    return is_leaf() ? 0 : next_delimiter() - node_.pos - 1;
  }
  /** @brief Return child i, requiring i in [0, degree()). */
  TreeNode child(size_t i) const {
    return i == 0 ? first_child() : tree_->child(node_, i);
  }
  /** @brief Return the first child; requires a non-leaf current node. */
  TreeNode first_child() const {
    if (first_.number == kUncached) {
      first_ = tree_->first_child(node_);
    }
    return first_;
  }
  /** @brief Return the next sibling; requires a non-root, non-last child. */
  TreeNode next_sibling() const {
    return TreeNode(node_.number + 1, next_delimiter());
  }

  /** @brief Test for the last sibling; requires a non-root current node. */
  bool is_last_child() const { return !tree_->bv[incoming_edge() + 1]; }

  /** @brief Return the parent, or root itself; keep the current node. */
  TreeNode parent() const {
    if (is_root()) {
      return node_;
    }
    if (parent_.number == kUncached) {
      parent_ = tree_->parent_from_edge(node_, incoming_edge());
    }
    return parent_;
  }

  /** @brief Enumerate children with TreeBase::for_each_child semantics. */
  template <class Visitor>
  void for_each_child(Visitor&& visitor) const {
    tree_->for_each_child(node_, std::forward<Visitor>(visitor));
  }

  /** @brief Bind to a valid node of an immutable tree that outlives this view.
   */
  LoudsNodeView(const Tree& tree, TreeNode node) : tree_(&tree), node_(node) {}

 private:
  static constexpr size_t kUncached = std::numeric_limits<size_t>::max();

  size_t incoming_edge() const {
    // Bit zero can be a valid incoming edge.
    if (incoming_ == kUncached) {
      incoming_ = tree_->bv.select(node_.number);
    }
    return incoming_;
  }

  size_t next_delimiter() const {
    if (next_ == 0) {
      next_ = tree_->next_delimiter(node_.pos, node_.number + 1);
    }
    return next_;
  }

  mutable size_t next_ = 0;  // A next record always starts after bit zero.
  mutable TreeNode first_{kUncached, 0}, parent_{kUncached, 0};
  const Tree* tree_;
  TreeNode node_;
  mutable size_t incoming_ = kUncached;
};

}  // namespace pixie::tree::detail
