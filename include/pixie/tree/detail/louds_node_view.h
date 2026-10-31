#pragma once

#include <pixie/tree.h>

#include <cstddef>
#include <optional>
#include <utility>

namespace pixie::detail {

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

  /**
   * @brief Move to a valid node of the same tree, invalidating the cache.
   * @param node Valid handle from this view's tree.
   */
  void reset(TreeNode node) {
    node_ = node;
    incoming_ = 0;
    next_.reset();
    first_.reset();
    parent_.reset();
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
    if (!first_) {
      first_ = tree_->first_child(node_);
    }
    return *first_;
  }
  /** @brief Return the next sibling; requires a non-root, non-last child. */
  TreeNode next_sibling() const {
    return TreeNode(node_.number + 1, next_delimiter());
  }

  /** @brief Test for the last sibling; requires a non-root current node. */
  bool is_last_child() const { return tree_->bv[incoming_edge() + 1]; }

  /** @brief Return the parent, or root itself; keep the current node. */
  TreeNode parent() const {
    if (is_root()) {
      return node_;
    }
    if (!parent_) {
      parent_ = tree_->parent_from_edge(node_, incoming_edge());
    }
    return *parent_;
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
  size_t incoming_edge() const {
    // Bit zero is the root delimiter, so no valid incoming edge is zero.
    if (incoming_ == 0) {
      incoming_ = tree_->bv.select0(node_.number);
    }
    return incoming_;
  }

  size_t next_delimiter() const {
    if (!next_) {
      next_ = tree_->next_delimiter(node_.pos, node_.number + 2);
    }
    return *next_;
  }

  mutable std::optional<size_t> next_;
  mutable std::optional<TreeNode> first_, parent_;
  const Tree* tree_;
  TreeNode node_;
  mutable size_t incoming_ = 0;
};

}  // namespace pixie::detail
