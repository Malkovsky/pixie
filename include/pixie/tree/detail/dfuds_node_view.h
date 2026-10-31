#pragma once

#include <pixie/tree.h>

#include <cstddef>
#include <optional>
#include <utility>

namespace pixie::detail {

/**
 * @brief Non-owning navigation view with lazy per-node query caches.
 * @details The tree and source must outlive the view and remain stable.
 * Copies have independent caches. Shared access to one view requires external
 * synchronization. Queries keep the current node; only reset moves it.
 */
template <class Tree>
class DfudsNodeView {
 public:
  /** @brief Return the current valid node handle. */
  TreeNode node() const { return node_; }
  /** @brief Move to a valid node of the same tree and invalidate all caches.
   */
  void reset(TreeNode node) { *this = DfudsNodeView(*tree_, node); }
  /** @brief Test whether the current node is root. */
  bool is_root() const { return tree_->is_root(node_); }
  /** @brief Test whether the current node is a leaf. */
  bool is_leaf() const { return tree_->is_leaf(node_); }
  /** @brief Return child i, requiring i in [0, degree()). */
  TreeNode child(size_t i) const {
    auto result = first_child();
    while (i--) {
      result = tree_->next_sibling(result);
    }
    return result;
  }
  /** @brief Return the current node's child count. */
  size_t degree() const { return delimiter() - node_.pos; }
  /** @brief Return the first child; requires a non-leaf node. */
  TreeNode first_child() const {
    return TreeNode(node_.number + 1, delimiter() + 1);
  }
  /** @brief Return the next sibling; requires a non-root, non-last child. */
  TreeNode next_sibling() const {
    if (!sibling_) {
      if (is_leaf()) {
        sibling_ = TreeNode(node_.number + 1, node_.pos + 1);
      } else {
        const auto pos = end() + 1;
        // A k-node DFUDS subtree occupies 2*k-1 bits.
        sibling_ = TreeNode(node_.number + (pos - node_.pos + 1) / 2, pos);
      }
    }
    return *sibling_;
  }
  /** @brief Return the parent, or root itself. */
  TreeNode parent() const {
    if (is_root()) {
      return node_;
    }
    if (!parent_) {
      parent_ = tree_->parent_from_open(node_, open());
    }
    return *parent_;
  }
  /** @brief Test for the last sibling; requires a non-root node. */
  bool is_last_child() const {
    if (!last_) {
      last_ = tree_->starts_degree_run(open());
    }
    return *last_;
  }

  /** @brief Enumerate children with TreeBase::for_each_child semantics. */
  template <class Visitor>
  void for_each_child(Visitor&& visitor) const {
    // degree() and first_child() share the delimiter already retained by
    // this view, including when the caller queried degree before inspection.
    const auto count = degree();
    if (count == 0) {
      return;
    }
    auto current = DfudsNodeView(*tree_, first_child());
    for (size_t i = 0; i < count; ++i) {
      visitor(current.node());
      if (i + 1 < count) {
        current.reset(current.next_sibling());
      }
    }
  }

  /** @brief Bind to a valid node of an immutable tree that outlives this view.
   */
  DfudsNodeView(const Tree& tree, TreeNode node) : tree_(&tree), node_(node) {}

 private:
  size_t delimiter() const {
    if (!delimiter_) {
      delimiter_ = tree_->rmm_.select0(node_.number + 1);
    }
    return *delimiter_;
  }
  size_t end() const {
    if (!end_) {
      end_ = tree_->rmm_.fwdsearch(node_.pos, -1);
    }
    return *end_;
  }
  size_t open() const {
    if (!open_) {
      open_ = tree_->incoming_open(node_);
    }
    return *open_;
  }
  mutable std::optional<size_t> delimiter_, end_, open_;

  const Tree* tree_;
  TreeNode node_;
  mutable std::optional<TreeNode> parent_, sibling_;
  mutable std::optional<bool> last_;
};

}  // namespace pixie::detail
