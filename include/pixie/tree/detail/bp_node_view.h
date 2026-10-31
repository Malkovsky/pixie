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
class BpNodeView {
 public:
  /** @brief Return the current valid node handle. */
  TreeNode node() const { return node_; }
  /** @brief Move to a valid node of the same tree and invalidate all caches.
   */
  void reset(TreeNode node) { *this = BpNodeView(*tree_, node); }
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
  /** @brief Return the child count; the first call takes O(degree)
   * navigation. */
  size_t degree() const {
    if (!degree_) {
      degree_ = tree_->degree(node_);
    }
    return *degree_;
  }
  /** @brief Return the first child; requires a non-leaf node. */
  TreeNode first_child() const { return tree_->first_child(node_); }
  /** @brief Return the next sibling; requires a non-root, non-last child. */
  TreeNode next_sibling() const {
    if (!sibling_) {
      if (is_leaf()) {
        sibling_ = TreeNode(node_.number + 1, node_.pos + 2);
      } else {
        const auto pos = close() + 1;
        // A k-node BP subtree occupies exactly 2*k bits.
        sibling_ = TreeNode(node_.number + (pos - node_.pos) / 2, pos);
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
      parent_ = tree_->parent(node_);
    }
    return *parent_;
  }
  /** @brief Test for the last sibling; requires a non-root node. */
  bool is_last_child() const {
    if (!last_) {
      const auto end = close();
      last_ = tree_->rmm_.bit(end + 1) == 0;
    }
    return *last_;
  }

  /** @brief Enumerate children with TreeBase::for_each_child semantics. */
  template <class Visitor>
  void for_each_child(Visitor&& visitor) const {
    tree_->for_each_child(node_, std::forward<Visitor>(visitor));
  }

  /** @brief Bind to a valid node of an immutable tree that outlives this view.
   */
  BpNodeView(const Tree& tree, TreeNode node) : tree_(&tree), node_(node) {}

 private:
  size_t close() const {
    if (!close_) {
      close_ = is_leaf() ? node_.pos + 1 : tree_->rmm_.close(node_.pos);
    }
    return *close_;
  }
  mutable std::optional<size_t> close_, degree_;

  const Tree* tree_;
  TreeNode node_;
  mutable std::optional<TreeNode> parent_, sibling_;
  mutable std::optional<bool> last_;
};

}  // namespace pixie::detail
