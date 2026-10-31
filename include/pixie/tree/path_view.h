#pragma once

#include <pixie/tree.h>

#include <array>
#include <cstddef>

namespace pixie {

/**
 * @brief Navigation cursor caching the first Capacity ancestors from the root.
 * @tparam Tree Tree exposing root() and a NodeView returned by
 * navigate().
 * @tparam Capacity Maximum cached ancestor count; zero disables parent caching.
 * @details The tree must be nonempty and outlive the cursor. The tree and its
 * backing source must remain stable. The current node is held separately from
 * the ancestor cache. Beyond Capacity levels, ascent uses ordinary navigation;
 * returning to cached depths resumes cached ascent. Storage is bounded and
 * inline, with no allocations. Copies have independent traversal state.
 * Concurrent access to one cursor requires external synchronization.
 */
template <class Tree, std::size_t Capacity>
class TreePathView {
 public:
  /** @brief Start at the root of a valid nonempty tree. */
  explicit TreePathView(const Tree& tree) : view_(tree.navigate(tree.root())) {}

  /** @brief Return the current node handle by value. */
  TreeNode node() const { return view_.node(); }
  /** @brief Return the edge distance from the tree root. */
  std::size_t depth() const { return depth_; }
  /** @brief Test whether the current node is root. */
  bool is_root() const { return depth_ == 0; }
  /** @brief Test whether the current node is a leaf. */
  bool is_leaf() const { return view_.is_leaf(); }
  /** @brief Return the current node's child count. */
  std::size_t degree() const { return view_.degree(); }
  /** @brief Test for the last sibling; requires a non-root node. */
  bool is_last_child() const { return view_.is_last_child(); }

  /** @brief Descend to the first child; requires a non-leaf node. */
  void move_to_first_child() { descend(view_.first_child()); }
  /** @brief Descend to child i, requiring i in [0, degree()). */
  void move_to_child(std::size_t i) { descend(view_.child(i)); }
  /** @brief Move to the next sibling; requires a non-root, non-last child. */
  void move_to_next_sibling() { view_.reset(view_.next_sibling()); }
  /** @brief Ascend to the parent, or remain at root when already there. */
  void move_to_parent() {
    if (depth_ == 0) {
      return;
    }
    --depth_;
    if (depth_ < Capacity) {
      const auto& parent = ancestors_[depth_];
      view_.reset(TreeNode(parent.number, parent.pos));
    } else {
      view_.reset(view_.parent());
    }
  }

 private:
  struct Ancestor {
    std::size_t number = 0;
    std::size_t pos = 0;
  };

  void descend(TreeNode child) {
    if (depth_ < Capacity) {
      const auto node = view_.node();
      ancestors_[depth_] = {node.number, node.pos};
    }
    ++depth_;
    view_.reset(child);
  }

  typename Tree::NodeView view_;
  std::size_t depth_ = 0;
  std::array<Ancestor, Capacity> ancestors_{};
};

}  // namespace pixie
