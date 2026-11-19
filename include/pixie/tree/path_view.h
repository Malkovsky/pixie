#pragma once

#include <pixie/tree.h>

#include <array>
#include <cstddef>
#include <utility>

namespace pixie::tree {

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
  /** @brief Handle returned by this tree's navigation interface. */
  using Node = decltype(std::declval<const Tree&>().root());

  /** @brief Start at the root of a valid nonempty tree. */
  explicit TreePathView(const Tree& tree)
      : view_(tree.navigate(tree.root())),
        ancestors_(filled(view_.node(), std::make_index_sequence<Capacity>{})) {
  }

  /** @brief Reject a path cursor borrowing a temporary tree. */
  TreePathView(const Tree&&) = delete;

  /** @brief Return the current node handle by value. */
  Node node() const { return view_.node(); }
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
      view_.reset(ancestors_[depth_]);
    } else {
      view_.reset(view_.parent());
    }
  }

 private:
  // Handles need only be copyable, not default constructible.
  template <std::size_t... I>
  static std::array<Node, Capacity> filled(const Node& node,
                                           std::index_sequence<I...>) {
    return {((void)I, node)...};
  }

  void descend(Node child) {
    if (depth_ < Capacity) {
      ancestors_[depth_] = view_.node();
    }
    ++depth_;
    view_.reset(child);
  }

  typename Tree::NodeView view_;
  std::size_t depth_ = 0;
  std::array<Node, Capacity> ancestors_;
};

}  // namespace pixie::tree
