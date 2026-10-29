#pragma once

/**
 * @file tree.h
 * @brief Common interface and node handle for rooted ordered trees.
 *
 * Include `<pixie/tree/implementations.h>` to use Pixie's succinct tree
 * encodings.
 */

#include <concepts>
#include <cstddef>
#include <utility>

namespace pixie {

/**
 * @brief Logical node handle shared by succinct rooted-tree encodings.
 */
struct TreeNode {
  /** @brief Logical node number in the encoding's traversal order. */
  std::size_t number;

  /** @brief Bit position representing the node in the succinct encoding. */
  std::size_t pos;

  /**
   * @brief Construct a node handle.
   * @param node_number Logical node number.
   * @param representation_position Position in the succinct representation.
   */
  TreeNode(std::size_t node_number, std::size_t representation_position)
      : number(node_number), pos(representation_position) {}
};

/**
 * @brief Required operations for a lightweight, per-node navigation view.
 * @details Behavioral and lifetime requirements are specified by
 * TreeBase::navigate.
 */
template <class View>
concept TreeNavigationView =
    std::copyable<View> &&
    requires(View view, const View current, TreeNode node, std::size_t i) {
      { current.node() } -> std::same_as<TreeNode>;
      { view.reset(node) } -> std::same_as<void>;
      { current.is_root() } -> std::same_as<bool>;
      { current.is_leaf() } -> std::same_as<bool>;
      { current.degree() } -> std::same_as<std::size_t>;
      { current.child(i) } -> std::same_as<TreeNode>;
      { current.first_child() } -> std::same_as<TreeNode>;
      { current.parent() } -> std::same_as<TreeNode>;
      { current.next_sibling() } -> std::same_as<TreeNode>;
      { current.is_last_child() } -> std::same_as<bool>;
      {
        current.for_each_child([](TreeNode) {})
      } -> std::same_as<void>;
    };

/**
 * @brief CRTP facade for rooted ordered trees.
 *
 * Implementations may use different succinct encodings, but expose the same
 * logical navigation operations through this facade.
 * Implementations and their selected backends must support concurrent const
 * queries without mutable shared navigation state. Query memoization belongs
 * to individual views; traversal scratch and auxiliary indexes are external.
 *
 * @see `<pixie/tree/implementations.h>` for the available succinct encodings.
 */
template <class Impl>
class TreeBase {
 public:
  using Node = TreeNode;

  /**
   * @brief Create an independent memoizing navigation view.
   * @param node Valid handle from this tree; invalid handles are undefined
   * behavior.
   * @details The tree and its immutable backing data must outlive the view and
   * remain unmoved. All query memoization belongs to the view, never the tree.
   * Independent views may query the same const tree concurrently. Concurrent
   * access to one view requires external synchronization.
   * @par Implementation requirement
   * navigate_impl(node) returns a copyable NodeView with node(), reset(node),
   * is_root(), is_leaf(), degree(), child(i), first_child(), parent(),
   * next_sibling(), is_last_child() and for_each_child(visitor). Queries follow
   * the matching TreeBase semantics. Only reset changes the current node and
   * invalidates memoized results; copies have independent memoization.
   */
  auto navigate(Node node) const {
    static_assert(TreeNavigationView<typename Impl::NodeView>);
    return impl().navigate_impl(node);
  }

  /**
   * @brief Object bytes plus owned index allocations.
   * @details Excludes borrowed representation words, caller-owned views,
   * external traversal/index objects and allocator metadata.
   * @par Implementation requirement
   * memory_usage_bytes_impl() reports this total without allocating or mutating
   * the tree. Available only when the selected backend supports accounting.
   */
  std::size_t memory_usage_bytes() const
    requires requires(const Impl& tree) { tree.memory_usage_bytes_impl(); }
  {
    return impl().memory_usage_bytes_impl();
  }

  /**
   * @brief Return the root node.
   * @return Handle identifying the root.
   */
  Node root() const { return impl().root_impl(); }

  /**
   * @brief Return the number of nodes.
   * @return Number of nodes in the tree.
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Check whether the tree contains no nodes.
   * @return `true` when `size() == 0`.
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Check whether @p node has no children.
   * @param node Valid node handle from this tree.
   * @return `true` when the node is a leaf.
   */
  bool is_leaf(const Node& node) const { return impl().is_leaf_impl(node); }

  /**
   * @brief Check whether @p node is the root.
   * @param node Valid node handle from this tree.
   * @return `true` when the node is the root.
   */
  bool is_root(const Node& node) const { return impl().is_root_impl(node); }

  /**
   * @brief Return the number of children of @p node.
   * @param node Valid node handle from this tree.
   * @return The node's child count.
   */
  std::size_t degree(const Node& node) const {
    return impl().degree_impl(node);
  }

  /**
   * @brief Return the first child of @p node.
   * @param node Valid non-leaf node handle from this tree.
   * @return Handle identifying the node's first child.
   */
  Node first_child(const Node& node) const {
    return impl().first_child_impl(node);
  }

  /**
   * @brief Return a child of @p node by zero-based index.
   * @param node Valid node handle from this tree.
   * @param index Child index in `[0, degree(node))`.
   * @return Handle identifying the requested child.
   */
  Node child(const Node& node, std::size_t index) const {
    return impl().child_impl(node, index);
  }

  /**
   * @brief Visit every immediate child in zero-based child-index order.
   * @param node Valid node handle from this tree; invalid handles are undefined
   * behavior.
   * @param visitor Callable invoked as an lvalue with a Node value for each
   * child in `[0, degree(node))`. Return values are ignored; leaves invoke it
   * zero times. Exceptions propagate and stop enumeration.
   * @details The callback is borrowed synchronously and never copied or
   * retained. Handles may be retained under the same lifetime rules as child():
   * the tree and its backing data must remain alive and unchanged. Nested
   * enumeration is supported. No heap allocation is required by enumeration
   * itself.
   * @par Implementation requirement
   * `for_each_child_impl(node, visitor)` must provide these semantics without
   * retaining callback or traversal state in the tree.
   */
  template <class Visitor>
  void for_each_child(const Node& node, Visitor&& visitor) const {
    impl().for_each_child_impl(node, std::forward<Visitor>(visitor));
  }

  /**
   * @brief Return the parent of @p node.
   * @param node Valid node handle from this tree.
   * @return The parent handle, or the root itself when @p node is the root.
   */
  Node parent(const Node& node) const { return impl().parent_impl(node); }

  /**
   * @brief Check whether @p node is its parent's last child.
   * @param node Valid non-root node handle from this tree.
   * @return `true` when no sibling follows @p node.
   */
  bool is_last_child(const Node& node) const {
    return impl().is_last_child_impl(node);
  }

  /**
   * @brief Return the next sibling of @p node.
   * @param node Valid non-root node that is not its parent's last child.
   * @return Handle identifying the next sibling.
   */
  Node next_sibling(const Node& node) const {
    return impl().next_sibling_impl(node);
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
