#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmm/tree.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/dfuds_node_view.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief A tree class based on the depth-first unary degree sequence (DFUDS)
 * representation
 */
template <typename RMMTree>
class DFUDSTree : public TreeBase<DFUDSTree<RMMTree>> {
 private:
  const size_t num_bits_;
  RMMTree rmm_;

  // Sibling edges match the parent's unary ones from right to left. The last
  // child's edge is therefore the first one of that run, preceded by zero
  // (or at bit zero for the root's degree run).
  bool starts_degree_run(size_t open) const {
    return open == 0 || rmm_.bit(open - 1) == 0;
  }

  // For a non-root first child, the parent's degree run ends in the matching
  // adjacent "10" pair immediately before the child's record.
  size_t incoming_open(const TreeNode& node) const {
    return rmm_.bit(node.pos - 2) ? node.pos - 2 : rmm_.open(node.pos - 1);
  }

  TreeNode parent_from_open(const TreeNode& node, size_t open) const {
    // [open, node.pos) is balanced: half its bits are zero delimiters.
    // Subtract those records from the child's preorder number.
    const auto number = node.number - (node.pos - open) / 2;
    // For a last child, the matching one is already the parent's record start.
    const auto pos = number == 0               ? 0
                     : starts_degree_run(open) ? open
                                               : rmm_.select0(number) + 1;
    return TreeNode(number, pos);
  }

 public:
  using Node = TreeNode;

  /** @brief Independent, lightweight per-node navigation memo. */
  using NodeView = detail::DfudsNodeView<DFUDSTree>;

 private:
  friend class detail::DfudsNodeView<DFUDSTree>;

 public:
  /** @brief Implement TreeBase::navigate without mutating the tree. */
  NodeView navigate_impl(TreeNode node) const { return NodeView(*this, node); }

  /** @brief Object and owned backend allocations; excludes borrowed words. */
  size_t memory_usage_bytes_impl() const
    requires HasMemoryUsageBytes<RMMTree>
  {
    return sizeof(*this) + nested_owned_memory_bytes(rmm_);
  }

  /**
   * @brief Constructor from an external array of uint64_t
   * @details The source words must remain alive and immutable for this tree
   * and its traversals. No padding beyond the logical final word is required.
   *
   * @param dfuds_sequence parenthesis sequence in dfuds representation
   */
  explicit DFUDSTree(std::span<const std::uint64_t> dfuds_sequence,
                     size_t tree_size)
      : num_bits_(2 * tree_size - 1), rmm_(dfuds_sequence, 2 * tree_size - 1) {}

  /** @brief Reject a temporary owning vector; this tree borrows its words. */
  DFUDSTree(std::vector<std::uint64_t>&&, size_t) = delete;
  /** @brief Reject a temporary const owning vector for the same reason. */
  DFUDSTree(const std::vector<std::uint64_t>&&, size_t) = delete;

  /**
   * @brief Returns the root node
   */
  Node root_impl() const { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (num_bits_ + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const { return rmm_.bit(node.pos) == 0; }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const Node& node) const {
    return rmm_.select0(node.number + 1) - node.pos;
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = rmm_.select0(node.number + 1);
    size_t num = node.number + 1;
    return Node(num, pos + 1);
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    auto parent = NodeView(*this, node);
    const auto count = parent.degree();
    if (count == 0) {
      return;
    }
    auto child = NodeView(*this, parent.first_child());
    for (size_t i = 0; i < count; ++i) {
      visitor(child.node());
      if (i + 1 < count) {
        child.reset(child.next_sibling());
      }
    }
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    // A leaf occupies exactly one zero and the next record is its sibling.
    if (is_leaf_impl(node)) {
      return Node(node.number + 1, node.pos + 1);
    }
    size_t end = rmm_.fwdsearch(node.pos, -1);
    size_t pos = end + 1;
    // Preorder numbering advances by the subtree's node count.
    size_t num = node.number + (pos - node.pos + 1) / 2;
    return Node(num, pos);
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_open(node, incoming_open(node));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    return starts_degree_run(incoming_open(node));
  }
};

std::vector<uint64_t> adj_to_dfuds(
    size_t tree_size,
    const std::vector<std::vector<size_t>>& adj) {
  size_t dfuds_size = tree_size * 2 - 1;
  std::vector<uint64_t> dfuds((dfuds_size + 63) / 64, 0);
  std::vector<size_t> stack;
  stack.push_back(0);
  size_t pos = 0;
  while (!stack.empty()) {
    auto v = stack.back();
    stack.pop_back();
    size_t edge_count = adj[v].size();
    for (size_t i = 0; i < edge_count - 1; ++i) {  // edge 0 goes to parent
      dfuds[pos >> 6] = dfuds[pos >> 6] | (1ULL << (pos & 63));
      pos++;
      stack.push_back(adj[v][edge_count - 1 - i]);
    }
    pos++;
  }
  return dfuds;
}

}  // namespace pixie
