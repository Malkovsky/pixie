#pragma once

#include <cstddef>
#include <vector>

namespace pixie::references {

/**
 * @brief Owning ordered-tree reference with exactly three pointers per node.
 * @details Nodes occupy one contiguous allocation in the input numbering order.
 * Each record stores first-child, next-sibling, and parent pointers; no degree
 * or node number is stored. Navigation requires a valid handle from this tree.
 * Handles remain valid until their allocation is destroyed or move-assigned
 * away; moving the owner transfers the allocation without rebasing pointers.
 * This reference uses pointer handles and does not inherit TreeBase.
 */
class PointerTree {
  struct Record {
    const Record* first = nullptr;
    const Record* next = nullptr;
    const Record* parent = nullptr;
  };

 public:
  /** @brief Borrowed handle into this tree's contiguous node allocation. */
  using Node = const Record*;

  /**
   * @brief Copy a well-formed parent-first adjacency list rooted at node zero.
   * @details adjacency[v][0] is v's parent; subsequent entries are its children
   * in order. The root's parent is zero. Empty input creates an empty tree.
   * Malformed topology or out-of-range node numbers violate the precondition.
   * The source may be destroyed immediately after construction.
   */
  explicit PointerTree(const std::vector<std::vector<std::size_t>>& adjacency)
      : nodes_(adjacency.size()) {
    for (std::size_t i = 0; i < adjacency.size(); ++i) {
      const auto& edges = adjacency[i];
      nodes_[i].parent = &nodes_[edges[0]];
      if (edges.size() > 1) {
        nodes_[i].first = &nodes_[edges[1]];
      }
      for (std::size_t j = 1; j + 1 < edges.size(); ++j) {
        nodes_[edges[j]].next = &nodes_[edges[j + 1]];
      }
    }
  }

  PointerTree(const PointerTree&) = delete;
  PointerTree& operator=(const PointerTree&) = delete;
  /** @brief Transfer the node allocation, preserving its handles. */
  PointerTree(PointerTree&&) noexcept = default;
  /** @brief Replace this allocation; invalidate its previous handles. */
  PointerTree& operator=(PointerTree&&) noexcept = default;

  /** @brief Number of owned nodes. */
  std::size_t size() const { return nodes_.size(); }
  /** @brief Whether there is no root. */
  bool empty() const { return nodes_.empty(); }
  /** @brief Root handle, or nullptr for an empty tree. */
  Node root() const { return empty() ? nullptr : nodes_.data(); }
  /** @brief Input node number in [0, size()); requires a valid handle. */
  std::size_t node_number(Node node) const { return node - nodes_.data(); }
  /** @brief Whether a valid node has no children. */
  bool is_leaf(Node node) const { return node->first == nullptr; }
  /** @brief Whether a valid node is the root. */
  bool is_root(Node node) const { return node == root(); }
  /** @brief Parent handle; the root is its own parent. */
  Node parent(Node node) const { return node->parent; }
  /** @brief First child, or the node itself for a leaf. */
  Node first_child(Node node) const { return node->first ? node->first : node; }
  /** @brief Whether no next sibling exists; true for the root. */
  bool is_last_child(Node node) const { return node->next == nullptr; }
  /** @brief Next sibling, or the node itself if it has no next sibling. */
  Node next_sibling(Node node) const { return node->next ? node->next : node; }
  /** @brief Count children by scanning sibling pointers: O(degree(node)). */
  std::size_t degree(Node node) const {
    std::size_t count = 0;
    for (auto child = node->first; child; child = child->next) {
      ++count;
    }
    return count;
  }
  /**
   * @brief Child at zero-based index i in O(i + 1) time.
   * @details Requires i in [0, degree(node)); invalid indices are undefined.
   */
  Node child(Node node, std::size_t i) const {
    auto result = node->first;
    while (i--) {
      result = result->next;
    }
    return result;
  }
  /**
   * @brief Visit immediate children in order, borrowing the visitor.
   * @details Passes handles by value and ignores callback return values.
   * Leaves invoke no callbacks. Exceptions propagate; nested calls are valid.
   * The allocation must remain alive and unchanged during callbacks and while
   * any retained handles are used.
   */
  template <class Visitor>
  void for_each_child(Node node, Visitor&& visitor) const {
    for (auto child = node->first; child; child = child->next) {
      visitor(child);
    }
  }
  /** @brief Object plus reserved node bytes, excluding allocator overhead. */
  std::size_t memory_usage_bytes() const {
    return sizeof(*this) + nodes_.capacity() * sizeof(Record);
  }

 private:
  std::vector<Record> nodes_;
};

}  // namespace pixie::references
