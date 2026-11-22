#pragma once

#include <pixie/integer_vector/packed.h>
#include <pixie/tree/detail/construction.h>
#include <pixie/tree/detail/tree_path.h>

#include <bit>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::references {

/**
 * @brief Owning ordered tree numbered in breadth-first order.
 * @details Each of two packed arrays holds n entries of bit_width(n-1) bits
 * (zero bits for n <= 1). Children are consecutive BFS ordinals. Unlike
 * PointerTree, degree and indexed child access do not scan siblings; no child
 * IDs, input numbering, or auxiliary path index are retained. All navigation
 * requires a valid handle from this tree.
 */
class BfsCsrTree {
 public:
  /// Zero-based BFS ordinal, independent of the source handle type.
  using Node = std::size_t;

  /// Copy the packed allocations independently.
  BfsCsrTree(const BfsCsrTree&) = default;
  /// Transfer ownership and leave the source empty.
  BfsCsrTree(BfsCsrTree&& other) noexcept
      : parents_(std::exchange(other.parents_, {})),
        offsets_(std::exchange(other.offsets_, {})) {}
  /// Replace both arrays together, including on self-assignment.
  BfsCsrTree& operator=(BfsCsrTree other) noexcept {
    std::swap(parents_, other.parents_);
    std::swap(offsets_, other.offsets_);
    return *this;
  }

  /// Copy valid adjacency rooted at zero: each row lists parent, then children
  /// in sibling order. IDs need not be in BFS or topological order.
  explicit BfsCsrTree(const std::vector<std::vector<Node>>& adjacency) {
    build(
        adjacency.size(), [] { return Node{0}; },
        [&](Node v, auto append) {
          const auto& row = adjacency[v];
          for (std::size_t i = 1; i < row.size(); ++i) {
            append(row[i]);
          }
        });
  }

  /// Copy a valid ordered tree. Empty sources are never navigated.
  template <tree::detail::OrderedTreeSource Source>
  explicit BfsCsrTree(const Source& source) {
    build(
        source.size(), [&] { return source.root(); },
        [&](auto v, auto append) {
          std::size_t remaining = source.degree(v);
          if (remaining) {
            auto child = source.first_child(v);
            append(child);
            while (--remaining) {
              child = source.next_sibling(child);
              append(child);
            }
          }
        });
  }

  /// Number of nodes.
  std::size_t size() const { return offsets_.size(); }
  /// Whether there is no root.
  bool empty() const { return offsets_.empty(); }
  /// Ordinal zero; equals the past-the-end ordinal size() on an empty tree.
  Node root() const { return 0; }
  /// BFS ordinal of a valid handle.
  std::size_t node_number(Node v) const { return v; }
  /// Whether this node is the root.
  bool is_root(Node v) const { return v == 0; }
  /// Parent ordinal, or the root itself for the root.
  Node parent(Node v) const { return is_root(v) ? v : parents_[v]; }
  /// Number of immediate children, in O(1).
  std::size_t degree(Node v) const { return offset_end(v) - offsets_[v]; }
  /// Whether the child interval is empty.
  bool is_leaf(Node v) const { return degree(v) == 0; }
  /// First child or v itself when v is a leaf.
  Node first_child(Node v) const {
    const auto start = offsets_[v];
    return start == offset_end(v) ? v : start + 1;
  }
  /// Child at index i in O(1); requires i < degree(v).
  Node child(Node v, std::size_t i) const { return offsets_[v] + 1 + i; }
  /// Whether v has no next sibling. The root has none.
  bool is_last_child(Node v) const {
    return is_root(v) || v == offset_end(parents_[v]);
  }
  /// Next sibling or v itself if none exists.
  Node next_sibling(Node v) const { return is_last_child(v) ? v : v + 1; }
  /// Visit immediate children by value, in sibling order. Exceptions propagate.
  template <class Visitor>
  void for_each_child(Node v, Visitor&& visitor) const {
    const auto end = offset_end(v);
    for (auto i = offsets_[v]; i < end; ++i) {
      visitor(Node{i + 1});
    }
  }

  /// Inclusive path from u to v, in traversal order.
  std::vector<Node> tree_path(Node u, Node v) const {
    std::vector<Node> out, scratch;
    tree_path(u, v, out, scratch);
    return out;
  }
  /// Replace out with the path and clear scratch, retaining capacities.
  /// Aliasing the buffers throws std::invalid_argument.
  void tree_path(Node u,
                 Node v,
                 std::vector<Node>& out,
                 std::vector<Node>& scratch) const {
    tree::detail::collect_tree_path(
        u, v, out, scratch, [](Node n) { return n; },
        [this](Node n) { return parent(n); });
  }

  /// Object and both padded allocations, excluding allocator overhead.
  std::size_t memory_usage_bytes() const {
    return sizeof(*this) + (parents_.memory_usage_bytes() - sizeof(parents_)) +
           (offsets_.memory_usage_bytes() - sizeof(offsets_));
  }

 private:
  using Packed = PackedIntegerVector<Node>;
  Packed parents_, offsets_;

  // offsets_[v] counts children of all earlier nodes. The last child of v
  // has ordinal offset_end(v); the final prefix sum n-1 is implicit.
  Node offset_end(Node v) const {
    return v + 1 < size() ? offsets_[v + 1] : size() - 1;
  }

  template <class Root, class Children>
  void build(std::size_t n, Root root_handle, Children children) {
    if (!n) {
      return;
    }
    const auto width = std::bit_width(n - 1);
    // Packed storage rounds the bit count up to a complete 64-bit word.
    constexpr auto max_bits = std::numeric_limits<Node>::max() - 63;
    if (width && n > max_bits / width) {
      throw std::length_error("BFS tree packed dimensions overflow");
    }

    std::vector<Node> entries(n);
    {
      std::vector<decltype(root_handle())> queue;
      queue.reserve(n);
      queue.push_back(root_handle());
      for (Node v = 0; v < n; ++v) {
        entries[v] = queue.size() - 1;
        children(queue[v], [&](auto child) { queue.push_back(child); });
      }
    }
    offsets_ = Packed(entries, width);
    // Fill parents backwards so the same construction buffer can supply the
    // offsets and then the parents. All child ordinals are greater than v.
    Node end = n - 1;
    for (Node v = n; v-- > 0;) {
      const Node start = entries[v];
      while (end > start) {
        entries[end--] = v;
      }
    }
    // The root needs no parent entry. Use its slot to keep the requested width
    // even for a star (Packed canonicalizes all-zero vectors to width zero).
    entries[0] = n - 1;
    parents_ = Packed(entries, width);
  }
};

}  // namespace pixie::references
