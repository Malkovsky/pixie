#pragma once

#include <pixie/tree/detail/tree_parentheses.h>
#include <pixie/tree/detail/tree_path.h>

#include <cstddef>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

namespace pixie::tree {

/**
 * @brief CRTP navigation interface for the represented portion of a tree.
 * @details Concrete parts own their topology. Local handles belong to one part
 * and must be in [0, size()). The roots occupy the leading handles in source
 * order. Navigation works after the producer and source are destroyed.
 * Only source_tree() requires the source to remain alive, immutable and
 * unmoved. Degree, leaf status, sibling order and paths describe represented
 * edges. A represented leaf is a source leaf if children_complete() is true.
 * Implementations provide the corresponding *_impl hooks for size, root count,
 * degree, child, child range, completeness, parent, next sibling and memory
 * accounting. children_impl returns an allocation-free range in source order.
 * parent_impl returns the given node for local roots. next_sibling_impl returns
 * the given node at the end of its represented sibling group.
 */
template <class Impl, class SourceTree>
class TreePart {
 public:
  /** @brief Stable zero-based node number within this part. */
  using Node = std::size_t;

  /** @brief Positioned navigation borrowing this part. */
  class NodeView {
   public:
    /** @brief Return the current local handle. */
    Node node() const { return node_; }
    /** @brief Return the current local number. */
    std::size_t number() const { return node_; }
    /** @brief Move to a valid local handle. */
    void reset(Node node) { node_ = node; }
    /** @brief Whether the current node is a local root. */
    bool is_root() const { return part_->is_root(node_); }
    /** @brief Whether the current node has zero represented children. */
    bool is_leaf() const { return part_->is_leaf(node_); }
    /** @brief Whether every immediate source child is represented. */
    bool children_complete() const { return part_->children_complete(node_); }
    /** @brief Number of represented children. */
    std::size_t degree() const { return part_->degree(node_); }
    /** @brief Return the represented child at index i in [0, degree()). */
    Node child(std::size_t i) const { return part_->child(node_, i); }
    /** @brief Return the first represented child of a non-leaf. */
    Node first_child() const { return part_->first_child(node_); }
    /** @brief Return the represented parent, or the current local root. */
    Node parent() const { return part_->parent(node_); }
    /** @brief Return the successor of a non-last, non-root sibling. */
    Node next_sibling() const { return part_->next_sibling(node_); }
    /** @brief Whether a non-root node has no following represented sibling. */
    bool is_last_child() const { return part_->is_last_child(node_); }
    /** @brief Visit represented children in source order. */
    template <class Visitor>
    void for_each_child(Visitor&& visitor) const {
      part_->for_each_child(node_, visitor);
    }

   private:
    friend class TreePart;
    NodeView(const TreePart& part, Node node) : part_(&part), node_(node) {}
    const TreePart* part_;
    Node node_;
  };

  /** @brief Return an independent view that borrows this live, unmoved part. */
  NodeView navigate(Node node) const& { return NodeView(*this, node); }
  NodeView navigate(Node) const&& = delete;

  /** @brief Borrow the original tree. The caller must keep it alive and
   * unmoved. */
  const SourceTree& source_tree() const { return *source_; }
  /** @brief Number of represented nodes. */
  std::size_t size() const { return impl().size_impl(); }
  /** @brief Whether the part represents zero nodes. */
  bool empty() const { return size() == 0; }
  /** @brief Number of leading local roots. */
  std::size_t root_count() const { return impl().root_count_impl(); }
  /** @brief Return root i in source order. Throws for i outside [0,
   * root_count()). */
  Node root(std::size_t i = 0) const {
    if (i >= root_count()) {
      throw std::out_of_range("tree part root index out of range");
    }
    return i;
  }
  /** @brief Return the stable zero-based number of a local handle. */
  std::size_t number(Node node) const { return node; }
  /** @brief Return the handle for i. Throws for i outside [0, size()). */
  Node find_by_number(std::size_t i) const {
    if (i >= size()) {
      throw std::out_of_range("tree part node number out of range");
    }
    return i;
  }
  /** @brief Enumerate local handles in number order without allocation. */
  auto nodes() const { return std::views::iota(Node{0}, size()); }
  /** @brief Whether the node is a local root. */
  bool is_root(Node node) const { return node < root_count(); }
  /** @brief Whether the node has zero represented children. */
  bool is_leaf(Node node) const { return degree(node) == 0; }
  /** @brief Whether all immediate source children are represented. */
  bool children_complete(Node node) const {
    return impl().children_complete_impl(node);
  }
  /** @brief Number of represented immediate children. */
  std::size_t degree(Node node) const { return impl().degree_impl(node); }
  /** @brief Return the represented child at index i in [0, degree(node)). */
  Node child(Node node, std::size_t i) const {
    return impl().child_impl(node, i);
  }
  /** @brief Return the first represented child of a non-leaf. */
  Node first_child(Node node) const { return child(node, 0); }
  /** @brief Visit represented children in source order without allocation. */
  template <class Visitor>
  void for_each_child(Node node, Visitor&& visitor) const {
    for (Node child : impl().children_impl(node)) {
      visitor(Node(child));
    }
  }
  /** @brief Return the represented parent, or the node itself for a local root.
   */
  Node parent(Node node) const { return impl().parent_impl(node); }
  /** @brief Whether a non-root node has no following represented sibling. */
  bool is_last_child(Node node) const {
    return impl().next_sibling_impl(node) == node;
  }
  /** @brief Return the successor of a non-last, non-root sibling. */
  Node next_sibling(Node node) const { return impl().next_sibling_impl(node); }

  /** @brief Return the inclusive local path. Different local roots give no
   * path. */
  std::vector<Node> tree_path(Node u, Node v) const {
    std::vector<Node> out, scratch;
    tree_path(u, v, out, scratch);
    return out;
  }
  /** @brief Replace out with the local path using two distinct reusable
   * vectors. Scratch is empty on return. Parent handles precede their children.
   */
  void tree_path(Node u,
                 Node v,
                 std::vector<Node>& out,
                 std::vector<Node>& scratch) const {
    detail::collect_tree_path<true>(
        u, v, out, scratch, [](Node node) { return node; },
        [this](Node node) { return parent(node); });
  }
  /** @brief Stream a zero-based slice of the represented balanced parentheses.
   * The encodings of local roots are concatenated in order.
   */
  std::string parentheses(
      std::size_t start = 0,
      std::size_t count = std::numeric_limits<std::size_t>::max()) const {
    return detail::tree_parentheses(
        *this, root_count(), [this](std::size_t i) { return root(i); }, start,
        count);
  }
  /** @brief Object and owned topology storage in bytes, excluding the source.
   */
  std::size_t memory_usage_bytes() const {
    return impl().memory_usage_bytes_impl();
  }

 protected:
  /** @brief Borrow the source independently of the owned topology. */
  explicit TreePart(const SourceTree& source) : source_(&source) {}
  TreePart(const SourceTree&&) = delete;

 private:
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
  const SourceTree* source_;
};

}  // namespace pixie::tree
