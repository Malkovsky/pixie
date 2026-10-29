#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rank_select/support.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/louds_children.h>
#include <pixie/tree/detail/louds_node_view.h>

#include <bit>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief A node class of LOUDS tree
 */
using LoudsNode = TreeNode;

/**
 * @brief A tree class based on the level order unary degree sequence (LOUDS)
 * representation
 * @tparam RankSelect Rank/select backend constructible from packed source words
 * and a bit count, with both select directions enabled by default. It must
 * provide size(), operator[], select() and select0() with RankSelectBase
 * semantics. Source words must remain alive and immutable for the tree
 * lifetime, even when the selected backend owns a copy: the word-level
 * navigation shortcuts read the original source directly.
 * @details Each BFS-ordered node is encoded as one followed by one zero per
 * child. A handle stores its zero-based BFS number and the position of its one.
 * Siblings therefore occupy consecutive one delimiters. The parent delimiter
 * is the nearest preceding one before the node's incoming zero edge. Degree,
 * sibling and parent navigation use a single-word search for those delimiters,
 * falling back to indexed select when the delimiter lies outside that word.
 */
template <class RankSelect = RankSelectSupport<>>
class LoudsTree : public TreeBase<LoudsTree<RankSelect>> {
 private:
  RankSelect bv;
  const std::uint64_t* words_;

  // Node delimiters are ones. Search only the remainder of this word so a
  // high-degree node cannot turn navigation into an unbounded linear scan.
  size_t next_delimiter(size_t position, size_t rank) const {
    const size_t start = position + 1;
    const auto word = words_[start / 64] >> (start % 64);
    return word ? start + std::countr_zero(word) : bv.select(rank);
  }

  size_t parent_delimiter(size_t zero_position, size_t rank) const {
    // zero_position itself is zero, so retaining it avoids a shift by 64.
    const auto word = words_[zero_position / 64] &
                      (~std::uint64_t{0} >> (63 - zero_position % 64));
    return word ? (zero_position / 64) * 64 + 63 - std::countl_zero(word)
                : bv.select(rank);
  }

  LoudsNode parent_from_edge(const LoudsNode& node, size_t zero_pos) const {
    const size_t parent_number = zero_pos - node.number;
    return LoudsNode(parent_number,
                     parent_delimiter(zero_pos, parent_number + 1));
  }

 public:
  /** @brief Independent, lightweight per-node navigation memo. */
  using NodeView = detail::LoudsNodeView<LoudsTree>;

 private:
  friend class detail::LoudsNodeView<LoudsTree>;

 public:
  /** @brief Implement TreeBase::navigate without mutating the tree. */
  NodeView navigate_impl(TreeNode node) const { return NodeView(*this, node); }

  /** @brief Object and owned backend allocations; excludes borrowed words. */
  size_t memory_usage_bytes_impl() const
    requires HasMemoryUsageBytes<RankSelect>
  {
    return sizeof(*this) + nested_owned_memory_bytes(bv);
  }

  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit LoudsTree(std::span<const uint64_t> louds, size_t tree_size)
      : bv(louds, 2 * tree_size - 1), words_(louds.data()) {}

  /** @brief Reject a temporary owning vector; this tree borrows its words. */
  LoudsTree(std::vector<std::uint64_t>&&, size_t) = delete;
  /** @brief Reject a temporary const owning vector for the same reason. */
  LoudsTree(const std::vector<std::uint64_t>&&, size_t) = delete;

  /**
   * @brief Returns the root node
   */
  LoudsNode root_impl() const { return LoudsNode(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (bv.size() + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const LoudsNode& node) const {
    return (node.pos + 1 == bv.size()) or bv[node.pos + 1];
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const LoudsNode& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const LoudsNode& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    return next_delimiter(node.pos, node.number + 2) - node.pos - 1;
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    detail::tree_unpack::louds_children(
        *this,
        std::span<const uint64_t>(words_,
                                  bv.size() / 64 + (bv.size() % 64 != 0)),
        bv.size(), node, 0, degree_impl(node), visitor);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  LoudsNode child_impl(const LoudsNode& node, size_t i) const {
    size_t zeros = node.pos + i + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns first child of a @p node
   */
  LoudsNode first_child_impl(const LoudsNode& node) const {
    size_t zeros = node.pos + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  LoudsNode parent_impl(const LoudsNode& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_edge(node, bv.select0(node.number));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const LoudsNode& node) const {
    size_t zero_pos = bv.select0(node.number);
    return bv[zero_pos + 1];
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  LoudsNode next_sibling_impl(const LoudsNode& node) const {
    size_t sibling_number = node.number + 1;
    return LoudsNode(sibling_number,
                     next_delimiter(node.pos, sibling_number + 1));
  }
};

}  // namespace pixie
