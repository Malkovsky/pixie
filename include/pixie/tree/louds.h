#pragma once

#include <pixie/rank_select/support.h>
#include <pixie/tree.h>

#include <cstdint>
#include <span>

namespace pixie {

/**
 * @brief A node class of LOUDS tree
 */
using LoudsNode = TreeNode;

/**
 * @brief A tree class based on the level order unary degree sequence (LOUDS)
 * representation
 */
class LoudsTree : public TreeBase<LoudsTree> {
 private:
  RankSelectSupport<> bv;

 public:
  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit LoudsTree(std::span<const uint64_t> louds, size_t tree_size)
      : bv(louds, 2 * tree_size - 1) {}

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
    return bv.select(node.number + 2) - node.pos - 1;
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
    size_t zero_pos = bv.select0(node.number);
    size_t parent_number = zero_pos - node.number;
    return LoudsNode(parent_number, bv.select(parent_number + 1));
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
    return LoudsNode(sibling_number, bv.select(sibling_number + 1));
  }
};

}  // namespace pixie
