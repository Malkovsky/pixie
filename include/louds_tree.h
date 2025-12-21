#pragma once

#include <cstdint>
#include <span>

#include "bitvector.h"

namespace pixie {

/**
 * @brief A node class of LOUDS tree
 */
struct LoudsNode {
  size_t number;
  size_t pos;

  LoudsNode(size_t node_number, size_t louds_pos)
      : number(node_number), pos(louds_pos) {}
};

/**
 * @brief A tree class based on the level order unary degree sequence (LOUDS)
 * representation
 */
class LoudsTree {
 private:
  BitVector bv;

 public:
  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit LoudsTree(std::span<const uint64_t> louds, size_t tree_size)
      : bv(louds, 2 * tree_size - 1) {}

  /**
   * @brief Returns the root node
   */
  LoudsNode root() const { return LoudsNode(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size() const { return (bv.size() + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf(const LoudsNode& node) const {
    return (node.pos + 1 == bv.size()) or bv[node.pos + 1];
  }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree(const LoudsNode& node) const {
    if (is_leaf(node)) {
      return 0;
    }
    return bv.select(node.number + 2) - node.pos - 1;
  }

  /**
   * @brief Returns the i-child of a @p node
   * Indexing starts at 0
   */
  LoudsNode child(const LoudsNode& node, size_t i) const {
    size_t zeros = node.pos + i + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  LoudsNode parent(const LoudsNode& node) const {
    if (node.number == 0) {
      return root();
    }
    size_t zero_pos = bv.select0(node.number);
    size_t parent_number = zero_pos - node.number;
    return LoudsNode(parent_number, bv.select(parent_number + 1));
  }
};

}  // namespace pixie
