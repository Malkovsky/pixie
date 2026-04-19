#pragma once

#include <pixie/rmm_tree.h>

#include <cstdint>

namespace pixie {
/**
 * @brief A node class of BP tree
 */
struct BpNode {
  size_t number;
  size_t pos;

  BpNode(size_t node_number, size_t bp_pos)
      : number(node_number), pos(bp_pos) {}
};

/**
 * @brief A tree class based on the balances parentheses (BP)
 * representation
 */
class BpTree {
 private:
  const size_t num_bits_;
  RmMTree rmm;

 public:
  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit BpTree(const std::vector<std::uint64_t>& words, size_t tree_size)
      : num_bits_(2 * tree_size), rmm(words, 2 * tree_size) {}

  /**
   * @brief Returns the root node
   */
  static BpNode root() { return BpNode(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size() const { return num_bits_ / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf(const BpNode& node) const {
    return (node.pos + 1 == num_bits_) or rmm.bit(node.pos + 1) == 0;
  }

  /**
   * @brief Indicates if @p node is a root
   */
  static bool is_root(const BpNode& node) { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   *    this method has O(d) time complexity!
   *
   *    TODO try make this faster
   */
  size_t degree(const BpNode& node) const {
    if (is_leaf(node)) {
      return 0;
    }
    BpNode child = first_child(node);
    size_t child_count = 1;
    while (true) {
      if (is_last_child(child)) {
        return child_count;
      }
      child = next_sibling(child);
      child_count++;
    }
  }

  /**
   * @brief Returns first child of a @p node
   */
  static BpNode first_child(const BpNode& node) {
    size_t pos = node.pos + 1;
    size_t num = node.number + 1;
    return BpNode(num, pos);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   *    this method has O(i) time complexity!
   *
   *    TODO try make this faster
   */
  BpNode child(const BpNode& node, size_t i) const {
    BpNode child = first_child(node);
    while (i--) {
      child = next_sibling(child);
    }
    return child;
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  BpNode parent(const BpNode& node) const {
    if (node.number == 0) {
      return root();
    }
    size_t pos = rmm.enclose(node.pos + 1);
    size_t num = rmm.rank1(pos) - 1;
    return BpNode(num, pos);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child(const BpNode& node) const {
    size_t end = rmm.close(node.pos + 1);
    if (end + 1 >= num_bits_) {
      return true;
    }
    return rmm.bit(end + 1) == 0;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  BpNode next_sibling(const BpNode& node) const {
    size_t pos = rmm.close(node.pos + 1) + 1;
    size_t num = rmm.rank1(pos + 1) - 1;
    return BpNode(num, pos);
  }
};
}  // namespace pixie
