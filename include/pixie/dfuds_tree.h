#pragma once

#include <pixie/rmm_tree.h>

#include <cstdint>

namespace pixie {

/**
 * @brief A tree class based on the depth-first unary degree sequence (DFUDS)
 * representation
 */
class DFUDSTree {
 private:
  const size_t num_bits_;
  RmMTree rmm_;

 public:
  struct Node {
    size_t number;

    size_t pos;

    /**
     * @brief A node class of DFUDS tree
     */
    Node(size_t node_number, size_t dfuds_pos)
        : number(node_number), pos(dfuds_pos) {}
  };
  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit DFUDSTree(const std::vector<std::uint64_t>& words, size_t tree_size)
      : num_bits_(2 * tree_size - 1), rmm_(words, 2 * tree_size - 1) {}

  /**
   * @brief Returns the root node
   */
  static Node root() { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size() const { return (num_bits_ + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf(const Node& node) const {
    return (node.pos + 1 == num_bits_) or rmm_.bit(node.pos) == 0;
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree(const Node& node) const {
    return rmm_.select0(node.number + 1) - node.pos;
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child(const Node& node) {
    size_t pos = rmm_.select0(node.number + 1);
    size_t num = node.number + 1;
    return Node(num, pos + 1);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  Node child(const Node& node, size_t i) const {
    size_t pos = rmm_.close(rmm_.select0(node.number + 1) - i) + 1;
    size_t num = rmm_.rank0(pos);
    return Node(num, pos);
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling(const Node& node) const {
    size_t end = rmm_.fwdsearch(node.pos, -1);
    size_t pos = end + 1;
    size_t num = rmm_.rank0(pos);
    return Node(num, pos);
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent(const Node& node) const {
    if (node.number == 0) {
      return root();
    }
    size_t open = rmm_.open(node.pos);
    size_t rank = rmm_.rank0(open);
    size_t pos = rmm_.select0(rank) + 1;  // some overflow-related magic here
    size_t num = rmm_.rank0(pos);
    return Node(num, pos);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child(const Node& node) const {
    size_t end = rmm_.fwdsearch(node.pos, -1);
    size_t pos = end + 1;
    size_t op = rmm_.open(node.pos);
    size_t op2 = rmm_.open(pos);
    return pos == num_bits_ || op != op2 + 1;
  }
};
}  // namespace pixie
