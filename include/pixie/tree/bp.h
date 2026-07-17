#pragma once

#include <pixie/rmm/tree.h>
#include <pixie/tree.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief A tree class based on the balances parentheses (BP)
 * representation
 */
template <typename RMMTree>
class BPTree : public TreeBase<BPTree<RMMTree>> {
 private:
  const size_t num_bits_;
  RMMTree rmm_;

 public:
  /**
   * @brief A node class of BP tree
   */
  using Node = TreeNode;

  /**
   * @brief Constructor from an external array of uint64_t
   */
  explicit BPTree(const std::vector<std::uint64_t>& words, size_t tree_size)
      : num_bits_(2 * tree_size), rmm_(words, 2 * tree_size) {}

  /**
   * @brief Returns the root node
   */
  Node root_impl() const { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return num_bits_ / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const {
    return (node.pos + 2 == num_bits_) or rmm_.bit(node.pos + 1) == 0;
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   *    this method has O(d) time complexity!
   *
   *    TODO try make this faster
   */
  size_t degree_impl(const Node& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    Node child = first_child_impl(node);
    size_t child_count = 1;
    while (true) {
      if (is_last_child_impl(child)) {
        return child_count;
      }
      child = next_sibling_impl(child);
      child_count++;
    }
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = node.pos + 1;
    size_t num = node.number + 1;
    return Node(num, pos);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   *    this method has O(i) time complexity!
   *
   *    TODO try make this faster
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    size_t pos = rmm_.enclose(node.pos);
    size_t num = rmm_.rank1(pos);
    return Node(num, pos);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    size_t end = rmm_.close(node.pos);

    return end + 2 >= num_bits_ or rmm_.bit(end + 1) == 0;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    size_t pos = rmm_.close(node.pos) + 1;
    size_t num = rmm_.rank1(pos + 1) - 1;
    return Node(num, pos);
  }
};

std::vector<uint64_t> adj_to_bp(size_t tree_size,
                                const std::vector<std::vector<size_t>>& adj) {
  size_t bp_size = tree_size * 2;
  std::vector<uint64_t> bp((bp_size + 63) / 64, 0);
  std::vector<std::pair<size_t, size_t>> stack;
  stack.push_back(std::make_pair(0, 0));
  size_t pos = 0;
  bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
  while (!stack.empty()) {
    auto& [v, p] = stack.back();
    p++;
    if (p >= adj[v].size()) {
      pos++;
      stack.pop_back();
      continue;
    }
    pos++;
    bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
    stack.push_back(std::make_pair(adj[v][p], 0));
  }
  return bp;
}
}  // namespace pixie
