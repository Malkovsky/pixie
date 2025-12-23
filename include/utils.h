#pragma once

#include <queue>
#include <random>
#include <vector>

#include "louds_tree.h"

using pixie::LoudsNode;

std::vector<std::vector<size_t>> generate_random_tree(size_t tree_size,
                                                      std::mt19937_64& rng) {
  if (tree_size == 0) {
    return {};
  }
  std::vector<std::vector<size_t>> adj(tree_size);
  adj[0].push_back(0);
  for (size_t i = 1; i < tree_size; i++) {
    size_t parent = rng() % i;
    adj[i].push_back(parent);
    adj[parent].push_back(i);
  }
  return adj;
}

std::vector<std::vector<size_t>> bfs_order(
    size_t tree_size,
    const std::vector<std::vector<size_t>>& adj) {
  std::vector<std::vector<size_t>> bfs_adj(tree_size);
  std::queue<std::pair<size_t, size_t>> q;
  bfs_adj[0].push_back(0);
  q.push({0, 0});
  int cnt = 1;
  while (!q.empty()) {
    size_t old_v = q.front().first;
    size_t cur_v = q.front().second;
    q.pop();
    for (size_t i = 1; i < adj[old_v].size(); i++) {
      size_t old_u = adj[old_v][i];
      size_t cur_u = cnt++;
      q.push({old_u, cur_u});
      bfs_adj[cur_u].push_back(cur_v);
      bfs_adj[cur_v].push_back(cur_u);
    }
  }
  return bfs_adj;
}

std::vector<uint64_t> adj_to_louds(
    size_t tree_size,
    const std::vector<std::vector<size_t>>& adj) {
  size_t louds_size = tree_size * 2 - 1;
  std::vector<uint64_t> louds((louds_size + 63) / 64, 0);
  size_t pos = 0;
  for (size_t i = 0; i < tree_size; i++) {
    louds[pos >> 6] = louds[pos >> 6] | (1ULL << (pos & 63));
    pos += adj[i].size();
  }
  return louds;
}

struct AdjListNode {
  size_t number;
};

bool operator==(const AdjListNode& a, const LoudsNode& b) {
  return a.number == b.number;
}

bool operator==(const LoudsNode& b, const AdjListNode& a) {
  return a.number == b.number;
}

class AdjListTree {
 private:
  std::vector<std::vector<size_t>> adj;

 public:
  /**
   * @brief Constructor from adjacency list (root is 0)
   */
  explicit AdjListTree(const std::vector<std::vector<size_t>>& adjacency_list)
      : adj(adjacency_list) {}

  /**
   * @brief Returns the root node
   */
  AdjListNode root() const { return AdjListNode(0); }

  /**
   * @brief Checks if @p node is a leaf
   */
  bool is_leaf(const AdjListNode& node) const {
    return adj[node.number].size() <= 1;
  }

  /**
   * @brief Checks if @p node is a root
   */
  bool is_root(const AdjListNode& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of @p node
   */
  size_t degree(const AdjListNode& node) const {
    return adj[node.number].size() - 1;
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  AdjListNode child(const AdjListNode& node, size_t i) const {
    return AdjListNode(adj[node.number][i + 1]);
  }

  /**
   * @brief Returns the first child of @p node
   */
  AdjListNode first_child(const AdjListNode& node) const {
    return AdjListNode(adj[node.number][1]);
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  AdjListNode parent(const AdjListNode& node) const {
    return AdjListNode(adj[node.number][0]);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child(const AdjListNode& node) const {
    size_t p = parent(node).number;
    return adj[p].back() == node.number;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  AdjListNode next_sibling(const AdjListNode& node) const {
    return AdjListNode(node.number + 1);
  }
};
