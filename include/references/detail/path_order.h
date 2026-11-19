#pragma once

#include <cstddef>
#include <vector>

namespace pixie::references::detail {

// A well-formed adjacency row starts with the parent, then ordered children.
// Keep input IDs as keys unless an edge violates parent-before-child order.
inline std::vector<std::size_t> optional_preorder(
    const std::vector<std::vector<std::size_t>>& adj) {
  for (std::size_t i = 1; i < adj.size(); ++i) {
    if (adj[i][0] < i) {
      continue;
    }
    std::vector<std::size_t> rank(adj.size()), pending{0};
    std::size_t next = 0;
    while (!pending.empty()) {
      const auto v = pending.back();
      pending.pop_back();
      rank[v] = next++;
      const auto& row = adj[v];
      pending.insert(pending.end(), row.rbegin(), row.rend() - 1);
    }
    return rank;
  }
  return {};
}

}  // namespace pixie::references::detail
