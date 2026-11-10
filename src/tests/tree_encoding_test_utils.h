#pragma once

#include <pixie/tree/detail/construction.h>
#include <pixie/utils.h>

// LOUDS fixtures use an independent BFS encoder. Other formats use reference
// navigation to exercise the generic construction path.
template <pixie::detail::TreeEncodingKind Kind>
std::vector<uint64_t> EncodeReference(
    size_t nodes,
    const std::vector<std::vector<size_t>>& adjacency) {
  if (nodes != adjacency.size()) {
    throw std::invalid_argument("Reference size mismatch");
  }
  if constexpr (Kind == pixie::detail::TreeEncodingKind::kLOUDS) {
    if (!nodes) {
      return {};
    }
    std::vector<uint64_t> words((2 * nodes - 1 + 63) / 64, 0);
    std::vector<size_t> queue{0};
    size_t position = 0;
    for (size_t i = 0; i < queue.size(); ++i) {
      const auto& edges = adjacency[queue[i]];
      for (size_t j = 1; j < edges.size(); ++j) {
        queue.push_back(edges[j]);
        words[position / 64] |= uint64_t{1} << (position % 64);
        ++position;
      }
      ++position;
    }
    return words;
  }
  const auto storage =
      pixie::detail::construct_tree_encoding<Kind>(AdjListTree(adjacency));
  const auto words = storage.as_words64();
  return {words.begin(), words.end()};
}
