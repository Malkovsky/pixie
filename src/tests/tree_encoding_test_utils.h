#pragma once

#include <pixie/tree/detail/construction.h>
#include <pixie/utils.h>

// Test fixtures use the common builder over reference navigation. Succinct
// conversion tests compare their decoded producers against this fallback.
template <pixie::detail::TreeEncodingKind Kind>
std::vector<uint64_t> EncodeReference(
    size_t nodes,
    const std::vector<std::vector<size_t>>& adjacency) {
  if (nodes != adjacency.size()) {
    throw std::invalid_argument("Reference size mismatch");
  }
  const auto storage =
      pixie::detail::construct_tree_encoding<Kind>(AdjListTree(adjacency));
  const auto words = storage.as_words64();
  return {words.begin(), words.end()};
}
