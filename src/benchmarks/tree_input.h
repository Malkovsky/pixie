#pragma once

#include <pixie/tree/bp.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace pixie::benchmarks {

// Fixture adapter for raw generator output. Build a borrowed BP source index,
// then use the library's common conversion producer and destination writer.
// This setup includes source-index construction and an output vector copy;
// encoding-only conversion benchmarks call construct_tree_encoding directly.
template <detail::TreeEncodingKind Kind>
std::vector<uint64_t> ConvertBpEncoding(std::span<const uint64_t> bp,
                                        size_t n) {
  const typename BPTree<RmMBTree<>>::TreeView source(bp, n);
  const auto storage = detail::construct_tree_encoding<Kind>(source);
  const auto words = storage.as_words64();
  return {words.begin(), words.end()};
}

// Benchmark-only generators for the registered sizes (at most 2^26 nodes).
// Reproduce generate_random_tree's parent[i] = rng() % i and ascending
// insertion order exactly. Two temporary uint32_t arrays replace the
// vector-of-vectors graph and its separately renumbered copy.
inline std::vector<uint64_t> GenerateRandomBp(size_t n, uint64_t seed = 42) {
  constexpr uint32_t none = std::numeric_limits<uint32_t>::max();
  if (n >= none) {
    throw std::length_error("Tree benchmark generator requires n < 2^32-1");
  }
  if (!n) {
    return {};
  }
  std::vector<uint32_t> first(n, none), next(n);
  std::mt19937_64 rng(seed);
  for (size_t i = 1; i < n; ++i) {
    next[i] = rng() % i;
  }
  // Reuse the parent array for sibling links. Descending insertion preserves
  // ascending child order; parent[i] < i means its parent slot is still intact.
  for (size_t i = n - 1; i > 0; --i) {
    const auto parent = next[i];
    next[i] = first[parent];
    first[parent] = i;
  }
  std::vector<uint64_t> bp((2 * n + 63) / 64, 0);
  std::vector<uint32_t> stack{0};
  size_t position = 0;
  bp[0] = 1;
  while (!stack.empty()) {
    const auto node = stack.back();
    const auto child = first[node];
    if (child == none) {
      ++position;
      stack.pop_back();
    } else {
      first[node] = next[child];
      ++position;
      bp[position / 64] |= uint64_t{1} << (position % 64);
      stack.push_back(child);
    }
  }
  return bp;
}

}  // namespace pixie::benchmarks
