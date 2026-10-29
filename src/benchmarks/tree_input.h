#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace pixie::benchmarks {

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

// Requires a valid BP encoding of n nodes, n < 2^32-1. Scratch space is one
// uint32_t degree per node plus an ancestor stack. Output uses Pixie's DFUDS
// convention: degree ones followed by zero, without an extra leading one.
inline std::vector<uint64_t> DfudsFromBp(std::span<const uint64_t> bp,
                                         size_t n) {
  if (!n) {
    return {};
  }
  std::vector<uint32_t> degree(n, 0), stack;
  size_t number = 0;
  for (size_t position = 0; position < 2 * n; ++position) {
    if ((bp[position / 64] >> (position % 64)) & 1) {
      if (!stack.empty()) {
        ++degree[stack.back()];
      }
      stack.push_back(number++);
    } else {
      stack.pop_back();
    }
  }
  std::vector<uint64_t> result((2 * n - 1 + 63) / 64, 0);
  size_t position = 0;
  for (auto count : degree) {
    for (size_t i = 0; i < count; ++i, ++position) {
      result[position / 64] |= uint64_t{1} << (position % 64);
    }
    ++position;
  }
  return result;
}

// Requires valid BP of n nodes. At a fixed depth, closing order equals opening
// order, so completed degree records can be appended to one packed stream per
// level. Concatenating those streams gives BFS order without a node-sized
// queue. Pixie's LOUDS convention is one followed by degree zeros per node.
inline std::vector<uint64_t> LoudsFromBp(std::span<const uint64_t> bp,
                                         size_t n) {
  if (!n) {
    return {};
  }
  struct Level {
    std::vector<uint64_t> words;
    size_t bits = 0;

    void append(size_t degree) {
      const size_t position = bits;
      bits += degree + 1;
      words.resize((bits + 63) / 64, 0);
      words[position / 64] |= uint64_t{1} << (position % 64);
    }
  };
  std::vector<Level> levels;
  std::vector<size_t> degree;
  for (size_t position = 0; position < 2 * n; ++position) {
    if ((bp[position / 64] >> (position % 64)) & 1) {
      if (!degree.empty()) {
        ++degree.back();
      }
      degree.push_back(0);
      if (levels.size() < degree.size()) {
        levels.emplace_back();
      }
    } else {
      const auto count = degree.back();
      degree.pop_back();
      levels[degree.size()].append(count);
    }
  }
  std::vector<uint64_t> result((2 * n - 1 + 63) / 64, 0);
  size_t base = 0;
  for (const auto& level : levels) {
    for (size_t word = 0; word < level.words.size(); ++word) {
      const size_t position = base + word * 64;
      result[position / 64] |= level.words[word] << (position % 64);
      if (position % 64 && position / 64 + 1 < result.size()) {
        result[position / 64 + 1] |= level.words[word] >> (64 - position % 64);
      }
    }
    base += level.bits;
  }
  return result;
}

}  // namespace pixie::benchmarks
