#pragma once

#include <benchmark/benchmark.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/traversal_context.h>

#include <vector>

namespace pixie::benchmarks {
using LoudsRegionTree = tree::LoudsTree<RankSelectSupportCompressed<>>;

static LoudsRegionTree region_tree() {
  const size_t n = 255;
  std::vector<uint64_t> words((2 * n - 1 + 63) / 64, 0);
  size_t pos = 0;
  for (size_t i = 0; i < n; ++i) {
    const size_t degree = i < n / 2 ? 2 : 0;
    for (size_t j = 0; j < degree; ++j, ++pos) {
      words[pos / 64] |= uint64_t{1} << (pos % 64);
    }
    ++pos;
  }
  return LoudsRegionTree(words, n);
}

static void LoudsFetch(benchmark::State& state) {
  auto tree = region_tree();
  pixie::tree::TreeTraversalContext context(tree);
  for (auto _ : state) {
    auto result = context.fetch(tree.root(), state.range(0));
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(LoudsFetch)->Arg(64)->Arg(256)->Arg(1024);

static void LoudsInspectFetchedPart(benchmark::State& state) {
  auto tree = region_tree();
  pixie::tree::TreeTraversalContext context(tree);
  auto result = context.fetch(tree.root(), 1024);
  for (auto _ : state) {
    // Consume actual node and child IDs. Per-ID barriers can leave a
    // materialized ID in memory while forcing a computed ID to be produced.
    size_t checksum = 0;
    for (auto node : result.nodes()) {
      checksum += result.source_node(node).number;
      result.for_each_child(node, [&](auto child) { checksum += child; });
    }
    benchmark::DoNotOptimize(checksum);
  }
}
BENCHMARK(LoudsInspectFetchedPart);

// Whole-source traversal; unlike snapshot inspection this emits DFS events.
static void LoudsDirectEvents(benchmark::State& state) {
  auto tree = region_tree();
  for (auto _ : state) {
    auto node = tree.root();
    bool entering = true;
    while (true) {
      benchmark::DoNotOptimize(node);
      benchmark::DoNotOptimize(entering);
      if (entering && !tree.is_leaf(node)) {
        node = tree.first_child(node);
        continue;
      }
      if (entering) {
        entering = false;
        continue;
      }
      if (tree.is_root(node)) {
        break;
      }
      if (!tree.is_last_child(node)) {
        node = tree.next_sibling(node);
        entering = true;
      } else {
        node = tree.parent(node);
      }
    }
  }
}
BENCHMARK(LoudsDirectEvents);

}  // namespace pixie::benchmarks
