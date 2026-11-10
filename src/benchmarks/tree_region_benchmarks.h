#pragma once

#include <benchmark/benchmark.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/traversal_context.h>

#include <vector>

namespace pixie::benchmarks {
using LoudsRegionTree = LoudsTree<RankSelectSupportCompressed<>>;

// Bounded region materialization versus replay of an already retained result.
static void LoudsScan(benchmark::State& state) {
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
  LoudsRegionTree tree(words, n);
  pixie::TreeTraversalContext context(tree);
  for (auto _ : state) {
    auto result = context.scan(tree.root(), state.range(0));
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(LoudsScan)->Arg(64)->Arg(256)->Arg(1024);

static void LoudsDecodedReplay(benchmark::State& state) {
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
  LoudsRegionTree tree(words, n);
  pixie::TreeTraversalContext context(tree);
  auto result = context.fetch(tree.root(), 1024);
  for (auto _ : state) {
    result->events([](auto node, bool entering) {
      benchmark::DoNotOptimize(node);
      benchmark::DoNotOptimize(entering);
    });
  }
}
BENCHMARK(LoudsDecodedReplay);

// Same full enter/exit event contract as DecodedTree::events().
static void LoudsDirectEvents(benchmark::State& state) {
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
  LoudsRegionTree tree(words, n);
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
