#include <pixie/rmm/btree.h>
#include <pixie/tree/implementations.h>

#include "tree_family_benchmarks.h"

using pixie::benchmarks::RunBestTreeBenchmark;
using pixie::benchmarks::RunDirectTreeBenchmark;
using pixie::benchmarks::TreeLayout;
using pixie::benchmarks::TreeSizes;
using pixie::benchmarks::TreeWorkload;
using Tree = pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>;
using PointerTree = pixie::references::PointerTree;

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfs)
    ->Name("louds_best_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_best_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPath)
    ->Name("louds_best_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_best_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfs)
    ->Name("louds_direct_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_direct_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPath)
    ->Name("louds_direct_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_direct_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfs)
    ->Name("adjlist_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfsChildren)
    ->Name("adjlist_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPath)
    ->Name("adjlist_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("adjlist_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfs)
    ->Name("pointer_bfs_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kDfsChildren)
    ->Name("pointer_bfs_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPath)
    ->Name("pointer_bfs_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kLouds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("pointer_bfs_random_path_children")
    ->Apply(TreeSizes);

// Bounded region materialization versus replay of an already retained result.
static void LoudsScan(benchmark::State& state) {
  const size_t n = 255;
  std::vector<uint64_t> words((2 * n - 1 + 63) / 64, 0);
  size_t pos = 0;
  for (size_t i = 0; i < n; ++i) {
    words[pos / 64] |= uint64_t{1} << (pos % 64);
    pos += i < n / 2 ? 3 : 1;
  }
  Tree tree(words, n);
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
    words[pos / 64] |= uint64_t{1} << (pos % 64);
    pos += i < n / 2 ? 3 : 1;
  }
  Tree tree(words, n);
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
    words[pos / 64] |= uint64_t{1} << (pos % 64);
    pos += i < n / 2 ? 3 : 1;
  }
  Tree tree(words, n);
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
