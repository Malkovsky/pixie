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
