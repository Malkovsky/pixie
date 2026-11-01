#include <pixie/rmm/btree.h>
#include <pixie/tree/implementations.h>

#include "tree_conversion_benchmarks.h"
#include "tree_family_benchmarks.h"

using pixie::benchmarks::RunBestTreeBenchmark;
using pixie::benchmarks::RunDirectTreeBenchmark;
using pixie::benchmarks::TreeLayout;
using pixie::benchmarks::TreeSizes;
using pixie::benchmarks::TreeWorkload;
using Tree = pixie::BPTree<pixie::RmMBTree<>>;
using PointerTree = pixie::references::PointerTree;

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfs)
    ->Name("bp_best_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_best_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPath)
    ->Name("bp_best_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_best_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfs)
    ->Name("bp_direct_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_direct_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPath)
    ->Name("bp_direct_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_direct_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfs)
    ->Name("adjlist_preorder_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfsChildren)
    ->Name("adjlist_preorder_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPath)
    ->Name("adjlist_preorder_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   AdjListTree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPathChildren)
    ->Name("adjlist_preorder_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfs)
    ->Name("pointer_preorder_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kBp,
                   TreeWorkload::kDfsChildren)
    ->Name("pointer_preorder_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPath)
    ->Name("pointer_preorder_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   PointerTree,
                   TreeLayout::kBp,
                   TreeWorkload::kRandomPathChildren)
    ->Name("pointer_preorder_random_path_children")
    ->Apply(TreeSizes);
