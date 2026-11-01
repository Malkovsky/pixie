#include <pixie/rmm/btree.h>
#include <pixie/tree/implementations.h>

#include "tree_family_benchmarks.h"

using pixie::benchmarks::RunBestTreeBenchmark;
using pixie::benchmarks::RunDirectTreeBenchmark;
using pixie::benchmarks::TreeLayout;
using pixie::benchmarks::TreeSizes;
using pixie::benchmarks::TreeWorkload;
using Tree = pixie::DFUDSTree<pixie::RmMBTree<>>;
using PointerTree = pixie::references::PointerTree;

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kDfs)
    ->Name("dfuds_best_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_best_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_best_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunBestTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_best_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kDfs)
    ->Name("dfuds_direct_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_direct_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_direct_random_path")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(RunDirectTreeBenchmark,
                   Tree,
                   TreeLayout::kDfuds,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_direct_random_path_children")
    ->Apply(TreeSizes);
