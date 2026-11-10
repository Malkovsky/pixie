#pragma once

#include <pixie/rmm/btree.h>
#include <pixie/tree/implementations.h>

#include <string>

#include "tree_family_benchmarks.h"

namespace pixie::benchmarks {
// Both executables register the same implementations for each workload.
template <TreeWorkload Workload>
void RegisterTreeWorkload(const char* family) {
  using Bp = BPTree<RmMBTree<>>;
  using Dfuds = DFUDSTree<RmMBTree<>>;
  using Louds = LoudsTree<RankSelectSupportCompressed<>>;
  using Pointer = references::PointerTree;
  auto add = [&](const char* implementation, auto run) {
    benchmark::RegisterBenchmark(
        (std::string(family) + "/" + implementation).c_str(), run)
        ->Apply(TreeSizes);
  };
  add("bp_best", &RunBestTreeBenchmark<Bp, TreeLayout::kBp, Workload>);
  add("bp_direct", &RunDirectTreeBenchmark<Bp, TreeLayout::kBp, Workload>);
  add("dfuds_best", &RunBestTreeBenchmark<Dfuds, TreeLayout::kDfuds, Workload>);
  add("dfuds_direct",
      &RunDirectTreeBenchmark<Dfuds, TreeLayout::kDfuds, Workload>);
  add("louds_best", &RunBestTreeBenchmark<Louds, TreeLayout::kLouds, Workload>);
  add("louds_direct",
      &RunDirectTreeBenchmark<Louds, TreeLayout::kLouds, Workload>);
  add("adjlist_preorder",
      &RunDirectTreeBenchmark<AdjListTree, TreeLayout::kBp, Workload>);
  add("adjlist_bfs",
      &RunDirectTreeBenchmark<AdjListTree, TreeLayout::kLouds, Workload>);
  add("pointer_preorder",
      &RunDirectTreeBenchmark<Pointer, TreeLayout::kBp, Workload>);
  add("pointer_bfs",
      &RunDirectTreeBenchmark<Pointer, TreeLayout::kLouds, Workload>);
}
}  // namespace pixie::benchmarks
