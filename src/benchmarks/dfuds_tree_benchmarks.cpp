#include <benchmark/benchmark.h>
#include <pixie/dfuds_tree.h>
#include <pixie/utils.h>

#include <random>

using Node = pixie::DFUDSTree::Node;
using pixie::adj_to_dfuds;
using pixie::DFUDSTree;

/**
 * DFS with O(1) extra memory
 */
static void BM_DfudsTreeDFS(benchmark::State& state) {
  size_t tree_size = state.range(0);
  std::mt19937_64 rng(42);

  for (auto _ : state) {
    state.PauseTiming();

    std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
    adj = dfs_order(tree_size, adj);
    std::vector<uint64_t> dfuds = adj_to_dfuds(tree_size, adj);
    DFUDSTree tree(dfuds, tree_size);

    Node cur = tree.root();
    bool above = 1;

    state.ResumeTiming();

    benchmark::DoNotOptimize(cur);

    while (true) {
      if (above) {
        if (tree.is_leaf(cur)) {
          above = 0;
        } else {
          cur = tree.first_child(cur);
        }
        benchmark::DoNotOptimize(cur);
      } else {
        if (tree.is_last_child(cur)) {
          cur = tree.parent(cur);
          if (tree.is_root(cur)) {
            break;
          }
          benchmark::DoNotOptimize(cur);
        } else {
          cur = tree.next_sibling(cur);
          above = 1;
          benchmark::DoNotOptimize(cur);
        }
      }
    }
  }
}

BENCHMARK(BM_DfudsTreeDFS)
    ->ArgNames({"tree_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 8, 1ull << 18)
    ->Iterations(100);

BENCHMARK(BM_DfudsTreeDFS)
    ->ArgNames({"tree_size"})
    ->RangeMultiplier(2)
    ->Range(1ull << 18, 1ull << 26)
    ->Iterations(10);
