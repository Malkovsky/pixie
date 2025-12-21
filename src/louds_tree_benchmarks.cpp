#include <benchmark/benchmark.h>

#include <numeric>
#include <random>
#include <stack>

#include "louds_tree.h"
#include "utils.h"

using pixie::LoudsNode;
using pixie::LoudsTree;

/**
 * DFS with O(1) extra memory
 */
static void BM_LoudsTreeDFS(benchmark::State& state) {
  size_t tree_size = state.range(0);
  std::mt19937_64 rng(42);
  std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
  adj = bfs_order(tree_size, adj);
  std::vector<uint64_t> louds = adj_to_louds(tree_size, adj);
  LoudsTree tree(louds, tree_size);

  for (auto _ : state) {
    LoudsNode cur = tree.root();
    bool above = 1;

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

/**
 * DFS with O(1) extra memory
 */
static void BM_AdjListTreeDFS(benchmark::State& state) {
  size_t tree_size = state.range(0);
  std::mt19937_64 rng(42);
  std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
  adj = bfs_order(tree_size, adj);
  AdjListTree tree(adj);

  for (auto _ : state) {
    AdjListNode cur = tree.root();
    bool above = 1;

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

BENCHMARK(BM_LoudsTreeDFS)
    ->ArgNames({"tree_size"})
    ->RangeMultiplier(2)
    ->Range(8, 1ull << 22);

BENCHMARK(BM_AdjListTreeDFS)
    ->ArgNames({"tree_size"})
    ->RangeMultiplier(2)
    ->Range(8, 1ull << 22);
