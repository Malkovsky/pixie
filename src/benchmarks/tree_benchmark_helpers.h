#pragma once

#include <benchmark/benchmark.h>
#include <pixie/tree/path_view.h>

#include <cstddef>
#include <limits>
#include <optional>

#include "perf_query_window.h"
#include "shifted_query_pool.h"
#include "tree_traversal_workloads.h"

namespace pixie::benchmarks {

// Traverse through the public navigation API with O(1) auxiliary memory.
template <class Tree>
std::size_t TraverseDepthFirst(const Tree& tree) {
  auto node = tree.root();
  bool descending = true;
  std::size_t visited = 1;
  while (true) {
    auto observed = node;
    benchmark::DoNotOptimize(observed);
    if (descending && !tree.is_leaf(node)) {
      node = tree.first_child(node);
      ++visited;
      continue;
    }
    descending = false;
    if (tree.is_root(node)) {
      break;
    }
    if (tree.is_last_child(node)) {
      node = tree.parent(node);
    } else {
      node = tree.next_sibling(node);
      descending = true;
      ++visited;
    }
  }
  return visited;
}

// Same DFS decisions, with a view retaining query results until reset.
template <class Tree>
std::size_t TraverseDepthFirstWithView(const Tree& tree) {
  auto view = tree.navigate(tree.root());
  bool descending = true;
  std::size_t visited = 1;
  while (true) {
    auto observed = view.node();
    benchmark::DoNotOptimize(observed);
    if (descending && !view.is_leaf()) {
      view.reset(view.first_child());
      ++visited;
      continue;
    }
    descending = false;
    if (view.is_root()) {
      break;
    }
    if (view.is_last_child()) {
      view.reset(view.parent());
    } else {
      view.reset(view.next_sibling());
      descending = true;
      ++visited;
    }
  }
  return visited;
}

template <bool UseView = false,
          bool UsePath = false,
          bool Unpack = false,
          class Tree>
void RunDepthFirstBenchmark(benchmark::State& state,
                            const Tree& tree,
                            std::size_t tree_size,
                            std::span<const std::uint64_t> words = {}) {
  const auto traverse = [&] {
    if constexpr (UsePath) {
      DepthFirst<Unpack>(tree, words,
                         [](auto node) { benchmark::DoNotOptimize(node); });
      return tree.size();
    } else if constexpr (UseView) {
      return TraverseDepthFirstWithView(tree);
    } else {
      return TraverseDepthFirst(tree);
    }
  };
  auto& perf = PerfQueryWindow::instance();
  if (perf.active() && traverse() != tree_size) {
    state.SkipWithError("Warmup DFS did not visit every node");
    return;
  }
  perf.start();
  for (auto _ : state) {
    auto visited = traverse();
    benchmark::DoNotOptimize(visited);
    if (visited != tree_size) {
      state.SkipWithError("DFS did not visit every node");
      break;
    }
  }
  perf.stop(state.iterations() * tree_size);
  state.SetItemsProcessed(state.iterations() * tree_size);
  state.counters["nodes"] = tree_size;
  if (perf.active()) {
    state.counters["perf_nodes"] = perf.queries();
    state.counters["perf_phases"] = perf.phases();
    state.counters["perf_traversals"] = state.iterations();
    state.counters["perf_warmup_traversals"] = 1;
  }
}

enum class TreeWorkload {
  kDfs,
  kDfsChildren,
  kRandomPath,
  kRandomPathChildren
};

// Each random-path iteration is one complete root-to-leaf walk. Report both
// paths and visited nodes, since different paths have different lengths.
template <TreeWorkload Workload,
          bool Unpack = false,
          bool Indexed = false,
          class Tree>
void RunExpandedTreeBenchmark(benchmark::State& state,
                              const Tree& tree,
                              std::size_t tree_size,
                              std::span<const uint64_t> words = {}) {
  std::optional<experimental::TreeNavigationIndex<Tree>> index;
  if constexpr (Indexed) {
    static_assert(Workload == TreeWorkload::kRandomPath ||
                  Workload == TreeWorkload::kRandomPathChildren);
    index.emplace(tree);
    state.counters["navigation_index_bytes"] = index->memory_usage_bytes();
    state.counters["navigation_index_bytes_per_node"] =
        double(index->memory_usage_bytes()) / tree_size;
  }
  std::optional<ShiftedQueryPool> queries;
  if constexpr (Workload != TreeWorkload::kDfsChildren) {
    queries.emplace(0, std::numeric_limits<std::size_t>::max(), 43);
  }
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  const auto traverse = [&] {
    if constexpr (Workload == TreeWorkload::kDfsChildren) {
      if constexpr (Unpack) {
        return UnpackedDepthFirstInspect(tree, words, observe, observe);
      } else {
        return DepthFirstInspect(tree, observe, observe);
      }
    } else {
      if constexpr (Indexed) {
        return IndexedRandomRootToLeaf<Workload ==
                                       TreeWorkload::kRandomPathChildren>(
            *index, *queries, observe, observe);
      } else if constexpr (Unpack) {
        return UnpackedRandomRootToLeaf<Workload ==
                                        TreeWorkload::kRandomPathChildren>(
            tree, words, *queries, observe, observe);
      } else {
        return RandomRootToLeaf<Workload == TreeWorkload::kRandomPathChildren>(
            tree, *queries, observe, observe);
      }
    }
  };
  auto& perf = PerfQueryWindow::instance();
  if (perf.active()) {
    traverse();
    // Keep the measured query stream identical with and without perf warmup.
    if constexpr (Workload != TreeWorkload::kDfsChildren) {
      queries.emplace(0, std::numeric_limits<std::size_t>::max(), 43);
    }
  }
  std::size_t nodes = 0, children = 0;
  perf.start();
  for (auto _ : state) {
    auto work = traverse();
    nodes += work.nodes;
    children += work.children;
    benchmark::DoNotOptimize(work);
    if constexpr (Workload == TreeWorkload::kDfsChildren) {
      if (work.nodes != tree_size || work.children != tree_size - 1) {
        state.SkipWithError("DFS child inspection did not cover every edge");
        break;
      }
    }
  }
  perf.stop(nodes);
  state.SetItemsProcessed(nodes);
  state.counters["nodes"] = tree_size;
  state.counters["visited_nodes_per_iteration"] =
      benchmark::Counter(nodes, benchmark::Counter::kAvgIterations);
  state.counters["inspected_children_per_iteration"] =
      benchmark::Counter(children, benchmark::Counter::kAvgIterations);
  if constexpr (Workload != TreeWorkload::kDfsChildren) {
    state.counters["paths_per_second"] =
        benchmark::Counter(state.iterations(), benchmark::Counter::kIsRate);
    state.counters["query_pool_bytes"] = ShiftedQueryPool::kPoolBytes;
  }
  if (perf.active()) {
    state.counters["perf_nodes"] = perf.queries();
    state.counters["perf_phases"] = perf.phases();
  }
}

inline void TreeSizes(benchmark::internal::Benchmark* benchmark) {
  benchmark->ArgNames({"tree_size"})
      ->RangeMultiplier(2)
      ->Range(1ull << 8, 1ull << 26)
      ->Arg(10000000)
      ->Unit(benchmark::kMillisecond);
}

}  // namespace pixie::benchmarks
