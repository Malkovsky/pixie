#pragma once

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "tree_dfs/optimized.h"
#include "tree_traversal_workloads.h"

namespace pixie::benchmarks {
// The immutable index gets the part of the budget not reserved for the walk.
template <class Tree>
auto MakeWalkIndex(const Tree& tree,
                   size_t budget = tree::kDefaultTraversalBytes) {
  const auto reserved =
      tree::TreeTraversalContext(tree, budget).budget_usage_bytes();
  return tree::RootNavigationIndex<Tree>(tree, budget - reserved);
}

enum class CoverShape {
  kRandomRecursive,
  kStar,
  kChain,
  kComb,
  kBalanced,
  kCatalan,
  kBinaryCatalan
};

// Source IDs are insertion IDs, not preorder ranks. Registered sizes fit in
// uint32_t. The source is released before any query timing begins.
class CoverSource {
 public:
  CoverSource(size_t n, CoverShape shape, uint64_t seed)
      : first_(n, kNone), next_(n, kNone), degree_(n, 0) {
    std::mt19937_64 random(seed);
    if (shape == CoverShape::kBinaryCatalan && n) {
      // Uniform full-binary preorder with n internal nodes; discard nulls.
      std::vector<int8_t> steps(2 * n + 1, -1);
      std::fill_n(steps.begin(), n, 1);
      std::shuffle(steps.begin(), steps.end(), random);
      int level = 0, minimum = 0;
      size_t start = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        level += steps[i];
        if (level < minimum) {
          minimum = level;
          start = i + 1;
        }
      }
      std::vector<uint32_t> slots{kNone};
      uint32_t id = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        auto p = slots.back();
        slots.pop_back();
        if (steps[(start + i) % steps.size()] > 0) {
          next_[id] = p;
          if (p != kNone) {
            ++degree_[p];
          }
          slots.push_back(id);
          slots.push_back(id++);
        }
      }
    } else if (shape == CoverShape::kCatalan && n) {
      // Cycle lemma: rotate n-1 up and n down steps after the first
      // global minimum. Dropping the last down step gives a uniform Dyck path.
      std::vector<int8_t> steps(2 * n - 1, -1);
      std::fill_n(steps.begin(), n - 1, 1);
      std::shuffle(steps.begin(), steps.end(), random);
      int level = 0, minimum = 0;
      size_t start = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        level += steps[i];
        if (level < minimum) {
          minimum = level;
          start = i + 1;
        }
      }
      uint32_t node = 0, id = 0;
      for (size_t i = 0; i + 1 < steps.size(); ++i) {
        if (steps[(start + i) % steps.size()] == 1) {
          next_[++id] = node;
          ++degree_[node];
          node = id;
        } else {
          node = next_[node];
        }
      }
    }
    for (size_t i = 1; i < n; ++i) {
      size_t parent = 0;
      switch (shape) {
        case CoverShape::kCatalan:
        case CoverShape::kBinaryCatalan:
          continue;
        case CoverShape::kRandomRecursive:
          parent = random() % i;
          break;
        case CoverShape::kStar:
          break;
        case CoverShape::kChain:
          parent = i - 1;
          break;
        case CoverShape::kComb:
          parent = (i - 1) & ~size_t{1};
          break;
        case CoverShape::kBalanced:
          parent = (i - 1) / 2;
          break;
      }
      next_[i] = static_cast<uint32_t>(parent);
      ++degree_[parent];
    }
    for (size_t i = n; i-- > 1;) {
      const auto parent = next_[i];
      next_[i] = first_[parent];
      first_[parent] = static_cast<uint32_t>(i);
    }
  }
  size_t size() const { return first_.size(); }
  uint32_t root() const { return 0; }
  size_t degree(uint32_t node) const { return degree_[node]; }
  uint32_t first_child(uint32_t node) const { return first_[node]; }
  uint32_t next_sibling(uint32_t node) const { return next_[node]; }

 private:
  static constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();
  std::vector<uint32_t> first_, next_, degree_;
};

template <class Tree>
constexpr bool kDecomposed =
    requires(const Tree& tree) { tree.component_count(); };

template <class Tree>
void CoverMemory(benchmark::State& state, const Tree& tree) {
  const auto bytes = tree.memory_usage_bytes();
  state.counters["nodes"] = tree.size();
  state.counters["total_bytes"] = bytes;
  state.counters["bits_per_node"] = 8.0 * bytes / tree.size();
  // One caller-owned cursor, not peak stack usage during decoding/enumeration.
  state.counters["cursor_bytes"] = sizeof(tree.navigate(tree.root()));
  if constexpr (kDecomposed<Tree>) {
    const auto memory = tree.memory_breakdown();
    state.counters["object_bytes"] = memory.object_bytes;
    state.counters["bottom_topology_bytes"] = memory.bottom_topology_bytes;
    state.counters["component_directory_bytes"] =
        memory.component_directory_bytes;
    state.counters["boundary_bytes"] = memory.boundary_bytes;
    state.counters["external_parent_bytes"] = memory.external_parent_bytes;
    state.counters["child_directory_bytes"] = memory.child_directory_bytes;
    state.counters["preorder_bytes"] = memory.preorder_bytes;
    state.counters["macro_bytes"] = memory.macro_bytes;
    state.counters["topology_bits"] = tree.topology_bits();
    state.counters["components"] = tree.component_count();
    state.counters["exits"] = tree.exit_count();
    state.counters["preorder_runs"] = tree.preorder_run_count();
    state.counters["macro_height"] = tree.macro_height();
    state.counters["routing_height"] = tree.routing_height();
    state.counters["cursor_bytes"] = Tree::navigation_scratch_bytes();
    size_t minimum = tree.size(), maximum = 0;
    for (size_t c = 0; c < tree.component_count(); ++c) {
      minimum = std::min<size_t>(minimum, tree.component_size(c));
      maximum = std::max<size_t>(maximum, tree.component_size(c));
    }
    state.counters["component_min_nodes"] = minimum;
    state.counters["component_max_nodes"] = maximum;
    state.counters["component_mean_nodes"] =
        static_cast<double>(tree.size()) / tree.component_count();
  } else {
    state.counters["encoding_bytes"] = tree.encoding_owned_bytes();
  }
}

template <bool Rank, class Tree, class Node>
void CoverObserve(const Tree& tree, Node node) {
  benchmark::DoNotOptimize(node);
  if constexpr (Rank) {
    if constexpr (kDecomposed<Tree>) {
      benchmark::DoNotOptimize(tree.preorder_number(node));
    } else {
      // BP and DFUDS handles already carry preorder rank. LOUDS rank rows are
      // deliberately not registered: its node number is breadth-first order.
      benchmark::DoNotOptimize(node.number);
    }
  }
}

// Inspected walks report entry, children, and nonleaf return events. Plain DFS
// reports entries only, matching the optimized streaming traversal API.
// No original-height-dependent call stack or cache.
template <bool Inspect, class Tree, class Visit, class Child>
void CoverWalk(const Tree& tree, Visit& visit, Child& child) {
  if constexpr (Inspect &&
                requires { tree.traverse_inspected(visit, child); }) {
    return tree.traverse_inspected(visit, child);
  } else {
    auto cursor = tree.navigate(tree.root());
    bool entering = true;
    for (;;) {
      if constexpr (Inspect) {
        visit(cursor.node(), entering);
      } else if (entering) {
        visit(cursor.node(), true);
      }
      if (entering) {
        if constexpr (Inspect) {
          cursor.for_each_child(child);
        }
        if (!cursor.is_leaf()) {
          cursor.reset(cursor.first_child());
          continue;
        }
      }
      if (cursor.is_root()) {
        return;
      }
      if (cursor.is_last_child()) {
        cursor.reset(cursor.parent());
        entering = false;
      } else {
        cursor.reset(cursor.next_sibling());
        entering = true;
      }
    }
  }
}

template <bool Inspect, bool Rank, bool Optimized, class Tree>
void CoverDfs(benchmark::State& state,
              const Tree& tree,
              size_t budget = tree::kDefaultTraversalBytes) {
  auto visit = [&](auto node, bool entry) {
    CoverObserve<Rank>(tree, node);
    if constexpr (Inspect) {
      benchmark::DoNotOptimize(entry);
    }
  };
  auto child = [&](auto node) { CoverObserve<Rank>(tree, node); };
  if constexpr (Optimized) {
    // The optimized baseline's immutable root index is setup, not scratch.
    // Every walk still creates its own decoder buffers inside the timed call.
    auto start = std::chrono::steady_clock::now();
    auto index = MakeWalkIndex(tree, budget);
    pixie::tree::TreeTraversalContext context(tree, budget, &index);
    state.counters["prepare_seconds"] =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    pixie::benchmarks::TreeTraversalWorkloads walk(context);
    for (auto _ : state) {
      if constexpr (Inspect) {
        benchmark::DoNotOptimize(walk.inspect(visit, child));
      } else {
        auto entry = [&](auto node) { CoverObserve<Rank>(tree, node); };
        walk.traverse(entry);
      }
    }
    const auto retained = context.memory_usage_bytes() -
                          context.scratch_bytes() + index.memory_usage_bytes();
    state.counters["traversal_index_bytes"] = retained;
    state.counters["total_bytes"] = tree.memory_usage_bytes() + retained;
    state.counters["bits_per_node"] =
        8.0 * (tree.memory_usage_bytes() + retained) / tree.size();
    state.counters["scratch_bytes"] = context.scratch_bytes();
    state.counters["budget_usage_bytes"] =
        context.budget_usage_bytes() + index.budget_usage_bytes();
    state.counters.erase("cursor_bytes");
  } else {
    if constexpr (Inspect && kDecomposed<Tree>) {
      state.counters["enumeration_cursor_bytes"] =
          Tree::navigation_scratch_bytes();
    }
    for (auto _ : state) {
      CoverWalk<Inspect>(tree, visit, child);
    }
  }
  state.SetItemsProcessed(state.iterations() * tree.size());
  state.counters["visited_nodes_per_walk"] = tree.size();
  if constexpr (Inspect) {
    size_t entered = 0, returned = 0, inspected = 0;
    auto count = [&](auto, bool entry) {
      entered += entry;
      returned += !entry;
    };
    auto child_count = [&](auto) { ++inspected; };
    if constexpr (Optimized) {
      auto index = MakeWalkIndex(tree, budget);
      pixie::tree::TreeTraversalContext context(tree, budget, &index);
      pixie::benchmarks::TreeTraversalWorkloads(context).inspect(count,
                                                                 child_count);
    } else {
      CoverWalk<true>(tree, count, child_count);
    }
    state.counters["visited_nodes_per_walk"] = entered;
    state.counters["return_events_per_walk"] = returned;
    state.counters["inspected_children_per_walk"] = inspected;
    if (entered != tree.size() || inspected + 1 != tree.size()) {
      state.SkipWithError("Incorrect DFS coverage");
    }
  }
}

// Match RandomRootToLeaf<true>: one pooled random value modulo degree per
// nonleaf (sampling child choices, not leaves uniformly). Inspect every child
// in order and retain the selected handle; never fetch it again via
// child(index). Each path owns a fresh cursor. Leaves have no child
// observations.
template <class Tree, class Queries, class Visit, class Inspect>
void CoverInspectedPath(const Tree& tree,
                        const Queries& queries,
                        size_t& position,
                        Visit& visit,
                        Inspect& inspect) {
  auto choose = [&](size_t degree) {
    return queries[position++ % queries.size()] % degree;
  };
  if constexpr (requires { tree.inspect_path(visit, inspect, choose); }) {
    return tree.inspect_path(visit, inspect, choose);
  } else {
    auto cursor = tree.navigate(tree.root());
    for (;;) {
      visit(cursor.node());
      if (cursor.is_leaf()) {
        return;
      }
      const auto index = choose(cursor.degree());
      auto chosen = cursor.node();
      size_t i = 0;
      auto child = [&](auto node) {
        inspect(node);
        if (i++ == index) {
          chosen = node;
        }
      };
      cursor.for_each_child(child);
      cursor.reset(chosen);
    }
  }
}

template <bool Indexed = false, bool FixedBatch = false, class Tree>
void CoverInspectedRootPath(benchmark::State& state,
                            const Tree& tree,
                            size_t batch = 1) {
  std::array<uint64_t, 8192> queries;
  std::mt19937_64 random(state.range(1) ^ 0x9e3779b97f4a7c15ULL);
  for (auto& query : queries) {
    query = random();
  }
  size_t position = 0;
  auto observe = [&](auto node) { CoverObserve<false>(tree, node); };
  auto path = [&](auto& visit, auto& child) {
    if constexpr (Indexed) {
      struct Stream {
        const std::array<uint64_t, 8192>& queries;
        size_t& position;
        auto next() { return queries[position++ % queries.size()]; }
      } stream{queries, position};
      pixie::benchmarks::IndexedRandomRootToLeaf<true>(tree, stream, visit,
                                                       child);
    } else {
      CoverInspectedPath(tree, queries, position, visit, child);
    }
  };
  for (auto _ : state) {
    if constexpr (FixedBatch) {
      position = 0;
    }
    for (size_t i = 0; i < batch; ++i) {
      path(observe, observe);
    }
  }
  // Replay exactly the measured query stream outside timing, without adding
  // per-node accounting to any implementation's timed path.
  position = 0;
  size_t visited = 0, inspected = 0;
  auto visit = [&](auto) { ++visited; };
  auto inspect = [&](auto) { ++inspected; };
  const auto replay =
      FixedBatch ? int64_t(batch) : state.iterations() * int64_t(batch);
  for (int64_t i = 0; i < replay; ++i) {
    path(visit, inspect);
  }
  state.counters["nodes_per_path"] = static_cast<double>(visited) / replay;
  state.counters["inspected_children_per_path"] =
      static_cast<double>(inspected) / replay;
  state.counters["paths_per_iteration"] = batch;
  state.SetItemsProcessed(state.iterations() * batch);
}

}  // namespace pixie::benchmarks
