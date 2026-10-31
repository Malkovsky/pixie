#pragma once

#include <pixie/utils.h>
#include <references/pointer_tree.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <type_traits>
#include <vector>

#include "tree_benchmark_helpers.h"
#include "tree_comparison.h"
#include "tree_input.h"

namespace pixie::benchmarks {
enum class TreeLayout { kBp, kDfuds, kLouds };

template <class Tree>
inline constexpr bool kReferenceTree =
    std::is_same_v<Tree, AdjListTree> ||
    std::is_same_v<Tree, references::PointerTree>;

// One immutable fixture per type/layout and size, shared by all workload rows.
// Replacement releases the previous tree before generating the next size.
template <class Tree, TreeLayout Layout>
const auto& GetTreeFixture(std::size_t size) {
  struct Fixture {
    std::size_t size;
    std::unique_ptr<Tree> tree;
    explicit Fixture(std::size_t n) : size(n) {
      if constexpr (kReferenceTree<Tree>) {
        std::mt19937_64 rng(42);
        const auto adjacency = Layout == TreeLayout::kLouds
                                   ? bfs_order(n, generate_random_tree(n, rng))
                                   : dfs_order(n, generate_random_tree(n, rng));
        tree = std::make_unique<Tree>(adjacency);
      } else {
        auto words = GenerateRandomBp(n);
        if constexpr (Layout == TreeLayout::kDfuds) {
          words = DfudsFromBp(words, n);
        } else if constexpr (Layout == TreeLayout::kLouds) {
          words = LoudsFromBp(words, n);
        }
        tree = std::make_unique<Tree>(words, n);
      }
    }
  };
  static std::unique_ptr<Fixture> cached;
  if (!cached || cached->size != size) {
    cached.reset();
    cached = std::make_unique<Fixture>(size);
  }
  return *cached;
}

template <class Tree, TreeLayout Layout, TreeWorkload Workload>
void RunDirectTreeBenchmark(benchmark::State& state) {
  const auto& fixture = GetTreeFixture<Tree, Layout>(state.range(0));
  const auto& tree = *fixture.tree;
  if constexpr (std::is_same_v<Tree, references::PointerTree>) {
    state.counters["reference_bytes"] = tree.memory_usage_bytes();
    state.counters["reference_bytes_per_node"] =
        double(tree.memory_usage_bytes()) / fixture.size;
  }
  if constexpr (Workload == TreeWorkload::kDfs) {
    RunDepthFirstBenchmark<!kReferenceTree<Tree>>(state, tree, fixture.size);
  } else {
    RunExpandedTreeBenchmark<Workload>(state, tree, fixture.size);
  }
  if constexpr (!kReferenceTree<Tree>) {
    state.counters["input_bytes"] = tree.encoding_owned_bytes();
  }
}

template <class Tree, TreeLayout Layout, TreeWorkload Workload>
void RunBestTreeBenchmark(benchmark::State& state) {
  const auto& fixture = GetTreeFixture<Tree, Layout>(state.range(0));
  RunTreeComparison<Workload>(state, *fixture.tree);
}
}  // namespace pixie::benchmarks
