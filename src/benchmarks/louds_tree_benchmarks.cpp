#include <benchmark/benchmark.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/tree/experimental/tree_louds_iterative.h>
#include <pixie/tree/experimental/tree_louds_prefetch.h>
#include <pixie/tree/implementations.h>
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

using pixie::benchmarks::RunDepthFirstBenchmark;
using pixie::benchmarks::RunExpandedTreeBenchmark;
using pixie::benchmarks::TreeSizes;
using pixie::benchmarks::TreeWorkload;

namespace {

using ClockTree =
    pixie::LoudsTree<pixie::RankSelectSupportCompressed<>,
                     pixie::experimental::SecondChanceTreeCache<>>;
using SelectiveUnpackTree =
    pixie::LoudsTree<pixie::RankSelectSupportCompressed<>,
                     pixie::experimental::UnpackedTreeCache<false>>;
using ChunkUnpackTree =
    pixie::LoudsTree<pixie::RankSelectSupportCompressed<>,
                     pixie::experimental::UnpackedTreeCache<true>>;

using PointerTree = pixie::references::PointerTree;

// Retain one immutable fixture per concrete type and size across repetitions
// and workload rows. Generation, construction, and replacement remain outside
// timing; this avoids repeating the expensive 10-million-node generator.
template <class Tree>
const auto& GetFixture(std::size_t tree_size) {
  struct Fixture {
    std::size_t size;
    std::vector<std::uint64_t> words;
    std::unique_ptr<Tree> tree;
    explicit Fixture(std::size_t tree_size) : size(tree_size) {
      // Every backend and repetition receives the same BFS-numbered tree.
      if constexpr (std::is_same_v<Tree, AdjListTree> ||
                    std::is_same_v<Tree, PointerTree>) {
        std::mt19937_64 rng(42);
        const auto adjacency =
            bfs_order(tree_size, generate_random_tree(tree_size, rng));
        tree = std::make_unique<Tree>(adjacency);
      } else {
        words = pixie::benchmarks::LoudsFromBp(
            pixie::benchmarks::GenerateRandomBp(tree_size), tree_size);
        tree = std::make_unique<Tree>(words, tree_size,
                                      pixie::benchmarks::TreeCacheBudget());
      }
    }
  };
  static std::unique_ptr<Fixture> cached;
  if (!cached || cached->size != tree_size) {
    cached.reset();
    cached = std::make_unique<Fixture>(tree_size);
  }
  return *cached;
}

template <class Tree,
          bool UseView = false,
          bool UsePath = false,
          bool Unpack = false,
          TreeWorkload Workload = TreeWorkload::kDfs,
          bool Indexed = false>
void BM_TreeDFS(benchmark::State& state) {
  const std::size_t tree_size = state.range(0);
  const auto& fixture = GetFixture<Tree>(tree_size);
  const auto& tree = *fixture.tree;
  const auto& words = fixture.words;

  if constexpr (std::is_same_v<Tree, PointerTree>) {
    state.counters["reference_bytes"] = tree.memory_usage_bytes();
    state.counters["reference_bytes_per_node"] =
        double(tree.memory_usage_bytes()) / tree_size;
  }

  // Construction and destruction of the tree and its source are outside timing.
  if constexpr (Workload == TreeWorkload::kDfs) {
    RunDepthFirstBenchmark<UseView, UsePath, Unpack>(state, tree, tree_size,
                                                     words);
  } else {
    RunExpandedTreeBenchmark<Workload, Unpack, Indexed>(state, tree, tree_size,
                                                        words);
  }
  if constexpr (requires {
                  tree.navigate(tree.root()).cache_memory_usage_bytes();
                }) {
    state.counters["cache_bytes"] =
        tree.navigate(tree.root()).cache_memory_usage_bytes();
  }
}

}  // namespace

template <TreeWorkload Workload>
void BM_LoudsComparison(benchmark::State& state) {
  const auto& f =
      GetFixture<pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>>(
          state.range(0));
  pixie::benchmarks::RunTreeComparison<Workload>(state, *f.tree, f.words);
}
BENCHMARK_TEMPLATE(BM_LoudsComparison, TreeWorkload::kDfs)
    ->Name("louds_best_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_LoudsComparison, TreeWorkload::kDfsChildren)
    ->Name("louds_best_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_LoudsComparison, TreeWorkload::kRandomPath)
    ->Name("louds_best_random_path")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_LoudsComparison, TreeWorkload::kRandomPathChildren)
    ->Name("louds_best_random_path_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, pixie::LoudsTree<pixie::RankSelectSupport<>>)
    ->Name("louds_baseline_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>)
    ->Name("louds_compressed_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, AdjListTree)
    ->Name("adjlist_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupport<>>,
                   true)
    ->Name("louds_baseline_view_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   true)
    ->Name("louds_compressed_view_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupport<>>,
                   true,
                   true)
    ->Name("louds_baseline_path_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   true,
                   true)
    ->Name("louds_compressed_path_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   true,
                   true,
                   true)
    ->Name("louds_compressed_unpacked_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_compressed_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("louds_compressed_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_compressed_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("adjlist_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("adjlist_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("adjlist_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_compressed_unpacked_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPath)
    ->Name("louds_compressed_unpacked_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_compressed_unpacked_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath,
                   true)
    ->Name("louds_compressed_indexed_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren,
                   true)
    ->Name("louds_compressed_indexed_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS, PointerTree)
    ->Name("pointer_bfs_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("pointer_bfs_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("pointer_bfs_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("pointer_bfs_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfs)
    ->Name("louds_clock_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_clock_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("louds_clock_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_clock_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("louds_selective_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_selective_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("louds_selective_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_selective_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("louds_chunk_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_chunk_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("louds_chunk_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_chunk_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

namespace {
template <size_t Capacity, bool Decoded = false>
void BM_LoudsPrefetchInspect(benchmark::State& state) {
  using Tree = pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  const auto& fixture = GetFixture<Tree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  if (perf.active()) {
    pixie::experimental::LoudsPrefetchTraversal<
        Tree, decltype(observe), decltype(observe), Capacity, Decoded>
        warmup(*fixture.tree, fixture.words, observe, observe);
    warmup.walk(fixture.tree->root(), 0);
  }
  perf.start();
  for (auto _ : state) {
    pixie::experimental::LoudsPrefetchTraversal<
        Tree, decltype(observe), decltype(observe), Capacity, Decoded>
        walk(*fixture.tree, fixture.words, observe, observe);
    walk.walk(fixture.tree->root(), 0);
    if (walk.node_count != fixture.size ||
        walk.child_count != fixture.size - 1) {
      state.SkipWithError("Incorrect prefetched LOUDS coverage");
      break;
    }
  }
  perf.stop(state.iterations() * fixture.size);
  state.counters["nodes"] = fixture.size;
  state.SetItemsProcessed(state.iterations() * fixture.size);
  if (perf.active()) {
    state.counters["perf_nodes"] = perf.queries();
    state.counters["perf_phases"] = perf.phases();
  }
}
}  // namespace
BENCHMARK_TEMPLATE(BM_LoudsPrefetchInspect, 8)
    ->Name("louds_prefetch8_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_LoudsPrefetchInspect, 32)
    ->Name("louds_prefetch32_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_LoudsPrefetchInspect, 32, true)
    ->Name("louds_prefetch_decoded32_dfs_children")
    ->Apply(TreeSizes);

namespace {
template <bool Decoded>
void BM_LoudsIterativeInspect(benchmark::State& state) {
  using Tree = pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  const auto& fixture = GetFixture<Tree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  if (perf.active()) {
    pixie::experimental::LoudsIterativeTraversal<Tree, decltype(observe),
                                                 decltype(observe), Decoded>
        warmup(*fixture.tree, fixture.words, observe, observe);
    warmup.walk(fixture.tree->root(), 0);
  }
  perf.start();
  for (auto _ : state) {
    pixie::experimental::LoudsIterativeTraversal<Tree, decltype(observe),
                                                 decltype(observe), Decoded>
        walk(*fixture.tree, fixture.words, observe, observe);
    walk.walk(fixture.tree->root(), 0);
    if (walk.node_count != fixture.size ||
        walk.child_count != fixture.size - 1) {
      state.SkipWithError("Incorrect prefetched LOUDS coverage");
      break;
    }
  }
  perf.stop(state.iterations() * fixture.size);
  state.counters["nodes"] = fixture.size;
  state.SetItemsProcessed(state.iterations() * fixture.size);
  if (perf.active()) {
    state.counters["perf_nodes"] = perf.queries();
    state.counters["perf_phases"] = perf.phases();
  }
}
}  // namespace

BENCHMARK_TEMPLATE(BM_LoudsIterativeInspect, true)
    ->Name("louds_iterative_decoded_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_LoudsIterativeInspect, false)
    ->Name("louds_iterative_scalar_dfs_children")
    ->Apply(TreeSizes);

namespace {
using Policy =
    pixie::experimental::SecondChanceTreeCache<512 * 1024, 32, false, true>;
using SnapshotClockTree =
    pixie::LoudsTree<pixie::RankSelectSupportCompressed<>, Policy>;
}  // namespace
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfs)
    ->Name("louds_clock32_snapshot_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("louds_clock32_snapshot_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("louds_clock32_snapshot_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("louds_clock32_snapshot_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
