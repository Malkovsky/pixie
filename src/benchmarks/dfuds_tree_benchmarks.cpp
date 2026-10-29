#include <benchmark/benchmark.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/experimental/tree_dfuds_micro.h>
#include <pixie/tree/experimental/tree_prepared_events.h>
#include <pixie/tree/experimental/tree_prepared_root.h>
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
    pixie::DFUDSTree<pixie::experimental::RmMBTree<>,
                     pixie::experimental::SecondChanceTreeCache<>>;
using SelectiveUnpackTree =
    pixie::DFUDSTree<pixie::experimental::RmMBTree<>,
                     pixie::experimental::UnpackedTreeCache<false>>;
using ChunkUnpackTree =
    pixie::DFUDSTree<pixie::experimental::RmMBTree<>,
                     pixie::experimental::UnpackedTreeCache<true>>;

using PointerTree = pixie::references::PointerTree;

using NativeDfudsTree = pixie::DFUDSTree<pixie::RmMTree>;
using CompactDfudsTree = pixie::DFUDSTree<pixie::experimental::RmMBTree<>>;

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
      // Same topology as LOUDS, numbered in the encoding's natural DFS order.
      if constexpr (std::is_same_v<Tree, AdjListTree> ||
                    std::is_same_v<Tree, PointerTree>) {
        std::mt19937_64 rng(42);
        const auto adjacency =
            dfs_order(tree_size, generate_random_tree(tree_size, rng));
        tree = std::make_unique<Tree>(adjacency);
      } else {
        words = pixie::benchmarks::DfudsFromBp(
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
  if constexpr (!std::is_same_v<Tree, AdjListTree> &&
                !std::is_same_v<Tree, PointerTree>) {
    state.counters["input_bytes"] = words.capacity() * sizeof(std::uint64_t);
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
void BM_DfudsComparison(benchmark::State& state) {
  const auto& f = GetFixture<CompactDfudsTree>(state.range(0));
  pixie::benchmarks::RunTreeComparison<Workload>(state, *f.tree, f.words);
}
BENCHMARK_TEMPLATE(BM_DfudsComparison, TreeWorkload::kDfs)
    ->Name("dfuds_best_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_DfudsComparison, TreeWorkload::kDfsChildren)
    ->Name("dfuds_best_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_DfudsComparison, TreeWorkload::kRandomPath)
    ->Name("dfuds_best_random_path")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_DfudsComparison, TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_best_random_path_children")
    ->Apply(TreeSizes);

// SdslRmMTree does not implement the select0 operation required by DFUDS.
BENCHMARK_TEMPLATE(BM_TreeDFS, NativeDfudsTree)
    ->Name("dfuds_native_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, CompactDfudsTree)
    ->Name("dfuds_rmmbtree_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, AdjListTree)
    ->Name("adjlist_preorder_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, NativeDfudsTree, true)
    ->Name("dfuds_native_view_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, CompactDfudsTree, true)
    ->Name("dfuds_rmmbtree_view_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, CompactDfudsTree, true, true)
    ->Name("dfuds_rmmbtree_path_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, CompactDfudsTree, true, true, true)
    ->Name("dfuds_rmmbtree_unpacked_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_rmmbtree_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_rmmbtree_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_rmmbtree_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("adjlist_preorder_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("adjlist_preorder_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   AdjListTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("adjlist_preorder_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_rmmbtree_unpacked_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_rmmbtree_unpacked_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_rmmbtree_unpacked_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath,
                   true)
    ->Name("dfuds_rmmbtree_indexed_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactDfudsTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren,
                   true)
    ->Name("dfuds_rmmbtree_indexed_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS, PointerTree)
    ->Name("pointer_preorder_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("pointer_preorder_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("pointer_preorder_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   PointerTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("pointer_preorder_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfs)
    ->Name("dfuds_clock_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_clock_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_clock_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_clock_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("dfuds_selective_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_selective_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_selective_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_selective_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("dfuds_chunk_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_chunk_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_chunk_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_chunk_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

namespace {
template <size_t Bits, unsigned Compiled = 0>
void BM_MicroInspect(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactDfudsTree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  if (perf.active()) {
    pixie::experimental::dfuds_detail::DfudsMicroTraversal<
        Bits, CompactDfudsTree, decltype(observe), decltype(observe), Compiled>
        warmup(*fixture.tree, fixture.words, observe, observe);
    warmup.walk(fixture.tree->root());
  }
  perf.start();
  for (auto _ : state) {
    pixie::experimental::dfuds_detail::DfudsMicroTraversal<
        Bits, CompactDfudsTree, decltype(observe), decltype(observe), Compiled>
        walk(*fixture.tree, fixture.words, observe, observe);
    walk.walk(fixture.tree->root());
    if (walk.nodes != fixture.size || walk.children != fixture.size - 1) {
      state.SkipWithError("Incorrect microtree coverage");
      break;
    }
  }
  perf.stop(state.iterations() * fixture.size);
  if (perf.active()) {
    state.counters["perf_nodes"] = perf.queries();
    state.counters["perf_phases"] = perf.phases();
  }
  state.counters["nodes"] = fixture.size;
  state.SetItemsProcessed(state.iterations() * fixture.size);
}
}  // namespace
BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 8)
    ->Name("dfuds_micro_absolute_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 8)
    ->Name("dfuds_micro_absolute_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
// Root-connected prefetch plus the microtree shape table must share a 512 KiB
// persistent-data allowance. Query semantics and pool match the existing path
// benchmark; setup stays outside timing, as for the earlier fixed index.
template <bool Inspect>
void BM_BoundedPrefetchPath(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactDfudsTree>(state.range(0));
  pixie::experimental::TreeNavigationIndex<CompactDfudsTree, 15600> index(
      *fixture.tree);
  constexpr size_t shape_bytes =
      sizeof(pixie::experimental::dfuds_detail::dfuds_micro_tables) +
      197 * sizeof(void (*)());
  const size_t bytes = index.memory_usage_bytes() + shape_bytes;
  if (bytes > 512 * 1024) {
    state.SkipWithError("Prefetch budget exceeded");
    return;
  }
  pixie::benchmarks::ShiftedQueryPool queries(
      0, std::numeric_limits<size_t>::max(), 43);
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  size_t nodes = 0, children = 0;
  for (auto _ : state) {
    auto work = pixie::benchmarks::IndexedRandomRootToLeaf<Inspect>(
        index, queries, observe, observe);
    nodes += work.nodes;
    children += work.children;
    benchmark::DoNotOptimize(work);
  }
  state.counters["nodes"] = fixture.size;
  state.counters["cache_bytes"] = bytes;
  state.counters["visited_nodes_per_iteration"] =
      double(nodes) / state.iterations();
  state.counters["inspected_children_per_iteration"] =
      double(children) / state.iterations();
  state.counters["query_pool_bytes"] = 65536;
}
}  // namespace
BENCHMARK_TEMPLATE(BM_BoundedPrefetchPath, false)
    ->Name("dfuds_prefetch512_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_BoundedPrefetchPath, true)
    ->Name("dfuds_prefetch512_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

namespace {
template <size_t Bits>
void BM_PlainChunkInspect(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactDfudsTree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  for (auto _ : state) {
    pixie::detail::tree_unpack::InspectTraversal<
        false, CompactDfudsTree, decltype(observe), decltype(observe), Bits>
        walk(*fixture.tree, fixture.words, observe, observe);
    walk.walk(fixture.tree->root());
    if (walk.node_count != fixture.size ||
        walk.child_count != fixture.size - 1) {
      state.SkipWithError("Incorrect plain chunk coverage");
      break;
    }
  }
  state.counters["nodes"] = fixture.size;
  state.SetItemsProcessed(state.iterations() * fixture.size);
}
}  // namespace
BENCHMARK_TEMPLATE(BM_PlainChunkInspect, 8192)
    ->Name("dfuds_plain_chunk_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
template <bool Inspect>
void BM_PreparedRoot(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactDfudsTree>(state.range(0));
  using Decoder =
      pixie::experimental::dfuds_detail::DfudsMicroSubtree<8192, 1, false,
                                                           true>;
  pixie::experimental::PreparedRootSubtree<Decoder, 8192> cache(
      fixture.words, 2 * fixture.size - 1, fixture.tree->root());
  pixie::experimental::TreeNavigationIndex<CompactDfudsTree, 15500> index(
      *fixture.tree);
  const size_t bytes =
      cache.memory_usage_bytes() + index.memory_usage_bytes() +
      sizeof(pixie::experimental::dfuds_detail::dfuds_micro_tables) +
      2 * 197 * sizeof(void (*)());
  if (bytes > 512 * 1024) {
    state.SkipWithError("Prepared-root budget exceeded");
    return;
  }
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto inspect = [&](auto node) {
    if constexpr (Inspect) {
      observe(node);
    }
  };
  for (auto _ : state) {
    if (!cache.try_walk(observe, inspect)) {
      if constexpr (Inspect) {
        pixie::experimental::dfuds_detail::DfudsMicroTraversal<
            8192, CompactDfudsTree, decltype(observe), decltype(inspect), 8>
            walk(*fixture.tree, fixture.words, observe, inspect);
        walk.walk(fixture.tree->root());
      } else {
        pixie::benchmarks::UnpackedDepthFirst(*fixture.tree, fixture.words,
                                              observe);
      }
    }
  }
  state.counters["nodes"] = fixture.size;
  state.counters["cache_bytes"] = bytes;
  state.counters["prepared_nodes"] = cache.size();
  state.SetItemsProcessed(state.iterations() * fixture.size);
}
}  // namespace
BENCHMARK_TEMPLATE(BM_PreparedRoot, false)
    ->Name("dfuds_prepared_root_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_PreparedRoot, true)
    ->Name("dfuds_prepared_root_dfs_children")
    ->Apply(TreeSizes);

namespace {
template <bool Inspect, bool Compiled = true>
void BM_PreparedEvents(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactDfudsTree>(state.range(0));
  pixie::experimental::PreparedTreeEvents<> cache(*fixture.tree);
  pixie::experimental::TreeNavigationIndex<CompactDfudsTree, 15400> index(
      *fixture.tree);
  const size_t bytes =
      cache.memory_usage_bytes() + index.memory_usage_bytes() +
      sizeof(pixie::experimental::dfuds_detail::dfuds_micro_tables) +
      (197 + 2 * 256) * sizeof(void (*)());
  if (bytes > 512 * 1024) {
    state.SkipWithError("Prepared-root budget exceeded");
    return;
  }
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto inspect = [&](auto node) {
    if constexpr (Inspect) {
      observe(node);
    }
  };
  for (auto _ : state) {
    if (!cache.template try_walk<Compiled, Inspect>(observe, inspect)) {
      if constexpr (Inspect) {
        pixie::experimental::dfuds_detail::DfudsMicroTraversal<
            8192, CompactDfudsTree, decltype(observe), decltype(inspect), 8>
            walk(*fixture.tree, fixture.words, observe, inspect);
        walk.walk(fixture.tree->root());
      } else {
        pixie::benchmarks::UnpackedDepthFirst(*fixture.tree, fixture.words,
                                              observe);
      }
    }
  }
  state.counters["nodes"] = fixture.size;
  state.counters["cache_bytes"] = bytes;
  state.counters["prepared_nodes"] = cache.size();
  state.SetItemsProcessed(state.iterations() * fixture.size);
}
}  // namespace
BENCHMARK_TEMPLATE(BM_PreparedEvents, false)
    ->Name("dfuds_prepared_events_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_PreparedEvents, true)
    ->Name("dfuds_prepared_events_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_PreparedEvents, true, false)
    ->Name("dfuds_prepared_events_interpreted_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_PreparedEvents, false, false)
    ->Name("dfuds_prepared_events_interpreted_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 9)
    ->Name("dfuds_micro_excess16_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 9)
    ->Name("dfuds_micro_excess16_8192_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 10)
    ->Name("dfuds_micro_record24_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
using Policy =
    pixie::experimental::SecondChanceTreeCache<512 * 1024, 32, false, true>;
using SnapshotClockTree =
    pixie::DFUDSTree<pixie::experimental::RmMBTree<>, Policy>;
}  // namespace
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfs)
    ->Name("dfuds_clock32_snapshot_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("dfuds_clock32_snapshot_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("dfuds_clock32_snapshot_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("dfuds_clock32_snapshot_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
