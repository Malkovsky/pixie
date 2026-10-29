#include <benchmark/benchmark.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/experimental/tree_micro_compact.h>
#include <pixie/tree/experimental/tree_micro_prefetch.h>
#include <pixie/tree/experimental/tree_window_prefetch.h>
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

#ifdef SDSL_SUPPORT
#include <pixie/rmm/implementations.h>
#endif

using pixie::benchmarks::RunDepthFirstBenchmark;
using pixie::benchmarks::RunExpandedTreeBenchmark;
using pixie::benchmarks::TreeSizes;
using pixie::benchmarks::TreeWorkload;

namespace {

using ClockTree = pixie::BPTree<pixie::experimental::RmMBTree<>,
                                pixie::experimental::SecondChanceTreeCache<>>;
using SelectiveUnpackTree =
    pixie::BPTree<pixie::experimental::RmMBTree<>,
                  pixie::experimental::UnpackedTreeCache<false>>;
using ChunkUnpackTree =
    pixie::BPTree<pixie::experimental::RmMBTree<>,
                  pixie::experimental::UnpackedTreeCache<true>>;

using PointerTree = pixie::references::PointerTree;

using NativeBpTree = pixie::BPTree<pixie::RmMTree>;
using CompactBpTree = pixie::BPTree<pixie::experimental::RmMBTree<>>;

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
        words = pixie::benchmarks::GenerateRandomBp(tree_size);
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
void BM_BpComparison(benchmark::State& state) {
  const auto& f = GetFixture<CompactBpTree>(state.range(0));
  pixie::benchmarks::RunTreeComparison<Workload>(state, *f.tree, f.words);
}
BENCHMARK_TEMPLATE(BM_BpComparison, TreeWorkload::kDfs)
    ->Name("bp_best_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_BpComparison, TreeWorkload::kDfsChildren)
    ->Name("bp_best_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_BpComparison, TreeWorkload::kRandomPath)
    ->Name("bp_best_random_path")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_BpComparison, TreeWorkload::kRandomPathChildren)
    ->Name("bp_best_random_path_children")
    ->Apply(TreeSizes);

// Native and compact RmM backends use the same traversal workload.
BENCHMARK_TEMPLATE(BM_TreeDFS, NativeBpTree)
    ->Name("bp_native_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, CompactBpTree)
    ->Name("bp_rmmbtree_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, AdjListTree)
    ->Name("adjlist_preorder_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, NativeBpTree, true)
    ->Name("bp_native_view_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, CompactBpTree, true)
    ->Name("bp_rmmbtree_view_dfs")
    ->Apply(TreeSizes);

#ifdef SDSL_SUPPORT
BENCHMARK_TEMPLATE(BM_TreeDFS, pixie::BPTree<pixie::SdslRmMTree>)
    ->Name("bp_sdsl_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS, pixie::BPTree<pixie::SdslRmMTree>, true)
    ->Name("bp_sdsl_view_dfs")
    ->Apply(TreeSizes);
#endif

BENCHMARK_TEMPLATE(BM_TreeDFS, CompactBpTree, true, true)
    ->Name("bp_rmmbtree_path_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS, CompactBpTree, true, true, true)
    ->Name("bp_rmmbtree_unpacked_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_rmmbtree_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("bp_rmmbtree_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_rmmbtree_random_path_children")
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
                   CompactBpTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_rmmbtree_unpacked_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPath)
    ->Name("bp_rmmbtree_unpacked_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   true,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_rmmbtree_unpacked_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPath,
                   true)
    ->Name("bp_rmmbtree_indexed_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   CompactBpTree,
                   false,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren,
                   true)
    ->Name("bp_rmmbtree_indexed_random_path_children")
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
    ->Name("bp_clock_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_clock_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("bp_clock_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_clock_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("bp_selective_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_selective_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("bp_selective_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SelectiveUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_selective_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfs)
    ->Name("bp_chunk_unpack_dfs")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   true,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_chunk_unpack_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("bp_chunk_unpack_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_TreeDFS,
                   ChunkUnpackTree,
                   true,
                   true,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_chunk_unpack_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

namespace {
template <size_t Bits>
void BM_WindowInspect(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactBpTree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  for (auto _ : state) {
    pixie::experimental::BPWindowInspect<Bits>(*fixture.tree, fixture.words,
                                               observe, observe);
  }
  state.counters["nodes"] = fixture.size;
  state.SetItemsProcessed(state.iterations() * fixture.size);
}
}  // namespace
BENCHMARK_TEMPLATE(BM_WindowInspect, 512)
    ->Name("bp_window_512_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_WindowInspect, 2048)
    ->Name("bp_window_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_WindowInspect, 8192)
    ->Name("bp_window_8192_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_WindowInspect, 32768)
    ->Name("bp_window_32768_dfs_children")
    ->Apply(TreeSizes);

namespace {
template <size_t Bits,
          unsigned Compiled = 0,
          bool Compact = false,
          unsigned VectorPolicy = 0>
void BM_MicroInspect(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactBpTree>(state.range(0));
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  using Traversal = std::conditional_t<
      Compact,
      pixie::experimental::compact_bp::BpMicroTraversal<
          Bits, CompactBpTree, decltype(observe), decltype(observe), Compiled,
          true, VectorPolicy>,
      pixie::experimental::detail::BpMicroTraversal<
          Bits, CompactBpTree, decltype(observe), decltype(observe), Compiled>>;
  auto& perf = pixie::benchmarks::PerfQueryWindow::instance();
  if (perf.active()) {
    Traversal warmup(*fixture.tree, fixture.words, observe, observe);
    warmup.walk(fixture.tree->root());
  }
  perf.start();
  for (auto _ : state) {
    Traversal walk(*fixture.tree, fixture.words, observe, observe);
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
BENCHMARK_TEMPLATE(BM_MicroInspect, 512)
    ->Name("bp_micro_512_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 2048)
    ->Name("bp_micro_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192)
    ->Name("bp_micro_8192_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 512, true)
    ->Name("bp_micro_compiled_512_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, true)
    ->Name("bp_micro_compiled_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, true)
    ->Name("bp_micro_compiled_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
using Selective2048 = pixie::BPTree<
    pixie::experimental::RmMBTree<>,
    pixie::experimental::UnpackedTreeCache<false, 512 * 1024, 2048>>;
using Chunks2048 = pixie::BPTree<
    pixie::experimental::RmMBTree<>,
    pixie::experimental::UnpackedTreeCache<true, 512 * 1024, 2048>>;
template <class Tree, size_t Bits = 2048>
void BM_RetentionInspect(benchmark::State& state) {
  const auto& f = GetFixture<Tree>(state.range(0));
  auto observe = [](auto n) { benchmark::DoNotOptimize(n); };
  for (auto _ : state) {
    pixie::detail::tree_unpack::InspectTraversal<true, Tree, decltype(observe),
                                                 decltype(observe), Bits>
        walk(*f.tree, f.words, observe, observe);
    walk.walk(f.tree->root());
    if (walk.node_count != f.size || walk.child_count != f.size - 1) {
      state.SkipWithError("Bad coverage");
      break;
    }
  }
  state.counters["nodes"] = f.size;
  if constexpr (requires {
                  f.tree->navigate(f.tree->root()).cache_memory_usage_bytes();
                }) {
    state.counters["cache_bytes"] =
        f.tree->navigate(f.tree->root()).cache_memory_usage_bytes();
  }
}
}  // namespace
BENCHMARK_TEMPLATE(BM_RetentionInspect, CompactBpTree, 256)
    ->Name("bp_discard_256_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_RetentionInspect, CompactBpTree, 2048)
    ->Name("bp_discard_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_RetentionInspect, Selective2048)
    ->Name("bp_selective_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_RetentionInspect, Chunks2048)
    ->Name("bp_chunks_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 3)
    ->Name("bp_micro_hybrid_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 3)
    ->Name("bp_micro_hybrid_8192_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 4)
    ->Name("bp_micro_simd_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 4)
    ->Name("bp_micro_simd_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
// Root-connected prefetch plus the microtree shape table must share a 512 KiB
// persistent-data allowance. Query semantics and pool match the existing path
// benchmark; setup stays outside timing, as for the earlier fixed index.
template <bool Inspect>
void BM_BoundedPrefetchPath(benchmark::State& state) {
  const auto& fixture = GetFixture<CompactBpTree>(state.range(0));
  pixie::experimental::TreeNavigationIndex<CompactBpTree, 15600> index(
      *fixture.tree);
  constexpr size_t shape_bytes =
      sizeof(pixie::experimental::detail::bp_micro_tables) +
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
    ->Name("bp_prefetch512_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_BoundedPrefetchPath, true)
    ->Name("bp_prefetch512_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 7)
    ->Name("bp_micro_reservoir_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 7)
    ->Name("bp_micro_reservoir_8192_dfs_children")
    ->Apply(TreeSizes);

BENCHMARK_TEMPLATE(BM_MicroInspect, 2048, 8)
    ->Name("bp_micro_absolute_2048_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 8)
    ->Name("bp_micro_absolute_8192_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 8, true)
    ->Name("bp_micro_compact24_8192_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_MicroInspect, 8192, 8, true, 1)
    ->Name("bp_micro_vector24_8192_dfs_children")
    ->Apply(TreeSizes);

namespace {
using Policy =
    pixie::experimental::SecondChanceTreeCache<512 * 1024, 32, false, true>;
using SnapshotClockTree =
    pixie::BPTree<pixie::experimental::RmMBTree<>, Policy>;
}  // namespace
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfs)
    ->Name("bp_clock32_snapshot_dfs")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kDfsChildren)
    ->Name("bp_clock32_snapshot_dfs_children")
    ->Apply(TreeSizes);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPath)
    ->Name("bp_clock32_snapshot_random_path")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_TreeDFS,
                   SnapshotClockTree,
                   true,
                   false,
                   false,
                   TreeWorkload::kRandomPathChildren)
    ->Name("bp_clock32_snapshot_random_path_children")
    ->Apply(TreeSizes)
    ->Unit(benchmark::kNanosecond);
