#include <benchmark/benchmark.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/detail/forest256_cache.h>
#include <pixie/tree/experimental/decomposed_tree.h>
#include <pixie/tree/experimental/layered_dfuds.h>
#include <pixie/tree/experimental/simd_micro_traversal.h>
#include <pixie/tree/experimental/tree_prepared_events.h>
#include <pixie/tree/experimental/tree_prepared_root.h>
#include <pixie/tree/forest_cover.h>
#include <pixie/tree/traversal_context.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "forest256_twin.h"
#include "tree_conversion_benchmarks.h"
#include "tree_dfs/optimized.h"
#include "tree_input.h"
#include "tree_region_benchmarks.h"
#include "tree_shape_workloads.h"
#include "tree_traversal_workloads.h"

namespace {
template <bool Dfuds, bool Events>
void Prepared(benchmark::State& state) {
  using Tree =
      std::conditional_t<Dfuds, pixie::tree::DFUDSTree<pixie::RmMBTree<>>,
                         pixie::tree::BPTree<pixie::RmMBTree<>>>;
  using Decoder = std::conditional_t<
      Dfuds, pixie::tree::detail::dfuds_micro::DfudsMicroSubtree<8192>,
      pixie::tree::detail::bp_micro::BpMicroSubtree<8192>>;
  const size_t n = state.range(0);
  auto words = pixie::benchmarks::GenerateRandomBp(n);
  if constexpr (Dfuds) {
    words = pixie::benchmarks::ConvertBpEncoding<
        pixie::tree::detail::TreeEncodingKind::kDFUDS>(words, n);
  }
  Tree tree(words, n);
  pixie::benchmarks::TreeTraversalWorkloads fallback(tree);
  auto prepared = [&] {
    if constexpr (Events) {
      return pixie::tree::experimental::PreparedTreeEvents<>(tree);
    } else {
      return pixie::tree::experimental::PreparedRootSubtree<Decoder, 8192>(
          words, 2 * n - (Dfuds ? 1 : 0), tree.root());
    }
  }();
  auto visit = [](auto node, bool entry) {
    benchmark::DoNotOptimize(node);
    benchmark::DoNotOptimize(entry);
  };
  auto inspect = [](auto node) { benchmark::DoNotOptimize(node); };
  for (auto _ : state) {
    if (!prepared.try_walk(visit, inspect)) {
      fallback.inspect(visit, inspect);
    }
  }
  state.counters["prepared_nodes"] = prepared.size();
  state.counters["prepared_bytes"] = prepared.memory_usage_bytes();
  state.SetItemsProcessed(state.iterations() * n);
}

enum class ExpShape { kRandom, kStar, kClusteredTwigs };

std::vector<uint64_t> GenerateShapeBp(size_t n, ExpShape shape) {
  if (shape == ExpShape::kRandom) {
    return pixie::benchmarks::GenerateRandomBp(n);
  }
  std::vector<uint64_t> words((2 * n + 63) / 64, 0);
  size_t pos = 0;
  auto emit1 = [&] {
    words[pos >> 6] |= 1ULL << (pos & 63);
    ++pos;
  };
  auto emit0 = [&] { ++pos; };

  if (shape == ExpShape::kStar) {
    emit1();
    for (size_t i = 1; i < n; ++i) {
      emit1();
      emit0();
    }
    emit0();
    return words;
  }

  // kClusteredTwigs: root -> hubs of ~20 mixed 1..3-node children (~40-80
  // bits/hub)
  emit1();
  size_t rem = n - 1;
  size_t seed = 1;
  while (rem > 0) {
    emit1();
    --rem;
    const size_t hub_budget = std::min<size_t>(rem, 36);
    size_t used = 0;
    while (used < hub_budget) {
      const size_t rem_hub = hub_budget - used;
      const size_t sz = rem_hub >= 3 ? (1 + (seed % 3))
                                     : (rem_hub >= 2 ? (1 + (seed & 1)) : 1);
      seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
      if (sz == 1) {
        emit1();
        emit0();
      } else if (sz == 2) {
        emit1();
        emit1();
        emit0();
        emit0();
      } else {
        emit1();
        emit1();
        emit0();
        emit1();
        emit0();
        emit0();
      }
      used += sz;
    }
    rem -= used;
    emit0();
  }
  emit0();
  return words;
}

template <bool Dfuds, bool UseSimd, size_t Bits, ExpShape Shape>
void MicroDfsChildren(benchmark::State& state) {
  using Tree =
      std::conditional_t<Dfuds, pixie::tree::DFUDSTree<pixie::RmMBTree<>>,
                         pixie::tree::BPTree<pixie::RmMBTree<>>>;
  const size_t n = state.range(0);
  auto words = GenerateShapeBp(n, Shape);
  if constexpr (Dfuds) {
    words = pixie::benchmarks::ConvertBpEncoding<
        pixie::tree::detail::TreeEncodingKind::kDFUDS>(words, n);
  }
  Tree tree(words, n);
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };

  constexpr size_t kMaxDepth = 128;
  size_t total_nodes = 0;
  for (auto _ : state) {
    if constexpr (!UseSimd) {
      using RefWalk = std::conditional_t<
          Dfuds,
          pixie::tree::detail::dfuds_micro::DfudsMicroTraversal<
              Bits, Tree, decltype(observe), decltype(observe)>,
          pixie::tree::detail::bp_micro::BpMicroTraversal<
              Bits, Tree, decltype(observe), decltype(observe)>>;
      RefWalk walk(tree, words, observe, observe, kMaxDepth);
      walk.walk(tree.root());
      if (walk.nodes != n || walk.children != n - 1) {
        state.SkipWithError("Coverage mismatch");
        break;
      }
      total_nodes += walk.nodes;
    } else {
      using SimdWalk = std::conditional_t<
          Dfuds,
          pixie::tree::experimental::DfudsSimdMicroTraversal<
              Bits, Tree, decltype(observe), decltype(observe)>,
          pixie::tree::experimental::BpSimdMicroTraversal<
              Bits, Tree, decltype(observe), decltype(observe)>>;
      SimdWalk walk(tree, words, observe, observe, kMaxDepth);
      walk.walk(tree.root());
      if (walk.nodes != n || walk.children != n - 1) {
        state.SkipWithError("Coverage mismatch");
        break;
      }
      total_nodes += walk.nodes;
    }
  }
  state.SetItemsProcessed(total_nodes);
  state.counters["nodes"] = n;
  if constexpr (!UseSimd) {
    using RefWalk = std::conditional_t<
        Dfuds,
        pixie::tree::detail::dfuds_micro::DfudsMicroTraversal<
            Bits, Tree, decltype(observe), decltype(observe)>,
        pixie::tree::detail::bp_micro::BpMicroTraversal<
            Bits, Tree, decltype(observe), decltype(observe)>>;
    state.counters["scratch_bytes"] = RefWalk::scratch_bytes(kMaxDepth);
  } else {
    using SimdWalk = std::conditional_t<
        Dfuds,
        pixie::tree::experimental::DfudsSimdMicroTraversal<
            Bits, Tree, decltype(observe), decltype(observe)>,
        pixie::tree::experimental::BpSimdMicroTraversal<
            Bits, Tree, decltype(observe), decltype(observe)>>;
    state.counters["scratch_bytes"] = SimdWalk::scratch_bytes(kMaxDepth);
  }
}

void SimdTreeSizes(benchmark::internal::Benchmark* b) {
  b->Arg(4096)->Arg(65536)->Arg(1048576)->Arg(16777216);
}
void ShapeSizes(benchmark::internal::Benchmark* b) {
  b->Arg(128)->Arg(1024)->Arg(65536);
}

using pixie::benchmarks::CoverDfs;
using pixie::benchmarks::CoverInspectedRootPath;
using pixie::benchmarks::CoverMemory;
using pixie::benchmarks::CoverObserve;
using pixie::benchmarks::CoverShape;
using pixie::benchmarks::CoverSource;
using pixie::benchmarks::kDecomposed;

enum class CoverWork {
  kBuild,
  kDfs,
  kDfsRank,
  kInspect,
  kInspectRank,
  kLookup,
  kCursorRank,
  kRootPath,
  kFetch,
  kLocalEdge,
  kBoundaryEdge,
  kInspectedRootPath
};

template <class Tree, class... Options>
Tree BuildCoverTree(benchmark::State& state,
                    CoverShape shape,
                    Options... options) {
  const CoverSource source(state.range(0), shape, state.range(1));
  struct Timer {
    benchmark::State& state;
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    ~Timer() {
      state.counters["build_seconds"] =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        start)
              .count();
    }
  } timer{state};
  return Tree(source, options...);
}

template <class Tree, bool Optimized = false>
void CoverBenchmark(benchmark::State& state, CoverShape shape, CoverWork work) {
  const size_t n = state.range(0);
  const uint64_t seed = state.range(1);
  if (work == CoverWork::kBuild) {
    const CoverSource source(n, shape, seed);
    for (auto _ : state) {
      Tree built(source);
      benchmark::DoNotOptimize(built);
    }
    const Tree measured(source);
    CoverMemory(state, measured);
    // Construction peak is not available from retained-memory accounting.
    state.counters.erase("cursor_bytes");
    state.SetItemsProcessed(state.iterations() * n);
    return;
  }
  const Tree tree = BuildCoverTree<Tree>(state, shape);
  CoverMemory(state, tree);
  if (work == CoverWork::kInspectedRootPath) {
    if constexpr (Optimized) {
      const auto start = std::chrono::steady_clock::now();
      const auto budget = pixie::tree::kDefaultTraversalBytes;
      pixie::benchmarks::TreeTraversalWorkloads walk(tree, budget);
      pixie::tree::RootNavigationIndex index(tree,
                                             budget - walk.reserved_bytes());
      pixie::tree::TreeTraversalContext context(tree, 0, &index);
      state.counters["prepare_seconds"] =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        start)
              .count();
      const auto bytes =
          walk.fixture_metadata_bytes() + index.memory_usage_bytes();
      state.counters["traversal_index_bytes"] = bytes;
      state.counters["scratch_reserved_bytes"] = walk.scratch_bytes();
      state.counters["cursor_bytes"] = sizeof(context.cursor());
      state.counters["total_bytes"] = tree.memory_usage_bytes() + bytes;
      state.counters["bits_per_node"] =
          8.0 * (tree.memory_usage_bytes() + bytes) / n;
      return CoverInspectedRootPath<true>(state, context);
    }
    return CoverInspectedRootPath(state, tree);
  }
  switch (work) {
    case CoverWork::kDfs:
      return CoverDfs<false, false, Optimized>(state, tree);
    case CoverWork::kDfsRank:
      return CoverDfs<false, true, Optimized>(state, tree);
    case CoverWork::kInspect:
      return CoverDfs<true, false, Optimized>(state, tree);
    case CoverWork::kInspectRank:
      return CoverDfs<true, true, Optimized>(state, tree);
    default:
      break;
  }

  constexpr size_t kQueries = 8192;
  std::array<uint64_t, kQueries> queries;
  std::mt19937_64 random(seed ^ 0x9e3779b97f4a7c15ULL);
  for (auto& query : queries) {
    query = random();
  }
  if (work == CoverWork::kRootPath) {
    size_t position = 0, visited = 0;
    for (auto _ : state) {
      auto cursor = tree.navigate(tree.root());
      for (;;) {
        CoverObserve<false>(tree, cursor.node());
        ++visited;
        if (cursor.is_leaf()) {
          break;
        }
        const auto index = queries[position++ % kQueries] % cursor.degree();
        cursor.reset(cursor.child(index));
      }
    }
    state.counters["nodes_per_path"] =
        static_cast<double>(visited) / state.iterations();
    state.SetItemsProcessed(visited);
    return;
  }

  if constexpr (kDecomposed<Tree>) {
    if (work == CoverWork::kLookup) {
      for (auto& query : queries) {
        query %= n;
      }
      for (auto _ : state) {
        for (auto query : queries) {
          benchmark::DoNotOptimize(tree.find_by_preorder(query));
        }
      }
      state.counters["cursor_bytes"] = 0;
      state.SetItemsProcessed(state.iterations() * kQueries);
      return;
    }
    using Node = decltype(tree.root());
    std::vector<Node> nodes;
    nodes.reserve(kQueries);
    for (auto query : queries) {
      nodes.push_back(tree.find_by_preorder(query % n));
    }
    if (work == CoverWork::kCursorRank || work == CoverWork::kFetch) {
      if (work == CoverWork::kCursorRank) {
        for (auto _ : state) {
          for (auto node : nodes) {
            benchmark::DoNotOptimize(tree.preorder_number(node));
          }
        }
        state.counters["cursor_bytes"] = 0;
      } else {
        for (auto _ : state) {
          for (auto node : nodes) {
            auto cursor = tree.navigate(node);
            benchmark::DoNotOptimize(cursor.is_leaf());
          }
        }
      }
      state.SetItemsProcessed(state.iterations() * kQueries);
      return;
    }

    // Select an actual FCNS edge outside timing. The local row retains one
    // warm component; the boundary row decodes both ends on each round trip.
    // These are navigation probes, not isolated codec or cold-cache timings.
    Node source = tree.root(), destination = source;
    bool first_child = false, found = false;
    const bool boundary = work == CoverWork::kBoundaryEdge;
    auto choose = [&](Node node) {
      auto cursor = tree.navigate(node);
      auto edge = [&](Node target, bool first) {
        if ((node.component != target.component) == boundary) {
          source = node;
          destination = target;
          first_child = first;
          found = true;
        }
      };
      if (!cursor.is_leaf()) {
        edge(cursor.first_child(), true);
      }
      if (!found && !cursor.is_root() && !cursor.is_last_child()) {
        edge(cursor.next_sibling(), false);
      }
    };
    for (auto node : nodes) {
      choose(node);
      if (found) {
        break;
      }
    }
    for (size_t rank = 0; !found && rank < n; ++rank) {
      choose(tree.find_by_preorder(rank));
    }
    if (!found) {
      state.SkipWithError("Shape/size has no eligible component edge");
      return;
    }
    auto cursor = tree.navigate(source);
    benchmark::DoNotOptimize(cursor.is_leaf());
    for (auto _ : state) {
      cursor.reset(source);
      benchmark::DoNotOptimize(cursor.is_leaf());
      // A mutable result uses Benchmark's read/write optimization barrier.
      auto target = first_child ? cursor.first_child() : cursor.next_sibling();
      benchmark::DoNotOptimize(target);
      cursor.reset(target);
      benchmark::DoNotOptimize(cursor.is_leaf());
    }
    state.counters["edge_source_rank"] = tree.preorder_number(source);
    state.counters["edge_destination_rank"] = tree.preorder_number(destination);
    state.counters["edge_is_first_child"] = first_child;
    state.counters["cursor_component_changes_per_iteration"] =
        boundary ? 2.0 - 1.0 / state.iterations() : 0;
    state.SetItemsProcessed(state.iterations());
  }
}

template <class Tree, bool Optimized = false, bool Preorder = true>
void RegisterCoverBenchmarks(const char* implementation) {
  constexpr std::array shapes{"random_recursive", "star", "chain", "comb",
                              "balanced"};
  constexpr std::array operations{"build_destroy",      "dfs",
                                  "dfs_rank",           "inspected_dfs",
                                  "inspected_dfs_rank", "lookup",
                                  "cursor_rank",        "root_path",
                                  "fetch_plus_leaf",    "local_edge_warm",
                                  "boundary_edge_pair"};
  for (size_t shape = 0; shape < shapes.size(); ++shape) {
    for (size_t operation = 0; operation < operations.size(); ++operation) {
      const auto work = static_cast<CoverWork>(operation);
      const bool traversal =
          work >= CoverWork::kDfs && work <= CoverWork::kInspectRank;
      if (Optimized && !traversal) {
        continue;
      }
      if (!Preorder &&
          (work == CoverWork::kDfsRank || work == CoverWork::kInspectRank)) {
        continue;
      }
      if (!kDecomposed<Tree> && !traversal && work != CoverWork::kBuild &&
          work != CoverWork::kRootPath) {
        continue;
      }
      const auto name = std::string("cover/") + implementation + "/" +
                        operations[operation] + "/" + shapes[shape];
      auto* row = benchmark::RegisterBenchmark(
          name.c_str(), [=](benchmark::State& state) {
            CoverBenchmark<Tree, Optimized>(
                state, static_cast<CoverShape>(shape), work);
          });
      row->ArgNames({"N", "seed"});
      for (const int64_t size : {128, 1024, 65536, 1048576, 4194304}) {
        if (work == CoverWork::kBoundaryEdge && size == 128) {
          continue;  // A threshold-sized tree has no component boundary.
        }
        row->Args({size, 42});
        if (shape == 0) {
          row->Args({size, 1729});
        }
      }
    }
  }
}

const bool kCoverRegistered = [] {
  RegisterCoverBenchmarks<pixie::tree::experimental::BTreeDecomposedTree<>>(
      "btree");
  RegisterCoverBenchmarks<pixie::tree::experimental::MacroDecomposedTree<>>(
      "macro");
  using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
  using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  RegisterCoverBenchmarks<Bp>("bp_cursor");
  RegisterCoverBenchmarks<Dfuds>("dfuds_cursor");
  RegisterCoverBenchmarks<Bp, true>("bp_optimized");
  RegisterCoverBenchmarks<Dfuds, true>("dfuds_optimized");
  RegisterCoverBenchmarks<pixie::tree::LoudsTree<>, false, false>(
      "louds_cursor");
  return true;
}();

// These rows use layer-major numbering, never original-tree preorder ranks.
// The final forest has no depth bound once the boundary budget stops splitting.
using PortalMode = pixie::tree::experimental::LayeredDFUDSTree::PortalMode;
void LayeredMemory(benchmark::State& state,
                   const pixie::tree::experimental::LayeredDFUDSTree& tree) {
  const auto m = tree.memory_breakdown();
  state.counters["nodes"] = tree.size();
  state.counters["total_bytes"] = m.total_bytes;
  state.counters["allocated_bits"] = 8.0 * m.total_bytes;
  state.counters["bits_per_node"] = 8.0 * m.total_bytes / tree.size();
  state.counters["object_bits"] = 8.0 * m.object_bytes;
  state.counters["forest_allocated_bits"] = 8.0 * m.forest_bytes;
  state.counters["connector_allocated_bits"] = 8.0 * m.connector_bytes;
  state.counters["boundary_allocated_bits"] = 8.0 * m.boundary_bytes;
  state.counters["forest_index_allocated_bits"] = 8.0 * m.forest_index_bytes;
  state.counters["connector_index_allocated_bits"] =
      8.0 * m.connector_index_bytes;
  state.counters["boundary_index_allocated_bits"] =
      8.0 * m.boundary_index_bytes;
  state.counters["directory_allocated_bits"] = 8.0 * m.directory_bytes;
  state.counters["logical_payload_bits"] = tree.topology_bits();
  state.counters["logical_with_portals_bits"] =
      tree.topology_bits() + tree.boundary_bits();
  state.counters["sparse_portals"] = tree.portal_mode() == PortalMode::Sparse;
  state.counters["forest_logical_bits"] = tree.forest_bits();
  state.counters["connector_logical_bits"] = tree.connector_bits();
  state.counters["boundary_logical_bits"] = tree.boundary_bits();
  state.counters["boundary_charge_bits"] = tree.boundary_overhead_bits();
  state.counters["boundary_budget_bits"] =
      tree.boundary_budget_bits_per_node() * tree.size();
  state.counters["boundary_budget_bits_per_node"] =
      tree.boundary_budget_bits_per_node();
  state.counters["band_height"] = tree.band_height();
  state.counters["layers"] = tree.layer_count();
  state.counters["boundary_nodes"] = tree.boundary_count();
  state.counters["cross_edges"] = tree.cross_edge_count();
  const auto final = tree.layer(tree.layer_count() - 1);
  state.counters["final_layer_nodes"] = final.node_end - final.node_begin;
  state.counters["final_layer_depths"] = final.depth_end - final.depth_begin;
  state.counters["final_layer_depth_unbounded"] = 1;
  // Requested retained bytes plus object, not allocator overhead/RSS or peak
  // construction/decoder stack. The outgoing budget excludes the forest index.
  state.counters["cursor_bytes"] = tree.navigation_scratch_bytes();
  state.SetLabel("layer-major numbering; final forest depth-unbounded");
}

void LayeredBenchmark(benchmark::State& state,
                      CoverShape shape,
                      CoverWork work,
                      size_t band_height,
                      double budget,
                      PortalMode mode) {
  using Tree = pixie::tree::experimental::LayeredDFUDSTree;
  const size_t n = state.range(0);
  const uint64_t seed = state.range(1);
  if (work == CoverWork::kBuild) {
    const CoverSource source(n, shape, seed);
    for (auto _ : state) {
      Tree built(source, band_height, budget, mode);
      benchmark::DoNotOptimize(built);
    }
    const Tree measured(source, band_height, budget, mode);
    LayeredMemory(state, measured);
    state.counters.erase("cursor_bytes");
    state.SetItemsProcessed(state.iterations() * n);
    return;
  }
  const Tree tree =
      BuildCoverTree<Tree>(state, shape, band_height, budget, mode);
  LayeredMemory(state, tree);
  if (work == CoverWork::kInspectedRootPath) {
    return CoverInspectedRootPath(state, tree);
  }
  if (work == CoverWork::kDfs) {
    return CoverDfs<false, false, false>(state, tree);
  }
  if (work == CoverWork::kInspect) {
    return CoverDfs<true, false, false>(state, tree);
  }

  constexpr size_t kQueries = 8192;
  std::array<uint64_t, kQueries> queries;
  std::mt19937_64 random(seed ^ 0x9e3779b97f4a7c15ULL);
  for (auto& query : queries) {
    query = random();
  }
  if (work == CoverWork::kLookup) {
    for (auto& query : queries) {
      query %= n;
    }
    for (auto _ : state) {
      for (auto query : queries) {
        benchmark::DoNotOptimize(tree.find_by_number(query));
      }
    }
    state.counters["cursor_bytes"] = 0;
    state.SetItemsProcessed(state.iterations() * kQueries);
    return;
  }

  // Same query stream, indexed child selection, and fresh path cursor as the
  // DFUDS cursor baseline. No lookup or numbering conversion is timed here.
  size_t position = 0, visited = 0;
  for (auto _ : state) {
    auto cursor = tree.navigate(tree.root());
    for (;;) {
      CoverObserve<false>(tree, cursor.node());
      ++visited;
      if (cursor.is_leaf()) {
        break;
      }
      const auto index = queries[position++ % kQueries] % cursor.degree();
      cursor.reset(cursor.child(index));
    }
  }
  state.counters["nodes_per_path"] =
      static_cast<double>(visited) / state.iterations();
  state.SetItemsProcessed(visited);
}

template <bool Optimized>
void LayeredBaseline(benchmark::State& state,
                     CoverShape shape,
                     CoverWork work) {
  using Tree = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  if constexpr (Optimized) {
    if (work == CoverWork::kDfs) {
      const Tree tree = [&] {
        const CoverSource source(state.range(0), shape, state.range(1));
        return Tree(source);
      }();
      CoverMemory(state, tree);
      const auto budget = pixie::tree::kDefaultTraversalBytes;
      pixie::benchmarks::TreeTraversalWorkloads walk(tree, budget);
      pixie::tree::RootNavigationIndex index(tree,
                                             budget - walk.reserved_bytes());
      for (auto _ : state) {
        // Streaming traversal emits unary entries AND returns. Preorder entry
        // numbers increase by one; returns revisit an already entered node.
        size_t next_entry = 0;
        auto entry = [&](auto node) {
          if (node.number == next_entry) {
            ++next_entry;
            CoverObserve<false>(tree, node);
          }
        };
        walk.traverse(entry);
        benchmark::DoNotOptimize(next_entry);
      }
      const auto retained =
          walk.fixture_metadata_bytes() + index.memory_usage_bytes();
      state.counters["traversal_index_bytes"] = retained;
      state.counters["total_bytes"] = tree.memory_usage_bytes() + retained;
      state.counters["bits_per_node"] =
          8.0 * (tree.memory_usage_bytes() + retained) / tree.size();
      state.counters["scratch_bytes"] = walk.scratch_bytes();
      state.counters.erase("cursor_bytes");
      state.SetItemsProcessed(state.iterations() * tree.size());
    } else {
      CoverBenchmark<Tree, true>(state, shape, work);
    }
  } else {
    CoverBenchmark<Tree>(state, shape, work);
  }
  const double traversal =
      Optimized ? state.counters["traversal_index_bytes"].value : 0;
  const double total = state.counters["total_bytes"].value;
  const double encoding = state.counters["encoding_bytes"].value;
  state.counters["allocated_bits"] = 8 * total;
  state.counters["object_bits"] = 8 * sizeof(Tree);
  state.counters["encoding_allocated_bits"] = 8 * encoding;
  state.counters["tree_index_allocated_bits"] =
      8 * (total - sizeof(Tree) - encoding - traversal);
  state.counters["traversal_index_allocated_bits"] = 8 * traversal;
  state.counters["logical_payload_bits"] = 2.0 * state.range(0) - 1;
}

namespace f256 = pixie::tree::detail::forest256;
namespace twin = pixie::benchmarks::forest256_twin;

// Two concatenated logical units, with no per-unit padding. The topology and
// local observations are identical for the two codecs and all cache budgets.
template <bool Twin>
std::vector<uint64_t> ForestStream(size_t n, size_t offset) {
  std::vector<uint64_t> stream((offset + 4 * n + 63) / 64, ~uint64_t{});
  for (size_t unit = 0; unit < 2; ++unit) {
    const size_t roots = unit + 1;
    std::vector<uint16_t> degrees(n);
    for (size_t child = roots; child < n; ++child) {
      ++degrees[(child - roots) / (roots + 1)];
    }
    const auto packed = [&] {
      if constexpr (Twin) {
        return twin::encode_twin_forest256(roots, degrees);
      } else {
        return f256::encode_louds_forest256(roots, degrees);
      }
    }();
    for (size_t i = 0; i < 2 * n; ++i) {
      const size_t src = Twin && i >= n ? 256 + i - n : i;
      if (!((packed.words[src / 64] >> (src % 64)) & 1)) {
        const size_t dst = offset + unit * 2 * n + i;
        stream[dst / 64] &= ~(uint64_t{1} << (dst % 64));
      }
    }
  }
  return stream;
}

enum class ForestProbe { Scalar, Native, Cache };

template <ForestProbe Probe, bool Twin>
void ForestBenchmark(benchmark::State& state, size_t peer) {
  const size_t n = state.range(0), offset = state.range(1);
  const auto words = ForestStream<Twin>(n, offset);
  using Cache = std::conditional_t<Twin, twin::Cache, f256::Forest256Cache>;
  const size_t budget =
      Probe == ForestProbe::Cache ? state.range(2) * Cache::slot_bytes() : 0;
  Cache cache(words, budget);
  constexpr size_t batch = 256;
  auto replay = [&](auto&& lookup) {
    for (size_t i = 0; i < batch; ++i) {
      const size_t second = peer != 0 && (i & 1);
      lookup(second * peer, offset + second * 2 * n, uint8_t((73 * i) % n));
    }
  };
  typename Cache::Stats stats;
  if constexpr (Probe == ForestProbe::Cache) {
    // Warm once and count one steady-state batch outside the timed path.
    for (int pass = 0; pass < 2; ++pass) {
      stats = {};
      replay([&](size_t unit, size_t bit, uint8_t u) {
        cache.lookup_counted(unit, bit, n, u, stats);
      });
    }
  }
  for (auto _ : state) {
    replay([&](size_t unit, size_t bit, uint8_t u) {
      f256::NodeInfo node;
      if constexpr (Probe == ForestProbe::Cache) {
        node = cache.lookup(unit, bit, n, u);
      } else {
        f256::LevelCsr256 csr;
        constexpr auto decode =
            Probe == ForestProbe::Scalar
                ? (Twin ? twin::unpack_dense_scalar : f256::unpack_dense_scalar)
                : (Twin ? twin::unpack_dense : f256::unpack_dense);
        decode(words, bit, n, csr);
        node = {uint16_t(csr.root_count()), uint16_t(csr.offset_m1[u] + 1),
                uint16_t(csr.degree(u)), csr.parent_of(u),
                csr.is_last_child(u)};
      }
      benchmark::DoNotOptimize(node);
    });
  }
  state.counters["cache_budget_bytes"] = budget;
  state.counters["cache_allocated_bytes"] = cache.allocated_bytes();
  state.counters["cache_reader_bytes"] =
      Probe == ForestProbe::Cache ? sizeof(Cache) : 0;
  state.counters["cache_slots"] = cache.slot_count();
  state.counters["csr_bytes"] = sizeof(f256::LevelCsr256);
  if constexpr (Probe == ForestProbe::Cache) {
    state.counters["hits_per_query"] = double(stats.hits) / batch;
    state.counters["misses_per_query"] = double(stats.misses) / batch;
    state.counters["evictions_per_query"] = double(stats.evictions) / batch;
  }
  state.counters["packed_bits_per_query"] =
      Probe == ForestProbe::Cache ? double(stats.packed_bits) / batch : 2 * n;
  state.SetItemsProcessed(state.iterations() * batch);
  state.SetLabel("dense unit lookup; no tree controller");
}

const bool kForestRegistered = [] {
  auto register_codec = []<bool Twin>(const char* name) {
    for (size_t peer : {0, 1, 4}) {
      const auto prefix = std::string("forest256/dense/") + name + "/" +
                          (peer == 0   ? "repeat/"
                           : peer == 1 ? "pair/"
                                       : "conflict/");
      auto add = [&]<ForestProbe Probe>(const char* suffix) {
        auto* row = benchmark::RegisterBenchmark(
            (prefix + suffix).c_str(), [=](benchmark::State& s) {
              ForestBenchmark<Probe, Twin>(s, peer);
            });
        if constexpr (Probe == ForestProbe::Cache) {
          row->ArgNames({"M", "offset", "K"});
        } else {
          row->ArgNames({"M", "offset"});
        }
        for (int n : {127, 256}) {
          for (int offset : {0, 1, 63}) {
            if constexpr (Probe == ForestProbe::Cache) {
              for (int k : {0, 1, 4}) {
                row->Args({n, offset, k});
              }
            } else {
              row->Args({n, offset});
            }
          }
        }
      };
      add.template operator()<ForestProbe::Scalar>("scalar");
      add.template operator()<ForestProbe::Native>("native");
      add.template operator()<ForestProbe::Cache>("cache");
    }
  };
  register_codec.template operator()<false>("louds");
  register_codec.template operator()<true>("twin");
  return true;
}();

using Ordinal = pixie::tree::ForestCoverTree;
using ForestContext = pixie::tree::TreeTraversalContext<Ordinal>;
using ForestIndex = pixie::tree::RootNavigationIndex<Ordinal>;

// Adapt only the handle-oriented reader, without rechecking returned handles.
// Each call to navigate constructs a new, initially empty reader for this walk.
struct OrdinalWalk {
  const Ordinal& tree;
  size_t budget;
  auto root() const { return tree.root(); }
  size_t size() const { return tree.size(); }
  struct Cursor {
    ForestContext reader;
    Ordinal::Node current;
    auto node() const { return current; }
    void reset(Ordinal::Node h) { current = h; }
    bool is_root() const { return current == Ordinal::Node{}; }
    bool is_leaf() { return reader.is_leaf(current); }
    bool is_last_child() { return reader.is_last_child(current); }
    auto first_child() { return reader.first_child(current); }
    auto next_sibling() { return reader.next_sibling(current); }
    auto parent() { return reader.parent(current); }
    auto degree() { return reader.degree(current); }
    void for_each_child(auto& f) { reader.for_each_child(current, f); }
  };
  Cursor navigate(Ordinal::Node h) const {
    return {ForestContext(tree, budget), h};
  }
  void traverse_inspected(auto& visit, auto& child) const {
    ForestContext(tree, budget).traverse_inspected(visit, child);
  }
  void inspect_path(auto& visit, auto& child, auto& choose) const {
    ForestContext(tree, budget).inspect_path(visit, child, choose);
  }
};

void OrdinalMemory(benchmark::State& state, const Ordinal& tree) {
  state.counters["nodes"] = tree.size();
  state.counters["units"] = tree.unit_count();
  state.counters["total_bytes"] = tree.memory_usage_bytes();
  state.counters["bits_per_node"] =
      8.0 * tree.memory_usage_bytes() / tree.size();
}

struct PersistentOrdinalWalk : OrdinalWalk {
  ForestContext& reader;
  void traverse_inspected(auto& visit, auto& child) const {
    reader.traverse_inspected(visit, child);
  }
  void inspect_path(auto& visit, auto& child, auto& choose) const {
    reader.inspect_path(visit, child, choose);
  }
};

struct UpperOrdinalWalk : PersistentOrdinalWalk {
  const ForestIndex& index;
  bool fresh;
  void inspect_path(auto& visit, auto& child, auto& choose) const {
    if (fresh) {
      ForestContext r(tree, budget, &index);
      r.inspect_path(visit, child, choose);
    } else {
      reader.inspect_path(visit, child, choose);
    }
  }
};

const bool kOrdinalUpperRegistered = [] {
  for (auto [shape, name] :
       {std::pair{CoverShape::kRandomRecursive, "random_recursive"},
        std::pair{CoverShape::kBinaryCatalan, "binary_catalan"},
        std::pair{CoverShape::kBalanced, "balanced"}}) {
    for (bool fresh : {false, true}) {
      auto row = std::string("ordinal/upper/inspected_root_path/") +
                 (fresh ? "fresh/" : "persistent/") + name;
      benchmark::RegisterBenchmark(
          row.c_str(),
          [=](benchmark::State& s) {
            const auto tree = BuildCoverTree<Ordinal>(s, shape);
            const size_t budget =
                s.range(2) * f256::Forest256Cache::slot_bytes();
            const auto start = std::chrono::steady_clock::now();
            ForestIndex index(tree, s.range(3));
            s.counters["index_prepare_seconds"] =
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              start)
                    .count();
            ForestContext reader(tree, fresh ? 0 : budget, &index);
            s.counters["index_charge_bytes"] = index.budget_usage_bytes();
            s.counters["index_entries"] = index.size();
            s.counters["reader_slot_bytes"] = budget;
            const UpperOrdinalWalk adapter{
                {{tree, budget}, reader}, index, fresh};
            CoverInspectedRootPath(s, adapter);
          })
          ->ArgNames({"N", "seed", "K", "index_bytes"})
          ->Args({1048576, 42, 16, 0})
          ->Args({1048576, 42, 16, 16 * f256::Forest256Cache::slot_bytes()})
          ->Args({1048576, 42, 32, 0});
    }
  }
  return true;
}();

void OrdinalBenchmark(benchmark::State& state,
                      CoverShape shape,
                      CoverWork work,
                      bool persistent = false) {
  const size_t n = state.range(0), seed = state.range(1);
  if (work == CoverWork::kBuild) {
    const CoverSource source(n, shape, seed);
    for (auto _ : state) {
      Ordinal built(source);
      benchmark::DoNotOptimize(built);
    }
    OrdinalMemory(state, Ordinal(source));
    state.SetItemsProcessed(state.iterations() * n);
    return;
  }
  const auto tree = BuildCoverTree<Ordinal>(state, shape);
  OrdinalMemory(state, tree);
  const size_t budget =
      work == CoverWork::kLookup
          ? 0
          : state.range(2) * f256::Forest256Cache::slot_bytes();
  const bool walk = work == CoverWork::kDfs || work == CoverWork::kInspect ||
                    work == CoverWork::kInspectedRootPath;
  state.counters["cache_budget_bytes"] = budget;
  state.counters["cache_allocated_bytes"] =
      work == CoverWork::kLookup
          ? 0
          : ForestContext(tree, budget).allocated_bytes();
  state.counters["reader_fixed_bytes"] =
      work == CoverWork::kLookup ? 0 : sizeof(ForestContext);
  state.counters["cursor_bytes"] =
      work == CoverWork::kDfs ? sizeof(OrdinalWalk::Cursor) : 0;
  if (walk) {
    if (persistent) {
      const auto start = std::chrono::steady_clock::now();
      ForestContext reader(tree, budget);
      state.counters["reader_prepare_seconds"] =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        start)
              .count();
      state.SetLabel("reader allocated before timing; initially empty");
      const PersistentOrdinalWalk adapter{{tree, budget}, reader};
      if (work == CoverWork::kInspect) {
        return CoverDfs<true, false, false>(state, adapter);
      }
      return CoverInspectedRootPath(state, adapter);
    }
    state.SetLabel("fresh reader per walk, allocation timed");
    const OrdinalWalk adapter{tree, budget};
    if (work == CoverWork::kDfs) {
      return CoverDfs<false, false, false>(state, adapter);
    }
    if (work == CoverWork::kInspect) {
      return CoverDfs<true, false, false>(state, adapter);
    }
    return CoverInspectedRootPath(state, adapter);
  }
  std::array<size_t, 8192> ids;
  std::mt19937_64 random(seed ^ 0x9e3779b97f4a7c15ULL);
  for (auto& id : ids) {
    id = random() % n;
  }
  if (work == CoverWork::kLookup) {
    state.SetLabel("uniform stable ID to handle; no reader or decode");
    for (auto _ : state) {
      for (auto id : ids) {
        benchmark::DoNotOptimize(tree.find_by_number(id));
      }
    }
  } else {
    // Alternate units 0 and 1; K=1 misses every query, K=4 hits every query
    // after priming. Both observe the same varying IDs and original degrees.
    const size_t b = tree.number({1, 0}), e = tree.number({2, 0});
    for (size_t i = 0; i < ids.size(); ++i) {
      ids[i] = i % 2 ? b + ids[i] % (e - b) : ids[i] % b;
    }
    ForestContext reader(tree, budget);
    benchmark::DoNotOptimize(reader.degree(tree.find_by_number(ids[0])));
    benchmark::DoNotOptimize(reader.degree(tree.find_by_number(ids[1])));
    state.SetLabel("two alternating units; primed persistent reader");
    for (auto _ : state) {
      for (auto id : ids) {
        auto h = tree.find_by_number(id);
        benchmark::DoNotOptimize(h);
        benchmark::DoNotOptimize(reader.degree(h));
      }
    }
  }
  state.SetItemsProcessed(state.iterations() * ids.size());
}

const bool kLayeredRegistered = [] {
  constexpr std::array shapes{"random_recursive", "star",   "chain", "comb",
                              "balanced",         "catalan"};
  constexpr std::array operations{"build_destroy",
                                  "dfs",
                                  "inspected_dfs",
                                  "indexed_child_root_path",
                                  "layer_major_number_lookup",
                                  "inspected_root_path"};
  constexpr std::array work{CoverWork::kBuild,   CoverWork::kDfs,
                            CoverWork::kInspect, CoverWork::kRootPath,
                            CoverWork::kLookup,  CoverWork::kInspectedRootPath};
  struct Configuration {
    size_t height;
    double budget;
    const char* name;
    bool small_only = false;
    PortalMode mode = PortalMode::Dense;
  };
  // One-factor controls, not a height x budget cross product. The generous
  // small-input control exercises splits when tight budgets retain one forest.
  constexpr std::array configs{
      Configuration{8, 0.05, "h8_b0.05"},
      Configuration{1, 0.05, "h1_b0.05"},
      Configuration{4, 0.05, "h4_b0.05"},
      Configuration{16, 0.05, "h16_b0.05"},
      Configuration{8, 0, "h8_b0"},
      Configuration{8, 0.25, "h8_b0.25"},
      Configuration{8, 64, "h8_b64", true},
      Configuration{6, 0.05, "h6_b0.05"},
      Configuration{6, 0.05, "positive_sparse/h6_b0.05", false,
                    PortalMode::Sparse},
      Configuration{8, 0.05, "positive_sparse/h8_b0.05", false,
                    PortalMode::Sparse}};
  for (size_t shape = 0; shape < shapes.size(); ++shape) {
    auto ordinal = [&](const char* op, CoverWork kind, bool reader,
                       bool persistent = false) {
      const auto name =
          std::string("ordinal/forest256/") + op + "/" + shapes[shape];
      auto* row = benchmark::RegisterBenchmark(
          name.c_str(), [=](benchmark::State& state) {
            OrdinalBenchmark(state, static_cast<CoverShape>(shape), kind,
                             persistent);
          });
      if (reader) {
        row->ArgNames({"N", "seed", "K"});
        for (int n : {65536, 1048576}) {
          for (int k : {0, 1, 4}) {
            row->Args({n, 42, k});
            if ((shape == 0 || shape == 5) &&
                (kind == CoverWork::kInspect ||
                 kind == CoverWork::kInspectedRootPath)) {
              row->Args({n, 1729, k});
            }
          }
        }
      } else {
        row->ArgNames({"N", "seed"})->Args({65536, 42})->Args({1048576, 42});
      }
    };
    ordinal("build_destroy", CoverWork::kBuild, false);
    ordinal("dfs", CoverWork::kDfs, true);
    ordinal("inspected_dfs", CoverWork::kInspect, true);
    ordinal("inspected_root_path", CoverWork::kInspectedRootPath, true);
    ordinal("persistent/inspected_dfs", CoverWork::kInspect, true, true);
    ordinal("persistent/inspected_root_path", CoverWork::kInspectedRootPath,
            true, true);
    ordinal("unit_major_number_lookup", CoverWork::kLookup, false);
    ordinal("lookup_degree_pair", CoverWork::kFetch, true);
    auto optimized = [&]<class Tree>(const char* impl) {
      for (size_t op : {2, 5}) {
        auto name = std::string("ordinal/") + impl + "/" + operations[op] +
                    "/" + shapes[shape];
        auto* row = benchmark::RegisterBenchmark(
            name.c_str(), [=](benchmark::State& s) {
              CoverBenchmark<Tree, true>(s, static_cast<CoverShape>(shape),
                                         work[op]);
            });
        row->ArgNames({"N", "seed"});
        for (int n : {65536, 1048576}) {
          row->Args({n, 42});
          if (shape == 0 || shape == 5) {
            row->Args({n, 1729});
          }
        }
      }
    };
    optimized.operator()<pixie::tree::BPTree<pixie::RmMBTree<>>>(
        "bp_optimized");
    optimized.operator()<pixie::tree::DFUDSTree<pixie::RmMBTree<>>>(
        "dfuds_optimized");
    for (size_t op = 0; op < work.size(); ++op) {
      const bool inspected = work[op] == CoverWork::kInspect ||
                             work[op] == CoverWork::kInspectedRootPath;
      if (shape == 5 && !inspected) {
        continue;
      }
      const auto suffix =
          std::string("/") + operations[op] + "/" + shapes[shape];
      for (auto config : configs) {
        const auto name = std::string("layered/dfuds/") + config.name + suffix;
        auto* row = benchmark::RegisterBenchmark(
            name.c_str(), [=](benchmark::State& state) {
              LayeredBenchmark(state, static_cast<CoverShape>(shape), work[op],
                               config.height, config.budget, config.mode);
            });
        row->ArgNames({"N", "seed"});
        for (int64_t n : {4096, 65536, 1048576}) {
          if (!config.small_only || n == 4096) {
            row->Args({n, 42});
            if (inspected && (shape == 0 || shape == 5) && n == 1048576) {
              row->Args({n, 1729});
            }
          }
        }
      }
      if (shape == 5 || work[op] == CoverWork::kLookup) {
        continue;  // Ordinary DFUDS has no layer-major number lookup contract.
      }
      const auto louds_name = std::string("layered/louds_cursor") + suffix;
      benchmark::RegisterBenchmark(louds_name.c_str(),
                                   [=](benchmark::State& s) {
                                     CoverBenchmark<pixie::tree::LoudsTree<>>(
                                         s, static_cast<CoverShape>(shape),
                                         work[op]);
                                   })
          ->ArgNames({"N", "seed"})
          ->Args({4096, 42})
          ->Args({65536, 42})
          ->Args({1048576, 42});
      const auto cursor_name = std::string("layered/dfuds_cursor") + suffix;
      benchmark::RegisterBenchmark(
          cursor_name.c_str(),
          [=](benchmark::State& s) {
            LayeredBaseline<false>(s, static_cast<CoverShape>(shape), work[op]);
          })
          ->ArgNames({"N", "seed"})
          ->Args({4096, 42})
          ->Args({65536, 42})
          ->Args({1048576, 42});
      if (work[op] == CoverWork::kDfs || work[op] == CoverWork::kInspect) {
        const auto optimized_name =
            std::string(work[op] == CoverWork::kDfs
                            ? "layered/dfuds_optimized_entry_filtered"
                            : "layered/dfuds_optimized") +
            suffix;
        benchmark::RegisterBenchmark(optimized_name.c_str(),
                                     [=](benchmark::State& s) {
                                       LayeredBaseline<true>(
                                           s, static_cast<CoverShape>(shape),
                                           work[op]);
                                     })
            ->ArgNames({"N", "seed"})
            ->Args({4096, 42})
            ->Args({65536, 42})
            ->Args({1048576, 42});
      }
    }
  }
  return true;
}();
}  // namespace
BENCHMARK_TEMPLATE(Prepared, false, false)
    ->Name("bp_prepared_root")
    ->Arg(1024)
    ->Arg(4096)
    ->Arg(4097);
BENCHMARK_TEMPLATE(Prepared, false, true)
    ->Name("bp_prepared_events")
    ->Arg(1024)
    ->Arg(4096)
    ->Arg(4097);
BENCHMARK_TEMPLATE(Prepared, true, false)
    ->Name("dfuds_prepared_root")
    ->Arg(1024)
    ->Arg(4096)
    ->Arg(4097);
BENCHMARK_TEMPLATE(Prepared, true, true)
    ->Name("dfuds_prepared_events")
    ->Arg(1024)
    ->Arg(4096)
    ->Arg(4097);

// Random trees: Four-Russians (8192-bit chunk + 197-entry LUT) vs SIMD
// Microtree
BENCHMARK_TEMPLATE(MicroDfsChildren, false, false, 8192, ExpShape::kRandom)
    ->Name("bp_4russians_8192_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, false, true, 256, ExpShape::kRandom)
    ->Name("bp_simd_256_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, false, true, 512, ExpShape::kRandom)
    ->Name("bp_simd_512_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, false, true, 2048, ExpShape::kRandom)
    ->Name("bp_simd_2048_random")
    ->Apply(SimdTreeSizes);

BENCHMARK_TEMPLATE(MicroDfsChildren, true, false, 8192, ExpShape::kRandom)
    ->Name("dfuds_4russians_8192_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, true, true, 256, ExpShape::kRandom)
    ->Name("dfuds_simd_256_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, true, true, 512, ExpShape::kRandom)
    ->Name("dfuds_simd_512_random")
    ->Apply(SimdTreeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, true, true, 2048, ExpShape::kRandom)
    ->Name("dfuds_simd_2048_random")
    ->Apply(SimdTreeSizes);

// Star & clustered-hub shapes
BENCHMARK_TEMPLATE(MicroDfsChildren, false, false, 8192, ExpShape::kStar)
    ->Name("bp_4russians_8192_star")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, false, true, 512, ExpShape::kStar)
    ->Name("bp_simd_512_star")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren,
                   false,
                   false,
                   8192,
                   ExpShape::kClusteredTwigs)
    ->Name("bp_4russians_8192_hubs")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren,
                   false,
                   true,
                   512,
                   ExpShape::kClusteredTwigs)
    ->Name("bp_simd_512_hubs")
    ->Apply(ShapeSizes);

BENCHMARK_TEMPLATE(MicroDfsChildren, true, false, 8192, ExpShape::kStar)
    ->Name("dfuds_4russians_8192_star")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, true, true, 512, ExpShape::kStar)
    ->Name("dfuds_simd_512_star")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren,
                   true,
                   false,
                   8192,
                   ExpShape::kClusteredTwigs)
    ->Name("dfuds_4russians_8192_hubs")
    ->Apply(ShapeSizes);
BENCHMARK_TEMPLATE(MicroDfsChildren, true, true, 512, ExpShape::kClusteredTwigs)
    ->Name("dfuds_simd_512_hubs")
    ->Apply(ShapeSizes);
