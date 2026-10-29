#pragma once

#include <pixie/tree/experimental/tree_traversal_cache.h>

#include <charconv>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

#include "tree_benchmark_helpers.h"

namespace pixie::benchmarks {
// Process-wide fixture configuration; parsed before any timed traversal.
inline size_t TreeCacheBudget() {
  static const size_t budget = [] {
    const char* raw = std::getenv("PIXIE_TREE_CACHE_BYTES");
    if (!raw) {
      return kDefaultTreeCacheBytes;
    }
    const std::string_view text(raw);
    size_t bytes = 0;
    const auto result =
        std::from_chars(text.data(), text.data() + text.size(), bytes);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
      throw std::invalid_argument("Invalid PIXIE_TREE_CACHE_BYTES");
    }
    return bytes;
  }();
  return budget;
}
// One configuration per encoding, with the same index/scratch reservation for
// every workload. No per-size or per-workload selection of policy winners.
template <TreeWorkload Workload, class Tree>
void RunTreeComparison(benchmark::State& state,
                       const Tree& tree,
                       std::span<const uint64_t> words) {
  experimental::TreeTraversalCache<Tree> cache(tree, words);
  if (cache.memory_usage_bytes() > tree.cache_limit_bytes()) {
    state.SkipWithError("Combined cache and scratch budget exceeded");
    return;
  }
  std::optional<ShiftedQueryPool> queries;
  if constexpr (Workload == TreeWorkload::kRandomPath ||
                Workload == TreeWorkload::kRandomPathChildren) {
    queries.emplace(0, std::numeric_limits<size_t>::max(), 43);
  }
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  size_t nodes = 0, children = 0;
  for (auto _ : state) {
    if constexpr (Workload == TreeWorkload::kDfs) {
      cache.traverse(observe);
      nodes += tree.size();
    } else if constexpr (Workload == TreeWorkload::kDfsChildren) {
      auto [n, c] = cache.inspect(observe, observe);
      if (n != tree.size() || c != tree.size() - 1) {
        state.SkipWithError("Incorrect traversal coverage");
        break;
      }
      nodes += n;
      children += c;
    } else {
      auto work = IndexedRandomRootToLeaf<Workload ==
                                          TreeWorkload::kRandomPathChildren>(
          cache, *queries, observe, observe);
      nodes += work.nodes;
      children += work.children;
      benchmark::DoNotOptimize(work);
    }
  }
  state.counters["nodes"] = tree.size();
  state.counters["cache_limit_bytes"] = tree.cache_limit_bytes();
  state.counters["cache_and_scratch_bytes"] = cache.memory_usage_bytes();
  state.counters["scratch_reserved_bytes"] = cache.scratch_bytes();
  state.counters["streaming_dfs"] = cache.streaming_enabled();
  state.counters["decoded_inspection"] = cache.optimized_inspection_enabled();
  state.counters["decode_chunk_bits"] = cache.decoded_chunk_bits();
  state.counters["decode_depth"] = cache.cached_depth();
  state.counters["visited_nodes_per_iteration"] =
      double(nodes) / state.iterations();
  state.counters["inspected_children_per_iteration"] =
      double(children) / state.iterations();
}
}  // namespace pixie::benchmarks
