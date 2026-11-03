#include <benchmark/benchmark.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/experimental/tree_prepared_events.h>
#include <pixie/tree/experimental/tree_prepared_root.h>
#include <pixie/tree/traversal_context.h>

#include "tree_dfs/optimized.h"
#include "tree_input.h"

namespace {
template <bool Dfuds, bool Events>
void Prepared(benchmark::State& state) {
  using Tree = std::conditional_t<Dfuds, pixie::DFUDSTree<pixie::RmMBTree<>>,
                                  pixie::BPTree<pixie::RmMBTree<>>>;
  using Decoder =
      std::conditional_t<Dfuds,
                         pixie::detail::dfuds_micro::DfudsMicroSubtree<8192>,
                         pixie::detail::bp_micro::BpMicroSubtree<8192>>;
  const size_t n = state.range(0);
  auto words = pixie::benchmarks::GenerateRandomBp(n);
  if constexpr (Dfuds) {
    words = pixie::benchmarks::ConvertBpEncoding<
        pixie::detail::TreeEncodingKind::kDFUDS>(words, n);
  }
  Tree tree(words, n);
  pixie::TreeTraversalContext context(tree);
  pixie::benchmarks::TreeTraversalWorkloads fallback(context);
  auto prepared = [&] {
    if constexpr (Events) {
      return pixie::experimental::PreparedTreeEvents<>(tree);
    } else {
      return pixie::experimental::PreparedRootSubtree<Decoder, 8192>(
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
