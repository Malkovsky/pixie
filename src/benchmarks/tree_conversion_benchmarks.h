#pragma once

#include <benchmark/benchmark.h>
#include <pixie/tree/implementations.h>

#include "tree_input.h"

namespace pixie::benchmarks {

template <tree::detail::TreeEncodingKind Kind>
using ConversionTree = std::conditional_t<
    Kind == tree::detail::TreeEncodingKind::kBP,
    tree::BPTree<RmMBTree<>>,
    std::conditional_t<Kind == tree::detail::TreeEncodingKind::kDFUDS,
                       tree::DFUDSTree<RmMBTree<>>,
                       tree::LoudsTree<RankSelectSupportCompressed<>>>>;

// Source preparation is excluded. Encoding-only rows include output allocation
// and all conversion scratch; complete rows also build the destination index.
template <tree::detail::TreeEncodingKind From,
          tree::detail::TreeEncodingKind To,
          bool Complete>
void RunTreeConversion(benchmark::State& state) {
  const size_t n = state.range(0);
  const auto bp_words = GenerateRandomBp(n);
  const tree::BPTree<RmMBTree<>>::TreeView bp(bp_words, n);
  const ConversionTree<From> source(bp);
  {
    const ConversionTree<To> converted(source);
    const tree::BPTree<RmMBTree<>> roundtrip(converted);
    if (!std::ranges::equal(roundtrip.encoding_words(), bp_words)) {
      state.SkipWithError("Conversion changes ordered topology");
      return;
    }
  }
  for (auto _ : state) {
    if constexpr (Complete) {
      ConversionTree<To> output(source);
      benchmark::DoNotOptimize(output.encoding_words().data());
      benchmark::ClobberMemory();
    } else {
      auto output = tree::detail::construct_tree_encoding<To>(source);
      benchmark::DoNotOptimize(output.as_bytes().data());
      benchmark::ClobberMemory();
    }
  }
  state.SetItemsProcessed(state.iterations() * n);
  state.counters["decoder_scratch_bytes"] = tree::detail::ConstructionTraversal<
      ConversionTree<From>>::scratch_bytes();
}

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kBP,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   false)
    ->Name("convert_bp_dfuds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kBP,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   true)
    ->Name("convert_bp_dfuds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kBP,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   false)
    ->Name("convert_bp_louds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kBP,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   true)
    ->Name("convert_bp_louds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   tree::detail::TreeEncodingKind::kBP,
                   false)
    ->Name("convert_dfuds_bp_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   tree::detail::TreeEncodingKind::kBP,
                   true)
    ->Name("convert_dfuds_bp_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   false)
    ->Name("convert_dfuds_louds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   true)
    ->Name("convert_dfuds_louds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   tree::detail::TreeEncodingKind::kBP,
                   false)
    ->Name("convert_louds_bp_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   tree::detail::TreeEncodingKind::kBP,
                   true)
    ->Name("convert_louds_bp_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   false)
    ->Name("convert_louds_dfuds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   tree::detail::TreeEncodingKind::kLOUDS,
                   tree::detail::TreeEncodingKind::kDFUDS,
                   true)
    ->Name("convert_louds_dfuds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

}  // namespace pixie::benchmarks
