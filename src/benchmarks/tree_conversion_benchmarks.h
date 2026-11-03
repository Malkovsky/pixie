#pragma once

#include <benchmark/benchmark.h>
#include <pixie/tree/implementations.h>

#include "tree_input.h"

namespace pixie::benchmarks {

template <detail::TreeEncodingKind Kind>
using ConversionTree = std::conditional_t<
    Kind == detail::TreeEncodingKind::kBP,
    BPTree<RmMBTree<>>,
    std::conditional_t<Kind == detail::TreeEncodingKind::kDFUDS,
                       DFUDSTree<RmMBTree<>>,
                       LoudsTree<RankSelectSupportCompressed<>>>>;

// Source preparation is excluded. Encoding-only rows include output allocation
// and all conversion scratch; complete rows also build the destination index.
template <detail::TreeEncodingKind From,
          detail::TreeEncodingKind To,
          bool Complete>
void RunTreeConversion(benchmark::State& state) {
  const size_t n = state.range(0);
  const auto bp_words = GenerateRandomBp(n);
  const BPTree<RmMBTree<>>::TreeView bp(bp_words, n);
  const ConversionTree<From> source(bp);
  {
    const ConversionTree<To> converted(source);
    const BPTree<RmMBTree<>> roundtrip(converted);
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
      auto output = detail::construct_tree_encoding<To>(source);
      benchmark::DoNotOptimize(output.as_bytes().data());
      benchmark::ClobberMemory();
    }
  }
  state.SetItemsProcessed(state.iterations() * n);
  state.counters["decoder_scratch_bytes"] =
      detail::ConstructionTraversal<ConversionTree<From>>::scratch_bytes();
}

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kBP,
                   detail::TreeEncodingKind::kDFUDS,
                   false)
    ->Name("convert_bp_dfuds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kBP,
                   detail::TreeEncodingKind::kDFUDS,
                   true)
    ->Name("convert_bp_dfuds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kBP,
                   detail::TreeEncodingKind::kLOUDS,
                   false)
    ->Name("convert_bp_louds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kBP,
                   detail::TreeEncodingKind::kLOUDS,
                   true)
    ->Name("convert_bp_louds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kDFUDS,
                   detail::TreeEncodingKind::kBP,
                   false)
    ->Name("convert_dfuds_bp_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kDFUDS,
                   detail::TreeEncodingKind::kBP,
                   true)
    ->Name("convert_dfuds_bp_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kDFUDS,
                   detail::TreeEncodingKind::kLOUDS,
                   false)
    ->Name("convert_dfuds_louds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kDFUDS,
                   detail::TreeEncodingKind::kLOUDS,
                   true)
    ->Name("convert_dfuds_louds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kLOUDS,
                   detail::TreeEncodingKind::kBP,
                   false)
    ->Name("convert_louds_bp_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kLOUDS,
                   detail::TreeEncodingKind::kBP,
                   true)
    ->Name("convert_louds_bp_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kLOUDS,
                   detail::TreeEncodingKind::kDFUDS,
                   false)
    ->Name("convert_louds_dfuds_encoding")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(RunTreeConversion,
                   detail::TreeEncodingKind::kLOUDS,
                   detail::TreeEncodingKind::kDFUDS,
                   true)
    ->Name("convert_louds_dfuds_complete")
    ->Arg(1 << 20)
    ->Arg(1 << 24)
    ->Unit(benchmark::kMillisecond);

}  // namespace pixie::benchmarks
