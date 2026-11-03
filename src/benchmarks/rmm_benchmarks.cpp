#include <pixie/rank_select/support.h>
#include <pixie/rmm/implementations.h>

#include "rmm_benchmark_base.h"

namespace pixie_bench {

template <class RankSelect>
struct RmMBenchmarkTraits<pixie::RmMBTree<4, 32, RankSelect>> {
  static constexpr std::size_t DefaultBlockBits = pixie::RmMBTree<>::kBlockBits;

  static bool SupportsOp(std::string_view) { return true; }
};

}  // namespace pixie_bench

int main(int argc, char** argv) {
  // Preserve our custom flag across Google Benchmark's possible re-exec.
  benchmark::MaybeReenterWithoutASLR(argc, argv);
  // Keep backend comparisons in one executable with identical fixture options.
  bool baseline = false;
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string_view option(argv[i]);
    constexpr std::string_view prefix = "--rank_select_backend=";
    if (option.starts_with(prefix)) {
      const auto backend = option.substr(prefix.size());
      if (backend != "baseline" && backend != "compressed") {
        std::cerr << "Expected baseline or compressed rank/select backend\n";
        return 1;
      }
      baseline = backend == "baseline";
    } else {
      argv[kept++] = argv[i];
    }
  }
  argc = kept;
  argv[argc] = nullptr;
  benchmark::AddCustomContext("rank_select_backend",
                              baseline ? "baseline" : "compressed");
  if (baseline) {
    using Tree = pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>;
    pixie_bench::RmMBenchmark<Tree> benchmark;
    return benchmark.Run(argc, argv);
  }
  pixie_bench::RmMBenchmark<pixie::RmMBTree<>> benchmark;
  return benchmark.Run(argc, argv);
}
