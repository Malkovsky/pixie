#include <pixie/rmm/implementations.h>

#include "rmm_benchmark_base.h"

namespace pixie_bench {

template <>
struct RmMBenchmarkTraits<pixie::experimental::RmMBTree<>> {
  static constexpr std::size_t DefaultBlockBits =
      pixie::experimental::RmMBTree<>::kBlockBits;

  static bool SupportsOp(std::string_view) { return true; }
};

}  // namespace pixie_bench

int main(int argc, char** argv) {
  pixie_bench::RmMBenchmark<pixie::experimental::RmMBTree<>> benchmark;
  return benchmark.Run(argc, argv);
}
