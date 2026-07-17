#include <pixie/rmm/implementations.h>

#include <string_view>

#include "rmm_benchmark_base.h"

namespace pixie_bench {

template <>
struct RmMBenchmarkTraits<pixie::SdslRmMTree> {
  static constexpr std::size_t DefaultBlockBits = 64;

  static bool SupportsOp(std::string_view op) {
    return op == "rank1" || op == "rank0" || op == "select1" ||
           op == "excess" || op == "fwdsearch" || op == "bwdsearch" ||
           op == "range_min_query_pos" || op == "range_min_query_val" ||
           op == "close" || op == "open" || op == "enclose";
  }
};

}  // namespace pixie_bench

int main(int argc, char** argv) {
  pixie_bench::RmMBenchmark<pixie::SdslRmMTree> benchmark;
  return benchmark.Run(argc, argv);
}
