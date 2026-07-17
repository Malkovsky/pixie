#include <pixie/rmm/implementations.h>

#include "rmm_benchmark_base.h"

int main(int argc, char** argv) {
  pixie_bench::RmMBenchmark<pixie::RmMTree> benchmark;
  return benchmark.Run(argc, argv);
}
