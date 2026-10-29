#include <benchmark/benchmark.h>

#include "micro256.h"
#include "src/benchmarks/tree_input.h"

template <unsigned Mode>
void kernel(benchmark::State& state) {
  const size_t n = state.range(0);
  const unsigned shape = state.range(1);
  std::vector<std::vector<uint64_t>> pool;
  for (unsigned seed = 0; seed < 256; ++seed) {
    auto words = pixie::benchmarks::GenerateRandomBp(n, seed + 42);
    if (shape) {
      std::fill(words.begin(), words.end(), 0);
      for (size_t k = 0; k < n; ++k) {
        size_t p = shape == 1 ? k : (k ? 2 * k - 1 : 0);
        words[p / 64] |= uint64_t{1} << (p % 64);
      }
    }
    pool.push_back(std::move(words));
  }
  auto observer = [](auto node) { benchmark::DoNotOptimize(node); };
  using Compiled =
      pixie::experimental::detail::BpMicroSubtree<256, 1, true, true>;
  using Block =
      std::conditional_t<Mode == 0, Compiled, micro256::Block<Mode - 1>>;
  Block block;
  size_t cursor = 0;
  for (auto _ : state) {
    auto& words = pool[cursor++ & 255];
    if (!block.load(words, 2 * n, {0, 0})) {
      state.SkipWithError("complete microtree rejected");
      break;
    }
    block.walk(observer, observer);
  }
  state.SetItemsProcessed(state.iterations() * n);
}

int main(int argc, char** argv) {
  auto add = [](const char* name, auto function) {
    auto* b = benchmark::RegisterBenchmark(name, function);
    for (int n : {16, 32, 64, 128}) {
      for (int shape : {0, 1, 2}) {
        b->Args({n, shape});
      }
    }
  };
  add("compiled", kernel<0>);
  add("stack", kernel<1>);
  add("simd_matching", kernel<2>);
  add("simd_direct", kernel<3>);
  add("simd_fused", kernel<4>);
  add("simd_depth_cache", kernel<5>);
  add("simd_levels", kernel<6>);
  add("simd_levels_hybrid", kernel<7>);
  add("simd_levels_register", kernel<8>);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}
