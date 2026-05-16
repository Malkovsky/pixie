#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/experimental/excess.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

using pixie::experimental::excess_positions_512_branching_lut;
using pixie::experimental::excess_positions_512_expand;
using pixie::experimental::excess_positions_512_expand8;
using pixie::experimental::excess_positions_512_expand_avx512;
using pixie::experimental::excess_positions_512_lut_avx512;

static std::vector<std::array<uint64_t, 8>> make_blocks(
    size_t num_blocks = 4096) {
  std::mt19937_64 rng(42);
  std::vector<std::array<uint64_t, 8>> blocks(num_blocks);
  for (auto& b : blocks) {
    for (auto& w : b) {
      w = rng();
    }
  }
  return blocks;
}

static void BM_ExcessPositions512(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_BranchingLUT(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_branching_lut(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_BranchingLUT)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_LUTAVX512(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_lut_avx512(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_LUTAVX512)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_Expand(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_expand(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_Expand)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_Expand8(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_expand8(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_Expand8)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_ExpandAVX512(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_expand_avx512(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_ExpandAVX512)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

static void BM_ExcessPositions512_Scalar(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    for (int w = 0; w < 8; ++w) {
      out[w] = 0;
    }
    int cur = 0;
    for (size_t i = 0; i < 512; ++i) {
      const int bit = int((s[i >> 6] >> (i & 63)) & 1ull);
      cur += bit ? +1 : -1;
      if (cur == target_x) {
        out[i >> 6] |= (uint64_t{1} << (i & 63));
      }
    }
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_Scalar)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});
