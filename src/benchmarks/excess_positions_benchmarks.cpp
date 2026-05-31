#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/experimental/excess.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

using pixie::experimental::excess_min_128_byte_lut;
using pixie::experimental::excess_min_128_hybrid_lut;
using pixie::experimental::excess_positions_512_branching_lut;
using pixie::experimental::excess_positions_512_byte_lut;
using pixie::experimental::excess_positions_512_expand;
using pixie::experimental::excess_positions_512_expand8;
using pixie::experimental::excess_positions_512_expand_avx512;
using pixie::experimental::excess_positions_512_lut_avx512;
#ifdef PIXIE_AVX2_SUPPORT
using pixie::experimental::excess_min_128_expand16_avx2;
using pixie::experimental::excess_min_128_lane64_sse;
using pixie::experimental::excess_min_128_short_skip;
using pixie::experimental::excess_min_128_split64_sse;
#endif
using pixie::experimental::excess_min_128_nibble_lut;
using pixie::experimental::excess_min_128_scalar_bits;

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

static std::vector<std::array<uint64_t, 2>> make_128_blocks(
    size_t num_blocks = 4096) {
  std::mt19937_64 rng(42);
  std::vector<std::array<uint64_t, 2>> blocks(num_blocks);
  for (auto& b : blocks) {
    b = {rng(), rng()};
  }
  return blocks;
}

static std::vector<std::pair<size_t, size_t>> make_128_ranges(
    size_t num_ranges = 4096) {
  std::mt19937_64 rng(43);
  std::uniform_int_distribution<size_t> offset_dist(0, 128);
  std::vector<std::pair<size_t, size_t>> ranges(num_ranges);
  for (auto& range : ranges) {
    size_t left = offset_dist(rng);
    size_t right = offset_dist(rng);
    if (left > right) {
      std::swap(left, right);
    }
    range = {left, right};
  }
  return ranges;
}

static std::vector<std::pair<size_t, size_t>> make_disjoint_boundary_ranges(
    size_t num_ranges = 4096) {
  std::mt19937_64 rng(45);
  std::uniform_int_distribution<size_t> prefix_dist(0, 126);
  std::vector<std::pair<size_t, size_t>> ranges(num_ranges);
  for (auto& [suffix_left, prefix_right] : ranges) {
    prefix_right = prefix_dist(rng);
    std::uniform_int_distribution<size_t> suffix_dist(prefix_right + 1, 127);
    suffix_left = suffix_dist(rng);
  }
  return ranges;
}

static std::vector<int> make_512_targets(size_t num_targets = 4096) {
  std::mt19937 rng(44);
  std::uniform_int_distribution<int> target_dist(-128, 128);
  std::vector<int> targets(num_targets);
  for (int& target : targets) {
    target = target_dist(rng);
  }
  return targets;
}

static void BM_ExcessMin128(benchmark::State& state) {
  const size_t left = static_cast<size_t>(state.range(0));
  const size_t right = static_cast<size_t>(state.range(1));
  const auto blocks = make_128_blocks();
  const size_t num_blocks = blocks.size();

  size_t idx = 0;
  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    ExcessResult result = excess_min_128(s.data(), left, right);
    benchmark::DoNotOptimize(result.min_excess);
    benchmark::DoNotOptimize(result.offset);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessMin128)
    ->ArgNames({"left", "right"})
    ->Args({0, 128})
    ->Args({0, 127})
    ->Args({0, 16})
    ->Args({0, 32})
    ->Args({0, 48})
    ->Args({0, 64})
    ->Args({0, 31})
    ->Args({1, 17})
    ->Args({3, 35})
    ->Args({5, 37})
    ->Args({32, 64})
    ->Args({33, 65})
    ->Args({64, 96})
    ->Args({61, 93})
    ->Args({96, 128})
    ->Args({56, 72})
    ->Args({60, 68})
    ->Args({63, 64})
    ->Args({17, 17});

static void BM_ExcessMin128BoundaryPairIndependent(benchmark::State& state) {
  const auto blocks = make_128_blocks();
  const auto ranges = make_disjoint_boundary_ranges();
  const size_t num_blocks = blocks.size();
  const size_t num_ranges = ranges.size();

  size_t idx = 0;
  for (auto _ : state) {
    const auto& suffix = blocks[idx % num_blocks];
    const auto& prefix = blocks[(idx + 1) % num_blocks];
    const auto [suffix_left, prefix_right] = ranges[idx % num_ranges];
    ExcessResult suffix_result =
        excess_min_128(suffix.data(), suffix_left, 127);
    ExcessResult prefix_result = excess_min_128(prefix.data(), 0, prefix_right);
    benchmark::DoNotOptimize(suffix_result.min_excess);
    benchmark::DoNotOptimize(suffix_result.offset);
    benchmark::DoNotOptimize(prefix_result.min_excess);
    benchmark::DoNotOptimize(prefix_result.offset);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessMin128BoundaryPairIndependent);

static void BM_ExcessMin128BoundaryPairFused(benchmark::State& state) {
  const auto blocks = make_128_blocks();
  const auto ranges = make_disjoint_boundary_ranges();
  const size_t num_blocks = blocks.size();
  const size_t num_ranges = ranges.size();

  size_t idx = 0;
  for (auto _ : state) {
    const auto& suffix = blocks[idx % num_blocks];
    const auto& prefix = blocks[(idx + 1) % num_blocks];
    const auto [suffix_left, prefix_right] = ranges[idx % num_ranges];
    ExcessBoundaryPairResult result = excess_min_128_disjoint_suffix_prefix(
        suffix.data(), suffix_left, prefix.data(), prefix_right);
    benchmark::DoNotOptimize(result.suffix.min_excess);
    benchmark::DoNotOptimize(result.suffix.offset);
    benchmark::DoNotOptimize(result.prefix.min_excess);
    benchmark::DoNotOptimize(result.prefix.offset);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessMin128BoundaryPairFused);

template <ExcessResult (*Fn)(const uint64_t*, size_t, size_t)>
static void BM_ExcessMin128Variant(benchmark::State& state) {
  const size_t left = static_cast<size_t>(state.range(0));
  const size_t right = static_cast<size_t>(state.range(1));
  const auto blocks = make_128_blocks();
  const size_t num_blocks = blocks.size();

  size_t idx = 0;
  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    ExcessResult result = Fn(s.data(), left, right);
    benchmark::DoNotOptimize(result.min_excess);
    benchmark::DoNotOptimize(result.offset);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

#define PIXIE_BENCH_EXCESS_MIN_VARIANT(name, fn) \
  BENCHMARK_TEMPLATE(BM_ExcessMin128Variant, fn) \
      ->Name(name)                               \
      ->ArgNames({"left", "right"})              \
      ->Args({0, 128})                           \
      ->Args({0, 127})                           \
      ->Args({0, 16})                            \
      ->Args({0, 32})                            \
      ->Args({0, 48})                            \
      ->Args({0, 64})                            \
      ->Args({0, 31})                            \
      ->Args({1, 17})                            \
      ->Args({3, 35})                            \
      ->Args({5, 37})                            \
      ->Args({32, 64})                           \
      ->Args({33, 65})                           \
      ->Args({64, 96})                           \
      ->Args({61, 93})                           \
      ->Args({96, 128})                          \
      ->Args({56, 72})                           \
      ->Args({60, 68})                           \
      ->Args({63, 64})                           \
      ->Args({17, 17})

PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_ScalarBits",
                               excess_min_128_scalar_bits);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_NibbleLUT",
                               excess_min_128_nibble_lut);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_ByteLUT",
                               excess_min_128_byte_lut);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_HybridLUT",
                               excess_min_128_hybrid_lut);
#ifdef PIXIE_AVX2_SUPPORT
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_Expand16AVX2",
                               excess_min_128_expand16_avx2);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_Lane64SSE",
                               excess_min_128_lane64_sse);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_Split64SSE",
                               excess_min_128_split64_sse);
PIXIE_BENCH_EXCESS_MIN_VARIANT("BM_ExcessMin128_ShortSkip",
                               excess_min_128_short_skip);
#endif

#undef PIXIE_BENCH_EXCESS_MIN_VARIANT

template <ExcessResult (*Fn)(const uint64_t*, size_t, size_t)>
static void BM_ExcessMin128RandomRange(benchmark::State& state) {
  const auto blocks = make_128_blocks();
  const auto ranges = make_128_ranges();
  const size_t num_blocks = blocks.size();
  const size_t num_ranges = ranges.size();

  size_t idx = 0;
  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    const auto [left, right] = ranges[idx % num_ranges];
    ExcessResult result = Fn(s.data(), left, right);
    benchmark::DoNotOptimize(result.min_excess);
    benchmark::DoNotOptimize(result.offset);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

#define PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT(name, fn) \
  BENCHMARK_TEMPLATE(BM_ExcessMin128RandomRange, fn)->Name(name)

PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_RandomRange",
                                      excess_min_128);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_ScalarBits_RandomRange",
                                      excess_min_128_scalar_bits);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_NibbleLUT_RandomRange",
                                      excess_min_128_nibble_lut);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_ByteLUT_RandomRange",
                                      excess_min_128_byte_lut);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_HybridLUT_RandomRange",
                                      excess_min_128_hybrid_lut);
#ifdef PIXIE_AVX2_SUPPORT
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT(
    "BM_ExcessMin128_Expand16AVX2_RandomRange",
    excess_min_128_expand16_avx2);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_Lane64SSE_RandomRange",
                                      excess_min_128_lane64_sse);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_Split64SSE_RandomRange",
                                      excess_min_128_split64_sse);
PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT("BM_ExcessMin128_ShortSkip_RandomRange",
                                      excess_min_128_short_skip);
#endif

#undef PIXIE_BENCH_EXCESS_MIN_RANDOM_VARIANT

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

static void excess_positions_512_scalar_benchmark(const uint64_t* s,
                                                  int target_x,
                                                  uint64_t* out) noexcept {
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
}

static void BM_ExcessPositions512_Scalar(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_scalar_benchmark(s.data(), target_x, out);
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

static void BM_ExcessPositions512_ByteLUT(benchmark::State& state) {
  const int target_x = state.range(0);
  const auto blocks = make_blocks();
  const size_t num_blocks = blocks.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    excess_positions_512_byte_lut(s.data(), target_x, out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_ExcessPositions512_ByteLUT)
    ->ArgNames({"X"})
    ->Args({-64})
    ->Args({-8})
    ->Args({0})
    ->Args({8})
    ->Args({64});

template <void (*Fn)(const uint64_t*, int, uint64_t*)>
static void BM_ExcessPositions512RandomTarget(benchmark::State& state) {
  const auto blocks = make_blocks();
  const auto targets = make_512_targets();
  const size_t num_blocks = blocks.size();
  const size_t num_targets = targets.size();

  alignas(64) uint64_t out[8];
  size_t idx = 0;

  for (auto _ : state) {
    const auto& s = blocks[idx % num_blocks];
    Fn(s.data(), targets[idx % num_targets], out);
    benchmark::DoNotOptimize(out);
    ++idx;
  }

  state.SetItemsProcessed(state.iterations());
}

#define PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(name, fn) \
  BENCHMARK_TEMPLATE(BM_ExcessPositions512RandomTarget, fn)->Name(name)

PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM("BM_ExcessPositions512_RandomTarget",
                                        excess_positions_512);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_BranchingLUT_RandomTarget",
    excess_positions_512_branching_lut);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_LUTAVX512_RandomTarget",
    excess_positions_512_lut_avx512);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_Expand_RandomTarget",
    excess_positions_512_expand);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_Expand8_RandomTarget",
    excess_positions_512_expand8);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_ExpandAVX512_RandomTarget",
    excess_positions_512_expand_avx512);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_Scalar_RandomTarget",
    excess_positions_512_scalar_benchmark);
PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM(
    "BM_ExcessPositions512_ByteLUT_RandomTarget",
    excess_positions_512_byte_lut);

#undef PIXIE_BENCH_EXCESS_POSITIONS_512_RANDOM
