#include <benchmark/benchmark.h>
#include <pixie/bits.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

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


namespace {
// clang-format off
static inline const __m256i excess_lut_em4 = _mm256_setr_epi8(
    0x08, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x08, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_lut_em3 = _mm256_setr_epi8(
    0x04, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_lut_em2 = _mm256_setr_epi8(
    0x02, 0x08, 0x08, 0x00,
    0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00,
    0x02, 0x00, 0x00, 0x00,
    0x02, 0x08, 0x08, 0x00,
    0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00,
    0x02, 0x00, 0x00, 0x00);

static inline const __m256i excess_lut_em1 = _mm256_setr_epi8(
    0x01, 0x04, 0x05, 0x00,
    0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00,
    0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00,
    0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00,
    0x05, 0x00, 0x01, 0x00);

static inline const __m256i excess_lut_e0 = _mm256_setr_epi8(
    0x00, 0x02, 0x02, 0x08,
    0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00,
    0x08, 0x02, 0x02, 0x00,
    0x00, 0x02, 0x02, 0x08,
    0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00,
    0x08, 0x02, 0x02, 0x00);

static inline const __m256i excess_lut_e1 = _mm256_setr_epi8(
    0x00, 0x01, 0x00, 0x05,
    0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05,
    0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05,
    0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05,
    0x00, 0x05, 0x04, 0x01);

static inline const __m256i excess_lut_e2 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x02,
    0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A,
    0x00, 0x08, 0x08, 0x02,
    0x00, 0x00, 0x00, 0x02,
    0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A,
    0x00, 0x08, 0x08, 0x02);

static inline const __m256i excess_lut_e3 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x04);

static inline const __m256i excess_lut_e4 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x08,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x08);

// clang-format on

}


#ifdef PIXIE_AVX2_SUPPORT
static inline void excess_positions_512_branching_lut(const uint64_t* s, int target_x, uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  int cur = 0;
  // Use a local copy since pixie::excess_lut_delta seems inaccessible or causes issues
  const __m256i vdelta = _mm256_setr_epi8(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4,
                                          -4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4);
  const __m256i vmult = _mm256_set1_epi16(0x1001);

  for (int k = 0; k < 4; ++k) {
    __m128i word_vec = _mm_loadu_si128((const __m128i*)&s[2 * k]);
    __m128i lo_nibbles = _mm_and_si128(word_vec, _mm_set1_epi8(0x0F));
    __m128i hi_nibbles =
        _mm_and_si128(_mm_srli_epi16(word_vec, 4), _mm_set1_epi8(0x0F));

    __m128i unpack_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
    __m128i unpack_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);

    __m256i nibbles = _mm256_inserti128_si256(_mm256_castsi128_si256(unpack_lo),
                                              unpack_hi, 1);

    __m256i ps = _mm256_shuffle_epi8(vdelta, nibbles);
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 1));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 2));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 4));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 8));

    __m128i ps_lo = _mm256_castsi256_si128(ps);
    __m128i ps_hi = _mm256_extracti128_si256(ps, 1);
    __m128i carry = _mm_set1_epi8((int8_t)_mm_extract_epi8(ps_lo, 15));
    ps_hi = _mm_add_epi8(ps_hi, carry);
    ps = _mm256_inserti128_si256(_mm256_castsi128_si256(ps_lo), ps_hi, 1);

    __m256i b = _mm256_permute2x128_si256(ps, ps, 0x08);
    __m256i excl_ps = _mm256_alignr_epi8(ps, b, 15);

    int target_rel = target_x - cur;
    int block_delta =
        2 * (std::popcount(s[2 * k]) + std::popcount(s[2 * k + 1])) - 128;

    const int d = 2 * target_rel - block_delta;
    if (d < -128 || d > 128) {
      cur += block_delta;
      continue;
    }

    if (target_rel == 128 || target_rel == -128) {
      out[2 * k + 1] |= (1ULL << 63);
      cur += block_delta;
      continue;
    }

    __m256i vtgt = _mm256_set1_epi8((int8_t)target_rel);
    __m256i t = _mm256_sub_epi8(vtgt, excl_ps);

    __m256i total_match = _mm256_setzero_si256();
    __m256i t_eq;
    
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-4));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_em4, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-3));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_em3, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-2));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_em2, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-1));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_em1, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(0));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_e0, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(1));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_e1, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(2));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_e2, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(3));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_e3, nibbles)));

    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(4));
    total_match = _mm256_or_si256(total_match, _mm256_and_si256(t_eq, _mm256_shuffle_epi8(excess_lut_e4, nibbles)));

    __m256i res = _mm256_maddubs_epi16(total_match, vmult);
    __m128i res_lo = _mm256_castsi256_si128(res);
    __m128i res_hi = _mm256_extracti128_si256(res, 1);
    __m128i packed = _mm_packus_epi16(res_lo, res_hi);

    _mm_storeu_si128((__m128i*)&out[2 * k], packed);

    cur += block_delta;
  }
}

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
#endif

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


