#include <gtest/gtest.h>
#include <pixie/bits.h>
#include <pixie/experimental/excess.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <random>
#include <utility>

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

static void naive_excess_positions_512(const uint64_t* s,
                                       int target_x,
                                       uint64_t* out) {
  for (int w = 0; w < 8; ++w) {
    out[w] = 0;
  }
  if (target_x < -512 || target_x > 512) {
    return;
  }
  int cur = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
}

static int naive_excess_positions_128(const uint64_t* s,
                                      int target_x,
                                      uint64_t* out) {
  out[0] = out[1] = 0;
  const int block_excess =
      2 * (std::popcount(s[0]) + std::popcount(s[1])) - 128;
  if (target_x < -128 || target_x > 128) {
    return block_excess;
  }
  int cur = 0;
  for (size_t i = 0; i < 128; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
  return block_excess;
}

static int naive_prefix_excess_128(const uint64_t* s, size_t end_offset) {
  end_offset = std::min<size_t>(end_offset, 128);
  int cur = 0;
  for (size_t i = 0; i < end_offset; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
  }
  return cur;
}

static ExcessResult naive_excess_min_128(const uint64_t* s,
                                         size_t left,
                                         size_t right) {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int cur = 0;
  int best = 0;
  size_t best_offset = 0;
  bool found = false;
  if (left == 0) {
    found = true;
  }
  for (size_t bit = 0; bit < right; ++bit) {
    cur += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
    const size_t offset = bit + 1;
    if (offset < left) {
      continue;
    }
    if (!found || cur < best) {
      best = cur;
      best_offset = offset;
      found = true;
    }
  }
  if (!found) {
    best = naive_prefix_excess_128(s, left);
    best_offset = left;
  }
  return {best, best_offset};
}

static size_t naive_forward_search_128(const uint64_t* s,
                                       int target_x,
                                       size_t start_offset) {
  if (start_offset >= 128) {
    return 128;
  }
  int cur = 0;
  for (size_t i = 0; i < 128; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (i >= start_offset && cur == target_x) {
      return i;
    }
  }
  return 128;
}

static size_t naive_backward_search_128(const uint64_t* s,
                                        int target_x,
                                        size_t end_offset) {
  if (end_offset == 0) {
    return 128;
  }
  const size_t max_prefix_length = end_offset - 1;
  for (size_t prefix_length = max_prefix_length; prefix_length > 0;
       --prefix_length) {
    int cur = 0;
    for (size_t i = 0; i < prefix_length; ++i) {
      const uint64_t w = s[i >> 6];
      const int bit = int((w >> (i & 63)) & 1ull);
      cur += bit ? +1 : -1;
    }
    if (cur == target_x) {
      return prefix_length;
    }
  }
  return target_x == 0 ? 0 : 128;
}

static size_t count_matches(const uint64_t* out) {
  size_t cnt = 0;
  for (int w = 0; w < 8; ++w) {
    cnt += std::popcount(out[w]);
  }
  return cnt;
}

static void check_boundary_pair_matches_independent(
    const std::array<uint64_t, 2>& suffix,
    size_t suffix_left,
    const std::array<uint64_t, 2>& prefix,
    size_t prefix_right) {
  const ExcessBoundaryPairResult result = excess_min_128_disjoint_suffix_prefix(
      suffix.data(), suffix_left, prefix.data(), prefix_right);
  const ExcessResult expected_suffix =
      excess_min_128(suffix.data(), suffix_left, 127);
  const ExcessResult expected_prefix =
      excess_min_128(prefix.data(), 0, prefix_right);

  ASSERT_EQ(result.suffix.min_excess, expected_suffix.min_excess)
      << "suffix_left=" << suffix_left << " prefix_right=" << prefix_right;
  ASSERT_EQ(result.suffix.offset, expected_suffix.offset)
      << "suffix_left=" << suffix_left << " prefix_right=" << prefix_right;
  ASSERT_EQ(result.prefix.min_excess, expected_prefix.min_excess)
      << "suffix_left=" << suffix_left << " prefix_right=" << prefix_right;
  ASSERT_EQ(result.prefix.offset, expected_prefix.offset)
      << "suffix_left=" << suffix_left << " prefix_right=" << prefix_right;
}

template <typename Fn>
static void check_matches_naive(Fn fn,
                                const char* fn_name,
                                const uint64_t* s,
                                int target_x,
                                int case_id = 0) {
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];
  fn(s, target_x, out);
  naive_excess_positions_512(s, target_x, ref);
  for (int w = 0; w < 8; ++w) {
    ASSERT_EQ(out[w], ref[w])
        << fn_name << " case=" << case_id << " x=" << target_x << " word=" << w;
  }
  ASSERT_EQ(count_matches(out), count_matches(ref))
      << fn_name << " case=" << case_id << " x=" << target_x;
}

template <typename Fn>
static void check_min_matches_naive(Fn fn,
                                    const char* fn_name,
                                    const uint64_t* s,
                                    size_t left,
                                    size_t right,
                                    int case_id = 0) {
  const ExcessResult result = fn(s, left, right);
  const ExcessResult expected = naive_excess_min_128(s, left, right);
  ASSERT_EQ(result.min_excess, expected.min_excess)
      << fn_name << " case=" << case_id << " left=" << left
      << " right=" << right;
  ASSERT_EQ(result.offset, expected.offset)
      << fn_name << " case=" << case_id << " left=" << left
      << " right=" << right;
}

TEST(ExcessPositions128, MatchesNaiveMasksAndDelta) {
  const std::array<std::array<uint64_t, 2>, 4> cases = {{
      {0, 0},
      {UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
  }};

  for (const auto& s : cases) {
    for (int x = -130; x <= 130; ++x) {
      uint64_t out[2];
      uint64_t ref[2];
      const int delta = excess_positions_128(s.data(), x, out);
      const int ref_delta = naive_excess_positions_128(s.data(), x, ref);
      EXPECT_EQ(delta, ref_delta) << "x=" << x;
      EXPECT_EQ(out[0], ref[0]) << "x=" << x;
      EXPECT_EQ(out[1], ref[1]) << "x=" << x;
    }
  }
}

TEST(ExcessPositions128, PrefixExcessMatchesNaive) {
  std::mt19937_64 rng(42);
  const std::array<size_t, 13> offsets = {0,  1,  2,  31,  32,  63, 64,
                                          65, 95, 96, 127, 128, 129};

  for (int t = 0; t < 1000; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    for (size_t offset : offsets) {
      EXPECT_EQ(prefix_excess_128(s.data(), offset),
                naive_prefix_excess_128(s.data(), offset))
          << "case=" << t << " offset=" << offset;
    }
  }
}

TEST(ExcessPositions128, MinMatchesNaiveFixedCases) {
  const std::array<std::array<uint64_t, 2>, 5> cases = {{
      {0, 0},
      {UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
      {0x0000FFFF0000FFFFull, 0xFFFF0000FFFF0000ull},
  }};
  const std::array<std::pair<size_t, size_t>, 12> ranges = {{
      {0, 128},
      {0, 0},
      {1, 1},
      {63, 65},
      {64, 64},
      {64, 128},
      {3, 6},
      {5, 5},
      {127, 128},
      {128, 128},
      {120, 127},
      {129, 140},
  }};

  for (const auto& s : cases) {
    for (const auto [left, right] : ranges) {
      const ExcessResult result = excess_min_128(s.data(), left, right);
      const ExcessResult expected = naive_excess_min_128(s.data(), left, right);
      EXPECT_EQ(result.min_excess, expected.min_excess)
          << "left=" << left << " right=" << right;
      EXPECT_EQ(result.offset, expected.offset)
          << "left=" << left << " right=" << right;
    }
  }
}

TEST(ExcessPositions128, MinReturnsFirstTie) {
  const std::array<uint64_t, 2> s = {0x5555555555555555ull,
                                     0x5555555555555555ull};
  const ExcessResult result = excess_min_128(s.data(), 0, 128);
  EXPECT_EQ(result.min_excess, 0);
  EXPECT_EQ(result.offset, 0u);

  const ExcessResult shifted = excess_min_128(s.data(), 1, 128);
  EXPECT_EQ(shifted.min_excess, 0);
  EXPECT_EQ(shifted.offset, 2u);

  const ExcessResult left_tie = excess_min_128(s.data(), 2, 128);
  EXPECT_EQ(left_tie.min_excess, 0);
  EXPECT_EQ(left_tie.offset, 2u);
}

TEST(ExcessPositions128, MinHandlesRightBoundary) {
  const std::array<uint64_t, 2> s = {0, 0};

  const ExcessResult without_last = excess_min_128(s.data(), 0, 127);
  EXPECT_EQ(without_last.min_excess, -127);
  EXPECT_EQ(without_last.offset, 127u);

  const ExcessResult with_last = excess_min_128(s.data(), 0, 128);
  EXPECT_EQ(with_last.min_excess, -128);
  EXPECT_EQ(with_last.offset, 128u);
}

TEST(ExcessPositions128, MinPartialNibbleBoundsExcludeOuterMin) {
  const std::array<uint64_t, 2> s = {0, 0};

  const ExcessResult short_prefix = excess_min_128(s.data(), 1, 2);
  EXPECT_EQ(short_prefix.min_excess, -2);
  EXPECT_EQ(short_prefix.offset, 2u);

  const ExcessResult short_suffix = excess_min_128(s.data(), 2, 3);
  EXPECT_EQ(short_suffix.min_excess, -3);
  EXPECT_EQ(short_suffix.offset, 3u);
}

TEST(ExcessPositions128, MinPositiveRangeKeepsLeftBoundary) {
  const std::array<uint64_t, 2> s = {UINT64_MAX, UINT64_MAX};

  const ExcessResult result = excess_min_128(s.data(), 64, 128);
  EXPECT_EQ(result.min_excess, 64);
  EXPECT_EQ(result.offset, 64u);
}

TEST(ExcessPositions128, MinInvalidRangeUsesSentinel) {
  const std::array<uint64_t, 2> s = {0, 0};
  const ExcessResult result = excess_min_128(s.data(), 17, 16);
  EXPECT_EQ(result.min_excess, 0);
  EXPECT_EQ(result.offset, 128u);
}

TEST(ExcessPositions128, MinMatchesNaiveRandom) {
  std::mt19937_64 rng(43);
  std::uniform_int_distribution<size_t> offset_dist(0, 128);

  for (int t = 0; t < 1000; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    for (int q = 0; q < 32; ++q) {
      size_t left = offset_dist(rng);
      size_t right = offset_dist(rng);
      if (left > right) {
        std::swap(left, right);
      }
      const ExcessResult result = excess_min_128(s.data(), left, right);
      const ExcessResult expected = naive_excess_min_128(s.data(), left, right);
      ASSERT_EQ(result.min_excess, expected.min_excess)
          << "case=" << t << " left=" << left << " right=" << right;
      ASSERT_EQ(result.offset, expected.offset)
          << "case=" << t << " left=" << left << " right=" << right;
    }
  }
}

TEST(ExcessPositions128, DisjointBoundaryPairMatchesIndependentFixedCases) {
  const std::array<std::array<uint64_t, 2>, 5> cases = {{
      {0, 0},
      {UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
      {0x0000FFFF0000FFFFull, 0xFFFF0000FFFF0000ull},
  }};
  const std::array<std::pair<size_t, size_t>, 10> ranges = {{
      {1, 0},
      {4, 3},
      {17, 16},
      {32, 31},
      {63, 62},
      {64, 32},
      {65, 63},
      {96, 31},
      {127, 0},
      {120, 119},
  }};

  for (const auto& suffix : cases) {
    for (const auto& prefix : cases) {
      for (const auto [suffix_left, prefix_right] : ranges) {
        check_boundary_pair_matches_independent(suffix, suffix_left, prefix,
                                                prefix_right);
      }
    }
  }
}

TEST(ExcessPositions128, DisjointBoundaryPairMatchesIndependentRandom) {
  std::mt19937_64 rng(45);
  std::uniform_int_distribution<size_t> prefix_dist(0, 126);

  for (int t = 0; t < 1000; ++t) {
    const std::array<uint64_t, 2> suffix = {rng(), rng()};
    const std::array<uint64_t, 2> prefix = {rng(), rng()};
    for (int q = 0; q < 16; ++q) {
      const size_t prefix_right = prefix_dist(rng);
      std::uniform_int_distribution<size_t> suffix_dist(prefix_right + 1, 127);
      const size_t suffix_left = suffix_dist(rng);
      check_boundary_pair_matches_independent(suffix, suffix_left, prefix,
                                              prefix_right);
    }
  }
}

TEST(ExcessPositions128, BoundaryPairFallbackMatchesIndependent) {
  const std::array<uint64_t, 2> suffix = {0x0123456789ABCDEFull,
                                          0xFEDCBA9876543210ull};
  const std::array<uint64_t, 2> prefix = {0x0000FFFF0000FFFFull,
                                          0xFFFF0000FFFF0000ull};

  check_boundary_pair_matches_independent(suffix, 32, prefix, 32);
  check_boundary_pair_matches_independent(suffix, 0, prefix, 127);
  check_boundary_pair_matches_independent(suffix, 128, prefix, 0);
}

TEST(ExcessPositions128Experimental, MinVariantsMatchNaive) {
  const std::array<std::array<uint64_t, 2>, 6> cases = {{
      {0, 0},
      {UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
      {0x0000FFFF0000FFFFull, 0xFFFF0000FFFF0000ull},
      {0x1111222233334444ull, 0x8888777766665555ull},
  }};
  const std::array<std::pair<size_t, size_t>, 25> ranges = {{
      {0, 128},   {0, 127},   {0, 16},    {0, 31},    {0, 32},
      {0, 48},    {0, 64},    {32, 64},   {64, 96},   {96, 128},
      {56, 72},   {60, 68},   {63, 64},   {17, 17},   {1, 2},
      {2, 3},     {3, 6},     {5, 5},     {64, 64},   {64, 128},
      {127, 128}, {128, 128}, {120, 127}, {129, 140}, {17, 16},
  }};

  int case_id = 0;
  for (const auto& s : cases) {
    for (const auto [left, right] : ranges) {
      check_min_matches_naive(excess_min_128_scalar_bits, "scalar_bits",
                              s.data(), left, right, case_id);
      check_min_matches_naive(excess_min_128_nibble_lut, "nibble_lut", s.data(),
                              left, right, case_id);
      check_min_matches_naive(excess_min_128_byte_lut, "byte_lut", s.data(),
                              left, right, case_id);
      check_min_matches_naive(excess_min_128_hybrid_lut, "hybrid_lut", s.data(),
                              left, right, case_id);
#ifdef PIXIE_AVX2_SUPPORT
      check_min_matches_naive(excess_min_128_expand16_avx2, "expand16_avx2",
                              s.data(), left, right, case_id);
      check_min_matches_naive(excess_min_128_lane64_sse, "lane64_sse", s.data(),
                              left, right, case_id);
      check_min_matches_naive(excess_min_128_split64_sse, "split64_sse",
                              s.data(), left, right, case_id);
      check_min_matches_naive(excess_min_128_short_skip, "short_skip", s.data(),
                              left, right, case_id);
#endif
      ++case_id;
    }
  }
}

TEST(ExcessPositions128Experimental, MinVariantsMatchNaiveRandom) {
  std::mt19937_64 rng(44);
  std::uniform_int_distribution<size_t> offset_dist(0, 140);

  for (int t = 0; t < 500; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    for (int q = 0; q < 16; ++q) {
      const size_t left = offset_dist(rng);
      const size_t right = offset_dist(rng);
      check_min_matches_naive(excess_min_128_scalar_bits, "scalar_bits",
                              s.data(), left, right, t);
      check_min_matches_naive(excess_min_128_nibble_lut, "nibble_lut", s.data(),
                              left, right, t);
      check_min_matches_naive(excess_min_128_byte_lut, "byte_lut", s.data(),
                              left, right, t);
      check_min_matches_naive(excess_min_128_hybrid_lut, "hybrid_lut", s.data(),
                              left, right, t);
#ifdef PIXIE_AVX2_SUPPORT
      check_min_matches_naive(excess_min_128_expand16_avx2, "expand16_avx2",
                              s.data(), left, right, t);
      check_min_matches_naive(excess_min_128_lane64_sse, "lane64_sse", s.data(),
                              left, right, t);
      check_min_matches_naive(excess_min_128_split64_sse, "split64_sse",
                              s.data(), left, right, t);
      check_min_matches_naive(excess_min_128_short_skip, "short_skip", s.data(),
                              left, right, t);
#endif
    }
  }
}

TEST(ExcessPositions128, ForwardAndBackwardSearchMatchNaive) {
  std::mt19937_64 rng(42);
  const std::array<size_t, 8> offsets = {0, 1, 63, 64, 65, 126, 127, 128};

  for (int t = 0; t < 1000; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    for (int x = -128; x <= 128; x += 7) {
      for (size_t offset : offsets) {
        int block_excess = 0;
        EXPECT_EQ(forward_search_128(s.data(), x, offset, &block_excess),
                  naive_forward_search_128(s.data(), x, offset))
            << "case=" << t << " x=" << x << " offset=" << offset;
        EXPECT_EQ(block_excess,
                  2 * (std::popcount(s[0]) + std::popcount(s[1])) - 128)
            << "case=" << t;

        EXPECT_EQ(backward_search_128(s.data(), x, offset),
                  naive_backward_search_128(s.data(), x, offset))
            << "case=" << t << " x=" << x << " offset=" << offset;
      }
    }
  }
}

TEST(ExcessPositions512, AllZeros) {
  alignas(64) uint64_t s[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (int x = -8; x <= 0; ++x) {
    excess_positions_512(s, x, out);
    naive_excess_positions_512(s, x, ref);
    for (int w = 0; w < 8; ++w) {
      EXPECT_EQ(out[w], ref[w]) << "x=" << x << " word=" << w;
    }
  }

  excess_positions_512(s, 1, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], 0u);
  }

  for (int x = -8; x <= 1; ++x) {
    check_matches_naive(excess_positions_512_expand, "expand", s, x);
    check_matches_naive(excess_positions_512_expand8, "expand8", s, x);
    check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512", s,
                        x);
    check_matches_naive(excess_positions_512_branching_lut, "branching_lut", s,
                        x);
    check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, x);
    check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, x);
  }
}

TEST(ExcessPositions512, AllOnes) {
  alignas(64) uint64_t s[8] = {UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX,
                               UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX};
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (int x = 0; x <= 8; ++x) {
    excess_positions_512(s, x, out);
    naive_excess_positions_512(s, x, ref);
    for (int w = 0; w < 8; ++w) {
      EXPECT_EQ(out[w], ref[w]) << "x=" << x << " word=" << w;
    }
  }

  excess_positions_512(s, -1, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], 0u);
  }

  for (int x = -1; x <= 8; ++x) {
    check_matches_naive(excess_positions_512_expand, "expand", s, x);
    check_matches_naive(excess_positions_512_expand8, "expand8", s, x);
    check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512", s,
                        x);
    check_matches_naive(excess_positions_512_branching_lut, "branching_lut", s,
                        x);
    check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, x);
    check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, x);
  }
}

TEST(ExcessPositions512, Alternating) {
  alignas(64) uint64_t s[8];
  for (int w = 0; w < 8; ++w) {
    s[w] = 0xAAAAAAAAAAAAAAAAull;
  }
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (int x = -2; x <= 2; ++x) {
    excess_positions_512(s, x, out);
    naive_excess_positions_512(s, x, ref);
    for (int w = 0; w < 8; ++w) {
      EXPECT_EQ(out[w], ref[w]) << "x=" << x << " word=" << w;
    }
    check_matches_naive(excess_positions_512_expand, "expand", s, x);
    check_matches_naive(excess_positions_512_expand8, "expand8", s, x);
    check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512", s,
                        x);
    check_matches_naive(excess_positions_512_branching_lut, "branching_lut", s,
                        x);
    check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, x);
    check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, x);
  }
}

TEST(ExcessPositions512, OutOfRange) {
  alignas(64) uint64_t s[8] = {UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX,
                               UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX};
  alignas(64) uint64_t out[8];
  excess_positions_512(s, 513, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], 0u);
  }
  excess_positions_512(s, -513, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], 0u);
  }
}

TEST(ExcessPositions512, ExhaustiveSmall16) {
  alignas(64) uint64_t s[8];
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (uint64_t pattern = 0; pattern < (1ull << 16); ++pattern) {
    s[0] = pattern;
    for (int w = 1; w < 8; ++w) {
      s[w] = 0;
    }
    for (int x = -20; x <= 20; ++x) {
      excess_positions_512(s, x, out);
      naive_excess_positions_512(s, x, ref);
      check_matches_naive(excess_positions_512_expand, "expand", s, x,
                          static_cast<int>(pattern));
      check_matches_naive(excess_positions_512_expand8, "expand8", s, x,
                          static_cast<int>(pattern));
      check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512",
                          s, x, static_cast<int>(pattern));
      check_matches_naive(excess_positions_512_branching_lut, "branching_lut",
                          s, x, static_cast<int>(pattern));
      check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, x,
                          static_cast<int>(pattern));
      check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, x,
                          static_cast<int>(pattern));
      for (int w = 0; w < 8; ++w) {
        ASSERT_EQ(out[w], ref[w])
            << "pattern=" << pattern << " x=" << x << " word=" << w;
      }
    }
  }
}

TEST(ExcessPositions512, Random) {
  const int cases = [] {
    const char* env = std::getenv("EXCESS_POS_CASES");
    return env ? std::atoi(env) : 1000;
  }();
  const uint64_t seed = [] {
    const char* env = std::getenv("EXCESS_POS_SEED");
    return env ? std::stoull(env) : 42ull;
  }();

  std::mt19937_64 rng(static_cast<uint64_t>(seed));
  std::uniform_int_distribution<int> x_dist(-512, 512);

  alignas(64) uint64_t s[8];
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (int t = 0; t < cases; ++t) {
    for (int w = 0; w < 8; ++w) {
      s[w] = rng();
    }
    const int x = x_dist(rng);

    excess_positions_512(s, x, out);
    naive_excess_positions_512(s, x, ref);
    check_matches_naive(excess_positions_512_expand, "expand", s, x, t);
    check_matches_naive(excess_positions_512_expand8, "expand8", s, x, t);
    check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512", s,
                        x, t);
    check_matches_naive(excess_positions_512_branching_lut, "branching_lut", s,
                        x, t);
    check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, x, t);
    check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, x, t);

    for (int w = 0; w < 8; ++w) {
      ASSERT_EQ(out[w], ref[w]) << "case=" << t << " x=" << x << " word=" << w;
    }

    ASSERT_EQ(count_matches(out), count_matches(ref))
        << "case=" << t << " x=" << x;
  }
}

TEST(ExcessPositions512, TargetZero) {
  const uint64_t seed = 12345;
  std::mt19937_64 rng(seed);

  alignas(64) uint64_t s[8];
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (int t = 0; t < 500; ++t) {
    for (int w = 0; w < 8; ++w) {
      s[w] = rng();
    }
    excess_positions_512(s, 0, out);
    naive_excess_positions_512(s, 0, ref);
    check_matches_naive(excess_positions_512_expand, "expand", s, 0, t);
    check_matches_naive(excess_positions_512_expand8, "expand8", s, 0, t);
    check_matches_naive(excess_positions_512_expand_avx512, "expand_avx512", s,
                        0, t);
    check_matches_naive(excess_positions_512_branching_lut, "branching_lut", s,
                        0, t);
    check_matches_naive(excess_positions_512_lut_avx512, "lut_avx512", s, 0, t);
    check_matches_naive(excess_positions_512_byte_lut, "byte_lut", s, 0, t);
    for (int w = 0; w < 8; ++w) {
      ASSERT_EQ(out[w], ref[w]) << "case=" << t << " word=" << w;
    }
  }
}
