#include <gtest/gtest.h>
#include <pixie/bits.h>
#include <pixie/experimental/excess.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <random>

using pixie::experimental::excess_positions_512_branching_lut;
using pixie::experimental::excess_positions_512_expand;
using pixie::experimental::excess_positions_512_expand8;
using pixie::experimental::excess_positions_512_expand_avx512;
using pixie::experimental::excess_positions_512_lut_avx512;

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
    for (int w = 0; w < 8; ++w) {
      ASSERT_EQ(out[w], ref[w]) << "case=" << t << " word=" << w;
    }
  }
}
