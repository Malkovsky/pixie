#include <gtest/gtest.h>
#include <pixie/bits.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <random>

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
    for (int w = 0; w < 8; ++w) {
      ASSERT_EQ(out[w], ref[w]) << "case=" << t << " word=" << w;
    }
  }
}

