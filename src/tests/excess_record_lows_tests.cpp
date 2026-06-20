#include <gtest/gtest.h>
#include <pixie/bits.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <random>

// ---------------------------------------------------------------------------
// Naive reference implementations
// ---------------------------------------------------------------------------

static void naive_excess_record_lows_128(const uint64_t* s, uint64_t* out) {
  out[0] = out[1] = 0;
  int cur = 0;
  int best = 0;
  for (size_t i = 0; i < 128; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = static_cast<int>((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur < best) {
      best = cur;
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
}

static void naive_excess_record_lows_512(const uint64_t* s, uint64_t* out) {
  for (int w = 0; w < 8; ++w) {
    out[w] = 0;
  }
  int cur = 0;
  int best = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = static_cast<int>((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur < best) {
      best = cur;
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void check_128_matches_naive(const uint64_t* s, int case_id = 0) {
  alignas(16) uint64_t out[2];
  alignas(16) uint64_t ref[2];
  excess_record_lows_128(s, out);
  naive_excess_record_lows_128(s, ref);
  EXPECT_EQ(out[0], ref[0]) << "case=" << case_id;
  EXPECT_EQ(out[1], ref[1]) << "case=" << case_id;
}

static void check_512_matches_naive(const uint64_t* s, int case_id = 0) {
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];
  excess_record_lows_512(s, out);
  naive_excess_record_lows_512(s, ref);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], ref[w]) << "case=" << case_id << " word=" << w;
  }
}

static size_t count_bits(const uint64_t* words, int num_words) {
  size_t cnt = 0;
  for (int i = 0; i < num_words; ++i) {
    cnt += static_cast<size_t>(std::popcount(words[i]));
  }
  return cnt;
}

// ---------------------------------------------------------------------------
// 128-bit tests
// ---------------------------------------------------------------------------

TEST(ExcessRecordLows128, AllZeros) {
  const std::array<uint64_t, 2> s = {0, 0};
  uint64_t out[2];
  excess_record_lows_128(s.data(), out);
  // Excess goes 0, -1, -2, ..., -128. Every position is a new record low.
  EXPECT_EQ(out[0], UINT64_MAX);
  EXPECT_EQ(out[1], UINT64_MAX);
}

TEST(ExcessRecordLows128, AllOnes) {
  const std::array<uint64_t, 2> s = {UINT64_MAX, UINT64_MAX};
  uint64_t out[2];
  excess_record_lows_128(s.data(), out);
  // Excess goes 0, 1, 2, ..., 128. No new record lows after start.
  EXPECT_EQ(out[0], 0u);
  EXPECT_EQ(out[1], 0u);
}

TEST(ExcessRecordLows128, Alternating) {
  const std::array<uint64_t, 2> s = {0xAAAAAAAAAAAAAAAAull,
                                     0x5555555555555555ull};
  uint64_t out[2];
  excess_record_lows_128(s.data(), out);
  // Pattern 1010... -> excess: 0, -1, 0, -1, 0, ...
  // Only the first position (i=0) is a record low (excess=-1).
  // Actually i=0 -> first bit=1 -> excess=+1, i=1 -> bit=0 -> excess=0, wait.
  // 0xAA = 10101010, 0x55 = 01010101.
  // Let's just trust the naive for correctness.
  alignas(16) uint64_t ref[2];
  naive_excess_record_lows_128(s.data(), ref);
  EXPECT_EQ(out[0], ref[0]);
  EXPECT_EQ(out[1], ref[1]);
}

TEST(ExcessRecordLows128, FixedCases) {
  const std::array<std::array<uint64_t, 2>, 6> cases = {{
      {0, 0},
      {UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull},
      {0x0000FFFF0000FFFFull, 0xFFFF0000FFFF0000ull},
      {0x1111222233334444ull, 0x8888777766665555ull},
  }};

  int case_id = 0;
  for (const auto& s : cases) {
    check_128_matches_naive(s.data(), case_id++);
  }
}

TEST(ExcessRecordLows128, ExhaustiveSmall16) {
  alignas(16) uint64_t s[2];
  alignas(16) uint64_t out[2];
  alignas(16) uint64_t ref[2];

  for (uint64_t pattern = 0; pattern < (1ull << 16); ++pattern) {
    s[0] = pattern;
    s[1] = 0;
    excess_record_lows_128(s, out);
    naive_excess_record_lows_128(s, ref);
    EXPECT_EQ(out[0], ref[0]) << "pattern=" << pattern;
    EXPECT_EQ(out[1], ref[1]) << "pattern=" << pattern;
  }
}

TEST(ExcessRecordLows128, Random) {
  const int cases = [] {
    const char* env = std::getenv("RECORD_LOWS_CASES");
    return env ? std::atoi(env) : 10000;
  }();
  const uint64_t seed = [] {
    const char* env = std::getenv("RECORD_LOWS_SEED");
    return env ? std::stoull(env) : 42ull;
  }();

  std::mt19937_64 rng(static_cast<uint64_t>(seed));

  for (int t = 0; t < cases; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    check_128_matches_naive(s.data(), t);
  }
}

#ifdef PIXIE_AVX2_SUPPORT
static void check_nibble_lut_matches_naive(const uint64_t* s, int case_id = 0) {
  alignas(16) uint64_t out[2];
  alignas(16) uint64_t ref[2];
  excess_record_lows_128_nibble_lut(s, out);
  naive_excess_record_lows_128(s, ref);
  EXPECT_EQ(out[0], ref[0]) << "nibble case=" << case_id;
  EXPECT_EQ(out[1], ref[1]) << "nibble case=" << case_id;
}

TEST(ExcessRecordLows128, NibbleLUTExhaustiveSmall16) {
  alignas(16) uint64_t s[2];
  for (uint64_t pattern = 0; pattern < (1ull << 16); ++pattern) {
    s[0] = pattern;
    s[1] = 0;
    check_nibble_lut_matches_naive(s, static_cast<int>(pattern));
  }
}

TEST(ExcessRecordLows128, NibbleLUTRandom) {
  std::mt19937_64 rng(12345);
  for (int t = 0; t < 10000; ++t) {
    const std::array<uint64_t, 2> s = {rng(), rng()};
    check_nibble_lut_matches_naive(s.data(), t);
  }
}
#endif

// ---------------------------------------------------------------------------
// 512-bit tests
// ---------------------------------------------------------------------------

TEST(ExcessRecordLows512, AllZeros) {
  alignas(64) uint64_t s[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  alignas(64) uint64_t out[8];
  excess_record_lows_512(s, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], UINT64_MAX) << "word=" << w;
  }
}

TEST(ExcessRecordLows512, AllOnes) {
  alignas(64) uint64_t s[8] = {UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX,
                               UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX};
  alignas(64) uint64_t out[8];
  excess_record_lows_512(s, out);
  for (int w = 0; w < 8; ++w) {
    EXPECT_EQ(out[w], 0u) << "word=" << w;
  }
}

TEST(ExcessRecordLows512, FixedCases) {
  const std::array<std::array<uint64_t, 8>, 4> cases = {{
      {0, 0, 0, 0, 0, 0, 0, 0},
      {UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX,
       UINT64_MAX, UINT64_MAX},
      {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull, 0xAAAAAAAAAAAAAAAAull,
       0x5555555555555555ull, 0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
       0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull},
      {0x0123456789ABCDEFull, 0xFEDCBA9876543210ull, 0x0000FFFF0000FFFFull,
       0xFFFF0000FFFF0000ull, 0x1111222233334444ull, 0x8888777766665555ull,
       0x5555555555555555ull, 0xAAAAAAAAAAAAAAAAull},
  }};

  int case_id = 0;
  for (const auto& s : cases) {
    check_512_matches_naive(s.data(), case_id++);
  }
}

TEST(ExcessRecordLows512, ExhaustiveSmall16) {
  alignas(64) uint64_t s[8];
  alignas(64) uint64_t out[8];
  alignas(64) uint64_t ref[8];

  for (uint64_t pattern = 0; pattern < (1ull << 16); ++pattern) {
    s[0] = pattern;
    for (int w = 1; w < 8; ++w) {
      s[w] = 0;
    }
    excess_record_lows_512(s, out);
    naive_excess_record_lows_512(s, ref);
    for (int w = 0; w < 8; ++w) {
      EXPECT_EQ(out[w], ref[w]) << "pattern=" << pattern << " word=" << w;
    }
  }
}

TEST(ExcessRecordLows512, Random) {
  const int cases = [] {
    const char* env = std::getenv("RECORD_LOWS_CASES");
    return env ? std::atoi(env) : 10000;
  }();
  const uint64_t seed = [] {
    const char* env = std::getenv("RECORD_LOWS_SEED");
    return env ? std::stoull(env) : 42ull;
  }();

  std::mt19937_64 rng(static_cast<uint64_t>(seed));
  alignas(64) uint64_t s[8];

  for (int t = 0; t < cases; ++t) {
    for (int w = 0; w < 8; ++w) {
      s[w] = rng();
    }
    check_512_matches_naive(s, t);
  }
}

TEST(ExcessRecordLows512, CountMonotonic) {
  // Record-low count in a 512-bit all-zeros should be 512.
  alignas(64) uint64_t s[8] = {0};
  alignas(64) uint64_t out[8];
  excess_record_lows_512(s, out);
  EXPECT_EQ(count_bits(out, 8), 512u);

  // Record-low count in a 512-bit all-ones should be 0.
  for (int w = 0; w < 8; ++w) {
    s[w] = UINT64_MAX;
  }
  excess_record_lows_512(s, out);
  EXPECT_EQ(count_bits(out, 8), 0u);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
