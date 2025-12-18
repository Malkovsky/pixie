#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "bits.h"
#include "bitvector.h"

using pixie::BitVector;
using pixie::BitVectorInterleaved;

TEST(Rank512, Ones) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};
  for (size_t i = 0; i < 512; ++i) {
    auto p = rank_512(a.data(), i);
    EXPECT_EQ(p, i);
  }
}

TEST(Rank512, Random) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};

  std::mt19937_64 rng(42);
  for (size_t t = 0; t < 1000; ++t) {
    for (size_t i = 0; i < 8; ++i) {
      a[i] = rng();
    }
    size_t rank = 0;
    for (size_t i = 0; i < 512; ++i) {
      auto p = rank_512(a.data(), i);
      ASSERT_EQ(p, rank);
      rank += 1 & (a[i >> 6] >> (i & 63));
    }
    auto p = rank_512(a.data(), 512);
    ASSERT_EQ(p, rank);
  }
}

TEST(Select64, Ones) {
  uint64_t x = std::numeric_limits<uint64_t>::max();
  for (size_t i = 0; i < 64; ++i) {
    auto p = select_64(x, i);
    EXPECT_EQ(p, i);
  }
}

TEST(Select64, Random) {
  uint64_t a;

  std::mt19937_64 rng(42);
  for (size_t t = 0; t < 1000; ++t) {
    a = rng();
    size_t rank = 0;
    for (size_t i = 0; i < 64; ++i) {
      if (1 & (a >> i)) {
        auto p = select_64(a, rank++);
        ASSERT_EQ(p, i);
      }
    }
  }
}

TEST(Select512, Ones) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};
  for (size_t i = 0; i < 512; ++i) {
    auto p = select_512(a.data(), i);
    EXPECT_EQ(p, i);
  }
}

TEST(SelectZero512, Zeros) {
  std::array<uint64_t, 8> a{0, 0, 0, 0, 0, 0, 0, 0};
  for (size_t i = 0; i < 512; ++i) {
    auto p = select0_512(a.data(), i);
    EXPECT_EQ(p, i);
  }
}

TEST(Select512, Random) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};

  std::mt19937_64 rng(42);
  for (size_t t = 0; t < 1000; ++t) {
    for (auto& x : a) {
      x = rng();
    }
    size_t rank = 0;
    for (size_t i = 0; i < 512; ++i) {
      if (1 & (a[i >> 6] >> (i & 63))) {
        auto p = select_512(a.data(), rank++);
        ASSERT_EQ(p, i);
      }
    }
  }
}

TEST(SelectZero512, Random) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};

  std::mt19937_64 rng(42);
  for (size_t t = 0; t < 1000; ++t) {
    for (auto& x : a) {
      x = rng();
    }
    size_t rank = 0;
    for (size_t i = 0; i < 512; ++i) {
      if ((1 & (a[i >> 6] >> (i & 63))) == 0) {
        auto p = select0_512(a.data(), rank++);
        ASSERT_EQ(p, i);
      }
    }
  }
}

TEST(Select512, RankCompativility) {
  std::array<uint64_t, 8> a{std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max(),
                            std::numeric_limits<uint64_t>::max()};

  std::mt19937_64 rng(42);
  for (size_t t = 0; t < 1000; ++t) {
    for (auto& x : a) {
      x = rng();
    }
    size_t rank = 0;
    for (size_t i = 0; i < 512; ++i) {
      if (1 & (a[i >> 6] >> (i & 63))) {
        auto p = select_512(a.data(), rank++);
        ASSERT_EQ(p, i);
        auto r = rank_512(a.data(), p);
        ASSERT_EQ(r + 1, rank);
        r = rank_512(a.data(), p + 1);
        ASSERT_EQ(r, rank);
      }
    }
  }
}

TEST(BitVectorTest, Basic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b101010101;
  BitVector bv(bits, 9);
  EXPECT_EQ(bv.size(), 9);
  EXPECT_EQ(bv[0], 1);
  EXPECT_EQ(bv[1], 0);
  EXPECT_EQ(bv[2], 1);
  EXPECT_EQ(bv[3], 0);
  EXPECT_EQ(bv[4], 1);
  EXPECT_EQ(bv[5], 0);
  EXPECT_EQ(bv[6], 1);
  EXPECT_EQ(bv[7], 0);
  EXPECT_EQ(bv[8], 1);
}

TEST(BitVectorTest, MultipleWords) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0xAAAAAAAAAAAAAAAA;
  bits[1] = 0x5555555555555555;
  BitVector bv(bits, 128);
  EXPECT_EQ(bv.size(), 128);

  // First word: 0xAAAAAAAAAAAAAAAA (alternating 1s and 0s, starting with 0)
  for (size_t i = 0; i < 64; i++) {
    EXPECT_EQ(bv[i], (i % 2 == 1));
  }

  // Second word: 0x5555555555555555 (alternating 0s and 1s, starting with 1)
  for (size_t i = 64; i < 128; i++) {
    EXPECT_EQ(bv[i], (i % 2 == 0));
  }
}

TEST(BitVectorTest, ToString) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b10101010101;
  BitVector bv(bits, 11);
  EXPECT_EQ(bv.to_string(), "10101010101");
}

TEST(BitVectorTest, RankBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b10110;
  BitVector bv(bits, 5);

  EXPECT_EQ(bv.rank(0), 0);  // No bits
  EXPECT_EQ(bv.rank(1), 0);  // 0
  EXPECT_EQ(bv.rank(2), 1);  // 10
  EXPECT_EQ(bv.rank(3), 2);  // 110
  EXPECT_EQ(bv.rank(4), 2);  // 0110
  EXPECT_EQ(bv.rank(5), 3);  // 10110
}

TEST(BitVectorTest, RankWithZeros) {
  std::vector<uint64_t> bits(8, 0);
  BitVector bv(bits, 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(bv.rank(i), 0);
  }
}

TEST(BitVectorTest, SelectBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b1100010110010110;
  BitVector bv(bits, 5);

  EXPECT_EQ(bv.select(1), 1);
  EXPECT_EQ(bv.select(2), 2);
  EXPECT_EQ(bv.select(3), 4);
  EXPECT_EQ(bv.select(4), 7);
  EXPECT_EQ(bv.select(5), 8);
  EXPECT_EQ(bv.select(6), 10);
  EXPECT_EQ(bv.select(7), 14);
  EXPECT_EQ(bv.select(8), 15);
}

TEST(BitVectorTest, MainRankTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  BitVector bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;
  for (size_t i = 0; i < bv.size(); ++i) {
    ASSERT_EQ(rank, bv.rank(i));
    rank += bv[i];
  }
}

TEST(BitVectorTest, MainSelectTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  BitVector bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;

  for (size_t i = 0; i < bv.size(); ++i) {
    if (bv[i]) {
      ASSERT_EQ(bv.select(++rank), i);
      ASSERT_EQ(bv.rank(i), rank - 1);
      ASSERT_EQ(bv.rank(i + 1), rank);
    }
  }
}

TEST(BitVectorInterleavedTest, AtTest) {
  std::mt19937_64 rng(42);
  /*
   * TODO: need to check edge cases, at least
   *       non-multiple of 64 size
   */
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  BitVectorInterleaved bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;
  for (size_t i = 0; i < 65536 * 32; ++i) {
    for (size_t j = 0; j < 64; ++j) {
      ASSERT_EQ((bits[i] >> j) & 1, bv[i * 64 + j]);
    }
  }
}

TEST(BitVectorInterleavedTest, MainTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  BitVectorInterleaved bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;
  for (size_t i = 0; i < bv.size(); ++i) {
    ASSERT_EQ(rank, bv.rank(i));
    rank += bv[i];
  }
}

TEST(LowerBound4x64, Random) {
  std::vector<uint64_t> x(8);
  std::mt19937_64 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    uint64_t y = rng();
    uint16_t cnt = 0;
    bool fl = 1;
    for (size_t j = 0; j < 4; j++) {
      x[j] = rng();
      fl &= x[j] < y;
      cnt += fl;
    }
    if (cnt < 4) {
      ASSERT_EQ(lower_bound_4x64(x.data(), y), cnt);
    } else {
      ASSERT_GE(lower_bound_4x64(x.data(), y), cnt);
    }
  }
}

TEST(LowerBound8x64, Random) {
  std::vector<uint64_t> x(8);
  std::mt19937_64 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    uint64_t y = rng();
    uint16_t cnt = 0;
    bool fl = 1;
    for (size_t j = 0; j < 8; j++) {
      x[j] = rng();
      fl &= x[j] < y;
      cnt += fl;
    }
    if (cnt < 8) {
      ASSERT_EQ(lower_bound_8x64(x.data(), y), cnt);
    } else {
      ASSERT_GE(lower_bound_8x64(x.data(), y), cnt);
    }
  }
}

TEST(LowerBoundDlt4x64, Random) {
  std::vector<uint64_t> x(4);
  uint64_t dlt_array[4];
  std::mt19937_64 rng(42);
  for (size_t i = 0; i < 100000; i++) {
    uint64_t y = rng();
    uint64_t dlt_scalar = rng();
    uint16_t cnt = 0;
    bool fl = 1;
    for (size_t j = 0; j < 4; j++) {
      dlt_array[j] = rng();
      x[j] = rng();
      fl &= dlt_scalar + dlt_array[j] - x[j] < y;
      cnt += fl;
    }
    if (cnt < 4) {
      ASSERT_EQ(lower_bound_dlt_4x64(x.data(), y, dlt_array, dlt_scalar), cnt);
    } else {
      ASSERT_GE(lower_bound_dlt_4x64(x.data(), y, dlt_array, dlt_scalar), cnt);
    }
  }
}

TEST(LowerBoundDlt8x64, Random) {
  std::vector<uint64_t> x(8);
  uint64_t dlt_array[8];
  std::mt19937_64 rng(42);
  for (size_t i = 0; i < 100000; i++) {
    uint64_t y = rng();
    uint64_t dlt_scalar = rng();
    uint16_t cnt = 0;
    bool fl = 1;
    for (size_t j = 0; j < 8; j++) {
      dlt_array[j] = rng();
      x[j] = rng();
      fl &= dlt_scalar + dlt_array[j] - x[j] < y;
      cnt += fl;
    }
    if (cnt < 8) {
      ASSERT_EQ(lower_bound_dlt_8x64(x.data(), y, dlt_array, dlt_scalar), cnt);
    } else {
      ASSERT_GE(lower_bound_dlt_8x64(x.data(), y, dlt_array, dlt_scalar), cnt);
    }
  }
}

TEST(LowerBound32x16, Random) {
  std::vector<uint16_t> x(32);
  std::mt19937 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    uint16_t y = rng();
    uint16_t cnt = 0;
    for (size_t j = 0; j < 32; j++) {
      x[j] = rng();
      cnt += x[j] < y;
    }
    ASSERT_EQ(lower_bound_32x16(x.data(), y), cnt);
  }
}

TEST(LowerBoundDlt32x16, Random) {
  std::vector<uint16_t> x(32);
  uint16_t dlt_array[32];
  std::mt19937 rng(42);
  for (size_t i = 0; i < 100000; i++) {
    uint16_t y = rng();
    uint16_t dlt_scalar = rng();
    uint16_t cnt = 0;
    for (size_t j = 0; j < 32; j++) {
      x[j] = rng();
      if (dlt_scalar < x[j]) {
        dlt_array[j] =
            x[j] - dlt_scalar + rng() % ((1 << 16) - x[j] + dlt_scalar);
      } else {
        dlt_array[j] = rng() % ((1 << 16) - dlt_scalar);
      }
      cnt += dlt_array[j] + dlt_scalar - x[j] < y;
    }
    ASSERT_EQ(lower_bound_dlt_32x16(x.data(), y, dlt_array, dlt_scalar), cnt);
  }
}

TEST(Popcount64x4, Random) {
  std::vector<uint8_t> x(32);
  std::vector<uint8_t> y(32);
  std::mt19937 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    for (size_t j = 0; j < 32; j++) {
      x[j] = rng();
    }
    popcount_64x4(x.data(), y.data());
    for (size_t j = 0; j < 32; ++j) {
      uint8_t a = x[j] & 0x0F;
      ASSERT_EQ(y[j] & 0x0F, std::popcount(a));
      a = x[j] >> 4;
      ASSERT_EQ(y[j] >> 4, std::popcount(a));
    }
  }
}

TEST(Popcount32x8, Random) {
  std::vector<uint8_t> x(32);
  std::vector<uint8_t> y(32);
  std::mt19937 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    for (size_t j = 0; j < 32; j++) {
      x[j] = rng();
    }
    popcount_32x8(x.data(), y.data());
    for (size_t j = 0; j < 32; ++j) {
      ASSERT_EQ(y[j], std::popcount(x[j]));
    }
  }
}

TEST(RankParallel32x8, Random) {
  std::vector<uint8_t> x(32);
  std::vector<uint8_t> sum(32);
  std::vector<uint8_t> y(32);
  std::mt19937 rng(42);
  for (size_t i = 0; i < 1000; i++) {
    x[0] = rng();
    sum[0] = std::popcount(x[0]);
    for (size_t j = 1; j < 32; j++) {
      x[j] = rng();
      sum[j] = sum[j - 1] + std::popcount(x[j]);
    }
    rank_32x8(x.data(), y.data());
    ASSERT_EQ(std::equal(y.begin(), y.end(), sum.begin()), true);
  }
}
