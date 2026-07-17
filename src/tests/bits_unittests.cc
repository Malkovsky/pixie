#include <gtest/gtest.h>
#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

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

TEST(LowerBoundDelta4x64, Random) {
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
      ASSERT_EQ(lower_bound_delta_4x64(x.data(), y, dlt_array, dlt_scalar),
                cnt);
    } else {
      ASSERT_GE(lower_bound_delta_4x64(x.data(), y, dlt_array, dlt_scalar),
                cnt);
    }
  }
}

TEST(LowerBoundDelta8x64, Random) {
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
      ASSERT_EQ(lower_bound_delta_8x64(x.data(), y, dlt_array, dlt_scalar),
                cnt);
    } else {
      ASSERT_GE(lower_bound_delta_8x64(x.data(), y, dlt_array, dlt_scalar),
                cnt);
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

TEST(LowerBoundDelta32x16, Random) {
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
    ASSERT_EQ(lower_bound_delta_32x16(x.data(), y, dlt_array, dlt_scalar), cnt);
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
