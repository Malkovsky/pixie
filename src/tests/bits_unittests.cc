#include <gtest/gtest.h>
#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace {

using SelectBlock = std::array<uint64_t, 8>;

uint64_t naive_select_512(const uint64_t* bits, uint64_t rank, bool value) {
  uint64_t seen = 0;
  for (uint64_t i = 0; i < 512; ++i) {
    const bool bit = ((bits[i >> 6] >> (i & 63)) & 1ull) != 0;
    if (bit == value) {
      if (seen == rank) {
        return i;
      }
      ++seen;
    }
  }
  return 512;
}

uint64_t count_ones(const SelectBlock& block) {
  uint64_t count = 0;
  for (uint64_t word : block) {
    count += std::popcount(word);
  }
  return count;
}

std::vector<uint64_t> sample_ranks(uint64_t count) {
  std::vector<uint64_t> ranks;
  if (count == 0) {
    return ranks;
  }

  ranks.push_back(0);
  ranks.push_back(count / 2);
  ranks.push_back(count - 1);
  std::sort(ranks.begin(), ranks.end());
  ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
  return ranks;
}

void expect_select_samples(const SelectBlock& block,
                           const char* label,
                           int case_id) {
  const uint64_t ones = count_ones(block);
  const uint64_t zeros = 512 - ones;

  for (uint64_t rank : sample_ranks(ones)) {
    ASSERT_EQ(select_512(block.data(), rank),
              naive_select_512(block.data(), rank, true))
        << "select1 rank=" << rank << " case=" << case_id << " label=" << label;
  }
  for (uint64_t rank : sample_ranks(zeros)) {
    ASSERT_EQ(select0_512(block.data(), rank),
              naive_select_512(block.data(), rank, false))
        << "select0 rank=" << rank << " case=" << case_id << " label=" << label;
  }
}

}  // namespace

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
    EXPECT_EQ(select_64(x, i), i);
    EXPECT_EQ(select_64_no_bmi2(x, i), i);
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
        ASSERT_EQ(select_64(a, rank), i);
        ASSERT_EQ(select_64_no_bmi2(a, rank), i);
        ++rank;
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

TEST(Select512, Alternating) {
  const SelectBlock block = {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                             0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                             0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                             0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull};
  for (uint64_t rank = 0; rank < 256; ++rank) {
    ASSERT_EQ(select_512(block.data(), rank),
              naive_select_512(block.data(), rank, true));
    ASSERT_EQ(select0_512(block.data(), rank),
              naive_select_512(block.data(), rank, false));
  }
}

TEST(Select512, SingleBitPositions) {
  for (int position = 0; position < 512; ++position) {
    SelectBlock block = {0, 0, 0, 0, 0, 0, 0, 0};
    block[position >> 6] = uint64_t{1} << (position & 63);

    ASSERT_EQ(select_512(block.data(), 0), position);
    expect_select_samples(block, "single_bit", position);
  }
}

TEST(SelectZero512, SingleZeroPositions) {
  for (int position = 0; position < 512; ++position) {
    SelectBlock block = {std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max(),
                         std::numeric_limits<uint64_t>::max()};
    block[position >> 6] &= ~(uint64_t{1} << (position & 63));

    ASSERT_EQ(select0_512(block.data(), 0), position);
    expect_select_samples(block, "single_zero", position);
  }
}

TEST(Select512, EmptyWordGaps) {
  const SelectBlock block = {0, 0x8000000000000000ull,
                             0, 0x0123456789ABCDEFull,
                             0, std::numeric_limits<uint64_t>::max(),
                             0, 0x0000000000000001ull};
  const uint64_t ones = count_ones(block);
  for (uint64_t rank = 0; rank < ones; ++rank) {
    ASSERT_EQ(select_512(block.data(), rank),
              naive_select_512(block.data(), rank, true))
        << "rank=" << rank;
  }
  expect_select_samples(block, "empty_word_gaps", 0);
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
