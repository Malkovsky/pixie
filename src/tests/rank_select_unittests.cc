#include <gtest/gtest.h>
#include <pixie/rank_select/implementations.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

using RankSelectSupport = pixie::RankSelectSupport<>;

TEST(RankSelectSupportTest, Basic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b101010101;
  RankSelectSupport bv(bits, 9);
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

TEST(RankSelectSupportTest, MultipleWords) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0xAAAAAAAAAAAAAAAA;
  bits[1] = 0x5555555555555555;
  RankSelectSupport bv(bits, 128);
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

TEST(RankSelectSupportTest, ToString) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b10101010101;
  RankSelectSupport bv(bits, 11);
  EXPECT_EQ(bv.to_string(), "10101010101");
}

TEST(RankSelectSupportTest, RankBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b10110;
  RankSelectSupport bv(bits, 5);

  EXPECT_EQ(bv.rank(0), 0);  // No bits
  EXPECT_EQ(bv.rank(1), 0);  // 0
  EXPECT_EQ(bv.rank(2), 1);  // 10
  EXPECT_EQ(bv.rank(3), 2);  // 110
  EXPECT_EQ(bv.rank(4), 2);  // 0110
  EXPECT_EQ(bv.rank(5), 3);  // 10110
}

TEST(RankSelectSupportTest, RankWithZeros) {
  std::vector<uint64_t> bits(8, 0);
  RankSelectSupport bv(bits, 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(bv.rank(i), 0);
  }
}

TEST(RankSelectSupportTest, SelectBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b1100010110010110;
  RankSelectSupport bv(bits, 16);

  EXPECT_EQ(bv.select(1), 1);
  EXPECT_EQ(bv.select(2), 2);
  EXPECT_EQ(bv.select(3), 4);
  EXPECT_EQ(bv.select(4), 7);
  EXPECT_EQ(bv.select(5), 8);
  EXPECT_EQ(bv.select(6), 10);
  EXPECT_EQ(bv.select(7), 14);
  EXPECT_EQ(bv.select(8), 15);
}

TEST(RankSelectSupportTest, MainRankTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  RankSelectSupport bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;
  for (size_t i = 0; i < bv.size(); ++i) {
    ASSERT_EQ(rank, bv.rank(i));
    rank += bv[i];
  }
}

TEST(RankSelectSupportTest, MainSelectTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  RankSelectSupport bv(bits, 65536 * 32 * 64);
  uint64_t rank = 0;

  for (size_t i = 0; i < bv.size(); ++i) {
    if (bv[i]) {
      ASSERT_EQ(bv.select(++rank), i);
      ASSERT_EQ(bv.rank(i), rank - 1);
      ASSERT_EQ(bv.rank(i + 1), rank);
    }
  }
}

TEST(RankSelectSupportTest, RankZeroBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b10110;
  RankSelectSupport bv(bits, 10);

  EXPECT_EQ(bv.rank0(0), 0);
  EXPECT_EQ(bv.rank0(1), 1);
  EXPECT_EQ(bv.rank0(2), 1);
  EXPECT_EQ(bv.rank0(3), 1);
  EXPECT_EQ(bv.rank0(4), 2);
  EXPECT_EQ(bv.rank0(5), 2);
}

TEST(RankSelectSupportTest, RankZeroWithOnes) {
  std::vector<uint64_t> bits(8, 0);
  RankSelectSupport bv(bits, 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(bv.rank0(i), i);
  }
}

TEST(RankSelectSupportTest, SelectZeroBasic) {
  std::vector<uint64_t> bits(8, 0);
  bits[0] = 0b1100010110010110;
  RankSelectSupport bv(bits, 16);

  EXPECT_EQ(bv.select0(1), 0);
  EXPECT_EQ(bv.select0(2), 3);
  EXPECT_EQ(bv.select0(3), 5);
  EXPECT_EQ(bv.select0(4), 6);
  EXPECT_EQ(bv.select0(5), 9);
  EXPECT_EQ(bv.select0(6), 11);
  EXPECT_EQ(bv.select0(7), 12);
  EXPECT_EQ(bv.select0(8), 13);
}

TEST(RankSelectSupportTest, EmptyAndZeroRankGuards) {
  RankSelectSupport default_bv;
  EXPECT_EQ(default_bv.size(), 0u);
  EXPECT_EQ(default_bv.rank(0), 0u);
  EXPECT_EQ(default_bv.rank0(0), 0u);
  EXPECT_EQ(default_bv.select(0), 0u);
  EXPECT_EQ(default_bv.select(1), 0u);
  EXPECT_EQ(default_bv.select0(0), 0u);
  EXPECT_EQ(default_bv.select0(1), 0u);
  EXPECT_EQ(default_bv.to_string(), "");

  std::vector<uint64_t> empty;
  RankSelectSupport bv(std::span<const uint64_t>(empty), 0);

  EXPECT_EQ(bv.size(), 0u);
  EXPECT_EQ(bv.rank(0), 0u);
  EXPECT_EQ(bv.rank0(0), 0u);
  EXPECT_EQ(bv.select(0), 0u);
  EXPECT_EQ(bv.select(1), 0u);
  EXPECT_EQ(bv.select0(0), 0u);
  EXPECT_EQ(bv.select0(1), 0u);
  EXPECT_EQ(bv.to_string(), "");

  std::vector<uint64_t> bits(1, 0b10110);
  RankSelectSupport non_empty(std::span<const uint64_t>(bits), 5);
  EXPECT_EQ(non_empty.select(0), 0u);
  EXPECT_EQ(non_empty.select0(0), 0u);
}

TEST(RankSelectSupportTest, ExactShortSpanRankSelect) {
  std::vector<uint64_t> bits = {0b1100010110010110};
  RankSelectSupport bv(bits, 16);

  EXPECT_EQ(bv.rank(16), 8);
  EXPECT_EQ(bv.rank0(16), 8);
  EXPECT_EQ(bv.select(1), 1);
  EXPECT_EQ(bv.select(8), 15);
  EXPECT_EQ(bv.select(9), bv.size());
  EXPECT_EQ(bv.select0(1), 0);
  EXPECT_EQ(bv.select0(8), 13);
  EXPECT_EQ(bv.select0(9), bv.size());
}

TEST(RankSelectSupportTest, OptionalSelectSupport) {
  std::vector<uint64_t> bits = {0b1100010110010110};

  RankSelectSupport select0_only(bits, 16,
                                 RankSelectSupport::SelectSupport::kSelect0);
  EXPECT_FALSE(select0_only.supports_select1());
  EXPECT_TRUE(select0_only.supports_select0());
  EXPECT_EQ(select0_only.rank(16), 8);
  EXPECT_EQ(select0_only.rank0(16), 8);
  EXPECT_EQ(select0_only.select(1), select0_only.size());
  EXPECT_EQ(select0_only.select0(1), 0);
  EXPECT_EQ(select0_only.select0(8), 13);

  RankSelectSupport hinted_both(bits, 16,
                                RankSelectSupport::SelectSupport::kBoth, 8);
  EXPECT_TRUE(hinted_both.supports_select1());
  EXPECT_TRUE(hinted_both.supports_select0());
  EXPECT_EQ(hinted_both.select(1), 1);
  EXPECT_EQ(hinted_both.select(8), 15);
  EXPECT_EQ(hinted_both.select0(1), 0);
  EXPECT_EQ(hinted_both.select0(8), 13);
  EXPECT_THROW((RankSelectSupport(bits, 16,
                                  RankSelectSupport::SelectSupport::kBoth, 17)),
               std::invalid_argument);

  RankSelectSupport select1_only(bits, 16,
                                 RankSelectSupport::SelectSupport::kSelect1);
  EXPECT_TRUE(select1_only.supports_select1());
  EXPECT_FALSE(select1_only.supports_select0());
  EXPECT_EQ(select1_only.rank(16), 8);
  EXPECT_EQ(select1_only.rank0(16), 8);
  EXPECT_EQ(select1_only.select(1), 1);
  EXPECT_EQ(select1_only.select(8), 15);
  EXPECT_EQ(select1_only.select0(1), select1_only.size());

  RankSelectSupport rank_only(bits, 16,
                              RankSelectSupport::SelectSupport::kNone);
  EXPECT_FALSE(rank_only.supports_select1());
  EXPECT_FALSE(rank_only.supports_select0());
  EXPECT_EQ(rank_only.rank(16), 8);
  EXPECT_EQ(rank_only.rank0(16), 8);
  EXPECT_EQ(rank_only.select(0), 0);
  EXPECT_EQ(rank_only.select0(0), 0);
  EXPECT_EQ(rank_only.select(1), rank_only.size());
  EXPECT_EQ(rank_only.select0(1), rank_only.size());

  std::vector<uint64_t> zero_words(512, 0);
  RankSelectSupport large_select0_only(
      zero_words, zero_words.size() * 64,
      RankSelectSupport::SelectSupport::kSelect0);
  EXPECT_EQ(large_select0_only.select(1), large_select0_only.size());
  EXPECT_EQ(large_select0_only.select0(1), 0);
  EXPECT_EQ(large_select0_only.select0(16384), 16383);
  EXPECT_EQ(large_select0_only.select0(16385), 16384);
  EXPECT_EQ(large_select0_only.select0(32768), 32767);

  RankSelectSupport large_hinted_select0(
      zero_words, zero_words.size() * 64,
      RankSelectSupport::SelectSupport::kSelect0, 0);
  EXPECT_EQ(large_hinted_select0.select(1), large_hinted_select0.size());
  EXPECT_EQ(large_hinted_select0.select0(1), 0);
  EXPECT_EQ(large_hinted_select0.select0(32768), 32767);

  std::vector<uint64_t> one_words(512, ~uint64_t{0});
  RankSelectSupport large_select1_only(
      one_words, one_words.size() * 64,
      RankSelectSupport::SelectSupport::kSelect1);
  EXPECT_EQ(large_select1_only.select0(1), large_select1_only.size());
  EXPECT_EQ(large_select1_only.select(1), 0);
  EXPECT_EQ(large_select1_only.select(16384), 16383);
  EXPECT_EQ(large_select1_only.select(16385), 16384);
  EXPECT_EQ(large_select1_only.select(32768), 32767);
}

TEST(RankSelectSupportTest, ThrowsWhenOneCountHintUnderestimatesSamples) {
  std::vector<uint64_t> one_words(512, ~uint64_t{0});

  EXPECT_THROW(
      (RankSelectSupport(one_words, one_words.size() * 64,
                         RankSelectSupport::SelectSupport::kSelect1, 0)),
      std::invalid_argument);
}

TEST(RankSelectSupportTest, ShortSpanUsesScalarRankSelectFallbacks) {
  std::vector<uint64_t> one_in_second_word = {0, 1};
  RankSelectSupport ones(std::span<const uint64_t>(one_in_second_word), 65);
  EXPECT_EQ(ones.rank(64), 0u);
  EXPECT_EQ(ones.rank(65), 1u);
  EXPECT_EQ(ones.rank0(65), 64u);
  EXPECT_EQ(ones.select(1), 64u);

  std::vector<uint64_t> zero_in_second_word = {
      std::numeric_limits<uint64_t>::max(), 0};
  RankSelectSupport zeros(std::span<const uint64_t>(zero_in_second_word), 65);
  EXPECT_EQ(zeros.rank(65), 64u);
  EXPECT_EQ(zeros.rank0(64), 0u);
  EXPECT_EQ(zeros.rank0(65), 1u);
  EXPECT_EQ(zeros.select0(1), 64u);
}

TEST(RankSelectSupportTest, ShortSpanScalarRankUsesPartialWordTail) {
  std::vector<uint64_t> bits = {~uint64_t{0}, 0, 0b101};
  RankSelectSupport bv(std::span<const uint64_t>(bits), 130);

  EXPECT_EQ(bv.rank(128), 64u);
  EXPECT_EQ(bv.rank(129), 65u);
  EXPECT_EQ(bv.rank(130), 65u);
  EXPECT_EQ(bv.rank0(129), 64u);
}

TEST(RankSelectSupportTest, IgnoresDirtyPaddingAndTrailingWords) {
  std::vector<uint64_t> bits = {~uint64_t{0}, ~uint64_t{0}};
  RankSelectSupport bv(bits, 3);

  EXPECT_EQ(bv.size(), 3);
  EXPECT_EQ(bv.rank(3), 3);
  EXPECT_EQ(bv.rank(128), 3);
  EXPECT_EQ(bv.rank0(3), 0);
  EXPECT_EQ(bv.rank0(128), 0);
  EXPECT_EQ(bv.select(1), 0);
  EXPECT_EQ(bv.select(3), 2);
  EXPECT_EQ(bv.select(4), bv.size());
  EXPECT_EQ(bv.select0(1), bv.size());
}

TEST(RankSelectSupportTest, SelectZeroSkewedDistributionFallback) {
  constexpr size_t kWordsPerBasicBlock = 8;
  constexpr size_t kBasicBlocksPerSuperBlock = 128;
  constexpr size_t kZeroBasicBlocks = 40;

  std::vector<uint64_t> bits(kWordsPerBasicBlock * kBasicBlocksPerSuperBlock,
                             ~uint64_t{0});
  std::fill(bits.begin(), bits.begin() + kZeroBasicBlocks * kWordsPerBasicBlock,
            0);

  RankSelectSupport bv(std::span<const uint64_t>(bits), bits.size() * 64);
  EXPECT_EQ(bv.rank0(bv.size()), kZeroBasicBlocks * 512u);
  EXPECT_EQ(bv.select0(1), 0u);
  EXPECT_EQ(bv.select0(10000), 9999u);
  EXPECT_EQ(bv.select0(kZeroBasicBlocks * 512u), kZeroBasicBlocks * 512u - 1);
  EXPECT_EQ(bv.select0(kZeroBasicBlocks * 512u + 1), bv.size());
}

TEST(RankSelectSupportTest, MainRankZeroTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < bits.size(); i++) {
    bits[i] = rng();
  }

  RankSelectSupport bv(bits, bits.size() * 64);
  uint64_t zero_count = 0;
  for (size_t i = 0; i < bv.size(); ++i) {
    ASSERT_EQ(zero_count, bv.rank0(i));
    zero_count += (bv[i] == 0) ? 1 : 0;
  }
}

TEST(RankSelectSupportTest, MainSelectZeroTest) {
  std::mt19937_64 rng(42);
  std::vector<uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < bits.size(); i++) {
    bits[i] = rng();
  }

  RankSelectSupport bv(bits, bits.size() * 64);
  uint64_t zero_rank = 0;

  for (size_t i = 0; i < bv.size(); ++i) {
    if (bv[i] == 0) {
      zero_rank++;
      EXPECT_EQ(bv.select0(zero_rank), i);
      EXPECT_EQ(bv.rank0(i), zero_rank - 1);
      EXPECT_EQ(bv.rank0(i + 1), zero_rank);
    }
  }
}
