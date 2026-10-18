#include <gtest/gtest.h>
#include <pixie/rank_select/three_star.h>

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace {

using Support = pixie::ThreeStarRankSelectSupport;
static_assert(std::derived_from<Support, pixie::RankSelectBase<Support>>);

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9E3779B97F4A7C15ull;
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ull;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBull;
  return value ^ (value >> 31);
}

TEST(ThreeStarBenchmarkAdapterTest, MatchesOneRankAndSelect) {
  constexpr std::size_t kBits = 1 << 16;
  std::vector<std::uint64_t> words(kBits / 64);
  for (std::size_t word_index = 0; word_index < words.size(); ++word_index) {
    words[word_index] = splitmix64(word_index);
  }

  pixie::ThreeStarRankSelectSupport support(words, kBits);
  std::uint64_t rank = 0;
  for (std::size_t position = 0; position < kBits; ++position) {
    const bool bit = (words[position >> 6] >> (position & 63)) & 1;
    EXPECT_EQ(support.rank(position), rank) << position;
    if (bit) {
      ++rank;
      EXPECT_EQ(support.select(rank), position) << rank;
    }
  }
  EXPECT_EQ(support.rank(kBits), rank);
}

TEST(ThreeStarBenchmarkAdapterTest, HandlesSmallSources) {
  constexpr std::array<std::uint64_t, 16> kWords{
      0x8000000000000001ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0x8000000000000000ull,
  };
  pixie::ThreeStarRankSelectSupport support(kWords, 1024);

  EXPECT_EQ(support.rank(0), 0);
  EXPECT_EQ(support.rank(1), 1);
  EXPECT_EQ(support.rank(64), 2);
  EXPECT_EQ(support.rank(1024), 3);
  EXPECT_EQ(support.select(1), 0);
  EXPECT_EQ(support.select(2), 63);
  EXPECT_EQ(support.select(3), 1023);
  EXPECT_EQ(support.select(4), 1024);
}

TEST(ThreeStarBenchmarkAdapterTest, HandlesPartialBlocksAndRepeatedBuilds) {
  constexpr std::array<std::size_t, 10> kSizes{
      1, 63, 65, 1024, 2047, 2048, 2049, 65535, 65536, 65537,
  };
  for (const auto num_bits : kSizes) {
    for (unsigned pattern = 0; pattern < 4; ++pattern) {
      SCOPED_TRACE(num_bits);
      SCOPED_TRACE(pattern);
      std::vector<std::uint64_t> words((num_bits + 63) / 64);
      for (std::size_t word = 0; word < words.size(); ++word) {
        words[word] = pattern == 0   ? 0
                      : pattern == 1 ? ~std::uint64_t{0}
                      : pattern == 2 ? 0x8000000000000001ull
                                     : splitmix64(word);
      }
      pixie::ThreeStarRankSelectSupport support(words, num_bits);
      EXPECT_EQ(support.size(), num_bits);
      EXPECT_FALSE(support.empty());
      EXPECT_TRUE(support.supports_select1());
      EXPECT_FALSE(support.supports_select0());
      std::string expected;
      std::uint64_t rank = 0;
      for (std::size_t position = 0; position < num_bits; ++position) {
        ASSERT_EQ(support.rank(position), rank) << position;
        ASSERT_EQ(support.rank0(position), position - rank) << position;
        const bool bit = (words[position >> 6] >> (position & 63)) & 1;
        ASSERT_EQ(support[position], bit) << position;
        expected.push_back(bit ? '1' : '0');
        if (bit) {
          ASSERT_EQ(support.select(++rank), position) << rank;
        }
      }
      EXPECT_EQ(support.rank(num_bits), rank);
      EXPECT_EQ(support.rank0(num_bits), num_bits - rank);
      EXPECT_EQ(support.select(0), 0);
      EXPECT_EQ(support.select(rank + 1), num_bits);
      EXPECT_EQ(support.rank(std::numeric_limits<std::size_t>::max()), rank);
      EXPECT_EQ(support.rank0(std::numeric_limits<std::size_t>::max()),
                num_bits - rank);
      EXPECT_EQ(support.select0(0), 0u);
      EXPECT_EQ(support.select0(1), num_bits);
      EXPECT_EQ(support.select0(std::numeric_limits<std::size_t>::max()),
                num_bits);
      EXPECT_EQ(support.to_string(), expected);
    }
  }
}

}  // namespace
