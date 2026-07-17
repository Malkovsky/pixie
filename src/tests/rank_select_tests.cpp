#include <gtest/gtest.h>
#include <pixie/rank_select/implementations.h>

#include <bit>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace {

template <class Support>
class RankSelectSpecificationTest : public testing::Test {};

using RankSelectImplementations = testing::Types<pixie::RankSelectSupport<>>;
TYPED_TEST_SUITE(RankSelectSpecificationTest, RankSelectImplementations);

TYPED_TEST(RankSelectSpecificationTest, EmptySequence) {
  const std::vector<std::uint64_t> words;
  const TypeParam bits(words, 0);

  EXPECT_TRUE(bits.empty());
  EXPECT_EQ(bits.size(), 0u);
  EXPECT_EQ(bits.rank(0), 0u);
  EXPECT_EQ(bits.rank0(0), 0u);
  EXPECT_EQ(bits.select(1), 0u);
  EXPECT_EQ(bits.select0(1), 0u);
  EXPECT_EQ(bits.to_string(), "");
}

TYPED_TEST(RankSelectSpecificationTest, RankSelectAndAccessMatchPackedBits) {
  constexpr std::size_t sizes[] = {1,    2,    63,    64,    65,    127,
                                   495,  496,  497,   511,   512,   513,
                                   1024, 4097, 65535, 65536, 65537, 131073};
  std::mt19937_64 rng(42);

  for (const std::size_t size : sizes) {
    std::vector<std::uint64_t> words((size + 63) / 64);
    for (auto& word : words) {
      word = rng();
    }
    const TypeParam bits(words, size);
    ASSERT_FALSE(bits.empty());
    ASSERT_EQ(bits.size(), size);
    ASSERT_TRUE(bits.supports_select1());
    ASSERT_TRUE(bits.supports_select0());

    std::size_t ones = 0;
    std::size_t zeros = 0;
    std::string expected;
    expected.reserve(size);
    for (std::size_t position = 0; position < size; ++position) {
      const bool value = ((words[position / 64] >> (position % 64)) & 1u) != 0;
      EXPECT_EQ(bits.rank(position), ones) << "size=" << size;
      EXPECT_EQ(bits.rank0(position), zeros) << "size=" << size;
      EXPECT_EQ(bits[position], value) << "size=" << size;
      expected.push_back(value ? '1' : '0');
      if (value) {
        ++ones;
        EXPECT_EQ(bits.select(ones), position) << "size=" << size;
      } else {
        ++zeros;
        EXPECT_EQ(bits.select0(zeros), position) << "size=" << size;
      }
    }
    EXPECT_EQ(bits.rank(size), ones);
    EXPECT_EQ(bits.rank0(size), zeros);
    EXPECT_EQ(bits.select(ones + 1), size);
    EXPECT_EQ(bits.select0(zeros + 1), size);
    EXPECT_EQ(bits.to_string(), expected);
  }
}

TEST(RankSelectSupportTest, AcceptsStorageSourceWithoutCopying) {
  pixie::AlignedStorage source(64);
  source.writable_words64()[0] = 0b101101;
  const pixie::RankSelectSupport<> support(source, 6);

  EXPECT_EQ(support.rank(6), 4u);
  EXPECT_EQ(support.select(4), 5u);
  source.writable_words64()[0] = 0b000001;
  EXPECT_EQ(support[2], 0);
}

}  // namespace
