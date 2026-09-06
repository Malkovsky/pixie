#include <gtest/gtest.h>
#include <pixie/huffman/implementations.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {

template <class Codec>
class HuffmanTest : public testing::Test {};

using HuffmanImplementations = testing::Types<pixie::PivCoHuffman>;
TYPED_TEST_SUITE(HuffmanTest, HuffmanImplementations);

static_assert(std::is_base_of_v<pixie::HuffmanBase<pixie::PivCoHuffman>,
                                pixie::PivCoHuffman>);

template <class Codec>
void expect_round_trip(const std::vector<std::uint8_t>& input) {
  const Codec codec(input);
  EXPECT_EQ(codec.uncompressed_size(), input.size());
  EXPECT_EQ(codec.empty(), input.empty());
  EXPECT_EQ(codec.compressed_data().size(), codec.compressed_size());
  EXPECT_EQ(codec.decode(), input);

  const Codec loaded(codec.compressed_data());
  EXPECT_EQ(loaded.uncompressed_size(), input.size());
  EXPECT_EQ(loaded.empty(), input.empty());
  EXPECT_EQ(loaded.decode(), input);
}

TYPED_TEST(HuffmanTest, EmptyInputRoundTrips) {
  expect_round_trip<TypeParam>({});
}

TYPED_TEST(HuffmanTest, KnownInputRoundTrips) {
  const std::vector<std::uint8_t> input{'p', 'i', 'v', 'c', 'o', ' ', 'h',
                                        'u', 'f', 'f', 'm', 'a', 'n'};
  expect_round_trip<TypeParam>(input);
}

TYPED_TEST(HuffmanTest, FullAlphabetRoundTrips) {
  std::vector<std::uint8_t> input;
  for (std::size_t repeat = 0; repeat < 32; ++repeat) {
    for (std::size_t symbol = 0; symbol < 256; ++symbol) {
      input.push_back(static_cast<std::uint8_t>(symbol));
    }
  }
  expect_round_trip<TypeParam>(input);
}

TYPED_TEST(HuffmanTest, MultipleBlocksRoundTrip) {
  std::mt19937_64 rng(42);
  std::vector<std::uint8_t> input(TypeParam::kBlockSize * 3 + 97);
  for (std::uint8_t& symbol : input) {
    symbol = static_cast<std::uint8_t>(rng());
  }
  expect_round_trip<TypeParam>(input);
}

TEST(PivCoHuffmanTest, RejectsMalformedSerializedStream) {
  const std::array<std::byte, 32> malformed{};
  EXPECT_THROW(pixie::PivCoHuffman(std::span<const std::byte>(malformed)),
               std::runtime_error);
}

}  // namespace
