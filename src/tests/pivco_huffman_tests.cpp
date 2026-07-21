#include <gtest/gtest.h>
#include <pixie/huffman/implementations.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <type_traits>
#include <vector>

using pixie::PivCoHuffman;

static_assert(
    std::is_base_of_v<pixie::HuffmanBase<PivCoHuffman>, PivCoHuffman>);
static_assert(std::is_same_v<PivCoHuffman::symbol_type, std::uint8_t>);

namespace {
// Round-trip through the in-memory tree of the same codec instance.
std::vector<std::uint8_t> decode_same(const PivCoHuffman& codec) {
  return codec.decode();
}

// Round-trip through the serialized byte stream: simulate writing the
// compressed form to disk and loading it into a fresh codec instance.
std::vector<std::uint8_t> decode_via_bytes(const PivCoHuffman& codec) {
  const auto bytes = codec.compressed_data();
  PivCoHuffman copy(bytes);
  return copy.decode();
}
}  // namespace

TEST(PivCoHuffmanSmoke, EmptyInputRoundTrips) {
  const std::vector<std::uint8_t> data;
  const PivCoHuffman codec(data);
  EXPECT_TRUE(codec.empty());
  EXPECT_EQ(codec.uncompressed_size(), 0u);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, SingleSymbolRoundTrips) {
  const std::vector<std::uint8_t> data(1000, 42);
  const PivCoHuffman codec(data);
  EXPECT_EQ(codec.uncompressed_size(), data.size());
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, TwoSymbolRoundTrips) {
  std::vector<std::uint8_t> data;
  for (std::size_t i = 0; i < 500; i++) {
    data.push_back(static_cast<std::uint8_t>(i % 2));
  }
  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, KnownInputRoundTrips) {
  const std::vector<std::uint8_t> data{3, 2, 0, 3, 1, 1, 2, 0, 3, 1};
  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, CompressedDataIsNonEmpty) {
  const std::vector<std::uint8_t> data{3, 2, 0, 3, 1, 1, 2, 0, 3, 1};
  const PivCoHuffman codec(data);
  EXPECT_GT(codec.compressed_size(), 0u);
  EXPECT_EQ(codec.compressed_data().size(), codec.compressed_size());
}

TEST(PivCoHuffmanSmoke, FullAlphabetRoundTrips) {
  std::vector<std::uint8_t> data(256);
  for (std::size_t i = 0; i < 256; i++) {
    data[i] = static_cast<std::uint8_t>(i);
  }
  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, RandomUniformRoundTrips) {
  std::mt19937 rng(239);
  std::uniform_int_distribution<int> dist(0, 255);
  std::vector<std::uint8_t> data(10000);
  for (auto& b : data) {
    b = static_cast<std::uint8_t>(dist(rng));
  }
  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, SkewedInputRoundTrips) {
  // Most symbols are one value, a few are outliers: exercises a deep, lopsided
  // Huffman tree without any single-symbol shortcut.
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> coin(0, 99);
  std::vector<std::uint8_t> data(10000, 7);
  for (std::size_t i = 0; i < data.size(); i++) {
    if (coin(rng) < 5) {
      data[i] = static_cast<std::uint8_t>(rng() % 256);
    }
  }
  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}
