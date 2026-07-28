#include <gtest/gtest.h>
#include <pixie/huffman/implementations.h>

#include <algorithm>
#include <array>
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

TEST(PivCoHuffmanSmoke, FlatPowerOfTwoAlphabetsRoundTrip) {
  std::mt19937 rng(1871);
  for (std::uint8_t depth = 1; depth <= 8; ++depth) {
    const std::size_t alphabet_size = std::size_t{1} << depth;
    std::vector<std::uint8_t> data;
    data.reserve(alphabet_size * 16);
    for (std::size_t repetition = 0; repetition < 16; ++repetition) {
      for (std::size_t symbol = 0; symbol < alphabet_size; ++symbol) {
        data.push_back(static_cast<std::uint8_t>(symbol));
      }
    }
    std::shuffle(data.begin(), data.end(), rng);

    const PivCoHuffman codec(data);
    EXPECT_EQ(decode_same(codec), data) << "flat depth " << +depth;
    EXPECT_EQ(decode_via_bytes(codec), data) << "flat depth " << +depth;
  }
}

TEST(PivCoSimd, FlatKernelsDecodeFullGroupsAndTails) {
  std::mt19937 rng(4189);
  for (std::uint8_t depth = 1; depth <= 8; ++depth) {
    const std::size_t alphabet_size = std::size_t{1} << depth;
    const std::size_t count = 129 + depth;
    std::vector<std::uint8_t> table(alphabet_size);
    for (std::size_t i = 0; i < alphabet_size; ++i) {
      table[i] = static_cast<std::uint8_t>(i);
    }
    std::shuffle(table.begin(), table.end(), rng);

    std::vector<std::uint64_t> bits((count * depth + 63) / 64, 0);
    std::vector<std::uint8_t> expected(count);
    for (std::size_t i = 0; i < count; ++i) {
      const std::uint64_t code = rng() % alphabet_size;
      const std::size_t bit_position = i * depth;
      const std::size_t word_index = bit_position / 64;
      const std::size_t bit_offset = bit_position % 64;
      bits[word_index] |= code << bit_offset;
      if (bit_offset + depth > 64) {
        bits[word_index + 1] |= code >> (64 - bit_offset);
      }
      expected[i] = table[code];
    }

    std::vector<std::uint8_t> decoded(count);
    pixie::pivco_flat_decode(decoded.data(), bits.data(), table.data(), count,
                             depth);
    EXPECT_EQ(decoded, expected) << "flat depth " << +depth;
  }
}

TEST(PivCoSimd, MergeVariantsMatchReferenceAcrossFullGroupsAndTail) {
  constexpr std::size_t kCount = 263;
  constexpr std::size_t kPadding = 64;
  std::vector<std::uint64_t> bits((kCount + 63) / 64, 0);
  std::vector<std::uint8_t> left(kCount + kPadding, 0);
  std::vector<std::uint8_t> right(kCount + kPadding, 0);
  std::vector<std::uint8_t> expected(kCount);
  std::size_t left_count = 0;
  std::size_t right_count = 0;
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool goes_right = ((i * 37 + i / 5 + 11) % 13) < 6;
    if (goes_right) {
      bits[i / 64] |= std::uint64_t{1} << (i % 64);
      right[right_count] = static_cast<std::uint8_t>(0x80u + right_count % 97);
      expected[i] = right[right_count++];
    } else {
      left[left_count] = static_cast<std::uint8_t>(left_count % 113);
      expected[i] = left[left_count++];
    }
  }

  std::vector<std::uint8_t> decoded(kCount + kPadding, 0);
  pixie::pivco_merge_decode(decoded.data(), left.data(), right.data(),
                            bits.data(), kCount);
  EXPECT_TRUE(std::equal(expected.begin(), expected.end(), decoded.begin()));

  constexpr std::uint8_t kLeftSymbol = 17;
  constexpr std::uint8_t kRightSymbol = 231;
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool goes_right = ((bits[i / 64] >> (i % 64)) & 1u) != 0;
    expected[i] = goes_right ? kRightSymbol : kLeftSymbol;
  }
  std::fill(decoded.begin(), decoded.end(), 0);
  pixie::pivco_merge_decode_cst_cst(decoded.data(), kLeftSymbol, kRightSymbol,
                                    bits.data(), kCount);
  EXPECT_TRUE(std::equal(expected.begin(), expected.end(), decoded.begin()));

  right_count = 0;
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool goes_right = ((bits[i / 64] >> (i % 64)) & 1u) != 0;
    expected[i] = goes_right ? right[right_count++] : kLeftSymbol;
  }
  std::fill(decoded.begin(), decoded.end(), 0);
  pixie::pivco_merge_decode_cst_vec(decoded.data(), kLeftSymbol, right.data(),
                                    bits.data(), kCount);
  EXPECT_TRUE(std::equal(expected.begin(), expected.end(), decoded.begin()));

  left_count = 0;
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool goes_right = ((bits[i / 64] >> (i % 64)) & 1u) != 0;
    expected[i] = goes_right ? kRightSymbol : left[left_count++];
  }
  std::fill(decoded.begin(), decoded.end(), 0);
  pixie::pivco_merge_decode_vec_cst(decoded.data(), left.data(), kRightSymbol,
                                    bits.data(), kCount);
  EXPECT_TRUE(std::equal(expected.begin(), expected.end(), decoded.begin()));
}

TEST(PivCoSimd, PartitionVariantsBuildTheSameDirectionMask) {
  constexpr std::size_t kCount = 257;
  std::array<std::uint8_t, 256> directions{};
  for (std::size_t symbol = 0; symbol < directions.size(); ++symbol) {
    directions[symbol] =
        static_cast<std::uint8_t>(((symbol * 73) ^ (symbol >> 2)) & 1u);
  }

  std::mt19937 rng(7127);
  std::vector<std::uint8_t> input(kCount);
  std::vector<std::uint8_t> expected_left;
  std::vector<std::uint8_t> expected_right;
  std::vector<std::uint64_t> expected_bits((kCount + 63) / 64, 0);
  for (std::size_t i = 0; i < input.size(); ++i) {
    input[i] = static_cast<std::uint8_t>(rng());
    if (directions[input[i]]) {
      expected_bits[i / 64] |= std::uint64_t{1} << (i % 64);
      expected_right.push_back(input[i]);
    } else {
      expected_left.push_back(input[i]);
    }
  }

  std::vector<std::uint8_t> left(expected_left.size());
  std::vector<std::uint8_t> right(expected_right.size());
  std::vector<std::uint64_t> full_bits(expected_bits.size(), 0);
  pixie::pivco_partition_encode(left.data(), right.data(), full_bits.data(),
                                input.data(), directions.data(), input.size(),
                                left.size(), right.size());
  EXPECT_EQ(left, expected_left);
  EXPECT_EQ(right, expected_right);
  EXPECT_EQ(full_bits, expected_bits);

  std::vector<std::uint64_t> left_bits(expected_bits.size(), 0);
  std::fill(left.begin(), left.end(), 0);
  pixie::pivco_partition_encode_left(left.data(), left_bits.data(),
                                     input.data(), directions.data(),
                                     input.size(), left.size());
  EXPECT_EQ(left, expected_left);
  EXPECT_EQ(left_bits, expected_bits);

  std::vector<std::uint64_t> right_bits(expected_bits.size(), 0);
  std::fill(right.begin(), right.end(), 0);
  pixie::pivco_partition_encode_right(right.data(), right_bits.data(),
                                      input.data(), directions.data(),
                                      input.size(), right.size());
  EXPECT_EQ(right, expected_right);
  EXPECT_EQ(right_bits, expected_bits);

  std::vector<std::uint64_t> bitmap_only(expected_bits.size(), 0);
  pixie::pivco_partition_encode_none(bitmap_only.data(), input.data(),
                                     directions.data(), input.size());
  EXPECT_EQ(bitmap_only, expected_bits);
}

TEST(PivCoSimd, ConstantLeafBitmapMatchesReferenceAcrossTail) {
  constexpr std::size_t kCount = 267;
  constexpr std::uint8_t kRightSymbol = 193;
  std::vector<std::uint8_t> input(kCount);
  std::vector<std::uint64_t> expected_bits((kCount + 63) / 64, 0);
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool is_right = ((i * 19 + i / 3) % 11) < 5;
    input[i] = is_right ? kRightSymbol : static_cast<std::uint8_t>(i % 127);
    if (is_right) {
      expected_bits[i / 64] |= std::uint64_t{1} << (i % 64);
    }
  }

  std::vector<std::uint64_t> bits(expected_bits.size(), 0);
  pixie::pivco_partition_encode_cst_cst(bits.data(), input.data(), kRightSymbol,
                                        input.size());
  EXPECT_EQ(bits, expected_bits);
}

#if (defined(__AVX2__) && defined(__SSE4_1__)) || \
    defined(PIXIE_PIVCO_NEON_SUPPORT)
TEST(PivCoSimd, OneConstantPartitionsMatchReferenceAcrossTail) {
  constexpr std::size_t kCount = 269;
  constexpr std::uint8_t kConstant = 241;
  std::vector<std::uint8_t> input(kCount);
  std::vector<std::uint8_t> expected_dense;
  std::vector<std::uint64_t> expected_bits((kCount + 63) / 64, 0);
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool is_constant = ((i * 23 + i / 7) % 17) < 8;
    input[i] = is_constant ? kConstant : static_cast<std::uint8_t>(i % 191);
    if (!is_constant) {
      expected_bits[i / 64] |= std::uint64_t{1} << (i % 64);
      expected_dense.push_back(input[i]);
    }
  }

  std::vector<std::uint8_t> dense(expected_dense.size());
  std::vector<std::uint64_t> bits(expected_bits.size(), 0);
  pixie::pivco_partition_encode_cst_vec(dense.data(), bits.data(), input.data(),
                                        kConstant, input.size(), dense.size());
  EXPECT_EQ(dense, expected_dense);
  EXPECT_EQ(bits, expected_bits);

  std::fill(expected_bits.begin(), expected_bits.end(), 0);
  expected_dense.clear();
  for (std::size_t i = 0; i < kCount; ++i) {
    const bool is_constant = input[i] == kConstant;
    if (is_constant) {
      expected_bits[i / 64] |= std::uint64_t{1} << (i % 64);
    } else {
      expected_dense.push_back(input[i]);
    }
  }
  std::fill(dense.begin(), dense.end(), 0);
  std::fill(bits.begin(), bits.end(), 0);
  pixie::pivco_partition_encode_vec_cst(dense.data(), bits.data(), input.data(),
                                        kConstant, input.size(), dense.size());
  EXPECT_EQ(dense, expected_dense);
  EXPECT_EQ(bits, expected_bits);
}
#endif

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

TEST(PivCoHuffmanSmoke, LengthLimitedTreeRoundTrips) {
  // Fibonacci weights create a maximally deep unconstrained Huffman tree and
  // exercise the 15-bit length correction.
  std::vector<std::uint8_t> data;
  std::size_t previous = 1;
  std::size_t current = 1;
  for (std::uint8_t symbol = 0; symbol < 17; ++symbol) {
    const std::size_t frequency = symbol < 2 ? 1 : previous + current;
    if (symbol >= 2) {
      previous = current;
      current = frequency;
    }
    data.insert(data.end(), frequency, symbol);
  }

  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, MultipleBlocksRoundTrip) {
  std::mt19937 rng(9127);
  std::uniform_int_distribution<int> dist(0, 255);
  std::vector<std::uint8_t> data(3 * 64 * 1024 + 137);
  for (auto& byte : data) {
    byte = static_cast<std::uint8_t>(dist(rng));
  }

  const PivCoHuffman codec(data);
  EXPECT_EQ(decode_same(codec), data);
  EXPECT_EQ(decode_via_bytes(codec), data);
}

TEST(PivCoHuffmanSmoke, SharedTreeHandlesBlockLocalMissingSymbols) {
  constexpr std::size_t kBlockSize = 64 * 1024;
  std::vector<std::uint8_t> data;
  data.reserve(3 * kBlockSize);
  data.insert(data.end(), kBlockSize, 17);
  data.insert(data.end(), kBlockSize, 93);
  data.insert(data.end(), kBlockSize, 241);

  const PivCoHuffman codec(data);
  EXPECT_EQ(codec.decode(), data);
  EXPECT_EQ(codec.decode(), data);

  const PivCoHuffman copy(codec.compressed_data());
  EXPECT_EQ(copy.decode(), data);
  EXPECT_EQ(copy.decode(), data);
}

TEST(PivCoHuffmanSmoke, RandomizedAlphabetAndBlockSizesRoundTrip) {
  std::mt19937 rng(3319);
  for (std::size_t test_case = 0; test_case < 32; ++test_case) {
    const std::size_t alphabet_size = 1 + rng() % 256;
    const std::size_t size = 1 + rng() % (2 * 64 * 1024);
    std::vector<std::uint8_t> data(size);
    for (auto& byte : data) {
      byte = static_cast<std::uint8_t>(rng() % alphabet_size);
    }

    const PivCoHuffman codec(data);
    EXPECT_EQ(decode_same(codec), data) << "test case " << test_case;
    EXPECT_EQ(decode_via_bytes(codec), data) << "test case " << test_case;
  }
}
