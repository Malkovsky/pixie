#include <numeric>
#include <random>

#include "bits.h"
#include "bitvector.h"
#include <gtest/gtest.h>

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
    for (size_t i = 0; i <= 512; ++i) {
      auto p = rank_512(a.data(), i);
      ASSERT_EQ(p, rank);
      rank += 1 & (a[i >> 6] >> (i & 63));
    }
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
    for (auto &x : a) {
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
    for (auto &x : a) {
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
  std::vector<uint64_t> bits = {0b101010101};
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
  std::vector<uint64_t> bits = {0xAAAAAAAAAAAAAAAA, 0x5555555555555555};
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
  std::vector<uint64_t> bits = {0b10101010101};
  BitVector bv(bits, 11);
  EXPECT_EQ(bv.to_string(), "10101010101");
}

TEST(BitVectorTest, RankBasic) {
  std::vector<uint64_t> bits = {0b10110};
  BitVector bv(bits, 5);

  EXPECT_EQ(bv.rank(0), 0); // No bits
  EXPECT_EQ(bv.rank(1), 0); // 0
  EXPECT_EQ(bv.rank(2), 1); // 10
  EXPECT_EQ(bv.rank(3), 2); // 110
  EXPECT_EQ(bv.rank(4), 2); // 0110
  EXPECT_EQ(bv.rank(5), 3); // 10110
}

TEST(BitVectorTest, RankWithZeros) {
  std::vector<uint64_t> bits = {0};
  BitVector bv(bits, 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(bv.rank(i), 0);
  }
}

TEST(BitVectorTest, SelectBasic) {
  std::vector<uint64_t> bits = {0b1100010110010110};
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

TEST(BitVectorTest, BenchmarkSelectTest) {
  for (size_t n = 4; n <= (1ull << 34); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(1 + n / 64);
    for (auto &x : bits) {
      x = rng();
    }
    pixie::BitVector bv(bits, n);

    auto max_rank = bv.rank(bv.size());
    // size_t rank = 0;
    auto r = bv.rank(bv.size());
    auto s = bv.select(r);

    // for (size_t i = 0; i < 10000 && r - i > 0; ++i) {
    //   s = bv.select(r - i);
    // }

    for (size_t i = 0; i < 20000000; ++i) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
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

// TEST(BitVectorTest, SelectBasic) {
//     std::vector<uint64_t> bits = {0b10110};
//     BitVector bv(bits, 5);

//     EXPECT_EQ(bv.select(0), 0);  // Position of the 1st 1-bit
//     EXPECT_EQ(bv.select(1), 2);  // Position of the 2nd 1-bit
//     EXPECT_EQ(bv.select(2), 3);  // Position of the 3rd 1-bit

//     EXPECT_THROW(bv.select(3), std::out_of_range);  // Only 3 bits set
// }

// TEST(BitVectorTest, SelectWithZeros) {
//     std::vector<uint64_t> bits = {0};
//     BitVector bv(bits, 5);

//     EXPECT_THROW(bv.select(0), std::out_of_range);  // No bits set
// }

// TEST(BitVectorTest, SelectWithAllOnes) {
//     std::vector<uint64_t> bits = {0b11111};
//     BitVector bv(bits, 5);

//     for (size_t i = 0; i < 5; i++) {
//         EXPECT_EQ(bv.select(i), i);  // i-th 1-bit is at position i
//     }
// }

// TEST(BitVectorTest, SelectMultipleWords) {
//     // Create a bitvector with bits set at every third position
//     std::vector<uint64_t> bits(10, 0);  // 10 words, 640 bits total
//     size_t num_bits = 10 * 64;

//     for (size_t i = 0; i < num_bits; i++) {
//         if (i % 3 == 0) {
//             size_t word_idx = i / 64;
//             size_t bit_offset = i % 64;
//             bits[word_idx] |= (1ULL << bit_offset);
//         }
//     }

//     BitVector bv(bits, num_bits);

//     // Test select for positions where bits are set
//     for (size_t i = 0; i < num_bits / 3; i++) {
//         EXPECT_EQ(bv.select(i), i * 3);
//     }
// }

// TEST(BitVectorTest, CrossBlockBoundaries) {
//     // Create a bitvector with specific bits set to test block boundaries
//     std::vector<uint64_t> bits(32, 0);  // 32 words = 2048 bits
//     size_t num_bits = 32 * 64;

//     // Set bits at basic block boundaries (multiples of 512 bits)
//     for (size_t i = 0; i < num_bits; i += 512) {
//         size_t word_idx = i / 64;
//         size_t bit_offset = i % 64;
//         bits[word_idx] |= (1ULL << bit_offset);
//     }

//     // Set bits at the end of each basic block (511, 1023, etc.)
//     for (size_t i = 511; i < num_bits; i += 512) {
//         size_t word_idx = i / 64;
//         size_t bit_offset = i % 64;
//         bits[word_idx] |= (1ULL << bit_offset);
//     }

//     BitVector bv(bits, num_bits);

//     // Test rank at basic block boundaries
//     for (size_t i = 0; i <= 4; i++) {
//         EXPECT_EQ(bv.rank(i * 512), i * 2);
//     }

//     // Test select for these specific bits
//     for (size_t i = 0; i < 8; i++) {
//         if (i % 2 == 0) {
//             EXPECT_EQ(bv.select(i), (i / 2) * 512);
//         } else {
//             EXPECT_EQ(bv.select(i), (i / 2) * 512 + 511);
//         }
//     }
// }

// TEST(BitVectorTest, RandomPattern) {
//     // Create a bitvector with a random pattern to test rank and select
//     std::vector<uint64_t> bits(16, 0);  // 16 words = 1024 bits
//     size_t num_bits = 1000;  // Use 1000 bits (not the full 1024)

//     std::mt19937 rng(42);  // Fixed seed for reproducibility

//     // Set random bits
//     for (size_t i = 0; i < num_bits; i++) {
//         if (rng() % 2 == 1) {  // ~50% chance of setting each bit
//             size_t word_idx = i / 64;
//             size_t bit_offset = i % 64;
//             bits[word_idx] |= (1ULL << bit_offset);
//         }
//     }

//     BitVector bv(bits, num_bits);

//     // Count the number of set bits up to each position
//     std::vector<size_t> expected_ranks(num_bits + 1, 0);
//     for (size_t i = 0; i < num_bits; i++) {
//         expected_ranks[i+1] = expected_ranks[i] + (bv[i] ? 1 : 0);
//     }

//     // Test rank for all positions
//     for (size_t i = 0; i <= num_bits; i++) {
//         EXPECT_EQ(bv.rank(i), expected_ranks[i]);
//     }

//     // Test select for all set bits
//     size_t total_ones = expected_ranks[num_bits];
//     std::vector<size_t> one_positions;
//     for (size_t i = 0; i < num_bits; i++) {
//         if (bv[i]) {
//             one_positions.push_back(i);
//         }
//     }

//     for (size_t i = 0; i < total_ones; i++) {
//         EXPECT_EQ(bv.select(i), one_positions[i]);
//     }
// }

// TEST(BitVectorTest, LargeBitVector) {
//     // Test with a larger bitvector to ensure the two-level structure works
//     correctly size_t num_words = 1024;  // 1024 words = 65536 bits (matches
//     kSuperBlockSize) std::vector<uint64_t> bits(num_words, 0); size_t
//     num_bits = num_words * 64;

//     // Set every 100th bit
//     for (size_t i = 0; i < num_bits; i += 100) {
//         size_t word_idx = i / 64;
//         size_t bit_offset = i % 64;
//         bits[word_idx] |= (1ULL << bit_offset);
//     }

//     BitVector bv(bits, num_bits);

//     // Test rank at multiples of 100
//     for (size_t i = 0; i <= 650; i++) {
//         EXPECT_EQ(bv.rank(i * 100), i);
//     }

//     // Test select
//     for (size_t i = 0; i < 650; i++) {
//         EXPECT_EQ(bv.select(i), i * 100);
//     }
// }

// // Tests for our BitVector implementation
// TEST(BitVectorTest2, EmptyConstructor) {
//     std::vector<uint64_t> bits = {0b101010101};
//     BitVector bv(bits, 8);
//     EXPECT_EQ(bv.size(), 8);
// }
