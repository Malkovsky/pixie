#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "bits.h"
#include "bitvector.h"

using pixie::BitVector;
using pixie::BitVectorInterleaved;

TEST(BitVectorBenchmarkTest, SelectNonInterleaved10PercentFill) {
  for (size_t n = 4; n <= (1ull << 34); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(1 + n / 64);
    size_t num_ones = n * 0.1;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::BitVector bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}

TEST(BitVectorBenchmarkTest, SelectNonInterleaved90PercentFill) {
  for (size_t n = 4; n <= (1ull << 34); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(1 + n / 64);
    size_t num_ones = n * 0.9;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::BitVector bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}

TEST(BitVectorBenchmarkTest, SelectNonInterleaved) {
  for (size_t n = 4; n <= (1ull << 34); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(1 + n / 64);
    for (auto& x : bits) {
      x = rng();
    }
    pixie::BitVector bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;

    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}
