#include <gtest/gtest.h>
#include <pixie/bits.h>
#include <pixie/rank_select/implementations.h>

#include <numeric>
#include <random>

TEST(RankSelectBenchmarkTest, Select10PercentFill) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    size_t num_ones = n * 0.1;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}

TEST(RankSelectBenchmarkTest, SelectZero10PercentFill) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    size_t num_ones = n * 0.1;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank0 = bv.rank0(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank0 = rng() % max_rank0;
      bv.select0(rank0);
    }
  }
}

TEST(RankSelectBenchmarkTest, Select90PercentFill) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    size_t num_ones = n * 0.9;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}

TEST(RankSelectBenchmarkTest, SelectZero90PercentFill) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    size_t num_ones = n * 0.9;
    for (int i = 0; i < num_ones; i++) {
      uint64_t pos = rng() % n;
      bits[pos / 64] |= (1ULL << pos % 64);
    }

    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank0 = bv.rank0(bv.size()) + 1;
    for (int i = 0; i < 100000; i++) {
      uint64_t rank0 = rng() % max_rank0;
      bv.select0(rank0);
    }
  }
}

TEST(RankSelectBenchmarkTest, Select) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    for (auto& x : bits) {
      x = rng();
    }
    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank = bv.rank(bv.size()) + 1;

    for (int i = 0; i < 100000; i++) {
      uint64_t rank = rng() % max_rank;
      bv.select(rank);
    }
  }
}

TEST(RankSelectBenchmarkTest, SelectZero) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    for (auto& x : bits) {
      x = rng();
    }
    pixie::RankSelectSupport<> bv(bits, n);

    auto max_rank0 = bv.rank0(bv.size()) + 1;

    for (int i = 0; i < 100000; i++) {
      uint64_t rank0 = rng() % max_rank0;
      bv.select0(rank0);
    }
  }
}
