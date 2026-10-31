#include <gtest/gtest.h>
#include <pixie/bits.h>
#include <pixie/rank_select/support.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/utils.h>

#include <limits>
#include <numeric>
#include <random>

#include "../benchmarks/shifted_query_pool.h"
#include "../benchmarks/tree_input.h"

TEST(TreeBenchmarkInputTest, RandomTreesMatchExistingEncodersExactly) {
  for (size_t n :
       {1, 2, 3, 7, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1024, 4097}) {
    for (uint64_t seed : {0, 1, 42, 12345}) {
      SCOPED_TRACE(testing::Message() << "n=" << n << " seed=" << seed);
      std::mt19937_64 rng(seed);
      const auto adjacency = generate_random_tree(n, rng);
      const auto ordered = dfs_order(n, adjacency);
      const auto bp = pixie::benchmarks::GenerateRandomBp(n, seed);
      EXPECT_EQ(bp, pixie::adj_to_bp(n, ordered));
      EXPECT_EQ(pixie::benchmarks::DfudsFromBp(bp, n),
                pixie::adj_to_dfuds(n, ordered));
      EXPECT_EQ(pixie::benchmarks::LoudsFromBp(bp, n),
                adj_to_louds(n, bfs_order(n, adjacency)));
    }
  }
}

TEST(TreeBenchmarkInputTest, ConversionsHandleDepthWidthAndWordBoundaries) {
  for (size_t n : {1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 4097}) {
    for (unsigned shape = 0; shape < 3; ++shape) {
      SCOPED_TRACE(testing::Message() << "n=" << n << " shape=" << shape);
      std::vector<std::vector<size_t>> adjacency(n);
      adjacency[0].push_back(0);
      for (size_t node = 1; node < n; ++node) {
        const size_t parent = shape == 0   ? node - 1
                              : shape == 1 ? 0
                                           : (node - 1) / 2;
        adjacency[node].push_back(parent);
        adjacency[parent].push_back(node);
      }
      const auto bp = pixie::adj_to_bp(n, adjacency);
      EXPECT_EQ(pixie::benchmarks::DfudsFromBp(bp, n),
                pixie::adj_to_dfuds(n, adjacency));
      EXPECT_EQ(pixie::benchmarks::LoudsFromBp(bp, n),
                adj_to_louds(n, bfs_order(n, adjacency)));
    }
  }
}

TEST(TreeBenchmarkInputTest, EmptyInputsAndGeneratorLimit) {
  EXPECT_TRUE(pixie::benchmarks::GenerateRandomBp(0).empty());
  EXPECT_TRUE(pixie::benchmarks::DfudsFromBp({}, 0).empty());
  EXPECT_TRUE(pixie::benchmarks::LoudsFromBp({}, 0).empty());
  EXPECT_THROW(
      pixie::benchmarks::GenerateRandomBp(std::numeric_limits<uint32_t>::max()),
      std::length_error);
}

TEST(RankSelectBenchmarkTest, Select10PercentFill) {
  for (size_t n = 8; n <= (1ull << 28); n <<= 2) {
    std::mt19937_64 rng(42);

    std::vector<uint64_t> bits(((8 + n / 64) / 8) * 8);
    size_t num_ones = n * 0.1;
    for (size_t i = 0; i < num_ones; i++) {
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
    for (size_t i = 0; i < num_ones; i++) {
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
    for (size_t i = 0; i < num_ones; i++) {
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
    for (size_t i = 0; i < num_ones; i++) {
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

TEST(ShiftedQueryPoolTest, MatchesWideArithmeticAcrossPoolWraps) {
  constexpr auto max = std::numeric_limits<std::size_t>::max();
  for (const auto& [first, last] :
       {std::pair<std::size_t, std::size_t>{0, 1024},
        {1, 1009},
        {7, 7},
        {max - 101, max},
        {0, max},
        {1, max}}) {
    constexpr std::uint64_t seed = 42;
    pixie::benchmarks::ShiftedQueryPool queries(first, last, seed);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::size_t> dist(0, last - first);
    std::vector<std::size_t> pool(
        pixie::benchmarks::ShiftedQueryPool::kQueryCount);
    for (auto& value : pool) {
      value = dist(rng);
    }
    std::uint64_t state = seed;
    std::size_t shift = 0;
    const auto range = static_cast<unsigned __int128>(last) - first + 1;
    for (std::size_t i = 0; i < 3 * pool.size() + 1; ++i) {
      if (i % pool.size() == 0) {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        shift = static_cast<std::size_t>(z ^ (z >> 31));
      }
      const auto expected =
          first + static_cast<std::size_t>(
                      (static_cast<unsigned __int128>(pool[i % pool.size()]) +
                       shift + i) %
                      range);
      ASSERT_EQ(queries.next(), expected) << first << ", " << last << ", " << i;
    }
  }
}

TEST(ShiftedQueryPoolTest, SmallPoolDoesNotReplayAndRejectsEmptyRange) {
  using Pool = pixie::benchmarks::ShiftedQueryPool;
  static_assert(Pool::kPoolBytes == 64 * 1024);
  Pool a(1, std::size_t{1} << 30, 42);
  Pool b(1, std::size_t{1} << 30, 42);
  std::vector<std::size_t> first_pass(Pool::kQueryCount);
  for (auto& query : first_pass) {
    query = a.next();
    ASSERT_EQ(query, b.next());
  }
  std::size_t changed = 0;
  for (const auto previous : first_pass) {
    const auto next = a.next();
    ASSERT_EQ(next, b.next());
    changed += next != previous;
  }
  EXPECT_GT(changed, 0u);
  EXPECT_THROW((Pool(2, 1, 42)), std::invalid_argument);
}
