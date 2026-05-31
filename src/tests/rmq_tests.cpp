#include <gtest/gtest.h>
#include <pixie/rmq.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <random>
#include <span>
#include <vector>

namespace {

template <class T, class Compare>
std::size_t naive_arg_min(std::span<const T> values,
                          std::size_t left,
                          std::size_t right,
                          Compare compare) {
  if (left > right || right >= values.size()) {
    return pixie::rmq::SparseTable<T, Compare>::npos;
  }
  std::size_t best = left;
  for (std::size_t i = left + 1; i <= right; ++i) {
    if (compare(values[i], values[best])) {
      best = i;
    }
  }
  return best;
}

template <class Rmq, class T, class Compare>
void check_all_ranges(const Rmq& rmq,
                      std::span<const T> values,
                      Compare compare) {
  ASSERT_EQ(rmq.size(), values.size());
  for (std::size_t left = 0; left < values.size(); ++left) {
    for (std::size_t right = left; right < values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << "]";
      EXPECT_EQ(rmq.range_min(left, right), values[expected])
          << "range=[" << left << "," << right << "]";
    }
  }
}

template <class Rmq, class T, class Compare>
void check_all_arg_min_ranges(const Rmq& rmq,
                              std::span<const T> values,
                              Compare compare) {
  ASSERT_EQ(rmq.size(), values.size());
  for (std::size_t left = 0; left < values.size(); ++left) {
    for (std::size_t right = left; right < values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << "]";
    }
  }
}

std::vector<std::uint64_t> pack_depth_deltas(
    std::span<const std::int64_t> depths) {
  std::vector<std::uint64_t> bits((depths.size() - 1 + 63) / 64, 0);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    if (depths[i] - depths[i - 1] == 1) {
      bits[(i - 1) >> 6] |= std::uint64_t{1} << ((i - 1) & 63);
    }
  }
  return bits;
}

}  // namespace

TEST(RmqSparseTable, ExhaustiveSmallArray) {
  const std::vector<int> values = {4, 1, 3, 1, 5, 0, 0, 2};
  const pixie::rmq::SparseTable<int> rmq(values);
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

TEST(RmqSegmentTree, ExhaustiveSmallArray) {
  const std::vector<int> values = {4, 1, 3, 1, 5, 0, 0, 2};
  const pixie::rmq::SegmentTree<int> rmq(values);
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

TEST(RmqBpPlusMinusOne, ExhaustiveSmallDepthArray) {
  const std::vector<std::int64_t> depths = {0, 1, 0, 1, 2, 1, 2, 1, 0};
  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  check_all_arg_min_ranges(rmq, std::span<const std::int64_t>(depths),
                           std::less<std::int64_t>());
}

TEST(RmqBpPlusMinusOne, CrossBlockRanges) {
  std::vector<std::int64_t> depths(385);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 7 == 0) || (i % 7 == 1) || (i % 7 == 4);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  check_all_arg_min_ranges(rmq, std::span<const std::int64_t>(depths),
                           std::less<std::int64_t>());
}

TEST(RmqBpPlusMinusOne, BoundaryRangesAround128PositionBlocks) {
  std::vector<std::int64_t> depths(260);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] = depths[i - 1] + ((i % 5 == 0 || i % 11 == 0) ? 1 : -1);
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 0},     {126, 127}, {127, 128}, {128, 128},
      {128, 255}, {129, 255}, {255, 256}, {0, 259},
  };

  for (const auto [left, right] : ranges) {
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                            std::less<std::int64_t>()))
        << "range=[" << left << "," << right << "]";
  }
}

TEST(RmqBpPlusMinusOne, LongSequenceRangesNearEnd) {
  std::vector<std::int64_t> depths(8193);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 13 == 0) || (i % 17 == 0) || (i % 19 == 0);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {depths.size() - 1, depths.size() - 1},
      {depths.size() - 128, depths.size() - 1},
      {depths.size() - 513, depths.size() - 3},
      {depths.size() - 4097, depths.size() - 7},
  };

  for (const auto [left, right] : ranges) {
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                            std::less<std::int64_t>()))
        << "range=[" << left << "," << right << "]";
  }
}

TEST(RmqBpPlusMinusOne, CrossBlockTieKeepsFirstPosition) {
  std::vector<std::int64_t> depths(384);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] = (i % 2 == 0) ? 0 : 1;
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  const std::size_t left = 120;
  const std::size_t right = 260;

  EXPECT_EQ(rmq.arg_min(left, right), left);
  EXPECT_EQ(rmq.arg_min(left, right),
            naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                          std::less<std::int64_t>()));
}

TEST(RmqBpPlusMinusOne, DisjointBoundaryRangeMatchesNaive) {
  std::vector<std::int64_t> depths(384);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 9 == 0) || (i % 9 == 3) || (i % 11 == 0);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {96, 160}, {120, 140}, {127, 129}, {190, 258}, {250, 260},
  };

  for (const auto [left, right] : ranges) {
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                            std::less<std::int64_t>()))
        << "range=[" << left << "," << right << "]";
  }
}

TEST(RmqBpPlusMinusOne, RejectsTooSmallBitSpan) {
  const std::vector<std::uint64_t> bits;
  EXPECT_THROW((pixie::rmq::BpPlusMinusOneRmq<>(bits, 2)),
               std::invalid_argument);
}

TEST(RmqCartesianTree, ExhaustiveSmallArray) {
  const std::vector<int> values = {4, 1, 3, 1, 5, 0, 0, 2};
  const pixie::rmq::CartesianTreeRmq<int> rmq(values);
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

TEST(Rmq, FirstMinimumTieBreaking) {
  const std::vector<int> values = {7, 2, 2, 3, 2};
  const pixie::rmq::SparseTable<int> sparse(values);
  const pixie::rmq::SegmentTree<int> segment(values);
  const pixie::rmq::CartesianTreeRmq<int> cartesian(values);

  EXPECT_EQ(sparse.arg_min(0, 4), 1u);
  EXPECT_EQ(segment.arg_min(0, 4), 1u);
  EXPECT_EQ(cartesian.arg_min(0, 4), 1u);
  EXPECT_EQ(sparse.arg_min(2, 4), 2u);
  EXPECT_EQ(segment.arg_min(2, 4), 2u);
  EXPECT_EQ(cartesian.arg_min(2, 4), 2u);
}

TEST(Rmq, InvalidAndEmptyRanges) {
  const std::vector<int> values = {3, 1, 2};
  const pixie::rmq::SparseTable<int> sparse(values);
  const pixie::rmq::SegmentTree<int> segment(values);
  const pixie::rmq::CartesianTreeRmq<int> cartesian(values);

  EXPECT_EQ(sparse.arg_min(2, 1), pixie::rmq::SparseTable<int>::npos);
  EXPECT_EQ(segment.arg_min(2, 1), pixie::rmq::SegmentTree<int>::npos);
  EXPECT_EQ(cartesian.arg_min(2, 1), pixie::rmq::CartesianTreeRmq<int>::npos);
  EXPECT_EQ(sparse.arg_min(0, values.size()),
            pixie::rmq::SparseTable<int>::npos);
  EXPECT_EQ(segment.arg_min(0, values.size()),
            pixie::rmq::SegmentTree<int>::npos);
  EXPECT_EQ(cartesian.arg_min(0, values.size()),
            pixie::rmq::CartesianTreeRmq<int>::npos);
  EXPECT_EQ(sparse.range_min(2, 1), 0);
  EXPECT_EQ(segment.range_min(2, 1), 0);
  EXPECT_EQ(cartesian.range_min(2, 1), 0);

  const std::vector<int> empty;
  const pixie::rmq::SparseTable<int> empty_sparse(empty);
  const pixie::rmq::SegmentTree<int> empty_segment(empty);
  const pixie::rmq::CartesianTreeRmq<int> empty_cartesian(empty);
  EXPECT_TRUE(empty_sparse.empty());
  EXPECT_TRUE(empty_segment.empty());
  EXPECT_TRUE(empty_cartesian.empty());
  EXPECT_EQ(empty_sparse.arg_min(0, 0), pixie::rmq::SparseTable<int>::npos);
  EXPECT_EQ(empty_segment.arg_min(0, 0), pixie::rmq::SegmentTree<int>::npos);
  EXPECT_EQ(empty_cartesian.arg_min(0, 0),
            pixie::rmq::CartesianTreeRmq<int>::npos);
}

TEST(Rmq, ComparatorCanSelectMaximum) {
  const std::vector<int> values = {1, 8, 3, 8, 4};
  const pixie::rmq::SparseTable<int, std::greater<int>> sparse(values);
  const pixie::rmq::SegmentTree<int, std::greater<int>> segment(values);
  const pixie::rmq::CartesianTreeRmq<int, std::greater<int>> cartesian(values);

  check_all_ranges(sparse, std::span<const int>(values), std::greater<int>());
  check_all_ranges(segment, std::span<const int>(values), std::greater<int>());
  check_all_ranges(cartesian, std::span<const int>(values),
                   std::greater<int>());
  EXPECT_EQ(sparse.arg_min(0, 4), 1u);
  EXPECT_EQ(segment.arg_min(0, 4), 1u);
  EXPECT_EQ(cartesian.arg_min(0, 4), 1u);
}

TEST(RmqCartesianTree, MonotoneArrays) {
  const std::vector<int> increasing = {1, 2, 3, 4, 5, 6};
  const std::vector<int> decreasing = {6, 5, 4, 3, 2, 1};
  const pixie::rmq::CartesianTreeRmq<int> increasing_rmq(increasing);
  const pixie::rmq::CartesianTreeRmq<int> decreasing_rmq(decreasing);

  check_all_ranges(increasing_rmq, std::span<const int>(increasing),
                   std::less<int>());
  check_all_ranges(decreasing_rmq, std::span<const int>(decreasing),
                   std::less<int>());
}

TEST(RmqCartesianTree, CopyAndMoveRebuildInternalSpans) {
  const std::vector<int> values = {5, 4, 3, 2, 1, 2, 3};
  const pixie::rmq::CartesianTreeRmq<int> original(values);
  pixie::rmq::CartesianTreeRmq<int> copied(original);
  pixie::rmq::CartesianTreeRmq<int> assigned;
  assigned = copied;
  pixie::rmq::CartesianTreeRmq<int> moved(std::move(copied));

  check_all_ranges(original, std::span<const int>(values), std::less<int>());
  check_all_ranges(assigned, std::span<const int>(values), std::less<int>());
  check_all_ranges(moved, std::span<const int>(values), std::less<int>());
}

TEST(RmqCartesianTree, EulerDepthsArePlusMinusOne) {
  const std::vector<int> values = {4, 1, 3, 2, 5};
  const pixie::rmq::CartesianTreeRmq<int> rmq(values);
  const auto depths = rmq.euler_depths();

  ASSERT_FALSE(depths.empty());
  for (std::size_t i = 1; i < depths.size(); ++i) {
    EXPECT_EQ(std::abs(depths[i] - depths[i - 1]), 1);
  }
}

TEST(Rmq, DifferentialRandom) {
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> value_dist(-50, 50);
  for (std::size_t size = 1; size <= 257; size += 17) {
    std::vector<int> values(size);
    std::generate(values.begin(), values.end(),
                  [&] { return value_dist(rng); });

    const pixie::rmq::SparseTable<int> sparse(values);
    const pixie::rmq::SegmentTree<int> segment(values);
    const pixie::rmq::CartesianTreeRmq<int> cartesian(values);
    check_all_ranges(sparse, std::span<const int>(values), std::less<int>());
    check_all_ranges(segment, std::span<const int>(values), std::less<int>());
    check_all_ranges(cartesian, std::span<const int>(values), std::less<int>());
  }
}

TEST(RmqBpPlusMinusOne, DifferentialRandomWalks) {
  std::mt19937_64 rng(77);
  for (std::size_t size = 1; size <= 257; size += 17) {
    std::vector<std::int64_t> depths(size);
    for (std::size_t i = 1; i < depths.size(); ++i) {
      depths[i] = depths[i - 1] + ((rng() & 1u) != 0 ? 1 : -1);
    }

    const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
    const pixie::rmq::BpPlusMinusOneRmq<> rmq(bits, depths.size());
    check_all_arg_min_ranges(rmq, std::span<const std::int64_t>(depths),
                             std::less<std::int64_t>());
  }
}
