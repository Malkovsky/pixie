#include <gtest/gtest.h>
#include <pixie/rmq.h>
#include <pixie/rmq/experimental/node_euler_btree_rmq.h>
#include <pixie/rmq/node_euler_btree_rmq.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {

template <class T, class Compare>
std::size_t naive_arg_min(std::span<const T> values,
                          std::size_t left,
                          std::size_t right,
                          Compare compare) {
  if (left >= right || right > values.size()) {
    return pixie::rmq::SparseTable<T, Compare>::npos;
  }
  std::size_t best = left;
  for (std::size_t i = left + 1; i < right; ++i) {
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
    for (std::size_t right = left + 1; right <= values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
      EXPECT_EQ(rmq.range_min(left, right), values[expected])
          << "range=[" << left << "," << right << ")";
    }
  }
}

template <class Rmq, class T, class Compare>
void check_all_arg_min_ranges(const Rmq& rmq,
                              std::span<const T> values,
                              Compare compare) {
  ASSERT_EQ(rmq.size(), values.size());
  for (std::size_t left = 0; left < values.size(); ++left) {
    for (std::size_t right = left + 1; right <= values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

std::vector<std::uint64_t> pack_depth_deltas(
    std::span<const std::int64_t> depths) {
  if (depths.size() <= 1) {
    return {};
  }
  std::vector<std::uint64_t> bits((depths.size() - 1 + 63) / 64, 0);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    if (depths[i] - depths[i - 1] == 1) {
      bits[(i - 1) >> 6] |= std::uint64_t{1} << ((i - 1) & 63);
    }
  }
  return bits;
}

bool packed_bit(std::span<const std::uint64_t> words, std::size_t position) {
  return ((words[position >> 6] >> (position & 63)) & 1u) != 0;
}

template <class Rmq>
void expect_valid_bp_encoding(const Rmq& rmq, std::size_t value_count) {
  const std::size_t bit_count = rmq.bp_bit_count();
  const std::span<const std::uint64_t> words = rmq.bp_words();
  EXPECT_EQ(bit_count, 2 * value_count);
  EXPECT_EQ(words.size(), (bit_count + 63) / 64);

  std::size_t ones = 0;
  std::size_t zeros = 0;
  std::int64_t excess = 0;
  for (std::size_t position = 0; position < bit_count; ++position) {
    if (packed_bit(words, position)) {
      ++ones;
      ++excess;
    } else {
      ++zeros;
      --excess;
    }
    EXPECT_GE(excess, 0) << "position=" << position;
  }

  EXPECT_EQ(ones, value_count);
  EXPECT_EQ(zeros, value_count);
  EXPECT_EQ(excess, 0);
}

struct SparseTableCase {
  using Rmq = pixie::rmq::SparseTable<int>;
  using MaxRmq = pixie::rmq::SparseTable<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct SegmentTreeCase {
  using Rmq = pixie::rmq::SegmentTree<int>;
  using MaxRmq = pixie::rmq::SegmentTree<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct CartesianTreeCase {
  using Rmq = pixie::rmq::CartesianTreeRmq<int>;
  using MaxRmq = pixie::rmq::CartesianTreeRmq<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct NodeEulerBTreeCase {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  using MaxRmq = pixie::rmq::NodeEulerBTreeRmq<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct ExperimentalNodeEulerBTreeCase {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  using MaxRmq =
      pixie::rmq::experimental::NodeEulerBTreeRmq<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct ExperimentalNodeEulerBTreeMaskLeafCase {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<
      int,
      std::less<int>,
      std::size_t,
      248,
      256,
      pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>;
  using MaxRmq = pixie::rmq::experimental::NodeEulerBTreeRmq<
      int,
      std::greater<int>,
      std::size_t,
      248,
      256,
      pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct SegmentBTreeXlCase {
  using Rmq = pixie::rmq::SegmentBTreeXl<int>;
  using MaxRmq = pixie::rmq::SegmentBTreeXl<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct BpPlusMinusOne128Case {
  using Rmq = pixie::rmq::BpPlusMinusOneRmq<>;
  static constexpr std::size_t kBlockSize = 128;

  static Rmq make(std::span<const std::uint64_t> bits,
                  std::size_t depth_count) {
    return Rmq(bits, depth_count);
  }
};

struct BpPlusMinusOne64Case {
  using Rmq = pixie::rmq::BpPlusMinusOneRmq<std::size_t, 64>;
  static constexpr std::size_t kBlockSize = 64;

  static Rmq make(std::span<const std::uint64_t> bits,
                  std::size_t depth_count) {
    return Rmq(bits, depth_count);
  }
};

struct OneIntervalBTreeCase {
  using Rmq = pixie::rmq::OneIntervalBTreeRmq<>;
  static constexpr std::size_t kBlockSize = Rmq::kBlockSize;

  static Rmq make(std::span<const std::uint64_t> bits,
                  std::size_t depth_count) {
    return Rmq(bits, depth_count);
  }
};

struct OneIntervalBTreeBoundaryRecordsCase {
  using Rmq = pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<>;
  static constexpr std::size_t kBlockSize = Rmq::kBlockSize;

  static Rmq make(std::span<const std::uint64_t> bits,
                  std::size_t depth_count) {
    return Rmq(bits, depth_count);
  }
};

template <class DepthCase>
typename DepthCase::Rmq make_depth_rmq(const std::vector<std::uint64_t>& bits,
                                       std::size_t depth_count) {
  return DepthCase::make(std::span<const std::uint64_t>(bits), depth_count);
}

template <class DepthCase>
void check_depths_all_arg_min(std::span<const std::int64_t> depths) {
  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const typename DepthCase::Rmq rmq =
      make_depth_rmq<DepthCase>(bits, depths.size());
  check_all_arg_min_ranges(rmq, depths, std::less<std::int64_t>());
}

template <class Rmq>
void check_depth_ranges(
    const Rmq& rmq,
    std::span<const std::int64_t> depths,
    std::span<const std::pair<std::size_t, std::size_t>> ranges) {
  for (const auto [left, right] : ranges) {
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(depths, left, right, std::less<std::int64_t>()))
        << "range=[" << left << "," << right << ")";
  }
}

template <class DepthCase>
void check_depth_case_ranges(
    std::span<const std::int64_t> depths,
    std::span<const std::pair<std::size_t, std::size_t>> ranges) {
  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const typename DepthCase::Rmq rmq =
      make_depth_rmq<DepthCase>(bits, depths.size());
  check_depth_ranges(rmq, depths, ranges);
}

void check_one_interval_btree_ranges(
    std::span<const std::int64_t> depths,
    std::span<const std::pair<std::size_t, std::size_t>> ranges) {
  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const pixie::rmq::OneIntervalBTreeRmq<> rmq(bits, depths.size());
  const pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<> records(bits,
                                                                 depths.size());
  const pixie::rmq::BpPlusMinusOneRmq<> baseline(bits, depths.size());
  for (const auto [left, right] : ranges) {
    const std::size_t expected =
        naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                      std::less<std::int64_t>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(rmq.arg_min(left, right), baseline.arg_min(left, right))
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(records.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(records.arg_min(left, right), baseline.arg_min(left, right))
        << "range=[" << left << "," << right << ")";
  }
}

}  // namespace

template <class Case>
class ValueRmqContractTest : public ::testing::Test {};

using ValueRmqCases = ::testing::Types<SparseTableCase,
                                       SegmentTreeCase,
                                       CartesianTreeCase,
                                       NodeEulerBTreeCase,
                                       ExperimentalNodeEulerBTreeCase,
                                       ExperimentalNodeEulerBTreeMaskLeafCase,
                                       SegmentBTreeXlCase>;
TYPED_TEST_SUITE(ValueRmqContractTest, ValueRmqCases);

TYPED_TEST(ValueRmqContractTest, ExhaustiveSmallArray) {
  const std::vector<int> values = {4, 1, 3, 1, 5, 0, 0, 2};
  const typename TypeParam::Rmq rmq =
      TypeParam::make(std::span<const int>(values));
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

TYPED_TEST(ValueRmqContractTest, FirstMinimumTieBreaking) {
  const std::vector<int> values = {7, 2, 2, 3, 2};
  const typename TypeParam::Rmq rmq =
      TypeParam::make(std::span<const int>(values));

  EXPECT_EQ(rmq.arg_min(0, 5), 1u);
  EXPECT_EQ(rmq.arg_min(2, 5), 2u);
}

TYPED_TEST(ValueRmqContractTest, InvalidAndEmptyRanges) {
  using Rmq = typename TypeParam::Rmq;

  const std::vector<int> values = {3, 1, 2};
  const Rmq rmq = TypeParam::make(std::span<const int>(values));

  EXPECT_EQ(rmq.arg_min(2, 2), Rmq::npos);
  EXPECT_EQ(rmq.arg_min(2, 1), Rmq::npos);
  EXPECT_EQ(rmq.arg_min(0, values.size() + 1), Rmq::npos);
  EXPECT_EQ(rmq.range_min(2, 2), 0);
  EXPECT_EQ(rmq.arg_min(0, values.size()), 1u);

  const std::vector<int> empty_values;
  const Rmq default_rmq;
  const Rmq empty_rmq = TypeParam::make(std::span<const int>(empty_values));
  EXPECT_TRUE(default_rmq.empty());
  EXPECT_TRUE(empty_rmq.empty());
  EXPECT_EQ(default_rmq.arg_min(0, 0), Rmq::npos);
  EXPECT_EQ(empty_rmq.arg_min(0, 0), Rmq::npos);
}

TYPED_TEST(ValueRmqContractTest, ComparatorCanSelectMaximum) {
  const std::vector<int> values = {1, 8, 3, 8, 4};
  const typename TypeParam::MaxRmq rmq =
      TypeParam::make_max(std::span<const int>(values));

  check_all_ranges(rmq, std::span<const int>(values), std::greater<int>());
  EXPECT_EQ(rmq.arg_min(0, 5), 1u);
}

TYPED_TEST(ValueRmqContractTest, MonotoneArrays) {
  const std::vector<int> increasing = {1, 2, 3, 4, 5, 6};
  const std::vector<int> decreasing = {6, 5, 4, 3, 2, 1};
  const typename TypeParam::Rmq increasing_rmq =
      TypeParam::make(std::span<const int>(increasing));
  const typename TypeParam::Rmq decreasing_rmq =
      TypeParam::make(std::span<const int>(decreasing));

  check_all_ranges(increasing_rmq, std::span<const int>(increasing),
                   std::less<int>());
  check_all_ranges(decreasing_rmq, std::span<const int>(decreasing),
                   std::less<int>());
}

TYPED_TEST(ValueRmqContractTest, DifferentialRandom) {
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> value_dist(-50, 50);
  for (std::size_t size = 1; size <= 257; size += 17) {
    std::vector<int> values(size);
    std::generate(values.begin(), values.end(),
                  [&] { return value_dist(rng); });

    const typename TypeParam::Rmq rmq =
        TypeParam::make(std::span<const int>(values));
    check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
  }
}

TEST(RmqNodeEulerBTree, BoundarySizesAroundLeavesAndInternalNodes) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  const std::vector<std::size_t> sizes = {
      kLeaf - 1,
      kLeaf,
      kLeaf + 1,
      kLeaf * kFanout - 1,
      kLeaf * kFanout,
      kLeaf * kFanout + 1,
      kLeaf * kFanout + kLeaf + 1,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 17 + i / 5) % 31);
      if (i % 9 == 0 || i % 32 == 3) {
        values[i] = 4;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    const auto check_range = [&](std::size_t left, std::size_t right) {
      ASSERT_LT(left, right);
      ASSERT_LE(right, values.size());
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
      EXPECT_EQ(rmq.range_min(left, right), values[expected])
          << "range=[" << left << "," << right << ")";
    };

    check_range(0, 1);
    check_range(0, size);
    check_range(size - 1, size);

    const std::vector<std::size_t> boundaries = {
        kLeaf,
        kLeaf * kFanout,
    };
    for (const std::size_t boundary : boundaries) {
      if (boundary >= size) {
        continue;
      }
      const std::size_t near_left = boundary > 19 ? boundary - 19 : 0;
      const std::size_t wide_left =
          boundary > kLeaf + 3 ? boundary - kLeaf - 3 : 0;
      check_range(boundary - 1, boundary);
      check_range(boundary, boundary + 1);
      check_range(boundary - 1, boundary + 1);
      check_range(near_left, std::min(size, boundary + 23));
      check_range(wide_left, std::min(size, boundary + kLeaf + 5));
      check_range(1, std::min(size, boundary + 1));
      check_range(boundary, size);
    }
  }
}

TEST(RmqNodeEulerBTree, CommonAncestorQueryShapes) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  const std::size_t first_level_span = kLeaf * kFanout;
  std::vector<int> values(first_level_span + 3 * kLeaf + 17);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 29 + i / 13) % 101);
    if (i % 257 == 0) {
      values[i] = -7;
    }
    if (i % 4099 == 3) {
      values[i] = -9;
    }
  }

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {7, 99},
      {kLeaf - 3, 2 * kLeaf + 5},
      {17, 18 * kLeaf + 111},
      {first_level_span - 4, first_level_span + kLeaf + 9},
      {kLeaf * (kFanout - 7) + 13, first_level_span + 2 * kLeaf + 19},
      {0, first_level_span},
      {kLeaf, first_level_span + kLeaf},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqNodeEulerBTree, SelectorFirstCorrectsInvalidBorderMinima) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  std::vector<int> values(3 * kLeaf, 1000);
  values[0] = -1000;
  values[kLeaf / 2] = 50;
  values[kLeaf + 44] = -10;
  values[2 * kLeaf + 10] = 60;
  values[2 * kLeaf + 200] = -900;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kLeaf / 2, 2 * kLeaf + 100},
      {kLeaf / 2, kLeaf + 100},
      {kLeaf - 50, 3 * kLeaf - 100},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqNodeEulerBTree, DuplicateHeavyRandomDifferentialTo8193) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  std::mt19937_64 rng(811);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 255, 256, 257, 1024, 8191, 8192, 8193, kLeaf * kFanout + 1,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 11) < 5) {
        values[i] = 0;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqNodeEulerBTree, TopDepthSparseSelectorLevelShapes) {
  using MinRmq =
      pixie::rmq::NodeEulerBTreeRmq<int, std::less<int>, std::size_t, 4, 4>;
  using MaxRmq =
      pixie::rmq::NodeEulerBTreeRmq<int, std::greater<int>, std::size_t, 4, 4>;
  const std::vector<std::size_t> sizes = {
      4,    // no internal level
      5,    // one internal level
      17,   // two internal levels
      65,   // three internal levels
      257,  // more than three internal levels
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 19 + i / 3) % 23);
      if (i % 7 == 0 || i % 11 == 3) {
        values[i] = -2;
      }
      if (i % 29 == 5) {
        values[i] = 31;
      }
    }

    const MinRmq min_rmq{std::span<const int>(values)};
    const MaxRmq max_rmq{std::span<const int>(values)};
    check_all_arg_min_ranges(min_rmq, std::span<const int>(values),
                             std::less<int>());
    check_all_arg_min_ranges(max_rmq, std::span<const int>(values),
                             std::greater<int>());
  }
}

TEST(RmqNodeEulerBTree, TopDepthSparseSelectorDefaultFanoutRandom) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  using MaxRmq = pixie::rmq::NodeEulerBTreeRmq<int, std::greater<int>>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  const std::size_t size = kLeaf * kFanout + kLeaf + 37;

  std::vector<int> values(size);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 31 + i / 9) % 127);
    if (i % 257 == 3 || i % 8191 == 17) {
      values[i] = -9;
    }
    if (i % 4099 == 19) {
      values[i] = 211;
    }
  }

  const Rmq rmq{std::span<const int>(values)};
  const MaxRmq max_rmq{std::span<const int>(values)};
  std::mt19937_64 rng(424242);
  std::uniform_int_distribution<std::size_t> width_dist(1, values.size());
  for (std::size_t query = 0; query < 3000; ++query) {
    const std::size_t width = width_dist(rng);
    std::uniform_int_distribution<std::size_t> left_dist(0,
                                                         values.size() - width);
    const std::size_t left = left_dist(rng);
    const std::size_t right = left + width;
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::less<int>()))
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(max_rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::greater<int>()))
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqNodeEulerBTree, CopyAndMovePreserveSelectors) {
  using Rmq = pixie::rmq::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  std::vector<int> values(kLeaf * kFanout + 512);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 43 + i / 7) % 97);
    if (i % 13 == 0) {
      values[i] = -5;
    }
  }

  const Rmq original{std::span<const int>(values)};
  Rmq copied(original);
  Rmq assigned;
  assigned = original;
  Rmq moved(std::move(copied));
  Rmq move_assigned;
  move_assigned = std::move(assigned);

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 1},       {0, values.size()},
      {255, 257},   {kLeaf * kFanout - 5, kLeaf * kFanout + 6},
      {1234, 5678}, {values.size() - 17, values.size()},
  };

  const std::array<const Rmq*, 3> rmqs = {&original, &moved, &move_assigned};
  for (const Rmq* rmq : rmqs) {
    ASSERT_EQ(rmq->size(), values.size());
    for (const auto [left, right] : ranges) {
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq->arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqExperimentalNodeEulerBTree, SplitLevelBoundarySizes) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  constexpr std::size_t kMediumSpan = kLeaf * kFanout;
  const std::vector<std::size_t> sizes = {
      kLeaf - 1,
      kLeaf,
      kLeaf + 1,
      kMediumSpan - 1,
      kMediumSpan,
      kMediumSpan + 1,
      kMediumSpan + 3 * kLeaf + 17,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 37 + i / 11) % 113);
      if (i % 17 == 0 || i % 257 == 3) {
        values[i] = -4;
      }
      if (i % 4099 == 19) {
        values[i] = -11;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    const auto check_range = [&](std::size_t left, std::size_t right) {
      ASSERT_LT(left, right);
      ASSERT_LE(right, values.size());
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
      EXPECT_EQ(rmq.range_min(left, right), values[expected])
          << "range=[" << left << "," << right << ")";
    };

    check_range(0, 1);
    check_range(0, size);
    check_range(size - 1, size);

    for (const std::size_t boundary : {kLeaf, kMediumSpan}) {
      if (boundary == 0 || boundary >= size) {
        continue;
      }
      check_range(boundary - 1, boundary);
      check_range(boundary, boundary + 1);
      check_range(boundary - 1, boundary + 1);
      check_range(boundary > 67 ? boundary - 67 : 0,
                  std::min(size, boundary + 71));
      check_range(boundary > kLeaf + 5 ? boundary - kLeaf - 5 : 0,
                  std::min(size, boundary + kLeaf + 9));
    }
  }
}

TEST(RmqExperimentalNodeEulerBTree, LeafEmbeddedOffsetBoundaryRanges) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;

  std::vector<int> values(2 * kLeaf + 5);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 29 + i / 7) % 97) + 10;
  }
  values[kLeaf - 1] = -7;
  values[kLeaf] = -7;
  values[2 * kLeaf] = -9;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kLeaf},         {0, kLeaf + 1},
      {kLeaf, kLeaf + 1}, {kLeaf - 3, kLeaf + 4},
      {kLeaf, kLeaf + 4}, {2 * kLeaf - 2, 2 * kLeaf + 3},
      {0, values.size()},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(rmq.range_min(left, right), values[expected])
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqExperimentalNodeEulerBTree, MaskLeafBoundaryAndFallbackRanges) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<
      int, std::less<int>, std::size_t, 248, 256,
      pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  static_assert(kLeaf == 248);

  std::vector<int> values(2 * kLeaf + 13, 1000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>(i % 17);
  }

  values[20] = 80;
  values[40] = 70;
  values[80] = 60;
  values[90] = 60;
  values[123] = -100;
  values[160] = 10;
  values[200] = 5;
  values[220] = 0;
  values[kLeaf - 1] = 0;

  values[kLeaf + 7] = -7;
  values[kLeaf + 91] = -11;
  values[2 * kLeaf + 3] = -13;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kLeaf},
      {0, 41},
      {0, 100},
      {0, 123},
      {90, 160},
      {130, kLeaf},
      {221, kLeaf},
      {130, 210},
      {kLeaf - 5, kLeaf + 8},
      {kLeaf, 2 * kLeaf},
      {2 * kLeaf - 3, values.size()},
      {0, values.size()},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(rmq.range_min(left, right), values[expected])
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqSegmentBTreeXl, BoundaryAndFallbackRanges) {
  using Rmq = pixie::rmq::SegmentBTreeXl<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  static_assert(kLeaf == 496);

  std::vector<int> values(2 * kLeaf + 17, 2000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>((i * 7 + i / 3) % 31);
  }

  values[32] = 120;
  values[96] = 70;
  values[160] = 55;
  values[240] = -100;
  values[241] = -100;
  values[320] = 30;
  values[440] = 5;
  values[kLeaf - 1] = 5;

  values[kLeaf + 11] = -9;
  values[kLeaf + 173] = -17;
  values[2 * kLeaf + 5] = -21;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kLeaf},         {0, 97},
      {0, 240},           {96, 241},
      {250, kLeaf},       {441, kLeaf},
      {300, 450},         {kLeaf - 7, kLeaf + 12},
      {kLeaf, 2 * kLeaf}, {2 * kLeaf - 5, values.size()},
      {0, values.size()},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(rmq.range_min(left, right), values[expected])
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqSegmentBTreeXl, LeafSelectorEnumVariants) {
  using MaskRmq = pixie::rmq::SegmentBTreeXL<
      int, std::less<int>, std::size_t, 248, 256,
      pixie::rmq::SegmentBTreeXLLeafSelector::PrefixSuffix>;
  using BpRmq =
      pixie::rmq::SegmentBTreeXL<int, std::less<int>, std::size_t, 252, 256,
                                 pixie::rmq::SegmentBTreeXLLeafSelector::BP>;

  std::vector<int> values(4099, 10000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 37 + i / 7) % 257);
  }
  values[13] = -50;
  values[251] = -70;
  values[252] = -70;
  values[747] = -90;
  values[2048] = -120;
  values[4098] = -110;

  const MaskRmq mask_rmq{std::span<const int>(values)};
  const BpRmq bp_rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 1},      {0, 252},     {1, 251},           {248, 253},   {251, 753},
      {700, 2100}, {2048, 2049}, {0, values.size()}, {3000, 4099},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(mask_rmq.arg_min(left, right), expected)
        << "mask range=[" << left << "," << right << ")";
    EXPECT_EQ(bp_rmq.arg_min(left, right), expected)
        << "bp range=[" << left << "," << right << ")";
  }
}

TEST(RmqExperimentalNodeEulerBTree, DuplicateHeavyRandomDifferentialTo8193) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  std::mt19937_64 rng(9127);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 255, 256, 257, 1024, 4096, 8191, 8192, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 13) < 7) {
        values[i] = 0;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqExperimentalNodeEulerBTree, MaskLeafDuplicateHeavyRandomTo8193) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<
      int, std::less<int>, std::size_t, 248, 256,
      pixie::rmq::experimental::PrefixSuffixMaskLeafSelectorTag>;
  std::mt19937_64 rng(901248);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 247, 248, 249, 1024, 4096, 8191, 8192, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 11) < 6) {
        values[i] = 0;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqSegmentBTreeXl, DuplicateHeavyRandomTo8193) {
  using Rmq = pixie::rmq::SegmentBTreeXl<int>;
  std::mt19937_64 rng(901496);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 495, 496, 497, 1024, 4096, 8191, 8192, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 11) < 6) {
        values[i] = 0;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqExperimentalNodeEulerBTree, HighSparseSelectorWithMiddleLevels) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  constexpr std::size_t kMiddleBoundary = kLeaf * kFanout * kFanout;
  const std::size_t size = kMiddleBoundary + 2 * kLeaf + 17;

  std::vector<int> values(size);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 53 + i / 17) % 251);
    if (i % 257 == 3 || i % 65537 == 11) {
      values[i] = -5;
    }
  }
  values[13] = -100;
  values[kMiddleBoundary - 29] = -90;
  values[kMiddleBoundary + kLeaf + 7] = -110;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, values.size()},
      {kMiddleBoundary - 2 * kLeaf, kMiddleBoundary + 2 * kLeaf},
      {kMiddleBoundary - 31, kMiddleBoundary + kLeaf + 31},
      {kMiddleBoundary + 1, values.size()},
      {kLeaf * (kFanout - 1) + 9, kMiddleBoundary + kLeaf + 13},
      {values.size() - 3 * kLeaf, values.size()},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqExperimentalNodeEulerBTree, CopyAndMovePreserveSplitStorage) {
  using Rmq = pixie::rmq::experimental::NodeEulerBTreeRmq<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kFanout = Rmq::kFanout;
  std::vector<int> values(kLeaf * kFanout + 2 * kLeaf + 9);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 41 + i / 3) % 127);
    if (i % 23 == 0) {
      values[i] = -8;
    }
  }

  const Rmq original{std::span<const int>(values)};
  Rmq copied(original);
  Rmq assigned;
  assigned = original;
  Rmq moved(std::move(copied));
  Rmq move_assigned;
  move_assigned = std::move(assigned);

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 1},
      {0, values.size()},
      {255, 257},
      {kLeaf * kFanout - 5, kLeaf * kFanout + 6},
      {values.size() - 513, values.size() - 3},
      {values.size() - 17, values.size()},
  };

  const std::array<const Rmq*, 3> rmqs = {&original, &moved, &move_assigned};
  for (const Rmq* rmq : rmqs) {
    ASSERT_EQ(rmq->size(), values.size());
    for (const auto [left, right] : ranges) {
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq->arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

template <class Case>
class DepthRmqContractTest : public ::testing::Test {};

using DepthRmqCases = ::testing::Types<BpPlusMinusOne128Case,
                                       BpPlusMinusOne64Case,
                                       OneIntervalBTreeCase,
                                       OneIntervalBTreeBoundaryRecordsCase>;
TYPED_TEST_SUITE(DepthRmqContractTest, DepthRmqCases);

TYPED_TEST(DepthRmqContractTest, EmptyAndSingleDepthInputs) {
  using Rmq = typename TypeParam::Rmq;

  const std::vector<std::uint64_t> bits;
  const Rmq empty = TypeParam::make(std::span<const std::uint64_t>(bits), 0);
  EXPECT_TRUE(empty.empty());
  EXPECT_EQ(empty.size(), 0u);
  EXPECT_EQ(empty.arg_min(0, 0), Rmq::npos);

  const Rmq single = TypeParam::make(std::span<const std::uint64_t>(bits), 1);
  EXPECT_FALSE(single.empty());
  EXPECT_EQ(single.size(), 1u);
  EXPECT_EQ(single.arg_min(0, 1), 0u);
  EXPECT_EQ(single.arg_min(1, 1), Rmq::npos);
  EXPECT_EQ(single.arg_min(0, 2), Rmq::npos);
}

TYPED_TEST(DepthRmqContractTest, ExhaustiveSmallDepthArray) {
  const std::vector<std::int64_t> depths = {0, 1, 0, 1, 2, 1, 2, 1, 0};
  check_depths_all_arg_min<TypeParam>(std::span<const std::int64_t>(depths));
}

TYPED_TEST(DepthRmqContractTest, CrossBlockRanges) {
  constexpr std::size_t kBlock = TypeParam::kBlockSize;
  std::vector<std::int64_t> depths(3 * kBlock + 1);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 7 == 0) || (i % 7 == 1) || (i % 7 == 4);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  check_depths_all_arg_min<TypeParam>(std::span<const std::int64_t>(depths));
}

TYPED_TEST(DepthRmqContractTest, BoundaryRangesAroundBlocks) {
  constexpr std::size_t kBlock = TypeParam::kBlockSize;
  std::vector<std::int64_t> depths(3 * kBlock + 6);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] =
        depths[i - 1] + ((i % 5 == 0 || i % 11 == 0 || i % 17 == 0) ? 1 : -1);
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 1},
      {kBlock - 2, kBlock},
      {kBlock - 1, kBlock + 1},
      {kBlock, kBlock + 1},
      {kBlock, 2 * kBlock},
      {kBlock + 1, 2 * kBlock},
      {2 * kBlock - 1, 2 * kBlock + 1},
      {0, depths.size()},
      {kBlock - 31, 2 * kBlock + 3},
      {kBlock - 12, 2 * kBlock + 6},
      {2 * kBlock, depths.size()},
  };
  check_depth_case_ranges<TypeParam>(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TYPED_TEST(DepthRmqContractTest, LongSequenceRangesNearEnd) {
  std::vector<std::int64_t> depths(8193);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 13 == 0) || (i % 17 == 0) || (i % 19 == 0);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {depths.size() - 1, depths.size()},
      {depths.size() - 128, depths.size()},
      {depths.size() - 513, depths.size() - 2},
      {depths.size() - 4097, depths.size() - 6},
  };
  check_depth_case_ranges<TypeParam>(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TYPED_TEST(DepthRmqContractTest, CrossBlockTieKeepsFirstPosition) {
  constexpr std::size_t kBlock = TypeParam::kBlockSize;
  std::vector<std::int64_t> depths(3 * kBlock);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] = (i % 2 == 0) ? 0 : 1;
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const typename TypeParam::Rmq rmq =
      make_depth_rmq<TypeParam>(bits, depths.size());
  const std::size_t left = kBlock - 8;
  const std::size_t right = 2 * kBlock + 5;

  EXPECT_EQ(rmq.arg_min(left, right), left);
  EXPECT_EQ(rmq.arg_min(left, right),
            naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                          std::less<std::int64_t>()));
}

TYPED_TEST(DepthRmqContractTest, DisjointBoundaryRangeMatchesNaive) {
  constexpr std::size_t kBlock = TypeParam::kBlockSize;
  std::vector<std::int64_t> depths(3 * kBlock + 7);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 9 == 0) || (i % 9 == 3) || (i % 11 == 0);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kBlock - 32, kBlock + 33},       {kBlock - 8, kBlock + 13},
      {kBlock - 1, kBlock + 2},         {2 * kBlock - 10, 2 * kBlock + 5},
      {2 * kBlock - 1, 2 * kBlock + 4},
  };
  check_depth_case_ranges<TypeParam>(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TYPED_TEST(DepthRmqContractTest, RejectsTooSmallBitSpan) {
  const std::vector<std::uint64_t> bits;
  EXPECT_THROW(
      {
        const typename TypeParam::Rmq rmq =
            TypeParam::make(std::span<const std::uint64_t>(bits), 2);
        (void)rmq;
      },
      std::invalid_argument);
}

TYPED_TEST(DepthRmqContractTest, DifferentialRandomWalks) {
  std::mt19937_64 rng(77);
  for (std::size_t size = 1; size <= 257; size += 17) {
    std::vector<std::int64_t> depths(size);
    for (std::size_t i = 1; i < depths.size(); ++i) {
      depths[i] = depths[i - 1] + ((rng() & 1u) != 0 ? 1 : -1);
    }

    check_depths_all_arg_min<TypeParam>(std::span<const std::int64_t>(depths));
  }
}

TEST(RmqOneIntervalBTree, FusedBoundaryTiesKeepFirstPosition) {
  constexpr std::size_t kBlock = pixie::rmq::OneIntervalBTreeRmq<>::kBlockSize;
  std::vector<std::int64_t> depths(4 * kBlock + 9);
  for (std::size_t i = 0; i < depths.size(); ++i) {
    depths[i] = (i & 1u) == 0 ? 0 : 1;
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kBlock - 8, 2 * kBlock + 8},
      {kBlock - 1, 3 * kBlock + 4},
      {kBlock + 3, 4 * kBlock + 2},
      {2 * kBlock - 5, 4 * kBlock + 6},
  };
  check_one_interval_btree_ranges(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TEST(RmqOneIntervalBTreeBoundaryRecords, CrossBlockAndTailRanges) {
  using Rmq = pixie::rmq::OneIntervalBTreeRmqBoundaryRecords<>;
  constexpr std::size_t kBlock = Rmq::kBlockSize;
  std::vector<std::int64_t> depths(2 * kBlock + 17);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    const bool up = (i % 6 == 0) || (i % 10 == 1) || (i % 17 == 3);
    depths[i] = depths[i - 1] + (up ? 1 : -1);
  }

  const std::vector<std::uint64_t> bits = pack_depth_deltas(depths);
  const Rmq records(bits, depths.size());
  const pixie::rmq::OneIntervalBTreeRmq<> baseline(bits, depths.size());
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kBlock - 6, kBlock + 1},     {kBlock - 3, 2 * kBlock + 2},
      {kBlock - 7, 2 * kBlock},     {kBlock, 2 * kBlock + 5},
      {kBlock - 5, 2 * kBlock + 6}, {2 * kBlock - 3, depths.size()},
      {0, depths.size()},
  };

  for (const auto [left, right] : ranges) {
    const std::size_t expected =
        naive_arg_min(std::span<const std::int64_t>(depths), left, right,
                      std::less<std::int64_t>());
    EXPECT_EQ(records.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(records.arg_min(left, right), baseline.arg_min(left, right))
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqOneIntervalBTree, LowNodeBoundaryRanges) {
  constexpr std::size_t kBlock = pixie::rmq::OneIntervalBTreeRmq<>::kBlockSize;
  constexpr std::size_t kLowBoundary =
      kBlock * pixie::rmq::OneIntervalBTreeRmq<>::kLowFanout;
  std::vector<std::int64_t> depths(kLowBoundary + 4 * kBlock + 1);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] =
        depths[i - 1] + ((i % 7 == 0 || i % 19 == 0 || i % 23 == 0) ? 1 : -1);
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kLowBoundary - 2, kLowBoundary + 2},
      {kLowBoundary - kBlock - 1, kLowBoundary + kBlock + 1},
      {1, kLowBoundary + 1},
      {kBlock, kLowBoundary + kBlock},
      {kLowBoundary, depths.size()},
      {0, depths.size()},
  };
  check_one_interval_btree_ranges(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TEST(RmqOneIntervalBTree, HighNodeBoundaryRanges) {
  constexpr std::size_t kBlock = pixie::rmq::OneIntervalBTreeRmq<>::kBlockSize;
  constexpr std::size_t kHighBoundary =
      kBlock * pixie::rmq::OneIntervalBTreeRmq<>::kLowFanout *
      pixie::rmq::OneIntervalBTreeRmq<>::kHighFanout;
  std::vector<std::int64_t> depths(kHighBoundary + 32 * kBlock + 1);
  for (std::size_t i = 1; i < depths.size(); ++i) {
    depths[i] =
        depths[i - 1] + ((i % 13 == 0 || i % 29 == 0 || i % 31 == 0) ? 1 : -1);
  }

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {kHighBoundary - 2, kHighBoundary + 2},
      {kHighBoundary - 16 * kBlock - 1, kHighBoundary + 16 * kBlock + 1},
      {kBlock, kHighBoundary + kBlock},
      {kHighBoundary, depths.size()},
      {0, depths.size()},
  };
  check_one_interval_btree_ranges(
      std::span<const std::int64_t>(depths),
      std::span<const std::pair<std::size_t, std::size_t>>(ranges));
}

TEST(RmqCartesianTree, DuplicateHeavyArrays) {
  const std::vector<std::vector<int>> cases = {
      {5, 5, 5, 5, 5, 5, 5},
      {4, 1, 4, 1, 4, 1, 4, 1},
      {3, 2, 2, 2, 1, 1, 2, 1, 3},
      {9, 8, 9, 8, 9, 7, 7, 7, 8, 9},
  };

  for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
    SCOPED_TRACE(case_index);
    const std::vector<int>& values = cases[case_index];
    const pixie::rmq::CartesianTreeRmq<int> cartesian(values);
    expect_valid_bp_encoding(cartesian, values.size());
    check_all_ranges(cartesian, std::span<const int>(values), std::less<int>());
  }
}

TEST(RmqCartesianTree, ComparatorTiesKeepFirstMaximum) {
  const std::vector<int> values = {1, 8, 8, 7, 8, 3, 8};
  const pixie::rmq::CartesianTreeRmq<int, std::greater<int>> cartesian(values);

  expect_valid_bp_encoding(cartesian, values.size());
  check_all_ranges(cartesian, std::span<const int>(values),
                   std::greater<int>());
  EXPECT_EQ(cartesian.arg_min(0, values.size()), 1u);
  EXPECT_EQ(cartesian.arg_min(2, values.size()), 2u);
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

TEST(RmqCartesianTree, BpEncodingIsBalanced) {
  const std::vector<int> values = {4, 1, 3, 2, 5};
  const pixie::rmq::CartesianTreeRmq<int> rmq(values);
  expect_valid_bp_encoding(rmq, values.size());
}

TEST(RmqCartesianTree, BoundarySizesExerciseBpBlocks) {
  for (std::size_t size : {1u, 63u, 64u, 65u, 127u, 128u, 129u, 255u, 256u,
                           257u, 511u, 512u, 513u}) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 37 + i / 3) % 23);
      if (i % 11 == 0) {
        values[i] = 7;
      }
    }

    const pixie::rmq::CartesianTreeRmq<int> rmq(values);
    expect_valid_bp_encoding(rmq, values.size());

    const auto check_range = [&](std::size_t left, std::size_t right) {
      EXPECT_EQ(rmq.arg_min(left, right),
                naive_arg_min(std::span<const int>(values), left, right,
                              std::less<int>()))
          << "range=[" << left << "," << right << ")";
    };

    check_range(0, 1);
    check_range(0, size);
    check_range(size - 1, size);
    check_range(size / 3, std::min(size, size / 3 + 9));
    check_range(size / 2, std::min(size, size / 2 + 17));

    for (std::size_t boundary : {64u, 128u, 256u, 512u}) {
      if (boundary < size) {
        check_range(boundary - 1, std::min(size, boundary + 2));
        check_range(boundary / 2, std::min(size, boundary + 1));
      }
    }
  }
}
