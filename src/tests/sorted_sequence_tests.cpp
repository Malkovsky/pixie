#include <gtest/gtest.h>
#include <pixie/permutations/sorted_sequence.h>
#include <pixie/permutations/sorted_vector.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using namespace pixie;

template <class Sorted>
void check(const Sorted& sorted, const std::vector<std::uint64_t>& expected) {
  ASSERT_EQ(sorted.size(), expected.size());
  EXPECT_EQ(sorted.empty(), expected.empty());
  EXPECT_GE(sorted.memory_usage_bytes(), sizeof(sorted));
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ(sorted[i], expected[i]) << "index=" << i;
  }
  for (std::uint64_t probe :
       {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{17},
        std::uint64_t{42}, std::uint64_t{999},
        std::numeric_limits<std::uint64_t>::max()}) {
    const auto lower =
        std::lower_bound(expected.begin(), expected.end(), probe);
    const auto upper =
        std::upper_bound(expected.begin(), expected.end(), probe);
    EXPECT_EQ(sorted.lower_bound_index(probe),
              static_cast<std::size_t>(lower - expected.begin()));
    EXPECT_EQ(sorted.upper_bound_index(probe),
              static_cast<std::size_t>(upper - expected.begin()));
  }
}

template <class Sorted>
void random_insert_differential() {
  Sorted sorted;
  std::vector<std::uint64_t> expected;
  std::mt19937_64 random(713);
  for (std::size_t i = 0; i < 1000; ++i) {
    const std::uint64_t value = random() % 113;
    const auto expected_position = static_cast<std::size_t>(
        std::lower_bound(expected.begin(), expected.end(), value) -
        expected.begin());
    EXPECT_EQ(sorted.insert(value), expected_position);
    expected.insert(expected.begin() + expected_position, value);
    if (i % 31 == 0) {
      check(sorted, expected);
    }
  }
  check(sorted, expected);
}

using PackedSorted = SortedPermutableSequence<std::uint64_t,
                                              std::less<>,
                                              ElementStorage::packed,
                                              512,
                                              4,
                                              LengthLayout::individual>;
using IndirectSorted = SortedPermutableSequence<std::uint64_t,
                                                std::less<>,
                                                ElementStorage::indirect,
                                                512,
                                                4,
                                                LengthLayout::cumulative,
                                                64>;
using PermutationSorted = SortedPermutationVector<std::uint64_t,
                                                  std::less<>,
                                                  std::uint16_t,
                                                  512,
                                                  4,
                                                  LengthLayout::individual>;

TEST(SortedPermutableSequence, PackedRandomInsertDifferential) {
  random_insert_differential<PackedSorted>();
}

TEST(SortedPermutableSequence, IndirectReferencesSurviveInsertion) {
  IndirectSorted sorted;
  sorted.insert(3);
  const auto* three = &sorted[0];
  sorted.insert(2);
  ASSERT_EQ(sorted.size(), 2);
  EXPECT_EQ(sorted[0], 2);
  EXPECT_EQ(sorted[1], 3);
  EXPECT_EQ(&sorted[1], three);
  random_insert_differential<IndirectSorted>();
}

TEST(SortedPermutationVector, RandomInsertDifferential) {
  random_insert_differential<PermutationSorted>();
}

TEST(SortedPermutationVector, ReserveRetainsPayloadAddresses) {
  PermutationSorted sorted;
  sorted.reserve(8);
  sorted.insert(3);
  const auto* three = &sorted[0];
  sorted.insert(2);
  ASSERT_EQ(sorted.size(), 2);
  EXPECT_EQ(sorted[0], 2);
  EXPECT_EQ(sorted[1], 3);
  EXPECT_EQ(&sorted[1], three);
}

TEST(SortedWrappers, FromRangeConsumesAndSorts) {
  std::vector<std::uint64_t> input{9, 1, 8, 1, 4, 7};
  auto packed = PackedSorted::from_range(input);
  auto permutation = PermutationSorted::from_range(std::move(input));
  const std::vector<std::uint64_t> expected{1, 1, 4, 7, 8, 9};
  check(packed, expected);
  check(permutation, expected);
}

TEST(SortedWrappers, FromRangeAcceptsDistinctSentinel) {
  auto input = std::views::iota(std::uint64_t{0}) | std::views::take(6);
  static_assert(!std::ranges::common_range<decltype(input)>);
  const std::vector<std::uint64_t> expected{0, 1, 2, 3, 4, 5};
  check(PackedSorted::from_range(input), expected);
  check(IndirectSorted::from_range(input), expected);
  check(PermutationSorted::from_range(input), expected);
}

TEST(SortedPermutationVector, FailedIndexInsertionDiscardsAppendedPayload) {
  using Sorted = SortedPermutationVector<std::uint64_t, std::less<>,
                                         std::uint64_t, 512, 4>;
  using Permutation = Sorted::permutation_type;
  using Tree = std::remove_cvref_t<
      decltype(std::declval<const Permutation&>().test_tree())>;
  for (const auto n : {std::size_t{0}, Tree::block_capacity}) {
    SCOPED_TRACE(n);
    Sorted sorted;
    std::vector<std::uint64_t> expected;
    for (std::size_t i = 0; i < n; ++i) {
      sorted.insert(i);
      expected.push_back(i);
    }
    sorted.reserve(n + 1);
    Tree::test_fail_after(0);
    EXPECT_THROW(sorted.insert(999), std::bad_alloc);
    Tree::test_fail_after(-1);
    check(sorted, expected);
    sorted.insert(777);
    expected.push_back(777);
    check(sorted, expected);
  }
}

TEST(SortedPermutationVector, EnforcesPermutationIndexDomainBeforeAppend) {
  using SmallIndex =
      SortedPermutationVector<std::uint64_t, std::less<>, std::uint8_t, 512, 4>;
  SmallIndex sorted;
  for (std::size_t i = 0; i <= std::numeric_limits<std::uint8_t>::max(); ++i) {
    sorted.insert(static_cast<std::uint64_t>(i));
  }
  EXPECT_EQ(
      sorted.size(),
      static_cast<std::size_t>(std::numeric_limits<std::uint8_t>::max()) + 1);
  EXPECT_THROW(sorted.insert(std::numeric_limits<std::uint64_t>::max()),
               std::length_error);
  EXPECT_EQ(
      sorted.size(),
      static_cast<std::size_t>(std::numeric_limits<std::uint8_t>::max()) + 1);
}

struct KeyedValue {
  std::uint64_t key;
  std::uint64_t id;
  bool operator==(const KeyedValue&) const = default;
};

struct ByKey {
  bool operator()(const KeyedValue& left, const KeyedValue& right) const {
    return left.key < right.key;
  }
};

TEST(SortedWrappers, EquivalentValuesInsertBeforeExistingValues) {
  using Tree = SortedPermutableSequence<KeyedValue, ByKey,
                                        ElementStorage::indirect, 512, 4>;
  using Vector =
      SortedPermutationVector<KeyedValue, ByKey, std::uint16_t, 512, 4>;
  auto exercise = []<class Sorted>() {
    Sorted sorted;
    sorted.insert({2, 1});
    sorted.insert({2, 2});
    sorted.insert({1, 3});
    ASSERT_EQ(sorted.size(), 3);
    EXPECT_EQ(sorted[0], (KeyedValue{1, 3}));
    EXPECT_EQ(sorted[1], (KeyedValue{2, 2}));
    EXPECT_EQ(sorted[2], (KeyedValue{2, 1}));
  };
  exercise.template operator()<Tree>();
  exercise.template operator()<Vector>();
}

}  // namespace
