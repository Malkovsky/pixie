#include <gtest/gtest.h>
#include <pixie/detail/sequence/packed_value_block.h>
#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/permutations/permutation.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using namespace pixie;
using namespace pixie::detail::sequence;

template <class T>
[[gnu::noinline]] void rotate_values(std::vector<T>& values,
                                     std::size_t left,
                                     std::size_t right,
                                     std::size_t distance) {
  if (left != right) {
    std::rotate(values.begin() + left,
                values.begin() + left + distance % (right - left),
                values.begin() + right);
  }
}

template <class Container>
auto contents(const Container& container) {
  std::vector<typename Container::value_type> result;
  for (std::size_t i = 0; i < container.size(); ++i) {
    result.push_back(container[i]);
  }
  return result;
}

template <class P>
void check_permutation(const P& permutation,
                       const std::vector<typename P::value_type>& expected) {
  ASSERT_TRUE(permutation.test_tree().test_validate());
  ASSERT_EQ(permutation.size(), expected.size());
  ASSERT_EQ(permutation.empty(), expected.empty());
  EXPECT_EQ(contents(permutation), expected);
  std::vector<bool> seen(permutation.size());
  for (std::size_t i = 0; i < permutation.size(); ++i) {
    const auto value = permutation[i];
    ASSERT_LT(value, permutation.size());
    ASSERT_FALSE(seen[value]);
    seen[value] = true;
  }
  const auto memory = permutation.memory_usage();
  EXPECT_EQ(memory.total_bytes,
            memory.facade_bytes + memory.payload_capacity_bytes +
                memory.block_metadata_bytes + memory.ordering_tree_bytes +
                memory.tag_padding_bytes);
  EXPECT_EQ(memory.total_bytes, permutation.memory_usage_bytes());
  EXPECT_EQ(memory.total_bytes, permutation.test_tree().memory_usage_bytes());
}

template <class T>
class PermutationSpec : public ::testing::Test {};
using Permutations = ::testing::Types<
    Permutation<std::uint16_t, 512, 4, LengthLayout::cumulative>,
    Permutation<std::uint16_t, 512, 4, LengthLayout::individual>,
    Permutation<std::uint32_t, 512, 8, LengthLayout::cumulative>,
    Permutation<std::uint32_t, 512, 8, LengthLayout::individual>,
    Permutation<std::uint64_t, 512, 16, LengthLayout::cumulative>,
    Permutation<std::uint64_t, 512, 16, LengthLayout::individual>>;
TYPED_TEST_SUITE(PermutationSpec, Permutations);

TYPED_TEST(PermutationSpec, PublicContractAliasesNoexceptAndLayout) {
  using P = TypeParam;
  using T = typename P::value_type;
  using Base = PermutationBase<P, T>;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  static_assert(std::is_base_of_v<Base, P>);
  static_assert(std::is_empty_v<Base>);
  static_assert(std::same_as<typename P::size_type, std::size_t>);
  static_assert(std::same_as<typename P::const_reference, T>);
  static_assert(std::same_as<decltype(P::identity(0)), P>);
  static_assert(std::same_as<decltype(std::declval<const Base&>()[0]), T>);
  static_assert(noexcept(std::declval<const Base&>().size()));
  static_assert(noexcept(std::declval<const Base&>().empty()));
  static_assert(noexcept(std::declval<const Base&>().memory_usage_bytes()));
  static_assert(!noexcept(std::declval<const Base&>()[0]));
  static_assert(sizeof(P) == sizeof(Tree));
  static_assert(alignof(P) == alignof(Tree));
  auto owner = Base::identity(3);
  Base& contract = owner;
  contract.rotate_left(0, 3, 1);
  auto donor = Base::identity(2);
  contract.merge(donor);
  EXPECT_EQ(contents(owner), (std::vector<T>{1, 2, 0, 3, 4}));
  EXPECT_THROW(contract[contract.size()], std::out_of_range);
  EXPECT_EQ(contract.memory_usage_bytes(), owner.memory_usage().total_bytes);
}

TYPED_TEST(PermutationSpec, IdentityOwnershipExactMergeAndRanges) {
  using P = TypeParam;
  using T = typename P::value_type;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  static_assert(!std::is_copy_constructible_v<P>);
  static_assert(std::is_nothrow_move_constructible_v<P>);
  static_assert(std::is_nothrow_move_assignable_v<P>);
  P empty;
  check_permutation(empty, {});
  auto zero = P::identity(0);
  check_permutation(zero, {});
  empty.rotate_left(0, 0, -1);
  EXPECT_THROW(empty[0], std::out_of_range);
  EXPECT_THROW(empty.rotate_left(0, 1, 0), std::out_of_range);
  EXPECT_THROW(empty.rotate_left(1, 0, 0), std::out_of_range);
  auto a = P::identity(3);
  auto b = P::identity(2);
  a.rotate_left(1, 3, 1);
  b.rotate_left(0, 2, 1);
  Tree::test_reset_counters();
  Tree::test_fail_after(0);
  a.merge(b);
  Tree::test_fail_after(-1);
  EXPECT_EQ(Tree::test_counters.allocations, 0);
  check_permutation(a, {0, 2, 1, 4, 3});
  check_permutation(b, {});
  a.merge(a);
  a.merge(empty);
  EXPECT_THROW(a.rotate_left(0, 6, 0), std::out_of_range);
  EXPECT_THROW(a.rotate_left(3, 2, 0), std::out_of_range);
  EXPECT_THROW(a[5], std::out_of_range);
  auto& alias = a;
  a = std::move(alias);
  P moved(std::move(a));
  check_permutation(a, {});
  b = std::move(moved);
  check_permutation(moved, {});
  empty.merge(b);
  check_permutation(b, {});
  check_permutation(empty, {0, 2, 1, 4, 3});
  const auto n = Tree::block_capacity * 23 + 1;
  auto identity = P::identity(n);
  std::vector<T> expected(n);
  std::iota(expected.begin(), expected.end(), T{0});
  check_permutation(identity, expected);
}

TYPED_TEST(PermutationSpec, EveryShortRotation) {
  using P = TypeParam;
  using T = typename P::value_type;
  for (std::size_t n = 0; n <= 8; ++n) {
    for (std::size_t left = 0; left <= n; ++left) {
      for (std::size_t right = left; right <= n; ++right) {
        for (std::size_t d = 0; d <= right - left + 1; ++d) {
          auto permutation = P::identity(n);
          std::vector<T> expected(n);
          std::iota(expected.begin(), expected.end(), T{0});
          permutation.rotate_left(left, right, d);
          rotate_values(expected, left, right, d);
          check_permutation(permutation, expected);
        }
      }
    }
  }
}

TYPED_TEST(PermutationSpec, TwoLeafSpillPreflightsBeforeRebasing) {
  using P = TypeParam;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  for (std::size_t fail = 0; fail <= 1; ++fail) {
    auto a = P::identity(Tree::block_capacity);
    auto b = P::identity(1);
    const auto before = contents(a);
    const auto biases_a = a.test_tree().test_bias_snapshot();
    const auto biases_b = b.test_tree().test_bias_snapshot();
    Tree::test_reset_counters();
    Tree::test_fail_after(fail);
    if (fail == 0) {
      EXPECT_THROW(a.merge(b), std::bad_alloc);
    } else {
      EXPECT_NO_THROW(a.merge(b));
    }
    Tree::test_fail_after(-1);
    if (fail == 0) {
      EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
      EXPECT_EQ(a.test_tree().test_bias_snapshot(), biases_a);
      EXPECT_EQ(b.test_tree().test_bias_snapshot(), biases_b);
      check_permutation(a, before);
      check_permutation(b, {0});
    } else {
      auto expected = before;
      expected.push_back(static_cast<typename P::value_type>(before.size()));
      check_permutation(a, expected);
      check_permutation(b, {});
      EXPECT_EQ(Tree::test_counters.allocations, 1);
    }
  }
}

TYPED_TEST(PermutationSpec, RebasedChildRotationAndFailedPreflightKeepRawTags) {
  using P = TypeParam;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  const auto make = [] {
    auto a = P::identity(83 * Tree::block_capacity);
    auto b = P::identity(131 * Tree::block_capacity);
    a.merge(b);
    auto prefix = P::identity(47 * Tree::block_capacity);
    prefix.merge(a);
    return prefix;
  };
  auto permutation = make();
  auto expected = contents(permutation);
  auto biases = permutation.test_tree().test_bias_snapshot();
  ASSERT_TRUE(std::ranges::any_of(
      biases, [](const auto& entry) { return entry.second != 0; }));
  const auto by_address = [](const auto& a, const auto& b) {
    return std::less<const void*>{}(a.first, b.first);
  };
  std::ranges::sort(biases, by_address);
  const auto shapes = permutation.test_tree().test_child_boundaries();
  for (auto index : {std::size_t{0}, shapes.size() / 2, shapes.size() - 1}) {
    const auto& boundaries = shapes[index];
    const auto l = boundaries.front(), r = boundaries.back();
    const auto d = boundaries[1] - l;
    Tree::test_reset_counters();
    Tree::test_fail_after(0);
    EXPECT_NO_THROW(permutation.rotate_left(l, r, d));
    Tree::test_fail_after(-1);
    EXPECT_EQ(Tree::test_counters.allocations, 0);
    EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
    auto after = permutation.test_tree().test_bias_snapshot();
    std::ranges::sort(after, by_address);
    EXPECT_EQ(after, biases);
    rotate_values(expected, l, r, d);
    check_permutation(permutation, expected);
    Tree::test_fail_after(0);
    EXPECT_NO_THROW(permutation.rotate_left(l, r, r - l - d));
    Tree::test_fail_after(-1);
    rotate_values(expected, l, r, r - l - d);
  }
  std::size_t allocations;
  {
    auto probe = make();
    Tree::test_reset_counters();
    probe.rotate_left(1, probe.size() - 1, probe.size() / 3);
    allocations = Tree::test_counters.allocations;
  }
  const auto raw = permutation.test_tree().test_bias_snapshot();
  const auto leaves = permutation.test_tree().test_leaf_identities();
  const auto nodes = permutation.test_tree().test_internal_node_identities();
  for (std::size_t fail = 0; fail < allocations; ++fail) {
    Tree::test_reset_counters();
    const auto bytes = Tree::test_counters.live_bytes;
    Tree::test_fail_after(fail);
    EXPECT_THROW(permutation.rotate_left(1, permutation.size() - 1,
                                         permutation.size() / 3),
                 std::bad_alloc);
    Tree::test_fail_after(-1);
    EXPECT_EQ(Tree::test_counters.live_bytes, bytes);
    EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
    EXPECT_EQ(permutation.test_tree().test_bias_snapshot(), raw);
    EXPECT_EQ(permutation.test_tree().test_leaf_identities(), leaves);
    EXPECT_EQ(permutation.test_tree().test_internal_node_identities(), nodes);
    check_permutation(permutation, expected);
  }
}

TYPED_TEST(PermutationSpec, NestedRebasesAssociationsAndDeepRotations) {
  using P = TypeParam;
  using T = typename P::value_type;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  constexpr auto capacity = Tree::block_capacity;
  std::mt19937_64 random(91723);
  auto permutation = P::identity(capacity * 9 + 1);
  auto expected = contents(permutation);
  for (std::size_t iteration = 0; iteration < 90; ++iteration) {
    auto donor = P::identity(1 + random() % (capacity * 3));
    auto extra = P::identity(1 + random() % (capacity * 2));
    donor.rotate_left(0, donor.size(), random());
    donor.merge(extra);
    donor.rotate_left(0, donor.size(), random());
    const auto incoming = contents(donor);
    if (iteration % 3 == 0) {
      auto next = incoming;
      for (auto value : expected) {
        next.push_back(static_cast<T>(value + incoming.size()));
      }
      donor.merge(permutation);
      permutation = std::move(donor);
      expected = std::move(next);
    } else {
      const auto offset = expected.size();
      for (auto value : incoming) {
        expected.push_back(static_cast<T>(value + offset));
      }
      permutation.merge(donor);
      check_permutation(donor, {});
    }
    check_permutation(permutation, expected);
    for (std::size_t j = 0; j < 3; ++j) {
      const auto left = j == 0 ? 0 : random() % permutation.size();
      const auto right =
          j == 0 ? permutation.size()
                 : left + random() % (permutation.size() - left + 1);
      const auto distance = random();
      const auto h = permutation.test_tree().height();
      Tree::test_reset_counters();
      permutation.rotate_left(left, right, distance);
      const auto counts = Tree::test_counters;
      EXPECT_LE(counts.allocations, 21 * h + 27);
      EXPECT_LE(counts.node_visits, 160 * (h + 4));
      EXPECT_LE(counts.child_transfers, 160 * 16 * (h + 4));
      // Three splits + three seams: <=15 redistributions. Each can normalize
      // at most two leaves, so <=45 payload calls, independent of n and h.
      EXPECT_LE(counts.payload_mutations, 45);
      rotate_values(expected, left, right, distance);
      check_permutation(permutation, expected);
    }
  }
}

TYPED_TEST(PermutationSpec,
           AllocationFailureEveryIdentityMergeAndRotationPoint) {
  using P = TypeParam;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  constexpr auto capacity = Tree::block_capacity;
  Tree::test_reset_counters();
  { auto probe = P::identity(capacity * 19 + 1); }
  const auto construction_allocations = Tree::test_counters.allocations;
  for (std::size_t fail = 0; fail <= construction_allocations; ++fail) {
    Tree::test_fail_after(fail);
    if (fail < construction_allocations) {
      EXPECT_THROW(P::identity(capacity * 19 + 1), std::bad_alloc);
    } else {
      EXPECT_NO_THROW(P::identity(capacity * 19 + 1));
    }
    Tree::test_fail_after(-1);
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
    EXPECT_EQ(Tree::test_counters.live_bytes, 0);
  }
  auto make = [](std::size_t n) {
    auto result = P::identity(n);
    // Rebase already rebased subtrees, then rotate without normalizing every
    // off-path tag. Both transaction participants must retain nonzero biases.
    for (int round = 0; round < 3; ++round) {
      auto prefix = P::identity(capacity * 17);
      prefix.merge(result);
      result = std::move(prefix);
    }
    result.rotate_left(0, result.size(), 1);
    return result;
  };
  for (int operation = 0; operation < 5; ++operation) {
    auto apply = [operation](P& a, P& b) {
      if (operation < 3) {
        a.merge(b);
      } else if (operation == 3) {
        a.rotate_left(0, a.size(), a.size() / 2 + 1);
      } else {
        a.rotate_left(1, a.size() - 1, a.size() / 3);
      }
    };
    const auto left_size = operation == 0 ? capacity : capacity * 29 + 1;
    const auto right_size = operation == 1 ? capacity * 31 : capacity + 1;
    std::size_t allocations;
    {
      auto a = make(left_size);
      auto b = make(right_size);
      Tree::test_reset_counters();
      apply(a, b);
      allocations = Tree::test_counters.allocations;
    }
    for (std::size_t fail = 0; fail <= allocations; ++fail) {
      auto a = make(left_size);
      auto b = make(right_size);
      const auto before_a = contents(a);
      const auto before_b = contents(b);
      const auto leaves_a = a.test_tree().test_leaf_identities();
      const auto nodes_a = a.test_tree().test_internal_node_identities();
      const auto leaves_b = b.test_tree().test_leaf_identities();
      const auto nodes_b = b.test_tree().test_internal_node_identities();
      const auto biases_a = a.test_tree().test_bias_snapshot();
      const auto biases_b = b.test_tree().test_bias_snapshot();
      ASSERT_TRUE(std::ranges::any_of(
          biases_a, [](const auto& entry) { return entry.second != 0; }));
      ASSERT_TRUE(std::ranges::any_of(
          biases_b, [](const auto& entry) { return entry.second != 0; }));
      const auto live = Tree::test_counters.live_allocations;
      const auto bytes = Tree::test_counters.live_bytes;
      Tree::test_reset_counters();
      Tree::test_fail_after(fail);
      if (fail < allocations) {
        EXPECT_THROW(apply(a, b), std::bad_alloc);
      } else {
        EXPECT_NO_THROW(apply(a, b));
      }
      Tree::test_fail_after(-1);
      if (fail < allocations) {
        EXPECT_EQ(Tree::test_counters.live_allocations, live);
        EXPECT_EQ(Tree::test_counters.live_bytes, bytes);
        EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
        EXPECT_EQ(a.test_tree().test_leaf_identities(), leaves_a);
        EXPECT_EQ(a.test_tree().test_internal_node_identities(), nodes_a);
        EXPECT_EQ(b.test_tree().test_leaf_identities(), leaves_b);
        EXPECT_EQ(b.test_tree().test_internal_node_identities(), nodes_b);
        EXPECT_EQ(a.test_tree().test_bias_snapshot(), biases_a);
        EXPECT_EQ(b.test_tree().test_bias_snapshot(), biases_b);
        check_permutation(a, before_a);
        check_permutation(b, before_b);
      } else {
        auto expected = before_a;
        if (operation < 3) {
          for (auto value : before_b) {
            expected.push_back(
                static_cast<typename P::value_type>(value + before_a.size()));
          }
          check_permutation(b, {});
        } else {
          const auto left = operation == 3 ? 0 : 1;
          const auto right =
              operation == 3 ? expected.size() : expected.size() - 1;
          const auto distance =
              operation == 3 ? expected.size() / 2 + 1 : expected.size() / 3;
          rotate_values(expected, left, right, distance);
          check_permutation(b, before_b);
        }
        check_permutation(a, expected);
      }
    }
  }
}

TEST(PermutationBoundary, FullUint16DomainAndOverflowBeforeAllocation) {
  using P = Permutation<std::uint16_t>;
  using Tree =
      std::remove_cvref_t<decltype(std::declval<const P&>().test_tree())>;
  auto a = P::identity(65535);
  auto b = P::identity(1);
  a.merge(b);
  EXPECT_TRUE(b.empty());
  ASSERT_EQ(a.size(), 65536);
  EXPECT_EQ(a[65535], 65535);
  a.rotate_left(0, a.size(), 32769);
  std::vector<std::uint16_t> expected(65536);
  std::iota(expected.begin(), expected.end(), std::uint16_t{0});
  rotate_values(expected, 0, expected.size(), 32769);
  check_permutation(a, expected);
  b = P::identity(1);
  const auto biases_a = a.test_tree().test_bias_snapshot();
  const auto biases_b = b.test_tree().test_bias_snapshot();
  Tree::test_fail_after(0);
  EXPECT_THROW(a.merge(b), std::length_error);
  EXPECT_THROW(P::identity(65537), std::length_error);
  Tree::test_fail_after(-1);
  EXPECT_EQ(a.test_tree().test_bias_snapshot(), biases_a);
  EXPECT_EQ(b.test_tree().test_bias_snapshot(), biases_b);
  check_permutation(a, expected);
  check_permutation(b, {0});
  auto entire = P::identity(65536);
  EXPECT_EQ(entire[65535], 65535);
}

TEST(PermutationLayout, OptionalBiasHasNoUntaggedLayoutTax) {
  using Block = PackedValueBlock<std::uint64_t>;
  static_assert(sizeof(Block) == 256);
  static_assert(SequenceTree<Block, 4>::node_storage_bytes == 64);
  static_assert(SequenceTree<Block, 8>::node_storage_bytes == 128);
  static_assert(SequenceTree<Block, 16>::node_storage_bytes == 256);
  static_assert(SequenceTree<Block>::block_storage_bytes == sizeof(Block));
  static_assert(Block::capacity * 64 == 1920);
  static_assert(SequenceTree<Block, 4, LengthLayout::cumulative,
                             true>::node_storage_bytes == 128);
  static_assert(SequenceTree<Block, 16, LengthLayout::cumulative,
                             true>::node_storage_bytes == 320);
  static_assert(
      SequenceTree<PackedValueBlock<std::uint64_t>, 8, LengthLayout::cumulative,
                   true>::node_storage_bytes == 192);
  static_assert(
      SequenceTree<PackedValueBlock<std::uint64_t>, 8, LengthLayout::cumulative,
                   true>::block_storage_bytes == 320);
  // F6 uses existing alignment slack rather than another aligned node line.
  // Accounting must also handle a zero incremental node-byte category.
  auto permutation = Permutation<std::uint32_t, 512, 6>::identity(100);
  check_permutation(permutation, contents(permutation));
}

template <class Tree>
concept RawBlockTraversal =
    requires(const Tree& tree) { tree.for_each_block([](const auto&) {}); };
static_assert(RawBlockTraversal<SequenceTree<PackedValueBlock<std::uint64_t>>>);
static_assert(!RawBlockTraversal<SequenceTree<PackedValueBlock<std::uint64_t>,
                                              8,
                                              LengthLayout::cumulative,
                                              true>>);

TEST(PermutationStress, ThousandsOfSingletonMerges) {
  using P = Permutation<std::uint16_t, 512>;
  P permutation;
  for (std::size_t i = 0; i < 3000; ++i) {
    auto singleton = P::identity(1);
    permutation.merge(singleton);
    ASSERT_TRUE(singleton.empty());
    ASSERT_EQ(permutation[i], i);
  }
  permutation.rotate_left(1, permutation.size() - 1, 701);
  std::vector<std::uint16_t> expected(permutation.size());
  std::iota(expected.begin(), expected.end(), std::uint16_t{0});
  rotate_values(expected, 1, expected.size() - 1, 701);
  check_permutation(permutation, expected);
}

}  // namespace
