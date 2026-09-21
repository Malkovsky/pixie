#include <gtest/gtest.h>
#include <pixie/detail/sequence/bit_block.h>
#include <pixie/detail/sequence/packed_value_block.h>
#include <pixie/detail/sequence/sequence_tree.h>

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

// Keep GCC 13's loop specialization from diagnosing impossible huge memcpy
// lengths in the exhaustive short-range tests (also affects std::rotate).
template <class Block>
[[gnu::noinline]] void rotate_block(Block& block,
                                    std::size_t left,
                                    std::size_t right,
                                    std::size_t distance) {
  block.rotate_left(left, right, distance);
}

template <class Container>
auto contents(const Container& container) {
  std::vector<typename Container::value_type> result;
  for (std::size_t i = 0; i < container.size(); ++i) {
    result.push_back(container[i]);
  }
  return result;
}

TEST(PackedFieldRead, EveryShortFieldOriginWordAndRingCrossing) {
  std::mt19937_64 random(271828);
  for (const auto n : {0u, 1u, 3u, 63u, 64u, 65u, 127u, 128u, 129u, 192u}) {
    std::array<std::uint64_t, 3> words{random(), random(), random()};
    if (n % 64) {
      words[n / 64] |= ~std::uint64_t{0} << (n % 64);
    }
    BitBlock<192> block(words, n);
    EXPECT_EQ(block.read_bits(n, 0), 0);
    for (std::size_t origin = 0; origin < n; ++origin) {
      for (std::size_t start = 0; start <= n; ++start) {
        std::uint64_t expected = 0;
        for (std::size_t width = 0;
             width <= std::min<std::size_t>(64, n - start); ++width) {
          if (width != 0) {
            const auto bit = (origin + start + width - 1) % n;
            expected |= ((words[bit / 64] >> (bit % 64)) & 1) << (width - 1);
          }
          ASSERT_EQ(block.read_bits(start, width), expected)
              << "n=" << n << " origin=" << origin << " start=" << start
              << " width=" << width;
        }
      }
      block.rotate_left(0, n, 1);
    }
  }
}

template <class T>
class PackedValueSpec : public ::testing::Test {};
using ValueBlocks = ::testing::Types<PackedValueBlock<bool, 512>,
                                     PackedValueBlock<std::uint64_t, 512, 1>,
                                     PackedValueBlock<std::uint64_t, 512, 3>,
                                     PackedValueBlock<std::uint64_t, 512, 5>,
                                     PackedValueBlock<std::uint8_t, 512>,
                                     PackedValueBlock<std::uint16_t, 512>,
                                     PackedValueBlock<std::uint32_t, 512>,
                                     PackedValueBlock<std::uint64_t, 512>>;
TYPED_TEST_SUITE(PackedValueSpec, ValueBlocks);

TYPED_TEST(PackedValueSpec, LayoutBoundsEveryShortRotationAndWrappedReads) {
  using Block = TypeParam;
  using T = typename Block::value_type;
  static_assert(SequenceBlock<Block>);
  static_assert(sizeof(Block) == 64);
  static_assert(alignof(Block) == 64);
  constexpr auto maximum = ~std::uint64_t{0} >> (64 - Block::width);
  std::mt19937_64 random(42);
  std::array<T, Block::capacity> fields{};
  for (auto& field : fields) {
    field = static_cast<T>(random() & maximum);
  }
  fields.back() = static_cast<T>(maximum);
  Block empty;
  EXPECT_TRUE(empty.empty());
  EXPECT_THROW(empty[0], std::out_of_range);
  EXPECT_EQ(empty.payload_capacity_bytes() + empty.metadata_bytes(),
            sizeof(Block));
  for (std::size_t n = 0; n <= std::min<std::size_t>(8, Block::capacity); ++n) {
    for (std::size_t left = 0; left <= n; ++left) {
      for (std::size_t right = left; right <= n; ++right) {
        for (std::size_t d = 0; d <= right - left + 1; ++d) {
          Block block(std::span<const T>(fields.data(), n));
          std::vector<T> expected(fields.begin(), fields.begin() + n);
          rotate_block(block, left, right, d);
          rotate_values(expected, left, right, d);
          EXPECT_EQ(contents(block), expected);
        }
      }
    }
  }
  Block block{std::span<const T>(fields)};
  std::vector<T> expected(fields.begin(), fields.end());
  for (std::size_t i = 0; i < Block::capacity; ++i) {
    block.rotate_left(0, block.size(), 1);
    rotate_values(expected, 0, expected.size(), 1);
    EXPECT_EQ(contents(block), expected);
    auto partial = block;
    auto oracle = expected;
    partial.rotate_left(1, partial.size() - 1,
                        std::numeric_limits<std::size_t>::max());
    rotate_values(oracle, 1, oracle.size() - 1,
                  std::numeric_limits<std::size_t>::max());
    EXPECT_EQ(contents(partial), oracle);
  }
  EXPECT_THROW(block[block.size()], std::out_of_range);
  std::array<T, Block::capacity + 1> oversized{};
  EXPECT_THROW((Block{std::span<const T>(oversized)}), std::invalid_argument);
  if constexpr (Block::width < std::numeric_limits<T>::digits) {
    const std::array<T, 1> invalid{static_cast<T>(maximum + 1)};
    EXPECT_THROW((Block{std::span<const T>(invalid)}), std::invalid_argument);
  }
}

TYPED_TEST(PackedValueSpec, WrappedRedistributionAndBiasReencoding) {
  using Block = TypeParam;
  using T = typename Block::value_type;
  constexpr auto maximum = ~std::uint64_t{0} >> (64 - Block::width);
  std::array<T, Block::capacity> fields{};
  for (std::size_t i = 0; i < fields.size(); ++i) {
    fields[i] = static_cast<T>(i & maximum);
  }
  Block source{std::span<const T>(fields)};
  source.rotate_left(0, source.size(), 3);
  const auto original = contents(source);
  for (std::size_t cut = 0; cut <= source.size(); ++cut) {
    auto left = source;
    Block right;
    left.redistribute(right, cut);
    EXPECT_EQ(contents(left),
              (std::vector<T>(original.begin(), original.begin() + cut)));
    EXPECT_EQ(contents(right),
              (std::vector<T>(original.begin() + cut, original.end())));
    left.redistribute(right, source.size());
    EXPECT_TRUE(right.empty());
    EXPECT_EQ(contents(left), original);
  }
  const auto half = Block::capacity / 2;
  Block a(std::span<const T>(fields.data(), half));
  Block b(std::span<const T>(fields.data(), half + 1));
  a.rotate_left(0, a.size(), 1);
  b.rotate_left(0, b.size(), 2);
  auto expected = contents(a);
  const auto second = contents(b);
  expected.insert(expected.end(), second.begin(), second.end());
  a.redistribute(b, Block::capacity);
  auto joined = contents(a);
  const auto tail = contents(b);
  joined.insert(joined.end(), tail.begin(), tail.end());
  EXPECT_EQ(joined, expected);

  fields.fill(T{0});
  Block biased{std::span<const T>(fields)};
  biased.rotate_left(0, biased.size(), 3);
  biased.add_bias(maximum);
  for (std::size_t i = 0; i < biased.size(); ++i) {
    EXPECT_EQ(biased[i], static_cast<T>(maximum));
  }
  biased.add_bias(0);
  EXPECT_EQ(biased[0], static_cast<T>(maximum));
}

template <class T>
class TaggedTreeSpec : public ::testing::Test {};
using TaggedTrees =
    ::testing::Types<SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  4,
                                  LengthLayout::cumulative,
                                  true>,
                     SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  4,
                                  LengthLayout::individual,
                                  true>,
                     SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  8,
                                  LengthLayout::cumulative,
                                  true>,
                     SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  8,
                                  LengthLayout::individual,
                                  true>,
                     SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  16,
                                  LengthLayout::cumulative,
                                  true>,
                     SequenceTree<PackedValueBlock<std::uint64_t, 512>,
                                  16,
                                  LengthLayout::individual,
                                  true>>;
TYPED_TEST_SUITE(TaggedTreeSpec, TaggedTrees);

template <class Tree>
Tree make_tagged(std::size_t n) {
  using Block = typename Tree::block_type;
  auto blocks =
      std::views::iota(std::size_t{0},
                       n / Block::capacity + (n % Block::capacity != 0)) |
      std::views::transform([n](std::size_t index) {
        std::array<std::uint64_t, Block::capacity> fields;
        const auto start = index * Block::capacity;
        const auto count = std::min(Block::capacity, n - start);
        std::iota(fields.begin(), fields.end(), start);
        return Block(std::span<const std::uint64_t>(fields.data(), count));
      });
  return Tree::from_blocks(blocks);
}

TYPED_TEST(TaggedTreeSpec, FullWidthBiasConstReadsNestedSplitsAndRejoins) {
  using Tree = TypeParam;
  constexpr auto maximum = std::numeric_limits<std::uint64_t>::max();
  constexpr std::size_t n = Tree::block_capacity * 103 + 1;
  auto tree = make_tagged<Tree>(n);
  tree.test_add_bias(std::uint64_t{1} << 63);
  tree.test_add_bias(maximum - (std::uint64_t{1} << 63) - (n - 1));
  std::vector<std::uint64_t> expected(n);
  std::iota(expected.begin(), expected.end(), maximum - (n - 1));
  const auto leaves = tree.test_leaf_identities();
  const auto nodes = tree.test_internal_node_identities();
  const auto biases = tree.test_bias_snapshot();
  Tree::test_reset_counters();
  EXPECT_EQ(contents(std::as_const(tree)), expected);
  EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
  EXPECT_EQ(Tree::test_counters.allocations, 0);
  EXPECT_EQ(tree.test_leaf_identities(), leaves);
  EXPECT_EQ(tree.test_internal_node_identities(), nodes);
  EXPECT_EQ(tree.test_bias_snapshot(), biases);
  ASSERT_TRUE(tree.test_validate());
  std::mt19937_64 random(173);
  for (std::size_t iteration = 0; iteration < 70; ++iteration) {
    const auto cut =
        iteration < Tree::block_capacity ? iteration : random() % n;
    const auto h = tree.height();
    Tree::test_reset_counters();
    auto right = tree.split_off(cut);
    EXPECT_LE(Tree::test_counters.allocations, 3 * h + 1);
    EXPECT_LE(Tree::test_counters.payload_mutations, 3);
    ASSERT_TRUE(tree.test_validate());
    ASSERT_TRUE(right.test_validate());
    EXPECT_EQ(contents(tree), (std::vector<std::uint64_t>(
                                  expected.begin(), expected.begin() + cut)));
    EXPECT_EQ(contents(right), (std::vector<std::uint64_t>(
                                   expected.begin() + cut, expected.end())));
    tree.merge(right);
    EXPECT_TRUE(right.empty());
    const auto left = random() % n;
    const auto end = left + random() % (n - left + 1);
    const auto distance = random();
    tree.rotate_left(left, end, distance);
    rotate_values(expected, left, end, distance);
    ASSERT_TRUE(tree.test_validate());
    EXPECT_EQ(contents(tree), expected);
  }
}

TYPED_TEST(TaggedTreeSpec, ChildRotationRetainsUniformParentBias) {
  using Tree = TypeParam;
  auto tree = make_tagged<Tree>(257 * Tree::block_capacity);
  tree.test_add_bias(std::uint64_t{1} << 63);
  const auto raw = tree.test_bias_snapshot();
  const auto boundaries = tree.test_child_boundaries().front();
  auto expected = contents(tree);
  Tree::test_reset_counters();
  Tree::test_fail_after(0);
  EXPECT_NO_THROW(tree.rotate_left(0, tree.size(), boundaries[1]));
  Tree::test_fail_after(-1);
  EXPECT_EQ(Tree::test_counters.allocations, 0);
  EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
  const auto after = tree.test_bias_snapshot();
  ASSERT_EQ(after.size(), raw.size());
  EXPECT_EQ(after.front(), raw.front());
  for (const auto& entry : after) {
    EXPECT_NE(std::ranges::find(raw, entry), raw.end());
  }
  rotate_values(expected, 0, expected.size(), boundaries[1]);
  EXPECT_EQ(contents(tree), expected);
  EXPECT_TRUE(tree.test_validate());
}

TYPED_TEST(TaggedTreeSpec, TaggedLeafFittingMergeAndOffPathOwnership) {
  using Tree = TypeParam;
  auto a = make_tagged<Tree>(2);
  auto b = make_tagged<Tree>(3);
  a.test_add_bias(17);
  b.test_add_bias(std::numeric_limits<std::uint64_t>::max() - 2);
  auto expected = contents(a);
  const auto donor = contents(b);
  expected.insert(expected.end(), donor.begin(), donor.end());
  Tree::test_fail_after(0);
  a.rotate_left(0, 2, 1);
  a.merge(b);
  Tree::test_fail_after(-1);
  std::swap(expected[0], expected[1]);
  EXPECT_EQ(contents(a), expected);
  ASSERT_TRUE(a.test_validate());
  auto deep = make_tagged<Tree>(Tree::block_capacity * 1000);
  deep.test_add_bias(12345);
  const auto leaves = deep.test_leaf_identities();
  const auto nodes = deep.test_internal_node_identities();
  const auto h = deep.height();
  Tree::test_reset_counters();
  deep.rotate_left(1, deep.size() - 1, deep.size() / 3);
  const auto counts = Tree::test_counters;
  EXPECT_EQ(counts.allocations, 15 * h + 27);
  EXPECT_LE(counts.node_visits, 160 * (h + 4));
  EXPECT_LE(counts.payload_mutations, 45);
  const auto after = deep.test_leaf_identities();
  std::size_t surviving = 0;
  for (auto p : leaves) {
    surviving += std::find(after.begin(), after.end(), p) != after.end();
  }
  EXPECT_GE(surviving + 12, leaves.size());
  const auto after_nodes = deep.test_internal_node_identities();
  surviving = 0;
  for (auto p : nodes) {
    surviving += std::find(after_nodes.begin(), after_nodes.end(), p) !=
                 after_nodes.end();
  }
  EXPECT_GE(surviving + counts.node_visits, nodes.size());
  ASSERT_TRUE(deep.test_validate());
}

TYPED_TEST(TaggedTreeSpec, SplitFailureEveryPointAndMultipleCutsInsideOneLeaf) {
  using Tree = TypeParam;
  constexpr auto capacity = Tree::block_capacity;
  constexpr auto n = capacity * 43 + 1;
  auto make = [] {
    auto tree = make_tagged<Tree>(n);
    tree.test_add_bias(std::uint64_t{1} << 62);
    tree.rotate_left(0, n, 1);
    tree.test_add_bias(std::uint64_t{1} << 62);
    return tree;
  };
  std::size_t allocations;
  {
    auto tree = make();
    Tree::test_reset_counters();
    auto right = tree.split_off(capacity + 1);
    allocations = Tree::test_counters.allocations;
  }
  for (std::size_t fail = 0; fail <= allocations; ++fail) {
    auto tree = make();
    const auto expected = contents(tree);
    const auto leaves = tree.test_leaf_identities();
    const auto nodes = tree.test_internal_node_identities();
    const auto biases = tree.test_bias_snapshot();
    Tree right;
    const auto right_biases = right.test_bias_snapshot();
    const auto bytes = Tree::test_counters.live_bytes;
    Tree::test_reset_counters();
    Tree::test_fail_after(fail);
    if (fail < allocations) {
      EXPECT_THROW(right = tree.split_off(capacity + 1), std::bad_alloc);
    } else {
      EXPECT_NO_THROW(right = tree.split_off(capacity + 1));
    }
    Tree::test_fail_after(-1);
    if (fail < allocations) {
      EXPECT_EQ(Tree::test_counters.live_bytes, bytes);
      EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
      EXPECT_EQ(tree.test_leaf_identities(), leaves);
      EXPECT_EQ(tree.test_internal_node_identities(), nodes);
      EXPECT_EQ(tree.test_bias_snapshot(), biases);
      EXPECT_EQ(right.test_bias_snapshot(), right_biases);
    } else {
      tree.merge(right);
    }
    EXPECT_EQ(contents(tree), expected);
    ASSERT_TRUE(tree.test_validate());
  }
  for (auto distance : {std::size_t{1}, capacity - 2}) {
    auto tree = make();
    auto expected = contents(tree);
    tree.rotate_left(1, n - 1, distance);
    rotate_values(expected, 1, n - 1, distance);
    ASSERT_TRUE(tree.test_validate());
    EXPECT_EQ(contents(tree), expected);
  }
}

}  // namespace
