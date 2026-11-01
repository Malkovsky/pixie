#include <gtest/gtest.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/bp_micro_tables.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/experimental/tree_prepared_events.h>
#include <pixie/tree/experimental/tree_prepared_root.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/path_view.h>
#include <pixie/tree/traversal_context.h>
#include <pixie/utils.h>
#include <references/pointer_tree.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <future>
#include <limits>
#include <memory>
#include <optional>
#include <queue>
#include <random>
#include <stack>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../benchmarks/shifted_query_pool.h"
#include "../benchmarks/tree_dfs/optimized.h"
#include "../benchmarks/tree_traversal_workloads.h"

namespace {

template <class RankSelect>
struct LoudsCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes - 1; }
  using Tree = pixie::LoudsTree<RankSelect>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return bfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return adj_to_louds(tree_size, adjacency);
  }
};

template <class RmM>
struct BpCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes; }
  using Tree = pixie::BPTree<RmM>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return dfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return pixie::adj_to_bp(tree_size, adjacency);
  }
};

template <class RmM>
struct DfudsCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes - 1; }
  using Tree = pixie::DFUDSTree<RmM>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return dfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return pixie::adj_to_dfuds(tree_size, adjacency);
  }
};

template <class Case>
struct ViewCase : Case {
  using Tree = typename Case::Tree::TreeView;
};

template <class Case>
class TreeSpecificationTest : public testing::Test {};

using TreeCases = testing::Types<
    ViewCase<LoudsCase<pixie::RankSelectSupport<>>>,
    ViewCase<LoudsCase<pixie::RankSelectSupportCompressed<>>>,
    ViewCase<BpCase<pixie::RmMBTree<>>>,
    ViewCase<DfudsCase<pixie::RmMBTree<>>>,
    ViewCase<BpCase<pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>>>,
    ViewCase<DfudsCase<pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>>>,
    LoudsCase<pixie::RankSelectSupport<>>,
    LoudsCase<pixie::RankSelectSupportCompressed<>>,
    BpCase<pixie::RmMBTree<>>,
    DfudsCase<pixie::RmMBTree<>>,
    BpCase<pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>>,
    DfudsCase<pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>>>;
TYPED_TEST_SUITE(TreeSpecificationTest, TreeCases);

TYPED_TEST(TreeSpecificationTest, IndependentConcurrentViews) {
  using Tree = typename TypeParam::Tree;
  static_assert(pixie::TreeNavigationView<typename Tree::NodeView>);
  static_assert(
      std::is_constructible_v<Tree, std::vector<uint64_t>&&, size_t> ==
      !std::is_same_v<Tree, typename Tree::TreeView>);
  static_assert(
      std::is_constructible_v<Tree, const std::vector<uint64_t>&&, size_t> ==
      !std::is_same_v<Tree, typename Tree::TreeView>);
  std::mt19937_64 rng(42);
  constexpr size_t n = 513;
  const auto adjacency = TypeParam::order(n, generate_random_tree(n, rng));
  const auto words = TypeParam::encode(n, adjacency);
  const Tree tree(std::span<const uint64_t>(words), n);
  auto inspect = [&] {
    bool correct = true;
    auto visit = [&](pixie::TreeNode node, bool entering) {
      if (!entering) {
        return;
      }
      auto view = tree.navigate(node);
      auto independent = view;
      const auto degree = tree.degree(node);
      correct &= view.degree() == degree;
      size_t index = 0;
      view.for_each_child([&](auto child) {
        const auto expected = tree.child(node, index++);
        correct &= child.number == expected.number && child.pos == expected.pos;
      });
      correct &= index == degree;
      if (degree) {
        view.reset(view.first_child());
      }
      correct &= independent.node().number == node.number;
      correct &= independent.degree() == degree;
      correct &= independent.parent().pos == tree.parent(node).pos;
    };
    pixie::detail::tree_unpack::subtree(tree, tree.root(), visit);
    return correct;
  };
  auto first = std::async(std::launch::async, inspect);
  auto second = std::async(std::launch::async, inspect);
  EXPECT_TRUE(first.get());
  EXPECT_TRUE(second.get());
}

TYPED_TEST(TreeSpecificationTest, EncodingOwnershipCopyAndMove) {
  using Tree = typename TypeParam::Tree;
  using View = typename Tree::TreeView;
  constexpr bool owns = !std::is_same_v<Tree, View>;
  static_assert(std::is_same_v<typename View::TreeView, View>);
  static_assert(
      std::is_same_v<typename Tree::PathView, pixie::TreePathView<Tree>>);
  static_assert(std::is_same_v<typename Tree::template BoundedPathView<2>,
                               pixie::TreePathView<Tree, 2>>);
  std::mt19937_64 rng(91);
  constexpr size_t n = 257;
  const auto adjacency = TypeParam::order(n, generate_random_tree(n, rng));
  auto words = TypeParam::encode(n, adjacency);
  const auto expected_words = words;
  auto check = [&](const auto& tree) {
    EXPECT_EQ(tree.size(), n);
    EXPECT_TRUE(std::ranges::equal(tree.encoding_words(), expected_words));
    std::vector<pixie::TreeNode> pending{tree.root()};
    size_t visited = 0;
    while (!pending.empty()) {
      const auto node = pending.back();
      pending.pop_back();
      ++visited;
      EXPECT_EQ(tree.degree(node), adjacency[node.number].size() - 1);
      if (node.number) {
        EXPECT_EQ(tree.parent(node).number, adjacency[node.number][0]);
      }
      tree.for_each_child(node, [&](auto child) { pending.push_back(child); });
    }
    EXPECT_EQ(visited, n);
    typename std::remove_cvref_t<decltype(tree)>::PathView path(tree);
    path.move_to_first_child();
    path.move_to_parent();
    EXPECT_TRUE(path.is_root());
  };
  auto tree = std::make_unique<Tree>(words, n);
  EXPECT_EQ(tree->encoding_words().data() == words.data(), !owns);
  EXPECT_EQ(tree->encoding_owned_bytes() == 0, !owns);
  if constexpr (pixie::HasMemoryUsageBytes<Tree>) {
    EXPECT_GE(tree->memory_usage_bytes(),
              sizeof(Tree) + tree->encoding_owned_bytes());
  }
  Tree copy(*tree);
  EXPECT_EQ(copy.encoding_words().data() == tree->encoding_words().data(),
            !owns);
  Tree assigned(words, n);
  assigned = *tree;
  tree.reset();
  if constexpr (owns) {
    words.assign(words.size(), 0);
    words.clear();
    words.shrink_to_fit();
  }
  check(copy);
  check(assigned);
  const auto address = copy.encoding_words().data();
  Tree moved(std::move(copy));
  EXPECT_EQ(moved.encoding_words().data(), address);
  check(moved);
  // Moved-from trees support assignment; restore both ownership modes.
  copy = assigned;
  check(copy);
  assigned = std::move(moved);
  EXPECT_EQ(assigned.encoding_words().data(), address);
  check(assigned);
  const Tree& self = assigned;
  assigned = self;
  EXPECT_EQ(assigned.encoding_words().data(), address);
  check(assigned);
  View view(assigned.encoding_words(), n);
  EXPECT_EQ(view.encoding_owned_bytes(), 0u);
  check(view);
  pixie::TreeTraversalContext context(view, 32768);
  EXPECT_EQ(context.source_words().data(), view.encoding_words().data());
  if constexpr (owns) {
    Tree from_lvalue = [&] {
      pixie::AlignedStorage source{std::span<const uint64_t>(expected_words)};
      Tree result(source, n);
      EXPECT_NE(result.encoding_words().data(), source.as_words64().data());
      std::ranges::fill(source.writable_bytes(), std::byte{0});
      return result;
    }();
    check(from_lvalue);
    pixie::AlignedStorage storage(TypeParam::bit_size(n));
    auto padded = storage.as_lines();
    std::memcpy(std::as_writable_bytes(padded).data(), expected_words.data(),
                expected_words.size() * 8);
    const auto data = storage.padded_view().as_words64().data();
    Tree consumed(std::move(storage), n);
    EXPECT_EQ(consumed.encoding_words().data(), data);
    check(consumed);
    // Temporary word vectors are safe for the copying owning constructor.
    Tree temporary(std::vector<uint64_t>(expected_words), n);
    check(temporary);
  }
}

TYPED_TEST(TreeSpecificationTest, RejectsInvalidEncodingStorage) {
  using Tree = typename TypeParam::Tree;
  using View = typename Tree::TreeView;
  const std::vector<uint64_t> empty;
  EXPECT_THROW((Tree(empty, 2)), std::invalid_argument);
  EXPECT_THROW((Tree(empty, std::numeric_limits<size_t>::max())),
               std::length_error);
  EXPECT_THROW((View(pixie::ReadOnlyStorageView{}, 1)), std::invalid_argument);
  if constexpr (!std::is_same_v<Tree, View>) {
    EXPECT_THROW((Tree(pixie::AlignedStorage{}, 1)), std::invalid_argument);
    // Allocation padding is not additional logical encoding capacity.
    pixie::AlignedStorage short_source(8);
    ASSERT_GE(short_source.padded_size_bytes(), 16u);
    EXPECT_THROW((Tree(std::move(short_source), 33)), std::invalid_argument);
  }
  alignas(8) std::array<std::byte, 24> bytes{};
  EXPECT_THROW((View(pixie::ReadOnlyStorageView(
                         std::span<const std::byte>(bytes).subspan(1, 8)),
                     1)),
               std::invalid_argument);
  EXPECT_THROW((View(pixie::ReadOnlyStorageView(
                         std::span<const std::byte>(bytes).first(1)),
                     1)),
               std::invalid_argument);
}

TYPED_TEST(TreeSpecificationTest, EncodingStorageWordBoundaries) {
  using Tree = typename TypeParam::Tree;
  constexpr bool owns = !std::is_same_v<Tree, typename Tree::TreeView>;
  std::mt19937_64 rng(942);
  for (size_t n : {0u, 1u, 31u, 32u, 33u, 63u, 64u, 65u}) {
    SCOPED_TRACE(n);
    std::vector<uint64_t> words;
    std::vector<std::vector<size_t>> adjacency;
    if (n) {
      adjacency = TypeParam::order(n, generate_random_tree(n, rng));
      words = TypeParam::encode(n, adjacency);
    }
    const size_t bits = n ? TypeParam::bit_size(n) : 0;
    pixie::AlignedStorage storage(bits);
    if (!words.empty()) {
      std::ranges::copy(std::as_bytes(std::span(words)),
                        std::as_writable_bytes(storage.as_lines()).begin());
    }
    const auto address = storage.padded_view().as_words64().data();
    const Tree tree = [&] {
      if constexpr (owns) {
        return Tree(std::move(storage), n);
      } else {
        return Tree(storage.padded_view(), n);
      }
    }();
    EXPECT_EQ(tree.size(), n);
    EXPECT_EQ(tree.encoding_words().data(), address);
    EXPECT_EQ(tree.encoding_words().size(), bits / 64 + (bits % 64 != 0));
    EXPECT_TRUE(std::ranges::equal(tree.encoding_words(), words));
    EXPECT_EQ(tree.encoding_owned_bytes(),
              owns ? ((bits + 511) / 512) * 64 : 0);
    if (n) {
      EXPECT_EQ(tree.degree(tree.root()), adjacency[0].size() - 1);
      size_t i = 1;
      tree.for_each_child(tree.root(), [&](auto child) {
        ASSERT_LT(i, adjacency[0].size());
        EXPECT_EQ(child.number, adjacency[0][i++]);
        EXPECT_EQ(tree.parent(child).number, 0u);
      });
      EXPECT_EQ(i, adjacency[0].size());
    }
    // Trailing source words must neither become part of the encoding nor
    // increase an owning tree's allocation.
    const auto expected = words;
    words.resize(words.size() + 16, ~uint64_t{0});
    const Tree padded_input(words, n);
    EXPECT_TRUE(std::ranges::equal(padded_input.encoding_words(), expected));
    EXPECT_EQ(padded_input.encoding_owned_bytes(), tree.encoding_owned_bytes());
  }
}

struct RejectZeroSelectRmM : pixie::RmMBTree<> {
  using Base = pixie::RmMBTree<>;
  using Base::Base;
  size_t select0(size_t rank) const {
    if (rank == 0) {
      throw std::logic_error("zero select rank");
    }
    return Base::select0(rank);
  }
};

TEST(DfudsNavigationTest, RootParentDoesNotQueryZeroSelectRank) {
  const std::vector<std::vector<size_t>> adjacency = {
      {0, 1, 2, 3}, {0}, {0}, {0}};
  const auto words = pixie::adj_to_dfuds(4, adjacency);
  const pixie::DFUDSTree<RejectZeroSelectRmM> tree(words, 4);
  for (size_t i = 0; i < 3; ++i) {
    const auto child = tree.child(tree.root(), i);
    EXPECT_EQ(tree.parent(child).number, 0u);
    EXPECT_EQ(tree.parent(child).pos, 0u);
    EXPECT_EQ(tree.navigate(child).parent().pos, 0u);
  }
}

template <class Case>
class BoundedTraversalTest : public testing::Test {};
using BoundedCases =
    testing::Types<BpCase<pixie::RmMBTree<>>,
                   DfudsCase<pixie::RmMBTree<>>,
                   LoudsCase<pixie::RankSelectSupportCompressed<>>>;
TYPED_TEST_SUITE(BoundedTraversalTest, BoundedCases);

TYPED_TEST(BoundedTraversalTest, ContextOwnershipAndAccounting) {
  using Tree = typename TypeParam::Tree;
  using Context = pixie::TreeTraversalContext<Tree>;
  static_assert(!std::is_constructible_v<Context, Tree&&>);
  static_assert(
      !std::is_constructible_v<Context, const Tree&, std::vector<uint64_t>&&>);
  std::mt19937_64 rng(812);
  const auto adjacency = TypeParam::order(32, generate_random_tree(32, rng));
  const auto words = TypeParam::encode(32, adjacency);
  Tree tree(words, 32);
  Context disabled(tree, 0);
  EXPECT_EQ(disabled.budget_usage_bytes(), 0u);
  EXPECT_GE(disabled.memory_usage_bytes(), sizeof(Context));
  EXPECT_EQ(disabled.source_words().data(), tree.encoding_words().data());
  Context context(tree, 32768);
  const auto copy = context;
  auto first = context.cursor();
  auto second = context.cursor();
  auto copied = copy.cursor();
  const auto root = second.node();
  ASSERT_GT(first.degree(), 0u);
  first.move_to_child(0);
  EXPECT_EQ(second.node().number, root.number);
  EXPECT_EQ(second.node().pos, root.pos);
  EXPECT_EQ(copied.node().number, root.number);
  EXPECT_EQ(copied.node().pos, root.pos);
  EXPECT_NE(first.node().number, root.number);
}

TYPED_TEST(BoundedTraversalTest, IndexRejectsEmptyTrees) {
  using Tree = typename TypeParam::Tree;
  const std::vector<uint64_t> words;
  const Tree tree(words, 0);
  for (size_t budget : {size_t{0}, size_t{32768}}) {
    EXPECT_THROW((pixie::RootNavigationIndex<Tree>(tree, budget)),
                 std::invalid_argument);
    EXPECT_THROW((pixie::TreeTraversalContext(tree, budget)),
                 std::invalid_argument);
  }
}

TYPED_TEST(BoundedTraversalTest, CursorContractsAcrossAdmissionBoundaries) {
  using Tree = typename TypeParam::Tree;
  using Index = pixie::RootNavigationIndex<Tree, 8>;
  // Unequal degrees and tied priorities; more nodes than the index can hold.
  const std::vector<std::vector<size_t>> source = {
      {0, 1, 2, 3}, {0, 4, 5}, {0, 6}, {0}, {1, 7, 8}, {1},
      {2, 9, 10},   {4},       {4},    {6}, {6}};
  const auto adjacency = TypeParam::order(source.size(), source);
  const auto words = TypeParam::encode(source.size(), adjacency);
  const Tree tree(words, source.size());
  for (size_t budget : {size_t{0}, sizeof(Index), sizeof(Index) + 1,
                        sizeof(Index) + 32, size_t{32768}}) {
    SCOPED_TRACE(budget);
    const Index index(tree, budget);
    EXPECT_LE(index.budget_usage_bytes(), budget);
    auto copied_index = index;
    auto moved_index = std::move(copied_index);
    std::vector<typename Index::Cursor> pending{moved_index.cursor()};
    size_t visited = 0;
    while (!pending.empty()) {
      auto cursor = pending.back();
      pending.pop_back();
      ++visited;
      const auto node = cursor.node();
      using Handle = std::pair<size_t, size_t>;
      const Handle origin{node.number, node.pos};
      std::vector<Handle> expected, actual;
      tree.for_each_child(node, [&](auto child) {
        expected.emplace_back(child.number, child.pos);
      });
      auto collect = [&, token = std::make_unique<int>(1)](auto child) {
        actual.emplace_back(child.number, child.pos);
      };
      cursor.for_each_child(collect);
      EXPECT_EQ(actual, expected);
      ASSERT_EQ(cursor.degree(), expected.size());
      for (size_t i = 0; i < expected.size(); ++i) {
        SCOPED_TRACE(i);
        // Every exception position must preserve the cursor, on hits and
        // misses.
        size_t calls = 0;
        auto throwing = [&](auto) {
          if (calls++ == i) {
            throw std::runtime_error("stop");
          }
        };
        EXPECT_THROW(cursor.for_each_child(throwing), std::runtime_error);
        EXPECT_EQ(calls, i + 1);
        calls = 0;
        EXPECT_THROW(cursor.inspect_children_and_move_to_child(i, throwing),
                     std::runtime_error);
        EXPECT_EQ(calls, i + 1);
        EXPECT_EQ(Handle(cursor.node().number, cursor.node().pos), origin);
        actual.clear();
        auto inspected = cursor;
        inspected.inspect_children_and_move_to_child(i, collect);
        EXPECT_EQ(actual, expected);
        EXPECT_EQ(Handle(inspected.node().number, inspected.node().pos),
                  expected[i]);
        auto selected = cursor;
        selected.move_to_child(i);
        EXPECT_EQ(Handle(selected.node().number, selected.node().pos),
                  expected[i]);
        auto supplied = cursor;
        supplied.move_to_child(i, tree.child(node, i));
        EXPECT_EQ(Handle(supplied.node().number, supplied.node().pos),
                  expected[i]);
        pending.push_back(supplied);
      }
    }
    EXPECT_EQ(visited, source.size());
  }
}

TEST(TreeUnpackingTest, WordReaderZeroExtendsUnalignedSourceTails) {
  const std::vector<uint64_t> words{0x8123456789abcdef, 0xfedcba9876543210};
  for (size_t count = 0; count <= words.size(); ++count) {
    const std::span<const uint64_t> source(words.data(), count);
    for (size_t position = 0; position <= 192; ++position) {
      uint64_t expected = 0;
      for (size_t bit = 0; bit < 64; ++bit) {
        const size_t at = position + bit;
        if (at < count * 64) {
          expected |= ((words[at / 64] >> (at % 64)) & 1) << bit;
        }
      }
      EXPECT_EQ(pixie::detail::tree_unpack::word_at(source, position), expected)
          << "words=" << count << " position=" << position;
    }
  }
}

TYPED_TEST(BoundedTraversalTest, AllWorkloadsShareBudgetAndExactEvents) {
  using Tree = typename TypeParam::Tree;
  std::mt19937_64 rng(812);
  for (size_t scenario = 0; scenario < 7; ++scenario) {
    const size_t n = scenario == 0 ? 1 : scenario == 1 ? 130 : 10001;
    auto source = generate_random_tree(n, rng);
    if (scenario >= 3) {
      source.assign(n, {});
      source[0].push_back(0);
      for (size_t i = 1; i < n; ++i) {
        size_t parent = scenario == 3 ? i - 1 : 0;
        if (scenario >= 5) {
          // Deep combs force sibling transitions after falling back from
          // bounded decoder stacks, rather than only ascending a pure chain.
          parent = i <= 2 ? 0 : 2 * ((i - 1) / 2) - 1;
        }
        source[parent].push_back(i);
        source[i].push_back(parent);
      }
    }
    if (scenario == 6) {
      for (auto& edges : source) {
        std::reverse(edges.begin() + 1, edges.end());
      }
    }
    const auto adjacency = TypeParam::order(n, source);
    const auto words = TypeParam::encode(n, adjacency);
    for (size_t budget :
         {size_t(0), size_t(64), size_t(4096), size_t(8192), size_t(16384),
          size_t(32767), size_t(32768), size_t(32769), size_t(65536),
          size_t(81920), size_t(196608), size_t(524288)}) {
      Tree tree(words, n);
      pixie::TreeTraversalContext cache(tree, budget);
      pixie::benchmarks::TreeTraversalWorkloads workload(cache);
      EXPECT_LE(cache.budget_usage_bytes(), budget);
      if (budget >= 32768) {
        EXPECT_TRUE(cache.streaming_enabled());
        EXPECT_TRUE(cache.optimized_inspection_enabled());
      }
      std::vector<std::pair<size_t, size_t>> expected, actual;
      auto collect = [&](auto node) {
        actual.emplace_back(node.number, node.pos);
      };
      auto direct = [&](auto node) {
        expected.emplace_back(node.number, node.pos);
      };
      pixie::detail::tree_unpack::subtree(tree, tree.root(), direct);
      workload.traverse(collect);
      EXPECT_EQ(actual, expected);
      actual.clear();
      expected.clear();
      auto ref_events = [&](auto node, bool entering) {
        direct(node);
        if (entering) {
          tree.for_each_child(node, direct);
        }
      };
      pixie::detail::tree_unpack::subtree(tree, tree.root(), ref_events);
      const auto work = workload.inspect(collect, collect);
      EXPECT_EQ(work.first, n);
      EXPECT_EQ(work.second, n - 1);
      EXPECT_EQ(actual, expected);
      for (size_t trial = 0; trial < 32; ++trial) {
        auto cursor = cache.cursor();
        auto node = tree.root();
        while (!tree.is_leaf(node)) {
          ASSERT_EQ(cursor.node().number, node.number);
          ASSERT_EQ(cursor.node().pos, node.pos);
          ASSERT_EQ(cursor.degree(), tree.degree(node));
          actual.clear();
          expected.clear();
          tree.for_each_child(node, direct);
          const size_t choice = rng() % tree.degree(node);
          if (trial % 2) {
            cursor.inspect_children_and_move_to_child(choice, collect);
            EXPECT_EQ(actual, expected);
          } else {
            cursor.move_to_child(choice);
          }
          node = tree.child(node, choice);
        }
        EXPECT_EQ(cursor.node().number, node.number);
        EXPECT_EQ(cursor.node().pos, node.pos);
        EXPECT_EQ(cursor.degree(), 0);
      }
    }
  }
}

TYPED_TEST(TreeSpecificationTest,
           ForEachChildBorrowsVisitorAndSupportsNesting) {
  using Tree = typename TypeParam::Tree;
  const std::vector<std::vector<size_t>> source = {
      {0, 1, 2}, {0, 3, 4}, {0}, {1}, {1}};
  const auto adjacency = TypeParam::order(source.size(), source);
  const auto words = TypeParam::encode(source.size(), adjacency);
  const Tree tree(words, source.size());
  const pixie::TreeBase<Tree>& facade = tree;
  const AdjListTree reference(adjacency);
  struct Visitor {
    Visitor() = default;
    Visitor(const Visitor&) = delete;
    std::vector<pixie::TreeNode> nodes;
    void operator()(pixie::TreeNode node) { nodes.push_back(node); }
  } visitor;
  facade.for_each_child(facade.root(), visitor);
  ASSERT_EQ(visitor.nodes.size(), 2u);
  size_t index = 0;
  reference.for_each_child(reference.root(), [&](auto expected) {
    const auto actual = visitor.nodes[index];
    EXPECT_EQ(actual.number, expected.number);
    EXPECT_EQ(actual.pos, tree.child(tree.root(), index).pos);
    ++index;
  });
  size_t nested_count = 0;
  facade.for_each_child(facade.root(), [&](auto child) {
    facade.for_each_child(child, [&](auto grandchild) {
      ++nested_count;
      facade.for_each_child(grandchild, [](auto) { ADD_FAILURE(); });
    });
  });
  EXPECT_EQ(nested_count, 2u);
  size_t calls = 0;
  EXPECT_THROW(facade.for_each_child(facade.root(),
                                     [&](auto) {
                                       ++calls;
                                       throw std::runtime_error("stop");
                                     }),
               std::runtime_error);
  EXPECT_EQ(calls, 1u);
  visitor.nodes.clear();
  facade.for_each_child(facade.root(), visitor);
  EXPECT_EQ(visitor.nodes.size(), 2u);
}

TYPED_TEST(TreeSpecificationTest,
           IndexedPathsMatchNavigationWithBoundedStorage) {
  using Tree = typename TypeParam::Tree;
  for (int shape = 0; shape < 3; ++shape) {
    std::mt19937_64 rng(42);
    constexpr size_t n = 521;
    auto source = generate_random_tree(n, rng);
    if (shape) {
      for (auto& edges : source) {
        edges.clear();
      }
      source[0].push_back(0);
      for (size_t i = 1; i < n; ++i) {
        const size_t parent = shape == 1 ? i - 1 : 0;
        source[i].push_back(parent);
        source[parent].push_back(i);
      }
    }
    const auto adjacency = TypeParam::order(n, source);
    const auto words = TypeParam::encode(n, adjacency);
    const Tree tree(words, n);
    // Deliberately tiny: chain and star force fallback; random trees mix paths.
    const pixie::RootNavigationIndex<Tree, 8> index(tree);
    EXPECT_LE(index.memory_usage_bytes(), sizeof(index) + 8 * 32u);
    auto copy = index.cursor();
    auto original = copy;
    original.move_to_child(0);
    EXPECT_EQ(copy.node().number, tree.root().number);
    EXPECT_EQ(original.node().pos, tree.child(tree.root(), 0).pos);
    size_t child_index = 0;
    copy.for_each_child([&](auto child) {
      EXPECT_EQ(child.pos, tree.child(tree.root(), child_index++).pos);
    });
    EXPECT_EQ(child_index, tree.degree(tree.root()));
    EXPECT_THROW(copy.inspect_children_and_move_to_child(
                     0, [](auto) { throw std::runtime_error("stop"); }),
                 std::runtime_error);
    EXPECT_EQ(copy.node().pos, tree.root().pos);
    while (copy.degree()) {
      copy.move_to_child(0);
    }
    copy.for_each_child([](auto) { ADD_FAILURE(); });
    pixie::benchmarks::ShiftedQueryPool a(0, SIZE_MAX, 43), b(0, SIZE_MAX, 43),
        c(0, SIZE_MAX, 43);
    for (size_t trial = 0; trial < 128; ++trial) {
      std::vector<size_t> expected, actual, inspected, expected_children,
          actual_children;
      const auto plain = pixie::benchmarks::RandomRootToLeaf<true>(
          tree, a, [&](auto node) { expected.push_back(node.pos); },
          [&](auto node) { expected_children.push_back(node.pos); });
      const auto indexed = pixie::benchmarks::IndexedRandomRootToLeaf<false>(
          index, b, [&](auto node) { actual.push_back(node.pos); },
          [](auto) { ADD_FAILURE(); });
      const auto indexed_inspect =
          pixie::benchmarks::IndexedRandomRootToLeaf<true>(
              index, c, [&](auto node) { inspected.push_back(node.pos); },
              [&](auto node) { actual_children.push_back(node.pos); });
      EXPECT_EQ(actual, expected);
      EXPECT_EQ(inspected, expected);
      EXPECT_EQ(actual_children, expected_children);
      EXPECT_EQ(indexed.nodes, plain.nodes);
      EXPECT_EQ(indexed_inspect.children, plain.children);
    }
  }
}

TYPED_TEST(TreeSpecificationTest, BasicNavigation) {
  using Tree = typename TypeParam::Tree;
  constexpr std::size_t tree_size = 5;
  const std::vector<std::vector<std::size_t>> source = {
      {0, 1}, {0, 2}, {1, 3}, {2, 4}, {3}};
  const auto adjacency = TypeParam::order(tree_size, source);
  const auto words = TypeParam::encode(tree_size, adjacency);
  const Tree tree(words, tree_size);
  const AdjListTree reference(adjacency);

  EXPECT_EQ(tree.size(), tree_size);
  EXPECT_FALSE(tree.empty());
  auto node = tree.root();
  auto expected = reference.root();
  EXPECT_TRUE(tree.is_root(node));
  for (std::size_t i = 0; i < tree_size; ++i) {
    EXPECT_EQ(node.number, expected.number);
    EXPECT_EQ(tree.degree(node), reference.degree(expected));
    EXPECT_EQ(tree.is_leaf(node), reference.is_leaf(expected));
    if (tree.is_leaf(node)) {
      break;
    }
    node = tree.child(node, 0);
    expected = reference.child(expected, 0);
  }
}

TYPED_TEST(TreeSpecificationTest, RandomTreesSupportDepthFirstNavigation) {
  using Tree = typename TypeParam::Tree;
  for (std::size_t tree_size = 8; tree_size < (1 << 18); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    auto adjacency =
        TypeParam::order(tree_size, generate_random_tree(tree_size, rng));
    const auto words = TypeParam::encode(tree_size, adjacency);
    const Tree tree(words, tree_size);
    const AdjListTree reference(adjacency);

    using NodePair = std::pair<typename Tree::Node, AdjListNode>;
    std::stack<NodePair> pending;
    pending.push({tree.root(), reference.root()});
    std::size_t visited = 0;

    while (!pending.empty()) {
      const auto [node, expected] = pending.top();
      pending.pop();
      ++visited;
      EXPECT_EQ(node.number, expected.number);
      EXPECT_EQ(tree.parent(node).number, reference.parent(expected).number);
      EXPECT_EQ(tree.is_leaf(node), reference.is_leaf(expected));
      EXPECT_EQ(tree.is_root(node), reference.is_root(expected));
      if (!tree.is_root(node)) {
        EXPECT_EQ(tree.is_last_child(node), reference.is_last_child(expected));
        if (!tree.is_last_child(node)) {
          EXPECT_EQ(tree.next_sibling(node).number,
                    reference.next_sibling(expected).number);
        }
      }

      const std::size_t degree = tree.degree(node);
      ASSERT_EQ(degree, reference.degree(expected));
      for (std::size_t i = 0; i < degree; ++i) {
        pending.push({tree.child(node, i), reference.child(expected, i)});
      }
      if (degree != 0) {
        EXPECT_EQ(tree.first_child(node).number,
                  reference.first_child(expected).number);
      }
    }
    EXPECT_EQ(visited, tree_size);
  }
}

TYPED_TEST(TreeSpecificationTest, RandomTreesSupportBreadthFirstNavigation) {
  using Tree = typename TypeParam::Tree;
  for (std::size_t tree_size = 8; tree_size < (1 << 18); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    auto adjacency =
        TypeParam::order(tree_size, generate_random_tree(tree_size, rng));
    const auto words = TypeParam::encode(tree_size, adjacency);
    const Tree tree(words, tree_size);
    const AdjListTree reference(adjacency);

    using NodePair = std::pair<typename Tree::Node, AdjListNode>;
    std::queue<NodePair> pending;
    pending.push({tree.root(), reference.root()});

    while (!pending.empty()) {
      const auto [node, expected] = pending.front();
      pending.pop();
      const std::size_t degree = tree.degree(node);
      ASSERT_EQ(degree, reference.degree(expected));
      for (std::size_t i = 0; i < degree; ++i) {
        pending.push({tree.child(node, i), reference.child(expected, i)});
      }
    }
  }
}

TYPED_TEST(TreeSpecificationTest, TraversesWithoutAuxiliaryNodeStorage) {
  using Tree = typename TypeParam::Tree;
  constexpr std::size_t tree_size = 4096;
  std::mt19937_64 rng(42);
  auto adjacency =
      TypeParam::order(tree_size, generate_random_tree(tree_size, rng));
  const auto words = TypeParam::encode(tree_size, adjacency);
  const Tree tree(words, tree_size);

  auto node = tree.root();
  bool descending = true;
  std::size_t visited = 1;
  while (true) {
    if (descending && !tree.is_leaf(node)) {
      node = tree.first_child(node);
      ++visited;
      continue;
    }
    descending = false;
    if (tree.is_root(node)) {
      break;
    }
    if (tree.is_last_child(node)) {
      node = tree.parent(node);
    } else {
      node = tree.next_sibling(node);
      descending = true;
      ++visited;
    }
  }
  EXPECT_EQ(visited, tree_size);
}

TYPED_TEST(TreeSpecificationTest, NavigationAcrossLongRunsAndWordBoundaries) {
  using Tree = typename TypeParam::Tree;
  for (const std::size_t size : {1, 2, 63, 64, 65, 129, 513}) {
    for (int shape = 0; shape < 3; ++shape) {
      SCOPED_TRACE(testing::Message() << "size=" << size << " shape=" << shape);
      std::vector<std::vector<std::size_t>> source(size);
      source[0].push_back(0);
      for (std::size_t i = 1; i < size; ++i) {
        // Chain, star, and a chain ending in a high-degree node.
        const auto parent = shape == 0   ? i - 1
                            : shape == 1 ? 0
                                         : std::min(i - 1, size / 2);
        source[i].push_back(parent);
        source[parent].push_back(i);
      }
      const auto adjacency = TypeParam::order(size, source);
      const auto words = TypeParam::encode(size, adjacency);
      const Tree tree(words, size);
      const AdjListTree reference(adjacency);
      std::queue<std::pair<typename Tree::Node, AdjListNode>> pending;
      pending.push({tree.root(), reference.root()});
      std::size_t visited = 0;
      EXPECT_EQ(tree.parent(tree.root()).pos, tree.root().pos);
      while (!pending.empty()) {
        const auto [node, expected] = pending.front();
        pending.pop();
        ++visited;
        ASSERT_EQ(node.number, expected.number);
        ASSERT_EQ(tree.degree(node), reference.degree(expected));
        EXPECT_EQ(tree.is_leaf(node), reference.is_leaf(expected));
        for (std::size_t i = 0; i < tree.degree(node); ++i) {
          const auto child = tree.child(node, i);
          EXPECT_EQ(tree.parent(child).number, node.number);
          EXPECT_EQ(tree.parent(child).pos, node.pos);
          EXPECT_EQ(tree.is_last_child(child), i + 1 == tree.degree(node));
          if (i == 0) {
            EXPECT_EQ(tree.first_child(node).pos, child.pos);
          }
          if (i + 1 < tree.degree(node)) {
            const auto sibling = tree.next_sibling(child);
            const auto indexed = tree.child(node, i + 1);
            EXPECT_EQ(sibling.number, indexed.number);
            EXPECT_EQ(sibling.pos, indexed.pos);
          }
          pending.push({child, reference.child(expected, i)});
        }
      }
      EXPECT_EQ(visited, size);
    }
  }
}

template <class Backend>
struct CountingNavigationRankSelect : Backend {
  using Backend::Backend;
  inline static std::size_t zero_selects = 0, one_selects = 0;
  std::size_t select(std::size_t rank) const {
    ++one_selects;
    return Backend::select(rank);
  }
  std::size_t select0(std::size_t rank) const {
    ++zero_selects;
    return Backend::select0(rank);
  }
};

template <class Backend>
class LoudsNavigationViewTest : public testing::Test {};
using NavigationBackends = testing::Types<pixie::RankSelectSupport<>,
                                          pixie::RankSelectSupportCompressed<>>;
TYPED_TEST_SUITE(LoudsNavigationViewTest, NavigationBackends);

TYPED_TEST(LoudsNavigationViewTest, TraversalReusesSelectAndInvalidatesOnMove) {
  using Backend = CountingNavigationRankSelect<TypeParam>;
  using Tree = pixie::LoudsTree<Backend>;
  for (const std::size_t size : {1, 2, 63, 64, 65, 513, 4096}) {
    for (int shape = 0; shape < 4; ++shape) {
      SCOPED_TRACE(testing::Message() << "size=" << size << " shape=" << shape);
      std::mt19937_64 rng(42);
      auto source = generate_random_tree(size, rng);
      if (shape != 0) {
        source.assign(size, {});
        source[0].push_back(0);
        for (std::size_t i = 1; i < size; ++i) {
          const auto parent = shape == 1   ? i - 1
                              : shape == 2 ? 0
                                           : std::min(i - 1, size / 2);
          source[i].push_back(parent);
          source[parent].push_back(i);
        }
      }
      const auto adjacency = bfs_order(size, source);
      const auto words = adj_to_louds(size, adjacency);
      const Tree tree(words, size);
      const AdjListTree reference(adjacency);
      auto view = tree.navigate(tree.root());
      auto expected = reference.root();
      Backend::zero_selects = 0;
      EXPECT_EQ(view.parent().number, 0u);
      EXPECT_EQ(Backend::zero_selects, 0u);
      bool descending = true;
      std::size_t visited = 1, checks = 0, steps = 0;
      while (true) {
        ASSERT_LT(++steps, 3 * size);
        ASSERT_EQ(view.node().number, expected.number);
        ASSERT_EQ(view.degree(), reference.degree(expected));
        ASSERT_EQ(view.is_leaf(), reference.is_leaf(expected));
        if (descending && !view.is_leaf()) {
          EXPECT_EQ(view.child(0).pos, view.first_child().pos);
          view.reset(view.first_child());
          expected = reference.first_child(expected);
          ++visited;
          continue;
        }
        descending = false;
        if (view.is_root()) {
          break;
        }
        const auto last = view.is_last_child();
        ++checks;
        ASSERT_EQ(Backend::zero_selects, checks);
        EXPECT_EQ(view.is_last_child(), last);
        EXPECT_EQ(last, reference.is_last_child(expected));
        // Parent queries reuse the cached edge even on non-last siblings.
        EXPECT_EQ(view.parent().number, reference.parent(expected).number);
        ASSERT_EQ(Backend::zero_selects, checks);
        if (last) {
          view.reset(view.parent());
          expected = reference.parent(expected);
        } else {
          view.reset(view.next_sibling());
          expected = reference.next_sibling(expected);
          descending = true;
          ++visited;
        }
      }
      EXPECT_EQ(visited, size);
      EXPECT_EQ(checks, size - 1);
    }
  }
}

TYPED_TEST(LoudsNavigationViewTest, ViewsAndCopiesHaveIndependentCaches) {
  using Backend = CountingNavigationRankSelect<TypeParam>;
  const std::vector<std::vector<std::size_t>> adjacency = {{0, 1, 2}, {0}, {0}};
  const auto words = adj_to_louds(3, adjacency);
  const pixie::LoudsTree<Backend> tree(words, 3);
  const auto first = tree.first_child(tree.root());
  auto a = tree.navigate(first);
  auto b = tree.navigate(first);
  Backend::zero_selects = 0;
  EXPECT_FALSE(a.is_last_child());
  auto copy = a;
  EXPECT_EQ(copy.parent().number, 0u);
  EXPECT_EQ(Backend::zero_selects, 1u);
  EXPECT_FALSE(b.is_last_child());
  EXPECT_EQ(Backend::zero_selects, 2u);
  a.reset(a.next_sibling());
  EXPECT_TRUE(a.is_last_child());
  EXPECT_EQ(Backend::zero_selects, 3u);
  EXPECT_FALSE(copy.is_last_child());
  EXPECT_FALSE(b.is_last_child());
  EXPECT_EQ(Backend::zero_selects, 3u);
  copy.reset(first);
  EXPECT_EQ(copy.parent().number,
            0u);  // Also supports parent before predicate.
  EXPECT_FALSE(copy.is_last_child());
  EXPECT_EQ(Backend::zero_selects, 4u);
}

TYPED_TEST(TreeSpecificationTest, NavigationViewsMatchDirectQueries) {
  using Tree = typename TypeParam::Tree;
  for (const std::size_t size : {1, 2, 63, 64, 65, 513, 4096}) {
    for (int shape = 0; shape < 3; ++shape) {
      SCOPED_TRACE(testing::Message() << "size=" << size << " shape=" << shape);
      std::mt19937_64 rng(42);
      auto source = generate_random_tree(size, rng);
      if (shape != 0) {
        source.assign(size, {});
        source[0].push_back(0);
        for (std::size_t i = 1; i < size; ++i) {
          const auto parent = shape == 1 ? i - 1 : 0;
          source[i].push_back(parent);
          source[parent].push_back(i);
        }
      }
      const auto adjacency = TypeParam::order(size, source);
      const auto words = TypeParam::encode(size, adjacency);
      const Tree tree(words, size);
      auto view = tree.navigate(tree.root());
      std::queue<typename Tree::Node> pending;
      pending.push(tree.root());
      std::size_t visited = 0;
      while (!pending.empty()) {
        const auto node = pending.front();
        pending.pop();
        view.reset(node);
        ++visited;
        for (int repeat = 0; repeat < 2; ++repeat) {
          EXPECT_EQ(view.node().number, node.number);
          EXPECT_EQ(view.is_root(), tree.is_root(node));
          EXPECT_EQ(view.is_leaf(), tree.is_leaf(node));
          EXPECT_EQ(view.parent().number, tree.parent(node).number);
          EXPECT_EQ(view.parent().pos, tree.parent(node).pos);
          EXPECT_EQ(view.degree(), tree.degree(node));
          if (!view.is_root()) {
            EXPECT_EQ(view.is_last_child(), tree.is_last_child(node));
            if (!view.is_last_child()) {
              EXPECT_EQ(view.next_sibling().number,
                        tree.next_sibling(node).number);
              EXPECT_EQ(view.next_sibling().pos, tree.next_sibling(node).pos);
            }
          }
          if (!view.is_leaf()) {
            EXPECT_EQ(view.first_child().number, tree.first_child(node).number);
            EXPECT_EQ(view.first_child().pos, tree.first_child(node).pos);
          }
        }
        std::size_t child_index = 0;
        view.for_each_child([&](pixie::TreeNode child) {
          const auto expected = tree.child(node, child_index++);
          EXPECT_EQ(child.number, expected.number);
          EXPECT_EQ(child.pos, expected.pos);
        });
        EXPECT_EQ(child_index, view.degree());
        auto copy = view;
        auto independent = tree.navigate(node);
        view.reset(tree.root());
        EXPECT_EQ(copy.node().pos, node.pos);
        EXPECT_EQ(copy.parent().pos, independent.parent().pos);
        for (std::size_t i = 0; i < copy.degree(); ++i) {
          const auto child = copy.child(i);
          EXPECT_EQ(child.number, tree.child(node, i).number);
          EXPECT_EQ(child.pos, tree.child(node, i).pos);
          pending.push(child);
        }
      }
      EXPECT_EQ(visited, size);
      // Exercise invalidation during an actual stackless DFS as well.
      view.reset(tree.root());
      bool descending = true;
      visited = 1;
      std::size_t steps = 0;
      while (true) {
        ASSERT_LT(++steps, 3 * size);
        if (descending && !view.is_leaf()) {
          view.reset(view.first_child());
          ++visited;
          continue;
        }
        descending = false;
        if (view.is_root()) {
          break;
        }
        if (view.is_last_child()) {
          view.reset(view.parent());
        } else {
          view.reset(view.next_sibling());
          descending = true;
          ++visited;
        }
      }
      EXPECT_EQ(visited, size);
    }
  }
}

template <class Backend>
struct CountingNavigationRmM : Backend {
  using Backend::Backend;
  inline static std::size_t forward_calls = 0, open_calls = 0, close_calls = 0;
  inline static std::size_t rank_calls = 0, select0_calls = 0;
  size_t select0(size_t rank) const {
    ++select0_calls;
    return Backend::select0(rank);
  }
  size_t rank0(size_t pos) const {
    ++rank_calls;
    return Backend::rank0(pos);
  }
  size_t rank1(size_t pos) const {
    ++rank_calls;
    return Backend::rank1(pos);
  }

  size_t fwdsearch(size_t pos, int delta) const {
    ++forward_calls;
    return Backend::fwdsearch(pos, delta);
  }
  size_t open(size_t pos) const {
    ++open_calls;
    return Backend::open(pos);
  }
  size_t close(size_t pos) const {
    ++close_calls;
    return Backend::close(pos);
  }
};
template <class Backend>
class ParenthesisNavigationViewTest : public testing::Test {};
using NavigationRmMBackends =
    testing::Types<pixie::RmMBTree<>,
                   pixie::RmMBTree<4, 32, pixie::RankSelectSupport<>>>;
TYPED_TEST_SUITE(ParenthesisNavigationViewTest, NavigationRmMBackends);

TYPED_TEST(ParenthesisNavigationViewTest,
           DfudsInspectionReusesTheDegreeDelimiter) {
  using Backend = CountingNavigationRmM<TypeParam>;
  constexpr size_t size = 66;
  std::vector<std::vector<size_t>> adjacency(size);
  adjacency[0].push_back(0);
  for (size_t i = 1; i < size; ++i) {
    adjacency[0].push_back(i);
    adjacency[i].push_back(0);
  }
  const auto words = pixie::adj_to_dfuds(size, adjacency);
  const auto check = [&]() {
    const pixie::DFUDSTree<Backend> tree(words, size);
    auto view = tree.navigate(tree.root());
    Backend::select0_calls = 0;
    EXPECT_EQ(view.degree(), size - 1);
    for (int pass = 0; pass < 2; ++pass) {
      size_t count = 0;
      view.for_each_child([&](auto child) {
        EXPECT_EQ(child.number, ++count);
        EXPECT_TRUE(tree.is_leaf(child));
      });
      EXPECT_EQ(count, size - 1);
    }
    EXPECT_EQ(Backend::select0_calls, 1u);
  };
  check();
}

TYPED_TEST(ParenthesisNavigationViewTest, DfsReusesParenthesisSearches) {
  using Backend = CountingNavigationRmM<TypeParam>;
  constexpr std::size_t size = 4096;
  std::mt19937_64 rng(42);
  const auto adjacency = dfs_order(size, generate_random_tree(size, rng));
  std::size_t nonleaf_nonroot = 0, nonleaf_nonlast = 0;
  std::size_t nonleaf = adjacency[0].size() > 1;
  for (std::size_t node = 1; node < size; ++node) {
    if (adjacency[node].size() > 1) {
      ++nonleaf_nonroot;
      ++nonleaf;
      nonleaf_nonlast += adjacency[adjacency[node][0]].back() != node;
    }
  }
  const auto check = [&]<class Tree>(const Tree& tree, bool dfuds) {
    Backend::forward_calls = Backend::open_calls = Backend::close_calls = 0;
    Backend::rank_calls = 0;
    auto view = tree.navigate(tree.root());
    bool descending = true;
    std::size_t visited = 1, steps = 0;
    while (true) {
      ASSERT_LT(++steps, 3 * size);
      if (descending && !view.is_leaf()) {
        view.reset(view.first_child());
        ++visited;
        continue;
      }
      descending = false;
      if (view.is_root()) {
        break;
      }
      const bool last = view.is_last_child();
      EXPECT_EQ(view.is_last_child(), last);
      if (last) {
        view.reset(view.parent());
      } else {
        view.reset(view.next_sibling());
        descending = true;
        ++visited;
      }
    }
    EXPECT_EQ(visited, size);
    EXPECT_EQ(Backend::forward_calls, dfuds ? nonleaf_nonlast : 0);
    EXPECT_EQ(Backend::open_calls, dfuds ? size - 1 - nonleaf : 0);
    EXPECT_EQ(Backend::close_calls, dfuds ? 0 : nonleaf_nonroot);
    EXPECT_EQ(Backend::rank_calls, 0u);
  };
  const auto dfuds = pixie::adj_to_dfuds(size, adjacency);
  check(pixie::DFUDSTree<Backend>(dfuds, size), true);
  const auto bp = pixie::adj_to_bp(size, adjacency);
  check(pixie::BPTree<Backend>(bp, size), false);
}

TYPED_TEST(LoudsNavigationViewTest, DegreeSiblingAndFirstChildShareSelects) {
  using Backend = CountingNavigationRankSelect<TypeParam>;
  // Node 1's long child run forces the next-delimiter select fallback.
  constexpr std::size_t size = 133;
  std::vector<std::vector<std::size_t>> adjacency(size);
  adjacency[0] = {0, 1, 2};
  adjacency[1] = {0};
  adjacency[2] = {0};
  for (std::size_t i = 3; i < size; ++i) {
    adjacency[1].push_back(i);
    adjacency[i] = {1};
  }
  const auto words = adj_to_louds(size, adjacency);
  const pixie::LoudsTree<Backend> tree(words, size);
  auto view = tree.navigate(tree.first_child(tree.root()));
  Backend::one_selects = Backend::zero_selects = 0;
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::one_selects, 1u);
  EXPECT_EQ(view.next_sibling().number, 2u);
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::one_selects, 1u);
  EXPECT_EQ(view.first_child().number, 3u);
  EXPECT_EQ(view.child(0).number, 3u);
  EXPECT_EQ(view.first_child().number, 3u);
  EXPECT_EQ(Backend::one_selects, 2u);
  view.reset(view.node());
  EXPECT_EQ(view.next_sibling().number, 2u);
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::one_selects, 3u);
}

TYPED_TEST(TreeSpecificationTest, ExhaustiveSmallOrderedTreeNavigation) {
  using Tree = typename TypeParam::Tree;
  std::size_t shapes = 0;
  for (std::size_t size = 1; size <= 9; ++size) {
    const unsigned inner_bits = 2 * (size - 1);
    for (unsigned mask = 0; mask < (1u << inner_bits); ++mask) {
      if (std::popcount(mask) != static_cast<int>(size - 1)) {
        continue;
      }
      std::vector<std::vector<std::size_t>> source(1, {0});
      std::vector<std::size_t> stack = {0};
      bool valid = true;
      for (unsigned bit = 0; bit < inner_bits; ++bit) {
        if ((mask >> bit) & 1u) {
          const auto node = source.size();
          source[stack.back()].push_back(node);
          source.push_back({stack.back()});
          stack.push_back(node);
        } else {
          if (stack.size() == 1) {
            valid = false;
            break;
          }
          stack.pop_back();
        }
      }
      if (!valid) {
        continue;
      }
      ASSERT_EQ(stack.size(), 1u);
      ASSERT_EQ(source.size(), size);
      ++shapes;
      SCOPED_TRACE(testing::Message() << "size=" << size << " mask=" << mask);
      const auto adjacency = TypeParam::order(size, source);
      const auto words = TypeParam::encode(size, adjacency);
      const Tree tree(words, size);
      const AdjListTree reference(adjacency);
      std::queue<std::pair<typename Tree::Node, AdjListNode>> pending;
      pending.push({tree.root(), reference.root()});
      while (!pending.empty()) {
        const auto [node, expected] = pending.front();
        pending.pop();
        auto view = tree.navigate(node);
        ASSERT_EQ(node.number, expected.number);
        ASSERT_EQ(tree.degree(node), reference.degree(expected));
        ASSERT_EQ(view.degree(), reference.degree(expected));
        EXPECT_EQ(tree.is_leaf(node), reference.is_leaf(expected));
        if (!tree.is_root(node)) {
          EXPECT_EQ(tree.is_last_child(node),
                    reference.is_last_child(expected));
          EXPECT_EQ(view.is_last_child(), reference.is_last_child(expected));
          if (!reference.is_last_child(expected)) {
            EXPECT_EQ(tree.next_sibling(node).number,
                      reference.next_sibling(expected).number);
            EXPECT_EQ(view.next_sibling().number,
                      reference.next_sibling(expected).number);
          }
        }
        EXPECT_EQ(tree.parent(node).number, reference.parent(expected).number);
        EXPECT_EQ(view.parent().number, reference.parent(expected).number);
        for (std::size_t child = 0; child < reference.degree(expected);
             ++child) {
          const auto actual = tree.child(node, child);
          EXPECT_EQ(tree.parent(actual).pos, node.pos);
          EXPECT_EQ(tree.navigate(actual).parent().pos, node.pos);
          pending.push({actual, reference.child(expected, child)});
        }
      }
    }
  }
  EXPECT_EQ(shapes, 2056u);
}

template <class Tree>
struct ParentCountingTree {
  const Tree& tree;
  mutable std::size_t parent_calls = 0;
  struct NodeView {
    typename Tree::NodeView view;
    std::size_t* calls;
    auto node() const { return view.node(); }
    void reset(pixie::TreeNode node) { view.reset(node); }
    bool is_root() const { return view.is_root(); }
    bool is_leaf() const { return view.is_leaf(); }
    bool is_last_child() const { return view.is_last_child(); }
    auto degree() const { return view.degree(); }
    auto first_child() const { return view.first_child(); }
    auto child(std::size_t i) const { return view.child(i); }
    auto next_sibling() const { return view.next_sibling(); }
    auto parent() const {
      ++*calls;
      return view.parent();
    }
  };
  auto root() const { return tree.root(); }
  NodeView navigate(pixie::TreeNode node) const {
    return {tree.navigate(node), &parent_calls};
  }
};

TYPED_TEST(TreeSpecificationTest, BoundedPathCacheFallsBackAndResumes) {
  using Tree = typename TypeParam::Tree;
  for (int shape = 0; shape < 4; ++shape) {
    const std::size_t size = shape == 0 ? 1 : 521;
    std::mt19937_64 rng(42);
    auto source = generate_random_tree(size, rng);
    if (shape == 1 || shape == 2) {
      source.assign(size, {});
      source[0].push_back(0);
      for (std::size_t i = 1; i < size; ++i) {
        const auto parent = i == 1 || (shape == 2 && i == 261) ? 0 : i - 1;
        source[i].push_back(parent);
        source[parent].push_back(i);
      }
    }
    const auto adjacency = TypeParam::order(size, source);
    const auto words = TypeParam::encode(size, adjacency);
    const Tree tree(words, size);
    const AdjListTree reference(adjacency);
    ParentCountingTree<Tree> counted{tree};
    const auto check = [&]<std::size_t Capacity>() {
      SCOPED_TRACE(testing::Message()
                   << "shape=" << shape << " capacity=" << Capacity);
      counted.parent_calls = 0;
      pixie::TreePathView<ParentCountingTree<Tree>, Capacity> path(counted);
      auto expected = reference.root();
      path.move_to_parent();
      EXPECT_EQ(counted.parent_calls, 0u);
      std::size_t depth = 0, visited = 1, fallbacks = 0, steps = 0;
      bool descending = true;
      while (true) {
        ASSERT_LT(++steps, 3 * size);
        ASSERT_EQ(path.node().number, expected.number);
        EXPECT_EQ(path.depth(), depth);
        EXPECT_EQ(path.is_root(), reference.is_root(expected));
        EXPECT_EQ(path.is_leaf(), reference.is_leaf(expected));
        EXPECT_EQ(path.degree(), reference.degree(expected));
        if (descending && !path.is_leaf()) {
          // Exercise both movement entry points.
          if (depth % 2) {
            path.move_to_child(0);
          } else {
            path.move_to_first_child();
          }
          expected = reference.first_child(expected);
          ++depth;
          ++visited;
          continue;
        }
        descending = false;
        if (path.is_root()) {
          break;
        }
        ASSERT_EQ(path.is_last_child(), reference.is_last_child(expected));
        if (path.is_last_child()) {
          --depth;
          fallbacks += depth >= Capacity;
          path.move_to_parent();
          expected = reference.parent(expected);
        } else {
          path.move_to_next_sibling();
          expected = reference.next_sibling(expected);
          descending = true;
          ++visited;
        }
        EXPECT_EQ(counted.parent_calls, fallbacks);
      }
      EXPECT_EQ(visited, size);
      if (shape == 1) {
        EXPECT_EQ(fallbacks, size - 1 > Capacity ? size - 1 - Capacity : 0);
      }
      if (shape == 2) {
        // The first branch reaches its sibling without ascending to root.
        EXPECT_EQ(fallbacks,
                  2 * (260 > Capacity ? 260 - Capacity : 0) - (Capacity == 0));
      }
    };
    check.template operator()<0>();
    check.template operator()<1>();
    check.template operator()<2>();
    check.template operator()<128>();
  }
}

TYPED_TEST(TreeSpecificationTest, BoundedPathCopiesAndIndexedDescent) {
  using Tree = typename TypeParam::Tree;
  const std::vector<std::vector<std::size_t>> source = {
      {0, 1, 2}, {0, 3}, {0, 4}, {1}, {2}};
  const auto adjacency = TypeParam::order(5, source);
  const auto words = TypeParam::encode(5, adjacency);
  const Tree tree(words, 5);
  pixie::TreePathView<Tree, 1> path(tree);
  path.move_to_child(1);
  auto copy = path;
  const auto expected_parent = path.node();
  path.move_to_first_child();
  EXPECT_EQ(copy.depth(), 1u);
  EXPECT_EQ(copy.node().number, expected_parent.number);
  copy.move_to_parent();
  EXPECT_TRUE(copy.is_root());
  path.move_to_parent();  // Fallback outside the single cached level.
  EXPECT_EQ(path.node().number, expected_parent.number);
  EXPECT_EQ(path.node().pos, expected_parent.pos);
  path.move_to_parent();  // Cache is usable again.
  EXPECT_TRUE(path.is_root());
}

// Compare every callback (including internal-node returns) against navigation
// and the independent adjacency topology, with exact-sized, unpadded sources.
template <class Case>
void CheckDepthFirst(const std::vector<std::vector<size_t>>& source) {
  using Tree = typename Case::Tree;
  const size_t n = source.size();
  const auto adjacency = Case::order(n, source);
  auto words = Case::encode(n, adjacency);
  const size_t bits = Case::bit_size(n);
  // Deliberately dirty unused tail bits: neither decoder may consume them.
  if (bits % 64) {
    words.back() |= ~uint64_t{0} << (bits % 64);
  }
  words.shrink_to_fit();
  const Tree tree(words, n);
  const AdjListTree reference(adjacency);
  std::vector<std::pair<size_t, size_t>> expected;
  std::vector<size_t> positions(n), child_numbers;
  auto node = tree.root();
  auto ref = reference.root();
  bool down = true;
  while (true) {
    ASSERT_EQ(node.number, ref.number);
    expected.emplace_back(node.number, node.pos);
    if (down) {
      positions[node.number] = node.pos;
      const auto& edges = adjacency[node.number];
      child_numbers.insert(child_numbers.end(), edges.begin() + 1, edges.end());
    }
    if (down && !reference.is_leaf(ref)) {
      node = tree.first_child(node);
      ref = reference.first_child(ref);
      continue;
    }
    down = false;
    if (reference.is_root(ref)) {
      break;
    }
    if (reference.is_last_child(ref)) {
      node = tree.parent(node);
      ref = reference.parent(ref);
    } else {
      node = tree.next_sibling(node);
      ref = reference.next_sibling(ref);
      down = true;
    }
  }
  for (bool unpack : {false, true}) {
    std::vector<std::pair<size_t, size_t>> actual;
    // Move-only visitors must be borrowed, not copied.
    auto visit = [&,
                  token = std::make_unique<int>(1)](pixie::TreeNode current) {
      actual.emplace_back(current.number, current.pos);
    };
    if (unpack) {
      pixie::benchmarks::DepthFirst<true>(tree, words, visit);
    } else {
      pixie::benchmarks::DepthFirst(tree, words, visit);
    }
    ASSERT_EQ(actual, expected) << "unpack=" << unpack;
  }
  std::vector<std::pair<size_t, size_t>> inspected_events, inspected_children,
      expected_children;
  for (const auto child : child_numbers) {
    expected_children.emplace_back(child, positions[child]);
  }
  const auto work = pixie::benchmarks::UnpackedDepthFirstInspect(
      tree, words,
      [&](auto current) {
        inspected_events.emplace_back(current.number, current.pos);
      },
      [&](auto child) {
        inspected_children.emplace_back(child.number, child.pos);
      });
  ASSERT_EQ(inspected_events, expected);
  ASSERT_EQ(inspected_children, expected_children);
  EXPECT_EQ(work.nodes, n);
  EXPECT_EQ(work.children, n - 1);
}

TYPED_TEST(TreeSpecificationTest, DepthFirstDecodingBoundariesAndShapes) {
  for (size_t n : {1, 2, 3, 7, 31, 32, 33, 63, 64, 65, 127, 128, 129, 130, 257,
                   521, 4097}) {
    for (int shape = 0; shape < 7; ++shape) {
      SCOPED_TRACE(testing::Message() << "n=" << n << " shape=" << shape);
      std::mt19937_64 rng(7823 + n);
      auto source = generate_random_tree(n, rng);
      if (shape) {
        for (auto& v : source) {
          v.clear();
        }
        source[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          size_t parent = shape == 1   ? i - 1
                          : shape == 2 ? 0
                          : shape == 3 ? (i - 1) / 2
                          : shape == 4 ? ((i - 1) & ~size_t(1))
                          : shape == 5 ? std::min(i - 1, n / 2)
                                       : (i == n / 2 ? 0 : i - 1);
          source[i].push_back(parent);
          source[parent].push_back(i);
        }
      }
      CheckDepthFirst<TypeParam>(source);
    }
  }
}

TYPED_TEST(TreeSpecificationTest, DepthFirstExhaustiveOrderedTrees) {
  // Every primitive balanced-parentheses string is one ordered rooted tree.
  for (unsigned n = 1; n <= 9; ++n) {
    for (uint64_t bits = 0; bits < (uint64_t{1} << (2 * n)); ++bits) {
      if (std::popcount(bits) != static_cast<int>(n)) {
        continue;
      }
      std::vector<std::vector<size_t>> source(n);
      std::vector<size_t> stack;
      size_t count = 0;
      bool valid = true;
      for (unsigned i = 0; i < 2 * n; ++i) {
        if ((bits >> i) & 1) {
          if (count >= n) {
            valid = false;
            break;
          }
          size_t node = count++;
          size_t parent = stack.empty() ? 0 : stack.back();
          source[node].push_back(parent);
          if (node) {
            source[parent].push_back(node);
          }
          stack.push_back(node);
        } else {
          if (stack.empty()) {
            valid = false;
            break;
          }
          stack.pop_back();
          if (stack.empty() && i + 1 < 2 * n) {
            valid = false;
            break;
          }
        }
      }
      if (valid && stack.empty()) {
        CheckDepthFirst<TypeParam>(source);
      }
    }
  }
}

// Embed every nonzero byte plus lookahead at a byte-aligned position in
// valid BP. Exercise the same ordered shapes in every encoding/backend.
TYPED_TEST(TreeSpecificationTest, DepthFirstAllBytePatterns) {
  for (unsigned pattern = 0; pattern < 512; ++pattern) {
    if (!(pattern & 255)) {
      continue;
    }
    SCOPED_TRACE(pattern);
    std::vector<std::vector<size_t>> source;
    std::vector<size_t> stack;
    const auto open = [&] {
      const size_t node = source.size();
      const size_t parent = stack.empty() ? 0 : stack.back();
      source.push_back({parent});
      if (node) {
        source[parent].push_back(node);
      }
      stack.push_back(node);
    };
    for (unsigned i = 0; i < 8; ++i) {
      open();
    }
    for (unsigned i = 0; i < 9; ++i) {
      if ((pattern >> i) & 1) {
        open();
      } else {
        stack.pop_back();
      }
    }
    CheckDepthFirst<TypeParam>(source);
  }
}

TYPED_TEST(TreeSpecificationTest, DepthFirstVisitorExceptionsAndReentrancy) {
  using Tree = typename TypeParam::Tree;
  std::mt19937_64 rng(42);
  const auto adjacency = TypeParam::order(257, generate_random_tree(257, rng));
  const auto words = TypeParam::encode(257, adjacency);
  const Tree tree(words, 257);
  for (bool unpack : {false, true}) {
    size_t visits = 0;
    auto stop = [&](pixie::TreeNode) {
      if (++visits == 17) {
        throw std::runtime_error("stop");
      }
    };
    if (unpack) {
      EXPECT_THROW(pixie::benchmarks::DepthFirst<true>(tree, words, stop),
                   std::runtime_error);
    } else {
      EXPECT_THROW(pixie::benchmarks::DepthFirst(tree, words, stop),
                   std::runtime_error);
    }
    EXPECT_EQ(visits, 17u);
  }
  size_t outer = 0, inner = 0;
  pixie::benchmarks::DepthFirst<true>(tree, words, [&](pixie::TreeNode) {
    if (++outer == 1) {
      pixie::benchmarks::DepthFirst<true>(tree, words,
                                          [&](pixie::TreeNode) { ++inner; });
    }
  });
  EXPECT_EQ(outer, inner);
  EXPECT_GT(outer, 257u);
}

TYPED_TEST(TreeSpecificationTest, BenchmarkChildInspectionAndRandomPaths) {
  using namespace pixie::benchmarks;
  using Tree = typename TypeParam::Tree;
  for (size_t n : {1, 17, 127, 128, 129, 257, 521}) {
    for (int shape = 0; shape < 6; ++shape) {
      SCOPED_TRACE(testing::Message() << "n=" << n << " shape=" << shape);
      std::mt19937_64 rng(42);
      auto source = generate_random_tree(n, rng);
      if (shape) {
        for (auto& edges : source) {
          edges.clear();
        }
        source[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          const size_t parent = shape == 1   ? i - 1
                                : shape == 2 ? 0
                                : shape == 3 ? (i - 1) / 2
                                : shape == 4 ? (i <= 2 ? 0 : 1)
                                             : ((i - 1) & ~size_t(1));
          source[i].push_back(parent);
          source[parent].push_back(i);
        }
      }
      const auto adjacency = TypeParam::order(n, source);
      auto words = TypeParam::encode(n, adjacency);
      const size_t bits = TypeParam::bit_size(n);
      if (bits % 64) {
        words.back() |= ~uint64_t{0} << (bits % 64);
      }
      const Tree tree(words, n);
      const AdjListTree reference(adjacency);
      std::vector<size_t> expected_events, expected_children;
      const auto dfs = [&](auto&& self, size_t node) -> void {
        expected_events.push_back(node);
        for (size_t i = 1; i < adjacency[node].size(); ++i) {
          expected_children.push_back(adjacency[node][i]);
        }
        for (size_t i = 1; i < adjacency[node].size(); ++i) {
          self(self, adjacency[node][i]);
        }
        if (adjacency[node].size() > 1) {
          expected_events.push_back(node);
        }
      };
      dfs(dfs, 0);
      const auto check_dfs = [&](const auto& candidate) {
        std::vector<size_t> events, children;
        const auto work = DepthFirstInspect(
            candidate, [&](auto node) { events.push_back(node.number); },
            [&](auto node) { children.push_back(node.number); });
        EXPECT_EQ(events, expected_events);
        EXPECT_EQ(children, expected_children);
        EXPECT_EQ(work.nodes, n);
        EXPECT_EQ(work.children, n - 1);
      };
      check_dfs(tree);
      check_dfs(reference);
      std::vector<size_t> expected_positions(n);
      auto remember = [&](pixie::TreeNode node) {
        expected_positions[node.number] = node.pos;
      };
      pixie::benchmarks::tree_dfs::cached(tree, remember);
      {
        std::vector<size_t> events, children;
        const auto work = UnpackedDepthFirstInspect(
            tree, words,
            [&](auto node) {
              events.push_back(node.number);
              EXPECT_EQ(node.pos, expected_positions[node.number]);
            },
            [&](auto node) {
              children.push_back(node.number);
              EXPECT_EQ(node.pos, expected_positions[node.number]);
            });
        EXPECT_EQ(events, expected_events);
        EXPECT_EQ(children, expected_children);
        EXPECT_EQ(work.nodes, n);
        EXPECT_EQ(work.children, n - 1);
      }

      ShiftedQueryPool expected_queries(0, std::numeric_limits<size_t>::max(),
                                        43);
      ShiftedQueryPool plain_queries(0, std::numeric_limits<size_t>::max(), 43);
      ShiftedQueryPool inspect_queries(0, std::numeric_limits<size_t>::max(),
                                       43);
      ShiftedQueryPool reference_queries(0, std::numeric_limits<size_t>::max(),
                                         43);
      ShiftedQueryPool unpacked_queries(0, std::numeric_limits<size_t>::max(),
                                        43);
      ShiftedQueryPool unpacked_inspect_queries(
          0, std::numeric_limits<size_t>::max(), 43);
      // Chains cross several pool wraps; singletons consume no choices.
      for (size_t trial = 0; trial < 128; ++trial) {
        std::vector<size_t> expected, children;
        size_t current = 0;
        while (true) {
          expected.push_back(current);
          const auto& edges = adjacency[current];
          if (edges.size() == 1) {
            break;
          }
          children.insert(children.end(), edges.begin() + 1, edges.end());
          current = edges[1 + expected_queries.next() % (edges.size() - 1)];
        }
        std::vector<size_t> plain, inspected, actual_children, ref_path,
            ref_children;
        auto unexpected = [](auto) {
          ADD_FAILURE() << "Unexpected inspection";
        };
        const auto a = RandomRootToLeaf<false>(
            tree, plain_queries, [&](auto v) { plain.push_back(v.number); },
            unexpected);
        const auto b = RandomRootToLeaf<true>(
            tree, inspect_queries,
            [&](auto v) { inspected.push_back(v.number); },
            [&](auto v) { actual_children.push_back(v.number); });
        const auto c = RandomRootToLeaf<true>(
            reference, reference_queries,
            [&](auto v) { ref_path.push_back(v.number); },
            [&](auto v) { ref_children.push_back(v.number); });
        std::vector<size_t> unpacked, unpacked_inspected, unpacked_children;
        const auto d = UnpackedRandomRootToLeaf<false>(
            tree, words, unpacked_queries,
            [&](auto v) { unpacked.push_back(v.number); }, unexpected);
        const auto e = UnpackedRandomRootToLeaf<true>(
            tree, words, unpacked_inspect_queries,
            [&](auto v) { unpacked_inspected.push_back(v.number); },
            [&](auto v) { unpacked_children.push_back(v.number); });
        ASSERT_EQ(unpacked, expected);
        ASSERT_EQ(unpacked_inspected, expected);
        ASSERT_EQ(unpacked_children, children);
        EXPECT_EQ(d.nodes, expected.size());
        EXPECT_EQ(d.children, 0u);
        EXPECT_EQ(e.nodes, expected.size());
        EXPECT_EQ(e.children, children.size());
        ASSERT_EQ(plain, expected);
        ASSERT_EQ(inspected, expected);
        ASSERT_EQ(ref_path, expected);
        ASSERT_EQ(actual_children, children);
        ASSERT_EQ(ref_children, children);
        EXPECT_EQ(a.nodes, expected.size());
        EXPECT_EQ(a.children, 0u);
        EXPECT_EQ(b.nodes, a.nodes);
        EXPECT_EQ(c.nodes, a.nodes);
        EXPECT_EQ(b.children, children.size());
        EXPECT_EQ(c.children, children.size());
      }
    }
  }
}

TYPED_TEST(TreeSpecificationTest, InspectedDfsCountsDeepFallbackEvents) {
  using Tree = typename TypeParam::Tree;
  constexpr size_t n = 2305;
  std::vector<std::vector<size_t>> adjacency(n);
  adjacency[0].push_back(0);
  for (size_t i = 1; i < n; ++i) {
    adjacency[i].push_back(i - 1);
    adjacency[i - 1].push_back(i);
  }
  const auto words = TypeParam::encode(n, adjacency);
  const Tree tree(words, n);
  std::vector<size_t> positions(n);
  auto node = tree.root();
  for (size_t i = 0; i < n; ++i) {
    positions[i] = node.pos;
    if (i + 1 < n) {
      node = tree.first_child(node);
    }
  }
  size_t events = 0, children = 0;
  const auto work = pixie::benchmarks::UnpackedDepthFirstInspect(
      tree, words,
      [&](auto current) {
        const auto expected = events < n ? events : 2 * n - 2 - events;
        EXPECT_EQ(current.number, expected);
        EXPECT_EQ(current.pos, positions[current.number]);
        ++events;
      },
      [&](auto child) {
        EXPECT_EQ(child.number, ++children);
        EXPECT_EQ(child.pos, positions[child.number]);
      });
  EXPECT_EQ(events, 2 * n - 1);
  EXPECT_EQ(children, n - 1);
  EXPECT_EQ(work.nodes, n);
  EXPECT_EQ(work.children, children);
}

template <bool BP, size_t Bits = 256>
void CheckBoundedUnpackingCache() {
  using RmM = pixie::RmMBTree<>;
  using Tree =
      std::conditional_t<BP, pixie::BPTree<RmM>, pixie::DFUDSTree<RmM>>;
  for (size_t n : {size_t{1}, Bits / 2 - 1, Bits / 2, Bits / 2 + 1}) {
    for (bool star : {false, true}) {
      std::vector<std::vector<size_t>> adjacency(n);
      adjacency[0].push_back(0);
      for (size_t i = 1; i < n; ++i) {
        const auto parent = star ? 0 : i - 1;
        adjacency[i].push_back(parent);
        adjacency[parent].push_back(i);
      }
      auto words = BP ? pixie::adj_to_bp(n, adjacency)
                      : pixie::adj_to_dfuds(n, adjacency);
      const size_t bits = 2 * n - (BP ? 0 : 1);
      if (bits % 64) {
        words.back() |= ~uint64_t{0} << (bits % 64);
      }
      const Tree tree(words, n);
      pixie::detail::tree_unpack::BoundedSubtreeCache<Bits> cache;
      const bool loaded =
          cache.template load<BP>(tree, words, bits, tree.root());
      ASSERT_EQ(loaded, n <= Bits / 2);
      ASSERT_EQ(cache.contains(tree.root()), loaded);
      if (loaded) {
        auto check = [&](pixie::TreeNode node) {
          ASSERT_TRUE(cache.contains(node));
          ASSERT_EQ(cache.degree(node), tree.degree(node));
          size_t index = 0;
          cache.children(node, [&](pixie::TreeNode child) {
            const auto expected = tree.child(node, index);
            const auto indexed = cache.child(node, index++);
            EXPECT_EQ(child.number, expected.number);
            EXPECT_EQ(child.pos, expected.pos);
            EXPECT_EQ(indexed.number, expected.number);
            EXPECT_EQ(indexed.pos, expected.pos);
          });
          EXPECT_EQ(index, tree.degree(node));
        };
        pixie::benchmarks::tree_dfs::cached(tree, check);
      } else {
        // The first child fits even though its parent did not; its handles
        // must retain the original tree's nonzero number and bit offset.
        const auto child = tree.first_child(tree.root());
        ASSERT_TRUE(cache.template load<BP>(tree, words, bits, child));
        EXPECT_TRUE(cache.contains(child));
        EXPECT_EQ(cache.degree(child), tree.degree(child));
        if (!tree.is_leaf(child)) {
          EXPECT_EQ(cache.child(child, 0).pos, tree.first_child(child).pos);
          EXPECT_EQ(cache.child(child, 0).number,
                    tree.first_child(child).number);
        }
        ASSERT_FALSE(cache.template load<BP>(tree, words, bits, tree.root()));
        EXPECT_FALSE(cache.contains(child));
      }
    }
  }
}
TEST(TreeUnpackingTest, ParentSlotFillPreservesAdjacentElements) {
  for (unsigned offset = 0; offset < 33; ++offset) {
    for (unsigned count = 0; count <= 129; ++count) {
      std::array<uint16_t, 192> values;
      values.fill(0xaaaa);
      pixie::detail::tree_unpack::fill_parent_slots(values.data() + offset,
                                                    count, 0x1234);
      for (size_t i = 0; i < values.size(); ++i) {
        EXPECT_EQ(values[i],
                  i >= offset && i < offset + count ? 0x1234 : 0xaaaa);
      }
    }
  }
}

TEST(TreeUnpackingTest, BpCacheAcceptsOnlyCompleteBoundedSubtrees) {
  CheckBoundedUnpackingCache<true>();
  CheckBoundedUnpackingCache<true, 1024>();
  CheckBoundedUnpackingCache<true, 2048>();
  CheckBoundedUnpackingCache<true, 4096>();
}
TEST(TreeUnpackingTest, DfudsCacheAcceptsOnlyCompleteBoundedSubtrees) {
  CheckBoundedUnpackingCache<false>();
  CheckBoundedUnpackingCache<false, 1024>();
  CheckBoundedUnpackingCache<false, 2048>();
  CheckBoundedUnpackingCache<false, 4096>();
}

// References have representation-specific handles, but share navigation and
// workload semantics. Run the same specification against both layouts.
template <class Tree>
class ReferenceTreeSpecification : public testing::Test {};
using ReferenceTrees =
    testing::Types<AdjListTree, pixie::references::PointerTree>;
TYPED_TEST_SUITE(ReferenceTreeSpecification, ReferenceTrees);

template <class Tree, class Node>
size_t ReferenceNumber(const Tree& tree, Node node) {
  if constexpr (requires { node.number; }) {
    return node.number;
  } else {
    return tree.node_number(node);
  }
}

TYPED_TEST(ReferenceTreeSpecification, NavigationAndTraversalMatchTopology) {
  using namespace pixie::benchmarks;
  for (size_t n : {1, 17, 129, 521}) {
    for (unsigned shape = 0; shape < 4; ++shape) {
      std::mt19937_64 rng(42 + n);
      auto source = generate_random_tree(n, rng);
      if (shape) {
        source.assign(n, {});
        source[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          const size_t parent = shape == 1   ? i - 1
                                : shape == 2 ? 0
                                             : (i - 1) / 2;
          source[i].push_back(parent);
          source[parent].push_back(i);
        }
      }
      for (bool preorder : {false, true}) {
        SCOPED_TRACE(testing::Message() << "n=" << n << " shape=" << shape
                                        << " preorder=" << preorder);
        const auto adj = preorder ? dfs_order(n, source) : bfs_order(n, source);
        const TypeParam tree(adj);
        auto number = [&](auto node) { return ReferenceNumber(tree, node); };
        using Node = decltype(tree.root());
        std::vector<Node> handles(n);
        std::queue<Node> pending;
        handles[0] = tree.root();
        pending.push(tree.root());
        while (!pending.empty()) {
          const auto node = pending.front();
          pending.pop();
          const size_t id = number(node);
          const auto& edges = adj[id];
          EXPECT_EQ(tree.is_root(node), id == 0);
          EXPECT_EQ(number(tree.parent(node)), edges[0]);
          EXPECT_EQ(tree.degree(node), edges.size() - 1);
          EXPECT_EQ(tree.is_leaf(node), edges.size() == 1);
          std::vector<size_t> children;
          tree.for_each_child(node, [&](auto child) {
            children.push_back(number(child));
            handles[number(child)] = child;
            pending.push(child);
          });
          EXPECT_EQ(children,
                    std::vector<size_t>(edges.begin() + 1, edges.end()));
          for (size_t i = 0; i < children.size(); ++i) {
            const auto child = tree.child(node, i);
            EXPECT_EQ(number(child), children[i]);
            if (i == 0) {
              EXPECT_EQ(number(tree.first_child(node)), children[i]);
            }
            EXPECT_EQ(tree.is_last_child(child), i + 1 == children.size());
            const auto next =
                i + 1 == children.size() ? children[i] : children[i + 1];
            EXPECT_EQ(number(tree.next_sibling(child)), next);
          }
        }
        std::vector<size_t> expected_events, expected_children;
        const auto dfs = [&](auto&& self, size_t id) -> void {
          expected_events.push_back(id);
          const auto& edges = adj[id];
          expected_children.insert(expected_children.end(), edges.begin() + 1,
                                   edges.end());
          for (size_t i = 1; i < edges.size(); ++i) {
            self(self, edges[i]);
          }
          if (edges.size() > 1) {
            expected_events.push_back(id);
          }
        };
        dfs(dfs, 0);
        std::vector<size_t> events, children;
        const auto work = DepthFirstInspect(
            tree, [&](auto node) { events.push_back(number(node)); },
            [&](auto node) { children.push_back(number(node)); });
        EXPECT_EQ(events, expected_events);
        EXPECT_EQ(children, expected_children);
        EXPECT_EQ(work.nodes, n);
        EXPECT_EQ(work.children, n - 1);

        ShiftedQueryPool expected_queries(0, std::numeric_limits<size_t>::max(),
                                          43);
        ShiftedQueryPool plain_queries(0, std::numeric_limits<size_t>::max(),
                                       43);
        ShiftedQueryPool inspected_queries(
            0, std::numeric_limits<size_t>::max(), 43);
        // Long chains also exercise shifted-pool wraparound.
        for (size_t trial = 0; trial < 128; ++trial) {
          std::vector<size_t> expected_path, expected_inspection;
          size_t id = 0;
          while (true) {
            expected_path.push_back(id);
            const auto& edges = adj[id];
            if (edges.size() == 1) {
              break;
            }
            expected_inspection.insert(expected_inspection.end(),
                                       edges.begin() + 1, edges.end());
            id = edges[1 + expected_queries.next() % (edges.size() - 1)];
          }
          std::vector<size_t> plain_path, inspected_path, inspection;
          const auto plain = RandomRootToLeaf<false>(
              tree, plain_queries,
              [&](auto node) { plain_path.push_back(number(node)); },
              [](auto) {
                ADD_FAILURE() << "Plain paths must not inspect children";
              });
          const auto inspected = RandomRootToLeaf<true>(
              tree, inspected_queries,
              [&](auto node) { inspected_path.push_back(number(node)); },
              [&](auto node) { inspection.push_back(number(node)); });
          EXPECT_EQ(plain_path, expected_path);
          EXPECT_EQ(inspected_path, expected_path);
          EXPECT_EQ(inspection, expected_inspection);
          EXPECT_EQ(plain.nodes, expected_path.size());
          EXPECT_EQ(plain.children, 0u);
          EXPECT_EQ(inspected.nodes, expected_path.size());
          EXPECT_EQ(inspected.children, expected_inspection.size());
        }
      }
    }
  }
}

TYPED_TEST(ReferenceTreeSpecification, ChildVisitorIsBorrowedAndReentrant) {
  const TypeParam tree({{0, 1, 2}, {0, 3}, {0}, {1}});
  using Node = decltype(tree.root());
  struct Visitor {
    Visitor() = default;
    Visitor(const Visitor&) = delete;
    std::vector<Node> seen;
    void operator()(Node node) { seen.push_back(node); }
  } visitor;
  tree.for_each_child(tree.root(), visitor);
  ASSERT_EQ(visitor.seen.size(), 2u);
  size_t nested_count = 0;
  tree.for_each_child(tree.root(), [&](auto child) {
    tree.for_each_child(child, [&](auto) { ++nested_count; });
  });
  EXPECT_EQ(nested_count, 1u);
  EXPECT_THROW(
      tree.for_each_child(tree.root(),
                          [](auto) { throw std::runtime_error("visitor"); }),
      std::runtime_error);
  EXPECT_EQ(tree.degree(tree.root()), 2u);
}

TEST(PointerTreeTest, OwnershipMovesAndMissingLinks) {
  using Tree = pixie::references::PointerTree;
  static_assert(!std::is_copy_constructible_v<Tree>);
  Tree empty({});
  EXPECT_TRUE(empty.empty());
  EXPECT_EQ(empty.size(), 0u);
  EXPECT_EQ(empty.root(), nullptr);
  Tree singleton({{0}});
  EXPECT_EQ(singleton.parent(singleton.root()), singleton.root());
  EXPECT_EQ(singleton.first_child(singleton.root()), singleton.root());
  EXPECT_EQ(singleton.next_sibling(singleton.root()), singleton.root());
  EXPECT_TRUE(singleton.is_last_child(singleton.root()));
  auto original = [] {
    const std::vector<std::vector<size_t>> source = {{0, 1, 2}, {0}, {0}};
    return Tree(source);
  }();
  const auto root = original.root();
  const auto first = original.first_child(root);
  Tree moved(std::move(original));
  EXPECT_EQ(moved.root(), root);
  EXPECT_EQ(moved.parent(first), root);
  EXPECT_EQ(moved.node_number(moved.next_sibling(first)), 2u);
  singleton = std::move(moved);
  EXPECT_EQ(singleton.root(), root);
  EXPECT_EQ(singleton.first_child(root), first);
  EXPECT_EQ(singleton.degree(root), 2u);
  EXPECT_EQ(singleton.size(), 3u);
}

}  // namespace

namespace {
// Compare one interleaved callback stream: every child's inspection must happen
// before descent, not merely produce the right independent lists of handles.
void CheckMicroPrefetch(const std::vector<uint64_t>& words, size_t n) {
  using Tree = pixie::BPTree<pixie::RmMBTree<>>;
  const Tree tree(words, n);
  using Event = std::tuple<size_t, size_t, unsigned>;
  std::vector<Event> expected, actual;
  auto current = tree.root();
  bool entering = true;
  while (true) {
    expected.emplace_back(current.number, current.pos, entering ? 0 : 1);
    if (entering) {
      tree.for_each_child(current, [&](auto child) {
        expected.emplace_back(child.number, child.pos, 2);
      });
      if (!tree.is_leaf(current)) {
        current = tree.first_child(current);
        continue;
      }
    }
    entering = false;
    if (tree.is_root(current)) {
      break;
    }
    if (tree.is_last_child(current)) {
      current = tree.parent(current);
    } else {
      current = tree.next_sibling(current);
      entering = true;
    }
  }
  auto visit = [&, token = std::make_unique<int>(1)](auto node, bool entry) {
    actual.emplace_back(node.number, node.pos, entry ? 0 : 1);
  };
  auto inspect = [&](auto node) {
    actual.emplace_back(node.number, node.pos, 2);
  };
  pixie::detail::bp_micro::BpMicroTraversal<2048, Tree, decltype(visit),
                                            decltype(inspect)>
      small(tree, words, visit, inspect);
  small.walk(tree.root());
  ASSERT_EQ(actual, expected);
  EXPECT_EQ(small.nodes, n);
  EXPECT_EQ(small.children, n - 1);
  actual.clear();
  pixie::detail::bp_micro::BpMicroTraversal<8192, Tree, decltype(visit),
                                            decltype(inspect)>
      large(tree, words, visit, inspect);
  large.walk(tree.root());
  ASSERT_EQ(actual, expected);
  EXPECT_EQ(large.nodes, n);
  EXPECT_EQ(large.children, n - 1);
}
}  // namespace

TEST(TreeMicroPrefetchTest, ExhaustiveShortSubtrees) {
  size_t shapes = 0;
  for (unsigned key = 0; key < 8192; ++key) {
    const unsigned length =
        pixie::detail::bp_micro::bp_micro_tables.headers[key] & 15;
    if (!length || key != (key & ((1u << (length - 1)) - 1))) {
      continue;
    }
    ++shapes;
    CheckMicroPrefetch({(uint64_t(key) << 1) | 1}, length / 2);
  }
  EXPECT_EQ(shapes, 197u);
}
TEST(TreeMicroPrefetchTest, BoundariesAndDeepShapes) {
  for (size_t n : {1, 63, 64, 65, 127, 128, 129, 1023, 1024, 1025, 4097}) {
    for (int shape = 0; shape < 4; ++shape) {
      SCOPED_TRACE(testing::Message() << "n=" << n << " shape=" << shape);
      std::mt19937_64 rng(73 + n);
      auto source = generate_random_tree(n, rng);
      if (shape) {
        source.assign(n, {});
        source[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          size_t parent = shape == 1 ? 0 : shape == 2 ? i - 1 : (i - 1) / 2;
          source[parent].push_back(i);
          source[i].push_back(parent);
        }
      }
      auto adjacency = dfs_order(n, source);
      auto words = pixie::adj_to_bp(n, adjacency);
      if ((2 * n) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n) % 64);
      }
      CheckMicroPrefetch(words, n);
    }
  }
}

namespace {
void CheckLoudsPrefetch(const std::vector<uint64_t>& words, size_t n) {
  using Tree = pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  const Tree tree(words, n);
  using Event = std::tuple<size_t, size_t, unsigned>;
  std::vector<Event> expected, actual;
  auto current = tree.root();
  bool entering = true;
  while (true) {
    expected.emplace_back(current.number, current.pos, entering ? 0 : 1);
    if (entering) {
      tree.for_each_child(current, [&](auto child) {
        expected.emplace_back(child.number, child.pos, 2);
      });
      if (!tree.is_leaf(current)) {
        current = tree.first_child(current);
        continue;
      }
    }
    entering = false;
    if (tree.is_root(current)) {
      break;
    }
    if (tree.is_last_child(current)) {
      current = tree.parent(current);
    } else {
      current = tree.next_sibling(current);
      entering = true;
    }
  }
  auto visit = [&, token = std::make_unique<int>(1)](auto node, bool entry) {
    actual.emplace_back(node.number, node.pos, entry ? 0 : 1);
  };
  bool nested_checked = false;
  auto inspect = [&](auto node) {
    actual.emplace_back(node.number, node.pos, 2);
    if (!nested_checked) {
      nested_checked = true;
      auto ignore = [](auto) {};
      pixie::detail::tree_unpack::LoudsPrefetchTraversal<
          Tree, decltype(ignore), decltype(ignore), 8, true>
          nested(tree, words, ignore, ignore);
      nested.walk(tree.root(), 0);
      EXPECT_EQ(nested.node_count, n);
      EXPECT_EQ(nested.child_count, n - 1);
    }
  };
  pixie::detail::tree_unpack::LoudsPrefetchTraversal<Tree, decltype(visit),
                                                     decltype(inspect), 8, true>
      small(tree, words, visit, inspect);
  small.walk(tree.root(), 0);
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(small.node_count, n);
  ASSERT_EQ(small.child_count, n - 1);
  actual.clear();
  pixie::detail::tree_unpack::LoudsPrefetchTraversal<
      Tree, decltype(visit), decltype(inspect), 32, true>
      large(tree, words, visit, inspect);
  large.walk(tree.root(), 0);
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(large.node_count, n);
  ASSERT_EQ(large.child_count, n - 1);
  actual.clear();
  pixie::detail::tree_unpack::LoudsPrefetchTraversal<
      Tree, decltype(visit), decltype(inspect), 32, false>
      scalar_cursor(tree, words, visit, inspect);
  scalar_cursor.walk(tree.root(), 0);
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(scalar_cursor.node_count, n);
  ASSERT_EQ(scalar_cursor.child_count, n - 1);
}

}  // namespace

TEST(TreeLoudsPrefetchTest, BoundariesShapesAndNestedTraversal) {
  for (size_t n : {1, 7, 8, 9, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1023,
                   1024, 1025, 4097}) {
    for (int shape = 0; shape < 4; ++shape) {
      std::mt19937_64 rng(32 + n);
      auto adj = generate_random_tree(n, rng);
      if (shape) {
        adj.assign(n, {});
        adj[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          size_t p = shape == 1 ? 0 : shape == 2 ? i - 1 : (i - 1) / 2;
          adj[p].push_back(i);
          adj[i].push_back(p);
        }
      }
      adj = bfs_order(n, adj);
      auto words = adj_to_louds(n, adj);
      if ((2 * n - 1) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n - 1) % 64);
      }
      CheckLoudsPrefetch(words, n);
    }
  }
}

namespace {
void CheckDfudsMicroPrefetch(const std::vector<uint64_t>& words, size_t n) {
  using Tree = pixie::DFUDSTree<pixie::RmMBTree<>>;
  const Tree tree(words, n);
  using Event = std::tuple<size_t, size_t, unsigned>;
  std::vector<Event> expected, actual;
  auto current = tree.root();
  bool entering = true;
  while (true) {
    expected.emplace_back(current.number, current.pos, entering ? 0 : 1);
    if (entering) {
      tree.for_each_child(current, [&](auto child) {
        expected.emplace_back(child.number, child.pos, 2);
      });
      if (!tree.is_leaf(current)) {
        current = tree.first_child(current);
        continue;
      }
    }
    entering = false;
    if (tree.is_root(current)) {
      break;
    }
    if (tree.is_last_child(current)) {
      current = tree.parent(current);
    } else {
      current = tree.next_sibling(current);
      entering = true;
    }
  }
  auto visit = [&, token = std::make_unique<int>(1)](auto node, bool entry) {
    actual.emplace_back(node.number, node.pos, entry ? 0 : 1);
  };
  auto inspect = [&](auto node) {
    actual.emplace_back(node.number, node.pos, 2);
  };
  pixie::detail::dfuds_micro::DfudsMicroTraversal<2048, Tree, decltype(visit),
                                                  decltype(inspect)>
      small(tree, words, visit, inspect);
  small.walk(tree.root());
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(small.nodes, n);
  ASSERT_EQ(small.children, n - 1);
  actual.clear();
  pixie::detail::dfuds_micro::DfudsMicroTraversal<8192, Tree, decltype(visit),
                                                  decltype(inspect)>
      large(tree, words, visit, inspect);
  large.walk(tree.root());
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(large.nodes, n);
  ASSERT_EQ(large.children, n - 1);
  using Decoder = pixie::detail::dfuds_micro::DfudsMicroSubtree<8192>;
  pixie::experimental::PreparedRootSubtree<Decoder, 8192> prepared(
      words, 2 * n - 1, tree.root());
  for (int repeat = 0; repeat < 2; ++repeat) {
    actual.clear();
    const bool hit = prepared.try_walk(visit, inspect);
    ASSERT_EQ(hit, n <= 4096);
    if (hit) {
      ASSERT_EQ(actual, expected);
      ASSERT_EQ(prepared.size(), n);
    } else {
      ASSERT_TRUE(actual.empty());
      ASSERT_EQ(prepared.size(), 0u);
    }
  }
  pixie::experimental::PreparedTreeEvents<> events(tree);
  for (int repeat = 0; repeat < 2; ++repeat) {
    actual.clear();
    const bool hit = repeat ? events.template try_walk<true>(visit, inspect)
                            : events.try_walk(visit, inspect);
    ASSERT_EQ(hit, n <= 4096);
    ASSERT_EQ(events.size(), hit ? n : 0u);
    if (hit) {
      ASSERT_EQ(actual, expected);
    } else {
      ASSERT_TRUE(actual.empty());
    }
  }
  std::vector<Event> visits_only;
  for (const auto& event : expected) {
    if (std::get<2>(event) != 2) {
      visits_only.push_back(event);
    }
  }
  actual.clear();
  const bool plain_hit = events.template try_walk<false, false>(visit, inspect);
  ASSERT_EQ(plain_hit, n <= 4096);
  if (plain_hit) {
    ASSERT_EQ(actual, visits_only);
  } else {
    ASSERT_TRUE(actual.empty());
  }
  actual.clear();
  EXPECT_EQ((events.template try_walk<true, false>(visit, inspect)), plain_hit);
  EXPECT_EQ(actual, plain_hit ? visits_only : std::vector<Event>{});
}
}  // namespace

TEST(TreeDfudsMicroPrefetchTest, ExhaustiveShapesAndBoundaries) {
  for (unsigned key = 0; key < 8192; ++key) {
    auto length =
        pixie::detail::dfuds_micro::dfuds_micro_tables.headers[key] & 15;
    if (length && key < (1u << length)) {
      CheckDfudsMicroPrefetch({uint64_t(key)}, (length + 1) / 2);
    }
  }
  for (size_t n : {1, 63, 64, 65, 127, 128, 129, 1023, 1024, 1025, 4095, 4096,
                   4097, 16385}) {
    for (int shape = 0; shape < 4; ++shape) {
      std::mt19937_64 rng(32 + n);
      auto adj = generate_random_tree(n, rng);
      if (shape) {
        adj.assign(n, {});
        adj[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          size_t p = shape == 1 ? 0 : shape == 2 ? i - 1 : (i - 1) / 2;
          adj[p].push_back(i);
          adj[i].push_back(p);
        }
      }
      adj = dfs_order(n, adj);
      auto words = pixie::adj_to_dfuds(n, adj);
      CheckDfudsMicroPrefetch(words, n);
      if ((2 * n - 1) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n - 1) % 64);
      }
      CheckDfudsMicroPrefetch(words, n);
    }
  }
}

TEST(TreePreparedEventsTest, EmptyAndNestedReplay) {
  using Tree = pixie::DFUDSTree<pixie::RmMBTree<>>;
  const std::vector<uint64_t> empty_words;
  const Tree empty(empty_words, 0);
  const pixie::experimental::PreparedTreeEvents<> empty_cache(empty);
  size_t outer_visits = 0, outer_children = 0, inner_visits = 0,
         inner_children = 0;
  auto inner_visit = [&](auto) { ++inner_visits; };
  auto inner_inspect = [&](auto) { ++inner_children; };
  EXPECT_FALSE(empty_cache.try_walk(inner_visit, inner_inspect));
  EXPECT_EQ(inner_visits + inner_children, 0u);
  const std::vector<uint64_t> words{5};  // Three-node chain: 10 10 0.
  const Tree tree(words, 3);
  const pixie::experimental::PreparedTreeEvents<> cache(tree);
  bool nested = false;
  auto visit = [&](auto) {
    ++outer_visits;
    if (!nested) {
      nested = true;
      EXPECT_TRUE(cache.try_walk(inner_visit, inner_inspect));
    }
  };
  auto inspect = [&](auto) { ++outer_children; };
  EXPECT_TRUE(cache.try_walk(visit, inspect));
  EXPECT_EQ(outer_visits, 5u);
  EXPECT_EQ(inner_visits, 5u);
  EXPECT_EQ(outer_children, 2u);
  EXPECT_EQ(inner_children, 2u);
}
