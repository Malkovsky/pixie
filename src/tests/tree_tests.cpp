#include <gtest/gtest.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/bp_micro_tables.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/experimental/decomposed_tree.h>
#include <pixie/tree/experimental/layered_dfuds.h>
#include <pixie/tree/experimental/simd_micro_traversal.h>
#include <pixie/tree/experimental/tree_prepared_events.h>
#include <pixie/tree/experimental/tree_prepared_root.h>
#include <pixie/tree/forest_cover.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/path_view.h>
#include <pixie/tree/traversal_context.h>
#include <pixie/utils.h>
#include <references/bfs_csr_tree.h>
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
#include <span>
#include <stack>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../benchmarks/shifted_query_pool.h"
#include "../benchmarks/tree_dfs/optimized.h"
#include "../benchmarks/tree_traversal_workloads.h"
#include "tree_encoding_test_utils.h"

namespace {

template <class Tree>
class TreePathTest : public testing::Test {};

using PathTrees =
    testing::Types<AdjListTree,
                   pixie::references::PointerTree,
                   pixie::references::BfsCsrTree,
                   pixie::tree::BPTree<pixie::RmMBTree<>>,
                   pixie::tree::DFUDSTree<pixie::RmMBTree<>>,
                   pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>,
                   pixie::tree::ForestCoverTree>;
TYPED_TEST_SUITE(TreePathTest, PathTrees);

TYPED_TEST(TreePathTest, InclusivePathsInBothDirections) {
  // Input IDs: 0 -> {6, 3}; 6 -> {1, 2, 7, ..., 307}; 1 -> 4; 3 -> 5.
  std::vector<std::vector<size_t>> adj{{0, 6, 3}, {6, 4}, {6},      {0, 5},
                                       {1},       {3},    {0, 1, 2}};
  for (size_t i = 7; i <= 307; ++i) {
    adj[6].push_back(i);
    adj.push_back({6});
  }
  const TypeParam tree = [&] {
    if constexpr (std::is_constructible_v<TypeParam, decltype(adj)&>) {
      return TypeParam(adj);
    } else {
      return TypeParam(AdjListTree(adj));
    }
  }();
  // Resolve handles by child order, independent of each encoding's numbering.
  const auto root = tree.root();
  const auto six = tree.child(root, 0), three = tree.child(root, 1);
  const auto one = tree.child(six, 0), two = tree.child(six, 1);
  const auto four = tree.child(one, 0), five = tree.child(three, 0);
  const auto last = tree.child(six, adj[6].size() - 2);
  if constexpr (std::is_same_v<TypeParam, pixie::tree::ForestCoverTree>) {
    ASSERT_NE(one.unit, last.unit);
  }

  using Node = decltype(tree.root());
  std::vector<Node> out, scratch;
  out.reserve(16);
  scratch.reserve(16);
  const auto out_capacity = out.capacity(),
             scratch_capacity = scratch.capacity();
  auto check = [&](std::vector<Node> expected) {
    for (int direction = 0; direction < 2; ++direction) {
      auto compare = [&](const auto& actual) {
        ASSERT_EQ(actual.size(), expected.size());
        for (size_t i = 0; i < expected.size(); ++i) {
          if constexpr (std::is_same_v<Node, AdjListNode>) {
            EXPECT_EQ(actual[i].number, expected[i].number) << i;
          } else {
            EXPECT_EQ(actual[i], expected[i]) << i;
          }
        }
      };
      const auto u = expected.front(), v = expected.back();
      compare(tree.tree_path(u, v));
      out.push_back(root);
      scratch.push_back(root);
      tree.tree_path(u, v, out, scratch);
      compare(out);
      EXPECT_TRUE(scratch.empty());
      EXPECT_EQ(out.capacity(), out_capacity);
      EXPECT_EQ(scratch.capacity(), scratch_capacity);
      std::ranges::reverse(expected);
    }
  };
  check({four, one, six, root, three, five});
  check({four, one, six, last});
  check({one, six, two});
  check({six, one, four});
  check({root, three, five});
  check({four});
  check({root});
  EXPECT_THROW(tree.tree_path(four, five, out, out), std::invalid_argument);
}

TEST(TreeConversion, ExhaustiveSourceEncodingsAndNormalizedEvents) {
  using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
  using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  using Louds = pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  for (size_t n = 1; n <= 7; ++n) {
    for (uint64_t bits = 1; bits < (uint64_t{1} << (2 * n)); bits += 2) {
      int excess = 0;
      bool valid = true;
      for (size_t i = 0; i < 2 * n; ++i) {
        excess += (bits >> i & 1) ? 1 : -1;
        if (excess <= 0 && i + 1 != 2 * n) {
          valid = false;
          break;
        }
      }
      if (!valid || excess) {
        continue;
      }
      const Bp bp(std::span(&bits, 1), n);
      const Dfuds dfuds(bp);
      const Louds louds(dfuds);
      EXPECT_EQ(Bp(dfuds).encoding_words()[0], bits);
      EXPECT_EQ(Bp(louds).encoding_words()[0], bits);
      auto check_events = [&](const auto& tree) {
        std::vector<pixie::tree::TreeNode> path;
        size_t entries = 0, exits = 0, edges = 0;
        pixie::tree::detail::construction_dfs(
            tree,
            [&](auto node, size_t degree) {
              EXPECT_EQ(degree, tree.degree(node));
              if (!path.empty()) {
                EXPECT_EQ(tree.parent(node).number, path.back().number);
                EXPECT_EQ(tree.parent(node).pos, path.back().pos);
              }
              path.push_back(node);
              ++entries;
              edges += degree;
            },
            [&](auto node) {
              ASSERT_FALSE(path.empty());
              EXPECT_EQ(path.back().number, node.number);
              EXPECT_EQ(path.back().pos, node.pos);
              path.pop_back();
              ++exits;
            });
        EXPECT_TRUE(path.empty());
        EXPECT_EQ(entries, n);
        EXPECT_EQ(exits, n);
        EXPECT_EQ(edges, n - 1);
      };
      check_events(bp);
      check_events(dfuds);
      check_events(louds);
    }
  }
}

TEST(TreeConversion, KnownOrderedEncodings) {
  using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
  using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  using Louds = pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  const uint64_t star = 0b001011, chain = 0b000111;
  const Bp star_tree(std::span(&star, 1), 3),
      chain_tree(std::span(&chain, 1), 3);
  EXPECT_EQ(Dfuds(star_tree).encoding_words()[0], 0b00011);
  EXPECT_EQ(Louds(star_tree).encoding_words()[0], 0b00011);
  EXPECT_EQ(Dfuds(chain_tree).encoding_words()[0], 0b00101);
  EXPECT_EQ(Louds(chain_tree).encoding_words()[0], 0b00101);
}

TEST(TreeConversion, OrderedTopologyAcrossSourcesAndDestinations) {
  using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
  using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  using Louds = pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>;
  static_assert(!std::is_constructible_v<Bp::TreeView, const Louds&>);
  static_assert(!std::is_constructible_v<Dfuds::TreeView, const Bp&>);
  static_assert(!std::is_constructible_v<Louds::TreeView, const Dfuds&>);
  std::mt19937_64 rng(42);
  for (size_t n : {size_t{0}, size_t{1}, size_t{2}, size_t{33}, size_t{257},
                   size_t{4097}, size_t{8193}}) {
    for (int shape = 0; shape < 3; ++shape) {
      std::vector<std::vector<size_t>> adjacency(n);
      if (n) {
        adjacency[0].push_back(0);
      }
      for (size_t i = 1; i < n; ++i) {
        const size_t parent = shape == 0 ? i - 1 : shape == 1 ? 0 : rng() % i;
        adjacency[i].push_back(parent);
        adjacency[parent].push_back(i);
      }
      const auto expected_bp =
          n ? EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(
                  n, adjacency)
            : std::vector<uint64_t>{};
      const auto expected_dfuds =
          n ? EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(
                  n, adjacency)
            : std::vector<uint64_t>{};
      const auto expected_louds =
          n ? EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(
                  n, bfs_order(n, adjacency))
            : std::vector<uint64_t>{};
      auto check = [&](const auto& source) {
        Bp bp(source);
        Dfuds dfuds(source);
        Louds louds(source);
        EXPECT_EQ(bp.size(), n);
        EXPECT_EQ(dfuds.size(), n);
        EXPECT_EQ(louds.size(), n);
        EXPECT_TRUE(std::ranges::equal(bp.encoding_words(), expected_bp));
        EXPECT_TRUE(std::ranges::equal(dfuds.encoding_words(), expected_dfuds));
        EXPECT_TRUE(std::ranges::equal(louds.encoding_words(), expected_louds));
      };
      AdjListTree reference(adjacency);
      pixie::references::PointerTree pointers(adjacency);
      check(reference);
      check(pointers);
      check(pixie::references::BfsCsrTree(pointers));
      Bp bp(reference);
      Dfuds dfuds(reference);
      Louds louds(reference);
      check(bp);
      check(dfuds);
      check(louds);
      check(Bp::TreeView(bp.encoding_words(), n));
      check(Dfuds::TreeView(dfuds.encoding_words(), n));
      check(Louds::TreeView(louds.encoding_words(), n));
      // Conversion owns its encoding and does not borrow the temporary source.
      Bp independent{Louds(reference)};
      check(independent);
    }
  }
}

template <class RankSelect>
struct LoudsCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes - 1; }
  using Tree = pixie::tree::LoudsTree<RankSelect>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return bfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(
        tree_size, adjacency);
  }
};

template <class RmM>
struct BpCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes; }
  using Tree = pixie::tree::BPTree<RmM>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return dfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(
        tree_size, adjacency);
  }
};

template <class RmM>
struct DfudsCase {
  static size_t bit_size(size_t nodes) { return 2 * nodes - 1; }
  using Tree = pixie::tree::DFUDSTree<RmM>;

  static std::vector<std::vector<std::size_t>> order(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return dfs_order(tree_size, adjacency);
  }

  static std::vector<std::uint64_t> encode(
      std::size_t tree_size,
      const std::vector<std::vector<std::size_t>>& adjacency) {
    return EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(
        tree_size, adjacency);
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

TYPED_TEST(TreeSpecificationTest, StableNumberLookup) {
  using Tree = typename TypeParam::Tree;
  using Base = pixie::tree::TreeBase<Tree>;
  static_assert(!requires(const Base& tree) { tree.encoding_words(); });
  static_assert(!requires(const Base& tree) { tree.encoding_owned_bytes(); });
  static_assert(!requires { typename Base::TreeView; });
  static_assert(!pixie::Serializable<Base>);
  const Tree empty(std::span<const uint64_t>{}, 0);
  EXPECT_THROW(empty.find_by_number(0), std::out_of_range);

  std::mt19937_64 rng(42);
  constexpr size_t n = 257;
  const auto adjacency = TypeParam::order(n, generate_random_tree(n, rng));
  const auto words = TypeParam::encode(n, adjacency);
  const Tree tree(words, n);
  const Base& facade = tree;
  EXPECT_THROW(facade.find_by_number(n), std::out_of_range);
  EXPECT_THROW(facade.find_by_number(std::numeric_limits<size_t>::max()),
               std::out_of_range);
  std::vector<std::pair<size_t, typename Tree::Node>> stack{{0, facade.root()}};
  auto view = facade.navigate(facade.root());
  size_t visited = 0;
  while (!stack.empty()) {
    const auto [expected, node] = stack.back();
    stack.pop_back();
    ++visited;
    EXPECT_EQ(facade.number(node), expected);
    EXPECT_EQ(facade.find_by_number(expected), node);
    view.reset(node);
    EXPECT_EQ(view.number(), expected);
    size_t j = 1;  // Adjacency entries begin with the parent.
    facade.for_each_child(node, [&](auto child) {
      ASSERT_LT(j, adjacency[expected].size());
      stack.emplace_back(adjacency[expected][j++], child);
    });
    EXPECT_EQ(j, adjacency[expected].size());
  }
  EXPECT_EQ(visited, n);
}

TYPED_TEST(TreeSpecificationTest, IndependentConcurrentViews) {
  using Tree = typename TypeParam::Tree;
  static_assert(pixie::tree::TreeNavigationView<typename Tree::NodeView>);
  static_assert(std::is_constructible_v<typename Tree::PathView, Tree&>);
  static_assert(std::is_constructible_v<typename Tree::PathView, const Tree&>);
  static_assert(!std::is_constructible_v<typename Tree::PathView, Tree&&>);
  static_assert(
      !std::is_constructible_v<typename Tree::PathView, const Tree&&>);
  static_assert(requires(const Tree& tree) { tree.navigate({0, 0}); });
  static_assert(!requires(Tree&& tree) { std::move(tree).navigate({0, 0}); });
  static_assert(
      !requires(const Tree&& tree) { std::move(tree).navigate({0, 0}); });
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
    auto visit = [&](pixie::tree::TreeNode node, bool entering) {
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
    pixie::tree::detail::tree_unpack::subtree(tree, tree.root(), visit);
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
      std::is_same_v<typename Tree::PathView, pixie::tree::TreePathView<Tree>>);
  static_assert(std::is_same_v<typename Tree::template BoundedPathView<2>,
                               pixie::tree::TreePathView<Tree, 2>>);
  std::mt19937_64 rng(91);
  constexpr size_t n = 257;
  const auto adjacency = TypeParam::order(n, generate_random_tree(n, rng));
  auto words = TypeParam::encode(n, adjacency);
  const auto expected_words = words;
  auto check = [&](const auto& tree) {
    EXPECT_EQ(tree.size(), n);
    EXPECT_TRUE(std::ranges::equal(tree.encoding_words(), expected_words));
    std::vector<pixie::tree::TreeNode> pending{tree.root()};
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
  pixie::tree::TreeTraversalContext context(view, 32768);
  EXPECT_EQ(context.source_words().data(), view.encoding_words().data());
  if constexpr (owns) {
    Tree from_lvalue = [&] {
      pixie::AlignedStorage source{std::span<const uint64_t>(expected_words)};
      Tree result(source, n);
      EXPECT_NE(result.encoding_words().data(), source.as_words64().data());
      source = pixie::AlignedStorage(source.size_bits());
      return result;
    }();
    check(from_lvalue);
    pixie::AlignedStorage storage(TypeParam::bit_size(n));
    std::ranges::copy(std::as_bytes(std::span(expected_words)),
                      std::as_writable_bytes(storage.as_lines()).begin());
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
    std::ranges::copy(std::as_bytes(std::span(words)),
                      std::as_writable_bytes(storage.as_lines()).begin());
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
  const auto words =
      EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(4,
                                                                     adjacency);
  const pixie::tree::DFUDSTree<RejectZeroSelectRmM> tree(words, 4);
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
  using Context = pixie::tree::TreeTraversalContext<Tree>;
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

TYPED_TEST(BoundedTraversalTest, FetchMatchesSubtreeBitPrefixes) {
  using Tree = typename TypeParam::Tree;
  std::mt19937_64 rng(614);
  for (size_t n : {1u, 2u, 33u, 257u}) {
    auto adjacency = TypeParam::order(n, generate_random_tree(n, rng));
    Tree tree(TypeParam::encode(n, adjacency), n);
    pixie::tree::RootNavigationIndex index(tree, 32768);
    pixie::tree::TreeTraversalContext context(tree, 524288, &index);
    Tree other(TypeParam::encode(n, adjacency), n);
    EXPECT_THROW((pixie::tree::TreeTraversalContext(other, 524288, &index)),
                 std::invalid_argument);
    // Independent oracle: enumerate descendants from the adjacency list. Source
    // order is preorder for BP/DFUDS and level order for LOUDS.
    const auto children = [&](size_t id) {
      return std::span<const size_t>(adjacency[id]).subspan(1);
    };
    for (size_t root : {size_t{0}, n / 2, n > 1 ? size_t{1} : size_t{0}}) {
      std::vector<size_t> ids{root}, subtree(n, 1), local(n), start(n), end(n);
      for (size_t i = 0; i < ids.size(); ++i) {
        const auto row = children(ids[i]);
        ASSERT_LE(row.size(), n - ids.size());
        ids.insert(ids.end(), row.begin(), row.end());
      }
      std::sort(ids.begin(), ids.end());
      for (auto it = ids.rbegin(); it != ids.rend(); ++it) {
        for (auto c : children(*it)) {
          subtree[*it] += subtree[c];
        }
      }
      const auto handle = tree.find_by_number(root);
      const bool bp = TypeParam::bit_size(1) == 2;
      size_t unary_offset = 0;
      for (size_t i = 0; i < ids.size(); ++i) {
        const auto id = ids[i];
        local[id] = i;
        start[id] =
            bp ? tree.find_by_number(id).pos - handle.pos : unary_offset;
        unary_offset += children(id).size() + 1;
        end[id] = bp ? start[id] + 2 * subtree[id] : unary_offset;
      }
      for (size_t bits = 0; bits <= TypeParam::bit_size(ids.size()) + 1;
           ++bits) {
        auto part = context.fetch(handle, bits);
        const size_t count = std::count_if(
            ids.begin(), ids.end(), [&](auto id) { return start[id] < bits; });
        ASSERT_EQ(part.size(), count);
        ASSERT_EQ(part.root_count(), count ? 1u : 0u);
        for (const auto node : part.nodes()) {
          const auto i = part.number(node);
          const auto id = ids[i];
          EXPECT_EQ(part.find_by_number(i), node);
          EXPECT_EQ(part.source_node(node), tree.find_by_number(id));
          std::vector<size_t> expected;
          for (auto c : children(id)) {
            if (start[c] < bits) {
              expected.push_back(local[c]);
            }
          }
          std::vector<size_t> represented;
          part.for_each_child(node, [&](size_t child) {
            EXPECT_EQ(part.parent(child), node);
            represented.push_back(child);
          });
          EXPECT_EQ(represented, expected);
          for (size_t j = 0; j < represented.size(); ++j) {
            const auto child = represented[j];
            const bool last = j + 1 == represented.size();
            EXPECT_EQ(part.is_last_child(child), last);
            if (!last) {
              EXPECT_EQ(part.next_sibling(child), represented[j + 1]);
            }
          }
          EXPECT_EQ(part.children_complete(node),
                    end[id] <= bits && expected.size() == children(id).size());
        }
      }
    }
    for (size_t id = 0; id < n; ++id) {
      auto node = tree.find_by_number(id);
      EXPECT_TRUE(index.contains(node));
      EXPECT_TRUE(index.children_complete(node));
      EXPECT_EQ(index.represented_children(node).size(), tree.degree(node));
    }
    auto text = tree.parentheses();
    EXPECT_EQ(text.size(), 2 * n);
    int balance = 0;
    for (char c : text) {
      balance += c == '(' ? 1 : -1;
      EXPECT_GE(balance, 0);
    }
    EXPECT_EQ(balance, 0);
    for (size_t start = 0; start <= text.size() + 1; ++start) {
      EXPECT_EQ(tree.parentheses(start, 7),
                start < text.size() ? text.substr(start, 7) : "");
    }
  }
}

TYPED_TEST(BoundedTraversalTest, FetchOwnsOutputAndBoundsMemory) {
  using Tree = typename TypeParam::Tree;
  const auto detached = [] {
    std::mt19937_64 rng(91);
    const auto adjacency =
        TypeParam::order(128, generate_random_tree(128, rng));
    Tree tree(TypeParam::encode(128, adjacency), 128);
    pixie::tree::TreeTraversalContext disabled(tree, 0);
    EXPECT_TRUE(disabled.fetch(tree.root(), 256).empty());
    constexpr size_t limit = 1024;
    pixie::tree::TreeTraversalContext context(tree, limit);
    auto part = context.fetch(tree.root(), 256);
    EXPECT_FALSE(part.empty());
    EXPECT_LT(part.size(), 128u);
    EXPECT_EQ(&part.source_tree(), &tree);
    EXPECT_EQ(context.memory_usage_bytes(),
              sizeof(context) + context.scan_scratch_bytes());
    EXPECT_LE(part.memory_usage_bytes() + context.memory_usage_bytes(), limit);
    (void)context.fetch(tree.child(tree.root(), 0), 1);
    return part;
  }();
  ASSERT_FALSE(detached.empty());
  EXPECT_EQ(detached.source_node(detached.root()).number, 0u);
  ASSERT_FALSE(detached.is_leaf(detached.root()));
  EXPECT_EQ(detached.first_child(detached.root()), 1u);
}

TYPED_TEST(BoundedTraversalTest, FetchBudgetIsMonotoneAcross4KiB) {
  using Tree = typename TypeParam::Tree;
  constexpr size_t n = 256;
  std::vector<std::vector<size_t>> chain(n);
  chain[0].push_back(0);
  for (size_t i = 1; i < n; ++i) {
    chain[i - 1].push_back(i);
    chain[i].push_back(i - 1);
  }
  const Tree tree(TypeParam::encode(n, TypeParam::order(n, chain)), n);
  size_t previous = 0;
  // Include the former admission threshold, which also charged context bytes.
  for (size_t budget = 4096 - 128; budget <= 4096 + 128; ++budget) {
    SCOPED_TRACE(budget);
    pixie::tree::TreeTraversalContext context(tree, budget);
    const auto part = context.fetch(tree.root(), TypeParam::bit_size(n));
    ASSERT_FALSE(part.empty());
    EXPECT_GE(part.size(), previous);
    EXPECT_LE(context.memory_usage_bytes() + part.memory_usage_bytes(), budget);
    previous = part.size();
  }
  EXPECT_LT(previous, n);
}

TYPED_TEST(BoundedTraversalTest, IndexRejectsEmptyTrees) {
  using Tree = typename TypeParam::Tree;
  const std::vector<uint64_t> words;
  const Tree tree(words, 0);
  for (size_t budget : {size_t{0}, size_t{32768}}) {
    EXPECT_THROW((pixie::tree::RootNavigationIndex<Tree>(tree, budget)),
                 std::invalid_argument);
    EXPECT_THROW((pixie::tree::TreeTraversalContext(tree, budget)),
                 std::invalid_argument);
  }
}

TYPED_TEST(BoundedTraversalTest, CursorContractsAcrossAdmissionBoundaries) {
  using Tree = typename TypeParam::Tree;
  using Index = pixie::tree::RootNavigationIndex<Tree, 8>;
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
      EXPECT_EQ(pixie::tree::detail::tree_unpack::word_at(source, position),
                expected)
          << "words=" << count << " position=" << position;
    }
  }
}

template <class Case>
void CheckBoundedTraversalWorkloads(size_t large_size, size_t trials) {
  using Tree = typename Case::Tree;
  std::mt19937_64 rng(812);
  for (size_t scenario = 0; scenario < 7; ++scenario) {
    const size_t n = scenario == 0 ? 1 : scenario == 1 ? 130 : large_size;
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
    const auto adjacency = Case::order(n, source);
    const auto words = Case::encode(n, adjacency);
    for (size_t budget :
         {size_t(0), size_t(64), size_t(4096), size_t(8192), size_t(16384),
          size_t(32767), size_t(32768), size_t(32769), size_t(65536),
          size_t(81920), size_t(196608), size_t(524288)}) {
      Tree tree(words, n);
      pixie::tree::TreeTraversalContext context(tree, 0);
      pixie::benchmarks::TreeTraversalWorkloads workload(tree, budget);
      EXPECT_LE(workload.reserved_bytes(), budget);
      if (budget >= 32768) {
        EXPECT_TRUE(workload.streaming_enabled());
        EXPECT_TRUE(workload.optimized_inspection_enabled());
      }
      std::vector<std::pair<size_t, size_t>> expected, actual;
      auto collect = [&](auto node) {
        actual.emplace_back(node.number, node.pos);
      };
      auto direct = [&](auto node) {
        expected.emplace_back(node.number, node.pos);
      };
      pixie::tree::detail::tree_unpack::subtree(tree, tree.root(), direct);
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
      pixie::tree::detail::tree_unpack::subtree(tree, tree.root(), ref_events);
      const auto work = workload.inspect(collect, collect);
      EXPECT_EQ(work.first, n);
      EXPECT_EQ(work.second, n - 1);
      EXPECT_EQ(actual, expected);
      for (size_t trial = 0; trial < trials; ++trial) {
        auto cursor = context.cursor();
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

TYPED_TEST(BoundedTraversalTest, AllWorkloadsShareBudgetAndExactEvents) {
  CheckBoundedTraversalWorkloads<TypeParam>(1025, 2);
}

TYPED_TEST(BoundedTraversalTest, HeavyAllWorkloadsShareBudgetAndExactEvents) {
  CheckBoundedTraversalWorkloads<TypeParam>(10001, 32);
}

TYPED_TEST(TreeSpecificationTest,
           ForEachChildBorrowsVisitorAndSupportsNesting) {
  using Tree = typename TypeParam::Tree;
  const std::vector<std::vector<size_t>> source = {
      {0, 1, 2}, {0, 3, 4}, {0}, {1}, {1}};
  const auto adjacency = TypeParam::order(source.size(), source);
  const auto words = TypeParam::encode(source.size(), adjacency);
  const Tree tree(words, source.size());
  const pixie::tree::TreeBase<Tree>& facade = tree;
  const AdjListTree reference(adjacency);
  struct Visitor {
    Visitor() = default;
    Visitor(const Visitor&) = delete;
    std::vector<pixie::tree::TreeNode> nodes;
    void operator()(pixie::tree::TreeNode node) { nodes.push_back(node); }
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
    const pixie::tree::RootNavigationIndex<Tree, 8> index(tree);
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
  using Tree = pixie::tree::LoudsTree<Backend>;
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
      const auto words =
          EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(
              size, adjacency);
      const Tree tree(words, size);
      const AdjListTree reference(adjacency);
      auto view = tree.navigate(tree.root());
      auto expected = reference.root();
      Backend::one_selects = 0;
      EXPECT_EQ(view.parent().number, 0u);
      EXPECT_EQ(Backend::one_selects, 0u);
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
        ASSERT_EQ(Backend::one_selects, checks);
        EXPECT_EQ(view.is_last_child(), last);
        EXPECT_EQ(last, reference.is_last_child(expected));
        // Parent queries reuse the cached edge even on non-last siblings.
        EXPECT_EQ(view.parent().number, reference.parent(expected).number);
        ASSERT_EQ(Backend::one_selects, checks);
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
  const auto words =
      EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(3,
                                                                     adjacency);
  const pixie::tree::LoudsTree<Backend> tree(words, 3);
  const auto first = tree.first_child(tree.root());
  auto a = tree.navigate(first);
  auto b = tree.navigate(first);
  Backend::one_selects = 0;
  EXPECT_FALSE(a.is_last_child());
  auto copy = a;
  EXPECT_EQ(copy.parent().number, 0u);
  EXPECT_EQ(Backend::one_selects, 1u);
  EXPECT_FALSE(b.is_last_child());
  EXPECT_EQ(Backend::one_selects, 2u);
  a.reset(a.next_sibling());
  EXPECT_TRUE(a.is_last_child());
  EXPECT_EQ(Backend::one_selects, 3u);
  EXPECT_FALSE(copy.is_last_child());
  EXPECT_FALSE(b.is_last_child());
  EXPECT_EQ(Backend::one_selects, 3u);
  copy.reset(first);
  EXPECT_EQ(copy.parent().number,
            0u);  // Also supports parent before predicate.
  EXPECT_FALSE(copy.is_last_child());
  EXPECT_EQ(Backend::one_selects, 4u);
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
        view.for_each_child([&](pixie::tree::TreeNode child) {
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
  const auto words =
      EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(size,
                                                                     adjacency);
  const auto check = [&]() {
    const pixie::tree::DFUDSTree<Backend> tree(words, size);
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
  const auto dfuds =
      EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(size,
                                                                     adjacency);
  check(pixie::tree::DFUDSTree<Backend>(dfuds, size), true);
  const auto bp = EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(
      size, adjacency);
  check(pixie::tree::BPTree<Backend>(bp, size), false);
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
  const auto words =
      EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(size,
                                                                     adjacency);
  const pixie::tree::LoudsTree<Backend> tree(words, size);
  auto view = tree.navigate(tree.first_child(tree.root()));
  Backend::one_selects = Backend::zero_selects = 0;
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::zero_selects, 1u);
  EXPECT_EQ(view.next_sibling().number, 2u);
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::zero_selects, 1u);
  EXPECT_EQ(view.first_child().number, 3u);
  EXPECT_EQ(view.child(0).number, 3u);
  EXPECT_EQ(view.first_child().number, 3u);
  EXPECT_EQ(Backend::zero_selects, 2u);
  view.reset(view.node());
  EXPECT_EQ(view.next_sibling().number, 2u);
  EXPECT_EQ(view.degree(), size - 3);
  EXPECT_EQ(Backend::zero_selects, 3u);
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
    void reset(pixie::tree::TreeNode node) { view.reset(node); }
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
  NodeView navigate(pixie::tree::TreeNode node) const {
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
      pixie::tree::TreePathView<ParentCountingTree<Tree>, Capacity> path(
          counted);
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
  pixie::tree::TreePathView<Tree, 1> path(tree);
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
    auto visit =
        [&, token = std::make_unique<int>(1)](pixie::tree::TreeNode current) {
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
    auto stop = [&](pixie::tree::TreeNode) {
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
  pixie::benchmarks::DepthFirst<true>(tree, words, [&](pixie::tree::TreeNode) {
    if (++outer == 1) {
      pixie::benchmarks::DepthFirst<true>(
          tree, words, [&](pixie::tree::TreeNode) { ++inner; });
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
      auto remember = [&](pixie::tree::TreeNode node) {
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
  using Tree = std::conditional_t<BP, pixie::tree::BPTree<RmM>,
                                  pixie::tree::DFUDSTree<RmM>>;
  for (size_t n : {size_t{1}, Bits / 2 - 1, Bits / 2, Bits / 2 + 1}) {
    for (bool star : {false, true}) {
      std::vector<std::vector<size_t>> adjacency(n);
      adjacency[0].push_back(0);
      for (size_t i = 1; i < n; ++i) {
        const auto parent = star ? 0 : i - 1;
        adjacency[i].push_back(parent);
        adjacency[parent].push_back(i);
      }
      auto words =
          BP ? EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(
                   n, adjacency)
             : EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(
                   n, adjacency);
      const size_t bits = 2 * n - (BP ? 0 : 1);
      if (bits % 64) {
        words.back() |= ~uint64_t{0} << (bits % 64);
      }
      const Tree tree(words, n);
      pixie::tree::detail::tree_unpack::BoundedSubtreeCache<Bits> cache;
      const bool loaded =
          cache.template load<BP>(tree, words, bits, tree.root());
      ASSERT_EQ(loaded, n <= Bits / 2);
      ASSERT_EQ(cache.contains(tree.root()), loaded);
      if (loaded) {
        auto check = [&](pixie::tree::TreeNode node) {
          ASSERT_TRUE(cache.contains(node));
          ASSERT_EQ(cache.degree(node), tree.degree(node));
          size_t index = 0;
          cache.children(node, [&](pixie::tree::TreeNode child) {
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
      pixie::tree::detail::tree_unpack::fill_parent_slots(
          values.data() + offset, count, 0x1234);
      for (size_t i = 0; i < values.size(); ++i) {
        EXPECT_EQ(values[i],
                  i >= offset && i < offset + count ? 0x1234 : 0xaaaa);
      }
    }
  }
}

TEST(TreeUnpackingTest, ByteDecoderMatchesBitwiseTraversal) {
  using namespace pixie::tree::detail::tree_unpack;
  for (size_t initial_depth : {8u, 120u}) {
    for (uint64_t initial_number :
         {uint64_t{0}, (uint64_t{1} << 32) - 3, uint64_t{1} << 48}) {
      for (unsigned pattern = 1; pattern < 512; ++pattern) {
        if (!(pattern & 255)) {
          continue;  // All-closing bytes use the traversal's scalar path.
        }
        SCOPED_TRACE(::testing::Message()
                     << pattern << ", depth=" << initial_depth
                     << ", number=" << initial_number);
        const uint64_t base = (uint64_t{1} << 40) + 123;
        std::array<uint64_t, 136> numbers, positions;
        for (size_t i = 0; i < numbers.size(); ++i) {
          numbers[i] = (uint64_t{1} << 45) + i * 13;
          positions[i] = (uint64_t{1} << 46) + i * 19;
        }
        auto expected_numbers = numbers, expected_positions = positions;
        size_t depth = initial_depth, number = initial_number;
        std::vector<std::pair<uint64_t, uint64_t>> expected;
        for (unsigned bit = 0; bit < 8; ++bit) {
          if ((pattern >> bit) & 1) {
            expected_numbers[depth + 8] = number;
            expected_positions[depth + 8] = base + bit;
            expected.emplace_back(number++, base + bit);
            ++depth;
          } else {
            --depth;
            if (!((pattern >> (bit + 1)) & 1)) {
              expected.emplace_back(expected_numbers[depth + 7],
                                    expected_positions[depth + 7]);
            }
          }
        }
        const auto expected_depth = depth, expected_number = number;
        depth = initial_depth;
        number = initial_number;
        alignas(64) std::array<uint64_t, 16> events_n, events_p;
        events_n.fill(UINT64_MAX);
        events_p.fill(UINT64_MAX);
        const auto& program = byte_programs[pattern];
        decode_byte(program, base, number, depth, numbers, positions,
                    events_n.data(), events_p.data());
        ASSERT_EQ(program.count, expected.size());
        ASSERT_EQ(depth, expected_depth);
        ASSERT_EQ(number, expected_number);
        for (size_t i = 0; i < expected.size(); ++i) {
          EXPECT_EQ(events_n[i], expected[i].first);
          EXPECT_EQ(events_p[i], expected[i].second);
        }
        for (size_t i = 0; i < depth + 8; ++i) {
          EXPECT_EQ(numbers[i], expected_numbers[i]);
          EXPECT_EQ(positions[i], expected_positions[i]);
        }
        for (size_t i = 8; i < events_n.size(); ++i) {
          EXPECT_EQ(events_n[i], UINT64_MAX);
          EXPECT_EQ(events_p[i], UINT64_MAX);
        }
      }
    }
  }
}

TEST(TreeUnpackingTest, ZeroOffsetsMatchBoundedBitScan) {
  auto check = [](uint64_t word, unsigned valid, unsigned offset) {
    SCOPED_TRACE(::testing::Message()
                 << word << ", valid=" << valid << ", offset=" << offset);
    std::array<uint8_t, 96> output;
    output.fill(0xaa);
    std::vector<uint8_t> expected;
    for (unsigned bit = 0; bit < valid; ++bit) {
      if (!((word >> bit) & 1)) {
        expected.push_back(bit);
      }
    }
    const auto count = pixie::tree::detail::tree_unpack::zero_offsets(
        word, valid, output.data() + offset);
    ASSERT_EQ(count, expected.size());
    EXPECT_EQ(std::vector<uint8_t>(output.begin() + offset,
                                   output.begin() + offset + count),
              expected);
    for (unsigned i = 0; i < output.size(); ++i) {
      if (i < offset || i >= offset + 64) {
        ASSERT_EQ(output[i], 0xaa);
      }
    }
  };
  for (unsigned bits = 0; bits < 65536; ++bits) {
    for (unsigned shift : {0u, 16u, 32u, 48u}) {
      check(~(uint64_t(bits) << shift), 64, bits % 32);
    }
  }
  std::mt19937_64 rng(42);
  for (unsigned valid = 0; valid <= 64; ++valid) {
    check(0, valid, valid % 32);
    check(UINT64_MAX, valid, valid % 32);
    for (unsigned i = 0; i < 64; ++i) {
      check(rng(), valid, i % 32);
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
// workload semantics. Run the same specification against each layout.
template <class Tree>
class ReferenceTreeSpecification : public testing::Test {};
using ReferenceTrees = testing::Types<AdjListTree,
                                      pixie::references::PointerTree,
                                      pixie::references::BfsCsrTree>;
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
        const auto input =
            preorder ? dfs_order(n, source) : bfs_order(n, source);
        const TypeParam tree(input);
        const auto adj =
            std::is_same_v<TypeParam, pixie::references::BfsCsrTree>
                ? bfs_order(n, input)
                : input;
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

TEST(BfsCsrTreeTest, NumbersPermutedSourceInOrderedBreadthFirstOrder) {
  // BFS visits source IDs 0, 5, 2, 6, 1, 3, 4. The unequal subtrees
  // distinguish sibling order even though the source IDs are not retained.
  const pixie::references::BfsCsrTree tree(
      {{0, 5, 2}, {5}, {0, 3}, {2}, {6}, {0, 6, 1}, {5, 4}});
  const std::vector<std::vector<size_t>> expected = {
      {0, 1, 2}, {0, 3, 4}, {0, 5}, {1, 6}, {1}, {2}, {3}};
  std::vector<decltype(tree.root())> pending{tree.root()};
  for (size_t i = 0; i < expected.size(); ++i) {
    ASSERT_LT(i, pending.size());
    const auto node = pending[i];
    ASSERT_EQ(tree.node_number(node), i);
    EXPECT_EQ(tree.parent(node), expected[i][0]);
    ASSERT_EQ(tree.degree(node), expected[i].size() - 1);
    for (size_t j = 1; j < expected[i].size(); ++j) {
      const auto child = tree.child(node, j - 1);
      EXPECT_EQ(child, expected[i][j]);
      pending.push_back(child);
    }
  }
  EXPECT_EQ(pending.size(), expected.size());
}

TEST(BfsCsrTreeTest, PackedFieldsCrossWordBoundaries) {
  // With seven bits per entry, row 45's child offset (90) and row 100's
  // parent ID (49) each have nonzero bits on both sides of a word boundary.
  constexpr size_t n = 127;
  std::vector<std::vector<size_t>> source(n);
  source[0].push_back(0);
  for (size_t i = 1; i < n; ++i) {
    const size_t parent = (i - 1) / 2;
    source[i].push_back(parent);
    source[parent].push_back(i);
  }
  const pixie::references::BfsCsrTree tree(source);
  for (size_t i = 0; i < n; ++i) {
    EXPECT_EQ(tree.parent(i), i ? (i - 1) / 2 : 0);
    ASSERT_EQ(tree.degree(i), i < n / 2 ? 2u : 0u);
    for (size_t j = 0; j < tree.degree(i); ++j) {
      EXPECT_EQ(tree.child(i, j), 2 * i + 1 + j);
    }
  }
}

TEST(BfsCsrTreeTest, EmptyAndSingleton) {
  using Tree = pixie::references::BfsCsrTree;
  const Tree empty(std::vector<std::vector<size_t>>{});
  EXPECT_TRUE(empty.empty());
  EXPECT_EQ(empty.size(), 0u);
  EXPECT_EQ(empty.root(), empty.size());
  EXPECT_EQ(empty.memory_usage_bytes(), sizeof(Tree));

  const Tree singleton(std::vector<std::vector<size_t>>{{0}});
  EXPECT_FALSE(singleton.empty());
  EXPECT_EQ(singleton.size(), 1u);
  const auto root = singleton.root();
  EXPECT_EQ(root, 0u);
  EXPECT_EQ(singleton.parent(root), root);
  EXPECT_EQ(singleton.degree(root), 0u);
  EXPECT_EQ(singleton.first_child(root), root);
  EXPECT_TRUE(singleton.is_last_child(root));
  EXPECT_EQ(singleton.next_sibling(root), root);
  EXPECT_EQ(singleton.memory_usage_bytes(), sizeof(Tree));
}

TEST(BfsCsrTreeTest, TwoPackedArraysAtWidthTransitions) {
  using Tree = pixie::references::BfsCsrTree;
  constexpr size_t line = pixie::kAlignedStorageLineBytes;
  for (size_t n : {31, 32, 33, 63, 64, 65, 511, 512, 513, 1023, 1024, 1025}) {
    SCOPED_TRACE(n);
    // Stars force the largest child offset but have only zero parent IDs.
    std::vector<std::vector<size_t>> source(n, {0});
    for (size_t i = 1; i < n; ++i) {
      source[0].push_back(i);
    }
    const Tree tree(source);
    ASSERT_EQ(tree.degree(tree.root()), n - 1);
    EXPECT_EQ(tree.child(tree.root(), n - 2), n - 1);
    EXPECT_EQ(tree.parent(n - 1), tree.root());
    EXPECT_TRUE(tree.is_leaf(n - 1));
    const size_t bits = n * std::bit_width(n - 1);
    const size_t padded_bytes = ((bits + 8 * line - 1) / (8 * line)) * line;
    EXPECT_EQ(tree.memory_usage_bytes(), sizeof(Tree) + 2 * padded_bytes);
  }
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
  using Tree = pixie::tree::BPTree<pixie::RmMBTree<>>;
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
  pixie::tree::detail::bp_micro::BpMicroTraversal<2048, Tree, decltype(visit),
                                                  decltype(inspect)>
      small(tree, words, visit, inspect);
  small.walk(tree.root());
  ASSERT_EQ(actual, expected);
  EXPECT_EQ(small.nodes, n);
  EXPECT_EQ(small.children, n - 1);
  actual.clear();
  pixie::tree::detail::bp_micro::BpMicroTraversal<8192, Tree, decltype(visit),
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
        pixie::tree::detail::bp_micro::bp_micro_tables.headers[key] & 15;
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
      auto words = EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(
          n, adjacency);
      if ((2 * n) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n) % 64);
      }
      CheckMicroPrefetch(words, n);
    }
  }
}

namespace {
void CheckLoudsPrefetch(const std::vector<uint64_t>& words, size_t n) {
  using Tree = pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>;
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
      pixie::tree::detail::tree_unpack::LoudsPrefetchTraversal<
          Tree, decltype(ignore), decltype(ignore), 8, true>
          nested(tree, words, ignore, ignore);
      nested.walk(tree.root(), 0);
      EXPECT_EQ(nested.node_count, n);
      EXPECT_EQ(nested.child_count, n - 1);
    }
  };
  pixie::tree::detail::tree_unpack::LoudsPrefetchTraversal<
      Tree, decltype(visit), decltype(inspect), 8, true>
      small(tree, words, visit, inspect);
  small.walk(tree.root(), 0);
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(small.node_count, n);
  ASSERT_EQ(small.child_count, n - 1);
  actual.clear();
  pixie::tree::detail::tree_unpack::LoudsPrefetchTraversal<
      Tree, decltype(visit), decltype(inspect), 32, true>
      large(tree, words, visit, inspect);
  large.walk(tree.root(), 0);
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(large.node_count, n);
  ASSERT_EQ(large.child_count, n - 1);
  actual.clear();
  pixie::tree::detail::tree_unpack::LoudsPrefetchTraversal<
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
      auto words =
          EncodeReference<pixie::tree::detail::TreeEncodingKind::kLOUDS>(n,
                                                                         adj);
      if ((2 * n - 1) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n - 1) % 64);
      }
      CheckLoudsPrefetch(words, n);
    }
  }
}

namespace {
void CheckDfudsMicroPrefetch(const std::vector<uint64_t>& words, size_t n) {
  using Tree = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
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
  pixie::tree::detail::dfuds_micro::DfudsMicroTraversal<
      2048, Tree, decltype(visit), decltype(inspect)>
      small(tree, words, visit, inspect);
  small.walk(tree.root());
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(small.nodes, n);
  ASSERT_EQ(small.children, n - 1);
  actual.clear();
  pixie::tree::detail::dfuds_micro::DfudsMicroTraversal<
      8192, Tree, decltype(visit), decltype(inspect)>
      large(tree, words, visit, inspect);
  large.walk(tree.root());
  ASSERT_EQ(actual, expected);
  ASSERT_EQ(large.nodes, n);
  ASSERT_EQ(large.children, n - 1);
  using Decoder = pixie::tree::detail::dfuds_micro::DfudsMicroSubtree<8192>;
  pixie::tree::experimental::PreparedRootSubtree<Decoder, 8192> prepared(
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
  pixie::tree::experimental::PreparedTreeEvents<> events(tree);
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
        pixie::tree::detail::dfuds_micro::dfuds_micro_tables.headers[key] & 15;
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
      auto words =
          EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(n,
                                                                         adj);
      CheckDfudsMicroPrefetch(words, n);
      if ((2 * n - 1) % 64) {
        words.back() |= ~uint64_t{0} << ((2 * n - 1) % 64);
      }
      CheckDfudsMicroPrefetch(words, n);
    }
  }
}

TEST(TreePreparedEventsTest, EmptyAndNestedReplay) {
  using Tree = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  const std::vector<uint64_t> empty_words;
  const Tree empty(empty_words, 0);
  const pixie::tree::experimental::PreparedTreeEvents<> empty_cache(empty);
  size_t outer_visits = 0, outer_children = 0, inner_visits = 0,
         inner_children = 0;
  auto inner_visit = [&](auto) { ++inner_visits; };
  auto inner_inspect = [&](auto) { ++inner_children; };
  EXPECT_FALSE(empty_cache.try_walk(inner_visit, inner_inspect));
  EXPECT_EQ(inner_visits + inner_children, 0u);
  const std::vector<uint64_t> words{5};  // Three-node chain: 10 10 0.
  const Tree tree(words, 3);
  const pixie::tree::experimental::PreparedTreeEvents<> cache(tree);
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

namespace {

using SimdTestEvent = std::tuple<size_t, size_t, int>;

template <class Tree>
std::vector<SimdTestEvent> CollectExpectedSimdEvents(
    const Tree& tree,
    pixie::tree::TreeNode root) {
  std::vector<SimdTestEvent> expected;
  if (tree.empty()) {
    return expected;
  }
  auto current = root;
  bool entering = true;
  while (true) {
    expected.emplace_back(current.number, current.pos, entering ? 1 : 0);
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
    if (current.number == root.number) {
      break;
    }
    if (tree.is_last_child(current)) {
      current = tree.parent(current);
    } else {
      current = tree.next_sibling(current);
      entering = true;
    }
  }
  return expected;
}

template <template <size_t, class, class, class> class Walk,
          size_t MaxSimdBits,
          class Tree>
void CheckSimdWalk(const Tree& tree,
                   const std::vector<uint64_t>& words,
                   pixie::tree::TreeNode root,
                   size_t nodes,
                   size_t max_depth = 128) {
  const auto expected = CollectExpectedSimdEvents(tree, root);
  std::vector<SimdTestEvent> actual;
  auto act_visit = [&, token = std::make_unique<int>(1)](
                       pixie::tree::TreeNode node, bool entering) {
    actual.emplace_back(node.number, node.pos, entering ? 1 : 0);
  };
  auto act_inspect = [&](pixie::tree::TreeNode node) {
    actual.emplace_back(node.number, node.pos, 2);
  };
  Walk<MaxSimdBits, Tree, decltype(act_visit), decltype(act_inspect)> simd(
      tree, words, act_visit, act_inspect, max_depth);
  simd.walk(root);
  ASSERT_EQ(simd.nodes, nodes);
  ASSERT_EQ(simd.children, nodes > 0 ? nodes - 1 : 0u);
  ASSERT_EQ(actual, expected)
      << "subtree nodes=" << nodes << " root=" << root.number
      << " MaxSimdBits=" << MaxSimdBits << " max_depth=" << max_depth;
}

template <size_t MaxSimdBits>
void CheckBpSimdEquivalence(const std::vector<uint64_t>& words,
                            size_t n,
                            size_t max_depth = 128) {
  using Tree = pixie::tree::BPTree<pixie::RmMBTree<>>;
  Tree tree(words, n);
  const auto root = n ? tree.root() : pixie::tree::TreeNode(0, 0);
  CheckSimdWalk<pixie::tree::experimental::BpSimdMicroTraversal, MaxSimdBits>(
      tree, words, root, n, max_depth);
}

template <size_t MaxSimdBits>
void CheckDfudsSimdEquivalence(const std::vector<uint64_t>& words,
                               size_t n,
                               size_t max_depth = 128) {
  using Tree = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
  Tree tree(words, n);
  const auto root = n ? tree.root() : pixie::tree::TreeNode(0, 0);
  CheckSimdWalk<pixie::tree::experimental::DfudsSimdMicroTraversal,
                MaxSimdBits>(tree, words, root, n, max_depth);
}

}  // namespace

TEST(TreeSimdMicroPrimitivesTest, LowLevelKernelsMatchScalarReference) {
  using namespace pixie::tree::experimental::microtree_simd;
  std::mt19937_64 rng(20261111);

  for (int trial = 0; trial < 2000; ++trial) {
    uint64_t s[2] = {rng(), rng()};
    if (trial == 0) {
      s[0] = s[1] = 0;
    } else if (trial == 1) {
      s[0] = s[1] = ~uint64_t{0};
    } else if (trial == 2) {
      s[0] = s[1] = 0x5555555555555555ULL;
    }

    // Check bp_pair_zero_mask_64
    uint32_t exp64 = 0;
    int h = 0;
    for (unsigned j = 0; j < 32; ++j) {
      const unsigned p = static_cast<unsigned>((s[0] >> (2 * j)) & 3u);
      h += static_cast<int>(p & 1u) + static_cast<int>((p >> 1) & 1u) - 1;
      if (h == 0) {
        exp64 |= (1u << j);
      }
    }
    ASSERT_EQ(bp_pair_zero_mask_64(s[0]), exp64);

    // Check bp_pair_zero_mask_128 and bp_children_128 across initial_h
    for (int init_h : {-70, -64, -5, -1, 0, 1, 5, 32, 64, 70}) {
      uint64_t exp128 = 0;
      uint64_t exp_child = 0;
      size_t exp_close = 0;
      bool exp_closed = false;
      int cur_h = init_h;
      for (unsigned w = 0; w < 2; ++w) {
        for (unsigned j = 0; j < 32; ++j) {
          const size_t idx = 32 * w + j;
          const unsigned p = static_cast<unsigned>((s[w] >> (2 * j)) & 3u);
          if (!exp_closed && cur_h == 0 && (p & 1u) == 0) {
            exp_closed = true;
            exp_close = idx;
            exp_child = exp128;
          }
          cur_h +=
              static_cast<int>(p & 1u) + static_cast<int>((p >> 1) & 1u) - 1;
          if (cur_h == 0) {
            exp128 |= (1ULL << idx);
          }
        }
      }
      int next_h = 0;
      ASSERT_EQ(bp_pair_zero_mask_128(s, init_h, &next_h), exp128);
      ASSERT_EQ(next_h, cur_h);

      uint64_t act_child = 0;
      size_t act_close = 999;
      int act_next_h = 0;
      const bool act_closed =
          bp_children_128(s, init_h, act_child, act_close, act_next_h);
      ASSERT_EQ(act_closed, exp_closed);
      ASSERT_EQ(act_next_h, cur_h);
      if (exp_closed) {
        ASSERT_EQ(act_close, exp_close);
        ASSERT_EQ(act_child, exp_child);
      } else {
        ASSERT_EQ(act_child, exp128);
      }
    }

    // Check dfuds_record_lows_64 and dfuds_record_lows_128
    uint64_t exp_rl64 = 0;
    int cur = 0, best = 0;
    for (unsigned b = 0; b < 64; ++b) {
      cur += ((s[0] >> b) & 1ULL) ? 1 : -1;
      if (cur < best) {
        best = cur;
        exp_rl64 |= (1ULL << b);
      }
    }
    ASSERT_EQ(dfuds_record_lows_64(s[0]), exp_rl64);

    for (int init_gap : {0, 1, 2, 3, 4, 15, 63, 127, 128, 200}) {
      uint64_t exp_lo = 0, exp_hi = 0;
      cur = 0;
      best = -init_gap;
      for (unsigned b = 0; b < 64; ++b) {
        cur += ((s[0] >> b) & 1ULL) ? 1 : -1;
        if (cur < best) {
          best = cur;
          exp_lo |= (1ULL << b);
        }
      }
      for (unsigned b = 0; b < 64; ++b) {
        cur += ((s[1] >> b) & 1ULL) ? 1 : -1;
        if (cur < best) {
          best = cur;
          exp_hi |= (1ULL << b);
        }
      }
      int next_gap = -1;
      const Mask128 act = dfuds_record_lows_128(s, init_gap, &next_gap);
      ASSERT_EQ(act.lo, exp_lo) << "init_gap=" << init_gap;
      ASSERT_EQ(act.hi, exp_hi) << "init_gap=" << init_gap;
      ASSERT_EQ(next_gap, cur - best);
    }
  }
}

TEST(TreeSimdMicroPrimitivesTest, LoadsZeroExtendAcrossSpanEnd) {
  using namespace pixie::tree::experimental::microtree_simd;
  const uint64_t words[]{~uint64_t{}, 3};
  for (size_t pos : {63, 64, 127, 128, 129}) {
    uint64_t want[2]{};
    for (size_t i = 0; pos + i < 128; ++i) {
      want[i / 64] |= ((words[(pos + i) / 64] >> ((pos + i) % 64)) & 1)
                      << (i % 64);
    }
    uint64_t got[]{~uint64_t{}, ~uint64_t{}};
    load_bits_128(words, pos, got);
    EXPECT_EQ(load_bits_64(words, pos), want[0]);
    EXPECT_EQ(got[0], want[0]);
    EXPECT_EQ(got[1], want[1]);
  }
}

TEST(TreeSimdMicroTraversalTest, ExhaustiveAndBoundaryEquivalence) {
  CheckBpSimdEquivalence<512>({}, 0);
  CheckDfudsSimdEquivalence<512>({}, 0);

  // Exhaustive small shapes from 8192-entry microtree tables
  for (unsigned key = 0; key < 8192; ++key) {
    const auto bp_len =
        pixie::tree::detail::bp_micro::bp_micro_tables.headers[key] & 15;
    if (bp_len) {
      const uint64_t raw = (uint64_t(key) << 1) | 1ULL;
      if (raw < (1ULL << bp_len)) {
        CheckBpSimdEquivalence<128>({raw}, bp_len / 2);
        CheckBpSimdEquivalence<512>({raw | (~uint64_t{0} << bp_len)},
                                    bp_len / 2);
      }
    }
    const auto df_len =
        pixie::tree::detail::dfuds_micro::dfuds_micro_tables.headers[key] & 15;
    if (df_len && key < (1u << df_len)) {
      CheckDfudsSimdEquivalence<128>({uint64_t(key)}, (df_len + 1) / 2);
      CheckDfudsSimdEquivalence<512>({uint64_t(key) | (~uint64_t{0} << df_len)},
                                     (df_len + 1) / 2);
    }
  }

  // Random, star, chain, and complete binary shapes across all SIMD boundaries
  for (size_t n :
       {1,   2,   3,   4,   5,   31,  32,  33,  34,  63,   64,   65,  66,
        127, 128, 129, 255, 256, 257, 258, 512, 513, 1024, 4096, 4097}) {
    for (int shape = 0; shape < 4; ++shape) {
      std::mt19937_64 rng(99 + n * 7 + shape);
      auto adj = generate_random_tree(n, rng);
      if (shape) {
        adj.assign(n, {});
        adj[0].push_back(0);
        for (size_t i = 1; i < n; ++i) {
          const size_t p = shape == 1 ? 0 : shape == 2 ? i - 1 : (i - 1) / 2;
          adj[p].push_back(i);
          adj[i].push_back(p);
        }
      }
      adj = dfs_order(n, adj);
      auto bp_words =
          EncodeReference<pixie::tree::detail::TreeEncodingKind::kBP>(n, adj);
      auto df_words =
          EncodeReference<pixie::tree::detail::TreeEncodingKind::kDFUDS>(n,
                                                                         adj);

      CheckBpSimdEquivalence<128>(bp_words, n, 4);
      CheckBpSimdEquivalence<512>(bp_words, n, 128);
      CheckBpSimdEquivalence<1024>(bp_words, n, 128);

      CheckDfudsSimdEquivalence<128>(df_words, n, 4);
      CheckDfudsSimdEquivalence<512>(df_words, n, 128);
      CheckDfudsSimdEquivalence<1024>(df_words, n, 128);

      // Dirty trailing bits past logical encoding end
      if ((2 * n) % 64) {
        bp_words.back() |= ~uint64_t{0} << ((2 * n) % 64);
        CheckBpSimdEquivalence<512>(bp_words, n, 128);
      }
      if ((2 * n - 1) % 64) {
        df_words.back() |= ~uint64_t{0} << ((2 * n - 1) % 64);
        CheckDfudsSimdEquivalence<512>(df_words, n, 128);
      }
    }
  }
}

TEST(TreeSimdMicroTraversalTest, NonrootMixedFanoutAcrossChildBuffer) {
  // An 80-child internal node with non-leaf children forces the overflow walk
  // rather than the all-leaf star shortcut. Its sibling bounds the subtree.
  std::vector<std::vector<size_t>> adj{{0, 1, 2}, {0}, {0}};
  for (size_t j = 0; j < 80; ++j) {
    const size_t child = adj.size();
    adj[1].push_back(child);
    adj.push_back({1, child + 1});
    adj.push_back({child});
  }
  adj = dfs_order(adj.size(), adj);
  auto check = [&]<template <size_t, class, class, class> class Walk>(
                   const auto& tree, const auto& words) {
    const auto inner = tree.first_child(tree.root());
    const auto leaf = tree.child(tree.root(), 1);
    CheckSimdWalk<Walk, 128>(tree, words, inner, 161);
    CheckSimdWalk<Walk, 128>(tree, words, leaf, 1);
    CheckSimdWalk<Walk, 512>(tree, words, inner, 161);
    CheckSimdWalk<Walk, 512>(tree, words, leaf, 1);
    CheckSimdWalk<Walk, 1024>(tree, words, inner, 161);
    CheckSimdWalk<Walk, 1024>(tree, words, leaf, 1);
  };
  using namespace pixie::tree;
  auto bp = EncodeReference<detail::TreeEncodingKind::kBP>(adj.size(), adj);
  check.template operator()<experimental::BpSimdMicroTraversal>(
      BPTree<pixie::RmMBTree<>>(bp, adj.size()), bp);
  auto df = EncodeReference<detail::TreeEncodingKind::kDFUDS>(adj.size(), adj);
  check.template operator()<experimental::DfudsSimdMicroTraversal>(
      DFUDSTree<pixie::RmMBTree<>>(df, adj.size()), df);
}

namespace {

// Independent ordered source: its integer handles are deliberately not assumed
// to be preorder numbers, and it supplies no packed encoding to the
// constructor.
struct DecompositionSource {
  std::vector<std::vector<size_t>> children;
  std::vector<size_t> parents;
  std::vector<size_t> next;

  explicit DecompositionSource(const std::vector<size_t>& parent)
      : children(parent.size()), parents(parent), next(parent.size()) {
    for (size_t v = 1; v < parent.size(); ++v) {
      auto& siblings = children[parent[v]];
      if (!siblings.empty()) {
        next[siblings.back()] = v;
      }
      siblings.push_back(v);
    }
  }
  size_t size() const { return children.size(); }
  size_t root() const { return 0; }
  size_t degree(size_t v) const { return children[v].size(); }
  size_t first_child(size_t v) const { return children[v].front(); }
  size_t next_sibling(size_t v) const { return next[v]; }
};

std::vector<size_t> DecompositionParents(size_t n, int shape, uint64_t seed) {
  std::vector<size_t> parent(n);
  std::mt19937_64 rng(seed);
  for (size_t v = 1; v < n; ++v) {
    switch (shape) {
      case 0:  // Random recursive, not uniform Catalan.
        parent[v] = rng() % v;
        break;
      case 1:
        parent[v] = v - 1;  // Chain.
        break;
      case 2:
        parent[v] = 0;  // Star.
        break;
      case 3:
        parent[v] = (v - 1) / 2;  // Balanced binary.
        break;
      case 4:
        parent[v] = v < 3 ? 0 : ((v - 1) / 2) * 2;  // Comb.
        break;
      default:
        parent[v] = (v - 1) / 257;  // High degree at multiple depths.
        break;
    }
  }
  return parent;
}

template <class Tree>
void CheckDecomposition(const std::vector<size_t>& parents) {
  const DecompositionSource source(parents);
  const Tree tree(source);
  const size_t n = source.size();
  ASSERT_EQ(tree.size(), n);
  ASSERT_EQ(tree.empty(), n == 0);
  ASSERT_EQ(tree.topology_bits(), 2 * n);
  ASSERT_TRUE(tree.validate_invariants());
  EXPECT_GE(tree.memory_usage_bytes(), sizeof(Tree) + (2 * n + 7) / 8);
  const auto memory = tree.memory_breakdown();
  EXPECT_EQ(memory.total_bytes, tree.memory_usage_bytes());
  EXPECT_EQ(memory.total_bytes,
            memory.object_bytes + memory.bottom_topology_bytes +
                memory.component_directory_bytes + memory.boundary_bytes +
                memory.external_parent_bytes + memory.child_directory_bytes +
                memory.preorder_bytes + memory.macro_bytes);
  size_t covered = 0, exits = 0;
  for (size_t c = 0; c < tree.component_count(); ++c) {
    const auto nodes = tree.component_size(c);
    const auto cuts = tree.component_exit_count(c);
    ASSERT_GE(nodes, 1);
    ASSERT_LE(nodes, 256);
    ASSERT_LE(cuts, nodes == 1 ? 2 : 1);
    covered += nodes;
    exits += cuts;
  }
  EXPECT_EQ(covered, n);
  EXPECT_EQ(exits, tree.exit_count());
  EXPECT_LE(tree.preorder_run_count(), 2 * tree.component_count());
  EXPECT_THROW(tree.find_by_preorder(n), std::out_of_range);
  EXPECT_THROW(tree.find_by_preorder(std::numeric_limits<size_t>::max()),
               std::out_of_range);
  if (!n) {
    EXPECT_EQ(tree.component_count(), 0);
    EXPECT_THROW(tree.root(), std::out_of_range);
    return;
  }
  EXPECT_EQ(exits, tree.component_count() - 1);

  std::vector<size_t> order, pending{0}, rank(n);
  while (!pending.empty()) {
    const auto v = pending.back();
    pending.pop_back();
    rank[v] = order.size();
    order.push_back(v);
    for (auto it = source.children[v].rbegin(); it != source.children[v].rend();
         ++it) {
      pending.push_back(*it);
    }
  }
  ASSERT_EQ(order.size(), n);
  EXPECT_EQ(tree.find_by_preorder(0), tree.root());
  auto view = tree.navigate(tree.root());
  for (size_t p = 0; p < n; ++p) {
    const auto v = order[p];
    const auto node = tree.find_by_preorder(p);
    ASSERT_EQ(tree.preorder_number(node), p);
    view.reset(node);
    ASSERT_EQ(view.node(), node);
    ASSERT_EQ(view.is_root(), p == 0);
    ASSERT_EQ(view.is_leaf(), source.children[v].empty());
    ASSERT_EQ(view.degree(), source.children[v].size());
    if (n < 300) {
      EXPECT_THROW(view.child(view.degree()), std::out_of_range);
      if (view.is_leaf()) {
        EXPECT_THROW(view.first_child(), std::out_of_range);
      }
    }
    ASSERT_EQ(tree.preorder_number(view.parent()), rank[source.parents[v]]);
    if (!source.children[v].empty()) {
      ASSERT_EQ(view.first_child(), view.child(0));
    }
    size_t i = 0;
    view.for_each_child([&](auto child) {
      ASSERT_LT(i, source.children[v].size());
      ASSERT_EQ(tree.preorder_number(child), rank[source.children[v][i]]);
      ASSERT_EQ(view.child(i), child);
      auto child_view = tree.navigate(child);
      ASSERT_EQ(child_view.parent(), node);
      ASSERT_EQ(child_view.is_last_child(), i + 1 == source.children[v].size());
      if (i + 1 < source.children[v].size()) {
        ASSERT_EQ(tree.preorder_number(child_view.next_sibling()),
                  rank[source.children[v][i + 1]]);
      }
      ++i;
    });
    ASSERT_EQ(i, source.children[v].size());
    if (n < 300 || p % 997 == 0) {
      auto independent = view;
      independent.reset(tree.find_by_preorder(n - 1));
      EXPECT_EQ(independent.degree(), source.children[order.back()].size());
      EXPECT_EQ(view.node(), node);
      EXPECT_EQ(view.degree(), source.children[v].size());
      EXPECT_EQ(independent.node(), tree.find_by_preorder(n - 1));
    }
  }

  // Ordered enter/inspect/leave events, including leaf exits. Neither expected
  // events nor their IDs use a Pixie encoding or navigation implementation.
  using Event = std::pair<unsigned, size_t>;
  std::vector<Event> expected, actual;
  std::vector<std::pair<size_t, bool>> source_tasks{{0, false}};
  while (!source_tasks.empty()) {
    const auto [v, leaving] = source_tasks.back();
    source_tasks.pop_back();
    if (leaving) {
      expected.emplace_back(2, rank[v]);
      continue;
    }
    expected.emplace_back(0, rank[v]);
    for (auto child : source.children[v]) {
      expected.emplace_back(1, rank[child]);
    }
    source_tasks.emplace_back(v, true);
    for (auto it = source.children[v].rbegin(); it != source.children[v].rend();
         ++it) {
      source_tasks.emplace_back(*it, false);
    }
  }
  using Node = decltype(tree.root());
  std::vector<std::pair<Node, bool>> tasks{{tree.root(), false}};
  while (!tasks.empty()) {
    const auto [node, leaving] = tasks.back();
    tasks.pop_back();
    const auto p = tree.preorder_number(node);
    if (leaving) {
      actual.emplace_back(2, p);
      continue;
    }
    actual.emplace_back(0, p);
    view.reset(node);
    std::vector<Node> children;
    view.for_each_child([&](auto child) {
      actual.emplace_back(1, tree.preorder_number(child));
      children.push_back(child);
    });
    tasks.emplace_back(node, true);
    for (auto it = children.rbegin(); it != children.rend(); ++it) {
      tasks.emplace_back(*it, false);
    }
  }
  EXPECT_EQ(actual, expected);
}

template <class Tree>
class DecomposedTreeTest : public ::testing::Test {};
using DecomposedTreeTypes =
    ::testing::Types<pixie::tree::experimental::BTreeDecomposedTree<2>,
                     pixie::tree::experimental::MacroDecomposedTree<2>,
                     pixie::tree::experimental::BTreeDecomposedTree<3>,
                     pixie::tree::experimental::MacroDecomposedTree<3>,
                     pixie::tree::experimental::BTreeDecomposedTree<>,
                     pixie::tree::experimental::MacroDecomposedTree<>>;
TYPED_TEST_SUITE(DecomposedTreeTest, DecomposedTreeTypes);

TYPED_TEST(DecomposedTreeTest, ExhaustiveSmallOrderedShapes) {
  CheckDecomposition<TypeParam>({});
  size_t shapes = 0;
  for (size_t n = 1; n <= 9; ++n) {
    for (uint64_t bits = 1; bits < (uint64_t{1} << (2 * n)); bits += 2) {
      std::vector<size_t> parents, stack;
      bool valid = true;
      for (size_t b = 0; b < 2 * n; ++b) {
        if ((bits >> b) & 1) {
          if (stack.empty() && !parents.empty()) {
            valid = false;
            break;
          }
          const size_t v = parents.size();
          parents.push_back(stack.empty() ? 0 : stack.back());
          stack.push_back(v);
        } else {
          if (stack.empty()) {
            valid = false;
            break;
          }
          stack.pop_back();
        }
      }
      if (valid && stack.empty() && parents.size() == n) {
        CheckDecomposition<TypeParam>(parents);
        ++shapes;
      }
    }
  }
  EXPECT_EQ(shapes, 2056);
}

TYPED_TEST(DecomposedTreeTest, BoundariesAndNonPreorderSources) {
  for (size_t n : {1, 2, 3, 31, 63, 127, 128, 129, 255, 256, 257, 511, 1025}) {
    for (int shape = 0; shape < 6; ++shape) {
      SCOPED_TRACE(::testing::Message() << "n=" << n << " shape=" << shape);
      CheckDecomposition<TypeParam>(DecompositionParents(n, shape, 313 + n));
    }
  }
  for (size_t seed = 0; seed < 20; ++seed) {
    CheckDecomposition<TypeParam>(DecompositionParents(513, 0, seed));
  }
}

TYPED_TEST(DecomposedTreeTest, IndependentContextsAndReentrantEnumeration) {
  const TypeParam tree(DecompositionSource(DecompositionParents(1025, 3, 1)));
  auto view = tree.navigate(tree.root());
  size_t children = 0;
  view.for_each_child([&](auto child) {
    auto nested = tree.navigate(child);
    nested.for_each_child([&](auto grandchild) {
      EXPECT_EQ(tree.navigate(grandchild).parent(), child);
      EXPECT_EQ(view.node(), tree.root());
    });
    ++children;
  });
  EXPECT_EQ(children, 2);
  EXPECT_THROW(
      view.for_each_child([](auto) { throw std::runtime_error("stop"); }),
      std::runtime_error);
  EXPECT_EQ(view.degree(), 2);
  auto check = [&tree] {
    size_t sum = 0;
    auto cursor = tree.navigate(tree.root());
    for (size_t p = 0; p < tree.size(); ++p) {
      cursor.reset(tree.find_by_preorder(p));
      sum += cursor.degree();
      if (tree.preorder_number(cursor.node()) != p) {
        return size_t{0};
      }
    }
    return sum;
  };
  auto a = std::async(std::launch::async, check);
  auto b = std::async(std::launch::async, check);
  EXPECT_EQ(a.get(), tree.size() - 1);
  EXPECT_EQ(b.get(), tree.size() - 1);
}

TYPED_TEST(DecomposedTreeTest, OwningCopiesSurviveSourceAndOriginal) {
  auto make_copy = [] {
    const TypeParam original(
        DecompositionSource(DecompositionParents(2049, 0, 19)));
    TypeParam copy = original;
    return copy;
  };
  TypeParam copy = make_copy();
  TypeParam moved = std::move(copy);
  ASSERT_TRUE(moved.validate_invariants());
  size_t edges = 0;
  auto cursor = moved.navigate(moved.root());
  for (size_t p = 0; p < moved.size(); ++p) {
    cursor.reset(moved.find_by_preorder(p));
    ASSERT_EQ(moved.preorder_number(cursor.node()), p);
    edges += cursor.degree();
  }
  EXPECT_EQ(edges, moved.size() - 1);
}

TEST(DecomposedTreeLargeTest, HeavyDeepAndWideWithoutRecursiveStack) {
  for (size_t n : {size_t{65535}, size_t{65536}, size_t{65537}}) {
    for (int shape : {1, 2}) {
      SCOPED_TRACE(::testing::Message() << "n=" << n << " shape=" << shape);
      const auto parents = DecompositionParents(n, shape, 1);
      CheckDecomposition<pixie::tree::experimental::BTreeDecomposedTree<>>(
          parents);
      CheckDecomposition<pixie::tree::experimental::MacroDecomposedTree<>>(
          parents);
    }
  }
}

TEST(DecomposedTreePrimitiveTest, FullByteDomainAndUnalignedForests) {
  using pixie::tree::experimental::decomposition_detail::DecodedForest;
  DecodedForest decoded;
  for (size_t offset = 0; offset < 128; ++offset) {
    for (int shape = 0; shape < 3; ++shape) {
      std::vector<uint64_t> words((offset + 512 + 63) / 64, ~uint64_t{0});
      for (size_t b = 0; b < 512; ++b) {
        const bool open = shape == 0   ? b % 2 == 0
                          : shape == 1 ? b < 256
                                       : b == 0 || (b < 511 && b % 2 == 1);
        if (!open) {
          words[(offset + b) / 64] &= ~(uint64_t{1} << ((offset + b) % 64));
        }
      }
      decoded.decode(words, offset, 256);
      ASSERT_EQ(decoded.size, 256);
      for (size_t i = 0; i < 256; ++i) {
        const bool has_parent = shape != 0 && i != 0;
        const bool has_child = shape == 1 ? i < 255 : shape == 2 && i == 0;
        const bool has_next =
            shape == 0 ? i < 255 : shape == 2 && i > 0 && i < 255;
        ASSERT_EQ(decoded.has_parent[i], has_parent);
        ASSERT_EQ(decoded.has_first_child[i], has_child);
        ASSERT_EQ(decoded.has_next_sibling[i], has_next);
        if (has_parent) {
          ASSERT_EQ(decoded.parent[i], shape == 1 ? i - 1 : 0);
        }
        if (has_child) {
          ASSERT_EQ(decoded.first_child[i], i + 1);
        }
        if (has_next) {
          ASSERT_EQ(decoded.next_sibling[i], i + 1);
        }
      }
      words.pop_back();
      EXPECT_THROW(decoded.decode(words, offset, 256), std::invalid_argument);
    }
  }
  decoded.decode({}, 0, 0);
  EXPECT_EQ(decoded.size, 0);
  EXPECT_TRUE(decoded.has_parent.none());
  EXPECT_TRUE(decoded.has_first_child.none());
  EXPECT_TRUE(decoded.has_next_sibling.none());
  EXPECT_THROW(decoded.decode({}, 0, 257), std::invalid_argument);
  EXPECT_THROW(decoded.decode({}, std::numeric_limits<size_t>::max(), 1),
               std::invalid_argument);
}

TEST(DecomposedTreePrimitiveTest, FullBottomGroupsAndWideRoutingBases) {
  namespace detail = pixie::tree::experimental::decomposition_detail;
  std::vector<detail::PreorderRun> runs(65537);
  for (size_t i = 0; i < runs.size(); ++i) {
    runs[i] = {i, 256, 0};
  }
  detail::BalancedRouting balanced;
  balanced.build(runs, {}, {});
  detail::RunDirectory sampled;
  sampled.build(runs);
  ASSERT_GT(balanced.height(), 2);
  for (size_t run : {size_t{0}, size_t{254}, size_t{255}, size_t{256},
                     size_t{65535}, size_t{65536}}) {
    for (size_t local : {size_t{0}, size_t{1}, size_t{255}}) {
      const auto rank = run * 256 + local;
      EXPECT_EQ(balanced.find(rank), run);
      EXPECT_EQ(sampled.find(rank), run);
      EXPECT_EQ(balanced.start(run), run * 256);
    }
  }
  // Exercise full-width directory comparisons without a multi-billion-node
  // fixture. The compact bottom offsets remain unchanged.
  const uint64_t base = uint64_t{1} << 40;
  for (auto& value : balanced.bottom.samples) {
    value += base;
  }
  for (auto& level : balanced.levels) {
    for (auto& value : level) {
      value += base;
    }
  }
  for (size_t run : {size_t{0}, size_t{255}, size_t{256}, size_t{65536}}) {
    EXPECT_EQ(balanced.find(base + run * 256 + 255), run);
    EXPECT_EQ(balanced.start(run), base + run * 256);
  }
  EXPECT_EQ(detail::checked_add(base, 65536), base + 65536);
  EXPECT_THROW(detail::checked_add(std::numeric_limits<size_t>::max(), 1),
               std::length_error);
}

TYPED_TEST(DecomposedTreeTest, RejectsInconsistentSourceSize) {
  struct Source : DecompositionSource {
    using DecompositionSource::DecompositionSource;
    size_t reported = 1;
    size_t size() const { return reported; }
  } source({0, 0});
  EXPECT_THROW((TypeParam{source}), std::invalid_argument);
  source.reported = 3;
  EXPECT_THROW((TypeParam{source}), std::invalid_argument);
}

TYPED_TEST(DecomposedTreeTest, MissingSiblingLeavesViewUsable) {
  const TypeParam tree(DecompositionSource({0, 0}));
  auto view = tree.navigate(tree.root());
  EXPECT_THROW(view.next_sibling(), std::out_of_range);
  view.reset(view.first_child());
  EXPECT_THROW(view.next_sibling(), std::out_of_range);
  EXPECT_EQ(view.parent(), tree.root());
  const auto leaf = view.node();
  EXPECT_THROW(view.reset({tree.component_count(), 0}), std::out_of_range);
  EXPECT_THROW(view.reset({leaf.component, 255}), std::out_of_range);
  EXPECT_EQ(view.node(), leaf);
}

TYPED_TEST(DecomposedTreeTest, RejectsUnrepresentableTopologyBeforeTraversal) {
  struct OversizedSource {
    size_t size() const { return std::numeric_limits<size_t>::max(); }
    size_t root() const { throw std::runtime_error("must not traverse"); }
    size_t degree(size_t) const { return 0; }
    size_t first_child(size_t) const { return 0; }
    size_t next_sibling(size_t) const { return 0; }
  };
  EXPECT_THROW((TypeParam{OversizedSource{}}), std::length_error);
}

TEST(BfsCsrTreeTest, RejectsOverflowBeforeTraversal) {
  struct Source : DecompositionSource {
    using DecompositionSource::DecompositionSource;
    size_t size() const { return std::numeric_limits<size_t>::max(); }
    size_t root() const { throw std::runtime_error("unexpected traversal"); }
  } source(std::vector<size_t>{});
  EXPECT_THROW((pixie::references::BfsCsrTree{source}), std::length_error);
}

}  // namespace

namespace {

// Expected identities come from the source adjacency, never packed-tree
// queries. Filtering original preorder by actual depth bands gives layer-major
// preorder.
struct LayeredTreeOracle {
  explicit LayeredTreeOracle(const std::vector<size_t>& parents)
      : source(parents), depth(parents.size()), number(parents.size()) {
    std::vector<size_t> pending;
    if (!parents.empty()) {
      pending.push_back(0);
    }
    while (!pending.empty()) {
      const size_t v = pending.back();
      pending.pop_back();
      preorder.push_back(v);
      for (auto it = source.children[v].rbegin();
           it != source.children[v].rend(); ++it) {
        depth[*it] = depth[v] + 1;
        pending.push_back(*it);
      }
    }
  }

  void SetLayers(const std::vector<size_t>& depth_begins) {
    order.clear();
    std::vector<std::vector<size_t>> layers(depth_begins.size());
    for (size_t v : preorder) {
      const auto end =
          std::upper_bound(depth_begins.begin(), depth_begins.end(), depth[v]);
      layers[static_cast<size_t>(end - depth_begins.begin()) - 1].push_back(v);
    }
    for (const auto& layer : layers) {
      for (size_t v : layer) {
        number[v] = order.size();
        order.push_back(v);
      }
    }
  }

  DecompositionSource source;
  std::vector<size_t> depth;
  std::vector<size_t> preorder;
  std::vector<size_t> order;
  std::vector<size_t> number;
};

using LayeredTree = pixie::tree::experimental::LayeredDFUDSTree;

template <class Owner, bool Expected>
constexpr bool LayeredWordAccessMatches =
    (requires { std::declval<Owner>().forest_words(); } == Expected) &&
    (requires { std::declval<Owner>().connector_words(); } == Expected) &&
    (requires { std::declval<Owner>().boundary_words(); } == Expected);
static_assert(LayeredWordAccessMatches<LayeredTree&, true>);
static_assert(LayeredWordAccessMatches<const LayeredTree&, true>);
static_assert(LayeredWordAccessMatches<LayeredTree&&, false>);
static_assert(LayeredWordAccessMatches<const LayeredTree&&, false>);

void CheckLayeredNavigation(const LayeredTree& tree,
                            const LayeredTreeOracle& oracle) {
  const auto& source = oracle.source;
  const size_t n = source.size();
  ASSERT_EQ(tree.size(), n);
  ASSERT_EQ(tree.empty(), n == 0);
  ASSERT_EQ(oracle.order.size(), n);
  EXPECT_THROW(tree.find_by_number(n), std::out_of_range);
  EXPECT_THROW(tree.find_by_number(std::numeric_limits<size_t>::max()),
               std::out_of_range);
  if (!n) {
    EXPECT_THROW(tree.root(), std::out_of_range);
    return;
  }
  ASSERT_EQ(tree.find_by_number(0), tree.root());
  EXPECT_EQ(tree.parent(tree.root()), tree.root());
  EXPECT_TRUE(tree.is_last_child(tree.root()));
  auto view = tree.navigate(tree.root());
  for (size_t p = 0; p < n; ++p) {
    SCOPED_TRACE(::testing::Message() << "number=" << p);
    const size_t v = oracle.order[p];
    const auto node = tree.find_by_number(p);
    view.reset(node);
    ASSERT_EQ(tree.number(node), p);
    ASSERT_EQ(view.number(), p);
    ASSERT_EQ(view.node(), node);
    ASSERT_EQ(tree.is_root(node), v == 0);
    ASSERT_EQ(view.is_root(), v == 0);
    ASSERT_EQ(tree.is_leaf(node), source.children[v].empty());
    ASSERT_EQ(view.is_leaf(), source.children[v].empty());
    ASSERT_EQ(tree.degree(node), source.children[v].size());
    ASSERT_EQ(view.degree(), source.children[v].size());
    ASSERT_EQ(tree.number(tree.parent(node)), oracle.number[source.parents[v]]);
    ASSERT_EQ(tree.parent(node), view.parent());
    const bool last = v == 0 || source.children[source.parents[v]].back() == v;
    ASSERT_EQ(tree.is_last_child(node), last);
    ASSERT_EQ(view.is_last_child(), last);
    if (last) {
      if (n < 300 || p % 997 == 0) {
        EXPECT_THROW(view.next_sibling(), std::out_of_range);
        EXPECT_THROW(tree.next_sibling(node), std::out_of_range);
      }
    } else {
      ASSERT_EQ(tree.number(view.next_sibling()),
                oracle.number[source.next[v]]);
      ASSERT_EQ(view.next_sibling(), tree.next_sibling(node));
    }
    if (source.children[v].empty()) {
      if (n < 300 || p % 997 == 0) {
        EXPECT_THROW(view.first_child(), std::out_of_range);
        EXPECT_THROW(tree.first_child(node), std::out_of_range);
      }
    } else {
      ASSERT_EQ(tree.number(view.first_child()),
                oracle.number[source.children[v].front()]);
      ASSERT_EQ(view.first_child(), tree.first_child(node));
    }
    if (n < 300 || p % 997 == 0) {
      EXPECT_THROW(view.child(source.children[v].size()), std::out_of_range);
      EXPECT_THROW(tree.child(node, source.children[v].size()),
                   std::out_of_range);
      EXPECT_THROW(view.child(std::numeric_limits<size_t>::max()),
                   std::out_of_range);
    }
    size_t i = 0;
    view.for_each_child([&](auto child) {
      ASSERT_LT(i, source.children[v].size());
      ASSERT_EQ(tree.number(child), oracle.number[source.children[v][i]]);
      ASSERT_EQ(view.child(i), child);
      ASSERT_EQ(tree.child(node, i), child);
      ASSERT_EQ(tree.parent(child), node);
      ++i;
    });
    ASSERT_EQ(i, source.children[v].size());
    ASSERT_EQ(view.node(), node);
    if (n < 300 || p % 997 == 0) {
      auto independent = view;
      independent.reset(tree.find_by_number(n - 1));
      EXPECT_EQ(independent.degree(),
                source.children[oracle.order.back()].size());
      EXPECT_EQ(view.node(), node);
      EXPECT_EQ(view.degree(), source.children[v].size());
    }
  }

  // Compare original source IDs, so a wrong layer-major permutation cannot hide
  // behind matching tree-generated numbers. Both walks are iterative.
  using Event = std::pair<unsigned, size_t>;
  std::vector<Event> expected, actual;
  std::vector<std::pair<size_t, bool>> pending{{0, false}};
  while (!pending.empty()) {
    const auto [v, leaving] = pending.back();
    pending.pop_back();
    if (leaving) {
      expected.emplace_back(2, v);
      continue;
    }
    expected.emplace_back(0, v);
    for (size_t child : source.children[v]) {
      expected.emplace_back(1, child);
    }
    pending.emplace_back(v, true);
    for (auto it = source.children[v].rbegin(); it != source.children[v].rend();
         ++it) {
      pending.emplace_back(*it, false);
    }
  }
  auto node = tree.root();
  bool entering = true;
  for (size_t step = 0; step < 2 * n; ++step) {
    view.reset(node);
    const size_t p = tree.number(node);
    ASSERT_LT(p, n);
    if (entering) {
      actual.emplace_back(0, oracle.order[p]);
      view.for_each_child([&](auto child) {
        const size_t child_number = tree.number(child);
        ASSERT_LT(child_number, n);
        actual.emplace_back(1, oracle.order[child_number]);
      });
      if (!view.is_leaf()) {
        node = view.first_child();
        continue;
      }
    }
    actual.emplace_back(2, oracle.order[p]);
    if (view.is_root()) {
      break;
    }
    entering = !view.is_last_child();
    node = entering ? view.next_sibling() : view.parent();
  }
  ASSERT_EQ(expected.size(), 3 * n - 1);
  EXPECT_EQ(actual, expected);
}

void CheckLayeredTree(const LayeredTree& tree,
                      const std::vector<size_t>& parents) {
  LayeredTreeOracle oracle(parents);
  const size_t n = parents.size();
  ASSERT_EQ(tree.size(), n);
  const auto memory = tree.memory_breakdown();
  EXPECT_EQ(memory.object_bytes, sizeof(LayeredTree));
  EXPECT_EQ(memory.total_bytes, tree.memory_usage_bytes());
  EXPECT_EQ(memory.total_bytes,
            memory.object_bytes + memory.forest_bytes + memory.connector_bytes +
                memory.boundary_bytes + memory.forest_index_bytes +
                memory.connector_index_bytes + memory.boundary_index_bytes +
                memory.directory_bytes);
  EXPECT_LE(
      tree.boundary_overhead_bits(),
      static_cast<size_t>(
          static_cast<long double>(tree.boundary_budget_bits_per_node()) * n));
  EXPECT_THROW(tree.layer(tree.layer_count()), std::out_of_range);
  EXPECT_THROW(tree.layer(std::numeric_limits<size_t>::max()),
               std::out_of_range);
  std::vector<size_t> depth_begins;
  size_t depth_end = 0;
  for (size_t l = 0; l < tree.layer_count(); ++l) {
    const auto layer = tree.layer(l);
    ASSERT_EQ(layer.depth_begin, depth_end);
    ASSERT_GT(layer.depth_end, layer.depth_begin);
    if (l + 1 < tree.layer_count()) {
      EXPECT_EQ(layer.depth_end - layer.depth_begin, tree.band_height());
    }
    depth_begins.push_back(layer.depth_begin);
    depth_end = layer.depth_end;
  }
  if (n) {
    ASSERT_FALSE(depth_begins.empty());
    EXPECT_EQ(depth_end,
              *std::max_element(oracle.depth.begin(), oracle.depth.end()) + 1);
    oracle.SetLayers(depth_begins);
  } else {
    EXPECT_TRUE(depth_begins.empty());
  }
  ASSERT_EQ(oracle.order.size(), n);
  std::vector<bool> forest, connector, boundary;
  const bool sparse = tree.portal_mode() == LayeredTree::PortalMode::Sparse;
  const size_t covered =
      tree.layer_count() ? tree.layer(tree.layer_count() - 1).node_begin : 0;
  const size_t width = covered > 1 ? std::bit_width(covered - 1) : 1;
  size_t node_end = 0, boundaries = 0, crossings = 0;
  for (size_t l = 0; l < tree.layer_count(); ++l) {
    SCOPED_TRACE(::testing::Message() << "layer=" << l);
    const auto layer = tree.layer(l);
    ASSERT_EQ(layer.node_begin, node_end);
    ASSERT_GT(layer.node_end, layer.node_begin);
    ASSERT_LE(layer.node_end, n);
    EXPECT_EQ(layer.forest_begin, forest.size());
    EXPECT_EQ(layer.connector_begin, connector.size());
    EXPECT_EQ(layer.boundary_begin, boundary.size());
    for (size_t p = layer.node_begin; p < layer.node_end; ++p) {
      const size_t v = oracle.order[p];
      ASSERT_GE(oracle.depth[v], layer.depth_begin);
      ASSERT_LT(oracle.depth[v], layer.depth_end);
      const bool outgoing = l + 1 < tree.layer_count();
      const bool edge_depth = oracle.depth[v] + 1 == layer.depth_end;
      const size_t degree = oracle.source.children[v].size();
      const size_t local_degree = edge_depth ? 0 : degree;
      forest.insert(forest.end(), local_degree, true);
      forest.push_back(false);
      if (outgoing) {
        if (!sparse) {
          boundary.push_back(edge_depth && degree);
        }
        if (edge_depth && degree) {
          if (sparse) {
            for (size_t i = 0; i < width; ++i) {
              boundary.push_back((p >> i) & 1);
            }
          }
          ++boundaries;
          crossings += degree;
          connector.insert(connector.end(), degree - 1, true);
          connector.push_back(false);
        }
      }
      const auto node = tree.find_by_number(p);
      EXPECT_EQ(node.layer, l);
      EXPECT_EQ(node.local, p - layer.node_begin);
    }
    EXPECT_EQ(layer.forest_end, forest.size());
    EXPECT_EQ(layer.connector_end, connector.size());
    EXPECT_EQ(layer.boundary_end, boundary.size());
    if (l + 1 == tree.layer_count()) {
      EXPECT_EQ(layer.connector_begin, layer.connector_end);
      EXPECT_EQ(layer.boundary_begin, layer.boundary_end);
    }
    node_end = layer.node_end;
  }
  EXPECT_EQ(node_end, n);
  EXPECT_EQ(tree.boundary_count(), boundaries);
  EXPECT_EQ(tree.cross_edge_count(), crossings);
  EXPECT_EQ(tree.forest_bits(), n ? 2 * n - 1 - crossings : 0);
  EXPECT_EQ(tree.connector_bits(), crossings);
  EXPECT_EQ(tree.topology_bits(), n ? 2 * n - 1 : 0);
  EXPECT_EQ(tree.boundary_bits(), boundary.size());
  if (!sparse) {
    EXPECT_EQ(
        static_cast<size_t>(std::count(boundary.begin(), boundary.end(), true)),
        boundaries);
  }
  const auto check_bits = [](std::span<const uint64_t> words, size_t bits,
                             const std::vector<bool>& expected) {
    ASSERT_EQ(bits, expected.size());
    ASSERT_EQ(words.size(), (bits + 63) / 64);
    for (size_t bit = 0; bit < bits; ++bit) {
      ASSERT_EQ((words[bit / 64] >> (bit % 64)) & 1, expected[bit])
          << "bit=" << bit;
    }
    if (bits % 64) {
      EXPECT_EQ(words.back() >> (bits % 64), 0);
    }
  };
  check_bits(tree.forest_words(), tree.forest_bits(), forest);
  check_bits(tree.connector_words(), tree.connector_bits(), connector);
  check_bits(tree.boundary_words(), tree.boundary_bits(), boundary);
  EXPECT_GE(memory.forest_bytes, tree.forest_words().size_bytes());
  EXPECT_GE(memory.connector_bytes, tree.connector_words().size_bytes());
  EXPECT_GE(memory.boundary_bytes, tree.boundary_words().size_bytes());
  EXPECT_GE(tree.boundary_overhead_bits(),
            8 * memory.connector_bytes - tree.connector_bits() +
                8 * (memory.boundary_bytes + memory.connector_index_bytes +
                     memory.boundary_index_bytes));
  if (tree.layer_count() <= 1) {
    EXPECT_EQ(tree.boundary_overhead_bits(), 0);
    EXPECT_EQ(memory.connector_bytes, 0);
    EXPECT_EQ(memory.boundary_bytes, 0);
    EXPECT_EQ(memory.connector_index_bytes, 0);
    EXPECT_EQ(memory.boundary_index_bytes, 0);
  }
  CheckLayeredNavigation(tree, oracle);
}

static void CheckLayeredOrderedShapes(size_t max_nodes,
                                      size_t expected_shapes) {
  size_t shapes = 0;
  for (size_t n = 1; n <= max_nodes; ++n) {
    for (uint64_t bits = 1; bits < (uint64_t{1} << (2 * n)); bits += 2) {
      std::vector<size_t> parents, pending;
      bool valid = true;
      for (size_t bit = 0; bit < 2 * n; ++bit) {
        if ((bits >> bit) & 1) {
          if (pending.empty() && !parents.empty()) {
            valid = false;
            break;
          }
          parents.push_back(pending.empty() ? 0 : pending.back());
          pending.push_back(parents.size() - 1);
        } else {
          if (pending.empty()) {
            valid = false;
            break;
          }
          pending.pop_back();
        }
      }
      if (!valid || !pending.empty() || parents.size() != n) {
        continue;
      }
      ++shapes;
      for (size_t height : {1, 2, 3}) {
        for (double budget : {0.0, 0.05, 1000000.0}) {
          SCOPED_TRACE(::testing::Message()
                       << "n=" << n << " bits=" << bits << " h=" << height
                       << " budget=" << budget);
          const LayeredTree tree(DecompositionSource(parents), height, budget);
          CheckLayeredTree(tree, parents);
          if (budget == 1000000.0) {
            CheckLayeredTree(
                LayeredTree(DecompositionSource(parents), height, budget,
                            LayeredTree::PortalMode::Sparse),
                parents);
            LayeredTreeOracle oracle(parents);
            EXPECT_EQ(
                tree.layer_count(),
                *std::max_element(oracle.depth.begin(), oracle.depth.end()) /
                        height +
                    1);
          }
        }
      }
    }
  }
  EXPECT_EQ(shapes, expected_shapes);
}

TEST(LayeredDFUDSTest, ExhaustiveSmallOrderedShapes) {
  CheckLayeredOrderedShapes(6, 65);
}

TEST(LayeredDFUDSTest, HeavyExhaustiveSmallOrderedShapes) {
  CheckLayeredOrderedShapes(8, 626);
}

static void CheckLayeredShapeBoundaries(std::span<const size_t> sizes,
                                        size_t seeds) {
  for (const size_t n : sizes) {
    for (int shape = 0; shape < 6; ++shape) {
      const auto parents = DecompositionParents(n, shape, 313 + n);
      for (size_t height : {1, 3, 8}) {
        for (double budget : {0.0, 0.05, 1000000.0}) {
          SCOPED_TRACE(::testing::Message()
                       << "n=" << n << " shape=" << shape << " h=" << height
                       << " budget=" << budget);
          CheckLayeredTree(
              LayeredTree(DecompositionSource(parents), height, budget),
              parents);
        }
      }
    }
  }
  for (size_t seed = 0; seed < seeds; ++seed) {
    const auto parents = DecompositionParents(513, 0, seed);
    CheckLayeredTree(LayeredTree(DecompositionSource(parents), 2, 1000000.0),
                     parents);
  }
}

TEST(LayeredDFUDSTest, ShapeBoundariesAndNonPreorderSources) {
  constexpr size_t sizes[] = {1, 63, 64, 65, 127, 128, 129, 255, 256, 257};
  CheckLayeredShapeBoundaries(sizes, 2);
}

TEST(LayeredDFUDSTest, HeavyShapeBoundariesAndNonPreorderSources) {
  constexpr size_t sizes[] = {1,   2,   3,   31,  32,  33,  63,  64,  65,
                              127, 128, 129, 255, 256, 257, 511, 1025};
  CheckLayeredShapeBoundaries(sizes, 12);
}

TEST(LayeredDFUDSTest, EmptyDefaultsAndInvalidOptions) {
  const DecompositionSource empty({});
  const LayeredTree tree(empty);
  EXPECT_TRUE(tree.empty());
  EXPECT_EQ(tree.size(), 0);
  EXPECT_EQ(tree.layer_count(), 0);
  EXPECT_EQ(tree.band_height(), 8);
  EXPECT_DOUBLE_EQ(tree.boundary_budget_bits_per_node(), 0.05);
  EXPECT_EQ(tree.topology_bits(), 0);
  EXPECT_EQ(tree.boundary_overhead_bits(), 0);
  CheckLayeredTree(tree, {});
  CheckLayeredTree(LayeredTree(), {});
  const auto invalid = static_cast<LayeredTree::PortalMode>(2);
  EXPECT_THROW((LayeredTree(empty, 8, 0.05, invalid)), std::invalid_argument);
  EXPECT_THROW((LayeredTree(DecompositionSource({0}), 8, 0.05, invalid)),
               std::invalid_argument);
  for (double budget : {-1.0, -std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}) {
    EXPECT_THROW((LayeredTree(empty, 8, budget)), std::invalid_argument);
    EXPECT_THROW((LayeredTree(DecompositionSource({0}), 8, budget)),
                 std::invalid_argument);
  }
  EXPECT_THROW((LayeredTree(empty, 0, 0.05)), std::invalid_argument);
  EXPECT_THROW((LayeredTree(DecompositionSource({0}), 0, 0.05)),
               std::invalid_argument);
}

TEST(LayeredDFUDSTest, IndependentViewsAndReentrantEnumeration) {
  const LayeredTree tree(DecompositionSource(DecompositionParents(1025, 3, 7)),
                         2, 1000000.0);
  ASSERT_GT(tree.layer_count(), 1);
  auto view = tree.navigate(tree.root());
  std::vector<LayeredTree::Node> children;
  tree.for_each_child(tree.root(),
                      [&](auto child) { children.push_back(child); });
  ASSERT_EQ(children.size(), 2);
  size_t calls = 0;
  view.for_each_child([&](auto child) {
    ASSERT_LT(calls, children.size());
    EXPECT_EQ(child, children[calls]);
    size_t nested_calls = 0;
    view.for_each_child([&](auto sibling) {
      ASSERT_LT(nested_calls, children.size());
      EXPECT_EQ(sibling, children[nested_calls++]);
    });
    EXPECT_EQ(nested_calls, children.size());
    auto other = view;
    other.reset(child);
    other.for_each_child([&](auto grandchild) {
      EXPECT_EQ(tree.parent(grandchild), child);
      EXPECT_EQ(view.node(), tree.root());
    });
    ++calls;
  });
  EXPECT_EQ(calls, children.size());
  EXPECT_THROW(view.for_each_child([&](auto child) {
    auto nested = tree.navigate(child);
    nested.for_each_child([](auto) { throw std::runtime_error("stop"); });
  }),
               std::runtime_error);
  EXPECT_EQ(view.node(), tree.root());
  EXPECT_EQ(view.degree(), 2);
  calls = 0;
  view.for_each_child([&](auto child) { EXPECT_EQ(child, children[calls++]); });
  EXPECT_EQ(calls, children.size());
  auto sum_degrees = [&tree] {
    size_t edges = 0;
    auto cursor = tree.navigate(tree.root());
    for (size_t number = 0; number < tree.size(); ++number) {
      cursor.reset(tree.find_by_number(number));
      if (tree.number(cursor.node()) != number) {
        return size_t{0};
      }
      edges += cursor.degree();
    }
    return edges;
  };
  auto first = std::async(std::launch::async, sum_degrees);
  auto second = std::async(std::launch::async, sum_degrees);
  EXPECT_EQ(first.get(), tree.size() - 1);
  EXPECT_EQ(second.get(), tree.size() - 1);
}

TEST(LayeredDFUDSTest, ZeroChildBoundaryNodesAndDistinctLowerRootParents) {
  // At depth one, child 2 is a leaf between two parents of lower roots.
  // At depth two, leaf 5 separates two lower subtrees in the packed forest.
  const std::vector<size_t> parents{0, 0, 0, 0, 1, 1, 3, 3, 4, 6, 6};
  const LayeredTree tree(DecompositionSource(parents), 2, 1000000.0);
  ASSERT_EQ(tree.layer_count(), 2);
  EXPECT_EQ(tree.boundary_count(), 2);
  EXPECT_EQ(tree.cross_edge_count(), 4);
  CheckLayeredTree(tree, parents);
  LayeredTreeOracle oracle(parents);
  oracle.SetLayers({0, 2});
  EXPECT_NE(oracle.order, oracle.preorder);
  EXPECT_EQ(oracle.order,
            (std::vector<size_t>{0, 1, 2, 3, 4, 8, 5, 6, 9, 10, 7}));
  for (size_t p = 0; p < oracle.order.size(); ++p) {
    EXPECT_EQ(tree.number(tree.find_by_number(p)), p);
  }
}

TEST(LayeredDFUDSTest, BudgetFloorIncludesPaddingIndexesAndExtraDirectory) {
  const auto parents = DecompositionParents(256, 1, 1);
  const DecompositionSource source(parents);
  for (auto mode :
       {LayeredTree::PortalMode::Dense, LayeredTree::PortalMode::Sparse}) {
    const LayeredTree unsplit(source, 128, 0.0, mode);
    CheckLayeredTree(unsplit, parents);
    const LayeredTree forced(source, 128, 1000000.0, mode);
    ASSERT_EQ(unsplit.layer_count(), 1);
    ASSERT_EQ(forced.layer_count(), 2);
    const auto memory = forced.memory_breakdown();
    const size_t extra_directory =
        memory.directory_bytes - unsplit.memory_breakdown().directory_bytes;
    const size_t charged =
        8 * memory.connector_bytes - forced.connector_bits() +
        8 * (memory.boundary_bytes + memory.connector_index_bytes +
             memory.boundary_index_bytes + extra_directory);
    ASSERT_EQ(charged, forced.boundary_overhead_bits());
    ASSERT_GT(charged, 0);
    // n=256 makes both budgets exactly representable in binary; the second
    // permits floor(charged - 0.5), so the only possible cut must be rejected.
    const double exact_budget = static_cast<double>(charged) / 256;
    const double below_budget = (static_cast<double>(charged) - 0.5) / 256;
    const LayeredTree exact(source, 128, exact_budget, mode);
    const LayeredTree below(source, 128, below_budget, mode);
    EXPECT_EQ(exact.layer_count(), 2);
    EXPECT_EQ(exact.boundary_overhead_bits(), charged);
    EXPECT_EQ(below.layer_count(), 1);
    CheckLayeredTree(exact, parents);
    CheckLayeredTree(below, parents);
  }

  const auto deep_parents = DecompositionParents(4096, 1, 1);
  const DecompositionSource deep_source(deep_parents);
  const LayeredTree all_cuts(deep_source, 8, 1000000.0);
  ASSERT_EQ(all_cuts.layer_count(), 512);
  const double limited_budget =
      static_cast<double>(all_cuts.boundary_overhead_bits()) / (2 * 4096);
  const LayeredTree limited(deep_source, 8, limited_budget);
  ASSERT_GT(limited.layer_count(), 1);
  ASSERT_LT(limited.layer_count(), all_cuts.layer_count());
  const auto final = limited.layer(limited.layer_count() - 1);
  EXPECT_GT(final.depth_end - final.depth_begin, 8);
  EXPECT_EQ(final.depth_end, 4096);
  CheckLayeredTree(limited, deep_parents);
}

TEST(LayeredDFUDSTest, SparseMulticutBudgetAtPackedWidthChange) {
  constexpr auto mode = LayeredTree::PortalMode::Sparse;
  // A 66-node chain has exactly 65 available cuts at height one. Its
  // completed allocation gives the charge for that prefix of a longer chain.
  const LayeredTree prefix(DecompositionSource(DecompositionParents(66, 1, 1)),
                           1, 1000000.0, mode);
  ASSERT_EQ(prefix.layer_count(), 66);
  const size_t charge = prefix.boundary_overhead_bits();
  const auto parents = DecompositionParents(128, 1, 1);
  const DecompositionSource source(parents);
  const LayeredTree exact(source, 1, static_cast<double>(charge) / 128, mode);
  const LayeredTree below(source, 1, static_cast<double>(charge - 1) / 128,
                          mode);
  ASSERT_EQ(exact.layer_count(), 66);
  EXPECT_EQ(exact.boundary_overhead_bits(), charge);
  EXPECT_EQ(exact.boundary_bits(), 65 * 7);
  ASSERT_EQ(below.layer_count(), 65);
  EXPECT_EQ(below.boundary_bits(), 64 * 6);
  CheckLayeredTree(exact, parents);
  CheckLayeredTree(below, parents);
}

TEST(LayeredDFUDSTest, OpaqueScrambledSourceHandles) {
  struct Source {
    struct Handle {
      size_t scrambled;
    };
    DecompositionSource source;
    size_t size() const { return source.size(); }
    Handle root() const { return {size_t{0} ^ size_t{0x9e3779b9}}; }
    size_t degree(Handle v) const {
      return source.degree(v.scrambled ^ size_t{0x9e3779b9});
    }
    Handle first_child(Handle v) const {
      return {source.first_child(v.scrambled ^ size_t{0x9e3779b9}) ^
              size_t{0x9e3779b9}};
    }
    Handle next_sibling(Handle v) const {
      return {source.next_sibling(v.scrambled ^ size_t{0x9e3779b9}) ^
              size_t{0x9e3779b9}};
    }
  };
  for (int shape : {0, 3, 4}) {
    const auto parents = DecompositionParents(257, shape, 1729);
    const Source source{DecompositionSource(parents)};
    for (size_t height : {1, 2, 8}) {
      CheckLayeredTree(LayeredTree(source, height, 1000000.0), parents);
    }
  }
}

TEST(LayeredDFUDSTest, OwningCopiesSurviveSourceAndOriginal) {
  const auto parents = DecompositionParents(2049, 0, 19);
  for (auto mode :
       {LayeredTree::PortalMode::Dense, LayeredTree::PortalMode::Sparse}) {
    auto copy = [&] {
      const DecompositionSource source(parents);
      const LayeredTree original(source, 2, 1000000.0, mode);
      return LayeredTree(original);
    }();
    ASSERT_GT(copy.layer_count(), 1);
    CheckLayeredTree(copy, parents);
    LayeredTree moved(std::move(copy));
    CheckLayeredTree(moved, parents);
    LayeredTree assigned(DecompositionSource({0}));
    {
      const LayeredTree temporary = moved;
      assigned = temporary;
    }
    CheckLayeredTree(assigned, parents);
    LayeredTree move_assigned(DecompositionSource({}));
    move_assigned = std::move(assigned);
    CheckLayeredTree(move_assigned, parents);
  }
}

TEST(LayeredDFUDSTest, InvalidNodeCoordinates) {
  const LayeredTree tree(DecompositionSource(DecompositionParents(17, 3, 1)), 1,
                         1000000.0);
  for (LayeredTree::Node node :
       {LayeredTree::Node{tree.layer_count(), 0}, LayeredTree::Node{0, 1},
        LayeredTree::Node{std::numeric_limits<size_t>::max(), 0},
        LayeredTree::Node{0, std::numeric_limits<size_t>::max()}}) {
    EXPECT_THROW(tree.number(node), std::out_of_range);
    EXPECT_THROW(tree.navigate(node), std::out_of_range);
    EXPECT_THROW(tree.parent(node), std::out_of_range);
    EXPECT_THROW(tree.degree(node), std::out_of_range);
    EXPECT_THROW(tree.child(node, 0), std::out_of_range);
    EXPECT_THROW(tree.first_child(node), std::out_of_range);
    EXPECT_THROW(tree.next_sibling(node), std::out_of_range);
    EXPECT_THROW(tree.is_root(node), std::out_of_range);
    EXPECT_THROW(tree.is_leaf(node), std::out_of_range);
    EXPECT_THROW(tree.is_last_child(node), std::out_of_range);
    auto cursor = tree.navigate(tree.root());
    EXPECT_THROW(cursor.reset(node), std::out_of_range);
    EXPECT_EQ(cursor.node(), tree.root());
    EXPECT_EQ(cursor.degree(), tree.degree(tree.root()));
  }
}

TEST(LayeredDFUDSTest, PackedInterlayerSearchesAtEveryBitOffset) {
  std::vector<bool> offsets(64);
  for (size_t height : {1, 2, 3, 8, 31, 32, 33, 64}) {
    const auto parents = DecompositionParents(64 * height + 1, 1, 1);
    const LayeredTree tree(DecompositionSource(parents), height, 1000000.0);
    ASSERT_EQ(tree.layer_count(), 65);
    for (size_t layer = 0; layer < tree.layer_count(); ++layer) {
      offsets[tree.layer(layer).forest_begin % 64] = true;
    }
    // Neighboring forests occupy the same words. Searches must not use their
    // real one/zero bits as if they belonged to the current forest.
    CheckLayeredTree(tree, parents);
  }
  EXPECT_EQ(std::count(offsets.begin(), offsets.end(), true), 64);
}

TEST(LayeredDFUDSTest, RejectsUnrepresentableSizeBeforeTraversal) {
  struct Source {
    size_t size() const { return std::numeric_limits<size_t>::max(); }
    size_t root() const { throw std::runtime_error("must not traverse"); }
    size_t degree(size_t) const { return 0; }
    size_t first_child(size_t) const { return 0; }
    size_t next_sibling(size_t) const { return 0; }
  };
  EXPECT_THROW((LayeredTree(Source{})), std::length_error);
}

TEST(LayeredDFUDSLargeTest, HeavyDeepAndWideWithoutRecursiveStack) {
  const auto chain = DecompositionParents(65537, 1, 1);
  for (double budget : {0.0, 0.05, 1000000.0}) {
    const LayeredTree tree(DecompositionSource(chain), 8, budget);
    if (budget == 1000000.0) {
      EXPECT_EQ(tree.layer_count(), 8193);
    }
    CheckLayeredTree(tree, chain);
  }
  const auto star = DecompositionParents(65537, 2, 1);
  const LayeredTree tree(DecompositionSource(star), 1, 1000000.0);
  ASSERT_EQ(tree.layer_count(), 2);
  EXPECT_EQ(tree.cross_edge_count(), 65536);
  CheckLayeredTree(tree, star);
}

TEST(LayeredDFUDSTest, CachedQueriesResetCopiesAndExceptions) {
  // Internal nodes, boundary leaves, distinct lower roots, and local siblings.
  const std::vector<size_t> parents{0, 0, 0, 0, 1, 1, 3, 3, 4, 6, 6};
  for (size_t height : {1, 2, 8}) {
    for (double budget : {0.0, 1000000.0}) {
      const auto tree = [&] {
        const DecompositionSource source(parents);
        const LayeredTree original(source, height, budget);
        return LayeredTree(original);
      }();
      LayeredTreeOracle oracle(parents);
      std::vector<size_t> depths;
      for (size_t i = 0; i < tree.layer_count(); ++i) {
        depths.push_back(tree.layer(i).depth_begin);
      }
      oracle.SetLayers(depths);
      auto handle = [&](size_t v) {
        return tree.find_by_number(oracle.number[v]);
      };
      auto check = [&](const auto& view, size_t v) {
        const auto& children = oracle.source.children[v];
        const bool last =
            v == 0 || oracle.source.children[parents[v]].back() == v;
        for (size_t repeat = 0; repeat < 3; ++repeat) {
          EXPECT_EQ(view.node(), handle(v));
          EXPECT_EQ(view.number(), oracle.number[v]);
          // Exercise leaf before degree, and incoming before child queries.
          EXPECT_EQ(view.is_leaf(), children.empty());
          EXPECT_EQ(view.is_last_child(), last);
          EXPECT_EQ(view.parent(), handle(parents[v]));
          EXPECT_EQ(view.degree(), children.size());
          EXPECT_EQ(view.is_root(), v == 0);
          if (last) {
            EXPECT_THROW(view.next_sibling(), std::out_of_range);
          } else {
            EXPECT_EQ(view.next_sibling(), handle(oracle.source.next[v]));
          }
          if (children.empty()) {
            EXPECT_THROW(view.first_child(), std::out_of_range);
          } else {
            EXPECT_EQ(view.first_child(), handle(children.front()));
          }
          for (size_t i = children.size(); i > 0; --i) {
            EXPECT_EQ(view.child(i - 1), handle(children[i - 1]));
          }
          EXPECT_THROW(view.child(children.size()), std::out_of_range);
          size_t count = 0;
          view.for_each_child([&](auto child) {
            ASSERT_LT(count, children.size());
            EXPECT_EQ(child, handle(children[count++]));
          });
          EXPECT_EQ(count, children.size());
        }
      };
      auto view = tree.navigate(tree.root());
      for (size_t v : oracle.preorder) {
        view.reset(handle(v));
        check(view, v);
        auto independent = view;  // Copy populated caches, then diverge both.
        const size_t other = (v + 1) % parents.size();
        view.reset(handle(other));
        check(view, other);
        check(independent, v);
        independent.reset(handle(other));
        view.reset(handle(v));
        check(independent, other);
        check(view, v);
        EXPECT_THROW(view.reset((LayeredTree::Node{tree.layer_count(), 0})),
                     std::out_of_range);
        EXPECT_THROW(view.reset((LayeredTree::Node{
                         0, std::numeric_limits<size_t>::max()})),
                     std::out_of_range);
        check(view, v);
        view.reset(view.node());
        check(view, v);
        size_t count = 0;
        view.for_each_child([&](auto child) {
          EXPECT_EQ(child, handle(oracle.source.children[v][count++]));
          size_t nested = 0;
          view.for_each_child([&](auto sibling) {
            EXPECT_EQ(sibling, handle(oracle.source.children[v][nested++]));
            check(view, v);
          });
          EXPECT_EQ(nested, oracle.source.children[v].size());
        });
        EXPECT_EQ(count, oracle.source.children[v].size());
        if (count) {
          EXPECT_THROW(view.for_each_child([&](auto) {
            view.for_each_child([](auto) { throw std::runtime_error("stop"); });
          }),
                       std::runtime_error);
        }
        check(view, v);
      }
    }
  }
}

TEST(LayeredDFUDSTest, ComputedDestinationPositionsAndLocalNumbers) {
  // Mixed nonleaf siblings, local leaves with external children, and many
  // roots.
  std::vector<std::vector<size_t>> shapes{
      {0, 0, 0, 0, 1, 1, 3, 3, 4, 6, 6, 8, 8, 9, 10}};
  shapes.push_back(DecompositionParents(4097, 2, 42));
  shapes.push_back(DecompositionParents(257, 3, 42));
  for (const auto& parents : shapes) {
    for (size_t height : {1, 2, 3, 512}) {
      const LayeredTree tree(DecompositionSource(parents), height, 1000000.0);
      LayeredTreeOracle oracle(parents);
      std::vector<size_t> depths;
      for (size_t layer = 0; layer < tree.layer_count(); ++layer) {
        depths.push_back(tree.layer(layer).depth_begin);
      }
      oracle.SetLayers(depths);
      auto handle = [&](size_t v) {
        return tree.find_by_number(oracle.number[v]);
      };
      auto check = [&](const auto& view, size_t v) {
        EXPECT_EQ(view.node(), handle(v));
        EXPECT_EQ(view.number(), oracle.number[v]);
        EXPECT_EQ(view.degree(), oracle.source.children[v].size());
        EXPECT_EQ(view.is_leaf(), oracle.source.children[v].empty());
        EXPECT_EQ(view.is_last_child(),
                  v == 0 || oracle.source.children[parents[v]].back() == v);
        EXPECT_EQ(view.parent(), handle(parents[v]));
      };
      for (size_t v : oracle.preorder) {
        auto view = tree.navigate(handle(v));
        const auto& children = oracle.source.children[v];
        // Reverse indexed selection exercises skips over multiple raw subtrees.
        for (size_t i = children.size(); i > 0; --i) {
          view.reset(handle(v));
          view.reset(view.child(i - 1));
          check(view, children[i - 1]);
          view.reset(view.parent());
          check(view, v);
        }
        if (children.empty()) {
          continue;
        }
        view.reset(handle(v));
        view.reset(view.first_child());
        for (size_t i = 0; i < children.size(); ++i) {
          check(view, children[i]);
          if (i + 1 < children.size()) {
            const auto next = view.next_sibling();
            auto copy = view;  // Copies carry the computed destination too.
            EXPECT_THROW(view.reset((LayeredTree::Node{tree.layer_count(), 0})),
                         std::out_of_range);
            view.reset(next);
            copy.reset(next);
            check(copy, children[i + 1]);
          } else {
            EXPECT_THROW(view.next_sibling(), std::out_of_range);
          }
        }
        view.reset(handle(v));
        size_t count = 0;
        view.for_each_child([&](auto child) {
          ASSERT_LT(count, children.size());
          EXPECT_EQ(child, handle(children[count++]));
          // Overwrite the caller's pending target during enumeration.
          EXPECT_EQ(view.child(children.size() - 1), handle(children.back()));
        });
        EXPECT_EQ(count, children.size());
        check(view, v);
      }
    }
  }
}

TEST(LayeredDFUDSTest, PendingDestinationManualResetAndNestedExceptions) {
  const LayeredTree tree(DecompositionSource({0, 0, 0, 1, 1, 2, 3, 5}), 2,
                         1000000.0);
  auto view = tree.navigate(tree.root());
  const auto first = view.first_child();
  const auto second = view.child(1);  // Replaces the first hint.
  view.reset(first);  // Must reconstruct, not use second's hint.
  EXPECT_EQ(view.degree(), 2);
  EXPECT_EQ(view.parent(), tree.root());
  view.reset(view.next_sibling());
  EXPECT_EQ(view.node(), second);
  EXPECT_EQ(view.degree(), 1);

  for (size_t number = 0; number < tree.size(); ++number) {
    const auto expected = tree.find_by_number(number);
    view.reset(tree.root());
    (void)view.first_child();
    // Caller-created coordinates carry no position that can disagree with them.
    view.reset({expected.layer, expected.local});
    EXPECT_EQ(view.number(), number);
    EXPECT_EQ(view.degree(), tree.degree(expected));
    EXPECT_EQ(view.parent(), tree.parent(expected));
  }

  view.reset(tree.root());
  const auto target = view.first_child();
  auto copy = view;
  EXPECT_THROW(view.for_each_child([&](auto) {
    view.for_each_child([&](auto) {
      (void)view.child(1);
      throw std::runtime_error("nested stop");
    });
  }),
               std::runtime_error);
  EXPECT_EQ(view.node(), tree.root());
  EXPECT_EQ(view.degree(), 2);
  // The callback replaced only view's hint; copy retains its independent hint.
  view.reset(target);
  copy.reset(target);
  EXPECT_EQ(view.node(), copy.node());
  EXPECT_EQ(view.degree(), 2);
  EXPECT_EQ(copy.degree(), 2);
  EXPECT_EQ(view.parent(), copy.parent());
}

TEST(LayeredDFUDSTest, IncomingEqualExcessBoundaries) {
  // Three bands, several roots per lower forest, and zero-child boundary nodes.
  const std::vector<size_t> parents{0, 0, 0, 1, 1, 2,  2,  3,
                                    3, 5, 7, 7, 9, 10, 10, 5};
  const LayeredTree tree(DecompositionSource(parents), 2, 1000000.0);
  ASSERT_GE(tree.layer_count(), 3);
  LayeredTreeOracle oracle(parents);
  std::vector<size_t> depths;
  for (size_t l = 0; l < tree.layer_count(); ++l) {
    depths.push_back(tree.layer(l).depth_begin);
  }
  oracle.SetLayers(depths);

  // Independently encode the forest and record prefix boundaries from source
  // degrees, rather than obtaining positions through the tree's index.
  std::vector<uint64_t> words((tree.forest_bits() + 63) / 64);
  std::vector<size_t> positions(tree.size());
  std::vector<int64_t> excess{0};
  size_t bits = 0;
  for (size_t l = 0; l < tree.layer_count(); ++l) {
    const auto layer = tree.layer(l);
    ASSERT_EQ(bits, layer.forest_begin);
    for (size_t n = layer.node_begin; n < layer.node_end; ++n) {
      const size_t v = oracle.order[n];
      positions[n] = bits;
      const size_t degree = oracle.depth[v] + 1 == layer.depth_end
                                ? 0
                                : oracle.source.children[v].size();
      for (size_t i = 0; i < degree; ++i) {
        words[bits / 64] |= uint64_t{1} << (bits % 64);
        ++bits;
        excess.push_back(excess.back() + 1);
      }
      ++bits;
      excess.push_back(excess.back() - 1);
    }
  }
  pixie::RmMBTree<> index(words, bits);
  size_t first_roots = 0, later_roots = 0, adjacent = 0;
  size_t at_begin = 0, inside = 0, zero_child_boundaries = 0;
  for (size_t n = 0; n < tree.size(); ++n) {
    const size_t v = oracle.order[n];
    const auto node = tree.find_by_number(n);
    const auto layer = tree.layer(node.layer);
    const size_t p = positions[n], begin = layer.forest_begin;
    size_t expected = index.npos;
    for (size_t q = p; q > 0;) {
      --q;
      if (excess[q] == excess[p]) {
        expected = q;
        break;
      }
    }
    const size_t open = index.bwdsearch(p, 0);
    ASSERT_EQ(open, expected);
    const bool root = oracle.depth[v] == layer.depth_begin;
    if (p == begin) {
      EXPECT_TRUE(root);
      ++first_roots;
    } else if (root) {
      // Every completed raw DFUDS tree ends at a strict new global minimum,
      // so a real later root cannot have a same-excess match in an earlier
      // layer.
      EXPECT_EQ(open, index.npos);
      ++later_roots;
    } else {
      ASSERT_NE(open, index.npos);
      EXPECT_GE(open, begin);
      if (p - begin >= 2 && index.bit(p - 2)) {
        EXPECT_EQ(open, p - 2);
        ++adjacent;
      } else if (open == begin) {
        ++at_begin;  // Equality must survive the non-adjacent search path too.
      } else {
        ++inside;
      }
    }
    if (node.layer + 1 < tree.layer_count() &&
        oracle.depth[v] + 1 == layer.depth_end &&
        oracle.source.children[v].empty()) {
      ++zero_child_boundaries;
    }
  }
  EXPECT_GT(first_roots, 1);
  EXPECT_GT(later_roots, 0);
  EXPECT_GT(adjacent, 0);
  EXPECT_GT(at_begin, 0);
  EXPECT_GT(inside, 0);
  EXPECT_GT(zero_child_boundaries, 0);
  CheckLayeredTree(tree, parents);  // Parent, last-child and sibling queries.
}

TEST(LayeredDFUDSTest, BackwardMatchBeforeSyntheticForestBegin) {
  // An arbitrary prefix can share a root's excess. This cannot be produced by
  // concatenating complete DFUDS forests, whose ends are strict prefix minima.
  // Prefix 11 precedes forest 00: its second root starts at p=3, begin=2.
  const std::vector<uint64_t> words{0b0011};
  const pixie::RmMBTree<> index(words, 4);
  const size_t begin = 2, p = 3;
  const size_t open = index.bwdsearch(p, 0);
  ASSERT_EQ(open, 1);
  EXPECT_LT(open, begin);
  EXPECT_EQ(index.bwdsearch(0, 0), index.npos);
}

}  // namespace
