#include <gtest/gtest.h>
#include <pixie/tree/implementations.h>
#include <pixie/utils.h>

#include <cstddef>
#include <cstdint>
#include <queue>
#include <random>
#include <stack>
#include <utility>
#include <vector>

namespace {

struct LoudsCase {
  using Tree = pixie::LoudsTree;

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

struct BpCase {
  using Tree = pixie::BPTree<pixie::RmMTree>;

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

struct DfudsCase {
  using Tree = pixie::DFUDSTree<pixie::RmMTree>;

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
class TreeSpecificationTest : public testing::Test {};

using TreeCases = testing::Types<LoudsCase, BpCase, DfudsCase>;
TYPED_TEST_SUITE(TreeSpecificationTest, TreeCases);

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

}  // namespace
