#include <gtest/gtest.h>
#include <pixie/tree/event_part.h>
#include <pixie/tree/forest_cover.h>
#include <pixie/tree/root_navigation_index.h>

namespace {
using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;

Bp star(size_t n) {
  std::vector<uint64_t> words((2 * n + 63) / 64);
  if (n) {
    words[0] = 1;
  }
  for (size_t i = 1; i < n; ++i) {
    const size_t pos = 2 * i - 1;
    words[pos / 64] |= uint64_t{1} << (pos % 64);
  }
  return Bp(std::span<const uint64_t>(words), n);
}

TEST(TreePartInspection, DisabledIndexFrontierAndLeafAreDifferent) {
  const auto tree = star(3);
  const auto disabled =
      pixie::tree::RootNavigationIndex(tree, 0).inspect_part();
  EXPECT_TRUE(disabled.empty());
  EXPECT_EQ(disabled.root_count(), 0u);
  EXPECT_THROW(disabled.root(), std::out_of_range);
  EXPECT_TRUE(disabled.parentheses().empty());
  const auto partial =
      pixie::tree::RootNavigationIndex<Bp, 1>(tree).inspect_part();
  ASSERT_EQ(partial.size(), 1u);
  EXPECT_EQ(partial.root_count(), 1u);
  EXPECT_EQ(&partial.source_tree(), &tree);
  const auto frontier = partial.navigate(partial.root());
  EXPECT_TRUE(frontier.is_leaf());
  EXPECT_FALSE(frontier.children_complete());

  const auto full = [] {
    const auto owner = star(3);
    return pixie::tree::RootNavigationIndex(owner).inspect_part();
  }();
  ASSERT_EQ(full.root_count(), 1u);
  ASSERT_EQ(full.size(), 3u);
  EXPECT_EQ(full.find_by_number(full.size() - 1), 2u);
  EXPECT_THROW(full.find_by_number(full.size()), std::out_of_range);
  const auto root = full.root();
  EXPECT_EQ(full.parent(root), root);
  EXPECT_EQ(full.first_child(root), 1u);
  EXPECT_EQ(full.child(root, 1), 2u);
  EXPECT_TRUE(full.children_complete(root));
  EXPECT_TRUE(full.children_complete(1));
  EXPECT_TRUE(full.is_leaf(1));
  EXPECT_FALSE(full.is_last_child(1));
  EXPECT_EQ(full.next_sibling(1), 2u);
  EXPECT_TRUE(full.is_last_child(2));
  EXPECT_EQ(full.source_node(1), pixie::tree::TreeNode(1, 1));
  EXPECT_EQ(full.tree_path(2, 1), (std::vector<size_t>{2, 0, 1}));
  EXPECT_EQ(full.parentheses(), "(()())");
  EXPECT_TRUE(full.parentheses(1, 0).empty());
  EXPECT_EQ(full.parentheses(1, 4), "()()");
}

TEST(TreePartInspection, RootIndexQueriesDistinguishFrontierAndAbsentNodes) {
  const auto tree = star(5);
  const auto root = tree.root();
  const auto child = tree.first_child(root);
  const pixie::tree::RootNavigationIndex<Bp, 1> frontier(tree);
  EXPECT_TRUE(frontier.contains(root));
  EXPECT_FALSE(frontier.contains(child));
  EXPECT_FALSE(frontier.children_complete(root));
  EXPECT_TRUE(frontier.represented_children(root).empty());
  EXPECT_THROW(frontier.children_complete(child), std::out_of_range);
  EXPECT_THROW(frontier.represented_children(child), std::out_of_range);
  EXPECT_THROW(frontier.represented_parent(child), std::out_of_range);

  const pixie::tree::RootNavigationIndex full(tree);
  EXPECT_TRUE(full.contains(child));
  EXPECT_TRUE(full.children_complete(child));
  EXPECT_TRUE(full.represented_children(child).empty());
  EXPECT_EQ(full.represented_parent(root), root);
  EXPECT_EQ(full.represented_parent(child), root);
}

template <class Tree>
void check_subtree_topology(const Tree& tree) {
  // Root with two children. The first child has two leaves, and the next
  // sibling's encoding immediately follows that complete three-node subtree.
  const auto child = tree.first_child(tree.root());
  const auto part = pixie::tree::inspect_event_part<128>(tree, child);
  ASSERT_TRUE(part);
  EXPECT_EQ(part->root_count(), 1u);
  ASSERT_EQ(part->size(), 3u);
  EXPECT_EQ(&part->source_tree(), &tree);
  EXPECT_EQ(part->first_child(part->root()), 1u);
  EXPECT_EQ(part->child(part->root(), 1), 2u);
  for (size_t i = 0; i < part->size(); ++i) {
    const auto source = i ? tree.child(child, i - 1) : child;
    EXPECT_EQ(part->source_node(i), source);
    EXPECT_TRUE(part->children_complete(i));
    EXPECT_EQ(part->degree(i), tree.degree(source));
  }
}

TEST(TreePartInspection, DecoderProbeCopiesOnlyTheSelectedSubtree) {
  // BP bits in source order: 1110100100.
  const uint64_t bits = 0b0010010111;
  const Bp bp(std::span(&bits, 1), 5);
  check_subtree_topology(bp);
  check_subtree_topology(Dfuds(bp));
}

template <class Tree>
void check_probe_limit(const Tree& empty,
                       const Tree& fits,
                       const Tree& misses) {
  EXPECT_FALSE(pixie::tree::inspect_event_part<128>(empty, {0, 0}));
  const auto hit = pixie::tree::inspect_event_part<128>(fits, fits.root());
  ASSERT_TRUE(hit);
  ASSERT_EQ(hit->size(), 64u);
  EXPECT_EQ(hit->degree(hit->root()), 63u);
  EXPECT_FALSE(pixie::tree::inspect_event_part<128>(misses, misses.root()));
  const auto last = misses.child(misses.root(), 63);
  const auto leaf = pixie::tree::inspect_event_part<128>(misses, last);
  ASSERT_TRUE(leaf);
  ASSERT_EQ(leaf->size(), 1u);
  EXPECT_EQ(leaf->source_node(leaf->root()), last);
  EXPECT_TRUE(leaf->children_complete(leaf->root()));
  EXPECT_TRUE(leaf->is_leaf(leaf->root()));
}

TEST(TreePartInspection, EventProbeAcceptsOnlyCompleteSubtrees) {
  const auto empty = star(0), fits = star(64), misses = star(65);
  check_probe_limit(empty, fits, misses);
  check_probe_limit(Dfuds(empty), Dfuds(fits), Dfuds(misses));
}

TEST(TreePartInspection,
     ForestPortalOmitsChildrenAndChildUnitsHaveLeadingRoots) {
  const auto parts = [] {
    const pixie::tree::ForestCoverTree tree(star(700));
    EXPECT_THROW(tree.inspect_unit(tree.unit_count()), std::out_of_range);
    std::vector<pixie::tree::ForestUnitPart> parts;
    for (size_t u = 0; u < tree.unit_count(); ++u) {
      parts.push_back(tree.inspect_unit(u));
    }
    return parts;
  }();
  ASSERT_GT(parts.size(), 1u);
  ASSERT_EQ(parts[0].root_count(), 1u);
  const auto& portal_part = parts[0];
  EXPECT_FALSE(portal_part.children_complete(portal_part.root()));
  EXPECT_TRUE(portal_part.is_leaf(portal_part.root()));
  EXPECT_LT(portal_part.memory_usage_bytes(), 640u);
  ASSERT_GT(parts[1].root_count(), 1u);
  EXPECT_TRUE(parts[1].tree_path(parts[1].root(0), parts[1].root(1)).empty());
  EXPECT_TRUE(parts[1].tree_path(parts[1].root(1), parts[1].root(0)).empty());
  EXPECT_EQ(parts[1].parentheses(1, 3), ")()");
  size_t leaves = 0;
  for (size_t u = 1; u < parts.size(); ++u) {
    const auto& part = parts[u];
    ASSERT_EQ(part.root_count(), part.size());
    for (size_t i = 0; i < part.root_count(); ++i) {
      const auto root = part.root(i);
      EXPECT_EQ(part.parent(root), root);
      EXPECT_TRUE(part.children_complete(root));
      EXPECT_TRUE(part.is_leaf(root));
      ++leaves;
    }
  }
  EXPECT_EQ(leaves, 699u);
}

TEST(TreePartInspection, ForestIndexAndUnitUseTheirOwnLocalCoordinates) {
  const pixie::tree::ForestCoverTree tree(star(5));
  const auto unit = tree.inspect_unit(0);
  const auto index =
      pixie::tree::RootNavigationIndex<pixie::tree::ForestCoverTree>(tree)
          .inspect_part();
  const auto check_star = [&](const auto& part) {
    EXPECT_EQ(&part.source_tree(), &tree);
    ASSERT_EQ(part.root_count(), 1u);
    ASSERT_EQ(part.size(), 5u);
    const auto root = part.root();
    EXPECT_TRUE(part.children_complete(root));
    ASSERT_EQ(part.degree(root), tree.degree(tree.root()));
    for (size_t i = 0; i < part.degree(root); ++i) {
      const auto leaf = part.child(root, i);
      EXPECT_EQ(part.number(leaf), i + 1);
      EXPECT_EQ(part.parent(leaf), root);
      EXPECT_TRUE(part.children_complete(leaf));
      EXPECT_TRUE(part.is_leaf(leaf));
      EXPECT_EQ(part.is_last_child(leaf), i + 1 == part.degree(root));
    }
  };
  check_star(index);
  check_star(unit);
  EXPECT_EQ(unit.tree_path(1, 4), (std::vector<size_t>{1, 0, 4}));
}
}  // namespace
