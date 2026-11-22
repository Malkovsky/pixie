#include <gtest/gtest.h>
#include <pixie/tree/event_part.h>
#include <pixie/tree/forest_cover.h>
#include <pixie/tree/root_navigation_index.h>

namespace {
using Bp = pixie::tree::BPTree<pixie::RmMBTree<>>;
using Dfuds = pixie::tree::DFUDSTree<pixie::RmMBTree<>>;
using Part = pixie::tree::TreePart;

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
  const auto partial =
      pixie::tree::RootNavigationIndex<Bp, 1>(tree).inspect_part();
  ASSERT_EQ(partial.nodes().size(), 1u);
  EXPECT_EQ(partial.root_count(), 1u);
  EXPECT_TRUE(partial.nodes()[0].children.empty());
  EXPECT_FALSE(partial.nodes()[0].children_complete);

  const auto full = [] {
    const auto owner = star(3);
    return pixie::tree::RootNavigationIndex(owner).inspect_part();
  }();
  ASSERT_EQ(full.root_count(), 1u);
  ASSERT_EQ(full.nodes().size(), 3u);
  EXPECT_EQ(full.nodes()[0].children, (std::vector<size_t>{1, 2}));
  EXPECT_TRUE(full.nodes()[0].children_complete);
  EXPECT_TRUE(full.nodes()[1].children_complete);
  EXPECT_TRUE(full.nodes()[1].children.empty());
  EXPECT_EQ(full.nodes()[1].source_id, 1u);
  EXPECT_EQ(full.nodes()[1].bit_position, 1u);
}

template <class Tree>
void check_subtree_topology(const Tree& tree) {
  // Root with two children. The first child has two leaves, and the next
  // sibling's encoding immediately follows that complete three-node subtree.
  const auto child = tree.first_child(tree.root());
  const auto part = pixie::tree::inspect_event_part<128>(tree, child);
  ASSERT_TRUE(part);
  EXPECT_EQ(part->root_count(), 1u);
  ASSERT_EQ(part->nodes().size(), 3u);
  EXPECT_EQ(part->nodes()[0].children, (std::vector<size_t>{1, 2}));
  for (size_t i = 0; i < 3; ++i) {
    const auto& node = part->nodes()[i];
    EXPECT_EQ(node.source_id, child.number + i);
    const auto source = i ? tree.child(child, i - 1) : child;
    EXPECT_EQ(node.bit_position, source.pos);
    EXPECT_TRUE(node.children_complete);
    EXPECT_EQ(node.children.size(), tree.degree(source));
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
  ASSERT_EQ(hit->nodes().size(), 64u);
  EXPECT_EQ(hit->nodes()[0].children.size(), 63u);
  EXPECT_FALSE(pixie::tree::inspect_event_part<128>(misses, misses.root()));
  const auto last = misses.child(misses.root(), 63);
  const auto leaf = pixie::tree::inspect_event_part<128>(misses, last);
  ASSERT_TRUE(leaf);
  ASSERT_EQ(leaf->nodes().size(), 1u);
  EXPECT_EQ(leaf->nodes()[0].bit_position, last.pos);
  EXPECT_TRUE(leaf->nodes()[0].children_complete);
  EXPECT_TRUE(leaf->nodes()[0].children.empty());
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
    std::vector<Part> parts;
    for (size_t u = 0; u < tree.unit_count(); ++u) {
      parts.push_back(tree.inspect_unit(u));
    }
    return parts;
  }();
  ASSERT_GT(parts.size(), 1u);
  ASSERT_EQ(parts[0].root_count(), 1u);
  const auto& portal = parts[0].nodes().at(0);
  EXPECT_FALSE(portal.children_complete);
  EXPECT_TRUE(portal.children.empty());
  EXPECT_GT(parts[1].root_count(), 1u);
  size_t next = 1;
  for (size_t u = 1; u < parts.size(); ++u) {
    ASSERT_EQ(parts[u].root_count(), parts[u].nodes().size());
    for (size_t root = 0; root < parts[u].root_count(); ++root) {
      const auto& leaf = parts[u].nodes().at(root);
      EXPECT_EQ(leaf.source_id, next++);
      EXPECT_TRUE(leaf.children_complete);
      EXPECT_TRUE(leaf.children.empty());
      EXPECT_FALSE(leaf.bit_position);
    }
  }
  EXPECT_EQ(next, 700u);
}

TEST(TreePartInspection, ForestIndexAndUnitUseTheirOwnLocalCoordinates) {
  const pixie::tree::ForestCoverTree tree(star(5));
  const auto unit = tree.inspect_unit(0);
  const auto index =
      pixie::tree::RootNavigationIndex<pixie::tree::ForestCoverTree>(tree)
          .inspect_part();
  for (const Part* part : {&index, &unit}) {
    ASSERT_EQ(part->root_count(), 1u);
    ASSERT_EQ(part->nodes().size(), 5u);
    const auto& parent = part->nodes()[0];
    EXPECT_TRUE(parent.children_complete);
    ASSERT_EQ(parent.children.size(), tree.degree(tree.root()));
    for (size_t i = 0; i < parent.children.size(); ++i) {
      const auto& leaf = part->nodes().at(parent.children[i]);
      EXPECT_EQ(leaf.source_id, tree.number(tree.child(tree.root(), i)));
      EXPECT_FALSE(leaf.bit_position);
      EXPECT_TRUE(leaf.children_complete);
      EXPECT_TRUE(leaf.children.empty());
    }
  }
}
}  // namespace
