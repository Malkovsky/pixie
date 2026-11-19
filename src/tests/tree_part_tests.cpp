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
  EXPECT_TRUE(disabled.roots().empty());
  const auto partial =
      pixie::tree::RootNavigationIndex<Bp, 1>(tree).inspect_part();
  ASSERT_EQ(partial.nodes().size(), 1u);
  EXPECT_EQ(partial.nodes()[0].degree, 2u);
  EXPECT_TRUE(partial.nodes()[0].children.empty());
  EXPECT_FALSE(partial.nodes()[0].children_complete);

  const auto full = [] {
    const auto owner = star(3);
    return pixie::tree::RootNavigationIndex(owner).inspect_part();
  }();
  EXPECT_EQ(full.kind(), Part::Kind::kRootIndex);
  ASSERT_EQ(full.roots(), std::vector<size_t>{0});
  ASSERT_EQ(full.nodes().size(), 3u);
  EXPECT_EQ(full.nodes()[0].children, (std::vector<size_t>{1, 2}));
  EXPECT_TRUE(full.nodes()[0].children_complete);
  EXPECT_TRUE(full.nodes()[1].children_complete);
  EXPECT_TRUE(full.nodes()[1].children.empty());
  EXPECT_EQ(full.nodes()[1].source_id, 1u);
  EXPECT_EQ(full.nodes()[1].bit_position, 1u);
  EXPECT_FALSE(full.source_bits());
  EXPECT_TRUE(full.events().empty());
}

template <class Tree>
void check_events(const Tree& tree, Part::Kind kind, size_t leaf_bits) {
  // Root with two children. The first child has two leaves, and the next
  // sibling's encoding immediately follows that complete three-node subtree.
  const auto child = tree.first_child(tree.root());
  const auto part = pixie::tree::inspect_event_part<128>(tree, child);
  ASSERT_TRUE(part);
  EXPECT_EQ(part->kind(), kind);
  EXPECT_EQ(part->source_bits(),
            (Part::Range{child.pos, child.pos + 4 + leaf_bits}));
  EXPECT_EQ(part->roots(), std::vector<size_t>{0});
  ASSERT_EQ(part->nodes().size(), 3u);
  EXPECT_EQ(part->nodes()[0].children, (std::vector<size_t>{1, 2}));
  for (size_t i = 0; i < 3; ++i) {
    const auto& node = part->nodes()[i];
    EXPECT_EQ(node.source_id, child.number + i);
    EXPECT_EQ(node.bit_position, i ? tree.child(child, i - 1).pos : child.pos);
    EXPECT_TRUE(node.children_complete);
    EXPECT_EQ(node.degree, node.children.size());
  }
  using E = Part::Event::Type;
  EXPECT_EQ(part->events(), (std::vector<Part::Event>{{E::kEnter, 0},
                                                      {E::kInspect, 1},
                                                      {E::kInspect, 2},
                                                      {E::kEnter, 1},
                                                      {E::kEnter, 2},
                                                      {E::kReturn, 0}}));
}

TEST(TreePartInspection, EventPartsHaveExactOrderedSubtreeBoundaries) {
  // BP bits in source order: 1110100100.
  const uint64_t bits = 0b0010010111;
  const Bp bp(std::span(&bits, 1), 5);
  check_events(bp, Part::Kind::kBpEvents, 2);
  check_events(Dfuds(bp), Part::Kind::kDfudsEvents, 1);
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

TEST(TreePartInspection, ForestPortalNamesUnitsWithoutMaterializingChildren) {
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
  ASSERT_EQ(parts[0].roots(), std::vector<size_t>{0});
  const auto& portal = parts[0].nodes().at(0);
  EXPECT_FALSE(portal.children_complete);
  EXPECT_TRUE(portal.children.empty());
  EXPECT_EQ(portal.degree, 699u);
  ASSERT_TRUE(portal.child_units);
  EXPECT_EQ(*portal.child_units, (Part::Range{1, parts.size()}));
  size_t next = 1;
  for (size_t u = 1; u < parts.size(); ++u) {
    EXPECT_EQ(parts[u].kind(), Part::Kind::kForestUnit);
    ASSERT_EQ(parts[u].roots().size(), parts[u].nodes().size());
    for (auto root : parts[u].roots()) {
      const auto& leaf = parts[u].nodes().at(root);
      EXPECT_EQ(leaf.source_id, next++);
      EXPECT_TRUE(leaf.children_complete);
      EXPECT_TRUE(leaf.children.empty());
      EXPECT_FALSE(leaf.child_units);
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
  EXPECT_EQ(index.kind(), Part::Kind::kRootIndex);
  EXPECT_EQ(unit.kind(), Part::Kind::kForestUnit);
  for (const Part* part : {&index, &unit}) {
    ASSERT_EQ(part->roots(), std::vector<size_t>{0});
    ASSERT_EQ(part->nodes().size(), 5u);
    const auto& parent = part->nodes()[0];
    EXPECT_TRUE(parent.children_complete);
    EXPECT_EQ(parent.degree, 4u);
    ASSERT_EQ(parent.children.size(), parent.degree);
    for (size_t i = 0; i < parent.degree; ++i) {
      const auto& leaf = part->nodes().at(parent.children[i]);
      EXPECT_EQ(leaf.source_id, tree.number(tree.child(tree.root(), i)));
      EXPECT_FALSE(leaf.bit_position);
      EXPECT_TRUE(leaf.children_complete);
      EXPECT_TRUE(leaf.children.empty());
    }
  }
}
}  // namespace
