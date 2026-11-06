#include <gtest/gtest.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rank_select/support.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>
#include <pixie/tree/traversal_context.h>
#include <pixie/utils.h>

#include <algorithm>
#include <array>
#include <random>

#include "tree_encoding_test_utils.h"

namespace {
using pixie::DeserializationValidation;
using Compressed = pixie::RankSelectSupportCompressed<>;
using Baseline = pixie::RankSelectSupport<>;
using Rmm = pixie::RmMBTree<>;
using BaselineRmm = pixie::RmMBTree<4, 32, Baseline>;

template <class T>
std::vector<std::byte> save(const T& value) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  value.serialize(writer);
  writer.finish();
  return sink.take();
}

void set_u64(std::vector<std::byte>& bytes, size_t offset, uint64_t value) {
  ASSERT_LE(offset + 8, bytes.size());
  for (size_t i = 0; i < 8; ++i) {
    bytes[offset + i] = std::byte(value >> (i * 8));
  }
}

uint64_t get_u64(const std::vector<std::byte>& bytes, size_t offset) {
  pixie::BinaryReader reader(std::span<const std::byte>(bytes).subspan(offset));
  return reader.read_u64();
}

TEST(CompressedSerialization, MixedCoarseDirectoriesRoundTripAndRejectDamage) {
  constexpr size_t super = 16384;
  constexpr size_t bits = 80 * super + 37;
  std::vector<uint64_t> words((bits + 63) / 64);
  for (size_t s : {0u, 40u, 79u}) {
    std::fill_n(words.begin() + s * super / 64, super / 64, UINT64_MAX);
  }
  words[3] &= ~(uint64_t{1} << 17);
  for (bool complement : {false, true}) {
    if (complement) {
      for (auto& word : words) {
        word = ~word;
      }
    }
    const Compressed original(words, bits);
    const Baseline reference(words, bits);
    const auto bytes = save(original);
    for (auto validation : {DeserializationValidation::kQuick,
                            DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(bytes);
      const auto restored = Compressed::deserialize(reader, words, validation);
      EXPECT_TRUE(reader.empty());
      EXPECT_EQ(save(restored), bytes);
      for (size_t rank = 1; rank <= bits + 1; rank += 4095) {
        EXPECT_EQ(restored.select(rank), reference.select(rank));
        EXPECT_EQ(restored.select0(rank), reference.select0(rank));
      }
    }

    // Walk the version-1 wire arrays, retaining their count-field offsets.
    // Mutate fields inside the selected directory without breaking framing,
    // so validation reaches the sample records rather than rejecting a header.
    size_t offset = 64;
    auto array = [&](size_t element_bytes) {
      const size_t begin = offset;
      offset += 8 + get_u64(bytes, offset) * element_bytes;
      return begin;
    };
    const auto ranks = array(8);
    std::vector<std::pair<size_t, uint64_t>> mutations{
        {ranks + 8, 1}, {48, original.rank(bits) - 1}};
    bool found_mixed = false;
    for (unsigned direction = 0; direction < 2; ++direction) {
      const auto shape = get_u64(bytes, offset + 8);
      EXPECT_EQ(get_u64(bytes, offset), 0u);  // Coarse, not predicted.
      offset += 16;
      const auto samples = array(8);
      const auto directory = array(8);
      const auto map = array(16);
      const auto position_header = offset;
      offset += 24;
      array(8);  // Empty position samples.
      array(8);  // Occupied-superblock sample records.
      array(8);  // Occupied-superblock directory.
      const auto wide = array(16);
      array(2);  // Rank endpoints for the default sampling frequency.
      mutations.emplace_back(position_header, 1);  // Noncanonical empty layout.
      if (shape != 0) {
        continue;
      }
      found_mixed = true;
      ASSERT_GT(get_u64(bytes, samples), 1u);
      ASSERT_GT(get_u64(bytes, directory), 0u);
      ASSERT_GT(get_u64(bytes, map), 0u);
      ASSERT_GT(get_u64(bytes, wide), 0u);
      mutations.insert(
          mutations.end(),
          {{samples + 8, UINT64_MAX},  // Base outside the source domain.
           {samples + 16, 0},          // Missing unary highs.
           {directory + 8, UINT64_MAX},
           {map + 8, 1},            // Incorrect bitmap prefix.
           {map + 16, UINT64_MAX},  // Bits beyond the final sample.
           {map + 16, 0},           // Bitmap population mismatch.
           {wide + 8, 1},           // Noncontiguous endpoint directory.
           {wide + 16, 0}});        // Empty wide interval.
    }
    ASSERT_TRUE(found_mixed);
    for (const auto& [field, value] : mutations) {
      SCOPED_TRACE(testing::Message() << complement << "/" << field);
      auto damaged = bytes;
      set_u64(damaged, field, value);
      pixie::BinaryReader reader(damaged);
      EXPECT_THROW(Compressed::deserialize(reader, words),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  }
}

template <class Tree>
struct TreeFormat;
template <class B, class S>
struct TreeFormat<pixie::BPTree<B, S>> {
  using Backend = B;
  static constexpr uint64_t invalid_two_node_encoding = 0b0101;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return EncodeReference<pixie::detail::TreeEncodingKind::kBP>(n, adj);
  }
};
template <class B, class S>
struct TreeFormat<pixie::DFUDSTree<B, S>> {
  using Backend = B;
  static constexpr uint64_t invalid_two_node_encoding = 0b010;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return EncodeReference<pixie::detail::TreeEncodingKind::kDFUDS>(n, adj);
  }
};
template <class B, class S>
struct TreeFormat<pixie::LoudsTree<B, S>> {
  using Backend = B;
  static constexpr uint64_t invalid_two_node_encoding = 0b011;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return EncodeReference<pixie::detail::TreeEncodingKind::kLOUDS>(
        n, bfs_order(n, adj));
  }
};

template <class Tree>
class TreeSerializationTest : public testing::Test {};
using Trees = testing::Types<pixie::BPTree<Rmm>,
                             pixie::BPTree<BaselineRmm>,
                             pixie::DFUDSTree<Rmm>,
                             pixie::DFUDSTree<BaselineRmm>,
                             pixie::LoudsTree<Compressed>,
                             pixie::LoudsTree<Baseline>>;
TYPED_TEST_SUITE(TreeSerializationTest, Trees);

template <class A, class B>
void compare_trees(const A& expected, const B& actual) {
  ASSERT_EQ(actual.size(), expected.size());
  if (actual.empty()) {
    return;
  }
  std::vector<std::pair<pixie::TreeNode, pixie::TreeNode>> pending;
  pending.push_back({expected.root(), actual.root()});
  size_t visited = 0;
  while (!pending.empty()) {
    const auto [a, b] = pending.back();
    pending.pop_back();
    ASSERT_EQ(a.number, b.number);
    ASSERT_EQ(a.pos, b.pos);
    ASSERT_EQ(expected.degree(a), actual.degree(b));
    ASSERT_EQ(expected.parent(a).pos, actual.parent(b).pos);
    std::vector<pixie::TreeNode> left, right;
    expected.for_each_child(a, [&](auto child) { left.push_back(child); });
    actual.for_each_child(b, [&](auto child) { right.push_back(child); });
    ASSERT_EQ(left.size(), right.size());
    for (size_t i = 0; i < left.size(); ++i) {
      pending.push_back({left[i], right[i]});
    }
    ASSERT_LE(++visited, actual.size());
  }
  EXPECT_EQ(visited, actual.size());
}

TYPED_TEST(TreeSerializationTest,
           OwningAndViewRoundTripsAndConcatenatedArtifacts) {
  std::mt19937_64 rng(42);
  for (size_t n :
       std::initializer_list<size_t>{0, 1, 2, 31, 32, 33, 513, 17003}) {
    SCOPED_TRACE(n);
    const auto words =
        n ? TreeFormat<TypeParam>::encode(n, generate_random_tree(n, rng))
          : std::vector<uint64_t>{};
    const TypeParam tree(words, n);
    auto bytes = save(tree);
    const auto single = bytes.size();
    const auto second = bytes;
    bytes.insert(bytes.end(), second.begin(), second.end());
    pixie::BinaryReader reader(bytes);
    auto restored =
        TypeParam::deserialize(reader, DeserializationValidation::kFull);
    EXPECT_EQ(reader.position(), single);
    compare_trees(tree, restored);
    EXPECT_EQ(save(restored), second);
    // Aligned backing is required by TreeView, independent of reader lifetime.
    pixie::AlignedStorage backing(
        std::span<const std::byte>(bytes).subspan(single));
    pixie::BinaryReader mapped(backing.as_bytes());
    auto view = TypeParam::TreeView::deserialize(mapped);
    EXPECT_TRUE(mapped.empty());
    compare_trees(tree, view);
    if (n) {
      EXPECT_NE(restored.encoding_words().data(), tree.encoding_words().data());
      EXPECT_GE(
          reinterpret_cast<const std::byte*>(view.encoding_words().data()),
          backing.as_bytes().data());
      EXPECT_EQ(view.encoding_owned_bytes(), 0u);
    }
    auto moved = std::move(restored);
    compare_trees(tree, moved);
    auto next = TypeParam::deserialize(reader);
    EXPECT_TRUE(reader.empty());
    compare_trees(tree, next);
  }
}

TYPED_TEST(TreeSerializationTest, RootIndexRestoresUnderRuntimeBudget) {
  std::mt19937_64 rng(71);
  constexpr size_t n = 4097;
  const auto words =
      TreeFormat<TypeParam>::encode(n, generate_random_tree(n, rng));
  const TypeParam tree(words, n);
  using Index = pixie::RootNavigationIndex<TypeParam>;
  for (size_t budget : std::initializer_list<size_t>{0, 1024, 32768, 524288}) {
    const Index index(tree, budget);
    const auto bytes = save(index);
    pixie::BinaryReader reader(bytes);
    const auto restored = Index::deserialize(reader, tree, budget,
                                             DeserializationValidation::kFull);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(save(restored), bytes);
    EXPECT_LE(restored.budget_usage_bytes(), budget);
    for (size_t path = 0; path < 50; ++path) {
      auto a = index.cursor();
      auto b = restored.cursor();
      while (a.degree()) {
        ASSERT_EQ(a.degree(), b.degree());
        const auto child = rng() % a.degree();
        a.move_to_child(child);
        b.move_to_child(child);
        EXPECT_EQ(a.node().number, b.node().number);
        EXPECT_EQ(a.node().pos, b.node().pos);
      }
      EXPECT_EQ(b.degree(), 0u);
    }
    if (index.budget_usage_bytes()) {
      pixie::BinaryReader too_small(bytes);
      EXPECT_THROW(Index::deserialize(too_small, tree, 0),
                   std::invalid_argument);
      EXPECT_EQ(too_small.position(), 0u);
    }
  }
}

TYPED_TEST(TreeSerializationTest,
           RejectsTruncationAndWrongVersionTransactionally) {
  std::mt19937_64 rng(73);
  const auto words =
      TreeFormat<TypeParam>::encode(5, generate_random_tree(5, rng));
  const TypeParam tree(words, 5);
  auto bytes = save(tree);
  for (size_t n = 0; n < bytes.size(); ++n) {
    pixie::BinaryReader reader(std::span<const std::byte>(bytes).first(n));
    EXPECT_THROW(TypeParam::deserialize(reader), std::invalid_argument);
    ASSERT_EQ(reader.position(), 0u);
  }
  set_u64(bytes, 8, 99);
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW(TypeParam::deserialize(reader), std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TYPED_TEST(TreeSerializationTest,
           RestoredRootIndexCanBeBorrowedByTraversalContext) {
  std::mt19937_64 rng(98);
  const auto words =
      TreeFormat<TypeParam>::encode(1025, generate_random_tree(1025, rng));
  const TypeParam tree(words, 1025);
  for (size_t budget : std::initializer_list<size_t>{0, 32768, 524288}) {
    const pixie::RootNavigationIndex<TypeParam> index(tree, budget);
    const pixie::TreeTraversalContext<TypeParam> context(tree, budget, &index);
    const auto bytes = save(*context.root_index());
    pixie::BinaryReader reader(bytes);
    auto restored = pixie::RootNavigationIndex<TypeParam>::deserialize(
        reader, tree, budget);
    const pixie::TreeTraversalContext<TypeParam> loaded(tree, budget,
                                                        &restored);
    EXPECT_EQ(loaded.cached_depth(), context.cached_depth());
    EXPECT_EQ(loaded.decoded_chunk_bits(), context.decoded_chunk_bits());
    EXPECT_EQ(loaded.scratch_bytes(), context.scratch_bytes());
    EXPECT_LE(loaded.budget_usage_bytes(), budget);
    EXPECT_EQ(save(*loaded.root_index()), bytes);
  }
}

TEST(CompressedSerialization, AllSelectModesAndSamplingShapes) {
  std::mt19937_64 rng(74);
  // Dense prediction, sparse wide intervals, homogeneous blocks and tiny tails.
  for (size_t n : std::initializer_list<size_t>{0, 1, 511, 512, 1024, 1025,
                                                16384, 1000003}) {
    for (unsigned shape = 0; shape < 4; ++shape) {
      std::vector<uint64_t> words(n / 64 + (n % 64 != 0));
      for (size_t i = 0; i < words.size(); ++i) {
        words[i] = shape == 0   ? rng()
                   : shape == 1 ? (i % 2000 == 0 ? 1 : 0)
                   : shape == 2 ? (i / 256 % 2 ? UINT64_MAX : 0)
                                : UINT64_MAX;
      }
      for (unsigned mode = 0; mode < 4; ++mode) {
        SCOPED_TRACE(testing::Message() << n << "/" << shape << "/" << mode);
        const Compressed original(words, n,
                                  static_cast<Compressed::SelectSupport>(mode));
        const auto bytes = save(original);
        pixie::BinaryReader reader(bytes);
        const auto restored = Compressed::deserialize(
            reader, words, DeserializationValidation::kFull);
        EXPECT_TRUE(reader.empty());
        EXPECT_EQ(save(restored), bytes);
        for (size_t i = 0; i < 200; ++i) {
          const auto end = n ? rng() % (n + 1) : 0;
          EXPECT_EQ(restored.rank(end), original.rank(end));
          const auto r1 = rng() % (original.rank(n) + 2);
          const auto r0 = rng() % (original.rank0(n) + 2);
          EXPECT_EQ(restored.select(r1), original.select(r1));
          EXPECT_EQ(restored.select0(r0), original.select0(r0));
        }
      }
    }
  }
}

TEST(RmMBTreeSerialization, StandaloneEmptyPartialAndMultipleLevels) {
  std::mt19937_64 rng(75);
  for (size_t n : std::initializer_list<size_t>{0, 1, 511, 512, 513, 65539}) {
    std::vector<uint64_t> words(n / 64 + (n % 64 != 0));
    for (auto& word : words) {
      word = rng();
    }
    const Rmm original(words, n);
    const auto bytes = save(original);
    pixie::BinaryReader reader(bytes);
    const auto restored =
        Rmm::deserialize(reader, words, DeserializationValidation::kFull);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(save(restored), bytes);
    for (size_t i = 0; i < 100 && n; ++i) {
      const auto left = rng() % n;
      const auto right = left + rng() % (n - left);
      EXPECT_EQ(restored.rank1(left), original.rank1(left));
      EXPECT_EQ(restored.range_min_query_pos(left, right),
                original.range_min_query_pos(left, right));
      EXPECT_EQ(restored.fwdsearch(left, 1), original.fwdsearch(left, 1));
      EXPECT_EQ(restored.bwdsearch(right, -1), original.bwdsearch(right, -1));
    }
  }
}

TYPED_TEST(TreeSerializationTest, FullRejectsInvalidTopologyWithValidBackend) {
  const std::vector<uint64_t> words{
      TreeFormat<TypeParam>::invalid_two_node_encoding};
  const TypeParam invalid(words, 2);
  const auto bytes = save(invalid);
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW(TypeParam::deserialize(reader, DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TYPED_TEST(TreeSerializationTest,
           RootIndexRejectsMalformedEntriesAndWrongTree) {
  constexpr size_t n = 5;
  std::vector<std::vector<size_t>> star(n), chain(n);
  star[0] = {0, 1, 2, 3, 4};
  chain[0] = {0, 1};
  for (size_t i = 1; i < n; ++i) {
    star[i] = {0};
    chain[i] = {i - 1};
    if (i + 1 < n) {
      chain[i].push_back(i + 1);
    }
  }
  const auto words = TreeFormat<TypeParam>::encode(n, star);
  const auto other_words = TreeFormat<TypeParam>::encode(n, chain);
  const TypeParam tree(words, n), other(other_words, n);
  using Index = pixie::RootNavigationIndex<TypeParam>;
  const auto bytes = save(Index(tree));
  pixie::BinaryReader wrong_tree(bytes);
  EXPECT_THROW(Index::deserialize(wrong_tree, other, 524288,
                                  DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(wrong_tree.position(), 0u);
  for (size_t offset :
       std::initializer_list<size_t>{0, 8, 16, 24, 32, 40, 48, 56, 64}) {
    auto bad = bytes;
    set_u64(bad, offset, UINT64_MAX);
    pixie::BinaryReader reader(bad);
    EXPECT_THROW(Index::deserialize(reader, tree), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TYPED_TEST(TreeSerializationTest, RootIndexRejectsBrokenGraphAndContextBudget) {
  constexpr size_t n = 7;
  const std::vector<std::vector<size_t>> adj{
      {0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {1}, {1}, {2}, {2}};
  const auto words = TreeFormat<TypeParam>::encode(n, adj);
  const TypeParam tree(words, n);
  using Index = pixie::RootNavigationIndex<TypeParam>;
  const auto bytes = save(Index(tree));
  // Version 1: 24-byte frame, tree size, entry count, then 32-byte entries.
  constexpr size_t entries = 40;
  ASSERT_EQ(get_u64(bytes, 32), n);
  const std::vector<std::pair<size_t, uint64_t>> mutations{
      {entries + 16, 0},           // Children attached to a zero-degree node.
      {entries + 24, 0},           // Backward/self child reference.
      {entries + 16, 1},           // An unreferenced cached child.
      {entries, 1},                // Wrong root identity.
      {entries + 2 * 32 + 24, 3},  // Overlap the first child's child list.
      {entries + 32, get_u64(bytes, entries + 2 * 32)},
      {entries + 3 * 32 + 16, 1}};  // Wrong leaf degree.
  for (const auto& [field, value] : mutations) {
    SCOPED_TRACE(field);
    auto damaged = bytes;
    set_u64(damaged, field, value);
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(Index::deserialize(reader, tree, 524288,
                                    DeserializationValidation::kFull),
                 std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
  pixie::BinaryReader reader(bytes);
  auto index = Index::deserialize(reader, tree);
  const pixie::TreeTraversalContext<TypeParam> context(tree, 0, &index);
  EXPECT_EQ(context.budget_usage_bytes(), 0u);
  EXPECT_EQ(context.root_index(), &index);
}

TEST(RmMBTreeSerialization, RejectsDamagedSummaryAndTruncatedLevel) {
  const std::vector<uint64_t> words(9, 0x5555555555555555ULL);
  const Rmm original(words, 513);
  const auto bytes = save(original);
  // Header plus four scalar fields precede the nested rank/select frame.
  const size_t summary = 56 + get_u64(bytes, 56 + 16);
  for (const auto& [field, value] : std::vector<std::pair<size_t, uint64_t>>{
           {40, 512},  // Nested rank index still describes 513 bits.
           {summary + 8, original.rank1(513) + 1}}) {
    auto damaged = bytes;
    set_u64(damaged, field, value);
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(Rmm::deserialize(reader, words), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
  auto truncated = bytes;
  truncated.resize(summary + 48);
  set_u64(truncated, 16, truncated.size());
  pixie::BinaryReader short_level(truncated);
  EXPECT_THROW(Rmm::deserialize(short_level, words), std::invalid_argument);
  EXPECT_EQ(short_level.position(), 0u);

  auto empty = save(Rmm{});
  set_u64(empty, 56 + 24, 1);  // Nonzero minimum in an empty summary.
  pixie::BinaryReader empty_reader(empty);
  EXPECT_THROW(Rmm::deserialize(empty_reader, std::span<const uint64_t>{}),
               std::invalid_argument);
  EXPECT_EQ(empty_reader.position(), 0u);
}

TYPED_TEST(TreeSerializationTest, RejectsBackendSizeAndPopulationMismatch) {
  using Backend = typename TreeFormat<TypeParam>::Backend;
  const std::vector<std::vector<size_t>> adj{{0, 1}, {0}};
  const auto words = TreeFormat<TypeParam>::encode(2, adj);
  const TypeParam tree(words, 2);
  const auto bytes = save(tree);
  const size_t bits = std::same_as<TypeParam, pixie::BPTree<Backend>> ? 4 : 3;
  const std::vector<uint64_t> zeros(1, 0);
  // Replace a complete nested artifact: it is independently valid, but does
  // not belong to the enclosing tree. Keep the outer framing valid as well.
  for (const auto& backend :
       {save(Backend(words, bits + 1)), save(Backend(zeros, bits))}) {
    auto damaged = bytes;
    damaged.resize(48);  // Frame, node count, encoding length and one word.
    damaged.insert(damaged.end(), backend.begin(), backend.end());
    set_u64(damaged, 16, damaged.size());
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(TypeParam::deserialize(reader), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(TreeSerialization, RejectsLoudsBackendWithoutBothSelectDirections) {
  const std::vector<uint64_t> words{0b101};
  const auto bytes = save(pixie::LoudsTree<Compressed>(words, 2));
  for (auto mode :
       {Compressed::SelectSupport::kNone, Compressed::SelectSupport::kSelect1,
        Compressed::SelectSupport::kSelect0}) {
    auto damaged = bytes;
    damaged.resize(48);
    const auto backend = save(Compressed(words, 3, mode));
    damaged.insert(damaged.end(), backend.begin(), backend.end());
    set_u64(damaged, 16, damaged.size());
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(pixie::LoudsTree<Compressed>::deserialize(reader),
                 std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(RmMBTreeSerialization, ValidatesEachLeafSummaryFieldAndPadding) {
  const std::vector<uint64_t> words{0b0101};
  const auto bytes = save(Rmm(words, 4));
  const size_t leaf = 56 + get_u64(bytes, 72) + 48;
  // Each low node stores four arrays of 32 little-endian 16-bit fields.
  // The source's summary is (excess=0, min=0, max=1, min_count=2).
  auto set_field = [](auto& data, size_t offset, int value) {
    const auto encoded = static_cast<uint16_t>(value);
    data.at(offset) = std::byte(encoded);
    data.at(offset + 1) = std::byte(encoded >> 8);
  };
  const std::vector<std::array<int, 4>> invalid{
      {-6, 0, 1, 2}, {6, 0, 1, 2},  {1, 0, 1, 2}, {0, -5, 1, 2}, {0, 0, 5, 2},
      {0, 2, 1, 2},  {-2, 0, 1, 2}, {2, 0, 1, 2}, {0, 0, 1, 0},  {0, 0, 1, 5}};
  const std::vector<std::array<int, 4>> inconsistent{
      {2, 0, 2, 2}, {0, -1, 1, 2}, {0, 0, 2, 2}, {0, 0, 1, 1}};
  for (auto validation :
       {DeserializationValidation::kQuick, DeserializationValidation::kFull}) {
    const auto& cases = validation == DeserializationValidation::kQuick
                            ? invalid
                            : inconsistent;
    for (const auto& values : cases) {
      SCOPED_TRACE(testing::PrintToString(values));
      auto damaged = bytes;
      for (size_t field = 0; field < 4; ++field) {
        set_field(damaged, leaf + field * 64, values[field]);
      }
      pixie::BinaryReader reader(damaged);
      EXPECT_THROW(Rmm::deserialize(reader, words, validation),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  }
  for (size_t field = 0; field < 4; ++field) {
    auto damaged = bytes;
    set_field(damaged, leaf + field * 64 + 2, 1);
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(Rmm::deserialize(reader, words), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(CompressedSerialization, DenseMixedPositionMapsRejectDamage) {
  constexpr size_t super = 16384;
  constexpr size_t bits = 80 * super + 37;
  std::vector<uint64_t> words((bits + 63) / 64);
  for (size_t s : {0u, 1u, 2u, 3u, 4u, 5u, 40u, 79u}) {
    std::fill_n(words.begin() + s * super / 64, super / 64, UINT64_MAX);
  }
  words[3] &= ~(uint64_t{1} << 17);
  const auto bytes = save(Compressed(words, bits));
  size_t offset = 64;
  auto array = [&](size_t element_bytes) {
    const size_t begin = offset;
    offset += 8 + get_u64(bytes, offset) * element_bytes;
    return begin;
  };
  array(8);                                   // Rank records.
  ASSERT_EQ(get_u64(bytes, offset), 1u);      // Predicted positions.
  ASSERT_EQ(get_u64(bytes, offset + 8), 0u);  // Mixed shape.
  offset += 16;
  array(8);   // Empty coarse records.
  array(8);   // Empty coarse directory.
  array(16);  // Empty coarse bitmap.
  const size_t header = offset;
  ASSERT_EQ(get_u64(bytes, header + 8), 4u);  // Inline wide bitmap.
  offset += 24;
  const size_t records = array(8) + 8;
  const std::vector<std::pair<size_t, uint64_t>> mutations{
      {header, 0},
      {header, get_u64(bytes, header) + 1},
      {header + 8, 3},
      {header + 16, 60},
      {records, UINT64_MAX},
      {records + 8, 0},
      {records + 16, 1},                   // Incorrect wide prefix.
      {records + 24, 0},                   // Incorrect wide population.
      {records + 24, uint64_t{1} << 63}};  // Out-of-range bitmap bit.
  for (const auto& [field, value] : mutations) {
    SCOPED_TRACE(testing::Message() << field << "/" << value);
    auto damaged = bytes;
    set_u64(damaged, field, value);
    pixie::BinaryReader reader(damaged);
    EXPECT_THROW(Compressed::deserialize(reader, words), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(CompressedSerialization, NondefaultParametersAndSourceMismatch) {
  using Custom =
      pixie::RankSelectSupportCompressed<pixie::AlignedStorage, 64, 2>;
  std::mt19937_64 rng(92);
  std::vector<uint64_t> words(20000);
  for (size_t i = 0; i < words.size(); ++i) {
    words[i] = i < 1000 || i > 15000 ? rng() : (i % 100 == 0 ? 1 : 0);
  }
  const Custom original(words, words.size() * 64 - 11);
  const auto bytes = save(original);
  pixie::BinaryReader reader(bytes);
  const auto restored =
      Custom::deserialize(reader, words, DeserializationValidation::kFull);
  EXPECT_EQ(save(restored), bytes);
  for (size_t i = 1; i <= original.rank(original.size()); i += 37) {
    EXPECT_EQ(restored.select(i), original.select(i));
  }
  for (size_t i = 1; i <= original.rank0(original.size()); i += 37) {
    EXPECT_EQ(restored.select0(i), original.select0(i));
  }
  pixie::BinaryReader wrong_config(bytes);
  EXPECT_THROW(Compressed::deserialize(wrong_config, words),
               std::invalid_argument);
  words[0] ^= 1;
  pixie::BinaryReader wrong_source(bytes);
  EXPECT_THROW(Custom::deserialize(wrong_source, words,
                                   DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(wrong_source.position(), 0u);
}

TEST(RmMBTreeSerialization, RejectsConfigurationSourceAndMalformedMetadata) {
  const std::vector<uint64_t> words(9, 0x5555555555555555ULL);
  const Rmm tree(words, 513);
  const auto bytes = save(tree);
  pixie::BinaryReader short_source(bytes);
  EXPECT_THROW(
      Rmm::deserialize(short_source, std::span<const uint64_t>(words).first(1)),
      std::invalid_argument);
  using OtherFanout = pixie::RmMBTree<2, 16>;
  pixie::BinaryReader wrong_config(bytes);
  EXPECT_THROW(OtherFanout::deserialize(wrong_config, words),
               std::invalid_argument);
  auto different = words;
  different[0] ^= 1;
  pixie::BinaryReader wrong_source(bytes);
  EXPECT_THROW(Rmm::deserialize(wrong_source, different,
                                DeserializationValidation::kFull),
               std::invalid_argument);
  for (size_t offset :
       std::initializer_list<size_t>{0, 8, 16, 24, 32, 40, 48}) {
    auto bad = bytes;
    set_u64(bad, offset, UINT64_MAX);
    pixie::BinaryReader reader(bad);
    EXPECT_THROW(Rmm::deserialize(reader, words), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
  const auto empty = save(Rmm{});
  pixie::BinaryReader reader(empty);
  auto restored = Rmm::deserialize(reader, std::span<const uint64_t>{});
  EXPECT_EQ(restored.size(), 0u);
  EXPECT_EQ(save(restored), empty);
}

TYPED_TEST(TreeSerializationTest, FullValidationHandlesMutatedNestedSections) {
  std::mt19937_64 rng(91);
  const auto words =
      TreeFormat<TypeParam>::encode(33, generate_random_tree(33, rng));
  const TypeParam tree(words, 33);
  const auto bytes = save(tree);
  for (size_t offset = 0; offset < bytes.size(); ++offset) {
    SCOPED_TRACE(offset);
    auto damaged = bytes;
    damaged[offset] ^= std::byte{0xff};
    pixie::BinaryReader reader(damaged);
    try {
      const auto restored =
          TypeParam::deserialize(reader, DeserializationValidation::kFull);
      // Without a hash, a changed encoding can still describe a valid tree
      // with identical rank/min-max metadata. Verify consistency with the
      // loaded encoding, not identity with the original tree.
      const TypeParam rebuilt(restored.encoding_words(), restored.size());
      compare_trees(rebuilt, restored);
    } catch (const std::invalid_argument&) {
      EXPECT_EQ(reader.position(), 0u);
    } catch (const std::length_error&) {
      EXPECT_EQ(reader.position(), 0u);
    }
  }
}

struct CountingRmm : Rmm {
  static inline size_t builds = 0;
  CountingRmm(std::span<const uint64_t> words, size_t bits) : Rmm(words, bits) {
    ++builds;
  }
  static CountingRmm deserialize(pixie::BinaryReader& reader,
                                 std::span<const uint64_t> words,
                                 DeserializationValidation validation) {
    return CountingRmm(Rmm::deserialize(reader, words, validation));
  }

 private:
  explicit CountingRmm(Rmm value) : Rmm(std::move(value)) {}
};

TEST(TreeSerialization, QuickLoadDoesNotReconstructBackend) {
  const std::vector<uint64_t> words{0b0011};
  using Tree = pixie::BPTree<CountingRmm>;
  CountingRmm::builds = 0;
  const Tree tree(words, 2);
  const auto bytes = save(tree);
  pixie::BinaryReader reader(bytes);
  auto restored = Tree::deserialize(reader);
  EXPECT_EQ(CountingRmm::builds, 1u);
  compare_trees(tree, restored);
}
}  // namespace
