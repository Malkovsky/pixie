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

#include <random>

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

template <class Tree>
struct TreeFormat;
template <class B, class S>
struct TreeFormat<pixie::BPTree<B, S>> {
  static constexpr uint64_t invalid_two_node_encoding = 0b0101;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return pixie::adj_to_bp(n, adj);
  }
};
template <class B, class S>
struct TreeFormat<pixie::DFUDSTree<B, S>> {
  static constexpr uint64_t invalid_two_node_encoding = 0b010;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return pixie::adj_to_dfuds(n, adj);
  }
};
template <class B, class S>
struct TreeFormat<pixie::LoudsTree<B, S>> {
  static constexpr uint64_t invalid_two_node_encoding = 0b011;
  static auto encode(size_t n, const std::vector<std::vector<size_t>>& adj) {
    return adj_to_louds(n, bfs_order(n, adj));
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
           RestoredRootIndexCanBeAdoptedByTraversalContext) {
  std::mt19937_64 rng(98);
  const auto words =
      TreeFormat<TypeParam>::encode(1025, generate_random_tree(1025, rng));
  const TypeParam tree(words, 1025);
  for (size_t budget : std::initializer_list<size_t>{0, 32768, 524288}) {
    const pixie::TreeTraversalContext<TypeParam> context(tree, budget);
    const auto bytes = save(context.root_index());
    pixie::BinaryReader reader(bytes);
    auto restored = pixie::RootNavigationIndex<TypeParam>::deserialize(
        reader, tree, budget);
    const pixie::TreeTraversalContext<TypeParam> loaded(std::move(restored),
                                                        budget);
    EXPECT_EQ(loaded.cached_depth(), context.cached_depth());
    EXPECT_EQ(loaded.decoded_chunk_bits(), context.decoded_chunk_bits());
    EXPECT_EQ(loaded.scratch_bytes(), context.scratch_bytes());
    EXPECT_LE(loaded.budget_usage_bytes(), budget);
    EXPECT_EQ(save(loaded.root_index()), bytes);
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
