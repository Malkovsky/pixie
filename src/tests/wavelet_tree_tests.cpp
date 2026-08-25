#include <gtest/gtest.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/utils.h>
#include <pixie/wavelet_tree/implementations.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <span>
#include <vector>

using WaveletTree = pixie::WaveletTree<std::uint64_t>;

namespace {

void overwrite_u64(std::vector<std::byte>& bytes,
                   std::size_t offset,
                   std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
  }
}

struct WaveletNodeOffsets {
  std::size_t parent;
  std::size_t left_child;
  std::size_t right_child;
  std::size_t middle;
  std::size_t rank_num_bits;
};

struct WaveletArtifactOffsets {
  std::size_t alphabet_size;
  std::size_t data_size;
  std::size_t root;
  std::size_t node_count;
  std::vector<WaveletNodeOffsets> nodes;
  std::size_t leaves;
  std::size_t permutation;
};

WaveletArtifactOffsets locate_wavelet_artifact(
    std::span<const std::byte> artifact) {
  pixie::BinaryReader reader(artifact);
  reader.skip(3 * sizeof(std::uint64_t));
  const std::size_t alphabet_size = reader.position();
  const std::size_t alphabet_count = reader.read_size();
  const std::size_t data_size = reader.position();
  reader.skip(sizeof(std::uint64_t));
  const std::size_t root = reader.position();
  reader.skip(sizeof(std::uint64_t));
  const std::size_t node_count_offset = reader.position();
  const std::size_t node_count = reader.read_size();

  const auto skip_storage = [&reader] { reader.skip(reader.read_size()); };
  std::vector<WaveletNodeOffsets> nodes;
  nodes.reserve(node_count);
  for (std::size_t node = 0; node < node_count; ++node) {
    const std::size_t parent = reader.position();
    reader.skip(sizeof(std::uint64_t));
    const std::size_t left_child = reader.position();
    reader.skip(sizeof(std::uint64_t));
    const std::size_t right_child = reader.position();
    reader.skip(sizeof(std::uint64_t));
    const std::size_t middle = reader.position();
    reader.skip(sizeof(std::uint64_t));
    nodes.push_back({parent, left_child, right_child, middle, 0});

    skip_storage();
    nodes.back().rank_num_bits = reader.position();
    reader.skip(7 * sizeof(std::uint64_t) + 2 * sizeof(std::uint32_t));
    for (std::size_t storage = 0; storage < 3; ++storage) {
      skip_storage();
    }
  }
  const std::size_t leaves = reader.position();
  const std::size_t permutation =
      leaves + alphabet_count * sizeof(std::uint64_t);
  return {alphabet_size,    data_size, root,       node_count_offset,
          std::move(nodes), leaves,    permutation};
}

}  // namespace

TEST(WaveletTreeTest, BasicSelect) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  std::vector<std::vector<size_t>> rank(alphabet_size);
  for (size_t i = 0; i < data_size; i++) {
    rank[data[i]].push_back(i);
  }

  WaveletTree wavelet_tree(alphabet_size, data);

  for (uint64_t symb = 0; symb < alphabet_size; symb++) {
    for (size_t i = 0; i <= rank[symb].size(); i++) {
      uint64_t exp = i == rank[symb].size() ? data_size : rank[symb][i];
      uint64_t act = wavelet_tree.select(symb, i + 1);
      EXPECT_EQ(act, exp);
    }
  }
}

TEST(WaveletTreeTest, BasicRank) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  std::vector<size_t> count(alphabet_size);

  WaveletTree wavelet_tree(alphabet_size, data);
  for (size_t i = 0; i <= data_size; i++) {
    for (uint64_t symb = 0; symb < alphabet_size; symb++) {
      uint64_t exp = count[symb];
      uint64_t act = wavelet_tree.rank(symb, i);
      EXPECT_EQ(act, exp);
    }

    if (i == data_size) {
      break;
    }
    count[data[i]]++;
  }
}

TEST(WaveletTreeTest, BasicSegment) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  WaveletTree wavelet_tree(alphabet_size, data);

  for (size_t begin = 0; begin <= data_size; begin++) {
    for (size_t end = begin; end <= data_size; end++) {
      auto segment = wavelet_tree.get_segment(begin, end);
      EXPECT_EQ(segment.size(), end - begin);
      for (size_t i = 0; i < end - begin; i++) {
        EXPECT_EQ(segment[i], data[begin + i]);
      }
    }
  }
}

TEST(WaveletTreeTest, SingleSymbolSupportsQueriesAndSerialization) {
  const std::vector<std::uint8_t> data = {0, 0, 0, 0};
  for (const auto build_type : {pixie::WaveletTreeBuildType::Standard,
                                pixie::WaveletTreeBuildType::Huffman}) {
    const pixie::WaveletTree<std::uint8_t> tree(1, data, build_type);
    EXPECT_EQ(tree.rank(0, 3), 3u);
    EXPECT_EQ(tree.select(0, 4), 3u);
    EXPECT_EQ(tree.select(0, 5), tree.size());
    EXPECT_EQ(tree.get_segment(1, 4), (std::vector<std::uint8_t>{0, 0, 0}));

    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    tree.serialize(writer);
    writer.finish();
    const std::vector<std::byte> artifact = output.take();
    pixie::BinaryReader reader(artifact);
    const auto restored =
        pixie::WaveletTreeView<std::uint8_t>::deserialize(reader);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(restored.get_segment(0, data.size()), data);
  }
}

TEST(WaveletTreeTest, OwningCopiesRetainIndependentRankSources) {
  const std::vector<std::uint8_t> data = {0, 1, 0, 1};
  const std::vector<std::uint8_t> initial_data = {1, 1, 0, 0};
  std::unique_ptr<pixie::WaveletTree<std::uint8_t>> constructed_copy;
  pixie::WaveletTree<std::uint8_t> assigned_copy(2, initial_data);
  {
    const pixie::WaveletTree<std::uint8_t> original(2, data);
    constructed_copy =
        std::make_unique<pixie::WaveletTree<std::uint8_t>>(original);
    assigned_copy = original;
  }

  for (const auto* copy : {constructed_copy.get(), &assigned_copy}) {
    EXPECT_EQ(copy->rank(1, 3), 1u);
    EXPECT_EQ(copy->select(1, 2), 3u);
    EXPECT_EQ(copy->get_segment(0, data.size()), data);
  }
}

TEST(WaveletTreeTest, TypedByteSymbolsRoundTripWithoutWidening) {
  std::vector<std::uint8_t> data(256);
  for (std::size_t symbol = 0; symbol < data.size(); ++symbol) {
    data[symbol] = static_cast<std::uint8_t>(symbol);
  }
  const pixie::WaveletTree<std::uint8_t> tree(
      256, data, pixie::WaveletTreeBuildType::Huffman);
  EXPECT_EQ(tree.get_segment(0, data.size()), data);
  EXPECT_EQ(tree.select(255, 1), 255u);

  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  tree.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();

  pixie::BinaryReader reader(artifact);
  const auto restored =
      pixie::WaveletTreeView<std::uint8_t>::deserialize(reader);
  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(restored.get_segment(0, data.size()), data);

  pixie::BinaryReader wrong_symbol_reader(artifact);
  EXPECT_THROW((void)pixie::WaveletTreeView<std::uint16_t>::deserialize(
                   wrong_symbol_reader),
               std::invalid_argument);
  EXPECT_EQ(wrong_symbol_reader.position(), 0u);
}

TEST(WaveletTreeTest, BuildsFromCountsAndOneStreamedPass) {
  const std::vector<std::uint8_t> data = {3, 0, 1, 3, 2, 1, 0};
  const std::array<std::size_t, 4> counts = {2, 2, 1, 2};
  std::size_t passes = 0;
  const pixie::WaveletTree<std::uint8_t> tree(
      4, counts,
      [&](auto&& emit) {
        ++passes;
        for (const std::uint8_t symbol : data) {
          emit(symbol);
        }
      },
      pixie::WaveletTreeBuildType::Huffman);
  EXPECT_EQ(passes, 1u);
  EXPECT_EQ(tree.get_segment(0, data.size()), data);

  const std::array<std::size_t, 4> wrong_counts = {2, 2, 2, 1};
  EXPECT_THROW((pixie::WaveletTree<std::uint8_t>(
                   4, wrong_counts,
                   [&](auto&& emit) {
                     for (const std::uint8_t symbol : data) {
                       emit(symbol);
                     }
                   })),
               std::invalid_argument);
  EXPECT_THROW((pixie::WaveletTree<std::uint8_t>(257, data)),
               std::invalid_argument);
}

TEST(WaveletTreeTest, SmokeSelect) {
  std::vector<std::vector<size_t>> rank;
  for (size_t data_size = 8; data_size < (1 << 22); data_size <<= 1) {
    size_t alphabet_size = 1024;
    std::mt19937_64 rng(239);
    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);

    rank.assign(alphabet_size, {});
    for (size_t i = 0; i < data_size; i++) {
      rank[data[i]].push_back(i);
    }

    for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                            pixie::WaveletTreeBuildType::Huffman}) {
      WaveletTree wavelet_tree(alphabet_size, data, build_type);

      for (uint64_t symb = 0; symb < alphabet_size; symb++) {
        for (size_t i = 0; i <= rank[symb].size(); i++) {
          uint64_t exp = i == rank[symb].size() ? data_size : rank[symb][i];
          uint64_t act = wavelet_tree.select(symb, i + 1);
          EXPECT_EQ(act, exp);
        }
      }
    }
  }
}

TEST(WaveletTreeTest, SmokeRank) {
  std::vector<size_t> count;
  for (size_t data_size = 8; data_size < (1 << 22); data_size <<= 1) {
    size_t alphabet_size = 1024;
    std::mt19937_64 rng(239);
    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);
    std::vector<uint64_t> query =
        generate_random_data(data_size + 1, alphabet_size, rng);

    for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                            pixie::WaveletTreeBuildType::Huffman}) {
      count.assign(alphabet_size, 0);
      WaveletTree wavelet_tree(alphabet_size, data, build_type);

      for (size_t i = 0; i <= data_size; i++) {
        uint64_t symb = query[i];
        uint64_t exp = count[symb];
        uint64_t act = wavelet_tree.rank(symb, i);
        EXPECT_EQ(act, exp);

        if (i == data_size) {
          break;
        }
        count[data[i]]++;
      }
    }
  }
}

TEST(WaveletTreeTest, SmokeSegment) {
  size_t data_size = 256, alphabet_size = 100;

  std::mt19937_64 rng(239);
  std::vector<uint64_t> data =
      generate_random_data(data_size, alphabet_size, rng);

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree wavelet_tree(alphabet_size, data, build_type);

    for (size_t begin = 0; begin <= data_size; begin++) {
      for (size_t end = begin; end <= data_size; end++) {
        auto segment = wavelet_tree.get_segment(begin, end);
        EXPECT_EQ(segment.size(), end - begin);
        for (size_t i = 0; i < end - begin; i++) {
          EXPECT_EQ(segment[i], data[begin + i]);
        }
      }
    }
  }
}

TEST(WaveletTreeTest, SerializationSmoke) {
  size_t data_size = 4096, alphabet_size = 100;

  std::mt19937_64 rng(239);
  std::vector<uint64_t> data =
      generate_random_data(data_size, alphabet_size, rng);

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree orig_tree(alphabet_size, data, build_type);

    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    orig_tree.serialize(writer);
    writer.finish();
    std::vector<std::byte> serialized_data = output.take();
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(serialized_data);
      auto view_tree = pixie::WaveletTreeView<std::uint64_t>::deserialize(
          reader, validation);
      EXPECT_TRUE(reader.empty());

      for (size_t i = 0; i <= data_size; i += 16) {
        uint64_t symb = data[i == data_size ? 0 : i];
        EXPECT_EQ(orig_tree.rank(symb, i), view_tree.rank(symb, i));
      }

      std::vector<size_t> count(alphabet_size, 0);
      for (auto symb : data) {
        count[symb]++;
      }

      for (uint64_t symb = 0; symb < alphabet_size; symb++) {
        for (uint64_t rank = 1; rank <= count[symb]; rank++) {
          EXPECT_EQ(orig_tree.select(symb, rank), view_tree.select(symb, rank));
        }
      }

      auto orig_segment = orig_tree.get_segment(0, data_size);
      auto view_segment = view_tree.get_segment(0, data_size);
      EXPECT_EQ(orig_segment, view_segment);
    }
  }
}

TEST(WaveletTreeTest, OwningRestorationSurvivesArtifactDestruction) {
  static_assert(std::same_as<
                pixie::WaveletTree<std::uint8_t>,
                pixie::WaveletTreeIndex<std::uint8_t, pixie::AlignedStorage>>);
  static_assert(
      std::same_as<
          pixie::WaveletTreeView<std::uint8_t>,
          pixie::WaveletTreeIndex<std::uint8_t, pixie::ReadOnlyStorageView>>);
  static_assert(pixie::Deserializable<pixie::WaveletTree<std::uint8_t>>);
  static_assert(pixie::Deserializable<pixie::WaveletTreeView<std::uint8_t>>);

  const std::vector<std::uint8_t> data = {3, 2, 0, 3, 1, 1, 2};
  std::optional<pixie::WaveletTree<std::uint8_t>> restored;
  std::vector<std::byte> canonical;
  {
    const pixie::WaveletTree<std::uint8_t> original(4, data);
    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    original.serialize(writer);
    writer.finish();
    canonical = output.take();
    pixie::BinaryReader reader(canonical);
    restored.emplace(pixie::WaveletTree<std::uint8_t>::deserialize(
        reader, pixie::DeserializationValidation::kFull));
    EXPECT_TRUE(reader.empty());
  }
  const std::vector<std::byte> expected = canonical;
  canonical.clear();
  canonical.shrink_to_fit();
  EXPECT_EQ(restored->get_segment(0, data.size()), data);

  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  restored->serialize(writer);
  writer.finish();
  EXPECT_EQ(output.take(), expected);
}

TEST(WaveletTreeTest, OwningRestorationAcceptsUnalignedArtifact) {
  const std::vector<std::uint8_t> data = {3, 2, 0, 3, 1, 1, 2};
  const pixie::WaveletTree<std::uint8_t> original(4, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();
  std::vector<std::byte> unaligned(artifact.size() + 1);
  std::ranges::copy(artifact, unaligned.begin() + 1);
  pixie::BinaryReader reader(std::span<const std::byte>(unaligned).subspan(1));
  const auto restored = pixie::WaveletTree<std::uint8_t>::deserialize(reader);
  EXPECT_EQ(restored.get_segment(0, data.size()), data);
}

TEST(WaveletTreeTest, SerializationAdvancesAcrossFramedArtifacts) {
  const std::vector<std::uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  const WaveletTree original(4, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifacts = output.take();
  pixie::BinaryReader reader(artifacts);

  const auto first = pixie::WaveletTreeView<std::uint64_t>::deserialize(reader);
  EXPECT_FALSE(reader.empty());
  const auto second =
      pixie::WaveletTreeView<std::uint64_t>::deserialize(reader);
  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(first.get_segment(0, data.size()), data);
  EXPECT_EQ(second.get_segment(0, data.size()), data);
}

TEST(WaveletTreeTest, SerializationRoundTripsAnEmptyTree) {
  const std::vector<std::uint64_t> data;
  const WaveletTree original(0, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();
  pixie::BinaryReader reader(artifact);
  const auto restored =
      pixie::WaveletTreeView<std::uint64_t>::deserialize(reader);

  EXPECT_TRUE(reader.empty());
  EXPECT_TRUE(restored.empty());
  EXPECT_EQ(restored.rank(0, 0), 0u);
  EXPECT_EQ(restored.select(0, 1), 0u);
  EXPECT_TRUE(restored.get_segment(0, 0).empty());
}

TEST(WaveletTreeTest, SerializesDirectlyToMappedFile) {
  const std::vector<std::uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  const WaveletTree original(4, data);
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_wavelet_serialization_test.bin";
  std::filesystem::remove(path);
  {
    pixie::io::FileOutputSink output(path);
    std::array<std::byte, 13> staging{};
    pixie::BinaryWriter writer(output, staging);
    original.serialize(writer);
    writer.finish();
  }

  pixie::io::MappedFile file(path);
  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader reader(file.as_bytes());
    const auto restored =
        pixie::WaveletTreeView<std::uint64_t>::deserialize(reader, validation);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(restored.get_segment(0, data.size()), data);
  }
  std::filesystem::remove(path);
}

TEST(WaveletTreeTest,
     FullValidationRejectsLeafAndChildLengthCorruptionTransactionally) {
  const std::vector<std::uint64_t> data = {0, 1, 2, 3, 4, 5, 6, 7};
  const WaveletTree original(8, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();
  const WaveletArtifactOffsets layout = locate_wavelet_artifact(valid);

  std::vector<std::byte> bad_leaf = valid;
  overwrite_u64(bad_leaf, layout.leaves, 3);
  pixie::BinaryReader leaf_quick_reader(bad_leaf);
  EXPECT_NO_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(
      leaf_quick_reader, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(leaf_quick_reader.empty());
  pixie::BinaryReader leaf_full_reader(bad_leaf);
  EXPECT_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(
                   leaf_full_reader, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(leaf_full_reader.position(), 0u);

  std::vector<std::byte> bad_child_length = valid;
  overwrite_u64(bad_child_length, layout.nodes[1].rank_num_bits, 3);
  pixie::BinaryReader child_full_reader(bad_child_length);
  EXPECT_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(
                   child_full_reader, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(child_full_reader.position(), 0u);
}

TEST(WaveletTreeTest, SerializationRejectsEveryTruncatedPrefixTransactionally) {
  const std::vector<std::uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  const WaveletTree original(4, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();

  for (std::size_t size = 0; size < artifact.size(); ++size) {
    SCOPED_TRACE(::testing::Message() << "size=" << size);
    pixie::BinaryReader reader(
        std::span<const std::byte>(artifact).first(size));
    EXPECT_THROW(
        (void)pixie::WaveletTreeView<std::uint64_t>::deserialize(reader),
        std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(WaveletTreeTest, SerializationRejectsUnalignedZeroCopyArtifacts) {
  const std::vector<std::uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  const WaveletTree original(4, data);

  pixie::VectorOutputSink unaligned_output;
  pixie::BinaryWriter unaligned_writer(unaligned_output);
  unaligned_writer.write_u8(0);
  EXPECT_THROW(original.serialize(unaligned_writer), std::invalid_argument);
  EXPECT_EQ(unaligned_writer.size_bytes(), 1u);

  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();
  std::vector<std::byte> unaligned_artifact(artifact.size() + 1);
  std::ranges::copy(artifact, unaligned_artifact.begin() + 1);
  pixie::BinaryReader reader(
      std::span<const std::byte>(unaligned_artifact).subspan(1));
  EXPECT_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(reader),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(WaveletTreeTest, SerializationRejectsMalformedTopologyTransactionally) {
  constexpr std::uint64_t kNoNode = std::numeric_limits<std::uint64_t>::max();
  const std::vector<std::uint64_t> data = {0, 1, 2, 3, 4, 5, 6, 7};
  const WaveletTree original(8, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();
  const auto layout = locate_wavelet_artifact(valid);
  const auto& nodes = layout.nodes;
  ASSERT_EQ(nodes.size(), 7u);

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(
                       reader, validation),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  };

  auto root_with_parent = valid;
  overwrite_u64(root_with_parent, nodes[0].parent, 0);
  expect_rejected(std::move(root_with_parent));

  auto self_loop = valid;
  overwrite_u64(self_loop, nodes[2].left_child, 2);
  expect_rejected(std::move(self_loop));

  auto inconsistent_parent = valid;
  overwrite_u64(inconsistent_parent, nodes[1].parent, kNoNode);
  expect_rejected(std::move(inconsistent_parent));

  auto duplicate_parent = valid;
  overwrite_u64(duplicate_parent, nodes[0].right_child, 1);
  expect_rejected(std::move(duplicate_parent));

  auto detached_node = valid;
  overwrite_u64(detached_node, nodes[0].left_child, kNoNode);
  overwrite_u64(detached_node, nodes[1].parent, kNoNode);
  expect_rejected(std::move(detached_node));

  auto disconnected_cycle = valid;
  overwrite_u64(disconnected_cycle, nodes[0].left_child, kNoNode);
  overwrite_u64(disconnected_cycle, nodes[1].parent, 2);
  overwrite_u64(disconnected_cycle, nodes[2].left_child, 1);
  expect_rejected(std::move(disconnected_cycle));

  auto invalid_middle = valid;
  overwrite_u64(invalid_middle, nodes[0].middle, 8);
  expect_rejected(std::move(invalid_middle));

  auto zero_middle = valid;
  overwrite_u64(zero_middle, nodes[0].middle, 0);
  expect_rejected(std::move(zero_middle));

  auto invalid_parent = valid;
  overwrite_u64(invalid_parent, nodes[1].parent, nodes.size());
  expect_rejected(std::move(invalid_parent));

  auto invalid_left_child = valid;
  overwrite_u64(invalid_left_child, nodes[0].left_child, nodes.size());
  expect_rejected(std::move(invalid_left_child));

  auto invalid_right_child = valid;
  overwrite_u64(invalid_right_child, nodes[0].right_child, nodes.size());
  expect_rejected(std::move(invalid_right_child));
}

TEST(WaveletTreeTest,
     SerializationRejectsMalformedFrameMetadataTransactionally) {
  constexpr std::size_t kReservedOffset =
      sizeof(std::uint64_t) + sizeof(std::uint32_t);
  const std::vector<std::uint64_t> data = {0, 1, 2, 3, 4, 5, 6, 7};
  const WaveletTree original(8, data);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();
  const auto layout = locate_wavelet_artifact(valid);
  ASSERT_EQ(layout.nodes.size(), 7u);

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW((void)pixie::WaveletTreeView<std::uint64_t>::deserialize(
                       reader, validation),
                   std::exception);
      EXPECT_EQ(reader.position(), 0u);
    }
  };

  auto incompatible_header = valid;
  incompatible_header[kReservedOffset] = std::byte{1};
  expect_rejected(std::move(incompatible_header));

  auto unpadded_size = valid;
  overwrite_u64(unpadded_size, 2 * sizeof(std::uint64_t),
                3 * sizeof(std::uint64_t) + 1);
  expect_rejected(std::move(unpadded_size));

  auto excessive_node_count = valid;
  overwrite_u64(excessive_node_count, layout.node_count,
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(excessive_node_count));

  auto truncated_nodes = valid;
  overwrite_u64(truncated_nodes, layout.node_count, 1000);
  expect_rejected(std::move(truncated_nodes));

  auto excessive_alphabet = valid;
  overwrite_u64(excessive_alphabet, layout.alphabet_size,
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(excessive_alphabet));

  auto duplicate_permutation = valid;
  overwrite_u64(duplicate_permutation, layout.permutation, 0);
  overwrite_u64(duplicate_permutation,
                layout.permutation + sizeof(std::uint64_t), 0);
  expect_rejected(std::move(duplicate_permutation));

  auto out_of_range_permutation = valid;
  overwrite_u64(out_of_range_permutation, layout.permutation, 8);
  expect_rejected(std::move(out_of_range_permutation));

  auto invalid_root = valid;
  overwrite_u64(invalid_root, layout.root, layout.nodes.size());
  expect_rejected(std::move(invalid_root));

  auto invalid_leaf = valid;
  overwrite_u64(invalid_leaf, layout.leaves, layout.nodes.size());
  expect_rejected(std::move(invalid_leaf));

  auto wrong_root_length = valid;
  overwrite_u64(wrong_root_length, layout.data_size, data.size() + 1);
  expect_rejected(std::move(wrong_root_length));
}
