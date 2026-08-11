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
#include <random>
#include <span>
#include <vector>

using pixie::WaveletTree;

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
    pixie::BinaryReader reader(serialized_data);
    auto view_tree = pixie::WaveletTreeView::deserialize(reader);
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

  const auto first = pixie::WaveletTreeView::deserialize(reader);
  EXPECT_FALSE(reader.empty());
  const auto second = pixie::WaveletTreeView::deserialize(reader);
  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(first.get_segment(0, data.size()), data);
  EXPECT_EQ(second.get_segment(0, data.size()), data);
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
  pixie::BinaryReader reader(file.as_bytes());
  const auto restored = pixie::WaveletTreeView::deserialize(reader);
  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(restored.get_segment(0, data.size()), data);
  std::filesystem::remove(path);
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
    EXPECT_THROW((void)pixie::WaveletTreeView::deserialize(reader),
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
  EXPECT_THROW((void)pixie::WaveletTreeView::deserialize(reader),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}
