#include <gtest/gtest.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/rank_select/implementations.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <random>
#include <span>
#include <string>
#include <vector>

namespace {

void overwrite_u64(std::vector<std::byte>& bytes,
                   std::size_t offset,
                   std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
  }
}

struct RankSelectSampleLayout {
  std::size_t select1_begin;
  std::size_t select1_count;
  std::size_t select0_begin;
  std::size_t select0_count;
  std::size_t storage_offset;
};

RankSelectSampleLayout locate_rank_select_samples(
    std::span<const std::byte> artifact) {
  pixie::BinaryReader reader(artifact);
  reader.skip(3 * sizeof(std::uint64_t));
  RankSelectSampleLayout result;
  result.select1_begin = reader.read_size();
  result.select1_count = reader.read_size();
  result.select0_begin = reader.read_size();
  result.select0_count = reader.read_size();
  reader.skip(2 * sizeof(std::uint32_t) + 8 * sizeof(std::uint64_t) +
              32 * sizeof(std::uint16_t));
  for (std::size_t storage = 0; storage < 2; ++storage) {
    reader.skip(reader.read_size());
  }
  const std::size_t sample_bytes = reader.read_size();
  result.storage_offset = reader.position();
  reader.skip(sample_bytes);
  return result;
}

template <class Support>
class RankSelectSpecificationTest : public testing::Test {};

using RankSelectImplementations = testing::Types<pixie::RankSelectSupport<>>;
TYPED_TEST_SUITE(RankSelectSpecificationTest, RankSelectImplementations);

TYPED_TEST(RankSelectSpecificationTest, EmptySequence) {
  const std::vector<std::uint64_t> words;
  const TypeParam bits(words, 0);

  EXPECT_TRUE(bits.empty());
  EXPECT_EQ(bits.size(), 0u);
  EXPECT_EQ(bits.rank(0), 0u);
  EXPECT_EQ(bits.rank0(0), 0u);
  EXPECT_EQ(bits.select(1), 0u);
  EXPECT_EQ(bits.select0(1), 0u);
  EXPECT_EQ(bits.to_string(), "");
}

TYPED_TEST(RankSelectSpecificationTest, RankSelectAndAccessMatchPackedBits) {
  constexpr std::size_t sizes[] = {1,    2,    63,    64,    65,    127,
                                   495,  496,  497,   511,   512,   513,
                                   1024, 4097, 65535, 65536, 65537, 131073};
  std::mt19937_64 rng(42);

  for (const std::size_t size : sizes) {
    std::vector<std::uint64_t> words((size + 63) / 64);
    for (auto& word : words) {
      word = rng();
    }
    const TypeParam bits(words, size);
    ASSERT_FALSE(bits.empty());
    ASSERT_EQ(bits.size(), size);
    ASSERT_TRUE(bits.supports_select1());
    ASSERT_TRUE(bits.supports_select0());

    std::size_t ones = 0;
    std::size_t zeros = 0;
    std::string expected;
    expected.reserve(size);
    for (std::size_t position = 0; position < size; ++position) {
      const bool value = ((words[position / 64] >> (position % 64)) & 1u) != 0;
      EXPECT_EQ(bits.rank(position), ones) << "size=" << size;
      EXPECT_EQ(bits.rank0(position), zeros) << "size=" << size;
      EXPECT_EQ(bits[position], value) << "size=" << size;
      expected.push_back(value ? '1' : '0');
      if (value) {
        ++ones;
        EXPECT_EQ(bits.select(ones), position) << "size=" << size;
      } else {
        ++zeros;
        EXPECT_EQ(bits.select0(zeros), position) << "size=" << size;
      }
    }
    EXPECT_EQ(bits.rank(size), ones);
    EXPECT_EQ(bits.rank0(size), zeros);
    EXPECT_EQ(bits.select(ones + 1), size);
    EXPECT_EQ(bits.select0(zeros + 1), size);
    EXPECT_EQ(bits.to_string(), expected);
  }
}

TEST(RankSelectSupportTest, AcceptsStorageSourceWithoutCopying) {
  pixie::AlignedStorage source(64);
  source.writable_words64()[0] = 0b101101;
  const pixie::RankSelectSupport<> support(source, 6);

  EXPECT_EQ(support.rank(6), 4u);
  EXPECT_EQ(support.select(4), 5u);
  source.writable_words64()[0] = 0b000001;
  EXPECT_EQ(support[2], 0);
}

TEST(RankSelectSupportTest, OwningMetadataDeserializationRoundTrips) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260718);
  for (std::uint64_t& word : words) {
    word = rng();
  }

  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  std::vector<std::byte> artifact = output.take();
  pixie::BinaryReader reader(artifact);
  const pixie::RankSelectSupport<> restored =
      pixie::RankSelectSupport<>::deserialize(words, reader);
  EXPECT_TRUE(reader.empty());

  artifact.clear();
  artifact.shrink_to_fit();
  for (std::size_t position = 0; position <= kBitCount; ++position) {
    EXPECT_EQ(restored.rank(position), original.rank(position));
    EXPECT_EQ(restored.rank0(position), original.rank0(position));
  }
  const std::size_t ones = original.rank(kBitCount);
  const std::size_t zeros = original.rank0(kBitCount);
  for (std::size_t rank = 1; rank <= ones + 1; ++rank) {
    EXPECT_EQ(restored.select(rank), original.select(rank));
  }
  for (std::size_t rank = 1; rank <= zeros + 1; ++rank) {
    EXPECT_EQ(restored.select0(rank), original.select0(rank));
  }
}

TEST(RankSelectSupportTest, SerializesDirectlyToMappedFile) {
  constexpr std::size_t kBitCount = 1025;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260721);
  for (std::uint64_t& word : words) {
    word = rng();
  }
  const pixie::RankSelectSupport<> original(words, kBitCount);
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_rank_select_serialization_test.bin";
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
  const auto restored =
      pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(words,
                                                                        reader);
  EXPECT_TRUE(reader.empty());
  for (std::size_t position = 0; position <= kBitCount; position += 17) {
    EXPECT_EQ(restored.rank(position), original.rank(position));
    EXPECT_EQ(restored.rank0(position), original.rank0(position));
  }
  std::filesystem::remove(path);
}

TEST(RankSelectSupportTest,
     DeserializationRejectsOutOfRangeSelectSamplesTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();
  const RankSelectSampleLayout layout = locate_rank_select_samples(valid);
  ASSERT_NE(layout.select1_count, 0u);
  ASSERT_NE(layout.select0_count, 0u);

  const auto expect_rejected = [&](std::size_t sample) {
    std::vector<std::byte> artifact = valid;
    overwrite_u64(artifact,
                  layout.storage_offset + sample * sizeof(std::uint64_t),
                  std::numeric_limits<std::uint64_t>::max());

    pixie::BinaryReader owning_reader(artifact);
    EXPECT_THROW(
        (void)pixie::RankSelectSupport<>::deserialize(words, owning_reader),
        std::invalid_argument);
    EXPECT_EQ(owning_reader.position(), 0u);

    pixie::BinaryReader view_reader(artifact);
    EXPECT_THROW(
        (void)pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(
            words, view_reader),
        std::invalid_argument);
    EXPECT_EQ(view_reader.position(), 0u);
  };

  expect_rejected(layout.select1_begin);
  expect_rejected(layout.select0_begin);
}

}  // namespace
