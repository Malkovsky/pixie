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

void overwrite_u32(std::vector<std::byte>& bytes,
                   std::size_t offset,
                   std::uint32_t value) {
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
  std::size_t super_bytes;
  std::size_t basic_bytes;
  std::size_t sample_bytes;
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
  reader.skip(2 * sizeof(std::uint32_t));
  result.super_bytes = reader.read_size();
  reader.skip(result.super_bytes);
  result.basic_bytes = reader.read_size();
  reader.skip(result.basic_bytes);
  result.sample_bytes = reader.read_size();
  result.storage_offset = reader.position();
  reader.skip(result.sample_bytes);
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

TEST(RankSelectSupportTest, AcceptsPartialWordAlignedStorageSource) {
  pixie::AlignedStorage storage(1);
  storage.writable_bytes()[0] = std::byte{1};

  const pixie::RankSelectSupport support(storage, 1);

  EXPECT_EQ(support.size(), 1u);
  EXPECT_EQ(support.rank(1), 1u);
  EXPECT_EQ(support.select(1), 0u);

  const pixie::RankSelectSupport clamped(storage, 65);
  EXPECT_EQ(clamped.size(), 8u);
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

TEST(RankSelectSupportTest, SerializesOnlySimdRoundedBasicBlockMetadata) {
  constexpr std::array<std::size_t, 6> bit_counts = {
      0, 1, 32 * 512, 32 * 512 + 1, 65536, 65537};
  for (const std::size_t bit_count : bit_counts) {
    SCOPED_TRACE(bit_count);
    const std::vector<std::uint64_t> words((bit_count + 63) / 64);
    const pixie::RankSelectSupport<> support(
        words, bit_count, pixie::RankSelectSupport<>::SelectSupport::kNone);
    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    support.serialize(writer);
    writer.finish();
    const std::vector<std::byte> artifact = output.take();
    const RankSelectSampleLayout layout = locate_rank_select_samples(artifact);
    const std::size_t data_superblocks =
        bit_count == 0 ? 0 : 1 + (bit_count - 1) / 65536;
    const std::size_t data_basicblocks =
        bit_count == 0 ? 0 : 1 + (bit_count - 1) / 512;
    const std::size_t stored_basicblocks = (data_basicblocks + 31) / 32 * 32;
    EXPECT_EQ(layout.super_bytes,
              (data_superblocks + 1) * sizeof(std::uint64_t));
    EXPECT_EQ(layout.basic_bytes, stored_basicblocks * sizeof(std::uint16_t));
    EXPECT_EQ(layout.sample_bytes, 0u);
  }
}

TEST(RankSelectSupportTest, QuickAndFullValidationAcceptValidMetadata) {
  constexpr std::size_t kBitCount = 65537;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260816);
  for (std::uint64_t& word : words) {
    word = rng();
  }
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();

  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader reader(artifact);
    const auto restored =
        pixie::RankSelectSupport<>::deserialize(words, reader, validation);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(restored.rank(kBitCount), original.rank(kBitCount));
  }
}

TEST(RankSelectSupportTest,
     FullValidationRejectsSourceAndMetadataMismatchesTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();

  std::vector<std::uint64_t> different_words = words;
  different_words[0] ^= 1;
  pixie::BinaryReader quick_reader(valid);
  EXPECT_NO_THROW((void)pixie::RankSelectSupport<>::deserialize(
      different_words, quick_reader, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(quick_reader.empty());

  pixie::BinaryReader full_reader(valid);
  EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(
                   different_words, full_reader,
                   pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(full_reader.position(), 0u);

  std::vector<std::byte> bad_rank = valid;
  overwrite_u64(bad_rank, 2 * sizeof(std::uint64_t),
                original.rank(kBitCount) + 1);
  pixie::BinaryReader metadata_quick_reader(bad_rank);
  EXPECT_NO_THROW((void)pixie::RankSelectSupport<>::deserialize(
      words, metadata_quick_reader, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(metadata_quick_reader.empty());

  pixie::BinaryReader metadata_full_reader(bad_rank);
  EXPECT_THROW(
      (void)pixie::RankSelectSupport<>::deserialize(
          words, metadata_full_reader, pixie::DeserializationValidation::kFull),
      std::invalid_argument);
  EXPECT_EQ(metadata_full_reader.position(), 0u);
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
  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader reader(file.as_bytes());
    const auto restored =
        pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(
            words, reader, validation);
    EXPECT_TRUE(reader.empty());
    for (std::size_t position = 0; position <= kBitCount; position += 17) {
      EXPECT_EQ(restored.rank(position), original.rank(position));
      EXPECT_EQ(restored.rank0(position), original.rank0(position));
    }
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

    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader owning_reader(artifact);
      EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(
                       words, owning_reader, validation),
                   std::invalid_argument);
      EXPECT_EQ(owning_reader.position(), 0u);

      pixie::BinaryReader view_reader(artifact);
      EXPECT_THROW(
          (void)
              pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(
                  words, view_reader, validation),
          std::invalid_argument);
      EXPECT_EQ(view_reader.position(), 0u);
    }
  };

  expect_rejected(layout.select1_begin);
  expect_rejected(layout.select0_begin);
}

TEST(RankSelectSupportTest,
     DeserializationRejectsMalformedFixedMetadataTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(words, reader,
                                                                 validation),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  };

  auto bad_padded_size = valid;
  overwrite_u64(bad_padded_size, sizeof(std::uint64_t), 0);
  expect_rejected(std::move(bad_padded_size));

  auto excessive_rank = valid;
  overwrite_u64(excessive_rank, 2 * sizeof(std::uint64_t), kBitCount + 1);
  expect_rejected(std::move(excessive_rank));

  auto invalid_select1_begin = valid;
  overwrite_u64(invalid_select1_begin, 3 * sizeof(std::uint64_t),
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(invalid_select1_begin));

  auto missing_select1_samples = valid;
  overwrite_u64(missing_select1_samples, 4 * sizeof(std::uint64_t), 0);
  expect_rejected(std::move(missing_select1_samples));

  auto invalid_select0_begin = valid;
  overwrite_u64(invalid_select0_begin, 5 * sizeof(std::uint64_t),
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(invalid_select0_begin));

  auto missing_select0_samples = valid;
  overwrite_u64(missing_select0_samples, 6 * sizeof(std::uint64_t), 0);
  expect_rejected(std::move(missing_select0_samples));

  auto invalid_support = valid;
  overwrite_u32(invalid_support, 7 * sizeof(std::uint64_t), 4);
  expect_rejected(std::move(invalid_support));

  auto invalid_boolean = valid;
  overwrite_u32(invalid_boolean,
                7 * sizeof(std::uint64_t) + sizeof(std::uint32_t), 2);
  expect_rejected(std::move(invalid_boolean));
}

}  // namespace
