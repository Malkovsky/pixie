#include <gtest/gtest.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/packed_bit_builder.h>
#include <pixie/serialization.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <vector>

namespace {

TEST(BinarySerializationTest, UsesCanonicalLittleEndianIntegerEncoding) {
  pixie::VectorOutputSink output;
  std::array<std::byte, 3> staging{};
  pixie::BinaryWriter writer(output, staging);
  writer.write_u8(0x12);
  writer.write_u16(0x3456);
  writer.write_u32(0x789abcde);
  writer.write_u64(0x0123456789abcdef);
  writer.write_i32(-2);
  writer.finish();

  const std::array expected = {
      std::byte{0x12}, std::byte{0x56}, std::byte{0x34}, std::byte{0xde},
      std::byte{0xbc}, std::byte{0x9a}, std::byte{0x78}, std::byte{0xef},
      std::byte{0xcd}, std::byte{0xab}, std::byte{0x89}, std::byte{0x67},
      std::byte{0x45}, std::byte{0x23}, std::byte{0x01}, std::byte{0xfe},
      std::byte{0xff}, std::byte{0xff}, std::byte{0xff}};
  EXPECT_TRUE(std::ranges::equal(output.bytes(), expected));

  pixie::BinaryReader reader(output.bytes());
  EXPECT_EQ(reader.read_u8(), 0x12);
  EXPECT_EQ(reader.read_u16(), 0x3456);
  EXPECT_EQ(reader.read_u32(), 0x789abcdeu);
  EXPECT_EQ(reader.read_u64(), 0x0123456789abcdefu);
  EXPECT_EQ(reader.read_i32(), -2);
  EXPECT_TRUE(reader.empty());
}

TEST(BinarySerializationTest, RoundTripsEverySignedIntegerWidth) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_i8(-2);
  writer.write_i16(-3);
  writer.write_i32(-4);
  writer.write_i64(-5);
  writer.finish();

  pixie::BinaryReader reader(output.bytes());
  EXPECT_EQ(reader.read_i8(), -2);
  EXPECT_EQ(reader.read_i16(), -3);
  EXPECT_EQ(reader.read_i32(), -4);
  EXPECT_EQ(reader.read_i64(), -5);
  EXPECT_TRUE(reader.empty());
}

TEST(BinarySerializationTest, WritesBytesFramesAndPadding) {
  const std::array payload = {std::byte{1}, std::byte{2}, std::byte{3}};
  std::array<std::byte, 16> storage{};
  std::array<std::byte, 5> staging{};
  pixie::SpanOutputSink output(storage);
  pixie::BinaryWriter writer(output, staging);
  const std::size_t size_position = writer.write_u64_placeholder();
  writer.write_bytes(payload);
  writer.align_to(8);
  writer.patch_u64(size_position, writer.size_bytes());
  writer.finish();

  ASSERT_EQ(writer.size_bytes(), 16u);
  EXPECT_EQ(output.size_bytes(), 16u);
  pixie::BinaryReader reader(output.bytes());
  EXPECT_EQ(reader.read_u64(), 16u);
  EXPECT_TRUE(std::ranges::equal(reader.read_bytes(payload.size()), payload));
  reader.require_zero_padding(5);
  EXPECT_TRUE(reader.empty());
}

TEST(BinarySerializationTest, BackpatchesAFieldSplitByABufferBoundary) {
  pixie::VectorOutputSink output;
  std::array<std::byte, 7> staging{};
  pixie::BinaryWriter writer(output, staging);
  writer.write_u8(0xaa);
  const std::size_t size_position = writer.write_u64_placeholder();
  writer.write_u16(0x1234);
  writer.patch_u64(size_position, writer.size_bytes());
  writer.finish();

  pixie::BinaryReader reader(output.bytes());
  EXPECT_EQ(reader.read_u8(), 0xaa);
  EXPECT_EQ(reader.read_u64(), 11u);
  EXPECT_EQ(reader.read_u16(), 0x1234);
  EXPECT_TRUE(reader.empty());
}

TEST(BinarySerializationTest, ExplicitVectorSinkCanBeTakenAndReused) {
  pixie::VectorOutputSink output;
  {
    pixie::BinaryWriter writer(output);
    writer.write_u32(42);
    writer.finish();
  }
  const std::vector<std::byte> first = output.take();
  EXPECT_EQ(first.size(), sizeof(std::uint32_t));
  EXPECT_TRUE(output.empty());

  {
    pixie::BinaryWriter writer(output);
    writer.write_u8(7);
    writer.finish();
  }
  const std::vector<std::byte> second = output.take();
  ASSERT_EQ(second.size(), 1u);
  EXPECT_EQ(second[0], std::byte{7});
}

TEST(BinarySerializationTest, ReaderReportsOffsetsWithoutPartialPrimitiveRead) {
  const std::array bytes = {std::byte{0}, std::byte{0}, std::byte{0},
                            std::byte{0}, std::byte{1}, std::byte{2}};
  pixie::BinaryReader reader(bytes);
  reader.skip(4);
  try {
    (void)reader.read_u32();
    FAIL() << "Expected a truncated-input error";
  } catch (const pixie::SerializationError& error) {
    EXPECT_EQ(error.byte_offset(), 4u);
  }
  EXPECT_EQ(reader.position(), 4u);
}

TEST(BinarySerializationTest, SubreaderIsBoundedAndTracksAbsoluteOffsets) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_u32(1);
  writer.write_u16(2);
  writer.write_u16(3);
  writer.finish();
  pixie::BinaryReader reader(output.bytes());
  reader.skip(sizeof(std::uint32_t));
  pixie::BinaryReader child = reader.read_subreader(sizeof(std::uint16_t));

  EXPECT_EQ(child.read_u16(), 2);
  EXPECT_TRUE(child.empty());
  EXPECT_EQ(reader.position(), 6u);
  EXPECT_EQ(reader.read_u16(), 3);
  EXPECT_TRUE(reader.empty());
}

TEST(BinarySerializationTest, RejectsBadPaddingAndInvalidPatch) {
  const std::array bytes = {std::byte{0}, std::byte{1}};
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW(reader.require_zero_padding(2), pixie::SerializationError);
  EXPECT_EQ(reader.position(), 0u);

  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_u32(0);
  EXPECT_THROW(writer.patch_u64(0, 1), std::out_of_range);
  EXPECT_THROW(writer.align_to(0), std::invalid_argument);
  writer.finish();
}

TEST(BinarySerializationTest, RejectsExcessivePaddingTransactionally) {
  const std::array bytes = {std::byte{0}, std::byte{0}};
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW(reader.require_zero_padding(1), pixie::SerializationError);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(BinarySerializationTest, OutputSinksCheckDirectWritesAndPatches) {
  const std::array initial = {std::byte{1}, std::byte{2}, std::byte{3}};
  const std::array patch = {std::byte{7}, std::byte{8}};

  pixie::VectorOutputSink vector_output;
  vector_output.reserve_bytes(16);
  vector_output.write(initial);
  vector_output.write_at(1, patch);
  vector_output.write_at(vector_output.size_bytes(), {});
  EXPECT_THROW(vector_output.write_at(2, patch), std::out_of_range);
  EXPECT_THROW(vector_output.write_at(4, {}), std::out_of_range);
  vector_output.finish();
  EXPECT_EQ(vector_output.bytes()[0], std::byte{1});
  EXPECT_EQ(vector_output.bytes()[1], std::byte{7});
  EXPECT_EQ(vector_output.bytes()[2], std::byte{8});

  std::array<std::byte, 4> storage{};
  pixie::SpanOutputSink span_output(storage);
  EXPECT_TRUE(span_output.empty());
  EXPECT_EQ(span_output.capacity_bytes(), storage.size());
  span_output.write(initial);
  span_output.write_at(1, patch);
  EXPECT_THROW(span_output.write(patch), std::length_error);
  EXPECT_THROW(span_output.write_at(2, patch), std::out_of_range);
  span_output.finish();
  EXPECT_EQ(span_output.bytes()[1], std::byte{7});
  EXPECT_EQ(span_output.bytes()[2], std::byte{8});
}

TEST(BinarySerializationTest, RejectsEmptyStagingBuffers) {
  pixie::VectorOutputSink output;
  EXPECT_THROW((void)pixie::BinaryWriter(output, std::size_t{0}),
               std::invalid_argument);
  EXPECT_THROW((void)pixie::BinaryWriter(output, std::span<std::byte>{}),
               std::invalid_argument);
}

TEST(BinarySerializationTest, FixedSpanReportsCapacityExhaustion) {
  std::array<std::byte, sizeof(std::uint32_t)> storage{};
  std::array<std::byte, 2> staging{};
  pixie::SpanOutputSink output(storage);
  pixie::BinaryWriter writer(output, staging);
  writer.write_u32(42);
  writer.write_u8(7);
  EXPECT_THROW(writer.finish(), std::length_error);
  EXPECT_EQ(output.size_bytes(), sizeof(std::uint32_t));
  EXPECT_THROW(writer.write_u8(1), std::logic_error);
}

class RecordingOutputSink {
 public:
  void write(std::span<const std::byte> bytes) {
    maximum_write_size = std::max(maximum_write_size, bytes.size());
    size_bytes += bytes.size();
  }

  void write_at(std::size_t position, std::span<const std::byte> bytes) {
    if (position > size_bytes || bytes.size() > size_bytes - position) {
      throw std::out_of_range("Recording sink patch is outside the output");
    }
  }

  void finish() { finished = true; }

  std::size_t size_bytes = 0;
  std::size_t maximum_write_size = 0;
  bool finished = false;
};

static_assert(pixie::SeekableBinaryOutputSink<RecordingOutputSink>);

TEST(BinarySerializationTest, UsesOnlyTheCallerProvidedBufferForGeneratedData) {
  RecordingOutputSink output;
  std::array<std::byte, 7> staging{};
  pixie::BinaryWriter writer(output, staging);
  writer.write_zeros(1024 * 1024);
  writer.finish();

  EXPECT_EQ(writer.buffer_size_bytes(), staging.size());
  EXPECT_EQ(output.size_bytes, 1024u * 1024u);
  EXPECT_LE(output.maximum_write_size, staging.size());
  EXPECT_TRUE(output.finished);
}

class FailingOutputSink {
 public:
  void write(std::span<const std::byte>) {
    throw std::runtime_error("injected write failure");
  }

  void write_at(std::size_t, std::span<const std::byte>) {
    throw std::runtime_error("injected patch failure");
  }

  void finish() {}
};

class PatchFailingOutputSink {
 public:
  void write(std::span<const std::byte> bytes) { size_bytes += bytes.size(); }

  void write_at(std::size_t, std::span<const std::byte>) {
    throw std::runtime_error("injected patch failure");
  }

  void finish() {}

  std::size_t size_bytes = 0;
};

class FinishFailingOutputSink {
 public:
  void write(std::span<const std::byte> bytes) { size_bytes += bytes.size(); }

  void write_at(std::size_t, std::span<const std::byte>) {}

  void finish() { throw std::runtime_error("injected finish failure"); }

  std::size_t size_bytes = 0;
};

TEST(BinarySerializationTest, BecomesUnusableAfterSinkFailure) {
  FailingOutputSink output;
  std::array<std::byte, 1> staging{};
  pixie::BinaryWriter writer(output, staging);
  EXPECT_THROW(writer.write_u8(1), std::runtime_error);
  EXPECT_THROW(writer.write_u8(2), std::logic_error);
  EXPECT_THROW(writer.finish(), std::logic_error);
}

TEST(BinarySerializationTest, BecomesUnusableAfterPatchOrFinishFailure) {
  PatchFailingOutputSink patch_output;
  std::array<std::byte, 4> patch_staging{};
  pixie::BinaryWriter patch_writer(patch_output, patch_staging);
  patch_writer.write_u64(0);
  ASSERT_EQ(patch_output.size_bytes, sizeof(std::uint64_t));
  EXPECT_THROW(patch_writer.patch_u64(0, 1), std::runtime_error);
  EXPECT_THROW(patch_writer.flush(), std::logic_error);

  FinishFailingOutputSink finish_output;
  pixie::BinaryWriter finish_writer(finish_output);
  finish_writer.write_u8(1);
  EXPECT_THROW(finish_writer.finish(), std::runtime_error);
  EXPECT_EQ(finish_output.size_bytes, 1u);
  EXPECT_THROW(finish_writer.finish(), std::logic_error);
}

TEST(BinarySerializationTest, FlushDeliversBufferedDataAndKeepsWriterOpen) {
  pixie::VectorOutputSink output;
  std::array<std::byte, 4> staging{};
  pixie::BinaryWriter writer(output, staging);
  EXPECT_TRUE(writer.empty());
  writer.write_u16(0x1234);
  EXPECT_TRUE(output.empty());
  writer.flush();
  EXPECT_EQ(output.size_bytes(), sizeof(std::uint16_t));
  writer.flush();
  writer.write_u8(5);
  writer.finish();
  EXPECT_EQ(output.size_bytes(), 3u);
}

TEST(BinarySerializationTest, FinishIsIdempotentAndClosesTheWriter) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_u8(1);
  writer.finish();
  EXPECT_NO_THROW(writer.finish());
  EXPECT_THROW(writer.write_u8(2), std::logic_error);
}

TEST(BinarySerializationTest, WritesAndBackpatchesAFileWithBoundedMemory) {
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_file_output_sink_test.bin";
  std::filesystem::remove(path);
  {
    pixie::io::FileOutputSink output(path);
    std::array<std::byte, 5> staging{};
    pixie::BinaryWriter writer(output, staging);
    const std::size_t size_position = writer.write_u64_placeholder();
    writer.write_u32(0x12345678);
    writer.align_to(sizeof(std::uint64_t));
    writer.patch_u64(size_position, writer.size_bytes());
    writer.finish();
    EXPECT_EQ(output.size_bytes(), 16u);
    EXPECT_THROW(output.write(std::span<const std::byte>{}), std::logic_error);
  }

  pixie::io::MappedFile file(path);
  pixie::BinaryReader reader(file.as_bytes());
  EXPECT_EQ(reader.read_u64(), 16u);
  EXPECT_EQ(reader.read_u32(), 0x12345678u);
  reader.require_zero_padding(4);
  std::filesystem::remove(path);
}

TEST(BinarySerializationTest, FileSinkSupportsMovesAndChecksItsState) {
  const auto first_path = std::filesystem::temp_directory_path() /
                          "pixie_file_output_sink_move_test.bin";
  const auto second_path = std::filesystem::temp_directory_path() /
                           "pixie_file_output_sink_assignment_test.bin";
  std::filesystem::remove(first_path);
  std::filesystem::remove(second_path);
  const std::array bytes = {std::byte{1}, std::byte{2}, std::byte{3}};
  const std::array patch = {std::byte{9}};
  {
    pixie::io::FileOutputSink first(first_path);
    first.write(bytes);
    pixie::io::FileOutputSink moved(std::move(first));
    EXPECT_THROW(first.write({}), std::logic_error);

    pixie::io::FileOutputSink assigned(second_path);
    assigned = std::move(moved);
    EXPECT_EQ(assigned.size_bytes(), bytes.size());
    assigned.write_at(1, patch);
    EXPECT_THROW(assigned.write_at(3, patch), std::out_of_range);
    assigned.finish();
    EXPECT_NO_THROW(assigned.finish());
    EXPECT_THROW(assigned.write({}), std::logic_error);
  }

  pixie::io::MappedFile file(first_path);
  ASSERT_EQ(file.as_bytes().size(), bytes.size());
  EXPECT_EQ(file.as_bytes()[1], std::byte{9});
  std::filesystem::remove(first_path);
  std::filesystem::remove(second_path);
}

TEST(BinarySerializationTest, FileSinkReportsOpenFailure) {
  const auto missing_directory =
      std::filesystem::temp_directory_path() / "pixie_missing_output_directory";
  std::filesystem::remove_all(missing_directory);
  EXPECT_THROW(
      (void)pixie::io::FileOutputSink(missing_directory / "artifact.bin"),
      std::system_error);
}

bool packed_bit(std::span<const std::uint64_t> words, std::size_t position) {
  return ((words[position / 64] >> (position % 64)) & 1u) != 0;
}

TEST(PackedBitBuilderTest, WritesAcrossEveryWordOffset) {
  constexpr std::uint64_t kPattern = 0xfedcba9876543210;
  for (std::size_t offset = 0; offset < 64; ++offset) {
    for (const std::size_t width : std::array<std::size_t, 4>{8, 16, 32, 64}) {
      SCOPED_TRACE(::testing::Message()
                   << "offset=" << offset << ", width=" << width);
      pixie::PackedBitBuilder builder;
      for (std::size_t bit = 0; bit < offset; ++bit) {
        builder.write_bit(false);
      }
      builder.write_bits(kPattern, width);
      EXPECT_EQ(builder.size_bits(), offset + width);
      const std::vector<std::uint64_t> words = builder.take_words();
      for (std::size_t bit = 0; bit < offset; ++bit) {
        EXPECT_FALSE(packed_bit(words, bit));
      }
      for (std::size_t bit = 0; bit < width; ++bit) {
        EXPECT_EQ(packed_bit(words, offset + bit),
                  ((kPattern >> bit) & 1u) != 0);
      }
      EXPECT_EQ(builder.size_bits(), 0u);
    }
  }
}

TEST(PackedBitBuilderTest, SupportsPartialFieldsAndSafeReuse) {
  pixie::PackedBitBuilder builder;
  builder.reserve_bits(129);
  builder.write_bits(0xff, 4);
  builder.write_bit(false);
  const auto first = builder.take_words();
  ASSERT_EQ(first.size(), 1u);
  EXPECT_EQ(first[0], 0xfu);

  builder.write_bit(true);
  builder.write_bits(42, 0);
  const auto second = builder.take_words();
  ASSERT_EQ(second.size(), 1u);
  EXPECT_EQ(second[0], 1u);
  EXPECT_THROW(builder.write_bits(0, 65), std::invalid_argument);
}

}  // namespace
