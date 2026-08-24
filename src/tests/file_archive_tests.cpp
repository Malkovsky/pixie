#include <gtest/gtest.h>
#include <pixie/file_archive/implementations.h>
#include <pixie/serialization.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<std::byte> Bytes(std::string_view value) {
  const auto bytes = std::as_bytes(std::span(value.data(), value.size()));
  return {bytes.begin(), bytes.end()};
}

std::string String(std::span<const std::byte> bytes) {
  return {reinterpret_cast<const char*>(bytes.data()), bytes.size()};
}

std::vector<std::byte> Serialize(const pixie::FileArchive& archive) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  archive.serialize(writer);
  writer.finish();
  return sink.take();
}

void OverwriteU32(std::vector<std::byte>& bytes,
                  std::size_t offset,
                  std::uint32_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (8 * byte)) & 0xffU);
  }
}

void OverwriteU64(std::vector<std::byte>& bytes,
                  std::size_t offset,
                  std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (8 * byte)) & 0xffU);
  }
}

std::vector<pixie::FileArchiveSource> Sources() {
  return {{.path = "src/main.cpp",
           .content = Bytes("zero\none\n\ntwo"),
           .type = pixie::FileArchiveEntryType::kRegular,
           .mode = 0755},
          {.path = "assets/data.bin",
           .content = {std::byte{0xff}, std::byte{0}, std::byte{'\n'}},
           .type = pixie::FileArchiveEntryType::kRegular,
           .mode = 0644},
          {.path = "current",
           .content = Bytes("src/main.cpp"),
           .type = pixie::FileArchiveEntryType::kSymlink,
           .mode = 0777}};
}

}  // namespace

TEST(FileArchiveTest, FindsAndExtractsSortedEntries) {
  for (const auto build_type : {pixie::WaveletTreeBuildType::Standard,
                                pixie::WaveletTreeBuildType::Huffman}) {
    const pixie::FileArchive archive(Sources(), build_type);
    ASSERT_EQ(archive.size(), 3u);
    ASSERT_TRUE(archive.find("src/main.cpp").has_value());
    EXPECT_FALSE(archive.find("missing").has_value());

    const std::size_t source = *archive.find("src/main.cpp");
    const auto source_metadata = archive.entry(source);
    EXPECT_EQ(source_metadata.line_count, 4u);
    EXPECT_TRUE(source_metadata.is_text);
    EXPECT_EQ(source_metadata.mode, 0755u);
    EXPECT_EQ(String(archive.extract(source)), "zero\none\n\ntwo");
    EXPECT_EQ(String(archive.extract_lines(source, 0, 2)), "zero\none\n");
    EXPECT_EQ(String(archive.extract_lines(source, 1, 1)), "");
    EXPECT_EQ(String(archive.extract_lines(source, 2, 3)), "\n");
    EXPECT_EQ(String(archive.extract_lines(source, 4, 4)), "");

    const std::size_t binary = *archive.find("assets/data.bin");
    EXPECT_FALSE(archive.entry(binary).is_text);
    EXPECT_THROW(archive.extract_lines(binary, 0, 1), std::invalid_argument);
    EXPECT_EQ(archive.extract(binary), Sources()[1].content);

    const std::size_t symlink = *archive.find("current");
    EXPECT_EQ(archive.entry(symlink).type,
              pixie::FileArchiveEntryType::kSymlink);
    EXPECT_EQ(String(archive.extract(symlink)), "src/main.cpp");
    EXPECT_THROW(archive.extract_lines(symlink, 0, 0), std::invalid_argument);
  }
}

TEST(FileArchiveTest, SerializesOwningAndReadOnlyQueries) {
  for (const auto build_type : {pixie::WaveletTreeBuildType::Standard,
                                pixie::WaveletTreeBuildType::Huffman}) {
    const pixie::FileArchive archive(Sources(), build_type);
    const std::vector<std::byte> bytes = Serialize(archive);
    pixie::BinaryReader reader(bytes);
    const pixie::FileArchiveView quick = pixie::FileArchiveView::deserialize(
        reader, pixie::DeserializationValidation::kQuick);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(quick.build_type(), build_type);
    const std::size_t source = *quick.find("src/main.cpp");
    EXPECT_EQ(String(quick.extract_lines(source, 1, 4)), "one\n\ntwo");

    pixie::BinaryReader full_reader(bytes);
    const pixie::FileArchiveView full = pixie::FileArchiveView::deserialize(
        full_reader, pixie::DeserializationValidation::kFull);
    EXPECT_EQ(String(full.extract(*full.find("current"))), "src/main.cpp");
  }
}

TEST(FileArchiveTest, SupportsEmptyAndAllZeroContents) {
  std::vector<pixie::FileArchiveSource> sources = {
      {.path = "empty", .content = {}},
      {.path = "zero", .content = std::vector<std::byte>(3, std::byte{0})}};
  const pixie::FileArchive archive(std::move(sources));
  EXPECT_TRUE(archive.extract(*archive.find("empty")).empty());
  EXPECT_EQ(archive.extract(*archive.find("zero")),
            std::vector<std::byte>(3, std::byte{0}));
}

TEST(FileArchiveTest, RejectsInvalidSourcesAndRanges) {
  EXPECT_THROW(pixie::FileArchive({{.path = "", .content = {}}}),
               std::invalid_argument);
  EXPECT_THROW(pixie::FileArchive({{.path = "same", .content = {}},
                                   {.path = "same", .content = {}}}),
               std::invalid_argument);

  const pixie::FileArchive archive(Sources());
  const std::size_t source = *archive.find("src/main.cpp");
  EXPECT_THROW(archive.extract_lines(source, 3, 2), std::out_of_range);
  EXPECT_THROW(archive.extract_lines(source, 0, 5), std::out_of_range);
  EXPECT_THROW(archive.entry(archive.size()), std::out_of_range);
}

TEST(FileArchiveTest, ValidatesUtf8PathsAndEntryTypes) {
  const auto make_archive = [](std::string path) {
    return pixie::FileArchive({{.path = std::move(path), .content = {}}});
  };

  EXPECT_NO_THROW(make_archive("ascii-\xc2\xa2-\xe2\x82\xac-\xf0\x9f\x98\x80"));
  EXPECT_THROW(make_archive(std::string("nul\0path", 8)),
               std::invalid_argument);
  EXPECT_THROW(make_archive("\x80"), std::invalid_argument);
  EXPECT_THROW(make_archive("\xc2x"), std::invalid_argument);
  EXPECT_THROW(make_archive("\xe0\x80\x80"), std::invalid_argument);
  EXPECT_THROW(make_archive("\xed\xa0\x80"), std::invalid_argument);
  EXPECT_THROW(make_archive("\xf4\x90\x80\x80"), std::invalid_argument);
  EXPECT_THROW(make_archive("\xf0"), std::invalid_argument);
  EXPECT_THROW(pixie::FileArchive(
                   {{.path = "invalid-type",
                     .content = {},
                     .type = static_cast<pixie::FileArchiveEntryType>(2)}}),
               std::invalid_argument);
}

TEST(FileArchiveTest, RejectsMalformedMetadataTransactionally) {
  const pixie::FileArchive archive(Sources());
  const std::vector<std::byte> valid = Serialize(archive);
  const auto expect_rejected =
      [](std::vector<std::byte> bytes,
         pixie::DeserializationValidation validation =
             pixie::DeserializationValidation::kQuick) {
        pixie::BinaryReader reader(bytes);
        EXPECT_THROW(pixie::FileArchiveView::deserialize(reader, validation),
                     std::invalid_argument);
        EXPECT_EQ(reader.position(), 0u);
      };

  auto bad_version = valid;
  OverwriteU32(bad_version, 8, 5);
  expect_rejected(std::move(bad_version));

  auto bad_header_reserved = valid;
  OverwriteU32(bad_header_reserved, 12, 1);
  expect_rejected(std::move(bad_header_reserved));

  auto bad_build_type = valid;
  OverwriteU32(bad_build_type, 24, 2);
  expect_rejected(std::move(bad_build_type));

  auto bad_payload_reserved = valid;
  OverwriteU32(bad_payload_reserved, 28, 1);
  expect_rejected(std::move(bad_payload_reserved));

  auto bad_section_size = valid;
  OverwriteU64(bad_section_size, 48, 0);
  expect_rejected(std::move(bad_section_size));

  constexpr std::size_t kFirstRecord = 72;
  auto bad_type = valid;
  bad_type[kFirstRecord + 52] = std::byte{2};
  expect_rejected(std::move(bad_type));

  auto bad_text_flag = valid;
  bad_text_flag[kFirstRecord + 53] = std::byte{2};
  expect_rejected(std::move(bad_text_flag));

  auto bad_record_reserved = valid;
  bad_record_reserved[kFirstRecord + 54] = std::byte{1};
  expect_rejected(std::move(bad_record_reserved));

  const std::size_t padding_begin =
      72 + archive.file_table_bytes() + archive.path_storage_bytes();
  ASSERT_LT(padding_begin, archive.metadata_bytes());
  auto bad_padding = valid;
  bad_padding[padding_begin] = std::byte{1};
  expect_rejected(std::move(bad_padding));

  auto bad_logical_size = valid;
  OverwriteU64(bad_logical_size, 40, std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(bad_logical_size));

  auto bad_path_bounds = valid;
  OverwriteU64(bad_path_bounds, kFirstRecord + 8,
               std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(bad_path_bounds));

  auto empty_path = valid;
  OverwriteU64(empty_path, kFirstRecord + 8, 0);
  expect_rejected(std::move(empty_path),
                  pixie::DeserializationValidation::kFull);

  auto wrong_derived_text = valid;
  wrong_derived_text[kFirstRecord + 53] = std::byte{1};
  expect_rejected(std::move(wrong_derived_text),
                  pixie::DeserializationValidation::kFull);
}

TEST(FileArchiveTest, RejectsUnalignedSerializationAndDeserialization) {
  const pixie::FileArchive archive(Sources());
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  writer.write_u8(0);
  EXPECT_THROW(archive.serialize(writer), std::invalid_argument);

  const std::vector<std::byte> valid = Serialize(archive);
  std::vector<std::byte> unaligned(valid.size() + 1);
  std::ranges::copy(valid, unaligned.begin() + 1);
  pixie::BinaryReader reader(std::span<const std::byte>(unaligned).subspan(1));
  EXPECT_THROW(pixie::FileArchiveView::deserialize(reader),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(FileArchiveTest, RejectsTrailingAndTruncatedArtifacts) {
  const pixie::FileArchive archive(Sources());
  std::vector<std::byte> bytes = Serialize(archive);
  bytes.push_back(std::byte{0});
  pixie::BinaryReader trailing(bytes);
  const auto view = pixie::FileArchiveView::deserialize(trailing);
  EXPECT_EQ(view.size(), 3u);
  EXPECT_FALSE(trailing.empty());

  bytes.resize(16);
  pixie::BinaryReader truncated(bytes);
  EXPECT_THROW(pixie::FileArchiveView::deserialize(truncated),
               std::invalid_argument);
}
