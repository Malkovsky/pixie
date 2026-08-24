#include <gtest/gtest.h>
#include <pixie/file_archive/implementations.h>
#include <pixie/serialization.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
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
