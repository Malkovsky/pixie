#include <gtest/gtest.h>
#include <pixie/io/mapped_file.h>
#include <pixie/storage/implementations.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <class Storage>
class StorageSpecificationTest : public testing::Test {};

using StorageTypes =
    testing::Types<pixie::AlignedStorage, pixie::ReadOnlyStorageView>;
TYPED_TEST_SUITE(StorageSpecificationTest, StorageTypes);

template <class Storage>
Storage make_storage(std::span<const std::byte> bytes) {
  if constexpr (std::same_as<Storage, pixie::AlignedStorage>) {
    Storage result(bytes.size() * 8);
    std::copy(bytes.begin(), bytes.end(), result.writable_bytes().begin());
    return result;
  } else {
    return Storage(bytes);
  }
}

template <class Storage>
concept HasWritableBytes = requires(Storage value) { value.writable_bytes(); };

template <class Storage>
concept HasResize = requires(Storage value) { value.resize(1); };

TYPED_TEST(StorageSpecificationTest, EmptyStorageHasConsistentViews) {
  const std::array<std::byte, 0> bytes{};
  auto storage = make_storage<TypeParam>(bytes);
  EXPECT_TRUE(storage.empty());
  EXPECT_EQ(storage.size_bytes(), 0u);
  EXPECT_EQ(storage.size_bits(), 0u);
  EXPECT_TRUE(storage.as_bytes().empty());
  EXPECT_TRUE(storage.view().empty());
}

TYPED_TEST(StorageSpecificationTest, SupportsByteAndNestedViews) {
  const std::array bytes = {std::byte{0}, std::byte{1}, std::byte{2},
                            std::byte{3}, std::byte{4}, std::byte{5},
                            std::byte{6}, std::byte{7}};
  auto storage = make_storage<TypeParam>(bytes);
  auto middle = storage.view(2, 4);
  auto nested = middle.view(1, 2);
  EXPECT_EQ(middle.as_bytes()[0], std::byte{2});
  EXPECT_EQ(nested.as_bytes()[0], std::byte{3});
  EXPECT_EQ(nested.as_bytes()[1], std::byte{4});
  EXPECT_THROW(storage.view(storage.size_bytes(), 1), std::out_of_range);
}

TYPED_TEST(StorageSpecificationTest, ProvidesAlignedWordViewsWhenValid) {
  alignas(8) const std::array<std::uint64_t, 8> words = {1, 2, 3, 4,
                                                         5, 6, 7, 8};
  const auto bytes = std::as_bytes(std::span(words));
  auto storage = make_storage<TypeParam>(bytes);
  EXPECT_EQ(storage.as_words16().size(), storage.size_bytes() / 2);
  EXPECT_EQ(storage.as_words64()[0], 1u);
  EXPECT_THROW(storage.view(1, 2).as_words16(), std::invalid_argument);
  EXPECT_THROW(storage.view(0, 3).as_words16(), std::invalid_argument);
}

TEST(StorageSerializationTest, OwningStorageAndViewSerializeIdentically) {
  pixie::AlignedStorage storage(1);
  storage.writable_bytes()[0] = std::byte{42};
  pixie::VectorOutputSink owning_output;
  pixie::VectorOutputSink view_output;
  pixie::BinaryWriter owning_writer(owning_output);
  pixie::BinaryWriter view_writer(view_output);
  storage.serialize(owning_writer);
  storage.view().serialize(view_writer);
  owning_writer.finish();
  view_writer.finish();
  const std::vector<std::byte> owning_bytes = owning_output.take();
  EXPECT_EQ(owning_bytes.size(), sizeof(std::uint64_t) + 1);
  EXPECT_EQ(owning_bytes, view_output.take());
}

TEST(StorageSerializationTest, ReadOnlyViewRoundTripsAndAdvancesInput) {
  pixie::AlignedStorage storage(1);
  storage.writable_bytes()[0] = std::byte{42};
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  storage.serialize(writer);
  writer.finish();
  const auto serialized_data = output.take();
  pixie::BinaryReader reader(serialized_data);
  const auto restored = pixie::ReadOnlyStorageView::deserialize(reader);
  EXPECT_TRUE(std::ranges::equal(restored.as_bytes(), storage.as_bytes()));
  EXPECT_TRUE(reader.empty());
}

TEST(StorageSerializationTest, OwningRestoreCopiesInputAndReturnsExactType) {
  static_assert(pixie::Serializable<pixie::AlignedStorage>);
  static_assert(pixie::Deserializable<pixie::AlignedStorage>);
  static_assert(pixie::Deserializable<pixie::ReadOnlyStorageView>);
  static_assert(std::same_as<decltype(pixie::AlignedStorage::deserialize(
                                 std::declval<pixie::BinaryReader&>())),
                             pixie::AlignedStorage>);

  pixie::AlignedStorage restored;
  {
    pixie::AlignedStorage original(24);
    original.writable_bytes()[0] = std::byte{1};
    original.writable_bytes()[1] = std::byte{2};
    original.writable_bytes()[2] = std::byte{3};
    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    original.serialize(writer);
    writer.finish();
    std::vector<std::byte> artifact = output.take();
    pixie::BinaryReader reader(artifact);
    restored = pixie::AlignedStorage::deserialize(reader);
    EXPECT_TRUE(reader.empty());
  }
  EXPECT_TRUE(
      std::ranges::equal(restored.as_bytes(),
                         std::array{std::byte{1}, std::byte{2}, std::byte{3}}));
}

TEST(StorageSerializationTest, BothRestorationsRollbackOnTruncation) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_size(1);
  writer.finish();
  const auto bytes = output.take();
  for (const auto restore :
       {+[](pixie::BinaryReader& reader) {
          (void)pixie::AlignedStorage::deserialize(reader);
        },
        +[](pixie::BinaryReader& reader) {
          (void)pixie::ReadOnlyStorageView::deserialize(reader);
        }}) {
    pixie::BinaryReader reader(bytes);
    EXPECT_THROW(restore(reader), std::invalid_argument);
    EXPECT_EQ(reader.position(), 0u);
  }
}

TEST(AlignedStorageTest, PadsResizesAndProvidesWritableStorage) {
  pixie::AlignedStorage storage(1);
  EXPECT_EQ(storage.size_bytes(), 1u);
  EXPECT_EQ(storage.logical_size_bytes(), 1u);
  EXPECT_EQ(storage.size_bits(), 8u);
  EXPECT_EQ(storage.padded_size_bytes(), pixie::kAlignedStorageLineBytes);
  EXPECT_EQ(storage.padded_view().size_bytes(),
            pixie::kAlignedStorageLineBytes);
  EXPECT_EQ(reinterpret_cast<std::uintptr_t>(storage.as_bytes().data()) % 64,
            0u);
  storage.writable_bytes()[0] = std::byte{42};
  EXPECT_EQ(storage.as_bytes()[0], std::byte{42});
  storage.resize(64);
  EXPECT_EQ(storage.size_bytes(), sizeof(std::uint64_t));
  storage.writable_words64()[0] = 42;
  EXPECT_EQ(storage.as_words64()[0], 42u);
  storage.resize(0);
  EXPECT_TRUE(storage.empty());
  EXPECT_TRUE(storage.padded_view().empty());
  EXPECT_GE(storage.allocated_bytes(), storage.size_bytes());
  storage.shrink_to_fit();
}

TEST(AlignedStorageTest, CopiesCompleteWordsIntoAlignedStorage) {
  std::array<std::uint64_t, 3> words = {1, 2, 3};
  const pixie::AlignedStorage storage{std::span<const std::uint64_t>(words)};
  words[0] = 4;

  EXPECT_EQ(storage.size_bytes(), 3 * sizeof(std::uint64_t));
  EXPECT_TRUE(std::ranges::equal(storage.as_words64(),
                                 std::array<std::uint64_t, 3>{1, 2, 3}));
  EXPECT_EQ(reinterpret_cast<std::uintptr_t>(storage.as_bytes().data()) % 64,
            0u);
}

TEST(ReadOnlyStorageViewTest, MutatingOperationsAreNotAvailable) {
  static_assert(!HasWritableBytes<pixie::ReadOnlyStorageView>);
  static_assert(!HasResize<pixie::ReadOnlyStorageView>);
}

TEST(ReadOnlyStorageViewTest, DeserializeRejectsTruncatedInput) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  writer.write_size(1);
  writer.finish();
  const auto bytes = output.take();
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW(pixie::ReadOnlyStorageView::deserialize(reader),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(MappedFileTest, MapsContentsAndIsMoveOnly) {
  static_assert(!std::is_copy_constructible_v<pixie::io::MappedFile>);
  const auto path =
      std::filesystem::temp_directory_path() / "pixie_mapped_file_test.bin";
  {
    std::ofstream output(path, std::ios::binary);
    output.write("pixie", 5);
  }
  pixie::io::MappedFile file(path);
  EXPECT_EQ(file.size_bytes(), 5u);
  EXPECT_EQ(file.as_bytes()[0], std::byte{'p'});
  std::filesystem::remove(path);
}

TEST(MappedFileTest, EmptyFileHasAnEmptyView) {
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_empty_mapped_file_test.bin";
  std::ofstream(path, std::ios::binary);

  pixie::io::MappedFile file(path);
  EXPECT_EQ(file.size_bytes(), 0u);
  EXPECT_TRUE(file.as_bytes().empty());
  std::filesystem::remove(path);
}

TEST(MappedFileTest, MoveTransfersMappingOwnership) {
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_mapped_file_move_test.bin";
  {
    std::ofstream output(path, std::ios::binary);
    output.write("pixie", 5);
  }

  pixie::io::MappedFile source(path);
  pixie::io::MappedFile moved(std::move(source));
  EXPECT_TRUE(source.as_bytes().empty());
  EXPECT_EQ(moved.as_bytes()[0], std::byte{'p'});

  pixie::io::MappedFile assigned(path);
  assigned = std::move(moved);
  EXPECT_TRUE(moved.as_bytes().empty());
  EXPECT_EQ(assigned.as_bytes()[0], std::byte{'p'});
  std::filesystem::remove(path);
}

TEST(MappedFileTest, RejectsMissingFile) {
  EXPECT_THROW(pixie::io::MappedFile(std::filesystem::temp_directory_path() /
                                     "pixie_missing_mapped_file.bin"),
               std::runtime_error);
}

}  // namespace
