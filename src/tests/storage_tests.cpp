#include <gtest/gtest.h>
#include <pixie/io/mapped_file.h>
#include <pixie/storage/implementations.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <span>
#include <type_traits>
#include <utility>

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
  pixie::OutputBitStream owning_stream;
  pixie::OutputBitStream view_stream;
  storage.serialize(owning_stream);
  storage.view().serialize(view_stream);
  EXPECT_EQ(owning_stream.extract(), view_stream.extract());
}

TEST(StorageSerializationTest, ReadOnlyViewRoundTripsAndAdvancesInput) {
  pixie::AlignedStorage storage(1);
  storage.writable_bytes()[0] = std::byte{42};
  pixie::OutputBitStream stream;
  storage.serialize(stream);
  const auto serialized_words = stream.extract();
  std::span<const std::byte> input =
      std::as_bytes(std::span<const std::uint64_t>(serialized_words));
  const auto restored = pixie::ReadOnlyStorageView::deserialize(input);
  EXPECT_TRUE(std::ranges::equal(restored.as_bytes(), storage.as_bytes()));
  EXPECT_TRUE(input.empty());
}

TEST(AlignedStorageTest, PadsResizesAndProvidesWritableStorage) {
  pixie::AlignedStorage storage(1);
  EXPECT_EQ(storage.size_bytes(), pixie::kAlignedStorageLineBytes);
  EXPECT_EQ(storage.size_bits(), pixie::kAlignedStorageLineBits);
  EXPECT_EQ(reinterpret_cast<std::uintptr_t>(storage.as_bytes().data()) % 64,
            0u);
  storage.writable_words64()[0] = 42;
  EXPECT_EQ(storage.as_words64()[0], 42u);
  storage.resize(0);
  EXPECT_TRUE(storage.empty());
  EXPECT_GE(storage.allocated_bytes(), storage.size_bytes());
  storage.shrink_to_fit();
}

TEST(ReadOnlyStorageViewTest, MutatingOperationsAreNotAvailable) {
  static_assert(!HasWritableBytes<pixie::ReadOnlyStorageView>);
  static_assert(!HasResize<pixie::ReadOnlyStorageView>);
}

TEST(ReadOnlyStorageViewTest, DeserializeRejectsTruncatedInput) {
  std::array<std::byte, sizeof(std::size_t)> bytes{};
  const std::size_t payload_size = 1;
  std::memcpy(bytes.data(), &payload_size, sizeof(payload_size));
  std::span<const std::byte> input(bytes);
  EXPECT_THROW(pixie::ReadOnlyStorageView::deserialize(input),
               std::invalid_argument);
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
