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

template <class Storage>
concept HasMutableSegments = requires(Storage& value) {
  { value.segments() } -> std::same_as<pixie::SplitSpan<std::byte>>;
};

template <class Storage>
concept HasContiguousBytes = requires(const Storage& value) {
  { value.as_bytes() } -> std::same_as<std::span<const std::byte>>;
};

template <class Byte>
std::vector<std::byte> collect_segments(
    const pixie::SplitSpan<Byte>& segments) {
  std::vector<std::byte> result;
  result.reserve(segments.size());
  for (const std::span<Byte> segment : segments) {
    result.insert(result.end(), segment.begin(), segment.end());
  }
  return result;
}

TYPED_TEST(StorageSpecificationTest, EmptyStorageHasConsistentViews) {
  const std::array<std::byte, 0> bytes{};
  auto storage = make_storage<TypeParam>(bytes);
  EXPECT_TRUE(storage.empty());
  EXPECT_EQ(storage.size_bytes(), 0u);
  EXPECT_EQ(storage.size_bits(), 0u);
  EXPECT_EQ(storage.begin_position(), 0u);
  EXPECT_EQ(storage.end_position(), 0u);
  EXPECT_TRUE(storage.contains(0, 0));
  EXPECT_TRUE(storage.as_bytes().empty());
  EXPECT_TRUE(storage.view().empty());
  EXPECT_TRUE(storage.segments().empty());
  EXPECT_EQ(storage.segments().segment_count(), 0u);
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
  EXPECT_EQ(storage.begin_position(), 0u);
  EXPECT_EQ(storage.end_position(), bytes.size());
  EXPECT_TRUE(storage.contains(2, 4));
  EXPECT_FALSE(storage.contains(7, 2));
  const auto segments = storage.segments(2, 4);
  ASSERT_EQ(segments.segment_count(), 1u);
  EXPECT_EQ(segments.size(), 4u);
  EXPECT_TRUE(std::ranges::equal(*segments.begin(), middle.as_bytes()));
  EXPECT_THROW(storage.view(storage.size_bytes(), 1), std::out_of_range);
  EXPECT_THROW(storage.segments(storage.size_bytes(), 1), std::out_of_range);
}

TEST(SplitSpanTest, CanonicalizesAndIteratesTwoSegmentsInLogicalOrder) {
  const std::array first = {1, 2};
  const std::array second = {3, 4, 5};
  const pixie::SplitSpan<const int> split(first, second);

  EXPECT_EQ(split.size(), 5u);
  EXPECT_EQ(split.segment_count(), 2u);
  auto segment = split.begin();
  EXPECT_TRUE(std::ranges::equal(*segment++, first));
  EXPECT_TRUE(std::ranges::equal(*segment++, second));
  EXPECT_EQ(segment, split.end());

  const pixie::SplitSpan<const int> canonical({}, second);
  ASSERT_EQ(canonical.segment_count(), 1u);
  EXPECT_TRUE(std::ranges::equal(*canonical.begin(), second));
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

TEST(AlignedStorageTest, MutableSegmentsModifyTheBackingStorage) {
  static_assert(HasMutableSegments<pixie::AlignedStorage>);
  static_assert(
      std::same_as<
          decltype(std::declval<const pixie::AlignedStorage&>().segments()),
          pixie::SplitSpan<const std::byte>>);

  pixie::AlignedStorage storage(4 * 8);
  auto segments = storage.segments(1, 2);
  ASSERT_EQ(segments.segment_count(), 1u);
  (*segments.begin())[0] = std::byte{42};
  (*segments.begin())[1] = std::byte{43};

  EXPECT_EQ(storage.as_bytes()[1], std::byte{42});
  EXPECT_EQ(storage.as_bytes()[2], std::byte{43});
  EXPECT_THROW(storage.segments(storage.size_bytes(), 1), std::out_of_range);
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
  static_assert(!HasMutableSegments<pixie::ReadOnlyStorageView>);
  static_assert(
      std::same_as<
          decltype(std::declval<pixie::ReadOnlyStorageView&>().segments()),
          pixie::SplitSpan<const std::byte>>);
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

TEST(SlidingWindowStorageTest, AppendsAndExposesAbsolutePositions) {
  pixie::SlidingWindowStorage storage(8);
  const std::array bytes = {std::byte{1}, std::byte{2}, std::byte{3}};
  storage.append(bytes);

  EXPECT_EQ(storage.capacity_bytes(), 8u);
  EXPECT_EQ(storage.allocated_bytes(), 8u);
  EXPECT_EQ(storage.begin_position(), 0u);
  EXPECT_EQ(storage.end_position(), 3u);
  EXPECT_EQ(storage.size_bytes(), 3u);
  EXPECT_TRUE(storage.contains(0, 3));
  EXPECT_TRUE(storage.contains(3, 0));
  EXPECT_FALSE(storage.contains(3, 1));
  EXPECT_EQ(collect_segments(std::as_const(storage).segments()),
            std::vector<std::byte>(bytes.begin(), bytes.end()));
}

TEST(SlidingWindowStorageTest, WrapsEvictsAndRejectsExpiredRanges) {
  pixie::SlidingWindowStorage storage(5);
  const std::array first = {std::byte{0}, std::byte{1}, std::byte{2}};
  const std::array second = {std::byte{3}, std::byte{4}, std::byte{5},
                             std::byte{6}};
  storage.append(first);
  storage.append(second);

  EXPECT_EQ(storage.begin_position(), 2u);
  EXPECT_EQ(storage.end_position(), 7u);
  EXPECT_EQ(storage.size_bytes(), 5u);
  EXPECT_EQ(storage.segments().segment_count(), 2u);
  EXPECT_EQ(collect_segments(std::as_const(storage).segments()),
            (std::vector{std::byte{2}, std::byte{3}, std::byte{4}, std::byte{5},
                         std::byte{6}}));
  EXPECT_THROW(storage.segments(1, 1), std::out_of_range);
  EXPECT_THROW(storage.segments(6, 2), std::out_of_range);
  EXPECT_NO_THROW(storage.segments(storage.end_position(), 0));
}

TEST(SlidingWindowStorageTest, SupportsMutableRangesAcrossTheWrapPoint) {
  pixie::SlidingWindowStorage storage(5);
  const std::array bytes = {std::byte{0}, std::byte{1}, std::byte{2},
                            std::byte{3}, std::byte{4}, std::byte{5},
                            std::byte{6}};
  storage.append(bytes);

  auto segments = storage.segments(4, 3);
  ASSERT_EQ(segments.segment_count(), 2u);
  std::byte replacement = std::byte{40};
  for (const std::span<std::byte> segment : segments) {
    for (std::byte& value : segment) {
      value = replacement;
      replacement =
          static_cast<std::byte>(static_cast<unsigned char>(replacement) + 10);
    }
  }

  EXPECT_EQ(collect_segments(std::as_const(storage).segments()),
            (std::vector{std::byte{2}, std::byte{3}, std::byte{40},
                         std::byte{50}, std::byte{60}}));
}

TEST(SlidingWindowStorageTest, RetainedSegmentsKeepTheirAddresses) {
  pixie::SlidingWindowStorage storage(6);
  const std::array initial = {std::byte{0}, std::byte{1}, std::byte{2},
                              std::byte{3}, std::byte{4}, std::byte{5}};
  storage.append(initial);
  const auto retained = std::as_const(storage).segments(1, 2);
  const std::byte* const address = retained.begin()->data();

  const std::array next = {std::byte{6}};
  storage.append(next);

  ASSERT_TRUE(storage.contains(1, 2));
  EXPECT_EQ(retained.begin()->data(), address);
  EXPECT_EQ((*retained.begin())[0], std::byte{1});
  EXPECT_EQ((*retained.begin())[1], std::byte{2});
}

TEST(SlidingWindowStorageTest, OversizedAppendRetainsOnlyItsSuffix) {
  pixie::SlidingWindowStorage storage(4);
  const std::array bytes = {std::byte{0}, std::byte{1}, std::byte{2},
                            std::byte{3}, std::byte{4}, std::byte{5}};
  storage.append(bytes);

  EXPECT_EQ(storage.begin_position(), 2u);
  EXPECT_EQ(storage.end_position(), 6u);
  EXPECT_EQ(
      collect_segments(std::as_const(storage).segments()),
      (std::vector{std::byte{2}, std::byte{3}, std::byte{4}, std::byte{5}}));
}

TEST(SlidingWindowStorageTest, SerializesRetainedBytesInLogicalOrder) {
  pixie::SlidingWindowStorage storage(4);
  const std::array bytes = {std::byte{0}, std::byte{1}, std::byte{2},
                            std::byte{3}, std::byte{4}, std::byte{5}};
  storage.append(bytes);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  storage.serialize(writer);
  writer.finish();

  const std::vector<std::byte> artifact = output.take();
  pixie::BinaryReader reader(artifact);
  EXPECT_EQ(reader.read_size(), 4u);
  EXPECT_TRUE(std::ranges::equal(
      reader.read_bytes(4),
      std::array{std::byte{2}, std::byte{3}, std::byte{4}, std::byte{5}}));
  EXPECT_TRUE(reader.empty());
}

TEST(SlidingWindowStorageTest, HasFixedCapacityAndIndependentCopies) {
  static_assert(pixie::StorageImplementation<pixie::SlidingWindowStorage>);
  static_assert(HasMutableSegments<pixie::SlidingWindowStorage>);
  static_assert(!HasContiguousBytes<pixie::SlidingWindowStorage>);
  static_assert(!HasWritableBytes<pixie::SlidingWindowStorage>);
  static_assert(!HasResize<pixie::SlidingWindowStorage>);
  static_assert(std::is_copy_constructible_v<pixie::SlidingWindowStorage>);
  static_assert(std::is_move_constructible_v<pixie::SlidingWindowStorage>);
  static_assert(!std::is_copy_assignable_v<pixie::SlidingWindowStorage>);
  static_assert(!std::is_move_assignable_v<pixie::SlidingWindowStorage>);

  pixie::SlidingWindowStorage original(3);
  const std::array bytes = {std::byte{1}, std::byte{2}};
  original.append(bytes);
  pixie::SlidingWindowStorage copy(original);
  (*copy.segments(0, 1).begin())[0] = std::byte{9};

  EXPECT_EQ((*std::as_const(original).segments(0, 1).begin())[0], std::byte{1});
  EXPECT_EQ((*std::as_const(copy).segments(0, 1).begin())[0], std::byte{9});
}

TEST(SlidingWindowStorageTest, ZeroCapacityAcceptsOnlyEmptyAppends) {
  pixie::SlidingWindowStorage storage(0);
  const std::array<std::byte, 0> empty{};
  storage.append(empty);
  EXPECT_TRUE(storage.empty());
  EXPECT_TRUE(storage.segments().empty());

  const std::array nonempty = {std::byte{1}};
  EXPECT_THROW(storage.append(nonempty), std::length_error);
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
