#include <gtest/gtest.h>
#include <pixie/integer_vector/packed_monotone.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/serialization.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

static_assert(pixie::IntegerVectorValue<std::uint8_t>);
static_assert(pixie::IntegerVectorValue<std::uint64_t>);
static_assert(!pixie::IntegerVectorValue<bool>);
static_assert(!pixie::IntegerVectorValue<const bool>);
static_assert(!pixie::IntegerVectorValue<std::int64_t>);
static_assert(!pixie::IntegerVectorValue<const std::uint64_t>);

struct MinimalIntegerVector
    : pixie::IntegerVectorBase<MinimalIntegerVector, std::uint64_t> {
  std::size_t size_impl() const { return 0; }
  std::uint64_t value_at_impl(std::size_t) const { return 0; }
  void copy_to_impl(std::size_t, std::span<std::uint64_t>) const {}
};

template <class Vector>
concept HasMemoryUsage =
    requires(const Vector& vector) { vector.memory_usage_bytes(); };

static_assert(
    !HasMemoryUsage<pixie::MonotoneIntegerVector<MinimalIntegerVector>>);

using PackedOwner = pixie::PackedIntegerVector<>;
using PackedView = pixie::PackedIntegerVectorView<>;
using PackedMonotoneOwner = pixie::PackedMonotoneIntegerVector<>;
using PackedMonotoneView = pixie::PackedMonotoneIntegerVectorView<>;
using MinimalMonotone = pixie::MonotoneIntegerVector<MinimalIntegerVector>;

static_assert(pixie::Serializable<PackedOwner>);
static_assert(pixie::Serializable<PackedView>);
static_assert(pixie::Serializable<PackedMonotoneOwner>);
static_assert(pixie::Serializable<PackedMonotoneView>);
static_assert(pixie::Deserializable<PackedOwner>);
static_assert(pixie::Deserializable<PackedView>);
static_assert(pixie::Deserializable<PackedMonotoneOwner>);
static_assert(pixie::Deserializable<PackedMonotoneView>);
static_assert(!pixie::Serializable<MinimalMonotone>);
static_assert(!pixie::Deserializable<MinimalMonotone>);
static_assert(std::same_as<decltype(PackedOwner::deserialize(
                               std::declval<pixie::BinaryReader&>())),
                           PackedOwner>);
static_assert(std::same_as<decltype(PackedMonotoneView::deserialize(
                               std::declval<pixie::BinaryReader&>())),
                           PackedMonotoneView>);

template <class Serializable>
std::vector<std::byte> serialize_to_bytes(const Serializable& value) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  value.serialize(writer);
  writer.finish();
  return sink.take();
}

class AlignedArtifact {
 public:
  explicit AlignedArtifact(std::span<const std::byte> bytes)
      : words_((bytes.size() + sizeof(std::uint64_t) - 1) /
               sizeof(std::uint64_t)),
        size_(bytes.size()) {
    std::ranges::copy(bytes, writable_bytes().begin());
  }

  std::span<const std::byte> bytes() const {
    return std::as_bytes(std::span<const std::uint64_t>(words_)).first(size_);
  }

  std::span<std::byte> writable_bytes() {
    return std::as_writable_bytes(std::span<std::uint64_t>(words_))
        .first(size_);
  }

 private:
  std::vector<std::uint64_t> words_;
  std::size_t size_;
};

template <class Integer>
void overwrite_little_endian(std::span<std::byte> bytes,
                             std::size_t offset,
                             Integer value) {
  using Unsigned = std::make_unsigned_t<Integer>;
  const Unsigned encoded = static_cast<Unsigned>(value);
  ASSERT_LE(offset + sizeof(Integer), bytes.size());
  for (std::size_t i = 0; i < sizeof(Integer); ++i) {
    bytes[offset + i] =
        static_cast<std::byte>((encoded >> (i * 8)) & Unsigned{0xff});
  }
}

template <class Vector>
struct VectorHolder {
  std::shared_ptr<AlignedArtifact> backing;
  Vector vector;
};

struct PackedOwnerCase {
  using Vector = pixie::PackedIntegerVector<>;
  static constexpr bool kOwnsPayload = true;

  static VectorHolder<Vector> make(std::span<const std::uint64_t> values) {
    return {nullptr, Vector(values)};
  }

  static VectorHolder<Vector> make(std::span<const std::uint64_t> values,
                                   std::size_t width) {
    return {nullptr, Vector(values, width)};
  }
};

struct PackedViewCase {
  using Vector = pixie::PackedIntegerVectorView<>;
  static constexpr bool kOwnsPayload = false;

  static VectorHolder<Vector> make(std::span<const std::uint64_t> values) {
    return make_from_owner(pixie::PackedIntegerVector<>(values));
  }

  static VectorHolder<Vector> make(std::span<const std::uint64_t> values,
                                   std::size_t width) {
    return make_from_owner(pixie::PackedIntegerVector<>(values, width));
  }

 private:
  static VectorHolder<Vector> make_from_owner(
      const pixie::PackedIntegerVector<>& owner) {
    const std::vector<std::byte> serialized = serialize_to_bytes(owner);
    auto backing = std::make_shared<AlignedArtifact>(
        std::span<const std::byte>(serialized));
    pixie::BinaryReader reader(backing->bytes());
    Vector vector =
        Vector::deserialize(reader, pixie::DeserializationValidation::kFull);
    EXPECT_TRUE(reader.empty());
    return {std::move(backing), std::move(vector)};
  }
};

template <class Case>
class PackedIntegerVectorSpecificationTest : public ::testing::Test {};

using PackedIntegerVectorCases =
    ::testing::Types<PackedOwnerCase, PackedViewCase>;
TYPED_TEST_SUITE(PackedIntegerVectorSpecificationTest,
                 PackedIntegerVectorCases);

TYPED_TEST(PackedIntegerVectorSpecificationTest, EmptyAndAllZeroUseWidthZero) {
  const typename TypeParam::Vector default_vector;
  EXPECT_TRUE(default_vector.empty());
  EXPECT_EQ(default_vector.width(), 0u);

  const std::vector<std::uint64_t> empty;
  const auto empty_result = TypeParam::make(empty);
  EXPECT_TRUE(empty_result.vector.empty());
  EXPECT_EQ(empty_result.vector.width(), 0u);

  const std::array<std::uint64_t, 129> zeros{};
  const auto zero_result = TypeParam::make(zeros, 37);
  EXPECT_EQ(zero_result.vector.size(), zeros.size());
  EXPECT_EQ(zero_result.vector.width(), 0u);
  for (std::size_t i = 0; i < zeros.size(); ++i) {
    EXPECT_EQ(zero_result.vector[i], 0u);
  }
}

TYPED_TEST(PackedIntegerVectorSpecificationTest, SupportsEveryRuntimeWidth) {
  for (std::size_t width = 1; width <= 64; ++width) {
    const std::uint64_t maximum =
        width == 64 ? std::numeric_limits<std::uint64_t>::max()
                    : (std::uint64_t{1} << width) - 1;
    const std::array values = {std::uint64_t{0}, maximum, maximum / 3,
                               std::uint64_t{1}};
    const auto inferred = TypeParam::make(values);
    ASSERT_EQ(inferred.vector.width(), width);
    const auto result = TypeParam::make(values, width);
    ASSERT_EQ(result.vector.width(), width);
    ASSERT_EQ(result.vector.size(), values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      EXPECT_EQ(result.vector[i], values[i]);
      EXPECT_EQ(result.vector.at(i), values[i]);
    }
  }
}

TYPED_TEST(PackedIntegerVectorSpecificationTest,
           InfersWidthsAndCrossesWordBoundaries) {
  for (const std::size_t width : {1u, 3u, 7u, 17u, 31u, 63u, 64u}) {
    const std::uint64_t mask = width == 64
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : (std::uint64_t{1} << width) - 1;
    std::vector<std::uint64_t> values(517);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = (i * 0x9e3779b97f4a7c15ULL) & mask;
    }
    values.back() = mask;
    const auto result = TypeParam::make(values);
    EXPECT_EQ(result.vector.width(), width);
    for (std::size_t i = 0; i < values.size(); ++i) {
      ASSERT_EQ(result.vector[i], values[i]) << "width=" << width << " i=" << i;
    }
  }
}

TYPED_TEST(PackedIntegerVectorSpecificationTest,
           CopyToChecksCompleteRangeBeforeWriting) {
  const std::array values = {std::uint64_t{3}, std::uint64_t{1},
                             std::uint64_t{4}, std::uint64_t{1},
                             std::uint64_t{5}, std::uint64_t{9}};
  const auto result = TypeParam::make(values);

  std::array<std::uint64_t, 3> middle{};
  result.vector.copy_to(2, middle);
  EXPECT_EQ(middle, (std::array<std::uint64_t, 3>{4, 1, 5}));
  std::span<std::uint64_t> empty;
  result.vector.copy_to(result.vector.size(), empty);

  std::array<std::uint64_t, 2> unchanged = {77, 88};
  EXPECT_THROW(result.vector.copy_to(result.vector.size(), unchanged),
               std::out_of_range);
  EXPECT_EQ(unchanged, (std::array<std::uint64_t, 2>{77, 88}));
  EXPECT_THROW(result.vector.copy_to(result.vector.size() + 1, empty),
               std::out_of_range);
  EXPECT_THROW((void)result.vector.at(result.vector.size()), std::out_of_range);
}

TYPED_TEST(PackedIntegerVectorSpecificationTest, MatchesDeterministicOracle) {
  std::mt19937_64 random(42);
  for (std::size_t count : {1u, 2u, 63u, 64u, 65u, 511u, 512u, 513u}) {
    std::vector<std::uint64_t> values(count);
    std::generate(values.begin(), values.end(), [&] { return random(); });
    const auto result = TypeParam::make(values);
    std::vector<std::uint64_t> copy(count);
    result.vector.copy_to(0, copy);
    EXPECT_EQ(copy, values);
  }
}

TYPED_TEST(PackedIntegerVectorSpecificationTest, ReportsOnlyOwnedMemory) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2},
                             std::uint64_t{3}};
  const auto result = TypeParam::make(values);
  if constexpr (TypeParam::kOwnsPayload) {
    EXPECT_GE(result.vector.memory_usage_bytes(),
              sizeof(typename TypeParam::Vector) + 64u);
  } else {
    EXPECT_EQ(result.vector.memory_usage_bytes(),
              sizeof(typename TypeParam::Vector));
  }
}

TEST(PackedIntegerVectorTest, RejectsInvalidExplicitWidths) {
  const std::array values = {std::uint64_t{0}, std::uint64_t{4}};
  EXPECT_THROW((void)pixie::PackedIntegerVector<>(values, 2),
               std::invalid_argument);
  EXPECT_THROW((void)pixie::PackedIntegerVector<>(values, 65),
               std::invalid_argument);
  const std::array zeros = {std::uint64_t{0}, std::uint64_t{0}};
  EXPECT_NO_THROW((void)pixie::PackedIntegerVector<>(zeros, 0));
}

template <class Monotone>
void check_monotone_queries(const Monotone& vector,
                            const std::vector<std::uint64_t>& values) {
  const std::array probes = {std::uint64_t{0},
                             std::uint64_t{1},
                             std::uint64_t{2},
                             std::uint64_t{3},
                             std::uint64_t{7},
                             std::uint64_t{8},
                             std::numeric_limits<std::uint64_t>::max() - 1,
                             std::numeric_limits<std::uint64_t>::max()};
  for (const std::uint64_t probe : probes) {
    const std::size_t lower = static_cast<std::size_t>(
        std::lower_bound(values.begin(), values.end(), probe) - values.begin());
    const std::size_t upper = static_cast<std::size_t>(
        std::upper_bound(values.begin(), values.end(), probe) - values.begin());
    EXPECT_EQ(vector.lower_bound_index(probe), lower);
    EXPECT_EQ(vector.upper_bound_index(probe), upper);
    EXPECT_EQ(vector.contains(probe), lower != upper);
    EXPECT_EQ(vector.count(probe), upper - lower);
  }
}

TEST(MonotoneIntegerVectorTest, PreservesDuplicatesAndMatchesStandardBounds) {
  const std::vector<std::vector<std::uint64_t>> cases = {
      {},
      {7},
      {0, 1, 2, 3, 4, 5},
      {7, 7, 7, 7},
      {0, 1, 1, 1, 2, 7, 7, 8},
      {0, 1000, 1000000, std::numeric_limits<std::uint64_t>::max()},
  };
  for (const auto& values : cases) {
    const pixie::PackedMonotoneIntegerVector<> vector{
        pixie::PackedIntegerVector<>(values)};
    EXPECT_EQ(vector.size(), values.size());
    check_monotone_queries(vector, values);
  }
}

TEST(MonotoneIntegerVectorTest, RejectsEveryDecreasingShape) {
  for (const std::vector<std::uint64_t>& values :
       {std::vector<std::uint64_t>{1, 0}, {0, 2, 1}, {3, 3, 2}, {9, 1, 8}}) {
    EXPECT_THROW((void)pixie::PackedMonotoneIntegerVector<>(
                     pixie::PackedIntegerVector<>(values)),
                 std::invalid_argument);
  }
}

TEST(IntegerVectorSerializationTest, PackedOwnerAndViewRoundTrip) {
  const std::vector<std::uint64_t> values = {
      0, 1, 2, 3, 127, 128, std::numeric_limits<std::uint64_t>::max()};
  const pixie::PackedIntegerVector<> original(values);
  std::vector<std::byte> bytes = serialize_to_bytes(original);
  AlignedArtifact aligned(bytes);

  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader owner_reader(bytes);
    auto owner =
        pixie::PackedIntegerVector<>::deserialize(owner_reader, validation);
    EXPECT_TRUE(owner_reader.empty());
    pixie::BinaryReader view_reader(aligned.bytes());
    auto view =
        pixie::PackedIntegerVectorView<>::deserialize(view_reader, validation);
    EXPECT_TRUE(view_reader.empty());
    for (std::size_t i = 0; i < values.size(); ++i) {
      EXPECT_EQ(owner[i], values[i]);
      EXPECT_EQ(view[i], values[i]);
    }
    EXPECT_EQ(serialize_to_bytes(owner), bytes);
    EXPECT_EQ(serialize_to_bytes(view), bytes);
  }

  pixie::BinaryReader owner_reader(bytes);
  auto independent = pixie::PackedIntegerVector<>::deserialize(owner_reader);
  bytes.clear();
  bytes.shrink_to_fit();
  EXPECT_EQ(independent[values.size() - 1], values.back());
}

TEST(IntegerVectorSerializationTest, SupportsConcatenatedArtifactsAndSpanApi) {
  const std::array first_values = {std::uint64_t{1}, std::uint64_t{2}};
  const std::array second_values = {std::uint64_t{7}, std::uint64_t{9},
                                    std::uint64_t{11}};
  const auto first =
      serialize_to_bytes(pixie::PackedIntegerVector<>(first_values));
  const auto second =
      serialize_to_bytes(pixie::PackedIntegerVector<>(second_values));
  std::vector<std::byte> joined = first;
  joined.insert(joined.end(), second.begin(), second.end());

  std::span<const std::byte> remaining(joined);
  const auto restored_first =
      pixie::PackedIntegerVector<>::deserialize(remaining);
  EXPECT_EQ(remaining.size(), second.size());
  const auto restored_second =
      pixie::PackedIntegerVector<>::deserialize(remaining);
  EXPECT_TRUE(remaining.empty());
  EXPECT_EQ(restored_first[1], 2u);
  EXPECT_EQ(restored_second[2], 11u);
}

TEST(IntegerVectorSerializationTest, SpanApiRollsBackOnFailure) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{17},
                             std::uint64_t{1023}};
  std::vector<std::byte> truncated =
      serialize_to_bytes(pixie::PackedIntegerVector<>(values));
  truncated.pop_back();
  std::span<const std::byte> remaining(truncated);

  EXPECT_THROW((void)pixie::PackedIntegerVector<>::deserialize(
                   remaining, pixie::DeserializationValidation::kFull),
               std::exception);
  EXPECT_EQ(remaining.size(), truncated.size());
}

TEST(IntegerVectorSerializationTest, RejectsEveryTruncatedPackedPrefix) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{17},
                             std::uint64_t{1023}};
  const auto bytes = serialize_to_bytes(pixie::PackedIntegerVector<>(values));
  for (std::size_t length = 0; length < bytes.size(); ++length) {
    pixie::BinaryReader reader(std::span<const std::byte>(bytes).first(length));
    EXPECT_THROW((void)pixie::PackedIntegerVector<>::deserialize(reader),
                 std::exception)
        << "length=" << length;
    EXPECT_EQ(reader.position(), 0u) << "length=" << length;
  }
}

TEST(IntegerVectorSerializationTest, RejectsEveryTruncatedMonotonePrefix) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{1},
                             std::uint64_t{1023}};
  const pixie::PackedMonotoneIntegerVector<> original{
      pixie::PackedIntegerVector<>(values)};
  const auto bytes = serialize_to_bytes(original);
  for (std::size_t length = 0; length < bytes.size(); ++length) {
    pixie::BinaryReader reader(std::span<const std::byte>(bytes).first(length));
    EXPECT_THROW(
        (void)pixie::PackedMonotoneIntegerVector<>::deserialize(reader),
        std::exception)
        << "length=" << length;
    EXPECT_EQ(reader.position(), 0u) << "length=" << length;
  }
}

TEST(IntegerVectorSerializationTest,
     RejectsPackedHeaderAndDimensionCorruptionTransactionally) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2},
                             std::uint64_t{3}};
  const auto original =
      serialize_to_bytes(pixie::PackedIntegerVector<>(values));
  const auto expect_rejected = [&](std::vector<std::byte> bytes) {
    pixie::BinaryReader reader(bytes);
    EXPECT_THROW((void)pixie::PackedIntegerVector<>::deserialize(reader),
                 std::exception);
    EXPECT_EQ(reader.position(), 0u);
  };

  const std::array<std::function<void(std::span<std::byte>)>, 8> mutations = {
      [](std::span<std::byte> bytes) { bytes[0] ^= std::byte{1}; },
      [](std::span<std::byte> bytes) {
        overwrite_little_endian(bytes, 8, std::uint32_t{2});
      },
      [](std::span<std::byte> bytes) { bytes[12] = std::byte{2}; },
      [](std::span<std::byte> bytes) { bytes[13] = std::byte{32}; },
      [](std::span<std::byte> bytes) { bytes[14] = std::byte{65}; },
      [](std::span<std::byte> bytes) { bytes[15] = std::byte{1}; },
      [](std::span<std::byte> bytes) {
        overwrite_little_endian(bytes, 16, std::uint64_t{24});
      },
      [](std::span<std::byte> bytes) {
        overwrite_little_endian(bytes, 24,
                                std::numeric_limits<std::uint64_t>::max());
      },
  };
  for (const auto& mutation : mutations) {
    auto corrupted = original;
    mutation(corrupted);
    expect_rejected(std::move(corrupted));
  }

  auto trailing = original;
  trailing.resize(trailing.size() + sizeof(std::uint64_t));
  overwrite_little_endian(std::span<std::byte>(trailing), 16,
                          static_cast<std::uint64_t>(trailing.size()));
  expect_rejected(std::move(trailing));
}

TEST(IntegerVectorSerializationTest, FullValidationRejectsUnusedFinalBits) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2},
                             std::uint64_t{3}};
  auto bytes = serialize_to_bytes(pixie::PackedIntegerVector<>(values, 2));
  bytes.back() |= std::byte{0x80};

  pixie::BinaryReader quick_reader(bytes);
  EXPECT_NO_THROW((void)pixie::PackedIntegerVector<>::deserialize(
      quick_reader, pixie::DeserializationValidation::kQuick));
  pixie::BinaryReader full_reader(bytes);
  EXPECT_THROW((void)pixie::PackedIntegerVector<>::deserialize(
                   full_reader, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(full_reader.position(), 0u);
}

TEST(IntegerVectorSerializationTest, FullValidationRejectsNoncanonicalZeros) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  auto bytes = serialize_to_bytes(pixie::PackedIntegerVector<>(values, 2));
  overwrite_little_endian(std::span<std::byte>(bytes), 32, std::uint64_t{0});

  pixie::BinaryReader quick_reader(bytes);
  const auto quick = pixie::PackedIntegerVector<>::deserialize(
      quick_reader, pixie::DeserializationValidation::kQuick);
  EXPECT_EQ(quick.width(), 2u);
  EXPECT_EQ(quick[0], 0u);
  pixie::BinaryReader full_reader(bytes);
  EXPECT_THROW((void)pixie::PackedIntegerVector<>::deserialize(
                   full_reader, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(full_reader.position(), 0u);
}

TEST(IntegerVectorSerializationTest, PackedViewRequiresAlignedPayload) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  const auto bytes = serialize_to_bytes(pixie::PackedIntegerVector<>(values));
  std::vector<std::byte> unaligned(bytes.size() + 1);
  std::ranges::copy(bytes, unaligned.begin() + 1);
  const std::span<const std::byte> artifact(unaligned.data() + 1, bytes.size());

  pixie::BinaryReader view_reader(artifact);
  EXPECT_THROW((void)pixie::PackedIntegerVectorView<>::deserialize(view_reader),
               std::invalid_argument);
  EXPECT_EQ(view_reader.position(), 0u);
  pixie::BinaryReader owner_reader(artifact);
  EXPECT_NO_THROW(
      (void)pixie::PackedIntegerVector<>::deserialize(owner_reader));
}

TEST(IntegerVectorSerializationTest, MonotoneWrapperRoundTripsOwnerAndView) {
  const std::vector<std::uint64_t> values = {0, 1, 1, 7, 7, 7, 1000};
  const pixie::PackedMonotoneIntegerVector<> original{
      pixie::PackedIntegerVector<>(values)};
  const auto bytes = serialize_to_bytes(original);
  AlignedArtifact aligned(bytes);

  pixie::BinaryReader owner_reader(bytes);
  const auto owner = pixie::PackedMonotoneIntegerVector<>::deserialize(
      owner_reader, pixie::DeserializationValidation::kFull);
  pixie::BinaryReader view_reader(aligned.bytes());
  const auto view = pixie::PackedMonotoneIntegerVectorView<>::deserialize(
      view_reader, pixie::DeserializationValidation::kFull);
  EXPECT_TRUE(owner_reader.empty());
  EXPECT_TRUE(view_reader.empty());
  check_monotone_queries(owner, values);
  check_monotone_queries(view, values);
  EXPECT_EQ(serialize_to_bytes(owner), bytes);
  EXPECT_EQ(serialize_to_bytes(view), bytes);
}

TEST(IntegerVectorSerializationTest,
     QuickMonotoneValidationTrustsOrderingAndFullRejectsIt) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  const pixie::PackedMonotoneIntegerVector<> original{
      pixie::PackedIntegerVector<>(values, 2)};
  auto bytes = serialize_to_bytes(original);
  constexpr std::size_t kWrapperHeaderBytes = 24;
  constexpr std::size_t kPackedHeaderBytes = 32;
  overwrite_little_endian(std::span<std::byte>(bytes),
                          kWrapperHeaderBytes + kPackedHeaderBytes,
                          std::uint64_t{6});

  pixie::BinaryReader quick_reader(bytes);
  EXPECT_NO_THROW((void)pixie::PackedMonotoneIntegerVector<>::deserialize(
      quick_reader, pixie::DeserializationValidation::kQuick));
  pixie::BinaryReader full_reader(bytes);
  EXPECT_THROW((void)pixie::PackedMonotoneIntegerVector<>::deserialize(
                   full_reader, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(full_reader.position(), 0u);
}

TEST(IntegerVectorSerializationTest,
     MonotoneWrapperRejectsPlainPackedArtifact) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  const auto bytes = serialize_to_bytes(pixie::PackedIntegerVector<>(values));
  pixie::BinaryReader reader(bytes);
  EXPECT_THROW((void)pixie::PackedMonotoneIntegerVector<>::deserialize(reader),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(IntegerVectorSerializationTest,
     RejectsMonotoneWrapperCorruptionTransactionally) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  const pixie::PackedMonotoneIntegerVector<> original{
      pixie::PackedIntegerVector<>(values)};
  const auto serialized = serialize_to_bytes(original);
  const auto expect_rejected = [](std::vector<std::byte> bytes) {
    pixie::BinaryReader reader(bytes);
    EXPECT_THROW(
        (void)pixie::PackedMonotoneIntegerVector<>::deserialize(reader),
        std::exception);
    EXPECT_EQ(reader.position(), 0u);
  };

  auto bad_magic = serialized;
  bad_magic[0] ^= std::byte{1};
  expect_rejected(std::move(bad_magic));
  auto bad_version = serialized;
  overwrite_little_endian(std::span<std::byte>(bad_version), 8,
                          std::uint32_t{2});
  expect_rejected(std::move(bad_version));
  auto bad_reserved = serialized;
  overwrite_little_endian(std::span<std::byte>(bad_reserved), 12,
                          std::uint32_t{1});
  expect_rejected(std::move(bad_reserved));
  auto bad_size = serialized;
  overwrite_little_endian(std::span<std::byte>(bad_size), 16,
                          std::uint64_t{16});
  expect_rejected(std::move(bad_size));
  auto trailing = serialized;
  trailing.resize(trailing.size() + sizeof(std::uint64_t));
  overwrite_little_endian(std::span<std::byte>(trailing), 16,
                          static_cast<std::uint64_t>(trailing.size()));
  expect_rejected(std::move(trailing));
}

TEST(IntegerVectorSerializationTest,
     SerializersRejectUnalignedOffsetsBeforeWriting) {
  const std::array values = {std::uint64_t{1}, std::uint64_t{2}};
  const pixie::PackedIntegerVector<> packed(values);
  const pixie::PackedMonotoneIntegerVector<> monotone{
      pixie::PackedIntegerVector<>(values)};
  {
    pixie::VectorOutputSink sink;
    pixie::BinaryWriter writer(sink);
    writer.write_u8(0);
    EXPECT_THROW(packed.serialize(writer), std::invalid_argument);
    EXPECT_EQ(writer.size_bytes(), 1u);
  }
  {
    pixie::VectorOutputSink sink;
    pixie::BinaryWriter writer(sink);
    writer.write_u8(0);
    EXPECT_THROW(monotone.serialize(writer), std::invalid_argument);
    EXPECT_EQ(writer.size_bytes(), 1u);
  }
}

TEST(IntegerVectorSerializationTest, MappedViewRetainsMappedPayload) {
  const std::vector<std::uint64_t> values = {0, 1, 7, 31, 1000, 1000000};
  const pixie::PackedMonotoneIntegerVector<> original{
      pixie::PackedIntegerVector<>(values)};
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_integer_vector_mapped_test.bin";
  std::filesystem::remove(path);
  {
    pixie::io::FileOutputSink sink(path);
    std::array<std::byte, 13> staging{};
    pixie::BinaryWriter writer(sink, staging);
    original.serialize(writer);
    writer.finish();
  }

  pixie::io::MappedFile mapped(path);
  pixie::BinaryReader reader(mapped.as_bytes());
  const auto view = pixie::PackedMonotoneIntegerVectorView<>::deserialize(
      reader, pixie::DeserializationValidation::kFull);
  std::filesystem::remove(path);
  EXPECT_TRUE(reader.empty());
  check_monotone_queries(view, values);
}

}  // namespace
