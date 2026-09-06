#include <gtest/gtest.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/rmq/cartesian_hybrid_btree.h>
#include <pixie/rmq/cartesian_rmm.h>
#include <pixie/rmq/hybrid_btree.h>
#include <pixie/rmq/segment_tree.h>
#include <pixie/rmq/sparse_table.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>
#include <pixie/storage/aligned.h>

#ifdef SDSL_SUPPORT
#include <pixie/rmq/sdsl_sct.h>
#endif

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <class T, class Compare>
std::size_t naive_arg_min(std::span<const T> values,
                          std::size_t left,
                          std::size_t right,
                          Compare compare) {
  if (left >= right || right > values.size()) {
    return pixie::rmq::SparseTable<T, Compare>::npos;
  }
  std::size_t best = left;
  for (std::size_t i = left + 1; i < right; ++i) {
    if (compare(values[i], values[best])) {
      best = i;
    }
  }
  return best;
}

template <class Rmq, class T, class Compare>
void check_all_ranges(const Rmq& rmq,
                      std::span<const T> values,
                      Compare compare) {
  ASSERT_EQ(rmq.size(), values.size());
  for (std::size_t left = 0; left < values.size(); ++left) {
    for (std::size_t right = left + 1; right <= values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
      EXPECT_EQ(rmq.range_min(left, right), values[expected])
          << "range=[" << left << "," << right << ")";
    }
  }
}

template <class Rmq, class T, class Compare>
void check_all_arg_min_ranges(const Rmq& rmq,
                              std::span<const T> values,
                              Compare compare) {
  ASSERT_EQ(rmq.size(), values.size());
  for (std::size_t left = 0; left < values.size(); ++left) {
    for (std::size_t right = left + 1; right <= values.size(); ++right) {
      const std::size_t expected = naive_arg_min(values, left, right, compare);
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

bool packed_bit(std::span<const std::uint64_t> words, std::size_t position) {
  return ((words[position >> 6] >> (position & 63)) & 1u) != 0;
}

std::vector<std::uint64_t> make_packed_delta_bits(std::size_t bit_count) {
  std::vector<std::uint64_t> words((bit_count + 63) / 64);
  for (std::size_t position = 0; position < bit_count; ++position) {
    const bool bit = ((position * 11 + position / 3) % 17) < 8;
    if (bit) {
      words[position >> 6] |= std::uint64_t{1} << (position & 63);
    }
  }
  return words;
}

std::vector<std::int64_t> depths_from_delta_bits(
    std::span<const std::uint64_t> words,
    std::size_t depth_count) {
  std::vector<std::int64_t> depths(depth_count);
  for (std::size_t position = 1; position < depth_count; ++position) {
    depths[position] =
        depths[position - 1] + (packed_bit(words, position - 1) ? 1 : -1);
  }
  return depths;
}

std::size_t naive_depth_arg_min(std::span<const std::int64_t> depths,
                                std::size_t left,
                                std::size_t right) {
  if (left >= right || right > depths.size()) {
    return std::numeric_limits<std::size_t>::max();
  }
  std::size_t best = left;
  for (std::size_t position = left + 1; position < right; ++position) {
    if (depths[position] < depths[best]) {
      best = position;
    }
  }
  return best;
}

template <class Rmq>
void expect_valid_bp_encoding(const Rmq& rmq, std::size_t value_count) {
  const std::size_t bit_count = rmq.bp_bit_count();
  const std::span<const std::uint64_t> words = rmq.bp_words();
  EXPECT_EQ(bit_count, 2 * value_count);
  EXPECT_EQ(words.size(), (bit_count + 63) / 64);

  std::size_t ones = 0;
  std::size_t zeros = 0;
  std::int64_t excess = 0;
  for (std::size_t position = 0; position < bit_count; ++position) {
    if (packed_bit(words, position)) {
      ++ones;
      ++excess;
    } else {
      ++zeros;
      --excess;
    }
    EXPECT_GE(excess, 0) << "position=" << position;
  }

  EXPECT_EQ(ones, value_count);
  EXPECT_EQ(zeros, value_count);
  EXPECT_EQ(excess, 0);
}

template <class Rmq>
std::string bp_string(const Rmq& rmq) {
  std::string bits;
  bits.reserve(rmq.bp_bit_count());
  const std::span<const std::uint64_t> words = rmq.bp_words();
  for (std::size_t position = 0; position < rmq.bp_bit_count(); ++position) {
    bits.push_back(packed_bit(words, position) ? '1' : '0');
  }
  return bits;
}

template <class Rmq>
void expect_cartesian_bp_shape(const std::vector<int>& values,
                               const std::string& expected) {
  const Rmq rmq{std::span<const int>(values)};
  EXPECT_EQ(bp_string(rmq), expected);
  check_all_arg_min_ranges(rmq, std::span<const int>(values), std::less<int>());
}

struct SparseTableCase {
  using Rmq = pixie::rmq::SparseTable<int>;
  using MaxRmq = pixie::rmq::SparseTable<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct SegmentTreeCase {
  using Rmq = pixie::rmq::SegmentTree<int>;
  using MaxRmq = pixie::rmq::SegmentTree<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct CartesianRmMCase {
  using Rmq = pixie::rmq::CartesianRmM<int>;
  using MaxRmq = pixie::rmq::CartesianRmM<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct CartesianHybridBTreeCase {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  using MaxRmq = pixie::rmq::CartesianHybridBTree<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct CartesianBTreeCase {
  using Rmq = pixie::rmq::CartesianBTree<int>;
  using MaxRmq = pixie::rmq::CartesianBTree<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

struct HybridBTreeCase {
  using Rmq = pixie::rmq::HybridBTree<int>;
  using MaxRmq = pixie::rmq::HybridBTree<int, std::greater<int>>;

  static Rmq make(std::span<const int> values) { return Rmq(values); }

  static MaxRmq make_max(std::span<const int> values) { return MaxRmq(values); }
};

}  // namespace

template <class Rmq>
concept HasRmqSerialization =
    requires(const Rmq& rmq, pixie::BinaryWriter& writer) {
      rmq.serialize(writer);
    };

template <class Rmq>
concept HasRmqDeserialization = requires(std::span<const std::int64_t> values,
                                         pixie::BinaryReader& reader) {
  { Rmq::deserialize(reader, values) } -> std::same_as<Rmq>;
};

using SerializableRmq = pixie::rmq::CartesianHybridBTree<std::int64_t>;
static_assert(HasRmqSerialization<SerializableRmq>);
static_assert(HasRmqDeserialization<SerializableRmq>);
static_assert(
    !HasRmqSerialization<pixie::rmq::CartesianHybridBTree<std::int32_t>>);
static_assert(!HasRmqSerialization<
              pixie::rmq::CartesianHybridBTree<std::int64_t,
                                               std::greater<std::int64_t>>>);
static_assert(!HasRmqSerialization<
              pixie::rmq::CartesianHybridBTree<std::int64_t,
                                               std::less<std::int64_t>,
                                               std::uint32_t>>);
static_assert(!HasRmqSerialization<
              pixie::rmq::CartesianHybridBTree<std::int64_t,
                                               std::less<std::int64_t>,
                                               std::size_t,
                                               1024>>);
static_assert(!HasRmqSerialization<pixie::rmq::CartesianBTree<std::int64_t>>);

static std::vector<std::byte> serialize_rmq(const SerializableRmq& rmq) {
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  rmq.serialize(writer);
  writer.finish();
  return output.take();
}

static void overwrite_u64(std::vector<std::byte>& bytes,
                          std::size_t offset,
                          std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
  }
}

struct RmqArtifactLayout {
  std::size_t bp_bit_count;
  std::size_t top_block_size;
  std::size_t top_block_count;
  std::size_t top_sparse_levels;
  std::size_t bp_storage_size;
  std::size_t candidate_count;
  std::size_t candidates;
  std::size_t has_rank_index;
  std::size_t depth_selectors;
  std::size_t depth_selector_count;
  std::size_t depth_min_depths;
  std::size_t candidate_size;
  std::size_t bp_storage_bytes;
};

static RmqArtifactLayout locate_rmq_artifact(
    std::span<const std::byte> artifact) {
  pixie::BinaryReader reader(artifact);
  reader.skip(6 * sizeof(std::uint64_t));
  RmqArtifactLayout result;
  result.bp_bit_count = reader.position();
  reader.skip(sizeof(std::uint64_t));
  result.top_block_size = reader.position();
  reader.skip(sizeof(std::uint64_t));
  result.top_block_count = reader.position();
  reader.skip(sizeof(std::uint64_t));
  result.top_sparse_levels = reader.position();
  reader.skip(sizeof(std::uint64_t));
  result.bp_storage_size = reader.position();
  result.bp_storage_bytes = reader.read_size();
  reader.skip(result.bp_storage_bytes);
  result.candidate_count = reader.position();
  result.candidate_size = reader.read_size();
  result.candidates = reader.position();
  reader.skip(result.candidate_size * sizeof(std::uint64_t));
  result.has_rank_index = reader.position();
  const bool has_rank_index = reader.read_u8() != 0;
  if (has_rank_index) {
    reader.skip(7 * sizeof(std::uint64_t) + 2 * sizeof(std::uint32_t));
    for (std::size_t storage = 0; storage < 3; ++storage) {
      reader.skip(reader.read_size());
    }
  }
  reader.skip(sizeof(std::uint64_t));
  result.depth_selector_count = reader.read_size();
  result.depth_selectors = reader.position();
  reader.skip(result.depth_selector_count * 8 * sizeof(std::uint64_t));
  reader.skip(reader.read_size() * sizeof(std::uint64_t));
  const std::size_t min_depth_count = reader.read_size();
  result.depth_min_depths = reader.position();
  reader.skip(min_depth_count * sizeof(std::int64_t));
  return result;
}

static void expect_serialized_rmq_ranges(const SerializableRmq& original,
                                         const SerializableRmq& restored,
                                         std::span<const std::int64_t> values,
                                         bool exhaustive) {
  if (exhaustive) {
    check_all_ranges(restored, values, std::less<std::int64_t>());
    return;
  }

  const std::array<std::pair<std::size_t, std::size_t>, 10> boundaries = {
      std::pair{std::size_t{0}, values.size()},
      std::pair{std::size_t{0}, std::min<std::size_t>(512, values.size())},
      std::pair{std::min<std::size_t>(511, values.size() - 1),
                std::min<std::size_t>(513, values.size())},
      std::pair{std::min<std::size_t>(4095, values.size() - 1),
                std::min<std::size_t>(4097, values.size())},
      std::pair{values.size() / 3, values.size()},
      std::pair{values.size() / 4, 3 * values.size() / 4},
      std::pair{values.size() - 1, values.size()},
      std::pair{std::size_t{1}, values.size()},
      std::pair{values.size() / 2, values.size() / 2 + 1},
      std::pair{std::size_t{0}, std::size_t{1}}};
  for (const auto& [left, right] : boundaries) {
    if (left >= right) {
      continue;
    }
    const std::size_t expected =
        naive_arg_min(values, left, right, std::less<std::int64_t>());
    EXPECT_EQ(original.arg_min(left, right), expected);
    EXPECT_EQ(restored.arg_min(left, right), expected);
    EXPECT_EQ(restored.range_min(left, right), values[expected]);
  }

  std::mt19937_64 rng(20260718);
  std::uniform_int_distribution<std::size_t> position(0, values.size() - 1);
  for (std::size_t query = 0; query < 2000; ++query) {
    std::size_t left = position(rng);
    std::size_t right = position(rng);
    if (left > right) {
      std::swap(left, right);
    }
    ++right;
    const std::size_t expected =
        naive_arg_min(values, left, right, std::less<std::int64_t>());
    EXPECT_EQ(restored.arg_min(left, right), expected);
  }
}

static void expect_rmq_serialization_round_trip(
    std::size_t size,
    pixie::DeserializationValidation validation) {
  std::vector<std::int64_t> values(size);
  for (std::size_t i = 0; i < size; ++i) {
    values[i] = static_cast<std::int64_t>((i * 37 + i / 5 + i / 97) % 31) - 15;
  }
  if (size > 3) {
    values[size / 3] = -100;
    values[size / 3 + 1] = -100;
  }

  const SerializableRmq original(values);
  std::vector<std::byte> artifact = serialize_rmq(original);
  pixie::BinaryReader reader(artifact);
  const SerializableRmq restored =
      SerializableRmq::deserialize(reader, values, validation);

  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(restored.size(), values.size());
  EXPECT_EQ(restored.bp_bit_count(), original.bp_bit_count());
  EXPECT_EQ(restored.bp_words().size(), original.bp_words().size());
  EXPECT_EQ(serialize_rmq(restored), artifact);

  artifact.clear();
  artifact.shrink_to_fit();
  if (values.empty()) {
    EXPECT_TRUE(restored.empty());
    EXPECT_EQ(restored.arg_min(0, 0), SerializableRmq::npos);
    return;
  }
  expect_serialized_rmq_ranges(original, restored, values, size <= 513);
}

TEST(RmqSerializationTest, RoundTripsOwningMetadataAtBoundaries) {
  for (const std::size_t size : std::array<std::size_t, 11>{
           0, 1, 255, 256, 257, 511, 512, 513, 4095, 4096, 4097}) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      SCOPED_TRACE(::testing::Message() << "size=" << size << " validation="
                                        << static_cast<int>(validation));
      expect_rmq_serialization_round_trip(size, validation);
    }
  }
}

TEST(RmqSerializationTest, AdvancesAcrossConcatenatedArtifacts) {
  std::vector<std::int64_t> values(33);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<std::int64_t>((i * 11) % 17);
  }
  const SerializableRmq original(values);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  original.serialize(writer);
  writer.finish();
  const auto artifacts = output.take();
  pixie::BinaryReader reader(artifacts);

  const auto first = SerializableRmq::deserialize(reader, values);
  EXPECT_FALSE(reader.empty());
  const auto second = SerializableRmq::deserialize(reader, values);
  EXPECT_TRUE(reader.empty());
  EXPECT_EQ(first.arg_min(0, values.size()), second.arg_min(0, values.size()));
}

TEST(RmqSerializationTest, SerializesDirectlyToMappedFile) {
  std::vector<std::int64_t> values(65);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<std::int64_t>((i * 19) % 23) - 11;
  }
  const SerializableRmq original(values);
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_rmq_serialization_test.bin";
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
    const SerializableRmq restored =
        SerializableRmq::deserialize(reader, values, validation);
    EXPECT_TRUE(reader.empty());
    check_all_ranges(restored, std::span<const std::int64_t>(values),
                     std::less<std::int64_t>());
  }
  std::filesystem::remove(path);
}

TEST(RmqSerializationTest, RejectsCorruptionWithoutAdvancingInput) {
  std::vector<std::int64_t> values(40);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<std::int64_t>((i * 13) % 23);
  }
  const SerializableRmq original(values);
  const auto valid = serialize_rmq(original);

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW(
          (void)SerializableRmq::deserialize(reader, values, validation),
          std::exception);
      EXPECT_EQ(reader.position(), 0u);
    }
  };
  const auto overwrite = [](std::vector<std::byte>& artifact,
                            std::size_t offset, auto value) {
    using Value = decltype(value);
    using Unsigned = std::make_unsigned_t<Value>;
    const Unsigned encoded = static_cast<Unsigned>(value);
    for (std::size_t byte = 0; byte < sizeof(Value); ++byte) {
      artifact[offset + byte] =
          static_cast<std::byte>((encoded >> (byte * 8)) & Unsigned{0xff});
    }
  };

  auto bad_magic = valid;
  bad_magic[0] ^= std::byte{1};
  expect_rejected(std::move(bad_magic));

  auto bad_version = valid;
  overwrite(bad_version, 8, std::uint32_t{5});
  expect_rejected(std::move(bad_version));

  auto bad_leaf_size = valid;
  overwrite(bad_leaf_size, 32, std::uint64_t{1024});
  expect_rejected(std::move(bad_leaf_size));

  auto bad_bp_count = valid;
  overwrite(bad_bp_count, 48, std::uint64_t{0});
  expect_rejected(std::move(bad_bp_count));

  auto impossible_storage = valid;
  overwrite(impossible_storage, 80, std::numeric_limits<std::size_t>::max());
  expect_rejected(std::move(impossible_storage));

  std::span<const std::byte> truncated =
      std::span<const std::byte>(valid).first(valid.size() - 1);
  pixie::BinaryReader truncated_reader(truncated);
  EXPECT_THROW((void)SerializableRmq::deserialize(truncated_reader, values),
               std::invalid_argument);
  EXPECT_EQ(truncated_reader.position(), 0u);

  std::span<const std::int64_t> short_values(values.data(), values.size() - 1);
  pixie::BinaryReader reader(valid);
  EXPECT_THROW((void)SerializableRmq::deserialize(reader, short_values),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);

  std::vector<std::int64_t> different_values(values.size(), -1);
  pixie::BinaryReader different_reader(valid);
  EXPECT_NO_THROW(
      (void)SerializableRmq::deserialize(different_reader, different_values));
  EXPECT_TRUE(different_reader.empty());

  pixie::BinaryReader different_full_reader(valid);
  EXPECT_THROW((void)SerializableRmq::deserialize(
                   different_full_reader, different_values,
                   pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(different_full_reader.position(), 0u);
}

TEST(RmqSerializationTest, RejectsMalformedBpEncodingAndPadding) {
  constexpr std::size_t kBpDataOffset = 88;
  const auto set_bit = [](std::vector<std::byte>& artifact,
                          std::size_t position) {
    artifact[kBpDataOffset + position / 8] |=
        static_cast<std::byte>(std::uint8_t{1} << (position % 8));
  };
  const auto clear_bit = [](std::vector<std::byte>& artifact,
                            std::size_t position) {
    artifact[kBpDataOffset + position / 8] &=
        static_cast<std::byte>(~(std::uint8_t{1} << (position % 8)));
  };

  for (const std::size_t size :
       {std::size_t{63}, std::size_t{64}, std::size_t{65}}) {
    SCOPED_TRACE(::testing::Message() << "size=" << size);
    std::vector<std::int64_t> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<std::int64_t>((i * 13 + i / 3) % 17);
    }
    const SerializableRmq original(values);
    const std::vector<std::byte> valid = serialize_rmq(original);
    const auto expect_rejected = [&](std::vector<std::byte> artifact) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW((void)SerializableRmq::deserialize(
                       reader, values, pixie::DeserializationValidation::kFull),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    };

    auto negative_first_prefix = valid;
    clear_bit(negative_first_prefix, 0);
    expect_rejected(std::move(negative_first_prefix));

    auto nonzero_final_excess = valid;
    set_bit(nonzero_final_excess, 2 * size - 1);
    expect_rejected(std::move(nonzero_final_excess));

    auto nonzero_padding = valid;
    set_bit(nonzero_padding, 2 * size);
    expect_rejected(std::move(nonzero_padding));
  }

  constexpr std::size_t kLateViolationSize = 129;
  std::vector<std::int64_t> values(kLateViolationSize);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<std::int64_t>((i * 19 + i / 7) % 23);
  }
  const SerializableRmq original(values);
  std::vector<std::byte> negative_later = serialize_rmq(original);
  for (std::size_t position = 128; position < 2 * values.size(); ++position) {
    clear_bit(negative_later, position);
  }
  pixie::BinaryReader reader(negative_later);
  EXPECT_THROW((void)SerializableRmq::deserialize(
                   reader, values, pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(reader.position(), 0u);
}

TEST(RmqSerializationTest, RejectsMalformedIndexMetadataTransactionally) {
  std::vector<std::int64_t> values(4097);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<std::int64_t>((i * 17 + i / 11) % 29);
  }
  const SerializableRmq original(values);
  const std::vector<std::byte> valid = serialize_rmq(original);
  const RmqArtifactLayout layout = locate_rmq_artifact(valid);
  ASSERT_EQ(layout.candidate_size, 4u);
  ASSERT_EQ(layout.depth_selector_count, 1u);

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW(
          (void)SerializableRmq::deserialize(reader, values, validation),
          std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  };

  auto misaligned_storage = valid;
  overwrite_u64(misaligned_storage, layout.bp_storage_size,
                layout.bp_storage_bytes - 1);
  expect_rejected(std::move(misaligned_storage));

  auto insufficient_storage = valid;
  overwrite_u64(insufficient_storage, layout.bp_bit_count,
                layout.bp_storage_bytes * 8 + 64);
  expect_rejected(std::move(insufficient_storage));

  auto wrong_bp_bit_count = valid;
  overwrite_u64(wrong_bp_bit_count, layout.bp_bit_count, 2 * values.size() - 1);
  expect_rejected(std::move(wrong_bp_bit_count));

  auto wrong_top_block_size = valid;
  overwrite_u64(wrong_top_block_size, layout.top_block_size, 1);
  expect_rejected(std::move(wrong_top_block_size));

  auto wrong_top_block_count = valid;
  overwrite_u64(wrong_top_block_count, layout.top_block_count, 0);
  expect_rejected(std::move(wrong_top_block_count));

  auto wrong_top_sparse_levels = valid;
  overwrite_u64(wrong_top_sparse_levels, layout.top_sparse_levels, 0);
  expect_rejected(std::move(wrong_top_sparse_levels));

  auto invalid_marker = valid;
  invalid_marker[layout.has_rank_index] = std::byte{2};
  expect_rejected(std::move(invalid_marker));

  auto invalid_candidate = valid;
  overwrite_u64(invalid_candidate, layout.candidates, values.size());
  expect_rejected(std::move(invalid_candidate));

  auto invalid_sparse_padding = valid;
  overwrite_u64(invalid_sparse_padding,
                layout.candidates + 3 * sizeof(std::uint64_t), 0);
  expect_rejected(std::move(invalid_sparse_padding));

  auto wrong_exact_block_minimum = valid;
  const std::size_t original_candidate =
      pixie::BinaryReader(std::span<const std::byte>(valid).subspan(
                              layout.candidates, sizeof(std::uint64_t)))
          .read_size();
  const std::size_t alternate_candidate =
      original_candidate == 0 ? 1 : original_candidate - 1;
  overwrite_u64(wrong_exact_block_minimum, layout.candidates,
                alternate_candidate);
  pixie::BinaryReader exact_quick_reader(wrong_exact_block_minimum);
  EXPECT_NO_THROW((void)SerializableRmq::deserialize(
      exact_quick_reader, values, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(exact_quick_reader.empty());
  pixie::BinaryReader exact_full_reader(wrong_exact_block_minimum);
  EXPECT_THROW(
      (void)SerializableRmq::deserialize(
          exact_full_reader, values, pixie::DeserializationValidation::kFull),
      std::invalid_argument);
  EXPECT_EQ(exact_full_reader.position(), 0u);

  auto wrong_depth_selector = valid;
  wrong_depth_selector[layout.depth_selectors] ^= std::byte{1};
  pixie::BinaryReader selector_quick_reader(wrong_depth_selector);
  EXPECT_NO_THROW((void)SerializableRmq::deserialize(
      selector_quick_reader, values, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(selector_quick_reader.empty());
  pixie::BinaryReader selector_full_reader(wrong_depth_selector);
  EXPECT_THROW((void)SerializableRmq::deserialize(
                   selector_full_reader, values,
                   pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(selector_full_reader.position(), 0u);

  auto wrong_depth_minimum = valid;
  overwrite_u64(wrong_depth_minimum, layout.depth_min_depths,
                std::numeric_limits<std::uint64_t>::max());
  pixie::BinaryReader depth_full_reader(wrong_depth_minimum);
  EXPECT_THROW(
      (void)SerializableRmq::deserialize(
          depth_full_reader, values, pixie::DeserializationValidation::kFull),
      std::invalid_argument);
  EXPECT_EQ(depth_full_reader.position(), 0u);

  const std::vector<std::int64_t> empty_values;
  const SerializableRmq empty(empty_values);
  std::vector<std::byte> invalid_empty = serialize_rmq(empty);
  const RmqArtifactLayout empty_layout = locate_rmq_artifact(invalid_empty);
  overwrite_u64(invalid_empty, empty_layout.top_block_size, 1);
  pixie::BinaryReader empty_reader(invalid_empty);
  EXPECT_THROW((void)SerializableRmq::deserialize(empty_reader, empty_values),
               std::invalid_argument);
  EXPECT_EQ(empty_reader.position(), 0u);
}

template <class Case>
class ValueRmqSpecificationTest : public ::testing::Test {};

using ValueRmqCases = ::testing::Types<SparseTableCase,
                                       SegmentTreeCase,
                                       CartesianRmMCase,
                                       CartesianHybridBTreeCase,
                                       CartesianBTreeCase,
                                       HybridBTreeCase>;
TYPED_TEST_SUITE(ValueRmqSpecificationTest, ValueRmqCases);

TYPED_TEST(ValueRmqSpecificationTest, ExhaustiveSmallArray) {
  const std::vector<int> values = {4, 1, 3, 1, 5, 0, 0, 2};
  const typename TypeParam::Rmq rmq =
      TypeParam::make(std::span<const int>(values));
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

TEST(RmqCartesianBuildShape, SuccinctStackPreservesBpEncoding) {
  using Hybrid = pixie::rmq::CartesianHybridBTree<int>;
  using RmM = pixie::rmq::CartesianRmM<int>;

  const std::vector<std::pair<std::vector<int>, std::string>> cases = {
      {{1}, "10"},
      {{1, 2}, "1010"},
      {{2, 1}, "1100"},
      {{1, 1}, "1010"},
      {{3, 1, 2}, "110010"},
      {{1, 3, 2}, "101100"},
      {{3, 2, 2, 2, 1, 1, 2, 1, 3}, "111001010010110010"},
      {{4, 1, 3, 1, 5, 0, 0, 2}, "1110011001001010"},
  };

  for (const auto& [values, expected] : cases) {
    SCOPED_TRACE(expected);
    expect_cartesian_bp_shape<Hybrid>(values, expected);
    expect_cartesian_bp_shape<RmM>(values, expected);
  }
}

TEST(RmqSuccinctIncreasingStack, PopsAcrossEmptyBlocks) {
  pixie::rmq::utils::SuccinctIncreasingStack stack(256);
  EXPECT_TRUE(stack.empty());

  stack.push(0);
  EXPECT_EQ(stack.top(), 0u);
  stack.push(62);
  stack.push(63);
  EXPECT_EQ(stack.top(), 63u);
  stack.pop();
  EXPECT_EQ(stack.top(), 62u);
  stack.pop();
  EXPECT_EQ(stack.top(), 0u);

  stack.push(130);
  EXPECT_EQ(stack.top(), 130u);
  stack.pop();
  EXPECT_EQ(stack.top(), 0u);
  stack.push(131);
  stack.pop();
  EXPECT_EQ(stack.top(), 0u);

  stack.push(190);
  stack.push(191);
  EXPECT_EQ(stack.top(), 191u);
  stack.pop();
  EXPECT_EQ(stack.top(), 190u);
  stack.pop();
  EXPECT_EQ(stack.top(), 0u);
  stack.pop();
  EXPECT_TRUE(stack.empty());
}

TYPED_TEST(ValueRmqSpecificationTest, FirstMinimumTieBreaking) {
  const std::vector<int> values = {7, 2, 2, 3, 2};
  const typename TypeParam::Rmq rmq =
      TypeParam::make(std::span<const int>(values));

  EXPECT_EQ(rmq.arg_min(0, 5), 1u);
  EXPECT_EQ(rmq.arg_min(2, 5), 2u);
}

TYPED_TEST(ValueRmqSpecificationTest, InvalidAndEmptyRanges) {
  using Rmq = typename TypeParam::Rmq;

  const std::vector<int> values = {3, 1, 2};
  const Rmq rmq = TypeParam::make(std::span<const int>(values));

  EXPECT_EQ(rmq.arg_min(2, 2), Rmq::npos);
  EXPECT_EQ(rmq.arg_min(2, 1), Rmq::npos);
  EXPECT_EQ(rmq.arg_min(0, values.size() + 1), Rmq::npos);
  EXPECT_EQ(rmq.range_min(2, 2), 0);
  EXPECT_EQ(rmq.arg_min(0, values.size()), 1u);

  const std::vector<int> empty_values;
  const Rmq default_rmq;
  const Rmq empty_rmq = TypeParam::make(std::span<const int>(empty_values));
  EXPECT_TRUE(default_rmq.empty());
  EXPECT_TRUE(empty_rmq.empty());
  EXPECT_EQ(default_rmq.arg_min(0, 0), Rmq::npos);
  EXPECT_EQ(empty_rmq.arg_min(0, 0), Rmq::npos);
  EXPECT_EQ(default_rmq.range_min(0, 0), 0);
  EXPECT_EQ(empty_rmq.range_min(0, 0), 0);
}

TYPED_TEST(ValueRmqSpecificationTest, MemoryUsageCountsOwnedIndexStorage) {
  using Rmq = typename TypeParam::Rmq;
  using MaxRmq = typename TypeParam::MaxRmq;

  const Rmq default_rmq;
  EXPECT_GE(default_rmq.memory_usage_bytes(), sizeof(Rmq));

  std::vector<int> values(1537);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 19 + i / 5) % 127);
  }

  const Rmq rmq = TypeParam::make(std::span<const int>(values));
  const MaxRmq max_rmq = TypeParam::make_max(std::span<const int>(values));
  EXPECT_GE(rmq.memory_usage_bytes(), sizeof(Rmq));
  EXPECT_GE(max_rmq.memory_usage_bytes(), sizeof(MaxRmq));
  EXPECT_GT(rmq.memory_usage_bytes(), default_rmq.memory_usage_bytes());
}

TYPED_TEST(ValueRmqSpecificationTest, ComparatorCanSelectMaximum) {
  using MaxRmq = typename TypeParam::MaxRmq;

  const std::vector<int> values = {1, 8, 3, 8, 4};
  const MaxRmq rmq = TypeParam::make_max(std::span<const int>(values));

  check_all_ranges(rmq, std::span<const int>(values), std::greater<int>());
  EXPECT_EQ(rmq.arg_min(0, 5), 1u);
  EXPECT_EQ(rmq.arg_min(2, 2), MaxRmq::npos);
  EXPECT_EQ(rmq.arg_min(4, 3), MaxRmq::npos);
  EXPECT_EQ(rmq.arg_min(0, values.size() + 1), MaxRmq::npos);
  EXPECT_EQ(rmq.range_min(2, 2), 0);
  EXPECT_EQ(rmq.range_min(4, 3), 0);
  EXPECT_EQ(rmq.range_min(0, values.size() + 1), 0);

  const std::vector<int> empty_values;
  const MaxRmq default_rmq;
  const MaxRmq empty_rmq =
      TypeParam::make_max(std::span<const int>(empty_values));
  EXPECT_TRUE(default_rmq.empty());
  EXPECT_TRUE(empty_rmq.empty());
  EXPECT_EQ(default_rmq.arg_min(0, 0), MaxRmq::npos);
  EXPECT_EQ(empty_rmq.arg_min(0, 0), MaxRmq::npos);
  EXPECT_EQ(default_rmq.range_min(0, 0), 0);
  EXPECT_EQ(empty_rmq.range_min(0, 0), 0);
}

TYPED_TEST(ValueRmqSpecificationTest, MonotoneArrays) {
  const std::vector<int> increasing = {1, 2, 3, 4, 5, 6};
  const std::vector<int> decreasing = {6, 5, 4, 3, 2, 1};
  const typename TypeParam::Rmq increasing_rmq =
      TypeParam::make(std::span<const int>(increasing));
  const typename TypeParam::Rmq decreasing_rmq =
      TypeParam::make(std::span<const int>(decreasing));

  check_all_ranges(increasing_rmq, std::span<const int>(increasing),
                   std::less<int>());
  check_all_ranges(decreasing_rmq, std::span<const int>(decreasing),
                   std::less<int>());
}

TYPED_TEST(ValueRmqSpecificationTest, DifferentialRandom) {
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> value_dist(-50, 50);
  for (std::size_t size = 1; size <= 257; size += 17) {
    std::vector<int> values(size);
    std::generate(values.begin(), values.end(),
                  [&] { return value_dist(rng); });

    const typename TypeParam::Rmq rmq =
        TypeParam::make(std::span<const int>(values));
    check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
  }
}

TEST(RmqSparseTable, OverlappingBlockCandidateDirectionsAndTies) {
  {
    const std::vector<int> values = {5, 4, 3, 2, 1, 0, 9};
    const pixie::rmq::SparseTable<int> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 5u);
    EXPECT_EQ(rmq.range_min(0, 6), 0);
  }

  {
    const std::vector<int> values = {0, 5, 6, 7, 8, 9, 1};
    const pixie::rmq::SparseTable<int> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 0u);
    EXPECT_EQ(rmq.range_min(0, 6), 0);
  }

  {
    const std::vector<int> values = {1, 0, 9, 9, 0, 2};
    const pixie::rmq::SparseTable<int> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 1u);
    EXPECT_EQ(rmq.range_min(0, 6), 0);
  }

  {
    const std::vector<int> values = {1, 2, 3, 4, 9, 10, 0};
    const pixie::rmq::SparseTable<int, std::greater<int>> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 5u);
    EXPECT_EQ(rmq.range_min(0, 6), 10);
  }

  {
    const std::vector<int> values = {10, 9, 8, 7, 6, 5, 20};
    const pixie::rmq::SparseTable<int, std::greater<int>> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 0u);
    EXPECT_EQ(rmq.range_min(0, 6), 10);
  }

  {
    const std::vector<int> values = {10, 1, 2, 3, 10, 4};
    const pixie::rmq::SparseTable<int, std::greater<int>> rmq(values);
    EXPECT_EQ(rmq.arg_min(0, 6), 0u);
    EXPECT_EQ(rmq.range_min(0, 6), 10);
  }
}

TEST(RmqSegmentTree, NonPowerOfTwoTailRanges) {
  const std::vector<int> values = {6, 5, 4, 3, 2};
  const pixie::rmq::SegmentTree<int> rmq(values);
  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
}

#ifdef SDSL_SUPPORT
TEST(RmqSdslSct,
     MatchesPixieValueRmqSpecificationForSignedValuesAndDuplicates) {
  const std::vector<int> values = {4, -3, 7, -3, 0, -8, -8, 2, 2};
  const pixie::rmq::SdslSct<int> rmq{std::span<const int>(values)};

  check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
  EXPECT_EQ(rmq.arg_min(0, values.size()), 5u);
  EXPECT_EQ(rmq.arg_min(5, 7), 5u);
  EXPECT_EQ(rmq.arg_min(6, 7), 6u);
  EXPECT_EQ(rmq.arg_min(3, 3), pixie::rmq::SdslSct<int>::npos);
  EXPECT_EQ(rmq.arg_min(0, values.size() + 1), pixie::rmq::SdslSct<int>::npos);
}

TEST(RmqSdslSct, DifferentialRandom) {
  std::mt19937_64 rng(123);
  std::uniform_int_distribution<int> value_dist(-20, 20);
  for (std::size_t size = 1; size <= 129; size += 16) {
    std::vector<int> values(size);
    std::generate(values.begin(), values.end(),
                  [&] { return value_dist(rng); });

    const pixie::rmq::SdslSct<int> rmq{std::span<const int>(values)};
    check_all_ranges(rmq, std::span<const int>(values), std::less<int>());
  }

  const std::vector<int> empty_values;
  const pixie::rmq::SdslSct<int> empty{std::span<const int>(empty_values)};
  EXPECT_TRUE(empty.empty());
  EXPECT_EQ(empty.arg_min(0, 0), pixie::rmq::SdslSct<int>::npos);
}
#endif

TEST(RmqHybridBTree, BoundaryAndFallbackRanges) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  static_assert(kLeaf == 496);

  std::vector<int> values(2 * kLeaf + 17, 2000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>((i * 7 + i / 3) % 31);
  }

  values[32] = 120;
  values[96] = 70;
  values[160] = 55;
  values[240] = -100;
  values[241] = -100;
  values[320] = 30;
  values[440] = 5;
  values[kLeaf - 1] = 5;

  values[kLeaf + 11] = -9;
  values[kLeaf + 173] = -17;
  values[2 * kLeaf + 5] = -21;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kLeaf},         {0, 97},
      {0, 240},           {96, 241},
      {250, kLeaf},       {441, kLeaf},
      {300, 450},         {kLeaf - 7, kLeaf + 12},
      {kLeaf, 2 * kLeaf}, {2 * kLeaf - 5, values.size()},
      {0, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(rmq.range_min(left, right), values[expected])
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqHybridBTree, MaskLeafSelectorPrefixSuffixAndInteriorFallback) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;

  {
    std::vector<int> values(kLeaf, 1000);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] += static_cast<int>(i % 17);
    }
    values[0] = -10;
    values[kLeaf - 1] = -100;

    const Rmq rmq{std::span<const int>(values)};
    const std::size_t right = kLeaf / 2;
    EXPECT_EQ(rmq.arg_min(0, right), naive_arg_min(std::span<const int>(values),
                                                   0, right, std::less<int>()));
  }

  {
    std::vector<int> values(kLeaf, 1000);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] += static_cast<int>((i * 5) % 23);
    }
    values[0] = -100;
    values[kLeaf - 1] = -10;

    const Rmq rmq{std::span<const int>(values)};
    const std::size_t left = kLeaf / 2;
    EXPECT_EQ(rmq.arg_min(left, kLeaf),
              naive_arg_min(std::span<const int>(values), left, kLeaf,
                            std::less<int>()));
  }

  {
    std::vector<int> values(kLeaf, 1000);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] += static_cast<int>((i * 7) % 31);
    }
    values[0] = -100;
    values[kLeaf / 2] = -50;
    values[kLeaf - 1] = -80;

    const Rmq rmq{std::span<const int>(values)};
    const std::size_t left = kLeaf / 2 - 20;
    const std::size_t right = kLeaf / 2 + 21;
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::less<int>()));
  }
}

TEST(RmqHybridBTree, LeafSelectorEnumVariants) {
  using MaskRmq = pixie::rmq::HybridBTree<
      int, std::less<int>, std::size_t, 248, 256,
      pixie::rmq::HybridBTreeLeafSelector::PrefixSuffix>;
  using BpRmq =
      pixie::rmq::HybridBTree<int, std::less<int>, std::size_t, 252, 256,
                              pixie::rmq::HybridBTreeLeafSelector::BP>;

  EXPECT_EQ(MaskRmq::top_sparse_block_size_for(0),
            MaskRmq::kMinTopSparseBlockSize);
  EXPECT_EQ(MaskRmq::top_sparse_block_count_for(0), 0u);
  EXPECT_EQ(BpRmq::top_sparse_block_size_for(0), BpRmq::kMinTopSparseBlockSize);
  EXPECT_EQ(BpRmq::top_sparse_block_count_for(0), 0u);

  std::vector<int> values(4099, 10000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = static_cast<int>((i * 37 + i / 7) % 257);
  }
  values[13] = -50;
  values[251] = -70;
  values[252] = -70;
  values[747] = -90;
  values[2048] = -120;
  values[4098] = -110;

  const MaskRmq mask_rmq{std::span<const int>(values)};
  const BpRmq bp_rmq{std::span<const int>(values)};
  const std::vector<int> empty_values;
  const MaskRmq empty_mask_rmq;
  const BpRmq empty_bp_rmq;
  const MaskRmq empty_mask_span_rmq{std::span<const int>(empty_values)};
  const BpRmq empty_bp_span_rmq{std::span<const int>(empty_values)};
  EXPECT_EQ(empty_mask_rmq.arg_min(0, 0), MaskRmq::npos);
  EXPECT_EQ(empty_bp_rmq.arg_min(0, 0), BpRmq::npos);
  EXPECT_EQ(empty_mask_span_rmq.arg_min(0, 0), MaskRmq::npos);
  EXPECT_EQ(empty_bp_span_rmq.arg_min(0, 0), BpRmq::npos);
  EXPECT_EQ(empty_mask_rmq.range_min(0, 0), 0);
  EXPECT_EQ(empty_bp_rmq.range_min(0, 0), 0);
  EXPECT_EQ(empty_mask_span_rmq.range_min(0, 0), 0);
  EXPECT_EQ(empty_bp_span_rmq.range_min(0, 0), 0);

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, 1},       {0, 13},     {0, 252},     {1, 251},
      {20, 40},     {100, 200},  {200, 248},   {248, 253},
      {251, 753},   {700, 2100}, {2048, 2049}, {0, values.size()},
      {3000, 4099},
  };

  for (const auto& [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(mask_rmq.arg_min(left, right), expected)
        << "mask range=[" << left << "," << right << ")";
    EXPECT_EQ(bp_rmq.arg_min(left, right), expected)
        << "bp range=[" << left << "," << right << ")";
    EXPECT_EQ(mask_rmq.range_min(left, right), values[expected])
        << "mask range=[" << left << "," << right << ")";
    EXPECT_EQ(bp_rmq.range_min(left, right), values[expected])
        << "bp range=[" << left << "," << right << ")";
  }

  const auto check_small_tree_paths = []<class Rmq>() {
    constexpr std::size_t kLeaf = Rmq::kLeafSize;
    std::vector<int> one_leaf_values(kLeaf / 2 + 3, 1000);
    for (std::size_t i = 0; i < one_leaf_values.size(); ++i) {
      one_leaf_values[i] += static_cast<int>((i * 17 + i / 5) % 71);
    }
    one_leaf_values[one_leaf_values.size() / 2] = -20;
    const Rmq one_leaf_rmq{std::span<const int>(one_leaf_values)};
    EXPECT_EQ(one_leaf_rmq.arg_min(0, one_leaf_values.size()),
              naive_arg_min(std::span<const int>(one_leaf_values), 0,
                            one_leaf_values.size(), std::less<int>()))
        << "one leaf path leaf=" << kLeaf;

    std::vector<int> small_values(3 * kLeaf + 5, 10000);
    for (std::size_t i = 0; i < small_values.size(); ++i) {
      small_values[i] += static_cast<int>((i * 41 + i / 9) % 263);
    }
    small_values[7] = -100;
    small_values[kLeaf + 3] = -500;
    small_values[2 * kLeaf + 11] = -300;

    const Rmq rmq{std::span<const int>(small_values)};
    const std::vector<std::pair<std::size_t, std::size_t>> small_ranges = {
        {0, small_values.size()},     {0, kLeaf},
        {kLeaf, 2 * kLeaf},           {2 * kLeaf, 3 * kLeaf},
        {1, small_values.size() - 1}, {kLeaf + 1, kLeaf + 2},
    };
    for (const auto& [left, right] : small_ranges) {
      const std::size_t expected = naive_arg_min(
          std::span<const int>(small_values), left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "small tree range=[" << left << "," << right << ") leaf=" << kLeaf;
    }
  };

  check_small_tree_paths.template operator()<MaskRmq>();
  check_small_tree_paths.template operator()<BpRmq>();
}

TEST(RmqHybridBTree, BorderCorrectionForLeafSelectorVariants) {
  using DefaultRmq = pixie::rmq::HybridBTree<int>;
  using SmallMaskRmq = pixie::rmq::HybridBTree<
      int, std::less<int>, std::size_t, 248, 256,
      pixie::rmq::HybridBTreeLeafSelector::PrefixSuffix>;
  using BpRmq =
      pixie::rmq::HybridBTree<int, std::less<int>, std::size_t, 252, 256,
                              pixie::rmq::HybridBTreeLeafSelector::BP>;

  const auto check = []<class Rmq>() {
    constexpr std::size_t kLeaf = Rmq::kLeafSize;

    {
      std::vector<int> values(3 * kLeaf, 100000);
      for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] += static_cast<int>((i * 13 + i / 5) % 97);
      }
      values[0] = -100000;
      values[17] = -100;
      values[kLeaf + 23] = -1000;
      values[2 * kLeaf + 31] = -500;

      const Rmq rmq{std::span<const int>(values)};
      const std::size_t left = 1;
      const std::size_t right = values.size() - 1;
      EXPECT_EQ(rmq.arg_min(left, right),
                naive_arg_min(std::span<const int>(values), left, right,
                              std::less<int>()))
          << "left-border correction leaf=" << kLeaf;
    }

    {
      std::vector<int> values(3 * kLeaf, 100000);
      for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] += static_cast<int>((i * 11 + i / 7) % 89);
      }
      values[11] = -500;
      values[kLeaf + 29] = -1000;
      values[2 * kLeaf + 41] = -100;
      values[values.size() - 1] = -100000;

      const Rmq rmq{std::span<const int>(values)};
      const std::size_t left = 1;
      const std::size_t right = values.size() - 1;
      EXPECT_EQ(rmq.arg_min(left, right),
                naive_arg_min(std::span<const int>(values), left, right,
                              std::less<int>()))
          << "right-border correction leaf=" << kLeaf;
    }
  };

  check.template operator()<DefaultRmq>();
  check.template operator()<SmallMaskRmq>();
  check.template operator()<BpRmq>();
}

TEST(RmqHybridBTree, MiddleFanoutBoundaryRanges) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;
  constexpr std::size_t kMiddleBoundary = kLeaf * Rmq::kMiddleFanout;

  std::vector<int> values(kMiddleBoundary + 2 * kLeaf + 17, 10000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>((i * 17 + i / 11) % 257);
  }
  values[13] = -100;
  values[kMiddleBoundary - 29] = -90;
  values[kMiddleBoundary + kLeaf + 7] = -110;

  const Rmq rmq{std::span<const int>(values)};
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, values.size()},
      {kMiddleBoundary - 2 * kLeaf, kMiddleBoundary + 2 * kLeaf},
      {kMiddleBoundary - 31, kMiddleBoundary + kLeaf + 31},
      {kMiddleBoundary + 1, values.size()},
      {kLeaf * (Rmq::kMiddleFanout - 1) + 9, kMiddleBoundary + kLeaf + 13},
      {values.size() - 3 * kLeaf, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqHybridBTree, TopSparseOverlayBoundaryRanges) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kTopBlock = Rmq::kMinTopSparseBlockSize;

  std::vector<int> values(3 * kTopBlock + 77, 100000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>((i * 23 + i / 7) % 311);
  }
  values[0] = -1000;
  values[kTopBlock + 17] = -900;
  values[kTopBlock + 18] = -900;
  values[2 * kTopBlock + 41] = -950;
  values[3 * kTopBlock + 40] = -800;

  const Rmq rmq{std::span<const int>(values)};
  EXPECT_EQ(rmq.top_sparse_block_size(), kTopBlock);
  EXPECT_EQ(rmq.top_sparse_block_count(), 4u);

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {1, 2 * kTopBlock + 100},
      {kTopBlock, kTopBlock + 100},
      {kTopBlock + 18, values.size()},
      {kTopBlock - 5, 2 * kTopBlock + 5},
      {0, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqHybridBTree, TopSparseOverlayCapsBlockCount) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kMinBlock = Rmq::kMinTopSparseBlockSize;
  constexpr std::size_t kMaxBlocks = Rmq::kMaxTopSparseBlocks;
  constexpr std::size_t kLargestFixedBlockInput = kMinBlock * kMaxBlocks;

  EXPECT_EQ(Rmq::top_sparse_block_size_for(0), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(0), 0u);
  EXPECT_EQ(Rmq::top_sparse_block_size_for(1), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(1), 1u);
  EXPECT_EQ(Rmq::top_sparse_block_size_for(kLargestFixedBlockInput), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(kLargestFixedBlockInput),
            kMaxBlocks);
  EXPECT_GT(Rmq::top_sparse_block_size_for(kLargestFixedBlockInput + 1),
            kMinBlock);
  EXPECT_LE(Rmq::top_sparse_block_count_for(kLargestFixedBlockInput + 1),
            kMaxBlocks);

  const std::size_t huge = std::numeric_limits<std::size_t>::max() / 4;
  EXPECT_LE(Rmq::top_sparse_block_count_for(huge), kMaxBlocks);

  std::vector<int> values(3 * kMinBlock + 7, 0);
  const Rmq rmq(values);
  EXPECT_EQ(rmq.top_sparse_block_size(), kMinBlock);
  EXPECT_EQ(rmq.top_sparse_block_count(), 4u);
}

TEST(RmqHybridBTree, TopSparseOverlayComparatorMaximum) {
  using Rmq = pixie::rmq::HybridBTree<int, std::greater<int>>;
  constexpr std::size_t kTopBlock = Rmq::kMinTopSparseBlockSize;
  EXPECT_EQ(Rmq::top_sparse_block_size_for(0), kTopBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(0), 0u);

  std::vector<int> values(2 * kTopBlock + 33, -1000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = -1000 + static_cast<int>((i * 13 + i / 5) % 71);
  }
  values[0] = 1000;
  values[kTopBlock + 17] = 900;
  values[kTopBlock + 18] = 900;

  const Rmq rmq(values);
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {1, values.size()},
      {kTopBlock, kTopBlock + 100},
      {kTopBlock + 18, values.size()},
      {0, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    const std::size_t expected = naive_arg_min(
        std::span<const int>(values), left, right, std::greater<int>());
    EXPECT_EQ(rmq.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqHybridBTree, InternalRootFullRangeShortcut) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  constexpr std::size_t kLeaf = Rmq::kLeafSize;

  std::vector<int> values(7 * kLeaf + 19, 10000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] += static_cast<int>((i * 29 + i / 13) % 257);
  }
  values[17] = -100;
  values[3 * kLeaf + 11] = -700;
  values[values.size() - 3] = -200;

  const Rmq rmq{std::span<const int>(values)};
  EXPECT_EQ(rmq.arg_min(0, values.size()), 3 * kLeaf + 11);
  EXPECT_EQ(rmq.range_min(0, values.size()), -700);
  EXPECT_EQ(rmq.arg_min(0, kLeaf), naive_arg_min(std::span<const int>(values),
                                                 0, kLeaf, std::less<int>()));
  EXPECT_EQ(rmq.arg_min(kLeaf - 3, 2 * kLeaf + 5),
            naive_arg_min(std::span<const int>(values), kLeaf - 3,
                          2 * kLeaf + 5, std::less<int>()));
}

TEST(RmqHybridBTree, DuplicateHeavyRandomTo8193) {
  using Rmq = pixie::rmq::HybridBTree<int>;
  std::mt19937_64 rng(901496);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 495, 496, 497, 1024, 4096, 8191, 8192, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 11) < 6) {
        values[i] = 0;
      }
    }

    const Rmq rmq{std::span<const int>(values)};
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqCartesianHybridBTree, BoundarySizesAndBpEncoding) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  const std::vector<std::size_t> sizes = {
      1u, 255u, 256u, 257u, 511u, 512u, 513u, 8191u, 8192u, 8193u,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 37 + i / 5) % 53);
      if (i % 17 == 0 || i % 257 == 3) {
        values[i] = -4;
      }
    }

    const Rmq rmq(values);
    expect_valid_bp_encoding(rmq, values.size());

    const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
        {0, 1},
        {0, size},
        {size - 1, size},
        {size / 3, std::min(size, size / 3 + 19)},
        {size / 2, std::min(size, size / 2 + 37)},
    };
    for (const auto& [left, right] : ranges) {
      ASSERT_LT(left, right);
      EXPECT_EQ(rmq.arg_min(left, right),
                naive_arg_min(std::span<const int>(values), left, right,
                              std::less<int>()))
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqCartesianHybridBTree, BorderMinimaOutsideQueryFallBackToMiddle) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  constexpr std::size_t kLeaf = 512;

  {
    std::vector<int> values(3 * kLeaf, 100000);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] += static_cast<int>((i * 17 + i / 9) % 101);
    }
    values[0] = -100000;
    values[13] = -100;
    values[kLeaf + 37] = -1000;
    values[2 * kLeaf + 53] = -500;

    const Rmq rmq{std::span<const int>(values)};
    const std::size_t left = 1;
    const std::size_t right = values.size() - 1;
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::less<int>()));
  }

  {
    std::vector<int> values(3 * kLeaf, 100000);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] += static_cast<int>((i * 19 + i / 11) % 107);
    }
    values[17] = -500;
    values[kLeaf + 41] = -1000;
    values[2 * kLeaf + 59] = -100;
    values[values.size() - 1] = -100000;

    const Rmq rmq{std::span<const int>(values)};
    const std::size_t left = 1;
    const std::size_t right = values.size() - 1;
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::less<int>()));
  }
}

TEST(RmqCartesianHybridBTree, DuplicateHeavyRandomDifferentialTo8193) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  std::mt19937_64 rng(20260610);
  std::uniform_int_distribution<int> value_dist(-3, 3);
  const std::vector<std::size_t> sizes = {
      1, 2, 17, 255, 256, 257, 1024, 4096, 8191, 8192, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = value_dist(rng);
      if ((i % 11) < 6) {
        values[i] = 0;
      }
    }

    const Rmq rmq(values);
    std::uniform_int_distribution<std::size_t> width_dist(1, size);
    for (std::size_t query = 0; query < 2000; ++query) {
      const std::size_t width = width_dist(rng);
      std::uniform_int_distribution<std::size_t> left_dist(0, size - width);
      const std::size_t left = left_dist(rng);
      const std::size_t right = left + width;
      const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                                 left, right, std::less<int>());
      EXPECT_EQ(rmq.arg_min(left, right), expected)
          << "range=[" << left << "," << right << ")";
    }
  }
}

TEST(RmqCartesianHybridBTree, DepthBackendDirectPaths) {
  using HighSparseBackend =
      pixie::rmq::detail::HybridBTreePlusMinusOne<std::size_t, 512, true, 1>;
  constexpr std::size_t kDepthCount = 3 * 512 + 17;

  const std::vector<std::uint64_t> words =
      make_packed_delta_bits(kDepthCount - 1);
  const std::vector<std::int64_t> depths = depths_from_delta_bits(
      std::span<const std::uint64_t>(words), kDepthCount);
  const HighSparseBackend high_sparse(std::span<const std::uint64_t>(words),
                                      kDepthCount);

  EXPECT_EQ(high_sparse.arg_min(kDepthCount, kDepthCount),
            HighSparseBackend::npos);
  EXPECT_EQ(high_sparse.arg_min(37, 38), 37u);
  EXPECT_GT(high_sparse.memory_usage_bytes(), sizeof(HighSparseBackend));

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kDepthCount},
      {0, 512},
      {9, 80},
      {120, 260},
      {384, 511},
      {511, 1027},
      {9, 1400},
      {700, kDepthCount - 3},
      {kDepthCount - 65, kDepthCount},
  };
  for (const auto& [left, right] : ranges) {
    const std::size_t expected =
        naive_depth_arg_min(std::span<const std::int64_t>(depths), left, right);
    EXPECT_EQ(high_sparse.arg_min(left, right), expected)
        << "high sparse depth range=[" << left << "," << right << ")";
  }

  EXPECT_THROW(HighSparseBackend(std::span<const std::uint64_t>(), 65),
               std::invalid_argument);
}

TEST(RmqCartesianHybridBTree, BpStorageIsCacheLineAligned) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  constexpr std::size_t kLeafWords = 512 / 64;
  const auto expect_aligned = [](const Rmq& rmq) {
    const std::span<const std::uint64_t> words = rmq.bp_words();
    ASSERT_FALSE(words.empty());
    const auto base = reinterpret_cast<std::uintptr_t>(
        static_cast<const void*>(words.data()));
    EXPECT_EQ(base % pixie::kAlignedStorageLineBytes, 0u);

    for (std::size_t word = 0; word < words.size(); word += kLeafWords) {
      EXPECT_EQ((base + word * sizeof(std::uint64_t)) %
                    pixie::kAlignedStorageLineBytes,
                0u)
          << "leaf_start_word=" << word;
    }
  };

  const std::vector<std::size_t> sizes = {
      1, 255, 256, 257, 4096, 8193,
  };

  for (const std::size_t size : sizes) {
    SCOPED_TRACE(size);
    std::vector<int> values(size);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<int>((i * 31 + i / 3) % 97);
    }

    const Rmq rmq(values);
    Rmq copied(rmq);
    Rmq assigned;
    assigned = copied;
    Rmq moved(std::move(copied));
    expect_aligned(rmq);
    expect_aligned(assigned);
    expect_aligned(moved);
  }
}

TEST(RmqCartesianHybridBTree, TopSparseOverlayBoundaryRanges) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  constexpr std::size_t kTopBlock = 4096;
  std::vector<int> values(3 * kTopBlock + 77, 1000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = 1000 + static_cast<int>((i * 17 + i / 7) % 89);
  }
  values[0] = -1000;
  values[211] = -700;
  values[kTopBlock + 123] = -900;
  values[2 * kTopBlock + 300] = -800;

  const Rmq original(values);
  Rmq copied(original);
  Rmq assigned;
  assigned = copied;
  Rmq moved(std::move(copied));

  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {0, kTopBlock + 100},
      {100, 2 * kTopBlock + 100},
      {kTopBlock - 20, 2 * kTopBlock + 20},
      {333, kTopBlock + 10},
      {kTopBlock + 1, kTopBlock + 2000},
      {2 * kTopBlock + 10, values.size()},
      {0, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    ASSERT_LT(left, right);
    const std::size_t expected = naive_arg_min(std::span<const int>(values),
                                               left, right, std::less<int>());
    EXPECT_EQ(original.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(assigned.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
    EXPECT_EQ(moved.arg_min(left, right), expected)
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqCartesianHybridBTree, TopSparseOverlayCapsBlockCount) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int>;
  constexpr std::size_t kMinBlock = Rmq::kMinTopSparseBlockSize;
  constexpr std::size_t kMaxBlocks = Rmq::kMaxTopSparseBlocks;
  constexpr std::size_t kLargestFixedBlockInput = kMinBlock * kMaxBlocks;

  EXPECT_EQ(Rmq::top_sparse_block_size_for(0), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(0), 0u);
  EXPECT_EQ(Rmq::top_sparse_block_size_for(1), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(1), 1u);
  EXPECT_EQ(Rmq::top_sparse_block_size_for(kLargestFixedBlockInput), kMinBlock);
  EXPECT_EQ(Rmq::top_sparse_block_count_for(kLargestFixedBlockInput),
            kMaxBlocks);
  EXPECT_GT(Rmq::top_sparse_block_size_for(kLargestFixedBlockInput + 1),
            kMinBlock);
  EXPECT_LE(Rmq::top_sparse_block_count_for(kLargestFixedBlockInput + 1),
            kMaxBlocks);

  const std::size_t huge = std::numeric_limits<std::size_t>::max() / 4;
  EXPECT_LE(Rmq::top_sparse_block_count_for(huge), kMaxBlocks);

  std::vector<int> values(3 * kMinBlock + 7, 0);
  const Rmq rmq(values);
  EXPECT_EQ(rmq.top_sparse_block_size(), kMinBlock);
  EXPECT_EQ(rmq.top_sparse_block_count(), 4u);
}

TEST(RmqCartesianHybridBTree, TopSparseOverlayComparatorMaximum) {
  using Rmq = pixie::rmq::CartesianHybridBTree<int, std::greater<int>>;
  constexpr std::size_t kTopBlock = 4096;
  std::vector<int> values(2 * kTopBlock + 33, -1000);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = -1000 + static_cast<int>((i * 13 + i / 5) % 71);
  }
  values[0] = 1000;
  values[kTopBlock + 17] = 900;
  values[kTopBlock + 18] = 900;

  const Rmq rmq(values);
  const std::vector<std::pair<std::size_t, std::size_t>> ranges = {
      {1, values.size()},
      {kTopBlock, kTopBlock + 100},
      {kTopBlock + 18, values.size()},
      {0, values.size()},
  };

  for (const auto& [left, right] : ranges) {
    ASSERT_LT(left, right);
    EXPECT_EQ(rmq.arg_min(left, right),
              naive_arg_min(std::span<const int>(values), left, right,
                            std::greater<int>()))
        << "range=[" << left << "," << right << ")";
  }
}

TEST(RmqCartesianHybridBTree, ComparatorCopyAndMove) {
  using MaxRmq = pixie::rmq::CartesianHybridBTree<int, std::greater<int>>;
  const std::vector<int> values = {1, 8, 8, 7, 8, 3, 8, 4, 4, 8};
  const MaxRmq original(values);
  MaxRmq copied(original);
  MaxRmq assigned;
  assigned = copied;
  MaxRmq moved(std::move(copied));

  check_all_ranges(original, std::span<const int>(values), std::greater<int>());
  check_all_ranges(assigned, std::span<const int>(values), std::greater<int>());
  check_all_ranges(moved, std::span<const int>(values), std::greater<int>());
  EXPECT_EQ(original.arg_min(0, values.size()), 1u);
}
