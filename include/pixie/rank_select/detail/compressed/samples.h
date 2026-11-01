#pragma once

#include <pixie/bits.h>
#include <pixie/rank_select/detail/compressed/metadata_array.h>
#include <pixie/storage/aligned.h>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::detail::compressed {
inline std::uint64_t read_monotone_low(const std::uint64_t* words,
                                       std::size_t index,
                                       unsigned width) {
  if (width == 0) {
    return 0;
  }
  const std::size_t bit = index * width;
  const unsigned shift = bit % 64;
  std::uint64_t value = words[bit / 64] >> shift;
  if (shift + width > 64) {
    value |= words[bit / 64 + 1] << (64 - shift);
  }
  return value & ((std::uint64_t{1} << width) - 1);
}

inline void write_monotone_low(std::uint64_t* words,
                               std::size_t index,
                               unsigned width,
                               std::uint64_t value) {
  if (width == 0) {
    return;
  }
  value &= (std::uint64_t{1} << width) - 1;
  const std::size_t bit = index * width;
  const unsigned shift = bit % 64;
  words[bit / 64] |= value << shift;
  if (shift + width > 64) {
    words[bit / 64 + 1] |= value >> (64 - shift);
  }
}

// Both packed sample layouts use a base, unary high word and packed low bits.
// Callers first validate their record bounds, low offset and width (<= 59).
// Return the last decoded value to carry monotonicity across record boundaries.
inline uint64_t validate_monotone_record(const uint64_t* record,
                                         size_t low_offset,
                                         size_t count,
                                         unsigned width,
                                         uint64_t limit,
                                         uint64_t previous) {
  if (static_cast<size_t>(std::popcount(record[1])) != count) {
    throw std::invalid_argument("Invalid compressed unary sample");
  }
  for (size_t i = 0; i < count; ++i) {
    const auto high = ::select_64(record[1], i) - i;
    const auto delta =
        (high << width) + read_monotone_low(record + low_offset, i, width);
    if (high > (UINT64_MAX >> width) || record[0] > limit ||
        delta > limit - record[0] || record[0] + delta < previous) {
      throw std::invalid_argument("Invalid compressed sample value");
    }
    previous = record[0] + delta;
  }
  return previous;
}

/**
 * @brief Internal block Elias--Fano array of non-decreasing 64-bit integers.
 * @details Each group of up to 32 values stores a base, a unary high word,
 * and packed low bits contiguously. A separate word-offset directory has one
 * entry per group, with the low width in its upper six bits. Widths can be
 * zero; duplicates are valid. There is no child pointer or exception flag per
 * integer.
 */
template <StorageImplementation MetadataStorage = AlignedStorage>
class BlockedMonotoneSamples {
 public:
  using metadata_storage_type = MetadataStorage;

  BlockedMonotoneSamples() = default;
  bool operator==(const BlockedMonotoneSamples&) const = default;

  explicit BlockedMonotoneSamples(std::span<const std::uint64_t> values)
    requires MutableMetadataStorage<MetadataStorage>
  {
    if (!std::is_sorted(values.begin(), values.end())) {
      throw std::invalid_argument("Monotone samples must be non-decreasing");
    }
    const std::size_t groups =
        values.empty() ? 0 : 1 + (values.size() - 1) / 32;
    offsets_ = MetadataArray<std::uint64_t, MetadataStorage>(groups);
    std::size_t words = 0;
    for (std::size_t group = 0; group < groups; ++group) {
      const std::size_t begin = group * 32;
      const std::size_t count =
          std::min<std::size_t>(32, values.size() - begin);
      const unsigned width =
          std::bit_width((values[begin + count - 1] - values[begin]) >> 5);
      offsets_[group] = words | (std::uint64_t{width} << 58);
      const std::size_t record_words = 2 + (count * width + 63) / 64;
      if (words > std::numeric_limits<std::size_t>::max() / 64 - record_words) {
        throw std::length_error("Monotone samples are too large");
      }
      words += record_words;
    }
    storage_ = MetadataArray<std::uint64_t, MetadataStorage>(words);
    auto* packed = storage_.data();
    for (std::size_t group = 0; group < groups; ++group) {
      const std::size_t begin = group * 32;
      const std::size_t count =
          std::min<std::size_t>(32, values.size() - begin);
      const auto descriptor = offsets_[group];
      auto* record = packed + (descriptor & kOffsetMask);
      record[0] = values[begin];
      const unsigned width = descriptor >> 58;
      for (std::size_t i = 0; i < count; ++i) {
        const std::uint64_t delta = values[begin + i] - record[0];
        record[1] |= std::uint64_t{1} << ((delta >> width) + i);
        write_monotone_low(record + 2, i, width, delta);
      }
    }
  }

  std::uint64_t operator[](std::size_t index) const {
    const std::uint64_t descriptor = offsets_[index / 32];
    const auto* record = storage_.data() + (descriptor & kOffsetMask);
    const std::size_t local = index % 32;
    const unsigned width = descriptor >> 58;
    const std::uint64_t high = ::select_64(record[1], local) - local;
    return record[0] + (high << width) +
           read_monotone_low(record + 2, local, width);
  }

  std::size_t allocated_bytes() const {
    return storage_.allocated_bytes() + offsets_.allocated_bytes();
  }

  void serialize(BinaryWriter& writer) const {
    storage_.serialize(writer);
    offsets_.serialize(writer);
  }

  static BlockedMonotoneSamples deserialize(BinaryReader& reader)
    requires MutableMetadataStorage<MetadataStorage>
  {
    BlockedMonotoneSamples result;
    result.storage_ = decltype(storage_)::deserialize(reader);
    result.offsets_ = decltype(offsets_)::deserialize(reader);
    return result;
  }

  void validate(size_t count, uint64_t limit) const {
    const size_t groups = count / 32 + (count % 32 != 0);
    if (offsets_.size() != groups) {
      throw std::invalid_argument("Invalid compressed sample directory");
    }
    size_t offset = 0;
    uint64_t previous = 0;
    for (size_t g = 0; g < groups; ++g) {
      const size_t n = std::min(size_t{32}, count - g * 32);
      const auto width = offsets_[g] >> 58;
      const size_t words = 2 + (n * width + 63) / 64;
      if (width > 59 || (offsets_[g] & kOffsetMask) != offset ||
          words > storage_.size() - offset) {
        throw std::invalid_argument("Invalid compressed sample record");
      }
      const auto* record = storage_.data() + offset;
      previous = validate_monotone_record(record, 2, n, width, limit, previous);
      offset += words;
    }
    if (offset != storage_.size()) {
      throw std::invalid_argument("Unexpected compressed sample payload");
    }
  }

  bool empty() const { return offsets_.empty() && storage_.empty(); }

 protected:
  static constexpr std::uint64_t kOffsetMask = (std::uint64_t{1} << 58) - 1;
  MetadataArray<std::uint64_t, MetadataStorage> storage_;
  MetadataArray<std::uint64_t, MetadataStorage> offsets_;
};

struct SampleWideGroup {
  std::uint64_t begin = 0;
  std::uint64_t mask = 0;
  bool operator==(const SampleWideGroup&) const = default;
};

// An empty key list allocates no bitmap. Callers only query the map for mixed
// shapes; homogeneous, narrow and wide-prefix shapes retain their shortcuts.
inline std::vector<SampleWideGroup> make_sample_wide_groups(
    std::size_t count,
    std::span<const std::uint64_t> keys) {
  if (keys.empty()) {
    return {};
  }
  std::vector<SampleWideGroup> groups(count == 0 ? 0 : 1 + (count - 1) / 32);
  std::size_t next = 0;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    groups[g].begin = next;
    while (next < keys.size() && keys[next] / 32 == g) {
      assert(keys[next] < count);
      assert(next == 0 || keys[next - 1] < keys[next]);
      groups[g].mask |= std::uint64_t{1} << (keys[next] % 32);
      ++next;
    }
  }
  assert(next == keys.size());
  return groups;
}

inline std::optional<std::size_t> sample_wide_entry(std::uint64_t begin,
                                                    std::uint64_t mask,
                                                    unsigned local) {
  const auto bit = std::uint64_t{1} << local;
  if (!(mask & bit)) {
    return std::nullopt;
  }
  return begin + std::popcount(mask & (bit - 1));
}

/** @brief Blocked samples with an exact wide-interval map. */
template <class Base>
class BitmapMonotoneSamples : public Base {
 public:
  BitmapMonotoneSamples() = default;
  bool operator==(const BitmapMonotoneSamples&) const = default;

  void serialize(BinaryWriter& writer) const {
    Base::serialize(writer);
    wide_.serialize(writer, [](BinaryWriter& out, SampleWideGroup value) {
      out.write_u64(value.begin);
      out.write_u64(value.mask);
    });
  }

  static BitmapMonotoneSamples deserialize(BinaryReader& reader) {
    BitmapMonotoneSamples result;
    static_cast<Base&>(result) = Base::deserialize(reader);
    result.wide_ =
        decltype(wide_)::deserialize(reader, 16, [](BinaryReader& in) {
          return SampleWideGroup{in.read_u64(), in.read_u64()};
        });
    return result;
  }

  void validate(size_t count, uint64_t limit, size_t wide_count) const {
    Base::validate(count, limit);
    if (wide_.empty()) {
      if (wide_count) {
        throw std::invalid_argument("Missing wide sample map");
      }
      return;
    }
    if (wide_.size() != count / 32 + (count % 32 != 0)) {
      throw std::invalid_argument("Invalid wide sample map size");
    }
    size_t next = 0;
    for (size_t g = 0; g < wide_.size(); ++g) {
      const auto valid = std::min(size_t{32}, count - g * 32);
      if (wide_[g].begin != next || (wide_[g].mask >> valid) != 0) {
        throw std::invalid_argument("Invalid wide sample map");
      }
      next += std::popcount(wide_[g].mask);
    }
    if (next != wide_count) {
      throw std::invalid_argument("Invalid wide sample count");
    }
  }
  BitmapMonotoneSamples(std::span<const std::uint64_t> values,
                        std::span<const std::uint64_t> wide = {})
      : Base(values),
        wide_(std::span<const SampleWideGroup>(
            make_sample_wide_groups(values.size(), wide))) {}

  std::optional<std::size_t> wide_entry(std::size_t index) const {
    if (wide_.empty()) {
      return std::nullopt;
    }
    const auto& group = wide_[index / 32];
    return sample_wide_entry(group.begin, group.mask, index % 32);
  }
  std::size_t allocated_bytes() const {
    return Base::allocated_bytes() + wide_.allocated_bytes();
  }

 private:
  MetadataArray<SampleWideGroup, typename Base::metadata_storage_type> wide_;
};

}  // namespace pixie::detail::compressed
