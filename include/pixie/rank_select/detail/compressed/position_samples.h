#pragma once

#include <pixie/rank_select/detail/compressed/samples.h>

namespace pixie::detail::compressed {

/**
 * @brief Exact select positions in fixed-stride Elias--Fano groups of 32.
 * @details Each group interleaves a base, unary high word, optional wide-index
 * prefix/mask, and packed low bits. A shared maximum low width pads every group
 * uniformly; there is no offset directory or duplicate coarse sample array.
 * Duplicates, partial final groups and all uint64_t values are representable.
 * Internal callers supply sorted wide keys below values.size(). Query indexes
 * must exist; bounds() additionally requires a following terminal sample.
 */
template <StorageImplementation MetadataStorage = AlignedStorage>
class PositionSamples {
  static constexpr std::size_t kHeaderWords = 2;

 public:
  PositionSamples() = default;
  bool operator==(const PositionSamples&) const = default;

  void serialize(BinaryWriter& writer) const {
    writer.write_size(stride_);
    writer.write_size(low_offset_);
    writer.write_size(width_);
    storage_.serialize(writer);
  }

  static PositionSamples deserialize(BinaryReader& reader) {
    PositionSamples result;
    result.stride_ = reader.read_size();
    result.low_offset_ = reader.read_size();
    const auto width = reader.read_size();
    if (width > 59) {
      throw std::invalid_argument("Invalid position sample width");
    }
    result.width_ = static_cast<unsigned>(width);
    result.storage_ = decltype(storage_)::deserialize(reader);
    return result;
  }

  void validate(size_t count, uint64_t limit, size_t wide_count) const {
    if (count == 0) {
      if (!storage_.empty() || stride_ != 0 || low_offset_ != 2 ||
          width_ != 0 || wide_count) {
        throw std::invalid_argument("Invalid empty position samples");
      }
      return;
    }
    const size_t groups = count / 32 + (count % 32 != 0);
    if ((low_offset_ != 2 && low_offset_ != 4) ||
        stride_ != low_offset_ + (32 * width_ + 63) / 64 ||
        storage_.size() / stride_ != groups || storage_.size() % stride_) {
      throw std::invalid_argument("Invalid position sample layout");
    }
    uint64_t previous = 0;
    size_t next = 0;
    for (size_t g = 0; g < groups; ++g) {
      const auto* record = group(g);
      const auto n = std::min(size_t{32}, count - g * 32);
      if (static_cast<size_t>(std::popcount(record[1])) != n) {
        throw std::invalid_argument("Invalid position unary sample");
      }
      if (low_offset_ == 4) {
        if (record[2] != next || (record[3] >> n)) {
          throw std::invalid_argument("Invalid position wide map");
        }
        next += std::popcount(record[3]);
      }
      for (size_t i = 0; i < n; ++i) {
        const auto high = ::select_64(record[1], i) - i;
        const auto delta = (high << width_) +
                           read_monotone_low(record + low_offset_, i, width_);
        if (high > (UINT64_MAX >> width_) || record[0] > limit ||
            delta > limit - record[0] || record[0] + delta < previous) {
          throw std::invalid_argument("Invalid position sample value");
        }
        previous = record[0] + delta;
      }
    }
    if (next != wide_count) {
      throw std::invalid_argument("Invalid position wide count");
    }
  }

  explicit PositionSamples(std::span<const std::uint64_t> values,
                           std::span<const std::uint64_t> wide = {})
    requires MutableMetadataStorage<MetadataStorage>
  {
    if (!std::is_sorted(values.begin(), values.end())) {
      throw std::invalid_argument("Monotone samples must be non-decreasing");
    }
    const std::size_t groups =
        values.empty() ? 0 : 1 + (values.size() - 1) / 32;
    for (std::size_t group = 0; group < groups; ++group) {
      width_ = std::max(width_, group_width(values, group));
    }
    const auto wide_groups = make_sample_wide_groups(values.size(), wide);
    low_offset_ = kHeaderWords + (wide_groups.empty() ? 0 : 2);
    stride_ = low_offset_ + (32 * width_ + 63) / 64;
    if (groups > std::numeric_limits<std::size_t>::max() / (stride_ * 64)) {
      throw std::length_error("Fixed sample allocation too large");
    }
    storage_ = MetadataArray<std::uint64_t, MetadataStorage>(groups * stride_);
    auto* data = storage_.data();
    for (std::size_t group = 0; group < groups; ++group) {
      auto* record = data + group * stride_;
      record[0] = values[group * 32];
      if (!wide_groups.empty()) {
        record[kHeaderWords] = wide_groups[group].begin;
        record[kHeaderWords + 1] = wide_groups[group].mask;
      }
      const auto count = std::min<std::size_t>(32, values.size() - group * 32);
      for (std::size_t i = 0; i < count; ++i) {
        const auto delta = values[group * 32 + i] - record[0];
        record[1] |= std::uint64_t{1} << ((delta >> width_) + i);
        write_monotone_low(record + low_offset_, i, width_, delta);
      }
    }
  }

  std::uint64_t operator[](std::size_t index) const {
    return decode(group(index / 32), index % 32);
  }

  std::pair<std::uint64_t, std::uint64_t> bounds(std::size_t index) const {
    if (index % 32 == 31) {
      return {(*this)[index], (*this)[index + 1]};
    }
    const auto* record = group(index / 32);
    return {decode(record, index % 32), decode(record, index % 32 + 1)};
  }

  std::optional<std::size_t> wide_entry(std::size_t index) const {
    if (low_offset_ == kHeaderWords) {
      return std::nullopt;
    }
    assert(low_offset_ == kHeaderWords + 2);
    const auto* record = group(index / 32);
    return sample_wide_entry(record[kHeaderWords], record[kHeaderWords + 1],
                             index % 32);
  }

  std::size_t allocated_bytes() const { return storage_.allocated_bytes(); }

 private:
  MetadataArray<std::uint64_t, MetadataStorage> storage_;
  std::size_t stride_ = 0;
  std::size_t low_offset_ = kHeaderWords;
  unsigned width_ = 0;

  static unsigned group_width(std::span<const std::uint64_t> values,
                              std::size_t group) {
    const auto end = std::min<std::size_t>(32, values.size() - group * 32);
    return std::bit_width((values[group * 32 + end - 1] - values[group * 32]) >>
                          5);
  }

  const std::uint64_t* group(std::size_t index) const {
    return storage_.data() + index * stride_;
  }

  std::uint64_t decode(const std::uint64_t* record, unsigned index) const {
    const auto high = ::select_64(record[1], index) - index;
    return record[0] + (high << width_) +
           read_monotone_low(record + low_offset_, index, width_);
  }
};

}  // namespace pixie::detail::compressed
