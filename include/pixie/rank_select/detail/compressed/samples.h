#pragma once

#include <pixie/bits.h>
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

/**
 * @brief Internal block Elias--Fano array of non-decreasing 64-bit integers.
 * @details Each group of up to 32 values stores a base, a unary high word,
 * and packed low bits contiguously. A separate word-offset directory has one
 * entry per group, with the low width in its upper six bits. Widths can be
 * zero; duplicates are valid. There is no child pointer or exception flag per
 * integer.
 */
class BlockedMonotoneSamples {
 public:
  BlockedMonotoneSamples() = default;

  explicit BlockedMonotoneSamples(std::span<const std::uint64_t> values) {
    if (!std::is_sorted(values.begin(), values.end())) {
      throw std::invalid_argument("Monotone samples must be non-decreasing");
    }
    const std::size_t groups =
        values.empty() ? 0 : 1 + (values.size() - 1) / 32;
    offsets_.resize(groups);
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
    storage_ = AlignedStorage(words * 64);
    auto packed = storage_.writable_words64();
    for (std::size_t group = 0; group < groups; ++group) {
      const std::size_t begin = group * 32;
      const std::size_t count =
          std::min<std::size_t>(32, values.size() - begin);
      auto* record = packed.data() + (offsets_[group] & kOffsetMask);
      record[0] = values[begin];
      const unsigned width =
          std::bit_width((values[begin + count - 1] - values[begin]) >> 5);
      for (std::size_t i = 0; i < count; ++i) {
        const std::uint64_t delta = values[begin + i] - record[0];
        record[1] |= std::uint64_t{1} << ((delta >> width) + i);
        write_monotone_low(record + 2, i, width, delta);
      }
    }
  }

  std::uint64_t operator[](std::size_t index) const {
    const std::uint64_t descriptor = offsets_[index / 32];
    const auto* record = words() + (descriptor & kOffsetMask);
    const std::size_t local = index % 32;
    const unsigned width = descriptor >> 58;
    const std::uint64_t high = ::select_64(record[1], local) - local;
    return record[0] + (high << width) +
           read_monotone_low(record + 2, local, width);
  }

  std::size_t allocated_bytes() const {
    return storage_.allocated_bytes() +
           offsets_.capacity() * sizeof(offsets_[0]);
  }

 protected:
  static constexpr std::uint64_t kOffsetMask = (std::uint64_t{1} << 58) - 1;
  AlignedStorage storage_;
  std::vector<std::uint64_t> offsets_;

  const std::uint64_t* words() const {
    // Construction allocates complete words in 64-byte-aligned storage.
    return reinterpret_cast<const std::uint64_t*>(storage_.as_lines().data());
  }
};

struct SampleWideGroup {
  std::uint64_t begin = 0;
  std::uint64_t mask = 0;
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
  BitmapMonotoneSamples(std::span<const std::uint64_t> values,
                        std::span<const std::uint64_t> wide = {})
      : Base(values), wide_(make_sample_wide_groups(values.size(), wide)) {}

  std::optional<std::size_t> wide_entry(std::size_t index) const {
    if (wide_.empty()) {
      return std::nullopt;
    }
    const auto& group = wide_[index / 32];
    return sample_wide_entry(group.begin, group.mask, index % 32);
  }
  std::size_t allocated_bytes() const {
    return Base::allocated_bytes() + wide_.capacity() * sizeof(SampleWideGroup);
  }

 private:
  std::vector<SampleWideGroup> wide_;
};

}  // namespace pixie::detail::compressed
