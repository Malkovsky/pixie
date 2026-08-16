#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Builder for a packed LSB-first bit sequence.
 *
 * @details This type is intended for constructing succinct bit vectors, not
 * for persistent serialization. `take_words()` transfers the packed words and
 * resets the builder to an empty, reusable state.
 */
class PackedBitBuilder {
 private:
  std::size_t size_ = 0;
  std::vector<std::uint64_t> data_;

 public:
  /**
   * @brief Append one bit.
   */
  void write_bit(bool bit) {
    if (size_ % 64 == 0) {
      data_.push_back(static_cast<std::uint64_t>(bit));
    } else if (bit) {
      data_.back() |= 1ull << (size_ % 64);
    }
    ++size_;
  }

  /**
   * @brief Append the low @p width bits of @p bits, least-significant bit
   * first.
   * @throws std::invalid_argument if @p width is greater than 64.
   * @throws std::length_error if the resulting bit count is not representable.
   */
  void write_bits(std::uint64_t bits, std::size_t width) {
    if (width > 64) {
      throw std::invalid_argument("Packed bit width is greater than 64");
    }
    if (width == 0) {
      return;
    }
    if (size_ > std::numeric_limits<std::size_t>::max() - width) {
      throw std::length_error("Packed bit sequence is too large");
    }

    const std::size_t offset = size_ % 64;
    if (offset == 0) {
      data_.push_back(width == 64 ? bits : bits & ((1ull << width) - 1));
    } else {
      const std::size_t prefix = std::min(width, 64 - offset);
      const std::uint64_t prefix_mask =
          prefix == 64 ? ~std::uint64_t{0} : (1ull << prefix) - 1;
      data_.back() |= (bits & prefix_mask) << offset;
      if (prefix < width) {
        data_.push_back(bits >> prefix);
      }
    }
    size_ += width;
  }

  /** @brief Return the number of appended bits. */
  std::size_t size_bits() const noexcept { return size_; }

  /**
   * @brief Reserve storage for at least @p size_bits bits.
   */
  void reserve_bits(std::size_t size_bits) {
    const std::size_t words =
        size_bits / 64 + static_cast<std::size_t>(size_bits % 64 != 0);
    data_.reserve(words);
  }

  /**
   * @brief Transfer the packed words and reset this builder.
   */
  std::vector<std::uint64_t> take_words() {
    size_ = 0;
    return std::exchange(data_, {});
  }
};

}  // namespace pixie
