#pragma once

#include <pixie/storage/aligned.h>

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
 * for persistent serialization. `take_storage()` transfers the aligned
 * backing without copying, while `take_words()` materializes a standard
 * vector. Both reset the builder to an empty, reusable state.
 */
class PackedBitBuilder {
 private:
  std::size_t size_ = 0;
  AlignedStorage data_;

  static std::size_t bits_for_words(std::size_t words) {
    if (words > std::numeric_limits<std::size_t>::max() / 64) {
      throw std::length_error("Packed bit sequence is too large");
    }
    return words * 64;
  }

  void ensure_word_capacity(std::size_t required_words) {
    const std::size_t current_words =
        data_.size_bytes() / sizeof(std::uint64_t);
    if (required_words <= current_words) {
      return;
    }
    if (current_words == 0) {
      data_ = AlignedStorage(
          bits_for_words(std::max(required_words, kAlignedStorageLineWords64)));
      return;
    }
    const std::size_t doubled_words =
        current_words > std::numeric_limits<std::size_t>::max() / 2
            ? required_words
            : current_words * 2;
    const std::size_t grown_words = std::max(required_words, doubled_words);
    data_.resize(bits_for_words(grown_words));
  }

 public:
  /**
   * @brief Append one bit.
   */
  void write_bit(bool bit) {
    if (size_ % 64 == 0) {
      const std::size_t word = size_ / 64;
      ensure_word_capacity(word + 1);
      data_.writable_words64()[word] = static_cast<std::uint64_t>(bit);
    } else if (bit) {
      data_.writable_words64()[size_ / 64] |= 1ull << (size_ % 64);
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
      const std::size_t word = size_ / 64;
      ensure_word_capacity(word + 1);
      data_.writable_words64()[word] =
          width == 64 ? bits : bits & ((1ull << width) - 1);
    } else {
      const std::size_t prefix = std::min(width, 64 - offset);
      const std::uint64_t prefix_mask =
          prefix == 64 ? ~std::uint64_t{0} : (1ull << prefix) - 1;
      auto words = data_.writable_words64();
      words[size_ / 64] |= (bits & prefix_mask) << offset;
      if (prefix < width) {
        const std::size_t word = size_ / 64 + 1;
        ensure_word_capacity(word + 1);
        data_.writable_words64()[word] = bits >> prefix;
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
    ensure_word_capacity(words);
  }

  /**
   * @brief Materialize the packed words and reset this builder.
   */
  std::vector<std::uint64_t> take_words() {
    const std::size_t words =
        size_ / 64 + static_cast<std::size_t>(size_ % 64 != 0);
    const auto source = data_.as_words64().first(words);
    std::vector<std::uint64_t> result(source.begin(), source.end());
    size_ = 0;
    data_ = {};
    return result;
  }

  /**
   * @brief Transfer the packed words in aligned storage and reset this builder.
   */
  AlignedStorage take_storage() {
    const std::size_t words =
        size_ / 64 + static_cast<std::size_t>(size_ % 64 != 0);
    data_.resize(bits_for_words(words));
    size_ = 0;
    return std::exchange(data_, {});
  }
};

}  // namespace pixie
