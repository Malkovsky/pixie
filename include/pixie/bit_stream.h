#pragma once

#include <cstdlib>
#include <vector>

namespace pixie {

class OutputBitStream {
 private:
  size_t size_;
  std::vector<uint64_t> data_;

 public:
  OutputBitStream() : size_(0) {}

  /**
   * @brief Writes one bit to the stream
   */
  OutputBitStream& operator<<(bool bit) {
    if (size_ % 64 == 0) {
      data_.push_back(bit);
    } else if (bit) {
      data_.back() |= 1ull << (size_ % 64);
    }
    size_++;
    return *this;
  }

  /**
   * @brief Writes bits of the integral number to the stream in little-endian
   */
  template <std::integral T>
  OutputBitStream& operator<<(T bits) {
    using UT = std::make_unsigned_t<T>;
    UT ubits = static_cast<UT>(bits);
    constexpr size_t length = sizeof(T) * 8;
    static_assert(length <= 64);
    if (size_ % 64 == 0) {
      data_.push_back(ubits);
    } else {
      const size_t prefix = std::min(length, 64 - (size_ % 64));
      data_.back() |= static_cast<uint64_t>(ubits & ((1ull << prefix) - 1))
                      << (size_ % 64);
      if (prefix < length) {
        data_.push_back(ubits >> prefix);
      }
    }
    size_ += length;
    return *this;
  }

  /**
   * @brief Returns the number of written bits
   *
   */
  size_t size() const { return size_; }

  /**
   * @brief Reserves memory for "size" bits
   *
   */
  void reserve(size_t size) { data_.reserve((size + 63) / 64); }

  /**
   * @brief Moves vector containing written bits. There must be no operator<<
   * after extract()
   */
  std::vector<uint64_t> extract() { return std::move(data_); }
};

}  // namespace pixie
