#pragma once

#include <cstddef>
#include <span>
#include <vector>

/**
 * @brief A simple struct to represent a aligned storage for a cache line.
 */
struct alignas(64) CacheLine {
  std::array<std::byte, 64> data;
};

/**
 * @brief A simple aligned storage for cache-line sized blocks.
 *
 * @details Provides typed views over the same underlying storage as cache
 * lines, 64-bit words, or bytes. All spans are contiguous and sized to the
 * total storage capacity.
 */
class AlignedStorage {
 private:
  std::vector<CacheLine> data_;

 public:
  AlignedStorage() = default;
  /**
   * @brief Construct storage for at least @p bits bytes, rounded up to 512
   * bits.
   */
  AlignedStorage(size_t bits) : data_((bits + 511) / 512) {}

  /**
   * @brief Resize storage to hold at least @p bits bits, rounded up to 512
   */
  void resize(size_t bits) { data_.resize((bits + 511) / 512); }
  /** @brief Mutable view as cache lines. */
  std::span<CacheLine> AsLines() { return data_; }
  /** @brief Const view as cache lines. */
  std::span<const CacheLine> AsConstLines() const { return data_; }

  /**
   * @brief Mutable view as 64-bit words.
   */
  std::span<uint64_t> As64BitInts() {
    return std::span<uint64_t>(reinterpret_cast<uint64_t*>(data_.data()),
                               data_.size() * 8);
  }

  /** @brief Const view as 64-bit words. */
  std::span<const uint64_t> AsConst64BitInts() const {
    return std::span<const uint64_t>(
        reinterpret_cast<const uint64_t*>(data_.data()), data_.size() * 8);
  }

  /**
   * @brief Mutable view as bytes.
   * @note Uses a byte pointer to the underlying storage.
   */
  std::span<std::byte> AsBytes() {
    return std::span<std::byte>(reinterpret_cast<std::byte*>(data_.data()),
                                data_.size() * 64);
  }

  /** @brief Const view as bytes. */
  std::span<const std::byte> AsConstBytes() const {
    return std::span<const std::byte>(
        reinterpret_cast<const std::byte*>(data_.data()), data_.size() * 64);
  }

  /**
   * @brief Mutable view as bytes.
   * @note Uses a byte pointer to the underlying storage.
   */
  std::span<std::uint16_t> As16BitInts() {
    return std::span<std::uint16_t>(
        reinterpret_cast<std::uint16_t*>(data_.data()), data_.size() * 32);
  }

  /** @brief Const view as bytes. */
  std::span<const std::uint16_t> AsConst16BitInts() const {
    return std::span<const std::uint16_t>(
        reinterpret_cast<const std::uint16_t*>(data_.data()),
        data_.size() * 32);
  }
};
