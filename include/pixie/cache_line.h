#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

inline constexpr std::size_t kAlignedStorageLineBytes = 64;
inline constexpr std::size_t kAlignedStorageLineBits =
    kAlignedStorageLineBytes * 8;
inline constexpr std::size_t kAlignedStorageLineWords64 =
    kAlignedStorageLineBytes / sizeof(std::uint64_t);
inline constexpr std::size_t kAlignedStorageLineWords16 =
    kAlignedStorageLineBytes / sizeof(std::uint16_t);

/**
 * @brief A 64-byte aligned storage block.
 */
struct alignas(kAlignedStorageLineBytes) CacheLine {
  std::array<std::byte, kAlignedStorageLineBytes> data{};
};

static_assert(alignof(CacheLine) == kAlignedStorageLineBytes);
static_assert(sizeof(CacheLine) == kAlignedStorageLineBytes);
static_assert(kAlignedStorageLineBytes % sizeof(std::uint64_t) == 0);
static_assert(kAlignedStorageLineBytes % sizeof(std::uint16_t) == 0);

/**
 * @brief Aligned storage for 64-byte blocks.
 *
 * @details The constructor and resize accept a logical size in bits. Storage is
 * rounded up to a full 64-byte block, and all views expose the padded capacity.
 */
class AlignedStorage {
 private:
  std::vector<CacheLine> data_;

  static constexpr std::size_t LinesForBits(std::size_t bits) {
    return bits / kAlignedStorageLineBits +
           (bits % kAlignedStorageLineBits != 0);
  }

 public:
  AlignedStorage() = default;
  /**
   * @brief Construct storage for at least @p bits bits.
   */
  explicit AlignedStorage(std::size_t bits) : data_(LinesForBits(bits)) {}

  /**
   * @brief Resize storage to hold at least @p bits bits.
   */
  void resize(std::size_t bits) { data_.resize(LinesForBits(bits)); }

  /** @brief Padded storage capacity in bits. */
  std::size_t capacity_bits() const {
    return data_.size() * kAlignedStorageLineBits;
  }

  /** @brief Padded storage capacity in bytes. */
  std::size_t capacity_bytes() const {
    return data_.size() * kAlignedStorageLineBytes;
  }

  /** @brief Mutable view as cache lines. */
  std::span<CacheLine> AsLines() { return data_; }
  /** @brief Const view as cache lines. */
  std::span<const CacheLine> AsConstLines() const { return data_; }

  /**
   * @brief Mutable view as 64-bit words.
   */
  std::span<std::uint64_t> As64BitInts() {
    return std::span<std::uint64_t>(
        reinterpret_cast<std::uint64_t*>(data_.data()),
        data_.size() * kAlignedStorageLineWords64);
  }

  /** @brief Const view as 64-bit words. */
  std::span<const std::uint64_t> AsConst64BitInts() const {
    return std::span<const std::uint64_t>(
        reinterpret_cast<const std::uint64_t*>(data_.data()),
        data_.size() * kAlignedStorageLineWords64);
  }

  /**
   * @brief Mutable view as bytes.
   */
  std::span<std::byte> AsBytes() { return std::as_writable_bytes(AsLines()); }

  /** @brief Const view as bytes. */
  std::span<const std::byte> AsConstBytes() const {
    return std::as_bytes(AsConstLines());
  }

  /**
   * @brief Mutable view as 16-bit words.
   */
  std::span<std::uint16_t> As16BitInts() {
    return std::span<std::uint16_t>(
        reinterpret_cast<std::uint16_t*>(data_.data()),
        data_.size() * kAlignedStorageLineWords16);
  }

  /** @brief Const view as 16-bit words. */
  std::span<const std::uint16_t> AsConst16BitInts() const {
    return std::span<const std::uint16_t>(
        reinterpret_cast<const std::uint16_t*>(data_.data()),
        data_.size() * kAlignedStorageLineWords16);
  }
};

}  // namespace pixie
