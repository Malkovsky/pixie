#pragma once

#include <pixie/storage.h>
#include <pixie/storage/read_only_view.h>

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

/** @brief A 64-byte aligned storage block. */
struct alignas(kAlignedStorageLineBytes) CacheLine {
  std::array<std::byte, kAlignedStorageLineBytes> data{};
};

static_assert(alignof(CacheLine) == kAlignedStorageLineBytes);
static_assert(sizeof(CacheLine) == kAlignedStorageLineBytes);

/**
 * @brief Owning storage rounded up to 64-byte blocks.
 *
 * @details Construction and resize accept a logical bit count. All exposed
 * views cover the padded allocation. Resizing or destroying this object
 * invalidates its read-only views.
 */
class AlignedStorage : public StorageBase<AlignedStorage> {
 public:
  AlignedStorage() = default;

  /** @brief Construct storage for at least @p size_bits bits. */
  explicit AlignedStorage(std::size_t size_bits)
      : data_(lines_for_bits(size_bits)) {}

  /** @brief Return the padded allocation size in bytes. */
  std::size_t size_bytes_impl() const {
    return data_.size() * kAlignedStorageLineBytes;
  }

  /** @brief Return the padded allocation as read-only bytes. */
  std::span<const std::byte> as_bytes_impl() const {
    return std::as_bytes(std::span<const CacheLine>(data_));
  }

  /** @brief Return a checked read-only byte subrange. */
  ReadOnlyStorageView view_impl(std::size_t offset_bytes,
                                std::size_t count_bytes) const {
    if (offset_bytes > size_bytes_impl() ||
        count_bytes > size_bytes_impl() - offset_bytes) {
      throw std::out_of_range("Storage view is outside the allocation");
    }
    return ReadOnlyStorageView(
        as_bytes_impl().subspan(offset_bytes, count_bytes));
  }

  /** @brief Resize to hold at least @p size_bits bits. */
  void resize_impl(std::size_t size_bits) {
    data_.resize(lines_for_bits(size_bits));
  }

  /** @brief Return writable allocation bytes. */
  std::span<std::byte> writable_bytes_impl() {
    return std::as_writable_bytes(std::span<CacheLine>(data_));
  }

  /** @brief Return writable allocation as 16-bit words. */
  std::span<std::uint16_t> writable_words16_impl() {
    return {reinterpret_cast<std::uint16_t*>(data_.data()),
            data_.size() * kAlignedStorageLineWords16};
  }

  /** @brief Return writable allocation as 64-bit words. */
  std::span<std::uint64_t> writable_words64_impl() {
    return {reinterpret_cast<std::uint64_t*>(data_.data()),
            data_.size() * kAlignedStorageLineWords64};
  }

  /** @brief Return bytes reserved by the underlying vector. */
  std::size_t allocated_bytes_impl() const {
    return data_.capacity() * kAlignedStorageLineBytes;
  }

  /** @brief Request release of unused vector capacity. */
  void shrink_to_fit_impl() { data_.shrink_to_fit(); }

  /** @brief Return mutable cache-line blocks. */
  std::span<CacheLine> as_lines() { return data_; }

  /** @brief Return read-only cache-line blocks. */
  std::span<const CacheLine> as_lines() const { return data_; }

 private:
  static constexpr std::size_t lines_for_bits(std::size_t size_bits) {
    return size_bits / kAlignedStorageLineBits +
           (size_bits % kAlignedStorageLineBits != 0);
  }

  std::vector<CacheLine> data_;
};

}  // namespace pixie
