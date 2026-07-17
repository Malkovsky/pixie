#pragma once

#include <pixie/storage.h>

#include <cstring>
#include <span>
#include <stdexcept>

namespace pixie {

/**
 * @brief A non-owning, read-only view of a byte sequence.
 *
 * @details The caller must keep the backing storage alive and must not resize
 * it while this view or a view derived from it is in use.
 */
class ReadOnlyStorageView : public StorageBase<ReadOnlyStorageView> {
 public:
  ReadOnlyStorageView() = default;

  /** @brief Construct a view over @p data. */
  explicit ReadOnlyStorageView(std::span<const std::byte> data) : data_(data) {}

  /** @brief Return the number of viewed bytes. */
  std::size_t size_bytes_impl() const { return data_.size(); }

  /** @brief Return the viewed bytes. */
  std::span<const std::byte> as_bytes_impl() const { return data_; }

  /** @brief Return a checked read-only byte subrange. */
  ReadOnlyStorageView view_impl(std::size_t offset_bytes,
                                std::size_t count_bytes) const {
    if (offset_bytes > data_.size() ||
        count_bytes > data_.size() - offset_bytes) {
      throw std::out_of_range("Storage view is outside the backing storage");
    }
    return ReadOnlyStorageView(data_.subspan(offset_bytes, count_bytes));
  }

  /**
   * @brief Deserialize a size-prefixed view and advance @p data.
   * @throws std::invalid_argument if the size prefix or payload is truncated.
   */
  static ReadOnlyStorageView deserialize(std::span<const std::byte>& data) {
    if (data.size() < sizeof(std::size_t)) {
      throw std::invalid_argument("Truncated storage size prefix");
    }
    std::size_t size = 0;
    std::memcpy(&size, data.data(), sizeof(size));
    if (size > data.size() - sizeof(size)) {
      throw std::invalid_argument("Truncated storage payload");
    }
    ReadOnlyStorageView result(data.subspan(sizeof(size), size));
    data = data.subspan(sizeof(size) + size);
    return result;
  }

 private:
  std::span<const std::byte> data_;
};

}  // namespace pixie
