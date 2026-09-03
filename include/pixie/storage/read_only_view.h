#pragma once

#include <pixie/storage.h>

#include <span>

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

  /** @brief Return the first logical byte position. */
  position_type begin_position_impl() const { return 0; }

  /** @brief Return the position one past the final logical byte. */
  position_type end_position_impl() const { return data_.size(); }

  /** @brief Return the viewed bytes. */
  std::span<const std::byte> as_bytes_impl() const { return data_; }

  /** @brief Return a checked viewed byte range as one physical segment. */
  SplitSpan<const std::byte> segments_impl(std::size_t offset_bytes,
                                           std::size_t count_bytes) const {
    return SplitSpan(data_.subspan(offset_bytes, count_bytes));
  }

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
   * @brief Deserialize a size-prefixed view and advance @p reader.
   *
   * @details The returned view references the reader's backing byte sequence,
   * which must remain alive and stable for the view's lifetime. The reader is
   * unchanged on failure.
   *
   * @throws std::invalid_argument if the size prefix or payload is truncated.
   * @throws std::length_error if the encoded size is not representable.
   */
  static ReadOnlyStorageView deserialize_impl(BinaryReader& reader) {
    const std::size_t size = reader.read_size();
    return ReadOnlyStorageView(reader.read_bytes(size));
  }

 private:
  std::span<const std::byte> data_;
};

}  // namespace pixie
