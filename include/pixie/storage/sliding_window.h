#pragma once

/**
 * @file sliding_window.h
 * @brief Fixed-capacity ring-backed sliding byte storage.
 */

#include <pixie/storage.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Owning fixed-capacity storage over a monotonically positioned window.
 *
 * @details The complete backing ring is allocated during construction and is
 * never resized. Extending advances `end_position()` and evicts the oldest
 * bytes as necessary without initializing the newly exposed storage. Mutable
 * and const access use absolute logical positions and return at most two
 * physical spans. Requesting a future mutable range extends the window through
 * that range automatically; the caller remains responsible for assigning its
 * contents.
 *
 * Segments remain valid while their complete logical range remains inside the
 * working window. Moving or destroying the storage invalidates all segments.
 */
class SlidingWindowStorage : public StorageBase<SlidingWindowStorage> {
 public:
  /** @brief Construct an empty window with immutable byte capacity. */
  explicit SlidingWindowStorage(std::size_t capacity_bytes)
      : capacity_bytes_(validate_capacity(capacity_bytes)),
        data_(capacity_bytes_) {}

  SlidingWindowStorage(const SlidingWindowStorage&) = default;

  /** @brief Transfer the fixed backing allocation. */
  SlidingWindowStorage(SlidingWindowStorage&& other) noexcept
      : capacity_bytes_(std::exchange(other.capacity_bytes_, 0)),
        data_(std::move(other.data_)),
        begin_position_(std::exchange(other.begin_position_, 0)),
        end_position_(std::exchange(other.end_position_, 0)) {}

  SlidingWindowStorage& operator=(const SlidingWindowStorage&) = delete;
  SlidingWindowStorage& operator=(SlidingWindowStorage&&) = delete;

  /** @brief Return the immutable logical ring capacity in bytes. */
  std::size_t capacity_bytes() const { return capacity_bytes_; }

  /** @brief Return the number of bytes currently retained. */
  std::size_t size_bytes_impl() const {
    return static_cast<std::size_t>(end_position_ - begin_position_);
  }

  /** @brief Return the oldest retained logical byte position. */
  position_type begin_position_impl() const { return begin_position_; }

  /** @brief Return the position one past the newest retained byte. */
  position_type end_position_impl() const { return end_position_; }

  /** @brief Return the bytes in the fixed backing allocation. */
  std::size_t allocated_bytes_impl() const { return data_.capacity(); }

  /**
   * @brief Advance the logical end by a number of bytes.
   * @details Newly exposed bytes are not initialized. Advancing beyond the
   * capacity evicts the oldest logical bytes while retaining the fixed backing
   * allocation.
   * @throws std::length_error if the logical position would overflow or if a
   * nonzero extension is requested for a zero-capacity window.
   */
  void extend(std::size_t count_bytes) {
    if (count_bytes == 0) {
      return;
    }
    if (capacity_bytes_ == 0) {
      throw std::length_error("Cannot extend a zero-capacity sliding window");
    }
    if (count_bytes >
        std::numeric_limits<position_type>::max() - end_position_) {
      throw std::length_error("Sliding-window position overflow");
    }

    const position_type new_end = end_position_ + count_bytes;
    begin_position_ = new_end > capacity_bytes_ ? new_end - capacity_bytes_ : 0;
    end_position_ = new_end;
  }

  /**
   * @brief Extend through a requested future writable range when necessary.
   * @throws std::out_of_range if the requested range cannot fit in the window.
   * @throws std::length_error if the requested logical range overflows.
   */
  void prepare_segments_impl(position_type position, std::size_t count_bytes) {
    if (count_bytes == 0 || (position <= end_position_ &&
                             count_bytes <= end_position_ - position)) {
      return;
    }
    if (count_bytes > std::numeric_limits<position_type>::max() - position) {
      throw std::length_error("Sliding-window range overflow");
    }

    const position_type requested_end = position + count_bytes;
    if (count_bytes > capacity_bytes_ || position < begin_position_) {
      throw std::out_of_range(
          "Storage range cannot fit inside the working window");
    }
    const position_type extension = requested_end - end_position_;
    if constexpr (sizeof(position_type) > sizeof(std::size_t)) {
      if (extension > std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("Sliding-window extension is too large");
      }
    }
    extend(static_cast<std::size_t>(extension));
  }

  /** @brief Return a retained logical range as writable physical segments. */
  SplitSpan<std::byte> segments_impl(position_type position,
                                     std::size_t count_bytes) {
    return physical_segments(position, count_bytes,
                             std::span<std::byte>(data_));
  }

  /** @brief Return a retained logical range as read-only physical segments. */
  SplitSpan<const std::byte> segments_impl(position_type position,
                                           std::size_t count_bytes) const {
    return physical_segments(position, count_bytes,
                             std::span<const std::byte>(data_));
  }

 private:
  static std::size_t validate_capacity(std::size_t capacity_bytes) {
    if constexpr (sizeof(std::size_t) > sizeof(position_type)) {
      if (capacity_bytes > std::numeric_limits<position_type>::max()) {
        throw std::length_error("Sliding-window capacity is too large");
      }
    }
    return capacity_bytes;
  }

  template <class Byte>
  SplitSpan<Byte> physical_segments(position_type position,
                                    std::size_t count_bytes,
                                    std::span<Byte> data) const {
    if (count_bytes == 0) {
      return {};
    }
    const std::size_t physical_begin =
        static_cast<std::size_t>(position % capacity_bytes_);
    const std::size_t first_size =
        std::min(count_bytes, capacity_bytes_ - physical_begin);
    return SplitSpan<Byte>(data.subspan(physical_begin, first_size),
                           data.first(count_bytes - first_size));
  }

  std::size_t capacity_bytes_;
  std::vector<std::byte> data_;
  position_type begin_position_ = 0;
  position_type end_position_ = 0;
};

}  // namespace pixie
