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
 * never resized. Appending advances `end_position()` and evicts the oldest
 * bytes as necessary. Mutable and const access use absolute logical positions
 * and return at most two physical spans.
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
   * @brief Append bytes and evict the oldest bytes beyond capacity.
   * @throws std::length_error if the logical position would overflow or if a
   * nonempty sequence is appended to a zero-capacity window.
   */
  void append(std::span<const std::byte> bytes) {
    if (bytes.empty()) {
      return;
    }
    if (capacity_bytes_ == 0) {
      throw std::length_error(
          "Cannot append to a zero-capacity sliding window");
    }
    if (bytes.size() >
        std::numeric_limits<position_type>::max() - end_position_) {
      throw std::length_error("Sliding-window position overflow");
    }

    const position_type new_end = end_position_ + bytes.size();
    if (bytes.size() >= capacity_bytes_) {
      const std::span<const std::byte> suffix = bytes.last(capacity_bytes_);
      const position_type suffix_position = new_end - capacity_bytes_;
      copy_into_ring(suffix_position, suffix);
      begin_position_ = suffix_position;
    } else {
      copy_into_ring(end_position_, bytes);
      const position_type retained_size =
          std::min<position_type>(new_end, capacity_bytes_);
      begin_position_ = std::max(begin_position_, new_end - retained_size);
    }
    end_position_ = new_end;
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

  void copy_into_ring(position_type position,
                      std::span<const std::byte> bytes) {
    SplitSpan<std::byte> destination =
        physical_segments(position, bytes.size(), std::span<std::byte>(data_));
    std::size_t copied = 0;
    for (const std::span<std::byte> segment : destination) {
      std::ranges::copy(bytes.subspan(copied, segment.size()), segment.begin());
      copied += segment.size();
    }
  }

  std::size_t capacity_bytes_;
  std::vector<std::byte> data_;
  position_type begin_position_ = 0;
  position_type end_position_ = 0;
};

}  // namespace pixie
