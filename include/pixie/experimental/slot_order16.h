#pragma once

/**
 * @file slot_order16.h
 * @brief Packed ordering primitive for the intermediate-node experiment.
 * @details Primitive probes exclude subtree measures, allocation, and tree
 * traversal; they do not establish sequence-level performance on their own.
 * Current measurements and reproduction commands are retained in
 * benchmarks/sequence_snapshot.md.
 */

#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace pixie::experimental {

/**
 * @brief Experimental logical ordering of sixteen stable physical slots.
 * @details Owns no payload. Four-bit entries map zero-based logical positions
 * to physical slots; a bitmap tracks occupancy independently of order. Callers
 * own the slots and must initialize a newly inserted slot before reading it.
 * All operations are allocation-free. Invalid inputs violate preconditions
 * checked by assertions. This primitive does not maintain subtree lengths.
 *
 * The tree experiment places this representation above leaves. The separate
 * GroupedSlotOrder256 experiment combines sixteen groups inside one node,
 * with boundary regrouping for arbitrary cross-group orderings.
 */
class SlotOrder16 {
 public:
  /** @brief Maximum occupied child count. */
  static constexpr std::size_t capacity = 16;
  /** @brief Construct empty, or identity order of count slots; count <= 16. */
  explicit SlotOrder16(std::size_t count = 0) noexcept
      : order_(0xfedcba9876543210ULL & mask(unsigned(count * 4))),
        occupied_(std::uint16_t(mask(unsigned(count)))),
        count_(std::uint8_t(count)) {
    assert(count <= 16);
  }

  /** @brief Number of occupied slots, in [0, 16]. */
  std::size_t size() const noexcept { return count_; }

  /** @brief Physical slot for logical index; requires index < size(). */
  unsigned operator[](std::size_t index) const noexcept {
    assert(index < size());
    return unsigned((order_ >> (4 * index)) & 15);
  }

  /** @brief Bit s is set exactly when physical slot s is occupied. */
  std::uint16_t occupied() const noexcept { return occupied_; }

  /** @brief Check logical mapping uniqueness against physical occupancy. */
  bool valid() const noexcept {
    if (size() > capacity) {
      return false;
    }
    unsigned seen = 0;
    for (std::size_t i = 0; i < size(); ++i) {
      const auto bit = 1u << (*this)[i];
      if (seen & bit) {
        return false;
      }
      seen |= bit;
    }
    return seen == occupied_;
  }

  /**
   * @brief Rotate logical [left,right) without moving any supplied pointers.
   * @details Preconditions and distance reduction match rotate_left().
   * @return Zero physical pointer writes; only ordering metadata changes.
   */
  std::size_t rotate_children(std::array<void*, capacity>&,
                              std::size_t left,
                              std::size_t right,
                              std::size_t distance) noexcept {
    rotate_left(left, right, distance);
    return 0;
  }

  /**
   * @brief Rotate logical [left, right) left by distance modulo its length.
   * @details Empty ranges are unchanged. Requires left <= right <= size().
   */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return;
    }
    const auto bits = unsigned(4 * length);
    const auto shift = unsigned(4 * left);
    const auto region_mask = mask(bits);
    const auto region = (order_ >> shift) & region_mask;
    const auto rotated =
        (region >> (4 * distance)) | (region << (4 * (length - distance)));
    order_ =
        (order_ & ~(region_mask << shift)) | ((rotated & region_mask) << shift);
  }

  /**
   * @brief Insert a free slot at logical position and return its physical
   * index.
   * @details Requires position <= size() < 16. Existing slot indices remain
   * stable; the caller initializes the returned physical slot.
   */
  unsigned insert(std::size_t position) noexcept {
    assert(position <= size() && size() < 16);
    const auto slot =
        unsigned(std::countr_zero(unsigned(std::uint16_t(~occupied_))));
    order_ |= std::uint64_t(slot) << (4 * count_);
    occupied_ |= std::uint16_t(1u << slot);
    ++count_;
    rotate_left(position, size(), size() - position - 1);
    return slot;
  }

  /**
   * @brief Remove a logical entry and return its now-free physical slot.
   * @details Requires position < size(). Does not destroy the caller's payload.
   */
  unsigned erase(std::size_t position) noexcept {
    assert(position < size());
    const auto slot = (*this)[position];
    rotate_left(position, size(), 1);
    --count_;
    order_ &= mask(unsigned(4 * count_));
    occupied_ &= std::uint16_t(~(1u << slot));
    return slot;
  }

 private:
  static constexpr std::uint64_t mask(unsigned bits) noexcept {
    return bits >= 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << bits) - 1;
  }
  std::uint64_t order_;
  std::uint16_t occupied_;
  std::uint8_t count_;
};

}  // namespace pixie::experimental
