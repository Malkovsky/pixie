#pragma once

/**
 * @file grouped_slot_order256.h
 * @brief Single-node 16x16 child-order experiment.
 * @details The integrated tree uses fanout 256 throughout, with grouped
 * addressing only at leaf parents. This is not a mixed-fanout tree. Grouped
 * metadata occupies 176 bytes; node allocations are 4272 bytes versus 4096
 * for direct F256. Fresh nodes use direct addressing until a mapped rotation.
 * Current measurements and workload trade-offs are retained in
 * benchmarks/sequence_snapshot.md.
 */

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace pixie::experimental {

/**
 * @brief Experimental ordering of 256 pointer slots in one 16x16 node.
 * @details Seventeen packed 64-bit permutations select the physical group and
 * the physical slot within it. Per-group occupancy bitmaps describe live
 * physical slots. Only the final logical group may be partial. This is one
 * node allocation, not sixteen separately allocated intermediate nodes.
 *
 * Whole-group and within-group rotations move only ordering metadata. A
 * rotation of full groups with a non-group-aligned distance exchanges at most
 * eight pointers per group. Other boundary combinations use bounded scratch
 * pointer storage. No operation allocates, owns, or destroys pointees; callers
 * supply live slots and retain ownership. Subtree sizes are maintained by the
 * caller in logical order. Invalid arguments violate asserted preconditions.
 */
class GroupedSlotOrder256 {
 public:
  /** @brief Maximum occupied child count. */
  static constexpr std::size_t capacity = 256;

  /** @brief Identity ordering of count live slots; requires count <= 256. */
  explicit GroupedSlotOrder256(std::size_t count = 0) noexcept
      : groups_(identity & mask(unsigned(4 * ((count + 15) / 16)))),
        count_(std::uint16_t(count)) {
    assert(count <= capacity);
    for (std::size_t group = 0; group < (count + 15) / 16; ++group) {
      const auto n = count - group * 16 < 16 ? count - group * 16 : 16;
      slots_[group] = identity & mask(unsigned(n * 4));
      occupied_[group] = std::uint16_t(mask(unsigned(n)));
    }
  }

  /** @brief Number of live logical slots in [0, 256]. */
  std::size_t size() const noexcept { return count_; }

  /** @brief Physical slot for zero-based index; requires index < size(). */
  unsigned operator[](std::size_t index) const noexcept {
    assert(index < size());
    const auto group = physical_group(index / 16);
    return group * 16 + unsigned((slots_[group] >> (4 * (index % 16))) & 15);
  }

  /** @brief Whether physical slot is live; requires slot < capacity. */
  bool occupied(std::size_t slot) const noexcept {
    assert(slot < capacity);
    return (occupied_[slot / 16] & (1u << (slot % 16))) != 0;
  }

  /** @brief Check mapping uniqueness and occupancy without allocation. */
  bool valid() const noexcept {
    if (size() > capacity) {
      return false;
    }
    std::array<std::uint16_t, 16> seen{};
    for (std::size_t i = 0; i < size(); ++i) {
      const auto slot = (*this)[i];
      const auto bit = std::uint16_t(1u << (slot % 16));
      if (seen[slot / 16] & bit) {
        return false;
      }
      seen[slot / 16] |= bit;
    }
    return seen == occupied_;
  }

  /**
   * @brief Rotate logical [left,right) left, updating order and pointer slots.
   * @details Requires left <= right <= size(). Distance is reduced modulo
   * length; empty ranges are unchanged. Inactive pointers are never read or
   * written. No pointee moves. Scratch storage is bounded by 256 pointers.
   * @return Number of writes to physical pointer slots (excludes metadata).
   */
  std::size_t rotate_children(std::array<void*, capacity>& pointers,
                              std::size_t left,
                              std::size_t right,
                              std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return 0;
    }
    if (left / 16 == (right - 1) / 16) {
      const auto group = physical_group(left / 16);
      slots_[group] =
          rotate_word(slots_[group], left % 16, (right - 1) % 16 + 1, distance);
      return 0;
    }
    if (left % 16 == 0 && right % 16 == 0) {
      const auto begin = left / 16;
      const auto end = right / 16;
      const auto remainder = distance % 16;
      const bool backwards = remainder > 8;
      groups_ = rotate_word(groups_, begin, end,
                            distance / 16 + std::size_t(backwards));
      if (remainder == 0) {
        return 0;
      }
      for (auto i = begin; i < end; ++i) {
        auto& word = slots_[physical_group(i)];
        word = rotate_word(word, 0, 16, remainder);
      }
      // Exchange the short spill between neighboring groups. The outer
      // permutation has already handled the whole-group part of the shift.
      const auto spill = backwards ? 16 - remainder : remainder;
      for (std::size_t j = 0; j < spill; ++j) {
        if (backwards) {
          auto* saved = pointers[(*this)[(end - 1) * 16 + j]];
          for (auto i = end - 1; i > begin; --i) {
            pointers[(*this)[i * 16 + j]] = pointers[(*this)[(i - 1) * 16 + j]];
          }
          pointers[(*this)[begin * 16 + j]] = saved;
        } else {
          const auto offset = 16 - spill + j;
          auto* saved = pointers[(*this)[begin * 16 + offset]];
          for (auto i = begin; i + 1 < end; ++i) {
            pointers[(*this)[i * 16 + offset]] =
                pointers[(*this)[(i + 1) * 16 + offset]];
          }
          pointers[(*this)[(end - 1) * 16 + offset]] = saved;
        }
      }
      return (end - begin) * spill;
    }
    std::array<void*, capacity> scratch;
    for (std::size_t i = 0; i < length; ++i) {
      scratch[i] = pointers[(*this)[left + i]];
    }
    auto source = distance;
    for (std::size_t i = 0; i < length; ++i) {
      pointers[(*this)[left + i]] = scratch[source];
      if (++source == length) {
        source = 0;
      }
    }
    return length;
  }

 private:
  static constexpr std::uint64_t identity = 0xfedcba9876543210ULL;
  static constexpr std::uint64_t mask(unsigned bits) noexcept {
    return bits >= 64 ? ~std::uint64_t{0} : (std::uint64_t{1} << bits) - 1;
  }
  static std::uint64_t rotate_word(std::uint64_t word,
                                   std::size_t left,
                                   std::size_t right,
                                   std::size_t distance) noexcept {
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return word;
    }
    const auto region_mask = mask(unsigned(4 * length));
    const auto shift = unsigned(4 * left);
    const auto region = (word >> shift) & region_mask;
    const auto rotated =
        (region >> (4 * distance)) | (region << (4 * (length - distance)));
    return (word & ~(region_mask << shift)) |
           ((rotated & region_mask) << shift);
  }
  unsigned physical_group(std::size_t index) const noexcept {
    return unsigned((groups_ >> (4 * index)) & 15);
  }
  std::uint64_t groups_;
  std::array<std::uint64_t, 16> slots_{};
  std::array<std::uint16_t, 16> occupied_{};
  std::uint16_t count_;
};

}  // namespace pixie::experimental
