#pragma once

#include <pixie/storage/aligned.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

/**
 * @brief Bounded, nonowning block of native pointers for a SequenceTree.
 * @details StorageBits budgets the complete inline object, including length,
 * circular origin and alignment padding. Pointers are never encoded as integer
 * fields. Local mutations move pointers only, never pointed-to objects, and use
 * bounded stack scratch without allocation. The caller owns all pointees.
 * @tparam T Pointee type; reads expose const T* without dereferencing it.
 * @tparam StorageBits Complete inline budget in bits, a positive
 * CacheLine-sized multiple large enough for both metadata fields and at least
 * one pointer.
 */
template <class T, std::size_t StorageBits = 2048>
class PointerBlock {
  static_assert(StorageBits % (sizeof(CacheLine) * CHAR_BIT) == 0);
  static_assert(StorageBits / CHAR_BIT >=
                2 * sizeof(std::size_t) + sizeof(const T*));

 public:
  /**
   * @brief Immutable-pointee native pointer returned by value.
   * @details Does not own or extend the lifetime of the pointed-to object.
   */
  using value_type = const T*;
  /**
   * @brief Actual slot count after paying for both local metadata fields.
   * @details Includes only native pointer slots; the complete aligned storage
   * object is statically checked against StorageBits.
   */
  static constexpr std::size_t capacity =
      (StorageBits / CHAR_BIT - 2 * sizeof(std::size_t)) / sizeof(value_type);
  static_assert(capacity <= std::numeric_limits<std::size_t>::max() / 2);

 private:
  struct alignas(CacheLine) Storage {
    std::size_t length = 0;
    std::size_t origin = 0;
    std::array<value_type, capacity> values{};
  };
  static_assert(sizeof(Storage) == StorageBits / CHAR_BIT);
  Storage storage_;

  value_type& slot(std::size_t i) noexcept {
    return storage_.values[(storage_.origin + i) % storage_.length];
  }

 public:
  /**
   * @brief Construct an empty block.
   * @details Initializes length and origin to zero without allocating.
   */
  PointerBlock() noexcept = default;
  /**
   * @brief Copy pointer values from input, not the pointees.
   * @details No input ownership is acquired, and no allocation is performed.
   * Copies pointer values in their input order with a zero circular origin.
   * @param input Pointer values to copy; the span need not outlive
   * construction.
   * @throws std::length_error If input.size() > capacity.
   */
  explicit PointerBlock(std::span<const value_type> input) {
    if (input.size() > capacity) {
      throw std::length_error("PointerBlock: capacity");
    }
    std::copy(input.begin(), input.end(), storage_.values.begin());
    storage_.length = input.size();
  }
  /**
   * @brief Return the number of logical pointers.
   * @details Constant time and nonthrowing; excludes unused slots.
   * @return Pointer count in [0,capacity].
   */
  std::size_t size() const noexcept { return storage_.length; }
  /**
   * @brief Read a pointer by value at a zero-based position.
   * @details Never dereferences the pointer. This unchecked local operation
   * asserts its precondition in debug builds; it does not throw.
   * @param i Zero-based logical position in [0,size()).
   * @return Native pointer value at position i, without acquiring ownership.
   * @pre i < size().
   */
  value_type operator[](std::size_t i) const noexcept {
    assert(i < size());
    return storage_.values[(storage_.origin + i) % size()];
  }
  /**
   * @brief Replace a logical pointer without touching either pointee.
   * @pre i < size(). No ownership is acquired or released.
   * @param i Zero-based position. @param value Replacement pointer.
   */
  void set_at(std::size_t i, value_type value) noexcept {
    assert(i < size());
    slot(i) = value;
  }
  /**
   * @brief Rotate the valid half-open range [left,right) without allocation.
   * @details Distance is reduced modulo a nonempty range's length. Empty valid
   * ranges are no-ops. Whole-block rotations update only the origin. Other
   * elements retain their order, and pointees are never moved or dereferenced.
   * @param left Inclusive start of the half-open range.
   * @param right Exclusive end of the half-open range.
   * @param distance Left-rotation distance in pointer elements.
   * @pre left <= right && right <= size(); debug builds assert this condition.
   */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    const auto n = right - left;
    if (n == 0 || (distance %= n) == 0) {
      return;
    }
    if (left == 0 && right == size()) {
      storage_.origin = (storage_.origin + distance) % size();
      return;
    }
    auto reverse = [&](std::size_t begin, std::size_t end) {
      while (begin < end && begin < --end) {
        std::swap(slot(begin++), slot(end));
      }
    };
    reverse(left, left + distance);
    reverse(left + distance, right);
    reverse(left, right);
  }
  /**
   * @brief Insert one pointer at a local offset without allocation.
   * @details Requires offset <= size() < capacity. Shifts [offset,size()) one
   * position right using bounded stack scratch and places value at offset.
   * Resets both circular origins. Does not allocate or throw.
   * @param offset Local insertion position in [0,size()].
   * @param value Pointer to insert.
   * @pre offset <= size() && size() < capacity; debug builds assert.
   */
  void insert_at(std::size_t offset, value_type value) noexcept {
    assert(offset <= size() && size() < capacity);
    std::array<value_type, capacity> scratch;
    for (std::size_t i = 0; i < offset; ++i) {
      scratch[i] = (*this)[i];
    }
    scratch[offset] = value;
    const auto n = size();
    for (std::size_t i = offset; i < n; ++i) {
      scratch[i + 1] = (*this)[i];
    }
    std::copy_n(scratch.begin(), n + 1, storage_.values.begin());
    storage_.length = n + 1;
    storage_.origin = 0;
  }
  /**
   * @brief Redistribute concatenated pointers into two distinct blocks.
   * @details Preserves concatenated order, leaving left_count pointers in this
   * block and the remaining pointers in right. Resets both circular origins;
   * uses at most 2*capacity scratch pointers and cannot allocate or throw.
   * @param right Distinct block holding the suffix of the input concatenation.
   * @param left_count Desired pointer count in this block after redistribution.
   * @pre This block and right are distinct; left_count <= capacity;
   * left_count <= size()+right.size(); and size()+right.size()-left_count <=
   * capacity. Debug builds assert these conditions.
   */
  void redistribute(PointerBlock& right, std::size_t left_count) noexcept {
    const auto total = size() + right.size();
    assert(this != &right && left_count <= capacity && left_count <= total &&
           total - left_count <= capacity);
    std::array<value_type, 2 * capacity> scratch;
    for (std::size_t i = 0; i < size(); ++i) {
      scratch[i] = (*this)[i];
    }
    for (std::size_t i = 0; i < right.size(); ++i) {
      scratch[size() + i] = right[i];
    }
    std::copy_n(scratch.begin(), left_count, storage_.values.begin());
    std::copy_n(scratch.begin() + left_count, total - left_count,
                right.storage_.values.begin());
    storage_.length = left_count;
    right.storage_.length = total - left_count;
    storage_.origin = right.storage_.origin = 0;
  }
};

}  // namespace pixie::detail::sequence
/// @endcond
