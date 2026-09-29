#pragma once

/**
 * @file sorted_vector.h
 * @brief Sorted logical vector backed by append-only values and a permutation.
 */

#include <pixie/permutations/permutation.h>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Duplicate-preserving sorted view over owned vector payload and order.
 * @details Values append physically to a std::vector while Permutation stores
 * their sorted logical order. Insertion moves only order indices after the
 * vector append, so it does not shift existing payloads. Reads return const T&
 * into the vector; every vector reallocation invalidates those references.
 * There is no contiguous storage promise in sorted order and no iterator
 * interface. Bounds use predicate-guided search through the permutation tree.
 * @tparam T Non-bool, unqualified stored object type.
 * @tparam Compare Strict weak ordering for T and lookup keys.
 * @tparam Index Unsigned permutation index type with at most 64 value bits.
 * @tparam StorageBits Local permutation block budget.
 * @tparam Fanout Even permutation-tree fanout of at least four.
 * @tparam Layout Permutation-tree child-length representation.
 */
template <class T = std::uint64_t,
          class Compare = std::less<>,
          class Index = std::uint64_t,
          std::size_t StorageBits = 2048,
          std::size_t Fanout = 8,
          LengthLayout Layout = LengthLayout::cumulative>
class SortedPermutationVector {
  static_assert(std::is_object_v<T> && !std::is_array_v<T> &&
                    std::same_as<T, std::remove_cv_t<T>> &&
                    !std::same_as<T, bool>,
                "SortedPermutationVector requires a non-bool object type");

 public:
  /** @brief Underlying order-index type. */
  using permutation_type = Permutation<Index, StorageBits, Fanout, Layout>;
  /** @brief Stored value type. */
  using value_type = T;
  /** @brief Immutable read result borrowing the owned vector. */
  using const_reference = const T&;
  /** @brief Count and position type. */
  using size_type = std::size_t;

  /**
   * @brief Construct empty with a comparator.
   * @details Performs no value or permutation allocation.
   * @param compare Ordering relation for values and lookup keys.
   */
  explicit SortedPermutationVector(Compare compare = {})
      : compare_(std::move(compare)) {}
  /** @brief Exclusive ownership disallows copying. */
  SortedPermutationVector(const SortedPermutationVector&) = delete;
  /** @brief Exclusive ownership disallows copy assignment. */
  SortedPermutationVector& operator=(const SortedPermutationVector&) = delete;
  /** @brief Move vector payload, permutation, and comparator state. */
  SortedPermutationVector(SortedPermutationVector&&) = default;
  /** @brief Replace vector payload, permutation, and comparator by move. */
  SortedPermutationVector& operator=(SortedPermutationVector&&) = default;

  /**
   * @brief Build a sorted owner by consuming values one at a time.
   * @details Equal incoming values are inserted before existing equivalent
   * values. Input values append physically in input order; only the permutation
   * describes sorted order.
   * @tparam Range Input range whose iterator-move result constructs T.
   * @param range Source consumed in iteration order.
   * @param compare Ordering relation retained by the result.
   * @return Owner containing every input value in nondecreasing order.
   * @throws Any Exception from input iteration, value construction, comparison,
   * vector growth, or permutation operations.
   */
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static SortedPermutationVector from_range(Range&& range,
                                                Compare compare = {}) {
    SortedPermutationVector result(std::move(compare));
    auto it = std::ranges::begin(range);
    const auto end = std::ranges::end(range);
    for (; it != end; ++it) {
      result.insert(T(std::ranges::iter_move(it)));
    }
    return result;
  }

  /** @brief Return the number of sorted values. */
  size_type size() const noexcept { return permutation_.size(); }
  /** @brief Return whether this owner has no values. */
  bool empty() const noexcept { return permutation_.empty(); }
  /**
   * @brief Read one sorted zero-based value.
   * @details The result references the append-only vector payload. It remains
   * valid until a vector reallocation or this wrapper's destruction.
   * @param position Position in [0,size()).
   * @return Const reference to the value at the sorted position.
   * @throws std::out_of_range If position >= size().
   */
  const_reference operator[](size_type position) const {
    return values_[permutation_[position]];
  }

  /**
   * @brief Reserve append-only payload capacity.
   * @details Does not alter sorted order. A vector reallocation invalidates all
   * references previously returned by operator[].
   * @param count Minimum physical value capacity.
   * @throws std::length_error If count exceeds vector capacity limits.
   * @throws std::bad_alloc If vector allocation fails.
   */
  void reserve(size_type count) { values_.reserve(count); }

  /**
   * @brief Return the first position whose value is not ordered before key.
   * @details Uses Compare(value,key) in the underlying permutation search. The
   * result is in [0,size()]; size() means every stored value is before key.
   * @tparam Key Lookup type accepted by Compare.
   * @param key Lookup value; it is not stored or retained.
   * @return Lower-bound insertion position.
   */
  template <class Key>
  size_type lower_bound_index(const Key& key) const {
    return permutation_.lower_bound(
        [&](Index idx) { return !compare_(values_[idx], key); });
  }

  /**
   * @brief Return the first position whose value is ordered after key.
   * @details Uses Compare(key,value) in the underlying permutation search. The
   * result is in [0,size()]; size() means no stored value is after key.
   * @tparam Key Lookup type accepted by Compare.
   * @param key Lookup value; it is not stored or retained.
   * @return Upper-bound position.
   */
  template <class Key>
  size_type upper_bound_index(const Key& key) const {
    return permutation_.lower_bound(
        [&](Index idx) { return compare_(key, values_[idx]); });
  }

  /**
   * @brief Insert a value before existing equivalent values.
   * @details The payload first appends at its next physical vector index. The
   * permutation's atomic insert_at then places the new index at the
   * lower-bound position. Existing payload objects are never shifted by
   * sorted-order maintenance. If permutation insertion fails, the appended
   * payload is removed; vector capacity and reference invalidation from its
   * preceding growth are not rolled back.
   * @param value Value to own and insert.
   * @return Sorted position occupied by the new value.
   * @throws std::length_error If the next physical index exceeds Index.
   * @throws Any Exception from comparison, value construction, vector growth,
   * or permutation insertion.
   */
  size_type insert(T value) {
    const size_type position = lower_bound_index(value);
    if (values_.size() > std::numeric_limits<Index>::max()) {
      throw std::length_error("SortedPermutationVector: index domain");
    }
    values_.push_back(std::move(value));
    try {
      permutation_.insert_at(position);
    } catch (...) {
      values_.pop_back();
      throw;
    }
    return position;
  }

  /**
   * @brief Return explicit requested bytes owned by this wrapper.
   * @details Includes wrapper/comparator bytes, vector payload capacity, and
   * permutation allocations. Excludes allocator bookkeeping, RSS, and memory
   * recursively owned by T. Arithmetic saturates at SIZE_MAX.
   * @return Total requested live bytes.
   */
  size_type memory_usage_bytes() const noexcept {
    const size_type maximum = std::numeric_limits<size_type>::max();
    const size_type payload = values_.capacity() > maximum / sizeof(T)
                                  ? maximum
                                  : values_.capacity() * sizeof(T);
    const size_type permutation =
        permutation_.memory_usage_bytes() - sizeof(permutation_);
    const size_type base =
        payload > maximum - sizeof(*this) ? maximum : sizeof(*this) + payload;
    return permutation > maximum - base ? maximum : base + permutation;
  }

 private:
  std::vector<T> values_;
  permutation_type permutation_;
  [[no_unique_address]] Compare compare_;
};

}  // namespace pixie
