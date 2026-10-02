#pragma once

/**
 * @file sorted_sequence.h
 * @brief Sorted owning wrapper over a permutable sequence.
 */

#include <pixie/permutations/sequence.h>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <ranges>
#include <utility>

namespace pixie {

/**
 * @brief Duplicate-preserving sorted sequence backed by PermutableSequence.
 * @details Values are in comparator order under Compare. insert() finds the
 * first equivalent-or-greater position and atomically inserts the owned value
 * there. Reads retain the
 * underlying storage policy: packed values return by value and indirect values
 * return stable const references. There is no contiguous storage or iterator
 * interface. Bounds use the underlying tree's predicate-guided search.
 * @tparam T Stored object type accepted by PermutableSequence.
 * @tparam Compare Strict weak ordering for T and lookup keys.
 * @tparam Storage Compile-time payload storage policy.
 * @tparam StorageBits Local sequence block budget.
 * @tparam Fanout Even sequence-tree fanout of at least four.
 * @tparam Layout Sequence-tree child-length representation.
 * @tparam ChunkBytes Indirect payload target bytes per chunk.
 */
template <class T = std::uint64_t,
          class Compare = std::less<>,
          ElementStorage Storage = ElementStorage::automatic,
          std::size_t StorageBits = 2048,
          std::size_t Fanout = 8,
          LengthLayout Layout = LengthLayout::cumulative,
          std::size_t ChunkBytes = 4096>
class SortedPermutableSequence {
 public:
  /** @brief Underlying owning sequence type. */
  using sequence_type =
      PermutableSequence<T, Storage, StorageBits, Fanout, Layout, ChunkBytes>;
  /** @brief Stored value type. */
  using value_type = T;
  /** @brief Immutable indexed result inherited from the sequence policy. */
  using const_reference = typename sequence_type::const_reference;
  /** @brief Count and position type. */
  using size_type = std::size_t;

  /**
   * @brief Construct empty with a comparator.
   * @details Performs no sequence allocation. The comparator is used for every
   * bound and insertion for this owner's lifetime.
   * @param compare Ordering relation for values and lookup keys.
   */
  explicit SortedPermutableSequence(Compare compare = {})
      : compare_(std::move(compare)) {}
  /** @brief Exclusive ownership disallows copying. */
  SortedPermutableSequence(const SortedPermutableSequence&) = delete;
  /** @brief Exclusive ownership disallows copy assignment. */
  SortedPermutableSequence& operator=(const SortedPermutableSequence&) = delete;
  /** @brief Move ownership and comparator state. */
  SortedPermutableSequence(SortedPermutableSequence&&) = default;
  /** @brief Replace ownership and comparator state by move. */
  SortedPermutableSequence& operator=(SortedPermutableSequence&&) = default;

  /**
   * @brief Build a sorted owner by consuming values one at a time.
   * @details This preserves input-range support and the underlying payload
   * requirements without materializing a sortable temporary vector. Equal
   * incoming values are inserted before existing equivalent values.
   * @tparam Range Input range whose iterator-move result constructs T.
   * @param range Source consumed in iteration order.
   * @param compare Ordering relation retained by the result.
   * @return Owner containing every input value in nondecreasing order.
   * @throws Any Exception from input iteration, value construction, comparison,
   * or the underlying sequence operations.
   */
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static SortedPermutableSequence from_range(Range&& range,
                                                 Compare compare = {}) {
    SortedPermutableSequence result(std::move(compare));
    auto it = std::ranges::begin(range);
    const auto end = std::ranges::end(range);
    for (; it != end; ++it) {
      result.insert(T(std::ranges::iter_move(it)));
    }
    return result;
  }

  /** @brief Return the number of sorted values. */
  size_type size() const noexcept { return sequence_.size(); }
  /** @brief Return whether this owner is empty. */
  bool empty() const noexcept { return sequence_.empty(); }
  /**
   * @brief Read one sorted zero-based value.
   * @details Packed storage returns T by value. Indirect storage returns const
   * T& whose lifetime follows the underlying sequence's documented rules.
   * @param position Position in [0,size()).
   * @return Value at the sorted position.
   * @throws std::out_of_range If position >= size().
   */
  const_reference operator[](size_type position) const {
    return sequence_[position];
  }

  /**
   * @brief Return the first position whose value is not ordered before key.
   * @details Uses Compare(value,key) in the underlying tree search. The
   * result is in [0,size()]; size() means every stored value is before key.
   * @tparam Key Lookup type accepted by Compare.
   * @param key Lookup value; it is not stored or retained.
   * @return Lower-bound insertion position.
   */
  template <class Key>
  size_type lower_bound_index(const Key& key) const {
    return sequence_.lower_bound([&](const T& v) { return !compare_(v, key); });
  }

  /**
   * @brief Return the first position whose value is ordered after key.
   * @details Uses Compare(key,value) in the underlying tree search. The
   * result is in [0,size()]; size() means no stored value is after key.
   * @tparam Key Lookup type accepted by Compare.
   * @param key Lookup value; it is not stored or retained.
   * @return Upper-bound position.
   */
  template <class Key>
  size_type upper_bound_index(const Key& key) const {
    return sequence_.lower_bound([&](const T& v) { return compare_(key, v); });
  }

  /**
   * @brief Insert a value before existing equivalent values.
   * @details Computes lower_bound_index(value) and delegates to the
   * underlying sequence's atomic insert_at. The resulting sequence remains
   * sorted on successful return.
   * @param value Value to own and insert.
   * @return Sorted position occupied by the new value.
   * @throws Any Exception from comparison, value construction, allocation, or
   * the underlying insert operation.
   */
  size_type insert(T value) {
    const size_type position = lower_bound_index(value);
    sequence_.insert_at(position, std::move(value));
    return position;
  }

  /**
   * @brief Return explicit requested bytes owned by this wrapper.
   * @details Includes wrapper and comparator storage, plus all underlying
   * sequence storage. Excludes allocator bookkeeping, RSS, and allocations
   * owned by stored values. Arithmetic saturates at SIZE_MAX.
   * @return Total requested live bytes.
   */
  size_type memory_usage_bytes() const noexcept {
    const size_type nested = sequence_.memory_usage_bytes() - sizeof(sequence_);
    return nested > std::numeric_limits<size_type>::max() - sizeof(*this)
               ? std::numeric_limits<size_type>::max()
               : sizeof(*this) + nested;
  }

 private:
  sequence_type sequence_;
  [[no_unique_address]] Compare compare_;
};

}  // namespace pixie
