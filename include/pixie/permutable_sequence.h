#pragma once

/**
 * @file permutable_sequence.h
 * @brief Lightweight immutable-read permutable-sequence CRTP contract.
 * @details Include pixie/permutations/sequence.h for the owning implementation.
 */

#include <pixie/sequence_options.h>

#include <concepts>
#include <cstddef>
#include <ranges>
#include <utility>

namespace pixie {
/**
 * @brief Compile-time element storage selection.
 * @details Automatic packs bool and unsigned integers of at most 64 value bits;
 * all other supported types use stable indirect ownership. Indirect may be
 * requested explicitly for otherwise packable types.
 */
enum class ElementStorage { automatic, packed, indirect };

/**
 * @brief CRTP contract for exclusively owned sequences with immutable reads.
 * @details Merge appends unchanged values, unlike permutation rebasing.
 * Implementations provide constrained from_range_impl(range), size_impl(),
 * value_at_impl(i), rotate_left_impl(left,right,distance), merge_impl(donor),
 * and memory_usage_bytes_impl() with the facade semantics documented below.
 * size_impl() and memory_usage_bytes_impl() must be noexcept. Other extension
 * points propagate the documented exceptions. Validation belongs to the
 * implementation; the facade does not repeat it. No virtual dispatch.
 * Required signatures are size_type size_impl() const noexcept,
 * ConstReference value_at_impl(size_type) const,
 * void rotate_left_impl(size_type,size_type,size_type), void merge_impl(Impl&),
 * and size_type memory_usage_bytes_impl() const noexcept. The static factory
 * returns Impl and accepts Range&& with the same input-range and iterator-move
 * constructibility constraints as from_range(). Private extension points grant
 * friendship to this base. Concrete owners default-construct canonical empty
 * and are move-only with nonthrowing move and destruction; moves leave the
 * source canonical empty and self-move is a no-op.
 * @tparam Impl Concrete owning implementation.
 * @tparam T Stored object type.
 * @tparam ConstReference Explicit result type: T or const T&, avoiding lookup
 * in the incomplete CRTP implementation and preventing reference copies.
 */
template <class Impl, class T, class ConstReference>
  requires(std::same_as<ConstReference, T> ||
           std::same_as<ConstReference, const T&>)
class PermutableSequenceBase {
 public:
  /** @brief Stored element type. @details Independent of representation. */
  using value_type = T;
  /** @brief Immutable result type. @details T by value or stable const T&. */
  using const_reference = ConstReference;
  /** @brief Count and position type. @details Full unsigned size_t range. */
  using size_type = std::size_t;
  /**
   * @brief Consume a possibly single-pass input range into an owning sequence.
   * @details from_range_impl consumes through ranges::iter_move, preserving
   * order and supporting non-common, unsized ranges. On allocation, iteration,
   * or element-construction failure all acquired ownership is reclaimed, but
   * visited inputs may already be consumed. Input references are not retained.
   * This is not the strong mutation guarantee for existing containers.
   * @tparam Range Input range whose iterator-move results construct T.
   * @param range Source consumed in iteration order.
   * @return Exclusive owner of the input values in order.
   * @throws std::length_error If count or configured chunk capacity is too
   * large.
   * @throws std::bad_alloc If storage allocation fails.
   * @throws Any Exception from input iteration or element construction.
   */
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static Impl from_range(Range&& range) {
    return Impl::from_range_impl(std::forward<Range>(range));
  }
  /**
   * @brief Return the logical element count.
   * @details size_impl() reads the count without allocation or mutation.
   * @return Count in [0,SIZE_MAX], measured in elements.
   */
  size_type size() const noexcept { return impl().size_impl(); }
  /**
   * @brief Test emptiness.
   * @details Equivalent to size()==0; canonical empty owns no heap storage.
   * @return Whether no elements remain.
   */
  bool empty() const noexcept { return size() == 0; }
  /**
   * @brief Read a checked zero-based element without exposing mutation.
   * @details value_at_impl(position) returns exactly const_reference. Indirect
   * references survive rotation, merge (including consumed donor references),
   * and ownership moves. They expire when the eventual owner destroys or
   * replaces their payload. No contiguous buffer or sentinel is promised.
   * @param position Position in [0,size()).
   * @return Immutable value or stable reference, according to storage policy.
   * @throws std::out_of_range If position >= size().
   */
  const_reference operator[](size_type position) const {
    return impl().value_at_impl(position);
  }
  /**
   * @brief Rotate [left,right) left while preserving outside values.
   * @details rotate_left_impl validates even at zero distance; empty intervals
   * are no-ops, otherwise distance is reduced modulo right-left. Allocation
   * failure leaves contents and indirect addresses unchanged.
   * @param left Inclusive start.
   * @param right Exclusive end, at most size().
   * @param distance Left rotation distance in elements.
   * @throws std::out_of_range If left > right or right > size().
   * @throws std::bad_alloc If preflight allocation fails.
   */
  void rotate_left(size_type left, size_type right, size_type distance) {
    impl().rotate_left_impl(left, right, distance);
  }
  /**
   * @brief Append unchanged donor values and consume donor ownership.
   * @details merge_impl accepts only identical configurations. Self merge and
   * empty donor are no-ops; success leaves donor canonical empty. Count or
   * allocation failure leaves both containers unchanged. Indirect objects
   * never move and references from either participant retain their addresses.
   * @param donor Owner whose unchanged values follow this sequence's values.
   * @throws std::length_error If the combined count exceeds SIZE_MAX.
   * @throws std::bad_alloc If preflight allocation fails.
   */
  void merge(Impl& donor) { impl().merge_impl(donor); }
  /**
   * @brief Explicitly account for requested live storage.
   * @details memory_usage_bytes_impl() traverses allocations and chunks in
   * linear time without allocating or throwing. Includes facade, metadata,
   * padding, and payload capacity including slack; saturates at SIZE_MAX.
   * Excludes allocations inside T, allocator overhead, RSS, stack scratch and
   * temporary peaks. Concrete memory_usage() has subset fields that must not
   * be added again to its disjoint total categories.
   * @return Requested live bytes including this owner exactly once.
   */
  size_type memory_usage_bytes() const noexcept {
    return impl().memory_usage_bytes_impl();
  }

 private:
  Impl& impl() noexcept { return static_cast<Impl&>(*this); }
  const Impl& impl() const noexcept { return static_cast<const Impl&>(*this); }
};
}  // namespace pixie
