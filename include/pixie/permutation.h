#pragma once

/**
 * @file permutation.h
 * @brief Lightweight permutation CRTP contract.
 * @details Include pixie/permutations/permutation.h for the owning
 * implementation.
 */

#include <pixie/sequence_options.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace pixie {
/**
 * @brief CRTP contract for exclusively owned permutations of [0,size()).
 * @details Construction is identity-only; reads are checked and immutable.
 * Rotation preserves the bijection and merge rebases the consumed donor.
 * Implementations provide identity_impl(n), size_impl(), value_at_impl(i),
 * insert_at_impl(position), rotate_left_impl(left,right,distance),
 * merge_impl(donor), and memory_usage_bytes_impl() with the corresponding
 * facade semantics below.
 * size_impl() and memory_usage_bytes_impl() must be noexcept; the other
 * extension points propagate the documented exceptions. Validation belongs
 * to the implementation, not a second facade check. No virtual dispatch.
 * Required signatures are static Impl identity_impl(size_type),
 * size_type size_impl() const noexcept, Index value_at_impl(size_type) const,
 * void insert_at_impl(size_type), void rotate_left_impl(size_type,size_type,
 * size_type), void merge_impl(Impl&), and size_type memory_usage_bytes_impl()
 * const noexcept. Private extension
 * points grant friendship to this base. Concrete owners default-construct
 * canonical empty and are move-only with nonthrowing move and destruction;
 * moves leave the source canonical empty and self-move is a no-op.
 * @tparam Impl Concrete owning implementation.
 * @tparam Index Unsigned integer other than bool with at most 64 value bits.
 */
template <class Impl, class Index = std::uint64_t>
class PermutationBase {
  static_assert(std::is_integral_v<Index> && std::is_unsigned_v<Index> &&
                !std::is_same_v<Index, bool> &&
                std::numeric_limits<Index>::digits <= 64);

 public:
  /** @brief Stored index type. @details Independent of the count width. */
  using value_type = Index;
  /** @brief Immutable read result. @details Returned by value, never writable.
   */
  using const_reference = Index;
  /** @brief Count and position type. @details Full unsigned size_t range. */
  using size_type = std::size_t;

  /**
   * @brief Construct an owning identity permutation.
   * @details identity_impl(n) returns Impl containing 0 through n-1, reclaiming
   * all acquired storage on failure. The domain is min(SIZE_MAX,max(Index)+1),
   * interpreted mathematically; fixed-width indices never widen.
   * @param n Number of entries, represented independently of Index.
   * @return Exclusive owner of the identity permutation.
   * @throws std::length_error If n exceeds the index domain.
   * @throws std::bad_alloc If construction allocation fails.
   */
  static Impl identity(size_type n) { return Impl::identity_impl(n); }
  /**
   * @brief Return the logical index count.
   * @details size_impl() returns the count without allocation or mutation.
   * @return Number of entries, in [0,SIZE_MAX].
   */
  size_type size() const noexcept { return impl().size_impl(); }
  /**
   * @brief Test emptiness.
   * @details Equivalent to size()==0; canonical empty owns no heap storage.
   * @return Whether no entries remain.
   */
  bool empty() const noexcept { return size() == 0; }
  /**
   * @brief Read a checked zero-based index by value.
   * @details value_at_impl(position) includes pending biases without pushing
   * tags or otherwise modifying storage. No sentinel value is used.
   * @param position Position in [0,size()).
   * @return Logical index in [0,size()).
   * @throws std::out_of_range If position >= size().
   */
  const_reference operator[](size_type position) const {
    return impl().value_at_impl(position);
  }
  /**
   * @brief Insert the next identity value at a zero-based position.
   * @details insert_at_impl(position) inserts the value size() at position,
   * shifting existing entries at or after position one to the right. The
   * result remains a valid permutation of [0,size()+1). Allocation failure
   * leaves all contents and pending biases unchanged.
   * @param position Insertion position in [0,size()].
   * @throws std::out_of_range If position > size().
   * @throws std::length_error If the new size exceeds the index domain.
   * @throws std::bad_alloc If preflight allocation fails.
   */
  void insert_at(size_type position) { impl().insert_at_impl(position); }
  /**
   * @brief Rotate [left,right) left, preserving outside indices.
   * @details rotate_left_impl validates even at zero distance. Empty intervals
   * are no-ops; otherwise distance is reduced modulo right-left. Allocation
   * failure leaves all contents and pending biases unchanged.
   * @param left Inclusive start.
   * @param right Exclusive end, at most size().
   * @param distance Left rotation distance in entries.
   * @throws std::out_of_range If left > right or right > size().
   * @throws std::bad_alloc If preflight allocation fails.
   */
  void rotate_left(size_type left, size_type right, size_type distance) {
    impl().rotate_left_impl(left, right, distance);
  }
  /**
   * @brief Append donor indices plus the old receiver size and consume donor.
   * @details merge_impl accepts only the same concrete configuration. Self
   * merge and empty donor are no-ops; success leaves donor canonical empty.
   * Size checks and allocating preflight precede rebasing or ownership edits;
   * either failure leaves both participants unchanged.
   * @param donor Owner of the permutation to append and rebase.
   * @throws std::length_error If the combined count or index domain overflows.
   * @throws std::bad_alloc If preflight allocation fails.
   */
  void merge(Impl& donor) { impl().merge_impl(donor); }
  /**
   * @brief Explicitly account for requested live storage.
   * @details memory_usage_bytes_impl() traverses allocations in linear time,
   * without allocating or throwing. Includes the facade, alignment, metadata,
   * and unused payload capacity; saturates at SIZE_MAX. Excludes allocator
   * overhead, RSS, stack scratch, and construction/preflight peaks. Concrete
   * memory_usage() supplies a structured breakdown, not a second total.
   * @return Total requested live bytes including this owner exactly once.
   */
  size_type memory_usage_bytes() const noexcept {
    return impl().memory_usage_bytes_impl();
  }

 private:
  Impl& impl() noexcept { return static_cast<Impl&>(*this); }
  const Impl& impl() const noexcept { return static_cast<const Impl&>(*this); }
};
}  // namespace pixie
