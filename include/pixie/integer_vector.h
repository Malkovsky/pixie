#pragma once

/**
 * @file integer_vector.h
 * @brief Common interfaces for immutable positional integer vectors.
 *
 * Include a concrete header under `<pixie/integer_vector/>` to use an
 * integer-vector implementation.
 */

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace pixie {

/** @brief An unsigned integer-vector value type no wider than 64 bits. */
template <class Value>
concept IntegerVectorValue = std::unsigned_integral<Value> &&
                             std::same_as<Value, std::remove_cv_t<Value>> &&
                             !std::same_as<std::remove_cv_t<Value>, bool> &&
                             (std::numeric_limits<Value>::digits <= 64);

/**
 * @brief CRTP facade for an immutable positional integer vector.
 *
 * @tparam Impl Concrete implementation providing `size_impl()`,
 * `value_at_impl()`, and `copy_to_impl()`. It may also provide
 * `memory_usage_bytes_impl()`.
 * @tparam Value Unsigned element type other than `bool`, no wider than 64 bits.
 *
 * @details Positions are zero-based. Implementations retain ownership and
 * lifetime semantics documented by their concrete types.
 */
template <class Impl, IntegerVectorValue Value = std::uint64_t>
class IntegerVectorBase {
 public:
  /** @brief Stored integer type. */
  using value_type = Value;

  /** @brief Type used for element counts and zero-based positions. */
  using size_type = std::size_t;

  /** @brief Return the logical number of elements. */
  size_type size() const { return impl().size_impl(); }

  /** @brief Return whether the vector has no elements. */
  bool empty() const { return size() == 0; }

  /**
   * @brief Read the element at zero-based @p position without bounds checking.
   * @param position Position in `[0, size())`.
   * @return The stored value.
   * @pre `position < size()`.
   */
  value_type operator[](size_type position) const {
    return impl().value_at_impl(position);
  }

  /**
   * @brief Read the element at zero-based @p position with bounds checking.
   * @throws std::out_of_range if `position >= size()`.
   */
  value_type at(size_type position) const {
    if (position >= size()) {
      throw std::out_of_range("Integer-vector position is out of range");
    }
    return (*this)[position];
  }

  /**
   * @brief Copy a checked contiguous source range into @p output.
   *
   * @details Copies `[begin, begin + output.size())`. The complete source
   * range is validated before the implementation is called, so an invalid
   * range leaves @p output unchanged.
   *
   * @throws std::out_of_range if the requested source range is invalid.
   */
  void copy_to(size_type begin, std::span<value_type> output) const {
    if (begin > size() || output.size() > size() - begin) {
      throw std::out_of_range("Integer-vector copy range is out of bounds");
    }
    impl().copy_to_impl(begin, output);
  }

  /**
   * @brief Return total memory owned by this vector when supported.
   * @return Inline object bytes plus storage owned below the object. Borrowed
   * backing bytes are excluded.
   */
  size_type memory_usage_bytes() const
    requires requires(const Impl& concrete) {
      { concrete.memory_usage_bytes_impl() } -> std::convertible_to<size_type>;
    }
  {
    return impl().memory_usage_bytes_impl();
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

/**
 * @brief CRTP refinement for duplicate-preserving monotone integer vectors.
 *
 * @tparam Impl Concrete implementation providing the integer-vector extension
 * points plus `lower_bound_index_impl()` and `upper_bound_index_impl()`.
 * @tparam Value Unsigned element type other than `bool`, no wider than 64 bits.
 *
 * @details Monotone means nondecreasing; equal adjacent values are retained.
 * This contract does not prescribe representation-independent query
 * complexity.
 */
template <class Impl, IntegerVectorValue Value = std::uint64_t>
class MonotoneIntegerVectorBase : public IntegerVectorBase<Impl, Value> {
 public:
  using typename IntegerVectorBase<Impl, Value>::size_type;
  using typename IntegerVectorBase<Impl, Value>::value_type;

  /**
   * @brief Return the first index whose value is greater than or equal to @p x.
   * @return An index in `[0, size()]`; `size()` means no such element exists.
   */
  size_type lower_bound_index(value_type x) const {
    return impl().lower_bound_index_impl(x);
  }

  /**
   * @brief Return the first index whose value is greater than @p x.
   * @return An index in `[0, size()]`; `size()` means no such element exists.
   */
  size_type upper_bound_index(value_type x) const {
    return impl().upper_bound_index_impl(x);
  }

  /** @brief Return whether at least one element equals @p x. */
  bool contains(value_type x) const {
    const size_type position = lower_bound_index(x);
    return position != this->size() && (*this)[position] == x;
  }

  /** @brief Return the number of elements equal to @p x. */
  size_type count(value_type x) const {
    return upper_bound_index(x) - lower_bound_index(x);
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
