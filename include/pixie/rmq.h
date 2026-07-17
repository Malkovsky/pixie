#pragma once

/**
 * @file rmq.h
 * @brief Common interface for static range-minimum-query indexes.
 *
 * Include `<pixie/rmq/implementations.h>` for the implementation catalog and
 * current benchmark comparison.
 */

#include <concepts>
#include <cstddef>
#include <limits>

namespace pixie::rmq {

/**
 * @brief CRTP facade for static range-minimum-query indexes.
 *
 * Implementations are non-owning indexes over an external random-access array.
 * Queries use half-open zero-based ranges `[left, right)` and return the first
 * position attaining the minimum. Invalid or empty ranges return `npos`.
 *
 * @see `<pixie/rmq/implementations.h>` for the complete implementation catalog
 * and current benchmark comparison.
 */
template <class Impl, class Value>
class RmqBase {
 public:
  /**
   * @brief Sentinel returned when no valid query answer exists.
   */
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

  /**
   * @brief Number of indexed values.
   *
   * @return The number of values covered by the underlying RMQ index.
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Whether the indexed array is empty.
   *
   * @return `true` when `size() == 0`.
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Return the first minimum position in [@p left, @p right).
   *
   * @details The query range is half-open. Ties are resolved by returning the
   * smallest position attaining the minimum. Invalid or empty ranges return
   * `npos`.
   *
   * @param left First position in the query range.
   * @param right One past the last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    return impl().arg_min_impl(left, right);
  }

  /**
   * @brief Return the minimum value in [@p left, @p right).
   *
   * @details Invalid or empty ranges return a default-constructed value.
   *
   * @param left First position in the query range.
   * @param right One past the last position in the query range.
   * @return The minimum value in the half-open range, or `Value{}`.
   */
  Value range_min(std::size_t left, std::size_t right) const {
    const std::size_t position = arg_min(left, right);
    if (position == npos) {
      return Value{};
    }
    return impl().value_at_impl(position);
  }

  /**
   * @brief Return owned auxiliary memory usage in bytes when implemented.
   *
   * @details Implementations opt in by exposing
   * `memory_usage_bytes_impl() const`. The reported value should include the
   * index object itself and heap buffers owned by the index, but not the
   * external value array indexed by the RMQ.
   */
  std::size_t memory_usage_bytes() const
    requires requires(const Impl& concrete) {
      {
        concrete.memory_usage_bytes_impl()
      } -> std::convertible_to<std::size_t>;
    }
  {
    return impl().memory_usage_bytes_impl();
  }

 private:
  /**
   * @brief Return this object as its concrete CRTP implementation.
   *
   * @return Reference to the derived RMQ implementation.
   */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie::rmq
