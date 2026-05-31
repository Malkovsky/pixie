#pragma once

#include <cstddef>
#include <limits>

namespace pixie::rmq {

/**
 * @brief CRTP facade for static range-minimum-query indexes.
 *
 * Implementations are non-owning indexes over an external random-access array.
 * Queries use inclusive zero-based ranges and return the first position
 * attaining the minimum. Invalid ranges return `npos`.
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
   * @brief Return the first minimum position in [@p left, @p right].
   *
   * @details The query range is inclusive. Ties are resolved by returning the
   * smallest position attaining the minimum. Invalid ranges return `npos`.
   *
   * @param left First position in the query range.
   * @param right Last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    return impl().arg_min_impl(left, right);
  }

  /**
   * @brief Return the minimum value in [@p left, @p right].
   *
   * @details Invalid ranges return a default-constructed value.
   *
   * @param left First position in the query range.
   * @param right Last position in the query range.
   * @return The minimum value in the inclusive range, or `Value{}`.
   */
  Value range_min(std::size_t left, std::size_t right) const {
    const std::size_t position = arg_min(left, right);
    if (position == npos) {
      return Value{};
    }
    return impl().value_at_impl(position);
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
