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
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Whether the indexed array is empty.
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Return the first minimum position in [@p left, @p right].
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    return impl().arg_min_impl(left, right);
  }

  /**
   * @brief Return the minimum value in [@p left, @p right].
   * @details Invalid ranges return a default-constructed value.
   */
  Value range_min(std::size_t left, std::size_t right) const {
    const std::size_t position = arg_min(left, right);
    if (position == npos) {
      return Value{};
    }
    return impl().value_at_impl(position);
  }

 private:
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie::rmq
