#pragma once

/**
 * @file rank_select.h
 * @brief Common interface for rank/select support over packed bit sequences.
 *
 * Include a concrete implementation header such as
 * `<pixie/rank_select/compressed.h>` or `<pixie/rank_select/support.h>`.
 * `<pixie/rank_select/implementations.h>` is the benchmark/reference catalog.
 */

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>

namespace pixie {

/**
 * @brief CRTP facade for immutable rank/select queries and bit access.
 * @tparam Impl Concrete type implementing the const `*_impl()` operations
 * documented on the facade methods. No virtual dispatch is used.
 * @details Source representation, ownership, and lifetime are defined by the
 * concrete implementation. Borrowed data must remain alive and immutable while
 * queries use it. Rank counts half-open prefixes; select takes a one-based
 * occurrence rank and returns a zero-based bit position.
 *
 * @see `<pixie/rank_select/implementations.h>` for the available concrete
 * implementations.
 */
template <class Impl>
class RankSelectBase {
 public:
  /**
   * @brief Return the number of valid bits.
   * @return Logical bit-vector length.
   * @details Requires `Impl::size_impl()` to return the logical length.
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Check whether the bit vector is empty.
   * @return `true` when `size() == 0`.
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Read a bit by zero-based position.
   * @param position Position in `[0, size())`.
   * @return Zero or one.
   * @details Requires `Impl::bit_impl(position)` to return the logical bit.
   * Positions at or beyond `size()` violate the precondition; no bounds check
   * is required.
   */
  int operator[](std::size_t position) const {
    return impl().bit_impl(position);
  }

  /**
   * @brief Count one bits in the prefix `[0, end_position)`.
   * @param end_position Prefix boundary, clamped to `size()`.
   * @return Number of one bits in the prefix.
   * @details Requires `Impl::rank_impl(end_position)` to count ones in
   * `[0, min(end_position, size()))`, including zero for an empty prefix.
   */
  std::uint64_t rank(std::size_t end_position) const {
    return impl().rank_impl(end_position);
  }

  /**
   * @brief Count zero bits in the prefix `[0, end_position)`.
   * @param end_position Prefix boundary, clamped to `size()`.
   * @return Number of zero bits in the prefix.
   * @details An optional `Impl::rank0_impl(end_position)` may compute this
   * prefix directly. It must clamp boundaries at or beyond `size()` to the
   * total zero count. Otherwise the facade derives zeros from one rank.
   */
  std::uint64_t rank0(std::size_t end_position) const {
    if constexpr (requires { impl().rank0_impl(end_position); }) {
      return impl().rank0_impl(end_position);
    } else {
      return end_position >= size() ? size() - rank(size())
                                    : end_position - rank(end_position);
    }
  }

  /**
   * @brief Return the position of the rank-th one bit.
   * @param rank One-based occurrence rank. Rank zero returns zero.
   * @return Zero-based bit position, or `size()` if unavailable or absent.
   * @details Requires `Impl::select_impl(rank)` to return zero for rank zero,
   * `size()` when one-select is disabled or rank exceeds the one population,
   * and otherwise the unique position of that occurrence.
   */
  std::uint64_t select(std::size_t rank) const {
    return impl().select_impl(rank);
  }

  /**
   * @brief Return the position of the rank-th zero bit.
   * @param rank One-based occurrence rank. Rank zero returns zero.
   * @return Zero-based bit position, or `size()` if unavailable or absent.
   * @details Requires `Impl::select0_impl(rank)` to return zero for rank zero,
   * `size()` when zero-select is disabled or rank exceeds the zero population,
   * and otherwise the unique position of that occurrence.
   */
  std::uint64_t select0(std::size_t rank) const {
    return impl().select0_impl(rank);
  }

  /**
   * @brief Check whether select queries for one bits are enabled.
   * @return `true` when `select()` is supported by this instance.
   * @details Requires `Impl::supports_select1_impl()` to report whether this
   * instance enables one-select, independently of its one population.
   */
  bool supports_select1() const { return impl().supports_select1_impl(); }

  /**
   * @brief Check whether select queries for zero bits are enabled.
   * @return `true` when `select0()` is supported by this instance.
   * @details Requires `Impl::supports_select0_impl()` to report whether this
   * instance enables zero-select, independently of its zero population.
   */
  bool supports_select0() const { return impl().supports_select0_impl(); }

  /**
   * @brief Return owned index memory usage when provided by the implementation.
   * @return Size of the object and its owned storage in bytes.
   * @details The optional `Impl::memory_usage_bytes_impl()` returns a value
   * convertible to `std::size_t`. Include owned source copies when present;
   * exclude borrowed storage. This reports allocation size, not resident memory
   * or allocator bookkeeping.
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

  /**
   * @brief Convert the logical bits to a binary string.
   * @return A string of `size()` characters containing only `0` and `1`.
   */
  std::string to_string() const {
    std::string result;
    result.reserve(size());
    for (std::size_t i = 0; i < size(); ++i) {
      result.push_back((*this)[i] ? '1' : '0');
    }
    return result;
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
