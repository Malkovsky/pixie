#pragma once

/**
 * @file rank_select.h
 * @brief Common interface for rank/select support over packed bit sequences.
 *
 * Include `<pixie/rank_select/implementations.h>` to use Pixie's concrete
 * rank/select implementations.
 */

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>

namespace pixie {

/**
 * @brief CRTP facade for immutable rank/select queries and bit access.
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
   */
  int operator[](std::size_t position) const {
    return impl().bit_impl(position);
  }

  /**
   * @brief Count one bits in the prefix `[0, end_position)`.
   * @param end_position Prefix boundary in `[0, size()]`.
   * @return Number of one bits in the prefix.
   */
  std::uint64_t rank(std::size_t end_position) const {
    return impl().rank_impl(end_position);
  }

  /**
   * @brief Count zero bits in the prefix `[0, end_position)`.
   * @param end_position Prefix boundary in `[0, size()]`.
   * @return Number of zero bits in the prefix.
   */
  std::uint64_t rank0(std::size_t end_position) const {
    return end_position >= size() ? size() - rank(size())
                                  : end_position - rank(end_position);
  }

  /**
   * @brief Return the position of the rank-th one bit.
   * @param rank One-based rank. Rank zero returns zero for compatibility.
   * @return Zero-based bit position, or `size()` if unavailable or absent.
   */
  std::uint64_t select(std::size_t rank) const {
    return impl().select_impl(rank);
  }

  /**
   * @brief Return the position of the rank-th zero bit.
   * @param rank One-based rank. Rank zero returns zero for compatibility.
   * @return Zero-based bit position, or `size()` if unavailable or absent.
   */
  std::uint64_t select0(std::size_t rank) const {
    return impl().select0_impl(rank);
  }

  /**
   * @brief Check whether select queries for one bits are enabled.
   * @return `true` when `select()` is supported by this instance.
   */
  bool supports_select1() const { return impl().supports_select1_impl(); }

  /**
   * @brief Check whether select queries for zero bits are enabled.
   * @return `true` when `select0()` is supported by this instance.
   */
  bool supports_select0() const { return impl().supports_select0_impl(); }

  /**
   * @brief Return owned index memory usage when provided by the implementation.
   * @return Size of the object and its owned auxiliary storage in bytes.
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
