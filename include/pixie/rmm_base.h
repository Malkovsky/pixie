#pragma once

#include <cstddef>
#include <limits>

namespace pixie {

/**
 * @brief CRTP facade for rank/select and range min-max tree operations.
 *
 * Positions are zero-based unless explicitly documented otherwise. Operations
 * returning positions use npos when the requested value does not exist or the
 * query is outside the valid domain. Balanced-parentheses navigation follows
 * SDSL-style zero-based parenthesis indexing: close/open/enclose accept and
 * return bit positions, not prefix-boundary positions.
 */
template <class Impl>
class RmMBase {
 public:
  /**
   * @brief Sentinel returned by position queries when no valid answer exists.
   */
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

  /**
   * @brief Number of bits in the represented sequence.
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Count 1 bits in the prefix [0, end_position).
   */
  std::size_t rank1(std::size_t end_position) const {
    return impl().rank1_impl(end_position);
  }

  /**
   * @brief Count 0 bits in the prefix [0, end_position).
   */
  std::size_t rank0(std::size_t end_position) const {
    return impl().rank0_impl(end_position);
  }

  /**
   * @brief Return the zero-based position of the rank-th 1 bit.
   * @param rank One-based rank of the requested 1 bit.
   */
  std::size_t select1(std::size_t rank) const {
    return impl().select1_impl(rank);
  }

  /**
   * @brief Return the zero-based position of the rank-th 0 bit.
   * @param rank One-based rank of the requested 0 bit.
   */
  std::size_t select0(std::size_t rank) const {
    return impl().select0_impl(rank);
  }

  /**
   * @brief Count "10" bit patterns fully contained in [0, end_position).
   *
   * The returned count is the number of start positions p such that
   * p + 1 < end_position, bit[p] == 1, and bit[p + 1] == 0.
   */
  std::size_t rank10(std::size_t end_position) const {
    return impl().rank10_impl(end_position);
  }

  /**
   * @brief Return the zero-based start position of the rank-th "10" pattern.
   * @param rank One-based rank of the requested pattern.
   */
  std::size_t select10(std::size_t rank) const {
    return impl().select10_impl(rank);
  }

  /**
   * @brief Prefix excess on [0, end_position): +1 for 1 bits, -1 for 0 bits.
   */
  int excess(std::size_t end_position) const {
    return impl().excess_impl(end_position);
  }

  /**
   * @brief Forward excess search from a prefix boundary.
   *
   * Finds the first bit position p >= start_position such that
   * excess(p + 1) == excess(start_position) + delta.
   */
  std::size_t fwdsearch(std::size_t start_position, int delta) const {
    return impl().fwdsearch_impl(start_position, delta);
  }

  /**
   * @brief Backward excess search from a prefix boundary.
   *
   * Finds the greatest prefix boundary p < start_position such that
   * excess(p) == excess(start_position) + delta.
   */
  std::size_t bwdsearch(std::size_t start_position, int delta) const {
    return impl().bwdsearch_impl(start_position, delta);
  }

  /**
   * @brief Position of the first minimum relative excess in [range_begin,
   * range_end].
   *
   * Relative excess is accumulated from zero over the inclusive bit range.
   */
  std::size_t range_min_query_pos(std::size_t range_begin,
                                  std::size_t range_end) const {
    return impl().range_min_query_pos_impl(range_begin, range_end);
  }

  /**
   * @brief Minimum relative excess value over [range_begin, range_end].
   *
   * Relative excess is accumulated from zero over the inclusive bit range.
   */
  int range_min_query_val(std::size_t range_begin,
                          std::size_t range_end) const {
    return impl().range_min_query_val_impl(range_begin, range_end);
  }

  /**
   * @brief Position of the first maximum relative excess in [range_begin,
   * range_end].
   *
   * Relative excess is accumulated from zero over the inclusive bit range.
   */
  std::size_t range_max_query_pos(std::size_t range_begin,
                                  std::size_t range_end) const {
    return impl().range_max_query_pos_impl(range_begin, range_end);
  }

  /**
   * @brief Maximum relative excess value over [range_begin, range_end].
   *
   * Relative excess is accumulated from zero over the inclusive bit range.
   */
  int range_max_query_val(std::size_t range_begin,
                          std::size_t range_end) const {
    return impl().range_max_query_val_impl(range_begin, range_end);
  }

  /**
   * @brief Count positions attaining the minimum relative excess in
   * [range_begin, range_end].
   */
  std::size_t mincount(std::size_t range_begin, std::size_t range_end) const {
    return impl().mincount_impl(range_begin, range_end);
  }

  /**
   * @brief Select a position attaining the minimum relative excess in
   * [range_begin, range_end].
   * @param rank One-based rank among positions attaining the range minimum.
   */
  std::size_t minselect(std::size_t range_begin,
                        std::size_t range_end,
                        std::size_t rank) const {
    return impl().minselect_impl(range_begin, range_end, rank);
  }

  /**
   * @brief Matching closing parenthesis for the parenthesis at open_position.
   *
   * Uses SDSL-style zero-based indexing. If open_position refers to a closing
   * parenthesis, the expected result is open_position.
   */
  std::size_t close(std::size_t open_position) const {
    return impl().close_impl(open_position);
  }

  /**
   * @brief Matching opening parenthesis for the parenthesis at close_position.
   *
   * Uses SDSL-style zero-based indexing. If close_position refers to an opening
   * parenthesis, the expected result is close_position.
   */
  std::size_t open(std::size_t close_position) const {
    return impl().open_impl(close_position);
  }

  /**
   * @brief Closest opening parenthesis strictly enclosing position.
   *
   * Uses SDSL-style zero-based indexing. If position refers to a closing
   * parenthesis, this is equivalent to open(position).
   */
  std::size_t enclose(std::size_t open_position) const {
    return impl().enclose_impl(open_position);
  }

  /**
   * @brief Read bit at position @p position (LSB-first across words).
   */
  int bit_impl(const size_t& position) const noexcept {
    return impl().bit_impl(position);
  }

 private:
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
