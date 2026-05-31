#pragma once

#include <pixie/rmq/rmq_base.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace pixie::rmq {

/**
 * @brief Static sparse-table RMQ baseline.
 *
 * @details Stores the index of the first minimum for each power-of-two range.
 * Query time is O(1), build time is O(n log n), and storage is O(n log n)
 * indices. The indexed values are not owned and must outlive this object.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class SparseTable : public RmqBase<SparseTable<T, Compare, Index>, T> {
 public:
  static constexpr std::size_t npos =
      RmqBase<SparseTable<T, Compare, Index>, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  /**
   * @brief Construct an empty sparse table.
   */
  SparseTable() = default;

  /**
   * @brief Build a sparse table over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller index as the RMQ answer.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit SparseTable(std::span<const T> values, Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  /**
   * @brief Return the number of indexed values.
   *
   * @return `values.size()` from construction.
   */
  std::size_t size_impl() const { return values_.size(); }

  /**
   * @brief Return the value at an indexed position.
   *
   * @param position Zero-based position in the indexed values.
   * @return Copy of the value at @p position.
   */
  T value_at_impl(std::size_t position) const { return values_[position]; }

  /**
   * @brief Return the first minimum position in [@p left, @p right].
   *
   * @details Answers in O(1) by comparing the two power-of-two ranges covering
   * the inclusive query interval. Ties return the smaller position.
   *
   * @param left First position in the query range.
   * @param right Last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left > right || right >= values_.size()) {
      return npos;
    }
    const std::size_t length = right - left + 1;
    const std::size_t level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << level;
    const std::size_t first = table_[level][left];
    const std::size_t second = table_[level][right + 1 - span];
    return better(first, second);
  }

 private:
  /**
   * @brief Choose the better of two candidate positions.
   *
   * @details `npos` is treated as missing. If both values compare equal, the
   * smaller position wins to preserve first-minimum semantics.
   *
   * @param left First candidate position, or `npos`.
   * @param right Second candidate position, or `npos`.
   * @return Position of the selected candidate.
   */
  std::size_t better(std::size_t left, std::size_t right) const {
    if (left == npos) {
      return right;
    }
    if (right == npos) {
      return left;
    }
    if (compare_(values_[right], values_[left])) {
      return right;
    }
    if (compare_(values_[left], values_[right])) {
      return left;
    }
    return std::min(left, right);
  }

  /**
   * @brief Build all sparse-table levels over the indexed values.
   *
   * @details Level 0 stores singleton positions. Each higher level stores the
   * first minimum of two adjacent half ranges from the previous level.
   *
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  void build() {
    table_.clear();
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ sparse table index type is too small");
    }

    table_.emplace_back(values_.size());
    for (std::size_t i = 0; i < values_.size(); ++i) {
      table_[0][i] = static_cast<Index>(i);
    }

    for (std::size_t span = 2, half = 1; span <= values_.size();
         half = span, span <<= 1) {
      const std::size_t level = table_.size();
      table_.emplace_back(values_.size() - span + 1);
      for (std::size_t i = 0; i < table_[level].size(); ++i) {
        table_[level][i] = static_cast<Index>(
            better(table_[level - 1][i], table_[level - 1][i + half]));
      }
    }
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<std::vector<Index>> table_;
};

}  // namespace pixie::rmq
