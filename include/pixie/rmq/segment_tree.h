#pragma once

#include <pixie/memory_usage.h>
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
 * @brief Static iterative segment-tree RMQ baseline.
 *
 * @details Stores the index of the first minimum for each segment in a flat
 * binary tree. Query time is O(log n), build time is O(n), and storage is O(n)
 * indices. The indexed values are not owned and must outlive this object.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class SegmentTree : public RmqBase<SegmentTree<T, Compare, Index>, T> {
 public:
  static constexpr std::size_t npos =
      RmqBase<SegmentTree<T, Compare, Index>, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  /**
   * @brief Construct an empty segment tree.
   */
  SegmentTree() = default;

  /**
   * @brief Build an iterative segment tree over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller index as the RMQ answer.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit SegmentTree(std::span<const T> values, Compare compare = Compare())
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
   * @brief Return the first minimum position in [@p left, @p right).
   *
   * @details Answers in O(log n) by walking the flat iterative segment tree.
   * Ties return the smaller position.
   *
   * @param left First position in the query range.
   * @param right One past the last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }

    left += leaf_base_;
    right += leaf_base_;
    std::size_t answer = npos;
    while (left < right) {
      if ((left & 1u) != 0) {
        answer = better(answer, tree_[left]);
        ++left;
      }
      if ((right & 1u) != 0) {
        --right;
        answer = better(answer, tree_[right]);
      }
      left >>= 1;
      right >>= 1;
    }
    return answer;
  }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this segment-tree object and its flat index buffer. The
   * external input values are not owned and are excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    return sizeof(*this) + pixie::vector_capacity_bytes(tree_);
  }

 private:
  /**
   * @brief Choose the better of two candidate positions.
   *
   * @details `npos` and `invalid_index` are treated as missing. If both values
   * compare equal, the smaller position wins to preserve first-minimum
   * semantics.
   *
   * @param left First candidate position, `npos`, or `invalid_index`.
   * @param right Second candidate position, `npos`, or `invalid_index`.
   * @return Position of the selected candidate.
   */
  std::size_t better(std::size_t left, std::size_t right) const {
    if (left == npos || left == invalid_index) {
      return right;
    }
    if (right == npos || right == invalid_index) {
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
   * @brief Choose the better child while building the tree.
   *
   * @details Child ranges are disjoint and ordered, so ties keep the left
   * child. Missing tail leaves use `invalid_index`.
   *
   * @param left Candidate from the left child.
   * @param right Candidate from the right child.
   * @return Candidate for the parent node, or `invalid_index`.
   */
  Index build_better(Index left, Index right) const {
    if (left == invalid_index) {
      return right;
    }
    if (right == invalid_index) {
      return left;
    }
    return compare_(values_[right], values_[left]) ? right : left;
  }

  /**
   * @brief Build the flat iterative segment tree.
   *
   * @details Leaves start at `leaf_base_`, which is the next power of two.
   * Unused leaves contain `invalid_index`, and internal nodes store the first
   * minimum position of their covered segment.
   *
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  void build() {
    tree_.clear();
    leaf_base_ = 0;
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ segment tree index type is too small");
    }

    leaf_base_ = std::bit_ceil(values_.size());
    tree_.clear();
    tree_.resize(2 * leaf_base_);
    for (std::size_t i = 0; i < values_.size(); ++i) {
      tree_[leaf_base_ + i] = static_cast<Index>(i);
    }
    std::fill(tree_.begin() + leaf_base_ + values_.size(), tree_.end(),
              invalid_index);
    for (std::size_t node = leaf_base_; node > 1;) {
      --node;
      tree_[node] = build_better(tree_[node << 1], tree_[(node << 1) | 1]);
    }
  }

  std::span<const T> values_;
  Compare compare_;
  std::size_t leaf_base_ = 0;
  std::vector<Index> tree_;
};

}  // namespace pixie::rmq
