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

  SegmentTree() = default;

  explicit SegmentTree(std::span<const T> values, Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  std::size_t size_impl() const { return values_.size(); }

  T value_at_impl(std::size_t position) const { return values_[position]; }

  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left > right || right >= values_.size()) {
      return npos;
    }

    left += leaf_base_;
    right += leaf_base_;
    std::size_t answer = npos;
    while (left <= right) {
      if ((left & 1u) != 0) {
        answer = better(answer, tree_[left]);
        ++left;
      }
      if ((right & 1u) == 0) {
        answer = better(answer, tree_[right]);
        if (right == 0) {
          break;
        }
        --right;
      }
      left >>= 1;
      right >>= 1;
    }
    return answer;
  }

 private:
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
    tree_.assign(2 * leaf_base_, invalid_index);
    for (std::size_t i = 0; i < values_.size(); ++i) {
      tree_[leaf_base_ + i] = static_cast<Index>(i);
    }
    for (std::size_t node = leaf_base_; node > 1;) {
      --node;
      tree_[node] =
          static_cast<Index>(better(tree_[node << 1], tree_[(node << 1) | 1]));
    }
  }

  std::span<const T> values_;
  Compare compare_;
  std::size_t leaf_base_ = 0;
  std::vector<Index> tree_;
};

}  // namespace pixie::rmq
