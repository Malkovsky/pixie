#pragma once

#include <pixie/rmq/bp_plus_minus_one_rmq.h>
#include <pixie/rmq/rmq_base.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::rmq {

/**
 * @brief General RMQ via Cartesian-tree reduction to ±1 RMQ.
 *
 * @details Builds a stable min Cartesian tree over the indexed values, takes
 * its Euler tour, and answers original RMQ queries as LCA queries implemented
 * by RMQ over the Euler depth sequence. The indexed values are not owned and
 * must outlive this object.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class CartesianTreeRmq
    : public RmqBase<CartesianTreeRmq<T, Compare, Index>, T> {
 public:
  static constexpr std::size_t npos =
      RmqBase<CartesianTreeRmq<T, Compare, Index>, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  /**
   * @brief Construct an empty Cartesian-tree RMQ index.
   */
  CartesianTreeRmq() = default;

  /**
   * @brief Build a Cartesian-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values stay stable: the smaller index remains the first minimum.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit CartesianTreeRmq(std::span<const T> values,
                            Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  /**
   * @brief Copy an RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   */
  CartesianTreeRmq(const CartesianTreeRmq& other)
      : values_(other.values_),
        compare_(other.compare_),
        left_child_(other.left_child_),
        right_child_(other.right_child_),
        first_occurrence_(other.first_occurrence_),
        euler_nodes_(other.euler_nodes_),
        depths_(other.depths_),
        euler_delta_bits_(other.euler_delta_bits_) {
    reset_depth_rmq();
  }

  /**
   * @brief Copy-assign an RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   * @return Reference to this object.
   */
  CartesianTreeRmq& operator=(const CartesianTreeRmq& other) {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = other.compare_;
    left_child_ = other.left_child_;
    right_child_ = other.right_child_;
    first_occurrence_ = other.first_occurrence_;
    euler_nodes_ = other.euler_nodes_;
    depths_ = other.depths_;
    euler_delta_bits_ = other.euler_delta_bits_;
    reset_depth_rmq();
    return *this;
  }

  /**
   * @brief Move an RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   */
  CartesianTreeRmq(CartesianTreeRmq&& other) noexcept
      : values_(other.values_),
        compare_(std::move(other.compare_)),
        left_child_(std::move(other.left_child_)),
        right_child_(std::move(other.right_child_)),
        first_occurrence_(std::move(other.first_occurrence_)),
        euler_nodes_(std::move(other.euler_nodes_)),
        depths_(std::move(other.depths_)),
        euler_delta_bits_(std::move(other.euler_delta_bits_)) {
    reset_depth_rmq();
  }

  /**
   * @brief Move-assign an RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   * @return Reference to this object.
   */
  CartesianTreeRmq& operator=(CartesianTreeRmq&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = std::move(other.compare_);
    left_child_ = std::move(other.left_child_);
    right_child_ = std::move(other.right_child_);
    first_occurrence_ = std::move(other.first_occurrence_);
    euler_nodes_ = std::move(other.euler_nodes_);
    depths_ = std::move(other.depths_);
    euler_delta_bits_ = std::move(other.euler_delta_bits_);
    reset_depth_rmq();
    return *this;
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
   * @details Converts the query to an LCA query over the Cartesian-tree Euler
   * tour and returns the corresponding original array position. Ties return the
   * smaller original position because the Cartesian tree is stable.
   *
   * @param left First position in the query range.
   * @param right Last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left > right || right >= values_.size()) {
      return npos;
    }
    std::size_t first = first_occurrence_[left];
    std::size_t second = first_occurrence_[right];
    if (first > second) {
      std::swap(first, second);
    }
    const std::size_t euler_position = depth_rmq_.arg_min(first, second);
    if (euler_position == npos) {
      return npos;
    }
    return euler_nodes_[euler_position];
  }

  /**
   * @brief Return the Euler-tour node sequence used by the reduction.
   *
   * @return Non-owning span of original array positions in Euler-tour order.
   */
  std::span<const Index> euler_nodes() const { return euler_nodes_; }

  /**
   * @brief Return the Euler-tour depth sequence used by the reduction.
   *
   * @return Non-owning span of depths corresponding to `euler_nodes()`.
   */
  std::span<const std::int64_t> euler_depths() const { return depths_; }

 private:
  /**
   * @brief Rebuild all Cartesian-tree and Euler-tour auxiliary data.
   *
   * @details Clears previous state, builds a stable Cartesian tree, records its
   * Euler tour, converts adjacent Euler-depth deltas to bits, and rebuilds the
   * ±1 RMQ backend.
   *
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  void build() {
    left_child_.clear();
    right_child_.clear();
    first_occurrence_.clear();
    euler_nodes_.clear();
    depths_.clear();
    euler_delta_bits_.clear();
    depth_rmq_ = BpPlusMinusOneRmq<Index>();

    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("Cartesian RMQ index type is too small");
    }

    left_child_.assign(values_.size(), invalid_index);
    right_child_.assign(values_.size(), invalid_index);

    const std::size_t root = build_cartesian_tree();
    first_occurrence_.assign(values_.size(), invalid_index);
    euler_nodes_.reserve(2 * values_.size() - 1);
    depths_.reserve(2 * values_.size() - 1);
    euler_tour(root, 0);
    build_euler_delta_bits();
    reset_depth_rmq();
  }

  /**
   * @brief Build the stable min Cartesian tree.
   *
   * @details Uses the standard monotone-stack construction. Strictly smaller
   * values become ancestors; equal values are not popped, preserving first
   * minimum tie-breaking.
   *
   * @return Root node position in the original value array.
   */
  std::size_t build_cartesian_tree() {
    std::vector<Index> stack;
    stack.reserve(values_.size());

    for (std::size_t i = 0; i < values_.size(); ++i) {
      Index last = invalid_index;
      while (!stack.empty() && compare_(values_[i], values_[stack.back()])) {
        last = stack.back();
        stack.pop_back();
      }
      if (last != invalid_index) {
        left_child_[i] = last;
      }
      if (!stack.empty()) {
        right_child_[stack.back()] = static_cast<Index>(i);
      }
      stack.push_back(static_cast<Index>(i));
    }

    return stack.front();
  }

  /**
   * @brief Append the Euler tour of a Cartesian-tree subtree.
   *
   * @details Visits @p node, recurses into each existing child, and appends
   * @p node again after returning from that child.
   *
   * @param node Current Cartesian-tree node.
   * @param depth Depth of @p node in the Cartesian tree.
   */
  void euler_tour(std::size_t node, std::int64_t depth) {
    append_euler(node, depth);
    if (left_child_[node] != invalid_index) {
      euler_tour(left_child_[node], depth + 1);
      append_euler(node, depth);
    }
    if (right_child_[node] != invalid_index) {
      euler_tour(right_child_[node], depth + 1);
      append_euler(node, depth);
    }
  }

  /**
   * @brief Append one node/depth pair to the Euler-tour arrays.
   *
   * @details Records the first Euler occurrence of @p node if this is the first
   * time the node is appended.
   *
   * @param node Cartesian-tree node, also an original value position.
   * @param depth Depth of @p node in the Cartesian tree.
   */
  void append_euler(std::size_t node, std::int64_t depth) {
    if (first_occurrence_[node] == invalid_index) {
      first_occurrence_[node] = static_cast<Index>(euler_nodes_.size());
    }
    euler_nodes_.push_back(static_cast<Index>(node));
    depths_.push_back(depth);
  }

  /**
   * @brief Rebuild the ±1 RMQ backend over the current Euler-depth deltas.
   *
   * @details Called after build, copy, and move operations because the backend
   * stores non-owning spans into this object's `euler_delta_bits_` storage.
   */
  void reset_depth_rmq() {
    depth_rmq_ = BpPlusMinusOneRmq<Index>(
        std::span<const std::uint64_t>(euler_delta_bits_), depths_.size());
  }

  /**
   * @brief Pack adjacent Euler-depth changes into BP-style delta bits.
   *
   * @details Bit `1` means the next Euler depth is current depth + 1; bit `0`
   * means current depth - 1. Cartesian-tree Euler tours have only ±1 adjacent
   * depth changes.
   */
  void build_euler_delta_bits() {
    euler_delta_bits_.assign((depths_.size() - 1 + 63) / 64, 0);
    for (std::size_t i = 1; i < depths_.size(); ++i) {
      const std::int64_t delta = depths_[i] - depths_[i - 1];
      if (delta == 1) {
        euler_delta_bits_[(i - 1) >> 6] |= std::uint64_t{1} << ((i - 1) & 63);
      }
    }
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<Index> left_child_;
  std::vector<Index> right_child_;
  std::vector<Index> first_occurrence_;
  std::vector<Index> euler_nodes_;
  std::vector<std::int64_t> depths_;
  std::vector<std::uint64_t> euler_delta_bits_;
  BpPlusMinusOneRmq<Index> depth_rmq_;
};

}  // namespace pixie::rmq
