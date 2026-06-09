#pragma once

#include <pixie/bitvector.h>
#include <pixie/rmq/bp_plus_minus_one_rmq.h>
#include <pixie/rmq/rmq_base.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::rmq {

/**
 * @brief General RMQ via the Ferrada-Navarro BP Cartesian-tree encoding.
 *
 * @details Builds the balanced-parentheses representation of the Cartesian-tree
 * RMQ information directly from the indexed values. The construction scans the
 * input right-to-left with a monotone stack and emits the BP sequence used by
 * the Ferrada-Navarro formulation. Queries map the half-open value interval to
 * positions in that BP sequence with `select0`, run a ±1 RMQ over BP prefix
 * excess, and map the winning BP position back to an array index with `rank0`.
 *
 * The indexed values are not owned and must outlive this object. `arg_min()`
 * answers from the BP indexes alone; `range_min()` still reads the external
 * value span after the position is known. Equal values are handled stably: the
 * smaller original position remains the first minimum.
 *
 * This implementation is a practical BP Cartesian-tree backend, not a fully
 * compressed 2n + o(n)-bit object: it stores the 2n BP bits plus the supporting
 * `BitVector` rank/select index and a `BpPlusMinusOneRmq` over BP excess.
 *
 * References:
 * - Johannes Fischer, "Optimal Succinctness for Range Minimum Queries",
 *   LATIN 2010; arXiv:0812.2775.
 * - Hector Ferrada and Gonzalo Navarro, "Improved Range Minimum Queries",
 *   Data Compression Conference 2016; Journal of Discrete Algorithms 43
 *   (2017), doi:10.1016/j.jda.2016.09.002.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class CartesianTreeRmq
    : public RmqBase<CartesianTreeRmq<T, Compare, Index>, T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "CartesianTreeRmq index type must be unsigned");

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
        bp_bits_(other.bp_bits_),
        bp_bit_count_(other.bp_bit_count_) {
    reset_bp_indexes();
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
    bp_bits_ = other.bp_bits_;
    bp_bit_count_ = other.bp_bit_count_;
    reset_bp_indexes();
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
        bp_bits_(std::move(other.bp_bits_)),
        bp_bit_count_(other.bp_bit_count_) {
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    reset_bp_indexes();
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
    bp_bits_ = std::move(other.bp_bits_);
    bp_bit_count_ = other.bp_bit_count_;
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    reset_bp_indexes();
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
   * @brief Return the first minimum position in [@p left, @p right).
   *
   * @details Adapts the Ferrada-Navarro BP formula to Pixie's zero-based
   * half-open ranges. Ties return the smaller original position.
   *
   * @param left First position in the query range.
   * @param right One past the last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    const BitVector& bp_index = *bp_index_;
    const std::size_t first_close = bp_index.select0(left + 1);
    const std::size_t last_close = bp_index.select0(right);
    if (first_close == bp_bit_count_ || last_close == bp_bit_count_ ||
        first_close > last_close) {
      return npos;
    }

    const std::size_t shifted_min =
        bp_depth_rmq_.arg_min(first_close + 1, last_close + 2);
    if (shifted_min == npos || shifted_min == 0) {
      return npos;
    }
    const std::size_t answer = bp_index.rank0(shifted_min) - 1;
    return answer < values_.size() ? answer : npos;
  }

  /**
   * @brief Return the number of BP bits in the Cartesian-tree RMQ encoding.
   *
   * @return Logical BP bit count, equal to `2 * size()`.
   */
  std::size_t bp_bit_count() const { return bp_bit_count_; }

  /**
   * @brief Return the packed BP words used by the RMQ encoding.
   *
   * @return Non-owning span of little-endian packed BP words.
   */
  std::span<const std::uint64_t> bp_words() const { return bp_bits_; }

 private:
  /**
   * @brief Rebuild the direct BP Cartesian-tree RMQ representation.
   *
   * @details Constructs the Ferrada-Navarro BP sequence with a right-to-left
   * monotone stack, then rebuilds select/rank and excess-min indexes over the
   * packed BP words.
   *
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  void build() {
    bp_bits_.clear();
    bp_bit_count_ = 0;
    reset_bp_indexes();

    if (values_.empty()) {
      return;
    }
    if (values_.size() > (static_cast<std::size_t>(invalid_index) - 1) / 2) {
      throw std::length_error("Cartesian RMQ index type is too small");
    }

    bp_bit_count_ = 2 * values_.size();
    bp_bits_.assign((bp_bit_count_ + 63) / 64, 0);
    build_bp_bits();
    reset_bp_indexes();
  }

  /**
   * @brief Build the BP bits with the Ferrada-Navarro stack construction.
   *
   * @details The paper describes prepending bits while scanning right-to-left.
   * This implementation fills the destination bitvector from right to left.
   */
  void build_bp_bits() {
    const auto stack = std::make_unique_for_overwrite<Index[]>(values_.size());
    std::size_t stack_size = 0;
    std::size_t write_position = bp_bit_count_;

    for (std::size_t i = values_.size(); i-- > 0;) {
      while (stack_size != 0 &&
             !compare_(values_[stack[stack_size - 1]], values_[i])) {
        --stack_size;
        prepend_bp_bit(write_position, true);
      }
      stack[stack_size++] = static_cast<Index>(i);
      prepend_bp_bit(write_position, false);
    }

    while (write_position != 0) {
      prepend_bp_bit(write_position, true);
    }
  }

  /**
   * @brief Prepend one BP bit into the right-to-left construction buffer.
   *
   * @param write_position Current first unwritten position; decremented here.
   * @param bit `true` for an opening parenthesis and `false` for closing.
   */
  void prepend_bp_bit(std::size_t& write_position, bool bit) {
    --write_position;
    if (bit) {
      bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                       << (write_position & 63);
    }
  }

  /**
   * @brief Rebuild indexes that store non-owning spans into `bp_bits_`.
   *
   * @details Called after build, copy, and move operations because both indexes
   * keep views into this object's packed BP word storage.
   */
  void reset_bp_indexes() {
    bp_index_.reset();
    bp_depth_rmq_ = BpPlusMinusOneRmq<Index>();
    if (bp_bit_count_ == 0) {
      return;
    }
    const std::span<const std::uint64_t> words(bp_bits_);
    bp_index_.emplace(words, bp_bit_count_);
    bp_depth_rmq_ = BpPlusMinusOneRmq<Index>(words, bp_bit_count_ + 1);
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<std::uint64_t> bp_bits_;
  std::size_t bp_bit_count_ = 0;
  std::optional<BitVector> bp_index_;
  BpPlusMinusOneRmq<Index> bp_depth_rmq_;
};

}  // namespace pixie::rmq
