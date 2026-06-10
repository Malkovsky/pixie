#pragma once

/**
 * Benchmark snapshot: experimental rmM-backed Cartesian/BP RMQ, 2026-06-10.
 *
 * Command shape:
 * taskset -c 0 ./build/release/bench_rmq
 *   --benchmark_filter='^(rmq_cartesian_tree|rmq_cartesian_tree_rmm_btree|
 *     rmq_bp_plus_minus_one|rmq_bp_rmm_btree)/...|^(rmq_build_...)/...'
 *   --benchmark_min_time=0.30s
 *   --benchmark_out=/tmp/rmq_rmm_btree_20260610.json
 *
 * Query CPU time, ns:
 * | N        | width    | cartesian | cartesian+rmm | bp +/-1 | bp+rmm  |
 * |----------|----------|-----------|---------------|---------|---------|
 * | 2^10     | 64       |      80.5 |         155.5 |    42.8 |   107.1 |
 * |----------|----------|-----------|---------------|---------|---------|
 * | 2^14     | 64       |      98.8 |         177.4 |    42.7 |   117.9 |
 * | 2^14     | 4096     |     123.9 |         722.2 |    69.0 |   634.9 |
 * |----------|----------|-----------|---------------|---------|---------|
 * | 2^18     | 64       |     110.8 |         176.7 |    42.5 |   118.1 |
 * | 2^18     | 4096     |     156.0 |         744.6 |    66.6 |   735.1 |
 * | 2^18     | 2^18     |     165.8 |        1222.1 |    67.6 |  1204.5 |
 * |----------|----------|-----------|---------------|---------|---------|
 * | 2^22     | 64       |     199.1 |         328.6 |    46.5 |   131.7 |
 * | 2^22     | 4096     |     327.0 |        1271.9 |    95.8 |   772.3 |
 * | 2^22     | 2^22     |     312.8 |        1601.3 |    88.6 |  1498.6 |
 * |----------|----------|-----------|---------------|---------|---------|
 * | 2^24     | 64       |     272.5 |         396.2 |    77.0 |   202.9 |
 * | 2^24     | 4096     |     500.1 |        1248.3 |   189.6 |  1038.2 |
 * | 2^24     | 2^24     |     483.5 |        1832.2 |   126.7 |  1687.9 |
 *
 * Build CPU time, ms:
 * | N        | cartesian | cartesian+rmm | bp +/-1 | bp+rmm |
 * |----------|-----------|---------------|---------|--------|
 * | 2^10     |     0.006 |         0.008 |   0.000 |  0.005 |
 * | 2^14     |     0.053 |         0.074 |   0.002 |  0.026 |
 * | 2^18     |     1.770 |         2.243 |   0.047 |  0.377 |
 * | 2^22     |    43.794 |        35.318 |   1.230 |  6.164 |
 * | 2^24     |   192.036 |       155.852 |  20.232 | 23.921 |
 */

#include <pixie/experimental/rmm_btree.h>
#include <pixie/rmq/rmq_base.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::rmq::experimental {

/**
 * @brief Depth RMQ adapter over the experimental rmM btree.
 *
 * @details The input words encode adjacent depth deltas: bit 1 means +1 and
 * bit 0 means -1. The indexed depth sequence has `depth_count` prefix
 * positions, so the bit sequence has `depth_count - 1` bits. Queries return the
 * first position of the minimum depth in a half-open depth-position range.
 *
 * `pixie::experimental::RmMBTree::range_min_query_pos()` works on inclusive
 * bit ranges and reports extrema after consuming each bit. This adapter adds
 * the missing left prefix boundary comparison, preserving first-minimum ties.
 */
template <class Index = std::size_t,
          std::size_t HighCacheLines = 4,
          std::size_t LowFanout = 32>
class RmMBTreePlusMinusOneRmq {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "RmMBTreePlusMinusOneRmq index type must be unsigned");

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kBlockSize =
      pixie::experimental::RmMBTree<HighCacheLines, LowFanout>::kBlockBits;

  RmMBTreePlusMinusOneRmq() = default;

  /**
   * @brief Build an rmM-backed ±1 RMQ over external packed delta bits.
   *
   * @param bits Little-endian packed delta bits.
   * @param depth_count Number of indexed depth positions.
   * @throws std::length_error if @p Index cannot represent all positions.
   * @throws std::invalid_argument if @p bits is shorter than required.
   */
  RmMBTreePlusMinusOneRmq(std::span<const std::uint64_t> bits,
                          std::size_t depth_count) {
    build(bits, depth_count);
  }

  /**
   * @brief Rebuild this adapter over external packed delta bits.
   */
  void build(std::span<const std::uint64_t> bits, std::size_t depth_count) {
    words_ = bits;
    depth_count_ = depth_count;
    rmm_ = RmM();

    if (depth_count_ == 0) {
      return;
    }
    if (depth_count_ > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ rmM btree index type is too small");
    }
    rmm_ = RmM(words_, depth_count_ - 1);
  }

  /**
   * @brief Return the number of indexed depth positions.
   */
  std::size_t size() const { return depth_count_; }

  /**
   * @brief Whether the indexed depth sequence is empty.
   */
  bool empty() const { return depth_count_ == 0; }

  /**
   * @brief Return the first minimum depth position in [@p left, @p right).
   *
   * @details Empty or invalid ranges return `npos`. Equal minima return the
   * smaller depth position.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    if (left >= right || right > depth_count_) {
      return npos;
    }
    if (right == left + 1) {
      return left;
    }

    const std::size_t bit_left = left;
    const std::size_t bit_right = right - 2;
    const int minimum_after_left =
        rmm_.range_min_query_val(bit_left, bit_right);
    if (minimum_after_left >= 0) {
      return left;
    }

    const std::size_t bit_position =
        rmm_.range_min_query_pos(bit_left, bit_right);
    return bit_position == RmM::npos ? npos : bit_position + 1;
  }

  /**
   * @brief Return the one-based rank-th closing parenthesis position.
   */
  std::size_t select0(std::size_t rank) const { return rmm_.select0(rank); }

  /**
   * @brief Count closing parentheses in the prefix [0, @p end_position).
   */
  std::size_t rank0(std::size_t end_position) const {
    return rmm_.rank0(end_position);
  }

 private:
  using RmM = pixie::experimental::RmMBTree<HighCacheLines, LowFanout>;

  std::span<const std::uint64_t> words_;
  std::size_t depth_count_ = 0;
  RmM rmm_;
};

/**
 * @brief Experimental Ferrada-Navarro Cartesian-tree RMQ using rmM support.
 *
 * @details This class keeps the same public value-RMQ contract as
 * `CartesianTreeRmq`, but replaces the `BitVector` + `BpPlusMinusOneRmq`
 * support layer with `RmMBTreePlusMinusOneRmq`, which in turn uses the
 * experimental range min-max btree over the BP excess sequence.
 *
 * This is intentionally not included from `pixie/rmq.h`; it exists to benchmark
 * the Ferrada-Navarro BP formula with an rmM-backed `rmq_D` primitive before
 * deciding whether to promote or optimize it.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t HighCacheLines = 4,
          std::size_t LowFanout = 32>
class CartesianTreeRmMBTreeRmq
    : public RmqBase<CartesianTreeRmMBTreeRmq<T,
                                              Compare,
                                              Index,
                                              HighCacheLines,
                                              LowFanout>,
                     T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "CartesianTreeRmMBTreeRmq index type must be unsigned");

  static constexpr std::size_t npos = RmqBase<
      CartesianTreeRmMBTreeRmq<T, Compare, Index, HighCacheLines, LowFanout>,
      T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  CartesianTreeRmMBTreeRmq() = default;

  /**
   * @brief Build an experimental Cartesian-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values stay stable: the smaller index remains the first minimum.
   */
  explicit CartesianTreeRmMBTreeRmq(std::span<const T> values,
                                    Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  CartesianTreeRmMBTreeRmq(const CartesianTreeRmMBTreeRmq& other)
      : values_(other.values_),
        compare_(other.compare_),
        bp_bits_(other.bp_bits_),
        bp_bit_count_(other.bp_bit_count_) {
    reset_bp_index();
  }

  CartesianTreeRmMBTreeRmq& operator=(const CartesianTreeRmMBTreeRmq& other) {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = other.compare_;
    bp_bits_ = other.bp_bits_;
    bp_bit_count_ = other.bp_bit_count_;
    reset_bp_index();
    return *this;
  }

  CartesianTreeRmMBTreeRmq(CartesianTreeRmMBTreeRmq&& other) noexcept
      : values_(other.values_),
        compare_(std::move(other.compare_)),
        bp_bits_(std::move(other.bp_bits_)),
        bp_bit_count_(other.bp_bit_count_) {
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    reset_bp_index();
  }

  CartesianTreeRmMBTreeRmq& operator=(
      CartesianTreeRmMBTreeRmq&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = std::move(other.compare_);
    bp_bits_ = std::move(other.bp_bits_);
    bp_bit_count_ = other.bp_bit_count_;
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    reset_bp_index();
    return *this;
  }

  /**
   * @brief Return the number of indexed values.
   */
  std::size_t size_impl() const { return values_.size(); }

  /**
   * @brief Return the value at an indexed position.
   */
  T value_at_impl(std::size_t position) const { return values_[position]; }

  /**
   * @brief Return the first minimum position in [@p left, @p right).
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    const std::size_t first_close = bp_support_.select0(left + 1);
    const std::size_t last_close = bp_support_.select0(right);
    if (first_close == BpSupport::npos || last_close == BpSupport::npos ||
        first_close > last_close) {
      return npos;
    }

    const std::size_t shifted_min =
        bp_support_.arg_min(first_close + 1, last_close + 2);
    if (shifted_min == npos || shifted_min == 0) {
      return npos;
    }
    const std::size_t answer = bp_support_.rank0(shifted_min) - 1;
    return answer < values_.size() ? answer : npos;
  }

  /**
   * @brief Return the number of BP bits in the Cartesian-tree RMQ encoding.
   */
  std::size_t bp_bit_count() const { return bp_bit_count_; }

  /**
   * @brief Return the packed BP words used by the RMQ encoding.
   */
  std::span<const std::uint64_t> bp_words() const { return bp_bits_; }

 private:
  using BpSupport = RmMBTreePlusMinusOneRmq<Index, HighCacheLines, LowFanout>;

  void build() {
    bp_bits_.clear();
    bp_bit_count_ = 0;
    bp_support_ = BpSupport();

    if (values_.empty()) {
      return;
    }
    if (values_.size() > (static_cast<std::size_t>(invalid_index) - 1) / 2) {
      throw std::length_error("Cartesian rmM RMQ index type is too small");
    }

    bp_bit_count_ = 2 * values_.size();
    bp_bits_.assign((bp_bit_count_ + 63) / 64, 0);
    build_bp_bits();
    reset_bp_index();
  }

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

  void prepend_bp_bit(std::size_t& write_position, bool bit) {
    --write_position;
    if (bit) {
      bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                       << (write_position & 63);
    }
  }

  void reset_bp_index() {
    bp_support_ = BpSupport();
    if (bp_bit_count_ == 0) {
      return;
    }
    bp_support_ =
        BpSupport(std::span<const std::uint64_t>(bp_bits_), bp_bit_count_ + 1);
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<std::uint64_t> bp_bits_;
  std::size_t bp_bit_count_ = 0;
  BpSupport bp_support_;
};

}  // namespace pixie::rmq::experimental
