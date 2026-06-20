#pragma once

/**
 * rmM-backed Cartesian RMQ.
 *
 * The benchmark harness registers this retained variant as
 * `rmq_cartesian_rmm` and `rmq_build_cartesian_rmm`.
 * It is included from `pixie/rmq.h` with the rest of the retained RMQ backends.
 *
 * Diagnostic snapshot, 2026-06-13:
 *
 * Focused value-RMQ rows, N=2^22 values, CPU mean across 5 repetitions.
 * Command shape:
 *   taskset -c 0 ./build/release/bench_rmq
 *     --benchmark_filter='^rmq_(cartesian_rmm|sdsl_sct)/4194304/'
 *
 * | max width | CartesianRmM (ns) | SdslSct (ns) |
 * | --------: | ----------------: | -----------: |
 * |        64 |             121.0 |        211.0 |
 * |      4096 |             310.0 |        757.0 |
 * |      2^18 |             529.0 |       1057.0 |
 * |      2^22 |             598.0 |        979.0 |
 *
 * Raw BP rmM rows at 2^23 bits, Q=32768, CPU mean across 5 repetitions.
 * Command shape:
 *   ./build/release/bench_rmm_{btree,sdsl}
 *     --ops=range_min_query_pos,range_min_query_val
 *     --explicit_sizes=8388608 --Q=32768
 *
 * | backend             | range_min_pos (ns) | range_min_val (ns) |
 * | ------------------- | -----------------: | -----------------: |
 * | RmMBTree            |              602.0 |              540.0 |
 * | SdslRmMTree         |              727.0 |              742.0 |
 *
 * The 2026-06-13 optimization pass streams cover nodes instead of materializing
 * a zero-initialized per-query cover, adds a combined min position/value query
 * for the Cartesian adapter, and uses a minimum-only boundary scanner for
 * range-min position/value queries. Current perf samples for raw
 * range_min_query_pos are concentrated in for_each_cover_node(),
 * scan_min_range(), and summary_at().
 *
 * Historical build snapshot after select0-only BP rank/select construction,
 * 2026-06-13, before the RmMBTree full-block summary fast path.
 * Command shape:
 *   taskset -c 0 ./build/release/bench_rmq
 *     --benchmark_filter='^rmq_build_(cartesian_rmm|cartesian_hybrid_btree|sdsl_sct)/(4194304|67108864)$'
 *     --benchmark_repetitions=5
 *
 * CPU mean, milliseconds.
 *
 * | N    | CartesianRmM | CartesianHybrid | SdslSct |
 * | ---: | -----------: | --------------: | ------: |
 * | 2^22 |       54.407 |          49.361 |  41.317 |
 * | 2^26 |      915.712 |         795.562 | 666.048 |
 *
 * Historical construction perf profile snapshot, 2026-06-13, before the
 * RmMBTree full-block summary fast path.
 * Profiling build: RelWithDebInfo, -O3, debug info, frame pointers. Rows are
 * N=2^26 build benchmarks sampled with perf at 999 Hz.
 *
 * CartesianHybridBTree: 867 ms CPU in the profiled run. BP construction is the
 * dominant cost: the succinct monotone-stack operations plus BP writes account
 * for roughly two thirds of samples. Top sparse-table construction is about 6%.
 * CartesianRmM: 1013 ms CPU in the profiled run. BP construction was still the
 * dominant cost, while rmM summary construction added a visible second cost:
 * summarize_bits()/bit() lines accounted for about 14%. Benchmark dataset
 * generation was about 5% in both rows.
 */

#include <pixie/experimental/rmm_btree.h>
#include <pixie/memory_usage.h>
#include <pixie/rmq/rmq_base.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::rmq {

namespace detail {

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
class RmMPlusMinusOne {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "RmMPlusMinusOne index type must be unsigned");

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kBlockSize =
      pixie::experimental::RmMBTree<HighCacheLines, LowFanout>::kBlockBits;

  RmMPlusMinusOne() = default;

  /**
   * @brief Build an rmM-backed ±1 RMQ over external packed delta bits.
   *
   * @param bits Little-endian packed delta bits.
   * @param depth_count Number of indexed depth positions.
   * @throws std::length_error if @p Index cannot represent all positions.
   * @throws std::invalid_argument if @p bits is shorter than required.
   */
  RmMPlusMinusOne(std::span<const std::uint64_t> bits,
                  std::size_t depth_count) {
    build(bits, depth_count);
  }

  /**
   * @brief Build with an exact count of +1 delta bits.
   *
   * @details Cartesian BP construction knows this count exactly, which lets the
   * nested BitVector allocate select samples without a second scan.
   */
  RmMPlusMinusOne(std::span<const std::uint64_t> bits,
                  std::size_t depth_count,
                  std::size_t one_count) {
    build(bits, depth_count, one_count);
  }

  /**
   * @brief Rebuild this adapter over external packed delta bits.
   */
  void build(std::span<const std::uint64_t> bits, std::size_t depth_count) {
    build(bits, depth_count, std::nullopt);
  }

  /**
   * @brief Rebuild this adapter with an optional exact count of +1 delta bits.
   */
  void build(std::span<const std::uint64_t> bits,
             std::size_t depth_count,
             std::optional<std::size_t> one_count) {
    words_ = bits;
    depth_count_ = depth_count;
    rmm_ = RmM();

    if (depth_count_ == 0) {
      return;
    }
    if (depth_count_ > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ rmM btree index type is too small");
    }
    rmm_ = RmM(words_, depth_count_ - 1, BitVector::SelectSupport::kSelect0,
               one_count);
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
    const auto minimum_after_left =
        rmm_.range_min_query_result(bit_left, bit_right);
    if (minimum_after_left.value >= 0) {
      return left;
    }

    const std::size_t bit_position = minimum_after_left.position;
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

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this ±1 RMQ adapter and the nested rmM support. The
   * external packed BP words are not owned and are excluded.
   */
  std::size_t memory_usage_bytes() const {
    return sizeof(*this) + pixie::nested_owned_memory_bytes(rmm_);
  }

 private:
  using RmM = pixie::experimental::RmMBTree<HighCacheLines, LowFanout>;

  std::span<const std::uint64_t> words_;
  std::size_t depth_count_ = 0;
  RmM rmm_;
};

}  // namespace detail

/**
 * @brief Ferrada-Navarro Cartesian-tree RMQ using rmM support.
 *
 * @details This class follows the same public value-RMQ specification as the
 * other value RMQ backends, but replaces the usual balanced-parentheses
 * rank/select and depth-RMQ support with `detail::RmMPlusMinusOne`, which in
 * turn uses the experimental range min-max btree over the BP excess sequence.
 * BP construction uses a succinct monotone bit-stack, preserving the same
 * stable Cartesian-tree shape without an n-entry index stack.
 *
 * This implementation is included from `pixie/rmq.h` as the rmM-backed
 * Cartesian-tree reduction.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t HighCacheLines = 4,
          std::size_t LowFanout = 32>
class CartesianRmM
    : public RmqBase<CartesianRmM<T, Compare, Index, HighCacheLines, LowFanout>,
                     T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "CartesianRmM index type must be unsigned");

  static constexpr std::size_t npos =
      RmqBase<CartesianRmM<T, Compare, Index, HighCacheLines, LowFanout>,
              T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  CartesianRmM() = default;

  /**
   * @brief Build a Cartesian-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values stay stable: the smaller index remains the first minimum.
   */
  explicit CartesianRmM(std::span<const T> values, Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  CartesianRmM(const CartesianRmM& other)
      : values_(other.values_),
        compare_(other.compare_),
        bp_bits_(other.bp_bits_),
        bp_bit_count_(other.bp_bit_count_) {
    reset_bp_index();
  }

  CartesianRmM& operator=(const CartesianRmM& other) {
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

  CartesianRmM(CartesianRmM&& other) noexcept
      : values_(other.values_),
        compare_(std::move(other.compare_)),
        bp_bits_(std::move(other.bp_bits_)),
        bp_bit_count_(other.bp_bit_count_) {
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    reset_bp_index();
  }

  CartesianRmM& operator=(CartesianRmM&& other) noexcept {
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

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this value-RMQ object, packed Cartesian BP words, and
   * nested BP rmM support. The external input values are not owned and are
   * excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    return sizeof(*this) + pixie::vector_capacity_bytes(bp_bits_) +
           pixie::nested_owned_memory_bytes(bp_support_);
  }

 private:
  using BpSupport = detail::RmMPlusMinusOne<Index, HighCacheLines, LowFanout>;

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
    utils::SuccinctIncreasingStack stack(values_.size());
    std::size_t write_position = bp_bit_count_;

    for (std::size_t i = values_.size(); i-- > 0;) {
      while (!stack.empty() &&
             !compare_(values_[stack_index(stack.top())], values_[i])) {
        stack.pop();
        prepend_bp_bit(write_position, true);
      }
      stack.push(stack_key(i));
      prepend_bp_bit(write_position, false);
    }

    while (write_position != 0) {
      prepend_bp_bit(write_position, true);
    }
  }

  /**
   * @brief Convert a value position to the increasing key used by the
   * construction stack.
   */
  std::size_t stack_key(std::size_t value_index) const {
    return values_.size() - 1 - value_index;
  }

  /**
   * @brief Convert a construction-stack key back to the original value
   * position.
   */
  std::size_t stack_index(std::size_t key) const {
    return values_.size() - 1 - key;
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
    bp_support_ = BpSupport(std::span<const std::uint64_t>(bp_bits_),
                            bp_bit_count_ + 1, values_.size());
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<std::uint64_t> bp_bits_;
  std::size_t bp_bit_count_ = 0;
  BpSupport bp_support_;
};

}  // namespace pixie::rmq
