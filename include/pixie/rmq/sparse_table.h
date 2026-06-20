#pragma once

#include <pixie/cache_line.h>
#include <pixie/memory_usage.h>
#include <pixie/rmq/rmq_base.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
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
 * @tparam Alignment Byte alignment for each sparse-table level allocation.
 * Use 0 to select the standard allocator.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t Alignment = pixie::kAlignedStorageLineBytes>
class SparseTable
    : public RmqBase<SparseTable<T, Compare, Index, Alignment>, T> {
 public:
  static constexpr std::size_t npos =
      RmqBase<SparseTable<T, Compare, Index, Alignment>, T>::npos;
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
   * @brief Return the first minimum position in [@p left, @p right).
   *
   * @details Answers in O(1) by comparing the two power-of-two ranges covering
   * the half-open query interval. Ties return the smaller position.
   *
   * @param left First position in the query range.
   * @param right One past the last position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    const std::size_t length = right - left;
    const std::size_t level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << level;
    const std::size_t first = table_[level][left];
    const std::size_t second = table_[level][right - span];
    return better(first, second);
  }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this sparse-table object and all stored index levels.
   * The external input values are not owned and are excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    std::size_t bytes = sizeof(*this);
    bytes += pixie::vector_capacity_bytes(table_);
    for (const TableLevel& level : table_) {
      bytes += pixie::vector_capacity_bytes(level);
    }
    return bytes;
  }

 private:
  /**
   * @brief Choose the better of two candidate positions.
   *
   * @details Both candidates are valid table entries. If both values compare
   * equal, the smaller position wins to preserve first-minimum semantics.
   *
   * @param left First candidate position.
   * @param right Second candidate position.
   * @return Position of the selected candidate.
   */
  std::size_t better(std::size_t left, std::size_t right) const {
    if (compare_(values_[right], values_[left])) {
      return right;
    }
    if (compare_(values_[left], values_[right])) {
      return left;
    }
    return std::min(left, right);
  }

  /**
   * @brief Choose the better candidate while building adjacent ranges.
   *
   * @details Build ranges are disjoint and ordered, so ties always keep the
   * left candidate. This avoids the second strict comparison needed by the
   * general query-time helper.
   *
   * @param left Candidate from the left half.
   * @param right Candidate from the right half.
   * @return First-minimum candidate for the combined range.
   */
  Index build_better(Index left, Index right) const {
    return compare_(values_[right], values_[left]) ? right : left;
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

    table_.reserve(std::bit_width(values_.size()));
    table_.emplace_back(values_.size());
    std::iota(table_[0].begin(), table_[0].end(), Index{0});

    for (std::size_t span = 2, half = 1; span <= values_.size();
         half = span, span <<= 1) {
      const std::size_t level = table_.size();
      table_.emplace_back(values_.size() - span + 1);
      for (std::size_t i = 0; i < table_[level].size(); ++i) {
        table_[level][i] = static_cast<Index>(
            build_better(table_[level - 1][i], table_[level - 1][i + half]));
      }
    }
  }

  template <class Value, std::size_t AllocationAlignment>
  class AlignedAllocator {
   public:
    static_assert(AllocationAlignment >= alignof(Value));
    static_assert((AllocationAlignment & (AllocationAlignment - 1)) == 0);

    using value_type = Value;

    AlignedAllocator() = default;

    template <class Other>
    AlignedAllocator(
        const AlignedAllocator<Other, AllocationAlignment>&) noexcept {}

    [[nodiscard]] Value* allocate(std::size_t count) {
      if (count == 0) {
        return nullptr;
      }
      if (count > std::numeric_limits<std::size_t>::max() / sizeof(Value)) {
        throw std::bad_array_new_length();
      }
      return static_cast<Value*>(::operator new(
          count * sizeof(Value), std::align_val_t{AllocationAlignment}));
    }

    void deallocate(Value* pointer, std::size_t) noexcept {
      ::operator delete(pointer, std::align_val_t{AllocationAlignment});
    }

    template <class Other>
    bool operator==(
        const AlignedAllocator<Other, AllocationAlignment>&) const noexcept {
      return true;
    }

    template <class Other>
    bool operator!=(
        const AlignedAllocator<Other, AllocationAlignment>&) const noexcept {
      return false;
    }

    template <class Other>
    struct rebind {
      using other = AlignedAllocator<Other, AllocationAlignment>;
    };
  };

  template <class Value, std::size_t AllocationAlignment>
  struct LevelAllocator {
    using type = AlignedAllocator<Value, AllocationAlignment>;
  };

  template <class Value>
  struct LevelAllocator<Value, 0> {
    using type = std::allocator<Value>;
  };

  using TableLevel =
      std::vector<Index, typename LevelAllocator<Index, Alignment>::type>;

  std::span<const T> values_;
  Compare compare_;
  std::vector<TableLevel> table_;
};

}  // namespace pixie::rmq
