#pragma once

/**
 * @file split_span.h
 * @brief Allocation-free logical ranges split across at most two spans.
 */

#include <array>
#include <cstddef>
#include <span>

namespace pixie {

/**
 * @brief A logical contiguous sequence stored in at most two physical spans.
 *
 * @details Empty input spans are canonicalized away. Iteration visits only
 * non-empty physical spans in logical order. The descriptor owns no elements;
 * its spans follow the lifetime and invalidation rules of their backing
 * storage.
 *
 * @tparam T Viewed element type, optionally const-qualified.
 */
template <class T>
class SplitSpan {
 public:
  using span_type = std::span<T>;
  using const_iterator = typename std::array<span_type, 2>::const_iterator;

  /** @brief Construct an empty logical range. */
  constexpr SplitSpan() = default;

  /** @brief Construct a range stored in one physical span. */
  constexpr explicit SplitSpan(span_type segment)
      : segments_{segment, {}}, segment_count_(segment.empty() ? 0 : 1) {}

  /** @brief Construct a range stored in up to two physical spans. */
  constexpr SplitSpan(span_type first, span_type second) {
    if (first.empty()) {
      first = second;
      second = {};
    }
    segments_ = {first, second};
    segment_count_ = static_cast<std::size_t>(!first.empty()) +
                     static_cast<std::size_t>(!second.empty());
  }

  /** @brief Return the total number of logical elements. */
  constexpr std::size_t size() const noexcept {
    return segments_[0].size() + segments_[1].size();
  }

  /** @brief Return whether the logical range is empty. */
  constexpr bool empty() const noexcept { return segment_count_ == 0; }

  /** @brief Return the number of non-empty physical segments. */
  constexpr std::size_t segment_count() const noexcept {
    return segment_count_;
  }

  /** @brief Iterate over the non-empty physical segments in logical order. */
  constexpr const_iterator begin() const noexcept { return segments_.begin(); }

  /** @brief Return the end iterator for the non-empty physical segments. */
  constexpr const_iterator end() const noexcept {
    return segments_.begin() + static_cast<std::ptrdiff_t>(segment_count_);
  }

 private:
  std::array<span_type, 2> segments_{};
  std::size_t segment_count_ = 0;
};

}  // namespace pixie
