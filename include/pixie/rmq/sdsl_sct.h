#pragma once

#ifndef SDSL_SUPPORT
#error "pixie/rmq/sdsl_sct.h requires SDSL_SUPPORT"
#endif

#include <pixie/rmq/rmq_base.h>

#include <cstddef>
#include <functional>
#include <sdsl/io.hpp>
#include <sdsl/rmq_succinct_sct.hpp>
#include <span>
#include <type_traits>

namespace pixie::rmq {

/**
 * @brief Optional value-RMQ adapter over `sdsl::rmq_succinct_sct`.
 *
 * @details SDSL's SCT implementation builds a global Super-Cartesian-tree
 * balanced-parentheses sequence (BPS-SCT) over the indexed values. It then
 * answers inclusive range queries with balanced-parentheses navigation
 * operations such as `select`, `find_close`, `rr_enclose`, and `rank`.
 *
 * SDSL's helper comments attribute the Super-Cartesian-tree BP construction to
 * Ohlebusch and Gog in a compressed enhanced suffix-array context. That is not
 * the canonical RMQ citation. For the succinct RMQ result behind this style of
 * representation, see:
 *
 * Johannes Fischer, "Optimal Succinctness for Range Minimum Queries",
 * LATIN 2010; arXiv:0812.2775.
 *
 * This adapter exposes the Pixie value-RMQ contract: half-open ranges
 * `[left, right)`, invalid ranges returning `npos`, and first-minimum tie
 * semantics.
 *
 * The SDSL structure does not retain values after construction, so this adapter
 * keeps a non-owning span for `range_min()`. The indexed values must outlive
 * the adapter. SDSL SCT supports min/max through a boolean template parameter,
 * not arbitrary comparators; this wrapper currently exposes only
 * `std::less<T>`.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Must be `std::less<T>`.
 * @tparam Index Interface compatibility parameter for benchmark templates.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class SdslSct : public RmqBase<SdslSct<T, Compare, Index>, T> {
 public:
  static_assert(std::is_same_v<Compare, std::less<T>>,
                "SDSL SCT RMQ wrapper supports only std::less");

  static constexpr std::size_t npos =
      RmqBase<SdslSct<T, Compare, Index>, T>::npos;

  /**
   * @brief Construct an empty SDSL SCT RMQ adapter.
   */
  SdslSct() = default;

  /**
   * @brief Build the SDSL SCT index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller index as the RMQ answer.
   *
   * @param values Values to index.
   * @param compare Must be the default `std::less<T>` comparator.
   */
  explicit SdslSct(std::span<const T> values, Compare compare = Compare())
      : values_(values), rmq_(&values_) {
    (void)compare;
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
   *
   * @details SDSL's query is inclusive, so the adapter translates the Pixie
   * half-open range to `[left, right - 1]`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    return static_cast<std::size_t>(rmq_(left, right - 1));
  }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts the wrapper object plus SDSL's serialized SCT byte count.
   * The external input values are not owned and are excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    sdsl::nullstream out;
    return sizeof(*this) + static_cast<std::size_t>(rmq_.serialize(out));
  }

 private:
  std::span<const T> values_;
  sdsl::rmq_succinct_sct<true> rmq_;
};

}  // namespace pixie::rmq
