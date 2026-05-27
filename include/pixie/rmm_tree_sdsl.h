#pragma once

#ifndef SDSL_SUPPORT
#error "pixie/rmm_tree_sdsl.h requires SDSL_SUPPORT"
#endif

#include <pixie/rmm_base.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <sdsl/bit_vectors.hpp>

// SDSL keeps the generic excess-search primitives private and exposes only
// navigation wrappers such as find_close/find_open. The Pixie comparison
// backend needs direct fwdsearch/bwdsearch, so expose them in this optional
// adapter instead of benchmarking a naive fallback.
#define private public
#include <sdsl/bp_support_sada.hpp>
#undef private

#include <span>
#include <utility>
#include <vector>

namespace pixie {

class SdslRmMTree : public RmMBase<SdslRmMTree> {
 public:
  using BpSupport = sdsl::bp_support_sada<>;
  static constexpr std::size_t npos = RmMBase<SdslRmMTree>::npos;

  SdslRmMTree() = default;

  SdslRmMTree(std::span<const std::uint64_t> words,
              std::size_t bit_count,
              std::size_t) {
    size_ = bit_count;
    const std::size_t valid_words = (bit_count + 63) / 64;
    for (std::size_t i = 0; i < valid_words; ++i) {
      std::uint64_t word = words[i];
      if (i + 1 == valid_words && (bit_count & 63) != 0) {
        word &= (std::uint64_t{1} << (bit_count & 63)) - 1;
      }
      ones_ += std::popcount(word);
    }
    zeros_ = size_ - ones_;

    bits_ = sdsl::bit_vector(size_);
    prefix_excess_.assign(size_ + 1, 0);
    int current_excess = 0;
    for (std::size_t i = 0; i < size_; ++i) {
      const bool bit = (words[i >> 6] >> (i & 63)) & 1u;
      bits_[i] = bit;
      current_excess += bit ? 1 : -1;
      prefix_excess_[i + 1] = current_excess;
      max_excess_ = std::max(max_excess_, current_excess);
    }
    build_excess_bounds();
    tree_ = BpSupport(&bits_);
  }

  SdslRmMTree(const SdslRmMTree& other)
      : size_(other.size_),
        ones_(other.ones_),
        zeros_(other.zeros_),
        max_excess_(other.max_excess_),
        prefix_excess_(other.prefix_excess_),
        prefix_min_excess_(other.prefix_min_excess_),
        prefix_max_excess_(other.prefix_max_excess_),
        suffix_min_excess_(other.suffix_min_excess_),
        suffix_max_excess_(other.suffix_max_excess_),
        bits_(other.bits_) {
    reset_support();
  }

  SdslRmMTree& operator=(const SdslRmMTree& other) {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    ones_ = other.ones_;
    zeros_ = other.zeros_;
    max_excess_ = other.max_excess_;
    prefix_excess_ = other.prefix_excess_;
    prefix_min_excess_ = other.prefix_min_excess_;
    prefix_max_excess_ = other.prefix_max_excess_;
    suffix_min_excess_ = other.suffix_min_excess_;
    suffix_max_excess_ = other.suffix_max_excess_;
    bits_ = other.bits_;
    reset_support();
    return *this;
  }

  SdslRmMTree(SdslRmMTree&& other) noexcept
      : size_(other.size_),
        ones_(other.ones_),
        zeros_(other.zeros_),
        max_excess_(other.max_excess_),
        prefix_excess_(std::move(other.prefix_excess_)),
        prefix_min_excess_(std::move(other.prefix_min_excess_)),
        prefix_max_excess_(std::move(other.prefix_max_excess_)),
        suffix_min_excess_(std::move(other.suffix_min_excess_)),
        suffix_max_excess_(std::move(other.suffix_max_excess_)),
        bits_(std::move(other.bits_)) {
    reset_support();
  }

  SdslRmMTree& operator=(SdslRmMTree&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    ones_ = other.ones_;
    zeros_ = other.zeros_;
    max_excess_ = other.max_excess_;
    prefix_excess_ = std::move(other.prefix_excess_);
    prefix_min_excess_ = std::move(other.prefix_min_excess_);
    prefix_max_excess_ = std::move(other.prefix_max_excess_);
    suffix_min_excess_ = std::move(other.suffix_min_excess_);
    suffix_max_excess_ = std::move(other.suffix_max_excess_);
    bits_ = std::move(other.bits_);
    reset_support();
    return *this;
  }

  std::size_t size_impl() const { return size_; }

  std::size_t rank1_impl(std::size_t end_position) const {
    if (size_ == 0 || end_position == 0) {
      return 0;
    }
    return tree_.rank(std::min(end_position, size_) - 1);
  }

  std::size_t rank0_impl(std::size_t end_position) const {
    if (size_ == 0 || end_position == 0) {
      return 0;
    }
    if (end_position >= size_) {
      return zeros_;
    }
    return tree_.preceding_closing_parentheses(end_position);
  }

  std::size_t select1_impl(std::size_t rank) const {
    if (rank == 0 || rank > ones_) {
      return npos;
    }
    return tree_.select(rank);
  }

  std::size_t select0_impl(std::size_t) const { return npos; }
  std::size_t rank10_impl(std::size_t) const { return 0; }
  std::size_t select10_impl(std::size_t) const { return npos; }

  int excess_impl(std::size_t end_position) const {
    if (size_ == 0 || end_position == 0) {
      return 0;
    }
    return tree_.excess(std::min(end_position, size_) - 1);
  }

  std::size_t fwdsearch_impl(std::size_t start_position, int delta) const {
    if (start_position >= size_) {
      return npos;
    }
    const int target = prefix_excess_[start_position] + delta;
    if (target > max_excess_) {
      return npos;
    }

    if (start_position == 0) {
      const int first_excess = bits_[0] ? 1 : -1;
      if (first_excess == delta) {
        return 0;
      }
      if (!suffix_contains(2, target)) {
        return npos;
      }
      const std::size_t position =
          tree_.fwd_excess(0, static_cast<typename BpSupport::difference_type>(
                                  delta - first_excess));
      return position < size_ ? position : npos;
    }

    if (!suffix_contains(start_position + 1, target)) {
      return npos;
    }
    const std::size_t position = tree_.fwd_excess(
        start_position - 1,
        static_cast<typename BpSupport::difference_type>(delta));
    return position < size_ ? position : npos;
  }

  std::size_t bwdsearch_impl(std::size_t start_position, int delta) const {
    if (start_position == 0 || start_position > size_) {
      return npos;
    }

    const std::size_t anchor = start_position - 1;
    const int target = prefix_excess_[start_position] + delta;
    if (target > max_excess_) {
      return npos;
    }
    if (prefix_excess_[anchor] == target) {
      return anchor;
    }
    if (anchor == 0) {
      return npos;
    }
    if (!prefix_contains(anchor - 1, target)) {
      return npos;
    }

    const std::size_t position = tree_.bwd_excess(
        anchor, static_cast<typename BpSupport::difference_type>(delta));
    if (position == static_cast<std::size_t>(-1)) {
      return target == 0 ? 0 : npos;
    }
    return position < size_ ? position + 1 : npos;
  }

  std::size_t range_min_query_pos_impl(std::size_t range_begin,
                                       std::size_t range_end) const {
    if (size_ == 0) {
      return 0;
    }
    return tree_.rmq(range_begin, range_end);
  }

  int range_min_query_val_impl(std::size_t range_begin,
                               std::size_t range_end) const {
    if (size_ == 0) {
      return 0;
    }
    const auto min_position = tree_.rmq(range_begin, range_end);
    if (min_position >= size_) {
      return 0;
    }
    const auto base_excess =
        (range_begin == 0 ? 0 : tree_.excess(range_begin - 1));
    return tree_.excess(min_position) - base_excess;
  }

  std::size_t range_max_query_pos_impl(std::size_t, std::size_t) const {
    return npos;
  }
  int range_max_query_val_impl(std::size_t, std::size_t) const { return 0; }
  std::size_t mincount_impl(std::size_t, std::size_t) const { return 0; }
  std::size_t minselect_impl(std::size_t, std::size_t, std::size_t) const {
    return npos;
  }

  std::size_t close_impl(std::size_t open_position) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t position = tree_.find_close(open_position);
    return position < size_ ? position : npos;
  }

  std::size_t open_impl(std::size_t close_position) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t position = tree_.find_open(close_position);
    return position < size_ ? position : npos;
  }

  std::size_t enclose_impl(std::size_t open_position) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t position = tree_.enclose(open_position);
    return position < size_ ? position : npos;
  }

  int bit_impl(const size_t& position) const noexcept {
    return (bits_[position >> 6] >> (position & 63)) & 1u;
  }

 private:
  void reset_support() { tree_ = BpSupport(&bits_); }

  void build_excess_bounds() {
    prefix_min_excess_.resize(size_ + 1);
    prefix_max_excess_.resize(size_ + 1);
    suffix_min_excess_.resize(size_ + 1);
    suffix_max_excess_.resize(size_ + 1);

    prefix_min_excess_[0] = prefix_excess_[0];
    prefix_max_excess_[0] = prefix_excess_[0];
    for (std::size_t i = 1; i <= size_; ++i) {
      prefix_min_excess_[i] =
          std::min(prefix_min_excess_[i - 1], prefix_excess_[i]);
      prefix_max_excess_[i] =
          std::max(prefix_max_excess_[i - 1], prefix_excess_[i]);
    }

    suffix_min_excess_[size_] = prefix_excess_[size_];
    suffix_max_excess_[size_] = prefix_excess_[size_];
    for (std::size_t i = size_; i > 0;) {
      --i;
      suffix_min_excess_[i] =
          std::min(prefix_excess_[i], suffix_min_excess_[i + 1]);
      suffix_max_excess_[i] =
          std::max(prefix_excess_[i], suffix_max_excess_[i + 1]);
    }
  }

  bool suffix_contains(std::size_t boundary_begin, int target) const {
    return boundary_begin <= size_ &&
           suffix_min_excess_[boundary_begin] <= target &&
           target <= suffix_max_excess_[boundary_begin];
  }

  bool prefix_contains(std::size_t boundary_end, int target) const {
    return prefix_min_excess_[boundary_end] <= target &&
           target <= prefix_max_excess_[boundary_end];
  }

  std::size_t size_{};
  std::size_t ones_{};
  std::size_t zeros_{};
  int max_excess_{};
  std::vector<int> prefix_excess_;
  std::vector<int> prefix_min_excess_;
  std::vector<int> prefix_max_excess_;
  std::vector<int> suffix_min_excess_;
  std::vector<int> suffix_max_excess_;
  sdsl::bit_vector bits_;
  BpSupport tree_;
};

}  // namespace pixie
