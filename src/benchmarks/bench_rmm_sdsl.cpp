#include <pixie/rmm_base.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "rmm_benchmark_base.h"

namespace {

class SdslRmMTree : public pixie::RmMBase<SdslRmMTree> {
 public:
  using BpSupport = sdsl::bp_support_sada<>;
  static constexpr std::size_t npos = pixie::RmMBase<SdslRmMTree>::npos;

  SdslRmMTree() = default;

  SdslRmMTree(std::span<const std::uint64_t> words,
              std::size_t bit_count,
              std::size_t) {
    size_ = bit_count;
    for (std::size_t i = 0; i < words.size(); ++i) {
      std::uint64_t word = words[i];
      if (i + 1 == words.size() && (bit_count & 63) != 0) {
        word &= (std::uint64_t(1) << (bit_count & 63)) - 1;
      }
      ones_ += std::popcount(word);
    }
    zeros_ = size_ - ones_;

    bits_ = sdsl::bit_vector(size_);
    for (std::size_t i = 0; i < size_; ++i) {
      bits_[i] = (words[i >> 6] >> (i & 63)) & 1u;
    }
    tree_ = BpSupport(&bits_);
  }

  SdslRmMTree(const SdslRmMTree& other)
      : size_(other.size_),
        ones_(other.ones_),
        zeros_(other.zeros_),
        bits_(other.bits_) {
    tree_ = BpSupport(&bits_);
  }

  SdslRmMTree& operator=(const SdslRmMTree& other) {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    ones_ = other.ones_;
    zeros_ = other.zeros_;
    bits_ = other.bits_;
    tree_ = BpSupport(&bits_);
    return *this;
  }

  SdslRmMTree(SdslRmMTree&& other) noexcept
      : size_(other.size_),
        ones_(other.ones_),
        zeros_(other.zeros_),
        bits_(std::move(other.bits_)) {
    tree_ = BpSupport(&bits_);
  }

  SdslRmMTree& operator=(SdslRmMTree&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    ones_ = other.ones_;
    zeros_ = other.zeros_;
    bits_ = std::move(other.bits_);
    tree_ = BpSupport(&bits_);
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

  std::size_t fwdsearch_impl(std::size_t, int) const { return npos; }
  std::size_t bwdsearch_impl(std::size_t, int) const { return npos; }

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
    return tree_.find_close(open_position);
  }

  std::size_t open_impl(std::size_t close_position_one_based) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t close_position_zero_based =
        close_position_one_based > 0 ? close_position_one_based - 1 : 0;
    return tree_.find_open(close_position_zero_based);
  }

  std::size_t enclose_impl(std::size_t open_position_one_based) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t open_position_zero_based =
        open_position_one_based > 0 ? open_position_one_based - 1 : 0;
    return tree_.enclose(open_position_zero_based);
  }

 private:
  std::size_t size_{};
  std::size_t ones_{};
  std::size_t zeros_{};
  sdsl::bit_vector bits_;
  BpSupport tree_;
};

}  // namespace

namespace pixie_bench {

template <>
struct RmMBenchmarkTraits<SdslRmMTree> {
  static bool SupportsOp(std::string_view op) {
    return op == "rank1" || op == "rank0" || op == "select1" ||
           op == "excess" || op == "range_min_query_pos" ||
           op == "range_min_query_val" || op == "close" || op == "open" ||
           op == "enclose";
  }
};

}  // namespace pixie_bench

int main(int argc, char** argv) {
  pixie_bench::RmMBenchmark<SdslRmMTree> benchmark;
  return benchmark.Run(argc, argv);
}
