#pragma once

#include <pixie/bits.h>
#include <pixie/rmq/sparse_table.h>

#include <algorithm>
#include <array>
#include <bit>
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
 * @brief FCB-style RMQ backend for arrays with adjacent differences ±1.
 *
 * @details The indexed depth sequence is represented by BP deltas: bit 1 means
 * the next depth is current + 1, and bit 0 means current - 1. A sequence with
 * @p depth_count depth positions has @p depth_count - 1 delta bits. Blocks
 * match the 128-bit excess primitives in `bits.h`; only the absolute minimum
 * value of each block is stored, and positions are recovered by rescanning the
 * selected block.
 *
 * @tparam Index Unsigned integer type used for stored positions.
 * @tparam BlockSize Number of depth positions per microblock.
 */
template <class Index = std::size_t, std::size_t BlockSize = 128>
class BpPlusMinusOneRmq {
 public:
  static_assert(BlockSize == 128);

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  BpPlusMinusOneRmq() = default;

  BpPlusMinusOneRmq(std::span<const std::uint64_t> bits,
                    std::size_t depth_count)
      : input_bits_(bits), depth_count_(depth_count) {
    build();
  }

  BpPlusMinusOneRmq(const BpPlusMinusOneRmq& other)
      : input_bits_(other.input_bits_),
        depth_count_(other.depth_count_),
        block_min_values_(other.block_min_values_) {
    reset_macro_rmq();
  }

  BpPlusMinusOneRmq& operator=(const BpPlusMinusOneRmq& other) {
    if (this == &other) {
      return *this;
    }
    input_bits_ = other.input_bits_;
    depth_count_ = other.depth_count_;
    block_min_values_ = other.block_min_values_;
    reset_macro_rmq();
    return *this;
  }

  BpPlusMinusOneRmq(BpPlusMinusOneRmq&& other) noexcept
      : input_bits_(other.input_bits_),
        depth_count_(other.depth_count_),
        block_min_values_(std::move(other.block_min_values_)) {
    reset_macro_rmq();
  }

  BpPlusMinusOneRmq& operator=(BpPlusMinusOneRmq&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    input_bits_ = other.input_bits_;
    depth_count_ = other.depth_count_;
    block_min_values_ = std::move(other.block_min_values_);
    reset_macro_rmq();
    return *this;
  }

  std::size_t size() const { return depth_count_; }

  bool empty() const { return depth_count_ == 0; }

  std::size_t arg_min(std::size_t left, std::size_t right) const {
    if (left > right || right >= depth_count_) {
      return npos;
    }

    const std::size_t left_block = left / BlockSize;
    const std::size_t right_block = right / BlockSize;
    if (left_block == right_block) {
      return scan_block_range(left_block, left % BlockSize, right % BlockSize)
          .position;
    }

    Candidate answer = scan_block_range(left_block, left % BlockSize,
                                        block_size(left_block) - 1);
    answer =
        better(answer, scan_block_range(right_block, 0, right % BlockSize));

    if (left_block + 1 < right_block) {
      const std::size_t block_position =
          macro_rmq_.arg_min(left_block + 1, right_block - 1);
      if (block_position !=
          SparseTable<std::int64_t, std::less<std::int64_t>, Index>::npos) {
        answer = better(answer, scan_full_block(block_position));
      }
    }

    return answer.position;
  }

 private:
  struct Candidate {
    std::size_t position = npos;
    std::int64_t value = std::numeric_limits<std::int64_t>::max();
  };

  std::size_t block_size(std::size_t block) const {
    const std::size_t begin = block * BlockSize;
    return std::min(BlockSize, depth_count_ - begin);
  }

  Candidate better(Candidate left, Candidate right) const {
    if (left.position == npos) {
      return right;
    }
    if (right.position == npos) {
      return left;
    }
    if (right.value < left.value) {
      return right;
    }
    if (left.value < right.value) {
      return left;
    }
    return right.position < left.position ? right : left;
  }

  void build() {
    block_min_values_.clear();
    macro_rmq_ = SparseTable<std::int64_t, std::less<std::int64_t>, Index>();

    if (depth_count_ == 0) {
      return;
    }
    if (depth_count_ > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ ±1 index type is too small");
    }
    if (input_bits_.size() < (depth_count_ - 1 + 63) / 64) {
      throw std::invalid_argument("RMQ ±1 bit span is too small");
    }

    const std::size_t block_count = (depth_count_ + BlockSize - 1) / BlockSize;
    block_min_values_.reserve(block_count);

    std::int64_t base_depth = 0;
    for (std::size_t block = 0; block < block_count; ++block) {
      const std::size_t begin = block * BlockSize;
      const std::size_t size = std::min(BlockSize, depth_count_ - begin);
      std::int64_t min_depth = base_depth;
      std::int64_t current_depth = base_depth;
      for (std::size_t offset = 1; offset < size; ++offset) {
        const std::size_t delta_position = begin + offset - 1;
        const bool up = bit(delta_position);
        current_depth += up ? 1 : -1;
        if (current_depth < min_depth) {
          min_depth = current_depth;
        }
      }

      block_min_values_.push_back(min_depth);
      if (block + 1 < block_count) {
        base_depth += block_excess(begin, next_block_delta_count(begin));
      }
    }

    reset_macro_rmq();
  }

  bool bit(std::size_t position) const {
    return ((input_bits_[position >> 6] >> (position & 63)) & 1u) != 0;
  }

  std::uint64_t word_or_zero(std::size_t word) const {
    return word < input_bits_.size() ? input_bits_[word] : 0;
  }

  std::array<std::uint64_t, 2> block_bits(std::size_t block) const {
    const std::size_t first_word = block * (BlockSize / 64);
    return {word_or_zero(first_word), word_or_zero(first_word + 1)};
  }

  std::size_t next_block_delta_count(std::size_t begin) const {
    const std::size_t next_begin = begin + BlockSize;
    return std::min(next_begin, depth_count_ - 1) - begin;
  }

  std::int64_t block_excess(std::size_t begin, std::size_t delta_count) const {
    std::int64_t excess = 0;
    for (std::size_t i = 0; i < delta_count; ++i) {
      excess += bit(begin + i) ? 1 : -1;
    }
    return excess;
  }

  std::int64_t block_base_depth(std::size_t block,
                                const std::array<std::uint64_t, 2>& bits,
                                std::size_t size) const {
    const ExcessMin128Result full_min =
        excess_min_128(bits.data(), 0, size - 1);
    return block_min_values_[block] - full_min.min_excess;
  }

  Candidate scan_block_range(std::size_t block,
                             std::size_t left_offset,
                             std::size_t right_offset) const {
    const std::size_t begin = block * BlockSize;
    const std::size_t size = block_size(block);
    right_offset = std::min(right_offset, size - 1);
    const auto bits = block_bits(block);
    const std::int64_t base_depth = block_base_depth(block, bits, size);
    const ExcessMin128Result result =
        excess_min_128(bits.data(), left_offset, right_offset);
    if (result.offset == npos || result.offset >= size) {
      return {};
    }
    return {begin + result.offset, base_depth + result.min_excess};
  }

  Candidate scan_full_block(std::size_t block) const {
    const std::size_t begin = block * BlockSize;
    const std::size_t size = block_size(block);
    const auto bits = block_bits(block);
    const ExcessMin128Result result = excess_min_128(bits.data(), 0, size - 1);
    if (result.offset == npos || result.offset >= size) {
      return {};
    }
    return {begin + result.offset, block_min_values_[block]};
  }

  void reset_macro_rmq() {
    macro_rmq_ = SparseTable<std::int64_t, std::less<std::int64_t>, Index>(
        std::span<const std::int64_t>(block_min_values_));
  }

  std::span<const std::uint64_t> input_bits_;
  std::size_t depth_count_ = 0;
  std::vector<std::int64_t> block_min_values_;
  SparseTable<std::int64_t, std::less<std::int64_t>, Index> macro_rmq_;
};

}  // namespace pixie::rmq
