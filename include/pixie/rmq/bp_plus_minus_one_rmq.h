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
 * match the 128-bit excess primitives in `bits.h`; each block stores a compact
 * 16-byte summary with its base depth, absolute minimum value, and first local
 * offset attaining that minimum.
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

  /**
   * @brief Construct an empty ±1 RMQ index.
   */
  BpPlusMinusOneRmq() = default;

  /**
   * @brief Build a ±1 RMQ index over a BP-encoded depth-delta sequence.
   *
   * @details The bit span is not copied and must outlive this object. Bit `1`
   * encodes a +1 step between adjacent depths, and bit `0` encodes a -1 step.
   * A sequence with @p depth_count positions requires @p depth_count - 1 bits.
   *
   * @param bits Packed delta bits in little-endian bit order within each word.
   * @param depth_count Number of depth positions indexed by the RMQ.
   * @throws std::length_error if @p Index cannot represent all positions or a
   * packed 48-bit block summary overflows.
   * @throws std::invalid_argument if @p bits does not contain enough words.
   */
  BpPlusMinusOneRmq(std::span<const std::uint64_t> bits,
                    std::size_t depth_count)
      : input_bits_(bits), depth_count_(depth_count) {
    build();
  }

  /**
   * @brief Copy a ±1 RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   */
  BpPlusMinusOneRmq(const BpPlusMinusOneRmq& other)
      : input_bits_(other.input_bits_),
        depth_count_(other.depth_count_),
        block_summaries_(other.block_summaries_) {
    reset_macro_rmq();
  }

  /**
   * @brief Copy-assign a ±1 RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   * @return Reference to this object.
   */
  BpPlusMinusOneRmq& operator=(const BpPlusMinusOneRmq& other) {
    if (this == &other) {
      return *this;
    }
    input_bits_ = other.input_bits_;
    depth_count_ = other.depth_count_;
    block_summaries_ = other.block_summaries_;
    reset_macro_rmq();
    return *this;
  }

  /**
   * @brief Move a ±1 RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   */
  BpPlusMinusOneRmq(BpPlusMinusOneRmq&& other) noexcept
      : input_bits_(other.input_bits_),
        depth_count_(other.depth_count_),
        block_summaries_(std::move(other.block_summaries_)) {
    reset_macro_rmq();
  }

  /**
   * @brief Move-assign a ±1 RMQ index and rebuild internal non-owning views.
   *
   * @param other Source index.
   * @return Reference to this object.
   */
  BpPlusMinusOneRmq& operator=(BpPlusMinusOneRmq&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    input_bits_ = other.input_bits_;
    depth_count_ = other.depth_count_;
    block_summaries_ = std::move(other.block_summaries_);
    reset_macro_rmq();
    return *this;
  }

  /**
   * @brief Return the number of indexed depth positions.
   *
   * @return The @p depth_count passed to construction.
   */
  std::size_t size() const { return depth_count_; }

  /**
   * @brief Whether the indexed depth sequence is empty.
   *
   * @return `true` when `size() == 0`.
   */
  bool empty() const { return depth_count_ == 0; }

  /**
   * @brief Return the first minimum depth position in [@p left, @p right].
   *
   * @details The query range is inclusive over depth positions, not delta-bit
   * positions. Ties return the smallest position attaining the minimum depth.
   *
   * @param left First depth position in the query range.
   * @param right Last depth position in the query range.
   * @return Zero-based position of the first range minimum, or `npos`.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    if (left > right || right >= depth_count_) {
      return npos;
    }

    const std::size_t left_block = left / BlockSize;
    const std::size_t right_block = right / BlockSize;
    if (left_block == right_block) {
      return scan_block_range_position(left_block, left % BlockSize,
                                       right % BlockSize);
    }

    const std::size_t left_offset = left % BlockSize;
    const std::size_t right_offset = right % BlockSize;
    Candidate answer =
        scan_block_range(left_block, left_offset, block_size(left_block) - 1);
    answer = better(answer, scan_block_range(right_block, 0, right_offset));

    if (left_block + 1 < right_block) {
      const std::size_t block_position =
          macro_rmq_.arg_min(left_block + 1, right_block - 1);
      if (block_position != MacroRmq::npos) {
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

  /**
   * @brief Packed summary for one 128-position depth block.
   *
   * @details Stores signed 48-bit base depth, signed 48-bit absolute block
   * minimum, and 16-bit first local offset of that minimum in two 64-bit words.
   * The hot macro-RMQ comparison keeps the sign-biased minimum in the low
   * 48 bits of `word0`, so signed ordering is available as one masked unsigned
   * comparison. The high 16 bits of `word0` store the local minimum offset, and
   * `word1` stores the base depth.
   */
  struct alignas(16) BlockSummary {
    static constexpr std::uint64_t kSigned48Mask = (std::uint64_t{1} << 48) - 1;
    static constexpr std::uint64_t kSigned48SignBit = std::uint64_t{1} << 47;
    static constexpr std::int64_t kSigned48Min = -(std::int64_t{1} << 47);
    static constexpr std::int64_t kSigned48Max = (std::int64_t{1} << 47) - 1;

    std::uint64_t word0 = 0;
    std::uint64_t word1 = 0;

    /**
     * @brief Pack block metadata into a 16-byte summary.
     *
     * @param base_depth Absolute depth at local offset 0.
     * @param min_value Absolute minimum depth in the block.
     * @param min_offset First local offset attaining @p min_value.
     * @return Packed block summary.
     * @throws std::length_error if a signed 48-bit field or 16-bit offset
     * overflows.
     */
    static BlockSummary make(std::int64_t base_depth,
                             std::int64_t min_value,
                             std::size_t min_offset) {
      if (!fits_signed48(base_depth) || !fits_signed48(min_value)) {
        throw std::length_error("RMQ +/-1 block summary depth overflow");
      }
      if (min_offset > std::numeric_limits<std::uint16_t>::max()) {
        throw std::length_error("RMQ +/-1 block summary offset overflow");
      }

      const std::uint64_t packed_base = pack_signed48(base_depth);
      const std::uint64_t ordered_min = pack_ordered48(min_value);
      return {ordered_min | (static_cast<std::uint64_t>(min_offset) << 48),
              packed_base};
    }

    /**
     * @brief Decode the absolute depth at local offset 0.
     *
     * @return Signed base depth for this block.
     */
    std::int64_t base_depth() const {
      return unpack_signed48(word1 & kSigned48Mask);
    }

    /**
     * @brief Decode the absolute minimum depth in this block.
     *
     * @return Signed block minimum value.
     */
    std::int64_t min_value() const {
      return unpack_ordered48(ordered_min_value());
    }

    /**
     * @brief Return the sign-biased minimum payload for fast comparisons.
     *
     * @details The returned low-48-bit value preserves signed order under
     * unsigned comparison: smaller signed depths have smaller payloads.
     *
     * @return Sign-biased 48-bit minimum depth payload.
     */
    std::uint64_t ordered_min_value() const { return word0 & kSigned48Mask; }

    /**
     * @brief Decode the first local minimum offset in this block.
     *
     * @return Local depth-position offset in [0, BlockSize).
     */
    std::size_t min_offset() const {
      return static_cast<std::size_t>(word0 >> 48);
    }

   private:
    /**
     * @brief Test whether a value fits in the signed 48-bit packed field.
     *
     * @param value Value to test.
     * @return `true` if @p value is representable.
     */
    static bool fits_signed48(std::int64_t value) {
      return value >= kSigned48Min && value <= kSigned48Max;
    }

    /**
     * @brief Truncate a signed value to its 48-bit two's-complement payload.
     *
     * @param value Signed 48-bit value.
     * @return Low 48 payload bits.
     */
    static std::uint64_t pack_signed48(std::int64_t value) {
      return static_cast<std::uint64_t>(value) & kSigned48Mask;
    }

    /**
     * @brief Encode a signed 48-bit value for unsigned ordered comparison.
     *
     * @param value Signed 48-bit value.
     * @return Low 48 payload bits with the sign bit flipped.
     */
    static std::uint64_t pack_ordered48(std::int64_t value) {
      return pack_signed48(value) ^ kSigned48SignBit;
    }

    /**
     * @brief Sign-extend a 48-bit two's-complement payload.
     *
     * @param value Low 48 payload bits.
     * @return Decoded signed value.
     */
    static std::int64_t unpack_signed48(std::uint64_t value) {
      if ((value & kSigned48SignBit) != 0) {
        value |= ~kSigned48Mask;
      }
      return static_cast<std::int64_t>(value);
    }

    /**
     * @brief Decode a sign-biased 48-bit ordered payload.
     *
     * @param value Low 48 payload bits with the sign bit flipped.
     * @return Decoded signed value.
     */
    static std::int64_t unpack_ordered48(std::uint64_t value) {
      return unpack_signed48(value ^ kSigned48SignBit);
    }
  };

  static_assert(sizeof(BlockSummary) == 16);
  static_assert(alignof(BlockSummary) == 16);

  struct BlockSummaryMinLess {
    /**
     * @brief Compare two block summaries by absolute minimum value.
     *
     * @param left First block summary.
     * @param right Second block summary.
     * @return `true` if @p left has a strictly smaller block minimum.
     */
    bool operator()(const BlockSummary& left, const BlockSummary& right) const {
      return left.ordered_min_value() < right.ordered_min_value();
    }
  };

  using MacroRmq = SparseTable<BlockSummary, BlockSummaryMinLess, Index>;

  /**
   * @brief Return the number of depth positions in a block.
   *
   * @param block Zero-based block index.
   * @return `BlockSize` for full blocks, or the tail size for the last block.
   */
  std::size_t block_size(std::size_t block) const {
    const std::size_t begin = block * BlockSize;
    return std::min(BlockSize, depth_count_ - begin);
  }

  /**
   * @brief Choose the better of two candidate minima.
   *
   * @details Missing candidates use `npos`. Smaller values win; equal values
   * return the smaller global position.
   *
   * @param left First candidate.
   * @param right Second candidate.
   * @return Selected candidate.
   */
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

  /**
   * @brief Build block summaries and the macro sparse table.
   *
   * @details Computes the absolute base depth, absolute block minimum, and
   * first local minimum offset for each 128-position block.
   *
   * @throws std::length_error if @p Index or a packed block summary overflows.
   * @throws std::invalid_argument if the input bit span is too small.
   */
  void build() {
    block_summaries_.clear();
    macro_rmq_ = MacroRmq();

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
    block_summaries_.reserve(block_count);

    std::int64_t base_depth = 0;
    for (std::size_t block = 0; block < block_count; ++block) {
      const std::size_t begin = block * BlockSize;
      const std::size_t size = std::min(BlockSize, depth_count_ - begin);
      std::size_t min_offset = 0;
      std::int64_t min_depth = base_depth;
      std::int64_t current_depth = base_depth;
      for (std::size_t offset = 1; offset < size; ++offset) {
        const std::size_t delta_position = begin + offset - 1;
        const bool up = bit(delta_position);
        current_depth += up ? 1 : -1;
        if (current_depth < min_depth) {
          min_depth = current_depth;
          min_offset = offset;
        }
      }

      block_summaries_.push_back(
          BlockSummary::make(base_depth, min_depth, min_offset));
      if (block + 1 < block_count) {
        base_depth += block_excess(begin, next_block_delta_count(begin));
      }
    }

    reset_macro_rmq();
  }

  /**
   * @brief Read one BP delta bit.
   *
   * @param position Zero-based delta-bit position.
   * @return `true` for a +1 step and `false` for a -1 step.
   */
  bool bit(std::size_t position) const {
    return ((input_bits_[position >> 6] >> (position & 63)) & 1u) != 0;
  }

  /**
   * @brief Read an input word or return zero past the available span.
   *
   * @param word Zero-based 64-bit word index.
   * @return Input word value, or zero when @p word is out of range.
   */
  std::uint64_t word_or_zero(std::size_t word) const {
    return word < input_bits_.size() ? input_bits_[word] : 0;
  }

  /**
   * @brief Load the two 64-bit words backing a 128-position block.
   *
   * @param block Zero-based block index.
   * @return Pair of words suitable for `excess_min_128`.
   */
  std::array<std::uint64_t, 2> block_bits(std::size_t block) const {
    const std::size_t first_word = block * (BlockSize / 64);
    return {word_or_zero(first_word), word_or_zero(first_word + 1)};
  }

  /**
   * @brief Count delta bits from a block start to the next block boundary.
   *
   * @param begin Global depth-position offset of the block start.
   * @return Number of delta bits contributing to the transition to the next
   * block base depth.
   */
  std::size_t next_block_delta_count(std::size_t begin) const {
    const std::size_t next_begin = begin + BlockSize;
    return std::min(next_begin, depth_count_ - 1) - begin;
  }

  /**
   * @brief Compute total excess over a contiguous delta-bit range.
   *
   * @param begin First delta-bit position.
   * @param delta_count Number of delta bits to scan.
   * @return Sum of +1/-1 deltas in the range.
   */
  std::int64_t block_excess(std::size_t begin, std::size_t delta_count) const {
    std::int64_t excess = 0;
    for (std::size_t i = 0; i < delta_count; ++i) {
      excess += bit(begin + i) ? 1 : -1;
    }
    return excess;
  }

  /**
   * @brief Return only the minimum position inside one block range.
   *
   * @details Used for same-block queries where the absolute minimum value is
   * not needed. The range is inclusive in local block offsets.
   *
   * @param block Zero-based block index.
   * @param left_offset First local depth-position offset.
   * @param right_offset Last local depth-position offset.
   * @return Global position of the first local minimum, or `npos`.
   */
  std::size_t scan_block_range_position(std::size_t block,
                                        std::size_t left_offset,
                                        std::size_t right_offset) const {
    const std::size_t begin = block * BlockSize;
    const std::size_t size = block_size(block);
    right_offset = std::min(right_offset, size - 1);
    const auto bits = block_bits(block);
    const ExcessResult result =
        excess_min_128(bits.data(), left_offset, right_offset);
    if (result.offset == npos || result.offset >= size) {
      return npos;
    }
    return begin + result.offset;
  }

  /**
   * @brief Scan an inclusive local range inside one block.
   *
   * @details Uses `excess_min_128` for the relative minimum and combines it
   * with the stored block base depth to return an absolute-depth candidate.
   *
   * @param block Zero-based block index.
   * @param left_offset First local depth-position offset.
   * @param right_offset Last local depth-position offset.
   * @return Candidate containing global position and absolute depth.
   */
  Candidate scan_block_range(std::size_t block,
                             std::size_t left_offset,
                             std::size_t right_offset) const {
    const std::size_t begin = block * BlockSize;
    const std::size_t size = block_size(block);
    right_offset = std::min(right_offset, size - 1);
    const auto bits = block_bits(block);
    const ExcessResult result =
        excess_min_128(bits.data(), left_offset, right_offset);
    if (result.offset == npos || result.offset >= size) {
      return {};
    }
    return {begin + result.offset,
            block_summaries_[block].base_depth() + result.min_excess};
  }

  /**
   * @brief Return the precomputed minimum candidate for a full block.
   *
   * @param block Zero-based block index.
   * @return Candidate containing global position and absolute depth.
   */
  Candidate scan_full_block(std::size_t block) const {
    const std::size_t begin = block * BlockSize;
    return {begin + block_summaries_[block].min_offset(),
            block_summaries_[block].min_value()};
  }

  /**
   * @brief Rebuild the macro sparse table over block summaries.
   *
   * @details Called after build, copy, and move operations because the sparse
   * table stores a non-owning span into this object's block-summary vector.
   */
  void reset_macro_rmq() {
    macro_rmq_ = MacroRmq(std::span<const BlockSummary>(block_summaries_));
  }

  std::span<const std::uint64_t> input_bits_;
  std::size_t depth_count_ = 0;
  std::vector<BlockSummary> block_summaries_;
  MacroRmq macro_rmq_;
};

}  // namespace pixie::rmq
