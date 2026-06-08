#pragma once

#include <pixie/bits.h>
#include <pixie/rmq/rmq_base.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::rmq::experimental {

/**
 * @brief Low-level selector tag for BP-based 252-value leaves.
 */
struct BpLeafSelectorTag {};

/**
 * @brief Experimental value RMQ with compact per-level BP selectors.
 *
 * @details Values are split into 252-value leaves. Low nodes store 504 BP bits
 * plus an 8-bit local-minimum offset in one 512-bit selector. Middle internal
 * nodes use 192-way fanout and store 384 BP bits, a 64-bit absolute
 * subtree-minimum position, and 64 bits of zero-rank prefix metadata in the
 * selector. High nodes use 256-way fanout and are only the root and its child
 * level, so there are at most 257 of them. High nodes add a sparse depth RMQ
 * side table and cache each child slot's value range and subtree minimum.
 *
 * This backend is intentionally not included by `pixie/rmq.h`; include this
 * header directly while evaluating the experiment.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 * @tparam LeafSize Number of original values per leaf. This experimental
 * backend currently supports only 252.
 * @tparam Fanout Number of children per high internal node. This experimental
 * backend currently supports only 256. Middle internal nodes use 192.
 * @tparam LowLevelSelector Compile-time low-level selector tag. The default
 * BP leaf selector is the only implemented tag in this cleanup pass.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 252,
          std::size_t Fanout = 256,
          class LowLevelSelector = BpLeafSelectorTag>
class NodeEulerBTreeRmq : public RmqBase<NodeEulerBTreeRmq<T,
                                                           Compare,
                                                           Index,
                                                           LeafSize,
                                                           Fanout,
                                                           LowLevelSelector>,
                                         T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "NodeEulerBTreeRmq index type must be unsigned");
  static_assert(LeafSize == 252,
                "experimental NodeEulerBTreeRmq currently requires 252-value "
                "leaves");
  static_assert(Fanout == 256,
                "experimental NodeEulerBTreeRmq currently requires 256-way "
                "internal nodes");
  static_assert(std::is_same_v<LowLevelSelector, BpLeafSelectorTag>,
                "only BP low-level leaves are implemented");

  using Self =
      NodeEulerBTreeRmq<T, Compare, Index, LeafSize, Fanout, LowLevelSelector>;

  static constexpr std::size_t npos = RmqBase<Self, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kLeafSize = LeafSize;
  static constexpr std::size_t kFanout = Fanout;
  static constexpr std::size_t kMiddleFanout = 192;

  /**
   * @brief Construct an empty experimental RMQ index.
   */
  NodeEulerBTreeRmq() = default;
  NodeEulerBTreeRmq(const NodeEulerBTreeRmq&) = default;
  NodeEulerBTreeRmq(NodeEulerBTreeRmq&&) noexcept = default;
  NodeEulerBTreeRmq& operator=(const NodeEulerBTreeRmq&) = default;
  NodeEulerBTreeRmq& operator=(NodeEulerBTreeRmq&&) noexcept = default;

  /**
   * @brief Build an experimental B-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller original position as the answer.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit NodeEulerBTreeRmq(std::span<const T> values,
                             Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
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
   * @details Empty or invalid ranges return `npos`. Ties are reduced by
   * comparing `(value, position)`, so traversal order cannot change first-min
   * semantics.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size() || level_sizes_.empty()) {
      return npos;
    }

    const std::size_t root_level = level_count() - 1;
    if (left == 0 && right == values_.size()) {
      return subtree_min_position(root_level, 0);
    }

    const std::size_t left_leaf = leaf_for_value(left);
    const std::size_t right_leaf = leaf_for_value(right - 1);
    if (left_leaf == right_leaf) {
      return leaf_range_min(left_leaf, left, right);
    }

    const auto [level, node_index] = covering_node(left_leaf, right_leaf);
    return query_node(level, node_index, left, right);
  }

 private:
  static constexpr std::size_t kSelectorEntries = 256;
  static constexpr std::size_t kSelectorBits = 2 * kSelectorEntries;
  static constexpr std::size_t kSelectorWords = kSelectorBits / 64;
  static constexpr std::size_t kEmbeddedOffsetEntries = 252;
  static constexpr std::size_t kEmbeddedPositionEntries = 192;
  static constexpr std::size_t kEmbeddedOffsetBit = 2 * kEmbeddedOffsetEntries;
  static constexpr std::size_t kEmbeddedPositionBit =
      2 * kEmbeddedPositionEntries;
  static constexpr std::size_t kMiddleZeroPrefixWord = 7;
  static constexpr std::size_t kHighSparseLevels =
      static_cast<std::size_t>(std::bit_width(Fanout));
  static constexpr std::size_t kHighSparseSlotsPerNode =
      kHighSparseLevels * Fanout;
  static constexpr std::size_t kLeafLinearScanThreshold = 64;
  static constexpr std::size_t kLeafAvx2ScanThreshold = 16;
  static constexpr std::uint64_t kEmbeddedOffsetMask =
      std::uint64_t{0xFF} << (kEmbeddedOffsetBit & 63);
  static constexpr bool kInvalidIndexEqualsNpos =
      static_cast<std::size_t>(invalid_index) == npos;

  static_assert(kEmbeddedOffsetBit + 8 == kSelectorBits);
  static_assert(kEmbeddedPositionBit + 128 == kSelectorBits);
  static_assert(LeafSize <= kEmbeddedOffsetEntries);
  static_assert(kMiddleFanout <= kEmbeddedPositionEntries);

  struct HighChildMetadata {
    std::size_t value_begin = 0;
    std::size_t value_end = 0;
    Index min_position = invalid_index;
  };

  class alignas(64) Bp512Selector {
   public:
    Bp512Selector() = default;

    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kSelectorEntries) {
        throw std::length_error("NodeEulerBTreeRmq local selector too large");
      }

      bp_bits_.fill(0);
      if (entry_count == 0) {
        return;
      }

      std::array<std::uint16_t, kSelectorEntries> stack{};
      std::size_t stack_size = 0;
      std::size_t write_position = 2 * entry_count;

      for (std::size_t i = entry_count; i-- > 0;) {
        while (stack_size != 0 && !entry_less(stack[stack_size - 1], i)) {
          --stack_size;
          prepend_bp_bit(write_position, true);
        }
        stack[stack_size++] = static_cast<std::uint16_t>(i);
        prepend_bp_bit(write_position, false);
      }

      while (write_position != 0) {
        prepend_bp_bit(write_position, true);
      }
    }

    std::size_t arg_min(std::size_t slot_left,
                        std::size_t slot_right,
                        std::size_t entry_count) const {
      if (slot_left >= slot_right || slot_right > entry_count ||
          entry_count > kSelectorEntries) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t bit_count = 2 * entry_count;
      const std::size_t first_close = close_position(slot_left);
      const std::size_t last_close = close_position(slot_right - 1);
      if (first_close > last_close || last_close >= bit_count) {
        return npos;
      }

      const std::size_t shifted_min =
          depth_arg_min(first_close + 1, last_close + 2, bit_count);
      if (shifted_min == npos || shifted_min == 0) {
        return npos;
      }

      const std::size_t zero_rank = rank0_at(shifted_min, bit_count);
      if (zero_rank == 0) {
        return npos;
      }
      const std::size_t entry = zero_rank - 1;
      return entry < entry_count ? entry : npos;
    }

    std::size_t arg_min_with_zero_prefix(std::size_t slot_left,
                                         std::size_t slot_right,
                                         std::size_t entry_count) const {
      if (slot_left >= slot_right || slot_right > entry_count ||
          entry_count > kEmbeddedPositionEntries) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t bit_count = 2 * entry_count;
      const std::size_t first_close =
          close_position_with_zero_prefix(slot_left, bit_count);
      const std::size_t last_close =
          close_position_with_zero_prefix(slot_right - 1, bit_count);
      if (first_close > last_close || last_close >= bit_count) {
        return npos;
      }

      const std::size_t shifted_min = depth_arg_min_with_zero_prefix(
          first_close + 1, last_close + 2, bit_count);
      if (shifted_min == npos || shifted_min == 0) {
        return npos;
      }

      const std::size_t zero_rank =
          rank0_at_with_zero_prefix(shifted_min, bit_count);
      if (zero_rank == 0) {
        return npos;
      }
      const std::size_t entry = zero_rank - 1;
      return entry < entry_count ? entry : npos;
    }

    std::size_t close_position(std::size_t slot) const {
      return select0_512(bp_bits_.data(), slot);
    }

    std::size_t rank0_at(std::size_t position, std::size_t bit_count) const {
      position = std::min(position, bit_count);
      return position - rank_512(bp_bits_.data(), position);
    }

    std::int16_t depth_at_position(std::size_t position) const {
      return static_cast<std::int16_t>(prefix_excess(position));
    }

    void set_embedded_min_offset(std::size_t offset) {
      bp_bits_[kEmbeddedOffsetBit >> 6] =
          (bp_bits_[kEmbeddedOffsetBit >> 6] & ~kEmbeddedOffsetMask) |
          ((static_cast<std::uint64_t>(offset) & 0xFFu)
           << (kEmbeddedOffsetBit & 63));
    }

    std::uint8_t embedded_min_offset() const {
      return static_cast<std::uint8_t>(
          (bp_bits_[kEmbeddedOffsetBit >> 6] & kEmbeddedOffsetMask) >>
          (kEmbeddedOffsetBit & 63));
    }

    void set_embedded_min_position(std::size_t position) {
      bp_bits_[kEmbeddedPositionBit >> 6] =
          static_cast<std::uint64_t>(position);
    }

    std::size_t embedded_min_position() const {
      return static_cast<std::size_t>(bp_bits_[kEmbeddedPositionBit >> 6]);
    }

    void build_zero_prefix_metadata(std::size_t bit_count) {
      const std::size_t word_count = (bit_count + 63) / 64;
      std::uint64_t packed = 0;
      std::size_t zeros = 0;
      for (std::size_t word = 0; word <= word_count; ++word) {
        packed |= (static_cast<std::uint64_t>(zeros) & 0xFFu) << (8 * word);
        if (word == word_count) {
          break;
        }
        const std::size_t word_begin = word * 64;
        const std::size_t word_bits =
            std::min<std::size_t>(64, bit_count - word_begin);
        const std::uint64_t bits = bp_bits_[word] & first_bits_mask(word_bits);
        zeros += word_bits - std::popcount(bits);
      }
      bp_bits_[kMiddleZeroPrefixWord] = packed;
    }

   private:
    std::size_t prepend_bp_bit(std::size_t& write_position, bool bit) {
      --write_position;
      if (bit) {
        bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                         << (write_position & 63);
      }
      return write_position;
    }

    int prefix_excess(std::size_t position) const {
      position = std::min(position, kSelectorBits);
      const std::size_t ones = rank_512(bp_bits_.data(), position);
      return static_cast<int>(2 * ones) - static_cast<int>(position);
    }

    std::size_t zero_prefix_at_word(std::size_t word) const {
      return static_cast<std::uint8_t>(bp_bits_[kMiddleZeroPrefixWord] >>
                                       (8 * word));
    }

    std::size_t rank0_at_with_zero_prefix(std::size_t position,
                                          std::size_t bit_count) const {
      position = std::min(position, bit_count);
      const std::size_t full_words = position >> 6;
      std::size_t zeros = zero_prefix_at_word(full_words);
      const std::size_t tail_bits = position & 63;
      if (tail_bits != 0) {
        const std::uint64_t tail =
            bp_bits_[full_words] & first_bits_mask(tail_bits);
        zeros += tail_bits - std::popcount(tail);
      }
      return zeros;
    }

    int prefix_excess_with_zero_prefix(std::size_t position,
                                       std::size_t bit_count) const {
      position = std::min(position, bit_count);
      return static_cast<int>(position) -
             2 * static_cast<int>(
                     rank0_at_with_zero_prefix(position, bit_count));
    }

    std::size_t close_position_with_zero_prefix(std::size_t slot,
                                                std::size_t bit_count) const {
      const std::size_t word_count = (bit_count + 63) / 64;
      for (std::size_t word = 0; word < word_count; ++word) {
        const std::size_t next_zero_prefix = zero_prefix_at_word(word + 1);
        if (next_zero_prefix <= slot) {
          continue;
        }
        const std::size_t local_rank = slot - zero_prefix_at_word(word);
        const std::size_t word_begin = word * 64;
        const std::size_t word_bits =
            std::min<std::size_t>(64, bit_count - word_begin);
        const std::uint64_t zeros =
            (~bp_bits_[word]) & first_bits_mask(word_bits);
        return word_begin + select_64(zeros, local_rank);
      }
      return npos;
    }

    std::size_t depth_arg_min(std::size_t left,
                              std::size_t right,
                              std::size_t bit_count) const {
      const std::size_t depth_count = bit_count + 1;
      if (left >= right || right > depth_count) {
        return npos;
      }

      std::size_t position = left;
      int best_depth = prefix_excess(position);
      std::size_t best_position = position;

      while (position < right) {
        const std::size_t chunk_begin = (position / 128) * 128;
        const std::size_t local_left = position - chunk_begin;
        const std::size_t local_right =
            std::min<std::size_t>(right - 1, chunk_begin + 128) - chunk_begin;

        int candidate_depth;
        std::size_t candidate_position;
        if (chunk_begin >= bit_count) {
          candidate_depth = prefix_excess(chunk_begin);
          candidate_position = chunk_begin;
        } else {
          const std::size_t word = chunk_begin >> 6;
          const ExcessResult result =
              excess_min_128(bp_bits_.data() + word, local_left, local_right);
          candidate_depth = prefix_excess(chunk_begin) + result.min_excess;
          candidate_position = chunk_begin + result.offset;
        }

        if (candidate_depth < best_depth) {
          best_depth = candidate_depth;
          best_position = candidate_position;
        }

        position = chunk_begin + local_right + 1;
      }

      return best_position;
    }

    std::size_t depth_arg_min_with_zero_prefix(std::size_t left,
                                               std::size_t right,
                                               std::size_t bit_count) const {
      const std::size_t depth_count = bit_count + 1;
      if (left >= right || right > depth_count) {
        return npos;
      }

      std::size_t position = left;
      int best_depth = prefix_excess_with_zero_prefix(position, bit_count);
      std::size_t best_position = position;

      while (position < right) {
        const std::size_t chunk_begin = (position / 128) * 128;
        const std::size_t local_left = position - chunk_begin;
        const std::size_t local_right =
            std::min<std::size_t>(right - 1, chunk_begin + 128) - chunk_begin;

        int candidate_depth;
        std::size_t candidate_position;
        if (chunk_begin >= bit_count) {
          candidate_depth =
              prefix_excess_with_zero_prefix(chunk_begin, bit_count);
          candidate_position = chunk_begin;
        } else {
          const std::size_t word = chunk_begin >> 6;
          const ExcessResult result =
              excess_min_128(bp_bits_.data() + word, local_left, local_right);
          candidate_depth =
              prefix_excess_with_zero_prefix(chunk_begin, bit_count) +
              result.min_excess;
          candidate_position = chunk_begin + result.offset;
        }

        if (candidate_depth < best_depth) {
          best_depth = candidate_depth;
          best_position = candidate_position;
        }

        position = chunk_begin + local_right + 1;
      }

      return best_position;
    }

    std::array<std::uint64_t, kSelectorWords> bp_bits_{};
  };

  static_assert(sizeof(Bp512Selector) == 64);

  class TopDepthSparseSelector {
   public:
    TopDepthSparseSelector() = default;

    void build(const Bp512Selector& selector, std::size_t entry_count) {
      const std::size_t depth_count = 2 * entry_count + 1;
      depth_count_ = depth_count;
      depths_.assign(depth_count, 0);
      if (depth_count == 0) {
        log_count_ = 0;
        sparse_positions_.clear();
        return;
      }

      for (std::size_t position = 0; position < depth_count; ++position) {
        depths_[position] = selector.depth_at_position(position);
      }

      log_count_ = std::bit_width(depth_count);
      sparse_positions_.assign(log_count_ * depth_count, 0);
      for (std::size_t position = 0; position < depth_count; ++position) {
        sparse_positions_[position] = static_cast<std::uint16_t>(position);
      }

      for (std::size_t level = 1; level < log_count_; ++level) {
        const std::size_t half_span = std::size_t{1} << (level - 1);
        const std::size_t span = half_span << 1;
        if (span > depth_count) {
          break;
        }
        const std::size_t previous_offset = (level - 1) * depth_count;
        const std::size_t current_offset = level * depth_count;
        for (std::size_t position = 0; position + span <= depth_count;
             ++position) {
          sparse_positions_[current_offset + position] = better_depth_position(
              sparse_positions_[previous_offset + position],
              sparse_positions_[previous_offset + position + half_span]);
        }
      }
    }

    std::size_t arg_min(const Bp512Selector& selector,
                        std::size_t slot_left,
                        std::size_t slot_right,
                        std::size_t entry_count) const {
      if (slot_left >= slot_right || slot_right > entry_count ||
          entry_count > kSelectorEntries) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t bit_count = 2 * entry_count;
      const std::size_t first_close = selector.close_position(slot_left);
      const std::size_t last_close = selector.close_position(slot_right - 1);
      if (first_close > last_close || last_close >= bit_count) {
        return npos;
      }

      const std::size_t shifted_min =
          depth_arg_min(first_close + 1, last_close + 2);
      if (shifted_min == npos || shifted_min == 0) {
        return npos;
      }

      const std::size_t zero_rank = selector.rank0_at(shifted_min, bit_count);
      if (zero_rank == 0) {
        return npos;
      }
      const std::size_t entry = zero_rank - 1;
      return entry < entry_count ? entry : npos;
    }

   private:
    std::uint16_t better_depth_position(std::uint16_t left,
                                        std::uint16_t right) const {
      const std::int16_t left_depth = depths_[left];
      const std::int16_t right_depth = depths_[right];
      if (right_depth < left_depth) {
        return right;
      }
      if (left_depth < right_depth) {
        return left;
      }
      return std::min(left, right);
    }

    std::size_t depth_arg_min(std::size_t left, std::size_t right) const {
      if (left >= right || right > depth_count_) {
        return npos;
      }
      const std::size_t length = right - left;
      const std::size_t level = std::bit_width(length) - 1;
      const std::size_t span = std::size_t{1} << level;
      const std::size_t offset = level * depth_count_;
      return better_depth_position(sparse_positions_[offset + left],
                                   sparse_positions_[offset + right - span]);
    }

    std::vector<std::int16_t> depths_;
    std::vector<std::uint16_t> sparse_positions_;
    std::size_t depth_count_ = 0;
    std::size_t log_count_ = 0;
  };

  bool missing_position(std::size_t position) const {
    if constexpr (kInvalidIndexEqualsNpos) {
      return position == npos;
    } else {
      return position == npos ||
             position == static_cast<std::size_t>(invalid_index);
    }
  }

  bool strictly_better_position(std::size_t left, std::size_t right) const {
    if (missing_position(left)) {
      return false;
    }
    if (missing_position(right)) {
      return true;
    }
    return compare_(values_[left], values_[right]);
  }

  std::size_t better_position(std::size_t left, std::size_t right) const {
    if (missing_position(left)) {
      return right;
    }
    if (missing_position(right)) {
      return left;
    }
    if (compare_(values_[right], values_[left])) {
      return right;
    }
    if (compare_(values_[left], values_[right])) {
      return left;
    }
    return std::min(left, right);
  }

  void build() {
    leaf_selectors_.clear();
    medium_selectors_.clear();
    high_selectors_.clear();
    high_depth_selectors_.clear();
    high_min_positions_.clear();
    high_child_metadata_.clear();
    high_sparse_min_slots_.clear();
    medium_level_offsets_.clear();
    high_level_offsets_.clear();
    level_sizes_.clear();
    level_value_spans_.clear();
    level_fanouts_.clear();
    high_level_begin_ = std::numeric_limits<std::size_t>::max();
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("NodeEulerBTreeRmq index type is too small");
    }

    initialize_layout((values_.size() + LeafSize - 1) / LeafSize);
    for (std::size_t leaf = 0; leaf < level_sizes_[0]; ++leaf) {
      build_leaf(leaf);
    }

    for (std::size_t level = 1; level < level_count(); ++level) {
      for (std::size_t node = 0; node < level_sizes_[level]; ++node) {
        build_internal_node(level, node);
      }
    }
  }

  void initialize_layout(std::size_t leaf_count) {
    level_sizes_.push_back(leaf_count);
    level_value_spans_.push_back(LeafSize);
    level_fanouts_.push_back(0);

    std::size_t current_count = leaf_count;
    std::size_t current_span = LeafSize;
    while (current_count > Fanout * Fanout) {
      level_fanouts_.push_back(kMiddleFanout);
      current_count = ceil_div(current_count, kMiddleFanout);
      current_span = saturating_product(current_span, kMiddleFanout);
      level_sizes_.push_back(current_count);
      level_value_spans_.push_back(current_span);
    }
    while (current_count > 1) {
      level_fanouts_.push_back(Fanout);
      current_count = ceil_div(current_count, Fanout);
      current_span = saturating_product(current_span, Fanout);
      level_sizes_.push_back(current_count);
      level_value_spans_.push_back(current_span);
    }

    leaf_selectors_.resize(level_sizes_[0]);

    medium_level_offsets_.assign(level_count(), 0);
    high_level_offsets_.assign(level_count(), 0);
    if (level_count() <= 1) {
      high_level_begin_ = std::numeric_limits<std::size_t>::max();
      return;
    }

    const std::size_t root_level = level_count() - 1;
    high_level_begin_ = root_level > 1 ? root_level - 1 : 1;

    std::size_t medium_node_count = 0;
    for (std::size_t level = 1; level < high_level_begin_; ++level) {
      medium_level_offsets_[level] = medium_node_count;
      medium_node_count += level_sizes_[level];
    }
    medium_selectors_.resize(medium_node_count);

    std::size_t high_node_count = 0;
    for (std::size_t level = high_level_begin_; level < level_count();
         ++level) {
      high_level_offsets_[level] = high_node_count;
      high_node_count += level_sizes_[level];
    }
    high_selectors_.resize(high_node_count);
    high_depth_selectors_.resize(high_node_count);
    high_min_positions_.resize(high_node_count, invalid_index);
    high_child_metadata_.resize(high_node_count * Fanout);
    high_sparse_min_slots_.resize(high_node_count * kHighSparseSlotsPerNode);
  }

  void build_leaf(std::size_t leaf) {
    Bp512Selector& selector = leaf_selectors_[leaf];
    const std::size_t begin = node_value_begin(0, leaf);
    const std::size_t count = entry_count(0, leaf);
    selector.build(count, [&](std::size_t left, std::size_t right) {
      return compare_(values_[begin + left], values_[begin + right]);
    });

    const std::size_t slot = selector.arg_min(0, count, count);
    selector.set_embedded_min_offset(slot);
  }

  void build_internal_node(std::size_t level, std::size_t node) {
    Bp512Selector& selector = mutable_selector_at(level, node);
    const std::size_t count = entry_count(level, node);
    const std::size_t first_child = node * fanout_at_level(level);
    const bool high_level = is_high_level(level);
    const std::size_t high_flat = high_level ? high_flat_index(level, node) : 0;
    if (high_level) {
      for (std::size_t slot = 0; slot < count; ++slot) {
        const std::size_t child = first_child + slot;
        HighChildMetadata& metadata =
            mutable_high_child_metadata_at(high_flat, slot);
        metadata.value_begin = node_value_begin(level - 1, child);
        metadata.value_end = node_value_end(level - 1, child);
        metadata.min_position =
            static_cast<Index>(subtree_min_position(level - 1, child));
      }
      build_high_sparse_min_slots(high_flat, count);
    }

    selector.build(count, [&](std::size_t left, std::size_t right) {
      if (high_level) {
        return strictly_better_position(
            high_child_metadata_at(high_flat, left).min_position,
            high_child_metadata_at(high_flat, right).min_position);
      }
      return strictly_better_position(
          subtree_min_position(level - 1, first_child + left),
          subtree_min_position(level - 1, first_child + right));
    });

    const std::size_t slot = selector.arg_min(0, count, count);
    const std::size_t min_position =
        high_level ? high_child_metadata_at(high_flat, slot).min_position
                   : subtree_min_position(level - 1, first_child + slot);
    if (high_level) {
      high_min_positions_[high_flat] = static_cast<Index>(min_position);
      mutable_high_depth_selector_at(level, node).build(selector, count);
    } else {
      selector.set_embedded_min_position(min_position);
      selector.build_zero_prefix_metadata(2 * count);
    }
  }

  static std::size_t saturating_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
      return std::numeric_limits<std::size_t>::max();
    }
    return left * right;
  }

  static std::size_t ceil_div(std::size_t value, std::size_t divisor) {
    return (value + divisor - 1) / divisor;
  }

  std::size_t leaf_range_min(std::size_t leaf,
                             std::size_t left,
                             std::size_t right) const {
    if (left >= right) {
      return npos;
    }

    const std::size_t begin = node_value_begin(0, leaf);
    const std::size_t end = node_value_end(0, leaf);
    if (left <= begin && end <= right) {
      return subtree_min_position(0, leaf);
    }
    if (right - left <= kLeafLinearScanThreshold) {
      return linear_range_min(left, right);
    }

    const std::size_t slot_left = left - begin;
    const std::size_t slot_right = right - begin;
    const std::size_t slot = leaf_selectors_[leaf].arg_min(
        slot_left, slot_right, entry_count(0, leaf));
    if (slot == npos) {
      return npos;
    }
    return begin + slot;
  }

  std::size_t linear_range_min(std::size_t left, std::size_t right) const {
    if (left >= right) {
      return npos;
    }
#ifdef PIXIE_AVX2_SUPPORT
    if constexpr (std::is_same_v<T, std::int64_t> &&
                  std::is_same_v<Compare, std::less<T>>) {
      if (right - left >= kLeafAvx2ScanThreshold) {
        return linear_range_min_i64_avx2(left, right);
      }
    }
#endif
    std::size_t best = left;
    for (std::size_t position = left + 1; position < right; ++position) {
      if (compare_(values_[position], values_[best])) {
        best = position;
      }
    }
    return best;
  }

#ifdef PIXIE_AVX2_SUPPORT
  std::size_t linear_range_min_i64_avx2(std::size_t left,
                                        std::size_t right) const {
    const std::int64_t* data = values_.data();
    std::size_t position = left;

    __m256i best_values =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + position));
    __m256i best_positions = _mm256_set_epi64x(
        static_cast<long long>(position + 3),
        static_cast<long long>(position + 2),
        static_cast<long long>(position + 1), static_cast<long long>(position));
    position += 4;

    for (; position + 4 <= right; position += 4) {
      const __m256i values =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + position));
      const __m256i positions =
          _mm256_set_epi64x(static_cast<long long>(position + 3),
                            static_cast<long long>(position + 2),
                            static_cast<long long>(position + 1),
                            static_cast<long long>(position));
      const __m256i take_new = _mm256_cmpgt_epi64(best_values, values);
      best_values = _mm256_blendv_epi8(best_values, values, take_new);
      best_positions = _mm256_blendv_epi8(best_positions, positions, take_new);
    }

    alignas(32) std::int64_t value_lanes[4];
    alignas(32) std::uint64_t position_lanes[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(value_lanes), best_values);
    _mm256_store_si256(reinterpret_cast<__m256i*>(position_lanes),
                       best_positions);

    std::int64_t best_value = value_lanes[0];
    std::size_t best_position = static_cast<std::size_t>(position_lanes[0]);
    for (std::size_t lane = 1; lane < 4; ++lane) {
      const std::size_t lane_position =
          static_cast<std::size_t>(position_lanes[lane]);
      if (value_lanes[lane] < best_value ||
          (value_lanes[lane] == best_value && lane_position < best_position)) {
        best_value = value_lanes[lane];
        best_position = lane_position;
      }
    }

    for (; position < right; ++position) {
      if (data[position] < best_value) {
        best_value = data[position];
        best_position = position;
      }
    }
    return best_position;
  }
#endif

  std::pair<std::size_t, std::size_t> covering_node(
      std::size_t left_leaf,
      std::size_t right_leaf) const {
    std::size_t level = 0;
    std::size_t left_node = left_leaf;
    std::size_t right_node = right_leaf;
    while (left_node != right_node) {
      ++level;
      const std::size_t fanout = fanout_at_level(level);
      left_node /= fanout;
      right_node /= fanout;
    }
    return {level, left_node};
  }

  std::size_t leaf_for_value(std::size_t position) const {
    return position / LeafSize;
  }

  std::size_t child_for_value(std::size_t child_level,
                              std::size_t position) const {
    return position / level_value_spans_[child_level];
  }

  std::size_t query_child_slots(std::size_t level,
                                std::size_t node,
                                std::size_t slot_left,
                                std::size_t slot_right,
                                std::size_t left,
                                std::size_t right) const {
    if (slot_left >= slot_right) {
      return npos;
    }

    const std::size_t count = entry_count(level, node);
    const std::size_t slot =
        slot_left + 1 == slot_right
            ? slot_left
            : selector_arg_min(level, node, slot_left, slot_right, count);
    if (slot == npos) {
      return npos;
    }

    const std::size_t child_level = level - 1;
    const std::size_t first_child = node * fanout_at_level(level);
    const std::size_t child = first_child + slot;
    const HighChildMetadata* high_children =
        is_high_level(level) ? high_child_metadata_begin(level, node) : nullptr;
    const std::size_t child_begin = high_children != nullptr
                                        ? high_children[slot].value_begin
                                        : node_value_begin(child_level, child);
    const std::size_t child_end = high_children != nullptr
                                      ? high_children[slot].value_end
                                      : node_value_end(child_level, child);
    const std::size_t child_min =
        high_children != nullptr ? high_children[slot].min_position
                                 : subtree_min_position(child_level, child);
    if ((left <= child_begin && child_end <= right) ||
        contains_position(left, right, child_min)) {
      return child_min;
    }

    const std::size_t last_slot = slot_right - 1;
    const std::size_t left_child_begin =
        high_children != nullptr
            ? high_children[slot_left].value_begin
            : node_value_begin(child_level, first_child + slot_left);
    const std::size_t left_child_end =
        high_children != nullptr
            ? high_children[slot_left].value_end
            : node_value_end(child_level, first_child + slot_left);
    std::size_t answer = query_node(child_level, first_child + slot_left,
                                    std::max(left, left_child_begin),
                                    std::min(right, left_child_end));

    if (slot_left != last_slot) {
      const std::size_t right_child_begin =
          high_children != nullptr
              ? high_children[last_slot].value_begin
              : node_value_begin(child_level, first_child + last_slot);
      const std::size_t right_child_end =
          high_children != nullptr
              ? high_children[last_slot].value_end
              : node_value_end(child_level, first_child + last_slot);
      answer = better_position(answer,
                               query_node(child_level, first_child + last_slot,
                                          std::max(left, right_child_begin),
                                          std::min(right, right_child_end)));
    }

    if (slot_left + 1 < last_slot) {
      answer = better_position(
          answer,
          full_child_slot_range_min(level, node, slot_left + 1, last_slot));
    }

    return answer;
  }

  std::size_t full_child_slot_range_min(std::size_t level,
                                        std::size_t node,
                                        std::size_t slot_left,
                                        std::size_t slot_right) const {
    if (slot_left >= slot_right) {
      return npos;
    }

    const std::size_t slot =
        slot_left + 1 == slot_right
            ? slot_left
            : selector_arg_min(level, node, slot_left, slot_right,
                               entry_count(level, node));
    if (slot == npos) {
      return npos;
    }

    if (is_high_level(level)) {
      return high_child_metadata_begin(level, node)[slot].min_position;
    }
    return subtree_min_position(level - 1,
                                node * fanout_at_level(level) + slot);
  }

  std::size_t query_node(std::size_t level,
                         std::size_t node,
                         std::size_t left,
                         std::size_t right) const {
    if (left >= right) {
      return npos;
    }
    const std::size_t begin = node_value_begin(level, node);
    const std::size_t end = node_value_end(level, node);
    if (left <= begin && end <= right) {
      return subtree_min_position(level, node);
    }
    if (level == 0) {
      return leaf_range_min(node, left, right);
    }

    const std::size_t child_level = level - 1;
    const std::size_t left_child = child_for_value(child_level, left);
    const std::size_t right_child = child_for_value(child_level, right - 1);
    const std::size_t first_child = node * fanout_at_level(level);
    const std::size_t left_slot = left_child - first_child;
    const std::size_t right_slot = right_child - first_child + 1;
    return query_child_slots(level, node, left_slot, right_slot, left, right);
  }

  bool contains_position(std::size_t left,
                         std::size_t right,
                         std::size_t position) const {
    return !missing_position(position) && left <= position && position < right;
  }

  std::size_t level_count() const { return level_sizes_.size(); }

  std::size_t entry_count(std::size_t level, std::size_t node) const {
    if (level == 0) {
      const std::size_t begin = node_value_begin(0, node);
      return std::min<std::size_t>(LeafSize, values_.size() - begin);
    }
    const std::size_t first_child = node * fanout_at_level(level);
    return std::min<std::size_t>(fanout_at_level(level),
                                 level_sizes_[level - 1] - first_child);
  }

  std::size_t node_value_begin(std::size_t level, std::size_t node) const {
    return node * level_value_spans_[level];
  }

  std::size_t node_value_end(std::size_t level, std::size_t node) const {
    return std::min(values_.size(),
                    node_value_begin(level, node) + level_value_spans_[level]);
  }

  std::size_t subtree_min_position(std::size_t level, std::size_t node) const {
    if (level == 0) {
      return node_value_begin(0, node) +
             leaf_selectors_[node].embedded_min_offset();
    }
    if (is_high_level(level)) {
      return high_min_positions_[high_flat_index(level, node)];
    }
    return selector_at(level, node).embedded_min_position();
  }

  const Bp512Selector& selector_at(std::size_t level, std::size_t node) const {
    if (level == 0) {
      return leaf_selectors_[node];
    }
    if (is_high_level(level)) {
      return high_selectors_[high_flat_index(level, node)];
    }
    return medium_selectors_[medium_flat_index(level, node)];
  }

  Bp512Selector& mutable_selector_at(std::size_t level, std::size_t node) {
    if (is_high_level(level)) {
      return high_selectors_[high_flat_index(level, node)];
    }
    return medium_selectors_[medium_flat_index(level, node)];
  }

  std::size_t selector_arg_min(std::size_t level,
                               std::size_t node,
                               std::size_t slot_left,
                               std::size_t slot_right,
                               std::size_t count) const {
    if (is_high_level(level)) {
      return high_sparse_arg_min(high_flat_index(level, node), slot_left,
                                 slot_right, count);
    }
    const Bp512Selector& selector = selector_at(level, node);
    return selector.arg_min_with_zero_prefix(slot_left, slot_right, count);
  }

  bool is_high_level(std::size_t level) const {
    return level > 0 && level >= high_level_begin_ && level < level_count();
  }

  std::size_t medium_flat_index(std::size_t level, std::size_t node) const {
    return medium_level_offsets_[level] + node;
  }

  std::size_t fanout_at_level(std::size_t level) const {
    return level_fanouts_[level];
  }

  std::size_t high_flat_index(std::size_t level, std::size_t node) const {
    return high_level_offsets_[level] + node;
  }

  const HighChildMetadata* high_child_metadata_begin(std::size_t level,
                                                     std::size_t node) const {
    return high_child_metadata_.data() + high_flat_index(level, node) * Fanout;
  }

  std::uint8_t* mutable_high_sparse_min_slots_begin(std::size_t high_flat) {
    return high_sparse_min_slots_.data() + high_flat * kHighSparseSlotsPerNode;
  }

  const std::uint8_t* high_sparse_min_slots_begin(std::size_t high_flat) const {
    return high_sparse_min_slots_.data() + high_flat * kHighSparseSlotsPerNode;
  }

  const HighChildMetadata& high_child_metadata_at(std::size_t high_flat,
                                                  std::size_t slot) const {
    return high_child_metadata_[high_flat * Fanout + slot];
  }

  HighChildMetadata& mutable_high_child_metadata_at(std::size_t high_flat,
                                                    std::size_t slot) {
    return high_child_metadata_[high_flat * Fanout + slot];
  }

  std::size_t better_high_child_slot(std::size_t high_flat,
                                     std::size_t left_slot,
                                     std::size_t right_slot) const {
    const std::size_t left_position =
        high_child_metadata_at(high_flat, left_slot).min_position;
    const std::size_t right_position =
        high_child_metadata_at(high_flat, right_slot).min_position;
    return better_position(left_position, right_position) == right_position
               ? right_slot
               : left_slot;
  }

  void build_high_sparse_min_slots(std::size_t high_flat, std::size_t count) {
    std::uint8_t* table = mutable_high_sparse_min_slots_begin(high_flat);
    for (std::size_t slot = 0; slot < count; ++slot) {
      table[slot] = static_cast<std::uint8_t>(slot);
    }

    for (std::size_t table_level = 1; table_level < kHighSparseLevels;
         ++table_level) {
      const std::size_t span = std::size_t{1} << table_level;
      if (span > count) {
        break;
      }
      const std::size_t half_span = span >> 1;
      const std::uint8_t* previous = table + (table_level - 1) * Fanout;
      std::uint8_t* current = table + table_level * Fanout;
      for (std::size_t slot = 0; slot + span <= count; ++slot) {
        current[slot] = static_cast<std::uint8_t>(better_high_child_slot(
            high_flat, previous[slot], previous[slot + half_span]));
      }
    }
  }

  std::size_t high_sparse_arg_min(std::size_t high_flat,
                                  std::size_t slot_left,
                                  std::size_t slot_right,
                                  std::size_t count) const {
    if (slot_left >= slot_right || slot_right > count) {
      return npos;
    }
    const std::size_t length = slot_right - slot_left;
    if (length == 1) {
      return slot_left;
    }

    const std::size_t table_level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << table_level;
    const std::uint8_t* table =
        high_sparse_min_slots_begin(high_flat) + table_level * Fanout;
    return better_high_child_slot(high_flat, table[slot_left],
                                  table[slot_right - span]);
  }

  const TopDepthSparseSelector& high_depth_selector_at(std::size_t level,
                                                       std::size_t node) const {
    return high_depth_selectors_[high_flat_index(level, node)];
  }

  TopDepthSparseSelector& mutable_high_depth_selector_at(std::size_t level,
                                                         std::size_t node) {
    return high_depth_selectors_[high_flat_index(level, node)];
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<Bp512Selector> leaf_selectors_;
  std::vector<Bp512Selector> medium_selectors_;
  std::vector<Bp512Selector> high_selectors_;
  std::vector<TopDepthSparseSelector> high_depth_selectors_;
  std::vector<Index> high_min_positions_;
  std::vector<HighChildMetadata> high_child_metadata_;
  std::vector<std::uint8_t> high_sparse_min_slots_;
  std::vector<std::size_t> medium_level_offsets_;
  std::vector<std::size_t> high_level_offsets_;
  std::vector<std::size_t> level_sizes_;
  std::vector<std::size_t> level_value_spans_;
  std::vector<std::size_t> level_fanouts_;
  std::size_t high_level_begin_ = std::numeric_limits<std::size_t>::max();
};

}  // namespace pixie::rmq::experimental
