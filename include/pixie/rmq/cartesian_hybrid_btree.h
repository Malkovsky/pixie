#pragma once

#include <pixie/bits.h>
#include <pixie/bitvector.h>
#include <pixie/cache_line.h>
#include <pixie/memory_usage.h>
#include <pixie/rmq/rmq_base.h>
#include <pixie/rmq/utils/succinct_monotone_stack.h>

#include <algorithm>
#include <array>
#include <bit>
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
 * @brief HybridBTree-style RMQ backend for ±1 depth sequences.
 *
 * @details The input words encode adjacent depth deltas: bit 1 means +1 and
 * bit 0 means -1. The indexed sequence has `depth_count` prefix positions, so
 * the bit sequence has `depth_count - 1` deltas. Public ranges are half-open
 * over depth positions, and ties return the smaller depth position.
 *
 * This backend mirrors the high-level query shape of `HybridBTree`: depth
 * positions are split into configurable-size leaves, leaves are grouped into a
 * B-tree, and every internal node stores a local Cartesian/BP selector over the
 * minima of its immediate children. Middle levels are fixed at 192 child slots:
 * their selector uses 384 BP bits and reserves the remaining 128 bits for the
 * embedded subtree minimum position/depth. A configurable number of top levels
 * additionally keep sparse tables over child-minimum slots. Unlike value
 * `HybridBTree`, leaves do not store another local Cartesian selector or
 * cached minima. Leaf minima are recomputed from the original BP delta bits
 * with the 128-bit excess primitives from `bits.h`, using rank support to
 * recover the absolute base depth at the leaf start.
 *
 * @tparam Index Unsigned integer type used for stored positions.
 * @tparam LeafSize Number of depth positions per low-level leaf; must be a
 * multiple of 512.
 * @tparam UseHighSparseLayout Whether top levels use sparse tables instead of
 * local BP selectors.
 * @tparam HighSparseLayoutLevels Number of top tree levels using that sparse
 * layout when enabled.
 */
template <class Index = std::size_t,
          std::size_t LeafSize = 512,
          bool UseHighSparseLayout = true,
          std::size_t HighSparseLayoutLevels = 2>
class HybridBTreePlusMinusOne {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "HybridBTreePlusMinusOne index type must be unsigned");
  static_assert(LeafSize != 0 && LeafSize % 512 == 0,
                "HybridBTreePlusMinusOne leaf size must be a positive "
                "multiple of 512");
  static_assert(!UseHighSparseLayout || HighSparseLayoutLevels > 0,
                "HybridBTreePlusMinusOne high sparse layout must cover "
                "at least one level when enabled");

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kLeafSize = LeafSize;
  static constexpr std::size_t kHighLevelFanout = 256;
  static constexpr std::size_t kMiddleFanout = 192;
  static constexpr bool kUseHighSparseLayout = UseHighSparseLayout;
  static constexpr std::size_t kHighSparseLayoutLevels = HighSparseLayoutLevels;

  /**
   * @brief Construct an empty ±1 RMQ index.
   */
  HybridBTreePlusMinusOne() = default;

  /**
   * @brief Build a ±1 RMQ index over external packed delta bits.
   *
   * @param bits Little-endian packed delta bits.
   * @param depth_count Number of indexed depth positions.
   * @throws std::length_error if @p Index cannot represent all positions.
   * @throws std::invalid_argument if @p bits is shorter than required.
   */
  HybridBTreePlusMinusOne(std::span<const std::uint64_t> bits,
                          std::size_t depth_count) {
    build(bits, depth_count);
  }

  /**
   * @brief Build a ±1 RMQ index over external packed bits and rank support.
   *
   * @details The supplied rank index is non-owning and must outlive this RMQ
   * object. This is used by CartesianHybridBTree to reuse its BP rank/select
   * index.
   */
  HybridBTreePlusMinusOne(std::span<const std::uint64_t> bits,
                          std::size_t depth_count,
                          const BitVector& rank_index) {
    build(bits, depth_count, rank_index);
  }

  /**
   * @brief Rebuild this index over external packed delta bits.
   */
  void build(std::span<const std::uint64_t> bits, std::size_t depth_count) {
    input_bits_ = bits;
    depth_count_ = depth_count;
    external_rank_index_ = nullptr;
    build();
  }

  /**
   * @brief Rebuild this index using non-owning rank support for the same bits.
   */
  void build(std::span<const std::uint64_t> bits,
             std::size_t depth_count,
             const BitVector& rank_index) {
    input_bits_ = bits;
    depth_count_ = depth_count;
    external_rank_index_ = &rank_index;
    build();
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
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this ±1 RMQ object and all selector/metadata buffers. The
   * external packed delta bits are not owned and are excluded.
   */
  std::size_t memory_usage_bytes() const {
    std::size_t bytes = sizeof(*this);
    if (external_rank_index_ == nullptr) {
      bytes += pixie::optional_nested_owned_memory_bytes(owned_rank_index_);
    }
    bytes += pixie::vector_capacity_bytes(internal_selectors_);
    bytes += pixie::vector_capacity_bytes(internal_min_positions_);
    bytes += pixie::vector_capacity_bytes(internal_min_depths_);
    bytes += pixie::vector_capacity_bytes(high_child_metadata_);
    bytes += pixie::vector_capacity_bytes(high_sparse_min_slots_);
    bytes += pixie::vector_capacity_bytes(internal_level_offsets_);
    bytes += pixie::vector_capacity_bytes(min_summary_level_offsets_);
    bytes += pixie::vector_capacity_bytes(high_level_offsets_);
    bytes += pixie::vector_capacity_bytes(level_sizes_);
    bytes += pixie::vector_capacity_bytes(level_position_spans_);
    bytes += pixie::vector_capacity_bytes(level_fanouts_);
    return bytes;
  }

  /**
   * @brief Return the first minimum depth position in [@p left, @p right).
   *
   * @details Empty or invalid ranges return `npos`. Equal depths return the
   * smaller position.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    if (left >= right || right > depth_count_ || level_sizes_.empty()) {
      return npos;
    }
    if (left + 1 == right) {
      return left;
    }

    const std::size_t root_level = level_count() - 1;
    if (left == 0 && right == depth_count_) {
      return subtree_min_candidate(root_level, 0).position;
    }

    const std::size_t left_leaf = leaf_for_position(left);
    const std::size_t right_leaf = leaf_for_position(right - 1);
    if (left_leaf == right_leaf) {
      return leaf_range_arg_min_relative(left_leaf, left, right);
    }

    const auto [level, node] = covering_node(left_leaf, right_leaf);
    return query_node(level, node, left, right).position;
  }

  /**
   * @brief Return the one-based rank-th zero delta bit position.
   *
   * @details This is the select operation needed by the Cartesian wrapper to
   * map value endpoints to close-parenthesis positions. The returned position
   * is a delta-bit position in `[0, size() - 1)`, or `npos` when @p rank is out
   * of range.
   */
  std::size_t select0(std::size_t rank) const {
    if (rank == 0 || depth_count_ == 0 || rank_index_or_null() == nullptr) {
      return npos;
    }
    const std::size_t delta_count = depth_count_ - 1;
    const std::size_t position = rank_index().select0(rank);
    return position < delta_count ? position : npos;
  }

 private:
  static constexpr std::size_t kLeafWords = LeafSize / 64;
  static constexpr std::size_t kLeafChunks = LeafSize / 128;
  static constexpr std::size_t kSelectorEntries = 256;
  static constexpr std::size_t kSelectorBits = 2 * kSelectorEntries;
  static constexpr std::size_t kSelectorWords = kSelectorBits / 64;
  static constexpr std::size_t kEmbeddedSummaryWords = 2;
  static constexpr std::size_t kEmbeddedSummaryBits =
      64 * kEmbeddedSummaryWords;
  static constexpr std::size_t kEmbeddedSummaryMaxEntries =
      (kSelectorBits - kEmbeddedSummaryBits) / 2;
  static constexpr std::size_t kEmbeddedSummaryPositionWord =
      kSelectorWords - 2;
  static constexpr std::size_t kEmbeddedSummaryDepthWord = kSelectorWords - 1;
  static constexpr std::size_t kHighSparseTableLevels =
      static_cast<std::size_t>(std::bit_width(kHighLevelFanout));
  static constexpr std::size_t kHighSparseSlotsPerNode =
      kHighSparseTableLevels * kHighLevelFanout;
  static constexpr bool kInvalidIndexEqualsNpos =
      static_cast<std::size_t>(invalid_index) == npos;
  static_assert(kEmbeddedSummaryMaxEntries == 192);
  static_assert(kMiddleFanout == kEmbeddedSummaryMaxEntries);
  static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t));

  struct DepthCandidate {
    std::size_t position = npos;
    std::int64_t depth = std::numeric_limits<std::int64_t>::max();
  };

  struct HighChildMetadata {
    std::size_t position_begin = 0;
    std::size_t position_end = 0;
    Index min_position = invalid_index;
    std::int64_t min_depth = std::numeric_limits<std::int64_t>::max();
  };

  class alignas(64) Bp512Selector {
   public:
    /**
     * @brief Construct an empty packed BP selector.
     */
    Bp512Selector() = default;

    /**
     * @brief Build a stable local Cartesian-tree BP selector over entries.
     *
     * @details The comparator must return whether the left slot is strictly
     * better than the right slot. Equal minima stay stable through the same
     * right-to-left construction rule used by the value Cartesian tree.
     */
    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kSelectorEntries) {
        throw std::length_error(
            "HybridBTreePlusMinusOne local selector too large");
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

    /**
     * @brief Return the first minimum slot in a local entry range.
     */
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

    /**
     * @brief Store a node's subtree minimum in the unused selector tail.
     *
     * @details This is valid for nodes with at most 192 entries. Their BP
     * sequence uses at most the first 384 bits, and all local selector queries
     * are bounded by `2 * entry_count`, so the final two words are available
     * for node metadata.
     */
    void set_embedded_min_summary(std::size_t position, std::int64_t depth) {
      bp_bits_[kEmbeddedSummaryPositionWord] =
          static_cast<std::uint64_t>(position);
      bp_bits_[kEmbeddedSummaryDepthWord] = std::bit_cast<std::uint64_t>(depth);
    }

    /**
     * @brief Return the embedded subtree-minimum position.
     */
    std::size_t embedded_min_position() const {
      return static_cast<std::size_t>(bp_bits_[kEmbeddedSummaryPositionWord]);
    }

    /**
     * @brief Return the embedded subtree-minimum depth.
     */
    std::int64_t embedded_min_depth() const {
      return std::bit_cast<std::int64_t>(bp_bits_[kEmbeddedSummaryDepthWord]);
    }

   private:
    /**
     * @brief Prepend one BP bit while building the sequence right-to-left.
     */
    void prepend_bp_bit(std::size_t& write_position, bool bit) {
      --write_position;
      if (bit) {
        bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                         << (write_position & 63);
      }
    }

    /**
     * @brief Return the BP close position for a local entry slot.
     */
    std::size_t close_position(std::size_t slot) const {
      return select0_512(bp_bits_.data(), slot);
    }

    /**
     * @brief Count zero bits before @p position in the raw BP sequence.
     */
    std::size_t rank0_at(std::size_t position, std::size_t bit_count) const {
      position = std::min(position, bit_count);
      return position - rank_512(bp_bits_.data(), position);
    }

    /**
     * @brief Return open-minus-close excess before a raw BP position.
     */
    int prefix_excess(std::size_t position) const {
      position = std::min(position, kSelectorBits);
      const std::size_t ones = rank_512(bp_bits_.data(), position);
      return static_cast<int>(2 * ones) - static_cast<int>(position);
    }

    /**
     * @brief Return the minimum BP-depth position in a selector depth range.
     */
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
          candidate_depth = prefix_excess(bit_count);
          candidate_position = bit_count;
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

    std::array<std::uint64_t, kSelectorWords> bp_bits_{};
  };

  static_assert(sizeof(Bp512Selector) == 64);

  /**
   * @brief Return whether a stored position is one of the missing sentinels.
   */
  bool missing_position(std::size_t position) const {
    if constexpr (kInvalidIndexEqualsNpos) {
      return position == npos;
    } else {
      return position == npos ||
             position == static_cast<std::size_t>(invalid_index);
    }
  }

  /**
   * @brief Choose the smaller-depth candidate, breaking ties by position.
   */
  DepthCandidate better_candidate(DepthCandidate left,
                                  DepthCandidate right) const {
    if (missing_position(left.position)) {
      return right;
    }
    if (missing_position(right.position)) {
      return left;
    }
    if (right.depth < left.depth) {
      return right;
    }
    if (left.depth < right.depth) {
      return left;
    }
    return right.position < left.position ? right : left;
  }

  /**
   * @brief Return whether @p left is strictly better than @p right.
   */
  bool strictly_better_candidate(DepthCandidate left,
                                 DepthCandidate right) const {
    if (missing_position(left.position)) {
      return false;
    }
    if (missing_position(right.position)) {
      return true;
    }
    if (left.depth != right.depth) {
      return left.depth < right.depth;
    }
    return left.position < right.position;
  }

  /**
   * @brief Build rank support, internal selectors, and top sparse tables.
   */
  void build() {
    owned_rank_index_.reset();
    internal_selectors_.clear();
    internal_min_positions_.clear();
    internal_min_depths_.clear();
    high_child_metadata_.clear();
    high_sparse_min_slots_.clear();
    internal_level_offsets_.clear();
    min_summary_level_offsets_.clear();
    high_level_offsets_.clear();
    level_sizes_.clear();
    level_position_spans_.clear();
    level_fanouts_.clear();
    high_level_begin_ = std::numeric_limits<std::size_t>::max();

    if (depth_count_ == 0) {
      return;
    }
    if (depth_count_ > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error(
          "HybridBTreePlusMinusOne index type is too small");
    }
    if (depth_count_ >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
      throw std::length_error(
          "HybridBTreePlusMinusOne depth range is too large");
    }
    if (input_bits_.size() < (depth_count_ - 1 + 63) / 64) {
      throw std::invalid_argument(
          "HybridBTreePlusMinusOne bit span is too small");
    }

    const std::size_t delta_count = depth_count_ - 1;
    if (external_rank_index_ == nullptr) {
      owned_rank_index_.emplace(input_bits_, delta_count,
                                BitVector::SelectSupport::kSelect0);
    } else if (external_rank_index_->size() < delta_count) {
      throw std::invalid_argument(
          "HybridBTreePlusMinusOne external rank index is too small");
    }

    initialize_layout((depth_count_ + LeafSize - 1) / LeafSize);
    for (std::size_t level = 1; level < level_count(); ++level) {
      for (std::size_t node = 0; node < level_sizes_[level]; ++node) {
        build_internal_node(level, node);
      }
    }
  }

  /**
   * @brief Compute B-tree level sizes, fanouts, and flat storage offsets.
   */
  void initialize_layout(std::size_t leaf_count) {
    level_sizes_.push_back(leaf_count);
    level_position_spans_.push_back(LeafSize);
    level_fanouts_.push_back(0);

    std::size_t current_count = leaf_count;
    std::size_t current_span = LeafSize;
    while (current_count > kHighLevelFanout * kHighLevelFanout) {
      level_fanouts_.push_back(kMiddleFanout);
      current_count = ceil_div(current_count, kMiddleFanout);
      current_span = saturating_product(current_span, kMiddleFanout);
      level_sizes_.push_back(current_count);
      level_position_spans_.push_back(current_span);
    }
    while (current_count > 1) {
      level_fanouts_.push_back(kHighLevelFanout);
      current_count = ceil_div(current_count, kHighLevelFanout);
      current_span = saturating_product(current_span, kHighLevelFanout);
      level_sizes_.push_back(current_count);
      level_position_spans_.push_back(current_span);
    }
    internal_level_offsets_.assign(level_count(), 0);
    min_summary_level_offsets_.assign(level_count(), npos);
    high_level_offsets_.assign(level_count(), 0);
    if (level_count() <= 1) {
      return;
    }

    std::size_t internal_count = 0;
    for (std::size_t level = 1; level < level_count(); ++level) {
      internal_level_offsets_[level] = internal_count;
      internal_count += level_sizes_[level];
    }
    internal_selectors_.resize(internal_count);

    if constexpr (UseHighSparseLayout) {
      const std::size_t root_level = level_count() - 1;
      std::size_t high_layout_levels = 0;
      for (std::size_t level = root_level;
           level > 0 && high_layout_levels < HighSparseLayoutLevels &&
           fanout_at_level(level) == kHighLevelFanout;
           --level) {
        ++high_layout_levels;
      }
      high_level_begin_ = high_layout_levels == 0
                              ? level_count()
                              : root_level + 1 - high_layout_levels;

      std::size_t high_node_count = 0;
      for (std::size_t level = high_level_begin_; level < level_count();
           ++level) {
        high_level_offsets_[level] = high_node_count;
        high_node_count += level_sizes_[level];
      }
      high_child_metadata_.resize(high_node_count * kHighLevelFanout);
      high_sparse_min_slots_.resize(high_node_count * kHighSparseSlotsPerNode);
    } else {
      high_level_begin_ = level_count();
    }

    std::size_t side_summary_count = 0;
    for (std::size_t level = 1; level < level_count(); ++level) {
      if (!level_embeds_min_summary(level)) {
        min_summary_level_offsets_[level] = side_summary_count;
        side_summary_count += level_sizes_[level];
      }
    }
    internal_min_positions_.resize(side_summary_count, invalid_index);
    internal_min_depths_.resize(side_summary_count,
                                std::numeric_limits<std::int64_t>::max());
  }

  /**
   * @brief Build one internal node selector and cached minimum.
   */
  void build_internal_node(std::size_t level, std::size_t node) {
    const std::size_t count = entry_count(level, node);
    const std::size_t first_child = node * fanout_at_level(level);
    const bool high_level = is_high_level(level);
    const std::size_t high_flat = high_level ? high_flat_index(level, node) : 0;

    std::array<DepthCandidate, kSelectorEntries> child_minima{};
    for (std::size_t slot = 0; slot < count; ++slot) {
      child_minima[slot] = subtree_min_candidate(level - 1, first_child + slot);
    }

    if (high_level) {
      for (std::size_t slot = 0; slot < count; ++slot) {
        const std::size_t child = first_child + slot;
        const DepthCandidate child_min = child_minima[slot];
        HighChildMetadata& metadata =
            mutable_high_child_metadata_at(high_flat, slot);
        metadata.position_begin = node_position_begin(level - 1, child);
        metadata.position_end = node_position_end(level - 1, child);
        metadata.min_position = static_cast<Index>(child_min.position);
        metadata.min_depth = child_min.depth;
      }
      build_high_sparse_min_slots(level, node, count);
    }

    Bp512Selector& selector = mutable_selector_at(level, node);
    selector.build(count, [&](std::size_t left, std::size_t right) {
      return strictly_better_candidate(child_minima[left], child_minima[right]);
    });

    const std::size_t slot = selector.arg_min(0, count, count);
    const DepthCandidate minimum = child_minima[slot];
    if (level_embeds_min_summary(level)) {
      selector.set_embedded_min_summary(minimum.position, minimum.depth);
    } else {
      const std::size_t flat = min_summary_flat_index(level, node);
      internal_min_positions_[flat] = static_cast<Index>(minimum.position);
      internal_min_depths_[flat] = minimum.depth;
    }
  }

  /**
   * @brief Return `ceil(value / divisor)` for positive divisors.
   */
  static std::size_t ceil_div(std::size_t value, std::size_t divisor) {
    return (value + divisor - 1) / divisor;
  }

  /**
   * @brief Multiply two extents, saturating at `std::size_t` maximum.
   */
  static std::size_t saturating_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
      return std::numeric_limits<std::size_t>::max();
    }
    return left * right;
  }

  /**
   * @brief Read an input word or return zero past the available span.
   */
  std::uint64_t word_or_zero(std::size_t word) const {
    return word < input_bits_.size() ? input_bits_[word] : 0;
  }

  /**
   * @brief Return two contiguous words for a 128-bit leaf chunk.
   *
   * @details Most Cartesian BP storage is padded, so the returned pointer can
   * usually address the original input span directly. The temporary storage is
   * used only for the final partial chunk of standalone depth-RMQ indexes.
   */
  const std::uint64_t* chunk_words_or_copy(
      std::size_t first_word,
      std::array<std::uint64_t, 2>& storage) const {
    if (first_word + 1 < input_bits_.size()) {
      return input_bits_.data() + first_word;
    }
    storage[0] = word_or_zero(first_word);
    storage[1] = word_or_zero(first_word + 1);
    return storage.data();
  }

  /**
   * @brief Return the active BP rank/select support, if one exists.
   */
  const BitVector* rank_index_or_null() const {
    return external_rank_index_ != nullptr
               ? external_rank_index_
               : (owned_rank_index_ ? &*owned_rank_index_ : nullptr);
  }

  /**
   * @brief Return the active BP rank/select support.
   */
  const BitVector& rank_index() const { return *rank_index_or_null(); }

  /**
   * @brief Return the absolute open-minus-close depth at a BP depth position.
   */
  std::int64_t depth_at_position(std::size_t position) const {
    const std::size_t delta_count = depth_count_ == 0 ? 0 : depth_count_ - 1;
    position = std::min(position, delta_count);
    const std::uint64_t ones = rank_index().rank(position);
    return static_cast<std::int64_t>(ones) -
           static_cast<std::int64_t>(position - ones);
  }

  /**
   * @brief Return the first minimum candidate inside one leaf local range.
   */
  DepthCandidate scan_leaf_range_with_base(std::size_t leaf,
                                           std::size_t left_offset,
                                           std::size_t right_offset,
                                           std::int64_t base_depth) const {
    const std::size_t begin = node_position_begin(0, leaf);
    const std::size_t count = entry_count(0, leaf);
    if (count == 0 || left_offset > right_offset || left_offset >= count) {
      return {};
    }
    right_offset = std::min(right_offset, count - 1);

    DepthCandidate answer;
    std::int64_t chunk_base_excess = 0;
    const std::size_t first_word = leaf * kLeafWords;
    for (std::size_t chunk = 0; chunk < kLeafChunks; ++chunk) {
      const std::size_t chunk_begin = chunk * 128;
      if (chunk_begin >= count || chunk_begin > right_offset) {
        break;
      }

      std::array<std::uint64_t, 2> chunk_storage{};
      const std::uint64_t* chunk_words =
          chunk_words_or_copy(first_word + 2 * chunk, chunk_storage);

      const std::size_t chunk_end =
          std::min<std::size_t>(count - 1, chunk_begin + 127);
      if (left_offset > chunk_end) {
        chunk_base_excess += prefix_excess_128(chunk_words, 128);
        continue;
      }

      const std::size_t local_left =
          std::max(left_offset, chunk_begin) - chunk_begin;
      const std::size_t local_right =
          std::min(right_offset, chunk_end) - chunk_begin;
      const ExcessResult result =
          excess_min_128(chunk_words, local_left, local_right);
      const std::size_t offset = chunk_begin + result.offset;
      if (result.offset != npos && offset < count) {
        answer = better_candidate(
            answer, {begin + offset,
                     base_depth + chunk_base_excess + result.min_excess});
      }

      chunk_base_excess += prefix_excess_128(chunk_words, 128);
    }
    return answer;
  }

  /**
   * @brief Return the first minimum position in one leaf without global rank.
   *
   * @details This is the same-leaf fast path. It compares prefix depths after
   * subtracting the depth at the query's left endpoint, which preserves the
   * arg-min position while avoiding an absolute `rank()` query.
   */
  std::size_t leaf_range_arg_min_relative(std::size_t leaf,
                                          std::size_t left,
                                          std::size_t right) const {
    if (left >= right) {
      return npos;
    }

    const std::size_t begin = node_position_begin(0, leaf);
    const std::size_t count = entry_count(0, leaf);
    std::size_t left_offset = left - begin;
    std::size_t right_offset = right - begin - 1;
    if (count == 0 || left_offset >= count || left_offset > right_offset) {
      return npos;
    }
    right_offset = std::min(right_offset, count - 1);

    std::size_t best_position = npos;
    std::int64_t best_depth = std::numeric_limits<std::int64_t>::max();
    std::int64_t chunk_base_excess = 0;
    bool first_chunk = true;

    const std::size_t first_word = leaf * kLeafWords;
    std::size_t chunk = left_offset / 128;
    for (; chunk < kLeafChunks; ++chunk) {
      const std::size_t chunk_begin = chunk * 128;
      if (chunk_begin >= count || chunk_begin > right_offset) {
        break;
      }

      std::array<std::uint64_t, 2> chunk_storage{};
      const std::uint64_t* chunk_words =
          chunk_words_or_copy(first_word + 2 * chunk, chunk_storage);

      const std::size_t chunk_end =
          std::min<std::size_t>(count - 1, chunk_begin + 127);
      const std::size_t local_left =
          std::max(left_offset, chunk_begin) - chunk_begin;
      const std::size_t local_right =
          std::min(right_offset, chunk_end) - chunk_begin;
      const int left_prefix =
          first_chunk ? prefix_excess_128(chunk_words, local_left) : 0;
      const std::int64_t local_base =
          first_chunk ? -static_cast<std::int64_t>(left_prefix)
                      : chunk_base_excess;
      const ExcessResult result =
          excess_min_128(chunk_words, local_left, local_right);
      const std::size_t offset = chunk_begin + result.offset;
      if (result.offset != npos && offset < count) {
        const std::int64_t candidate_depth = local_base + result.min_excess;
        if (best_position == npos || candidate_depth < best_depth) {
          best_position = begin + offset;
          best_depth = candidate_depth;
        }
      }

      const int chunk_excess = prefix_excess_128(chunk_words, 128);
      if (first_chunk) {
        chunk_base_excess =
            static_cast<std::int64_t>(chunk_excess) - left_prefix;
        first_chunk = false;
      } else {
        chunk_base_excess += chunk_excess;
      }
    }

    return best_position;
  }

  /**
   * @brief Return the minimum candidate in a possibly partial leaf.
   */
  DepthCandidate leaf_range_min(std::size_t leaf,
                                std::size_t left,
                                std::size_t right) const {
    if (left >= right) {
      return {};
    }

    const std::size_t begin = node_position_begin(0, leaf);
    const std::size_t slot_left = left - begin;
    const std::size_t slot_right = right - begin;
    return scan_leaf_range_with_base(leaf, slot_left, slot_right - 1,
                                     depth_at_position(begin));
  }

  /**
   * @brief Find the lowest tree node that contains both endpoint leaves.
   */
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

  /**
   * @brief Return the leaf index containing a depth position.
   */
  std::size_t leaf_for_position(std::size_t position) const {
    return position / LeafSize;
  }

  /**
   * @brief Return the child index at @p child_level containing a position.
   */
  std::size_t child_for_position(std::size_t child_level,
                                 std::size_t position) const {
    return position / level_position_spans_[child_level];
  }

  /**
   * @brief Query a child-slot range of an internal node.
   */
  DepthCandidate query_child_slots(std::size_t level,
                                   std::size_t node,
                                   std::size_t slot_left,
                                   std::size_t slot_right,
                                   std::size_t left,
                                   std::size_t right) const {
    if (slot_left >= slot_right) {
      return {};
    }

    const std::size_t count = entry_count(level, node);
    const std::size_t slot =
        slot_left + 1 == slot_right
            ? slot_left
            : selector_arg_min(level, node, slot_left, slot_right, count);
    if (slot == npos) {
      return {};
    }

    const std::size_t child_level = level - 1;
    const std::size_t first_child = node * fanout_at_level(level);
    const std::size_t child = first_child + slot;
    const HighChildMetadata* high_children =
        is_high_level(level) ? high_child_metadata_begin(level, node) : nullptr;
    const DepthCandidate child_min =
        high_children != nullptr ? high_child_min_candidate(level, node, slot)
                                 : subtree_min_candidate(child_level, child);
    const std::size_t child_begin =
        high_children != nullptr ? high_children[slot].position_begin
                                 : node_position_begin(child_level, child);
    const std::size_t child_end = high_children != nullptr
                                      ? high_children[slot].position_end
                                      : node_position_end(child_level, child);
    if ((left <= child_begin && child_end <= right) ||
        contains_position(left, right, child_min.position)) {
      return child_min;
    }

    const std::size_t last_slot = slot_right - 1;
    const std::size_t left_child_begin =
        high_children != nullptr
            ? high_children[slot_left].position_begin
            : node_position_begin(child_level, first_child + slot_left);
    const std::size_t left_child_end =
        high_children != nullptr
            ? high_children[slot_left].position_end
            : node_position_end(child_level, first_child + slot_left);
    DepthCandidate answer = query_node(child_level, first_child + slot_left,
                                       std::max(left, left_child_begin),
                                       std::min(right, left_child_end));

    if (slot_left != last_slot) {
      const std::size_t right_child_begin =
          high_children != nullptr
              ? high_children[last_slot].position_begin
              : node_position_begin(child_level, first_child + last_slot);
      const std::size_t right_child_end =
          high_children != nullptr
              ? high_children[last_slot].position_end
              : node_position_end(child_level, first_child + last_slot);
      answer = better_candidate(answer,
                                query_node(child_level, first_child + last_slot,
                                           std::max(left, right_child_begin),
                                           std::min(right, right_child_end)));
    }

    if (slot_left + 1 < last_slot) {
      answer = better_candidate(
          answer,
          full_child_slot_range_min(level, node, slot_left + 1, last_slot));
    }

    return answer;
  }

  /**
   * @brief Return the best candidate among fully covered child slots.
   */
  DepthCandidate full_child_slot_range_min(std::size_t level,
                                           std::size_t node,
                                           std::size_t slot_left,
                                           std::size_t slot_right) const {
    if (slot_left >= slot_right) {
      return {};
    }

    const std::size_t slot =
        slot_left + 1 == slot_right
            ? slot_left
            : selector_arg_min(level, node, slot_left, slot_right,
                               entry_count(level, node));
    if (slot == npos) {
      return {};
    }
    return child_min_candidate(level, node, slot);
  }

  /**
   * @brief Query a tree node for the minimum candidate in a depth range.
   */
  DepthCandidate query_node(std::size_t level,
                            std::size_t node,
                            std::size_t left,
                            std::size_t right) const {
    if (left >= right) {
      return {};
    }

    const std::size_t begin = node_position_begin(level, node);
    const std::size_t end = node_position_end(level, node);
    if (left <= begin && end <= right) {
      return subtree_min_candidate(level, node);
    }
    if (level == 0) {
      return leaf_range_min(node, left, right);
    }

    const std::size_t child_level = level - 1;
    const std::size_t left_child = child_for_position(child_level, left);
    const std::size_t right_child = child_for_position(child_level, right - 1);
    const std::size_t first_child = node * fanout_at_level(level);
    const std::size_t left_slot = left_child - first_child;
    const std::size_t right_slot = right_child - first_child + 1;
    return query_child_slots(level, node, left_slot, right_slot, left, right);
  }

  /**
   * @brief Return whether a position lies inside a half-open range.
   */
  bool contains_position(std::size_t left,
                         std::size_t right,
                         std::size_t position) const {
    return !missing_position(position) && left <= position && position < right;
  }

  /**
   * @brief Return the number of B-tree levels, including leaves.
   */
  std::size_t level_count() const { return level_sizes_.size(); }

  /**
   * @brief Return the number of entries in a leaf or child slots in a node.
   */
  std::size_t entry_count(std::size_t level, std::size_t node) const {
    if (level == 0) {
      const std::size_t begin = node_position_begin(0, node);
      return std::min<std::size_t>(LeafSize, depth_count_ - begin);
    }
    const std::size_t first_child = node * fanout_at_level(level);
    return std::min<std::size_t>(fanout_at_level(level),
                                 level_sizes_[level - 1] - first_child);
  }

  /**
   * @brief Return the first depth position covered by a node.
   */
  std::size_t node_position_begin(std::size_t level, std::size_t node) const {
    return node * level_position_spans_[level];
  }

  /**
   * @brief Return one past the last depth position covered by a node.
   */
  std::size_t node_position_end(std::size_t level, std::size_t node) const {
    return std::min(depth_count_, node_position_begin(level, node) +
                                      level_position_spans_[level]);
  }

  /**
   * @brief Return a node's cached subtree-minimum candidate.
   */
  DepthCandidate subtree_min_candidate(std::size_t level,
                                       std::size_t node) const {
    if (level == 0) {
      return leaf_range_min(node, node_position_begin(0, node),
                            node_position_end(0, node));
    }
    if (level_embeds_min_summary(level)) {
      const Bp512Selector& selector = selector_at(level, node);
      return {selector.embedded_min_position(), selector.embedded_min_depth()};
    }
    const std::size_t flat = min_summary_flat_index(level, node);
    return {static_cast<std::size_t>(internal_min_positions_[flat]),
            internal_min_depths_[flat]};
  }

  /**
   * @brief Return a child slot's subtree-minimum candidate.
   */
  DepthCandidate child_min_candidate(std::size_t level,
                                     std::size_t node,
                                     std::size_t slot) const {
    if (is_high_level(level)) {
      return high_child_min_candidate(level, node, slot);
    }
    return subtree_min_candidate(level - 1,
                                 node * fanout_at_level(level) + slot);
  }

  /**
   * @brief Return a high-node child candidate using cached metadata.
   */
  DepthCandidate high_child_min_candidate(std::size_t level,
                                          std::size_t node,
                                          std::size_t slot) const {
    const HighChildMetadata& metadata =
        high_child_metadata_at(high_flat_index(level, node), slot);
    return {static_cast<std::size_t>(metadata.min_position),
            metadata.min_depth};
  }

  /**
   * @brief Return an immutable internal-node BP selector.
   */
  const Bp512Selector& selector_at(std::size_t level, std::size_t node) const {
    return internal_selectors_[internal_flat_index(level, node)];
  }

  /**
   * @brief Return a mutable internal-node BP selector while building.
   */
  Bp512Selector& mutable_selector_at(std::size_t level, std::size_t node) {
    return internal_selectors_[internal_flat_index(level, node)];
  }

  /**
   * @brief Run the appropriate local selector for an internal node.
   */
  std::size_t selector_arg_min(std::size_t level,
                               std::size_t node,
                               std::size_t slot_left,
                               std::size_t slot_right,
                               std::size_t count) const {
    if constexpr (UseHighSparseLayout) {
      if (is_high_level(level)) {
        return high_sparse_arg_min(level, node, slot_left, slot_right, count);
      }
    }
    return selector_at(level, node).arg_min(slot_left, slot_right, count);
  }

  /**
   * @brief Return whether a level uses the high-node layout.
   */
  bool is_high_level(std::size_t level) const {
    if constexpr (!UseHighSparseLayout) {
      (void)level;
      return false;
    }
    return level > 0 && level >= high_level_begin_ && level < level_count();
  }

  /**
   * @brief Map an internal node to its flat storage index.
   */
  std::size_t internal_flat_index(std::size_t level, std::size_t node) const {
    return internal_level_offsets_[level] + node;
  }

  /**
   * @brief Return whether a level embeds subtree minima in selector tail bits.
   */
  bool level_embeds_min_summary(std::size_t level) const {
    return level > 0 && !is_high_level(level) &&
           fanout_at_level(level) <= kEmbeddedSummaryMaxEntries;
  }

  /**
   * @brief Map a non-embedded internal node to side summary storage.
   */
  std::size_t min_summary_flat_index(std::size_t level,
                                     std::size_t node) const {
    return min_summary_level_offsets_[level] + node;
  }

  /**
   * @brief Return the fanout used to group children at a level.
   */
  std::size_t fanout_at_level(std::size_t level) const {
    return level_fanouts_[level];
  }

  /**
   * @brief Map a high-level node to its flat high-node storage index.
   */
  std::size_t high_flat_index(std::size_t level, std::size_t node) const {
    return high_level_offsets_[level] + node;
  }

  /**
   * @brief Return the first child-metadata record for a high node.
   */
  const HighChildMetadata* high_child_metadata_begin(std::size_t level,
                                                     std::size_t node) const {
    return high_child_metadata_.data() +
           high_flat_index(level, node) * kHighLevelFanout;
  }

  /**
   * @brief Return high-child metadata by flat high-node index and slot.
   */
  const HighChildMetadata& high_child_metadata_at(std::size_t high_flat,
                                                  std::size_t slot) const {
    return high_child_metadata_[high_flat * kHighLevelFanout + slot];
  }

  /**
   * @brief Return mutable high-child metadata while building.
   */
  HighChildMetadata& mutable_high_child_metadata_at(std::size_t high_flat,
                                                    std::size_t slot) {
    return high_child_metadata_[high_flat * kHighLevelFanout + slot];
  }

  /**
   * @brief Return mutable sparse-slot storage for one high node.
   */
  std::uint8_t* mutable_high_sparse_min_slots_begin(std::size_t high_flat) {
    return high_sparse_min_slots_.data() + high_flat * kHighSparseSlotsPerNode;
  }

  /**
   * @brief Return sparse-slot storage for one high node.
   */
  const std::uint8_t* high_sparse_min_slots_begin(std::size_t high_flat) const {
    return high_sparse_min_slots_.data() + high_flat * kHighSparseSlotsPerNode;
  }

  /**
   * @brief Choose the better high-node child slot.
   */
  std::size_t better_high_child_slot(std::size_t level,
                                     std::size_t node,
                                     std::size_t left_slot,
                                     std::size_t right_slot) const {
    const DepthCandidate left =
        high_child_min_candidate(level, node, left_slot);
    const DepthCandidate right =
        high_child_min_candidate(level, node, right_slot);
    return better_candidate(left, right).position == right.position ? right_slot
                                                                    : left_slot;
  }

  /**
   * @brief Build sparse tables over high-node child minima.
   */
  void build_high_sparse_min_slots(std::size_t level,
                                   std::size_t node,
                                   std::size_t count) {
    const std::size_t high_flat = high_flat_index(level, node);
    std::uint8_t* table = mutable_high_sparse_min_slots_begin(high_flat);
    for (std::size_t slot = 0; slot < count; ++slot) {
      table[slot] = static_cast<std::uint8_t>(slot);
    }

    for (std::size_t table_level = 1; table_level < kHighSparseTableLevels;
         ++table_level) {
      const std::size_t span = std::size_t{1} << table_level;
      if (span > count) {
        break;
      }
      const std::size_t half_span = span >> 1;
      const std::uint8_t* previous =
          table + (table_level - 1) * kHighLevelFanout;
      std::uint8_t* current = table + table_level * kHighLevelFanout;
      for (std::size_t slot = 0; slot + span <= count; ++slot) {
        current[slot] = static_cast<std::uint8_t>(better_high_child_slot(
            level, node, previous[slot], previous[slot + half_span]));
      }
    }
  }

  /**
   * @brief Return the best high-node child slot in a slot range.
   */
  std::size_t high_sparse_arg_min(std::size_t level,
                                  std::size_t node,
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

    const std::size_t high_flat = high_flat_index(level, node);
    const std::size_t table_level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << table_level;
    const std::uint8_t* table =
        high_sparse_min_slots_begin(high_flat) + table_level * kHighLevelFanout;
    return better_high_child_slot(level, node, table[slot_left],
                                  table[slot_right - span]);
  }

  std::span<const std::uint64_t> input_bits_;
  std::size_t depth_count_ = 0;
  std::optional<BitVector> owned_rank_index_;
  const BitVector* external_rank_index_ = nullptr;
  std::vector<Bp512Selector> internal_selectors_;
  std::vector<Index> internal_min_positions_;
  std::vector<std::int64_t> internal_min_depths_;
  std::vector<HighChildMetadata> high_child_metadata_;
  std::vector<std::uint8_t> high_sparse_min_slots_;
  std::vector<std::size_t> internal_level_offsets_;
  std::vector<std::size_t> min_summary_level_offsets_;
  std::vector<std::size_t> high_level_offsets_;
  std::vector<std::size_t> level_sizes_;
  std::vector<std::size_t> level_position_spans_;
  std::vector<std::size_t> level_fanouts_;
  std::size_t high_level_begin_ = std::numeric_limits<std::size_t>::max();
};

}  // namespace detail

/**
 * @brief Cartesian-tree value RMQ using HybridBTree-style LCA.
 *
 * @details This class follows the same public value-RMQ specification as the
 * other value RMQ backends. It builds a stable Ferrada-Navarro BP
 * Cartesian-tree encoding, uses `BitVector` for close-parenthesis rank/select,
 * and delegates the BP-depth minimum query to
 * `detail::HybridBTreePlusMinusOne`. The BP-depth backend keeps a configurable
 * low-level leaf size, fixed 192-entry middle nodes with embedded minima, and
 * fixed 256-entry high nodes. A single coarse value-level sparse table is
 * checked first; it uses at least 4096-value blocks and grows the block width
 * when needed so the top layer has at most 2^14 blocks. Wide queries whose
 * padded block-cover minimum lies inside the requested range return from this
 * top table without touching the global BP rank/select path. BP construction
 * uses a succinct monotone bit-stack, preserving the same stable Cartesian-tree
 * shape without an n-entry index stack.
 *
 * This implementation is included from `pixie/rmq.h` as the compact
 * Cartesian-tree reduction backed by a HybridBTree-shaped ±1 RMQ index.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 512,
          bool UseTopSparseOverlay = true>
class CartesianHybridBTree
    : public RmqBase<CartesianHybridBTree<T,
                                          Compare,
                                          Index,
                                          LeafSize,
                                          UseTopSparseOverlay>,
                     T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "CartesianHybridBTree index type must be unsigned");
  static_assert(LeafSize != 0 && LeafSize % 512 == 0,
                "CartesianHybridBTree leaf size must be a positive "
                "multiple of 512");

  using Self =
      CartesianHybridBTree<T, Compare, Index, LeafSize, UseTopSparseOverlay>;

  static constexpr std::size_t npos = RmqBase<Self, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kMinTopSparseBlockSize = 4096;
  static constexpr std::size_t kMaxTopSparseBlocks = std::size_t{1} << 14;
  static constexpr bool kUseTopSparseOverlay = UseTopSparseOverlay;

  /**
   * @brief Construct an empty Cartesian-tree RMQ index.
   */
  CartesianHybridBTree() = default;

  /**
   * @brief Build a Cartesian-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values stay stable: the smaller index remains the first minimum.
   */
  explicit CartesianHybridBTree(std::span<const T> values,
                                Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  /**
   * @brief Copy an RMQ index and rebuild internal non-owning views.
   */
  CartesianHybridBTree(const CartesianHybridBTree& other)
      : values_(other.values_),
        compare_(other.compare_),
        bp_bits_(other.bp_bits_),
        bp_bit_count_(other.bp_bit_count_),
        top_sparse_candidates_(other.top_sparse_candidates_),
        top_block_size_(other.top_block_size_),
        top_block_count_(other.top_block_count_),
        top_sparse_levels_(other.top_sparse_levels_) {
    reset_bp_indexes();
  }

  /**
   * @brief Copy-assign an RMQ index and rebuild internal non-owning views.
   */
  CartesianHybridBTree& operator=(const CartesianHybridBTree& other) {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = other.compare_;
    bp_bits_ = other.bp_bits_;
    bp_bit_count_ = other.bp_bit_count_;
    top_sparse_candidates_ = other.top_sparse_candidates_;
    top_block_size_ = other.top_block_size_;
    top_block_count_ = other.top_block_count_;
    top_sparse_levels_ = other.top_sparse_levels_;
    reset_bp_indexes();
    return *this;
  }

  /**
   * @brief Move an RMQ index and rebuild internal non-owning views.
   */
  CartesianHybridBTree(CartesianHybridBTree&& other) noexcept
      : values_(other.values_),
        compare_(std::move(other.compare_)),
        bp_bits_(std::move(other.bp_bits_)),
        bp_bit_count_(other.bp_bit_count_),
        top_sparse_candidates_(std::move(other.top_sparse_candidates_)),
        top_block_size_(other.top_block_size_),
        top_block_count_(other.top_block_count_),
        top_sparse_levels_(other.top_sparse_levels_) {
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    other.top_block_size_ = kMinTopSparseBlockSize;
    other.top_block_count_ = 0;
    other.top_sparse_levels_ = 0;
    reset_bp_indexes();
  }

  /**
   * @brief Move-assign an RMQ index and rebuild internal non-owning views.
   */
  CartesianHybridBTree& operator=(CartesianHybridBTree&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    values_ = other.values_;
    compare_ = std::move(other.compare_);
    bp_bits_ = std::move(other.bp_bits_);
    bp_bit_count_ = other.bp_bit_count_;
    top_sparse_candidates_ = std::move(other.top_sparse_candidates_);
    top_block_size_ = other.top_block_size_;
    top_block_count_ = other.top_block_count_;
    top_sparse_levels_ = other.top_sparse_levels_;
    other.values_ = std::span<const T>();
    other.bp_bit_count_ = 0;
    other.top_block_size_ = kMinTopSparseBlockSize;
    other.top_block_count_ = 0;
    other.top_sparse_levels_ = 0;
    reset_bp_indexes();
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
    if constexpr (!UseTopSparseOverlay) {
      return cartesian_arg_min(left, right);
    }
    if (right - left <= top_block_size_) {
      return cartesian_arg_min(left, right);
    }
    const std::size_t top_answer = top_sparse_arg_min(left, right);
    if (top_answer != npos) {
      return top_answer;
    }
    return cartesian_arg_min(left, right);
  }

  /**
   * @brief Return the number of BP bits in the Cartesian-tree RMQ encoding.
   */
  std::size_t bp_bit_count() const { return bp_bit_count_; }

  /**
   * @brief Return the packed BP words used by the RMQ encoding.
   */
  std::span<const std::uint64_t> bp_words() const {
    return bp_storage_words().first(bp_word_count());
  }

  /**
   * @brief Return the top sparse-table block width chosen for a value count.
   */
  static std::size_t top_sparse_block_size_for(std::size_t value_count) {
    if (value_count == 0) {
      return kMinTopSparseBlockSize;
    }
    return std::max(kMinTopSparseBlockSize,
                    ceil_div(value_count, kMaxTopSparseBlocks));
  }

  /**
   * @brief Return the number of top sparse-table blocks for a value count.
   */
  static std::size_t top_sparse_block_count_for(std::size_t value_count) {
    if (value_count == 0) {
      return 0;
    }
    return ceil_div(value_count, top_sparse_block_size_for(value_count));
  }

  /**
   * @brief Return the current top sparse-table block width.
   */
  std::size_t top_sparse_block_size() const { return top_block_size_; }

  /**
   * @brief Return the current number of top sparse-table blocks.
   */
  std::size_t top_sparse_block_count() const { return top_block_count_; }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this value-RMQ object, packed Cartesian BP words, the top
   * sparse overlay, and nested BP rank/select and ±1 RMQ indexes. The external
   * input values are not owned and are excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    return sizeof(*this) + bp_bits_.allocated_bytes() +
           pixie::vector_capacity_bytes(top_sparse_candidates_) +
           pixie::optional_nested_owned_memory_bytes(bp_index_) +
           pixie::nested_owned_memory_bytes(bp_depth_rmq_);
  }

 private:
  using BpDepthRmq = detail::HybridBTreePlusMinusOne<Index, LeafSize, true, 1>;

  struct TopCandidate {
    Index position = invalid_index;
  };
  static_assert(sizeof(TopCandidate) == sizeof(Index));

  /**
   * @brief Return the first minimum position through the Cartesian BP
   * reduction.
   */
  std::size_t cartesian_arg_min(std::size_t left, std::size_t right) const {
    const std::size_t first_close = select_close_position(left + 1);
    const std::size_t last_close = select_close_position(right);
    if (first_close == npos || last_close == npos || first_close > last_close) {
      return npos;
    }

    const std::size_t shifted_min =
        bp_depth_rmq_.arg_min(first_close + 1, last_close + 2);
    if (shifted_min == npos || shifted_min == 0) {
      return npos;
    }
    const BitVector& bp_index = *bp_index_;
    const std::size_t answer = bp_index.rank0(shifted_min) - 1;
    return answer < values_.size() ? answer : npos;
  }

  /**
   * @brief Rebuild the BP Cartesian-tree representation and support indexes.
   */
  void build() {
    bp_bits_.resize(0);
    bp_bit_count_ = 0;
    top_sparse_candidates_.clear();
    top_block_size_ = kMinTopSparseBlockSize;
    top_block_count_ = 0;
    top_sparse_levels_ = 0;
    reset_bp_indexes();

    if (values_.empty()) {
      return;
    }
    if (values_.size() > (static_cast<std::size_t>(invalid_index) - 1) / 2) {
      throw std::length_error(
          "CartesianHybridBTree RMQ index type is too small");
    }

    bp_bit_count_ = 2 * values_.size();
    bp_bits_.resize(padded_bp_bit_capacity());
    std::ranges::fill(bp_storage_words(), std::uint64_t{0});
    build_bp_bits();
    if constexpr (UseTopSparseOverlay) {
      build_top_sparse_table();
    }
    reset_bp_indexes();
  }

  /**
   * @brief Build the Ferrada-Navarro BP bits with a monotone stack.
   */
  void build_bp_bits() {
    utils::SuccinctIncreasingStack stack(values_.size());
    std::size_t write_position = bp_bit_count_;
    std::span<std::uint64_t> words = bp_storage_words();
    const auto prepend_open = [&]() {
      --write_position;
      words[write_position >> 6] |= std::uint64_t{1} << (write_position & 63);
    };

    for (std::size_t i = values_.size(); i-- > 0;) {
      while (!stack.empty() &&
             !compare_(values_[stack_index(stack.top())], values_[i])) {
        stack.pop();
        prepend_open();
      }
      stack.push(stack_key(i));
      --write_position;
    }

    while (write_position != 0) {
      prepend_open();
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

  /**
   * @brief Rebuild indexes that store non-owning spans into `bp_bits_`.
   */
  void reset_bp_indexes() {
    bp_index_.reset();
    bp_depth_rmq_ = BpDepthRmq();
    if (bp_bit_count_ == 0) {
      return;
    }
    const std::span<const std::uint64_t> words = bp_words();
    const std::span<const std::uint64_t> padded_words = bp_storage_words();
    // TODO: try incorporating rank/select information into the tree.
    bp_index_.emplace(words, bp_bit_count_, BitVector::SelectSupport::kSelect0,
                      values_.size());
    bp_depth_rmq_ = BpDepthRmq(padded_words, bp_bit_count_ + 1, *bp_index_);
  }

  /**
   * @brief Return the exact number of 64-bit words in the BP encoding.
   */
  std::size_t bp_word_count() const { return ceil_div(bp_bit_count_, 64); }

  /**
   * @brief Return BP storage capacity padded to full low-level depth leaves.
   */
  std::size_t padded_bp_bit_capacity() const {
    if (bp_bit_count_ == 0) {
      return 0;
    }
    const std::size_t depth_count = bp_bit_count_ + 1;
    return ceil_div(depth_count, LeafSize) * LeafSize;
  }

  /**
   * @brief Return mutable padded BP storage as 64-bit words.
   */
  std::span<std::uint64_t> bp_storage_words() { return bp_bits_.As64BitInts(); }

  /**
   * @brief Return padded BP storage as 64-bit words.
   */
  std::span<const std::uint64_t> bp_storage_words() const {
    return bp_bits_.AsConst64BitInts();
  }

  /**
   * @brief Build a single top sparse table over original-value block minima.
   */
  void build_top_sparse_table() {
    top_sparse_candidates_.clear();
    top_block_size_ = top_sparse_block_size_for(values_.size());
    top_block_count_ = ceil_div(values_.size(), top_block_size_);
    top_sparse_levels_ =
        top_block_count_ == 0 ? 0 : std::bit_width(top_block_count_);
    if (top_block_count_ == 0) {
      return;
    }

    top_sparse_candidates_.assign(top_sparse_levels_ * top_block_count_,
                                  TopCandidate{});
    for (std::size_t block = 0; block < top_block_count_; ++block) {
      const std::size_t begin = block * top_block_size_;
      const std::size_t end = std::min(values_.size(), begin + top_block_size_);
      std::size_t minimum = begin;
      for (std::size_t position = begin + 1; position < end; ++position) {
        if (strictly_better_value_position(position, minimum)) {
          minimum = position;
        }
      }
      top_sparse_candidates_[block] = make_top_candidate(minimum);
    }

    for (std::size_t level = 1; level < top_sparse_levels_; ++level) {
      const std::size_t span = std::size_t{1} << level;
      const std::size_t half_span = span >> 1;
      TopCandidate* current =
          top_sparse_candidates_.data() + level * top_block_count_;
      const TopCandidate* previous =
          top_sparse_candidates_.data() + (level - 1) * top_block_count_;
      for (std::size_t block = 0; block + span <= top_block_count_; ++block) {
        current[block] =
            better_top_candidate(previous[block], previous[block + half_span]);
      }
    }
  }

  /**
   * @brief Return `ceil(value / divisor)` for positive divisors.
   */
  static std::size_t ceil_div(std::size_t value, std::size_t divisor) {
    return value == 0 ? 0 : 1 + (value - 1) / divisor;
  }

  /**
   * @brief Wrap an original value position as a top sparse-table candidate.
   */
  TopCandidate make_top_candidate(std::size_t position) const {
    if (position >= values_.size()) {
      return {};
    }
    return {static_cast<Index>(position)};
  }

  /**
   * @brief Return whether @p position is a valid stored original-value index.
   */
  bool valid_value_position(std::size_t position) const {
    return position != npos &&
           position != static_cast<std::size_t>(invalid_index) &&
           position < values_.size();
  }

  /**
   * @brief Return whether value position @p left is strictly better than @p
   * right.
   */
  bool strictly_better_value_position(std::size_t left,
                                      std::size_t right) const {
    if (!valid_value_position(left)) {
      return false;
    }
    if (!valid_value_position(right)) {
      return true;
    }
    if (compare_(values_[left], values_[right])) {
      return true;
    }
    if (compare_(values_[right], values_[left])) {
      return false;
    }
    return left < right;
  }

  /**
   * @brief Choose the better original-value candidate, preserving first ties.
   */
  TopCandidate better_top_candidate(TopCandidate left,
                                    TopCandidate right) const {
    const std::size_t left_position = static_cast<std::size_t>(left.position);
    const std::size_t right_position = static_cast<std::size_t>(right.position);
    return strictly_better_value_position(right_position, left_position) ? right
                                                                         : left;
  }

  /**
   * @brief Return the sparse-table candidate over a top-block range.
   */
  TopCandidate top_sparse_block_arg_min(std::size_t block_left,
                                        std::size_t block_right) const {
    if (block_left >= block_right || block_right > top_block_count_ ||
        top_sparse_levels_ == 0) {
      return {};
    }
    const std::size_t length = block_right - block_left;
    const std::size_t level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << level;
    const TopCandidate* table =
        top_sparse_candidates_.data() + level * top_block_count_;
    return better_top_candidate(table[block_left], table[block_right - span]);
  }

  /**
   * @brief Return whether a candidate lies inside a half-open value range.
   */
  bool top_candidate_inside(TopCandidate candidate,
                            std::size_t left,
                            std::size_t right) const {
    const std::size_t position = static_cast<std::size_t>(candidate.position);
    return valid_value_position(position) && left <= position &&
           position < right;
  }

  /**
   * @brief Return the top-overlay answer, or `npos` when the BP path should
   * run.
   */
  std::size_t top_sparse_arg_min(std::size_t left, std::size_t right) const {
    if (top_block_count_ <= 1) {
      return npos;
    }

    const std::size_t padded_block_left = left / top_block_size_;
    const std::size_t padded_block_right = (right - 1) / top_block_size_ + 1;
    if (padded_block_left + 1 >= padded_block_right) {
      return npos;
    }

    const TopCandidate padded =
        top_sparse_block_arg_min(padded_block_left, padded_block_right);
    if (top_candidate_inside(padded, left, right)) {
      return static_cast<std::size_t>(padded.position);
    }

    const std::size_t first_full_block =
        (left + top_block_size_ - 1) / top_block_size_;
    const std::size_t full_block_right = right / top_block_size_;
    if (first_full_block >= full_block_right) {
      return npos;
    }

    TopCandidate answer =
        top_sparse_block_arg_min(first_full_block, full_block_right);

    const std::size_t left_border_end = first_full_block * top_block_size_;
    if (left < left_border_end) {
      answer = better_top_candidate(
          answer, make_top_candidate(cartesian_arg_min(left, left_border_end)));
    }

    const std::size_t right_border_begin = full_block_right * top_block_size_;
    if (right_border_begin < right) {
      answer = better_top_candidate(
          answer,
          make_top_candidate(cartesian_arg_min(right_border_begin, right)));
    }

    return valid_value_position(static_cast<std::size_t>(answer.position))
               ? static_cast<std::size_t>(answer.position)
               : npos;
  }

  /**
   * @brief Return the one-based rank-th Cartesian close parenthesis.
   */
  std::size_t select_close_position(std::size_t rank) const {
    if (rank == 0 || rank > values_.size() || !bp_index_) {
      return npos;
    }
    const std::size_t position = bp_index_->select0(rank);
    return position < bp_bit_count_ ? position : npos;
  }

  std::span<const T> values_;
  Compare compare_;
  pixie::AlignedStorage bp_bits_;
  std::size_t bp_bit_count_ = 0;
  std::vector<TopCandidate> top_sparse_candidates_;
  std::size_t top_block_size_ = kMinTopSparseBlockSize;
  std::size_t top_block_count_ = 0;
  std::size_t top_sparse_levels_ = 0;
  std::optional<BitVector> bp_index_;
  BpDepthRmq bp_depth_rmq_;
};

/**
 * @brief Cartesian-tree RMQ variant without the value-level top sparse overlay.
 *
 * @details This alias keeps the same Cartesian BP encoding and BP-depth B-tree
 * backend as `CartesianHybridBTree`, but every query goes directly through the
 * Cartesian-tree reduction. It is useful as a stable benchmark reference for
 * measuring the value-level top sparse table separately.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 512>
using CartesianBTree = CartesianHybridBTree<T, Compare, Index, LeafSize, false>;

}  // namespace pixie::rmq
