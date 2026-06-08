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

namespace pixie::rmq {

/**
 * @brief Value RMQ with per-node Cartesian/Euler selectors.
 *
 * @details Values are split into fixed-size leaves, then grouped by a B-tree.
 * Every leaf stores a local Cartesian BP selector over original values. Every
 * internal node stores the same kind of selector over child subtree minima.
 * Queries first enter the lowest B-tree node covering both endpoints, then ask
 * that node's selector over all intersecting children, including partial border
 * children. If the selected child is fully covered, or its stored subtree
 * minimum is inside the query, the answer is known immediately. Otherwise the
 * selected border child is corrected recursively and compared with the
 * remaining child slots.
 *
 * This backend is available as the main per-node Euler B-tree RMQ
 * implementation.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 * @tparam LeafSize Number of original values per leaf.
 * @tparam Fanout Maximum number of children per internal node.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 256,
          std::size_t Fanout = 256>
class NodeEulerBTreeRmq
    : public RmqBase<NodeEulerBTreeRmq<T, Compare, Index, LeafSize, Fanout>,
                     T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "NodeEulerBTreeRmq index type must be unsigned");
  static_assert(LeafSize > 0);
  static_assert(Fanout > 1);

  static constexpr std::size_t npos =
      RmqBase<NodeEulerBTreeRmq<T, Compare, Index, LeafSize, Fanout>, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kLeafSize = LeafSize;
  static constexpr std::size_t kFanout = Fanout;

  /**
   * @brief Construct an empty RMQ index.
   */
  NodeEulerBTreeRmq() = default;
  NodeEulerBTreeRmq(const NodeEulerBTreeRmq&) = default;
  NodeEulerBTreeRmq(NodeEulerBTreeRmq&&) noexcept = default;
  NodeEulerBTreeRmq& operator=(const NodeEulerBTreeRmq&) = default;
  NodeEulerBTreeRmq& operator=(NodeEulerBTreeRmq&&) noexcept = default;

  /**
   * @brief Build a B-tree RMQ index over @p values.
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
    if (left >= right || right > values_.size() || nodes_.empty()) {
      return npos;
    }

    const NodeMetadata& root = nodes_.front();
    if (left == root.value_begin && right == root.value_end) {
      return root.subtree_min_position;
    }

    const bool allow_leaf_linear_scan =
        right - left <= kLeafLinearScanThreshold;
    const std::size_t left_leaf = leaf_for_value(left);
    const std::size_t right_leaf = leaf_for_value(right - 1);
    if (left_leaf == right_leaf) {
      const NodeMetadata& leaf = node_at(0, left_leaf);
      if (left == leaf.value_begin && right == leaf.value_end) {
        return leaf.subtree_min_position;
      }
      return leaf_range_min(leaf, selector_at(0, left_leaf), left, right,
                            allow_leaf_linear_scan);
    }

    const auto [level, node_index] = covering_node(left_leaf, right_leaf);
    return query_node(level, node_index, left, right, allow_leaf_linear_scan);
  }

 private:
  static constexpr std::size_t kMaxLocalEntries =
      LeafSize > Fanout ? LeafSize : Fanout;
  static constexpr std::size_t kMaxLocalBits = 2 * kMaxLocalEntries;
  static constexpr std::size_t kMaxLocalWords = (kMaxLocalBits + 63) / 64;
  static constexpr std::size_t kMaxDepthBlocks = (kMaxLocalBits + 1 + 63) / 64;
  static constexpr std::size_t kInternalScalarSlotThreshold = 4;
  static constexpr std::size_t kLeafLinearScanThreshold = 64;
  static constexpr std::size_t kLeafAvx2ScanThreshold = 16;
  static constexpr bool kLeafSizePowerOfTwo = (LeafSize & (LeafSize - 1)) == 0;
  static constexpr bool kFanoutPowerOfTwo = (Fanout & (Fanout - 1)) == 0;
  static constexpr bool kPowerOfTwoLayout =
      kLeafSizePowerOfTwo && kFanoutPowerOfTwo;
  static constexpr bool kInvalidIndexEqualsNpos =
      static_cast<std::size_t>(invalid_index) == npos;

  static_assert(
      kMaxLocalBits <=
      static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()));
  static_assert(
      kMaxLocalEntries <=
      static_cast<std::size_t>(std::numeric_limits<std::int16_t>::max()));

  struct DepthBlockSummary {
    std::int16_t base_depth = 0;
    std::int16_t min_depth = 0;
    std::uint16_t min_offset = 0;
    std::uint16_t size = 0;
  };

  class LocalSelector {
   public:
    LocalSelector() = default;

    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kMaxLocalEntries) {
        throw std::length_error("NodeEulerBTreeRmq local selector too large");
      }

      entry_count_ = static_cast<std::uint16_t>(entry_count);
      bit_count_ = static_cast<std::uint16_t>(2 * entry_count);
      word_count_ = static_cast<std::uint16_t>((bit_count_ + 63) / 64);
      depth_block_count_ =
          static_cast<std::uint16_t>((bit_count_ + 1 + 63) / 64);
      bp_bits_.fill(0);
      close_positions_.fill(0);
      depth_blocks_.fill({});
      zero_rank_prefix_.fill(0);

      if (entry_count == 0) {
        return;
      }

      std::array<std::uint16_t, kMaxLocalEntries> stack{};
      std::size_t stack_size = 0;
      std::size_t write_position = bit_count_;

      for (std::size_t i = entry_count; i-- > 0;) {
        while (stack_size != 0 && !entry_less(stack[stack_size - 1], i)) {
          --stack_size;
          prepend_bp_bit(write_position, true);
        }
        stack[stack_size++] = static_cast<std::uint16_t>(i);
        close_positions_[i] = prepend_bp_bit(write_position, false);
      }

      while (write_position != 0) {
        prepend_bp_bit(write_position, true);
      }
      build_depth_blocks();
      build_zero_rank_prefix();
    }

    std::size_t arg_min(std::size_t slot_left,
                        std::size_t slot_right,
                        bool force_scalar = false) const {
      if (slot_left >= slot_right || slot_right > entry_count_) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t first_close = close_positions_[slot_left];
      const std::size_t last_close = close_positions_[slot_right - 1];
      if (first_close > last_close) {
        return npos;
      }

      const std::size_t shifted_min =
          depth_arg_min(first_close + 1, last_close + 2, force_scalar);
      if (shifted_min == npos || shifted_min == 0) {
        return npos;
      }

      const std::size_t zero_rank = rank0(shifted_min);
      if (zero_rank == 0) {
        return npos;
      }
      const std::size_t entry = zero_rank - 1;
      return entry < entry_count_ ? entry : npos;
    }

    std::size_t entry_count() const { return entry_count_; }

    std::size_t depth_count() const { return bit_count_ + 1; }

    std::size_t close_position(std::size_t slot) const {
      return close_positions_[slot];
    }

    std::int16_t depth_at_position(std::size_t position) const {
      return depth_at(position);
    }

    std::size_t rank0_at(std::size_t position) const { return rank0(position); }

   private:
    std::size_t prepend_bp_bit(std::size_t& write_position, bool bit) {
      --write_position;
      if (bit) {
        bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                         << (write_position & 63);
      }
      return write_position;
    }

    bool bit(std::size_t position) const {
      return ((bp_bits_[position >> 6] >> (position & 63)) & 1u) != 0;
    }

    static std::uint64_t low_bits_mask(std::size_t count) {
      if (count == 0) {
        return 0;
      }
      if (count >= 64) {
        return ~std::uint64_t{0};
      }
      return (std::uint64_t{1} << count) - 1;
    }

    std::int16_t prefix_excess_in_word(std::size_t word,
                                       std::size_t count) const {
      if (count == 0 || word >= word_count_) {
        return 0;
      }
      const std::uint64_t bits = bp_bits_[word] & low_bits_mask(count);
      const std::size_t ones = std::popcount(bits);
      return static_cast<std::int16_t>(2 * static_cast<std::ptrdiff_t>(ones) -
                                       static_cast<std::ptrdiff_t>(count));
    }

    std::int16_t depth_at(std::size_t position) const {
      const std::size_t block = position >> 6;
      const std::size_t offset = position & 63;
      return static_cast<std::int16_t>(depth_blocks_[block].base_depth +
                                       prefix_excess_in_word(block, offset));
    }

    void build_depth_blocks() {
      const std::size_t depth_count = bit_count_ + 1;
      std::int16_t current_depth = 0;

      for (std::size_t block = 0; block < depth_block_count_; ++block) {
        const std::size_t block_begin = block * 64;
        const std::size_t block_size =
            std::min<std::size_t>(64, depth_count - block_begin);
        DepthBlockSummary summary;
        summary.base_depth = current_depth;
        summary.min_depth = current_depth;
        summary.size = static_cast<std::uint16_t>(block_size);

        for (std::size_t offset = 1; offset < block_size; ++offset) {
          current_depth = static_cast<std::int16_t>(
              current_depth + (bit(block_begin + offset - 1) ? 1 : -1));
          if (current_depth < summary.min_depth) {
            summary.min_depth = current_depth;
            summary.min_offset = static_cast<std::uint16_t>(offset);
          }
        }

        depth_blocks_[block] = summary;
        if (block_begin + block_size < depth_count) {
          current_depth = static_cast<std::int16_t>(
              current_depth + (bit(block_begin + block_size - 1) ? 1 : -1));
        }
      }
    }

    void build_zero_rank_prefix() {
      zero_rank_prefix_[0] = 0;
      for (std::size_t word = 0; word < word_count_; ++word) {
        const std::size_t word_begin = word * 64;
        const std::size_t word_bits =
            std::min<std::size_t>(64, bit_count_ - word_begin);
        const std::uint64_t logical_bits =
            bp_bits_[word] & low_bits_mask(word_bits);
        const std::size_t zeros = word_bits - std::popcount(logical_bits);
        zero_rank_prefix_[word + 1] =
            static_cast<std::uint16_t>(zero_rank_prefix_[word] + zeros);
      }
    }

    std::size_t depth_arg_min_scalar(std::size_t left,
                                     std::size_t right) const {
      const std::size_t depth_count = bit_count_ + 1;
      if (left >= right || right > depth_count) {
        return npos;
      }

      std::size_t position = left;
      std::int16_t best_depth = depth_at(position);
      std::size_t best_position = position;
      std::int16_t current_depth = best_depth;
      ++position;

      while (position < right && (position & 63) != 0) {
        current_depth = static_cast<std::int16_t>(current_depth +
                                                  (bit(position - 1) ? 1 : -1));
        if (current_depth < best_depth) {
          best_depth = current_depth;
          best_position = position;
        }
        ++position;
      }

      while (position + 64 <= right) {
        const DepthBlockSummary& summary = depth_blocks_[position >> 6];
        if (summary.min_depth < best_depth) {
          best_depth = summary.min_depth;
          best_position = position + summary.min_offset;
        }
        position += 64;
      }

      if (position < right) {
        current_depth = depth_at(position);
        if (current_depth < best_depth) {
          best_depth = current_depth;
          best_position = position;
        }
        ++position;
      }
      while (position < right) {
        current_depth = static_cast<std::int16_t>(current_depth +
                                                  (bit(position - 1) ? 1 : -1));
        if (current_depth < best_depth) {
          best_depth = current_depth;
          best_position = position;
        }
        ++position;
      }

      return best_position;
    }

    std::size_t depth_arg_min(std::size_t left,
                              std::size_t right,
                              bool force_scalar) const {
      const std::size_t depth_count = bit_count_ + 1;
      if (left >= right || right > depth_count) {
        return npos;
      }
      if (force_scalar) {
        return depth_arg_min_scalar(left, right);
      }

      std::size_t position = left;
      std::int16_t best_depth = depth_at(position);
      std::size_t best_position = position;

      while (position < right) {
        const std::size_t chunk_begin = (position / 128) * 128;
        const std::size_t local_left = position - chunk_begin;
        const std::size_t local_right =
            std::min<std::size_t>(right - 1, chunk_begin + 128) - chunk_begin;

        if (chunk_begin >= bit_count_) {
          const std::int16_t candidate_depth = depth_at(chunk_begin);
          if (candidate_depth < best_depth) {
            best_depth = candidate_depth;
            best_position = chunk_begin;
          }
        } else {
          const std::size_t word = chunk_begin >> 6;
          const ExcessResult candidate =
              excess_min_128(bp_bits_.data() + word, local_left, local_right);
          const std::int16_t candidate_depth = static_cast<std::int16_t>(
              depth_at(chunk_begin) + candidate.min_excess);
          const std::size_t candidate_position = chunk_begin + candidate.offset;
          if (candidate_depth < best_depth) {
            best_depth = candidate_depth;
            best_position = candidate_position;
          }
        }

        position = chunk_begin + local_right + 1;
      }

      return best_position;
    }

    std::size_t rank0(std::size_t position) const {
      position = std::min<std::size_t>(position, bit_count_);
      const std::size_t full_words = position >> 6;
      std::size_t zeros = zero_rank_prefix_[full_words];
      const std::size_t tail_bits = position & 63;
      if (tail_bits != 0) {
        const std::uint64_t tail =
            bp_bits_[full_words] & low_bits_mask(tail_bits);
        zeros += tail_bits - std::popcount(tail);
      }
      return zeros;
    }

    std::array<std::uint64_t, kMaxLocalWords> bp_bits_{};
    std::array<std::uint16_t, kMaxLocalEntries> close_positions_{};
    std::array<DepthBlockSummary, kMaxDepthBlocks> depth_blocks_{};
    std::array<std::uint16_t, kMaxLocalWords + 1> zero_rank_prefix_{};
    std::uint16_t entry_count_ = 0;
    std::uint16_t bit_count_ = 0;
    std::uint16_t word_count_ = 0;
    std::uint16_t depth_block_count_ = 0;
  };

  class TopDepthSparseSelector {
   public:
    TopDepthSparseSelector() = default;

    void build(const LocalSelector& selector) {
      const std::size_t depth_count = selector.depth_count();
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

    std::size_t arg_min(const LocalSelector& selector,
                        std::size_t slot_left,
                        std::size_t slot_right) const {
      if (slot_left >= slot_right || slot_right > selector.entry_count()) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t first_close = selector.close_position(slot_left);
      const std::size_t last_close = selector.close_position(slot_right - 1);
      if (first_close > last_close) {
        return npos;
      }

      const std::size_t shifted_min =
          depth_arg_min(first_close + 1, last_close + 2);
      if (shifted_min == npos || shifted_min == 0) {
        return npos;
      }

      const std::size_t zero_rank = selector.rank0_at(shifted_min);
      if (zero_rank == 0) {
        return npos;
      }
      const std::size_t entry = zero_rank - 1;
      return entry < selector.entry_count() ? entry : npos;
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

  static constexpr std::size_t log2_exact(std::size_t value) {
    std::size_t shift = 0;
    while (value > 1) {
      value >>= 1;
      ++shift;
    }
    return shift;
  }

  static constexpr std::size_t kLeafShift =
      kLeafSizePowerOfTwo ? log2_exact(LeafSize) : 0;
  static constexpr std::size_t kFanoutShift =
      kFanoutPowerOfTwo ? log2_exact(Fanout) : 0;

  struct alignas(64) NodeMetadata {
    std::size_t value_begin = 0;
    std::size_t value_end = 0;
    std::size_t first_entry = 0;
    Index subtree_min_position = invalid_index;
    std::uint16_t entry_count = 0;
    std::array<std::byte,
               64 - 3 * sizeof(std::size_t) - sizeof(Index) -
                   sizeof(std::uint16_t)>
        padding{};
  };

  static_assert(sizeof(NodeMetadata) == 64);

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
    nodes_.clear();
    selectors_.clear();
    top_depth_selectors_.clear();
    top_level_offsets_.clear();
    level_offsets_.clear();
    level_sizes_.clear();
    level_value_spans_.clear();
    top_level_begin_ = std::numeric_limits<std::size_t>::max();
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("NodeEulerBTreeRmq index type is too small");
    }

    initialize_layout((values_.size() + LeafSize - 1) / LeafSize);
    const std::size_t leaf_count = level_sizes_[0];
    for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
      build_leaf(leaf);
    }

    for (std::size_t child_level = 0; child_level + 1 < level_count();
         ++child_level) {
      const std::size_t parent_count = level_sizes_[child_level + 1];
      for (std::size_t parent = 0; parent < parent_count; ++parent) {
        build_internal_node(child_level, parent);
      }
    }
  }

  void initialize_layout(std::size_t leaf_count) {
    level_sizes_.push_back(leaf_count);
    level_value_spans_.push_back(LeafSize);
    while (level_sizes_.back() > 1) {
      level_sizes_.push_back((level_sizes_.back() + Fanout - 1) / Fanout);
      level_value_spans_.push_back(saturating_product(
          level_value_spans_.back(), static_cast<std::size_t>(Fanout)));
    }

    level_offsets_.assign(level_sizes_.size(), 0);
    std::size_t total_nodes = 0;
    for (std::size_t level = level_sizes_.size(); level-- > 0;) {
      level_offsets_[level] = total_nodes;
      total_nodes += level_sizes_[level];
    }

    nodes_.resize(total_nodes);
    selectors_.resize(total_nodes);
    initialize_top_depth_selector_layout();
  }

  void initialize_top_depth_selector_layout() {
    top_level_begin_ = std::numeric_limits<std::size_t>::max();
    top_level_offsets_.assign(level_sizes_.size(), 0);
    top_depth_selectors_.clear();
    if (level_sizes_.size() <= 1) {
      return;
    }

    const std::size_t root_level = level_sizes_.size() - 1;
    top_level_begin_ = root_level > 1 ? root_level - 1 : 1;

    std::size_t top_node_count = 0;
    for (std::size_t level = top_level_begin_; level < level_sizes_.size();
         ++level) {
      top_level_offsets_[level] = top_node_count;
      top_node_count += level_sizes_[level];
    }
    top_depth_selectors_.resize(top_node_count);
  }

  void build_leaf(std::size_t leaf) {
    NodeMetadata& node = mutable_node_at(0, leaf);
    LocalSelector& selector = mutable_selector_at(0, leaf);
    node.first_entry = leaf * LeafSize;
    node.value_begin = node.first_entry;
    node.entry_count = static_cast<std::uint16_t>(
        std::min<std::size_t>(LeafSize, values_.size() - node.first_entry));
    node.value_end = node.value_begin + node.entry_count;
    selector.build(node.entry_count, [&](std::size_t left, std::size_t right) {
      return compare_(values_[node.first_entry + left],
                      values_[node.first_entry + right]);
    });

    const std::size_t slot = selector.arg_min(0, node.entry_count);
    node.subtree_min_position = static_cast<Index>(node.first_entry + slot);
  }

  void build_internal_node(std::size_t child_level, std::size_t parent) {
    const std::size_t parent_level = child_level + 1;
    NodeMetadata& node = mutable_node_at(parent_level, parent);
    LocalSelector& selector = mutable_selector_at(parent_level, parent);
    node.first_entry = parent * Fanout;
    node.entry_count = static_cast<std::uint16_t>(std::min<std::size_t>(
        Fanout, level_sizes_[child_level] - node.first_entry));
    node.value_begin = node_at(child_level, node.first_entry).value_begin;
    node.value_end =
        node_at(child_level, node.first_entry + node.entry_count - 1).value_end;
    selector.build(node.entry_count, [&](std::size_t left, std::size_t right) {
      return strictly_better_position(
          node_at(child_level, node.first_entry + left).subtree_min_position,
          node_at(child_level, node.first_entry + right).subtree_min_position);
    });

    const std::size_t slot = selector.arg_min(0, node.entry_count);
    node.subtree_min_position =
        node_at(child_level, node.first_entry + slot).subtree_min_position;
    if (is_top_internal_level(parent_level)) {
      mutable_top_depth_selector_at(parent_level, parent).build(selector);
    }
  }

  static std::size_t saturating_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
      return std::numeric_limits<std::size_t>::max();
    }
    return left * right;
  }

  std::size_t leaf_range_min(const NodeMetadata& node,
                             const LocalSelector& selector,
                             std::size_t left,
                             std::size_t right,
                             bool allow_linear_scan) const {
    if (allow_linear_scan && right - left <= kLeafLinearScanThreshold) {
      return linear_range_min(left, right);
    }
    const std::size_t slot_left = left - node.value_begin;
    const std::size_t slot_right = right - node.value_begin;
    const std::size_t slot = selector.arg_min(slot_left, slot_right);
    if (slot == npos) {
      return npos;
    }
    return node.value_begin + slot;
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
      if constexpr (kFanoutPowerOfTwo) {
        left_node >>= kFanoutShift;
        right_node >>= kFanoutShift;
      } else {
        left_node /= Fanout;
        right_node /= Fanout;
      }
    }
    return {level, left_node};
  }

  std::size_t leaf_for_value(std::size_t position) const {
    if constexpr (kLeafSizePowerOfTwo) {
      return position >> kLeafShift;
    } else {
      return position / LeafSize;
    }
  }

  std::size_t child_for_value(std::size_t child_level,
                              std::size_t position) const {
    if constexpr (kPowerOfTwoLayout) {
      return position >> (kLeafShift + child_level * kFanoutShift);
    } else {
      return position / level_value_spans_[child_level];
    }
  }

  bool contains_position(std::size_t left,
                         std::size_t right,
                         std::size_t position) const {
    return !missing_position(position) && left <= position && position < right;
  }

  std::size_t query_child_slots(std::size_t level,
                                std::size_t node_index,
                                std::size_t slot_left,
                                std::size_t slot_right,
                                std::size_t left,
                                std::size_t right,
                                bool allow_leaf_linear_scan) const {
    if (slot_left >= slot_right) {
      return npos;
    }

    const NodeMetadata& node = node_at(level, node_index);
    const std::size_t child_level = level - 1;
    const std::size_t slot_count = slot_right - slot_left;
    const std::size_t slot =
        slot_count == 1
            ? slot_left
            : selector_arg_min(level, node_index, slot_left, slot_right,
                               slot_count <= kInternalScalarSlotThreshold);
    if (slot == npos) {
      return npos;
    }

    const std::size_t child_index = node.first_entry + slot;
    const NodeMetadata& child = node_at(child_level, child_index);
    if ((left <= child.value_begin && child.value_end <= right) ||
        contains_position(left, right, child.subtree_min_position)) {
      return child.subtree_min_position;
    }

    std::size_t answer =
        query_node(child_level, child_index, std::max(left, child.value_begin),
                   std::min(right, child.value_end), allow_leaf_linear_scan);
    answer = better_position(
        answer, query_child_slots(level, node_index, slot_left, slot, left,
                                  right, allow_leaf_linear_scan));
    answer = better_position(
        answer, query_child_slots(level, node_index, slot + 1, slot_right, left,
                                  right, allow_leaf_linear_scan));
    return answer;
  }

  std::size_t query_node(std::size_t level,
                         std::size_t node_index,
                         std::size_t left,
                         std::size_t right,
                         bool allow_leaf_linear_scan) const {
    if (left >= right) {
      return npos;
    }
    const NodeMetadata& node = node_at(level, node_index);
    if (left <= node.value_begin && node.value_end <= right) {
      return node.subtree_min_position;
    }
    if (level == 0) {
      return leaf_range_min(node, selector_at(level, node_index), left, right,
                            allow_leaf_linear_scan);
    }

    const std::size_t child_level = level - 1;
    const std::size_t left_child = child_for_value(child_level, left);
    const std::size_t right_child = child_for_value(child_level, right - 1);
    const std::size_t left_slot = left_child - node.first_entry;
    const std::size_t right_slot = right_child - node.first_entry + 1;
    return query_child_slots(level, node_index, left_slot, right_slot, left,
                             right, allow_leaf_linear_scan);
  }

  std::size_t level_count() const { return level_offsets_.size(); }

  std::size_t flat_index(std::size_t level, std::size_t node_index) const {
    return level_offsets_[level] + node_index;
  }

  bool is_top_internal_level(std::size_t level) const {
    return level > 0 && level >= top_level_begin_ && level < level_count();
  }

  std::size_t top_flat_index(std::size_t level, std::size_t node_index) const {
    return top_level_offsets_[level] + node_index;
  }

  std::size_t selector_arg_min(std::size_t level,
                               std::size_t node_index,
                               std::size_t slot_left,
                               std::size_t slot_right,
                               bool force_scalar) const {
    const LocalSelector& selector = selector_at(level, node_index);
    if (is_top_internal_level(level)) {
      return top_depth_selector_at(level, node_index)
          .arg_min(selector, slot_left, slot_right);
    }
    return selector.arg_min(slot_left, slot_right, force_scalar);
  }

  const NodeMetadata& node_at(std::size_t level, std::size_t node_index) const {
    return nodes_[flat_index(level, node_index)];
  }

  NodeMetadata& mutable_node_at(std::size_t level, std::size_t node_index) {
    return nodes_[flat_index(level, node_index)];
  }

  const LocalSelector& selector_at(std::size_t level,
                                   std::size_t node_index) const {
    return selectors_[flat_index(level, node_index)];
  }

  LocalSelector& mutable_selector_at(std::size_t level,
                                     std::size_t node_index) {
    return selectors_[flat_index(level, node_index)];
  }

  const TopDepthSparseSelector& top_depth_selector_at(
      std::size_t level,
      std::size_t node_index) const {
    return top_depth_selectors_[top_flat_index(level, node_index)];
  }

  TopDepthSparseSelector& mutable_top_depth_selector_at(
      std::size_t level,
      std::size_t node_index) {
    return top_depth_selectors_[top_flat_index(level, node_index)];
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<NodeMetadata> nodes_;
  std::vector<LocalSelector> selectors_;
  std::vector<TopDepthSparseSelector> top_depth_selectors_;
  std::vector<std::size_t> level_offsets_;
  std::vector<std::size_t> top_level_offsets_;
  std::vector<std::size_t> level_sizes_;
  std::vector<std::size_t> level_value_spans_;
  std::size_t top_level_begin_ = std::numeric_limits<std::size_t>::max();
};

}  // namespace pixie::rmq
