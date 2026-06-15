#pragma once

#include <pixie/bits.h>
#include <pixie/bitvector.h>
#include <pixie/rmm_base.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::experimental {

/**
 * @brief Cache-aligned btree implementation of the range min-max index.
 * @details `RmMBTree` is a non-owning succinct index over a little-endian bit
 * sequence. It provides the `RmMBase` operations by combining a rank/select
 * helper with a hierarchy of min-max summaries over fixed 512-bit blocks. Each
 * summary stores subtree size, number of ones, total excess, minimum/maximum
 * relative excess, and minimum multiplicity, which lets searches skip whole
 * subtrees when the requested excess cannot occur there.
 *
 * The internal layout has three layers:
 * - the original bit words, referenced by `std::span<const uint64_t>`;
 * - level 0 block summaries, reconstructed from parent nodes and scanned with
 *   128-bit chunk primitives inside each 512-bit block;
 * - higher btree levels stored as cache-line-aligned summary nodes. The first
 *   level above blocks uses compact `int16_t` excess fields because one low
 *   node covers at most `512 * LowFanout` bits. Higher levels use `int64_t`
 *   excess and count fields.
 *
 * Low nodes are one-cache-line-aligned arrays of compact per-child summaries:
 * @code
 * LowNode<LowFanout = 32>
 * +-----------------------+
 * | cumulative excess     |  int16_t  prefix_excess[32], prefix through child
 * +-----------------------+
 * | child relative min    |  int16_t  min_excess[32]
 * +-----------------------+
 * | child relative max    |  int16_t  max_excess[32]
 * +-----------------------+
 * | count of minima       |  uint16_t min_count[32]
 * +-----------------------+
 * @endcode
 *
 * High nodes have the same logical fields, but use wider lanes and a fanout
 * chosen from @p HighCacheLines:
 * @code
 * HighNode<kHighFanout>
 * +-----------------------+
 * | cumulative excess     |  int64_t  prefix_excess[kHighFanout], prefix
 * |                       |  through child
 * +-----------------------+
 * | child relative min    |  int64_t  min_excess[kHighFanout]
 * +-----------------------+
 * | child relative max    |  int64_t  max_excess[kHighFanout]
 * +-----------------------+
 * | count of minima       |  uint64_t min_count[kHighFanout]
 * +-----------------------+
 * @endcode
 *
 * Forward and backward excess searches ascend from the starting block until a
 * sibling summary can contain the target, then descend through matching child
 * ranges. Node scans use SIMD helpers from `bits.h` when available: low nodes
 * compare 16 `int16_t` lanes at a time, and high nodes compare four `int64_t`
 * lanes at a time. Scalar fallbacks are kept for partial chunks and targets
 * that cannot be represented in the low-node lane type.
 *
 * @tparam HighCacheLines Number of 64-byte cache lines assigned to one high
 * summary node.
 * @tparam LowFanout Number of child summaries stored in one low-level node.
 */
template <std::size_t HighCacheLines = 4, std::size_t LowFanout = 32>
class RmMBTree : public RmMBase<RmMBTree<HighCacheLines, LowFanout>> {
 public:
  static_assert(HighCacheLines > 0);
  static_assert(LowFanout > 0);

  static constexpr std::size_t npos =
      RmMBase<RmMBTree<HighCacheLines, LowFanout>>::npos;
  static constexpr std::size_t kBlockBits = 512;
  static constexpr std::size_t kBlockWords = kBlockBits / 64;
  static constexpr std::size_t kCacheLineBytes = 64;
  static constexpr std::size_t kLowFanout = LowFanout;
  static constexpr std::size_t kHighFanout =
      std::max<std::size_t>(2, (512 * HighCacheLines) / (4 * 64));
  static constexpr std::size_t kMaxFanout = std::max(kLowFanout, kHighFanout);
  static_assert(kMaxFanout <= 64);
  static_assert(
      kBlockBits * kLowFanout <=
      static_cast<std::size_t>(std::numeric_limits<std::int16_t>::max()));

  RmMBTree() = default;
  RmMBTree(const RmMBTree&) = default;
  RmMBTree(RmMBTree&&) noexcept = default;
  RmMBTree& operator=(const RmMBTree&) = default;
  RmMBTree& operator=(RmMBTree&&) noexcept = default;

  /**
   * @brief Minimum position and value returned by one range-min traversal.
   */
  struct RangeMinQueryResult {
    std::size_t position = npos;
    int value = 0;
  };

  /**
   * @brief Construct an RmM btree over an external bit-vector span.
   * @details The tree stores a non-owning view of @p words and builds its
   * rank/select helper plus min-max summaries over the first @p bit_count bits.
   * The optional block-size argument is accepted for API compatibility; this
   * implementation uses the fixed 512-bit block size.
   * @param words External little-endian 64-bit words that contain the bit
   * sequence.
   * @param bit_count Number of valid bits in @p words.
   */
  explicit RmMBTree(std::span<const std::uint64_t> words,
                    std::size_t bit_count,
                    std::size_t = kBlockBits) {
    build(words, bit_count);
  }

  /**
   * @brief Construct with explicit rank/select support options.
   *
   * @details This keeps the normal RmMBTree API full-featured by default while
   * allowing adapters that only need rank/select0 to skip select1 samples.
   * @param one_count Optional exact number of 1-bits for exact select-sample
   * allocation.
   */
  RmMBTree(std::span<const std::uint64_t> words,
           std::size_t bit_count,
           BitVector::SelectSupport select_support,
           std::optional<std::size_t> one_count = std::nullopt,
           std::size_t = kBlockBits) {
    build(words, bit_count, select_support, one_count);
  }

  std::size_t size_impl() const { return bit_count_; }

  std::size_t rank1_impl(std::size_t end_position) const {
    return rank_index_ ? rank_index_->rank(end_position) : 0;
  }

  std::size_t rank0_impl(std::size_t end_position) const {
    return rank_index_ ? rank_index_->rank0(end_position) : 0;
  }

  std::size_t select1_impl(std::size_t rank) const {
    if (!rank_index_ || rank == 0) {
      return npos;
    }
    const std::size_t position = rank_index_->select(rank);
    return position < bit_count_ ? position : npos;
  }

  std::size_t select0_impl(std::size_t rank) const {
    if (!rank_index_ || rank == 0) {
      return npos;
    }
    const std::size_t position = rank_index_->select0(rank);
    return position < bit_count_ ? position : npos;
  }

  std::size_t rank10_impl(std::size_t end_position) const {
    if (end_position <= 1 || bit_count_ == 0) {
      return 0;
    }
    end_position = std::min(end_position, bit_count_);
    std::size_t count = 0;
    for (std::size_t position = 0; position + 1 < end_position; ++position) {
      count += bit(position) == 1 && bit(position + 1) == 0;
    }
    return count;
  }

  std::size_t select10_impl(std::size_t rank) const {
    if (rank == 0) {
      return npos;
    }
    for (std::size_t position = 0; position + 1 < bit_count_; ++position) {
      if (bit(position) == 1 && bit(position + 1) == 0 && --rank == 0) {
        return position;
      }
    }
    return npos;
  }

  int excess_impl(std::size_t end_position) const {
    end_position = std::min(end_position, bit_count_);
    return static_cast<int>(
        2 * static_cast<std::int64_t>(rank1_impl(end_position)) -
        static_cast<std::int64_t>(end_position));
  }

  std::size_t fwdsearch_impl(std::size_t start_position, int delta) const {
    if (start_position >= bit_count_) {
      return npos;
    }

    const std::size_t block_index = start_position / kBlockBits;
    const std::size_t block_begin = block_index * kBlockBits;
    const std::size_t start_offset = start_position - block_begin;
    const std::size_t block_result =
        find_fwd_in_block(block_index, start_offset, delta);
    if (block_result != npos) {
      return block_result;
    }

    const std::int64_t relative_target =
        block_excess_at_local(block_index, start_offset) + delta;
    const std::int64_t block_start_excess = prefix_excess_impl(block_begin);
    const std::int64_t target = block_start_excess + relative_target;
    std::size_t level = 0;
    std::size_t index = block_index;
    while (has_parent_level(level)) {
      const std::size_t fanout = fanout_to_parent(level);
      const std::size_t parent = index / fanout;
      const std::size_t sibling_end =
          std::min(level_count(level), parent * fanout + fanout);
      const NodeScanResult scan = scan_children_fwd(
          level, parent, index % fanout + 1, sibling_end - parent * fanout,
          target, prefix_excess_impl(node_end_bit(level, index)));
      if (scan.found) {
        const std::size_t result =
            descend_fwd(level, scan.index, target, scan.node_start_excess);
        if (result != npos) {
          return result;
        }
      }
      level += 1;
      index = parent;
    }
    return npos;
  }

  std::size_t bwdsearch_impl(std::size_t start_position, int delta) const {
    if (start_position == 0 || start_position > bit_count_) {
      return npos;
    }

    const std::size_t block_index = (start_position - 1) / kBlockBits;
    const std::size_t block_begin = block_index * kBlockBits;
    const std::size_t end_offset = start_position - block_begin;
    const std::size_t block_result =
        find_bwd_in_block(block_index, end_offset, delta);
    if (block_result != npos) {
      return block_result;
    }

    const std::int64_t relative_target =
        block_excess_at_local(block_index, end_offset) + delta;
    const std::int64_t block_start_excess = prefix_excess_impl(block_begin);
    const std::int64_t target = block_start_excess + relative_target;
    std::size_t level = 0;
    std::size_t index = block_index;
    while (has_parent_level(level)) {
      const std::size_t fanout = fanout_to_parent(level);
      const std::size_t parent = index / fanout;
      const NodeScanResult scan =
          scan_children_bwd(level, parent, 0, index % fanout, target,
                            prefix_excess_impl(node_start_bit(level, index)));
      if (scan.found) {
        if (scan.boundary_only) {
          return node_start_bit(level, scan.index);
        }
        const std::size_t result =
            descend_bwd(level, scan.index, target, scan.node_start_excess);
        if (result != npos) {
          return result;
        }
      }
      level += 1;
      index = parent;
    }
    return npos;
  }

  std::size_t range_min_query_pos_impl(std::size_t range_begin,
                                       std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return npos;
    }
    return range_extreme_query_pos(range_begin, range_end, true);
  }

  int range_min_query_val_impl(std::size_t range_begin,
                               std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return 0;
    }
    return range_extreme_query_val(range_begin, range_end, true);
  }

  /**
   * @brief Return the first minimum position and value in one traversal.
   * @details This is a concrete-type helper for adapters that need to compare
   * a boundary prefix against the interior range minimum. Invalid ranges return
   * `{npos, 0}`.
   */
  RangeMinQueryResult range_min_query_result(std::size_t range_begin,
                                             std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return {};
    }
    const auto result = range_extreme_query(range_begin, range_end, true);
    return {result.position, static_cast<int>(result.value)};
  }

  std::size_t range_max_query_pos_impl(std::size_t range_begin,
                                       std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return npos;
    }
    return range_extreme_query_pos(range_begin, range_end, false);
  }

  int range_max_query_val_impl(std::size_t range_begin,
                               std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return 0;
    }
    return range_extreme_query_val(range_begin, range_end, false);
  }

  std::size_t mincount_impl(std::size_t range_begin,
                            std::size_t range_end) const {
    if (range_begin > range_end || range_end >= bit_count_) {
      return 0;
    }
    return range_min_stats(range_begin, range_end).count;
  }

  std::size_t minselect_impl(std::size_t range_begin,
                             std::size_t range_end,
                             std::size_t rank) const {
    if (range_begin > range_end || range_end >= bit_count_ || rank == 0) {
      return npos;
    }
    const RangeMinStats stats = range_min_stats(range_begin, range_end);
    if (rank > stats.count) {
      return npos;
    }
    return range_min_select(range_begin, range_end, stats.value, rank);
  }

  std::size_t close_impl(std::size_t open_position) const {
    if (open_position >= bit_count_) {
      return npos;
    }
    if (!bit(open_position)) {
      return open_position;
    }
    return fwd_excess_at(open_position, -1);
  }

  std::size_t open_impl(std::size_t close_position) const {
    if (close_position >= bit_count_) {
      return npos;
    }
    if (bit(close_position)) {
      return close_position;
    }
    return bwdsearch_impl(close_position + 1, 0);
  }

  std::size_t enclose_impl(std::size_t position) const {
    if (position >= bit_count_) {
      return npos;
    }
    if (!bit(position)) {
      return open_impl(position);
    }
    return bwdsearch_impl(position + 1, -2);
  }

 private:
  /**
   * @brief Search forward using SDSL-style inclusive-position excess semantics.
   * @details Public `fwdsearch_impl` starts from a prefix boundary. This helper
   * starts at @p position as an included bit position, so the equivalent btree
   * boundary search begins at `position + 1`. It is used by
   * balanced-parentheses operations such as `close`, where SDSL's
   * `fwd_excess(i, -1)` maps to this wrapper.
   * @param position Included bit position used as the search origin.
   * @param delta Desired excess delta from the prefix after @p position.
   * @return Matching bit position, or `npos` when no forward match exists.
   */
  std::size_t fwd_excess_at(std::size_t position, int delta) const {
    if (position >= bit_count_) {
      return npos;
    }
    if (position + 1 >= bit_count_) {
      return npos;
    }
    return fwdsearch_impl(position + 1, delta);
  }

  struct Summary {
    std::uint64_t size_bits = 0;
    std::uint64_t ones = 0;
    std::int64_t block_excess = 0;
    std::int64_t min_excess = 0;
    std::int64_t max_excess = 0;
    std::uint64_t min_count = 0;
  };

  template <class Excess, class Count, std::size_t Fanout>
  struct alignas(kCacheLineBytes) SummaryNode {
    using ExcessType = Excess;
    using CountType = Count;
    static constexpr std::size_t kFanout = Fanout;
    std::array<Excess, Fanout> prefix_excess{};
    std::array<Excess, Fanout> min_excess{};
    std::array<Excess, Fanout> max_excess{};
    std::array<Count, Fanout> min_count{};
  };

  using LowNode = SummaryNode<std::int16_t, std::uint16_t, kLowFanout>;
  using HighNode = SummaryNode<std::int64_t, std::uint64_t, kHighFanout>;
  static_assert(alignof(LowNode) == kCacheLineBytes);
  static_assert(alignof(HighNode) == kCacheLineBytes);
  static_assert(sizeof(LowNode) % kCacheLineBytes == 0);
  static_assert(sizeof(HighNode) % kCacheLineBytes == 0);

  struct ByteAgg {
    std::int8_t block_excess = 0;
    std::int8_t min_excess = 0;
    std::int8_t max_excess = 0;
    std::uint8_t min_count = 0;
    std::uint8_t pos_first_min = 0;
    std::uint8_t pos_first_max = 0;
  };

  static constexpr std::size_t kSearchChunkBits = 128;
  static constexpr std::size_t kSearchChunkWords = kSearchChunkBits / 64;
  static constexpr std::size_t kSearchChunkCount =
      kBlockBits / kSearchChunkBits;

  /**
   * @brief Byte-level summaries used by range scans.
   * @details Builds a process-local lookup table indexed by byte value. Each
   * entry records the byte excess, local min/max excess, minimum multiplicity,
   * and first positions of the min/max values.
   * @return Immutable 256-entry byte summary table.
   */
  static const std::array<ByteAgg, 256>& byte_lut() {
    static const std::array<ByteAgg, 256> table = [] {
      std::array<ByteAgg, 256> result{};
      for (int byte_value = 0; byte_value < 256; ++byte_value) {
        ByteAgg agg;
        int current = 0;
        int minimum = std::numeric_limits<int>::max();
        int maximum = std::numeric_limits<int>::min();
        const auto bit_at = [&](int bit_index) {
          return (byte_value >> bit_index) & 1;
        };
        for (int bit_index = 0; bit_index < 8; ++bit_index) {
          const int value = bit_at(bit_index);
          current += value ? 1 : -1;
          if (current < minimum) {
            minimum = current;
            agg.min_count = 1;
            agg.pos_first_min = static_cast<std::uint8_t>(bit_index);
          } else if (current == minimum) {
            ++agg.min_count;
          }
          if (current > maximum) {
            maximum = current;
            agg.pos_first_max = static_cast<std::uint8_t>(bit_index);
          }
        }
        agg.block_excess = static_cast<std::int8_t>(current);
        agg.min_excess = static_cast<std::int8_t>(minimum);
        agg.max_excess = static_cast<std::int8_t>(maximum);
        result[static_cast<std::size_t>(byte_value)] = agg;
      }
      return result;
    }();
    return table;
  }

  /**
   * @brief Build the rank index and min-max tree levels over the input bits.
   * @details Validates that @p words contains enough storage, records the
   * external span, constructs the internal rank/select helper, summarizes each
   * 512-bit block, and builds summary levels above the blocks.
   * @param words External little-endian words that remain owned by the caller.
   * @param bit_count Number of valid bits in @p words.
   * @throws std::invalid_argument if @p words is shorter than required for
   * @p bit_count.
   */
  void build(std::span<const std::uint64_t> words, std::size_t bit_count) {
    build(words, bit_count, BitVector::SelectSupport::kBoth, std::nullopt);
  }

  /**
   * @brief Build with explicit rank/select support options.
   */
  void build(std::span<const std::uint64_t> words,
             std::size_t bit_count,
             BitVector::SelectSupport select_support,
             std::optional<std::size_t> one_count = std::nullopt) {
    const std::size_t required_words = (bit_count + 63) / 64;
    if (words.size() < required_words) {
      throw std::invalid_argument(
          "RmMBTree input span is shorter than bit_count");
    }

    bits_ = words;
    bit_count_ = bit_count;
    rank_index_.emplace(words, bit_count, select_support, one_count);
    block_count_ = (bit_count_ + kBlockBits - 1) / kBlockBits;
    std::vector<Summary> block_summaries(block_count_);

    for (std::size_t block_index = 0; block_index < block_count_;
         ++block_index) {
      const std::size_t block_begin = block_index * kBlockBits;
      const std::size_t block_end =
          std::min(bit_count_, block_begin + kBlockBits);
      block_summaries[block_index] =
          summarize_bits(block_begin, block_end - block_begin);
    }

    build_levels(block_summaries);
  }

  /**
   * @brief Build low and high summary levels from 512-bit block summaries.
   * @details Level 0 is the block-summary level. The first parent level uses
   * compact 16-bit low nodes, and all higher levels use 64-bit high nodes until
   * a single top summary remains.
   * @param block_summaries Relative summaries for each 512-bit block.
   */
  void build_levels(const std::vector<Summary>& block_summaries) {
    level_counts_.clear();
    low_levels_.clear();
    high_levels_.clear();
    top_summary_ = Summary{};
    level_counts_.push_back(block_summaries.size());
    if (block_summaries.empty()) {
      return;
    }

    std::vector<Summary> current = block_summaries;
    current = build_parent_level<LowNode>(current, low_levels_.emplace_back());
    level_counts_.push_back(current.size());
    current =
        build_parent_level<HighNode>(current, high_levels_.emplace_back());
    level_counts_.push_back(current.size());
    while (current.size() > 1) {
      current =
          build_parent_level<HighNode>(current, high_levels_.emplace_back());
      level_counts_.push_back(current.size());
    }
    top_summary_ = current.front();
  }

  /**
   * @brief Build one parent level and return summaries for the produced nodes.
   * @details Groups @p in by the fanout of @p Node, stores each child summary
   * in the parent node, and returns relative summaries for the newly built
   * parent level.
   * @tparam Node LowNode or HighNode.
   * @param in Child summaries from the previous level.
   * @param nodes Destination storage for the constructed parent nodes.
   * @return Relative summaries for @p nodes.
   */
  template <class Node>
  static std::vector<Summary> build_parent_level(const std::vector<Summary>& in,
                                                 std::vector<Node>& nodes) {
    constexpr std::size_t fanout = Node::kFanout;
    std::vector<Summary> out((in.size() + fanout - 1) / fanout);
    nodes.resize(out.size());
    for (std::size_t parent = 0; parent < out.size(); ++parent) {
      const std::size_t begin = parent * fanout;
      const std::size_t end = std::min(in.size(), begin + fanout);
      Summary combined;
      for (std::size_t i = begin; i < end; ++i) {
        store_child_summary(nodes[parent], i - begin, combined.block_excess,
                            in[i]);
        combined = append(combined, in[i]);
      }
      out[parent] = combined;
    }
    return out;
  }

  /**
   * @brief Store one child summary in a node using inclusive prefix excess.
   * @details Stores cumulative excess through @p slot in
   * `prefix_excess[slot]`, along with the child min/max excess and minimum
   * count. The child total excess is later reconstructed from adjacent
   * cumulative prefixes.
   * @tparam Node LowNode or HighNode.
   * @param node Parent node being filled.
   * @param slot Child slot inside @p node.
   * @param prefix_excess Excess before this child within @p node.
   * @param summary Relative summary of the child subtree.
   */
  template <class Node>
  static void store_child_summary(Node& node,
                                  std::size_t slot,
                                  std::int64_t prefix_excess,
                                  const Summary& summary) {
    node.prefix_excess[slot] =
        static_cast<typename decltype(node.prefix_excess)::value_type>(
            prefix_excess + summary.block_excess);
    node.min_excess[slot] =
        static_cast<typename decltype(node.min_excess)::value_type>(
            summary.min_excess);
    node.max_excess[slot] =
        static_cast<typename decltype(node.max_excess)::value_type>(
            summary.max_excess);
    node.min_count[slot] =
        static_cast<typename decltype(node.min_count)::value_type>(
            summary.min_count);
  }

  /**
   * @brief Summarize a contiguous bit range relative to its beginning.
   * @details Scans the range bit by bit and computes total excess, min/max
   * relative excess, number of positions attaining the minimum, and number of 1
   * bits.
   * @param begin First bit position.
   * @param length Number of bits to summarize.
   * @return Relative summary of `[begin, begin + length)`.
   */
  Summary summarize_bits(std::size_t begin, std::size_t length) const {
    if (length == kBlockBits && (begin % kBlockBits) == 0) {
      const std::size_t block_index = begin / kBlockBits;
      if (full_block_has_words(block_index)) {
        return summarize_full_block(block_index);
      }
    }

    Summary summary;
    summary.size_bits = length;
    if (length == 0) {
      return summary;
    }
    int current = 0;
    int minimum = std::numeric_limits<int>::max();
    int maximum = std::numeric_limits<int>::min();
    for (std::size_t offset = 0; offset < length; ++offset) {
      const std::uint8_t value = bit(begin + offset);
      summary.ones += value;
      current += value ? 1 : -1;
      if (current < minimum) {
        minimum = current;
        summary.min_count = 1;
      } else if (current == minimum) {
        ++summary.min_count;
      }
      if (current > maximum) {
        maximum = current;
      }
    }
    summary.block_excess = current;
    summary.min_excess = minimum;
    summary.max_excess = maximum;
    return summary;
  }

  /**
   * @brief Summarize one complete 512-bit block with 128-bit excess kernels.
   * @details Full Cartesian-BP construction blocks are aligned and contain all
   * eight backing words. This fast path reuses the existing SIMD/minimum
   * primitives from `bits.h` instead of reading every bit individually.
   */
  Summary summarize_full_block(std::size_t block_index) const {
    const std::uint64_t* block = bits_.data() + block_index * kBlockWords;
    Summary summary;
    for (std::size_t chunk = 0; chunk < kSearchChunkCount; ++chunk) {
      summary = append(summary,
                       summarize_128_chunk(block + chunk * kSearchChunkWords));
    }
    return summary;
  }

  /**
   * @brief Summarize one full 128-bit chunk relative to its first bit.
   * @details Minima are found directly; maxima are minima of the bit-inverted
   * chunk with the sign flipped. Minimum multiplicity is computed by the
   * existing SIMD target-position kernel.
   */
  static Summary summarize_128_chunk(const std::uint64_t* chunk) {
    Summary summary;
    summary.size_bits = kSearchChunkBits;
    summary.ones = std::popcount(chunk[0]) + std::popcount(chunk[1]);
    summary.block_excess = 2 * static_cast<std::int64_t>(summary.ones) -
                           static_cast<std::int64_t>(kSearchChunkBits);

    const ExcessResult minimum = excess_min_128(chunk, 1, kSearchChunkBits);
    summary.min_excess = minimum.min_excess;

    const std::array<std::uint64_t, kSearchChunkWords> inverted = {~chunk[0],
                                                                   ~chunk[1]};
    const ExcessResult inverted_min =
        excess_min_128(inverted.data(), 1, kSearchChunkBits);
    summary.max_excess = -static_cast<std::int64_t>(inverted_min.min_excess);

    std::uint64_t minimum_positions[kSearchChunkWords];
    excess_positions_128(chunk, static_cast<int>(summary.min_excess),
                         minimum_positions);
    summary.min_count = std::popcount(minimum_positions[0]) +
                        std::popcount(minimum_positions[1]);
    return summary;
  }

  /**
   * @brief Concatenate two relative summaries.
   * @details Produces the summary for `left || right`, translating right-side
   * extrema by the total excess of @p left and merging minimum counts when both
   * sides attain the combined minimum.
   * @param left Summary of the first range.
   * @param right Summary of the following range.
   * @return Relative summary of the concatenation.
   */
  static Summary append(Summary left, const Summary& right) {
    if (left.size_bits == 0) {
      return right;
    }
    if (right.size_bits == 0) {
      return left;
    }
    Summary result;
    result.size_bits = left.size_bits + right.size_bits;
    result.ones = left.ones + right.ones;
    result.block_excess = left.block_excess + right.block_excess;
    result.min_excess =
        std::min(left.min_excess, left.block_excess + right.min_excess);
    result.max_excess =
        std::max(left.max_excess, left.block_excess + right.max_excess);
    result.min_count = 0;
    if (left.min_excess == result.min_excess) {
      result.min_count += left.min_count;
    }
    if (left.block_excess + right.min_excess == result.min_excess) {
      result.min_count += right.min_count;
    }
    return result;
  }

  /**
   * @brief Number of valid bits in a 512-bit block.
   * @details All blocks except the last have `kBlockBits` bits. The last block
   * may be partial when `size()` is not a multiple of 512.
   * @param block_index Zero-based block index.
   * @return Number of valid bits in the block.
   */
  std::size_t block_size(std::size_t block_index) const {
    const std::size_t begin = block_index * kBlockBits;
    return std::min(bit_count_ - begin, kBlockBits);
  }

  /**
   * @brief Whether a block is full and backed by all eight storage words.
   * @details Fast SIMD/chunk search paths require a full 512-bit block and all
   * eight words to be present in the external span.
   * @param block_index Zero-based block index.
   * @return True when the block can be processed as exactly eight words.
   */
  bool full_block_has_words(std::size_t block_index) const {
    return block_size(block_index) == kBlockBits &&
           (block_index + 1) * kBlockWords <= bits_.size();
  }

  /**
   * @brief Prefix excess at a block-local boundary using direct popcounts.
   * @details Reads only the words in the target block instead of using the rank
   * index. This is used on local-search miss paths before ascending the tree.
   * @param block_index Zero-based block index.
   * @param offset Exclusive block-local prefix boundary, clamped to the block
   * length.
   * @return Prefix excess relative to the beginning of the block.
   */
  std::int64_t block_excess_at_local(std::size_t block_index,
                                     std::size_t offset) const {
    const std::size_t length = block_size(block_index);
    offset = std::min(offset, length);
    if (offset == 0) {
      return 0;
    }

    const std::size_t first_word = block_index * kBlockWords;
    std::size_t remaining = offset;
    std::int64_t ones = 0;
    std::size_t word_offset = 0;
    while (remaining >= 64) {
      ones += std::popcount(bits_[first_word + word_offset]);
      remaining -= 64;
      ++word_offset;
    }
    if (remaining != 0) {
      ones += std::popcount(bits_[first_word + word_offset] &
                            first_bits_mask(remaining));
    }
    return 2 * ones - static_cast<std::int64_t>(offset);
  }

  /**
   * @brief Total excess of a 128-bit search chunk.
   * @details Computes `2 * popcount(chunk) - 128` for the two input words.
   * @param chunk Pointer to two little-endian words.
   * @return Total excess of the 128-bit chunk.
   */
  static int chunk_excess_128(const std::uint64_t* chunk) {
    return 2 * static_cast<int>(std::popcount(chunk[0]) +
                                std::popcount(chunk[1])) -
           static_cast<int>(kSearchChunkBits);
  }

  /**
   * @brief Find a forward excess target inside one block from a local boundary.
   * @details Searches for the first bit position `p >= start_offset` in the
   * block such that the prefix excess at `p + 1` equals the excess at
   * `start_offset` plus @p delta. Full blocks use 128-bit chunk primitives;
   * partial blocks fall back to direct scalar scanning.
   * @param block_index Zero-based block index.
   * @param start_offset Block-local start boundary.
   * @param delta Desired excess delta from `start_offset`.
   * @return Global bit position of the match, or `npos`.
   */
  std::size_t find_fwd_in_block(std::size_t block_index,
                                std::size_t start_offset,
                                std::int64_t delta) const {
    const std::size_t length = block_size(block_index);
    if (start_offset >= length) {
      return npos;
    }

    if (full_block_has_words(block_index)) {
      const std::uint64_t* block = bits_.data() + block_index * kBlockWords;
      const std::size_t first_chunk = start_offset / kSearchChunkBits;
      std::int64_t target = delta;
      for (std::size_t chunk = first_chunk; chunk < kSearchChunkCount;
           ++chunk) {
        const std::uint64_t* chunk_words = block + chunk * kSearchChunkWords;
        const std::size_t local_start =
            chunk == first_chunk ? start_offset - chunk * kSearchChunkBits : 0;
        if (chunk == first_chunk) {
          target += prefix_excess_128(chunk_words, local_start);
        }
        int block_excess = 0;
        const std::size_t offset = forward_search_128(
            chunk_words, static_cast<int>(target), local_start, &block_excess);
        if (offset != kSearchChunkBits) {
          return block_index * kBlockBits + chunk * kSearchChunkBits + offset;
        }
        target -= block_excess;
      }
      return npos;
    }

    std::int64_t current = block_excess_at_local(block_index, start_offset);
    const std::int64_t relative_target = current + delta;
    const std::size_t block_begin = block_index * kBlockBits;
    for (std::size_t offset = start_offset; offset < length; ++offset) {
      current += bit(block_begin + offset) ? 1 : -1;
      if (current == relative_target) {
        return block_index * kBlockBits + offset;
      }
    }
    return npos;
  }

  /**
   * @brief Find a backward excess target inside one block from a local
   * boundary.
   * @details Searches for the greatest prefix boundary `p < end_offset` in the
   * block such that the prefix excess at `p` equals the excess at @p end_offset
   * plus @p delta. The returned value is a global prefix-boundary/bit position.
   * @param block_index Zero-based block index.
   * @param end_offset Exclusive block-local boundary to search before.
   * @param delta Desired excess delta from `end_offset`.
   * @return Global position of the matching boundary, or `npos`.
   */
  std::size_t find_bwd_in_block(std::size_t block_index,
                                std::size_t end_offset,
                                std::int64_t delta) const {
    if (end_offset == 0) {
      return npos;
    }
    const std::size_t max_prefix_length = end_offset - 1;

    if (full_block_has_words(block_index)) {
      const std::uint64_t* block = bits_.data() + block_index * kBlockWords;
      const std::size_t last_chunk = max_prefix_length / kSearchChunkBits;
      std::int64_t target = delta;
      for (std::size_t chunk = last_chunk + 1; chunk > 0;) {
        --chunk;
        const std::uint64_t* chunk_words = block + chunk * kSearchChunkWords;
        const std::size_t local_end =
            chunk == last_chunk ? end_offset - chunk * kSearchChunkBits
                                : kSearchChunkBits;
        if (chunk == last_chunk) {
          target += prefix_excess_128(chunk_words, local_end);
        }
        int block_excess = 0;
        const std::size_t offset = backward_search_128(
            chunk_words, static_cast<int>(target), local_end, &block_excess);
        if (offset != kSearchChunkBits) {
          return block_index * kBlockBits + chunk * kSearchChunkBits + offset;
        }
        if (chunk > 0) {
          target += chunk_excess_128(block + (chunk - 1) * kSearchChunkWords);
        }
      }
      return npos;
    }

    const std::int64_t relative_target =
        block_excess_at_local(block_index, end_offset) + delta;
    std::int64_t current =
        block_excess_at_local(block_index, max_prefix_length);
    const std::size_t block_begin = block_index * kBlockBits;
    for (std::size_t prefix_length = max_prefix_length; prefix_length > 0;
         --prefix_length) {
      if (current == relative_target) {
        return block_index * kBlockBits + prefix_length;
      }
      current -= bit(block_begin + prefix_length - 1) ? 1 : -1;
    }
    return relative_target == 0 ? block_index * kBlockBits : npos;
  }

  /**
   * @brief Descend from a matching summary node to the first forward match.
   * @details Starting from a node whose min/max range contains @p target,
   * repeatedly chooses the leftmost matching child and finally scans the leaf
   * block.
   * @param level Summary level of the starting node.
   * @param index Node index at @p level.
   * @param target Absolute prefix excess target.
   * @param node_start_excess Absolute excess before the starting node.
   * @return Global bit position of the forward match, or `npos`.
   */
  std::size_t descend_fwd(std::size_t level,
                          std::size_t index,
                          std::int64_t target,
                          std::int64_t node_start_excess) const {
    while (level > 0) {
      const std::size_t child_level = level - 1;
      const std::size_t fanout = fanout_to_parent(child_level);
      const std::size_t child_begin = index * fanout;
      const std::size_t child_end =
          std::min(level_count(child_level), child_begin + fanout);
      const NodeScanResult scan =
          scan_children_fwd(child_level, index, 0, child_end - child_begin,
                            target, node_start_excess);
      if (!scan.found) {
        return npos;
      }
      index = scan.index;
      level = child_level;
      node_start_excess = scan.node_start_excess;
    }
    return find_fwd_in_block(index, 0, target - node_start_excess);
  }

  /**
   * @brief Descend from a matching summary node to the last backward match.
   * @details Starting from a node whose min/max range can contain @p target,
   * repeatedly chooses the rightmost matching child and finally scans the leaf
   * block backward. Boundary-only matches return the child start directly.
   * @param level Summary level of the starting node.
   * @param index Node index at @p level.
   * @param target Absolute prefix excess target.
   * @param node_start_excess Absolute excess before the starting node.
   * @return Global position of the backward match, or `npos`.
   */
  std::size_t descend_bwd(std::size_t level,
                          std::size_t index,
                          std::int64_t target,
                          std::int64_t node_start_excess) const {
    while (level > 0) {
      const std::size_t child_level = level - 1;
      const std::size_t fanout = fanout_to_parent(child_level);
      const std::size_t child_begin = index * fanout;
      const std::size_t child_end =
          std::min(level_count(child_level), child_begin + fanout);
      std::int64_t child_start_excess =
          node_start_excess + summary_at(level, index).block_excess;
      const NodeScanResult scan =
          scan_children_bwd(child_level, index, 0, child_end - child_begin,
                            target, child_start_excess);
      if (!scan.found) {
        return npos;
      }
      if (scan.boundary_only) {
        return node_start_bit(child_level, scan.index);
      }
      index = scan.index;
      level = child_level;
      node_start_excess = scan.node_start_excess;
    }
    const std::int64_t block_excess = summary_at(0, index).block_excess;
    return find_bwd_in_block(index, block_size(index),
                             target - node_start_excess - block_excess);
  }

  struct NodeScanResult {
    bool found = false;
    bool boundary_only = false;
    std::size_t index = 0;
    std::int64_t node_start_excess = 0;
  };

  /**
   * @brief Scan child summaries in increasing order for a forward match.
   * @details Dispatches to either low-node or high-node scanning depending on
   * @p child_level. The scan interval is `[begin_slot, end_slot)`.
   * @param child_level Level of the children being scanned.
   * @param parent Parent node index.
   * @param begin_slot First child slot to inspect.
   * @param end_slot One-past-last child slot to inspect.
   * @param target Absolute prefix excess target.
   * @param begin_excess Absolute excess at `begin_slot`.
   * @return Scan result describing the first matching child, if any.
   */
  NodeScanResult scan_children_fwd(std::size_t child_level,
                                   std::size_t parent,
                                   std::size_t begin_slot,
                                   std::size_t end_slot,
                                   std::int64_t target,
                                   std::int64_t begin_excess) const {
    if (begin_slot >= end_slot) {
      return {};
    }
    if (child_level == 0) {
      return scan_node_fwd(low_levels_[0][parent], child_level, parent,
                           begin_slot, end_slot, target, begin_excess);
    }
    return scan_node_fwd(high_levels_[child_level - 1][parent], child_level,
                         parent, begin_slot, end_slot, target, begin_excess);
  }

  /**
   * @brief Scan child summaries in decreasing order for a backward match.
   * @details Dispatches to either low-node or high-node scanning depending on
   * @p child_level. The scan interval is `[begin_slot, end_slot)`, inspected
   * from right to left.
   * @param child_level Level of the children being scanned.
   * @param parent Parent node index.
   * @param begin_slot First child slot in the scan interval.
   * @param end_slot One-past-last child slot in the scan interval.
   * @param target Absolute prefix excess target.
   * @param end_excess Absolute excess at `end_slot`.
   * @return Scan result describing the last matching child, if any.
   */
  NodeScanResult scan_children_bwd(std::size_t child_level,
                                   std::size_t parent,
                                   std::size_t begin_slot,
                                   std::size_t end_slot,
                                   std::int64_t target,
                                   std::int64_t end_excess) const {
    if (begin_slot >= end_slot) {
      return {};
    }
    if (child_level == 0) {
      return scan_node_bwd(low_levels_[0][parent], child_level, parent,
                           begin_slot, end_slot, target, end_excess);
    }
    return scan_node_bwd(high_levels_[child_level - 1][parent], child_level,
                         parent, begin_slot, end_slot, target, end_excess);
  }

  /**
   * @brief Scan one concrete node type in increasing slot order.
   * @details Converts the absolute target into each child's relative target and
   * checks the child min/max range in SIMD-sized groups.
   * @tparam Node LowNode or HighNode.
   * @param node Parent node being scanned.
   * @param child_level Level of @p node's children.
   * @param parent Parent node index.
   * @param begin_slot First child slot to inspect.
   * @param end_slot One-past-last child slot to inspect.
   * @param target Absolute prefix excess target.
   * @param begin_excess Absolute excess at `begin_slot`.
   * @return First matching child, if any.
   */
  template <class Node>
  NodeScanResult scan_node_fwd(const Node& node,
                               std::size_t child_level,
                               std::size_t parent,
                               std::size_t begin_slot,
                               std::size_t end_slot,
                               std::int64_t target,
                               std::int64_t begin_excess) const {
    const std::int64_t node_base_excess =
        begin_excess - prefix_excess_at(node, begin_slot);
    for (std::size_t slot = begin_slot; slot < end_slot;) {
      const std::size_t lane_count =
          std::min(vector_lane_count<Node>(), end_slot - slot);
      const std::uint32_t mask = matching_chunk_mask(
          node, slot, lane_count, target - node_base_excess, false);
      if (mask != 0) {
        const std::size_t lane = std::countr_zero(mask);
        const std::size_t matched_slot = slot + lane;
        return {true, false,
                parent * fanout_to_parent(child_level) + matched_slot,
                child_start_excess(node, node_base_excess, matched_slot)};
      }
      slot += lane_count;
    }
    return {};
  }

  /**
   * @brief Scan one concrete node type in decreasing slot order.
   * @details Converts the absolute target into each child's relative target and
   * checks the child min/max range in SIMD-sized groups from right to left.
   * @tparam Node LowNode or HighNode.
   * @param node Parent node being scanned.
   * @param child_level Level of @p node's children.
   * @param parent Parent node index.
   * @param begin_slot First child slot in the scan interval.
   * @param end_slot One-past-last child slot in the scan interval.
   * @param target Absolute prefix excess target.
   * @param end_excess Absolute excess at `end_slot`.
   * @return Last matching child, if any.
   */
  template <class Node>
  NodeScanResult scan_node_bwd(const Node& node,
                               std::size_t child_level,
                               std::size_t parent,
                               std::size_t begin_slot,
                               std::size_t end_slot,
                               std::int64_t target,
                               std::int64_t end_excess) const {
    const std::int64_t node_base_excess =
        end_excess - prefix_excess_at(node, end_slot);
    for (std::size_t slot_end = end_slot; slot_end > begin_slot;) {
      const std::size_t lane_count =
          std::min(vector_lane_count<Node>(), slot_end - begin_slot);
      const std::size_t slot = slot_end - lane_count;
      const std::uint32_t mask = matching_chunk_mask(
          node, slot, lane_count, target - node_base_excess, true);
      if (mask != 0) {
        const std::size_t lane =
            static_cast<std::size_t>(std::bit_width(mask) - 1);
        const std::size_t matched_slot = slot + lane;
        const std::int64_t relative_target =
            target - child_start_excess(node, node_base_excess, matched_slot);
        const bool interior_match =
            node.min_excess[matched_slot] <= relative_target &&
            relative_target <= node.max_excess[matched_slot];
        return {true, !interior_match,
                parent * fanout_to_parent(child_level) + matched_slot,
                child_start_excess(node, node_base_excess, matched_slot)};
      }
      slot_end = slot;
    }
    return {};
  }

  /**
   * @brief Inclusive prefix excess before child slot.
   * @details Node storage keeps inclusive cumulative excess through each child.
   * This helper converts that representation to the exclusive prefix before a
   * slot.
   * @tparam Node LowNode or HighNode.
   * @param node Node containing cumulative prefixes.
   * @param slot Child slot boundary in `[0, Node::kFanout]`.
   * @return Excess before @p slot within @p node.
   */
  template <class Node>
  static std::int64_t prefix_excess_at(const Node& node, std::size_t slot) {
    if (slot == 0) {
      return 0;
    }
    return prefix_through(node, slot - 1);
  }

  /**
   * @brief Absolute excess at the start of a child slot.
   * @details Adds the node's absolute base excess to the child-local exclusive
   * prefix before @p slot.
   * @tparam Node LowNode or HighNode.
   * @param node Node containing cumulative prefixes.
   * @param node_base_excess Absolute excess before child slot 0.
   * @param slot Child slot.
   * @return Absolute excess before the child at @p slot.
   */
  template <class Node>
  static std::int64_t child_start_excess(const Node& node,
                                         std::int64_t node_base_excess,
                                         std::size_t slot) {
    return node_base_excess + prefix_excess_at(node, slot);
  }

  /**
   * @brief Inclusive prefix excess through child slot.
   * @details Reads the stored cumulative prefix value for a child slot.
   * @tparam Node LowNode or HighNode.
   * @param node Node containing cumulative prefixes.
   * @param slot Child slot to read.
   * @return Excess through @p slot within @p node.
   */
  template <class Node>
  static std::int64_t prefix_through(const Node& node, std::size_t slot) {
    return static_cast<std::int64_t>(node.prefix_excess[slot]);
  }

  /**
   * @brief Total excess of one child reconstructed from inclusive prefixes.
   * @details Computes the difference between the inclusive prefix through
   * @p slot and the exclusive prefix before @p slot.
   * @tparam Node LowNode or HighNode.
   * @param node Node containing cumulative prefixes.
   * @param slot Child slot.
   * @return Total excess of the child subtree at @p slot.
   */
  template <class Node>
  static std::int64_t child_excess(const Node& node, std::size_t slot) {
    return prefix_through(node, slot) - prefix_excess_at(node, slot);
  }

  /**
   * @brief SIMD lane count used for one node scan chunk.
   * @details Low nodes use 16 int16 lanes, while high nodes use 4 int64 lanes.
   * @tparam Node LowNode or HighNode.
   * @return Number of child slots processed per vector chunk.
   */
  template <class Node>
  static constexpr std::size_t vector_lane_count() {
    if constexpr (std::is_same_v<typename Node::ExcessType, std::int16_t>) {
      return 16;
    } else {
      return 4;
    }
  }

  /**
   * @brief Return a bit mask of child lanes whose min/max range can match.
   * @details Uses AVX2 specializations when the lane count exactly matches the
   * node type; otherwise falls back to scalar range checks.
   * @tparam Node LowNode or HighNode.
   * @param node Node whose child ranges are checked.
   * @param slot First child slot represented by lane 0.
   * @param lane_count Number of valid lanes to check.
   * @param target_in_node Target excess relative to the start of @p node.
   * @param include_zero_boundary Whether a relative target of zero is accepted
   * as a boundary-only match for backward search.
   * @return Bit mask with bit `i` set when lane `i` can match.
   */
  template <class Node>
  static std::uint32_t matching_chunk_mask(const Node& node,
                                           std::size_t slot,
                                           std::size_t lane_count,
                                           std::int64_t target_in_node,
                                           bool include_zero_boundary) {
#ifdef PIXIE_AVX2_SUPPORT
    if constexpr (std::is_same_v<typename Node::ExcessType, std::int16_t>) {
      if (lane_count == 16 &&
          target_in_node >= std::numeric_limits<std::int16_t>::min() &&
          target_in_node <= std::numeric_limits<std::int16_t>::max()) {
        alignas(32) std::int16_t prefix_before[16]{};
        fill_prefix_before(node, slot, prefix_before);
        return rmm_btree_match_mask_i16x16(
            prefix_before, node.min_excess.data() + slot,
            node.max_excess.data() + slot,
            static_cast<std::int16_t>(target_in_node), include_zero_boundary);
      }
    } else if constexpr (std::is_same_v<typename Node::ExcessType,
                                        std::int64_t>) {
      if (lane_count == 4) {
        alignas(32) std::int64_t prefix_before[4]{};
        fill_prefix_before(node, slot, prefix_before);
        return rmm_btree_match_mask_i64x4(
            prefix_before, node.min_excess.data() + slot,
            node.max_excess.data() + slot, target_in_node,
            include_zero_boundary);
      }
    }
#endif
    std::uint32_t result = 0;
    for (std::size_t lane = 0; lane < lane_count; ++lane) {
      const std::int64_t rel =
          target_in_node - prefix_excess_at(node, slot + lane);
      const bool found = (include_zero_boundary && rel == 0) ||
                         (node.min_excess[slot + lane] <= rel &&
                          rel <= node.max_excess[slot + lane]);
      if (found) {
        result |= std::uint32_t{1} << lane;
      }
    }
    return result;
  }

  /**
   * @brief Fill SIMD lanes with child-start prefix excess values.
   * @details Summary nodes store inclusive prefix excess through each child.
   * This helper converts that representation into exclusive prefix-before-child
   * values for the vector chunk beginning at @p slot. Lane zero is zero when
   * the chunk starts at the first child; otherwise each lane reads the previous
   * stored inclusive prefix.
   * @tparam Node LowNode or HighNode.
   * @param node Summary node containing inclusive child prefixes.
   * @param slot First child slot represented by output lane zero.
   * @param out Destination array with `vector_lane_count<Node>()` entries.
   */
  template <class Node>
  static void fill_prefix_before(const Node& node,
                                 std::size_t slot,
                                 typename Node::ExcessType* out) {
    if (slot == 0) {
      out[0] = 0;
      for (std::size_t lane = 1; lane < vector_lane_count<Node>(); ++lane) {
        out[lane] = node.prefix_excess[lane - 1];
      }
      return;
    }
    for (std::size_t lane = 0; lane < vector_lane_count<Node>(); ++lane) {
      out[lane] = node.prefix_excess[slot + lane - 1];
    }
  }

  /**
   * @brief Whether a summary can contain a forward target.
   * @details Forward search needs an interior prefix whose relative excess lies
   * in the child's `[min_excess, max_excess]` interval.
   * @param summary Child summary to test.
   * @param relative_target Target excess relative to the child start.
   * @return True if the target may occur in the child.
   */
  static bool contains_fwd(const Summary& summary,
                           std::int64_t relative_target) {
    return summary.min_excess <= relative_target &&
           relative_target <= summary.max_excess;
  }

  /**
   * @brief Whether a summary can contain a backward target or left boundary.
   * @details Backward search accepts a relative target of zero as a match at
   * the child boundary; otherwise it uses the same min/max containment as
   * forward search.
   * @param summary Child summary to test.
   * @param relative_target Target excess relative to the child start.
   * @return True if the target may occur in or at the child boundary.
   */
  static bool contains_bwd(const Summary& summary,
                           std::int64_t relative_target) {
    return relative_target == 0 || contains_fwd(summary, relative_target);
  }

  std::int64_t prefix_excess_impl(std::size_t end_position) const {
    return 2 * static_cast<std::int64_t>(rank1_impl(end_position)) -
           static_cast<std::int64_t>(end_position);
  }

  /**
   * @brief Whether a level has a non-empty parent level.
   * @details Used while ascending the btree from a block toward the root.
   * @param level Current summary level.
   * @return True when `level + 1` exists and has at least one node.
   */
  bool has_parent_level(std::size_t level) const {
    return level + 1 < total_levels() && level_count(level + 1) != 0;
  }

  /**
   * @brief Number of summary levels including the block level.
   * @details Level 0 corresponds to 512-bit blocks; higher levels are stored in
   * low/high summary-node arrays.
   * @return Number of levels currently represented.
   */
  std::size_t total_levels() const { return level_counts_.size(); }

  /**
   * @brief Number of nodes at a summary level.
   * @details Returns zero for out-of-range levels.
   * @param level Summary level to query.
   * @return Node count at @p level.
   */
  std::size_t level_count(std::size_t level) const {
    return level < level_counts_.size() ? level_counts_[level] : 0;
  }

  /**
   * @brief Reconstruct the relative summary for one node or block.
   * @details Retrieves a child summary from its parent node. For the top level,
   * returns the cached top summary.
   * @param level Summary level of the requested node.
   * @param index Node index within @p level.
   * @return Relative summary of the requested subtree.
   */
  Summary summary_at(std::size_t level, std::size_t index) const {
    if (level + 1 >= total_levels()) {
      return top_summary_;
    }
    const std::size_t parent_level = level + 1;
    const std::size_t fanout = fanout_to_parent(level);
    const std::size_t parent = index / fanout;
    const std::size_t slot = index % fanout;
    Summary summary;
    if (parent_level == 1) {
      const LowNode& node = low_levels_[0][parent];
      summary.block_excess = child_excess(node, slot);
      summary.min_excess = node.min_excess[slot];
      summary.max_excess = node.max_excess[slot];
      summary.min_count = node.min_count[slot];
    } else {
      const HighNode& node = high_levels_[parent_level - 2][parent];
      summary.block_excess = child_excess(node, slot);
      summary.min_excess = node.min_excess[slot];
      summary.max_excess = node.max_excess[slot];
      summary.min_count = node.min_count[slot];
    }
    return summary;
  }

  /**
   * @brief Fanout from a level to its parent level.
   * @details Blocks are grouped by `kLowFanout` into low nodes; all higher
   * levels use `kHighFanout`.
   * @param level Child level.
   * @return Fanout used to compute the parent of @p level.
   */
  static std::size_t fanout_to_parent(std::size_t level) {
    return level == 0 ? kLowFanout : kHighFanout;
  }

  /**
   * @brief Multiply sizes, saturating on overflow.
   * @details Returns `max_size_t` instead of overflowing when the product
   * cannot be represented.
   * @param lhs Left operand.
   * @param rhs Right operand.
   * @return Saturating product.
   */
  static std::size_t mul_clamped(std::size_t lhs, std::size_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
      return std::numeric_limits<std::size_t>::max();
    }
    return lhs * rhs;
  }

  /**
   * @brief Maximum number of bits covered by one node at a level.
   * @details Level 0 spans one 512-bit block. Higher levels multiply by the low
   * fanout once and then by the high fanout for each additional level.
   * @param level Summary level.
   * @return Maximum bit span of one node at @p level, saturated on overflow.
   */
  static std::size_t level_span_bits(std::size_t level) {
    std::size_t span = kBlockBits;
    if (level >= 1) {
      span = mul_clamped(span, kLowFanout);
    }
    if (level >= 2) {
      for (std::size_t i = 2; i <= level; ++i) {
        span = mul_clamped(span, kHighFanout);
      }
    }
    return span;
  }

  /**
   * @brief First bit position covered by a node.
   * @details Computes `index * level_span_bits(level)` and clamps the result to
   * the sequence size on overflow or beyond-end nodes.
   * @param level Summary level.
   * @param index Node index within @p level.
   * @return First covered bit position, clamped to `bit_count_`.
   */
  std::size_t node_start_bit(std::size_t level, std::size_t index) const {
    const std::size_t span = level_span_bits(level);
    if (span != 0 && index > std::numeric_limits<std::size_t>::max() / span) {
      return bit_count_;
    }
    return std::min(bit_count_, index * span);
  }

  /**
   * @brief Number of valid bits covered by a node.
   * @details Returns the full level span for interior nodes and the remaining
   * sequence length for the last partial node.
   * @param level Summary level.
   * @param index Node index within @p level.
   * @return Number of valid bits covered by the node.
   */
  std::size_t node_size_bits(std::size_t level, std::size_t index) const {
    const std::size_t start = node_start_bit(level, index);
    if (start >= bit_count_) {
      return 0;
    }
    return std::min(level_span_bits(level), bit_count_ - start);
  }

  /**
   * @brief Exclusive end bit position covered by a node.
   * @details Adds `node_size_bits(level, index)` to `node_start_bit(level,
   * index)`, clamping to the sequence size on overflow.
   * @param level Summary level.
   * @param index Node index within @p level.
   * @return Exclusive end position of the covered bit interval.
   */
  std::size_t node_end_bit(std::size_t level, std::size_t index) const {
    const std::size_t start = node_start_bit(level, index);
    const std::uint64_t size = node_size_bits(level, index);
    if (size > std::numeric_limits<std::size_t>::max() - start) {
      return bit_count_;
    }
    return std::min(bit_count_, start + static_cast<std::size_t>(size));
  }

  struct NodeRef {
    std::size_t level;
    std::size_t index;
  };

  static constexpr std::size_t kMaxCoverNodes = 512;

  struct ScanResult {
    std::int64_t block_excess = 0;
    std::int64_t min_value = std::numeric_limits<std::int64_t>::max();
    std::int64_t max_value = std::numeric_limits<std::int64_t>::min();
    std::uint64_t min_count = 0;
    std::size_t min_position = npos;
    std::size_t max_position = npos;
  };

  struct MinScanResult {
    std::int64_t block_excess = 0;
    std::int64_t min_value = std::numeric_limits<std::int64_t>::max();
    std::size_t min_position = npos;
  };

  struct RangeExtremeResult {
    std::size_t position = npos;
    std::int64_t value = 0;
  };

  struct RangeMinStats {
    std::int64_t value = 0;
    std::uint64_t count = 0;
  };

  /**
   * @brief Return the position of the first range minimum or maximum.
   * @details Wrapper around `range_extreme_query` that extracts only the
   * position.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @param find_min True for minimum, false for maximum.
   * @return Position of the first extreme, or `npos`.
   */
  std::size_t range_extreme_query_pos(std::size_t range_begin,
                                      std::size_t range_end,
                                      bool find_min) const {
    return range_extreme_query(range_begin, range_end, find_min).position;
  }

  /**
   * @brief Return the value of the range minimum or maximum.
   * @details Wrapper around `range_extreme_query` that extracts only the
   * relative excess value.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @param find_min True for minimum, false for maximum.
   * @return Extreme relative excess value.
   */
  int range_extreme_query_val(std::size_t range_begin,
                              std::size_t range_end,
                              bool find_min) const {
    return static_cast<int>(
        range_extreme_query_value(range_begin, range_end, find_min));
  }

  /**
   * @brief Return only the minimum or maximum relative excess in a range.
   * @details Uses the same decomposition as `range_extreme_query`, but skips
   * the final descent needed to recover a bit position.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @param find_min True for minimum, false for maximum.
   * @return Extreme relative excess value.
   */
  std::int64_t range_extreme_query_value(std::size_t range_begin,
                                         std::size_t range_end,
                                         bool find_min) const {
    std::int64_t value = 0;
    std::int64_t best = find_min ? std::numeric_limits<std::int64_t>::max()
                                 : std::numeric_limits<std::int64_t>::min();

    auto consider_value = [&](std::int64_t candidate) {
      if ((find_min && candidate < best) || (!find_min && candidate > best)) {
        best = candidate;
      }
    };

    const std::size_t range_end_exclusive = range_end + 1;
    const std::size_t first_full_block =
        (range_begin + kBlockBits - 1) / kBlockBits;
    const std::size_t full_begin =
        std::min(range_end_exclusive, first_full_block * kBlockBits);
    if (range_begin < full_begin) {
      if (find_min) {
        const MinScanResult scan = scan_min_range(range_begin, full_begin);
        consider_value(scan.min_value);
        value += scan.block_excess;
      } else {
        const ScanResult scan = scan_range(range_begin, full_begin);
        consider_value(scan.max_value);
        value += scan.block_excess;
      }
    }

    const std::size_t last_full_block_exclusive =
        range_end_exclusive / kBlockBits;
    const std::size_t middle_begin = full_begin;
    const std::size_t middle_end =
        std::max(middle_begin, last_full_block_exclusive * kBlockBits);
    if (middle_begin < middle_end) {
      for_each_cover_node(middle_begin, middle_end, [&](NodeRef node) {
        const Summary summary = summary_at(node.level, node.index);
        consider_value(value +
                       (find_min ? summary.min_excess : summary.max_excess));
        value += summary.block_excess;
        return true;
      });
    }

    if (middle_end < range_end_exclusive) {
      if (find_min) {
        const MinScanResult scan =
            scan_min_range(middle_end, range_end_exclusive);
        consider_value(value + scan.min_value);
      } else {
        const ScanResult scan = scan_range(middle_end, range_end_exclusive);
        consider_value(value + scan.max_value);
      }
    }
    if (best == std::numeric_limits<std::int64_t>::max() ||
        best == std::numeric_limits<std::int64_t>::min()) {
      return 0;
    }
    return best;
  }

  /**
   * @brief Compute both position and value for a range minimum or maximum.
   * @details Scans partial edge blocks directly, covers the aligned middle with
   * summary nodes, and descends only if the best extreme is represented by a
   * summary node.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @param find_min True for minimum, false for maximum.
   * @return Extreme position and value.
   */
  RangeExtremeResult range_extreme_query(std::size_t range_begin,
                                         std::size_t range_end,
                                         bool find_min) const {
    std::int64_t value = 0;
    std::int64_t best = find_min ? std::numeric_limits<std::int64_t>::max()
                                 : std::numeric_limits<std::int64_t>::min();
    std::size_t best_position = npos;
    NodeRef best_node;
    std::int64_t prefix_at_best_node = 0;
    bool best_is_node = false;

    auto consider_point = [&](std::int64_t candidate, std::size_t position) {
      if ((find_min && candidate < best) || (!find_min && candidate > best)) {
        best = candidate;
        best_position = position;
        best_is_node = false;
      }
    };

    const std::size_t range_end_exclusive = range_end + 1;
    const std::size_t first_full_block =
        (range_begin + kBlockBits - 1) / kBlockBits;
    const std::size_t full_begin =
        std::min(range_end_exclusive, first_full_block * kBlockBits);
    if (range_begin < full_begin) {
      if (find_min) {
        const MinScanResult scan = scan_min_range(range_begin, full_begin);
        consider_point(scan.min_value, scan.min_position);
        value += scan.block_excess;
      } else {
        const ScanResult scan = scan_range(range_begin, full_begin);
        consider_point(scan.max_value, scan.max_position);
        value += scan.block_excess;
      }
    }

    const std::size_t last_full_block_exclusive =
        range_end_exclusive / kBlockBits;
    const std::size_t middle_begin = full_begin;
    const std::size_t middle_end =
        std::max(middle_begin, last_full_block_exclusive * kBlockBits);
    if (middle_begin < middle_end) {
      for_each_cover_node(middle_begin, middle_end, [&](NodeRef node) {
        const Summary summary = summary_at(node.level, node.index);
        const std::int64_t candidate =
            value + (find_min ? summary.min_excess : summary.max_excess);
        if ((find_min && candidate < best) || (!find_min && candidate > best)) {
          best = candidate;
          best_node = node;
          prefix_at_best_node = value;
          best_is_node = true;
        }
        value += summary.block_excess;
        return true;
      });
    }

    if (middle_end < range_end_exclusive) {
      if (find_min) {
        const MinScanResult scan =
            scan_min_range(middle_end, range_end_exclusive);
        consider_point(value + scan.min_value, scan.min_position);
      } else {
        const ScanResult scan = scan_range(middle_end, range_end_exclusive);
        consider_point(value + scan.max_value, scan.max_position);
      }
    }

    if (best_is_node) {
      best_position =
          descend_first_extreme(best_node.level, best_node.index,
                                best - prefix_at_best_node, find_min);
    }
    return {best_position,
            best == std::numeric_limits<std::int64_t>::max() ||
                    best == std::numeric_limits<std::int64_t>::min()
                ? 0
                : best};
  }

  /**
   * @brief Compute minimum value and multiplicity in an inclusive range.
   * @details Uses the same edge-block plus aligned-cover decomposition as range
   * extrema, accumulating the number of positions attaining the best minimum.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @return Minimum relative excess and its multiplicity.
   */
  RangeMinStats range_min_stats(std::size_t range_begin,
                                std::size_t range_end) const {
    std::int64_t value = 0;
    std::int64_t best = std::numeric_limits<std::int64_t>::max();
    std::uint64_t count = 0;

    auto consider = [&](std::int64_t candidate, std::uint64_t candidate_count) {
      if (candidate < best) {
        best = candidate;
        count = candidate_count;
      } else if (candidate == best) {
        count += candidate_count;
      }
    };

    const std::size_t range_end_exclusive = range_end + 1;
    const std::size_t first_full_block =
        (range_begin + kBlockBits - 1) / kBlockBits;
    const std::size_t full_begin =
        std::min(range_end_exclusive, first_full_block * kBlockBits);
    if (range_begin < full_begin) {
      const ScanResult scan = scan_range(range_begin, full_begin);
      consider(scan.min_value, scan.min_count);
      value += scan.block_excess;
    }

    const std::size_t last_full_block_exclusive =
        range_end_exclusive / kBlockBits;
    const std::size_t middle_begin = full_begin;
    const std::size_t middle_end =
        std::max(middle_begin, last_full_block_exclusive * kBlockBits);
    if (middle_begin < middle_end) {
      for_each_cover_node(middle_begin, middle_end, [&](NodeRef node) {
        const Summary summary = summary_at(node.level, node.index);
        consider(value + summary.min_excess, summary.min_count);
        value += summary.block_excess;
        return true;
      });
    }

    if (middle_end < range_end_exclusive) {
      const ScanResult scan = scan_range(middle_end, range_end_exclusive);
      consider(value + scan.min_value, scan.min_count);
    }
    return {best == std::numeric_limits<std::int64_t>::max() ? 0 : best, count};
  }

  /**
   * @brief Select the q-th position attaining the range minimum.
   * @details Assumes @p target is the range-minimum value. Skips whole scanned
   * pieces whose minimum is not @p target or whose minimum count is before
   * @p rank, then descends/scans the selected piece.
   * @param range_begin Inclusive range start.
   * @param range_end Inclusive range end.
   * @param target Minimum relative excess value to select.
   * @param rank One-based rank among positions attaining @p target.
   * @return Selected bit position, or `npos`.
   */
  std::size_t range_min_select(std::size_t range_begin,
                               std::size_t range_end,
                               std::int64_t target,
                               std::uint64_t rank) const {
    std::int64_t value = 0;
    const std::size_t range_end_exclusive = range_end + 1;
    const std::size_t first_full_block =
        (range_begin + kBlockBits - 1) / kBlockBits;
    const std::size_t full_begin =
        std::min(range_end_exclusive, first_full_block * kBlockBits);
    if (range_begin < full_begin) {
      const ScanResult scan = scan_range(range_begin, full_begin);
      if (scan.min_value == target) {
        if (rank <= scan.min_count) {
          return qth_min_in_range(range_begin, full_begin, target, rank);
        }
        rank -= scan.min_count;
      }
      value += scan.block_excess;
    }

    const std::size_t last_full_block_exclusive =
        range_end_exclusive / kBlockBits;
    const std::size_t middle_begin = full_begin;
    const std::size_t middle_end =
        std::max(middle_begin, last_full_block_exclusive * kBlockBits);
    if (middle_begin < middle_end) {
      bool found_node = false;
      std::size_t selected = npos;
      for_each_cover_node(middle_begin, middle_end, [&](NodeRef node) {
        const Summary summary = summary_at(node.level, node.index);
        const std::int64_t candidate = value + summary.min_excess;
        if (candidate == target) {
          if (rank <= summary.min_count) {
            selected =
                descend_qth_min(node.level, node.index, target - value, rank);
            found_node = true;
            return false;
          }
          rank -= summary.min_count;
        }
        value += summary.block_excess;
        return true;
      });
      if (found_node) {
        return selected;
      }
    }

    if (middle_end < range_end_exclusive) {
      const ScanResult scan = scan_range(middle_end, range_end_exclusive);
      if (value + scan.min_value == target) {
        return qth_min_in_range(middle_end, range_end_exclusive, target - value,
                                rank);
      }
    }
    return npos;
  }

  /**
   * @brief Visit an aligned block interval as summary nodes.
   * @details Produces the same left-to-right cover as the previous materialized
   * cover helper, but streams nodes into @p callback and only keeps the right
   * boundary stack needed to preserve bit order. The callback should return
   * false to stop early.
   */
  template <class Callback>
  bool for_each_cover_node(std::size_t begin,
                           std::size_t end,
                           Callback&& callback) const {
    if (begin >= end || total_levels() == 0 || (begin % kBlockBits) != 0 ||
        (end % kBlockBits) != 0) {
      return true;
    }

    std::array<NodeRef, kMaxCoverNodes> right_nodes;
    std::size_t right_size = 0;
    std::size_t emitted = 0;
    std::size_t level = 0;
    std::size_t left = begin / kBlockBits;
    std::size_t right = end / kBlockBits;

    auto emit = [&](NodeRef node) {
      if (emitted >= kMaxCoverNodes) {
        return true;
      }
      ++emitted;
      return callback(node);
    };

    auto push_right = [&](NodeRef node) {
      if (right_size < right_nodes.size()) {
        right_nodes[right_size++] = node;
      }
    };

    while (left < right) {
      if (!has_parent_level(level)) {
        for (std::size_t index = left; index < right; ++index) {
          if (!emit({level, index})) {
            return false;
          }
        }
        break;
      }

      const std::size_t fanout = fanout_to_parent(level);
      while (left < right && (left % fanout) != 0) {
        if (!emit({level, left})) {
          return false;
        }
        ++left;
      }
      while (left < right && (right % fanout) != 0) {
        --right;
        push_right({level, right});
      }
      left /= fanout;
      right /= fanout;
      ++level;
    }

    while (right_size > 0) {
      if (!emit(right_nodes[--right_size])) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Descend to the first bit position attaining an extreme value.
   * @details Starting from a node whose relative min/max equals @p target,
   * chooses the first child preserving that target until a block is reached.
   * @param level Starting summary level.
   * @param index Node index at @p level.
   * @param target Relative target excess inside the starting node.
   * @param find_min True for minimum, false for maximum.
   * @return First bit position attaining the target, or `npos`.
   */
  std::size_t descend_first_extreme(std::size_t level,
                                    std::size_t index,
                                    std::int64_t target,
                                    bool find_min) const {
    while (level > 0) {
      const std::size_t child_level = level - 1;
      const std::size_t fanout = fanout_to_parent(child_level);
      const std::size_t child_begin = index * fanout;
      const std::size_t child_end =
          std::min(level_count(child_level), child_begin + fanout);
      std::int64_t prefix = 0;
      bool found = false;
      for (std::size_t child = child_begin; child < child_end; ++child) {
        Summary summary = summary_at(child_level, child);
        const std::int64_t candidate =
            prefix + (find_min ? summary.min_excess : summary.max_excess);
        if (candidate == target) {
          index = child;
          level = child_level;
          target -= prefix;
          found = true;
          break;
        }
        prefix += summary.block_excess;
      }
      if (!found) {
        return npos;
      }
    }
    return first_prefix_in_block(index, target);
  }

  /**
   * @brief Find the first prefix in a 512-bit block with the target excess.
   * @details Uses `excess_positions_512` for full blocks and scalar scanning
   * for partial blocks.
   * @param block_index Zero-based block index.
   * @param target Target excess relative to the beginning of the block.
   * @return First bit position whose inclusive block prefix reaches @p target,
   * or `npos`.
   */
  std::size_t first_prefix_in_block(std::size_t block_index,
                                    std::int64_t target) const {
    if (full_block_has_words(block_index) && target >= -512 && target <= 512) {
      std::uint64_t out[kBlockWords];
      excess_positions_512(bits_.data() + block_index * kBlockWords,
                           static_cast<int>(target), out);
      for (std::size_t word = 0; word < kBlockWords; ++word) {
        const std::uint64_t mask = out[word];
        if (mask != 0) {
          return block_index * kBlockBits + word * 64 + std::countr_zero(mask);
        }
      }
      return npos;
    }

    const std::size_t begin = block_index * kBlockBits;
    const std::size_t length = block_size(block_index);
    std::int64_t current = 0;
    for (std::size_t offset = 0; offset < length; ++offset) {
      current += bit(begin + offset) ? 1 : -1;
      if (current == target) {
        return begin + offset;
      }
    }
    return npos;
  }

  /**
   * @brief Descend to the q-th position attaining a minimum value.
   * @details Starting from a node known to contain occurrences of @p target,
   * skips child subtrees by minimum count until the q-th occurrence is located.
   * @param level Starting summary level.
   * @param index Node index at @p level.
   * @param target Relative minimum value inside the starting node.
   * @param rank One-based rank among occurrences of @p target.
   * @return Selected bit position, or `npos`.
   */
  std::size_t descend_qth_min(std::size_t level,
                              std::size_t index,
                              std::int64_t target,
                              std::uint64_t rank) const {
    while (level > 0) {
      const std::size_t child_level = level - 1;
      const std::size_t fanout = fanout_to_parent(child_level);
      const std::size_t child_begin = index * fanout;
      const std::size_t child_end =
          std::min(level_count(child_level), child_begin + fanout);
      std::int64_t prefix = 0;
      bool found = false;
      for (std::size_t child = child_begin; child < child_end; ++child) {
        Summary summary = summary_at(child_level, child);
        const std::int64_t candidate = prefix + summary.min_excess;
        if (candidate == target) {
          if (rank <= summary.min_count) {
            index = child;
            level = child_level;
            target -= prefix;
            found = true;
            break;
          }
          rank -= summary.min_count;
        }
        prefix += summary.block_excess;
      }
      if (!found) {
        return npos;
      }
    }
    return qth_min_in_range(index * kBlockBits,
                            index * kBlockBits + block_size(index), target,
                            rank);
  }

  /**
   * @brief Scan bits to select the q-th occurrence of a target minimum.
   * @details Performs a scalar scan over `[begin, end)`, accumulating relative
   * excess from zero.
   * @param begin Inclusive bit start.
   * @param end Exclusive bit end.
   * @param target Target relative excess.
   * @param rank One-based occurrence rank.
   * @return Selected bit position, or `npos`.
   */
  std::size_t qth_min_in_range(std::size_t begin,
                               std::size_t end,
                               std::int64_t target,
                               std::uint64_t rank) const {
    std::int64_t current = 0;
    for (std::size_t position = begin; position < end; ++position) {
      current += bit(position) ? 1 : -1;
      if (current == target && --rank == 0) {
        return position;
      }
    }
    return npos;
  }

  /**
   * @brief Byte-accelerated scan of an arbitrary bit interval.
   * @details Handles unaligned edge bits individually and consumes aligned
   * bytes via `byte_lut`, producing total excess plus min/max positions and min
   * count.
   * @param begin Inclusive bit start.
   * @param end Exclusive bit end.
   * @return Relative scan summary for `[begin, end)`.
   */
  ScanResult scan_range(std::size_t begin, std::size_t end) const {
    ScanResult result;
    const auto& lut = byte_lut();
    while (begin < end && (begin & 7) != 0) {
      append_scanned_bit(result, begin);
      ++begin;
    }
    while (begin + 8 <= end) {
      const ByteAgg& byte = lut[get_byte(begin)];
      const std::int64_t min_candidate = result.block_excess + byte.min_excess;
      if (min_candidate < result.min_value) {
        result.min_value = min_candidate;
        result.min_count = byte.min_count;
        result.min_position = begin + byte.pos_first_min;
      } else if (min_candidate == result.min_value) {
        result.min_count += byte.min_count;
      }
      const std::int64_t max_candidate = result.block_excess + byte.max_excess;
      if (max_candidate > result.max_value) {
        result.max_value = max_candidate;
        result.max_position = begin + byte.pos_first_max;
      }
      result.block_excess += byte.block_excess;
      begin += 8;
    }
    while (begin < end) {
      append_scanned_bit(result, begin);
      ++begin;
    }
    return result;
  }

  /**
   * @brief Byte-accelerated minimum-only scan of an arbitrary bit interval.
   * @details Used by range-min position/value queries that do not need maximum
   * fields or minimum multiplicity.
   */
  MinScanResult scan_min_range(std::size_t begin, std::size_t end) const {
    MinScanResult result;
    const auto& lut = byte_lut();
    while (begin < end && (begin & 7) != 0) {
      append_scanned_min_bit(result, begin);
      ++begin;
    }
    while (begin + 8 <= end) {
      const ByteAgg& byte = lut[get_byte(begin)];
      const std::int64_t min_candidate = result.block_excess + byte.min_excess;
      if (min_candidate < result.min_value) {
        result.min_value = min_candidate;
        result.min_position = begin + byte.pos_first_min;
      }
      result.block_excess += byte.block_excess;
      begin += 8;
    }
    while (begin < end) {
      append_scanned_min_bit(result, begin);
      ++begin;
    }
    return result;
  }

  /**
   * @brief Append one bit to an incremental range scan.
   * @details Updates total excess, min/max values, first min/max positions, and
   * minimum count after consuming @p position.
   * @param result Scan accumulator to update.
   * @param position Bit position to consume.
   */
  void append_scanned_bit(ScanResult& result, std::size_t position) const {
    result.block_excess += bit(position) ? 1 : -1;
    if (result.block_excess < result.min_value) {
      result.min_value = result.block_excess;
      result.min_count = 1;
      result.min_position = position;
    } else if (result.block_excess == result.min_value) {
      ++result.min_count;
    }
    if (result.block_excess > result.max_value) {
      result.max_value = result.block_excess;
      result.max_position = position;
    }
  }

  /**
   * @brief Append one bit to an incremental minimum-only range scan.
   */
  void append_scanned_min_bit(MinScanResult& result,
                              std::size_t position) const {
    result.block_excess += bit(position) ? 1 : -1;
    if (result.block_excess < result.min_value) {
      result.min_value = result.block_excess;
      result.min_position = position;
    }
  }

  /**
   * @brief Return the byte starting at a byte-aligned bit position.
   * @details Reads eight bits from the backing span by shifting the containing
   * word. Callers only use this on byte-aligned positions that do not cross a
   * word boundary.
   * @param bit_position Byte-aligned bit position.
   * @return Low eight bits starting at @p bit_position.
   */
  std::uint8_t get_byte(std::size_t bit_position) const {
    return static_cast<std::uint8_t>(
        (bits_[bit_position >> 6] >> (bit_position & 63)) & 0xffu);
  }

  /**
   * @brief Read one bit from the backing span.
   * @details Interprets the backing words as little-endian bit order.
   * @param position Bit position to read.
   * @return 0 or 1.
   */
  std::uint8_t bit(std::size_t position) const {
    return static_cast<std::uint8_t>((bits_[position >> 6] >> (position & 63)) &
                                     1ull);
  }

  std::span<const std::uint64_t> bits_;
  std::optional<BitVector> rank_index_;
  std::size_t bit_count_ = 0;
  std::size_t block_count_ = 0;
  Summary top_summary_;
  std::vector<std::size_t> level_counts_;
  std::vector<std::vector<LowNode>> low_levels_;
  std::vector<std::vector<HighNode>> high_levels_;
};

}  // namespace pixie::experimental
