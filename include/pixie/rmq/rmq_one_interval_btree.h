#pragma once

#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pixie::rmq {

/**
 * @brief Interval B-tree backend for arrays with adjacent differences +/-1.
 *
 * @details The indexed depth sequence is represented by BP deltas: bit 1 means
 * the next depth is current + 1, and bit 0 means current - 1. The structure
 * keeps 128-depth-position lower blocks separate and builds a cache-aligned
 * interval hierarchy above those blocks. Queries scan partial edge blocks with
 * the existing 128-bit excess-min primitives and cover the aligned middle with
 * B-tree summary nodes.
 *
 * @tparam Index Unsigned integer type used to validate representable
 * positions.
 * @tparam HighCacheLines Number of 64-byte cache lines assigned to one high
 * summary node.
 * @tparam LowFanout Number of lower-block summaries stored in one low node.
 * @tparam UseBoundaryRecords Whether to precompute left-to-right and
 * right-to-left block minimum record positions for cross-block edge queries.
 */
template <class Index = std::size_t,
          std::size_t HighCacheLines = 2,
          std::size_t LowFanout = 32,
          bool UseBoundaryRecords = false>
class OneIntervalBTreeRmq {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "OneIntervalBTreeRmq index type must be unsigned");
  static_assert(HighCacheLines > 0);
  static_assert(LowFanout > 0);

  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kBlockSize = 128;
  static constexpr std::size_t kChunkSize = 128;
  static constexpr std::size_t kBlockChunks = kBlockSize / kChunkSize;
  static constexpr std::size_t kChunkWords = kChunkSize / 64;
  static constexpr bool kUseBoundaryRecords = UseBoundaryRecords;
  static constexpr std::size_t kCacheLineBytes = 64;
  static constexpr std::size_t kLowFanout = LowFanout;
  static constexpr std::size_t kHighFanout =
      std::max<std::size_t>(2, (512 * HighCacheLines) / (2 * 64));
  static constexpr std::size_t kMaxFanout = std::max(kLowFanout, kHighFanout);

  static_assert(kBlockSize % kChunkSize == 0);
  static_assert(kChunkSize == 128);
  static_assert(kMaxFanout <= 64);
  static_assert(
      kBlockSize * kLowFanout <=
      static_cast<std::size_t>(std::numeric_limits<std::int16_t>::max()));

  /**
   * @brief Construct an empty +/-1 RMQ index.
   */
  OneIntervalBTreeRmq() = default;
  OneIntervalBTreeRmq(const OneIntervalBTreeRmq&) = default;
  OneIntervalBTreeRmq(OneIntervalBTreeRmq&&) noexcept = default;
  OneIntervalBTreeRmq& operator=(const OneIntervalBTreeRmq&) = default;
  OneIntervalBTreeRmq& operator=(OneIntervalBTreeRmq&&) noexcept = default;

  /**
   * @brief Build an interval B-tree over a BP-encoded depth-delta sequence.
   *
   * @details The bit span is not copied and must outlive this object. A
   * sequence with @p depth_count depth positions requires
   * @p depth_count - 1 delta bits.
   *
   * @param bits Packed delta bits in little-endian bit order within each word.
   * @param depth_count Number of depth positions indexed by the RMQ.
   * @throws std::length_error if @p Index cannot represent all positions.
   * @throws std::invalid_argument if @p bits does not contain enough words.
   */
  OneIntervalBTreeRmq(std::span<const std::uint64_t> bits,
                      std::size_t depth_count)
      : input_bits_(bits), depth_count_(depth_count) {
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
   * @brief Return the first minimum depth position in [@p left, @p right).
   *
   * @details The query range is half-open over depth positions, not delta-bit
   * positions. Empty or invalid ranges return `npos`.
   */
  std::size_t arg_min(std::size_t left, std::size_t right) const {
    if (left >= right || right > depth_count_) {
      return npos;
    }

    const std::size_t last = right - 1;
    const std::size_t left_block = left / kBlockSize;
    const std::size_t right_block = last / kBlockSize;
    if (left_block == right_block) {
      return scan_block_range_position(left_block, left % kBlockSize,
                                       last % kBlockSize);
    }

    Candidate best;
    std::int64_t prefix = 0;
    NodeRef best_node;
    std::int64_t best_node_target = 0;
    bool best_is_node = false;

    const auto consider_position = [&](std::size_t position,
                                       std::int64_t value) {
      if (position == npos) {
        return;
      }
      if (value < best.value) {
        best = {position, value};
        best_is_node = false;
      }
    };

    const auto consider_node = [&](NodeRef node, std::int64_t target,
                                   std::int64_t value) {
      if (value < best.value) {
        best = {npos, value};
        best_node = node;
        best_node_target = target;
        best_is_node = true;
      }
    };

    const auto consider_excess_result =
        [&](std::size_t block, ExcessResult result, std::int64_t value) {
          if (result.offset >= block_size(block)) {
            return;
          }
          consider_position(block * kBlockSize + result.offset, value);
        };

    const std::size_t left_offset = left % kBlockSize;
    const std::size_t right_offset = last % kBlockSize;
    ExcessResult fused_right_boundary;
    bool has_fused_right_boundary = false;

    const std::size_t first_full_block = (left + kBlockSize - 1) / kBlockSize;
    const std::size_t full_begin =
        std::min(right, first_full_block * kBlockSize);
    if constexpr (UseBoundaryRecords) {
      if (left < full_begin) {
        const ScanResult scan =
            scan_suffix_boundary_record(left_block, left_offset);
        consider_position(scan.position, scan.min_value);
        prefix += scan.block_excess;
      }
    } else {
      const bool use_fused_boundaries = left_offset > right_offset &&
                                        left_offset < kBlockSize - 1 &&
                                        right_offset > 0;
      if (use_fused_boundaries) {
        const auto left_bits = block_bits(left_block);
        const auto right_bits = block_bits(right_block);
        const ExcessBoundaryPairResult result =
            excess_min_128_disjoint_suffix_prefix(
                left_bits.data(), left_offset, right_bits.data(), right_offset);
        const std::int64_t left_base =
            block_prefix_excess(left_block, left_offset);
        consider_excess_result(left_block, result.suffix,
                               result.suffix.min_excess - left_base);
        prefix += block_prefix_excess(left_block, kBlockSize) - left_base;
        fused_right_boundary = result.prefix;
        has_fused_right_boundary = true;
      } else if (left < full_begin) {
        const ScanResult scan = scan_range(left, full_begin);
        consider_position(scan.position, scan.min_value);
        prefix += scan.block_excess;
      }
    }

    const std::size_t last_full_block_exclusive = right / kBlockSize;
    const std::size_t middle_begin = full_begin;
    const std::size_t middle_end =
        std::max(middle_begin, last_full_block_exclusive * kBlockSize);
    if (middle_begin < middle_end) {
      Cover cover;
      collect_cover(middle_begin, middle_end, cover);
      for (std::size_t i = 0; i < cover.size; ++i) {
        const NodeRef node = cover.nodes[i];
        const Summary summary = summary_at(node.level, node.index);
        consider_node(node, summary.min_excess, prefix + summary.min_excess);
        prefix += summary.block_excess;
      }
    }

    if (middle_end < right) {
      if constexpr (UseBoundaryRecords) {
        const ExcessResult result =
            scan_prefix_boundary_record(right_block, right_offset);
        consider_excess_result(right_block, result, prefix + result.min_excess);
      } else {
        if (has_fused_right_boundary) {
          consider_excess_result(right_block, fused_right_boundary,
                                 prefix + fused_right_boundary.min_excess);
        } else {
          const ScanResult scan = scan_range(middle_end, right);
          consider_position(scan.position, prefix + scan.min_value);
        }
      }
    }

    if (best_is_node) {
      return descend_first_min(best_node.level, best_node.index,
                               best_node_target);
    }
    return best.position;
  }

 private:
  struct Summary {
    std::uint64_t size_positions = 0;
    std::int64_t block_excess = 0;
    std::int64_t min_excess = 0;
  };

  struct BlockSummary {
    Summary summary;
    std::uint16_t min_offset = 0;
  };

  struct BoundaryRecords {
    std::array<std::uint64_t, kChunkWords> prefix_records{};
    std::array<std::uint64_t, kChunkWords> suffix_records{};
  };

  static_assert(sizeof(BoundaryRecords) == 4 * sizeof(std::uint64_t));

  template <class Excess, std::size_t Fanout>
  struct alignas(kCacheLineBytes) SummaryNode {
    using ExcessType = Excess;
    static constexpr std::size_t kFanout = Fanout;

    std::array<Excess, Fanout> prefix_excess{};
    std::array<Excess, Fanout> min_excess{};
  };

  using LowNode = SummaryNode<std::int16_t, kLowFanout>;
  using HighNode = SummaryNode<std::int64_t, kHighFanout>;
  static_assert(alignof(LowNode) == kCacheLineBytes);
  static_assert(alignof(HighNode) == kCacheLineBytes);
  static_assert(sizeof(LowNode) % kCacheLineBytes == 0);
  static_assert(sizeof(HighNode) % kCacheLineBytes == 0);

  struct Candidate {
    std::size_t position = npos;
    std::int64_t value = std::numeric_limits<std::int64_t>::max();
  };

  struct ScanResult {
    std::size_t position = npos;
    std::int64_t min_value = std::numeric_limits<std::int64_t>::max();
    std::int64_t block_excess = 0;
  };

  struct NodeRef {
    std::size_t level = 0;
    std::size_t index = 0;
  };

  struct ChildSearchResult {
    bool found = false;
    std::size_t index = 0;
    std::int64_t target = 0;
  };

  static constexpr std::size_t kMaxCoverItems = 512;

  struct Cover {
    std::array<NodeRef, kMaxCoverItems> nodes{};
    std::size_t size = 0;

    void push(NodeRef node) {
      if (size < nodes.size()) {
        nodes[size++] = node;
      }
    }
  };

  /**
   * @brief Build lower block summaries and interval B-tree levels.
   */
  void build() {
    level_counts_.clear();
    low_levels_.clear();
    high_levels_.clear();
    block_summaries_.clear();
    boundary_records_.clear();
    top_summary_ = Summary{};
    block_count_ = 0;

    if (depth_count_ == 0) {
      return;
    }
    if (depth_count_ > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("RMQ +/-1 interval B-tree index is too small");
    }
    const std::size_t delta_count = depth_count_ - 1;
    if (input_bits_.size() < (delta_count + 63) / 64) {
      throw std::invalid_argument(
          "RMQ +/-1 interval B-tree bit span is too small");
    }

    block_count_ = (depth_count_ + kBlockSize - 1) / kBlockSize;
    block_summaries_.resize(block_count_);
    if constexpr (UseBoundaryRecords) {
      boundary_records_.resize(block_count_);
    }
    std::vector<Summary> block_summaries(block_count_);
    for (std::size_t block = 0; block < block_count_; ++block) {
      block_summaries_[block] = summarize_block(block);
      if constexpr (UseBoundaryRecords) {
        boundary_records_[block] = build_boundary_records(block);
      }
      block_summaries[block] = block_summaries_[block].summary;
    }
    build_levels(block_summaries);
  }

  /**
   * @brief Build low and high summary levels from 128-position blocks.
   */
  void build_levels(const std::vector<Summary>& block_summaries) {
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
  }

  BlockSummary summarize_block(std::size_t block) const {
    const std::size_t size = block_size(block);
    BlockSummary block_summary;
    Summary& summary = block_summary.summary;
    summary.size_positions = size;
    if (size == 0) {
      return block_summary;
    }
    const auto bits = block_bits(block);
    const ExcessResult scan = excess_min_128(bits.data(), 0, size - 1);
    summary.min_excess = scan.min_excess;
    block_summary.min_offset = static_cast<std::uint16_t>(scan.offset);
    summary.block_excess =
        block_prefix_excess(block, next_block_delta_count(block));
    return block_summary;
  }

  BoundaryRecords build_boundary_records(std::size_t block) const {
    const std::size_t size = block_size(block);
    BoundaryRecords records;
    if (size == 0) {
      return records;
    }

    const auto bits = block_bits(block);
    std::array<std::int16_t, kBlockSize> local_excess{};
    std::int16_t current = 0;
    std::int16_t prefix_best = 0;
    set_record_bit(records.prefix_records, 0);

    for (std::size_t offset = 1; offset < size; ++offset) {
      current += local_bit(bits, offset - 1) ? 1 : -1;
      local_excess[offset] = current;
      if (current < prefix_best) {
        prefix_best = current;
        set_record_bit(records.prefix_records, offset);
      }
    }

    std::int16_t suffix_best = std::numeric_limits<std::int16_t>::max();
    for (std::size_t offset = size; offset-- > 0;) {
      if (local_excess[offset] <= suffix_best) {
        suffix_best = local_excess[offset];
        set_record_bit(records.suffix_records, offset);
      }
    }
    return records;
  }

  static Summary append(Summary left, const Summary& right) {
    if (left.size_positions == 0) {
      return right;
    }
    if (right.size_positions == 0) {
      return left;
    }
    Summary result;
    result.size_positions = left.size_positions + right.size_positions;
    result.block_excess = left.block_excess + right.block_excess;
    result.min_excess =
        std::min(left.min_excess, left.block_excess + right.min_excess);
    return result;
  }

  std::size_t block_size(std::size_t block) const {
    const std::size_t begin = block * kBlockSize;
    return begin >= depth_count_ ? 0
                                 : std::min(kBlockSize, depth_count_ - begin);
  }

  std::size_t next_block_delta_count(std::size_t block) const {
    const std::size_t begin = block * kBlockSize;
    if (begin >= depth_count_ - 1) {
      return 0;
    }
    const std::size_t next_begin = begin + kBlockSize;
    return std::min(next_begin, depth_count_ - 1) - begin;
  }

  std::uint64_t word_or_zero(std::size_t word) const {
    return word < input_bits_.size() ? input_bits_[word] : 0;
  }

  std::array<std::uint64_t, 2> block_bits(std::size_t block) const {
    const std::size_t first_word = (block * kBlockSize) / 64;
    return {word_or_zero(first_word), word_or_zero(first_word + 1)};
  }

  std::array<std::uint64_t, 2> chunk_bits(std::size_t block,
                                          std::size_t chunk) const {
    const std::size_t first_word =
        (block * kBlockSize + chunk * kChunkSize) / 64;
    return {word_or_zero(first_word), word_or_zero(first_word + 1)};
  }

  std::int64_t block_prefix_excess(std::size_t block,
                                   std::size_t prefix_offset) const {
    prefix_offset = std::min(prefix_offset, next_block_delta_count(block));
    const auto bits = block_bits(block);
    return prefix_excess_128(bits.data(), prefix_offset);
  }

  static bool local_bit(const std::array<std::uint64_t, kChunkWords>& bits,
                        std::size_t offset) {
    return ((bits[offset >> 6] >> (offset & 63)) & 1u) != 0;
  }

  static void set_record_bit(std::array<std::uint64_t, kChunkWords>& records,
                             std::size_t offset) {
    records[offset >> 6] |= std::uint64_t{1} << (offset & 63);
  }

  ExcessResult scan_prefix_boundary_record(std::size_t block,
                                           std::size_t right_offset) const {
    const std::size_t size = block_size(block);
    if (size == 0 || block >= boundary_records_.size()) {
      return {};
    }
    right_offset = std::min(right_offset, size - 1);
    const BoundaryRecords& records = boundary_records_[block];

    std::uint64_t candidates = records.prefix_records[right_offset >> 6] &
                               first_bits_mask((right_offset & 63) + 1);
    if (candidates != 0) {
      const std::size_t offset = (right_offset & ~std::size_t{63}) +
                                 (63 - std::countl_zero(candidates));
      const auto bits = block_bits(block);
      return {prefix_excess_128(bits.data(), offset), offset};
    }
    if ((right_offset >> 6) != 0 && records.prefix_records[0] != 0) {
      const std::size_t offset =
          63 - std::countl_zero(records.prefix_records[0]);
      const auto bits = block_bits(block);
      return {prefix_excess_128(bits.data(), offset), offset};
    }
    return {};
  }

  ExcessResult scan_suffix_boundary_record_result(
      std::size_t block,
      std::size_t left_offset) const {
    const std::size_t size = block_size(block);
    if (size == 0 || block >= boundary_records_.size()) {
      return {};
    }
    left_offset = std::min(left_offset, size - 1);
    const BoundaryRecords& records = boundary_records_[block];

    const std::size_t word = left_offset >> 6;
    std::uint64_t candidates =
        records.suffix_records[word] & ~first_bits_mask(left_offset & 63);
    if (candidates != 0) {
      const std::size_t offset =
          word * 64 + static_cast<std::size_t>(std::countr_zero(candidates));
      const auto bits = block_bits(block);
      return {prefix_excess_128(bits.data(), offset), offset};
    }
    if (word == 0 && records.suffix_records[1] != 0) {
      const std::size_t offset =
          64 +
          static_cast<std::size_t>(std::countr_zero(records.suffix_records[1]));
      const auto bits = block_bits(block);
      return {prefix_excess_128(bits.data(), offset), offset};
    }
    return {};
  }

  ScanResult scan_suffix_boundary_record(std::size_t block,
                                         std::size_t left_offset) const {
    const std::size_t size = block_size(block);
    if (size == 0) {
      return {};
    }
    left_offset = std::min(left_offset, size - 1);
    const std::int64_t base = block_prefix_excess(block, left_offset);

    ScanResult scan;
    const ExcessResult result =
        scan_suffix_boundary_record_result(block, left_offset);
    if (result.offset == npos || result.offset >= size) {
      return scan;
    }

    scan.position = block * kBlockSize + result.offset;
    scan.min_value = result.min_excess - base;
    scan.block_excess = block_prefix_excess(block, kBlockSize) - base;
    return scan;
  }

  ScanResult scan_range(std::size_t begin, std::size_t end) const {
    const std::size_t block = begin / kBlockSize;
    const std::size_t block_begin = block * kBlockSize;
    return scan_block_range(block, begin - block_begin, end - block_begin - 1);
  }

  std::size_t scan_block_range_position(std::size_t block,
                                        std::size_t left_offset,
                                        std::size_t right_offset) const {
    const std::size_t size = block_size(block);
    right_offset = std::min(right_offset, size - 1);
    const auto bits = block_bits(block);
    const ExcessResult result =
        excess_min_128(bits.data(), left_offset, right_offset);
    if (result.offset == npos || result.offset >= size) {
      return npos;
    }
    return block * kBlockSize + result.offset;
  }

  ScanResult scan_block_range(std::size_t block,
                              std::size_t left_offset,
                              std::size_t right_offset) const {
    const std::size_t size = block_size(block);
    right_offset = std::min(right_offset, size - 1);
    const std::int64_t base = block_prefix_excess(block, left_offset);

    ScanResult scan;
    const auto bits = block_bits(block);
    const ExcessResult result =
        excess_min_128(bits.data(), left_offset, right_offset);
    if (result.offset == npos || result.offset >= size) {
      return scan;
    }

    scan.position = block * kBlockSize + result.offset;
    scan.min_value = result.min_excess - base;
    scan.block_excess = block_prefix_excess(block, right_offset + 1) - base;
    return scan;
  }

  std::size_t total_levels() const { return level_counts_.size(); }

  std::size_t level_count(std::size_t level) const {
    return level < level_counts_.size() ? level_counts_[level] : 0;
  }

  Summary summary_at(std::size_t level, std::size_t index) const {
    if (level == 0) {
      return block_summaries_[index].summary;
    }
    if (level + 1 >= total_levels()) {
      return top_summary_;
    }
    const std::size_t parent_level = level + 1;
    const std::size_t fanout = fanout_to_parent(level);
    const std::size_t parent = index / fanout;
    const std::size_t slot = index % fanout;
    Summary summary;
    summary.size_positions = node_size_positions(level, index);
    if (parent_level == 1) {
      const LowNode& node = low_levels_[0][parent];
      summary.block_excess = child_excess(node, slot);
      summary.min_excess = node.min_excess[slot];
    } else {
      const HighNode& node = high_levels_[parent_level - 2][parent];
      summary.block_excess = child_excess(node, slot);
      summary.min_excess = node.min_excess[slot];
    }
    return summary;
  }

  static std::size_t fanout_to_parent(std::size_t level) {
    return level == 0 ? kLowFanout : kHighFanout;
  }

  static std::size_t mul_clamped(std::size_t lhs, std::size_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs) {
      return std::numeric_limits<std::size_t>::max();
    }
    return lhs * rhs;
  }

  static std::size_t level_span_positions(std::size_t level) {
    std::size_t span = kBlockSize;
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

  std::size_t node_start_position(std::size_t level, std::size_t index) const {
    const std::size_t span = level_span_positions(level);
    if (span != 0 && index > std::numeric_limits<std::size_t>::max() / span) {
      return depth_count_;
    }
    return std::min(depth_count_, index * span);
  }

  std::size_t node_size_positions(std::size_t level, std::size_t index) const {
    const std::size_t start = node_start_position(level, index);
    if (start >= depth_count_) {
      return 0;
    }
    return std::min(level_span_positions(level), depth_count_ - start);
  }

  template <class Node>
  static std::int64_t prefix_excess_at(const Node& node, std::size_t slot) {
    if (slot == 0) {
      return 0;
    }
    return static_cast<std::int64_t>(node.prefix_excess[slot - 1]);
  }

  template <class Node>
  static std::int64_t child_excess(const Node& node, std::size_t slot) {
    return static_cast<std::int64_t>(node.prefix_excess[slot]) -
           prefix_excess_at(node, slot);
  }

  template <class Node>
  ChildSearchResult find_first_child_with_min(const Node& node,
                                              std::size_t child_level,
                                              std::size_t parent,
                                              std::size_t child_count,
                                              std::int64_t target) const {
    for (std::size_t slot = 0; slot < child_count; ++slot) {
      const std::int64_t candidate =
          prefix_excess_at(node, slot) +
          static_cast<std::int64_t>(node.min_excess[slot]);
      if (candidate == target) {
        return {true, parent * fanout_to_parent(child_level) + slot,
                static_cast<std::int64_t>(node.min_excess[slot])};
      }
    }
    return {};
  }

  ChildSearchResult find_first_child_with_min(std::size_t child_level,
                                              std::size_t parent,
                                              std::size_t child_count,
                                              std::int64_t target) const {
    if (child_level == 0) {
      return find_first_child_with_min(low_levels_[0][parent], child_level,
                                       parent, child_count, target);
    }
    return find_first_child_with_min(high_levels_[child_level - 1][parent],
                                     child_level, parent, child_count, target);
  }

  void collect_cover(std::size_t begin, std::size_t end, Cover& out) const {
    if (begin >= end || total_levels() == 0 || (begin % kBlockSize) != 0 ||
        (end % kBlockSize) != 0) {
      return;
    }

    Cover right_cover;
    std::size_t level = 0;
    std::size_t left = begin / kBlockSize;
    std::size_t right = end / kBlockSize;

    while (left < right) {
      if (!has_parent_level(level)) {
        for (std::size_t index = left; index < right; ++index) {
          out.push({level, index});
        }
        break;
      }

      const std::size_t fanout = fanout_to_parent(level);
      while (left < right && (left % fanout) != 0) {
        out.push({level, left});
        ++left;
      }
      while (left < right && (right % fanout) != 0) {
        --right;
        right_cover.push({level, right});
      }
      left /= fanout;
      right /= fanout;
      ++level;
    }

    while (right_cover.size > 0) {
      out.push(right_cover.nodes[--right_cover.size]);
    }
  }

  bool has_parent_level(std::size_t level) const {
    return level + 1 < total_levels() && level_count(level + 1) != 0;
  }

  std::size_t descend_first_min(std::size_t level,
                                std::size_t index,
                                std::int64_t target) const {
    while (level > 0) {
      const std::size_t child_level = level - 1;
      const std::size_t fanout = fanout_to_parent(child_level);
      const std::size_t child_begin = index * fanout;
      const std::size_t child_end =
          std::min(level_count(child_level), child_begin + fanout);
      const ChildSearchResult child = find_first_child_with_min(
          child_level, index, child_end - child_begin, target);
      if (!child.found) {
        return npos;
      }
      index = child.index;
      level = child_level;
      target = child.target;
    }

    if (index >= block_summaries_.size()) {
      return npos;
    }
    const BlockSummary& block = block_summaries_[index];
    return block.summary.min_excess == target
               ? index * kBlockSize + block.min_offset
               : npos;
  }

  std::span<const std::uint64_t> input_bits_;
  std::size_t depth_count_ = 0;
  std::size_t block_count_ = 0;
  Summary top_summary_;
  std::vector<BlockSummary> block_summaries_;
  std::vector<BoundaryRecords> boundary_records_;
  std::vector<std::size_t> level_counts_;
  std::vector<std::vector<LowNode>> low_levels_;
  std::vector<std::vector<HighNode>> high_levels_;
};

template <class Index = std::size_t,
          std::size_t HighCacheLines = 4,
          std::size_t LowFanout = 32>
using OneIntervalBTreeRmqBoundaryRecords =
    OneIntervalBTreeRmq<Index, HighCacheLines, LowFanout, true>;

}  // namespace pixie::rmq
