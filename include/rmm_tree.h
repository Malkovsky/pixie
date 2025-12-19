#pragma once
#include <immintrin.h>

#include <algorithm>
#include <array>
#include <bit>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "bits.h"

namespace pixie {
/**
 * @brief Range min–max tree over a bitvector (LSB-first) tailored for
 * balanced-parentheses (BP).
 * @details
 *  Implements the classic rmq/rMq family and BP navigation in O(log n)
 *  using a perfectly balanced binary tree of blocks. Supported queries:
 *  - rank1 / rank0
 *  - select1 / select0
 *  - rank10 / select10 (starts of "10")
 *  - excess (prefix +1 for '1', −1 for '0')
 *  - fwdsearch / bwdsearch (prefix–sum search)
 *  - range_min_query_pos / range_min_query_val (first minimum within a range)
 *  - range_max_query_pos / range_max_query_val (first maximum within a range)
 *  - mincount / minselect (count/selection of minima)
 *  - close / open / enclose (BP navigation)
 *
 *  The bitvector is LSB-first inside each 64-bit word.
 */
class RmMTree {
  // ------------ bitvector ------------
  std::vector<std::uint64_t> bits;  // LSB-first
  size_t num_bits = 0;              // number of bits

  // ------------ blocking ------------
  size_t block_bits = 64;  // block size (bits), leaf covers <= block_bits bits
  size_t leaf_count = 0;   // #leaves = ceil(num_bits/block_bits)

  // ------------ tree arrays (heap order: 1 is root) ------------
  // size of segment (in bits) covered by node
  // needed for: rank1/rank0, select1/select0, select10,
  //             excess, fwdsearch/bwdsearch/close/open/enclose,
  //             range_min_query/range_max_query, minselect.
  std::vector<uint32_t> segment_size_bits;

  // node_total_excess = total excess (+1 for '1', -1 for '0') on the node
  // needed for: rank1/rank0, select1/select0, excess,
  //             fwdsearch/bwdsearch/close/open/enclose,
  //             range_min_query/range_max_query, mincount/minselect.
  std::vector<int32_t> node_total_excess;

  // node_min_prefix_excess = minimum pref-excess on the node (from 0)
  // needed for: fwdsearch/bwdsearch/close/open/enclose, range_min_query,
  // mincount/minselect.
  std::vector<int32_t> node_min_prefix_excess;

  // node_max_prefix_excess = maximum pref-excess on the node (from 0)
  // needed for: fwdsearch/bwdsearch/close/open/enclose, range_max_query.
  std::vector<int32_t> node_max_prefix_excess;

  // node_min_count = number of positions where the minimum is attained
  // needed for: mincount/minselect.
  std::vector<uint32_t> node_min_count;

  // node_pattern10_count = # of "10" pattern occurrences inside the node
  // needed for: rank10, select10.
  std::vector<uint32_t> node_pattern10_count;

  // node_first_bit = first bit (0/1), node_last_bit = last bit (0/1) of the
  // segment (to handle "10" crossing) both needed for: rank10, select10.
  std::vector<uint8_t> node_first_bit, node_last_bit;

  // leaf_prefix[k * leaf_prefix_stride + j] =
  // excess(block_begin + j) - excess(block_begin), j in [0..leaf_size],
  // block_begin = k * block_bits. needed for:
  // fwdsearch/bwdsearch/close/open/enclose
  std::vector<int16_t> leaf_prefix;
  size_t leaf_prefix_stride = 0;

 public:
  /**
   * @brief Sentinel for "not found".
   */
  static constexpr size_t npos = std::numeric_limits<size_t>::max();

#ifdef DEBUG
  float built_overhead = 0.0;
#endif

  // --------- construction ----------

  /**
   * @brief Construct empty structure.
   */
  RmMTree() = default;

  /**
   * @brief Build from a '0'/'1' string.
   * @param bit_string Bitstring (characters '0' and '1').
   * @param leaf_block_bits Desired leaf size (power of two, 0 = auto).
   * @param max_overhead Max allowed overhead fraction (<0 to disable
   * constraint).
   * @details Block size priority: (1) respect @p max_overhead, (2) explicit @p
   * leaf_block_bits, (3) set to ceil_pow2(log2(num_bits)).
   */
  explicit RmMTree(const std::string& bit_string,
                   const size_t& leaf_block_bits /*0=auto*/ = 0,
                   const float& max_overhead /*<0=off*/ = -1.0) {
    build_from_string(bit_string, leaf_block_bits, max_overhead);
  }

  /**
   * @brief Build from 64-bit words (LSB-first).
   * @param words Array of words holding bits LSB-first.
   * @param bit_count Number of valid bits.
   * @param leaf_block_bits Desired leaf size (power of two, 0 = auto).
   * @param max_overhead Max allowed overhead fraction (<0 to disable
   * constraint).
   * @details Block size priority: (1) respect @p max_overhead, (2) explicit @p
   * leaf_block_bits, (3) set to ceil_pow2(log2(num_bits)).
   */
  explicit RmMTree(const std::vector<std::uint64_t>& words,
                   size_t bit_count,
                   const size_t& leaf_block_bits /*0=auto*/ = 0,
                   const float& max_overhead /*<0=off*/ = -1.0) {
    build_from_words(words, bit_count, leaf_block_bits, max_overhead);
  }

  // --------- queries: rank/select/excess ----------

  /**
   * @brief Number of ones in prefix [0, @p end_position).
   * @details Returns 0 for @p end_position == 0.
   */
  size_t rank1(const size_t& end_position) const {
    if (end_position == 0) {
      return 0;
    }
    const size_t block_index = block_of(end_position - 1);
    size_t ones_count = 0;
    if (block_index > 0) {
      size_t nodes_buffer[64];
      const size_t node_count =
          cover_blocks_collect(0, block_index - 1, nodes_buffer);
      for (size_t j = 0; j < node_count; ++j) {
        ones_count += ones_in_node(nodes_buffer[j]);
      }
    }
    const size_t block_begin = block_index * block_bits;
    const size_t block_end = std::min(num_bits, block_begin + block_bits);
    ones_count +=
        rank1_in_block(block_begin, std::min(end_position, block_end));
    return ones_count;
  }

  /**
   * @brief Number of zeros in prefix [0, @p end_position).
   * @details Computed as @p end_position - rank1(@p end_position).
   */
  size_t rank0(const size_t& end_position) const {
    return end_position - rank1(end_position);
  }

  /**
   * @brief 1-based select of the @p target_one_rank-th one.
   * @return Position of @p target_one_rank-th '1' or npos if not found.
   */
  size_t select1(size_t target_one_rank) const {
    if (target_one_rank == 0 || num_bits == 0) {
      return npos;
    }
    size_t node_index = 1;
    if (ones_in_node(node_index) < target_one_rank) {
      return npos;
    }
    size_t segment_base = 0;
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      const uint32_t ones_in_left_child = ones_in_node(left_child);
      if (ones_in_left_child >= target_one_rank) {
        node_index = left_child;
      } else {
        target_one_rank -= ones_in_left_child;
        segment_base += segment_size_bits[left_child];
        node_index = right_child;
      }
    }
    return select1_in_block(
        segment_base,
        std::min(segment_base + segment_size_bits[node_index], num_bits),
        target_one_rank);
  }

  /**
   * @brief 1-based select of the @p target_zero_rank-th zero.
   * @return Position of @p target_zero_rank-th '0' or npos if not found.
   */
  size_t select0(size_t target_zero_rank) const {
    if (target_zero_rank == 0 || num_bits == 0) {
      return npos;
    }
    size_t node_index = 1;
    const auto zeros_in_node = [&](const size_t& node) noexcept {
      return segment_size_bits[node] - ones_in_node(node);
    };
    if (zeros_in_node(node_index) < target_zero_rank) {
      return npos;
    }
    size_t segment_base = 0;
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      const size_t zeros_in_left_child = zeros_in_node(left_child);
      if (zeros_in_left_child >= target_zero_rank) {
        node_index = left_child;
      } else {
        target_zero_rank -= zeros_in_left_child;
        segment_base += segment_size_bits[left_child];
        node_index = right_child;
      }
    }
    return select0_in_block(
        segment_base,
        std::min(segment_base + segment_size_bits[node_index], num_bits),
        target_zero_rank);
  }

  /**
   * @brief Rank of the pattern "10" (starts) within [0, @p end_position).
   * @details Counts p where bit[p]==1 and bit[p+1]==0 with p+1<@p
   * end_position.
   */
  size_t rank10(const size_t& end_position) const {
    if (end_position <= 1) {
      return 0;
    }
    const size_t block_index = block_of(end_position - 1);
    size_t pattern_count = 0;
    int previous_last_bit = -1;

    if (block_index > 0) {
      const auto covered_nodes = cover_blocks(0, block_index - 1);
      for (const size_t& node_index : covered_nodes) {
        pattern_count += node_pattern10_count[node_index];
        if (previous_last_bit != -1 && previous_last_bit == 1 &&
            node_first_bit[node_index] == 0) {
          ++pattern_count;
        }
        previous_last_bit = node_last_bit[node_index];
      }
    }
    const size_t block_begin = block_index * block_bits;
    pattern_count += rr_in_block(block_begin, end_position);
    // boundary between the last full node and the leaf tail
    if (block_index > 0 && end_position > block_begin &&
        previous_last_bit == 1 && bit(block_begin) == 0) {
      ++pattern_count;
    }
    return pattern_count;
  }

  /**
   * @brief 1-based select of the @p target_pattern_rank-th "10" start.
   * @return Position p such that bits[p..p+1]=="10", or npos if not found.
   */
  size_t select10(size_t target_pattern_rank) const {
    if (target_pattern_rank == 0 || num_bits == 0) {
      return npos;
    }
    size_t node_index = 1;
    if (node_pattern10_count[node_index] < target_pattern_rank) {
      return npos;
    }
    size_t segment_base = 0;
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1, right_child = left_child | 1;
      const size_t crossing_pattern =
          (node_last_bit[left_child] == 1 && node_first_bit[right_child] == 0)
              ? 1u
              : 0u;
      if (node_pattern10_count[left_child] >= target_pattern_rank) {
        node_index = left_child;
        continue;
      }
      size_t remaining_rank =
          target_pattern_rank - node_pattern10_count[left_child];
      if (crossing_pattern) {
        if (remaining_rank == 1) {
          return segment_base + segment_size_bits[left_child] - 1;
        }
        --remaining_rank;
      }
      segment_base += segment_size_bits[left_child];
      node_index = right_child;
      target_pattern_rank = remaining_rank;
    }
    return select10_in_block(
        segment_base,
        std::min(segment_base + segment_size_bits[node_index], num_bits),
        target_pattern_rank);
  }

  /**
   * @brief Prefix excess on [0, @p end_position): +1 for '1', −1 for '0'.
   */
  inline int excess(const size_t& end_position) const {
    return int64_t(rank1(end_position)) * 2 - int64_t(end_position);
  }

  /**
   * @brief Forward search: first position p ≥ @p start_position where
   * excess(p) = excess(@p start_position) + @p delta.
   * @details Scans remainder of current leaf, then descends using precomputed
   * bounds. Returns npos if no such position exists.
   */
  size_t fwdsearch(const size_t& start_position, const int& delta) const {
    if (start_position >= num_bits) {
      return npos;
    }

    // 1) scan the remainder of the current leaf
    const size_t leaf_block_index = block_of(start_position);
    const size_t block_begin = leaf_block_index * block_bits;
    const size_t block_end = std::min(num_bits, block_begin + block_bits);
    int leaf_delta = 0;
    const size_t leaf_result = leaf_fwd_bp_simd(
        leaf_block_index, block_begin, start_position, delta, leaf_delta);
    if (leaf_result != npos) {
      return leaf_result;
    }

    int remaining_delta = delta - leaf_delta;
    if (leaf_block_index + 1 < leaf_count) {
      size_t nodes_buffer[64];
      const size_t node_count = cover_blocks_collect(
          leaf_block_index + 1, leaf_count - 1, nodes_buffer);
      size_t segment_base = (leaf_block_index + 1) * block_bits;
      for (size_t j = 0; j < node_count; ++j) {
        const size_t node_index = nodes_buffer[j];
        if (remaining_delta == 0) {
          return segment_base;
        }
        if (node_min_prefix_excess[node_index] <= remaining_delta &&
            remaining_delta <= node_max_prefix_excess[node_index]) {
          return descend_fwd(node_index, remaining_delta, segment_base);
        }
        remaining_delta -= node_total_excess[node_index];
        segment_base += segment_size_bits[node_index];
      }
    }
    return npos;
  }

  /**
   * @brief Backward search: last position p ≤ @p start_position where
   * excess(p) = excess(@p start_position) + @p delta.
   * @details Scans inside the leaf to the left, then climbs to examine left
   * siblings. Returns npos if no such position exists.
   */
  size_t bwdsearch(const size_t& start_position, const int& delta) const {
    if (start_position > num_bits || start_position == 0) {
      return npos;
    }

    // 1) scan inside the block
    const size_t leaf_block_index = block_of(start_position - 1);
    const size_t block_begin = leaf_block_index * block_bits;
    int leaf_delta = 0;  // excess(start_position) - excess(block_begin)
    const size_t leaf_result = leaf_bwd_bp_simd(
        leaf_block_index, block_begin, start_position, delta, leaf_delta);
    if (leaf_result != npos) {
      return leaf_result;
    }

    // need = target - excess(block_begin) = excess(start_position) + delta -
    // excess(block_begin) = leaf_delta + delta
    int remaining_delta = leaf_delta + delta;
    size_t node_index = leaf_index_of(block_begin);
    size_t segment_base = block_begin;
    while (node_index > 1) {
      if (node_index & 1) {  // node_index is the right child
        const size_t sibling_index = node_index ^ 1;  // left sibling
        const size_t sibling_border =
            segment_base;  // right border of the sibling (== start(node_index))
        const int needed_inside_sibling =
            remaining_delta +
            node_total_excess[sibling_index];  // target in coordinates relative
                                               // to the start of sibling
        const bool allow_right_border =
            (sibling_border != start_position);  // j must be < start_position

        // try inside the sibling, but return only if a position is found
        if (needed_inside_sibling == 0 ||
            (node_min_prefix_excess[sibling_index] <= needed_inside_sibling &&
             needed_inside_sibling <= node_max_prefix_excess[sibling_index])) {
          const size_t result = descend_bwd(
              sibling_index, sibling_border - segment_size_bits[sibling_index],
              needed_inside_sibling, sibling_border, allow_right_border);
          if (result != npos) {
            return result;
          }
        }
        // junction between children is a separate branch (allowed only if < i)
        if (needed_inside_sibling == node_total_excess[sibling_index] &&
            sibling_border < start_position) {
          return sibling_border;
        }

        // stepped over the sibling, shifted the zero point of the coordinates
        remaining_delta += node_total_excess[sibling_index];
        segment_base -= segment_size_bits[sibling_index];
      }
      node_index >>= 1;
    }
    return npos;
  }

  /**
   * @brief Position of the first minimum of excess on [@p range_begin, @p
   * range_end]
   * (inclusive).
   * @return Position of first occurrence of minimum, or npos on invalid range.
   */
  size_t range_min_query_pos(const size_t& range_begin,
                             const size_t& range_end) const {
    if (range_begin > range_end || range_end >= num_bits) {
      return npos;
    }

    const size_t begin_block_index = block_of(range_begin);
    const size_t begin_block_start = begin_block_index * block_bits;
    const size_t begin_block_end =
        std::min(num_bits, begin_block_start + block_bits);
    const size_t end_block_index = block_of(range_end);
    const size_t end_block_start = end_block_index * block_bits;

    int best_value = INT_MAX;
    size_t best_position = npos;
    size_t chosen_node_index = 0;
    int prefix_excess = 0, prefix_at_choice = 0;

    // prefix
    int min_prefix_first_chunk = INT_MAX;
    size_t first_chunk_position = npos;
    const size_t end_of_first_chunk = std::min(
        range_end, (size_t)(begin_block_end ? begin_block_end - 1 : 0));
    if (range_begin <= end_of_first_chunk) {
      first_min_value_pos8(range_begin, end_of_first_chunk,
                           min_prefix_first_chunk, first_chunk_position);
      prefix_excess =
          (int64_t)rank1_in_block(range_begin, end_of_first_chunk + 1) * 2 -
          int64_t(end_of_first_chunk + 1 - range_begin);
      best_value = min_prefix_first_chunk;
      best_position = first_chunk_position;
      chosen_node_index = 0;
    }

    // middle
    if (begin_block_index + 1 <= end_block_index - 1) {
      size_t left_index = first_leaf_index + (begin_block_index + 1);
      size_t right_index = first_leaf_index + (end_block_index - 1);
      size_t right_nodes[64];
      int right_nodes_count = 0;

      while (left_index <= right_index) {
        if (left_index & 1) {
          const size_t node_index = left_index++;
          const int candidate =
              prefix_excess + node_min_prefix_excess[node_index];
          if (candidate < best_value) {
            best_value = candidate;
            best_position = npos;
            chosen_node_index = node_index;
            prefix_at_choice = prefix_excess;
          }
          prefix_excess += node_total_excess[node_index];
        }
        if ((right_index & 1) == 0) {
          right_nodes[right_nodes_count++] = right_index--;
        }
        left_index >>= 1;
        right_index >>= 1;
      }
      while (right_nodes_count--) {
        const size_t node_index = right_nodes[right_nodes_count];
        const int candidate =
            prefix_excess + node_min_prefix_excess[node_index];
        if (candidate < best_value) {
          best_value = candidate;
          best_position = npos;
          chosen_node_index = node_index;
          prefix_at_choice = prefix_excess;
        }
        prefix_excess += node_total_excess[node_index];
      }
    }

    // tail
    if (end_block_index != begin_block_index) {
      int min_prefix_last_chunk;
      size_t last_chunk_position;
      first_min_value_pos8(end_block_start, range_end, min_prefix_last_chunk,
                           last_chunk_position);
      const int candidate = prefix_excess + min_prefix_last_chunk;
      if (candidate < best_value) {
        best_value = candidate;
        best_position = last_chunk_position;
        chosen_node_index = 0;
      }
    }

    if (best_position != npos) {
      return best_position;
    }

    return descend_first_min(chosen_node_index, best_value - prefix_at_choice,
                             node_base(chosen_node_index));
  }

  /**
   * @brief Value of the minimum prefix excess on [@p range_begin, @p range_end]
   * relative to @p range_begin.
   * @details Equivalent to min_{t in [@p range_begin..@p range_end]}
   * (excess(t+1) - excess(@p range_begin)).
   */
  int range_min_query_val(const size_t& range_begin,
                          const size_t& range_end) const {
    if (range_begin > range_end || range_end >= num_bits) {
      return 0;
    }
    size_t min_position = range_min_query_pos(range_begin, range_end);
    if (min_position == npos) {
      return 0;
    }
    return excess(min_position + 1) - excess(range_begin);
  }

  /**
   * @brief Position of the first maximum of excess on [@p range_begin, @p
   * range_end]
   * (inclusive).
   * @return Position of first occurrence of maximum, or npos on invalid range.
   */
  size_t range_max_query_pos(const size_t& range_begin,
                             const size_t& range_end) const {
    if (range_begin > range_end || range_end >= num_bits) {
      return npos;
    }

    const size_t begin_block_index = block_of(range_begin);
    const size_t begin_block_start = begin_block_index * block_bits;
    const size_t begin_block_end =
        std::min(num_bits, begin_block_start + block_bits);
    const size_t end_block_index = block_of(range_end);
    const size_t end_block_start = end_block_index * block_bits;

    int best_value = INT_MIN;
    size_t best_position = npos;
    size_t chosen_node_index = 0;
    int prefix_excess = 0, prefix_at_choice = 0;

    // prefix
    int max_prefix_first_chunk = INT_MIN;
    size_t first_chunk_position = npos;
    const size_t end_of_first_chunk = std::min(
        range_end, (size_t)(begin_block_end ? begin_block_end - 1 : 0));
    if (range_begin <= end_of_first_chunk) {
      first_max_value_pos8(range_begin, end_of_first_chunk,
                           max_prefix_first_chunk, first_chunk_position);
      prefix_excess =
          (int64_t)rank1_in_block(range_begin, end_of_first_chunk + 1) * 2 -
          int64_t(end_of_first_chunk + 1 - range_begin);
      best_value = max_prefix_first_chunk;
      best_position = first_chunk_position;
      chosen_node_index = 0;
    }

    // middle
    if (begin_block_index + 1 <= end_block_index - 1) {
      size_t left_index = first_leaf_index + (begin_block_index + 1);
      size_t right_index = first_leaf_index + (end_block_index - 1);
      size_t right_nodes[64];
      int right_nodes_count = 0;

      while (left_index <= right_index) {
        if (left_index & 1) {
          const size_t node_index = left_index++;
          const int candidate =
              prefix_excess + node_max_prefix_excess[node_index];
          if (candidate > best_value) {
            best_value = candidate;
            best_position = npos;
            chosen_node_index = node_index;
            prefix_at_choice = prefix_excess;
          }
          prefix_excess += node_total_excess[node_index];
        }
        if ((right_index & 1) == 0) {
          right_nodes[right_nodes_count++] = right_index--;
        }
        left_index >>= 1;
        right_index >>= 1;
      }
      while (right_nodes_count--) {
        const size_t node_index = right_nodes[right_nodes_count];
        const int candidate =
            prefix_excess + node_max_prefix_excess[node_index];
        if (candidate > best_value) {
          best_value = candidate;
          best_position = npos;
          chosen_node_index = node_index;
          prefix_at_choice = prefix_excess;
        }
        prefix_excess += node_total_excess[node_index];
      }
    }

    // tail
    if (end_block_index != begin_block_index) {
      int max_prefix_last_chunk;
      size_t last_chunk_position;
      first_max_value_pos8(end_block_start, range_end, max_prefix_last_chunk,
                           last_chunk_position);
      const int candidate = prefix_excess + max_prefix_last_chunk;
      if (candidate > best_value) {
        best_value = candidate;
        best_position = last_chunk_position;
        chosen_node_index = 0;
      }
    }

    if (best_position != npos) {
      return best_position;
    }

    return descend_first_max(chosen_node_index, best_value - prefix_at_choice,
                             node_base(chosen_node_index));
  }

  /**
   * @brief Value of the maximum prefix excess on [@p range_begin, @p
   * range_end] relative to @p range_begin.
   */
  int range_max_query_val(const size_t& range_begin,
                          const size_t& range_end) const {
    if (range_begin > range_end || range_end >= num_bits) {
      return 0;
    }
    size_t max_position = range_max_query_pos(range_begin, range_end);
    if (max_position == npos) {
      return 0;
    }
    return excess(max_position + 1) - excess(range_begin);
  }

  /**
   * @brief How many times the minimum prefix excess occurs on [@p range_begin,
   * @p range_end].
   */
  size_t mincount(const size_t& range_begin, const size_t& range_end) const {
    if (range_begin > range_end || range_end >= num_bits) {
      return 0;
    }

    const size_t begin_block_index = block_of(range_begin);
    const size_t begin_block_start = begin_block_index * block_bits;
    const size_t begin_block_end =
        std::min(num_bits, begin_block_start + block_bits);
    const size_t end_block_index = block_of(range_end);
    const size_t end_block_start = end_block_index * block_bits;

    int best_value = INT_MAX;
    size_t min_count = 0;
    int prefix_excess = 0;

    // first chunk
    {
      int current_excess = 0, min_value = INT_MAX, local_count = 0;
      const size_t end_of_first_chunk =
          std::min(range_end, begin_block_end - 1);
      for (size_t position = range_begin; position <= end_of_first_chunk;
           ++position) {
        current_excess += bit(position) ? +1 : -1;
        if (current_excess < min_value) {
          min_value = current_excess;
          local_count = 1;
        } else if (current_excess == min_value) {
          ++local_count;
        }
      }
      best_value = min_value;
      min_count = local_count;
      prefix_excess = current_excess;  // offset toward the middle
    }

    // middle
    if (begin_block_index + 1 <= end_block_index - 1) {
      const auto middle_nodes =
          cover_blocks(begin_block_index + 1, end_block_index - 1);
      for (const size_t& node_index : middle_nodes) {
        const int candidate =
            prefix_excess + node_min_prefix_excess[node_index];
        if (candidate < best_value) {
          best_value = candidate;
          min_count = node_min_count[node_index];
        } else if (candidate == best_value) {
          min_count += node_min_count[node_index];
        }
        prefix_excess += node_total_excess[node_index];
      }
    }

    // last chunk
    if (end_block_index != begin_block_index) {
      int current_excess = 0, min_value = INT_MAX, local_count = 0;
      for (size_t position = end_block_start; position <= range_end;
           ++position) {
        current_excess += bit(position) ? +1 : -1;
        if (current_excess < min_value) {
          min_value = current_excess;
          local_count = 1;
        } else if (current_excess == min_value) {
          ++local_count;
        }
      }
      const int candidate = prefix_excess + min_value;
      if (candidate < best_value) {
        best_value = candidate;
        min_count = local_count;
      } else if (candidate == best_value) {
        min_count += local_count;
      }
    }
    return min_count;
  }

  /**
   * @brief Position of the @p target_min_rank-th (1-based) occurrence of the
   * minimum on [@p range_begin, @p range_end].
   * @return Position or npos if @p target_min_rank exceeds the number of
   * minima.
   */
  size_t minselect(const size_t& range_begin,
                   const size_t& range_end,
                   size_t target_min_rank) const {
    if (range_begin > range_end || range_end >= num_bits ||
        target_min_rank == 0) {
      return npos;
    }

    const size_t begin_block_index = block_of(range_begin);
    const size_t begin_block_start = begin_block_index * block_bits;
    const size_t begin_block_end =
        std::min(num_bits, begin_block_start + block_bits);
    const size_t end_block_index = block_of(range_end);
    const size_t end_block_start = end_block_index * block_bits;

    // prefix
    const size_t end_of_first_chunk = std::min(range_end, begin_block_end - 1);
    int current_first_chunk_excess = 0, min_first_chunk = 0;
    uint32_t count_first_chunk = 0;

    if (range_begin <= end_of_first_chunk) {
      scan_range_min_count8(range_begin, end_of_first_chunk,
                            current_first_chunk_excess, min_first_chunk,
                            count_first_chunk);
    } else {
      current_first_chunk_excess = 0;
      min_first_chunk = INT_MAX;
      count_first_chunk = 0;
    }

    int best_value = (min_first_chunk == INT_MAX ? INT_MAX : min_first_chunk);
    size_t total_count =
        (min_first_chunk == INT_MAX ? 0u : (size_t)count_first_chunk);
    int prefix_excess = current_first_chunk_excess;  // offset for middle

    size_t left_index = first_leaf_index + begin_block_index + 1;
    size_t right_index = first_leaf_index + end_block_index - 1;
    size_t right_nodes[64];
    int right_nodes_count = 0;

    // middle
    if (begin_block_index + 1 <= end_block_index - 1) {
      while (left_index <= right_index) {
        if (left_index & 1) {
          const int candidate =
              prefix_excess + node_min_prefix_excess[left_index];
          if (candidate < best_value) {
            best_value = candidate;
            total_count = node_min_count[left_index];
          } else if (candidate == best_value) {
            total_count += node_min_count[left_index];
          }
          prefix_excess += node_total_excess[left_index++];
        }
        if ((right_index & 1) == 0) {
          right_nodes[right_nodes_count++] = right_index--;
        }
        left_index >>= 1;
        right_index >>= 1;
      }
      while (right_nodes_count--) {
        const size_t node_index = right_nodes[right_nodes_count];
        const int candidate =
            prefix_excess + node_min_prefix_excess[node_index];
        if (candidate < best_value) {
          best_value = candidate;
          total_count = node_min_count[node_index];
        } else if (candidate == best_value) {
          total_count += node_min_count[node_index];
        }
        prefix_excess += node_total_excess[node_index];
      }
    }

    // tail
    int current_last_chunk_excess = 0, min_last_chunk = INT_MAX;
    uint32_t count_last_chunk = 0;
    if (end_block_index != begin_block_index) {
      scan_range_min_count8(end_block_start, range_end,
                            current_last_chunk_excess, min_last_chunk,
                            count_last_chunk);
      const int candidate = prefix_excess + min_last_chunk;
      if (candidate < best_value) {
        best_value = candidate;
        total_count = count_last_chunk;
      } else if (candidate == best_value) {
        total_count += count_last_chunk;
      }
    }

    if (target_min_rank > total_count) {
      return npos;
    }

    // prefix
    if (min_first_chunk == best_value && count_first_chunk) {
      if (target_min_rank <= count_first_chunk) {
        return qth_min_in_block(range_begin, end_of_first_chunk,
                                target_min_rank);
      }
      target_min_rank -= count_first_chunk;
    }

    // middle
    prefix_excess = current_first_chunk_excess;
    if (begin_block_index + 1 <= end_block_index - 1) {
      left_index = first_leaf_index + (begin_block_index + 1);
      right_index = first_leaf_index + (end_block_index - 1);
      right_nodes_count = 0;
      while (left_index <= right_index) {
        if (left_index & 1) {
          const size_t node_index = left_index++;
          const int candidate =
              prefix_excess + node_min_prefix_excess[node_index];
          if (candidate == best_value) {
            if (target_min_rank <= node_min_count[node_index]) {
              return descend_qth_min(node_index, best_value - prefix_excess,
                                     target_min_rank, node_base(node_index));
            }
            target_min_rank -= node_min_count[node_index];
          }
          prefix_excess += node_total_excess[node_index];
        }
        if (!(right_index & 1)) {
          right_nodes[right_nodes_count++] = right_index--;
        }
        left_index >>= 1;
        right_index >>= 1;
      }
      while (right_nodes_count--) {
        const size_t node_index = right_nodes[right_nodes_count];
        const int candidate =
            prefix_excess + node_min_prefix_excess[node_index];
        if (candidate == best_value) {
          if (target_min_rank <= node_min_count[node_index]) {
            return descend_qth_min(node_index, best_value - prefix_excess,
                                   target_min_rank, node_base(node_index));
          }
          target_min_rank -= node_min_count[node_index];
        }
        prefix_excess += node_total_excess[node_index];
      }
    }

    // tail
    if (end_block_index != begin_block_index &&
        (prefix_excess + min_last_chunk) == best_value) {
      return qth_min_in_block(end_block_start, range_end, target_min_rank);
    }

    return npos;
  }

  // ----- parentheses navigation (BP) -----

  /**
   * @brief close(@p open_position): matching ')' for '(' at @p open_position.
   * @return Position of matching ')', or npos.
   */
  inline size_t close(const size_t& open_position) const {
    if (open_position >= num_bits) {
      return npos;
    }
    return fwdsearch(open_position, -1);
  }

  /**
   * @brief open(@p close_position): matching '(' for ')' at @p close_position.
   * @return Position of matching '(', or npos.
   */
  inline size_t open(const size_t& close_position) const {
    // bwdsearch allows i in [1..num_bits]
    if (close_position == 0 || close_position > num_bits) {
      return npos;
    }
    const size_t result = bwdsearch(close_position, 0);
    return (result == npos ? npos : result + 1);
  }

  /**
   * @brief enclose(@p position): opening '(' that strictly encloses @p
   * position.
   * @return Position of enclosing '(', or npos.
   */
  inline size_t enclose(const size_t& position) const {
    if (position == 0 || position > num_bits) {
      return npos;
    }
    const size_t result = bwdsearch(position, -2);
    return (result == npos ? npos : result + 1);
  }

 private:
  /**
   * @brief Count "10" occurrences inside a 64-bit slice of given logical
   * length @p length.
   * @details Only positions fully inside the slice @p slice are counted.
   */
  static inline size_t pop10_in_slice64(const std::uint64_t& slice,
                                        const int& length) noexcept {
    if (length <= 1) {
      return 0;
    }
    std::uint64_t pattern_mask = slice & ~(slice >> 1);  // candidates for "10"
    if (length < 64) {
      pattern_mask &= ((std::uint64_t(1) << (length - 1)) - 1);
    } else {
      pattern_mask &= 0x7FFFFFFFFFFFFFFFull;
    }
    return (size_t)std::popcount(pattern_mask);
  }

  /**
   * @brief Rank of ones within [@p block_begin, @p block_end).
   * @details Works on word boundaries; @p block_end may equal @p block_begin.
   */
  size_t rank1_in_block(const size_t& block_begin,
                        const size_t& block_end) const noexcept {
    if (block_end <= block_begin) {
      return 0;
    }
    size_t left_word_index = block_begin >> 6;
    const size_t right_word_index = block_end >> 6;
    size_t left_offset = block_begin & 63;
    const size_t right_offset = block_end & 63;
    size_t count = 0;
    if (left_word_index == right_word_index) {
      const std::uint64_t mask =
          ((right_offset == 0) ? 0 : ((std::uint64_t(1) << right_offset) - 1)) &
          (~std::uint64_t(0) << left_offset);
      return (size_t)std::popcount(bits[left_word_index] & mask);
    }
    if (left_offset) {
      count += (size_t)std::popcount(bits[left_word_index] &
                                     (~std::uint64_t(0) << left_offset));
      ++left_word_index;
    }
    while (left_word_index < right_word_index) {
      count += (size_t)std::popcount(bits[left_word_index]);
      ++left_word_index;
    }
    if (right_offset) {
      count += (size_t)std::popcount(bits[right_word_index] &
                                     ((std::uint64_t(1) << right_offset) - 1));
    }
    return count;
  }

  /**
   * @brief Count "10" starts within [@p block_begin, @p block_end).
   * @details Accounts for cross-word boundaries.
   */
  size_t rr_in_block(const size_t& block_begin,
                     const size_t& block_end) const noexcept {
    if (block_end <= block_begin + 1) {
      return 0;
    }
    size_t left_word_index = block_begin >> 6;
    const size_t right_word_index = (block_end - 1) >> 6;
    const int left_offset = block_begin & 63;
    const int right_offset = (block_end - 1) & 63;
    size_t count = 0;

    if (left_word_index == right_word_index) {
      const int length = right_offset - left_offset + 1;
      const std::uint64_t slice = bits[left_word_index] >> left_offset;
      return pop10_in_slice64(slice, length);
    }

    // prefix word
    {
      const int length = 64 - left_offset;
      const std::uint64_t slice = bits[left_word_index] >> left_offset;
      count += pop10_in_slice64(slice, length);
    }
    // full interior words
    for (size_t word_index = left_word_index + 1; word_index < right_word_index;
         ++word_index) {
      const std::uint64_t word = bits[word_index];
      count += pop10_in_slice64(word, 64);
    }
    // suffix word
    {
      const int length = right_offset + 1;
      const std::uint64_t mask = (length == 64)
                                     ? ~std::uint64_t(0)
                                     : ((std::uint64_t(1) << length) - 1);
      const std::uint64_t slice = bits[right_word_index] & mask;
      count += pop10_in_slice64(slice, length);
    }
    // cross-word boundaries (bit 63 of w and bit 0 of w+1)
    for (size_t word_index = left_word_index; word_index < right_word_index;
         ++word_index) {
      if (((bits[word_index] >> 63) & 1u) &&
          ((bits[word_index + 1] & 1u) == 0)) {
        ++count;
      }
    }
    return count;
  }

  /**
   * @brief 1-based select of @p target_pattern_rank-th "10" within [@p
   * block_begin, @p block_end).
   * @return Position or npos.
   */
  size_t select10_in_block(const size_t& block_begin,
                           const size_t& block_end,
                           size_t target_pattern_rank) const noexcept {
    if (block_end <= block_begin + 1) {
      return npos;
    }
    size_t left_word_index = block_begin >> 6;
    const size_t right_word_index = (block_end - 1) >> 6;
    const int left_offset = block_begin & 63;
    const int right_offset = (block_end - 1) & 63;

    const auto select_in_masked_slice =
        [&](const std::uint64_t& slice, const int& length,
            const size_t& target_index) noexcept -> int {
      if (length <= 1) {
        return -1;
      }
      std::uint64_t pattern_mask = slice & ~(slice >> 1);
      if (length < 64) {
        pattern_mask &= ((std::uint64_t(1) << (length - 1)) - 1);
      } else {
        pattern_mask &= 0x7FFFFFFFFFFFFFFFull;
      }
      return select_in_word(pattern_mask, target_index);
    };

    if (left_word_index == right_word_index) {
      const int length = right_offset - left_offset + 1;
      const std::uint64_t slice = bits[left_word_index] >> left_offset;
      const int offset =
          select_in_masked_slice(slice, length, target_pattern_rank);
      return offset >= 0 ? (block_begin + (size_t)offset) : npos;
    }

    // prefix word
    {
      const int length = 64 - left_offset;
      const std::uint64_t slice = bits[left_word_index] >> left_offset;
      std::uint64_t pattern_mask = slice & ~(slice >> 1);
      pattern_mask &= ((std::uint64_t(1) << (length - 1)) - 1);
      const int count = std::popcount(pattern_mask);
      if (target_pattern_rank <= (size_t)count) {
        const int offset =
            select_in_masked_slice(slice, length, target_pattern_rank);
        return block_begin + (size_t)offset;
      }
      target_pattern_rank -= count;
    }

    // walk interior boundaries and words
    for (size_t word_index = left_word_index; word_index + 1 < right_word_index;
         ++word_index) {
      // boundary between w and w+1
      if (((bits[word_index] >> 63) & 1u) &&
          ((bits[word_index + 1] & 1u) == 0)) {
        if (--target_pattern_rank == 0) {
          return (word_index << 6) + 63;
        }
      }
      // full word w+1 (positions 0..62)
      const std::uint64_t next_word = bits[word_index + 1];
      const std::uint64_t pattern_mask =
          (next_word & ~(next_word >> 1)) & 0x7FFFFFFFFFFFFFFFull;
      const int count = std::popcount(pattern_mask);
      if (target_pattern_rank <= (size_t)count) {
        const int offset = select_in_word(pattern_mask, target_pattern_rank);
        if (offset == -1) {
          return npos;
        }
        return ((word_index + 1) << 6) + (size_t)offset;
      }
      target_pattern_rank -= count;
    }

    // boundary (w_r-1, w_r)
    if (((bits[right_word_index - 1] >> 63) & 1u) &&
        ((bits[right_word_index] & 1u) == 0)) {
      if (--target_pattern_rank == 0) {
        return ((right_word_index - 1) << 6) + 63;
      }
    }

    // suffix word w_r: [0..off_r]
    {
      const int length = right_offset + 1;
      const std::uint64_t mask = (length == 64)
                                     ? ~std::uint64_t(0)
                                     : ((std::uint64_t(1) << length) - 1);
      const std::uint64_t slice = bits[right_word_index] & mask;
      const int offset =
          select_in_masked_slice(slice, length, target_pattern_rank);
      if (offset >= 0) {
        return (right_word_index << 6) + (size_t)offset;
      }
    }
    return npos;
  }

  struct ByteAgg {
    int8_t excess_total;  // total excess for the byte
    int8_t min_prefix;    // minimum prefix within the byte (from 0)
    int8_t max_prefix;    // maximum prefix within the byte (from 0)
    uint8_t min_count;  // number of positions attaining the minimum in the byte
    uint8_t pattern10_count;  // number of "10" patterns inside the byte
    uint8_t first_bit;        // first bit (LSB)
    uint8_t last_bit;         // last bit (MSB)
    uint8_t pos_first_min;    // pos of first minimum in this byte
    uint8_t pos_first_max;    // pos of first maximum in this byte
  };

  struct LUT8Tables {
    std::array<ByteAgg, 256> agg;
    /**
     * @brief For each byte and each delta ∈ [-8..+8] — the first position in
     * the byte (0..7) where the prefix equals delta; -1 if none
     */
    std::array<std::array<int8_t, 17>, 256> fwd_pos;
    /**
     * @brief For each byte and each delta ∈ [-8..+8] — the last position in the
     * byte (0..7) where the prefix equals delta; -1 if none
     */
    std::array<std::array<int8_t, 17>, 256> bwd_pos;
  };

  /**
   * @brief Unified initialization of all byte-based LUTs.
   */
  static inline const LUT8Tables& LUT8_ALL() noexcept {
    static const LUT8Tables tables = [] {
      LUT8Tables lookup_tables{};
      for (int byte_value = 0; byte_value < 256; ++byte_value) {
        int current_excess = 0, min_prefix = INT_MAX, max_prefix = INT_MIN,
            min_count = 0, pattern10_count = 0;
        int first_min_position = 0, first_max_position = 0;
        int prefixes[8];
        const auto bit_at = [&](const int& bit_index) {
          return (byte_value >> bit_index) & 1;
        };  // LSB-first
        for (int bit_index = 0; bit_index < 8; ++bit_index) {
          int bit_value = bit_at(bit_index);
          if (bit_index + 1 < 8 && bit_value && bit_at(bit_index + 1) == 0) {
            ++pattern10_count;
          }
          current_excess += bit_value ? +1 : -1;
          prefixes[bit_index] = current_excess;
          if (current_excess < min_prefix) {
            min_prefix = current_excess;
            min_count = 1;
            first_min_position = bit_index;
          } else if (current_excess == min_prefix) {
            ++min_count;
          }
          if (current_excess > max_prefix) {
            max_prefix = current_excess;
            first_max_position = bit_index;
          }
        }
        ByteAgg aggregates{};
        aggregates.excess_total = current_excess;
        aggregates.min_prefix = (min_prefix == INT_MAX ? 0 : min_prefix);
        aggregates.max_prefix = (max_prefix == INT_MIN ? 0 : max_prefix);
        aggregates.min_count = min_count;
        aggregates.pattern10_count = pattern10_count;
        aggregates.first_bit = bit_at(0);
        aggregates.last_bit = bit_at(7);
        aggregates.pos_first_min = first_min_position;
        aggregates.pos_first_max = first_max_position;
        lookup_tables.agg[byte_value] = aggregates;
        auto& forward_positions = lookup_tables.fwd_pos[byte_value];
        auto& backward_positions = lookup_tables.bwd_pos[byte_value];
        forward_positions.fill(-1);
        backward_positions.fill(-1);
        for (int delta = -8; delta <= 8; ++delta) {
          for (int bit_index = 0; bit_index < 8; ++bit_index) {
            if (prefixes[bit_index] == delta) {
              forward_positions[delta + 8] = bit_index;
              break;
            }
          }
          for (int bit_index = 7; bit_index >= 0; --bit_index) {
            if (prefixes[bit_index] == delta) {
              backward_positions[delta + 8] = bit_index;
              break;
            }
          }
        }
      }
      return lookup_tables;
    }();
    return tables;
  }

  /**
   * @brief Byte-level aggregates.
   */
  static inline const std::array<ByteAgg, 256>& LUT8() noexcept {
    return LUT8_ALL().agg;
  }

  /**
   * @brief Position of the first occurrence of delta within the byte.
   */
  static inline const std::array<std::array<int8_t, 17>, 256>&
  LUT8_FWD_POS() noexcept {
    return LUT8_ALL().fwd_pos;
  }

  /**
   * @brief Position of the last occurrence of delta within the byte.
   */
  static inline const std::array<std::array<int8_t, 17>, 256>&
  LUT8_BWD_POS() noexcept {
    return LUT8_ALL().bwd_pos;
  }

  /**
   * @brief Range search: first index i in [@p begin, @p end_excl) with
   * @p arr[i] == @p target.
   */
  static inline size_t find_first_equal_i16_range(
      const int16_t* arr,
      size_t begin,
      size_t end_excl,
      const int16_t& target) noexcept {
    if (begin >= end_excl) {
      return npos;
    }
#if defined(PIXIE_AVX512_SUPPORT)
    static constexpr size_t STEP = 32;
    size_t i = begin;
    for (; i + STEP <= end_excl; i += STEP) {
      const size_t offset =
          ::find_forward_equal_i16_avx512(arr + i, target, npos);
      if (offset != npos) {
        return i + offset;
      }
    }
    for (; i < end_excl; ++i) {
      if (arr[i] == target) {
        return i;
      }
    }
#elif defined(PIXIE_AVX2_SUPPORT)
    static constexpr size_t STEP = 16;
    size_t i = begin;
    for (; i + STEP <= end_excl; i += STEP) {
      const size_t offset =
          ::find_forward_equal_i16_avx2(arr + i, target, npos);
      if (offset != npos) {
        return i + offset;
      }
    }
    for (; i < end_excl; ++i) {
      if (arr[i] == target) {
        return i;
      }
    }
#else
    for (size_t i = begin; i < end_excl; ++i) {
      if (arr[i] == target) {
        return i;
      }
    }
#endif
    return npos;
  }

  /**
   * @brief Range search: last index i in [@p begin, @p end_excl) with
   * @p arr[i] == @p target.
   */
  static inline size_t find_last_equal_i16_range(
      const int16_t* arr,
      size_t begin,
      size_t end_excl,
      const int16_t& target) noexcept {
    if (begin >= end_excl) {
      return npos;
    }
#if defined(PIXIE_AVX512_SUPPORT)
    static constexpr size_t STEP = 32;
    size_t block_start = end_excl;
    while (block_start >= begin + STEP) {
      block_start -= STEP;
      const size_t offset =
          ::find_backward_equal_i16_avx512(arr + block_start, target, npos);
      if (offset != npos) {
        return block_start + offset;
      }
    }
    for (size_t i = block_start; i > begin;) {
      --i;
      if (arr[i] == target) {
        return i;
      }
    }
#elif defined(PIXIE_AVX2_SUPPORT)
    static constexpr size_t STEP = 16;
    size_t block_start = end_excl;
    while (block_start >= begin + STEP) {
      block_start -= STEP;
      const size_t offset =
          ::find_backward_equal_i16_avx2(arr + block_start, target, npos);
      if (offset != npos) {
        return block_start + offset;
      }
    }
    for (size_t i = block_start; i > begin;) {
      --i;
      if (arr[i] == target) {
        return i;
      }
    }
#else
    for (size_t i = end_excl; i > begin;) {
      --i;
      if (arr[i] == target) {
        return i;
      }
    }
#endif
    return npos;
  }

  /**
   * @brief Forward search within a single leaf over the range [@p search_start,
   * @p search_end).
   * @param required_delta required relative excess (relative to @p
   * search_start).
   * @return The position, or npos.
   */
  inline size_t scan_leaf_fwd(const size_t& search_start,
                              const size_t& search_end,
                              const int& required_delta) const noexcept {
    if (search_start >= search_end) {
      return npos;
    }
    const auto& aggregates_table = LUT8();
    const auto& forward_lookup = LUT8_FWD_POS();
    int current_excess = 0;
    size_t position = search_start;
    while (position + 8 <= search_end) {
      const uint8_t byte_value = get_byte(position);
      const auto& byte_aggregate = aggregates_table[byte_value];
      const int local_need = required_delta - current_excess;
      if (local_need >= byte_aggregate.min_prefix &&
          local_need <= byte_aggregate.max_prefix && local_need >= -8 &&
          local_need <= 8) {
        const int8_t offset = forward_lookup[byte_value][local_need + 8];
        if (offset >= 0) {
          return position + size_t(offset);
        }
      }
      current_excess += byte_aggregate.excess_total;
      position += 8;
    }

    while (position < search_end) {
      current_excess += bit(position) ? 1 : -1;
      if (current_excess == required_delta) {
        return position;
      }
      ++position;
    }

    return npos;
  }

  /**
   * @brief Backward search within a single leaf over [@p block_begin, @p
   * block_end) (does not look to the right of @p block_end).
   * @param required_delta relative excess from @p block_begin.
   * @param allow_right_boundary whether @p block_end may be used as a valid
   * answer.
   * @param global_right_border global right limit (as in descend_bwd).
   */
  inline size_t scan_leaf_bwd(
      const size_t& block_begin,
      const size_t& block_end,
      const int& required_delta,
      const bool& allow_right_boundary,
      const size_t& global_right_border) const noexcept {
    if (block_begin >= block_end) {
      if ((block_begin < global_right_border || allow_right_boundary) &&
          required_delta == 0) {
        return block_begin;
      }
      return npos;
    }

    const int16_t* leaf_prefix_ptr =
        &leaf_prefix[block_of(block_begin) * leaf_prefix_stride];
    const size_t end_inclusive = block_end - block_begin;
    if (end_inclusive > 0) {
      const size_t position = find_last_equal_i16_range(
          leaf_prefix_ptr, 1, end_inclusive + 1, required_delta);
      if (position != npos) {
        return block_begin + position;
      }
    }

    if ((block_begin < global_right_border || allow_right_boundary) &&
        required_delta == 0) {
      return block_begin;
    }

    return npos;
  }

  /**
   * @brief Extract 8 bits starting at @p position (LSB-first across words).
   */
  inline uint8_t get_byte(const size_t& position) const noexcept {
    const size_t word_index = position >> 6;
    const size_t offset = position & 63;
    const std::uint64_t lower_word = bits[word_index] >> offset;
    if (offset <= 56) {
      return uint8_t(lower_word & 0xFFu);
    }
    const std::uint64_t higher_word =
        (word_index + 1 < bits.size()) ? bits[word_index + 1] : 0;
    const std::uint64_t byte_value =
        (lower_word | (higher_word << (64 - offset))) & 0xFFu;
    return uint8_t(byte_value);
  }

  /**
   * @brief Descend to the first (leftmost) maximum with node-relative prefix
   * equal to @p target_prefix.
   * @param node_index Node index.
   * @param target_prefix Target prefix within node coordinates.
   * @param segment_base Starting global position of node.
   * @return Position or npos.
   */
  size_t descend_first_max(size_t node_index,
                           int target_prefix,
                           size_t segment_base) const noexcept {
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1, right_child = left_child | 1;
      const int left_max = node_max_prefix_excess[left_child];
      const int right_max =
          node_total_excess[left_child] + node_max_prefix_excess[right_child];
      if (left_max >= right_max && left_max == target_prefix) {
        node_index = left_child;
      } else if (right_max == target_prefix) {
        segment_base += segment_size_bits[left_child];
        target_prefix -= node_total_excess[left_child];
        node_index = right_child;
      } else {
        return npos;
      }
    }

    const size_t segment_begin = segment_base;
    const size_t segment_end =
        std::min(segment_base + segment_size_bits[node_index], num_bits);
    int max_value;
    size_t position;

    first_max_value_pos8(segment_begin,
                         segment_end ? (segment_end - 1) : segment_begin,
                         max_value, position);
    return (max_value == target_prefix ? position : npos);
  }

  /**
   * @brief Heap index of the first leaf node.
   */
  size_t first_leaf_index = 1;

  /**
   * @brief Block index containing position @p position.
   */
  size_t block_of(const size_t& position) const noexcept {
    return position / block_bits;
  }

  /**
   * @brief Heap index of the leaf whose segment starts at @p block_start.
   */
  size_t leaf_index_of(const size_t& block_start) const noexcept {
    return first_leaf_index + block_of(block_start);
  }

  /**
   * @brief Starting global position for node @p node_index (0-indexed).
   * @details Walk up to compute base for internal nodes; direct for leaves.
   */
  size_t node_base(size_t node_index) const noexcept {
    if (node_index >= first_leaf_index) {
      return (node_index - first_leaf_index) * block_bits;
    }

    size_t base = 0;
    for (; node_index > 1; node_index >>= 1) {
      if (node_index & 1) {
        base += segment_size_bits[node_index - 1];
      }
    }
    return base;
  }

  /**
   * @brief Cover a range of whole blocks [@p block_begin_index..@p
   * block_end_index] (inclusive) with O(log n) maximal nodes.
   * @details Returns node indices in left-to-right order.
   */
  std::vector<size_t> cover_blocks(const size_t& block_begin_index,
                                   const size_t& block_end_index) const {
    size_t left_index = first_leaf_index + block_begin_index;
    size_t right_index = first_leaf_index + block_end_index;
    std::vector<size_t> left_nodes, right_nodes;
    while (left_index <= right_index) {
      if ((left_index & 1) == 1) {
        left_nodes.push_back(left_index++);
      }
      if ((right_index & 1) == 0) {
        right_nodes.push_back(right_index--);
      }
      left_index >>= 1;
      right_index >>= 1;
    }
    std::reverse(right_nodes.begin(), right_nodes.end());
    left_nodes.insert(left_nodes.end(), right_nodes.begin(), right_nodes.end());
    return left_nodes;
  }

  /**
   * @brief Descend for fwdsearch to find first position where relative prefix
   * equals @p need.
   */
  size_t descend_fwd(size_t node_index,
                     int required_delta,
                     size_t segment_base) const noexcept {
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      if (node_min_prefix_excess[left_child] <= required_delta &&
          required_delta <= node_max_prefix_excess[left_child]) {
        node_index = left_child;
      } else {
        required_delta -= node_total_excess[left_child];
        segment_base += segment_size_bits[left_child];
        node_index = right_child;
      }
    }
    return scan_leaf_fwd(
        segment_base,
        std::min(segment_base + segment_size_bits[node_index], num_bits),
        required_delta);
  }

  /**
   * @brief Descend for bwdsearch to return the rightmost solution.
   * @param node_index Current node.
   * @param segment_base Left border of node @p node_index.
   * @param required_delta Target relative prefix inside @p node_index.
   * @param global_right_border Global right limit (exclusive if
   * !@p allow_right_boundary).
   * @param allow_right_boundary Whether right border is allowed to match.
   */
  size_t descend_bwd(size_t node_index,
                     const size_t& segment_base,
                     const int& required_delta,
                     const size_t& global_right_border,
                     const bool& allow_right_boundary) const noexcept {
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      const int required_in_right =
          required_delta - node_total_excess[left_child];

      // 1) try the right child first (to capture the rightmost j)
      if (node_min_prefix_excess[right_child] <= required_in_right &&
          required_in_right <= node_max_prefix_excess[right_child]) {
        const size_t result = descend_bwd(
            right_child, segment_base + segment_size_bits[left_child],
            required_in_right, global_right_border, allow_right_boundary);
        if (result != npos) {
          return result;
        }
      }

      // 2) junction between children (end of the left child)
      const size_t junction = segment_base + segment_size_bits[left_child];
      if (required_delta == node_total_excess[left_child] &&
          (junction < global_right_border || allow_right_boundary)) {
        return junction;
      }

      // 3) can we move left within the range?
      if (node_min_prefix_excess[left_child] <= required_delta &&
          required_delta <= node_max_prefix_excess[left_child]) {
        node_index = left_child;
        continue;
      }

      // None of (1)-(3) worked. The only possible point is the left border of
      // the node.
      if (required_delta == 0 &&
          (segment_base < global_right_border || allow_right_boundary)) {
        return segment_base;
      }

      return npos;
    }

    return scan_leaf_bwd(
        segment_base,
        std::min(
            global_right_border,
            std::min(segment_base + segment_size_bits[node_index], num_bits)),
        required_delta, allow_right_boundary, global_right_border);
  }

  /**
   * @brief Descend to find first position where node-relative prefix equals
   * @p target_prefix (minimum).
   */
  size_t descend_first_min(size_t node_index,
                           int target_prefix,
                           size_t segment_base) const noexcept {
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1, right_child = left_child | 1;
      const int left_min = node_min_prefix_excess[left_child];
      const int right_min =
          node_total_excess[left_child] + node_min_prefix_excess[right_child];
      if (left_min <= right_min && left_min == target_prefix) {
        node_index = left_child;
      } else if (right_min == target_prefix) {
        segment_base += segment_size_bits[left_child];
        target_prefix -= node_total_excess[left_child];
        node_index = right_child;
      } else {
        return npos;
      }
    }

    const size_t segment_begin = segment_base;
    const size_t segment_end =
        std::min(segment_base + segment_size_bits[node_index], num_bits);
    int min_value;
    size_t position;

    first_min_value_pos8(segment_begin,
                         segment_end ? (segment_end - 1) : segment_begin,
                         min_value, position);
    return (min_value == target_prefix ? position : npos);
  }

  /**
   * @brief Descend to find the @p target_min_rank-th minimum (1-based) where
   * node-relative prefix equals @p target_prefix.
   */
  size_t descend_qth_min(size_t node_index,
                         int target_prefix,
                         size_t target_min_rank,
                         size_t segment_base) const noexcept {
    while (node_index < first_leaf_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      const int left_min = node_min_prefix_excess[left_child];
      const int right_min =
          node_total_excess[left_child] + node_min_prefix_excess[right_child];
      if (left_min == target_prefix) {
        if (node_min_count[left_child] >= target_min_rank) {
          node_index = left_child;
          continue;
        }
        target_min_rank -= node_min_count[left_child];
      }
      if (right_min == target_prefix) {
        segment_base += segment_size_bits[left_child];
        target_prefix -= node_total_excess[left_child];
        node_index = right_child;
        continue;
      }
      return npos;
    }
    return qth_min_in_block(
        segment_base,
        std::min(segment_base + segment_size_bits[node_index], num_bits) - 1,
        target_min_rank);
  }

  /**
   * @brief 1-based select of @p target_one_rank-th '1' within [@p block_begin,
   * @p block_end).
   */
  size_t select1_in_block(const size_t& block_begin,
                          const size_t& block_end,
                          size_t target_one_rank) const noexcept {
    size_t left_word_index = block_begin >> 6;
    const size_t right_word_index = (block_end >> 6);
    const size_t left_offset = block_begin & 63;
    const std::uint64_t left_mask =
        (left_offset ? (~std::uint64_t(0) << left_offset) : ~std::uint64_t(0));
    if (left_word_index == right_word_index) {
      const std::uint64_t word =
          bits[left_word_index] & left_mask &
          ((block_end & 63) ? ((std::uint64_t(1) << (block_end & 63)) - 1)
                            : ~std::uint64_t(0));
      return block_begin + select_in_word(word, target_one_rank);
    }
    // prefix
    if (left_offset) {
      const std::uint64_t word = bits[left_word_index] & left_mask;
      const int count = std::popcount(word);
      if (target_one_rank <= (size_t)count) {
        return block_begin + select_in_word(word, target_one_rank);
      }
      target_one_rank -= count;
      left_word_index++;
    }
    // full words
    while (left_word_index < right_word_index) {
      const std::uint64_t word = bits[left_word_index];
      const int count = std::popcount(word);
      if (target_one_rank <= (size_t)count) {
        return (left_word_index << 6) + select_in_word(word, target_one_rank);
      }
      target_one_rank -= count;
      ++left_word_index;
    }
    // tail
    const size_t right_offset = block_end & 63;
    if (right_offset) {
      const std::uint64_t word =
          bits[left_word_index] & ((std::uint64_t(1) << right_offset) - 1);
      const int count = std::popcount(word);
      if (target_one_rank <= (size_t)count) {
        return (left_word_index << 6) + select_in_word(word, target_one_rank);
      }
    }
    return npos;
  }

  /**
   * @brief 1-based select of @p target_zero_rank-th '0' within [@p
   * block_begin, @p block_end).
   */
  size_t select0_in_block(const size_t& block_begin,
                          const size_t& block_end,
                          size_t target_zero_rank) const noexcept {
    if (block_end <= block_begin) {
      return npos;
    }

    size_t left_word_index = block_begin >> 6;
    const size_t right_word_index = block_end >> 6;
    const size_t left_offset = block_begin & 63;

    if (left_word_index == right_word_index) {
      const std::uint64_t left_mask =
          (left_offset ? (~std::uint64_t(0) << left_offset)
                       : ~std::uint64_t(0));
      const std::uint64_t right_mask =
          ((block_end & 63) ? ((std::uint64_t(1) << (block_end & 63)) - 1)
                            : ~std::uint64_t(0));
      const std::uint64_t word =
          (~bits[left_word_index]) & left_mask & right_mask;
      const int offset = select_in_word(word, target_zero_rank);
      return (offset >= 0) ? (block_begin + (size_t)offset) : npos;
    }

    // prefix
    if (left_offset) {
      const std::uint64_t word =
          (~bits[left_word_index]) & (~std::uint64_t(0) << left_offset);
      const int count = std::popcount(word);
      if (target_zero_rank <= (size_t)count) {
        const int offset = select_in_word(word, target_zero_rank);
        return (offset >= 0) ? (block_begin + (size_t)offset) : npos;
      }
      target_zero_rank -= count;
      ++left_word_index;
    }

    // full words
    while (left_word_index < right_word_index) {
      const std::uint64_t word = ~bits[left_word_index];
      const int count = std::popcount(word);
      if (target_zero_rank <= (size_t)count) {
        const int offset = select_in_word(word, target_zero_rank);
        return (offset >= 0) ? ((left_word_index << 6) + (size_t)offset) : npos;
      }
      target_zero_rank -= count;
      ++left_word_index;
    }

    // tail
    const size_t right_offset = block_end & 63;
    if (right_offset) {
      const std::uint64_t word =
          (~bits[left_word_index]) & ((std::uint64_t(1) << right_offset) - 1);
      const int count = std::popcount(word);
      if (target_zero_rank <= (size_t)count) {
        const int offset = select_in_word(word, target_zero_rank);
        return (offset >= 0) ? ((left_word_index << 6) + (size_t)offset) : npos;
      }
    }
    return npos;
  }

  /**
   * @brief 1-based select of @p target_rank-th set bit inside a 64-bit word.
   * @return Bit index [0..63] or −1 if not found.
   */
  static inline int select_in_word(std::uint64_t word,
                                   size_t target_rank) noexcept {
    while (word) {
      if (--target_rank == 0) {
        return std::countr_zero(word);
      }
      word &= (word - 1);
    }
    return -1;
  }

  /**
   * @brief Ceil division of positive integers.
   */
  static inline size_t ceil_div(const size_t& numerator,
                                const size_t& denominator) noexcept {
    return (numerator + denominator - 1) / denominator;
  }

  /**
   * @brief Number of node slots needed for @p bit_count with leaf size @p
   * block_size_pow2.
   * @details Includes internal node space (heap layout) plus leaves.
   */
  static inline size_t nodeslots_for(const size_t& bit_count,
                                     const size_t& block_size_pow2) noexcept {
    if (bit_count == 0) {
      return 0;
    }
    size_t leaf_node_count = ceil_div(bit_count, block_size_pow2);
    return std::bit_ceil(std::max<size_t>(1, leaf_node_count)) +
           leaf_node_count;
  }

  /**
   * @brief Auxiliary overhead in bytes/bitvector bytes for given parameters.
   */
  static inline float overhead_for(const size_t& bit_count,
                                   const size_t& block_size_pow2) noexcept {
    static constexpr size_t AUX_SLOT_BYTES =
        sizeof(uint32_t) + sizeof(int32_t) + sizeof(int32_t) + sizeof(int32_t) +
        sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint8_t) +
        sizeof(uint8_t) + sizeof(int16_t);

    size_t bitvector_bytes = ceil_div(bit_count, 64) * 8;
    if (bitvector_bytes == 0) {
      return 0;
    }
    size_t slot_count = nodeslots_for(bit_count, block_size_pow2);
    size_t aux_bytes = slot_count * AUX_SLOT_BYTES;
    return ((float)aux_bytes) / ((float)bitvector_bytes);
  }

  /**
   * @brief Choose minimal block size (power of two) keeping overhead ≤ @p
   * overhead_cap.
   * @details Returns 64 if @p overhead_cap < 0 (no constraint). Clamped to
   * ≤16384 or @p bit_count.
   */
  static inline size_t choose_block_bits_for_overhead(
      const size_t& bit_count,
      const float& overhead_cap) noexcept {
    if (overhead_cap < 0.f) {
      return 64;
    }

    const size_t max_block_bits = std::min<size_t>(bit_count, 16384);
    size_t candidate_block_bits = 64;
    while (candidate_block_bits < max_block_bits) {
      if (overhead_for(bit_count, candidate_block_bits) <= overhead_cap) {
        break;
      }
      candidate_block_bits <<= 1;
    }
    return candidate_block_bits;
  }

  /**
   * @brief Build internal structures from a 0/1 string.
   * @param leaf_block_bits Desired leaf size (0 = auto).
   * @param max_overhead Overhead cap (<0 = no cap).
   */
  void build_from_string(const std::string& bit_string_input,
                         const size_t& leaf_block_bits = 0,
                         const float& max_overhead = -1.0) {
    num_bits = bit_string_input.size();
    bits.assign((num_bits + 63) / 64, 0);
    for (size_t i = 0; i < num_bits; ++i) {
      if (bit_string_input[i] == '1') {
        set1(i);
      }
    }
    build(leaf_block_bits, max_overhead);
  }

  /**
   * @brief Build internal structures from 64-bit words.
   * @param words Words with LSB-first bits.
   * @param bit_count Number of valid bits.
   * @param leaf_block_bits Desired leaf size (0 = auto).
   * @param max_overhead Overhead cap (<0 = no cap).
   */
  void build_from_words(const std::vector<std::uint64_t>& words,
                        const size_t& bit_count,
                        const size_t& leaf_block_bits = 0,
                        const float& max_overhead = -1.0) {
    bits = words;
    num_bits = bit_count;
    if (bits.size() * 64 < num_bits) {
      bits.resize((num_bits + 63) / 64);
    }
    build(leaf_block_bits, max_overhead);
  }

  /**
   * @brief Read bit at position @p position (LSB-first across words).
   */
  inline int bit(const size_t& position) const noexcept {
    return (bits[position >> 6] >> (position & 63)) & 1u;
  }

  /**
   * @brief Set bit at position @p position to 1.
   */
  inline void set1(const size_t& position) noexcept {
    bits[position >> 6] |= (std::uint64_t(1) << (position & 63));
  }

  /**
   * @brief Number of ones in node @p node_index computed from size and total
   * excess.
   */
  inline uint32_t ones_in_node(const size_t& node_index) const noexcept {
    return ((int64_t)segment_size_bits[node_index] +
            (int64_t)node_total_excess[node_index]) >>
           1;
  }

  /**
   * @brief Scan [@p range_begin..@p range_end] inclusive, computing minimum
   * value and how many
   * times it is attained.
   * @details Uses 8-bit LUT for speed; outputs @p current_excess at end, @p
   * min_value, and @p count.
   */
  inline void scan_range_min_count8(size_t range_begin,
                                    const size_t& range_end,
                                    int& current_excess,
                                    int& min_value,
                                    uint32_t& count) const noexcept {
    current_excess = 0;
    min_value = INT_MAX;
    count = 0;
    if (range_end < range_begin) {
      min_value = 0;
      return;
    }
    // to byte alignment
    while (range_begin <= range_end && (range_begin & 7)) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
        count = 1;
      } else if (current_excess == min_value) {
        ++count;
      }
      ++range_begin;
    }
    // full bytes
    const auto& aggregates_table = LUT8();
    while (range_begin + 7 <= range_end) {
      const auto& byte_aggregate = aggregates_table[get_byte(range_begin)];
      const int candidate = current_excess + byte_aggregate.min_prefix;
      if (candidate < min_value) {
        min_value = candidate;
        count = byte_aggregate.min_count;
      } else if (candidate == min_value) {
        count += byte_aggregate.min_count;
      }
      current_excess += byte_aggregate.excess_total;
      range_begin += 8;
    }
    // tail
    while (range_begin <= range_end) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
        count = 1;
      } else if (current_excess == min_value) {
        ++count;
      }
      ++range_begin;
    }
    if (min_value == INT_MAX) {
      min_value = count = 0;
    }
  }

  /**
   * @brief Cover the range of whole blocks [@p block_begin_index..@p
   * block_end_index] (inclusive) with tree nodes, writing them to
   * @p out_nodes from left to right.
   * @return The number of nodes written.
   */
  inline size_t cover_blocks_collect(const size_t& block_begin_index,
                                     const size_t& block_end_index,
                                     size_t (&out_nodes)[64]) const noexcept {
    if (leaf_count == 0 || block_begin_index > block_end_index) {
      return 0;
    }
    size_t left_index = first_leaf_index + block_begin_index;
    size_t right_index = first_leaf_index + block_end_index;
    size_t left_nodes[32];
    size_t right_nodes[32];
    size_t left_count = 0, right_count = 0;
    while (left_index <= right_index) {
      if (left_index & 1) {
        left_nodes[left_count++] = left_index++;
      }
      if ((right_index & 1) == 0) {
        right_nodes[right_count++] = right_index--;
      }
      left_index >>= 1;
      right_index >>= 1;
    }
    size_t out_count = 0;
    for (size_t i = 0; i < left_count; ++i) {
      out_nodes[out_count++] = left_nodes[i];
    }
    while (right_count > 0) {
      out_nodes[out_count++] = right_nodes[--right_count];
    }
    return out_count;
  }

  /**
   * @brief Select @p target_min_rank-th minimum (1-based) inside [@p
   * range_begin..@p range_end] inclusive in two 8-bit passes.
   * @details First pass finds global minimum, second selects the requested
   * position.
   */
  inline size_t qth_min_in_block(const size_t& range_begin,
                                 const size_t& range_end,
                                 size_t target_min_rank) const noexcept {
    if (range_end < range_begin || target_min_rank == 0) {
      return npos;
    }

    const auto& aggregates_table = LUT8();

    int current_excess = 0, min_value = INT_MAX;
    size_t position = range_begin;

    while (position <= range_end && (position & 7)) {
      current_excess += bit(position) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
      }
      ++position;
    }
    while (position + 7 <= range_end) {
      const auto& byte_aggregate = aggregates_table[get_byte(position)];
      min_value =
          std::min(min_value, current_excess + byte_aggregate.min_prefix);
      current_excess += byte_aggregate.excess_total;
      position += 8;
    }
    while (position <= range_end) {
      current_excess += bit(position) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
      }
      ++position;
    }

    current_excess = 0;
    position = range_begin;

    // to byte alignment
    while (position <= range_end && (position & 7)) {
      current_excess += bit(position) ? +1 : -1;
      if (current_excess == min_value) {
        if (--target_min_rank == 0) {
          return position;
        }
      }
      ++position;
    }

    // full bytes
    while (position + 7 <= range_end) {
      const uint8_t byte_value = get_byte(position);
      const auto& byte_aggregate = aggregates_table[byte_value];
      const int candidate = current_excess + byte_aggregate.min_prefix;
      if (candidate == min_value) {
        int prefix_sum = 0;
        for (int k = 0; k < 8; ++k) {
          prefix_sum += ((byte_value >> k) & 1u) ? +1 : -1;
          if (prefix_sum == byte_aggregate.min_prefix) {
            if (--target_min_rank == 0) {
              return position + k;
            }
          }
        }
      }
      current_excess += byte_aggregate.excess_total;
      position += 8;
    }

    // tail
    while (position <= range_end) {
      current_excess += bit(position) ? +1 : -1;
      if (current_excess == min_value) {
        if (--target_min_rank == 0) {
          return position;
        }
      }
      ++position;
    }

    return npos;
  }

  /**
   * @brief Performs a forward search within a single leaf.
   * @details Starting from position @p start_position (global index) inside the
   * leaf that begins at @p leaf_block_begin, it looks for the nearest position
   * where the excess changes by @p delta relative to @p start_position. The
   * function uses a SIMD-accelerated search over the leaf’s prefix-sum array.
   * If a matching position is found, returns its global index. If not found,
   * returns npos and sets @p leaf_delta to the net excess change from @p
   * start_position to the end of the leaf.
   */
  inline size_t leaf_fwd_bp_simd(const size_t& leaf_index,
                                 const size_t& leaf_block_begin,
                                 const size_t& start_position,
                                 const int& delta,
                                 int& leaf_delta) const noexcept {
    const size_t leaf_length = segment_size_bits[first_leaf_index + leaf_index];
    const int16_t* leaf_prefix_ptr =
        &leaf_prefix[leaf_index * leaf_prefix_stride];
    const size_t offset_in_leaf = start_position - leaf_block_begin;
    if (offset_in_leaf > leaf_length) {
      leaf_delta = 0;
      return npos;
    }
    const int16_t start_prefix = leaf_prefix_ptr[offset_in_leaf];
    const size_t match_index =
        find_first_equal_i16_range(leaf_prefix_ptr, offset_in_leaf + 1,
                                   leaf_length + 1, start_prefix + delta);
    if (match_index != npos) {
      return leaf_block_begin + match_index - 1;
    }
    leaf_delta = leaf_prefix_ptr[leaf_length] - start_prefix;
    return npos;
  }

  /**
   * @brief Performs a backward search within a leaf segment for a position
   * earlier than @p start_position such that the excess difference from @p
   * start_position equals @p delta.
   * @details Uses SIMD (AVX2) to scan the leaf’s prefix-sum array of excess
   * values. Returns the absolute bit position if found; otherwise returns npos
   * and outputs the leaf-local excess delta at @p start_position.
   */
  inline size_t leaf_bwd_bp_simd(const size_t& leaf_index,
                                 const size_t& leaf_block_begin,
                                 const size_t& start_position,
                                 const int& delta,
                                 int& leaf_delta) const noexcept {
    const int16_t* leaf_prefix_ptr =
        &leaf_prefix[leaf_index * leaf_prefix_stride];
    const size_t offset_in_leaf = start_position - leaf_block_begin;
    if (offset_in_leaf > segment_size_bits[first_leaf_index + leaf_index]) {
      leaf_delta = 0;
      return npos;
    }
    const int16_t start_prefix = leaf_prefix_ptr[offset_in_leaf];
    if (offset_in_leaf > 0) {
      const size_t match_index = find_last_equal_i16_range(
          leaf_prefix_ptr, 0, offset_in_leaf, start_prefix + delta);
      if (match_index != npos) {
        return leaf_block_begin + match_index;
      }
    }
    leaf_delta = start_prefix;
    return npos;
  }

  /**
   * @brief Find first minimum value and its first position in [@p range_begin
   * ..@p range_end] inclusive using 8-bit LUT.
   * @param min_value_out Output minimum value (0 if empty).
   * @param first_position Output position of first minimum (npos if none).
   */
  inline void first_min_value_pos8(size_t range_begin,
                                   const size_t& range_end,
                                   int& min_value_out,
                                   size_t& first_position) const noexcept {
    const auto& aggregates_table = LUT8();
    int current_excess = 0;
    int min_value = INT_MAX;
    first_position = npos;

    // to byte allignment
    while (range_begin <= range_end && (range_begin & 7)) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
        first_position = range_begin;
      }
      ++range_begin;
    }

    // full bytes
    while (range_begin + 7 <= range_end) {
      const auto& byte_aggregate = aggregates_table[get_byte(range_begin)];
      const int candidate = current_excess + byte_aggregate.min_prefix;
      if (candidate < min_value) {
        min_value = candidate;
        first_position = range_begin + byte_aggregate.pos_first_min;
      }
      current_excess += byte_aggregate.excess_total;
      range_begin += 8;
    }

    // tail
    while (range_begin <= range_end) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess < min_value) {
        min_value = current_excess;
        first_position = range_begin;
      }
      ++range_begin;
    }

    min_value_out = (min_value == INT_MAX ? 0 : min_value);
  }

  /**
   * @brief Find first maximum value and its first position in [@p range_begin
   * ..@p range_end] inclusive using 8-bit LUT.
   * @param max_value_out Output maximum value (0 if empty).
   * @param first_position Output position of first maximum (npos if none).
   */
  inline void first_max_value_pos8(size_t range_begin,
                                   const size_t& range_end,
                                   int& max_value_out,
                                   size_t& first_position) const noexcept {
    const auto& aggregates_table = LUT8();
    int current_excess = 0;
    int max_value = INT_MIN;
    first_position = npos;

    while (range_begin <= range_end && (range_begin & 7)) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess > max_value) {
        max_value = current_excess;
        first_position = range_begin;
      }
      ++range_begin;
    }

    while (range_begin + 7 <= range_end) {
      const auto& byte_aggregate = aggregates_table[get_byte(range_begin)];
      const int candidate = current_excess + byte_aggregate.max_prefix;
      if (candidate > max_value) {
        max_value = candidate;
        first_position = range_begin + byte_aggregate.pos_first_max;
      }
      current_excess += byte_aggregate.excess_total;
      range_begin += 8;
    }

    while (range_begin <= range_end) {
      current_excess += bit(range_begin) ? +1 : -1;
      if (current_excess > max_value) {
        max_value = current_excess;
        first_position = range_begin;
      }
      ++range_begin;
    }

    max_value_out = (max_value == INT_MIN ? 0 : max_value);
  }

  /**
   * @brief Build the tree arrays and per-node aggregates.
   * @details Chooses block_bits honoring @p max_overhead or explicit @p
   * leaf_block_bits, allocates arrays, fills leaves via LUT, and builds
   * internal nodes bottom-up.
   */
  void build(const size_t& leaf_block_bits, const float& max_overhead) {
    // the lower clamp depends on the desired overhead fraction; otherwise use
    // 64
    const size_t clamp_by_overhead =
        (max_overhead >= 0.0
             ? choose_block_bits_for_overhead(num_bits, max_overhead)
             : size_t(64));

    // chosen block_bits: honor an explicit request, but not below
    // clamp_by_overhead
    if (leaf_block_bits == 0) {
      block_bits =
          std::max(clamp_by_overhead,
                   std::bit_ceil<size_t>(
                       (num_bits <= 1) ? 1 : std::bit_width(num_bits - 1)));
    } else {
      block_bits =
          std::max(clamp_by_overhead,
                   std::bit_ceil(std::max<size_t>(1, leaf_block_bits)));
    }

#ifdef DEBUG
    // finalizes the achieved overhead percentage
    built_overhead = overhead_for(num_bits, block_bits);
#endif

    leaf_count = ceil_div(num_bits, block_bits);
    first_leaf_index = std::bit_ceil(std::max<size_t>(1, leaf_count));
    const size_t tree_size = first_leaf_index + leaf_count - 1;
    segment_size_bits.assign(tree_size + 1, 0);
    node_total_excess.assign(tree_size + 1, 0);
    node_min_prefix_excess.assign(tree_size + 1, 0);
    node_max_prefix_excess.assign(tree_size + 1, 0);
    node_min_count.assign(tree_size + 1, 0);
    node_pattern10_count.assign(tree_size + 1, 0);
    node_first_bit.assign(tree_size + 1, 0);
    node_last_bit.assign(tree_size + 1, 0);
    leaf_prefix_stride = block_bits + 1;
    leaf_prefix.assign(leaf_count * leaf_prefix_stride, 0);

    // leaves
    for (size_t leaf_block_index = 0; leaf_block_index < leaf_count;
         ++leaf_block_index) {
      const size_t leaf_node_index = first_leaf_index + leaf_block_index;
      const size_t segment_begin = leaf_block_index * block_bits;
      const size_t segment_end = std::min(num_bits, segment_begin + block_bits);
      segment_size_bits[leaf_node_index] = segment_end - segment_begin;

      if (segment_begin < segment_end) {
        node_first_bit[leaf_node_index] = bit(segment_begin);
      }

      const auto& aggregates_table = LUT8();

      int current_excess = 0, min_value = INT_MAX, max_value = INT_MIN;
      uint32_t min_count = 0;
      uint32_t pattern10_count = 0;

      uint8_t previous_bit = 0;

      size_t position = segment_begin;

      // Full bytes
      while (position + 8 <= segment_end) {
        const uint8_t byte_value = get_byte(position);
        const auto& byte_aggregate = aggregates_table[byte_value];

        // internal "10" inside the byte
        pattern10_count += byte_aggregate.pattern10_count;
        // stitching across the boundary between the previous and current byte
        // (within the segment)
        if (previous_bit == 1 && byte_aggregate.first_bit == 0) {
          pattern10_count++;
        }

        // prefix min/max accounting for the current offset
        const int candidate_min = current_excess + byte_aggregate.min_prefix;
        if (candidate_min < min_value) {
          min_value = candidate_min;
          min_count = byte_aggregate.min_count;
        } else if (candidate_min == min_value) {
          min_count += byte_aggregate.min_count;
        }

        max_value =
            std::max(max_value, current_excess + byte_aggregate.max_prefix);
        current_excess += byte_aggregate.excess_total;
        previous_bit = byte_aggregate.last_bit;
        position += 8;
      }

      // Tail < 8 bits
      while (position < segment_end) {
        const uint8_t bit_value = bit(position);
        if (previous_bit == 1 && bit_value == 0) {
          pattern10_count++;
        }
        const int step = bit_value ? +1 : -1;
        current_excess += step;
        if (current_excess < min_value) {
          min_value = current_excess;
          min_count = 1;
        } else if (current_excess == min_value) {
          ++min_count;
        }
        if (current_excess > max_value) {
          max_value = current_excess;
        }

        previous_bit = bit_value;
        ++position;
      }

      if (segment_begin < segment_end) {
        node_last_bit[leaf_node_index] = previous_bit;
      }

      node_total_excess[leaf_node_index] = current_excess;
      node_min_prefix_excess[leaf_node_index] =
          (segment_size_bits[leaf_node_index] == 0 ? 0 : min_value);
      node_max_prefix_excess[leaf_node_index] =
          (segment_size_bits[leaf_node_index] == 0 ? 0 : max_value);
      node_min_count[leaf_node_index] = min_count;
      node_pattern10_count[leaf_node_index] = (uint32_t)pattern10_count;
      int16_t* leaf_prefix_ptr =
          &leaf_prefix[leaf_block_index * leaf_prefix_stride];
      int16_t prefix_accumulator = 0;
      leaf_prefix_ptr[0] = 0;
      for (size_t position_in_leaf = segment_begin, prefix_index = 1;
           position_in_leaf < segment_end; ++position_in_leaf, ++prefix_index) {
        prefix_accumulator += bit(position_in_leaf) ? 1 : -1;
        leaf_prefix_ptr[prefix_index] = prefix_accumulator;
      }
    }
    // internal nodes
    for (size_t node_index = first_leaf_index - 1; node_index >= 1;
         --node_index) {
      const size_t left_child = node_index << 1;
      const size_t right_child = left_child | 1;
      const bool has_left =
          (left_child <= tree_size) && segment_size_bits[left_child];
      const bool has_right =
          (right_child <= tree_size) && segment_size_bits[right_child];
      if (!has_left && !has_right) {
        segment_size_bits[node_index] = 0;
        continue;
      }
      if (has_left && !has_right) {
        segment_size_bits[node_index] = segment_size_bits[left_child];
        node_total_excess[node_index] = node_total_excess[left_child];
        node_min_prefix_excess[node_index] = node_min_prefix_excess[left_child];
        node_max_prefix_excess[node_index] = node_max_prefix_excess[left_child];
        node_min_count[node_index] = node_min_count[left_child];
        node_pattern10_count[node_index] = node_pattern10_count[left_child];
        node_first_bit[node_index] = node_first_bit[left_child];
        node_last_bit[node_index] = node_last_bit[left_child];
      } else if (!has_left && has_right) {
        segment_size_bits[node_index] = segment_size_bits[right_child];
        node_total_excess[node_index] = node_total_excess[right_child];
        node_min_prefix_excess[node_index] =
            node_min_prefix_excess[right_child];
        node_max_prefix_excess[node_index] =
            node_max_prefix_excess[right_child];
        node_min_count[node_index] = node_min_count[right_child];
        node_pattern10_count[node_index] = node_pattern10_count[right_child];
        node_first_bit[node_index] = node_first_bit[right_child];
        node_last_bit[node_index] = node_last_bit[right_child];
      } else {
        segment_size_bits[node_index] =
            segment_size_bits[left_child] + segment_size_bits[right_child];
        node_total_excess[node_index] =
            node_total_excess[left_child] + node_total_excess[right_child];
        const int right_min_candidate =
            node_total_excess[left_child] + node_min_prefix_excess[right_child];
        const int right_max_candidate =
            node_total_excess[left_child] + node_max_prefix_excess[right_child];
        node_min_prefix_excess[node_index] =
            std::min(node_min_prefix_excess[left_child], right_min_candidate);
        node_max_prefix_excess[node_index] =
            std::max(node_max_prefix_excess[left_child], right_max_candidate);
        node_min_count[node_index] =
            (node_min_prefix_excess[left_child] ==
                     node_min_prefix_excess[node_index]
                 ? node_min_count[left_child]
                 : 0) +
            (right_min_candidate == node_min_prefix_excess[node_index]
                 ? node_min_count[right_child]
                 : 0);
        node_pattern10_count[node_index] = node_pattern10_count[left_child] +
                                           node_pattern10_count[right_child] +
                                           ((node_last_bit[left_child] == 1 &&
                                             node_first_bit[right_child] == 0)
                                                ? 1u
                                                : 0u);
        node_first_bit[node_index] = node_first_bit[left_child];
        node_last_bit[node_index] = node_last_bit[right_child];
      }
      if (node_index == 1) {
        break;
      }
    }
  }
};

}  // namespace pixie
