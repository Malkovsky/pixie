#pragma once

#include <pixie/bits.h>
#include <pixie/memory_usage.h>
#include <pixie/rmq.h>

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
 * @brief Low-level selector implementation used by HybridBTree leaves.
 */
enum class HybridBTreeLeafSelector {
  PrefixSuffix,
  BP,
};

/**
 * @brief Hybrid B-tree RMQ with compact per-level selectors.
 *
 * @details This is a static, non-owning value-RMQ index over an external value
 * array. Values are split into leaves, leaves are grouped into a B-tree, and
 * each node stores a selector over the minima of its immediate children. Query
 * ranges use the library-wide half-open convention `[left, right)`, invalid
 * ranges return `npos`, and equal values resolve to the smaller original
 * position.
 *
 * The default low level uses 496-value prefix/suffix mask leaves. A leaf stores
 * one bit for each prefix/suffix record position plus a 16-bit local-minimum
 * offset in one 64-byte object. Prefix or suffix ranges can be answered from
 * this mask; interior partial ranges fall back to a linear scan of the original
 * values. The same implementation also supports 248-value mask leaves and
 * 252-value BP leaves for controlled experiments.
 *
 * Internal nodes use 192-way fanout. Their selector is a 512-bit local
 * Cartesian-tree balanced-parentheses encoding over child minima: 384 BP bits,
 * a 64-bit absolute subtree-minimum position, and 64 bits of zero-rank prefix
 * metadata. Node minimum values are cached in a side vector so comparing
 * candidates does not require repeatedly descending to leaves.
 *
 * Wide queries first try a single coarse sparse table over original-value block
 * minima. The sparse table uses at least 4096-value blocks and grows the block
 * width when needed so the top layer has at most 2^14 blocks. If the padded
 * block-cover candidate lies inside the query, it is returned immediately; on
 * a miss, the sparse table answers the fully covered middle block range and
 * the B-tree handles the two border ranges.
 *
 * @code
 * value array, n entries
 * |
 * +-- L0: leaves
 * |       ceil(n / 496) nodes by default
 * |       each leaf covers <=496 values and stores one 64-byte selector
 * |
 * +-- T: value-block sparse overlay
 * |       block width >=4096 values, at most 2^14 blocks
 * |       stores one original-value minimum position per sparse-table cell
 * |
 * +-- L1..Lm: internal levels
 * |       fanout 192
 * |       each node stores one 512-bit BP selector over child minima
 *
 * Examples with default 496-value leaves:
 *
 * n = 2^24:
 *   33,826 leaves -> 177 internal nodes -> 1 root
 *
 * n = 2^26:
 *   135,301 leaves -> 705 internal nodes -> 4 internal nodes -> 1 root
 * @endcode
 *
 * Querying starts at the lowest node that covers both endpoint leaves. At each
 * internal node, the selector is asked for the best child over the intersecting
 * child-slot range. If that child is fully covered, or its stored subtree
 * minimum lies inside the query, that candidate is final for the node.
 * Otherwise only the affected border child is corrected recursively and
 * compared with the middle full-child range. All comparisons use `(value,
 * position)` semantics so traversal order cannot change tie-breaking.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 * @tparam LeafSize Number of original values per leaf. This backend currently
 * supports 252 with BP leaves and 248 or 496 with mask leaves.
 * @tparam Fanout Template parameter kept fixed at 256 for the current layout;
 * internal tree nodes use fixed 192-way fanout.
 * @tparam LeafSelectorKind Compile-time low-level selector kind.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 496,
          std::size_t Fanout = 256,
          HybridBTreeLeafSelector LeafSelectorKind =
              HybridBTreeLeafSelector::PrefixSuffix>
class HybridBTree
    : public RmqBase<
          HybridBTree<T, Compare, Index, LeafSize, Fanout, LeafSelectorKind>,
          T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "HybridBTree index type must be unsigned");
  static constexpr bool kBpLeafSelector =
      LeafSelectorKind == HybridBTreeLeafSelector::BP;
  static constexpr bool kMaskLeafSelector =
      LeafSelectorKind == HybridBTreeLeafSelector::PrefixSuffix;
  static_assert((kBpLeafSelector && LeafSize == 252) ||
                    (kMaskLeafSelector && (LeafSize == 248 || LeafSize == 496)),
                "HybridBTree requires 252-value BP leaves "
                "or 248/496-value prefix/suffix mask leaves");
  static_assert(Fanout == 256,
                "HybridBTree currently requires a 256 fanout template "
                "argument");
  static_assert(kBpLeafSelector || kMaskLeafSelector,
                "unsupported HybridBTree low-level selector kind");

  using Self =
      HybridBTree<T, Compare, Index, LeafSize, Fanout, LeafSelectorKind>;

  static constexpr std::size_t npos = RmqBase<Self, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kLeafSize = LeafSize;
  static constexpr std::size_t kFanout = Fanout;
  static constexpr std::size_t kMiddleFanout = 192;
  static constexpr std::size_t kMinTopSparseBlockSize = 4096;
  static constexpr std::size_t kMaxTopSparseBlocks = std::size_t{1} << 14;

  /**
   * @brief Construct an empty RMQ index.
   */
  HybridBTree() = default;

  /**
   * @brief Copy an RMQ index while preserving its non-owning value span.
   */
  HybridBTree(const HybridBTree&) = default;

  /**
   * @brief Move an RMQ index while preserving selector and cache storage.
   */
  HybridBTree(HybridBTree&&) noexcept = default;

  /**
   * @brief Copy-assign an RMQ index and its cached metadata.
   */
  HybridBTree& operator=(const HybridBTree&) = default;

  /**
   * @brief Move-assign an RMQ index and its cached metadata.
   */
  HybridBTree& operator=(HybridBTree&&) noexcept = default;

  /**
   * @brief Build a hybrid B-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller original position as the answer.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit HybridBTree(std::span<const T> values, Compare compare = Compare())
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

    if (right - left <= top_block_size_) {
      return tree_arg_min(left, right);
    }
    const std::size_t top_answer = top_sparse_arg_min(left, right);
    return top_answer != npos ? top_answer : tree_arg_min(left, right);
  }

  /**
   * @brief Return the top sparse-table block width chosen for a value count.
   */
  static std::size_t top_sparse_block_size_for(std::size_t value_count) {
    if (value_count == 0) {
      return kMinTopSparseBlockSize;
    }
    return std::max(kMinTopSparseBlockSize,
                    1 + (value_count - 1) / kMaxTopSparseBlocks);
  }

  /**
   * @brief Return the number of top sparse-table blocks for a value count.
   */
  static std::size_t top_sparse_block_count_for(std::size_t value_count) {
    if (value_count == 0) {
      return 0;
    }
    const std::size_t block_size = top_sparse_block_size_for(value_count);
    return 1 + (value_count - 1) / block_size;
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
   * @details Counts this RMQ object and all selector/metadata buffers. The
   * external input values are not owned and are excluded.
   */
  std::size_t memory_usage_bytes_impl() const {
    std::size_t bytes = sizeof(*this);
    bytes += pixie::vector_capacity_bytes(leaf_selectors_);
    bytes += pixie::vector_capacity_bytes(medium_selectors_);
    bytes += pixie::vector_capacity_bytes(medium_min_values_);
    bytes += pixie::vector_capacity_bytes(top_sparse_candidates_);
    bytes += pixie::vector_capacity_bytes(medium_level_offsets_);
    bytes += pixie::vector_capacity_bytes(level_sizes_);
    bytes += pixie::vector_capacity_bytes(level_value_spans_);
    bytes += pixie::vector_capacity_bytes(level_fanouts_);
    return bytes;
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
  static constexpr std::size_t kLeafLinearScanThreshold = 64;
  static constexpr std::size_t kLeafAvx2ScanThreshold = 16;
  static constexpr std::uint64_t kEmbeddedOffsetMask =
      std::uint64_t{0xFF} << (kEmbeddedOffsetBit & 63);
  static constexpr bool kInvalidIndexEqualsNpos =
      static_cast<std::size_t>(invalid_index) == npos;

  static_assert(kEmbeddedOffsetBit + 8 == kSelectorBits);
  static_assert(kEmbeddedPositionBit + 128 == kSelectorBits);
  static_assert(!kBpLeafSelector || LeafSize <= kEmbeddedOffsetEntries);
  static_assert(kMiddleFanout <= kEmbeddedPositionEntries);

  struct MinCandidate {
    std::size_t position = npos;
    const T* value = nullptr;
  };

  struct TopCandidate {
    Index position = invalid_index;
  };
  static_assert(sizeof(TopCandidate) == sizeof(Index));

  class alignas(64) Bp512Selector {
   public:
    /**
     * @brief Construct an empty packed BP selector.
     */
    Bp512Selector() = default;

    /**
     * @brief Build a stable local Cartesian-tree BP selector.
     *
     * @details The selector indexes @p entry_count logical entries. The
     * supplied comparator returns whether the left entry is strictly better
     * than the right entry. Equal entries keep the smaller slot by using the
     * same stable rule as the value-level RMQ.
     */
    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kSelectorEntries) {
        throw std::length_error("HybridBTree local selector too large");
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
     *
     * @details This path is used by BP leaves and high nodes whose BP words do
     * not share storage with zero-rank prefix metadata.
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
     * @brief Return the first minimum slot using packed zero-rank metadata.
     *
     * @details Middle nodes overwrite the last selector word with zero-rank
     * prefixes, so rank/select operations must consult the packed metadata
     * instead of treating all eight words as raw BP bits.
     */
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
     * @brief Pack a leaf-local minimum offset into the spare BP bits.
     */
    void set_embedded_min_offset(std::size_t offset) {
      bp_bits_[kEmbeddedOffsetBit >> 6] =
          (bp_bits_[kEmbeddedOffsetBit >> 6] & ~kEmbeddedOffsetMask) |
          ((static_cast<std::uint64_t>(offset) & 0xFFu)
           << (kEmbeddedOffsetBit & 63));
    }

    /**
     * @brief Return the leaf-local minimum offset stored in the selector.
     */
    std::uint8_t embedded_min_offset() const {
      return static_cast<std::uint8_t>(
          (bp_bits_[kEmbeddedOffsetBit >> 6] & kEmbeddedOffsetMask) >>
          (kEmbeddedOffsetBit & 63));
    }

    /**
     * @brief Store a node's absolute subtree minimum position.
     */
    void set_embedded_min_position(std::size_t position) {
      bp_bits_[kEmbeddedPositionBit >> 6] =
          static_cast<std::uint64_t>(position);
    }

    /**
     * @brief Return the absolute subtree minimum position stored in a node.
     */
    std::size_t embedded_min_position() const {
      return static_cast<std::size_t>(bp_bits_[kEmbeddedPositionBit >> 6]);
    }

    /**
     * @brief Pack per-word zero-count prefixes for middle-node selectors.
     */
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
    /**
     * @brief Prepend one BP bit while building the sequence right-to-left.
     */
    std::size_t prepend_bp_bit(std::size_t& write_position, bool bit) {
      --write_position;
      if (bit) {
        bp_bits_[write_position >> 6] |= std::uint64_t{1}
                                         << (write_position & 63);
      }
      return write_position;
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
     * @brief Return the packed zero-count prefix at a 64-bit word boundary.
     */
    std::size_t zero_prefix_at_word(std::size_t word) const {
      return static_cast<std::uint8_t>(bp_bits_[kMiddleZeroPrefixWord] >>
                                       (8 * word));
    }

    /**
     * @brief Count close parentheses before @p position using packed prefixes.
     */
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

    /**
     * @brief Return BP excess when the last word stores zero-rank metadata.
     */
    int prefix_excess_with_zero_prefix(std::size_t position,
                                       std::size_t bit_count) const {
      position = std::min(position, bit_count);
      return static_cast<int>(position) -
             2 * static_cast<int>(
                     rank0_at_with_zero_prefix(position, bit_count));
    }

    /**
     * @brief Select a close parenthesis using packed per-word zero counts.
     */
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

    /**
     * @brief Return the position of the minimum BP depth in a depth range.
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

    /**
     * @brief Return the minimum-depth position for middle-node packed BP data.
     */
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

  class alignas(LeafSize == 496 ? 64 : 32) PrefixSuffixMaskLeafSelector {
   public:
    /**
     * @brief Construct an empty prefix/suffix mask selector.
     */
    PrefixSuffixMaskLeafSelector() = default;

    /**
     * @brief Build prefix/suffix record masks and store the local minimum.
     *
     * @details A bit is set for every prefix-minimum record and every
     * suffix-minimum record. Equal suffix values set the smaller slot so
     * first-minimum tie-breaking is preserved.
     */
    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kMaskEntries) {
        throw std::length_error(
            "HybridBTree prefix/suffix leaf selector too large");
      }

      words_.fill(0);
      if (entry_count == 0) {
        set_embedded_min_offset(0);
        return;
      }

      std::size_t prefix_best = 0;
      set_mask_bit(0);
      for (std::size_t slot = 1; slot < entry_count; ++slot) {
        if (entry_less(slot, prefix_best)) {
          prefix_best = slot;
          set_mask_bit(slot);
        }
      }

      std::size_t suffix_best = entry_count - 1;
      set_mask_bit(suffix_best);
      for (std::size_t slot = entry_count - 1; slot > 0;) {
        --slot;
        if (!entry_less(suffix_best, slot)) {
          suffix_best = slot;
          set_mask_bit(slot);
        }
      }

      set_embedded_min_offset(prefix_best);
    }

    /**
     * @brief Return a local minimum slot when the mask can answer directly.
     *
     * @details The mask answers full-leaf ranges, prefixes before the global
     * leaf minimum, and suffixes after the global leaf minimum. Interior ranges
     * that cannot be represented by prefix/suffix records return `npos` so the
     * caller can scan original values.
     */
    std::size_t arg_min(std::size_t slot_left,
                        std::size_t slot_right,
                        std::size_t entry_count) const {
      if (slot_left >= slot_right || slot_right > entry_count ||
          entry_count > kMaskEntries) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t min_offset = embedded_min_offset();
      if (slot_left <= min_offset && min_offset < slot_right) {
        return min_offset;
      }
      if (slot_left == 0 && slot_right <= min_offset) {
        return previous_set_bit_before(slot_right);
      }
      if (slot_right == entry_count && min_offset < slot_left) {
        return next_set_bit_at_or_after(slot_left);
      }
      return npos;
    }

    /**
     * @brief Pack the local minimum offset into the high bits of the mask.
     */
    void set_embedded_min_offset(std::size_t offset) {
      words_[kOffsetWord] =
          (words_[kOffsetWord] & kMaskWordMask) |
          ((static_cast<std::uint64_t>(offset) & kOffsetMask) << kOffsetShift);
    }

    /**
     * @brief Return the local minimum offset packed into the mask words.
     */
    std::size_t embedded_min_offset() const {
      return static_cast<std::size_t>((words_[kOffsetWord] >> kOffsetShift) &
                                      kOffsetMask);
    }

   private:
    static constexpr std::size_t kMaskEntries = LeafSize;
    static constexpr std::size_t kOffsetBits = LeafSize == 496 ? 16 : 8;
    static constexpr std::size_t kPackedBits = kMaskEntries + kOffsetBits;
    static constexpr std::size_t kWordCount = (kPackedBits + 63) / 64;
    static constexpr std::size_t kOffsetWord = kMaskEntries / 64;
    static constexpr std::size_t kOffsetShift = kMaskEntries & 63;

    /**
     * @brief Return a mask with the lowest @p bits set.
     */
    static constexpr std::uint64_t low_bits_mask(std::size_t bits) {
      if (bits == 0) {
        return 0;
      }
      if (bits >= 64) {
        return std::numeric_limits<std::uint64_t>::max();
      }
      return (std::uint64_t{1} << bits) - 1;
    }

    static constexpr std::uint64_t kOffsetMask = low_bits_mask(kOffsetBits);
    static constexpr std::uint64_t kMaskWordMask = low_bits_mask(kOffsetShift);

    static_assert(!kMaskLeafSelector || kPackedBits % 64 == 0);
    static_assert(!kMaskLeafSelector || kOffsetShift + kOffsetBits == 64);

    /**
     * @brief Set the prefix/suffix-record bit for a local leaf slot.
     */
    void set_mask_bit(std::size_t slot) {
      words_[slot >> 6] |= std::uint64_t{1} << (slot & 63);
    }

    /**
     * @brief Return a mask word with embedded offset bits hidden.
     */
    std::uint64_t mask_word(std::size_t word) const {
      return word == kOffsetWord ? words_[word] & kMaskWordMask : words_[word];
    }

    /**
     * @brief Return the last set record bit strictly before @p limit.
     */
    std::size_t previous_set_bit_before(std::size_t limit) const {
      if (limit == 0) {
        return npos;
      }

      std::size_t word = (limit - 1) >> 6;
      std::uint64_t bits =
          mask_word(word) & first_bits_mask(((limit - 1) & 63) + 1);
      while (true) {
        if (bits != 0) {
          return word * 64 + 63 - std::countl_zero(bits);
        }
        if (word == 0) {
          break;
        }
        --word;
        bits = mask_word(word);
      }
      return npos;
    }

    /**
     * @brief Return the first set record bit at or after @p slot.
     */
    std::size_t next_set_bit_at_or_after(std::size_t slot) const {
      std::size_t word = slot >> 6;
      std::uint64_t bits = mask_word(word) & ~first_bits_mask(slot & 63);
      while (word < kWordCount) {
        if (bits != 0) {
          const std::size_t result = word * 64 + std::countr_zero(bits);
          return result < kMaskEntries ? result : npos;
        }
        ++word;
        bits = word < kWordCount ? mask_word(word) : 0;
      }
      return npos;
    }

    std::array<std::uint64_t, kWordCount> words_{};
  };

  static_assert(!kMaskLeafSelector || sizeof(PrefixSuffixMaskLeafSelector) ==
                                          (LeafSize == 496 ? 64 : 32));
  static_assert(!kMaskLeafSelector || alignof(PrefixSuffixMaskLeafSelector) ==
                                          (LeafSize == 496 ? 64 : 32));

  using LeafSelector = std::conditional_t<kMaskLeafSelector,
                                          PrefixSuffixMaskLeafSelector,
                                          Bp512Selector>;

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
   * @brief Wrap an original value position as a comparable candidate.
   */
  MinCandidate value_candidate(std::size_t position) const {
    if (missing_position(position) || position >= values_.size()) {
      return {};
    }
    return {position, values_.data() + position};
  }

  /**
   * @brief Return a node's cached subtree-minimum candidate.
   */
  MinCandidate subtree_min_candidate(std::size_t level,
                                     std::size_t node) const {
    const std::size_t position = subtree_min_position(level, node);
    if (missing_position(position)) {
      return {};
    }
    return {position, &subtree_min_value(level, node)};
  }

  /**
   * @brief Return a child slot's subtree-minimum candidate.
   */
  MinCandidate subtree_child_min_candidate(std::size_t level,
                                           std::size_t node,
                                           std::size_t slot) const {
    const std::size_t child_level = level - 1;
    return subtree_min_candidate(child_level,
                                 node * fanout_at_level(level) + slot);
  }

  /**
   * @brief Return whether @p left is strictly better than @p right.
   */
  bool strictly_better_candidate(MinCandidate left, MinCandidate right) const {
    if (missing_position(left.position)) {
      return false;
    }
    if (missing_position(right.position)) {
      return true;
    }
    return compare_(*left.value, *right.value);
  }

  /**
   * @brief Choose the better candidate using value, then smaller position.
   */
  MinCandidate better_candidate(MinCandidate left, MinCandidate right) const {
    if (missing_position(left.position)) {
      return right;
    }
    if (missing_position(right.position)) {
      return left;
    }
    if (compare_(*right.value, *left.value)) {
      return right;
    }
    if (compare_(*left.value, *right.value)) {
      return left;
    }
    return right.position < left.position ? right : left;
  }

  /**
   * @brief Compare two regular child slots while building a local selector.
   */
  bool strictly_better_subtree_child_slot(std::size_t level,
                                          std::size_t node,
                                          std::size_t left_slot,
                                          std::size_t right_slot) const {
    return strictly_better_candidate(
        subtree_child_min_candidate(level, node, left_slot),
        subtree_child_min_candidate(level, node, right_slot));
  }

  /**
   * @brief Build all levels, selectors, and cached minimum metadata.
   */
  void build() {
    leaf_selectors_.clear();
    medium_selectors_.clear();
    medium_min_values_.clear();
    top_sparse_candidates_.clear();
    top_block_size_ = kMinTopSparseBlockSize;
    top_block_count_ = 0;
    top_sparse_levels_ = 0;
    medium_level_offsets_.clear();
    level_sizes_.clear();
    level_value_spans_.clear();
    level_fanouts_.clear();
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("HybridBTree index type is too small");
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
    build_top_sparse_table();
  }

  /**
   * @brief Compute level sizes, fanouts, flat offsets, and cache storage.
   *
   * @details Level zero always stores leaf selectors. Every internal level uses
   * the fixed 192-way middle-node layout.
   */
  void initialize_layout(std::size_t leaf_count) {
    level_sizes_.push_back(leaf_count);
    level_value_spans_.push_back(LeafSize);
    level_fanouts_.push_back(0);

    std::size_t current_count = leaf_count;
    std::size_t current_span = LeafSize;
    while (current_count > 1) {
      level_fanouts_.push_back(kMiddleFanout);
      current_count = ceil_div(current_count, kMiddleFanout);
      current_span = saturating_product(current_span, kMiddleFanout);
      level_sizes_.push_back(current_count);
      level_value_spans_.push_back(current_span);
    }

    leaf_selectors_.resize(level_sizes_[0]);

    medium_level_offsets_.assign(level_count(), 0);
    if (level_count() <= 1) {
      return;
    }

    std::size_t medium_node_count = 0;
    for (std::size_t level = 1; level < level_count(); ++level) {
      medium_level_offsets_[level] = medium_node_count;
      medium_node_count += level_sizes_[level];
    }
    medium_selectors_.resize(medium_node_count);
    medium_min_values_.reserve(medium_node_count);
  }

  /**
   * @brief Build one leaf selector over original values.
   */
  void build_leaf(std::size_t leaf) {
    LeafSelector& selector = leaf_selectors_[leaf];
    const std::size_t begin = node_value_begin(0, leaf);
    const std::size_t count = entry_count(0, leaf);
    selector.build(count, [&](std::size_t left, std::size_t right) {
      return compare_(values_[begin + left], values_[begin + right]);
    });

    if constexpr (kBpLeafSelector) {
      const std::size_t slot = selector.arg_min(0, count, count);
      selector.set_embedded_min_offset(slot);
    }
  }

  /**
   * @brief Build one internal node selector and its cached minima.
   */
  void build_internal_node(std::size_t level, std::size_t node) {
    Bp512Selector& selector = mutable_selector_at(level, node);
    const std::size_t count = entry_count(level, node);
    const std::size_t first_child = node * fanout_at_level(level);

    selector.build(count, [&](std::size_t left, std::size_t right) {
      return strictly_better_subtree_child_slot(level, node, left, right);
    });

    const std::size_t slot = selector.arg_min(0, count, count);
    const std::size_t min_position =
        subtree_min_position(level - 1, first_child + slot);
    selector.set_embedded_min_position(min_position);
    medium_min_values_.push_back(values_[min_position]);
    selector.build_zero_prefix_metadata(2 * count);
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
   * @brief Return `ceil(value / divisor)` for positive divisors.
   */
  static std::size_t ceil_div(std::size_t value, std::size_t divisor) {
    return value == 0 ? 0 : 1 + (value - 1) / divisor;
  }

  /**
   * @brief Return the first minimum position through the B-tree fallback.
   */
  std::size_t tree_arg_min(std::size_t left, std::size_t right) const {
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
    return query_node(level, node_index, left, right).position;
  }

  /**
   * @brief Return the first minimum position for a possibly partial leaf.
   *
   * @details Full leaves use the embedded minimum. Mask leaves answer prefix
   * and suffix ranges through their record bits and fall back to scanning only
   * for unsupported interior ranges.
   */
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
    const std::size_t slot_left = left - begin;
    const std::size_t slot_right = right - begin;

    if constexpr (kMaskLeafSelector) {
      const std::size_t slot = leaf_selectors_[leaf].arg_min(
          slot_left, slot_right, entry_count(0, leaf));
      if (slot != npos) {
        return begin + slot;
      }
      return linear_range_min(left, right);
    } else {
      if (right - left <= kLeafLinearScanThreshold) {
        return linear_range_min(left, right);
      }

      const std::size_t slot = leaf_selectors_[leaf].arg_min(
          slot_left, slot_right, entry_count(0, leaf));
      if (slot == npos) {
        return npos;
      }
      return begin + slot;
    }
  }

  /**
   * @brief Scan original values to answer a range inside a leaf.
   */
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
  /**
   * @brief AVX2 scan for signed 64-bit minimum ranges.
   *
   * @details Vector lanes track the best value and first position seen so far.
   * Equal values keep the earlier position by only replacing lanes on strict
   * value improvement.
   */
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
   * @brief Return the leaf index containing an original value position.
   */
  std::size_t leaf_for_value(std::size_t position) const {
    return position / LeafSize;
  }

  /**
   * @brief Return the child index at @p child_level containing a value.
   */
  std::size_t child_for_value(std::size_t child_level,
                              std::size_t position) const {
    return position / level_value_spans_[child_level];
  }

  /**
   * @brief Query a child-slot range of an internal node.
   *
   * @details The node selector first picks the best child in the requested slot
   * range. If that child cannot be accepted directly because only a border part
   * is queried, the two border children are queried recursively and the middle
   * full-child range is answered by the current node selector.
   */
  MinCandidate query_child_slots(std::size_t level,
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
    const MinCandidate child_min = subtree_min_candidate(child_level, child);
    const std::size_t child_begin = node_value_begin(child_level, child);
    const std::size_t child_end = node_value_end(child_level, child);
    if ((left <= child_begin && child_end <= right) ||
        contains_position(left, right, child_min.position)) {
      return child_min;
    }

    const std::size_t last_slot = slot_right - 1;
    const std::size_t left_child_begin =
        node_value_begin(child_level, first_child + slot_left);
    const std::size_t left_child_end =
        node_value_end(child_level, first_child + slot_left);
    MinCandidate answer = query_node(child_level, first_child + slot_left,
                                     std::max(left, left_child_begin),
                                     std::min(right, left_child_end));

    if (slot_left != last_slot) {
      const std::size_t right_child_begin =
          node_value_begin(child_level, first_child + last_slot);
      const std::size_t right_child_end =
          node_value_end(child_level, first_child + last_slot);
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
  MinCandidate full_child_slot_range_min(std::size_t level,
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

    return subtree_child_min_candidate(level, node, slot);
  }

  /**
   * @brief Query a tree node for the minimum candidate in a value range.
   */
  MinCandidate query_node(std::size_t level,
                          std::size_t node,
                          std::size_t left,
                          std::size_t right) const {
    if (left >= right) {
      return {};
    }
    const std::size_t begin = node_value_begin(level, node);
    const std::size_t end = node_value_end(level, node);
    if (left <= begin && end <= right) {
      return subtree_min_candidate(level, node);
    }
    if (level == 0) {
      return value_candidate(leaf_range_min(node, left, right));
    }

    const std::size_t child_level = level - 1;
    const std::size_t left_child = child_for_value(child_level, left);
    const std::size_t right_child = child_for_value(child_level, right - 1);
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
      const std::size_t begin = node_value_begin(0, node);
      return std::min<std::size_t>(LeafSize, values_.size() - begin);
    }
    const std::size_t first_child = node * fanout_at_level(level);
    return std::min<std::size_t>(fanout_at_level(level),
                                 level_sizes_[level - 1] - first_child);
  }

  /**
   * @brief Return the first original value position covered by a node.
   */
  std::size_t node_value_begin(std::size_t level, std::size_t node) const {
    return node * level_value_spans_[level];
  }

  /**
   * @brief Return one past the last original value position covered by a node.
   */
  std::size_t node_value_end(std::size_t level, std::size_t node) const {
    return std::min(values_.size(),
                    node_value_begin(level, node) + level_value_spans_[level]);
  }

  /**
   * @brief Return a node's absolute subtree-minimum position.
   */
  std::size_t subtree_min_position(std::size_t level, std::size_t node) const {
    if (level == 0) {
      return node_value_begin(0, node) +
             leaf_selectors_[node].embedded_min_offset();
    }
    return selector_at(level, node).embedded_min_position();
  }

  /**
   * @brief Return a cached reference to a node's subtree-minimum value.
   */
  const T& subtree_min_value(std::size_t level, std::size_t node) const {
    if (level == 0) {
      return values_[subtree_min_position(0, node)];
    }
    return medium_min_values_[medium_flat_index(level, node)];
  }

  /**
   * @brief Return an immutable internal-node BP selector.
   */
  const Bp512Selector& selector_at(std::size_t level, std::size_t node) const {
    return medium_selectors_[medium_flat_index(level, node)];
  }

  /**
   * @brief Return a mutable internal-node BP selector while building.
   */
  Bp512Selector& mutable_selector_at(std::size_t level, std::size_t node) {
    return medium_selectors_[medium_flat_index(level, node)];
  }

  /**
   * @brief Run the appropriate local selector for an internal node.
   */
  std::size_t selector_arg_min(std::size_t level,
                               std::size_t node,
                               std::size_t slot_left,
                               std::size_t slot_right,
                               std::size_t count) const {
    const Bp512Selector& selector = selector_at(level, node);
    return selector.arg_min_with_zero_prefix(slot_left, slot_right, count);
  }

  /**
   * @brief Map a medium-level node to its flat storage index.
   */
  std::size_t medium_flat_index(std::size_t level, std::size_t node) const {
    return medium_level_offsets_[level] + node;
  }

  /**
   * @brief Return the fanout used to group children at a level.
   */
  std::size_t fanout_at_level(std::size_t level) const {
    return level_fanouts_[level];
  }

  /**
   * @brief Build a single top sparse table over original-value block minima.
   */
  void build_top_sparse_table() {
    top_sparse_candidates_.clear();
    top_block_size_ = top_sparse_block_size_for(values_.size());
    top_block_count_ = top_sparse_block_count_for(values_.size());
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
   * @brief Wrap an original value position as a top sparse-table candidate.
   */
  TopCandidate make_top_candidate(std::size_t position) const {
    if (!valid_value_position(position)) {
      return {};
    }
    return {static_cast<Index>(position)};
  }

  /**
   * @brief Return whether @p position is a valid original-value index.
   */
  bool valid_value_position(std::size_t position) const {
    return !missing_position(position) && position < values_.size();
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
   * @brief Return the top-overlay answer, or `npos` when the tree should run.
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
          answer, make_top_candidate(tree_arg_min(left, left_border_end)));
    }

    const std::size_t right_border_begin = full_block_right * top_block_size_;
    if (right_border_begin < right) {
      answer = better_top_candidate(
          answer, make_top_candidate(tree_arg_min(right_border_begin, right)));
    }

    return valid_value_position(static_cast<std::size_t>(answer.position))
               ? static_cast<std::size_t>(answer.position)
               : npos;
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<LeafSelector> leaf_selectors_;
  std::vector<Bp512Selector> medium_selectors_;
  std::vector<T> medium_min_values_;
  std::vector<TopCandidate> top_sparse_candidates_;
  std::vector<std::size_t> medium_level_offsets_;
  std::vector<std::size_t> level_sizes_;
  std::vector<std::size_t> level_value_spans_;
  std::vector<std::size_t> level_fanouts_;
  std::size_t top_block_size_ = kMinTopSparseBlockSize;
  std::size_t top_block_count_ = 0;
  std::size_t top_sparse_levels_ = 0;
};

}  // namespace pixie::rmq
