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
 * @brief Low-level selector tag for BP-based 252-value leaves.
 */
struct BpLeafSelectorTag {};

/**
 * @brief Low-level selector tag for prefix/suffix mask leaves.
 */
struct PrefixSuffixMaskLeafSelectorTag {};

/**
 * @brief Segment B-tree RMQ with compact per-level selectors and XL leaves.
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
 * Middle internal nodes use 192-way fanout. Their selector is a 512-bit local
 * Cartesian-tree balanced-parentheses encoding over child minima: 384 BP bits,
 * a 64-bit absolute subtree-minimum position, and 64 bits of zero-rank prefix
 * metadata. Middle-node minimum values are cached in a side vector so comparing
 * candidates does not require repeatedly descending to leaves.
 *
 * High nodes use 256-way fanout and are limited to the root and its child
 * level, so there are at most 257 high nodes. They keep the same local BP
 * selector, cache child value ranges and child subtree minima, and add sparse
 * tables over child-minimum slots. This spends more space at the top of the
 * tree to reduce work on wide queries while keeping the high-level metadata
 * small enough to stay cache-resident.
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
 * @tparam Fanout Number of children per high internal node. This backend
 * currently supports only 256. Middle internal nodes use 192.
 * @tparam LowLevelSelector Compile-time low-level selector tag.
 */
template <class T,
          class Compare = std::less<T>,
          class Index = std::size_t,
          std::size_t LeafSize = 496,
          std::size_t Fanout = 256,
          class LowLevelSelector = PrefixSuffixMaskLeafSelectorTag>
class SegmentBTreeXl
    : public RmqBase<
          SegmentBTreeXl<T, Compare, Index, LeafSize, Fanout, LowLevelSelector>,
          T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "SegmentBTreeXl index type must be unsigned");
  static constexpr bool kBpLeafSelector =
      std::is_same_v<LowLevelSelector, BpLeafSelectorTag>;
  static constexpr bool kMaskLeafSelector =
      std::is_same_v<LowLevelSelector, PrefixSuffixMaskLeafSelectorTag>;
  static_assert((kBpLeafSelector && LeafSize == 252) ||
                    (kMaskLeafSelector && (LeafSize == 248 || LeafSize == 496)),
                "SegmentBTreeXl requires 252-value BP leaves "
                "or 248/496-value prefix/suffix mask leaves");
  static_assert(Fanout == 256,
                "SegmentBTreeXl currently requires 256-way internal nodes");
  static_assert(kBpLeafSelector || kMaskLeafSelector,
                "unsupported SegmentBTreeXl low-level selector tag");

  using Self =
      SegmentBTreeXl<T, Compare, Index, LeafSize, Fanout, LowLevelSelector>;

  static constexpr std::size_t npos = RmqBase<Self, T>::npos;
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();
  static constexpr std::size_t kLeafSize = LeafSize;
  static constexpr std::size_t kFanout = Fanout;
  static constexpr std::size_t kMiddleFanout = 192;

  /**
   * @brief Construct an empty RMQ index.
   */
  SegmentBTreeXl() = default;

  /**
   * @brief Copy an RMQ index while preserving its non-owning value span.
   */
  SegmentBTreeXl(const SegmentBTreeXl&) = default;

  /**
   * @brief Move an RMQ index while preserving selector and cache storage.
   */
  SegmentBTreeXl(SegmentBTreeXl&&) noexcept = default;

  /**
   * @brief Copy-assign an RMQ index and its cached metadata.
   */
  SegmentBTreeXl& operator=(const SegmentBTreeXl&) = default;

  /**
   * @brief Move-assign an RMQ index and its cached metadata.
   */
  SegmentBTreeXl& operator=(SegmentBTreeXl&&) noexcept = default;

  /**
   * @brief Build a segment B-tree RMQ index over @p values.
   *
   * @details The values are not copied and must outlive this object. Equal
   * values keep the smaller original position as the answer.
   *
   * @param values Values to index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if @p Index cannot represent all positions.
   */
  explicit SegmentBTreeXl(std::span<const T> values,
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
    return query_node(level, node_index, left, right).position;
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
  static_assert(!kBpLeafSelector || LeafSize <= kEmbeddedOffsetEntries);
  static_assert(kMiddleFanout <= kEmbeddedPositionEntries);

  struct HighChildMetadata {
    std::size_t value_begin = 0;
    std::size_t value_end = 0;
    Index min_position = invalid_index;
  };

  struct MinCandidate {
    std::size_t position = npos;
    const T* value = nullptr;
  };

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
        throw std::length_error("SegmentBTreeXl local selector too large");
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
            "SegmentBTreeXl prefix/suffix leaf selector too large");
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
    if (missing_position(position)) {
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
   * @brief Return a high-node child candidate using high-child metadata.
   */
  MinCandidate high_child_min_candidate(std::size_t level,
                                        std::size_t node,
                                        std::size_t slot) const {
    const std::size_t child_level = level - 1;
    const std::size_t child = node * fanout_at_level(level) + slot;
    const std::size_t position =
        high_child_metadata_at(high_flat_index(level, node), slot).min_position;
    if (missing_position(position)) {
      return {};
    }
    return {position, &subtree_min_value(child_level, child)};
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
   * @brief Compare two high-node child slots while building high metadata.
   */
  bool strictly_better_high_child_slot(std::size_t level,
                                       std::size_t node,
                                       std::size_t left_slot,
                                       std::size_t right_slot) const {
    return strictly_better_candidate(
        high_child_min_candidate(level, node, left_slot),
        high_child_min_candidate(level, node, right_slot));
  }

  /**
   * @brief Build all levels, selectors, and cached minimum metadata.
   */
  void build() {
    leaf_selectors_.clear();
    medium_selectors_.clear();
    medium_min_values_.clear();
    high_selectors_.clear();
    high_min_positions_.clear();
    high_min_values_.clear();
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
      throw std::length_error("SegmentBTreeXl index type is too small");
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

  /**
   * @brief Compute level sizes, fanouts, flat offsets, and cache storage.
   *
   * @details The top two tree levels are marked as high levels. Levels below
   * them use the middle-node layout. Level zero always stores leaf selectors.
   */
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
    medium_min_values_.reserve(medium_node_count);

    std::size_t high_node_count = 0;
    for (std::size_t level = high_level_begin_; level < level_count();
         ++level) {
      high_level_offsets_[level] = high_node_count;
      high_node_count += level_sizes_[level];
    }
    high_selectors_.resize(high_node_count);
    high_min_positions_.resize(high_node_count, invalid_index);
    high_min_values_.reserve(high_node_count);
    high_child_metadata_.resize(high_node_count * Fanout);
    high_sparse_min_slots_.resize(high_node_count * kHighSparseSlotsPerNode);
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
      build_high_sparse_min_slots(level, node, count);
    }

    selector.build(count, [&](std::size_t left, std::size_t right) {
      if (high_level) {
        return strictly_better_high_child_slot(level, node, left, right);
      }
      return strictly_better_subtree_child_slot(level, node, left, right);
    });

    const std::size_t slot = selector.arg_min(0, count, count);
    const std::size_t min_position =
        high_level ? high_child_metadata_at(high_flat, slot).min_position
                   : subtree_min_position(level - 1, first_child + slot);
    if (high_level) {
      high_min_positions_[high_flat] = static_cast<Index>(min_position);
      high_min_values_.push_back(values_[min_position]);
    } else {
      selector.set_embedded_min_position(min_position);
      medium_min_values_.push_back(values_[min_position]);
      selector.build_zero_prefix_metadata(2 * count);
    }
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
    return (value + divisor - 1) / divisor;
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
    const HighChildMetadata* high_children =
        is_high_level(level) ? high_child_metadata_begin(level, node) : nullptr;
    const MinCandidate child_min =
        high_children != nullptr ? high_child_min_candidate(level, node, slot)
                                 : subtree_min_candidate(child_level, child);
    const std::size_t child_begin = high_children != nullptr
                                        ? high_children[slot].value_begin
                                        : node_value_begin(child_level, child);
    const std::size_t child_end = high_children != nullptr
                                      ? high_children[slot].value_end
                                      : node_value_end(child_level, child);
    if ((left <= child_begin && child_end <= right) ||
        contains_position(left, right, child_min.position)) {
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
    MinCandidate answer = query_node(child_level, first_child + slot_left,
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

    if (is_high_level(level)) {
      return high_child_min_candidate(level, node, slot);
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
    if (is_high_level(level)) {
      return high_min_positions_[high_flat_index(level, node)];
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
    if (is_high_level(level)) {
      return high_min_values_[high_flat_index(level, node)];
    }
    return medium_min_values_[medium_flat_index(level, node)];
  }

  /**
   * @brief Return an immutable internal-node BP selector.
   */
  const Bp512Selector& selector_at(std::size_t level, std::size_t node) const {
    if (is_high_level(level)) {
      return high_selectors_[high_flat_index(level, node)];
    }
    return medium_selectors_[medium_flat_index(level, node)];
  }

  /**
   * @brief Return a mutable internal-node BP selector while building.
   */
  Bp512Selector& mutable_selector_at(std::size_t level, std::size_t node) {
    if (is_high_level(level)) {
      return high_selectors_[high_flat_index(level, node)];
    }
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
    if (is_high_level(level)) {
      return high_sparse_arg_min(level, node, slot_left, slot_right, count);
    }
    const Bp512Selector& selector = selector_at(level, node);
    return selector.arg_min_with_zero_prefix(slot_left, slot_right, count);
  }

  /**
   * @brief Return whether a level uses the high-node layout.
   */
  bool is_high_level(std::size_t level) const {
    return level > 0 && level >= high_level_begin_ && level < level_count();
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
   * @brief Map a high-level node to its flat storage index.
   */
  std::size_t high_flat_index(std::size_t level, std::size_t node) const {
    return high_level_offsets_[level] + node;
  }

  /**
   * @brief Return the first child-metadata record for a high node.
   */
  const HighChildMetadata* high_child_metadata_begin(std::size_t level,
                                                     std::size_t node) const {
    return high_child_metadata_.data() + high_flat_index(level, node) * Fanout;
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
   * @brief Return high-child metadata by flat high-node index and slot.
   */
  const HighChildMetadata& high_child_metadata_at(std::size_t high_flat,
                                                  std::size_t slot) const {
    return high_child_metadata_[high_flat * Fanout + slot];
  }

  /**
   * @brief Return mutable high-child metadata while building.
   */
  HighChildMetadata& mutable_high_child_metadata_at(std::size_t high_flat,
                                                    std::size_t slot) {
    return high_child_metadata_[high_flat * Fanout + slot];
  }

  /**
   * @brief Choose the better high-node child slot.
   */
  std::size_t better_high_child_slot(std::size_t level,
                                     std::size_t node,
                                     std::size_t left_slot,
                                     std::size_t right_slot) const {
    const MinCandidate left = high_child_min_candidate(level, node, left_slot);
    const MinCandidate right =
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
        high_sparse_min_slots_begin(high_flat) + table_level * Fanout;
    return better_high_child_slot(level, node, table[slot_left],
                                  table[slot_right - span]);
  }

  std::span<const T> values_;
  Compare compare_;
  std::vector<LeafSelector> leaf_selectors_;
  std::vector<Bp512Selector> medium_selectors_;
  std::vector<T> medium_min_values_;
  std::vector<Bp512Selector> high_selectors_;
  std::vector<Index> high_min_positions_;
  std::vector<T> high_min_values_;
  std::vector<HighChildMetadata> high_child_metadata_;
  std::vector<std::uint8_t> high_sparse_min_slots_;
  std::vector<std::size_t> medium_level_offsets_;
  std::vector<std::size_t> high_level_offsets_;
  std::vector<std::size_t> level_sizes_;
  std::vector<std::size_t> level_value_spans_;
  std::vector<std::size_t> level_fanouts_;
  std::size_t high_level_begin_ = std::numeric_limits<std::size_t>::max();
};

}  // namespace pixie::rmq
