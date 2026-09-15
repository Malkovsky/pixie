#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmq.h>
#include <pixie/storage/aligned.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::rmq {

/**
 * @brief Block RMQ with one cache-line mask per block and a sparse table.
 *
 * @details This is a static, non-owning RMQ index over an external value array.
 * Values are partitioned into fixed 496-value blocks. Each block stores one
 * prefix/suffix record bit per value and a 16-bit local-minimum offset in an
 * aligned 64-byte selector. A sparse table over block minima answers full-block
 * ranges.
 *
 * Queries first inspect the minimum of the padded block range covering the
 * request. If that speculative candidate lies inside the requested half-open
 * range, it is final. Otherwise, a multi-block query combines its two partial
 * boundary blocks with the sparse-table answer for fully covered middle
 * blocks. Prefixes and suffixes use the record mask; an unsupported range
 * internal to one block falls back to a linear scan of the original values.
 *
 * Equal values resolve to the smaller original position. The input values are
 * not copied and must remain alive and stable for this object's lifetime.
 *
 * @tparam T Value type in the indexed array.
 * @tparam Compare Strict weak ordering used to choose minima.
 * @tparam Index Unsigned integer type used for stored positions.
 */
template <class T, class Compare = std::less<T>, class Index = std::size_t>
class SimpleBlockRmq : public RmqBase<SimpleBlockRmq<T, Compare, Index>, T> {
 public:
  static_assert(std::is_unsigned_v<Index>,
                "SimpleBlockRmq index type must be unsigned");

  /**
   * @brief Concrete CRTP type used by `RmqBase`.
   */
  using Self = SimpleBlockRmq<T, Compare, Index>;

  /**
   * @brief Sentinel returned for invalid or empty query ranges.
   */
  static constexpr std::size_t npos = RmqBase<Self, T>::npos;

  /**
   * @brief Sentinel reserved in the stored index representation.
   */
  static constexpr Index invalid_index = std::numeric_limits<Index>::max();

  /**
   * @brief Number of original values represented by one block selector.
   */
  static constexpr std::size_t kBlockSize = 496;

  /**
   * @brief Construct an empty RMQ index.
   */
  SimpleBlockRmq() = default;

  /**
   * @brief Copy an RMQ index while preserving its non-owning value span.
   *
   * @param other Index and borrowed value span to copy.
   */
  SimpleBlockRmq(const SimpleBlockRmq& other) = default;

  /**
   * @brief Move an RMQ index and its owned metadata buffers.
   *
   * @details The moved-from index is reset to the empty state.
   *
   * @param other Index whose buffers and borrowed value span are moved.
   */
  SimpleBlockRmq(SimpleBlockRmq&& other) noexcept(
      std::is_nothrow_move_constructible_v<Compare>)
      : compare_(std::move(other.compare_)),
        block_selectors_(std::move(other.block_selectors_)),
        sparse_table_(std::move(other.sparse_table_)) {
    values_ = other.values_;
    other.reset_moved_from();
  }

  /**
   * @brief Copy-assign an RMQ index and its owned metadata.
   *
   * @param other Index and borrowed value span to copy.
   * @return Reference to this index.
   */
  SimpleBlockRmq& operator=(const SimpleBlockRmq& other) = default;

  /**
   * @brief Move-assign an RMQ index and its owned metadata.
   *
   * @details The moved-from index is reset to the empty state.
   *
   * @param other Index whose buffers and borrowed value span are moved.
   * @return Reference to this index.
   */
  SimpleBlockRmq& operator=(SimpleBlockRmq&& other) noexcept(
      std::is_nothrow_move_assignable_v<Compare>) {
    if (this != &other) {
      compare_ = std::move(other.compare_);
      block_selectors_ = std::move(other.block_selectors_);
      sparse_table_ = std::move(other.sparse_table_);
      values_ = other.values_;
      other.reset_moved_from();
    }
    return *this;
  }

  /**
   * @brief Build a simple blocked RMQ index over @p values.
   *
   * @details The values are borrowed rather than copied. Equal values keep the
   * smaller original position as the answer.
   *
   * @param values Values to index. They must outlive the constructed index.
   * @param compare Ordering used to choose minima.
   * @throws std::length_error if `Index` cannot represent every position
   * while reserving `invalid_index` as a sentinel.
   */
  explicit SimpleBlockRmq(std::span<const T> values,
                          Compare compare = Compare())
      : values_(values), compare_(compare) {
    build();
  }

  /**
   * @brief Return the number of indexed values.
   *
   * @return Number of values in the borrowed input span.
   */
  std::size_t size_impl() const { return values_.size(); }

  /**
   * @brief Return the value stored at @p position in the external array.
   *
   * @param position Valid zero-based input position.
   * @return Copy of the value at @p position.
   */
  T value_at_impl(std::size_t position) const { return values_[position]; }

  /**
   * @brief Return the first minimum position in [@p left, @p right).
   *
   * @details The query first tries the sparse-table minimum of the padded
   * covering block range. On a miss it combines the full-block middle and up to
   * two partial boundaries. Invalid or empty ranges return `npos`.
   *
   * @param left First position in the half-open query range.
   * @param right One past the last position in the query range.
   * @return Zero-based position of the first minimum, or `npos`.
   */
  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    if (left + 1 == right) {
      return left;
    }

    const std::size_t covering_block_left = left / kBlockSize;
    const std::size_t covering_block_right = (right - 1) / kBlockSize + 1;
    const std::size_t covering =
        sparse_block_arg_min(covering_block_left, covering_block_right);
    if (left <= covering && covering < right) {
      return covering;
    }

    if (covering_block_left + 1 == covering_block_right) {
      return block_range_arg_min(covering_block_left, left, right);
    }

    const std::size_t first_full_block =
        left / kBlockSize + (left % kBlockSize != 0);
    const std::size_t full_block_right = right / kBlockSize;

    std::size_t answer = npos;
    if (first_full_block < full_block_right) {
      answer = sparse_block_arg_min(first_full_block, full_block_right);
    }

    if (left % kBlockSize != 0) {
      answer = better_position(
          answer, block_range_arg_min(covering_block_left, left,
                                      block_value_end(covering_block_left)));
    }

    if (right % kBlockSize != 0) {
      const std::size_t right_block = right / kBlockSize;
      answer = better_position(
          answer, block_range_arg_min(right_block,
                                      block_value_begin(right_block), right));
    }
    return answer;
  }

  /**
   * @brief Return owned auxiliary memory usage in bytes.
   *
   * @details Counts this object and all selector and sparse-table capacities.
   * The external input values are borrowed and excluded.
   *
   * @return Total owned bytes, including the inline object.
   */
  std::size_t memory_usage_bytes_impl() const {
    std::size_t bytes = sizeof(*this);
    bytes += pixie::vector_capacity_bytes(block_selectors_);
    bytes += pixie::vector_capacity_bytes(sparse_table_);
    for (const TableLevel& level : sparse_table_) {
      bytes += pixie::vector_capacity_bytes(level);
    }
    return bytes;
  }

 private:
  static constexpr std::size_t kMaskWordCount = 8;
  static constexpr std::size_t kOffsetWord = kBlockSize / 64;
  static constexpr std::size_t kOffsetShift = kBlockSize & 63;
  static constexpr std::size_t kOffsetBits = 16;

  /**
   * @brief Prefix/suffix record mask for one value block.
   *
   * @details The first 496 bits identify prefix and suffix record minima. The
   * final 16 bits store the first local-minimum offset.
   */
  class alignas(pixie::kAlignedStorageLineBytes) BlockSelector {
   public:
    /**
     * @brief Construct an empty selector.
     */
    BlockSelector() = default;

    /**
     * @brief Build record masks and the first local-minimum offset.
     *
     * @details Prefix records update only on strict improvement. Suffix
     * records update on improvement or equivalence so smaller positions win
     * ties.
     *
     * @param entry_count Number of active values in the block.
     * @param entry_less Callable returning whether one local slot is strictly
     * better than another.
     * @throws std::length_error if @p entry_count exceeds `kBlockSize`.
     */
    template <class EntryLess>
    void build(std::size_t entry_count, EntryLess entry_less) {
      if (entry_count > kBlockSize) {
        throw std::length_error("SimpleBlockRmq block selector too large");
      }

      words_.fill(0);
      if (entry_count == 0) {
        set_min_offset(0);
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
      set_min_offset(prefix_best);
    }

    /**
     * @brief Return the first minimum slot when the mask can answer directly.
     *
     * @details Full ranges, ranges containing the block minimum, prefixes, and
     * suffixes are supported. Other block-internal ranges return `npos`.
     *
     * @param slot_left First local slot in the half-open query range.
     * @param slot_right One past the final local slot in the query range.
     * @param entry_count Number of active values in the block.
     * @return First local minimum slot, or `npos` for invalid or unsupported
     * ranges.
     */
    std::size_t arg_min(std::size_t slot_left,
                        std::size_t slot_right,
                        std::size_t entry_count) const {
      if (slot_left >= slot_right || slot_right > entry_count ||
          entry_count > kBlockSize) {
        return npos;
      }
      if (slot_left + 1 == slot_right) {
        return slot_left;
      }

      const std::size_t minimum = min_offset();
      if (slot_left <= minimum && minimum < slot_right) {
        return minimum;
      }
      if (slot_left == 0 && slot_right <= minimum) {
        return previous_set_bit_before(slot_right);
      }
      if (slot_right == entry_count && minimum < slot_left) {
        return next_set_bit_at_or_after(slot_left, entry_count);
      }
      return npos;
    }

    /**
     * @brief Return the first minimum offset stored in the selector tail.
     *
     * @return Zero-based local offset of the first block minimum.
     */
    std::size_t min_offset() const {
      return static_cast<std::size_t>((words_[kOffsetWord] >> kOffsetShift) &
                                      kOffsetMask);
    }

   private:
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

    void set_mask_bit(std::size_t slot) {
      words_[slot >> 6] |= std::uint64_t{1} << (slot & 63);
    }

    void set_min_offset(std::size_t offset) {
      words_[kOffsetWord] =
          (words_[kOffsetWord] & kMaskWordMask) |
          ((static_cast<std::uint64_t>(offset) & kOffsetMask) << kOffsetShift);
    }

    std::uint64_t mask_word(std::size_t word) const {
      return word == kOffsetWord ? words_[word] & kMaskWordMask : words_[word];
    }

    std::size_t previous_set_bit_before(std::size_t limit) const {
      if (limit == 0) {
        return npos;
      }

      std::size_t word = (limit - 1) >> 6;
      std::uint64_t bits =
          mask_word(word) & low_bits_mask(((limit - 1) & 63) + 1);
      while (true) {
        if (bits != 0) {
          return word * 64 + 63 - std::countl_zero(bits);
        }
        if (word == 0) {
          return npos;
        }
        --word;
        bits = mask_word(word);
      }
    }

    std::size_t next_set_bit_at_or_after(std::size_t slot,
                                         std::size_t entry_count) const {
      std::size_t word = slot >> 6;
      std::uint64_t bits = mask_word(word) & ~low_bits_mask(slot & 63);
      while (word < kMaskWordCount) {
        if (bits != 0) {
          const std::size_t result = word * 64 + std::countr_zero(bits);
          return result < entry_count ? result : npos;
        }
        ++word;
        bits = word < kMaskWordCount ? mask_word(word) : 0;
      }
      return npos;
    }

    std::array<std::uint64_t, kMaskWordCount> words_{};
  };

  static_assert(kBlockSize + kOffsetBits ==
                kMaskWordCount * std::numeric_limits<std::uint64_t>::digits);
  static_assert(kOffsetShift + kOffsetBits ==
                std::numeric_limits<std::uint64_t>::digits);
  static_assert(sizeof(BlockSelector) == pixie::kAlignedStorageLineBytes);
  static_assert(alignof(BlockSelector) == pixie::kAlignedStorageLineBytes);

  /**
   * @brief Allocator used to start every sparse-table level on a cache line.
   *
   * @tparam Value Element type stored by the allocated container.
   */
  template <class Value>
  class CacheLineAllocator {
   public:
    static_assert(pixie::kAlignedStorageLineBytes >= alignof(Value));

    using value_type = Value;

    CacheLineAllocator() = default;

    template <class Other>
    CacheLineAllocator(const CacheLineAllocator<Other>&) noexcept {}

    [[nodiscard]] Value* allocate(std::size_t count) {
      if (count == 0) {
        return nullptr;
      }
      if (count > std::numeric_limits<std::size_t>::max() / sizeof(Value)) {
        throw std::bad_array_new_length();
      }
      return static_cast<Value*>(
          ::operator new(count * sizeof(Value),
                         std::align_val_t{pixie::kAlignedStorageLineBytes}));
    }

    void deallocate(Value* pointer, std::size_t) noexcept {
      ::operator delete(pointer,
                        std::align_val_t{pixie::kAlignedStorageLineBytes});
    }

    template <class Other>
    bool operator==(const CacheLineAllocator<Other>&) const noexcept {
      return true;
    }

    template <class Other>
    bool operator!=(const CacheLineAllocator<Other>&) const noexcept {
      return false;
    }

    /**
     * @brief Rebind this cache-line allocator to another element type.
     *
     * @tparam Other Replacement element type.
     */
    template <class Other>
    struct rebind {
      /**
       * @brief Cache-line allocator specialized for `Other`.
       */
      using other = CacheLineAllocator<Other>;
    };
  };

  using TableLevel = std::vector<Index, CacheLineAllocator<Index>>;

  /**
   * @brief Build block selectors and all valid sparse-table levels.
   *
   * @details Level zero stores each block's absolute minimum position. Higher
   * levels combine adjacent power-of-two ranges.
   *
   * @throws std::length_error if `Index` cannot represent every input position
   * while reserving `invalid_index` as a sentinel.
   */
  void build() {
    block_selectors_.clear();
    sparse_table_.clear();
    if (values_.empty()) {
      return;
    }
    if (values_.size() > static_cast<std::size_t>(invalid_index)) {
      throw std::length_error("SimpleBlockRmq index type is too small");
    }

    const std::size_t block_count = 1 + (values_.size() - 1) / kBlockSize;
    block_selectors_.resize(block_count);
    sparse_table_.reserve(std::bit_width(block_count));
    sparse_table_.emplace_back(block_count);

    for (std::size_t block = 0; block < block_count; ++block) {
      const std::size_t begin = block_value_begin(block);
      const std::size_t count = block_entry_count(block);
      BlockSelector& selector = block_selectors_[block];
      selector.build(count, [&](std::size_t left, std::size_t right) {
        return compare_(values_[begin + left], values_[begin + right]);
      });
      sparse_table_[0][block] =
          static_cast<Index>(begin + selector.min_offset());
    }

    for (std::size_t span = 2, half = 1; span <= block_count;
         half = span, span <<= 1) {
      const std::size_t previous_level = sparse_table_.size() - 1;
      sparse_table_.emplace_back(block_count - span + 1);
      TableLevel& current = sparse_table_.back();
      const TableLevel& previous = sparse_table_[previous_level];
      for (std::size_t block = 0; block < current.size(); ++block) {
        current[block] = static_cast<Index>(build_better_position(
            static_cast<std::size_t>(previous[block]),
            static_cast<std::size_t>(previous[block + half])));
      }
    }
  }

  /**
   * @brief Return the first minimum over a range of whole blocks.
   *
   * @param block_left First block in the half-open block range.
   * @param block_right One past the final block in the range.
   * @return Absolute first-minimum position, or `npos` for an invalid or empty
   * block range.
   */
  std::size_t sparse_block_arg_min(std::size_t block_left,
                                   std::size_t block_right) const {
    if (block_left >= block_right || sparse_table_.empty() ||
        block_right > block_selectors_.size()) {
      return npos;
    }

    const std::size_t length = block_right - block_left;
    const std::size_t level = std::bit_width(length) - 1;
    const std::size_t span = std::size_t{1} << level;
    const TableLevel& table = sparse_table_[level];
    return better_position(static_cast<std::size_t>(table[block_left]),
                           static_cast<std::size_t>(table[block_right - span]));
  }

  /**
   * @brief Return the first minimum in a range contained in one block.
   *
   * @details The selector handles supported ranges; unsupported interiors are
   * scanned in the borrowed value array.
   *
   * @param block Block containing the complete query range.
   * @param left First absolute position in the query range.
   * @param right One past the final absolute position in the query range.
   * @return Absolute first-minimum position, or `npos` for an empty range.
   */
  std::size_t block_range_arg_min(std::size_t block,
                                  std::size_t left,
                                  std::size_t right) const {
    if (left >= right) {
      return npos;
    }

    const std::size_t begin = block_value_begin(block);
    const std::size_t slot = block_selectors_[block].arg_min(
        left - begin, right - begin, block_entry_count(block));
    return slot == npos ? linear_range_arg_min(left, right) : begin + slot;
  }

  /**
   * @brief Scan an unsupported block-internal range in the original values.
   *
   * @param left First absolute position in the nonempty range.
   * @param right One past the final absolute position in the range.
   * @return Absolute first-minimum position.
   */
  std::size_t linear_range_arg_min(std::size_t left, std::size_t right) const {
    std::size_t best = left;
    for (std::size_t position = left + 1; position < right; ++position) {
      if (compare_(values_[position], values_[best])) {
        best = position;
      }
    }
    return best;
  }

  /**
   * @brief Choose the better candidate from adjacent build ranges.
   *
   * @details The left range precedes the right range, so a tie keeps @p left.
   *
   * @param left Candidate from the left range.
   * @param right Candidate from the right range.
   * @return Absolute position of the selected candidate.
   */
  std::size_t build_better_position(std::size_t left, std::size_t right) const {
    return compare_(values_[right], values_[left]) ? right : left;
  }

  /**
   * @brief Choose the better valid candidate and preserve first-minimum ties.
   *
   * @param left First absolute candidate position, or `npos`.
   * @param right Second absolute candidate position, or `npos`.
   * @return Better valid position, the sole valid position, or `npos` if both
   * candidates are missing.
   */
  std::size_t better_position(std::size_t left, std::size_t right) const {
    if (left == npos) {
      return right;
    }
    if (right == npos) {
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

  /**
   * @brief Return the first original-value position represented by @p block.
   *
   * @param block Zero-based block index.
   * @return Absolute position at the start of @p block.
   */
  static std::size_t block_value_begin(std::size_t block) {
    return block * kBlockSize;
  }

  /**
   * @brief Return one past the final original-value position in @p block.
   *
   * @param block Zero-based block index.
   * @return Absolute exclusive end, clamped to the input size.
   */
  std::size_t block_value_end(std::size_t block) const {
    return std::min(values_.size(), block_value_begin(block) + kBlockSize);
  }

  /**
   * @brief Return the number of active original values in @p block.
   *
   * @param block Zero-based block index.
   * @return Number of represented values, including a partial final block.
   */
  std::size_t block_entry_count(std::size_t block) const {
    return block_value_end(block) - block_value_begin(block);
  }

  /**
   * @brief Restore the coherent empty state after moving owned metadata.
   */
  void reset_moved_from() {
    values_ = {};
    block_selectors_.clear();
    sparse_table_.clear();
  }

  std::span<const T> values_;
  Compare compare_{};
  std::vector<BlockSelector> block_selectors_;
  std::vector<TableLevel> sparse_table_;
};

}  // namespace pixie::rmq
