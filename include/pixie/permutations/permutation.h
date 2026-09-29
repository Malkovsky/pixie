#pragma once

#include <pixie/detail/sequence/packed_value_block.h>
#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/permutation.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace pixie {

/**
 * @brief Exclusively owned permutation of the integers in [0, size()).
 * @details Construction is identity-only. Rotations preserve the bijection;
 * consuming merge appends donor entries plus the old receiver size. Fixed-width
 * indices never widen. Reads return values, not references. Access, rotation,
 * and merge use logarithmic tree work plus bounded local re-encoding; identity
 * construction and explicit memory accounting are linear. No public arbitrary
 * import, split, element assignment, or value addition is provided. Allocation
 * failure leaves existing containers unchanged, including pending lazy index
 * biases. Single-leaf and eligible complete-child rotations are nonallocating;
 * general rotations use transactional split/join with bounded leaf work.
 * Default leaves use 256-byte blocks with 1920 payload bits, wrapped to 320
 * bytes for lazy tags/alignment; F8 tagged nodes use 192 bytes on 64-bit hosts
 * versus 128 untagged. Metadata is not promised strictly succinct. There is
 * no whole-container contiguous buffer. Internal detail types are unsupported.
 * @tparam Index Unsigned integer, excluding bool, with at most 64 value bits.
 * @tparam StorageBits Complete packed block budget, excluding tagged wrappers.
 * @tparam Fanout Even tree fanout of at least four.
 * @tparam Layout Cumulative or individual child measures.
 */
template <class Index = std::uint64_t,
          std::size_t StorageBits = 2048,
          std::size_t Fanout = 8,
          LengthLayout Layout = LengthLayout::cumulative>
class Permutation
    : public PermutationBase<Permutation<Index, StorageBits, Fanout, Layout>,
                             Index> {
  friend class PermutationBase<Permutation, Index>;
  using block_type = detail::sequence::PackedValueBlock<Index, StorageBits>;
  using tree_type =
      detail::sequence::SequenceTree<block_type, Fanout, Layout, true>;
  static_assert(std::is_integral_v<Index> && std::is_unsigned_v<Index> &&
                !std::is_same_v<Index, bool> &&
                std::numeric_limits<Index>::digits <= 64);

 public:
  /** @brief Construct empty. @details No allocation or implicit identity data.
   */
  Permutation() noexcept = default;
  /** @brief Disallow copying. @details Ownership is exclusive. */
  Permutation(const Permutation&) = delete;
  /** @brief Disallow copy assignment. @details Ownership is exclusive. */
  Permutation& operator=(const Permutation&) = delete;
  /**
   * @brief Transfer ownership without allocation.
   * @details Source becomes canonical empty; no stored index is rewritten.
   * @param other Source whose ownership is consumed.
   */
  Permutation(Permutation&& other) noexcept = default;
  /**
   * @brief Replace ownership without throwing.
   * @details Source becomes canonical empty; self-move is a no-op.
   * @param other Source whose ownership replaces this owner's storage.
   * @return This owner after replacement.
   */
  Permutation& operator=(Permutation&& other) noexcept = default;

  // Stream identity blocks without a full temporary index vector.
  /**
   * @brief Find the first position where pred(logical_index) is true.
   * @details Uses the tree's guided descent: at each node children are
   * binary-searched by their rightmost stored index, then the search descends
   * into the single child that can contain the first match. The predicate
   * receives the decoded logical index (pending biases applied). Caller must
   * guarantee monotonicity: all false results precede all true results.
   * @param pred Predicate on Index returning true at or past the search target.
   * @return Position in [0, size()]; size() means pred is false everywhere.
   */
  template <class Pred>
  std::size_t lower_bound(Pred pred) const {
    return tree_.lower_bound(pred);
  }
  /**
   * @brief Requested live storage breakdown, excluding allocator bookkeeping.
   * @details The five byte categories payload_capacity_bytes,
   * block_metadata_bytes, ordering_tree_bytes, tag_padding_bytes, and
   * facade_bytes are disjoint and sum to total_bytes. Payload includes unused
   * capacity. Tag/padding bytes are incremental over the same untagged layout.
   * Counts exclude preflight peaks and RSS; arithmetic saturates at SIZE_MAX.
   */
  struct MemoryUsage {
    /** @brief Leaf count. @details Counts live leaf allocations only. */
    std::size_t blocks = 0;
    /** @brief Internal node count. @details Excludes leaf allocations. */
    std::size_t nodes = 0;
    /** @brief Packed payload bytes. @details Includes unused block capacity. */
    std::size_t payload_capacity_bytes = 0;
    /** @brief Local bookkeeping bytes. @details Excludes payload and tags. */
    std::size_t block_metadata_bytes = 0;
    /** @brief Untagged node bytes. @details Includes their alignment padding.
     */
    std::size_t ordering_tree_bytes = 0;
    /** @brief Incremental tag bytes. @details Includes additional padding. */
    std::size_t tag_padding_bytes = 0;
    /** @brief Inline owner bytes. @details Includes the embedded tree handle.
     */
    std::size_t facade_bytes = sizeof(Permutation);
    /** @brief Total live bytes. @details Saturated sum of disjoint categories.
     */
    std::size_t total_bytes = sizeof(Permutation);
  };
  /**
   * @brief Explicit linear-time live allocation accounting.
   * @details Never called implicitly during queries or mutation; visits each
   * allocation once, with no allocation of its own.
   * @return Disjoint storage categories and total requested live bytes.
   */
  MemoryUsage memory_usage() const noexcept {
    using Untagged = detail::sequence::SequenceTree<block_type, Fanout, Layout>;
    const auto tree = tree_.memory_usage();
    MemoryUsage result;
    result.blocks = tree.blocks;
    result.nodes = tree.nodes;
    result.payload_capacity_bytes = tree_type::saturated_multiply(
        tree.blocks, block_type{}.payload_capacity_bytes());
    result.block_metadata_bytes = tree_type::saturated_multiply(
        tree.blocks, block_type{}.metadata_bytes());
    result.ordering_tree_bytes =
        tree_type::saturated_multiply(tree.nodes, Untagged::node_storage_bytes);
    result.tag_padding_bytes = tree_type::saturated_add(
        tree_type::saturated_multiply(
            tree.blocks, tree_type::block_storage_bytes - sizeof(block_type)),
        tree_type::saturated_multiply(
            tree.nodes,
            tree_type::node_storage_bytes - Untagged::node_storage_bytes));
    result.total_bytes = tree_type::saturated_add(
        sizeof(*this),
        tree_type::saturated_add(tree.block_bytes, tree.node_bytes));
    return result;
  }
#ifdef PIXIE_SEQUENCE_TREE_TESTING
  /**
   * @brief Inspect the tree in test builds without exposing mutable ownership.
   * @details Provides structural validation and stable identities, not
   * rebasing.
   * @return Read-only tree reference borrowing this owner's lifetime.
   */
  const tree_type& test_tree() const noexcept { return tree_; }
#endif

 private:
  static Permutation identity_impl(std::size_t n) {
    check_domain(n);
    const auto count =
        n / block_type::capacity + (n % block_type::capacity != 0);
    auto blocks =
        std::views::iota(std::size_t{0}, count) |
        std::views::transform([n](std::size_t block) {
          std::array<Index, block_type::capacity> fields;
          const auto start = block * block_type::capacity;
          const auto length = std::min(block_type::capacity, n - start);
          for (std::size_t i = 0; i < length; ++i) {
            fields[i] = static_cast<Index>(start + i);
          }
          return block_type(std::span<const Index>(fields.data(), length));
        });
    Permutation result;
    result.tree_ = tree_type::from_blocks(blocks);
    return result;
  }
  std::size_t size_impl() const noexcept { return tree_.size(); }
  Index value_at_impl(std::size_t position) const { return tree_[position]; }
  /**
   * @brief Insert the next identity value at a zero-based position.
   * @details Inserts the value size() at position, shifting later entries
   * right. Uses the tree's atomic insert: fast path is leaf-local with zero
   * allocation; slow path is a single split/concatenate transaction. On
   * success the result is a valid permutation of [0,size()+1).
   * @param position Insertion position in [0,size()].
   * @throws std::out_of_range If position > size().
   * @throws std::length_error If the new size exceeds the index domain.
   * @throws std::bad_alloc On preflight failure; contents remain unchanged.
   */
  void insert_at_impl(std::size_t position) {
    if (position > tree_.size()) {
      throw std::out_of_range("Permutation: insert position");
    }
    check_domain(tree_.size() + 1);
    tree_.insert_at(position, static_cast<Index>(tree_.size()));
  }
  void rotate_left_impl(std::size_t left,
                        std::size_t right,
                        std::size_t distance) {
    tree_.rotate_left(left, right, distance);
  }
  // Validate the resulting index domain before rebasing the donor.
  void merge_impl(Permutation& donor) {
    if (this == &donor || donor.empty()) {
      return;
    }
    if (donor.size() > std::numeric_limits<std::size_t>::max() - this->size()) {
      throw std::length_error("Permutation: concatenation size");
    }
    check_domain(this->size() + donor.size());
    tree_.merge_rebased(donor.tree_);
  }
  std::size_t memory_usage_bytes_impl() const noexcept {
    return memory_usage().total_bytes;
  }
  static void check_domain(std::size_t n) {
    if (n != 0 && n - 1 > std::numeric_limits<Index>::max()) {
      throw std::length_error("Permutation: index domain");
    }
  }
  tree_type tree_;
};

}  // namespace pixie
