#pragma once

/**
 * @file power_of_two_sequence.h
 * @brief Experimental packed sequence with configurable storage policies.
 * @details Power-of-two leaves, bounded node retention, cached root bounds,
 * and child ordering are opt-in research configurations. Public defaults are
 * unchanged. Current measurements, limitations, and reproduction commands are
 * retained in benchmarks/sequence_snapshot.md.
 */

#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/experimental/grouped_slot_order256.h>
#include <pixie/experimental/power_of_two_value_block.h>
#include <pixie/experimental/sequence_tree_policy.h>
#include <pixie/experimental/slot_order.h>
#include <pixie/experimental/slot_order16.h>
#include <pixie/permutable_sequence.h>

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>

namespace pixie::experimental {

/**
 * @brief Experimental owning packed sequence with power-of-two leaf capacity.
 * @details Implements PermutableSequenceBase with immutable by-value reads,
 * exclusive ownership, unchanged-value consuming merge, and checked half-open
 * rotations. Uses the same split/join engine as PermutableSequence, with
 * experimental leaf capacity, reserve, and navigation policies. Metadata and
 * alignment padding are additional to the Values*sizeof(T) leaf payload. No
 * contiguous storage or stable reference is promised. Construction consumes
 * input through iter_move; source references are not retained. Mutation
 * allocation failure leaves owners unchanged.
 * @tparam T Unsigned integral type with at most 64 value bits.
 * @tparam Values Power-of-two leaf capacity; payload must be word-sized.
 * @tparam Fanout Even branching factor, at least four.
 * @tparam ReservedNodes Maximum idle internal nodes retained after rotations.
 * Zero disables retention. At >=1 MiB of logical value storage, successful
 * rotations retain up to this many unused nodes for subsequent rotations.
 * Retained bytes are owned, bounded, and included in memory_usage().
 * @tparam CacheRegularPrefix Cache the leading regular-region size in one
 * owner word, refreshed after each successful mutation. Requires power-of-two
 * fanout. Pure radix reads in that region avoid loading root size metadata.
 * @tparam PermutedLeaves Map children through a slot permutation immediately
 * above leaves; requires Fanout == 16, 32, 64, or 256. Fanout 16 uses packed
 * nibbles, 32/64 use byte indices with primitive-local SIMD dispatch, and 256
 * uses one node containing sixteen groups of sixteen child slots. Higher levels
 * retain direct addressing. Fresh/repacked nodes also address children directly
 * until their first mapped rotation constructs the policy state without
 * allocation. All internal allocations reserve space for ordering state so the
 * existing allocation pool can recycle nodes between levels.
 */
template <class T = std::uint64_t,
          std::size_t Values = 256,
          std::size_t Fanout = 32,
          std::size_t ReservedNodes = 0,
          bool CacheRegularPrefix = false,
          bool PermutedLeaves = false>
class PowerOfTwoPermutableSequence
    : public PermutableSequenceBase<
          PowerOfTwoPermutableSequence<T,
                                       Values,
                                       Fanout,
                                       ReservedNodes,
                                       CacheRegularPrefix,
                                       PermutedLeaves>,
          T,
          T> {
  using Block = PowerOfTwoValueBlock<T, Values>;
  static_assert(!PermutedLeaves || Fanout == 16 || Fanout == 32 ||
                Fanout == 64 || Fanout == 256);
  using LeafOrder =
      std::conditional_t<Fanout == 16,
                         SlotOrder16,
                         std::conditional_t<Fanout == 256,
                                            GroupedSlotOrder256,
                                            ByteSlotOrder<Fanout>>>;
  using Tree = detail::sequence::SequenceTree<
      Block,
      Fanout,
      LengthLayout::cumulative,
      false,
      SequenceTreePolicy<ReservedNodes,
                         std::conditional_t<PermutedLeaves, LeafOrder, void>>>;
  friend class PermutableSequenceBase<PowerOfTwoPermutableSequence, T, T>;
  static_assert(!CacheRegularPrefix || std::has_single_bit(Fanout));
  struct NoCache {};
  struct RegularCache {
    std::size_t end = 0;
    RegularCache() noexcept = default;
    RegularCache(RegularCache&& other) noexcept
        : end(std::exchange(other.end, 0)) {}
    RegularCache& operator=(RegularCache&& other) noexcept {
      if (this != &other) {
        end = std::exchange(other.end, 0);
      }
      return *this;
    }
  };
  Tree tree_;
  [[no_unique_address]]
  std::conditional_t<CacheRegularPrefix, RegularCache, NoCache> regular_;
  void refresh_regular() noexcept {
    if constexpr (CacheRegularPrefix) {
      regular_.end = tree_.regular_prefix_size();
    }
  }

  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static PowerOfTwoPermutableSequence from_range_impl(Range&& range) {
    auto it = std::ranges::begin(range);
    const auto end = std::ranges::end(range);
    auto next = [&] {
      std::array<T, Values> values;
      std::size_t count = 0;
      while (count != Values && it != end) {
        values[count++] = T(std::ranges::iter_move(it));
        ++it;
      }
      return Block(std::span<const T>(values.data(), count));
    };
    // A consuming input cursor: dereferencing does not advance the source.
    struct Blocks {
      decltype(next)& read;
      Block current;
      struct Iterator {
        using value_type [[maybe_unused]] = Block;
        using difference_type [[maybe_unused]] = std::ptrdiff_t;
        using iterator_concept [[maybe_unused]] = std::input_iterator_tag;
        Blocks* source;
        Block& operator*() const noexcept { return source->current; }
        Iterator& operator++() {
          source->current = source->read();
          return *this;
        }
        void operator++(int) { ++*this; }
        bool operator==(std::default_sentinel_t) const noexcept {
          return source->current.empty();
        }
      };
      Iterator begin() {
        current = read();
        return {this};
      }
      std::default_sentinel_t end() const noexcept { return {}; }
    } blocks{next, {}};
    PowerOfTwoPermutableSequence result;
    result.tree_ = Tree::from_blocks(blocks);
    result.refresh_regular();
    return result;
  }

  std::size_t size_impl() const noexcept { return tree_.size(); }
  T value_at_impl(std::size_t i) const {
    if constexpr (CacheRegularPrefix) {
      if (i < regular_.end) {
        return tree_.read_regular_unchecked(i);
      }
    }
    return tree_[i];
  }
  void insert_at_impl(std::size_t position, T value) {
    tree_.insert_at(position, value);
    refresh_regular();
  }
  void rotate_left_impl(std::size_t left,
                        std::size_t right,
                        std::size_t distance) {
    tree_.rotate_left(left, right, distance);
    refresh_regular();
  }
  void merge_impl(PowerOfTwoPermutableSequence& donor) {
    tree_.merge(donor.tree_);
    refresh_regular();
    donor.refresh_regular();
  }
  std::size_t memory_usage_bytes_impl() const noexcept {
    return memory_usage().total_bytes;
  }

 public:
  /** @brief Extra owner bytes used by the optional regular-prefix cache. */
  static constexpr std::size_t regular_cache_bytes =
      CacheRegularPrefix ? sizeof(std::size_t) : 0;
  /** @brief Resolved storage strategy; elements are packed and read by value.
   */
  static constexpr ElementStorage storage = ElementStorage::packed;
  /** @brief Construct canonical empty without allocating. */
  PowerOfTwoPermutableSequence() noexcept = default;
  /** @brief Exclusive owners cannot be copied. */
  PowerOfTwoPermutableSequence(const PowerOfTwoPermutableSequence&) = delete;
  /** @brief Exclusive owners cannot be copy-assigned. */
  PowerOfTwoPermutableSequence& operator=(const PowerOfTwoPermutableSequence&) =
      delete;
  /** @brief Transfer ownership, leaving the source canonical empty. */
  PowerOfTwoPermutableSequence(PowerOfTwoPermutableSequence&&) noexcept =
      default;
  /** @brief Replace ownership; self-move is a no-op, source becomes empty. */
  PowerOfTwoPermutableSequence& operator=(
      PowerOfTwoPermutableSequence&&) noexcept = default;
  /**
   * @brief Account for this object and all live requested tree storage.
   * @details Includes unused leaf/node capacity and alignment padding;
   * excludes allocator bookkeeping. Traverses the tree without allocating or
   * throwing.
   * @return Tree memory breakdown, with total_bytes including this owner.
   */
  auto memory_usage() const noexcept {
    auto result = tree_.memory_usage();
    const auto extra = sizeof(*this) - sizeof(Tree);
    result.total_bytes += std::min(extra, SIZE_MAX - result.total_bytes);
    return result;
  }

#ifdef PIXIE_SEQUENCE_TREE_TESTING
  /** @brief Test-only tree type for structural and allocation-failure probes.
   */
  using test_tree_type = Tree;
  /** @brief Borrow a read-only tree view; lifetime is bounded by this owner.
   */
  const Tree& test_tree() const noexcept { return tree_; }
#endif
};

}  // namespace pixie::experimental
