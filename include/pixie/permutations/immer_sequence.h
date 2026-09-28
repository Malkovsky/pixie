#pragma once

/**
 * @file immer_sequence.h
 * @brief Optional by-value permutable sequence backed by pinned Immer.
 */

#include <pixie/permutable_sequence.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <immer/flex_vector.hpp>
#include <immer/flex_vector_transient.hpp>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace pixie {

/**
 * @brief Move-only, by-value adapter for immer::flex_vector<T>.
 * @details Implements PermutableSequenceBase, including unchanged-value
 * consuming merge. No snapshots or references escape the adapter. Mutations
 * keep the original roots alive until success to preserve the facade's strong
 * exception guarantees; temporary slices use Immer's consuming operations.
 * Build with PIXIE_THIRD_PARTY_BACKENDS. This is not the stable-reference
 * indirect storage policy and does not implement permutation rebasing.
 * @tparam T Copy-constructible value type with at most max_align_t alignment.
 * Reads and persistent updates may copy values. Move-only payloads and
 * stable-reference reads are not supported.
 * @tparam HeapPolicy Immer heap policy; defaults to its standard pooled heap.
 * Changing the heap does not change the reference-counting or value policy.
 */
template <std::copy_constructible T = std::uint64_t,
          class HeapPolicy = immer::default_heap_policy>
class ImmerPermutableSequence
    : public PermutableSequenceBase<ImmerPermutableSequence<T, HeapPolicy>,
                                    T,
                                    T> {
  using Base = PermutableSequenceBase<ImmerPermutableSequence, T, T>;
  friend Base;
  using Policy = immer::memory_policy<HeapPolicy,
                                      immer::default_refcount_policy,
                                      immer::default_lock_policy>;
  using Vector = immer::flex_vector<T, Policy>;
  using Tree = std::remove_cvref_t<decltype(std::declval<Vector>().impl())>;
  using Node = typename Tree::node_t;
  static_assert(Node::keep_headroom);
  static_assert(alignof(T) <= alignof(std::max_align_t));

 public:
  /** @brief Maximum number of child pointers per internal node. */
  static constexpr std::size_t fanout = std::size_t{1} << Vector::bits;
  /** @brief Maximum number of values per leaf. */
  static constexpr std::size_t leaf_capacity = std::size_t{1}
                                               << Vector::bits_leaf;

  /** @brief Construct canonical empty without heap allocation. */
  ImmerPermutableSequence() noexcept = default;
  /** @brief Ownership is exclusive; copying is disabled. */
  ImmerPermutableSequence(const ImmerPermutableSequence&) = delete;
  /** @brief Ownership is exclusive; copying is disabled. */
  ImmerPermutableSequence& operator=(const ImmerPermutableSequence&) = delete;
  /** @brief Transfer ownership, leaving other canonical empty. */
  ImmerPermutableSequence(ImmerPermutableSequence&& other) noexcept
      : values_(std::exchange(other.values_, Vector{})) {}
  /** @brief Transfer ownership; self-move is a no-op, other becomes empty. */
  ImmerPermutableSequence& operator=(ImmerPermutableSequence&& other) noexcept {
    if (this != &other) {
      values_ = std::exchange(other.values_, Vector{});
    }
    return *this;
  }

  /** @brief Requested live tree bytes, excluding heap free lists and RSS. */
  struct MemoryUsage {
    /** @brief Number of live leaf allocations, including the tail. */
    std::size_t blocks = 0;
    /** @brief Number of live internal nodes. */
    std::size_t nodes = 0;
    /** @brief Leaf allocation bytes including unused capacity. */
    std::size_t block_bytes = 0;
    /** @brief Internal node and relaxed-size-table allocation bytes. */
    std::size_t node_bytes = 0;
    /** @brief Facade plus block_bytes plus node_bytes; saturates at SIZE_MAX.
     */
    std::size_t total_bytes = sizeof(ImmerPermutableSequence);
  };

  /**
   * @brief Account for live storage in linear time, without allocations.
   * @details Uses the pinned Immer node layout and its fixed-capacity
   * reference-counted allocations. Excludes static empty nodes, allocator
   * pools/bookkeeping, RSS and temporary persistent snapshots. Returned owners
   * do not share payload subtrees with other adapter instances.
   * @return Disjoint allocation counts and bytes, including this facade once.
   */
  MemoryUsage memory_usage() const noexcept {
    MemoryUsage result;
    if (!values_.empty()) {
      const auto& tree = values_.impl();
      const auto tail_offset = tree.tail_offset();
      if (tail_offset != 0) {
        count_tree(tree.root, tree.shift, tail_offset, result);
      }
      ++result.blocks;
      add(result.block_bytes, Node::max_sizeof_leaf);
    }
    add(result.total_bytes, result.block_bytes);
    add(result.total_bytes, result.node_bytes);
    return result;
  }

 private:
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static ImmerPermutableSequence from_range_impl(Range&& range) {
    auto transient = Vector{}.transient();
    auto it = std::ranges::begin(range);
    const auto end = std::ranges::end(range);
    for (; it != end; ++it) {
      if (transient.size() == Vector::max_size()) {
        throw std::length_error("ImmerPermutableSequence: size limit");
      }
      transient.push_back(T(std::ranges::iter_move(it)));
    }
    ImmerPermutableSequence result;
    result.values_ = std::move(transient).persistent();
    return result;
  }

  std::size_t size_impl() const noexcept { return values_.size(); }
  T value_at_impl(std::size_t position) const {
    if (position >= values_.size()) {
      throw std::out_of_range("ImmerPermutableSequence: read position");
    }
    return values_[position];
  }
  void insert_at_impl(std::size_t position, T value) {
    if (position > values_.size()) {
      throw std::out_of_range("ImmerPermutableSequence: insert position");
    }
    check_growth(1);
    // Keep the original root for rollback; consume the temporary slices.
    values_ = values_.insert(position, std::move(value));
  }
  void rotate_left_impl(std::size_t left,
                        std::size_t right,
                        std::size_t distance) {
    if (left > right || right > values_.size()) {
      throw std::out_of_range("ImmerPermutableSequence: rotation range");
    }
    if (left == right || (distance %= right - left) == 0) {
      return;
    }
    auto prefix = values_.take(left);
    auto first = values_.drop(left).take(distance);
    auto second = values_.drop(left + distance).take(right - left - distance);
    auto suffix = values_.drop(right);
    auto result = std::move(prefix) + std::move(second) + std::move(first) +
                  std::move(suffix);
    values_ = std::move(result);
  }
  void merge_impl(ImmerPermutableSequence& donor) {
    if (this == &donor || donor.empty()) {
      return;
    }
    check_growth(donor.size());
    if (this->empty()) {
      *this = std::move(donor);
      return;
    }
    auto result = values_ + donor.values_;
    values_ = std::move(result);
    donor.values_ = Vector{};
  }
  std::size_t memory_usage_bytes_impl() const noexcept {
    return memory_usage().total_bytes;
  }
  void check_growth(std::size_t count) const {
    if (count > Vector::max_size() - values_.size()) {
      throw std::length_error("ImmerPermutableSequence: size limit");
    }
  }
  static void add(std::size_t& total, std::size_t bytes) noexcept {
    const auto maximum = std::numeric_limits<std::size_t>::max();
    total = bytes > maximum - total ? maximum : total + bytes;
  }
  static void count_tree(Node* node,
                         unsigned shift,
                         std::size_t size,
                         MemoryUsage& result) noexcept {
    if (shift < Vector::bits_leaf) {
      ++result.blocks;
      add(result.block_bytes, Node::max_sizeof_leaf);
      return;
    }
    ++result.nodes;
    const auto* relaxed = node->relaxed();
    add(result.node_bytes,
        relaxed ? Node::max_sizeof_inner_r : Node::max_sizeof_inner);
    if (relaxed && !Node::embed_relaxed) {
      add(result.node_bytes, Node::max_sizeof_relaxed);
    }
    const auto child_capacity = std::size_t{1} << shift;
    const auto count =
        relaxed ? relaxed->d.count : 1 + (size - 1) / child_capacity;
    std::size_t prefix = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const auto end = relaxed ? relaxed->d.sizes[i]
                               : std::min(size, prefix + child_capacity);
      if (shift == Vector::bits_leaf) {
        ++result.blocks;
        add(result.block_bytes, Node::max_sizeof_leaf);
      } else {
        count_tree(node->inner()[i], shift - Vector::bits, end - prefix,
                   result);
      }
      prefix = end;
    }
  }

  Vector values_;
};
}  // namespace pixie
