#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <utility>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

// Default child addressing. Policy headers retain the count and regular
// prefix in one word; state-free policies add no per-node storage.
template <std::size_t Fanout>
struct DirectChildOrder {
  static constexpr bool supports_bias = true;
  struct Header {
    std::size_t count : std::bit_width(Fanout) = 0;
    std::size_t regular_prefix : std::bit_width(Fanout) = 0;
  };
  template <class Node>
  static std::size_t slot(const Node&, std::size_t i, std::size_t) noexcept {
    return i;
  }
  template <class Node>
  static void reset(Node&) noexcept {}
  template <class Tree, class Node>
  static std::size_t rotate(Node& node,
                            std::size_t begin,
                            std::size_t middle,
                            std::size_t end,
                            std::size_t) noexcept {
    std::rotate(node.children.begin() + begin, node.children.begin() + middle,
                node.children.begin() + end);
    return end - begin;
  }
  template <class Node>
  static bool valid(const Node&, std::size_t) noexcept {
    return true;
  }
};

// The default owner retains no spare allocations between operations.
template <class Node>
struct NoNodeReserve {
  NoNodeReserve* for_rotation(std::size_t, std::size_t) noexcept {
    return nullptr;
  }
  void trim(std::size_t, std::size_t) noexcept {}
  void clear() noexcept {}
  std::size_t node_count() const noexcept { return 0; }
};

// Transaction-local preflight storage. Dismantled nodes are recycled through
// their unoccupied first pointer slot; occupied children belong to Owners.
// Allocator supplies instrumented allocation and node-header initialization.
template <class Node, class Owner, class Allocator>
struct SequenceNodePool {
  Node* nodes = nullptr;
  std::array<Owner, 3> leaves;
  std::size_t leaf_count = 0;

  explicit SequenceNodePool(NoNodeReserve<Node>* = nullptr) noexcept {}
  SequenceNodePool(const SequenceNodePool&) = delete;
  SequenceNodePool& operator=(const SequenceNodePool&) = delete;
  ~SequenceNodePool() {
    while (nodes) {
      auto* next = static_cast<Node*>(nodes->children[0]);
      Allocator::dispose(nodes);
      nodes = next;
    }
  }
  void commit() noexcept {}
  void reserve(std::size_t node_count, std::size_t leaf_spares) {
    for (std::size_t i = 0; i < node_count; ++i) {
      recycle(Allocator::node());
    }
    for (; leaf_count < leaf_spares; ++leaf_count) {
      leaves[leaf_count] = Allocator::leaf();
    }
  }
  void recycle(Node* node) noexcept {
    Allocator::reset(node);
    node->children[0] = nodes;
    nodes = node;
  }
  Node* node() noexcept {
    assert(nodes && nodes->count == 0);
    auto* result = nodes;
    nodes = static_cast<Node*>(nodes->children[0]);
    result->children[0] = nullptr;
    return result;
  }
  Owner leaf() noexcept {
    assert(leaf_count);
    return std::move(leaves[--leaf_count]);
  }
};

// Internal extension boundary: production uses direct pointers and strictly
// transaction-local allocation. Other policies must preserve noexcept commit,
// ownership, child occupancy, and preflight rollback contracts.
struct SequenceTreePolicy {
  template <std::size_t Fanout>
  using Order = DirectChildOrder<Fanout>;
  template <class Node, class Allocator>
  using Reserve = NoNodeReserve<Node>;
  template <class Node, class Owner, class Allocator>
  using Pool = SequenceNodePool<Node, Owner, Allocator>;
  struct TestCounters {};
};

}  // namespace pixie::detail::sequence
/// @endcond
