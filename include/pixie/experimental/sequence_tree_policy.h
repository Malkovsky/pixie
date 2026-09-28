#pragma once

#include <pixie/detail/sequence/tree_policy.h>

#include <bit>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::experimental {

// Lazy ordering is confined to leaf parents. Its activation flag shares the
// node's count word; inactive union storage is never read. Uniform node sizes
// allow the existing allocation pool to recycle across tree levels.
template <std::size_t Fanout, class Order>
class LeafChildOrder {
  static_assert(Order::capacity == Fanout);
  static_assert(std::is_nothrow_constructible_v<Order, std::size_t> &&
                std::is_trivially_destructible_v<Order>);
  union Storage {
    Order value;
    Storage() noexcept {}
  } storage_;

 public:
  static constexpr bool supports_bias = false;
  struct Header {
    std::size_t count : std::bit_width(Fanout) = 0;
    std::size_t regular_prefix : std::bit_width(Fanout) = 0;
    std::size_t mapped : 1 = 0;
  };
  template <class Node>
  std::size_t slot(const Node& node,
                   std::size_t i,
                   std::size_t height) const noexcept {
    return height == 1 && node.mapped ? storage_.value[i] : i;
  }
  template <class Node>
  void reset(Node& node) noexcept {
    node.mapped = false;
  }
  template <class Tree, class Node>
  std::size_t rotate(Node& node,
                     std::size_t begin,
                     std::size_t middle,
                     std::size_t end,
                     std::size_t height) noexcept {
    if (height != 1) {
      return detail::sequence::DirectChildOrder<Fanout>::template rotate<Tree>(
          node, begin, middle, end, height);
    }
    if (!node.mapped) {
      std::construct_at(&storage_.value, std::size_t(node.count));
      node.mapped = true;
#ifdef PIXIE_SEQUENCE_TREE_TESTING
      ++Tree::test_counters.order_materializations;
#endif
    }
    return storage_.value.rotate_children(node.children, begin, end,
                                          middle - begin);
  }
  template <class Node>
  bool valid(const Node& node, std::size_t height) const noexcept {
    return !node.mapped ||
           (height == 1 && storage_.value.size() == node.count &&
            storage_.value.valid());
  }
};

// Experimental bounded rotation reserve. Small owners retain nothing. Failed
// preflight restores the old reserve footprint; only successful commits can
// retain newly allocated spares up to Limit. No payload leaf is cached.
template <class Node, class Allocator, std::size_t Limit>
struct BoundedNodeReserve {
  Node* nodes = nullptr;
  std::size_t count = 0;
  BoundedNodeReserve() = default;
  BoundedNodeReserve(const BoundedNodeReserve&) = delete;
  BoundedNodeReserve& operator=(const BoundedNodeReserve&) = delete;
  BoundedNodeReserve(BoundedNodeReserve&& other) noexcept
      : nodes(std::exchange(other.nodes, nullptr)),
        count(std::exchange(other.count, 0)) {}
  BoundedNodeReserve& operator=(BoundedNodeReserve&& other) noexcept {
    if (this != &other) {
      clear();
      nodes = std::exchange(other.nodes, nullptr);
      count = std::exchange(other.count, 0);
    }
    return *this;
  }
  ~BoundedNodeReserve() { clear(); }
  void clear() noexcept {
    while (nodes) {
      auto* next = static_cast<Node*>(nodes->children[0]);
      Allocator::dispose(nodes);
      nodes = next;
    }
    count = 0;
  }
  BoundedNodeReserve* for_rotation(std::size_t size,
                                   std::size_t value_bytes) noexcept {
    return size >= (std::size_t{1} << 20) / value_bytes ? this : nullptr;
  }
  void trim(std::size_t size, std::size_t value_bytes) noexcept {
    if (!for_rotation(size, value_bytes)) {
      clear();
    }
  }
  std::size_t node_count() const noexcept { return count; }
};

template <class Node, class Owner, class Allocator, std::size_t Limit>
class RetainedNodePool
    : public detail::sequence::SequenceNodePool<Node, Owner, Allocator> {
  using Base = detail::sequence::SequenceNodePool<Node, Owner, Allocator>;
  using Reserve = BoundedNodeReserve<Node, Allocator, Limit>;
  Reserve* retained_ = nullptr;
  std::size_t borrowed_ = 0;
  std::size_t retain_limit_ = 0;

 public:
  explicit RetainedNodePool(Reserve* reserve = nullptr) noexcept
      : retained_(reserve) {
    if (retained_) {
      this->nodes = std::exchange(retained_->nodes, nullptr);
      borrowed_ = std::exchange(retained_->count, 0);
      retain_limit_ = borrowed_;
    }
  }
  void commit() noexcept { retain_limit_ = Limit; }
  ~RetainedNodePool() {
    if (retained_) {
      while (this->nodes && retained_->count < retain_limit_) {
        auto* next = static_cast<Node*>(this->nodes->children[0]);
        this->nodes->children[0] = retained_->nodes;
        retained_->nodes = this->nodes;
        ++retained_->count;
        this->nodes = next;
      }
    }
  }
  void reserve(std::size_t count, std::size_t leaf_spares) {
    const auto available = std::exchange(borrowed_, 0);
    Base::reserve(count > available ? count - available : 0, leaf_spares);
  }
};

template <std::size_t ReservedNodes, class LeafOrder = void>
struct SequenceTreePolicy {
  template <std::size_t Fanout>
  using Order = std::conditional_t<std::is_void_v<LeafOrder>,
                                   detail::sequence::DirectChildOrder<Fanout>,
                                   LeafChildOrder<Fanout, LeafOrder>>;
  template <class Node, class Allocator>
  using Reserve =
      std::conditional_t<ReservedNodes == 0,
                         detail::sequence::NoNodeReserve<Node>,
                         BoundedNodeReserve<Node, Allocator, ReservedNodes>>;
  template <class Node, class Owner, class Allocator>
  using Pool = std::conditional_t<
      ReservedNodes == 0,
      detail::sequence::SequenceNodePool<Node, Owner, Allocator>,
      RetainedNodePool<Node, Owner, Allocator, ReservedNodes>>;
  struct TestCounters {
    std::size_t order_materializations = 0;
  };
};

}  // namespace pixie::experimental
/// @endcond
