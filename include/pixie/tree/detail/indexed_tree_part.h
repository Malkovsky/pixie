#pragma once

#include <pixie/tree/tree_part.h>

#include <cstddef>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::tree::detail {

// Numbering follows the producer's order. Cached child runs and LOUDS scans
// have consecutive local child IDs. Preorder scans need explicit child IDs.
template <class Tree, bool ContiguousChildren>
class IndexedTreePart
    : public TreePart<IndexedTreePart<Tree, ContiguousChildren>, Tree> {
  using Base = TreePart<IndexedTreePart, Tree>;
  using SourceNode = typename Tree::Node;
  struct Unused {};
  static constexpr std::size_t kIncomplete =
      std::numeric_limits<std::size_t>::max() / 2 + 1;

 public:
  using Node = typename Base::Node;

  struct Record {
    Record(SourceNode source,
           std::size_t first,
           std::size_t count,
           bool complete)
        : source(source),
          first(first),
          count_and_status(complete ? count : count | kIncomplete) {}

    std::size_t count() const { return count_and_status & ~kIncomplete; }
    bool complete() const { return (count_and_status & kIncomplete) == 0; }
    void append_child() { ++count_and_status; }

    SourceNode source;
    std::size_t first;
    Node parent = 0;
    [[no_unique_address]] std::conditional_t<ContiguousChildren, Unused, Node>
        sibling_position{};

   private:
    // The count is bounded by the number of allocated records.
    std::size_t count_and_status;
  };

  using ChildStorage =
      std::conditional_t<ContiguousChildren, Unused, std::vector<Node>>;
  static constexpr std::size_t bytes_per_node =
      sizeof(Record) + (ContiguousChildren ? 0 : sizeof(Node));

  IndexedTreePart(const Tree& tree,
                  std::vector<Record> records,
                  ChildStorage children = {})
      : Base(tree),
        records_(std::move(records)),
        children_(std::move(children)) {
    for (Node node = 0; node < records_.size(); ++node) {
      const auto& record = records_[node];
      const auto count = record.count();
      for (std::size_t i = 0; i < count; ++i) {
        auto& child = records_[child_impl(node, i)];
        child.parent = node;
        if constexpr (!ContiguousChildren) {
          child.sibling_position = record.first + i;
        }
      }
    }
  }

  /** @brief Source handle copied into this snapshot for this local node. */
  SourceNode source_node(Node node) const { return records_[node].source; }

 private:
  friend Base;

  std::size_t size_impl() const { return records_.size(); }
  std::size_t root_count_impl() const { return !records_.empty(); }
  bool children_complete_impl(Node node) const {
    return records_[node].complete();
  }
  std::size_t degree_impl(Node node) const { return records_[node].count(); }
  Node child_impl(Node node, std::size_t i) const {
    const auto position = records_[node].first + i;
    if constexpr (ContiguousChildren) {
      return position;
    } else {
      return children_[position];
    }
  }
  auto children_impl(Node node) const {
    const auto& record = records_[node];
    if constexpr (ContiguousChildren) {
      return std::views::iota(record.first, record.first + record.count());
    } else {
      return std::span<const Node>(children_).subspan(record.first,
                                                      record.count());
    }
  }
  Node parent_impl(Node node) const { return records_[node].parent; }
  Node next_sibling_impl(Node node) const {
    const auto& parent = records_[parent_impl(node)];
    if constexpr (ContiguousChildren) {
      return node + 1 < parent.first + parent.count() ? node + 1 : node;
    } else {
      const auto next = records_[node].sibling_position + 1;
      return next < parent.first + parent.count() ? children_[next] : node;
    }
  }
  std::size_t memory_usage_bytes_impl() const {
    const auto records_bytes =
        sizeof(*this) + records_.capacity() * sizeof(Record);
    if constexpr (ContiguousChildren) {
      return records_bytes;
    } else {
      return records_bytes + children_.capacity() * sizeof(Node);
    }
  }

  std::vector<Record> records_;
  [[no_unique_address]] ChildStorage children_;
};

}  // namespace pixie::tree::detail
