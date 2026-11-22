#pragma once

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace pixie::tree {

/** @brief Owning inspection snapshot of the nodes and edges in one tree part.
 * @details Child links index nodes() from zero. The first
 * root_count() nodes are the local roots in source order: they have no parent
 * in the part, though their source parents may exist. Remaining node order is
 * producer-specific. No record borrows a tree, index or decoder.
 * Inspection allocations are separate from navigation and cache budgets.
 */
class TreePart {
 public:
  /** @brief A source node and the ordered subset of its children in this part.
   */
  struct Node {
    /// Number in the source encoding, when available; not necessarily preorder.
    std::optional<std::size_t> source_id;
    /// Original encoding bit position, when meaningful.
    std::optional<std::size_t> bit_position;
    /// True when the producer established that all immediate children are
    /// present. An empty list proves a leaf only when complete.
    bool children_complete = true;
    std::vector<std::size_t> children;
  };

  /** @brief Take ownership of producer records. The producer supplies valid
   * local indices and lists in source child order; no validation is performed.
   */
  TreePart(std::vector<Node> nodes, std::size_t root_count)
      : nodes_(std::move(nodes)), root_count_(root_count) {}

  /** @brief True when no nodes are represented, including a disabled index. */
  bool empty() const { return nodes_.empty(); }
  /** @brief Local node records, owned by this snapshot. */
  const std::vector<Node>& nodes() const { return nodes_; }
  /** @brief Number of leading local roots; zero for an empty part. */
  std::size_t root_count() const { return root_count_; }

 private:
  std::vector<Node> nodes_;
  std::size_t root_count_;
};
}  // namespace pixie::tree
