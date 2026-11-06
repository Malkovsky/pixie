#pragma once

#include <pixie/incomplete_tree.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace pixie {
template <class Tree>
class TreeTraversalContext;

/** @brief An owned, immutable decoded portion of a borrowed complete tree.
 * @details Records form a connected tree rooted at the requested scan node.
 * Membership lookup is O(log n); represented children are returned in original
 * order. An absent parent is outside the result, not necessarily nonexistent.
 * Returned records survive context reuse. The source tree must outlive them.
 * No operation performs decoding or source-tree navigation.
 */
template <class Tree>
class DecodedTree {
  friend class TreeTraversalContext<Tree>;
  static constexpr size_t missing = std::numeric_limits<size_t>::max();
  struct Record {
    TreeNode node;
    size_t parent = missing, first = missing, last = missing, next = missing;
    bool complete = false;
  };
  const Tree* tree_;
  std::vector<Record> records_;
  size_t bits_ = 0, regions_ = 0;

  size_t find(TreeNode node) const {
    const auto it = std::lower_bound(
        records_.begin(), records_.end(), node.number,
        [](const Record& r, size_t number) { return r.node.number < number; });
    return it != records_.end() && it->node.number == node.number &&
                   it->node.pos == node.pos
               ? size_t(it - records_.begin())
               : missing;
  }
  size_t require(TreeNode node) const {
    auto i = find(node);
    if (i == missing) {
      throw std::out_of_range("node is outside decoded tree");
    }
    return i;
  }
  size_t add(TreeNode node, size_t parent) {
    size_t i = records_.size();
    records_.push_back({node, parent});
    if (parent != missing) {
      auto& p = records_[parent];
      if (p.first == missing) {
        p.first = i;
      } else {
        records_[p.last].next = i;
      }
      p.last = i;
    }
    return i;
  }

 public:
  /** @brief Empty decoded result bound to tree; tree must outlive this object.
   */
  explicit DecodedTree(const Tree& tree) : tree_(&tree) {}
  /** @brief Reject a dangling source tree. */
  DecodedTree(Tree&&) = delete;
  /** @brief Reject a dangling const source tree. */
  DecodedTree(const Tree&&) = delete;
  /** @brief Original complete tree, without ownership transfer. */
  const Tree& tree_base() const { return *tree_; }
  /** @brief Number of represented node records. */
  size_t size() const { return records_.size(); }
  /** @brief Whether no node record was admitted. */
  bool empty() const { return records_.empty(); }
  /** @brief Logical bits examined, including repeated examinations. */
  size_t scanned_bits() const { return bits_; }
  /** @brief Number of nonempty contiguous scan regions. */
  size_t scanned_regions() const { return regions_; }
  /** @brief Object and allocated records, excluding backing tree. */
  size_t memory_usage_bytes() const {
    return sizeof(*this) + records_.capacity() * sizeof(Record);
  }
  /** @brief Membership for a handle from tree_base(), without fallback. */
  bool contains(TreeNode node) const { return find(node) != missing; }
  /** @brief Full immediate child list known; throws for an absent node. */
  bool children_complete(TreeNode node) const {
    return records_[require(node)].complete;
  }
  /** @brief Requested root if admitted, otherwise nullopt. */
  std::optional<TreeNode> root() const {
    return empty() ? std::nullopt : std::optional(records_[0].node);
  }
  /** @brief Represented parent, nullopt at an external boundary; complete root
   * returns itself. Throws for an absent node. */
  std::optional<TreeNode> represented_parent(TreeNode node) const {
    auto i = records_[require(node)].parent;
    if (i != missing) {
      return records_[i].node;
    }
    return node.number == 0 ? std::optional(node) : std::nullopt;
  }
  /** @brief Represented immediate children in order, possibly incomplete.
   * Throws for an absent node. */
  std::vector<TreeNode> represented_children(TreeNode node) const {
    std::vector<TreeNode> result;
    for (auto i = records_[require(node)].first; i != missing;
         i = records_[i].next) {
      result.push_back(records_[i].node);
    }
    return result;
  }
  /** @brief Visit all records in encoding node-number order. */
  template <class Visitor>
  void for_each_node(Visitor&& visit) const {
    for (const auto& record : records_) {
      visit(record.node);
    }
  }
  /** @brief DFS over represented topology, with enter and exit for every node.
   * This describes the partial result, not a complete-source DFS prefix.
   * Call children_complete() to interpret boundary nodes. Uses O(1) local
   * traversal state; independent calls may interleave or run concurrently. */
  template <class Visitor>
  void events(Visitor&& visit) const {
    if (empty()) {
      return;
    }
    size_t i = 0;
    while (true) {
      visit(records_[i].node, true);
      if (records_[i].first != missing) {
        i = records_[i].first;
        continue;
      }
      while (true) {
        visit(records_[i].node, false);
        if (records_[i].next != missing) {
          i = records_[i].next;
          break;
        }
        i = records_[i].parent;
        if (i == missing) {
          return;
        }
      }
    }
  }
};
}  // namespace pixie
