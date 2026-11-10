#pragma once

#include <pixie/tree.h>

#include <concepts>

namespace pixie {
/**
 * @brief Inspection of represented topology bound to an immutable complete
 * tree.
 * @details Handles must come from tree_base(). contains() never fetches data;
 * children_complete() requires a represented node and means all immediate
 * children are represented, not all descendants. represented_children() visits
 * only represented children in original order. represented_parent() returns an
 * optional handle: nullopt means outside the representation (the complete root
 * returns itself). No operation silently falls back to complete-tree
 * navigation. The complete tree must remain alive, immutable and unmoved.
 * Concrete objects document lookup complexity and how missing-node
 * preconditions are checked.
 */
template <class T>
concept IncompleteTree = requires(const T& tree, TreeNode node) {
  tree.tree_base();
  { tree.contains(node) } -> std::same_as<bool>;
  { tree.children_complete(node) } -> std::same_as<bool>;
  tree.represented_parent(node);
  tree.represented_children(node);
};
}  // namespace pixie
