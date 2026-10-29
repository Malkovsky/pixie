#pragma once

#include <pixie/tree.h>

#include <cstddef>
#include <utility>

namespace pixie {

/** @brief Default per-tree cache plus traversal-buffer allowance, in bytes. */
inline constexpr std::size_t kDefaultTreeCacheBytes = 512 * 1024;

/**
 * @brief Navigation cache policy preserving encoding-specific query
 * memoization.
 * @details A tree cache policy exposes `template<class Tree> class Cache`.
 * The context is constructed with a byte limit as part of the tree and must not
 * access the tree during construction. It supplies a `NodeView` type and
 * `NodeView navigate(const Tree&, TreeNode) const`. The supplied tree always
 * owns that context. A view may borrow the context and tree; both must remain
 * alive, unmoved, and backed by unchanged source data while a view is used.
 * Policies must not retain pointers into another context when copied or moved.
 * Contexts provide memory_usage_bytes(), bounded by the constructor limit.
 * Zero disables shared retention. The limit excludes source/index storage,
 * caller-owned views and allocator overhead. Traversal adapters must reserve
 * their temporary buffers from the same allowance before retaining entries.
 * Immutable constexpr lookup/dispatch tables are excluded.
 *
 * A NodeView is copyable and provides node(), reset(TreeNode), is_root(),
 * is_leaf(), degree(), child(size_t), first_child(), parent(), next_sibling(),
 * is_last_child(), and for_each_child(visitor). Queries obey TreeBase semantics
 * for the current node, including zero-based child indices and root's parent
 * being root itself. Invalid handles/indices have undefined behavior. Queries
 * do not move the view; reset moves only that view. Copies have independent
 * current nodes but may share cache storage according to the policy. Stable
 * TreeNode handles must never depend on cache-entry addresses or eviction.
 *
 * Only tree.navigate() uses this policy; direct TreeBase queries bypass it.
 * On a miss, policies resolve information through the tree's direct operations
 * or its encoding-specific LocalNodeView; they must not call tree.navigate(),
 * which would re-enter the policy. Decoding may return uncached results when
 * retention is not possible. Policies define cache storage, synchronization,
 * and retention rules; implementing partial child storage is optional.
 *
 * This policy owns no shared storage: all memoized results live in the view.
 * Copies have independent caches. Independent views can be used concurrently;
 * sharing a single view between threads requires external synchronization.
 */
struct LocalNavigationCache {
  /** @brief Stateless context selecting Tree's encoding-specific memoized view.
   */
  template <class Tree>
  class Cache {
   public:
    /** @brief View caching query results and reusable decoding intermediates.
     */
    using NodeView = typename Tree::LocalNodeView;

    /** @brief Accept a budget; this policy retains no shared storage. */
    explicit Cache(std::size_t = kDefaultTreeCacheBytes) {}
    /** @brief Shared persistent cache bytes (always zero). */
    std::size_t memory_usage_bytes() const { return 0; }

    /** @brief Bind to a valid node; lifetime and preconditions are stated
     * above. */
    NodeView navigate(const Tree& tree, TreeNode node) const {
      return NodeView(tree, node);
    }
  };
};

/**
 * @brief Navigation policy forwarding every query to the representation.
 * @details Implements LocalNavigationCache's policy contract without retaining
 * query results or allocating storage. Copies have independent current nodes.
 * Independent views support concurrent use; reset requires exclusive access
 * to the affected view. The tree and immutable backing data outlive its views.
 */
struct NoNavigationCache {
  /** @brief Stateless context providing direct-query node views. */
  template <class Tree>
  class Cache {
   public:
    /** @brief Accept a budget; this policy retains no shared storage. */
    explicit Cache(std::size_t = kDefaultTreeCacheBytes) {}
    /** @brief Shared persistent cache bytes (always zero). */
    std::size_t memory_usage_bytes() const { return 0; }
    /** @brief Borrowed node view with TreeBase query semantics and no cache. */
    class NodeView {
     public:
      /** @brief Return the current valid node handle by value. */
      TreeNode node() const { return node_; }
      /** @brief Move this view to a valid node of the same tree. */
      void reset(TreeNode node) { node_ = node; }
      /** @brief Test whether the current node is root. */
      bool is_root() const { return tree_->is_root(node_); }
      /** @brief Test whether the current node has no children. */
      bool is_leaf() const { return tree_->is_leaf(node_); }
      /** @brief Return the current node's number of children. */
      std::size_t degree() const { return tree_->degree(node_); }
      /** @brief Return child i; requires i in [0, degree()). */
      TreeNode child(std::size_t i) const { return tree_->child(node_, i); }
      /** @brief Return the first child; requires a non-leaf node. */
      TreeNode first_child() const { return tree_->first_child(node_); }
      /** @brief Return the parent, or root itself for root. */
      TreeNode parent() const { return tree_->parent(node_); }
      /** @brief Return the next sibling; requires a non-root, non-last child.
       */
      TreeNode next_sibling() const { return tree_->next_sibling(node_); }
      /** @brief Test for the last sibling; requires a non-root node. */
      bool is_last_child() const { return tree_->is_last_child(node_); }
      /** @brief Enumerate children with TreeBase::for_each_child semantics. */
      template <class Visitor>
      void for_each_child(Visitor&& visitor) const {
        tree_->for_each_child(node_, std::forward<Visitor>(visitor));
      }

     private:
      friend class Cache;
      NodeView(const Tree& tree, TreeNode node) : tree_(&tree), node_(node) {}
      const Tree* tree_;
      TreeNode node_;
    };

    /** @brief Bind to a valid node; lifetime and preconditions are stated
     * above. */
    NodeView navigate(const Tree& tree, TreeNode node) const {
      return NodeView(tree, node);
    }
  };
};

}  // namespace pixie
