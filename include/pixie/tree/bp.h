#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmm/tree.h>
#include <pixie/tree.h>
#include <pixie/tree/cache.h>
#include <pixie/tree/detail/bp_node_view.h>
#include <pixie/tree/detail/unpack_cache.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace pixie {
namespace detail {
struct NavigationChunkAccess;
}

/**
 * @brief A tree class based on the balances parentheses (BP)
 * representation
 * @tparam CachePolicy Navigation cache policy satisfying tree/cache.h.
 */
template <typename RMMTree, class CachePolicy = LocalNavigationCache>
class BPTree : public TreeBase<BPTree<RMMTree, CachePolicy>> {
 private:
  const size_t num_bits_;
  RMMTree rmm_;
  const uint64_t* words_;

 public:
  /** @brief Encoding-specific query memoization used by the local cache policy.
   */
  using LocalNodeView = detail::BpNodeView<BPTree>;
  /** @brief Cache context selected by CachePolicy; owned by this tree. */
  using Cache = typename CachePolicy::template Cache<BPTree>;
  /** @brief Node view supplied by the selected cache context. */
  using NodeView = typename Cache::NodeView;

 private:
  friend class detail::BpNodeView<BPTree>;
  friend struct detail::NavigationChunkAccess;
  [[no_unique_address]] Cache cache_;
  size_t cache_limit_;

 public:
  /**
   * @brief Create an independent view for a valid node of this tree.
   * @details The tree, its cache context, and borrowed source must outlive the
   * returned view and remain unmoved. CachePolicy defines cache sharing and
   * synchronization; node handles remain independent of cache entry locations.
   */
  NodeView navigate(TreeNode node) const {
    return cache_.navigate(*this, node);
  }

  /** @brief Cache and traversal-buffer allowance; excludes source/RmM storage.
   */
  size_t cache_limit_bytes() const { return cache_limit_; }
  /** @brief Object, owned index allocations and cache-reported storage.
   * @details Cache policies may include reserved scratch. Excludes borrowed
   * words, external traversal contexts and allocator metadata. */
  size_t memory_usage_bytes() const
    requires HasMemoryUsageBytes<RMMTree>
  {
    return sizeof(*this) + nested_owned_memory_bytes(rmm_) +
           nested_owned_memory_bytes(cache_);
  }
  /** @brief Navigation-policy storage plus reserved temporary buffers. */
  size_t cache_memory_usage_bytes() const {
    return cache_.memory_usage_bytes();
  }

  /**
   * @brief A node class of BP tree
   */
  using Node = TreeNode;

  /**
   * @brief Constructor from an external array of uint64_t
   * @details The source words must remain alive and immutable for this tree
   * and its traversals. No padding beyond the logical final word is required.
   * @param cache_bytes Cache plus traversal-buffer limit; zero disables shared
   * retention. Immutable constexpr tables and caller-owned views are excluded.
   */
  explicit BPTree(const std::vector<std::uint64_t>& words,
                  size_t tree_size,
                  size_t cache_bytes = kDefaultTreeCacheBytes)
      : num_bits_(2 * tree_size),
        rmm_(words, 2 * tree_size),
        words_(words.data()),
        cache_(cache_bytes),
        cache_limit_(cache_bytes) {}

  /**
   * @brief Returns the root node
   */
  Node root_impl() const { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return num_bits_ / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const {
    return rmm_.bit(node.pos + 1) == 0;
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   *    this method has O(d) time complexity!
   *
   *    TODO try make this faster
   */
  size_t degree_impl(const Node& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    Node child = first_child_impl(node);
    size_t child_count = 1;
    while (true) {
      if (is_last_child_impl(child)) {
        return child_count;
      }
      child = next_sibling_impl(child);
      child_count++;
    }
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = node.pos + 1;
    size_t num = node.number + 1;
    return Node(num, pos);
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    if (is_leaf_impl(node)) {
      return;
    }
    auto child = LocalNodeView(*this, first_child_impl(node));
    while (true) {
      visitor(child.node());
      if (child.is_last_child()) {
        return;
      }
      child.reset(child.next_sibling());
    }
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   *    this method has O(i) time complexity!
   *
   *    TODO try make this faster
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    // An opening immediately before this node identifies its parent. It is
    // also the preceding preorder node, so neither enclose nor rank is needed.
    if (rmm_.bit(node.pos - 1)) {
      return Node(node.number - 1, node.pos - 1);
    }
    size_t pos = rmm_.enclose(node.pos);
    // Between the parent's open and this child are complete sibling subtrees,
    // each occupying two bits per node, plus the parent's own opening bit.
    size_t num = node.number - (node.pos - pos + 1) / 2;
    return Node(num, pos);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    const size_t end = is_leaf_impl(node) ? node.pos + 1 : rmm_.close(node.pos);
    // A non-root node's close is followed by a sibling open or a parent close.
    return rmm_.bit(end + 1) == 0;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    // A leaf is exactly "10"; its sibling is the next preorder node.
    if (is_leaf_impl(node)) {
      return Node(node.number + 1, node.pos + 2);
    }
    size_t pos = rmm_.close(node.pos) + 1;
    size_t num = node.number + (pos - node.pos) / 2;
    return Node(num, pos);
  }
};

std::vector<uint64_t> adj_to_bp(size_t tree_size,
                                const std::vector<std::vector<size_t>>& adj) {
  size_t bp_size = tree_size * 2;
  std::vector<uint64_t> bp((bp_size + 63) / 64, 0);
  std::vector<std::pair<size_t, size_t>> stack;
  stack.push_back(std::make_pair(0, 0));
  size_t pos = 0;
  bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
  while (!stack.empty()) {
    auto& [v, p] = stack.back();
    p++;
    if (p >= adj[v].size()) {
      pos++;
      stack.pop_back();
      continue;
    }
    pos++;
    bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
    stack.push_back(std::make_pair(adj[v][p], 0));
  }
  return bp;
}
}  // namespace pixie
