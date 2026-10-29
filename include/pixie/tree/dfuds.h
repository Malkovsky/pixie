#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmm/tree.h>
#include <pixie/tree.h>
#include <pixie/tree/cache.h>
#include <pixie/tree/detail/dfuds_node_view.h>
#include <pixie/tree/detail/unpack_cache.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {
namespace detail {
struct NavigationChunkAccess;
}

/**
 * @brief A tree class based on the depth-first unary degree sequence (DFUDS)
 * representation
 * @tparam CachePolicy Navigation cache policy satisfying tree/cache.h.
 */
template <typename RMMTree, class CachePolicy = LocalNavigationCache>
class DFUDSTree : public TreeBase<DFUDSTree<RMMTree, CachePolicy>> {
 private:
  const size_t num_bits_;
  RMMTree rmm_;
  const uint64_t* words_;

  // Sibling edges match the parent's unary ones from right to left. The last
  // child's edge is therefore the first one of that run, preceded by zero
  // (or at bit zero for the root's degree run).
  bool starts_degree_run(size_t open) const {
    return open == 0 || rmm_.bit(open - 1) == 0;
  }

  // For a non-root first child, the parent's degree run ends in the matching
  // adjacent "10" pair immediately before the child's record.
  size_t incoming_open(const TreeNode& node) const {
    return rmm_.bit(node.pos - 2) ? node.pos - 2 : rmm_.open(node.pos - 1);
  }

  TreeNode parent_from_open(const TreeNode& node, size_t open) const {
    // [open, node.pos) is balanced: half its bits are zero delimiters.
    // Subtract those records from the child's preorder number.
    const auto number = node.number - (node.pos - open) / 2;
    // For a last child, the matching one is already the parent's record start.
    const auto pos = starts_degree_run(open) ? open : rmm_.select0(number) + 1;
    return TreeNode(number, pos);
  }

 public:
  using Node = TreeNode;

  /** @brief Encoding-specific query memoization used by the local cache policy.
   */
  using LocalNodeView = detail::DfudsNodeView<DFUDSTree>;
  /** @brief Cache context selected by CachePolicy; owned by this tree. */
  using Cache = typename CachePolicy::template Cache<DFUDSTree>;
  /** @brief Node view supplied by the selected cache context. */
  using NodeView = typename Cache::NodeView;

 private:
  friend class detail::DfudsNodeView<DFUDSTree>;
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
   * @brief Constructor from an external array of uint64_t
   * @details The source words must remain alive and immutable for this tree
   * and its traversals. No padding beyond the logical final word is required.
   *
   * @param dfuds_sequence parenthesis sequence in dfuds representation
   * @param cache_bytes Cache plus traversal-buffer limit; zero disables shared
   * retention. Immutable constexpr tables and caller-owned views are excluded.
   */
  explicit DFUDSTree(const std::vector<std::uint64_t>& dfuds_sequence,
                     size_t tree_size,
                     size_t cache_bytes = kDefaultTreeCacheBytes)
      : num_bits_(2 * tree_size - 1),
        rmm_(dfuds_sequence, 2 * tree_size - 1),
        words_(dfuds_sequence.data()),
        cache_(cache_bytes),
        cache_limit_(cache_bytes) {}

  /**
   * @brief Returns the root node
   */
  Node root_impl() const { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (num_bits_ + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const { return rmm_.bit(node.pos) == 0; }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const Node& node) const {
    return rmm_.select0(node.number + 1) - node.pos;
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = rmm_.select0(node.number + 1);
    size_t num = node.number + 1;
    return Node(num, pos + 1);
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    auto parent = LocalNodeView(*this, node);
    const auto count = parent.degree();
    if (count == 0) {
      return;
    }
    auto child = LocalNodeView(*this, parent.first_child());
    for (size_t i = 0; i < count; ++i) {
      visitor(child.node());
      if (i + 1 < count) {
        child.reset(child.next_sibling());
      }
    }
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    // A leaf occupies exactly one zero and the next record is its sibling.
    if (is_leaf_impl(node)) {
      return Node(node.number + 1, node.pos + 1);
    }
    size_t end = rmm_.fwdsearch(node.pos, -1);
    size_t pos = end + 1;
    // Preorder numbering advances by the subtree's node count.
    size_t num = node.number + (pos - node.pos + 1) / 2;
    return Node(num, pos);
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_open(node, incoming_open(node));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    return starts_degree_run(incoming_open(node));
  }
};

std::vector<uint64_t> adj_to_dfuds(
    size_t tree_size,
    const std::vector<std::vector<size_t>>& adj) {
  size_t dfuds_size = tree_size * 2 - 1;
  std::vector<uint64_t> dfuds((dfuds_size + 63) / 64, 0);
  std::vector<size_t> stack;
  stack.push_back(0);
  size_t pos = 0;
  while (!stack.empty()) {
    auto v = stack.back();
    stack.pop_back();
    size_t edge_count = adj[v].size();
    for (size_t i = 0; i < edge_count - 1; ++i) {  // edge 0 goes to parent
      dfuds[pos >> 6] = dfuds[pos >> 6] | (1ULL << (pos & 63));
      pos++;
      stack.push_back(adj[v][edge_count - 1 - i]);
    }
    pos++;
  }
  return dfuds;
}

}  // namespace pixie
