#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rank_select/support.h>
#include <pixie/tree.h>
#include <pixie/tree/cache.h>
#include <pixie/tree/detail/louds_node_view.h>
#include <pixie/tree/detail/unpack_louds.h>

#include <bit>
#include <cstdint>
#include <span>

namespace pixie {
namespace detail {
struct NavigationChunkAccess;
}

/**
 * @brief A node class of LOUDS tree
 */
using LoudsNode = TreeNode;

/**
 * @brief A tree class based on the level order unary degree sequence (LOUDS)
 * representation
 * @tparam CachePolicy Navigation cache policy satisfying tree/cache.h.
 * @tparam RankSelect Rank/select backend constructible from packed source words
 * and a bit count, with both select directions enabled by default. It must
 * provide size(), operator[], select() and select0() with RankSelectBase
 * semantics. Source words must remain alive and immutable for the tree
 * lifetime.
 * @details Each BFS-ordered node is encoded as one followed by one zero per
 * child. A handle stores its zero-based BFS number and the position of its one.
 * Siblings therefore occupy consecutive one delimiters. The parent delimiter
 * is the nearest preceding one before the node's incoming zero edge. Degree,
 * sibling and parent navigation use a single-word search for those delimiters,
 * falling back to indexed select when the delimiter lies outside that word.
 */
template <class RankSelect = RankSelectSupport<>,
          class CachePolicy = LocalNavigationCache>
class LoudsTree : public TreeBase<LoudsTree<RankSelect, CachePolicy>> {
 private:
  RankSelect bv;
  const std::uint64_t* words_;

  // Node delimiters are ones. Search only the remainder of this word so a
  // high-degree node cannot turn navigation into an unbounded linear scan.
  size_t next_delimiter(size_t position, size_t rank) const {
    const size_t start = position + 1;
    const auto word = words_[start / 64] >> (start % 64);
    return word ? start + std::countr_zero(word) : bv.select(rank);
  }

  size_t parent_delimiter(size_t zero_position, size_t rank) const {
    // zero_position itself is zero, so retaining it avoids a shift by 64.
    const auto word = words_[zero_position / 64] &
                      (~std::uint64_t{0} >> (63 - zero_position % 64));
    return word ? (zero_position / 64) * 64 + 63 - std::countl_zero(word)
                : bv.select(rank);
  }

  LoudsNode parent_from_edge(const LoudsNode& node, size_t zero_pos) const {
    const size_t parent_number = zero_pos - node.number;
    return LoudsNode(parent_number,
                     parent_delimiter(zero_pos, parent_number + 1));
  }

 public:
  /** @brief Encoding-specific query memoization used by the local cache policy.
   */
  using LocalNodeView = detail::LoudsNodeView<LoudsTree>;
  /** @brief Cache context selected by CachePolicy; owned by this tree. */
  using Cache = typename CachePolicy::template Cache<LoudsTree>;
  /** @brief Node view supplied by the selected cache context. */
  using NodeView = typename Cache::NodeView;

 private:
  friend class detail::LoudsNodeView<LoudsTree>;
  friend struct detail::NavigationChunkAccess;
  [[no_unique_address]] Cache cache_;
  size_t cache_limit_;

 public:
  /**
   * @brief Create an independent navigation view for a valid node of this tree.
   * @details The tree, its cache context, and borrowed source must outlive the
   * returned view and remain unmoved. CachePolicy defines cache sharing and
   * synchronization; node handles remain independent of cache entry locations.
   */
  NodeView navigate(LoudsNode node) const {
    return cache_.navigate(*this, node);
  }

  /** @brief Cache and traversal-buffer allowance; excludes source/rank-select.
   */
  size_t cache_limit_bytes() const { return cache_limit_; }
  /** @brief Object, owned index allocations and cache-reported storage.
   * @details Cache policies may include reserved scratch. Excludes borrowed
   * words, external traversal contexts and allocator metadata. */
  size_t memory_usage_bytes() const
    requires HasMemoryUsageBytes<RankSelect>
  {
    return sizeof(*this) + nested_owned_memory_bytes(bv) +
           nested_owned_memory_bytes(cache_);
  }
  /** @brief Navigation-policy storage plus reserved temporary buffers. */
  size_t cache_memory_usage_bytes() const {
    return cache_.memory_usage_bytes();
  }

  /**
   * @brief Constructor from an external array of uint64_t
   * @param cache_bytes Cache plus traversal-buffer limit; zero disables shared
   * retention. Immutable constexpr tables and caller-owned views are excluded.
   */
  explicit LoudsTree(std::span<const uint64_t> louds,
                     size_t tree_size,
                     size_t cache_bytes = kDefaultTreeCacheBytes)
      : bv(louds, 2 * tree_size - 1),
        words_(louds.data()),
        cache_(cache_bytes),
        cache_limit_(cache_bytes) {}

  /**
   * @brief Returns the root node
   */
  LoudsNode root_impl() const { return LoudsNode(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (bv.size() + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const LoudsNode& node) const {
    return (node.pos + 1 == bv.size()) or bv[node.pos + 1];
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const LoudsNode& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const LoudsNode& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    return next_delimiter(node.pos, node.number + 2) - node.pos - 1;
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    detail::tree_unpack::louds_children(
        *this,
        std::span<const uint64_t>(words_,
                                  bv.size() / 64 + (bv.size() % 64 != 0)),
        bv.size(), node, 0, degree_impl(node), visitor);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  LoudsNode child_impl(const LoudsNode& node, size_t i) const {
    size_t zeros = node.pos + i + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns first child of a @p node
   */
  LoudsNode first_child_impl(const LoudsNode& node) const {
    size_t zeros = node.pos + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  LoudsNode parent_impl(const LoudsNode& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_edge(node, bv.select0(node.number));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const LoudsNode& node) const {
    size_t zero_pos = bv.select0(node.number);
    return bv[zero_pos + 1];
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  LoudsNode next_sibling_impl(const LoudsNode& node) const {
    size_t sibling_number = node.number + 1;
    return LoudsNode(sibling_number,
                     next_delimiter(node.pos, sibling_number + 1));
  }
};

}  // namespace pixie
