#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection_iterative.h>
#include <pixie/tree/detail/traversal_profiles.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

namespace pixie {
namespace detail::tree_traversal {
template <class Tree>
struct Encoding;
template <class B, detail::TreeEncodingStorage S>
struct Encoding<BPTree<B, S>> {
  static constexpr int value = 0;
};
template <class B, detail::TreeEncodingStorage S>
struct Encoding<DFUDSTree<B, S>> {
  static constexpr int value = 1;
};
template <class B, detail::TreeEncodingStorage S>
struct Encoding<LoudsTree<B, S>> {
  static constexpr int value = 2;
};
struct SizingVisitor {
  void operator()(TreeNode) const {}
};
}  // namespace detail::tree_traversal

/**
 * @brief Runtime-budgeted root index and decoding configuration.
 * @details Borrows a nonempty immutable tree and obtains its encoded words;
 * the tree must outlive this object and its cursors. Its constructor's byte
 * limit covers this external context, the immutable root index and peak
 * explicit traversal buffers. The tree owns none of these resources. Up to half
 * the available bytes are assigned to decoding, capped at 192 KiB for BP/DFUDS
 * and 80 KiB for LOUDS. Plain streaming needs only 4 KiB and remains available
 * when inspected decoding cannot fit. Remaining bytes hold the root index. The
 * constructor selects compiled chunk/level sizes; BP/DFUDS cached recursion
 * depth is additionally computed at runtime.
 *
 * Zero/tiny budgets use direct navigation. Partial child lists and deeper
 * subtrees retain exact navigation fallbacks. Constexpr tables, source/index
 * storage, allocator metadata, compiler call-frame overhead, user callback
 * storage and construction's temporary priority queue are excluded. The limit
 * covers one active traversal; nested/concurrent walks require independent
 * scratch allowances. This context never owns mutable decoded data. It provides
 * root cursors and a decoding plan; DFS event orchestration remains in the
 * benchmark harness. Inactive contexts consume their fixed object bytes even
 * when the optimization budget charge is zero.
 */
template <class Tree>
class TreeTraversalContext {
  static constexpr int encoding = detail::tree_traversal::Encoding<Tree>::value;
  static constexpr size_t kStreamBytes = 4096;
  struct Plan {
    size_t scratch = 0;
    size_t bits = 0;
    size_t depth = 0;
  };
  template <size_t Bits, class V, class I>
  using ParenthesisWalk = std::conditional_t<
      encoding == 0,
      detail::bp_micro::BpMicroTraversal<Bits, Tree, V, I>,
      detail::dfuds_micro::DfudsMicroTraversal<Bits, Tree, V, I>>;
  using SizingVisitor = detail::tree_traversal::SizingVisitor;
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Plan plan_;
  RootNavigationIndex<Tree> index_;

  static constexpr size_t overhead() {
    return sizeof(TreeTraversalContext) - sizeof(RootNavigationIndex<Tree>);
  }
  template <size_t Bits>
  static Plan parentheses_plan(size_t allowance) {
    using Walk = ParenthesisWalk<Bits, SizingVisitor, SizingVisitor>;
    if (Walk::scratch_bytes(8) > allowance) {
      return {};
    }
    const size_t depth =
        std::min(size_t{128},
                 (allowance - Walk::scratch_bytes(0)) / (64 * sizeof(size_t)));
    return {std::max(kStreamBytes, Walk::scratch_bytes(depth)), Bits, depth};
  }
  template <size_t Depth, size_t Children>
  static Plan louds_plan(size_t allowance) {
    using Walk = detail::tree_unpack::LoudsIterativeTraversal<
        Tree, SizingVisitor, SizingVisitor, true, Depth, Children>;
    constexpr size_t bytes = std::max(kStreamBytes, Walk::scratch_bytes());
    return bytes <= allowance ? Plan{bytes, 0, Depth} : Plan{};
  }
  static Plan choose_plan(size_t budget) {
    if (budget < sizeof(TreeTraversalContext) + kStreamBytes) {
      return {};
    }
    const size_t available = budget - sizeof(TreeTraversalContext);
    const size_t allowance = std::min(size_t{encoding == 2 ? 80 : 192} * 1024,
                                      std::max(kStreamBytes, available / 2));
    using Profiles =
        std::conditional_t<encoding == 2, detail::tree_traversal::LoudsProfiles,
                           detail::tree_traversal::ParenthesisProfiles>;
    Plan plan;
    if (Profiles::visit_until([&](auto profile) {
          using Profile = decltype(profile);
          if constexpr (encoding == 2) {
            plan = louds_plan<Profile::depth, Profile::children>(allowance);
          } else {
            plan = parentheses_plan<Profile::bits>(allowance);
          }
          return plan.depth != 0;
        })) {
      return plan;
    }
    return {kStreamBytes, 0, 0};
  }

 public:
  /** @brief Select index and scratch sizes within byte_limit.
   * @details After construction this context is read-only. Each
   * traversal consumer owns its mutable scratch. The allowance covers one
   * active traversal; callers budget simultaneous or nested traversals
   * separately. The tree, context and immutable words must remain alive and
   * unmoved while borrowed cursors or decoding operations use them.
   * @throws std::invalid_argument if the tree is empty.
   */
  TreeTraversalContext(const Tree& tree,
                       size_t byte_limit = kDefaultTraversalBytes)
      : tree_(tree),
        words_(tree.encoding_words()),
        plan_(choose_plan(byte_limit)),
        index_(tree,
               byte_limit >= sizeof(*this) + plan_.scratch
                   ? byte_limit - overhead() - plan_.scratch
                   : 0) {}
  /** @brief Adopt a restored root index and choose a fresh decoding plan.
   * @details The index keeps borrowing its original tree. No nodes are
   * expanded. The budget includes this context, retained entries and traversal
   * scratch.
   * @throws std::invalid_argument if the chosen plan and index do not fit the
   * budget. Save root_index() to persist the immutable portion of a context.
   */
  explicit TreeTraversalContext(RootNavigationIndex<Tree> index,
                                size_t byte_limit = kDefaultTraversalBytes)
      : tree_(index.tree()),
        words_(tree_.encoding_words()),
        plan_(choose_plan(byte_limit)),
        index_(std::move(index)) {
    if (budget_usage_bytes() > byte_limit) {
      throw std::invalid_argument(
          "Restored root index exceeds traversal budget");
    }
  }

  /** @brief Peak explicit traversal-buffer allowance, including decoder loads.
   */
  size_t scratch_bytes() const { return plan_.scratch; }
  /** @brief Whether plain DFS uses its streaming decoder. */
  bool streaming_enabled() const { return plan_.scratch >= kStreamBytes; }
  /** @brief Whether inspected DFS can use decoded subtrees/level cursors. */
  bool optimized_inspection_enabled() const { return plan_.depth != 0; }
  /** @brief BP/DFUDS chunk size in bits; zero for LOUDS or direct inspection.
   */
  size_t decoded_chunk_bits() const { return plan_.bits; }
  /** @brief Cached traversal depth; BP/DFUDS keep decoding below this depth.
   */
  size_t cached_depth() const { return plan_.depth; }
  /** @brief Active index/context/scratch budget charge; zero when disabled. */
  size_t budget_usage_bytes() const {
    const size_t indexed = index_.budget_usage_bytes();
    return indexed || plan_.scratch
               ? overhead() + std::max(sizeof(index_), indexed) + plan_.scratch
               : 0;
  }
  /** @brief Actual context/index storage plus reserved per-traversal scratch.
   * @details Scratch is an allowance, not a persistent allocation. Includes
   * fixed object bytes even when disabled; excludes tree/source and metadata.
   */
  size_t memory_usage_bytes() const {
    return overhead() + index_.memory_usage_bytes() + plan_.scratch;
  }
  /** @brief Reject a temporary tree because this context borrows it. */
  TreeTraversalContext(Tree&&, size_t = kDefaultTraversalBytes) = delete;
  /** @brief Reject a temporary const tree for the same reason. */
  TreeTraversalContext(const Tree&&, size_t = kDefaultTraversalBytes) = delete;
  /** @brief Independent root cursor borrowing this index and its tree. */
  auto cursor() const& { return index_.cursor(); }
  /** @brief Reject cursors borrowing a temporary context. */
  auto cursor() const&& = delete;

  /** @brief Borrow the immutable tree; this context does not extend its
   * lifetime. */
  const Tree& tree() const { return tree_; }
  /** @brief Original encoded words borrowed by the context, without padding. */
  std::span<const uint64_t> source_words() const { return words_; }
  /** @brief Borrow the immutable root index; cursors must not outlive it. */
  const RootNavigationIndex<Tree>& root_index() const& { return index_; }
};
}  // namespace pixie
