#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection_iterative.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

namespace pixie {
namespace detail::tree_traversal {
template <class Tree>
struct Encoding;
template <class B>
struct Encoding<BPTree<B>> {
  static constexpr int value = 0;
};
template <class B>
struct Encoding<DFUDSTree<B>> {
  static constexpr int value = 1;
};
template <class B>
struct Encoding<LoudsTree<B>> {
  static constexpr int value = 2;
};
struct SizingVisitor {
  void operator()(TreeNode) const {}
};
}  // namespace detail::tree_traversal

/**
 * @brief Runtime-budgeted root index and decoding configuration.
 * @details Borrows a nonempty immutable tree and its matching encoded words;
 * both must outlive this object and its cursors. Its constructor's byte limit
 * covers this external context, the immutable root index and peak explicit
 * traversal buffers. The tree owns none of these resources. Up to half the
 * available bytes are assigned to decoding, capped at 192 KiB for BP/DFUDS and
 * 80 KiB for LOUDS. Plain streaming needs only 4 KiB and remains available when
 * inspected decoding cannot fit. Remaining bytes hold the root index. The
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
    if constexpr (encoding == 2) {
      for (auto p :
           {louds_plan<128, 32>(allowance), louds_plan<64, 16>(allowance),
            louds_plan<32, 16>(allowance), louds_plan<16, 8>(allowance),
            louds_plan<8, 8>(allowance)}) {
        if (p.depth) {
          return p;
        }
      }
    } else {
      for (auto p :
           {parentheses_plan<8192>(allowance),
            parentheses_plan<2048>(allowance), parentheses_plan<512>(allowance),
            parentheses_plan<128>(allowance)}) {
        if (p.depth) {
          return p;
        }
      }
    }
    return {kStreamBytes, 0, 0};
  }

 public:
  /** @brief Select index and scratch sizes within byte_limit; words match tree.
   * @details After construction this context is read-only. Each
   * traversal consumer owns its mutable scratch. The allowance covers one
   * active traversal; callers budget simultaneous or nested traversals
   * separately. The tree, context and immutable words must remain alive and
   * unmoved while borrowed cursors or decoding operations use them.
   * @throws std::invalid_argument if the tree is empty or words are too short.
   * A sufficiently long but mismatched encoding violates the precondition.
   */
  TreeTraversalContext(const Tree& tree,
                       std::span<const uint64_t> words,
                       size_t byte_limit = kDefaultTraversalBytes)
      : tree_(tree),
        words_(words),
        plan_(choose_plan(byte_limit)),
        index_(tree,
               byte_limit >= sizeof(*this) + plan_.scratch
                   ? byte_limit - overhead() - plan_.scratch
                   : 0) {
    const size_t bits = 2 * tree.size() - (encoding == 0 ? 0 : 1);
    if (words.size() < bits / 64 + (bits % 64 != 0)) {
      throw std::invalid_argument("Traversal context source is too short");
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
  TreeTraversalContext(Tree&&,
                       std::span<const uint64_t>,
                       size_t = kDefaultTraversalBytes) = delete;
  /** @brief Reject a temporary const tree for the same reason. */
  TreeTraversalContext(const Tree&&,
                       std::span<const uint64_t>,
                       size_t = kDefaultTraversalBytes) = delete;
  /** @brief Reject temporary source vectors; callers supply long-lived words.
   */
  TreeTraversalContext(const Tree&,
                       std::vector<uint64_t>&&,
                       size_t = kDefaultTraversalBytes) = delete;
  /** @brief Reject temporary const source vectors for the same reason. */
  TreeTraversalContext(const Tree&,
                       const std::vector<uint64_t>&&,
                       size_t = kDefaultTraversalBytes) = delete;

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
