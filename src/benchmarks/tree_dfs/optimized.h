#pragma once

#include <pixie/tree/traversal_context.h>

namespace pixie::benchmarks {
// Benchmark event contract: preorder entries and one return after the final
// child of every nonleaf. Child inspection occurs before descent.
template <class Tree>
class TreeTraversalWorkloads {
  static constexpr int encoding = detail::tree_traversal::Encoding<Tree>::value;
  template <size_t Bits, class V, class I>
  using ParenthesisWalk = std::conditional_t<
      encoding == 0,
      detail::bp_micro::BpMicroTraversal<Bits, Tree, V, I>,
      detail::dfuds_micro::DfudsMicroTraversal<Bits, Tree, V, I>>;
  const TreeTraversalContext<Tree>& context_;
  const Tree& tree_;
  std::span<const uint64_t> words_;

 public:
  explicit TreeTraversalWorkloads(const TreeTraversalContext<Tree>& context)
      : context_(context),
        tree_(context.tree()),
        words_(context.source_words()) {}
  /** @brief Emit DFS entry and final-child return events to a borrowed visitor.
   */
  template <class Visitor>
  void traverse(Visitor& visit) const {
    if (context_.streaming_enabled()) {
      traverse_decoded(visit);
    } else {
      pixie::detail::tree_unpack::subtree(tree_, tree_.root(), visit);
    }
  }
  /** @brief DFS events plus immediate children in order before descent.
   * @details Callbacks are borrowed; exceptions propagate. Returns
   * (entered nodes, inspected children). Return events do not inspect children.
   */
  template <class Visitor, class Inspector>
  std::pair<size_t, size_t> inspect(Visitor& visit, Inspector& child) const {
    if constexpr (encoding == 2) {
      switch (context_.cached_depth()) {
        case 128:
          return inspect_louds<128, 32>(visit, child);
        case 64:
          return inspect_louds<64, 16>(visit, child);
        case 32:
          return inspect_louds<32, 16>(visit, child);
        case 16:
          return inspect_louds<16, 8>(visit, child);
        case 8:
          return inspect_louds<8, 8>(visit, child);
      }
    } else {
      switch (context_.decoded_chunk_bits()) {
        case 8192:
          return inspect_parentheses<8192>(visit, child);
        case 2048:
          return inspect_parentheses<2048>(visit, child);
        case 512:
          return inspect_parentheses<512>(visit, child);
        case 128:
          return inspect_parentheses<128>(visit, child);
      }
    }
    size_t nodes = 0, children = 0;
    auto events = [&](TreeNode n, bool entering) {
      pixie::detail::tree_unpack::emit(visit, n, entering);
      if (entering) {
        ++nodes;
        tree_.for_each_child(n, [&](TreeNode c) {
          ++children;
          child(c);
        });
      }
    };
    pixie::detail::tree_unpack::subtree(tree_, tree_.root(), events);
    return {nodes, children};
  }

 private:
  // Separate frames prevent large-profile buffers appearing on smaller paths.
  template <class Visitor>
  [[gnu::noinline]] void traverse_decoded(Visitor& visit) const {
    if constexpr (encoding == 0) {
      pixie::detail::tree_unpack::bp_byte(tree_, words_, 2 * tree_.size(),
                                          visit);
    } else if constexpr (encoding == 1) {
      pixie::detail::tree_unpack::dfuds_excess_stream(
          tree_, words_, 2 * tree_.size() - 1, visit);
    } else {
      pixie::detail::tree_unpack::louds(tree_, words_, 2 * tree_.size() - 1,
                                        visit);
    }
  }
  template <size_t Bits, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_parentheses(
      Visitor& visit,
      Inspector& child) const {
    ParenthesisWalk<Bits, Visitor, Inspector> walk(tree_, words_, visit, child,
                                                   context_.cached_depth());
    walk.walk(tree_.root());
    return {walk.nodes, walk.children};
  }
  template <size_t Depth, size_t Children, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_louds(
      Visitor& visit,
      Inspector& child) const {
    detail::tree_unpack::LoudsIterativeTraversal<Tree, Visitor, Inspector, true,
                                                 Depth, Children>
        walk(tree_, words_, visit, child);
    walk.walk(tree_.root(), 0);
    return {walk.node_count, walk.child_count};
  }
};
}  // namespace pixie::benchmarks
