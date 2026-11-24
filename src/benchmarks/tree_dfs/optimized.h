#pragma once

#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection_iterative.h>
#include <pixie/tree/traversal_context.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

#include "traversal_profiles.h"

namespace pixie::benchmarks {
// Benchmark event contract: preorder entries and one return after the final
// child of every nonleaf. Child inspection occurs before descent.
template <class Tree>
class TreeTraversalWorkloads {
  static constexpr int encoding =
      tree::detail::tree_traversal::Encoding<Tree>::value;
  using Profiles = std::conditional_t<encoding == 2,
                                      detail::LoudsProfiles,
                                      detail::ParenthesisProfiles>;
  template <size_t Index>
  using Profile = typename Profiles::template At<Index>;
  template <size_t Bits, class V, class I>
  using ParenthesisWalk = std::conditional_t<
      encoding == 0,
      tree::detail::bp_micro::BpMicroTraversal<Bits, Tree, V, I>,
      tree::detail::dfuds_micro::DfudsMicroTraversal<Bits, Tree, V, I>>;
  static constexpr size_t kStreamBytes = 4096;
  struct Plan {
    size_t scratch = 0;
    size_t bits = 0;
    size_t depth = 0;
  };
  struct SizingVisitor {
    void operator()(tree::TreeNode) const {}
  };
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Plan plan_;

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
    using Walk = tree::detail::tree_unpack::LoudsIterativeTraversal<
        Tree, SizingVisitor, SizingVisitor, true, Depth, Children>;
    constexpr size_t bytes = std::max(kStreamBytes, Walk::scratch_bytes());
    return bytes <= allowance ? Plan{bytes, 0, Depth} : Plan{};
  }
  static Plan choose_plan(size_t budget) {
    if (budget < fixture_metadata_bytes() + kStreamBytes) {
      return {};
    }
    const size_t available = budget - fixture_metadata_bytes();
    const size_t allowance = std::min(size_t{encoding == 2 ? 80 : 192} * 1024,
                                      std::max(kStreamBytes, available / 2));
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
  explicit TreeTraversalWorkloads(const Tree& tree,
                                  size_t budget = tree::kDefaultTraversalBytes)
      : tree_(tree),
        words_(tree.encoding_words()),
        plan_(choose_plan(budget)) {}

  // The matched fixtures budget one navigation context and one decoder plan
  // for both DFS and root paths, so the index gets the same remaining budget.
  static constexpr size_t fixture_metadata_bytes() {
    return sizeof(tree::TreeTraversalContext<Tree>) + sizeof(Plan);
  }
  size_t reserved_bytes() const {
    return streaming_enabled() ? fixture_metadata_bytes() + scratch_bytes() : 0;
  }
  size_t scratch_bytes() const { return plan_.scratch; }
  bool streaming_enabled() const { return scratch_bytes() >= kStreamBytes; }
  bool optimized_inspection_enabled() const { return plan_.depth != 0; }
  size_t decoded_chunk_bits() const { return plan_.bits; }
  size_t cached_depth() const { return plan_.depth; }

  /** @brief Emit DFS entry and final-child return events to a borrowed visitor.
   */
  template <class Visitor>
  void traverse(Visitor& visit) const {
    if (streaming_enabled()) {
      traverse_decoded(visit);
    } else {
      pixie::tree::detail::tree_unpack::subtree(tree_, tree_.root(), visit);
    }
  }
  /** @brief DFS events plus immediate children in order before descent.
   * @details Callbacks are borrowed; exceptions propagate. Returns
   * (entered nodes, inspected children). Return events do not inspect children.
   */
  template <class Visitor, class Inspector>
  std::pair<size_t, size_t> inspect(Visitor& visit, Inspector& child) const {
    // Keep switch dispatch and separate execution frames. Profile additions
    // must extend dispatch explicitly rather than silently using the fallback.
    if constexpr (encoding == 2) {
      static_assert(Profiles::size == 5);
      switch (cached_depth()) {
        case Profile<0>::depth:
          return inspect_louds<Profile<0>::depth, Profile<0>::children>(visit,
                                                                        child);
        case Profile<1>::depth:
          return inspect_louds<Profile<1>::depth, Profile<1>::children>(visit,
                                                                        child);
        case Profile<2>::depth:
          return inspect_louds<Profile<2>::depth, Profile<2>::children>(visit,
                                                                        child);
        case Profile<3>::depth:
          return inspect_louds<Profile<3>::depth, Profile<3>::children>(visit,
                                                                        child);
        case Profile<4>::depth:
          return inspect_louds<Profile<4>::depth, Profile<4>::children>(visit,
                                                                        child);
      }
    } else {
      static_assert(Profiles::size == 4);
      switch (decoded_chunk_bits()) {
        case Profile<0>::bits:
          return inspect_parentheses<Profile<0>::bits>(visit, child);
        case Profile<1>::bits:
          return inspect_parentheses<Profile<1>::bits>(visit, child);
        case Profile<2>::bits:
          return inspect_parentheses<Profile<2>::bits>(visit, child);
        case Profile<3>::bits:
          return inspect_parentheses<Profile<3>::bits>(visit, child);
      }
    }
    size_t nodes = 0, children = 0;
    auto events = [&](tree::TreeNode n, bool entering) {
      pixie::tree::detail::tree_unpack::emit(visit, n, entering);
      if (entering) {
        ++nodes;
        tree_.for_each_child(n, [&](tree::TreeNode c) {
          ++children;
          child(c);
        });
      }
    };
    pixie::tree::detail::tree_unpack::subtree(tree_, tree_.root(), events);
    return {nodes, children};
  }

 private:
  // Separate frames prevent large-profile buffers appearing on smaller paths.
  template <class Visitor>
  [[gnu::noinline]] void traverse_decoded(Visitor& visit) const {
    if constexpr (encoding == 0) {
      pixie::tree::detail::tree_unpack::bp_byte(tree_, words_, 2 * tree_.size(),
                                                visit);
    } else if constexpr (encoding == 1) {
      pixie::tree::detail::tree_unpack::dfuds_excess_stream(
          tree_, words_, 2 * tree_.size() - 1, visit);
    } else {
      pixie::tree::detail::tree_unpack::louds(tree_, words_,
                                              2 * tree_.size() - 1, visit);
    }
  }
  template <size_t Bits, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_parentheses(
      Visitor& visit,
      Inspector& child) const {
    ParenthesisWalk<Bits, Visitor, Inspector> walk(tree_, words_, visit, child,
                                                   cached_depth());
    walk.walk(tree_.root());
    return {walk.nodes, walk.children};
  }
  template <size_t Depth, size_t Children, class Visitor, class Inspector>
  [[gnu::noinline]] std::pair<size_t, size_t> inspect_louds(
      Visitor& visit,
      Inspector& child) const {
    tree::detail::tree_unpack::LoudsIterativeTraversal<Tree, Visitor, Inspector,
                                                       true, Depth, Children>
        walk(tree_, words_, visit, child);
    walk.walk(tree_.root(), 0);
    return {walk.node_count, walk.child_count};
  }
};
}  // namespace pixie::benchmarks
