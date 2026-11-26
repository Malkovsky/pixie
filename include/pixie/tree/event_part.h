#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/detail/indexed_tree_part.h>
#include <pixie/tree/dfuds.h>

#include <optional>
#include <utility>
#include <vector>

namespace pixie::tree {
namespace detail {
template <class Decoder, class Tree>
std::optional<IndexedTreePart<Tree, false>> inspect_decoded_events(
    const Tree& tree,
    TreeNode root) {
  using Part = IndexedTreePart<Tree, false>;
  Decoder decoder;
  if (tree.empty() || !decoder.load(tree.encoding_words(),
                                    Decoder::bit_size(tree.size()), root)) {
    return std::nullopt;
  }
  std::vector<typename Part::Record> rows;
  rows.reserve(decoder.size());
  std::vector<typename Part::Node> children;
  children.reserve(decoder.size() - 1);
  auto visit = [&](TreeNode node, bool entering) {
    if (entering) {
      rows.emplace_back(node, children.size(), 0, true);
    }
  };
  // The decoder emits each child list immediately after its parent's entry
  // and before visiting its first child.
  auto inspect = [&](TreeNode child) {
    children.push_back(child.number - root.number);
    rows.back().append_child();
  };
  decoder.walk(visit, inspect);
  return Part(tree, std::move(rows), std::move(children));
}
}  // namespace detail

/** @brief Copy the subtree topology from one bounded BP decoder probe.
 * @details root must be a valid handle from tree (ignored on an empty tree).
 * Returns nullopt for an empty tree or when the subtree does not close within
 * Bits bits; no prefix or navigation fallback is returned. On a hit the probe
 * spans the opening bit through its matching closing bit, exactly 2*n bits
 * for n local nodes. The root's parent and siblings are outside the part.
 * Bits is a multiple of 128 in [128,32768]. */
template <size_t Bits = 256, class Backend, detail::TreeEncodingStorage Storage>
auto inspect_event_part(const BPTree<Backend, Storage>& tree, TreeNode root) {
  return detail::inspect_decoded_events<detail::bp_micro::BpMicroSubtree<Bits>>(
      tree, root);
}

/** @brief Inspect one bounded DFUDS decoder probe, with the same hit/miss and
 * handle contract as BP. The probe starts at the root's degree run and ends
 * just after the final descendant's zero delimiter (2*n-1 bits for n nodes).
 * Neither word padding nor any following sibling record belongs to the part. */
template <size_t Bits = 256, class Backend, detail::TreeEncodingStorage Storage>
auto inspect_event_part(const DFUDSTree<Backend, Storage>& tree,
                        TreeNode root) {
  return detail::inspect_decoded_events<
      detail::dfuds_micro::DfudsMicroSubtree<Bits>>(tree, root);
}
}  // namespace pixie::tree
