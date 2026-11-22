#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/tree_part.h>

namespace pixie::tree {
namespace detail {
template <class Decoder, class Tree>
std::optional<TreePart> inspect_decoded_events(const Tree& tree,
                                               TreeNode root) {
  Decoder decoder;
  if (tree.empty() || !decoder.load(tree.encoding_words(),
                                    Decoder::bit_size(tree.size()), root)) {
    return std::nullopt;
  }
  std::vector<TreePart::Node> nodes(decoder.size());
  size_t parent = 0;
  auto visit = [&](TreeNode n, bool entering) {
    if (entering) {
      parent = n.number - root.number;
      nodes[parent].source_id = n.number;
      nodes[parent].bit_position = n.pos;
    }
  };
  // All inspect callbacks immediately follow their parent's entry, before
  // the first child is visited. No source navigation or stack is needed.
  auto inspect = [&](TreeNode child) {
    const auto i = child.number - root.number;
    nodes[parent].children.push_back(i);
  };
  decoder.walk(visit, inspect);
  return TreePart(std::move(nodes), 1);
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
std::optional<TreePart> inspect_event_part(const BPTree<Backend, Storage>& tree,
                                           TreeNode root) {
  return detail::inspect_decoded_events<detail::bp_micro::BpMicroSubtree<Bits>>(
      tree, root);
}

/** @brief Inspect one bounded DFUDS decoder probe, with the same hit/miss and
 * handle contract as BP. The probe starts at the root's degree run and ends
 * just after the final descendant's zero delimiter (2*n-1 bits for n nodes).
 * Neither word padding nor any following sibling record belongs to the part. */
template <size_t Bits = 256, class Backend, detail::TreeEncodingStorage Storage>
std::optional<TreePart> inspect_event_part(
    const DFUDSTree<Backend, Storage>& tree,
    TreeNode root) {
  return detail::inspect_decoded_events<
      detail::dfuds_micro::DfudsMicroSubtree<Bits>>(tree, root);
}
}  // namespace pixie::tree
