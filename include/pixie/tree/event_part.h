#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/tree_part.h>

namespace pixie::tree {
namespace detail {
template <class Decoder, class Tree>
std::optional<TreePart> inspect_decoded_events(const Tree& tree,
                                               TreeNode root,
                                               TreePart::Kind kind) {
  Decoder decoder;
  if (tree.empty() || !decoder.load(tree.encoding_words(),
                                    Decoder::bit_size(tree.size()), root)) {
    return std::nullopt;
  }
  using Type = TreePart::Event::Type;
  std::vector<TreePart::Node> nodes(decoder.size());
  std::vector<TreePart::Event> events;
  events.reserve(3 * nodes.size() - 2);
  size_t parent = 0;
  auto visit = [&](TreeNode n, bool entering) {
    const auto i = n.number - root.number;
    events.push_back({entering ? Type::kEnter : Type::kReturn, i});
    if (entering) {
      nodes[i].source_id = n.number;
      nodes[i].bit_position = n.pos;
      parent = i;
    }
  };
  // All inspect callbacks immediately follow their parent's entry, before
  // the first child is visited. No source navigation or stack is needed.
  auto inspect = [&](TreeNode child) {
    const auto i = child.number - root.number;
    nodes[parent].children.push_back(i);
    ++nodes[parent].degree;
    events.push_back({Type::kInspect, i});
  };
  decoder.walk(visit, inspect);
  return TreePart(
      kind, std::move(nodes), {0}, std::move(events),
      TreePart::Range{root.pos, root.pos + Decoder::bit_size(decoder.size())});
}
}  // namespace detail

/** @brief Inspect one bounded BP decoder probe as a complete event part.
 * @details root must be a valid handle from tree (ignored on an empty tree).
 * Returns nullopt for an empty tree or when the subtree does not close within
 * Bits bits; no prefix or navigation fallback is returned. On a hit the part
 * contains the opening bit through its matching closing bit, exactly 2*n bits
 * for n local nodes. The root's parent and siblings are outside the part.
 * Bits is a multiple of 128 in [128,32768]. */
template <size_t Bits = 256, class Backend, detail::TreeEncodingStorage Storage>
std::optional<TreePart> inspect_event_part(const BPTree<Backend, Storage>& tree,
                                           TreeNode root) {
  return detail::inspect_decoded_events<detail::bp_micro::BpMicroSubtree<Bits>>(
      tree, root, TreePart::Kind::kBpEvents);
}

/** @brief Inspect one bounded DFUDS decoder probe, with the same hit/miss and
 * handle contract as BP. The part starts at the root's degree run and ends
 * just after the final descendant's zero delimiter (2*n-1 bits for n nodes).
 * Neither word padding nor any following sibling record belongs to the part. */
template <size_t Bits = 256, class Backend, detail::TreeEncodingStorage Storage>
std::optional<TreePart> inspect_event_part(
    const DFUDSTree<Backend, Storage>& tree,
    TreeNode root) {
  return detail::inspect_decoded_events<
      detail::dfuds_micro::DfudsMicroSubtree<Bits>>(
      tree, root, TreePart::Kind::kDfudsEvents);
}
}  // namespace pixie::tree
