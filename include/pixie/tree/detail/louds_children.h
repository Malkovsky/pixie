#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <algorithm>
#include <span>

namespace pixie::detail::tree_unpack {
// Decode consecutive child handles in windows of at most 256 source bits.
// A long unary run falls back to indexed child access at the next window.
// first/count must describe a valid consecutive child range of node.
template <class Tree, class Visitor>
void louds_children(const Tree& tree,
                    std::span<const uint64_t> words,
                    size_t nbits,
                    TreeNode node,
                    size_t first,
                    size_t count,
                    Visitor&& visit) {
  alignas(64) uint8_t offsets[64];
  size_t done = 0;
  while (done < count) {
    const auto start = tree.child(node, first + done);
    const size_t length = std::min<size_t>(256, nbits - start.pos);
    const size_t number = start.number;
    size_t decoded = 0;
    for (size_t offset = 0; offset < length && done < count; offset += 64) {
      const unsigned valid = std::min<size_t>(64, length - offset);
      const unsigned n =
          zero_offsets(~word_at(words, start.pos + offset), valid, offsets);
      for (unsigned i = 0; i < n && done < count; ++i, ++done, ++decoded) {
        visit(TreeNode(number + decoded, start.pos + offset + offsets[i]));
      }
    }
  }
}
}  // namespace pixie::detail::tree_unpack
