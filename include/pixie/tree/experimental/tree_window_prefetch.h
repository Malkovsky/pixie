#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/unpack.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>

namespace pixie::experimental {
/**
 * @brief Experimental fixed-bit-window BP traversal with child inspection.
 * @tparam Bits Window bits, a multiple of 64 in [128, 32768].
 * @details Source is the immutable BP encoding of tree, with exactly 2*size()
 * logical bits; no padding is required. The tree/source outlive the call.
 * Before descending from each node, inspect receives all children in order.
 * Visit receives entry and internal-node return events (optionally a bool
 * distinguishing them, as in the internal unpack kernels). Both callbacks are
 * borrowed, may throw, and may perform nested traversals. Storage is bounded
 * and local to this call. Complete child lists use decoded local links;
 * lists crossing a window boundary use representation navigation. Deep return
 * paths beyond 128 levels use representation parent queries. Empty trees emit
 * nothing. The source must match tree; a mismatch has undefined behavior.
 */
template <size_t Bits = 2048, class Tree, class Visitor, class Inspector>
void BPWindowInspect(const Tree& tree,
                     std::span<const uint64_t> words,
                     Visitor&& visit,
                     Inspector&& inspect) {
  static_assert(Bits >= 128 && Bits <= 32768 && Bits % 64 == 0);
  // A window may consist entirely of opens, so capacity is Bits, not Bits/2.
  constexpr uint16_t missing = Bits;
  std::array<uint16_t, Bits + 1> position, next, last;
  std::array<uint16_t, 2 * Bits + 1> ancestors;
  std::array<bool, Bits> crossing;
  struct Handle {
    size_t number, position;
  };
  std::array<Handle, 128> path;
  size_t active_depth = 0, number = 0;
  TreeNode current(0, 0);
  auto finish_to = [&](size_t depth) {
    while (active_depth > depth) {
      --active_depth;
      if (active_depth < path.size()) {
        current =
            TreeNode(path[active_depth].number, path[active_depth].position);
      } else {
        current = tree.parent(current);
      }
      pixie::detail::tree_unpack::emit(visit, current, false);
    }
  };
  const size_t nbits = 2 * tree.size();
  alignas(64) uint8_t offsets[64];
  for (size_t base = 0; base < nbits; base += Bits) {
    const size_t length = std::min(Bits, nbits - base);
    size_t count = 0;
    int minimum_depth = Bits;
    last[missing] = missing;
    for (size_t offset = 0; offset < length; offset += 64) {
      const auto n = pixie::detail::tree_unpack::zero_offsets(
          ~words[(base + offset) / 64], std::min<size_t>(64, length - offset),
          offsets);
      for (unsigned j = 0; j < n; ++j) {
        const auto bit = offset + offsets[j];
        const int depth = 2 * int(count) - int(bit);
        const auto parent =
            depth > minimum_depth ? ancestors[Bits + depth - 1] : missing;
        minimum_depth = std::min(minimum_depth, depth);
        position[count] = bit;
        next[count] = missing;
        last[count] = missing;
        next[last[parent]] = count;
        last[parent] = count;
        ancestors[Bits + depth] = count;
        crossing[count] = false;
        ++count;
      }
    }
    const int final_depth = 2 * int(count) - int(length);
    // Only still-open ancestors lack a complete child list. Everything else
    // can use local links, including subtrees rooted in the middle of a window.
    for (int depth = minimum_depth; depth < final_depth; ++depth) {
      crossing[ancestors[Bits + depth]] = true;
    }
    for (size_t i = 0; i < count; ++i) {
      const size_t pos = base + position[i];
      const size_t depth = 2 * (number + i) - pos;
      finish_to(depth);
      current = TreeNode(number + i, pos);
      pixie::detail::tree_unpack::emit(visit, current, true);
      const size_t following = pos + 1;
      const bool internal = (words[following / 64] >> (following % 64)) & 1;
      if (!internal) {
        continue;
      }
      if (active_depth < path.size()) {
        path[active_depth] = {current.number, current.pos};
      }
      ++active_depth;
      if (crossing[i]) {
        tree.for_each_child(current, inspect);
      } else {
        for (size_t child = i + 1; child != missing; child = next[child]) {
          inspect(TreeNode(number + child, base + position[child]));
        }
      }
    }
    number += count;
  }
  finish_to(0);
}
}  // namespace pixie::experimental
