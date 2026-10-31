#pragma once
#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <span>

namespace pixie::detail::tree_unpack {
// An internal node finishes when the global DFUDS excess reaches its
// start excess minus one. Targets on the ancestor stack are nondecreasing.
template <class Tree, class Visitor>
void dfuds_excess_stream(const Tree& tree,
                         std::span<const uint64_t> words,
                         size_t nbits,
                         Visitor& visit) {
  struct Entry {
    size_t number, pos;
    int64_t target;
  };
  std::array<Entry, 129> stack;
  stack[0] = {0, 0, std::numeric_limits<int64_t>::min()};
  size_t depth = 1, number = 0, start = 0;
  int64_t pending = 0, top = stack[0].target;
  alignas(64) uint8_t offsets[64];
  for (size_t base = 0; base < nbits;
       base += std::min<size_t>(64, nbits - base)) {
    unsigned count = zero_offsets(words[base / 64],
                                  std::min<size_t>(64, nbits - base), offsets);
    for (unsigned i = 0; i < count; ++i) {
      size_t delimiter = base + offsets[i], degree = delimiter - start;
      TreeNode current(number, start);
      int64_t target = pending - 1;
      pending += int64_t(degree) - 1;
      visit(current);
      if (depth <= 128) {
        stack[depth] = {number, start, target};
        depth += degree != 0;
        top = degree ? target : top;
        while (top == pending) {
          auto e = stack[--depth];
          current = TreeNode(e.number, e.pos);
          visit(current);
          top = stack[depth - 1].target;
        }
      } else {
        if (degree) {
          ++depth;
        } else {
          while (depth > 129 && tree.is_last_child(current)) {
            current = tree.parent(current);
            --depth;
            visit(current);
          }
          if (depth <= 129) {
            top = stack[depth - 1].target;
            while (top == pending) {
              auto e = stack[--depth];
              current = TreeNode(e.number, e.pos);
              visit(current);
              top = stack[depth - 1].target;
            }
          }
        }
      }
      ++number;
      start = delimiter + 1;
    }
  }
}

}  // namespace pixie::detail::tree_unpack
