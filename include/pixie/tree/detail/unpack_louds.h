#pragma once

#include <pixie/tree/detail/unpack.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <span>
#include <type_traits>

namespace pixie::detail::tree_unpack {
// DFS visits each level in LOUDS order. One monotone word cursor per level
// replaces all but the first select at that level. No per-tree directory.
struct LoudsLevel {
  size_t position = 0;
  size_t base = 0;
  uint64_t mask = 0;
};
inline size_t louds_next(std::span<const uint64_t> words,
                         size_t nbits,
                         LoudsLevel& cursor,
                         size_t position) {
  const size_t start = position + 1;
  if (start == nbits) {
    return nbits;
  }
  const size_t base = start / 64;
  if (cursor.base != base || !cursor.mask) {
    cursor.base = base;
    cursor.mask = words[base] & (~uint64_t{0} << (start % 64));
  }
  while (!cursor.mask) {
    if (++cursor.base >= words.size()) {
      return nbits;
    }
    cursor.mask = words[cursor.base];
  }
  const size_t result = cursor.base * 64 + std::countr_zero(cursor.mask);
  cursor.mask &= cursor.mask - 1;
  return std::min(result, nbits);
}
// Keep invariant arguments in one traversal object rather than passing the
// source span, tree, visitor and level array through every recursive call.
template <class Tree, class Visitor, class Inspector = void>
struct LoudsTraversal {
  const Tree& tree;
  std::span<const uint64_t> words;
  size_t nbits;
  Visitor& visit;
  std::array<LoudsLevel, 129> levels{};
  Inspector* inspect = nullptr;
  size_t node_count = 0, child_count = 0;

  void walk(TreeNode node, size_t depth) {
    const size_t end = louds_next(words, nbits, levels[depth], node.pos);
    levels[depth].position = end;
    walk_known(node, depth, end);
  }

  void walk_known(TreeNode node, size_t depth, size_t end) {
    if (depth == 128) {
      if constexpr (std::is_void_v<Inspector>) {
        subtree(tree, node, visit);
      } else {
        auto events = [&](TreeNode current, bool entering) {
          emit(visit, current, entering);
          if (entering) {
            ++node_count;
            tree.for_each_child(current, [&](TreeNode child) {
              (*inspect)(child);
              ++child_count;
            });
          }
        };
        subtree(tree, node, events);
      }
      return;
    }
    if constexpr (!std::is_void_v<Inspector>) {
      ++node_count;
    }
    emit(visit, node, true);
    const size_t degree = end - node.pos - 1;
    if (!degree) {
      return;
    }
    const size_t number = node.pos + 1 - node.number;
    size_t position = levels[depth + 1].position;
    if (!position) {
      position = tree.first_child(node).pos;
    }
    if constexpr (!std::is_void_v<Inspector>) {
      // Inspection must finish before descent. Borrow a copy of the next
      // level's monotone cursor, avoiding rank/select for every child list.
      auto cursor = levels[depth + 1];
      auto child_position = position;
      for (size_t i = 0; i < degree; ++i) {
        (*inspect)(TreeNode(number + i, child_position));
        if (i + 1 < degree) {
          child_position = louds_next(words, nbits, cursor, child_position);
        }
      }
    }
    if constexpr (!std::is_void_v<Inspector>) {
      child_count += degree;
    }
    for (size_t i = 0; i < degree; ++i) {
      const auto child = TreeNode(number + i, position);
      const size_t child_end =
          louds_next(words, nbits, levels[depth + 1], position);
      levels[depth + 1].position = child_end;
      if (child_end == position + 1) {
        // Leaves need no recursive frame, child inspection, or return event.
        if constexpr (!std::is_void_v<Inspector>) {
          ++node_count;
        }
        emit(visit, child, true);
      } else {
        walk_known(child, depth + 1, child_end);
      }
      position = child_end;
    }
    emit(visit, node, false);
  }
};
template <class Tree, class Visitor>
void louds(const Tree& tree,
           std::span<const uint64_t> words,
           size_t nbits,
           Visitor& visit) {
  LoudsTraversal<Tree, Visitor>{tree, words, nbits, visit}.walk(tree.root(), 0);
}
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
