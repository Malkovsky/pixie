#pragma once
#include <pixie/tree/detail/unpack_louds.h>

namespace pixie::detail::tree_unpack {
using pixie::detail::tree_unpack::emit;
using pixie::detail::tree_unpack::louds_next;
using pixie::detail::tree_unpack::LoudsLevel;
using pixie::detail::tree_unpack::subtree;
// Cache decoded marker positions for a monotone level cursor. The scalar/SIMD
// dispatch remains inside zero_offsets; callers use the same primitive API.
struct LoudsDecodedLevel {
  size_t position = 0, base = 0;
  unsigned index = 0, count = 0;
  alignas(64) uint8_t offsets[64];
};
inline size_t LoudsDecodedNext(std::span<const uint64_t> words,
                               size_t nbits,
                               LoudsDecodedLevel& cursor,
                               size_t position) {
  if (cursor.index < cursor.count) {
    return cursor.base + cursor.offsets[cursor.index++];
  }
  size_t start = position + 1;
  while (start < nbits) {
    const size_t base = start & ~size_t{63};
    const auto mask = words[base / 64] & (~uint64_t{0} << (start % 64));
    const unsigned valid = std::min<size_t>(64, nbits - base);
    cursor.count =
        pixie::detail::tree_unpack::zero_offsets(~mask, valid, cursor.offsets);
    if (cursor.count) {
      cursor.base = base;
      cursor.index = 1;
      return base + cursor.offsets[0];
    }
    start = base + 64;
  }
  return nbits;
}
// Traversal-local bounded buffers share child-position decoding between
// inspection and descent. Large degrees and depths retain the original path.
template <class Tree,
          class Visitor,
          class Inspector,
          size_t Capacity = 32,
          bool Decoded = false,
          size_t MaxDepth = 128>
struct LoudsPrefetchTraversal {
  static_assert(MaxDepth > 0 && MaxDepth <= 128);
  const Tree& tree;
  std::span<const uint64_t> words;
  size_t nbits;
  Visitor& visit;
  using Level = std::conditional_t<Decoded, LoudsDecodedLevel, LoudsLevel>;
  std::array<Level, MaxDepth + 1> levels{};
  Inspector* inspect = nullptr;
  size_t node_count = 0, child_count = 0;
  std::array<std::array<size_t, Capacity>, MaxDepth> child_ends;

  LoudsPrefetchTraversal(const Tree& t,
                         std::span<const uint64_t> w,
                         Visitor& v,
                         Inspector& i)
      : tree(t),
        words(w),
        nbits(2 * t.size() - (t.empty() ? 0 : 1)),
        visit(v),
        inspect(&i) {}

  size_t next(Level& level, size_t position) {
    if constexpr (Decoded) {
      return LoudsDecodedNext(words, nbits, level, position);
    } else {
      return louds_next(words, nbits, level, position);
    }
  }
  void walk(TreeNode node, size_t depth) {
    const size_t end = next(levels[depth], node.pos);
    levels[depth].position = end;
    walk_known(node, depth, end);
  }

  void walk_known(TreeNode node, size_t depth, size_t end) {
    if (depth == MaxDepth) {
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
    if (degree <= Capacity) {
      // A child's descent only uses buffers at deeper levels. This row stays
      // valid until all siblings finish, including across nested callbacks.
      auto& ends = child_ends[depth];
      const size_t first_position = position;
      for (size_t i = 0; i < degree; ++i) {
        (*inspect)(TreeNode(number + i, position));
        position = next(levels[depth + 1], position);
        ends[i] = position;
      }
      levels[depth + 1].position = position;
      child_count += degree;
      position = first_position;
      for (size_t i = 0; i < degree; ++i) {
        const auto child = TreeNode(number + i, position);
        const size_t child_end = ends[i];
        if (child_end == position + 1) {
          ++node_count;
          emit(visit, child, true);
        } else {
          walk_known(child, depth + 1, child_end);
        }
        position = child_end;
      }
      emit(visit, node, false);
      return;
    }
    if constexpr (!std::is_void_v<Inspector>) {
      // Inspection must finish before descent. Borrow a copy of the next
      // level's monotone cursor, avoiding rank/select for every child list.
      auto cursor = levels[depth + 1];
      auto child_position = position;
      for (size_t i = 0; i < degree; ++i) {
        (*inspect)(TreeNode(number + i, child_position));
        if (i + 1 < degree) {
          child_position = next(cursor, child_position);
        }
      }
    }
    if constexpr (!std::is_void_v<Inspector>) {
      child_count += degree;
    }
    for (size_t i = 0; i < degree; ++i) {
      const auto child = TreeNode(number + i, position);
      const size_t child_end = next(levels[depth + 1], position);
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
}  // namespace pixie::detail::tree_unpack
