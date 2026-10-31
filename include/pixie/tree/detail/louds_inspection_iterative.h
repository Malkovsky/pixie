#pragma once

#include <pixie/tree/detail/louds_inspection.h>

namespace pixie::detail::tree_unpack {
// Explicit traversal stack with the same per-level child buffers. Large child
// lists and depths retain the validated fallback, outside the hot loop.
template <class Tree,
          class Visitor,
          class Inspector,
          bool Decoded = true,
          size_t MaxDepth = 128,
          size_t Capacity = 32>
struct LoudsIterativeTraversal : LoudsPrefetchTraversal<Tree,
                                                        Visitor,
                                                        Inspector,
                                                        Capacity,
                                                        Decoded,
                                                        MaxDepth> {
  using Base = LoudsPrefetchTraversal<Tree,
                                      Visitor,
                                      Inspector,
                                      Capacity,
                                      Decoded,
                                      MaxDepth>;
  using Base::Base;
  struct Frame {
    size_t node_number, node_position, number, position, remaining, index;
  };
  std::array<Frame, MaxDepth> frames;

  static constexpr size_t scratch_bytes() {
    return sizeof(LoudsIterativeTraversal) +
           (MaxDepth + 1) * sizeof(typename Base::Level) + 1024;
  }

  [[gnu::cold, gnu::noinline]] void fallback(TreeNode node,
                                             size_t depth,
                                             size_t end) {
    Base::walk_known(node, depth, end);
  }

  void walk(TreeNode node, size_t depth) {
    const size_t root_depth = depth;
    size_t end = this->next(this->levels[depth], node.pos);
    this->levels[depth].position = end;
    Frame frame;
  enter: {
    const size_t degree = end - node.pos - 1;
    if (depth == MaxDepth || degree > Capacity) {
      fallback(node, depth, end);
      goto ascend;
    }
    ++this->node_count;
    emit(this->visit, node, true);
    if (!degree) {
      goto ascend;
    }
    const size_t number = node.pos + 1 - node.number;
    size_t position = this->levels[depth + 1].position;
    if (!position) {
      position = this->tree.first_child(node).pos;
    }
    frame = {node.number, node.pos, number, position, degree, 0};
    auto& ends = this->child_ends[depth];
    for (size_t i = 0; i < degree; ++i) {
      (*this->inspect)(TreeNode(number + i, position));
      position = this->next(this->levels[depth + 1], position);
      ends[i] = position;
    }
    this->levels[depth + 1].position = position;
    this->child_count += degree;
  }
  children:
    while (frame.remaining) {
      const auto child = TreeNode(frame.number++, frame.position);
      const size_t child_end = this->child_ends[depth][frame.index++];
      --frame.remaining;
      frame.position = child_end;
      if (child_end == child.pos + 1) {
        ++this->node_count;
        emit(this->visit, child, true);
      } else {
        frames[depth++] = frame;
        node = child;
        end = child_end;
        goto enter;
      }
    }
    emit(this->visit, TreeNode(frame.node_number, frame.node_position), false);
  ascend:
    if (depth == root_depth) {
      return;
    }
    frame = frames[--depth];
    goto children;
  }
};
}  // namespace pixie::detail::tree_unpack
