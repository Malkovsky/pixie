#pragma once

#include <pixie/tree/detail/unpack.h>

#include <algorithm>
#include <array>
#include <span>

namespace pixie::detail::tree_unpack {
// Decoder supplies bit_size(nodes), load_scratch_bytes(), load(), walk(), and
// size(). Keep encoding-specific decoding separate from the shared traversal.
template <class Decoder, class Tree, class Visit, class Inspect>
class MicroTraversal {
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Visit& visit_;
  Inspect& inspect_;
  size_t max_depth_;
  Decoder cache_;

 public:
  size_t nodes = 0, children = 0;
  // Includes load-local arrays and all simultaneously live recursive child
  // buffers. Compiler call-frame overhead and callback state are separate.
  static constexpr size_t scratch_bytes(size_t depth = 128) {
    return sizeof(MicroTraversal) +
           std::min(depth, size_t{128}) * 64 * sizeof(size_t) +
           Decoder::load_scratch_bytes() + 1024;
  }
  MicroTraversal(const Tree& tree,
                 std::span<const uint64_t> words,
                 Visit& visit,
                 Inspect& inspect,
                 size_t max_depth = 128)
      : tree_(tree),
        words_(words),
        visit_(visit),
        inspect_(inspect),
        max_depth_(std::min(max_depth, size_t{128})) {}
  void walk(TreeNode node, size_t depth = 0) {
    if (tree_.is_leaf(node)) {
      ++nodes;
      emit(visit_, node, true);
      return;
    }
    if (cache_.load(words_, Decoder::bit_size(tree_.size()), node)) {
      cache_.walk(visit_, inspect_);
      nodes += cache_.size();
      children += cache_.size() - 1;
      return;
    }
    if (depth == max_depth_) {
      walk_uncached(node);
      return;
    }
    ++nodes;
    emit(visit_, node, true);
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    tree_.for_each_child(node, [&](auto c) {
      inspect_(c);
      if (count < 32) {
        numbers[count] = c.number;
        positions[count] = c.pos;
      }
      ++count;
    });
    children += count;
    if (count <= 32) {
      for (size_t i = 0; i < count; ++i) {
        walk({numbers[i], positions[i]}, depth + 1);
      }
    } else {
      tree_.for_each_child(node, [&](auto c) { walk(c, depth + 1); });
    }
    emit(visit_, node, false);
  }

 private:
  // Keep decoding below the child-buffer depth limit, without another
  // recursive child array. Navigation views supply ascent/sibling state.
  void walk_uncached(TreeNode root) {
    auto view = tree_.navigate(root);
    bool descending = true;
    while (true) {
      if (descending && !view.is_leaf() &&
          cache_.load(words_, Decoder::bit_size(tree_.size()), view.node())) {
        cache_.walk(visit_, inspect_);
        nodes += cache_.size();
        children += cache_.size() - 1;
      } else {
        emit(visit_, view.node(), descending);
        if (descending) {
          ++nodes;
          tree_.for_each_child(view.node(), [&](auto c) {
            ++children;
            inspect_(c);
          });
          if (!view.is_leaf()) {
            view.reset(view.first_child());
            continue;
          }
        }
      }
      descending = false;
      if (view.node().number == root.number) {
        return;
      }
      if (view.is_last_child()) {
        view.reset(view.parent());
      } else {
        view.reset(view.next_sibling());
        descending = true;
      }
    }
  }
};
}  // namespace pixie::detail::tree_unpack
