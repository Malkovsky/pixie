#pragma once

#include <pixie/tree/detail/unpack_cache.h>

namespace pixie::benchmarks {
using pixie::detail::tree_unpack::BoundedSubtreeCache;
using pixie::detail::tree_unpack::emit;
using pixie::detail::tree_unpack::subtree;

// Decode a fitting subtree once, then share its child links between inspection
// and descent. One scratch cache belongs to the entire traversal, not each
// recursive frame. Deep non-fitting paths use the stackless navigation
// fallback.
template <bool BP,
          class Tree,
          class Visitor,
          class Inspector,
          size_t BitCapacity = 2048>
class InspectTraversal {
 public:
  InspectTraversal(const Tree& tree,
                   std::span<const uint64_t> words,
                   Visitor& visit,
                   Inspector& inspect)
      : tree_(tree), words_(words), visit_(visit), inspect_(inspect) {}

  // Counts are valid after a successful traversal. Aggregate complete decoded
  // chunks rather than adding a counter operation to every visitor callback.
  size_t node_count = 0, child_count = 0;

  void walk(TreeNode node, size_t depth = 0) {
    if (tree_.is_leaf(node)) {
      ++node_count;
      emit(visit_, node, true);
      return;
    }
    if (cache_.template load<BP>(tree_, words_, 2 * tree_.size() - (BP ? 0 : 1),
                                 node)) {
      cache_.inspect_depth_first(visit_, inspect_);
      node_count += cache_.size();
      child_count += cache_.size() - 1;
      return;
    }
    if (depth == 128) {
      auto events = [&](TreeNode current, bool entering) {
        emit(visit_, current, entering);
        if (entering) {
          ++node_count;
          tree_.for_each_child(current, [&](TreeNode child) {
            inspect_(child);
            ++child_count;
          });
        }
      };
      subtree(tree_, node, events);
      return;
    }
    ++node_count;
    emit(visit_, node, true);
    // Retain a small inspected child list across recursive descent. High-degree
    // nodes fall back to a second enumeration without unbounded stack storage.
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    tree_.for_each_child(node, [&](TreeNode child) {
      inspect_(child);
      if (count < numbers.size()) {
        numbers[count] = child.number;
        positions[count] = child.pos;
      }
      ++count;
    });
    child_count += count;
    if (count <= numbers.size()) {
      for (size_t i = 0; i < count; ++i) {
        walk(TreeNode(numbers[i], positions[i]), depth + 1);
      }
    } else {
      tree_.for_each_child(node,
                           [&](TreeNode child) { walk(child, depth + 1); });
    }
    emit(visit_, node, false);
  }

 private:
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Visitor& visit_;
  Inspector& inspect_;
  BoundedSubtreeCache<BitCapacity> cache_;
};

}  // namespace pixie::benchmarks
