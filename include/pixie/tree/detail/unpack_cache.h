#pragma once

#include <pixie/tree/detail/unpack.h>
#include <pixie/tree/detail/unpack_bp.h>
#include <pixie/tree/detail/unpack_dfuds.h>

#include <algorithm>

namespace pixie::detail::tree_unpack {
// A bounded, per-traversal decoded subtree. A failed BitCapacity-bit probe
// leaves the cache empty. Only complete subtrees are accepted. Source words
// are borrowed during load; cached handles refer to the original tree.
template <size_t BitCapacity = 256>
class BoundedSubtreeCache {
  static_assert(BitCapacity >= 128 && BitCapacity % 128 == 0);
  static_assert(BitCapacity / 2 < UINT16_MAX);
  static constexpr size_t kNodes = BitCapacity / 2;
  // Keep each handle beside its navigation metadata. Callback emission copies
  // the ready-made handle instead of reconstructing it from separate arrays.
  struct Record {
    size_t number, position;
    uint16_t parent, next, degree;
  };

 public:
  size_t size() const { return count_; }
  bool contains(TreeNode node) const {
    return node.number >= first_number_ && node.number - first_number_ < count_;
  }

  template <bool BP, class Tree>
  [[gnu::noinline]] bool load(const Tree&,
                              std::span<const uint64_t> words,
                              size_t nbits,
                              TreeNode root) {
    count_ = 0;
    std::array<uint64_t, BitCapacity / 64> local{};
    const size_t available = std::min<size_t>(BitCapacity, nbits - root.pos);
    for (size_t offset = 0; offset < available; offset += 64) {
      local[offset / 64] = word_at(words, root.pos + offset);
    }
    size_t length = 0;
    int target = BP ? 0 : -1;
    for (size_t offset = 0; offset < available; offset += 128) {
      int delta;
      const size_t found =
          forward_search_128(local.data() + offset / 64, target, 0, &delta);
      if (found < std::min<size_t>(128, available - offset)) {
        length = offset + found + 1;
        break;
      }
      target -= delta;
    }
    if (!length) {
      return false;
    }
    first_number_ = root.number;
    std::array<uint16_t, kNodes> ancestors, last;
    size_t record_start = 0;
    alignas(64) uint8_t offsets[64];
    for (size_t base = 0; base < length; base += 64) {
      const unsigned valid = std::min<size_t>(64, length - base);
      const unsigned records = zero_offsets(
          BP ? ~local[base / 64] : local[base / 64], valid, offsets);
      for (unsigned j = 0; j < records; ++j) {
        const size_t bit = base + offsets[j];
        const size_t position = BP ? bit : record_start;
        const size_t i = count_++;
        const uint16_t degree = BP ? 0 : bit - record_start;
        size_t parent = 0;
        if constexpr (BP) {
          const size_t depth = 2 * i - position;
          if (i) {
            parent = ancestors[depth - 1];
          }
          ancestors[depth] = static_cast<uint16_t>(i);
        } else {
          // DFUDS excess before this record is the top pending-child slot.
          // A degree-d record replaces it with d slots owned by this node.
          // Across the subtree, exactly n-1 slots are filled: no ancestor
          // completion loop or per-level remaining-child counters are needed.
          const size_t pending = position - 2 * i;
          if (i) {
            parent = ancestors[pending];
          }
          fill_parent_slots(ancestors.data() + pending, degree,
                            static_cast<uint16_t>(i));
          // The immediate child consumes the highest newly filled slot next.
          // Keep its scalar store separate from the masked fill.
          // For leaves the slot is discarded; it remains inside the scratch.
          ancestors[pending + degree - (degree != 0)] =
              static_cast<uint16_t>(i);
          record_start = bit + 1;
        }
        nodes_[i] = {root.number + i, root.pos + position,
                     static_cast<uint16_t>(parent), kNodes, degree};
        last[i] = kNodes;
        if (i) {
          // The dummy record absorbs the first-child write, avoiding a
          // data-dependent branch for every new sibling list.
          nodes_[last[parent]].next = static_cast<uint16_t>(i);
          last[parent] = static_cast<uint16_t>(i);
          if constexpr (BP) {
            ++nodes_[parent].degree;
          }
        }
      }
    }
    return true;
  }

  size_t degree(TreeNode node) const {
    return nodes_[node.number - first_number_].degree;
  }
  TreeNode child(TreeNode node, size_t index) const {
    size_t child = node.number - first_number_ + 1;
    while (index--) {
      child = nodes_[child].next;
    }
    return handle(child);
  }
  template <class Visitor>
  size_t children(TreeNode node, Visitor&& visit) const {
    const size_t index = node.number - first_number_;
    const auto count = nodes_[index].degree;
    if (!count) {
      return 0;
    }
    for (size_t child = index + 1; child != kNodes;
         child = nodes_[child].next) {
      visit(handle(child));
    }
    return count;
  }

  // Preorder records are contiguous, so inspection and DFS events need no
  // recursive calls or explicit path stack within a decoded subtree.
  template <class Visitor, class Inspector>
  [[gnu::noinline]] void inspect_depth_first(Visitor& visit,
                                             Inspector& inspect) const {
    const size_t count = count_;
    for (size_t i = 0; i < count; ++i) {
      emit(visit, handle(i), true);
      if (nodes_[i].degree) {
        for (size_t child = i + 1; child != kNodes;
             child = nodes_[child].next) {
          inspect(handle(child));
        }
      } else {
        size_t completed = i;
        while (completed && nodes_[completed].next == kNodes) {
          completed = nodes_[completed].parent;
          emit(visit, handle(completed), false);
        }
      }
    }
  }

 private:
  TreeNode handle(size_t index) const {
    return TreeNode(nodes_[index].number, nodes_[index].position);
  }
  size_t first_number_ = 0, count_ = 0;
  std::array<Record, kNodes + 1> nodes_;
};
using SubtreeCache = BoundedSubtreeCache<>;

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
    if constexpr (requires {
                    tree_.navigate(node).inspect_decoded(visit_, inspect_);
                  }) {
      const auto count = tree_.navigate(node).inspect_decoded(visit_, inspect_);
      if (count) {
        node_count += count;
        child_count += count - 1;
        return;
      }
    } else if (cache_.template load<BP>(
                   tree_, words_, 2 * tree_.size() - (BP ? 0 : 1), node)) {
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

}  // namespace pixie::detail::tree_unpack
