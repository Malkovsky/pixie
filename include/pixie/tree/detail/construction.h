#pragma once

#include <pixie/tree/detail/encoding_storage.h>

#include <concepts>
#include <type_traits>
#include <vector>

namespace pixie::detail {

// Sources are valid ordered trees. Handles remain valid throughout traversal;
// first_child and next_sibling are called only when that node exists.
template <class Tree>
concept OrderedTreeSource =
    requires(const Tree& tree, decltype(tree.root()) node) {
      { tree.size() } -> std::convertible_to<size_t>;
      { tree.degree(node) } -> std::convertible_to<size_t>;
      { tree.first_child(node) } -> std::same_as<decltype(node)>;
      { tree.next_sibling(node) } -> std::same_as<decltype(node)>;
    };

// Shared DFS event stream. Each node generates enter(node, degree), then its
// children in source order, then leave(node). Scratch is O(height) handles.
template <OrderedTreeSource Tree, class Enter, class Leave>
void generic_construction_dfs(const Tree& tree, Enter enter, Leave leave) {
  if (!tree.size()) {
    return;
  }
  using Node = decltype(tree.root());
  struct Frame {
    Node node;
    Node child;
    size_t remaining;
  };
  std::vector<Frame> stack;
  auto visit = [&](Node node) {
    const auto degree = tree.degree(node);
    enter(node, degree);
    if (degree) {
      stack.push_back({node, tree.first_child(node), degree});
    } else {
      leave(node);
    }
  };
  visit(tree.root());
  while (!stack.empty()) {
    auto& frame = stack.back();
    if (!frame.remaining) {
      leave(frame.node);
      stack.pop_back();
      continue;
    }
    const auto child = frame.child;
    if (--frame.remaining) {
      frame.child = tree.next_sibling(child);
    }
    visit(child);
  }
}

// Concrete encodings specialize this producer without depending on any
// destination encoding. References and user-defined sources use navigation.
template <class Tree>
struct ConstructionTraversal {
  template <class Enter, class Leave>
  static void walk(const Tree& tree, Enter& enter, Leave& leave) {
    generic_construction_dfs(tree, enter, leave);
  }
  static constexpr size_t scratch_bytes() { return 0; }
};

// Degree-only destinations do not observe exit events. Naming this callback
// lets decoded producers skip exit normalization instead of doing invisible
// work on every return from a subtree.
struct IgnoreConstructionExit {
  template <class Node>
  void operator()(Node) const {}
};

// Normalize inspected DFS kernels: entry is followed by immediate children;
// return events exist only for nonleaves. Delay entry until degree is known,
// and synthesize the missing leaf exit. No per-node degree array is needed.
template <class Tree, class Enter, class Leave, class InspectWalk>
void inspected_construction_dfs(const Tree& tree,
                                Enter& enter,
                                Leave& leave,
                                InspectWalk walk) {
  auto pending_node = tree.root();
  bool pending = false;
  size_t degree = 0;
  auto flush = [&] {
    if (pending) {
      enter(pending_node, degree);
      if (!degree) {
        leave(pending_node);
      }
      pending = false;
    }
  };
  auto visit = [&](auto node, bool entering) {
    if constexpr (std::same_as<Leave, IgnoreConstructionExit>) {
      if (!entering) {
        return;
      }
    }
    flush();
    if (entering) {
      pending_node = node;
      degree = 0;
      pending = true;
    } else {
      leave(node);
    }
  };
  auto child = [&](auto) { ++degree; };
  walk(visit, child);
  flush();
}

// Common contract: preorder enter(node, degree), ordered children, leave(node)
// for every node, including leaves. Never query root on an empty source.
template <OrderedTreeSource Tree, class Enter, class Leave>
void construction_dfs(const Tree& tree, Enter enter, Leave leave) {
  if (tree.size()) {
    ConstructionTraversal<Tree>::walk(tree, enter, leave);
  }
}

// DFS order restricted to one depth is BFS order. Staging packed records per
// depth avoids a queue of node handles while preserving LOUDS sibling order.
struct ConstructionLevel {
  std::vector<uint64_t> words;
  size_t bits = 0;
  void append(size_t degree) {
    const auto position = bits;
    bits += degree + 1;
    words.resize(bits / 64 + (bits % 64 != 0), 0);
    words[position / 64] |= uint64_t{1} << (position % 64);
  }
};

// Builds directly into the destination allocation. No source node number or
// bit position is copied. Source topology and size must agree.
template <TreeEncodingKind Kind, OrderedTreeSource Tree>
AlignedStorage construct_tree_encoding(const Tree& tree) {
  const auto bits =
      tree_encoding_bits(tree.size(), Kind == TreeEncodingKind::kBP);
  const auto words = bits / 64 + (bits % 64 != 0);
  if (words > std::numeric_limits<size_t>::max() / 64) {
    throw std::length_error("Tree encoding allocation is too large");
  }
  AlignedStorage output(words * 64);
  auto data = output.writable_words64();
  size_t position = 0;
  auto one = [&] {
    data[position / 64] |= uint64_t{1} << (position % 64);
    ++position;
  };
  if constexpr (Kind == TreeEncodingKind::kBP) {
    construction_dfs(
        tree, [&](auto, size_t) { one(); }, [&](auto) { ++position; });
  } else if constexpr (Kind == TreeEncodingKind::kDFUDS) {
    construction_dfs(
        tree,
        [&](auto, size_t degree) {
          for (size_t i = 0; i < degree; ++i) {
            one();
          }
          ++position;
        },
        IgnoreConstructionExit{});
  } else {
    std::vector<ConstructionLevel> levels;
    size_t depth = 0;
    construction_dfs(
        tree,
        [&](auto, size_t degree) {
          if (levels.size() <= depth) {
            levels.resize(depth + 1);
          }
          levels[depth++].append(degree);
        },
        [&](auto) { --depth; });
    size_t base = 0;
    for (const auto& level : levels) {
      for (size_t i = 0; i < level.words.size(); ++i) {
        const size_t pos = base + i * 64;
        data[pos / 64] |= level.words[i] << (pos % 64);
        if (pos % 64 && pos / 64 + 1 < data.size()) {
          data[pos / 64 + 1] |= level.words[i] >> (64 - pos % 64);
        }
      }
      base += level.bits;
    }
  }
  return output;
}

}  // namespace pixie::detail
