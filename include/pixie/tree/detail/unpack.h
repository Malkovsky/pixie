#pragma once

#include <pixie/tree.h>

namespace pixie::detail::tree_unpack {
// Existing one-argument visitors keep the DFS event sequence; consumers that
// inspect children can additionally distinguish entry from return events.
template <class Visitor>
void emit(Visitor& visit, TreeNode node, bool entering) {
  if constexpr (requires { visit(node, entering); }) {
    visit(node, entering);
  } else {
    visit(node);
  }
}
// Bounded-depth decoders use ordinary navigation for an uncached subtree.
template <class Tree, class Visitor>
void subtree(const Tree& tree, TreeNode root, Visitor& visit) {
  auto view = tree.navigate(root);
  bool descending = true;
  while (true) {
    emit(visit, view.node(), descending);
    if (descending && !view.is_leaf()) {
      view.reset(view.first_child());
      continue;
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
}  // namespace pixie::detail::tree_unpack
