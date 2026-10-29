#pragma once

#include <pixie/tree/path_view.h>

namespace pixie::benchmarks::tree_dfs {
// All traversal state is local; visitors are borrowed and never copied.
template <class Tree, class Visitor>
void cached(const Tree& tree, Visitor& visit) {
  TreePathView<Tree> view(tree);
  bool descending = true;
  while (true) {
    visit(view.node());
    if (descending && !view.is_leaf()) {
      view.move_to_first_child();
      continue;
    }
    descending = false;
    if (view.is_root()) {
      return;
    }
    if (view.is_last_child()) {
      view.move_to_parent();
    } else {
      view.move_to_next_sibling();
      descending = true;
    }
  }
}

}  // namespace pixie::benchmarks::tree_dfs
