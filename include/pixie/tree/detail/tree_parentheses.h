#pragma once

#include <cstddef>
#include <string>

namespace pixie::tree::detail {

template <class Tree, class Root>
std::string tree_parentheses(const Tree& tree,
                             std::size_t roots,
                             Root root,
                             std::size_t start,
                             std::size_t count) {
  std::string output;
  if (!count) {
    return output;
  }
  std::size_t offset = 0;
  for (std::size_t i = 0; i < roots; ++i) {
    auto node = root(i);
    bool entering = true;
    while (true) {
      if (offset++ >= start) {
        output.push_back(entering ? '(' : ')');
        if (output.size() == count) {
          return output;
        }
      }
      if (entering && !tree.is_leaf(node)) {
        node = tree.first_child(node);
        continue;
      }
      if (entering) {
        entering = false;
        continue;
      }
      if (tree.is_root(node)) {
        break;
      }
      if (!tree.is_last_child(node)) {
        node = tree.next_sibling(node);
        entering = true;
      } else {
        node = tree.parent(node);
      }
    }
  }
  return output;
}

}  // namespace pixie::tree::detail
