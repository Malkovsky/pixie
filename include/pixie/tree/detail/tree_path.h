#pragma once

#include <stdexcept>
#include <vector>

namespace pixie::tree::detail {

// Keys are unique and strictly decrease on each non-root parent step. The
// larger endpoint therefore cannot be an ancestor of the smaller one.
template <bool AllowDisconnected = false, class Node, class Key, class Parent>
void collect_tree_path(Node u,
                       Node v,
                       std::vector<Node>& out,
                       std::vector<Node>& scratch,
                       Key key,
                       Parent parent) {
  if (&out == &scratch) {
    throw std::invalid_argument("tree_path requires distinct vectors");
  }
  out.clear();
  scratch.clear();
  auto a = key(u), b = key(v);
  while (a != b) {
    if (a > b) {
      const Node next = parent(u);
      if constexpr (AllowDisconnected) {
        if (key(next) == a) {
          out.clear();
          scratch.clear();
          return;
        }
      }
      out.push_back(u);
      u = next;
      a = key(u);
    } else {
      const Node next = parent(v);
      if constexpr (AllowDisconnected) {
        if (key(next) == b) {
          out.clear();
          scratch.clear();
          return;
        }
      }
      scratch.push_back(v);
      v = next;
      b = key(v);
    }
  }
  out.push_back(u);  // Lowest common ancestor, included exactly once.
  out.insert(out.end(), scratch.rbegin(), scratch.rend());
  scratch.clear();
}

}  // namespace pixie::tree::detail
