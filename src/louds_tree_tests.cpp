#include "louds_tree.h"

#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "utils.h"

using pixie::LoudsNode;
using pixie::LoudsTree;

TEST(LoudsTreeTest, Basic) {
  std::vector<std::vector<size_t>> adj = {{0, 1}, {0, 2}, {1, 3}, {2, 4}, {3}};
  size_t tree_size = 5;

  std::vector<uint64_t> louds = adj_to_louds(tree_size, adj);

  LoudsTree louds_tree(louds, 5);
  AdjListTree debug_tree(adj);

  LoudsNode cur = louds_tree.root();
  AdjListNode debug = debug_tree.root();
  for (size_t i = 0; i < tree_size - 1; i++) {
    EXPECT_EQ(cur, debug);
    cur = louds_tree.child(cur, 0);
    debug = debug_tree.child(debug, 0);
  }
  EXPECT_EQ(cur, debug);
}
