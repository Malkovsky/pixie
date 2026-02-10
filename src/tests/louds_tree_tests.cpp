#include <gtest/gtest.h>
#include <pixie/louds_tree.h>
#include <pixie/utils.h>

#include <iostream>
#include <numeric>
#include <random>
#include <stack>

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

TEST(LoudsTreeTest, RandomTreeDFS) {
  for (size_t tree_size = 8; tree_size < (1 << 22); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
    adj = bfs_order(tree_size, adj);
    std::vector<uint64_t> louds = adj_to_louds(tree_size, adj);
    LoudsTree louds_tree(louds, tree_size);
    AdjListTree debug_tree(adj);

    std::stack<std::pair<LoudsNode, AdjListNode>> st;

    st.push({louds_tree.root(), debug_tree.root()});

    while (!st.empty()) {
      auto cur = st.top().first;
      auto debug = st.top().second;
      st.pop();
      EXPECT_EQ(cur, debug);
      EXPECT_EQ(louds_tree.parent(cur), debug_tree.parent(debug));
      EXPECT_EQ(louds_tree.is_last_child(cur), debug_tree.is_last_child(debug));

      if (!debug_tree.is_last_child(debug)) {
        EXPECT_EQ(louds_tree.next_sibling(cur), debug_tree.next_sibling(debug));
      }

      size_t deg = louds_tree.degree(cur);
      EXPECT_EQ(deg, debug_tree.degree(debug));

      for (size_t i = 0; i < deg; i++) {
        st.push({louds_tree.child(cur, i), debug_tree.child(debug, i)});
      }
    }
  }
}

TEST(LoudsTreeTest, RandomTreeDFSWithZeroExtraMemory) {
  for (size_t tree_size = 8; tree_size < (1 << 22); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
    adj = bfs_order(tree_size, adj);
    std::vector<uint64_t> louds = adj_to_louds(tree_size, adj);
    LoudsTree louds_tree(louds, tree_size);
    AdjListTree debug_tree(adj);

    LoudsNode cur = louds_tree.root();
    AdjListNode debug = debug_tree.root();

    bool above = 1;
    size_t cnt_visited = 1;
    while (true) {
      if (above) {
        ASSERT_EQ(louds_tree.is_leaf(cur), debug_tree.is_leaf(debug));
        if (louds_tree.is_leaf(cur)) {
          above = 0;
        } else {
          cur = louds_tree.first_child(cur);
          debug = debug_tree.first_child(debug);
          ASSERT_EQ(cur, debug);
          cnt_visited++;
        }
      } else {
        ASSERT_EQ(louds_tree.is_last_child(cur),
                  debug_tree.is_last_child(debug));
        if (louds_tree.is_last_child(cur)) {
          cur = louds_tree.parent(cur);
          debug = debug_tree.parent(debug);
          ASSERT_EQ(cur, debug);
          if (louds_tree.is_root(cur)) {
            break;
          }
        } else {
          cur = louds_tree.next_sibling(cur);
          debug = debug_tree.next_sibling(debug);
          ASSERT_EQ(cur, debug);
          above = 1;
          cnt_visited++;
        }
      }
    }
    ASSERT_EQ(tree_size, cnt_visited);
  }
}
