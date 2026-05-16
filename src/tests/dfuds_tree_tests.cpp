#include "pixie/dfuds_tree.h"

#include <gtest/gtest.h>
#include <pixie/utils.h>

#include <random>
#include <stack>

using Node = pixie::DFUDSTree::Node;
using pixie::DFUDSTree;

TEST(DfudsTreeTest, Basic) {
  std::vector<std::vector<size_t>> adj = {{0, 1}, {0, 2}, {1, 3}, {2, 4}, {3}};
  size_t tree_size = 5;

  std::vector<uint64_t> dfuds = adj_to_dfuds(tree_size, adj);

  DFUDSTree dfuds_tree(dfuds, 5);
  AdjListTree debug_tree(adj);

  Node cur = dfuds_tree.root();
  AdjListNode debug = debug_tree.root();
  for (size_t i = 0; i < tree_size - 1; i++) {
    EXPECT_EQ(cur, debug);
    cur = dfuds_tree.child(cur, 0);
    debug = debug_tree.child(debug, 0);
  }
  EXPECT_EQ(cur, debug);
}

TEST(DfudsTreeTest, RandomTreeDFS) {
  for (size_t tree_size = 8; tree_size < (1 << 22); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
    adj = dfs_order(tree_size, adj);
    std::vector<uint64_t> dfuds = adj_to_dfuds(tree_size, adj);
    DFUDSTree dfuds_tree(dfuds, tree_size);
    AdjListTree debug_tree(adj);

    std::stack<std::pair<Node, AdjListNode>> st;

    st.push({dfuds_tree.root(), debug_tree.root()});

    while (!st.empty()) {
      auto cur = st.top().first;
      auto debug = st.top().second;
      st.pop();
      EXPECT_EQ(cur, debug);
      EXPECT_EQ(dfuds_tree.parent(cur), debug_tree.parent(debug));

      if (cur.number > 0) {
        EXPECT_EQ(dfuds_tree.is_last_child(cur),
                  debug_tree.is_last_child(debug));
      }
      size_t deg = dfuds_tree.degree(cur);
      EXPECT_EQ(deg, debug_tree.degree(debug));

      if (deg == 0) {
        continue;
      }
      auto child = dfuds_tree.first_child(cur);
      auto debug_child = debug_tree.first_child(debug);
      st.push({child, debug_child});
      for (size_t i = 1; i < deg; i++) {
        child = dfuds_tree.next_sibling(child);
        st.push({child, debug_tree.child(debug, i)});
      }
    }
  }
}

TEST(DfudsTreeTest, RandomTreeBFS) {
  for (size_t tree_size = 8; tree_size < (1 << 22); tree_size <<= 1) {
    std::mt19937_64 rng(42);
    std::vector<std::vector<size_t>> adj = generate_random_tree(tree_size, rng);
    adj = dfs_order(tree_size, adj);
    std::vector<uint64_t> dfuds = adj_to_dfuds(tree_size, adj);
    DFUDSTree dfuds_tree(dfuds, tree_size);
    AdjListTree debug_tree(adj);

    std::queue<std::pair<Node, AdjListNode>> st;

    st.push({dfuds_tree.root(), debug_tree.root()});

    while (!st.empty()) {
      auto cur = st.front().first;
      auto debug = st.front().second;
      st.pop();
      EXPECT_EQ(dfuds_tree.parent(cur), debug_tree.parent(debug));

      if (cur.number > 0) {
        EXPECT_EQ(dfuds_tree.is_last_child(cur),
                  debug_tree.is_last_child(debug));
      }
      size_t deg = dfuds_tree.degree(cur);
      EXPECT_EQ(deg, debug_tree.degree(debug));

      if (deg == 0) {
        continue;
      }
      auto child = dfuds_tree.first_child(cur);
      auto debug_child = debug_tree.first_child(debug);
      st.push({child, debug_child});
      for (size_t i = 1; i < deg; i++) {
        child = dfuds_tree.next_sibling(child);
        st.push({child, debug_tree.child(debug, i)});
      }
    }
  }
}
