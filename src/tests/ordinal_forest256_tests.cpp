#include <gtest/gtest.h>
#include <pixie/tree/forest_cover.h>

#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using Tree = pixie::tree::ForestCoverTree;
using Context = pixie::tree::TreeTraversalContext<Tree>;
using Index = pixie::tree::RootNavigationIndex<Tree>;
using Node = Tree::Node;
constexpr size_t kSlot = Context::slot_bytes();

template <class T>
concept CanRead = std::is_constructible_v<Context, T, size_t>;
template <class T>
concept CanBind = requires(Context& r, T&& t) { r.bind(std::forward<T>(t)); };
static_assert(CanRead<const Tree&> && CanBind<const Tree&>);
static_assert(!CanRead<Tree> && !CanRead<const Tree> && !CanBind<Tree> &&
              !CanBind<const Tree>);
static_assert(!std::is_move_constructible_v<Tree> &&
              !std::is_copy_assignable_v<Tree>);
static_assert(pixie::tree::TreeNavigationView<Tree::NodeView, Node>);

// Opaque source handles; adjacency and insertion IDs are independent of unit
// IDs.
struct Source {
  struct Handle {
    size_t bits;
  };
  explicit Source(std::vector<size_t> parents)
      : n(parents.size()), children(n), next(n) {
    for (size_t i = 1; i < n; ++i) {
      auto& c = children[parents[i]];
      if (!c.empty()) {
        next[c.back()] = i;
      }
      c.push_back(i);
    }
  }
  size_t size() const { return n; }
  Handle root() const { return {~size_t{0}}; }
  size_t degree(Handle h) const { return children.at(~h.bits).size(); }
  Handle first_child(Handle h) const { return {~children.at(~h.bits).at(0)}; }
  Handle next_sibling(Handle h) const { return {~next.at(~h.bits)}; }
  size_t n;
  std::vector<std::vector<size_t>> children;
  std::vector<size_t> next;
};

void Check(const Source& s, const Tree& t, size_t budget) {
  ASSERT_EQ(t.size(), s.size());
  Context r(t, budget);
  std::vector<std::pair<size_t, Node>> pending{{0, t.root()}};
  std::vector<bool> seen(t.size());
  size_t visited = 0;
  while (!pending.empty()) {
    auto [v, h] = pending.back();
    pending.pop_back();
    const auto id = t.number(h);
    ASSERT_LT(id, t.size());
    ASSERT_FALSE(seen[id]);
    seen[id] = true;
    ++visited;
    ASSERT_EQ(t.find_by_number(id), h);
    const auto& expected = s.children[v];
    ASSERT_EQ(r.degree(h), expected.size());
    ASSERT_EQ(r.is_leaf(h), expected.empty());
    size_t j = 0;
    Node previous{};
    r.for_each_child(h, [&](Node c) {
      ASSERT_LT(j, expected.size());
      EXPECT_EQ(r.child(h, j), c);
      EXPECT_EQ(j ? r.next_sibling(previous) : r.first_child(h), c);
      EXPECT_EQ(r.is_last_child(c), j + 1 == expected.size());
      EXPECT_EQ(r.parent(c), h);
      pending.emplace_back(expected[j++], c);
      previous = c;
    });
    ASSERT_EQ(j, expected.size());
  }
  EXPECT_EQ(visited, t.size());
  EXPECT_EQ(r.parent(t.root()), t.root());
  EXPECT_TRUE(r.is_last_child(t.root()));
}

TEST(OrdinalForest256, EmptyAndInvalidExtents) {
  Source s({});
  Tree t(s);
  EXPECT_EQ(t.unit_count(), 0u);
  EXPECT_THROW(t.root(), std::out_of_range);
  EXPECT_THROW(t.find_by_number(0), std::out_of_range);
  s.n = SIZE_MAX;
  EXPECT_THROW({ Tree bad(s); }, std::length_error);
  s.n = UINT32_MAX;
  EXPECT_THROW({ Tree bad(s); }, std::length_error);
  Source edge({0, 0});
  edge.n = 1;
  EXPECT_THROW({ Tree bad(edge); }, std::invalid_argument);
  edge.n = 3;
  EXPECT_THROW({ Tree bad(edge); }, std::invalid_argument);
}

TEST(OrdinalForest256, ExhaustiveSmallOrderedTrees) {
  size_t shapes = 0;
  for (size_t n = 1; n <= 8; ++n) {
    const unsigned bits = 2 * n - 2;
    for (unsigned mask = 0; mask < (1u << bits); ++mask) {
      std::vector<size_t> p{0}, stack{0};
      bool valid = true;
      for (unsigned i = 0; i < bits; ++i) {
        if (mask >> i & 1) {
          p.push_back(stack.back());
          stack.push_back(p.size() - 1);
        } else if (stack.size() > 1) {
          stack.pop_back();
        } else {
          valid = false;
          break;
        }
      }
      if (valid && stack.size() == 1 && p.size() == n) {
        Source s(std::move(p));
        Check(s, Tree(s), kSlot);
        ++shapes;
      }
    }
  }
  EXPECT_EQ(shapes, 626u);
}

TEST(ForestCover, NodeViewEnumerationSurvivesReset) {
  // Small root: local run. Wide root: portal spanning several child units.
  for (size_t n : {7u, 1025u}) {
    Tree tree{Source(std::vector<size_t>(n))};
    auto view = tree.navigate(tree.root());
    auto original = view;
    size_t k = 0;
    view.for_each_child([&](Node c) {
      EXPECT_EQ(original.child(k++), c);
      view.reset(c);
      EXPECT_TRUE(view.is_leaf());
      EXPECT_EQ(view.parent(), tree.root());
      EXPECT_EQ(view.number(), tree.number(c));
    });
    EXPECT_EQ(k, n - 1);
    EXPECT_TRUE(original.is_root());
    EXPECT_EQ(original.degree(), n - 1);
  }
}

TEST(OrdinalForest256, ThresholdsDeepChainAndWideRoot) {
  for (size_t n : {127, 128, 129, 255, 256, 257, 65537, 65538}) {
    for (bool chain : {false, true}) {
      SCOPED_TRACE(::testing::Message() << n << " chain=" << chain);
      std::vector<size_t> p(n);
      if (chain) {
        for (size_t i = 1; i < n; ++i) {
          p[i] = i - 1;
        }
      }
      Source s(std::move(p));
      Tree t(s);
      Check(s, t, kSlot);
      EXPECT_LE(t.unit_count(), 1 + n / 128);
    }
  }
}

TEST(OrdinalForest256, RandomTopologyAndReaderBudgets) {
  std::mt19937_64 rng(1729);
  std::vector<size_t> p(8192);
  for (size_t i = 1; i < p.size(); ++i) {
    p[i] = rng() % i;
  }
  Source s(std::move(p));
  Tree t(s);
  const size_t retained = t.memory_usage_bytes();
  for (size_t b : {size_t{0}, kSlot - 1, kSlot, kSlot + 1, 2 * kSlot - 1,
                   2 * kSlot, 8 * kSlot}) {
    EXPECT_EQ(Context(t, b).allocated_bytes(), b / kSlot * kSlot);
    Check(s, t, b);
  }
  EXPECT_EQ(t.memory_usage_bytes(), retained);
}

TEST(OrdinalForest256, CopiesMovesRebindingAndThrowingCallbacks) {
  Tree t = []() {
    Source original(std::vector<size_t>(1025));
    Tree first(original);
    return Tree(first);
  }();
  Source s({0, 0, 1});
  Tree other(s);
  Context r(t, kSlot);
  auto root = t.root();
  auto last = r.child(root, 1023);
  EXPECT_EQ(r.parent(last), root);
  auto copy = r;
  copy.bind(other);
  EXPECT_EQ(copy.degree(other.root()), 1u);
  EXPECT_EQ(r.degree(root), 1024u);
  auto moved = std::move(copy);
  EXPECT_THROW(copy.degree({}), std::logic_error);
  copy.bind(t);
  EXPECT_EQ(copy.allocated_bytes(), 0u);
  EXPECT_EQ(copy.degree(root), 1024u);
  copy = moved;
  EXPECT_EQ(copy.degree({}), 1u);
  r = std::move(copy);
  EXPECT_EQ(r.degree({}), 1u);
  r.bind(t);
  size_t j = 0;
  EXPECT_THROW(r.for_each_child(root,
                                [&](Node c) {
                                  EXPECT_EQ(r.parent(c), root);
                                  r.for_each_child(c, [](Node) {
                                    ADD_FAILURE() << "star leaf";
                                  });
                                  if (++j == 257) {
                                    throw std::runtime_error("callback");
                                  }
                                }),
               std::runtime_error);
  EXPECT_EQ(j, 257u);
  EXPECT_EQ(r.child(root, 1023), last);
  EXPECT_THROW(r.child(root, SIZE_MAX), std::out_of_range);
  EXPECT_THROW(r.first_child(last), std::out_of_range);
  EXPECT_THROW(r.next_sibling(last), std::out_of_range);
  EXPECT_THROW(r.parent({UINT32_MAX, 0}), std::out_of_range);
  EXPECT_THROW(t.number({0, 256}), std::out_of_range);
  EXPECT_THROW(t.find_by_number(t.size()), std::out_of_range);
}

TEST(OrdinalForest256, ExactInspectedEventsAcrossUnits) {
  // The root spans several units; below it mix local branches and portals.
  std::mt19937_64 rng(42);
  std::vector<size_t> parents(8192);
  for (size_t i = 1025; i < parents.size(); ++i) {
    parents[i] = rng() % i;
  }
  Source s(std::move(parents));
  Tree tree(s);
  Tree other(Source({0}));
  Context reference(tree, kSlot);
  std::vector<Node> handles(s.n);
  for (size_t v = 0; v < s.n; ++v) {
    for (size_t i = 0; i < s.children[v].size(); ++i) {
      handles[s.children[v][i]] = reference.child(handles[v], i);
    }
  }
  using Event = std::pair<char, Node>;
  std::vector<Event> expected, actual;
  auto dfs = [&](auto&& self, size_t v) -> void {
    expected.emplace_back('E', handles[v]);
    for (auto c : s.children[v]) {
      expected.emplace_back('C', handles[c]);
    }
    for (auto c : s.children[v]) {
      self(self, c);
    }
    if (!s.children[v].empty()) {
      expected.emplace_back('R', handles[v]);
    }
  };
  dfs(dfs, 0);  // Shallow source oracle, independent of encoded parents.
  for (size_t budget : {size_t{0}, kSlot, 4 * kSlot}) {
    actual.clear();
    Context r(tree, budget);
    auto inspect = [&](Node h) {
      actual.emplace_back('C', h);
      if (actual.size() == 2) {
        r.inspect_path([](Node) {}, [](Node) {},
                       [](size_t d) { return d - 1; });
      }
      if (actual.size() == 5) {
        // Transfer and destroy the actual slots during the outer callback.
        // Neither copying a reader nor a throwing inner walk restores them.
        auto moved = std::move(r);
        auto independent = moved;
        r = independent;
        std::vector<Event> prefix;
        EXPECT_THROW(
            r.traverse_inspected(
                [&](Node n, bool e) { prefix.emplace_back(e ? 'E' : 'R', n); },
                [&](Node n) {
                  prefix.emplace_back('C', n);
                  r = Context(other, budget);
                  throw std::runtime_error("inner child");
                }),
            std::runtime_error);
        EXPECT_EQ(prefix, (std::vector<Event>{expected[0], expected[1]}));
        prefix.clear();
        r.traverse_inspected(
            [&](Node n, bool e) { prefix.emplace_back(e ? 'E' : 'R', n); },
            [&](Node n) { prefix.emplace_back('C', n); });
        EXPECT_EQ(prefix, (std::vector<Event>{{'E', other.root()}}));
      }
    };
    r.traverse_inspected(
        [&](Node h, bool entry) { actual.emplace_back(entry ? 'E' : 'R', h); },
        inspect);
    EXPECT_EQ(actual, expected);
    EXPECT_TRUE(
        r.is_leaf(other.root()));  // Only the replacement root is a leaf.
  }

  for (auto seed : {0u, 255u, 256u, 1023u, 1729u}) {
    expected.clear();
    actual.clear();
    size_t v = 0;
    std::vector<size_t> degrees, seen_degrees;
    for (;;) {
      expected.emplace_back('E', handles[v]);
      const auto& ch = s.children[v];
      if (ch.empty()) {
        break;
      }
      degrees.push_back(ch.size());
      for (auto c : ch) {
        expected.emplace_back('C', handles[c]);
      }
      v = ch[seed % ch.size()];
    }
    Context r(tree, kSlot);
    r.inspect_path([&](Node h) { actual.emplace_back('E', h); },
                   [&](Node h) {
                     actual.emplace_back('C', h);
                     r.bind(other);
                   },
                   [&](size_t d) {
                     seen_degrees.push_back(d);
                     return seed % d;
                   });
    EXPECT_EQ(degrees, seen_degrees);
    EXPECT_EQ(expected, actual);
    EXPECT_TRUE(r.is_leaf({}));  // Callback binding survives the outer walk.
  }
}

TEST(OrdinalForest256, InspectedThrowStopsAtEventAndReaderRemainsValid) {
  Tree tree(Source({0, 0, 0, 1}));
  Context r(tree, kSlot);
  const auto left = r.child({}, 0), right = r.child({}, 1),
             grand = r.child(left, 0);
  using Event = std::pair<char, Node>;
  const std::vector<Event> dfs{{'E', {}},   {'C', left},  {'C', right},
                               {'E', left}, {'C', grand}, {'E', grand},
                               {'R', left}, {'E', right}, {'R', {}}};
  const std::vector<Event> path{{'E', {}},    {'Q', {}},   {'C', left},
                                {'C', right}, {'E', left}, {'Q', {}},
                                {'C', grand}, {'E', grand}};
  for (bool is_path : {false, true}) {
    const auto& expected = is_path ? path : dfs;
    for (size_t stop = 0; stop < expected.size(); ++stop) {
      size_t pos = 0;
      auto emit = [&](char kind, Node h) {
        EXPECT_EQ((Event{kind, h}), expected.at(pos));
        if (pos++ == stop) {
          throw std::runtime_error("event");
        }
      };
      auto inspect = [&](Node h) { emit('C', h); };
      auto walk = [&] {
        if (is_path) {
          r.inspect_path([&](Node h) { emit('E', h); }, inspect,
                         [&](size_t) {
                           emit('Q', {});
                           return 0;
                         });
        } else {
          r.traverse_inspected([&](Node h, bool e) { emit(e ? 'E' : 'R', h); },
                               inspect);
        }
      };
      EXPECT_THROW(walk(), std::runtime_error);
      EXPECT_EQ(pos, stop + 1);
      EXPECT_EQ(r.parent(grand), left);
    }
  }
  auto noop = [](auto...) {};
  EXPECT_THROW(
      r.inspect_path(
          noop, [](Node) { ADD_FAILURE(); }, [](size_t d) { return d; }),
      std::out_of_range);
  Tree empty(Source({}));
  r.bind(empty);
  EXPECT_THROW(r.traverse_inspected(noop, noop), std::out_of_range);
  EXPECT_THROW(r.inspect_path(noop, noop, [](size_t) { return 0; }),
               std::out_of_range);
}

TEST(OrdinalForest256, DeepInspectedWalkIsIterative) {
  std::vector<size_t> parents(131073);
  for (size_t i = 1; i < parents.size(); ++i) {
    parents[i] = i - 1;
  }
  Tree tree{Source(std::move(parents))};
  Context r(tree, kSlot);
  size_t entries = 0, returns = 0, children = 0;
  r.traverse_inspected([&](Node, bool e) { ++(e ? entries : returns); },
                       [&](Node) { ++children; });
  EXPECT_EQ(entries, tree.size());
  EXPECT_EQ(returns, tree.size() - 1);
  EXPECT_EQ(children, tree.size() - 1);
}

TEST(OrdinalForest256, UpperPrefixBudgetOwnerAndExactFrontier) {
  std::vector<size_t> parents(1301);
  for (size_t i = 301; i < parents.size(); ++i) {
    parents[i] = i - 300;
  }
  Source s(std::move(parents));
  Tree tree(s), other(s);
  Context r(tree);
  EXPECT_EQ(r.allocated_bytes(), 16 * kSlot);
  EXPECT_EQ(Context(tree, kSlot - 1).allocated_bytes(), 0u);
  EXPECT_EQ(Context(tree, 2 * kSlot - 1).allocated_bytes(), kSlot);
  std::vector<Node> handles(s.n);
  for (size_t v = 0; v < s.n; ++v) {
    for (size_t j = 0; j < s.children[v].size(); ++j) {
      handles[s.children[v][j]] = r.child(handles[v], j);
    }
  }
  for (size_t budget : {0u, 40u, 5000u, 65536u}) {
    Index index(tree, budget);
    EXPECT_LE(index.budget_usage_bytes(), budget);
    if (!budget) {
      EXPECT_EQ(index.size(), 0u);
    }
    auto noop = [](Node) {};
    auto first = [](size_t) { return 0; };
    EXPECT_THROW(Context(other, 0, &index), std::invalid_argument);
    EXPECT_THROW(r.bind(other, &index), std::invalid_argument);
    for (size_t k : {0u, 255u, 256u, 299u}) {
      using Event = std::pair<char, size_t>;
      std::vector<Event> want, got;
      for (size_t v = 0;;) {
        want.emplace_back('E', tree.number(handles[v]));
        const auto& ch = s.children[v];
        if (ch.empty()) {
          break;
        }
        want.emplace_back('Q', ch.size());
        for (auto c : ch) {
          want.emplace_back('C', tree.number(handles[c]));
        }
        v = ch[k % ch.size()];
      }
      r.bind(tree, &index);
      auto hits =
          r.inspect_path([&](Node n) { got.emplace_back('E', tree.number(n)); },
                         [&](Node n) {
                           got.emplace_back('C', tree.number(n));
                           if (got.size() == 3) {
                             r.inspect_path(noop, noop, first);
                             r.bind(other);
                           }
                         },
                         [&](size_t d) {
                           got.emplace_back('Q', d);
                           return k % d;
                         });
      EXPECT_EQ(got, want);
      EXPECT_EQ(hits != 0, budget >= 5000);
    }
    r.bind(tree, &index);
    EXPECT_THROW(
        r.inspect_path(
            noop, [](Node) { ADD_FAILURE(); }, [](size_t d) { return d; }),
        std::out_of_range);
    EXPECT_THROW(
        r.inspect_path(
            noop, [](Node) { throw std::runtime_error("stop"); }, first),
        std::runtime_error);
  }
}
}  // namespace
