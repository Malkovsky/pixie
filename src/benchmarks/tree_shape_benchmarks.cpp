#include <pixie/rmm/btree.h>
#include <pixie/tree/implementations.h>
#include <pixie/utils.h>
#include <references/bfs_csr_tree.h>
#include <references/pointer_tree.h>

#include <any>
#include <array>
#include <memory>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "tree_benchmark_helpers.h"
#include "tree_shape_workloads.h"

namespace pixie::benchmarks {
namespace {
using Forest = tree::ForestCoverTree;
using Context = tree::TreeTraversalContext<Forest>;

struct ForestWalk {
  const Forest& tree;
  Context& context;
  size_t size() const { return tree.size(); }
  void traverse_inspected(auto& visit, auto& child) const {
    context.traverse_inspected(visit, child);
  }
  void inspect_path(auto& visit, auto& child, auto& choose) const {
    context.inspect_path(visit, child, choose);
  }
};

template <class Tree>
struct ReferenceWalk {
  using Node = decltype(std::declval<const Tree&>().root());
  const Tree& tree;
  struct Cursor {
    const Tree& tree;
    Node current;
    Node node() const { return current; }
    void reset(Node n) { current = n; }
    bool is_leaf() const { return tree.is_leaf(current); }
    bool is_root() const { return tree.is_root(current); }
    bool is_last_child() const { return tree.is_last_child(current); }
    Node parent() const { return tree.parent(current); }
    Node first_child() const { return tree.first_child(current); }
    Node next_sibling() const { return tree.next_sibling(current); }
    size_t degree() const { return tree.degree(current); }
    void for_each_child(auto& visit) const {
      tree.for_each_child(current, visit);
    }
  };
  Node root() const { return tree.root(); }
  Cursor navigate(Node n) const { return {tree, n}; }
  size_t size() const { return tree.size(); }
};

template <class Tree>
Tree BuildMatchedTree(const CoverSource& source);

// Reuse only the last immutable reference across calibration and repetitions.
// In particular, never retain two large reference layouts simultaneously.
std::any last_reference;

template <class Tree>
void MatchedReference(benchmark::State& state, CoverShape shape, bool dfs) {
  struct Fixture {
    CoverShape shape;
    uint64_t seed;
    Tree tree;
  };
  const size_t n = state.range(0);
  const uint64_t seed = state.range(1);
  auto* stored = std::any_cast<std::shared_ptr<Fixture>>(&last_reference);
  if (!stored || (*stored)->tree.size() != n || (*stored)->shape != shape ||
      (*stored)->seed != seed) {
    last_reference.reset();
    stored = &last_reference.emplace<std::shared_ptr<Fixture>>(
        std::make_shared<Fixture>(Fixture{
            shape, seed, BuildMatchedTree<Tree>(CoverSource(n, shape, seed))}));
  }
  const auto& tree = (*stored)->tree;
  state.counters["tree_bytes"] = tree.memory_usage_bytes();
  state.counters["nodes"] = n;
  state.counters["budget_limit_bytes"] = 0;
  state.counters["budget_usage_bytes"] = 0;
  ReferenceWalk<Tree> walk{tree};
  if (dfs) {
    return CoverDfs<true, false, false>(state, walk);
  }
  CoverInspectedRootPath<false, true>(state, walk, 4096);
}

template <class Tree>
void Matched(benchmark::State& state, CoverShape shape, bool dfs) {
  const size_t budget = state.range(2);
  // All families replay the same batch regardless of calibration. A single
  // star or chain path already inspects the entire tree.
  const size_t paths =
      shape == CoverShape::kStar || shape == CoverShape::kChain ? 1 : 4096;
  const Tree tree(CoverSource(state.range(0), shape, state.range(1)));
  state.counters["tree_bytes"] = tree.memory_usage_bytes();
  state.counters["nodes"] = tree.size();
  state.counters["budget_limit_bytes"] = budget;
  if constexpr (std::is_same_v<Tree, Forest>) {
    // Baseline is one cache slot without a root index. Larger profiles
    // reserve the context and a decoded unit before spending on the index.
    const size_t slots = budget == 0 ? 1 : budget == 32768 ? 4 : 16;
    const size_t cache = slots * Context::slot_bytes();
    const size_t scratch = sizeof(tree::detail::forest256::LevelCsr256);
    const size_t fixed = sizeof(Context) + cache + scratch;
    tree::RootNavigationIndex<Forest> index(tree, budget ? budget - fixed : 0);
    Context context(tree, cache, budget ? &index : nullptr);
    const auto charge =
        context.memory_usage_bytes() + scratch + index.budget_usage_bytes();
    if (budget && charge > budget) {
      state.SkipWithError("Forest traversal exceeds combined budget");
      return;
    }
    state.counters["budget_usage_bytes"] = charge;
    state.counters["cache_slots"] = slots;
    state.counters["cache_allocated_bytes"] = context.allocated_bytes();
    state.counters["scratch_bytes"] = scratch;
    state.counters["traversal_index_bytes"] = index.budget_usage_bytes();
    const ForestWalk walk{tree, context};
    if (dfs) {
      return CoverDfs<true, false, false>(state, walk);
    }
    return CoverInspectedRootPath<false, true>(state, walk, paths);
  } else {
    if (!budget) {
      state.counters["budget_usage_bytes"] = 0;
      if (dfs) {
        return CoverDfs<true, false, false>(state, tree);
      }
      return CoverInspectedRootPath<false, true>(state, tree, paths);
    }
    if (dfs) {
      return CoverDfs<true, false, true>(state, tree, budget);
    }
    auto index = MakeWalkIndex(tree, budget);
    tree::TreeTraversalContext context(tree, budget, &index);
    state.counters["budget_usage_bytes"] =
        context.budget_usage_bytes() + index.budget_usage_bytes();
    state.counters["scratch_bytes"] = context.scratch_bytes();
    state.counters["traversal_index_bytes"] = context.memory_usage_bytes() -
                                              context.scratch_bytes() +
                                              index.memory_usage_bytes();
    return CoverInspectedRootPath<true, true>(state, context, paths);
  }
}

constexpr std::array kShapes{
    std::pair{CoverShape::kRandomRecursive, "random_recursive"},
    std::pair{CoverShape::kBinaryCatalan, "binary_catalan"},
    std::pair{CoverShape::kCatalan, "ordered_plane_catalan"}};
template <class Tree>
void RegisterMatched(const char* family, bool dfs) {
  for (auto [shape, label] : kShapes) {
    const auto name = std::string("tree_shape/") + family + "/" +
                      (dfs ? "inspected_dfs/" : "inspected_root_path/") + label;
    constexpr bool reference = !requires(const Tree& t) {
      t.navigate(t.root());
    } && !std::is_same_v<Tree, Forest>;
    auto* row = benchmark::RegisterBenchmark(
        name.c_str(), [=](benchmark::State& state) {
          if constexpr (reference) {
            MatchedReference<Tree>(state, shape, dfs);
          } else {
            Matched<Tree>(state, shape, dfs);
          }
        });
    row->ArgNames({"N", "seed", "budget"});
    for (int budget : {0, 32768, 524288}) {
      if (reference && budget) {
        continue;
      }
      for (auto n = kMinTreeNodes; n <= kMaxTreeNodes;
           n *= kTreeSizeMultiplier) {
        row->Args({n, 42, budget});
      }
    }
  }
}

template <class Tree>
class PrefixParents {
  const Tree& tree_;
  std::vector<tree::TreeNode> parents_;

 public:
  PrefixParents(const Tree& tree, size_t budget) : tree_(tree) {
    const auto count = std::min(
        tree.size(), (budget - sizeof(*this)) / sizeof(tree::TreeNode));
    parents_.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      parents_.push_back(i ? tree.parent(tree.find_by_number(i)) : tree.root());
    }
  }
  size_t entries() const { return parents_.size(); }
  size_t bytes() const {
    return sizeof(*this) + parents_.capacity() * sizeof(tree::TreeNode);
  }
  void tree_path(tree::TreeNode u,
                 tree::TreeNode v,
                 std::vector<tree::TreeNode>& out,
                 std::vector<tree::TreeNode>& scratch) const {
    // Native preorder (BP/DFUDS) or BFS (LOUDS) decreases toward the root.
    tree::detail::collect_tree_path(
        u, v, out, scratch, [](tree::TreeNode h) { return h.number; },
        [this](tree::TreeNode h) {
          return h.number < parents_.size() ? parents_[h.number]
                                            : tree_.parent(h);
        });
  }
};

template <class Tree>
Tree BuildMatchedTree(const CoverSource& source) {
  if constexpr (std::is_constructible_v<Tree, const CoverSource&>) {
    return Tree(source);
  } else {
    std::vector<std::vector<size_t>> adjacency(source.size());
    for (uint32_t p = 0; p < source.size(); ++p) {
      adjacency[p].resize(1 + source.degree(p));
    }
    for (uint32_t p = 0; p < source.size(); ++p) {
      auto c = source.first_child(p);
      for (size_t i = 0; i < source.degree(p); ++i) {
        adjacency[p][1 + i] = c;
        adjacency[c][0] = p;
        c = source.next_sibling(c);
      }
    }
    return Tree(adjacency);
  }
}

template <class Tree>
struct PathFixture {
  using Node = decltype(std::declval<const Tree&>().root());
  CoverShape shape;
  uint64_t seed;
  bool siblings;
  const Tree tree;
  std::vector<std::pair<Node, Node>> pairs;
  std::vector<std::vector<Node>> expected;
  size_t max_depth = 0;

  PathFixture(const CoverSource& source,
              CoverShape shape,
              uint64_t seed,
              bool siblings)
      : shape(shape),
        seed(seed),
        siblings(siblings),
        tree(BuildMatchedTree<Tree>(source)),
        expected(256) {
    const auto n = source.size();
    std::vector<uint32_t> depth(n), parent(n), sibling_starts;
    for (uint32_t p = 0; p < n; ++p) {
      auto c = source.first_child(p);
      for (size_t j = 0; j < source.degree(p); ++j) {
        parent[c] = p;
        depth[c] = depth[p] + 1;
        max_depth = std::max<size_t>(max_depth, depth[c]);
        const auto next = source.next_sibling(c);
        if (siblings && next < n) {
          sibling_starts.push_back(c);
        }
        c = next;
      }
    }
    if (siblings && sibling_starts.empty()) {
      return;
    }
    std::mt19937_64 random(seed ^ 0x9e3779b97f4a7c15ULL);
    std::vector<std::pair<uint32_t, uint32_t>> ids;
    std::unordered_map<uint32_t, Node> handles;
    for (size_t i = 0; i < expected.size(); ++i) {
      uint32_t a, b;
      if (siblings) {
        a = sibling_starts[random() % sibling_starts.size()];
        b = source.next_sibling(a);
      } else {
        a = random() % n;
        b = random() % n;
      }
      ids.emplace_back(a, b);
      for (auto id : {a, b}) {
        while (handles.emplace(id, tree.root()).second && id) {
          id = parent[id];
        }
      }
    }
    // Translate only ancestors of selected endpoints, by matched child order.
    // Parent-before-child source IDs let us process these without recursion.
    std::vector<uint32_t> ancestors;
    for (const auto& [id, node] : handles) {
      ancestors.push_back(id);
    }
    std::sort(ancestors.begin(), ancestors.end());
    for (auto p : ancestors) {
      auto c = source.first_child(p);
      tree.for_each_child(handles.at(p), [&](Node node) {
        if (auto it = handles.find(c); it != handles.end()) {
          it->second = node;
        }
        c = source.next_sibling(c);
      });
    }
    for (size_t i = 0; i < ids.size(); ++i) {
      auto [a, b] = ids[i];
      pairs.emplace_back(handles.at(a), handles.at(b));
      std::vector<Node> suffix;
      while (a != b) {
        if (depth[a] >= depth[b]) {
          expected[i].push_back(handles.at(a));
          a = parent[a];
        } else {
          suffix.push_back(handles.at(b));
          b = parent[b];
        }
      }
      expected[i].push_back(handles.at(a));
      expected[i].insert(expected[i].end(), suffix.rbegin(), suffix.rend());
    }
  }
};

// Retain at most one immutable fixture across calibration/repetitions and
// budgets. Readers and their scratch remain fresh for every invocation.
std::any last_path_fixture;

template <class Tree>
void MatchedPaths(benchmark::State& state,
                  CoverShape shape,
                  bool siblings,
                  size_t budget) {
  using Fixture = PathFixture<Tree>;
  using Node = typename Fixture::Node;
  const size_t n = state.range(0);
  const uint64_t seed = state.range(1);
  auto* stored = std::any_cast<std::shared_ptr<Fixture>>(&last_path_fixture);
  if (!stored || (*stored)->tree.size() != n || (*stored)->shape != shape ||
      (*stored)->siblings != siblings || (*stored)->seed != seed) {
    last_path_fixture.reset();
    stored = &last_path_fixture.emplace<std::shared_ptr<Fixture>>(
        std::make_shared<Fixture>(CoverSource(n, shape, seed), shape, seed,
                                  siblings));
  }
  const auto& [fixture_shape, fixture_seed, fixture_siblings, tree, pairs,
               expected, max_depth] = **stored;
  if (pairs.empty()) {
    state.SkipWithError("Tree has no sibling pairs");
    return;
  }
  std::vector<Node> out, scratch;
  out.reserve(2 * max_depth + 2);
  scratch.reserve(2 * max_depth + 2);
  size_t nodes = 0, longest = 0;
  state.counters["budget_limit_bytes"] = budget;
  state.counters["tree_bytes"] = tree.memory_usage_bytes();
  state.counters["path_buffer_bytes"] =
      (out.capacity() + scratch.capacity()) * sizeof(Node);
  auto measure = [&](auto& reader, size_t charge) {
    state.counters["budget_usage_bytes"] = charge;
    if (budget && charge > budget) {
      state.SkipWithError("Parent navigation exceeds budget");
      return;
    }
    while (state.KeepRunningBatch(pairs.size())) {
      for (const auto& [u, v] : pairs) {
        reader.tree_path(u, v, out, scratch);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
      }
    }
    for (size_t i = 0; i < pairs.size(); ++i) {
      const auto& [u, v] = pairs[i];
      reader.tree_path(u, v, out, scratch);
      const auto equal = [](Node a, Node b) {
        if constexpr (std::is_same_v<Node, AdjListNode>) {
          return a.number == b.number;
        } else {
          return a == b;
        }
      };
      if (!std::equal(out.begin(), out.end(), expected[i].begin(),
                      expected[i].end(), equal)) {
        state.SkipWithError("Path differs from source parent/depth oracle");
        return;
      }
      nodes += out.size();
      longest = std::max(longest, out.size());
    }
  };
  if constexpr (std::is_same_v<Tree, Forest>) {
    using Reader = tree::TreeTraversalContext<Forest>;
    const auto slots =
        budget ? std::min(tree.unit_count(),
                          (budget - sizeof(Reader)) / Reader::slot_bytes())
               : 1;
    Reader reader(tree, slots * Reader::slot_bytes());
    measure(reader, reader.memory_usage_bytes());
    state.counters["decoder_slots"] = slots;
    state.counters["decoder_cache_bytes"] = reader.allocated_bytes();
  } else if constexpr (std::is_same_v<Node, tree::TreeNode>) {
    if (budget) {
      PrefixParents reader(tree, budget);
      measure(reader, reader.bytes());
      state.counters["parent_entries"] = reader.entries();
    } else {
      measure(tree, 0);
    }
  } else {
    measure(tree, 0);
  }
  state.counters["nodes_per_path"] = double(nodes) / pairs.size();
  state.counters["max_nodes_per_path"] = longest;
  state.counters["pairs_per_batch"] = pairs.size();
  state.SetItemsProcessed(state.iterations());
}

template <class Tree>
void RegisterPaths(const char* family) {
  for (auto [shape, label] : kShapes) {
    for (bool siblings : {false, true}) {
      const auto name = std::string("tree_path/") + family +
                        (siblings ? "/sibling/" : "/random/") + label;
      auto* row = benchmark::RegisterBenchmark(name.c_str(), [=](auto& state) {
        MatchedPaths<Tree>(state, shape, siblings, 0);
      });
      row->ArgNames({"N", "seed"});
      for (int power = 10; power <= 26; power += 2) {
        if (!siblings || power == 10 || power == 16 || power == 20) {
          row->Args({1 << power, 42});
        }
      }
      if constexpr (std::is_same_v<Tree, Forest> ||
                    std::is_same_v<typename PathFixture<Tree>::Node,
                                   tree::TreeNode>) {
        if (!siblings) {
          auto* indexed =
              benchmark::RegisterBenchmark(name.c_str(), [=](auto& state) {
                MatchedPaths<Tree>(state, shape, false, state.range(2));
              });
          indexed->ArgNames({"N", "seed", "budget"});
          for (int power = 10; power <= 26; power += 2) {
            for (int bytes : {32768, 524288}) {
              indexed->Args({1 << power, 42, bytes});
            }
          }
        }
      }
    }
  }
}

void RegisterMatchedFamilies(bool dfs) {
  RegisterMatched<AdjListTree>("adj_list", dfs);
  RegisterMatched<references::PointerTree>("pointer", dfs);
  RegisterMatched<references::BfsCsrTree>("bfs_csr", dfs);
  RegisterMatched<tree::BPTree<RmMBTree<>>>("bp", dfs);
  RegisterMatched<tree::DFUDSTree<RmMBTree<>>>("dfuds", dfs);
  RegisterMatched<tree::LoudsTree<>>("louds", dfs);
  RegisterMatched<Forest>("forest_cover", dfs);
}
}  // namespace

void RegisterShapeDfs() {
  RegisterMatchedFamilies(true);
}

void RegisterShapePaths() {
  RegisterMatchedFamilies(false);
  RegisterPaths<tree::BPTree<RmMBTree<>>>("bp");
  RegisterPaths<tree::DFUDSTree<RmMBTree<>>>("dfuds");
  RegisterPaths<tree::LoudsTree<>>("louds");
  RegisterPaths<Forest>("forest_cover");
  RegisterPaths<AdjListTree>("adj_list");
  RegisterPaths<references::PointerTree>("pointer");
  RegisterPaths<references::BfsCsrTree>("bfs_csr");
}
}  // namespace pixie::benchmarks
