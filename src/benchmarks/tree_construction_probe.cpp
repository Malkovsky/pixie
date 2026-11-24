#include <pixie/rmm/btree.h>
#include <pixie/tree/forest_cover.h>
#include <pixie/tree/traversal_context.h>
#include <pixie/utils.h>
#include <references/bfs_csr_tree.h>
#include <references/pointer_tree.h>
#include <sys/resource.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "tree_cover_source.h"
#include "tree_dfs/optimized.h"
#include "tree_input.h"

using Clock = std::chrono::steady_clock;
using Rmm = pixie::RmMBTree<>;
using Bp = pixie::tree::BPTree<Rmm>;
using Dfuds = pixie::tree::DFUDSTree<Rmm>;
using Louds = pixie::tree::LoudsTree<pixie::RankSelectSupportCompressed<>>;
using Pointer = pixie::references::PointerTree;

struct Result {
  double generate = 0, convert = 0, construct = 0, cache = 0;
  size_t source = 0, encoding = 0, structure = 0, persistent = 0, scratch = 0;
};
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

template <class Tree>
Result succinct(size_t n, size_t budget, bool copy_input) {
  Result r;
  auto start = Clock::now();
  auto words = pixie::benchmarks::GenerateRandomBp(n);
  r.generate = seconds(start);
  start = Clock::now();
  if constexpr (std::is_same_v<Tree, Dfuds>) {
    words = pixie::benchmarks::ConvertBpEncoding<
        pixie::tree::detail::TreeEncodingKind::kDFUDS>(words, n);
  } else if constexpr (std::is_same_v<Tree, Louds>) {
    words = pixie::benchmarks::ConvertBpEncoding<
        pixie::tree::detail::TreeEncodingKind::kLOUDS>(words, n);
  }
  pixie::AlignedStorage encoding;
  if (!copy_input) {
    encoding = pixie::AlignedStorage(std::span<const uint64_t>(words));
    std::vector<uint64_t>().swap(words);
  }
  r.convert = seconds(start);
  start = Clock::now();
  Tree tree = copy_input ? Tree(words, n) : Tree(std::move(encoding), n);
  r.construct = seconds(start);
  if (budget == 0) {
    // Owning trees include the encoding once; their intrinsic indexes borrow
    // it.
    r.encoding = tree.encoding_owned_bytes();
    r.structure = tree.memory_usage_bytes();
    if (tree.size() != n || !r.encoding || r.structure <= r.encoding ||
        (n > 1 && tree.degree(tree.root()) == 0)) {
      throw std::runtime_error("uncached construction validation failed");
    }
    return r;
  }
  start = Clock::now();
  pixie::benchmarks::TreeTraversalWorkloads walk(tree, budget);
  pixie::tree::RootNavigationIndex index(tree, budget - walk.reserved_bytes());
  r.cache = seconds(start);
  r.encoding = tree.encoding_owned_bytes();
  r.structure = tree.memory_usage_bytes();
  r.scratch = walk.scratch_bytes();
  r.persistent = walk.fixture_metadata_bytes() + index.memory_usage_bytes();
  if (tree.size() != n ||
      walk.reserved_bytes() + index.budget_usage_bytes() > budget ||
      !walk.optimized_inspection_enabled()) {
    throw std::runtime_error("construction validation failed");
  }
  return r;
}

Result forest_cover(size_t n, size_t budget) {
  using namespace pixie::benchmarks;
  using Tree = pixie::tree::ForestCoverTree;
  using Context = pixie::tree::TreeTraversalContext<Tree>;
  if (budget != 0 && budget != 32768 && budget != 524288) {
    throw std::invalid_argument("unsupported ForestCover budget");
  }
  Result r;
  auto start = Clock::now();
  CoverSource source(n, CoverShape::kRandomRecursive, 42);
  r.generate = seconds(start);
  start = Clock::now();
  Tree tree(source);
  r.construct = seconds(start);
  if (budget == 0) {
    r.structure = tree.memory_usage_bytes();
    if (tree.size() != n ||
        tree.degree(tree.root()) != source.degree(source.root())) {
      throw std::runtime_error(
          "uncached ForestCover construction validation failed");
    }
    return r;
  }
  // Match the charged profiles used for ForestCover inspected traversals.
  const size_t cache = (budget == 32768 ? 4 : 16) * Context::slot_bytes();
  r.scratch = sizeof(pixie::tree::detail::forest256::LevelCsr256);
  start = Clock::now();
  pixie::tree::RootNavigationIndex<Tree> index(
      tree, budget - sizeof(Context) - cache - r.scratch);
  Context context(tree, cache, &index);
  r.cache = seconds(start);
  r.structure = tree.memory_usage_bytes();
  r.persistent = context.memory_usage_bytes() + index.budget_usage_bytes();
  if (tree.size() != n || !index.size() || context.allocated_bytes() != cache ||
      r.persistent + r.scratch > budget) {
    throw std::runtime_error("ForestCover construction validation failed");
  }
  return r;
}

template <class Tree>
Result reference(size_t n, bool bfs) {
  Result r;
  std::mt19937_64 rng(42);
  auto start = Clock::now();
  auto adjacency = generate_random_tree(n, rng);
  r.generate = seconds(start);
  if constexpr (!std::is_same_v<Tree, pixie::references::BfsCsrTree>) {
    start = Clock::now();
    adjacency = bfs ? bfs_order(n, adjacency) : dfs_order(n, adjacency);
    r.convert = seconds(start);
  }
  start = Clock::now();
  Tree tree(adjacency);
  r.construct = seconds(start);
  r.structure = tree.memory_usage_bytes();
  if constexpr (std::is_same_v<Tree, pixie::references::BfsCsrTree>) {
    if (tree.size() != n ||
        (n && tree.degree(tree.root()) + 1 != adjacency[0].size())) {
      throw std::runtime_error("BFS CSR construction validation failed");
    }
  }
  return r;
}

int main(int argc, char** argv) {
  if (argc != 4 && argc != 5) {
    return 2;
  }
  const std::string name = argv[1];
  const size_t n = std::stoull(argv[2]), budget = std::stoull(argv[3]);
  const std::string mode = argc == 5 ? argv[4] : "consume";
  if (mode != "copy" && mode != "consume") {
    return 2;
  }
  const bool copy_input = mode == "copy";
  Result r;
  if (name == "forest_cover") {
    r = forest_cover(n, budget);
  } else if (name == "bp") {
    r = succinct<Bp>(n, budget, copy_input);
  } else if (name == "dfuds") {
    r = succinct<Dfuds>(n, budget, copy_input);
  } else if (name == "louds") {
    r = succinct<Louds>(n, budget, copy_input);
  } else if (name == "vectors_bfs") {
    r = reference<AdjListTree>(n, true);
  } else if (name == "vectors_preorder") {
    r = reference<AdjListTree>(n, false);
  } else if (name == "pointers_bfs") {
    r = reference<Pointer>(n, true);
  } else if (name == "pointers_preorder") {
    r = reference<Pointer>(n, false);
  } else if (name == "bfs_csr" && budget == 0) {
    r = reference<pixie::references::BfsCsrTree>(n, true);
  } else {
    return 2;
  }
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  std::cout << std::setprecision(12) << "{\"tree\":\"" << name
            << "\",\"construction_mode\":\""
            << (name == "forest_cover" ? "ordered_source"
                : name == "bfs_csr"    ? "adjacency"
                                       : mode)
            << "\",\"nodes\":" << n << ",\"budget\":" << budget
            << ",\"generate_s\":" << r.generate
            << ",\"convert_s\":" << r.convert
            << ",\"construct_s\":" << r.construct << ",\"cache_s\":" << r.cache
            << ",\"total_s\":" << r.generate + r.convert + r.construct + r.cache
            << ",\"source_bytes\":" << r.source
            << ",\"encoding_bytes\":" << r.encoding
            << ",\"structure_bytes\":" << r.structure
            << ",\"persistent_cache_bytes\":" << r.persistent
            << ",\"scratch_bytes\":" << r.scratch
            << ",\"retained_bytes\":" << r.source + r.structure + r.persistent
            << ",\"working_bytes\":"
            << r.source + r.structure + r.persistent + r.scratch
            << ",\"peak_process_rss_kib\":" << usage.ru_maxrss << "}\n";
}
