#include <pixie/rmm/btree.h>
#include <pixie/tree/traversal_context.h>
#include <pixie/utils.h>
#include <references/pointer_tree.h>
#include <sys/resource.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "tree_input.h"

using Clock = std::chrono::steady_clock;
using Rmm = pixie::experimental::RmMBTree<>;
using Bp = pixie::BPTree<Rmm>;
using Dfuds = pixie::DFUDSTree<Rmm>;
using Louds = pixie::LoudsTree<pixie::RankSelectSupportCompressed<>>;
using Pointer = pixie::references::PointerTree;

struct Result {
  double generate = 0, convert = 0, construct = 0, cache = 0;
  size_t source = 0, structure = 0, persistent = 0, scratch = 0;
};
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

template <class Tree>
Result succinct(size_t n, size_t budget) {
  Result r;
  auto start = Clock::now();
  auto words = pixie::benchmarks::GenerateRandomBp(n);
  r.generate = seconds(start);
  start = Clock::now();
  if constexpr (std::is_same_v<Tree, Dfuds>) {
    words = pixie::benchmarks::DfudsFromBp(words, n);
  } else if constexpr (std::is_same_v<Tree, Louds>) {
    words = pixie::benchmarks::LoudsFromBp(words, n);
  }
  r.convert = seconds(start);
  start = Clock::now();
  Tree tree(words, n);
  r.construct = seconds(start);
  start = Clock::now();
  pixie::TreeTraversalContext cache(tree, words, budget);
  r.cache = seconds(start);
  r.source = words.capacity() * sizeof(uint64_t);
  r.structure = tree.memory_usage_bytes();
  r.scratch = cache.scratch_bytes();
  r.persistent = cache.memory_usage_bytes() - r.scratch;
  if (tree.size() != n || cache.budget_usage_bytes() > budget ||
      !cache.optimized_inspection_enabled()) {
    throw std::runtime_error("construction validation failed");
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
  start = Clock::now();
  adjacency = bfs ? bfs_order(n, adjacency) : dfs_order(n, adjacency);
  r.convert = seconds(start);
  start = Clock::now();
  Tree tree(adjacency);
  r.construct = seconds(start);
  r.structure = tree.memory_usage_bytes();
  return r;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    return 2;
  }
  const std::string name = argv[1];
  const size_t n = std::stoull(argv[2]), budget = std::stoull(argv[3]);
  Result r;
  if (name == "bp") {
    r = succinct<Bp>(n, budget);
  } else if (name == "dfuds") {
    r = succinct<Dfuds>(n, budget);
  } else if (name == "louds") {
    r = succinct<Louds>(n, budget);
  } else if (name == "vectors_bfs") {
    r = reference<AdjListTree>(n, true);
  } else if (name == "vectors_preorder") {
    r = reference<AdjListTree>(n, false);
  } else if (name == "pointers_bfs") {
    r = reference<Pointer>(n, true);
  } else if (name == "pointers_preorder") {
    r = reference<Pointer>(n, false);
  } else {
    return 2;
  }
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  std::cout << std::setprecision(12) << "{\"tree\":\"" << name
            << "\",\"nodes\":" << n << ",\"budget\":" << budget
            << ",\"generate_s\":" << r.generate
            << ",\"convert_s\":" << r.convert
            << ",\"construct_s\":" << r.construct << ",\"cache_s\":" << r.cache
            << ",\"total_s\":" << r.generate + r.convert + r.construct + r.cache
            << ",\"source_bytes\":" << r.source
            << ",\"structure_bytes\":" << r.structure
            << ",\"persistent_cache_bytes\":" << r.persistent
            << ",\"scratch_bytes\":" << r.scratch
            << ",\"retained_bytes\":" << r.source + r.structure + r.persistent
            << ",\"working_bytes\":"
            << r.source + r.structure + r.persistent + r.scratch
            << ",\"peak_process_rss_kib\":" << usage.ru_maxrss << "}\n";
}
