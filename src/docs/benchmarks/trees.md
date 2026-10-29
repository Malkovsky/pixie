# Trees {#benchmarks_trees}

This page summarizes benchmarks for tree implementation: BP, DFUDS, LOUDS and reference implementations.

## Machine

AMD Ryzen 7 8845HS (8 cores, 16 threads, 16 MiB L3), 15 GiB RAM,
WSL2 (Linux 6.6.87.2-microsoft-standard-WSL2). GCC 13.3, C++20,
`-O3 -DNDEBUG -march=native`; measurements pinned to CPU 0.
Results captured October 28, 2026, before the ownership refactor.
The current API supplies the same combined cache/scratch budget to the external
`TreeTraversalContext(tree, words, byte_limit)`, rather than to the tree constructor.
The graphs are the saved measurements, not a new run of the refactored code.

## Implementations and inputs

Random recursive trees, seed 42. All measurements cover every even power from
\f$2^{10}\f$ through \f$2^{26}\f$, with **32 KiB and 512 KiB** cache budgets for succinct trees.

| Implementation | Representation and navigation | Space |
| --- | --- | --- |
| BP | Balanced parentheses, RmMBTree with compressed rank/select | \f$2n\f$ bits + navigation index + cache |
| DFUDS | Depth-first unary degree sequence, RmMBTree with compressed rank/select | \f$2n-1\f$ bits + navigation index + cache |
| LOUDS | Level-order unary degree sequence, compressed rank/select | \f$2n-1\f$ bits + navigation index + cache |
| Vectors | One adjacency vector per node, storing parent followed by children | \f$n\f$ instances of `std::vector<size_t>` with parent + children entries, approximately \f$40n\f$ bytes on this machine |
| Pointers | One contiguous array of first-child, next-sibling and parent pointers | \f$3n\f$ pointers, approximately \f$24n\f$ bytes on this machine |

The reference implementations are measured in both BFS and preorder numbering, succinct representation
use their native numbering and bounded cache.

## Measurement

- **DFS:** traverse the whole tree, consuming preorder entry events and a return
  event for each nonleaf after its children.
- **DFS + child inspection:** also enumerate every immediate child before
  descending.
- **Root-to-leaf path:** start from root, choose a random child at each level until reaching a leaf.
- **Path + child inspection:** also enumerate every immediate child along the path.
- **Construction:** build from prepared input (including indexes and cache),
  and full setup including input generation and conversion/reordering.
- **Memory usage:** bytes/node including the representation, indexes, cache and
  reserved traversal scratch; excludes allocator metadata and construction scratch.

Traversal timings are median CPU times over five repetitions; construction
timings are median wall times over three fresh-process repetitions.
Both cache sizes share each graph: dotted, hollow-marker curves denote 32 KiB;
solid colored curves denote 512 KiB. References are shown once. All succinct
curves use the current runtime-adaptive configuration.

Prepared input is a bitvector for succinct trees and an adjacency list for
references. Full setup generates BP and converts it to DFUDS/LOUDS for succinct
trees; references use the ordinary adjacency generator and BFS/preorder
reordering. These timings compare the current fixture pipelines, not a common
input format. Final-object destruction is excluded from construction timing.

## Reproduction

Run from the repository root with CMake and Python with Matplotlib:

```bash
cmake --preset benchmarks
cmake --build --preset benchmarks --target bp_tree_benchmarks dfuds_tree_benchmarks louds_tree_benchmarks -j 2
out=build/tree-reproduction
mkdir -p "$out"

python3 scripts/benchmark_trees.py --succinct-mode best --scope succinct \
  --cache-bytes 524288 --exponents 10 12 14 16 18 20 22 24 26 --random-interleaving \
  --repetitions 5 --min-time 0.2 --warmup 0.1 --output-dir "$out/512"
python3 scripts/benchmark_trees.py --succinct-mode best --scope succinct \
  --cache-bytes 32768 --exponents 10 12 14 16 18 20 22 24 26 --random-interleaving \
  --repetitions 5 --min-time 0.2 --warmup 0.1 --output-dir "$out/32"
python3 scripts/benchmark_trees.py --scope references \
  --exponents 10 12 14 16 18 20 22 24 26 --random-interleaving \
  --repetitions 5 --min-time 0.2 --warmup 0.1 --output-dir "$out/references"

c++ -O3 -DNDEBUG -march=native -std=c++20 -Wall -Wextra -Werror -Iinclude \
  src/benchmarks/tree_construction_probe.cpp -o "$out/construction-probe"
python3 scripts/benchmark_tree_construction.py --binary "$out/construction-probe" \
  --output-dir "$out/construction" --exponents 10 12 14 16 18 20 22 24 26

python3 scripts/plot_tree_docs.py --reference "$out/references/combined.json" \
  --cache32 "$out/32/combined.json" --cache512 "$out/512/combined.json" \
  --construction "$out/construction/results.json" --output-dir "$out/graphs"
```

## Results

<img src="tree_dfs.svg" alt="DFS at 32 KiB and 512 KiB" />

<img src="tree_dfs_children.svg" alt="DFS with child inspection" />

<img src="tree_random_path.svg" alt="Root-to-leaf paths" />

<img src="tree_random_path_children.svg" alt="Paths with child inspection" />

<img src="tree_build.svg" alt="Construction from prepared input and full setup" />

<img src="tree_memory.svg" alt="Memory usage including reserved traversal scratch" />
