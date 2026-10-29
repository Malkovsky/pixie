# Succinct-tree experiment notes

Consolidated 2026-10-28, before cleanup of the local profiling workspace.
These notes preserve decisions, algorithms and representative evidence from
the BP/DFUDS/LOUDS experiments. They are not an API stability promise or a claim
that discarded prototypes remain reproducible after their sources are deleted.
No experiment files were deleted while preparing these notes.

## Read first

- [Navigation and caching](navigation-and-caching.md): local query reuse,
  bounded root index, CLOCK, retained unpacking and runtime memory budgets.
- [Decoding and SIMD](decoding-and-simd.md): streaming DFS, compiled topology
  programs, retained optimizations, negative results and 128-node SIMD work.
- [Measurement and construction](measurement-and-construction.md): workload
  semantics, current rounded results, input generation, space accounting and
  validation requirements.
- [Archived SIMD prototype](microtree-simd/README.md): the 128-node BP kernels,
  exact-trace checks and optional timing probes extracted from profiling.

The public-facing snapshot is [the tree benchmark page](../../../../../src/docs/benchmarks/trees.md).
The notes below are deliberately self-contained: historical report names are
provenance labels, not links that must survive cleanup. Numerical comparisons
within each historical experiment use its own controls; do not combine them
into one synthetic speedup or compare their absolute times across runs.

## Decisions worth retaining

1. **Keep three layers separate.** A local node view memoizes navigation
   results; a bounded root-connected index accelerates repeated paths;
   temporary decoding accelerates scans. A general eviction cache is a fourth,
   optional mechanism, not a free substitute for those layers.
2. **Optimize whole workloads.** Fitting succinct storage in L3 does not remove
   dependent navigation, decoding, dispatch, callbacks or stack bookkeeping.
   Smaller structures and fewer cache misses did not guarantee lower time.
3. **Use the same configuration for all four workloads.** The current adapter
   chooses an encoding-specific configuration from the constructor's runtime
   byte limit. It does not pick a different memory budget for paths and DFS.
4. **Charge temporary decoder storage.** The 512 KiB/32 KiB allowances cover
   the navigation context/root index plus explicit traversal buffers. Immutable
   constexpr tables are excluded under the agreed accounting convention.
5. **Retain compiled seven-node programs inside larger chunks.** Enumerating
   more shapes, materializing generic event streams and replacing the decoder
   with arbitrary 128-node SIMD matching all failed to improve large random
   inspected DFS in the tested forms.
6. **Preserve the benchmarks' actual observations.** Every entry, nonleaf
   return and inspected-child handle matters; count-only or weakened observers
   can make a candidate look artificially fast.
7. **Keep both references and both layouts.** Vectors and three pointers win
   in different scenarios. Choosing one convenient layout can hide a loss.

The objective of beating the fastest reference on all four workloads is
**not achieved**. BP is particularly strong in plain DFS and is the fastest
succinct inspected DFS in the current large-tree snapshot. LOUDS is strongest
on the indexed paths. Inspected DFS remains behind BFS pointers; LOUDS plain
DFS also remains behind pointers. No single encoding wins every workload.

## Current code map

| Responsibility | Location relative to `include/pixie/` |
| --- | --- |
| Tree facade and child enumeration contract | `tree.h` |
| Encodings and local view policies | `tree/{bp,dfuds,louds,cache}.h` |
| Lazy per-node results | `tree/detail/{bp,dfuds,louds}_node_view.h` |
| Streaming kernels and bit primitives | `tree/detail/unpack*.h` |
| Bounded root index | `tree/experimental/tree_navigation_index.h` |
| Combined runtime budget selection | `tree/experimental/tree_traversal_cache.h` |
| BP compact records and vector program arguments | `tree/experimental/tree_micro_compact.h`, `tree_micro_vector.h` |
| Shared small BP programs | `tree/experimental/tree_micro_prefetch.h` |
| DFUDS programs/records | `tree/experimental/tree_dfuds_micro.h` |
| LOUDS level and child-end buffers | `tree/experimental/tree_louds_prefetch.h`, `tree_louds_iterative.h` |
| Optional general caches | `tree/experimental/second_chance_tree_cache.h`, `unpacked_tree_cache.h` |

Repository tools retained for reproduction: `scripts/benchmark_trees.py`,
`scripts/benchmark_tree_construction.py`, `scripts/plot_tree_docs.py`, and
`src/benchmarks/tree_construction_probe.cpp`. Input generation is in
`src/benchmarks/tree_input.h`; exact-event and navigation tests are in
`src/tests/tree_tests.cpp`. Some additional experimental headers remain from
earlier candidates; their presence does not mean the final adapter selects them.

## Next iteration: hypotheses, not established wins

- Consider sharing/reclaiming scratch and root-index capacity across operation
  lifetimes, without breaking the budget or reentrant-callback semantics.
  Fixed reservation currently costs 32 KiB path performance.
- A wide/shallow-shape SIMD specialization may be worth a cheap admission test
  inside the existing large-chunk decoder. Do not switch every subtree to the
  losing generic 256-bit decoder.
- A general cache might benefit from workload-aware admission or complete-child
  coverage, but must be compared with the root index, not only an expensive
  uncached baseline. The root index assumes uniform child choices.
- Continue reducing inspected-DFS decode/replay work. Extra buffers, smaller
  records, fewer instructions or wider SIMD are only hypotheses until integrated
  timings and exact-event checks support them.
- Recheck on other machines and tree shapes before generalizing. Primary
  performance data are one random-recursive topology seed on one x86 host.
