# Trees {#benchmarks_trees}

These benchmarks compare BP, DFUDS, LOUDS and ForestCover on three tasks:
walking a whole tree, descending from the root while inspecting children, and
returning the path between two random nodes. Adjacency lists and pointers
provide reference times for the last task.

## DFS and root-to-leaf benchmarks

`tree_dfs_benchmarks` contains the `dfs` and `dfs_children` families.
`tree_path_benchmarks` contains `random_path` and `random_path_children`.
Both start at the root and choose a child at each step until they reach a leaf.
These four families include BP, DFUDS, and LOUDS in optimized (`*_best`) and direct
(`*_direct`) modes, plus adjacency-list and pointer references in BFS and
preorder layouts. Names start with the workload, for example
`dfs/bp_best/tree_size:1024`.

To run one family at one size after building:

```bash
taskset -c 0 build/benchmarks/tree_dfs_benchmarks \
  --benchmark_filter='^dfs/[^/]+/tree_size:1024$'
taskset -c 0 build/benchmarks/tree_path_benchmarks \
  --benchmark_filter='^random_path_children/[^/]+/tree_size:1024$'
```

The collector (`scripts/benchmark_trees.py`) accepts `--workload dfs`, `--workload path`, or `--workload all`
(the default). `--succinct-mode both` collects optimized and direct modes in
one dataset. It runs each implementation and size in a separate process for
each binary, so large reference fixtures are released between groups.

Scan/replay, conversion, and prepared-tree diagnostics are available through
the optional `tree_experimental_benchmarks` target.

## Matched tree benchmarks

For the three-shape comparison, `tree_dfs_benchmarks` registers the inspected
DFS rows. `tree_path_benchmarks` registers the inspected root-to-leaf rows and
the paths between two random or sibling endpoints. To run one row of each kind:

```bash
taskset -c 0 build/benchmarks/tree_dfs_benchmarks \
  --benchmark_filter='^tree_shape/bp/inspected_dfs/random_recursive/N:1024/seed:42/budget:32768$'
taskset -c 0 build/benchmarks/tree_path_benchmarks \
  --benchmark_filter='^tree_shape/bp/inspected_root_path/random_recursive/N:1024/seed:42/budget:32768$'
taskset -c 0 build/benchmarks/tree_path_benchmarks \
  --benchmark_filter='^tree_path/forest_cover/random/random_recursive/N:1024/seed:42$'
```

\htmlonly
<h2>Interactive charts</h2>
<p>Select implementations and traversal budgets to compare. Each point gives
equal weight to the three tree families and appears only when all three have measurements.</p>
<ul>
<li><a href="inspected_dfs_observable.html">Inspected DFS</a>: CPU nanoseconds per node, linear scale.</li>
<li><a href="inspected_root_path_observable.html">Inspected root-to-leaf paths</a>: CPU nanoseconds per path, logarithmic scale.</li>
<li><a href="tree_path_observable.html">Random endpoint paths</a>: CPU nanoseconds per path, logarithmic scale.</li>
</ul>
\endhtmlonly

Use **Visible curves** to search and select implementations and budgets. Hover
over a point for its latency and, for either path chart, its mean node count.
Scroll to zoom, drag to pan and use **Reset** to restore the axes. The HTML
pages include their data and load their scripts locally. Each inspected chart
has 12 curves; the random-endpoint chart adds two references for a total of
14. All cover nine sizes, every even exponent from 2^10 through 2^26 nodes.

## Random endpoint paths

`tree.tree_path(u, v, out, scratch)` writes the unique path from `u` to `v`
in order, including both endpoints and their common ancestor exactly once.
The benchmark replays the same 256 seed-42 endpoint pairs for each encoding
of a given source tree. It times materializing the paths into reused vectors.
Tree construction, endpoint translation, buffer reservation and checks against
the source parent/depth oracle are outside the timed loop. Google Benchmark
counts each pair as one query, so the plotted time is already **ns/path**.
Sibling-only queries and root-to-leaf walks are separate workloads.

The budget is for extra navigation state, not for the tree or the two path
buffers. Baseline BP, DFUDS and LOUDS use direct parent queries; their 32 KiB
and 512 KiB curves cache parent handles for a prefix of the native numbering
and fall back to direct navigation elsewhere. This prefix cache is a benchmark
helper, not `pixie::tree::RootNavigationIndex`. ForestCover's baseline uses one
decoder slot; its budgeted readers fit the context and slots in the stated
allowance. Each reference has one curve and no extra navigation cache.

At 2^26 nodes the three-family mean path had 7,698 nodes. The pointer reference
took 117 microseconds per path, one-slot ForestCover 223, BP with 512 KiB 225,
and the adjacency list 302. The mean spans very different shapes: the same 256
pairs averaged 35 nodes on random recursive trees, 15,351 on binary Catalan
trees and 7,708 on ordered plane Catalan trees. The respective paths touched
22, 883 and 856 distinct forest units on average. A path may leave and re-enter
a unit, and the two endpoint climbs issue parent lookups in a different order
from the returned nodes. Distinct units and contiguous runs in the output are
structural counts, not decoder-cache misses.

To regenerate the three interactive pages from completed Google Benchmark
JSON, set `inspected` to the directory of `*-timing-*.json` shards containing
only the two `tree_shape/` workloads from both executables at all nine sizes,
`initial` to the initial random-endpoint JSON and `extensions` to the directory
of additional `p<exponent>.json` summaries or per-family shards:

```bash
python3 scripts/plot_tree_curves.py --inspected-dir "$inspected" \
  --random-initial "$initial" --random-shards "$extensions"
```

The generator writes `*_observable.html` under
`src/docs/images/benchmarks/trees/preview`. Doxygen copies them and their local
assets beside this page. Only complete three-family means produce a dot.

## Which objects are public?

The four measured encodings are `pixie::tree::BPTree`, `DFUDSTree`,
`LoudsTree` and `ForestCoverTree`, all in `pixie::tree`. Their shared
`TreeBase` contract exposes node handles, positioned navigation views and
endpoint paths. `RootNavigationIndex` and `TreeTraversalContext` hold reusable
navigation state outside the tree; `TreePart` owns inspection snapshots.
The decomposed and layered trees, prepared replay objects and SIMD traversal
prototypes in `pixie::tree::experimental` are opt-in experiments, not names
for the promoted forest cover. The adjacency-list and pointer references and
benchmark workloads are outside the public family. See \ref rooted_tree_public_objects
for the header and type inventory.

## Matched tree shape latencies

For each implementation, workload, node count, seed and traversal budget, the
table reports the arithmetic mean of three median CPU latencies: one each for
`random_recursive`, `binary_catalan` and `ordered_plane_catalan`. Each family
has equal weight regardless of its repetition count. Times are nanoseconds
per **whole inspected DFS** or per **inspected root-to-leaf path**. The latter
uses the measured iteration time divided by `paths_per_iteration`.

<!-- tree_shape results begin -->

\image html tree_shape_inspected_dfs.svg "Whole inspected DFS CPU latency by traversal budget"

\image html tree_shape_inspected_root_path.svg "Inspected root-to-leaf path CPU latency by traversal budget"

| Budget (bytes) | Nodes | Seed | Implementation | Workload | Mean CPU latency (ns) |
| ---: | ---: | ---: | --- | --- | ---: |
| 0 | 1024 | 42 | bp | inspected_dfs | 29132.1 |
| 0 | 1024 | 42 | bp | inspected_root_path | 503.093 |
| 0 | 1024 | 42 | dfuds | inspected_dfs | 56742.8 |
| 0 | 1024 | 42 | dfuds | inspected_root_path | 599.474 |
| 0 | 1024 | 42 | forest_cover | inspected_dfs | 3393.5 |
| 0 | 1024 | 42 | forest_cover | inspected_root_path | 35.6784 |
| 0 | 1024 | 42 | louds | inspected_dfs | 42769.9 |
| 0 | 1024 | 42 | louds | inspected_root_path | 101.264 |
| 0 | 4096 | 42 | bp | inspected_dfs | 129938 |
| 0 | 4096 | 42 | bp | inspected_root_path | 553.981 |
| 0 | 4096 | 42 | dfuds | inspected_dfs | 214030 |
| 0 | 4096 | 42 | dfuds | inspected_root_path | 492.352 |
| 0 | 4096 | 42 | forest_cover | inspected_dfs | 15339.4 |
| 0 | 4096 | 42 | forest_cover | inspected_root_path | 56.3353 |
| 0 | 4096 | 42 | louds | inspected_dfs | 199476 |
| 0 | 4096 | 42 | louds | inspected_root_path | 88.3361 |
| 0 | 16384 | 42 | bp | inspected_dfs | 721030 |
| 0 | 16384 | 42 | bp | inspected_root_path | 1306.52 |
| 0 | 16384 | 42 | dfuds | inspected_dfs | 951820 |
| 0 | 16384 | 42 | dfuds | inspected_root_path | 440.843 |
| 0 | 16384 | 42 | forest_cover | inspected_dfs | 77375.2 |
| 0 | 16384 | 42 | forest_cover | inspected_root_path | 75.7234 |
| 0 | 16384 | 42 | louds | inspected_dfs | 979340 |
| 0 | 16384 | 42 | louds | inspected_root_path | 167.435 |
| 0 | 65536 | 42 | bp | inspected_dfs | 3.0744e+06 |
| 0 | 65536 | 42 | bp | inspected_root_path | 1939.97 |
| 0 | 65536 | 42 | dfuds | inspected_dfs | 3.8389e+06 |
| 0 | 65536 | 42 | dfuds | inspected_root_path | 822.659 |
| 0 | 65536 | 42 | forest_cover | inspected_dfs | 700235 |
| 0 | 65536 | 42 | forest_cover | inspected_root_path | 106.781 |
| 0 | 65536 | 42 | louds | inspected_dfs | 3.97334e+06 |
| 0 | 65536 | 42 | louds | inspected_root_path | 222.241 |
| 0 | 262144 | 42 | bp | inspected_dfs | 1.23464e+07 |
| 0 | 262144 | 42 | bp | inspected_root_path | 1606.9 |
| 0 | 262144 | 42 | dfuds | inspected_dfs | 1.54532e+07 |
| 0 | 262144 | 42 | dfuds | inspected_root_path | 886.422 |
| 0 | 262144 | 42 | forest_cover | inspected_dfs | 2.93003e+06 |
| 0 | 262144 | 42 | forest_cover | inspected_root_path | 91.9272 |
| 0 | 262144 | 42 | louds | inspected_dfs | 1.74321e+07 |
| 0 | 262144 | 42 | louds | inspected_root_path | 186.703 |
| 0 | 1048576 | 42 | bp | inspected_dfs | 5.00449e+07 |
| 0 | 1048576 | 42 | bp | inspected_root_path | 2440.29 |
| 0 | 1048576 | 42 | dfuds | inspected_dfs | 6.22384e+07 |
| 0 | 1048576 | 42 | dfuds | inspected_root_path | 830.489 |
| 0 | 1048576 | 42 | forest_cover | inspected_dfs | 1.21293e+07 |
| 0 | 1048576 | 42 | forest_cover | inspected_root_path | 100.751 |
| 0 | 1048576 | 42 | louds | inspected_dfs | 7.16627e+07 |
| 0 | 1048576 | 42 | louds | inspected_root_path | 168.167 |
| 0 | 4194304 | 42 | bp | inspected_dfs | 2.00553e+08 |
| 0 | 4194304 | 42 | bp | inspected_root_path | 3044.1 |
| 0 | 4194304 | 42 | dfuds | inspected_dfs | 2.51445e+08 |
| 0 | 4194304 | 42 | dfuds | inspected_root_path | 960.693 |
| 0 | 4194304 | 42 | forest_cover | inspected_dfs | 5.22174e+07 |
| 0 | 4194304 | 42 | forest_cover | inspected_root_path | 108.066 |
| 0 | 4194304 | 42 | louds | inspected_dfs | 2.91167e+08 |
| 0 | 4194304 | 42 | louds | inspected_root_path | 192.755 |
| 0 | 16777216 | 42 | bp | inspected_dfs | 9.35543e+08 |
| 0 | 16777216 | 42 | bp | inspected_root_path | 3687.67 |
| 0 | 16777216 | 42 | dfuds | inspected_dfs | 1.14521e+09 |
| 0 | 16777216 | 42 | dfuds | inspected_root_path | 1421.54 |
| 0 | 16777216 | 42 | forest_cover | inspected_dfs | 2.3326e+08 |
| 0 | 16777216 | 42 | forest_cover | inspected_root_path | 166.943 |
| 0 | 16777216 | 42 | louds | inspected_dfs | 1.31988e+09 |
| 0 | 16777216 | 42 | louds | inspected_root_path | 212.903 |
| 0 | 67108864 | 42 | bp | inspected_dfs | 3.4061e+09 |
| 0 | 67108864 | 42 | bp | inspected_root_path | 4103.81 |
| 0 | 67108864 | 42 | dfuds | inspected_dfs | 3.97459e+09 |
| 0 | 67108864 | 42 | dfuds | inspected_root_path | 1506.86 |
| 0 | 67108864 | 42 | forest_cover | inspected_dfs | 8.75888e+08 |
| 0 | 67108864 | 42 | forest_cover | inspected_root_path | 206.332 |
| 0 | 67108864 | 42 | louds | inspected_dfs | 4.81125e+09 |
| 0 | 67108864 | 42 | louds | inspected_root_path | 276.579 |
| 32768 | 1024 | 42 | bp | inspected_dfs | 9224.42 |
| 32768 | 1024 | 42 | bp | inspected_root_path | 11.7421 |
| 32768 | 1024 | 42 | dfuds | inspected_dfs | 14386.5 |
| 32768 | 1024 | 42 | dfuds | inspected_root_path | 12.2639 |
| 32768 | 1024 | 42 | forest_cover | inspected_dfs | 3236.41 |
| 32768 | 1024 | 42 | forest_cover | inspected_root_path | 10.3625 |
| 32768 | 1024 | 42 | louds | inspected_dfs | 24404.5 |
| 32768 | 1024 | 42 | louds | inspected_root_path | 12.1214 |
| 32768 | 4096 | 42 | bp | inspected_dfs | 45960.9 |
| 32768 | 4096 | 42 | bp | inspected_root_path | 11.5294 |
| 32768 | 4096 | 42 | dfuds | inspected_dfs | 51075.9 |
| 32768 | 4096 | 42 | dfuds | inspected_root_path | 12.5133 |
| 32768 | 4096 | 42 | forest_cover | inspected_dfs | 14935 |
| 32768 | 4096 | 42 | forest_cover | inspected_root_path | 9.74907 |
| 32768 | 4096 | 42 | louds | inspected_dfs | 133056 |
| 32768 | 4096 | 42 | louds | inspected_root_path | 11.223 |
| 32768 | 16384 | 42 | bp | inspected_dfs | 230229 |
| 32768 | 16384 | 42 | bp | inspected_root_path | 26.6223 |
| 32768 | 16384 | 42 | dfuds | inspected_dfs | 211803 |
| 32768 | 16384 | 42 | dfuds | inspected_root_path | 27.9028 |
| 32768 | 16384 | 42 | forest_cover | inspected_dfs | 72668.8 |
| 32768 | 16384 | 42 | forest_cover | inspected_root_path | 21.1549 |
| 32768 | 16384 | 42 | louds | inspected_dfs | 668734 |
| 32768 | 16384 | 42 | louds | inspected_root_path | 25.7442 |
| 32768 | 65536 | 42 | bp | inspected_dfs | 1.22453e+06 |
| 32768 | 65536 | 42 | bp | inspected_root_path | 38.6548 |
| 32768 | 65536 | 42 | dfuds | inspected_dfs | 1.19656e+06 |
| 32768 | 65536 | 42 | dfuds | inspected_root_path | 40.5399 |
| 32768 | 65536 | 42 | forest_cover | inspected_dfs | 709812 |
| 32768 | 65536 | 42 | forest_cover | inspected_root_path | 29.8647 |
| 32768 | 65536 | 42 | louds | inspected_dfs | 2.79026e+06 |
| 32768 | 65536 | 42 | louds | inspected_root_path | 36.2386 |
| 32768 | 262144 | 42 | bp | inspected_dfs | 5.02259e+06 |
| 32768 | 262144 | 42 | bp | inspected_root_path | 36.0283 |
| 32768 | 262144 | 42 | dfuds | inspected_dfs | 4.88291e+06 |
| 32768 | 262144 | 42 | dfuds | inspected_root_path | 38.1649 |
| 32768 | 262144 | 42 | forest_cover | inspected_dfs | 2.9276e+06 |
| 32768 | 262144 | 42 | forest_cover | inspected_root_path | 23.7927 |
| 32768 | 262144 | 42 | louds | inspected_dfs | 1.2414e+07 |
| 32768 | 262144 | 42 | louds | inspected_root_path | 30.168 |
| 32768 | 1048576 | 42 | bp | inspected_dfs | 2.05222e+07 |
| 32768 | 1048576 | 42 | bp | inspected_root_path | 47.0818 |
| 32768 | 1048576 | 42 | dfuds | inspected_dfs | 1.99469e+07 |
| 32768 | 1048576 | 42 | dfuds | inspected_root_path | 44.3487 |
| 32768 | 1048576 | 42 | forest_cover | inspected_dfs | 1.21641e+07 |
| 32768 | 1048576 | 42 | forest_cover | inspected_root_path | 27.9044 |
| 32768 | 1048576 | 42 | louds | inspected_dfs | 5.14525e+07 |
| 32768 | 1048576 | 42 | louds | inspected_root_path | 30.0683 |
| 32768 | 4194304 | 42 | bp | inspected_dfs | 8.44862e+07 |
| 32768 | 4194304 | 42 | bp | inspected_root_path | 73.8318 |
| 32768 | 4194304 | 42 | dfuds | inspected_dfs | 8.17902e+07 |
| 32768 | 4194304 | 42 | dfuds | inspected_root_path | 62.6895 |
| 32768 | 4194304 | 42 | forest_cover | inspected_dfs | 5.21235e+07 |
| 32768 | 4194304 | 42 | forest_cover | inspected_root_path | 34.008 |
| 32768 | 4194304 | 42 | louds | inspected_dfs | 2.09446e+08 |
| 32768 | 4194304 | 42 | louds | inspected_root_path | 36.2962 |
| 32768 | 16777216 | 42 | bp | inspected_dfs | 3.90453e+08 |
| 32768 | 16777216 | 42 | bp | inspected_root_path | 124.827 |
| 32768 | 16777216 | 42 | dfuds | inspected_dfs | 3.59253e+08 |
| 32768 | 16777216 | 42 | dfuds | inspected_root_path | 93.2989 |
| 32768 | 16777216 | 42 | forest_cover | inspected_dfs | 2.33459e+08 |
| 32768 | 16777216 | 42 | forest_cover | inspected_root_path | 41.8304 |
| 32768 | 16777216 | 42 | louds | inspected_dfs | 9.16325e+08 |
| 32768 | 16777216 | 42 | louds | inspected_root_path | 42.9839 |
| 32768 | 67108864 | 42 | bp | inspected_dfs | 1.47661e+09 |
| 32768 | 67108864 | 42 | bp | inspected_root_path | 167.233 |
| 32768 | 67108864 | 42 | dfuds | inspected_dfs | 1.39675e+09 |
| 32768 | 67108864 | 42 | dfuds | inspected_root_path | 127.575 |
| 32768 | 67108864 | 42 | forest_cover | inspected_dfs | 8.66079e+08 |
| 32768 | 67108864 | 42 | forest_cover | inspected_root_path | 47.5249 |
| 32768 | 67108864 | 42 | louds | inspected_dfs | 3.81994e+09 |
| 32768 | 67108864 | 42 | louds | inspected_root_path | 48.9679 |
| 524288 | 1024 | 42 | bp | inspected_dfs | 3800.26 |
| 524288 | 1024 | 42 | bp | inspected_root_path | 11.8473 |
| 524288 | 1024 | 42 | dfuds | inspected_dfs | 4552.13 |
| 524288 | 1024 | 42 | dfuds | inspected_root_path | 12.1942 |
| 524288 | 1024 | 42 | forest_cover | inspected_dfs | 3057.33 |
| 524288 | 1024 | 42 | forest_cover | inspected_root_path | 10.5007 |
| 524288 | 1024 | 42 | louds | inspected_dfs | 5483.32 |
| 524288 | 1024 | 42 | louds | inspected_root_path | 11.7909 |
| 524288 | 4096 | 42 | bp | inspected_dfs | 15054.7 |
| 524288 | 4096 | 42 | bp | inspected_root_path | 10.4494 |
| 524288 | 4096 | 42 | dfuds | inspected_dfs | 17560.3 |
| 524288 | 4096 | 42 | dfuds | inspected_root_path | 11.3849 |
| 524288 | 4096 | 42 | forest_cover | inspected_dfs | 14657.6 |
| 524288 | 4096 | 42 | forest_cover | inspected_root_path | 9.8886 |
| 524288 | 4096 | 42 | louds | inspected_dfs | 44489.3 |
| 524288 | 4096 | 42 | louds | inspected_root_path | 10.549 |
| 524288 | 16384 | 42 | bp | inspected_dfs | 141519 |
| 524288 | 16384 | 42 | bp | inspected_root_path | 23.204 |
| 524288 | 16384 | 42 | dfuds | inspected_dfs | 158509 |
| 524288 | 16384 | 42 | dfuds | inspected_root_path | 23.3041 |
| 524288 | 16384 | 42 | forest_cover | inspected_dfs | 66664.1 |
| 524288 | 16384 | 42 | forest_cover | inspected_root_path | 19.2336 |
| 524288 | 16384 | 42 | louds | inspected_dfs | 378032 |
| 524288 | 16384 | 42 | louds | inspected_root_path | 23.065 |
| 524288 | 65536 | 42 | bp | inspected_dfs | 889953 |
| 524288 | 65536 | 42 | bp | inspected_root_path | 29.9941 |
| 524288 | 65536 | 42 | dfuds | inspected_dfs | 984392 |
| 524288 | 65536 | 42 | dfuds | inspected_root_path | 29.3732 |
| 524288 | 65536 | 42 | forest_cover | inspected_dfs | 700212 |
| 524288 | 65536 | 42 | forest_cover | inspected_root_path | 25.9267 |
| 524288 | 65536 | 42 | louds | inspected_dfs | 2.60852e+06 |
| 524288 | 65536 | 42 | louds | inspected_root_path | 29.9987 |
| 524288 | 262144 | 42 | bp | inspected_dfs | 3.75648e+06 |
| 524288 | 262144 | 42 | bp | inspected_root_path | 21.04 |
| 524288 | 262144 | 42 | dfuds | inspected_dfs | 4.12804e+06 |
| 524288 | 262144 | 42 | dfuds | inspected_root_path | 20.9579 |
| 524288 | 262144 | 42 | forest_cover | inspected_dfs | 2.92103e+06 |
| 524288 | 262144 | 42 | forest_cover | inspected_root_path | 18.4847 |
| 524288 | 262144 | 42 | louds | inspected_dfs | 1.18766e+07 |
| 524288 | 262144 | 42 | louds | inspected_root_path | 20.3468 |
| 524288 | 1048576 | 42 | bp | inspected_dfs | 1.49846e+07 |
| 524288 | 1048576 | 42 | bp | inspected_root_path | 18.8346 |
| 524288 | 1048576 | 42 | dfuds | inspected_dfs | 1.66857e+07 |
| 524288 | 1048576 | 42 | dfuds | inspected_root_path | 19.9204 |
| 524288 | 1048576 | 42 | forest_cover | inspected_dfs | 1.19687e+07 |
| 524288 | 1048576 | 42 | forest_cover | inspected_root_path | 17.5048 |
| 524288 | 1048576 | 42 | louds | inspected_dfs | 5.25008e+07 |
| 524288 | 1048576 | 42 | louds | inspected_root_path | 18.7672 |
| 524288 | 4194304 | 42 | bp | inspected_dfs | 6.24187e+07 |
| 524288 | 4194304 | 42 | bp | inspected_root_path | 25.2052 |
| 524288 | 4194304 | 42 | dfuds | inspected_dfs | 6.83533e+07 |
| 524288 | 4194304 | 42 | dfuds | inspected_root_path | 27.4901 |
| 524288 | 4194304 | 42 | forest_cover | inspected_dfs | 5.14995e+07 |
| 524288 | 4194304 | 42 | forest_cover | inspected_root_path | 22.0059 |
| 524288 | 4194304 | 42 | louds | inspected_dfs | 2.18765e+08 |
| 524288 | 4194304 | 42 | louds | inspected_root_path | 24.887 |
| 524288 | 16777216 | 42 | bp | inspected_dfs | 2.94751e+08 |
| 524288 | 16777216 | 42 | bp | inspected_root_path | 26.4437 |
| 524288 | 16777216 | 42 | dfuds | inspected_dfs | 3.12837e+08 |
| 524288 | 16777216 | 42 | dfuds | inspected_root_path | 28.8054 |
| 524288 | 16777216 | 42 | forest_cover | inspected_dfs | 2.35586e+08 |
| 524288 | 16777216 | 42 | forest_cover | inspected_root_path | 24.7512 |
| 524288 | 16777216 | 42 | louds | inspected_dfs | 9.07791e+08 |
| 524288 | 16777216 | 42 | louds | inspected_root_path | 25.5908 |
| 524288 | 67108864 | 42 | bp | inspected_dfs | 1.04496e+09 |
| 524288 | 67108864 | 42 | bp | inspected_root_path | 35.6297 |
| 524288 | 67108864 | 42 | dfuds | inspected_dfs | 1.07454e+09 |
| 524288 | 67108864 | 42 | dfuds | inspected_root_path | 35.456 |
| 524288 | 67108864 | 42 | forest_cover | inspected_dfs | 8.76486e+08 |
| 524288 | 67108864 | 42 | forest_cover | inspected_root_path | 30.0865 |
| 524288 | 67108864 | 42 | louds | inspected_dfs | 3.29095e+09 |
| 524288 | 67108864 | 42 | louds | inspected_root_path | 31.184 |

<!-- tree_shape results end -->

To publish the table and both SVGs, put the inspected timing shards for both
workloads from one campaign in `$out`, then run from the repository root with Matplotlib:

```bash
python3 scripts/report_tree_shapes.py --page src/docs/benchmarks/trees.md \
  "$out"/*-timing-*.json
```

Without `--page`, the script prints complete groups and identifies incomplete
ones on stderr. Publication requires every even power from 2^10 through 2^26
at seed 42. Each graph compares BP, DFUDS, LOUDS and forest cover in separate
panels for 0, 32 and 512 KiB. The plots use milliseconds per whole DFS and
nanoseconds per path. The historical single-family plots below are unchanged.

## Historical plot machine

AMD Ryzen 7 8845HS (8 cores, 16 threads, 16 MiB L3), 15 GiB RAM,
WSL2 (Linux 6.6.87.2-microsoft-standard-WSL2). GCC 13.3, C++20,
`-O3 -DNDEBUG -march=native`; measurements pinned to CPU 0.
Results captured October 30, 2026, at commit `8e72c47`, after the ownership
refactor, before the subsequent owning-encoding change. These graphs retain
that measured revision. BP, DFUDS and LOUDS now own their encodings by default.
Their `pixie::tree::TreeTraversalContext` takes a byte limit for the context;
an optional root index has a separate budget.

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
cmake --build --preset benchmarks --target tree_dfs_benchmarks tree_path_benchmarks -j 2
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
  --construction-mode consume --output-dir "$out/construction" --exponents 10 12 14 16 18 20 22 24 26

python3 scripts/plot_tree_docs.py --reference "$out/references/combined.json" \
  --cache32 "$out/32/combined.json" --cache512 "$out/512/combined.json" \
  --construction "$out/construction/results.json" --output-dir "$out/graphs"
```

Use `--construction-mode copy` in a separate output directory to include copying
prepared encoding words in construction timing. The consuming mode prepares aligned
storage before construction, and includes that preparation in full setup.

## Historical single-family results

<img src="tree_dfs.svg" alt="DFS at 32 KiB and 512 KiB" />

<img src="tree_dfs_children.svg" alt="DFS with child inspection" />

<img src="tree_random_path.svg" alt="Root-to-leaf paths" />

<img src="tree_random_path_children.svg" alt="Paths with child inspection" />

<img src="tree_build.svg" alt="Construction from prepared input and full setup" />

<img src="tree_memory.svg" alt="Memory usage including reserved traversal scratch" />
