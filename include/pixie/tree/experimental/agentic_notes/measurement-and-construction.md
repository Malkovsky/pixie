# Measurement, construction and validation

## Final measurement contract

The latest completed grid uses every even exponent from **2^10 to 2^26**;
10 million nodes is a historical experiment size, not an extra point in the
final sweep. Both 32 KiB and 512 KiB use the same current runtime-adaptive
implementation binaries across the grid. The earlier decoder-disabled 32 KiB
run and pre-adaptation 512 KiB small-size run are not the final succinct data.

Four traversal workloads must retain their semantics:

1. Plain DFS consumes every preorder entry and one return per nonleaf after
   its children. DFS units per node use the input node count, not event count.
2. Inspected DFS also consumes each immediate child handle before descent.
3. Random root-to-leaf paths choose a child uniformly at each nonleaf.
4. Inspected paths enumerate all immediate children before choosing/descent;
   observing the child list is part of timed work.

Topology is the seed-42 random-recursive tree; paths use a shifted 64 KiB
random pool with seed 43. BP/DFUDS number naturally in preorder, LOUDS in BFS.
Both vector and three-pointer references are measured in both layouts.
Index/cache construction is outside query timing. The bounded root index uses
topology probabilities, not sampled query frequencies.

Host: AMD Ryzen 7 8845HS, 16 MiB L3, 15 GiB RAM, WSL2/Linux 6.6.87.2,
GCC 13.3, C++20, Release `-O3 -DNDEBUG -march=native`, pinned CPU 0.
Traversal uses five randomly interleaved repetitions, minimum 0.2 s and
0.1 s warmup. Medians are CPU time. Construction uses three fresh-process
repetitions in shuffled order and steady-clock wall time. These clocks are
not interchangeable. Construction is not Google Benchmark query timing.

Completion counts: 216 succinct traversal cases / 1080 raw repetitions over
both budgets; 90 construction configurations / 270 samples including the four
reference/layout combinations. The reference query curves come from the
earlier full reference sweep; no claim of same-process paired measurement
against those reference curves is intended.

## Rounded current large-tree checkpoint

At 2^26 = 67,108,864 nodes. DFS columns are milliseconds per traversal; path
columns are nanoseconds per root-to-leaf path. References are reused from the
full 512 KiB comparison dataset; their own implementations have no cache budget.

| Implementation | Budget | DFS ms | DFS + children ms | Path ns | Path + children ns |
| --- | --- | ---: | ---: | ---: | ---: |
| BP | 32 KiB | 184.95 | 649.24 | 266.01 | 382.60 |
| DFUDS | 32 KiB | 318.13 | 773.05 | 214.86 | 328.38 |
| LOUDS | 32 KiB | 560.14 | 759.25 | 75.49 | 110.32 |
| BP | 512 KiB | 184.48 | 518.14 | 60.34 | 84.10 |
| DFUDS | 512 KiB | 320.72 | 582.52 | 58.85 | 85.54 |
| LOUDS | 512 KiB | 556.56 | 583.85 | 39.96 | 59.13 |
| Vectors preorder | — | 702.86 | 764.13 | 66.05 | 79.29 |
| Vectors BFS | — | 835.53 | 905.66 | 74.59 | 75.15 |
| Pointers preorder | — | 314.84 | 573.96 | 102.32 | 117.61 |
| Pointers BFS | — | 326.25 | 434.99 | 74.36 | 97.53 |

Interpretation: plain DFS needs little scratch, so 32 KiB retains essentially
the 512 KiB speed. Smaller inspected decoders are slower, and smaller root
indexes especially hurt paths. No whole-family/all-workload victory is implied.
For example, even BP's 512 KiB inspected DFS is slower than BFS pointers, and
its inspected path is slower than the best vector layout at this size.
Differences of only a few percent across separate runs need confirmation.

The final 1M/16M/67M adaptive-cache comparison had median repetition CVs about
0.91%/1.11% and maxima 2.67%/2.90% for 32/512 KiB. Do not transfer those
statistics to every earlier experiment or to the later small-size extension.
Some earlier controls had double-digit CVs and substantial drift.

## Construction: BP first, not LOUDS first

The compact benchmark generator reproduces `parent[i] = rng() % i` and the
same ascending child order as the adjacency generator:

1. Fill a uint32 parent array.
2. Insert nodes in descending order, reusing their parent slots as next-sibling
   links; a second uint32 array stores first children. Parent IDs are smaller
   than child IDs, so this in-place reuse does not overwrite an unread parent.
3. DFS emits the BP bitvector. Release both arrays afterward.
4. BP→DFUDS gathers preorder degrees in a uint32 array and writes degree ones
   followed by zero, with no extra initial one.
5. BP→LOUDS accumulates packed degree records per depth. Within a fixed depth,
   completion order equals opening order; concatenating levels yields BFS order
   without an O(n)-node BFS queue. LOUDS records are one then degree zeros.

This is **not** generation entirely in succinct space. Peak scratch is about
8 bytes/node plus depth state. Its supported benchmark domain is nonempty
sizes below 2^32−1, separate from the library's full-width coordinate goal.
Conversion needs no rank/select or RmM queries. Exact encoded-word differential
tests establish that the compact route preserves the earlier fixture topology.

Historical generation-only 10M medians were BP 10.71→0.76 s,
DFUDS 10.39→0.90 s and LOUDS 8.74→0.85 s; peak RSS fell from roughly
1.2 GiB to 82 MiB. These exclude tree-index construction and use their own
old adjacency-generation controls. They are not current query speedups.
Evidence: `tree-input-construction`, 2026-10-20.

References still use the ordinary adjacency generator, reorder to BFS or
preorder, then copy into the owning target. Thus full setup compares current
fixture pipelines, not equally optimized builders from a common input format.
Prepared-input construction starts from a bitvector for succinct trees and
an adjacency list for references. It includes the root-index/cache build.
Final-object/source destruction is excluded; intermediate cleanup inside
generation/conversion is included.

At 2^26 nodes, current median times and working storage are:

| Implementation | Budget | Full setup s | Prepared-input build s | Working bytes/node |
| --- | --- | ---: | ---: | ---: |
| BP | 32 KiB | 6.407 | 0.0857 | 0.2951 |
| DFUDS | 32 KiB | 7.780 | 0.0731 | 0.2951 |
| LOUDS | 32 KiB | 7.170 | 0.0044 | 0.2594 |
| BP | 512 KiB | 6.763 | 0.0951 | 0.3025 |
| DFUDS | 512 KiB | 7.616 | 0.0771 | 0.3025 |
| LOUDS | 512 KiB | 7.171 | 0.0064 | 0.2667 |
| Vectors BFS | — | 118.954 | 30.1384 | 40.0000 |
| Vectors preorder | — | 141.329 | 27.6468 | 40.0000 |
| Pointers BFS | — | 117.255 | 26.6148 | 24.0000 |
| Pointers preorder | — | 113.098 | 21.2792 | 24.0000 |

The largest source-generation cost explains most succinct full-setup time.
Large reference construction times reflect the current allocation/reordering
pipeline and should not be advertised as unavoidable costs of pointer trees.

## Memory accounting

- BP source length is 2n bits; the selected DFUDS/LOUDS conventions use 2n−1.
  Actual source buffers round to complete words and report vector capacity.
- Owned structure bytes include inline objects and reserved dynamic index
  storage. Nested objects are not counted twice. Under the selected local
  navigation policy there is no additional shared node cache inside the tree.
- Retained bytes add the external source buffer and persistent traversal
  context/root index. Working bytes additionally include reserved explicit
  traversal scratch. Scratch reservation is not the same as resident RSS.
- On the measured ABI, copied adjacency vectors account for about 40 bytes/node
  (24-byte vector control plus about two 8-byte entries). Three pointers account
  for 24 bytes/node. These are **bytes**, not bits, and omit allocator metadata.
- The combined cache allowance includes explicit child/ancestor/load buffers
  for one active traversal. It excludes source/RmM/rank-select storage,
  constexpr tables, allocator metadata, compiler scalar frame/spill overhead,
  caller callback state and the index-construction priority queue. Reentrant
  or concurrent walks need independent budget/lifetime reasoning.
- Peak process RSS includes code/runtime, intermediates, allocator retention
  and construction scratch. At 2^26 the medians were 531.2 MiB for succinct
  fixtures, 7976.0 MiB for BFS references and 8483.4 MiB for preorder references.
  Do not divide setup RSS by n and call it the final succinct bytes/node.

Old experiments sometimes counted only persistent cache, or also charged
shared program tables. Their 512 KiB labels are not interchangeable with the
final combined runtime budget. The current maximum reported cache+scratch in
the 1M/16M/67M check was 32,752 B / 524,272 B, within the requested limits.
Excluding constexpr tables from the byte allowance does not make their data
or the compiled handler code free of hardware cache/instruction-cache cost.

## Validation and experimental discipline

The runtime-budget change passed 473 Release tree CTests and a focused native
and Clang ASan/UBSan probe across twelve budgets, including values immediately
around 32 KiB. The probe covered random trees, a long chain and a wide star.
This is not a claim that the entire tree suite ran under every sanitizer/ISA.
Earlier individual SIMD experiments have their own validation scopes.

Keep these correctness gates for future kernel changes:

- Compare exact interleaved entry/return/inspection sequences, full node
  numbers and representation positions; totals alone are insufficient.
- Exhaust ordered small shapes where feasible: 197 through seven nodes,
  2056 through nine. Add random, chain, star and binary shapes, word/chunk
  boundaries, clean/dirty tail bits and deep navigation fallbacks.
- For low-level SIMD, test AVX-512, AVX2 and scalar paths, unaligned/bounded
  loads, output canaries and coordinates above 2^40. Avoid silent narrowing.
- For caches, test callback-triggered eviction, partial lists, copy/reset,
  exceptions and retained views. Check byte limits including fallback paths.

For performance, pin CPU, finish compilation/tests before timing, preserve
the observer, and compare same-binary/same-fixture candidates where possible.
Use alternating/rotating blocks to reduce shared drift. Bootstrap intervals
from these experiments were exploratory and did not model serial correlation.
An interval containing 1 does not prove equivalence. Confirm promising screens
in registered benchmarks and across the required workloads/sizes before
promotion. An unchanged streaming algorithm's timing drift is not a cache win.

PMU captures must exclude setup/warmup, verify the enabled phase and equal
logical work, reject multiplexed/not-counted events, and preserve matching
binary/source identity. Sampled cycle percentages are not source-line stopwatches;
fewer instructions or cache misses alone do not establish faster execution.

The published reproduction entry point is
[the tree benchmark page](../../../../../src/docs/benchmarks/trees.md). The
collector retains raw iterations separately from aggregates and validates
repetition counts, errors and byte limits. If raw experiment archives are
deleted, the historical rounded values in these notes remain evidence summaries,
not a substitute for rerunning a changed implementation.
