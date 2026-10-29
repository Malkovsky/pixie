# Traversal decoding and SIMD experiments

## What unpacking means here

Read a bounded region of the succinct bitstring and derive temporary navigation
or event information that can be traversed without repeated rank/select/RmM
queries. It is not necessarily a hardware prefetch, a persistent node cache,
or a conversion of the whole tree to pointers. The successful large-tree
implementations discard/reuse decoded scratch as the traversal progresses.

Plain DFS and DFS with child inspection are materially different workloads.
Preorder entry events can follow a stream, but observing all immediate children
before descent requires knowing child starts while skipping their subtrees.
Neither a faster plain-DFS stream nor a faster topology-only scanner proves
that inspected DFS improved.

## Plain DFS kernels

- **BP:** process byte patterns using precomputed event descriptions, derive
  node handles and maintain a bounded ancestor buffer. Handle deep regions
  and boundaries through exact fallback paths. A closing parenthesis does not
  itself always correspond to an emitted return event: emit the same nonleaf
  returns as ordinary DFS.
- **DFUDS:** enumerate zero delimiters in word-sized batches. The distance
  from the previous record gives degree. Update excess by `degree-1`; an
  internal node finishes when excess reaches its entry excess minus one.
  A bounded stack of completion targets emits returns; deeper cases retain
  navigation fallback. The hot loop no longer performs an RmM query per node.
- **LOUDS:** use monotone cursors for levels, because each node's children are
  consecutive in BFS order. A level cursor can retain source position/state
  while DFS revisits that depth, avoiding repeated global select operations.

A query-only 10M-node DFUDS plain-DFS profile measured about 45.3 ms unsampled
and 47.1 ms sampled. Almost all samples were in the excess-stream kernel;
completion-target comparisons and ancestor-stack pops received about 45% of
cycle weight. Rank/select was not a meaningful sampled hotspot. This is
evidence to inspect control/stack work, not a precise latency attribution to
those source lines. Evidence: the historical `dfuds-dfs-report` snapshot.

## Compiled microprograms inside bounded chunks

Two scales must not be confused:

1. An outer accepted complete subtree/chunk may occupy thousands of bits
   (8192 bits in the selected 512 KiB configuration; 512 at 32 KiB).
2. A compiled microprogram represents one of the **197 ordered rooted shapes
   with one through seven nodes**. Many such programs execute inside one
   outer chunk, with explicit macro records connecting them.

The loader probes whether a complete subtree fits in the bounded bit window.
On success it decodes a compact skeleton of macro nodes and recognized small
shapes. On failure it navigates downward and tries again. This makes use
common near the bottom of a large tree, but there is no fixed level at which
microtrees become legal: any fitting subtree can be decoded, including a root.
Trailing padding and word boundaries must not create a false complete subtree.

A small topology lookup maps the local bits to a program ID/length. Shared
constexpr metadata describes relative node-number/bit-position offsets,
inspection events and return flags. A separately compiled handler per shape
expands those offsets against the full-width root handle and invokes the
original callbacks. The topology program is reusable across trees; its
instance coordinates are transient. The compiler sees constant offsets and
event order, removing a general replay interpreter from this small unit.
BP and DFUDS have encoding-specific programs, not interchangeable bit keys.

Retained improvements, each with its own historical control:

| Change | Evidence and interpretation |
| --- | --- |
| Buffered bit reader/reservoir | About 4–5% at 2048-bit BP chunks, little 8192-bit gain; reduces refill/control work rather than instruction count alone. |
| Absolute handles in decoded records | Roughly 2% in a BP screen despite increasing record storage; avoids repeated coordinate reconstruction. |
| Compact 24-byte sibling-link records | BP registered 84.66→81.33 ms at 10M; preceding paired result was only about 1.4%. DFUDS registered gain about 2%. Saves about 32 KiB per 8192-bit decoder versus 32-byte records. |
| Initialize last-child state only for macro records | DFUDS paired gain about 0.4–0.5% at 10M; larger sequential registered deltas were not attributed entirely to this store removal. |
| Vector root argument for BP handlers | Registered compact control 78.86→76.26 ms at 10M, about 3.3%; preceding paired screen about 2.6%. Retained in the selected BP path. |

The vector-root optimization is small and specific. Keep the two 64-bit root
coordinates in a 128-bit integer vector across the internal handler call;
add compile-time offsets in two lanes, then produce the unchanged `TreeNode`
callback arguments. It avoids reconstruction/reloads at the call boundary.
It is not a new 128-node decoder or a wider tree format. ISA guards and scalar
fallback remain in the low-level helper. Do not add the percentages above:
controls, code layouts, accounting rules and measurement dates differ.

Historical evidence labels: `reservoir-report`, `absolute-report`,
`bp-compact-registered`, `dfuds-record24-retained-report`,
`dfuds-macrolast-report`, `bp-vector-registered`, 2026-10-24 through 2026-10-26.
Their old persistent-only cache counts are superseded by the combined runtime
accounting described in the navigation notes.

## LOUDS inspected traversal

Inspection and descent share bounded per-level buffers of child-end positions.
An explicit stack keeps traversal frames compact. Decoded delimiters let the
same handles be reconstructed without selecting them again. Deep or wide
branches have exact navigation fallbacks. This is an encoding-specific stream
and child-list optimization, not BP-style complete-subtree microprogramming.

The decoded child-position approach improved a historical 10M inspected-DFS
run from about 115 to 87 ms, while BFS pointers were about 68 ms. Attempts to
retain full child handles, update frames in place, add an all-leaf-group
shortcut or batch node-count accounting did not improve the retained path.
One profile put decoder refill at only about 4% of local sampled cycles;
optimizing refill alone could not explain or close the whole remaining gap.

## Why fitting in L3 was insufficient

Query-only inspected-DFS counters at 10M, per logical node, are useful evidence:

| Historical kernel | Instructions | Cycles | Branch misses | Generic cache misses |
| --- | ---: | ---: | ---: | ---: |
| DFUDS absolute records | 65.35 | 40.45 | 0.850 | 0.00978 |
| BFS pointers, same DFUDS comparison | 24.25 | 30.34 | 0.976 | 0.01191 |
| LOUDS iterative decoder | 82.63 | — | 0.962 | 0.00237 |
| BFS pointers, same LOUDS comparison | 24.25 | — | 0.997 | 0.01609 |

Succinct traversal performed much more work despite fewer measured cache
misses. These generic PMU events are **not** measured L3 miss counts. Separate
groups were 100% scheduled, but these are fixed-work captures, not statistical
counter distributions or additive latency breakdowns.

For compiled BP/DFUDS inspected traversal, another sampled capture attributed
about 37%/44% of cycles to validation/decoding and 43%/38% to microprograms,
with roughly 16%/15% in replay/wrapper code. The microprogram region includes
callbacks and function overhead. Both decode and replay deserve attention.
Evidence: `dfuds-cost-counters-report`, `louds-iterative-profile-report`,
`large-micro-profile-report`, 2026-10-24/25.

## Negative results worth remembering

- **Shared relative-event tables:** wide and compact byte tables lost to
  compiled handlers for both BP and DFUDS. At 8192 bits, one BP comparison was
  82.07→89.93 ms and DFUDS 94.70→101.38 ms. Smaller metadata did not remove
  ordered per-event loop/reconstruction costs.
- **Cache full handles inside each microprogram:** local arrays/packed pairs
  added preparation, spills and reloads; no repeatable improvement. Reuse of
  arithmetic does not guarantee reuse in registers.
- **Single-pass BP fit validation plus decode:** removing a separate fit pass
  regressed a screen from about 84 to 111 ms. Speculative decode of failed
  candidates can cost more than a fast preliminary check.
- **Subtree-end records instead of sibling links:** 24/32-byte variants lost
  by about 4–5% in paired BP tests. They retired fewer instructions but took
  more cycles. The compact sibling-link control was the useful candidate.
- **Larger chunks:** 16K/32K-bit BP chunks did not establish a paired advantage
  over 8192 bits and enlarged scratch. A lower noisy raw median was not enough.
- **Composite microtrees built from small child programs:** more coverage did
  not pay for admission probes, extra record types and replay dispatch.
- **Prepared complete-root event streams:** genuine wins at 1K/4K, but their
  storage grows with the tree. They do not solve the >=1M-node target under a
  fixed cache budget and are not the current large-tree solution.

The composite census is a useful caution for future SIMD work: 488,263 eligible
roots covered 5,151,943 nodes (51.52%) of the seed-42 10M-node BP tree, yet the
measured composite implementations regressed by roughly 8–11%. There were
only 32 eligible pure stars, covering 259 nodes. A fast isolated star kernel
can therefore have negligible opportunity on this particular random input.
Admission cost and structural frequency both need measurement.

## SIMD experiments through 128 nodes / 256 bits

Enumerating seven, eight and nine-node programs requires 197, 626 and 2056
shapes. SIMD expansion of compact/wide event descriptors into temporary handle
arrays was tested with AVX-512, AVX2 and generic C++ controls. It could improve
on generic replay yet still lose to compiled seven-node programs. One compact
eight-node repeated comparison was about 1.08x slower. This separates a
primitive-level SIMD win from an end-to-end algorithm win.

The arbitrary 128-node experiments avoided enumerating all shapes:

1. Copy up to four words and test whether the matching close fits in 256 bits.
2. Expand bits to +1/-1 byte lanes. Shift/add stages at 1, 2, 4 and 8 positions
   compute prefix excess inside 128-bit lanes; popcounts supply cross-lane and
   cross-word carries. Depth 128 must remain distinct under byte arithmetic.
3. Match each opening at excess `d` with the first later `d-1`; variants
   materialize links, replay directly from masks, fuse child inspection with
   matching, or retain temporary depth masks. Relative preorder number is
   `(opening_position + inclusive_excess - 1)/2` for local BP coordinates.
4. A level-wise variant forms opening/closing masks for each depth: the kth
   opening matches the kth closing at that level. It extracts positions with
   trailing-zero loops or AVX-512 VBMI2 byte compression. A height threshold
   falls back to scalar stack matching.

Representative final 10M-node BP inspected-DFS comparison:

| Kernel | CPU ms | Instructions/node | Cycles/node |
| --- | ---: | ---: | ---: |
| Compiled seven-node programs, 8192-bit outer chunks | 80.52 | 66.33 | 36.16 |
| Same programs, 256-bit chunks | 101.14 | 87.94 | 45.14 |
| Level masks, scalar extraction | 175.97 | 133.61 | 77.13 |
| Byte compression directly to memory | 393.07 | 118.09 | 176.96 |
| Byte compression to registers, then ordinary stores | 149.69 | 118.81 | 65.71 |

Register compression used padded temporary arrays for full-width stores and
was much better than the memory-destination compression instruction on this
host. It still lost to the retained integrated decoder. A separate per-opening
matching/direct/fused replay screen was around 161–166 ms versus 82.5 ms;
temporary depth-mask caching was still about twice the compiled control.

The interesting exception was an isolated complete 128-node **star**:
register compression measured about 544 ns/tree versus 891 ns for compiled
programs, approximately 1.64x faster. Random complete 128-node subtrees still
lost (about 1199 versus 567 ns). This supports a possible selective wide/shallow
specialization, not replacing every random subtree. Keep the large outer
chunks when testing such a specialization, so small-chunk fallback overhead
does not dominate the comparison.

These were BP-only experiments. No DFUDS/LOUDS gain follows from them. The
unfinished eight-node vector-program prototype was explicitly left unvalidated
when the optimization iteration closed; do not promote it during cleanup.
Evidence: `microtree-simd`, `microtree256-simd`, `microtree-level-simd`,
`bp-subtree-end`, `bp-chunk-capacity`, `bp-composite-micro`,
`compiled-handle-cache`, and `bp-vector-registered`, 2026-10-26/27.
