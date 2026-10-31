# Navigation and caching

## Node handles and local query reuse

A handle carries a logical node number and a representation bit position.
These are different coordinates. The extra position avoids rediscovering the
encoding location from the number at every operation; do not silently narrow
either coordinate to 32 bits in the library. Benchmark generators can have a
narrower, explicitly documented domain.

A local navigation view borrows the immutable tree and retains results for
one current node. Reset invalidates the memoized state; copying a local view
copies independent state. It is not a root-to-node path stack, a shared graph
cache, or an owner of the source bitvector.

Useful retained quantities:

- BP: matching close, degree, parent, sibling and last-child result.
- DFUDS: unary-record delimiter, subtree end, incoming opening, parent,
  sibling and last-child result.
- LOUDS: next node delimiter, incoming zero-edge position, first child and parent.

Formulas reduce work before caching is considered. With preorder number `v`
and start bit `p`, a BP subtree ending at close `q` contains `(q-p+1)/2`
nodes; the next sibling starts at `q+1`. A DFUDS subtree has `2k-1` bits,
so a next-sibling start `s` advances the node number by `(s-p+1)/2`.
DFUDS degree and first child share the select0 delimiter `d`: degree `d-p`,
first child `(v+1,d+1)`. These statements use Pixie's encoding conventions,
not interchangeable formulas for arbitrary BP/DFUDS definitions.

LOUDS stores one delimiter followed by degree zeros. Degree uses the next
delimiter distance; sibling numbers are consecutive. The adapter searches the
remainder of a machine word for a nearby delimiter and falls back to indexed
select when absent. This is bounded local scanning, not an unbounded linear
search replacing select. Parent and last-child operations can reuse the same
incoming zero-edge position. BP/DFUDS node numbers and bit distances similarly
remove explicit rank calls from the navigation adapters; their RmM operations
can still depend on the backend rank/select support.

One concrete remaining duplication was DFUDS `degree()` followed by child
inspection: enumeration created another local view and repeated select0.
Enumeration through the existing view fixed this with no extra field. Paired
ordinary inspected-path screens at 1M/10M found approximately 3% improvement.
Registered follow-ups supported modest path gains, but ordinary inspected DFS
was effectively unchanged. The optimized streaming/microtree kernels mostly
bypass these views, so the result is not a 3% improvement to every workload.
Evidence label: `node-view-query-reuse-{,integrated-}report`, 2026-10-25.

Caching a single rank/select result also cannot remove a whole repeated
child-list walk. Prefer an enumeration interface that allows implementation-
specific reuse. When inspection already obtains the chosen child handle,
pass that handle into descent instead of performing `child(i)` again.

## Bounded root-connected index

This is a prebuilt partial graph, not an LRU cache and not a current-path stack.
On the measured 64-bit ABI an entry occupies 32 bytes: a full `TreeNode`,
degree, and a 32-bit first-child entry index plus padding. Child entries are
contiguous. A sentinel marks a child list that was not expanded.

Construction starts at the root with probability 1. A max-priority queue
expands the most likely indexed node; each child receives probability
`parent_probability / parent_degree`. This models uniform child selection,
uses topology only and does not train on the measured random query pool.
Ties are resolved deterministically by index. A node is expanded only if its
**entire immediate child list** fits. If it does not fit, that node is skipped
and other pending nodes can still be expanded. The index may therefore leave
some reserved capacity unused; a high-degree root can prevent useful expansion.

A cursor uses direct entry indexing for cached degree, children and descent.
Inspecting a cached list is one contiguous loop. Once it descends beyond the
stored frontier, it carries an ordinary tree handle and uses representation
navigation for the rest of that root-to-leaf path. There is no address hash
lookup to rediscover unrelated cached descendants. The next path starts at
the indexed root again. The current public cursor is downward-oriented.

When capacity is reached, construction stops expanding useful pending nodes.
There is **no eviction, hit promotion, runtime refill or copying of the whole
tree**. Reading cached children is ordinary data access, not necessarily a CPU
prefetch instruction. The runtime byte limit caps reserved entries; construction's
priority queue is temporary and outside the traversal-cache allowance.

This design made repeated paths much faster than general CLOCK policies on the
measured distribution. It is distribution-sensitive: topology-based uniform
child probabilities may be a poor predictor for biased application queries.

## General cache experiments

CLOCK entries retain local navigation results, known degree and a child prefix.
A hit sets a recent bit. A persistent clock hand clears recent bits and evicts
the first unmarked entry. A fixed linear-probing directory uses backward-shift
deletion. This is approximate recency, not exact ordered LRU promotion.

Unconditional shared lookup/admission hurt one-pass DFS badly. Operation-based
admission avoided this: parent/sibling/leaf/first-child operations stayed local;
degree and indexed-child requests activated shared caching. A full historical
512 KiB acceptance sweep had a worst CLOCK/local DFS ratio of 1.0087, but this
mostly demonstrated that those DFS loops avoided the shared cache. It did not
make per-node hashing free. Evidence: `clock512`, 2026-10-24.

The eight-child prefix was a real weakness: the seed-42 10M-node root has 19
children. Increasing the prefix to 32 and snapshotting complete child lists
before callbacks traded entry count for useful child coverage. Historical
registered BP paths improved roughly 1342→331 ns and inspected paths
2795→516 ns; DFUDS 1195→238 ns and 2551→376 ns. These remained slower than
the bounded root index. LOUDS results were noisier and did not establish a
general win. Evidence: `cache-prefix`, 2026-10-25.

LOUDS is not an equal-cost cache control for BP/DFUDS. Its indexed children
and delimiter scans are cheaper; BP/DFUDS prefix misses can walk preceding
sibling subtrees. In retained-unpacking experiments, LOUDS also kept its
specialized level-cursor traversal and largely bypassed shared retention.

Two combinations were tried: decode temporary chunks and retain requested node
results selectively, or retain/evict whole decoded chunks. At 10M, early
inspected DFS measured BP about 197/204 ms versus standalone unpacking 105 ms;
DFUDS about 196 ms versus 113 ms. Initial chunk caps differed (256 versus
2048 bits), so that comparison did not isolate retention alone. Matching-cap
controls also lost (BP roughly 187/197 ms retained versus 132 ms discard).
The conclusion is about these implementations, not an impossibility theorem
for cache/unpacking combinations. Evidence: `unpack-cache-comparison` and
matching-size controls in `cache-prefetch-goal`, 2026-10-24.

### Callback and lifetime constraints

Never hold an evictable cache-entry reference across a user callback. Copy
handles or bounded snapshots first; nested queries can evict the entry.
Views must not pin slots or become dangling when entries move. Tests need
partial lists, nested eviction, exceptions and copy/reset behavior. Account
for callback-safe temporary copies in any combined budget. External
synchronization is required for shared mutable cache use; a const tree API
alone does not imply thread safety.

## Runtime budget adaptation

The selected adapter borrows a tree using local navigation and the matching
source words. The tree constructor carries a runtime byte limit (default
512 KiB). The adapter divides that allowance between its context, the root
index and explicit temporary traversal buffers. It does not add another full
cache budget on top of a shared tree policy.

| Encoding | Scratch at 32 KiB | Inspected decoder at 32 KiB | Scratch at 512 KiB | Inspected decoder at 512 KiB |
| --- | ---: | --- | ---: | --- |
| BP | 16,064 B | 512-bit chunks, cached depth 15 | 182,400 B | 8192-bit chunks, depth 128 |
| DFUDS | 16,064 B | 512-bit chunks, cached depth 14 | 190,592 B | 8192-bit chunks, depth 128 |
| LOUDS | 15,232 B | 32 levels, 16 children/level | 73,088 B | 128 levels, 32 children/level |

Up to half the available bytes go to scratch, capped at 192 KiB for BP/DFUDS
and 80 KiB for LOUDS. Plain streaming needs only 4 KiB. Parenthesis candidates
are 8192/2048/512/128 bits with runtime-derived cached depth; LOUDS uses a small
set of compiled depth/child-capacity choices. Dispatch occurs once per walk,
not per node. Separate noinline helper frames prevent a large-profile array
from living in the small-budget path's frame. Zero/tiny budgets remain correct
through ordinary navigation.

The first 32 KiB implementation disabled decoding because it reserved the
largest scratch profile unconditionally. That artificial cliff caused most
of its DFS slowdown. A second issue was falling back to ordinary traversal
for an entire deep subtree once cached depth was exhausted. The final
BP/DFUDS fallback navigates without another ancestor child array but continues
trying the same reusable microtree decoder at every depth.

At 2^26 nodes, the final 32 KiB versus decoder-disabled 32 KiB speedups were
8.73/4.06x for BP plain/inspected DFS, 8.64/6.00x for DFUDS and 4.72/5.44x
for LOUDS. However, reserving decoder scratch reduced root-index capacity:
paths became 1.13–1.63x slower across the 1M/16M/67M comparison sizes.
This is a real tradeoff, not a uniformly faster policy. The early 32 KiB
fallback data must not be mixed into the final runtime-adaptive curves.
