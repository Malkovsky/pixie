# Succinct tree traversal kernels and the bounded path index

> Historical experiment snapshot. Names, ownership, and candidate variants below
> describe the experiments, not the current API. See the maintained
> [experiment notes](README.md)
> and `include/pixie/tree/root_navigation_index.h` / `traversal_context.h`
> for the finalized implementation.

Implementation note, 20 October 2026.

This note describes Pixie's current LOUDS, balanced-parentheses (BP), and
depth-first unary degree sequence (DFUDS) traversal implementations. It records
the design, memory costs, and fallback behavior rather than a benchmark snapshot.

## Two kinds of retained state

| Mechanism         | Retained information                                     | Lifetime          | Purpose                                  |
|-------------------|----------------------------------------------------------|-------------------|------------------------------------------|
| Traversal kernels | Bit cursors, ancestors, optionally one decoded subtree    | One traversal     | Avoid navigation queries during DFS      |
| Bounded index     | Selected handles, degrees, and contiguous child lists     | Across many walks | Accelerate frequently used path steps    |

Both preserve the original representation and produce its normal node handles.
The kernels live under `include/pixie/tree/detail/`. The DFS workload and dispatch
remain in the benchmark harness, not in `TreeBase`. The bounded path index is an
explicit experimental auxiliary structure, not an implicit addition to every
tree.

Earlier navigation improvements remain useful on fallback paths: navigation
views cache query results, formulas avoid redundant rank/select or excess
queries, and a bounded ancestor path view reuses parent navigation state.
`for_each_child(node, visitor)` lets each implementation enumerate children
without allocating a result container.

## Traversal semantics

The DFS workload observes each node on entry and each internal node again after
its final child finishes. Leaves have no separate return event. For a root with
two leaves, the sequence is:

```text
enter root
enter child 0
enter child 1
return root
```

The inspection variant enumerates all immediate children on entry, before
descending. This ordering is important: a sequential encoding scan can emit
ordinary DFS events cheaply, but child inspection must identify sibling roots
across intervening subtrees.

BP and DFUDS consequently use different kernels for plain and inspected DFS.
All callbacks still receive complete handles, in the prescribed order.

## LOUDS: advancing cursors for each level

Pixie encodes a BFS-ordered node as a one followed by one zero per child. Its
handle stores a zero-based BFS number and the position of that one.

Although DFS jumps between depths, nodes encountered at any one depth appear in
LOUDS order. Each depth can therefore retain an independently advancing cursor,
containing its position, word index, and remaining delimiter mask.

Finding another delimiter within a word uses:

```cpp
position = word_base + countr_zero(mask);
mask &= mask - 1;
```

The kernel loads another word when the mask becomes empty. A known next
delimiter yields:

```text
degree = next_delimiter - current_position - 1
first_child_number = current_position + 1 - current_node_number
```

The first child position initially comes from normal navigation. Later positions
come from advancing the cursor for the next depth.

For inspection, the kernel copies the next level's cursor, enumerates the
children using the copy, and then descends with the original cursor. Inspection
adds another delimiter scan without a rank/select operation per child.

Child boundaries are computed once and passed into descent. Leaves are emitted
directly without recursive calls. There are 129 cursor slots, occupying roughly
3 KiB on the measured 64-bit platform, plus bounded recursive state. At depth
128, ordinary navigation handles the remaining subtree.

A separate LOUDS child-enumeration primitive works in windows of at most 256
source bits. Each window begins at an indexed child position and decodes
consecutive node delimiters. Long unary runs therefore do not force an
unbounded scan to reach a requested starting child.

Source: [unpack_louds.h](../../detail/unpack_louds.h).

## BP plain DFS: precomputed byte programs

BP stores the tree in DFS order. An opening one introduces a node; closing zeros
describe completed subtrees.

The kernel reads eight bits plus one lookahead bit. Those nine bits select one
of 512 precomputed programs describing:

- Visitor events and their order.
- Events referring to nodes opened within this byte.
- Events referring to existing ancestors.
- Changes in node numbering and depth.
- Newly opened nodes that must remain on the ancestor stack.

The lookahead distinguishes a close followed by another close, helping identify
the internal-node return events used by the workload.

The AVX-512 implementation builds event numbers and positions in vectors and
uses permutations to choose between local nodes and saved ancestors. Callback
invocations remain individual and ordered.

The universal table depends only on bit patterns, not on a particular tree. It
occupies 16 KiB shared across trees. Per-traversal ancestor storage is roughly
2 KiB. Tail fragments, special patterns, and deep paths use scalar handling.
Beyond the cached ancestor depth, normal parent navigation provides the
fallback.

Source: [unpack_bp.h](../../detail/unpack_bp.h).

## DFUDS plain DFS: degree streaming and completion excess

Pixie's DFUDS stream stores each preorder node as `degree` one bits followed by
a zero. The kernel extracts zero positions from each 64-bit source word:

```text
degree = zero_position - record_start
next_record_start = zero_position + 1
```

Node numbers increase sequentially. The running excess changes by:

```text
excess += degree - 1
```

An internal node entered at excess `e` finishes when excess reaches `e - 1`.
The kernel saves this target alongside the node handle. After consuming a node,
matching targets pop ancestors and emit return events. One leaf can finish
several ancestors.

This avoids parent and last-child queries on the normal path. The stack holds
128 cached ancestor entries, with a sentinel slot. Deeper traversal uses
navigation until it returns to the cached portion.

Source: [unpack_dfuds.h](../../detail/unpack_dfuds.h).

## BP and DFUDS inspected DFS: complete-subtree unpacking

The inspected traversal probes up to 2,048 source bits, accepting at most 1,024
nodes. It searches for the subtree endpoint using relative excess:

- BP: relative excess returns to zero.
- DFUDS: relative excess reaches minus one.

Only a complete subtree is accepted. The current implementation does not carry
a partially decoded subtree across chunk boundaries.

### Decoded representation

A successful probe produces preorder records:

```cpp
struct Record {
  size_t number;
  size_t position;
  uint16_t parent;
  uint16_t next_sibling;
  uint16_t degree;
};
```

Records occupy 24 bytes on the measured platform. The record buffer is roughly
24 KiB, including an additional dummy record; construction uses further
temporary arrays. A single buffer is reused throughout the traversal, rather
than allocated in each recursive frame. This expansion applies only to the
bounded scratch subtree, not to the whole tree.

A nonleaf's first child is implicit: its local preorder index is the parent's
index plus one. Sibling links skip entire child subtrees. A dummy record absorbs
the first sibling-link write, avoiding a special first-child branch.

### BP decoding

For an opening at local bit position `p`, with `i` preceding openings:

```text
depth = 2*i - p
```

An ancestor array indexed by depth supplies the parent. The decoder links
siblings and increments the parent's degree as each child is encountered.

### DFUDS decoding

For local node index `i` and record start `p`:

```text
pending = p - 2*i
```

This identifies the pending-child slot owned by the node's parent. A degree-`d`
node fills slots `[pending, pending + d)` with its own local index. Across a
complete subtree, these fills cover exactly `n - 1` child slots.

The implementation uses AVX-512 masked stores where supported and a scalar
fallback otherwise. It additionally writes the immediately consumed slot
`pending + d - 1` with a scalar store for nonleaves. For leaves, an extra scalar
write touches a discarded but valid scratch slot. The additional scalar store
improved measured performance; better handling of the following dependent load
is a working microarchitectural explanation, not a separately proven result.

### Traversing the decoded subtree

The kernel scans records in preorder:

1. Emit the current node's entry.
2. Inspect its child list through sibling links.
3. Advance to the next preorder record.
4. At leaves, follow completed ancestor links to emit return events.

Within the decoded subtree, these operations need no rank/select or RmM query.
Decoding and traversal are kept in separate non-inlined kernels to avoid
expanding the caller into a large function with excessive register spills.

If a subtree does not fit, normal navigation enumerates its children. Up to 32
handles are retained across descent. Larger child lists are enumerated again
for descent, rather than requiring an unbounded child buffer. At recursion
depth 128, a navigation fallback handles the subtree. The recursive child
buffers add bounded stack storage beyond the shared decoded-record buffer.

Sources: [unpack_cache.h](../../detail/unpack_cache.h),
[unpack_bits.h](../../detail/unpack_bits.h).

## SIMD scope

The low-level primitives own feature dispatch:

- BP byte programs use AVX-512 permutations to assemble event handles.
- Delimiter extraction uses AVX-512 VBMI2 byte compression when available,
  otherwise scalar trailing-zero extraction.
- Pending-parent filling uses AVX-512BW stores when available, otherwise a
  scalar loop.

These primitives do not currently have dedicated AVX2 implementations. A
platform lacking the required AVX-512 feature uses the corresponding scalar
fallback. SIMD accelerates decoding; it does not batch or omit callbacks.

## Bounded path index

`TreeNavigationIndex<Tree, Capacity>` is an explicit experimental structure
borrowing an immutable tree. Its default capacity is 32,768 entries:

```cpp
struct Entry {
  TreeNode node;
  size_t degree;
  uint32_t first_child_entry;  // Or the missing-children sentinel.
};
```

Entries occupy 32 bytes on the measured platform. Default reserved storage plus
the index object is 1,048,608 bytes, approximately 1 MiB or 0.105 bytes per node
at ten million nodes. This excludes allocator overhead and temporary build
state. The 32-bit index refers only to this bounded entry array; original node
handles retain full-width fields.

### Choosing which nodes to expand

Construction starts at the root and maintains a priority queue. Under uniform
child choices:

```text
P(root) = 1
P(child) = P(parent) / degree(parent)
```

The pending node with greatest reach probability is expanded first:

1. Check that its complete child list fits in the remaining capacity.
2. Append its children contiguously and cache their degrees.
3. Enqueue nonleaf children with their reach probabilities.

If the list does not fit, the node remains unexpanded; construction can continue
with other pending nodes. This uses topology only, with no sampled-query
training. Likely single-child chains may be expanded deeply while less likely
branches remain shallow.

The heuristic prioritizes reach probability, not a proven optimal benefit per
byte. Whole-child-list admission also matters: a root with degree greater than
the available budget cannot be expanded. Reserved bytes can exceed the bytes
actually populated with entries.

### Using a cursor

Each walk obtains a fresh cursor initially at the root. Within expanded nodes:

```text
degree()         -> read cached degree
move_to_child(i)  -> entry = first_child_entry + i
inspect children -> scan contiguous child entries
```

No succinct navigation query is needed for these steps. The combined
`inspect_children_and_move_to_child(i, visitor)` scans all child entries and
then chooses `first + i`, without a selection branch for every inspected child.
On fallback, it remembers the selected handle while enumerating children, so
it does not fetch that child again.

After leaving the indexed region, the cursor uses the original tree for the
rest of its downward walk. The index has no eviction, online population, or
per-walk mutation. Independent cursors share the same immutable index. The
borrowed tree and its backing data must outlive the index; the index must remain
alive and unmoved while cursors use it.

The index is constructed before benchmark timing and reused across walks.
Reported path timings therefore measure navigation, not construction. Its
benefit is intended for repeated walks over the same immutable tree.

Source:
[tree_navigation_index.h](../../root_navigation_index.h).

## Evaluation conventions

The four workloads are full DFS, DFS with child inspection, random root-to-leaf
paths, and those same paths with child inspection. Randomness comes from a
64 KiB shifted query pool. Each nonleaf consumes one choice, uniformly selecting
a child rather than uniformly selecting a leaf.

The optimized comparison uses unpacked DFS and indexed paths. These are
specific traversal strategies, not a claim that every uncached navigation
operation beats a reference implementation.

Correctness checks compare complete event and child sequences, handle
positions, counts, and identical path choices. Tests include short exhaustive
trees, boundary sizes, dirty unused tail bits, deep fallback paths, and
multiple shapes and backends. Timings use optimized native builds pinned to
one CPU, with construction outside the measured region and raw results kept
outside the repository.

## Additional pointer-tree reference

The comparison suite also provides
[`pixie::references::PointerTree`](../../../../references/pointer_tree.h).
Each node stores exactly three native pointers: first child, next sibling, and
parent. Records live in one contiguous allocation, in either preorder or BFS
order, to match the relevant succinct encoding's numbering.

There is no cached degree or node number. On a 64-bit platform, each record is
24 bytes. Parent, first-child, and next-sibling navigation take constant time.
Counting children scans the sibling chain, and selecting child `i` follows
`i` next-sibling links after accessing the first child. Child inspection walks
the same chain once.

This is a stronger DFS reference than the existing adjacency-list wrapper in
one respect: the adjacency wrapper searches the parent's child list to obtain
the next sibling, whereas the pointer tree has that link explicitly. Conversely,
the adjacency wrapper has constant-time degree and indexed-child access, while
the three-pointer representation does not. The four workloads expose both
trade-offs without adding a fourth field for cached degree.

The reference owns its records and does not retain its construction adjacency
list. It uses pointer handles rather than `TreeNode` and does not inherit
`TreeBase`. Benchmark rows are named `pointer_preorder_*` for BP/DFUDS and
`pointer_bfs_*` for LOUDS. The `reference_bytes` and
`reference_bytes_per_node` counters include the owner and reserved records,
excluding allocator overhead. Construction is outside timing, as for the
other implementations.
