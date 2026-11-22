# Tree navigation and traversal contexts

BP, DFUDS, LOUDS and ForestCover expose the same `pixie::tree::TreeBase` navigation
contract. Use `Tree::Node` for a handle and `tree.navigate(node)` for a
positioned view. A handle belongs to its tree: BP, DFUDS and LOUDS use a
node number and bit position; ForestCover uses a forest unit and a local
breadth-first number. `tree.number(node)` and `tree.find_by_number(id)`
convert between handles and the representation's stable numbers.

`auto path = tree.tree_path(u, v)` returns the nodes in order from `u` to `v`,
including both endpoints. Equal endpoints give one node. For repeated queries,
pass two distinct vectors of handles to `tree.tree_path(u, v, out, scratch)`.
The call replaces `out` and leaves `scratch` empty. Both keep their capacity.

Each query takes one parent step per path edge and stops where the endpoints
meet. BP and DFUDS compare preorder numbers, LOUDS compares breadth-first
numbers, and ForestCover compares `(unit, local)` coordinates. In each order,
a parent precedes its children. The adjacency-list and pointer references use
input IDs when parents precede children; otherwise they build a preorder index
once at construction. The original handles and numbering stay unchanged.

`ForestCoverTree` stores disjoint forests of up to 256 nodes. Its IDs are in
unit order, then breadth-first within each unit. It does not expose a single
BP/DFUDS/LOUDS encoding, a borrowing whole-tree view, or serialization. Direct
navigation needs no persistent cache. For repeated walks, a forest traversal
context owns decoded slots and optionally borrows a root index:

```cpp
#include <pixie/tree/forest_cover.h>

pixie::tree::ForestCoverTree forest(source);
pixie::tree::RootNavigationIndex<pixie::tree::ForestCoverTree> index(forest, 16 * 1024);
using Context = pixie::tree::TreeTraversalContext<pixie::tree::ForestCoverTree>;
Context context(forest, 4 * Context::slot_bytes(), &index);
context.traverse_inspected(visit, inspect_child);
context.inspect_path(visit_node, inspect_child, choose_child);
```

Keep the index and tree alive and unmoved while the context uses them. The
second constructor argument limits **slot bytes**, not total query memory.
Count the context object, its allocated slots, the index, and the walk's
scratch separately. The index accelerates inspected paths; DFS uses one
walk-local decoded forest so a nested callback cannot invalidate it.

Use `index.inspect_part()` or `forest.inspect_unit(unit)` to inspect an owning
snapshot without keeping a cache slot live. Node and edge indices inside the
snapshot are local. The first `root_count()` nodes are all the local roots in
source order. Each node's `children` lists only its children in the part.
`children_complete` states whether the producer verified the full child list. An
empty list denotes a leaf only when it is complete. `<pixie/tree/event_part.h>` also
provides `inspect_event_part<Bits>(tree, root)` for BP and DFUDS.

BP, DFUDS and LOUDS own their encoded words in `AlignedStorage` by default. Construction
from a word span copies; construction from an `AlignedStorage` rvalue transfers
its allocation. Use the concrete header for the chosen representation.

`TreeView` on these three contiguous encodings is the same
representation with `ReadOnlyStorageView`: it borrows aligned complete encoding
words, but owns its constructed backend metadata. Its backing bytes must remain
alive, immutable and at a stable address. It can also be built directly from a
read-only storage view. This is distinct from a per-node `NodeView`.

Owning copies duplicate the encoding and rebuild backend references; view copies
share encoding and copy backend metadata. Moves transfer storage and indexes.
Assignment, moving or destruction invalidates existing node views, indexes and
traversal contexts borrowing that tree; using invalidated borrowers is undefined
behavior. Moved-from trees may only be destroyed or assigned.

For these encodings, `encoding_words()` exposes a read-only span of complete encoding words.
`encoding_owned_bytes()` reports allocation bytes (zero for a view), and
`memory_usage_bytes()` includes owned encoding and backend allocations.
Trees serialize their encoding, node count and nested navigation backend as one
versioned artifact. `RmMBTree` and compressed rank/select also serialize
independently; their artifacts exclude the encoding they borrow. The original
baseline rank/select serialization remains supported.

```cpp
pixie::VectorOutputSink output;
pixie::BinaryWriter writer(output);
tree.serialize(writer);
writer.finish();
pixie::BinaryReader reader(output.bytes());
auto restored = Tree::deserialize(reader);
```

`deserialize(reader, DeserializationValidation::kQuick)` restores saved metadata
without rebuilding it from the encoding. It checks framing, configuration,
dimensions and internal metadata consistency. It is intended for trusted
artifacts: the caller must provide a valid encoding and its corresponding
metadata. Successful quick loading does not guarantee safe queries on arbitrary
malformed input. It may scan metadata, and owning loaders still copy the encoding.

`kFull` additionally verifies topology and source-derived metadata; this scans
the encoding and can allocate temporary metadata for validation. For example,
BP bits `1010` can have consistent navigation metadata but describe two roots,
so full tree validation rejects them. Neither mode computes a hash or establishes
artifact identity or authenticity: full validation checks consistency with the
supplied contents. Deserializers advance the input cursor only on success and
consume exactly one artifact.

Owning trees copy the encoding. `Tree::TreeView::deserialize` borrows the
artifact's word-aligned encoding bytes; the backing storage must stay alive,
immutable and unmoved. Compressed rank/select and RmMBTree metadata are restored
into owning allocations. A TreeView is therefore not a wholly zero-copy index.
Artifact offsets and backing addresses must keep embedded encoding words
aligned to eight bytes for view loading. File placement and bundling belong to
the application; multiple artifacts can be written to the same writer.

Owning trees also accept an ordered tree as a converting-constructor argument,
for example `pixie::tree::DFUDSTree<pixie::RmMBTree<>> dfuds(bp)`. The source provides `size()`,
`root()`, `degree(node)`, `first_child(node)` and `next_sibling(node)`. It must
describe a valid ordered tree, with stable handles during construction; sibling
navigation is called only when a next sibling exists. No CRTP inheritance is
required. Child order is preserved, but node numbers and bit positions belong
to the destination representation. The result owns its encoding and rebuilt
navigation metadata; it retains no source handles or caches. TreeView cannot
be a conversion destination. Same-type construction remains ordinary copying.

Conversion shares an internal DFS stream with entry, degree and exit events
for every node, including leaves. BP and DFUDS sources use their microtree
decoders; LOUDS sources use decoded child traversal. Generic sources use
navigation with an ancestor stack. The traversal producer is independent of
the destination encoding, so converters are not implemented pair by pair.

BP and DFUDS write directly into destination storage. LOUDS stages packed
degree records by depth and concatenates them: DFS restricted to one depth has
the same sibling order as BFS. This produces BFS-ordered output without a queue
of node handles. It requires O(nodes) packed staging bits plus O(height) level
records and vector capacity overhead. There is no per-node degree array.

Construction uses `size_t` counts. Its temporary memory is independent of the
query traversal-cache budget: decoded producers have fixed explicit scratch
allowances, generic traversal needs O(height) handles, and LOUDS additionally
needs the staging storage above. A root navigation cache is not built for
conversion. The constructor from encoding words remains available and avoids
traversal when the encoding already exists.

```cpp
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

// words contains a valid LOUDS encoding of node_count nodes.
pixie::tree::LoudsTree<> tree(words, node_count);  // Owns a copy.
using Tree = decltype(tree);
Tree::TreeView representation(tree.encoding_words(), node_count);
Tree::PathView path(tree);
auto view = tree.navigate(tree.root());
const auto degree = view.degree();

pixie::tree::RootNavigationIndex index(tree, 32 * 1024);
auto cursor = index.cursor();
if (cursor.degree() != 0) cursor.move_to_child(0);
```

A navigation view caches results for its current node. The root index instead
stores complete child lists for a bounded, root-connected part of the tree.
Construction prioritizes nodes by the probability of reaching them when each
child is equally likely. It does not learn from benchmark queries. There is no
eviction: navigation outside the indexed region falls back to the tree.
Cursors are independent and borrow the index; neither the tree nor index may
move or be destroyed while a cursor uses them.

`RootNavigationIndex::serialize(writer)` writes a separate artifact containing
only immutable entries. Restore it with
`RootNavigationIndex<Tree>::deserialize(reader, tree, byte_limit, validation)`.
The caller must provide the corresponding tree; quick validation checks entry
bounds and connected child ranges but cannot identify a different same-sized
tree. Full validation checks each cached node and child list against the tree;
it does not establish the identity of topology outside the indexed region.
Loading never expands the indexed region. Insufficient budgets are rejected,
and allocator spare capacity from the original index is not persisted.

For BP, DFUDS and LOUDS, `TreeTraversalContext<Tree>` borrows the tree and an
optional external root index. Its byte limit (512 KiB by default) determines
the decoding plan. The context exposes root cursors, chunk and depth sizes,
and scratch accounting. The whole-walk benchmarks in `src/benchmarks/tree_dfs`
use the internal decoding kernels; the context has no whole-tree DFS event API.

```cpp
#include <pixie/tree/traversal_context.h>

pixie::tree::RootNavigationIndex index(tree, 32768);
pixie::tree::TreeTraversalContext context(tree, 512 * 1024, &index);
auto part = context.fetch(tree.root(), 1024);
```

The index is independently owned and budgeted. `root_index()` returns a nullable
const pointer. Save/restore the index separately and pass the restored index to
another context; no constructor transfers index ownership. Cursors use native
navigation when no index is provided.

`fetch(node, max_bits)` returns a `TreePart` by value. A nonempty result has one
local root, `nodes()[0]`, and can outlive both the context and the source tree.
BP reads preorder parentheses: an opening bit admits the node and its matching
close confirms its children. DFUDS reads preorder unary degrees: starting a
record admits the node, and its zero terminator establishes its degree. Its list
is complete once every immediate child's record has started, even if the last
child is incomplete. LOUDS uses the same unary rule in breadth-first order.
Its native navigation to locate child runs is separate from the logical bits
examined across those runs. The reader stops at the requested subtree boundary.

Bit and byte limits can return partial or empty parts. The context byte limit
covers the context, its whole-walk scratch allowance, reusable fetch scratch,
and the nodes and child links of one output under construction. Previously
returned parts, the external index, allocator overhead and stack frames are
excluded. The context retains no output. Operations on a single context must
not overlap; sequential thread handoff is allowed.

`TreeBase::parentheses(start, count)` streams a canonical parenthesis slice for
all encodings without allocating a full intermediate string.

## Public objects and experiments {#rooted_tree_public_objects}

The following objects are in `pixie::tree`. Include the header for the object
you use; [implementations.h](implementations.h) lists the four promoted encodings.

| Header | Public objects | Role |
| --- | --- | --- |
| [tree.h](../tree.h) | `TreeBase`, `TreeTraits`, `TreeNavigationView`, `TreeNode` | Navigation contract, representation-specific handle/view types, view requirements, and the number/bit-position handle for contiguous encodings. Use `Tree::Node` and `Tree::NodeView` in generic code. |
| [bp.h](bp.h), [dfuds.h](dfuds.h), [louds.h](louds.h) | `BPTree`, `DFUDSTree`, `LoudsTree`, `LoudsNode` | Owning or borrowing contiguous encodings; `LoudsNode` names the same type as `TreeNode`. Each tree exposes its borrowing representation as `TreeView`. |
| [forest_cover.h](forest_cover.h) | `ForestCoverTree` and its `TreeTraversalContext` and `RootNavigationIndex` specializations | Owning bounded forests, caller-owned decoder slots and optional root index. Use `ForestCoverTree::Node` for handles. |
| [path_view.h](path_view.h) | `TreePathView` | Positioned navigation with bounded path state (`Tree::PathView`, `Tree::BoundedPathView<Capacity>`); use `tree_path` to materialize a path between endpoints. |
| [root_navigation_index.h](root_navigation_index.h), [traversal_context.h](traversal_context.h), [traversal_budget.h](traversal_budget.h) | `RootNavigationIndex`, `TreeTraversalContext`, `kDefaultTraversalBytes` | Bounded navigation and traversal for the contiguous encodings. The default allowance is 512 KiB. |
| [tree_part.h](tree_part.h), [event_part.h](event_part.h) | `TreePart`, `inspect_event_part` | Owning snapshot of nodes and ordered child links. A bounded BP/DFUDS probe returns the complete subtree topology or `std::nullopt`. |

The following headers are opt-in experiments in `pixie::tree::experimental`.
None is required to navigate a promoted tree.

| Header | Objects | Experiment |
| --- | --- | --- |
| [decomposed_tree.h](experimental/decomposed_tree.h) | `StaticDecomposedTree`, `BTreeDecomposedTree`, `MacroDecomposedTree`, `DecomposedNode`, `DecomposedMemoryBreakdown` | Two routing schemes over a first-child/next-sibling decomposition. |
| [layered_dfuds.h](experimental/layered_dfuds.h) | `LayeredDFUDSTree`, `LayeredDFUDSNode` | Depth bands with DFUDS forests and inter-layer connectors. |
| [tree_prepared_root.h](experimental/tree_prepared_root.h), [tree_prepared_events.h](experimental/tree_prepared_events.h) | `PreparedRootSubtree`, `PreparedTreeEvents` | Precomputed navigation and traversal replay for small trees. |
| [simd_micro_traversal.h](experimental/simd_micro_traversal.h) | `BpSimdMicroTraversal`, `DfudsSimdMicroTraversal` | Alternative SIMD traversal kernels. |
| [microtree_simd.h](../experimental/microtree_simd.h) | `microtree_simd::Mask128` and bit-mask operations | Low-level bounded loads, BP leaf/child scans and DFUDS record-low masks used by the experimental kernels. |

The forest decoder, its slot cache and the event decoders under `tree/detail`
are implementation details. The adjacency-list and pointer trees under
`include/references` are comparison implementations, not promoted encodings.
