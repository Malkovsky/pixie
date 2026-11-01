# Tree navigation and traversal contexts

BP, DFUDS and LOUDS expose the same `TreeBase` navigation contract.
Trees own their encoded words in `AlignedStorage` by default. Construction
from a word span copies; construction from an `AlignedStorage` rvalue transfers
its allocation. Use the concrete header for the chosen representation.

Each tree inherits `TreeView`, `NodeView`, `PathView` and
`BoundedPathView<Capacity>` aliases from `TreeBase`. `TreeView` is the same
representation with `ReadOnlyStorageView`: it borrows aligned complete encoding
words, but owns its constructed backend metadata. Its backing bytes must remain
alive, immutable and at a stable address. It can also be built directly from a
read-only storage view. This is distinct from a per-node `NodeView`.

Owning copies duplicate the encoding and rebuild backend references; view copies
share encoding and copy backend metadata. Moves transfer storage and indexes.
Assignment, moving or destruction invalidates existing node views, indexes and
traversal contexts borrowing that tree; using invalidated borrowers is undefined
behavior. Moved-from trees may only be destroyed or assigned.

`encoding_words()` exposes a read-only span of complete encoding words.
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
for example `DFUDSTree<RmMBTree<>> dfuds(bp)`. The source provides `size()`,
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
pixie::LoudsTree<> tree(words, node_count);  // Owns a copy.
using Tree = decltype(tree);
Tree::TreeView representation(tree.encoding_words(), node_count);
Tree::PathView path(tree);
auto view = tree.navigate(tree.root());
const auto degree = view.degree();

pixie::RootNavigationIndex index(tree, 32 * 1024);
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

`TreeTraversalContext<Tree>` in `pixie/tree/traversal_context.h` combines
this index with a decoding buffer plan. Construct it with
`(tree, byte_limit)`; the default is 512 KiB, and 32 KiB retains
optimized decoding. It exposes root cursors, selected chunk/depth sizes and
scratch accounting. It does not expose a public DFS event API. The selected
decoding kernels are internal implementation details under `tree/detail`;
the four benchmark workloads orchestrate them in `src/benchmarks/tree_dfs`.

To persist a context's static cache, serialize `context.root_index()`. After
loading that index, pass it to `TreeTraversalContext<Tree>(std::move(index),
byte_limit)`. This recomputes the small decoding plan and checks the combined
budget without reconstructing cached nodes. Cursors, per-node memoization,
temporary decoded chunks and constexpr microprogram tables are not serialized.

The limit covers active context/index storage and peak explicit traversal
buffers for one walk. Constexpr program tables, original encoding/backend,
allocator metadata, compiler stack frames, callback storage and construction
scratch are excluded. Concurrent or nested walks require separate scratch
allowances. The context itself is immutable and does not retain decoded nodes
between walks.

`budget_usage_bytes()` reports the charge against that limit.
`memory_usage_bytes()` also includes the inert object's fixed bytes when the
budget is zero; for a context it includes reserved scratch, not a persistent
scratch allocation. Zero disables optimization, not the C++ object's existence.

Experimental whole-tree replay remains in
`tree/experimental/tree_prepared_root.h` and `tree_prepared_events.h`.
It targets small trees and is not part of the default traversal strategy.
