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
Serialization is not implemented for these tree types.

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

`TreeTraversalContext<Tree>` in `pixie/tree/traversal_context.h` combines
this index with a decoding buffer plan. Construct it with
`(tree, byte_limit)`; the default is 512 KiB, and 32 KiB retains
optimized decoding. It exposes root cursors, selected chunk/depth sizes and
scratch accounting. It does not expose a public DFS event API. The selected
decoding kernels are internal implementation details under `tree/detail`;
the four benchmark workloads orchestrate them in `src/benchmarks/tree_dfs`.

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
