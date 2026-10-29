# Tree navigation and traversal contexts

BP, DFUDS and LOUDS expose the same `TreeBase` navigation contract.
Their encoded words are borrowed: keep the words alive and immutable for the
tree's lifetime. Use the concrete header for the chosen representation.

```cpp
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

// words contains a valid LOUDS encoding of node_count nodes.
pixie::LoudsTree<> tree(words, node_count);
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
`(tree, words, byte_limit)`; the default is 512 KiB, and 32 KiB retains
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
