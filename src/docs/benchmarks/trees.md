# Trees {#benchmarks_trees}

These charts compare four succinct ordered-tree representations and three reference implementations.

| Implementation | Representation and navigation | Storage |
| --- | --- | --- |
| BP | Balanced parentheses in depth-first order; a range min-max index locates matching and enclosing parentheses. | Two topology bits per node, plus the index. |
| DFUDS | Unary node degrees in depth-first order; a range min-max index supports navigation. | About two topology bits per node, plus the index. |
| LOUDS | Unary node degrees in breadth-first order; rank/select locates parent and child records. | About two topology bits per node, plus the index. |
| ForestCover | Small LOUDS forests of up to 256 nodes are decoded during navigation. | Packed topology plus unit directories and links between units. |
| Adjacency list (vectors) | Each node's vector stores its parent ID followed by its ordered child IDs. | About 320 bits per node on the benchmark machine, excluding allocator overhead. |
| Pointers | Each node stores first-child, next-sibling, and parent pointers. | Three pointers per node (192 bits on the benchmark machine). |
| BFS CSR | Children occupy consecutive BFS IDs. Two packed arrays store child offsets and parent IDs. | Two fields of `ceil(log2(n))` bits per node, plus word padding and the tree object. |

Each matched point is an equal-weight mean over three tree families at the same node count:

- **Random recursive:** each new node attaches to a uniformly chosen earlier node.
- **Binary Catalan:** sample uniformly from binary tree shapes of this size. Left and right child positions are distinct.
- **Ordered plane Catalan:** sample uniformly from ordered rooted tree shapes of this size. A node may have any number of children.

The baseline, 32 KiB and 512 KiB curves vary extra navigation state, not the tree itself. Use the curve selector to compare budgets and hover over a point for its value. The matched path charts also report the mean number of nodes per path.

\htmlonly
<link rel="stylesheet" href="tree-curves.css">
<script src="d3-7.9.0.min.js"></script>
<script src="plot-0.6.17.min.js"></script>
<script src="tree-matched.js"></script>
<script src="tree-history.js"></script>
<script src="tree-curves.js"></script>
<h2 id="matched-inspected_dfs">Inspected DFS</h2>
<p>Walk the whole tree, enumerating each node's children before descending. CPU time is per node.</p>
<div class="pixie-tree-chart" data-chart="matched-inspected_dfs"></div>
<h2 id="matched-inspected_root_path">Inspected root-to-leaf paths</h2>
<p>Inspect the children, then choose one and descend. CPU time is per completed root-to-leaf path.</p>
<div class="pixie-tree-chart" data-chart="matched-inspected_root_path"></div>
<h2 id="matched-tree_path">Random endpoint paths</h2>
<p>Return the ordered path between each of 256 endpoint pairs, including both endpoints. Timing includes writing the path into a reusable buffer.</p>
<p>All three references use parent-before-child numbering, a property shared by BFS order and DFS preorder. Each query moves the higher-numbered endpoint up one edge until the two meet at their lowest common ancestor. This avoids walking to the root when the path lies entirely below it.</p>
<div class="pixie-tree-chart" data-chart="matched-tree_path"></div>

<h2 id="build-time">Build time</h2>
<p>Build and storage use seed-42 random recursive trees. BFS CSR, all uncached succinct trees and all ForestCover configurations were measured on November 20. The remaining curves are from October 30 (revision <code>8e72c47</code>).</p>
<p>Build time includes making the tree ready for queries from prepared input: packed words for BP, DFUDS and LOUDS, or ordered adjacency for ForestCover and the references. Times are median wall milliseconds over three fresh processes, excluding input generation and final destruction.</p>
<div class="pixie-tree-chart" data-chart="history-build-prepared"></div>
<h2 id="storage">Storage per node</h2>
<p>Bits per node include the tree encoding and its intrinsic navigation indexes. Uncached variants have no traversal context or optional root index. Budgeted variants also count retained traversal state and reserved scratch. Allocator overhead and temporary construction memory are excluded.</p>
<div class="pixie-tree-chart" data-chart="history-memory"></div>

<details id="tree-history"><summary>Historical traversal measurements: random recursive trees</summary>
<p>The BFS CSR run uses the same seed-42 random recursive workload as the archived curves from revision <code>8e72c47</code>. Traversals report median CPU time over five repetitions on a Ryzen 7 8845HS, GCC 13.3, CPU 0.</p>
<h3>DFS</h3><div class="pixie-tree-chart" data-chart="history-dfs"></div>
<h3>DFS with child inspection</h3><div class="pixie-tree-chart" data-chart="history-dfs-children"></div>
<h3>Random root-to-leaf paths</h3><div class="pixie-tree-chart" data-chart="history-random-path"></div>
<h3>Random root-to-leaf paths with child inspection</h3><div class="pixie-tree-chart" data-chart="history-random-path-children"></div>
</details>
\endhtmlonly
