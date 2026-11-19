#pragma once

/**
 * @file decomposed_tree.h
 * @brief Static bounded FCNS covers with two experimental routing controllers.
 *
 * @par Experimental snapshot (2026-11-14)
 * Ryzen 7 8845HS under WSL, GCC 13.3, C++20 Release (-O3 -DNDEBUG
 * -march=native), Google Benchmark 1.9.4, pinned CPU 0. Median CPU times from
 * five randomly interleaved repetitions, minimum 0.1 s for large traversals.
 * Random-recursive inputs use parent[i] = rng() % i, not Catalan sampling.
 * The table uses 1,048,576 nodes, seed 42, and threshold 128. Times include
 * navigation, decoding, callbacks and fresh traversal scratch, but exclude
 * source generation, tree construction and optimized traversal-index creation.
 * Rank variants request original preorder rank at every callback. LOUDS uses
 * breadth-first numbering, so no preorder-rank variant is reported for it.
 *
 * | Implementation    | DFS ns/node | DFS+rank | Inspected | Inspected+rank |
 * Retained bits/node |
 * |-------------------|------------:|---------:|----------:|---------------:|-------------------:|
 * | BTree decomposed  |        28.3 |     27.6 |      50.9 |           51.8 |
 * 15.08 | | Macro decomposed  |        27.8 |     30.1 |      52.6
 * |           54.8 |              13.10 | | BP cursor         |        23.9 |
 * 23.7 |      37.9 |           37.9 |               2.37 | | BP optimized | 2.8
 * |      3.0 |       7.7 |            8.4 |               4.98 | | DFUDS cursor
 * |        41.2 |     41.3 |      60.5 |           61.2 |               2.37 |
 * | DFUDS optimized   |         4.8 |      4.9 |       8.6 |            9.2 |
 * 4.91 | | LOUDS cursor      |        49.9 |        - |      64.1 | - | 2.07 |
 *
 * Optimized baseline retained sizes include their immutable traversal indexes;
 * their separate scratch allowances are 182,384 B (BP) and 190,576 B (DFUDS).
 * Each decomposed navigation context is 904 B, with another 904 B context for
 * child enumeration. These are object sizes, not measured peak stack usage;
 * decoder locals and compiler spills are additional temporary storage.
 *
 * Disjoint retained-byte accounting for the same random tree includes vector
 * capacity and object padding, but not allocator bookkeeping or source data:
 *
 * | Category                   | BTree bytes | Macro bytes |
 * |----------------------------|------------:|------------:|
 * | Packed bottom topology     |     262,144 |     262,144 |
 * | Component size/offset data |     102,400 |     102,400 |
 * | Explicit boundary links    |     262,144 |           0 |
 * | External ordinal parents   |     200,176 |     200,176 |
 * | Sparse child-run directory |     655,360 |     655,360 |
 * | Preorder routing/reverse   |     493,544 |     493,544 |
 * | Macro topology/index       |           0 |       3,640 |
 * | Owning object              |         304 |         312 |
 * | Total                      |   1,976,072 |   1,717,576 |
 *
 * Both covers have 12,511 components, 12,510 exits and 15,368 preorder runs.
 * Component sizes are min/mean/max 1/83.8/253; structural macro height is 34,
 * routing height is 2. A separate public-cursor DFS trace counted 19,431
 * component transitions, not internal decoder calls. Allocation tracking
 * matched every retained allocation; constructor peak requested allocation
 * was 78,082,264 B (BTree) and 77,823,760 B (Macro), excluding preexisting
 * source storage, stack and allocator bookkeeping. These are not RSS figures.
 *
 * Across 1,048,576 and 4,194,304 random-recursive nodes and seeds 42/1729,
 * unranked DFS was 28.1-28.3 ns/node (BTree), 27.4-28.0 (Macro), and 2.8
 * (optimized BP); all large traversal CPU CVs were below 3.6%. At 1,048,576
 * nodes, random preorder lookup was 50/48 ns (BTree/Macro), handle-to-rank
 * 2.9/2.8 ns, and fresh cursor plus leaf test 1.15/1.16 us. The last probe
 * includes decoding and navigation checks, not an isolated decoder kernel.
 * At 65,536 random nodes, warm local-edge probes took about 11.5 ns for either
 * controller; the selected boundary-edge pair took 1.02/1.04 us per pair.
 *
 * Shape matters: at 65,536 nodes, inspected chain traversal was 40/41 ns/node
 * (BTree/Macro), versus 355 for optimized BP; inspected comb traversal was
 * 39/42 versus 186. On stars it was 33/35 versus 3.4. Root-path rankings also
 * differed: at that size LOUDS used about 17 ns/visited node on stars versus
 * about 400 for either controller. Small local wins do not imply a universal
 * traversal improvement. Metadata dominates the two-bit bottom topology.
 * Construction was measured separately, including destruction. Even after
 * a bounded 0.5 s follow-up, several construction CPU CVs exceeded 5%, so
 * no close construction ranking is retained here. These remain experiments.
 *
 * Reproduce through the local /bench wrapper with target
 * tree_experimental_benchmarks, five repetitions and 0.1s minimum, using:
 * @code{.text}
 * ^cover/(btree|macro|bp_cursor|bp_optimized|dfuds_cursor|dfuds_optimized|louds_cursor)/(dfs|dfs_rank|inspected_dfs|inspected_dfs_rank)/random_recursive/N:1048576/seed:42$
 * @endcode
 * Other registered cover/ workloads are build_destroy, lookup, cursor_rank,
 * fetch_plus_leaf, local_edge_warm, boundary_edge_pair and root_path. Probe
 * lookup/rank/fetch times cover 8,192 queries per iteration; root_path times
 * cover one path, with visited-node throughput reported separately. Raw
 * result files are intentionally not retained in the repository.
 */

#include <pixie/tree/detail/construction.h>
#include <pixie/tree/experimental/decomposed_tree_core.h>

#include <optional>
#include <type_traits>
#include <utility>

namespace pixie::tree::experimental {

/** @brief Canonical disjoint-ownership address, not a global BP bit position.
 */
struct DecomposedNode {
  size_t component;
  uint8_t local;
  bool operator==(const DecomposedNode&) const = default;
};

/** @brief Disjoint retained-byte categories, including allocation capacity. */
struct DecomposedMemoryBreakdown {
  size_t object_bytes = 0;
  size_t bottom_topology_bytes = 0;
  size_t component_directory_bytes = 0;
  size_t boundary_bytes = 0;
  size_t external_parent_bytes = 0;
  size_t child_directory_bytes = 0;
  size_t preorder_bytes = 0;
  size_t macro_bytes = 0;
  size_t total_bytes = 0;
};

/**
 * @brief Shared static cover; use BTreeDecomposedTree or MacroDecomposedTree.
 * @tparam Threshold Cover threshold in [2,128], not the maximum component size.
 * @details Components have at most 2*Threshold-1 nodes. Handles belong to one
 * tree and remain valid while it is unchanged. Numerically invalid handles and
 * nonexistent navigation targets throw std::out_of_range. Foreign handles with
 * coincidentally valid coordinates cannot be detected. Copying a tree preserves
 * its coordinates; contexts borrow the specific original object.
 *
 * This is an immutable experiment, without serialization or TreeBase encoding
 * compatibility. Construction uses O(n) temporary storage, released before the
 * constructor returns. All retained topology is the packed local BP forests:
 * exactly 2n logical bits, with final-word padding only. Metadata is O(n/T)
 * records, not a second n-node encoding or address table.
 */
template <size_t Threshold, bool Macro>
class StaticDecomposedTree {
  static_assert(Threshold >= 2 && Threshold <= 128);
  using Routing = std::conditional_t<Macro,
                                     decomposition_detail::MacroRouting,
                                     decomposition_detail::BalancedRouting>;
  static constexpr size_t kAbsent = decomposition_detail::kAbsent;
  static constexpr size_t kOffsetStride = 32;

  struct Component {
    uint16_t nodes = 0;
    uint8_t exits = 0;
    std::array<uint8_t, 2> source{};
    uint8_t next_edges = 0;
  };
  struct ChildRun {
    DecomposedNode first;
    size_t end;
  };
  struct ChildDirectory {
    DecomposedNode parent;
    size_t begin;
    size_t end;
  };

  size_t size_ = 0;
  size_t macro_height_ = 0;
  std::vector<uint64_t> topology_;
  std::vector<Component> components_;
  std::vector<size_t> offset_samples_;
  std::vector<DecomposedNode> external_parents_;
  std::vector<decomposition_detail::PreorderRun> runs_;
  std::vector<std::array<size_t, 2>> reverse_runs_;
  std::vector<ChildDirectory> child_directories_;
  std::vector<ChildRun> child_runs_;
  Routing routing_;

  static bool less_node(DecomposedNode a, DecomposedNode b) {
    return a.component < b.component ||
           (a.component == b.component && a.local < b.local);
  }

  void check_node(DecomposedNode node) const {
    if (node.component >= components_.size() ||
        node.local >= components_[node.component].nodes) {
      throw std::out_of_range("Invalid decomposed node");
    }
  }

  size_t bit_offset(size_t component) const {
    const size_t begin = component / kOffsetStride * kOffsetStride;
    size_t offset = offset_samples_[component / kOffsetStride];
    for (size_t i = begin; i < component; ++i) {
      offset += 2 * components_[i].nodes;
    }
    return offset;
  }

  const ChildDirectory* child_directory(DecomposedNode node) const {
    const auto found = std::lower_bound(
        child_directories_.begin(), child_directories_.end(), node,
        [](const ChildDirectory& entry, DecomposedNode key) {
          return less_node(entry.parent, key);
        });
    return found != child_directories_.end() && found->parent == node ? &*found
                                                                      : nullptr;
  }

  template <pixie::tree::detail::OrderedTreeSource Source>
  void construct(const Source& source) {
    using decomposition_detail::checked_add;
    size_ = source.size();
    if (size_ > std::numeric_limits<size_t>::max() / 2) {
      throw std::length_error("Decomposed topology bit count overflow");
    }
    if (!size_) {
      return;
    }
    struct Original {
      size_t first = kAbsent;
      size_t next = kAbsent;
      size_t parent = kAbsent;
      size_t last = kAbsent;
      size_t subtree = 0;
      size_t residual = 0;
      bool marked = false;
    };
    std::vector<Original> original;
    original.reserve(size_);
    std::vector<size_t> stack;
    pixie::tree::detail::construction_dfs(
        source,
        [&](auto, size_t) {
          if (original.size() == size_) {
            throw std::invalid_argument("Source size disagrees with topology");
          }
          const size_t node = original.size();
          original.emplace_back();
          if (!stack.empty()) {
            const size_t parent = stack.back();
            original[node].parent = parent;
            if (original[parent].first == kAbsent) {
              original[parent].first = node;
            } else {
              original[original[parent].last].next = node;
            }
            original[parent].last = node;
          }
          stack.push_back(node);
        },
        [&](auto) { stack.pop_back(); });
    if (original.size() != size_) {
      throw std::invalid_argument("Source size disagrees with topology");
    }

    // Binary preorder is ordinal preorder, so reverse preorder evaluates both
    // FCNS children before the parent. Original, not residual, sizes decide the
    // two-heavy case (Farzan--Munro's binary covering rule).
    //
    // At a non-fork node only one original child can contain a marked subtree:
    // the other has <T original nodes. Thus an uncut residual has at most one
    // exit. At a two-heavy node both children are cut and v is a singleton.
    // Residuals are <T; ordinary closed components have <=2T-1 nodes.
    // Components of size >=T can be charged disjointly to their nodes. Extra
    // small forced roots occur only at branch points of the pruned
    // heavy-subtree tree, whose branch count is bounded by its large leaves;
    // the root adds at most one.
    for (size_t i = size_; i-- > 0;) {
      auto& node = original[i];
      const auto subtree = [&](size_t child) {
        return child == kAbsent ? size_t{0} : original[child].subtree;
      };
      const auto residual = [&](size_t child) {
        return child == kAbsent ? size_t{0} : original[child].residual;
      };
      node.subtree =
          checked_add(1, checked_add(subtree(node.first), subtree(node.next)));
      const size_t count = 1 + residual(node.first) + residual(node.next);
      if (subtree(node.first) >= Threshold && subtree(node.next) >= Threshold) {
        node.marked = true;
        original[node.first].marked = true;
        original[node.next].marked = true;
        node.residual = 0;
      } else if (count >= Threshold) {
        node.marked = true;
        node.residual = 0;
      } else {
        node.residual = count;
      }
    }
    original[0].marked = true;

    std::vector<DecomposedNode> address(size_);
    std::vector<size_t> roots{0};
    std::vector<size_t> depths{1};
    std::vector<std::array<size_t, 2>> destinations;
    std::vector<uint8_t> macro_degrees;
    std::vector<size_t> affected_parents;
    topology_.resize(topology_bits() / 64 + (topology_bits() % 64 != 0));
    size_t bit = 0;
    struct EdgeTask {
      size_t node;
      uint8_t source;
      bool next;
    };
    std::vector<EdgeTask> pending;
    std::vector<size_t> local_nodes;
    // Queue roots in cut traversal order. This assigns BFS component IDs for
    // LOUDS directly, with identical coordinates in both controller variants.
    for (size_t component = 0; component < roots.size(); ++component) {
      Component descriptor;
      std::array<size_t, 2> links{kAbsent, kAbsent};
      local_nodes.clear();
      pending.push_back({roots[component], 0, false});
      while (!pending.empty()) {
        const auto task = pending.back();
        pending.pop_back();
        const auto& node = original[task.node];
        if (task.node != roots[component] && node.marked) {
          if (descriptor.exits == 2) {
            throw std::logic_error("FCNS cover has more than two exits");
          }
          const size_t ordinal = descriptor.exits++;
          descriptor.source[ordinal] = task.source;
          descriptor.next_edges |= static_cast<uint8_t>(task.next << ordinal);
          links[ordinal] = roots.size();
          roots.push_back(task.node);
          depths.push_back(depths[component] + 1);
          const size_t source_node = local_nodes[task.source];
          affected_parents.push_back(task.next ? original[source_node].parent
                                               : source_node);
          continue;
        }
        if (local_nodes.size() == 2 * Threshold - 1) {
          throw std::logic_error("FCNS component exceeds its bound");
        }
        const auto local = static_cast<uint8_t>(local_nodes.size());
        address[task.node] = {component, local};
        local_nodes.push_back(task.node);
        if (node.next != kAbsent) {
          pending.push_back({node.next, local, true});
        }
        if (node.first != kAbsent) {
          pending.push_back({node.first, local, false});
        }
      }
      descriptor.nodes = static_cast<uint16_t>(local_nodes.size());
      if (descriptor.exits > 1 && descriptor.nodes != 1) {
        throw std::logic_error("FCNS two-exit component is not a singleton");
      }
      components_.push_back(descriptor);
      destinations.push_back(links);
      macro_degrees.push_back(descriptor.exits);
      macro_height_ = std::max(macro_height_, depths[component]);
      if (component % kOffsetStride == 0) {
        offset_samples_.push_back(bit);
      }
      // Encode forest: open, local first-child forest, close, local siblings.
      // The temporary event stack is bounded by the local component size.
      struct Event {
        size_t node;
        bool close;
      };
      std::vector<Event> events{{roots[component], false}};
      while (!events.empty()) {
        const auto event = events.back();
        events.pop_back();
        if (event.close) {
          ++bit;
          continue;
        }
        topology_[bit / 64] |= uint64_t{1} << (bit % 64);
        ++bit;
        const auto& node = original[event.node];
        if (node.next != kAbsent && !original[node.next].marked) {
          events.push_back({node.next, false});
        }
        events.push_back({event.node, true});
        if (node.first != kAbsent && !original[node.first].marked) {
          events.push_back({node.first, false});
        }
      }
    }
    if (bit != topology_bits()) {
      throw std::logic_error("FCNS cover lost real nodes");
    }
    external_parents_.reserve(roots.size());
    for (size_t root : roots) {
      external_parents_.push_back(root ? address[original[root].parent]
                                       : DecomposedNode{0, 0});
    }
    reverse_runs_.resize(components_.size(), {kAbsent, kAbsent});
    for (size_t i = 0; i < size_;) {
      const auto first = address[i];
      size_t end = i + 1;
      while (end < size_ && address[end].component == first.component &&
             address[end].local == first.local + end - i) {
        ++end;
      }
      auto& reverse = reverse_runs_[first.component];
      if (reverse[0] == kAbsent) {
        reverse[0] = runs_.size();
      } else if (reverse[1] == kAbsent) {
        reverse[1] = runs_.size();
      } else {
        throw std::logic_error(
            "FCNS component has more than two preorder runs");
      }
      runs_.push_back(
          {first.component, static_cast<uint16_t>(end - i), first.local});
      i = end;
    }
    routing_.build(runs_, macro_degrees, std::move(destinations));

    std::sort(affected_parents.begin(), affected_parents.end());
    affected_parents.erase(
        std::unique(affected_parents.begin(), affected_parents.end()),
        affected_parents.end());
    for (size_t parent : affected_parents) {
      if (parent == kAbsent) {
        throw std::logic_error("Root has an external sibling");
      }
      const size_t begin = child_runs_.size();
      size_t count = 0;
      for (size_t child = original[parent].first; child != kAbsent;) {
        const auto first = address[child];
        do {
          ++count;
          child = original[child].next;
        } while (child != kAbsent &&
                 address[child].component == first.component);
        child_runs_.push_back({first, count});
      }
      child_directories_.push_back(
          {address[parent], begin, child_runs_.size()});
    }
    std::sort(child_directories_.begin(), child_directories_.end(),
              [](const auto& a, const auto& b) {
                return less_node(a.parent, b.parent);
              });
  }

 public:
  /** @brief Immutable canonical node handle. */
  using Node = DecomposedNode;

  /**
   * @brief Independent copyable navigation context with one decoded component.
   * @details Borrows an unmoved, unchanged tree that must outlive it. Const
   * operations may replace this context's scratch, never its node(). Separate
   * contexts may query concurrently; sharing one requires synchronization.
   */
  class NodeView {
    const StaticDecomposedTree* tree_;
    Node node_;
    mutable size_t cached_ = kAbsent;
    mutable decomposition_detail::DecodedForest decoded_;

    void fetch(size_t component) const {
      if (cached_ != component) {
        decoded_.decode(tree_->topology_, tree_->bit_offset(component),
                        tree_->components_[component].nodes);
        cached_ = component;
      }
    }

    std::optional<Node> edge(Node node, bool next) const {
      fetch(node.component);
      const auto& valid =
          next ? decoded_.has_next_sibling : decoded_.has_first_child;
      const auto& link = next ? decoded_.next_sibling : decoded_.first_child;
      if (valid[node.local]) {
        return Node{node.component, link[node.local]};
      }
      const auto& component = tree_->components_[node.component];
      for (size_t i = 0; i < component.exits; ++i) {
        if (component.source[i] == node.local &&
            bool((component.next_edges >> i) & 1) == next) {
          return Node{tree_->routing_.destination(node.component, i), 0};
        }
      }
      return std::nullopt;
    }

    friend class StaticDecomposedTree;
    NodeView(const StaticDecomposedTree& tree, Node node)
        : tree_(&tree), node_(node) {
      tree.check_node(node);
    }

   public:
    /** @brief Current handle; queries do not change it. */
    Node node() const { return node_; }
    /** @brief Set a checked handle, retaining scratch if its component matches.
     */
    void reset(Node node) {
      tree_->check_node(node);
      node_ = node;
    }
    /** @brief Whether the current node is the global root. */
    bool is_root() const { return node_ == Node{0, 0}; }
    /** @brief Whether the current node has no original children. */
    bool is_leaf() const { return !edge(node_, false); }
    /** @brief Original child count, using a sparse directory across cuts. */
    size_t degree() const {
      if (const auto* directory = tree_->child_directory(node_)) {
        return tree_->child_runs_[directory->end - 1].end;
      }
      size_t result = 0;
      for (auto child = edge(node_, false); child; child = edge(*child, true)) {
        ++result;
      }
      return result;
    }
    /** @brief Zero-based child in [0,degree()); otherwise throws out_of_range.
     */
    Node child(size_t index) const {
      if (const auto* directory = tree_->child_directory(node_)) {
        const auto first = tree_->child_runs_.begin() + directory->begin;
        const auto last = tree_->child_runs_.begin() + directory->end;
        const auto run = std::upper_bound(
            first, last, index,
            [](size_t key, const ChildRun& value) { return key < value.end; });
        if (run == last) {
          throw std::out_of_range("Child index out of range");
        }
        const size_t base = run == first ? 0 : (run - 1)->end;
        auto result = run->first;
        fetch(result.component);
        for (size_t i = base; i < index; ++i) {
          result.local = decoded_.next_sibling[result.local];
        }
        return result;
      }
      auto result = edge(node_, false);
      while (result && index) {
        result = edge(*result, true);
        --index;
      }
      if (!result) {
        throw std::out_of_range("Child index out of range");
      }
      return *result;
    }
    /** @brief First child; throws out_of_range on a leaf. */
    Node first_child() const {
      const auto result = edge(node_, false);
      if (!result) {
        throw std::out_of_range("Leaf has no first child");
      }
      return *result;
    }
    /** @brief Original ordinal parent, with parent(root)==root. */
    Node parent() const {
      fetch(node_.component);
      return decoded_.has_parent[node_.local]
                 ? Node{node_.component, decoded_.parent[node_.local]}
                 : tree_->external_parents_[node_.component];
    }
    /** @brief Next sibling; throws out_of_range for the last child or root. */
    Node next_sibling() const {
      const auto result = edge(node_, true);
      if (!result) {
        throw std::out_of_range("No next sibling");
      }
      return *result;
    }
    /** @brief Whether no next sibling exists; true for the root. */
    bool is_last_child() const { return !edge(node_, true); }
    /**
     * @brief Visit child handles in order, invoking visitor as a borrowed
     * lvalue.
     * @details Independent local scratch supports nested enumeration, including
     * callbacks using this same context. Exceptions propagate and stop visits.
     * No callback or traversal state is retained and no heap allocation occurs.
     */
    template <class Visitor>
    void for_each_child(Visitor&& visitor) const {
      auto walk = *this;
      for (auto child = walk.edge(node_, false); child;
           child = walk.edge(*child, true)) {
        visitor(*child);
      }
    }
  };

  /** @brief Construct an empty tree. */
  StaticDecomposedTree() = default;
  /**
   * @brief Own a static copy of a valid ordered source, preserving child order.
   * @details No source data or source handle is retained. Source traversal and
   * covering are iterative, including chains and high-degree sibling lists.
   * Unrepresentable sizes throw length_error; allocation exceptions propagate.
   */
  template <pixie::tree::detail::OrderedTreeSource Source>
    requires(!std::same_as<std::remove_cvref_t<Source>, StaticDecomposedTree>)
  explicit StaticDecomposedTree(const Source& source) {
    construct(source);
  }

  /** @brief Node count. */
  size_t size() const { return size_; }
  /** @brief Whether no root exists. */
  bool empty() const { return size_ == 0; }
  /** @brief Canonical root; throws out_of_range on an empty tree. */
  Node root() const {
    if (empty()) {
      throw std::out_of_range("Empty decomposed tree has no root");
    }
    return {0, 0};
  }
  /** @brief Checked zero-based original preorder lookup in O(log n). */
  Node find_by_preorder(size_t rank) const {
    if (rank >= size_) {
      throw std::out_of_range("Preorder rank out of range");
    }
    const size_t run = routing_.find(rank);
    return {runs_[run].component,
            static_cast<uint8_t>(runs_[run].local_start + rank -
                                 routing_.start(run))};
  }
  /** @brief Checked handle's zero-based original preorder rank in O(1). */
  size_t preorder_number(Node node) const {
    check_node(node);
    const auto& reverse = reverse_runs_[node.component];
    const size_t run =
        reverse[1] != kAbsent && node.local >= runs_[reverse[1]].local_start
            ? reverse[1]
            : reverse[0];
    return routing_.start(run) + node.local - runs_[run].local_start;
  }
  /** @brief Borrow this tree with independent one-component decode scratch. */
  NodeView navigate(Node node) const& { return NodeView(*this, node); }
  /** @brief Reject contexts borrowing a temporary tree. */
  NodeView navigate(Node) const&& = delete;
  /** @brief See NodeView::is_root. */
  bool is_root(Node node) const { return navigate(node).is_root(); }
  /** @brief See NodeView::is_leaf. */
  bool is_leaf(Node node) const { return navigate(node).is_leaf(); }
  /** @brief See NodeView::degree. */
  size_t degree(Node node) const { return navigate(node).degree(); }
  /** @brief See NodeView::child. */
  Node child(Node node, size_t index) const {
    return navigate(node).child(index);
  }
  /** @brief See NodeView::first_child. */
  Node first_child(Node node) const { return navigate(node).first_child(); }
  /** @brief See NodeView::parent. */
  Node parent(Node node) const { return navigate(node).parent(); }
  /** @brief See NodeView::next_sibling. */
  Node next_sibling(Node node) const { return navigate(node).next_sibling(); }
  /** @brief See NodeView::is_last_child. */
  bool is_last_child(Node node) const { return navigate(node).is_last_child(); }
  /** @brief See NodeView::for_each_child. */
  template <class Visitor>
  void for_each_child(Node node, Visitor&& visitor) const {
    navigate(node).for_each_child(std::forward<Visitor>(visitor));
  }

  /** @brief Logical bottom topology bits, excluding final-word padding. */
  size_t topology_bits() const { return 2 * size_; }
  /** @brief Number of disjoint local components. */
  size_t component_count() const { return components_.size(); }
  /** @brief Checked component occupancy, at most 2*Threshold-1. */
  size_t component_size(size_t component) const {
    return components_.at(component).nodes;
  }
  /** @brief Checked cut count; two exits imply a singleton fork. */
  size_t component_exit_count(size_t component) const {
    return components_.at(component).exits;
  }
  /** @brief Number of cut FCNS edges. */
  size_t exit_count() const { return empty() ? 0 : component_count() - 1; }
  /** @brief Number of global preorder runs, at most twice the component count.
   */
  size_t preorder_run_count() const { return runs_.size(); }
  /** @brief Maximum structural component depth, counting the root as one. */
  size_t macro_height() const { return macro_height_; }
  /** @brief Balanced levels or sampled-directory stages; zero on empty trees.
   */
  size_t routing_height() const { return routing_.height(); }
  /** @brief Retained bytes in one caller-owned context, excluded from tree
   * bytes. */
  static constexpr size_t navigation_scratch_bytes() {
    return sizeof(NodeView);
  }

  /**
   * @brief Disjoint object and allocation-capacity byte counts.
   * @details Final-word padding and unused vector capacity are included. Object
   * bytes include vector headers, scalar fields and struct padding; allocation
   * categories include all element padding. Allocator bookkeeping, construction
   * temporaries and caller-owned contexts are excluded. Boundary source markers
   * reside in the component-directory category; boundary_bytes holds explicit
   * destinations only. Reverse run references belong to preorder_bytes.
   */
  DecomposedMemoryBreakdown memory_breakdown() const {
    using decomposition_detail::allocated_bytes;
    DecomposedMemoryBreakdown result;
    result.object_bytes = sizeof(*this);
    result.bottom_topology_bytes = allocated_bytes(topology_);
    result.component_directory_bytes =
        allocated_bytes(components_) + allocated_bytes(offset_samples_);
    result.boundary_bytes = routing_.boundary_bytes();
    result.external_parent_bytes = allocated_bytes(external_parents_);
    result.child_directory_bytes =
        allocated_bytes(child_directories_) + allocated_bytes(child_runs_);
    result.preorder_bytes = allocated_bytes(runs_) +
                            allocated_bytes(reverse_runs_) +
                            routing_.preorder_bytes();
    result.macro_bytes = routing_.macro_bytes();
    result.total_bytes = result.object_bytes + result.bottom_topology_bytes +
                         result.component_directory_bytes +
                         result.boundary_bytes + result.external_parent_bytes +
                         result.child_directory_bytes + result.preorder_bytes +
                         result.macro_bytes;
    return result;
  }
  /** @brief Total retained object plus allocation-capacity bytes. */
  size_t memory_usage_bytes() const { return memory_breakdown().total_bytes; }

  /**
   * @brief Check cover, topology, routing and ordinal navigation consistency.
   * @details Diagnostic O(n log n) traversal with O(component count) temporary
   * counters. Returns false on an inconsistent representation; allocation
   * failures may propagate. Does not mutate retained state or cache scratch.
   */
  bool validate_invariants() const {
    if (empty()) {
      return components_.empty() && topology_.empty() && runs_.empty() &&
             child_directories_.empty() && child_runs_.empty();
    }
    if (components_.empty() || external_parents_.size() != component_count() ||
        reverse_runs_.size() != component_count() ||
        offset_samples_.size() !=
            (component_count() + kOffsetStride - 1) / kOffsetStride) {
      return false;
    }
    size_t nodes = 0;
    size_t exits = 0;
    std::vector<size_t> incoming(component_count());
    decomposition_detail::DecodedForest decoded;
    for (size_t c = 0; c < component_count(); ++c) {
      const auto& component = components_[c];
      if (!component.nodes || component.nodes > 2 * Threshold - 1 ||
          component.exits > 2 ||
          (component.exits == 2 && component.nodes != 1) ||
          bit_offset(c) != 2 * nodes) {
        return false;
      }
      decoded.decode(topology_, bit_offset(c), component.nodes);
      nodes += component.nodes;
      for (size_t e = 0; e < component.exits; ++e) {
        const size_t destination = routing_.destination(c, e);
        const auto source = component.source[e];
        if (source >= component.nodes || destination <= c ||
            destination >= component_count() || ++incoming[destination] != 1 ||
            ((component.next_edges >> e) & 1
                 ? decoded.has_next_sibling[source]
                 : decoded.has_first_child[source])) {
          return false;
        }
        ++exits;
      }
    }
    if (nodes != size_ || exits != component_count() - 1 ||
        topology_.size() !=
            topology_bits() / 64 + (topology_bits() % 64 != 0)) {
      return false;
    }
    std::vector<size_t> local_counts(component_count());
    size_t rank = 0;
    for (size_t i = 0; i < runs_.size(); ++i) {
      const auto& run = runs_[i];
      if (run.component >= component_count() || !run.length ||
          run.local_start != local_counts[run.component] ||
          run.length > components_[run.component].nodes - run.local_start ||
          routing_.start(i) != rank || routing_.find(rank) != i ||
          routing_.find(rank + run.length - 1) != i) {
        return false;
      }
      local_counts[run.component] += run.length;
      rank += run.length;
    }
    if (rank != size_) {
      return false;
    }
    for (size_t c = 0; c < component_count(); ++c) {
      if (local_counts[c] != components_[c].nodes) {
        return false;
      }
    }
    // Stackless ordinal DFS independently checks the global run stream, parent
    // overrides, child order and both kinds of boundary links.
    auto current = root();
    auto view = navigate(current);
    rank = 0;
    for (;;) {
      view.reset(current);
      if (preorder_number(current) != rank ||
          find_by_preorder(rank) != current) {
        return false;
      }
      ++rank;
      size_t children = 0;
      bool correct = true;
      view.for_each_child([&](Node child_node) {
        if (parent(child_node) != current ||
            view.child(children) != child_node) {
          correct = false;
        }
        ++children;
      });
      if (!correct || children != view.degree()) {
        return false;
      }
      if (!view.is_leaf()) {
        current = view.first_child();
        continue;
      }
      while (!view.is_root() && view.is_last_child()) {
        current = view.parent();
        view.reset(current);
      }
      if (view.is_root()) {
        break;
      }
      current = view.next_sibling();
      if (rank >= size_) {
        return false;
      }
    }
    return rank == size_ && parent(root()) == root();
  }
};

/** @brief Bounded FCNS cover with balanced weighted preorder routing and links.
 */
template <size_t Threshold = 128>
using BTreeDecomposedTree = StaticDecomposedTree<Threshold, false>;

/** @brief Same cover with LOUDS boundary routing and sampled preorder lookup.
 */
template <size_t Threshold = 128>
using MacroDecomposedTree = StaticDecomposedTree<Threshold, true>;

}  // namespace pixie::tree::experimental
