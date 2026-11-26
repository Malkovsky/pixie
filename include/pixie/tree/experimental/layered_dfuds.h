#pragma once

#include <pixie/rank_select/support.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree/detail/construction.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::tree::experimental {

/** @brief A layer and zero-based local forest-preorder ordinal, not a bit
 * offset.
 * @details Handles belong to their originating topology. Foreign handles are
 * unsupported; numerical bounds are checked. No original-preorder ID is kept.
 */
struct LayeredDFUDSNode {
  size_t layer = 0;
  size_t local = 0;
  bool operator==(const LayeredDFUDSNode&) const = default;
};

/** @brief Owning experimental depth-band DFUDS forests with unary connectors.
 * @details Each forest stores raw degree records 1^degree 0, without an initial
 * open. Every constituent tree has excess -1. Only cut-depth parents with
 * positive degree d contribute a portal and a group of d-1 ones and one zero.
 * Each group position denotes one ordered root in the next forest.
 * The final forest retains all remaining depths and has no outgoing records.
 * The three streams are densely concatenated, without per-layer word padding.
 *
 * Global numbers are layer-major, then local preorder, NOT original preorder.
 * This class intentionally does not inherit TreeBase or expose one contiguous
 * tree encoding. Query indexes borrow only this object's owned streams. There
 * are no retained source handles, per-node arrays, or original-number maps.
 *
 * A cut prefix is planned in O(n) using the current RankSelectSupport
 * allocation sizes, including cache-line rounding and exact population hints.
 * The first unaffordable cut stops splitting; later cuts are not tried. Actual
 * completed outgoing allocations are checked before accepting that prefix. If
 * allocator capacity or support layout differs from the sizing model, an
 * over-budget prefix is discarded and one unsplit forest is built instead.
 * There are at most two encoding builds, never a rebuild per candidate cut. A
 * zero budget always yields one forest for a nonempty source and no outgoing
 * allocations.
 *
 * The boundary budget counts outgoing stream and index allocations plus
 * directory capacity beyond the two entries needed by one forest, minus one
 * bit per crossing child (each relocates a forest edge). The baseline forest
 * RmM index and object are outside this budget, but included in full memory.
 * Consequently this is not a bound on total bits per node.
 *
 * Navigation uses rank/select and RmM searches, not forest decoding. Searches
 * have a proven answer within the layer, or explicitly bounded range queries;
 * excess in adjacent layers cannot supply an answer. Large root ordinals use
 * a bounded binary search rather than narrowing an excess delta to int.
 */
class LayeredDFUDSTree {
 public:
  using Node = LayeredDFUDSNode;
  /** @brief Dense membership bitmap or sorted, fixed-width absolute node IDs.
   */
  enum class PortalMode { Dense, Sparse };

  /** @brief Half-open bit, number, and original-depth ranges of one layer. */
  struct LayerRange {
    size_t forest_begin, forest_end;
    size_t connector_begin, connector_end;
    size_t boundary_begin, boundary_end;
    size_t node_begin, node_end;
    size_t depth_begin, depth_end;
  };

  /** @brief Disjoint retained bytes, including capacity and allocation padding.
   * @details Index fields exclude their embedded objects (in object_bytes).
   * Excludes allocator bookkeeping, source, construction scratch, and caller
   * views; reports allocated bytes, not resident memory or serialized size.
   */
  struct MemoryBreakdown {
    size_t object_bytes = 0;
    size_t forest_bytes = 0;
    size_t connector_bytes = 0;
    size_t boundary_bytes = 0;
    size_t forest_index_bytes = 0;
    size_t connector_index_bytes = 0;
    size_t boundary_index_bytes = 0;
    size_t directory_bytes = 0;
    size_t total_bytes = 0;
  };

  /** @brief Independent, copyable navigation state borrowing an immutable tree.
   * @details Tree destruction, assignment, or move invalidates views.
   * Navigation does not move the cursor; only reset does. Enumeration uses
   * local state, permits nested calls, and leaves the view unchanged if the
   * visitor throws. Lazy scalar caches and one computed destination-position
   * hint are copied independently. Only reset to that destination reuses the
   * hint; other handles reconstruct their position normally. Shared access to
   * one view requires synchronization.
   */
  class NodeView {
   public:
    /** @brief Construct a checked view; invalid bounds throw out_of_range. */
    NodeView(const LayeredDFUDSTree& tree, Node node)
        : tree_(&tree), node_(node) {
      tree_->check(node_);
    }
    /** @brief Return the current handle. */
    Node node() const { return node_; }
    /** @brief Replace the checked handle, retaining a matching position hint.
     * @details Other caches are invalidated. Failure leaves handle and caches
     * unchanged.
     */
    void reset(Node node) {
      NodeView next(*tree_, node);
      if (node == pending_node_) {
        next.position_ = pending_position_;
      }
      *this = next;
    }
    /** @brief Return the layer-major number of the current node. */
    size_t number() const {
      return tree_->directory_[node_.layer].nodes + node_.local;
    }
    /** @brief Return the original ordered-tree degree. */
    size_t degree() const {
      if (boundary()) {
        load_connector();
        return connector_end_ - connector_begin_;
      }
      return delimiter() - position();
    }
    /** @brief Test whether the node has no children. */
    bool is_leaf() const {
      return !tree_->forest_index_.bit(position()) && !boundary();
    }
    /** @brief Test whether this is the original tree root. */
    bool is_root() const { return node_.layer == 0 && node_.local == 0; }
    /** @brief Test whether this is the last child; true for the root. */
    bool is_last_child() const {
      if (is_root()) {
        return true;
      }
      const size_t open = incoming();
      if (open != Rmm::npos) {
        return open == tree_->directory_[node_.layer].forest ||
               !tree_->forest_index_.bit(open - 1);
      }
      return !(*tree_->connector_index_)[root_connector()];
    }
    /** @brief Return the parent, or the root itself for the root. */
    Node parent() const {
      if (is_root()) {
        return node_;
      }
      const size_t open = incoming();
      if (open != Rmm::npos) {
        // [open, position()) is balanced and contains one zero per node.
        const Node result{node_.layer, node_.local - (position() - open) / 2};
        const size_t begin = tree_->directory_[node_.layer].forest;
        const size_t p = !result.local ? begin
                         : open == begin || !tree_->forest_index_.bit(open - 1)
                             ? open
                             : Rmm::npos;
        return remember(result, p);
      }
      const size_t id =
          tree_->portal_id(tree_->connector_index_->rank0(root_connector()));
      const size_t layer = node_.layer - 1;
      return {layer, id - tree_->directory_[layer].nodes};
    }
    /** @brief Return child i (zero-based); missing children throw out_of_range.
     */
    Node child(size_t i) const {
      if (i >= degree()) {
        throw std::out_of_range("Layered DFUDS child");
      }
      if (boundary()) {
        const size_t ordinal =
            connector_begin_ - tree_->directory_[node_.layer].connector + i;
        const size_t layer = node_.layer + 1;
        const size_t p = tree_->forest_root_position(layer, ordinal);
        // The prefix contains ordinal complete raw trees: bits = 2*k-ordinal.
        return remember(
            {layer, (p - tree_->directory_[layer].forest + ordinal) / 2}, p);
      }
      if (!i) {
        // The next preorder record starts at delimiter()+1.
        return remember({node_.layer, node_.local + 1}, delimiter() + 1);
      }
      const size_t start = delimiter() + 1;
      const size_t p = tree_->forest_index_.close(delimiter() - 1 - i) + 1;
      // Skip i raw subtrees, each occupying 2*k-1 bits.
      return remember({node_.layer, node_.local + 1 + (p - start + i) / 2}, p);
    }
    /** @brief Return the first child; a leaf throws out_of_range. */
    Node first_child() const { return child(0); }
    /** @brief Return the next sibling; its absence throws out_of_range. */
    Node next_sibling() const {
      if (is_last_child()) {
        throw std::out_of_range("Layered DFUDS next sibling");
      }
      return next_sibling_unchecked();
    }
    /** @brief Visit each child once in order; propagate visitor exceptions. */
    template <class Visitor>
    void for_each_child(Visitor&& visitor) const {
      const size_t count = degree();
      if (!count) {
        return;
      }
      auto current = *this;
      current.reset(current.first_child());
      for (size_t i = 0; i < count; ++i) {
        visitor(current.node());
        if (i + 1 != count) {
          current.reset(current.next_sibling_unchecked());
        }
      }
    }

   private:
    const LayeredDFUDSTree* tree_;
    Node node_;
    mutable size_t position_ = Rmm::npos;
    mutable size_t delimiter_ = Rmm::npos;
    mutable size_t incoming_ = Rmm::npos;
    mutable size_t connector_begin_ = Rmm::npos;
    mutable size_t connector_end_ = Rmm::npos;
    mutable Node pending_node_;
    mutable size_t pending_position_ = Rmm::npos;
    mutable bool incoming_known_ = false;

    Node remember(Node node, size_t p) const {
      pending_node_ = node;
      pending_position_ = p;
      return node;
    }
    // The caller has established that another original sibling exists.
    Node next_sibling_unchecked() const {
      const size_t start = position();
      // A boundary node is a local leaf even when it has external children.
      const size_t p = !tree_->forest_index_.bit(start)
                           ? start + 1
                           : tree_->forest_index_.fwdsearch(start, -1) + 1;
      return remember({node_.layer, node_.local + (p - start + 1) / 2}, p);
    }
    bool boundary() const { return tree_->boundary(node_, number()); }
    size_t position() const {
      if (position_ == Rmm::npos) {
        position_ = tree_->position(node_, number());
      }
      return position_;
    }
    size_t delimiter() const {
      if (delimiter_ == Rmm::npos) {
        delimiter_ = tree_->forest_index_.select0(number() + 1);
      }
      return delimiter_;
    }
    size_t incoming() const {
      if (!incoming_known_) {
        incoming_ = tree_->incoming(node_, position());
        incoming_known_ = true;
      }
      return incoming_;
    }
    void load_connector() const {
      if (connector_begin_ == Rmm::npos) {
        const size_t rank = tree_->portal_rank(number()) + 1;
        connector_end_ = tree_->connector_index_->select0(rank) + 1;
        connector_begin_ = tree_->connector_start(rank);
      }
    }
    size_t root_connector() const {
      return tree_->root_connector(node_.layer, position(), number());
    }
  };

  /** @brief Construct a valid empty tree with the default policy. */
  LayeredDFUDSTree() = default;

  /** @brief Copy a valid OrderedTreeSource, preserving ordered topology.
   * @param source Valid immutable tree; handles need only survive construction.
   * @param band_height Positive number of depths per proposed nonfinal band.
   * @param boundary_budget_bits_per_node Finite nonnegative boundary allowance
   * per original node, default 0.05. See class accounting and stop policy.
   * @param mode Dense or Sparse portal storage; other values are invalid.
   * @details Construction uses O(n) temporary records and iterative DFS, also
   * for deep chains. Source storage is never retained. Empty sources are valid.
   * @throws std::invalid_argument for invalid policy parameters.
   * @throws std::length_error for unrepresentable sizes; allocation errors
   * propagate. Source topology and size must agree.
   */
  template <pixie::tree::detail::OrderedTreeSource Source>
  explicit LayeredDFUDSTree(const Source& source,
                            size_t band_height = 8,
                            double boundary_budget_bits_per_node = 0.05,
                            PortalMode mode = PortalMode::Dense)
      : size_(source.size()),
        band_height_(band_height),
        budget_(boundary_budget_bits_per_node),
        mode_(mode) {
    if (!band_height_ || !std::isfinite(budget_) || budget_ < 0 ||
        (mode_ != PortalMode::Dense && mode_ != PortalMode::Sparse)) {
      throw std::invalid_argument("Invalid layered DFUDS policy");
    }
    if (size_ > (std::numeric_limits<size_t>::max() - 64) / 4) {
      throw std::length_error("Layered DFUDS size is not representable");
    }
    if (!size_) {
      return;
    }
    std::vector<Record> records;
    records.reserve(size_);
    std::vector<Counts> bands;
    size_t depth = 0;
    pixie::tree::detail::generic_construction_dfs(
        source,
        [&](auto, size_t degree) {
          records.push_back({depth, degree});
          const size_t band = depth / band_height_;
          if (bands.size() <= band) {
            bands.resize(band + 1);
          }
          auto& counts = bands[band];
          ++counts.nodes;
          counts.edges = add(counts.edges, degree);
          if (degree && depth % band_height_ == band_height_ - 1) {
            ++counts.boundaries;
            counts.cross_edges = add(counts.cross_edges, degree);
          }
          depth_count_ = std::max(depth_count_, ++depth);
        },
        [&](auto) { --depth; });
    if (records.size() != size_) {
      throw std::invalid_argument("Layered DFUDS source size mismatch");
    }
    build_with_budget(records, bands);
  }

  /** @brief Deep-copy streams and rebuild indexes against the new storage. */
  LayeredDFUDSTree(const LayeredDFUDSTree& other)
      : size_(other.size_),
        band_height_(other.band_height_),
        depth_count_(other.depth_count_),
        budget_(other.budget_),
        mode_(other.mode_),
        width_(other.width_),
        boundary_count_(other.boundary_count_),
        cross_edges_(other.cross_edges_),
        directory_(other.directory_),
        forest_(other.forest_),
        connectors_(other.connectors_),
        boundaries_(other.boundaries_) {
    build_indexes();
  }
  /** @brief Independent copy; invalidates this tree's existing views. */
  LayeredDFUDSTree& operator=(const LayeredDFUDSTree& other) {
    if (this != &other) {
      LayeredDFUDSTree copy(other);
      *this = std::move(copy);
    }
    return *this;
  }
  /** @brief Transfer allocations and their index spans without rebinding.
   * @details Invalidates source views. Moved-from objects may only be assigned
   * or destroyed. Standard vector moves preserve the transferred addresses.
   */
  LayeredDFUDSTree(LayeredDFUDSTree&&) noexcept = default;
  /** @brief Transfer ownership; invalidates views of both trees. */
  LayeredDFUDSTree& operator=(LayeredDFUDSTree&&) noexcept = default;

  /** @brief Number of original nodes. */
  size_t size() const { return size_; }
  /** @brief Whether the tree has no nodes. */
  bool empty() const { return size_ == 0; }
  /** @brief Original root; an empty tree throws out_of_range. */
  Node root() const {
    check({0, 0});
    return {0, 0};
  }
  /** @brief Checked layer-major number in [0,size()), not original preorder. */
  size_t number(Node node) const {
    check(node);
    return directory_[node.layer].nodes + node.local;
  }
  /** @brief Invert number(); invalid numbers throw out_of_range. */
  Node find_by_number(size_t number) const {
    if (number >= size_) {
      throw std::out_of_range("Layered DFUDS node number");
    }
    const auto it = std::upper_bound(
        directory_.begin(), directory_.end(), number,
        [](size_t n, const Start& start) { return n < start.nodes; });
    const size_t layer = static_cast<size_t>(it - directory_.begin() - 1);
    return {layer, number - directory_[layer].nodes};
  }
  /** @brief Create a checked independent cursor borrowing this tree. */
  NodeView navigate(Node node) const { return NodeView(*this, node); }
  /** @brief True only for the original root; checks numerical bounds. */
  bool is_root(Node node) const { return number(node) == 0; }
  /** @brief Original ordered-tree degree; checks numerical bounds. */
  size_t degree(Node node) const { return navigate(node).degree(); }
  /** @brief True when the original node has no children. */
  bool is_leaf(Node node) const { return navigate(node).is_leaf(); }
  /** @brief Parent handle; the root's parent is the root itself. */
  Node parent(Node node) const { return navigate(node).parent(); }
  /** @brief Return child i in original sibling order, indexed from zero.
   * @throws std::out_of_range for i >= degree(node) or invalid node bounds.
   */
  Node child(Node node, size_t i) const { return navigate(node).child(i); }
  /** @brief Return child zero; leaves throw out_of_range. */
  Node first_child(Node node) const { return child(node, 0); }
  /** @brief True for a last sibling, including the original root. */
  bool is_last_child(Node node) const { return navigate(node).is_last_child(); }
  /** @brief Return the next original sibling; absence throws out_of_range. */
  Node next_sibling(Node node) const { return navigate(node).next_sibling(); }
  /** @brief Visit children once in order; reentrant and exception-transparent.
   * @details Does not mutate the tree or any caller cursor. No heap scratch.
   */
  template <class Visitor>
  void for_each_child(Node node, Visitor&& visitor) const {
    navigate(node).for_each_child(std::forward<Visitor>(visitor));
  }

  /** @brief Retained forest count (zero for an empty tree). */
  size_t layer_count() const {
    return directory_.empty() ? 0 : directory_.size() - 1;
  }
  /** @brief Layer ranges; i >= layer_count() throws out_of_range.
   * @details Final depth_end includes all remaining source depths. The internal
   * terminal directory sentinel is not itself a layer.
   */
  LayerRange layer(size_t i) const {
    if (i >= layer_count()) {
      throw std::out_of_range("Layered DFUDS layer");
    }
    const auto& a = directory_[i];
    const auto& b = directory_[i + 1];
    return {a.forest,
            b.forest,
            a.connector,
            b.connector,
            a.boundary,
            b.boundary,
            a.nodes,
            b.nodes,
            i * band_height_,
            i + 1 == layer_count() ? depth_count_ : (i + 1) * band_height_};
  }
  /** @brief Requested band height; stopping does not change this parameter. */
  size_t band_height() const { return band_height_; }
  /** @brief Requested finite nonnegative boundary allowance per original node.
   */
  double boundary_budget_bits_per_node() const { return budget_; }
  /** @brief Selected portal representation. */
  PortalMode portal_mode() const { return mode_; }
  /** @brief Connector zero count, one per positive-degree portal. */
  size_t boundary_count() const { return boundary_count_; }
  /** @brief Number of relocated edges, equal to the connector length. */
  size_t cross_edge_count() const { return cross_edges_; }
  /** @brief Logical raw forest bits, excluding word padding. */
  size_t forest_bits() const {
    return directory_.empty() ? 0 : directory_.back().forest;
  }
  /** @brief Logical connector bits, excluding word padding. */
  size_t connector_bits() const {
    return directory_.empty() ? 0 : directory_.back().connector;
  }
  /** @brief Logical bitmap bits or packed absolute portal ID bits. */
  size_t boundary_bits() const {
    return directory_.empty() ? 0 : directory_.back().boundary;
  }
  /** @brief Borrow packed LSB-first forest words; final padding is zero.
   * @details All three word spans require a live, unmoved, unassigned owner.
   */
  std::span<const uint64_t> forest_words() const& { return forest_; }
  std::span<const uint64_t> forest_words() const&& = delete;
  /** @brief Read-only packed LSB-first connector words; final padding is zero.
   */
  std::span<const uint64_t> connector_words() const& { return connectors_; }
  std::span<const uint64_t> connector_words() const&& = delete;
  /** @brief LSB-first bitmap or packed ID words, with zero final padding. */
  std::span<const uint64_t> boundary_words() const& { return boundaries_; }
  std::span<const uint64_t> boundary_words() const&& = delete;
  /** @brief Forest plus connector logical bits: 2*size()-1 if nonempty. */
  size_t topology_bits() const { return forest_bits() + connector_bits(); }
  /** @brief Full disjoint retained allocation accounting; see MemoryBreakdown.
   */
  MemoryBreakdown memory_breakdown() const {
    MemoryBreakdown result;
    result.object_bytes = sizeof(*this);
    result.forest_bytes = bytes(forest_);
    result.connector_bytes = bytes(connectors_);
    result.boundary_bytes = bytes(boundaries_);
    result.forest_index_bytes =
        forest_index_.memory_usage_bytes() - sizeof(Rmm);
    result.connector_index_bytes = owned(connector_index_);
    result.boundary_index_bytes = owned(boundary_index_);
    result.directory_bytes = bytes(directory_);
    result.total_bytes = result.object_bytes + result.forest_bytes +
                         result.connector_bytes + result.boundary_bytes +
                         result.forest_index_bytes +
                         result.connector_index_bytes +
                         result.boundary_index_bytes + result.directory_bytes;
    return result;
  }
  /** @brief Object plus all retained allocation capacities, in bytes. */
  size_t memory_usage_bytes() const { return memory_breakdown().total_bytes; }
  /** @brief Actual budgeted bits, including outgoing padding and support.
   * @details Zero for no accepted cuts. Subtracts each relocated edge once.
   */
  size_t boundary_overhead_bits() const {
    if (layer_count() < 2) {
      return 0;
    }
    const auto m = memory_breakdown();
    const size_t extra_directory = m.directory_bytes - 2 * sizeof(Start);
    return multiply(
               add(add(m.connector_bytes, m.boundary_bytes),
                   add(add(m.connector_index_bytes, m.boundary_index_bytes),
                       extra_directory)),
               8) -
           cross_edges_;
  }
  /** @brief Size of one caller-owned NodeView, excluding scalar call frames. */
  static constexpr size_t navigation_scratch_bytes() {
    return sizeof(NodeView);
  }

 private:
  using Rmm = RmMBTree<>;
  using Rank = RankSelectSupport<>;
  struct Start {
    size_t forest = 0, connector = 0, boundary = 0, nodes = 0;
  };
  struct Record {
    size_t depth, degree;
  };
  struct Counts {
    size_t nodes = 0, edges = 0, boundaries = 0, cross_edges = 0;
  };
  size_t size_ = 0;
  size_t band_height_ = 8;
  size_t depth_count_ = 0;
  double budget_ = 0.05;
  PortalMode mode_ = PortalMode::Dense;
  unsigned width_ = 1;
  size_t boundary_count_ = 0;
  size_t cross_edges_ = 0;
  std::vector<Start> directory_;
  std::vector<uint64_t> forest_, connectors_, boundaries_;
  Rmm forest_index_;
  std::optional<Rank> connector_index_, boundary_index_;

  static size_t add(size_t a, size_t b) {
    if (b > std::numeric_limits<size_t>::max() - a) {
      throw std::length_error("Layered DFUDS size overflow");
    }
    return a + b;
  }
  static size_t multiply(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max() / b) {
      throw std::length_error("Layered DFUDS size overflow");
    }
    return a * b;
  }
  static size_t ceil_div(size_t n, size_t d) { return n / d + (n % d != 0); }
  template <class T>
  static size_t bytes(const std::vector<T>& values) {
    return multiply(values.capacity(), sizeof(T));
  }
  static size_t owned(const std::optional<Rank>& index) {
    return index ? index->memory_usage_bytes() - sizeof(Rank) : 0;
  }
  bool fits(size_t bits) const {
    return static_cast<long double>(bits) <=
           static_cast<long double>(budget_) * size_;
  }
  static size_t rank_bytes(size_t bits, size_t selected_population) {
    // RankSelectSupport's three separately cache-line-rounded allocations.
    const size_t supers =
        multiply(ceil_div(add(ceil_div(bits, 65536), 1), 8), 64);
    const size_t basics = multiply(ceil_div(bits, 16384), 64);
    const size_t samples = add(1, selected_population / 16384);
    return add(add(supers, basics), multiply(ceil_div(samples, 8), 64));
  }
  size_t planned_overhead(size_t bitmap,
                          size_t boundaries,
                          size_t edges,
                          size_t cuts) const {
    const size_t portal_bits = mode_ == PortalMode::Dense
                                   ? bitmap
                                   : multiply(boundaries, id_width(bitmap));
    size_t allocation =
        multiply(add(ceil_div(portal_bits, 64), ceil_div(edges, 64)), 8);
    if (mode_ == PortalMode::Dense) {
      allocation = add(allocation, rank_bytes(bitmap, boundaries));
    }
    allocation = add(allocation, rank_bytes(edges, boundaries));
    allocation = add(allocation, multiply(cuts, sizeof(Start)));
    return multiply(allocation, 8) - edges;
  }
  void build_with_budget(const std::vector<Record>& records,
                         const std::vector<Counts>& bands) {
    size_t cuts = 0;
    size_t bitmap = 0, boundaries = 0, cross_edges = 0;
    for (size_t i = 0; budget_ > 0 && i + 1 < bands.size(); ++i) {
      bitmap = add(bitmap, bands[i].nodes);
      boundaries = add(boundaries, bands[i].boundaries);
      cross_edges = add(cross_edges, bands[i].cross_edges);
      if (!fits(planned_overhead(bitmap, boundaries, cross_edges, i + 1))) {
        break;
      }
      ++cuts;
    }
    build(records, bands, cuts);
    if (cuts && !fits(boundary_overhead_bits())) {
      build(records, bands, 0);
    }
  }
  static void set_one(std::vector<uint64_t>& words, size_t p) {
    words[p / 64] |= uint64_t{1} << (p % 64);
  }
  static void append_degree(std::vector<uint64_t>& words,
                            size_t& position,
                            size_t degree) {
    for (size_t i = 0; i < degree; ++i) {
      set_one(words, position++);
    }
    ++position;
  }
  void build(const std::vector<Record>& records,
             const std::vector<Counts>& bands,
             size_t cuts) {
    // Fresh vectors, not resize of old buffers: fallback must release capacity.
    directory_ = std::vector<Start>(cuts + 2);
    boundary_count_ = cross_edges_ = 0;
    for (size_t i = 0; i < bands.size(); ++i) {
      auto& end = directory_[std::min(i, cuts) + 1];
      end.nodes = add(end.nodes, bands[i].nodes);
      end.forest = add(end.forest, add(bands[i].nodes, bands[i].edges));
      if (i < cuts) {
        end.forest -= bands[i].cross_edges;
        end.connector = bands[i].cross_edges;
        end.boundary = bands[i].nodes;
        boundary_count_ = add(boundary_count_, bands[i].boundaries);
        cross_edges_ = add(cross_edges_, bands[i].cross_edges);
      }
    }
    size_t covered = 0;
    for (size_t i = 0; i < cuts; ++i) {
      covered += bands[i].nodes;
    }
    width_ = id_width(covered);
    for (size_t i = 1; i < directory_.size(); ++i) {
      const auto prev = directory_[i - 1];
      auto& current = directory_[i];
      if (mode_ == PortalMode::Sparse && i <= cuts) {
        current.boundary = multiply(bands[i - 1].boundaries, width_);
      }
      current.forest = add(current.forest, prev.forest);
      current.connector = add(current.connector, prev.connector);
      current.boundary = add(current.boundary, prev.boundary);
      current.nodes = add(current.nodes, prev.nodes);
    }
    forest_ = std::vector<uint64_t>(ceil_div(forest_bits(), 64));
    connectors_ = std::vector<uint64_t>(ceil_div(connector_bits(), 64));
    boundaries_ = std::vector<uint64_t>(ceil_div(boundary_bits(), 64));
    auto cursor = directory_;
    for (const auto& record : records) {
      const size_t band = std::min(record.depth / band_height_, cuts);
      auto& out = cursor[band];
      const bool boundary = record.degree && band < cuts &&
                            record.depth % band_height_ == band_height_ - 1;
      append_degree(forest_, out.forest, boundary ? 0 : record.degree);
      if (boundary) {
        if (mode_ == PortalMode::Dense) {
          set_one(boundaries_, out.nodes);
        } else {
          const size_t word = out.boundary / 64, shift = out.boundary % 64;
          const uint64_t id = out.nodes;
          boundaries_[word] |= id << shift;
          if (shift + width_ > 64) {
            boundaries_[word + 1] |= id >> (64 - shift);
          }
          out.boundary += width_;
        }
        append_degree(connectors_, out.connector, record.degree - 1);
      }
      ++out.nodes;
    }
    build_indexes();
  }
  void build_indexes() {
    forest_index_ =
        size_ ? Rmm(forest_, forest_bits(), Rmm::SelectSupport::kSelect0,
                    forest_bits() - size_)
              : Rmm();
    connector_index_.reset();
    boundary_index_.reset();
    if (layer_count() > 1) {
      connector_index_.emplace(std::span<const uint64_t>(connectors_),
                               connector_bits(), Rank::SelectSupport::kSelect0,
                               cross_edges_ - boundary_count_);
      if (mode_ == PortalMode::Dense) {
        boundary_index_.emplace(std::span<const uint64_t>(boundaries_),
                                boundary_bits(), Rank::SelectSupport::kSelect1,
                                boundary_count_);
      }
    }
  }
  void check(Node node) const {
    if (node.layer >= layer_count() ||
        node.local >=
            directory_[node.layer + 1].nodes - directory_[node.layer].nodes) {
      throw std::out_of_range("Layered DFUDS node");
    }
  }
  static unsigned id_width(size_t covered) {
    return covered > 1 ? std::bit_width(covered - 1) : 1;
  }
  size_t portal_id(size_t ordinal) const {
    if (mode_ == PortalMode::Dense) {
      return boundary_index_->select(ordinal + 1);
    }
    const size_t bit = ordinal * width_, shift = bit % 64;
    uint64_t id = boundaries_[bit / 64] >> shift;
    if (shift + width_ > 64) {
      id |= boundaries_[bit / 64 + 1] << (64 - shift);
    }
    return id & ((uint64_t{1} << width_) - 1);
  }
  size_t portal_rank(size_t number) const {
    if (mode_ == PortalMode::Dense) {
      return boundary_index_->rank(number);
    }
    size_t lo = 0, hi = boundary_count_;
    while (lo < hi) {
      const size_t mid = lo + (hi - lo) / 2;
      if (portal_id(mid) < number) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return lo;
  }
  bool boundary(Node node, size_t number) const {
    if (node.layer + 1 == layer_count()) {
      return false;
    }
    if (mode_ == PortalMode::Dense) {
      return (*boundary_index_)[number];
    }
    const size_t p = portal_rank(number);
    return p < boundary_count_ && portal_id(p) == number;
  }
  size_t position(Node node, size_t number) const {
    return node.local ? forest_index_.select0(number) + 1
                      : directory_[node.layer].forest;
  }
  int64_t excess(size_t layer, size_t p) const {
    const auto& start = directory_[layer];
    const size_t zeros = forest_index_.rank0(p) - start.nodes;
    return static_cast<int64_t>(p - start.forest) -
           static_cast<int64_t>(2 * zeros);
  }
  size_t incoming(Node node, size_t p) const {
    const size_t begin = directory_[node.layer].forest;
    if (p == begin) {
      return Rmm::npos;
    }
    if (p - begin >= 2 && forest_index_.bit(p - 2)) {
      return p - 2;
    }
    // The nearest equal-excess prefix is the incoming open boundary. A forest
    // root has no match in this layer; begin itself is a valid nonroot match.
    const size_t open = forest_index_.bwdsearch(p, 0);
    return open == Rmm::npos || open < begin ? Rmm::npos : open;
  }
  size_t connector_start(size_t one_based_zero) const {
    return one_based_zero == 1
               ? 0
               : connector_index_->select0(one_based_zero - 1) + 1;
  }
  size_t root_connector(size_t layer, size_t p, size_t number) const {
    const auto& start = directory_[layer];
    const size_t ordinal = 2 * (number - start.nodes) - (p - start.forest);
    return directory_[layer - 1].connector + ordinal;
  }
  size_t forest_root_position(size_t layer, size_t ordinal) const {
    const size_t begin = directory_[layer].forest;
    if (!ordinal) {
      return begin;
    }
    if (ordinal <= static_cast<size_t>(std::numeric_limits<int>::max())) {
      return forest_index_.fwdsearch(begin, -static_cast<int>(ordinal)) + 1;
    }
    // Locate the first prefix minimum at -ordinal without an int-sized delta.
    size_t lo = begin, hi = directory_[layer + 1].forest - 1;
    while (lo < hi) {
      const size_t mid = lo + (hi - lo) / 2;
      const size_t at = forest_index_.range_min_query_pos(begin, mid);
      if (excess(layer, at + 1) <= -static_cast<int64_t>(ordinal)) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo + 1;
  }
};

}  // namespace pixie::tree::experimental
