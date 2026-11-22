#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/construction.h>
#include <pixie/tree/detail/forest256_cache.h>
#include <pixie/tree/detail/forest_cover_node_view.h>
#include <pixie/tree/root_navigation_index.h>
#include <pixie/tree/tree_part.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::tree {
template <class Tree>
class TreeTraversalContext;
class ForestCoverTree;

template <>
struct TreeTraits<ForestCoverTree> {
  using Node = detail::ForestCoverNode;
  using NodeView = detail::ForestCoverNodeView<ForestCoverTree>;
};
/** @brief Store an ordered tree as disjoint LOUDS forests of at most 256 nodes.
 * @details The threshold is 128. Units and their nodes are each numbered in
 * breadth-first order. Global IDs concatenate units; they need not match the
 * source IDs or preorder. No source or address permutation survives
 * construction. Copies own their data. Readers borrow a live owner, which
 * cannot be assigned or moved. Separate readers may query the same tree
 * concurrently.
 *
 * On the 13-shape x86-64 space screen (including both combs, attachment
 * gadgets, random-recursive and Catalan trees), maximum retained tree storage
 * was 308,964 bytes at 1,048,576 nodes and 1,243,364 at 4,194,304 nodes.
 * These are measured capacities including the object, not a universal bound;
 * they exclude the caller's cache, construction scratch and allocator overhead.
 *
 * 2026-11-18, ordinal/upper/inspected_root_path/persistent, CPU 0, GCC 13.3
 * -O3 -march=native, 1,048,576 nodes, seed 42, 5 interleaved repetitions of
 * 32,768 paths. Median CPU ns/path; tree and index construction excluded:
 *
 * Shape              16 slots   32 slots   16 slots + 11,264-byte prefix
 * Random-recursive      154.1      149.4                            58.8
 * Binary Catalan         90.7       91.2                            26.1
 * Balanced              466.9      457.4                           393.0
 */
class ForestCoverTree : public TreeBase<ForestCoverTree> {
  friend class TreeBase<ForestCoverTree>;
  friend class TreeTraversalContext<ForestCoverTree>;
  friend class detail::ForestCoverNodeView<ForestCoverTree>;
  using U = uint32_t;
  using V = std::vector<U>;
  static constexpr U missing = UINT32_MAX;
  template <class T>
  static size_t bytes(const std::vector<T>& v) {
    return v.capacity() * sizeof(T);
  }
  struct Prefix {
    std::vector<uint8_t> delta;
    std::vector<uint64_t> samples;
    void build(const V& v) {
      delta.resize(v.size());
      samples.reserve(v.size() / 16 + 1);
      uint64_t sum = 0;
      for (size_t i = 0; i < v.size(); ++i) {
        if (i % 16 == 0) {
          samples.push_back(sum);
        }
        if (v[i] == 0 || v[i] > 256 || sum + v[i] >= missing) {
          throw std::logic_error("ordinal forest prefix out of range");
        }
        delta[i] = v[i] - 1;
        sum += v[i];
      }
      if (v.size() % 16 == 0) {
        samples.push_back(sum);
      }
    }
    U at(size_t i) const {
      uint64_t sum = samples[i / 16];
      for (size_t j = i / 16 * 16; j < i; ++j) {
        sum += 1u + delta[j];
      }
      return sum;
    }
    U find(U x) const {
      size_t s = std::upper_bound(samples.begin(), samples.end(), x) -
                 samples.begin() - 1,
             i = s * 16;
      uint64_t sum = samples[s];
      while (sum + 1 + delta[i] <= x) {
        sum += 1u + delta[i++];
      }
      return i;
    }
    size_t memory() const { return bytes(delta) + bytes(samples); }
  } sizes_, roots_;
  std::vector<uint64_t> stream_, keys_, ends_;
  V rank_;
  U size_ = 0, portals_ = 0, width_ = 0;
  U key(U p) const {
    if (!width_) {
      return 0;
    }
    const uint64_t b = uint64_t(p) * width_;
    unsigned shift = b % 64;
    uint64_t v = keys_[b / 64] >> shift;
    if (shift + width_ > 64) {
      v |= keys_[b / 64 + 1] << (64 - shift);
    }
    return v & ((1ull << width_) - 1);
  }
  U lower_portal(U id) const {
    U l = 0, r = portals_;
    while (l < r) {
      U m = l + (r - l) / 2;
      if (key(m) < id) {
        l = m + 1;
      } else {
        r = m;
      }
    }
    return l;
  }
  U portal(U id) const {
    U l = lower_portal(id);
    return l < portals_ && key(l) == id ? l : missing;
  }
  // One terminator per positive group of CHILD UNITS. Ordinal child counts
  // instead come from differences of cumulative root counts.
  U group(U u) const {
    U b = u - 1;
    return rank_[b / 64] +
           std::popcount(ends_[b / 64] & ((1ull << (b % 64)) - 1));
  }
  U end(U p) const {
    size_t i =
        std::upper_bound(rank_.begin(), rank_.end(), p) - rank_.begin() - 1;
    uint64_t w = ends_[i];
    for (U j = rank_[i]; j < p; ++j) {
      w &= w - 1;
    }
    return i * 64 + std::countr_zero(w) + 2;
  }
  U begin(U p) const { return p ? end(p - 1) : 1; }

 public:
  /// Deep copy; neither the original tree nor the construction source is
  /// borrowed.
  ForestCoverTree(const ForestCoverTree&) = default;
  ForestCoverTree(ForestCoverTree&&) = delete;
  ForestCoverTree& operator=(const ForestCoverTree&) = delete;
  /// Build from a valid ordered source with arbitrary stable handles in O(n)
  /// source operations and temporary space. A size >= UINT32_MAX or above
  /// SIZE_MAX/32 throws length_error before allocation/traversal. A mismatch
  /// between declared size and traversed node/edge counts throws
  /// invalid_argument.
  template <pixie::tree::detail::OrderedTreeSource S>
  explicit ForestCoverTree(const S& source) {
    const size_t n = source.size();
    if (n >= missing || n > SIZE_MAX / 32) {
      throw std::length_error("ordinal forest size");
    }
    size_ = n;
    if (!n) {
      return;
    }
    struct Tmp {
      U parent = missing, first = missing, last = missing, next = missing,
        owner = missing, res = 1;
    };
    std::vector<Tmp> a;
    a.reserve(n);
    V stack;
    size_t edges = 0;
    pixie::tree::detail::generic_construction_dfs(
        source,
        [&](auto, size_t d) {
          if (a.size() == n || d > n - 1 - edges) {
            throw std::invalid_argument("source size/degree");
          }
          edges += d;
          U v = a.size();
          a.emplace_back();
          if (!stack.empty()) {
            U p = stack.back();
            a[v].parent = p;
            if (a[p].first == missing) {
              a[p].first = v;
            } else {
              a[a[p].last].next = v;
            }
            a[p].last = v;
          }
          stack.push_back(v);
        },
        [&](auto) { stack.pop_back(); });
    if (a.size() != n || edges != n - 1) {
      throw std::invalid_argument("source size/topology");
    }
    // Reverse traversal order is postorder for residual sizes. A portal sends
    // all its children to consecutive buckets and returns one node to its
    // parent.
    U units = 1;
    for (size_t i = n; i-- > 0;) {
      U sum = 0;
      for (U c = a[i].first; c != missing; c = a[c].next) {
        sum += a[c].res;
      }
      if (sum < 128) {
        a[i].res += sum;
      } else {
        U fill = 256;
        for (U c = a[i].first; c != missing; c = a[c].next) {
          if (fill + a[c].res > 256) {
            ++units;
            fill = 0;
          }
          a[c].owner = units - 1;
          fill += a[c].res;
        }
      }
    }
    a[0].owner = 0;
    for (U i = 1; i < n; ++i) {
      if (a[i].owner == missing) {
        a[i].owner = a[a[i].parent].owner;
      }
    }
    // Number child units in encounter order while encoding their parent unit in
    // forest-wide BFS. A key therefore owns one contiguous range of whole
    // units.
    std::vector<V> buckets(units);
    buckets[0].push_back(0);
    V order{0}, sizes, roots, keys, ends;
    U base = 0;
    stream_.resize((uint64_t(n) * 2 + 63) / 64);
    for (size_t j = 0; j < order.size(); ++j) {
      U id = order[j];
      V nodes = std::move(buckets[id]);
      roots.push_back(nodes.size());
      std::vector<uint16_t> degrees;
      for (size_t l = 0; l < nodes.size(); ++l) {
        U degree = 0, last = missing;
        for (U c = a[nodes[l]].first; c != missing; c = a[c].next) {
          U target = a[c].owner;
          if (target == id) {
            ++degree;
            nodes.push_back(c);
          } else {
            if (last != target) {
              order.push_back(target);
              last = target;
            }
            buckets[target].push_back(c);
          }
        }
        if (last != missing) {
          keys.push_back(base + l);
          ends.push_back(order.size());
        }
        degrees.push_back(degree);
      }
      if (nodes.size() > size_ - base) {
        throw std::logic_error("ordinal forest exceeds source size");
      }
      auto encoded =
          detail::forest256::encode_louds_forest256(roots.back(), degrees);
      const size_t bits = 2 * nodes.size();
      const size_t start = 2ull * base;
      const size_t shift = start % 64;
      for (size_t b = 0; b < bits; b += 64) {
        uint64_t word = encoded.words[b / 64];
        const size_t remaining = bits - b;
        if (remaining < 64) {
          // A short forest's encoder padding is not part of the stream.
          word &= (uint64_t{1} << remaining) - 1;
        }
        const size_t dst = (start + b) / 64;
        stream_[dst] |= word << shift;
        if (shift && remaining > 64 - shift) {
          stream_[dst + 1] |= word >> (64 - shift);
        }
      }
      sizes.push_back(nodes.size());
      base += nodes.size();
    }
    if (base != size_ || order.size() != units || units > 1 + n / 128) {
      throw std::logic_error("ordinal cover");
    }
    sizes_.build(sizes);
    roots_.build(roots);
    portals_ = keys.size();
    width_ = keys.empty() ? 0 : std::bit_width(keys.back());
    keys_.assign((uint64_t(portals_) * width_ + 63) / 64, 0);
    for (size_t i = 0; width_ && i < keys.size(); ++i) {
      size_t b = i * width_;
      keys_[b / 64] |= uint64_t(keys[i]) << (b % 64);
      if (b % 64 + width_ > 64) {
        keys_[b / 64 + 1] |= keys[i] >> (64 - b % 64);
      }
    }
    ends_.assign((size_t(units) + 62) / 64, 0);
    for (U e : ends) {
      U b = e - 2;
      ends_[b / 64] |= 1ull << (b % 64);
    }
    rank_.resize(ends_.size() + 1);
    for (size_t i = 0; i < ends_.size(); ++i) {
      rank_[i + 1] = rank_[i] + std::popcount(ends_[i]);
    }
  }
  /// Number of disjoint bounded forests.
  size_t unit_count() const { return sizes_.delta.size(); }

  /// Copy one unit's local BFS rows. Portal child lists omit external
  /// children and are incomplete. Throws out_of_range for an invalid unit.
  TreePart inspect_unit(size_t unit) const {
    if (unit >= unit_count()) {
      throw std::out_of_range("forest unit");
    }
    const U base = sizes_.at(unit), count = 1u + sizes_.delta[unit];
    detail::forest256::LevelCsr256 decoded;
    unpack_dense(detail::forest256::Codec::louds, stream_, 2ull * base, count,
                 decoded);
    std::vector<TreePart::Node> rows(count);
    for (U i = 0; i < count; ++i) {
      auto& row = rows[i];
      row.source_id = base + i;
      if (const U degree = decoded.degree(i)) {
        const U first = decoded.offset_m1[i] + 1;
        for (U k = 0; k < degree; ++k) {
          row.children.push_back(first + k);
        }
      } else if (portal(base + i) != missing) {
        row.children_complete = false;
      }
    }
    return {std::move(rows), 1u + roots_.delta[unit]};
  }

 private:
  size_t size_impl() const { return size_; }
  /// Root handle, or out_of_range on empty input.
  Node root_impl() const {
    if (!size_) {
      throw std::out_of_range("empty ordinal forest");
    }
    return {};
  }
  // Cross-unit parents precede child units; local BFS puts parents first.
  auto path_key_impl(Node h) const { return std::pair{h.unit, h.local}; }
  /// Zero-based stable ID. Invalid numerical coordinates throw out_of_range.
  size_t number_impl(Node h) const {
    check(h);
    return sizes_.at(h.unit) + h.local;
  }
  /// Inverse of number. Out-of-range ID throws out_of_range.
  Node find_by_number_impl(size_t id) const {
    U u = sizes_.find(id);
    return {u, U(id) - sizes_.at(u)};
  }
  /// Owning object plus vector capacity bytes, excluding allocator bookkeeping.
  size_t memory_usage_bytes_impl() const {
    return sizeof(*this) + sizes_.memory() + roots_.memory() + bytes(stream_) +
           bytes(keys_) + bytes(ends_) + bytes(rank_);
  }

 private:
  bool contains_impl(Node h) const {
    return h.unit < unit_count() && h.local <= sizes_.delta[h.unit];
  }
  void check(Node h) const {
    if (!contains_impl(h)) {
      throw std::out_of_range("ordinal forest node");
    }
  }

  detail::forest256::NodeInfo info(
      Node h,
      detail::forest256::Forest256Cache& cache) const {
    check(h);
    return cache.lookup(h.unit, detail::forest256::Codec::louds,
                        2ull * sizes_.at(h.unit), 1u + sizes_.delta[h.unit],
                        h.local);
  }
  detail::forest256::NodeInfo info(Node h) const {
    detail::forest256::Forest256Cache transient(stream_, 0);
    return info(h, transient);
  }
  NodeView navigate_impl(Node h) const { return {*this, h, info(h)}; }
  NodeView cached_view(Node h, detail::forest256::Forest256Cache& cache) const {
    return {*this, h, info(h, cache)};
  }
  bool is_root_impl(Node h) const {
    check(h);
    return h == Node{};
  }
  size_t degree_impl(Node h) const { return navigate_impl(h).degree(); }
  bool is_leaf_impl(Node h) const { return navigate_impl(h).is_leaf(); }
  Node child_impl(Node h, size_t k) const { return navigate_impl(h).child(k); }
  Node first_child_impl(Node h) const { return navigate_impl(h).first_child(); }
  Node parent_impl(Node h) const { return navigate_impl(h).parent(); }
  Node next_sibling_impl(Node h) const {
    return navigate_impl(h).next_sibling();
  }
  bool is_last_child_impl(Node h) const {
    return navigate_impl(h).is_last_child();
  }
  template <class F>
  void for_each_child_impl(Node h, F&& f) const {
    navigate_impl(h).for_each_child(f);
  }
};

template <>
class RootNavigationIndex<ForestCoverTree>;

/// Caller-owned decoder cache borrowing a tree and an optional root index.
/// Keep both original owners alive throughout a walk, even if a callback
/// rebinds this context. Independent contexts may query one tree and index
/// concurrently. A moved-from context has no binding or slots; queries throw
/// logic_error until rebound.
template <>
class TreeTraversalContext<ForestCoverTree> {
  friend class RootNavigationIndex<ForestCoverTree>;
  using Index = RootNavigationIndex<ForestCoverTree>;
  using Node = ForestCoverTree::Node;
  using U = uint32_t;
  static constexpr U missing = UINT32_MAX;
  const ForestCoverTree* t_;
  const Index* index_;
  detail::forest256::Forest256Cache cache_;
  detail::forest256::NodeInfo info(Node h) {
    if (!t_) {
      throw std::logic_error("unbound reader");
    }
    return t_->info(h, cache_);
  }
  ForestCoverTree::NodeView view(Node h) {
    if (!t_) {
      throw std::logic_error("unbound context");
    }
    return t_->cached_view(h, cache_);
  }

  // Every handle in the walk is generated from the original tree. A visitor
  // can evict our slot or rebind the reader; never borrow a decoded row.
  template <bool Path, class Visit, class Inspect, class Choose>
  void walk(const ForestCoverTree* tree,
            Node h,
            Visit& visit,
            Inspect& inspect,
            Choose& choose) {
    std::conditional_t<Path, detail::forest256::NodeInfo,
                       detail::forest256::LevelCsr256>
        local;
    bool entering = true;
    U previous = missing, base = 0, first_portal = 0, first_key = missing,
      second_key = missing;
    for (;;) {
      if constexpr (Path) {
        visit(Node(h));
      } else {
        visit(Node(h), bool(entering));
      }
      if constexpr (Path) {
        base = tree->sizes_.at(h.unit);
        local = t_ == tree ? info(h) : tree->info(h);
      } else if (previous != h.unit) {
        // DFS usually touches many local leaves before changing units. No
        // packed-key search is needed for IDs below the second portal key.
        previous = h.unit;
        base = tree->sizes_.at(h.unit);
        const auto count = 1u + tree->sizes_.delta[h.unit];
        if (tree == t_) {
          cache_.snapshot(h.unit, detail::forest256::Codec::louds, 2ull * base,
                          count, local);
        } else {
          unpack_dense(detail::forest256::Codec::louds, tree->stream_,
                       2ull * base, count, local);
        }
        first_portal = tree->lower_portal(base);
        first_key =
            first_portal < tree->portals_ ? tree->key(first_portal) : missing;
        second_key = first_portal + 1 < tree->portals_
                         ? tree->key(first_portal + 1)
                         : missing;
      }
      if (entering) {
        U first, count, unit = h.unit, end = unit + 1;
        if constexpr (Path) {
          first = local.first;
          count = local.degree;
        } else {
          first = local.offset_m1[h.local] + 1;
          count = local.degree(h.local);
        }
        if (!count) {
          U id = base + h.local;
          U p = Path || id >= second_key ? tree->portal(id)
                : id == first_key        ? first_portal
                                         : missing;
          if (p != missing) {
            unit = tree->begin(p);
            end = tree->end(p);
            first = 0;
            count = 1u + tree->roots_.delta[unit];
          }
        }
        if (count) {
          Node selected{unit, first};
          size_t k = 0;
          U skipped = 0;
          if constexpr (Path) {
            U degree = end == unit + 1
                           ? count
                           : tree->roots_.at(end) - tree->roots_.at(unit);
            k = choose(size_t(degree));
            if (k >= degree) {
              throw std::out_of_range("child index");
            }
          }
          // A local sibling run is contiguous. External children are runs
          // of roots, each starting at local zero in the next unit.
          for (;;) {
            if constexpr (Path) {
              if (k >= skipped && k < skipped + count) {
                selected = {unit, first + U(k - skipped)};
              }
              skipped += count;
            }
            for (U i = 0; i < count; ++i) {
              inspect(Node{unit, first + i});
            }
            if (++unit == end) {
              break;
            }
            first = 0;
            count = 1u + tree->roots_.delta[unit];
          }
          h = selected;
          continue;
        }
      }
      if constexpr (Path) {
        return;
      } else {
        if (h == Node{}) {
          return;
        }
        entering = true;
        if (!local.is_last_child(h.local)) {
          ++h.local;
        } else if (!local.is_root(h.local)) {
          h.local = local.parent_of(h.local);
          entering = false;
        } else if ((tree->ends_[(h.unit - 1) / 64] >> ((h.unit - 1) % 64)) &
                   1) {
          h = tree->find_by_number(tree->key(tree->group(h.unit)));
          entering = false;
        } else {
          h = {h.unit + 1, 0};
        }
      }
    }
  }

 public:
  static constexpr size_t slot_bytes() {
    return detail::forest256::Forest256Cache::slot_bytes();
  }
  static constexpr size_t default_cache_bytes =
      16 * detail::forest256::Forest256Cache::slot_bytes();
  /// Borrow tree and optional matching index. Only whole decoder slots
  /// fit in cache_bytes; zero uses transient decoding without allocation.
  explicit TreeTraversalContext(const ForestCoverTree& tree,
                                size_t cache_bytes = default_cache_bytes,
                                const Index* index = nullptr);
  TreeTraversalContext(const ForestCoverTree&&,
                       size_t = default_cache_bytes,
                       const Index* = nullptr) = delete;
  /// Borrow the same tree and allocate independent, initially empty slots.
  TreeTraversalContext(const TreeTraversalContext&) = default;
  /// Copy the binding and byte allowance; a failed allocation leaves it
  /// unchanged.
  TreeTraversalContext& operator=(const TreeTraversalContext& other) {
    if (this != &other) {
      cache_ = other.cache_;
      t_ = other.t_;
      index_ = other.index_;
    }
    return *this;
  }
  /// Transfer the binding and slots without allocation.
  TreeTraversalContext(TreeTraversalContext&& other) noexcept
      : t_(std::exchange(other.t_, nullptr)),
        index_(std::exchange(other.index_, nullptr)),
        cache_(std::move(other.cache_)) {}
  /// Replace the binding and slots with those of the source reader.
  TreeTraversalContext& operator=(TreeTraversalContext&& other) noexcept {
    if (this != &other) {
      cache_ = std::move(other.cache_);
      t_ = std::exchange(other.t_, nullptr);
      index_ = std::exchange(other.index_, nullptr);
    }
    return *this;
  }
  /// Rebind to a live owner, invalidating all tags and preserving the budget.
  void bind(const ForestCoverTree& t, const Index* index = nullptr);
  void bind(const ForestCoverTree&&, const Index* = nullptr) = delete;
  /// Retained heap bytes for decoder slots, excluding the context object.
  size_t allocated_bytes() const { return cache_.allocated_bytes(); }
  /// Context object and retained cache slots, excluding the borrowed index.
  size_t memory_usage_bytes() const {
    return sizeof(*this) + allocated_bytes();
  }

 private:
  Node root() const {
    if (!t_) {
      throw std::logic_error("unbound reader");
    }
    return t_->root();
  }

 public:
  /// Visit (node,true), inspect every child in order, then traverse the
  /// children and visit (node,false) once on return from a nonleaf. The
  /// walk uses one 512-byte temporary decoded unit, independent of depth,
  /// with no extra heap slots or stack of ancestors.
  /// Empty trees throw out_of_range. Callback exceptions propagate. As for
  /// for_each_child, keep the original tree alive across callbacks/rebinding.
  template <class Visit, class Inspect>
  void traverse_inspected(Visit&& visit, Inspect&& inspect) {
    auto first = [](size_t) { return 0u; };
    walk<false>(t_, root(), visit, inspect, first);
  }
  /// Visit each node on one root-to-leaf path. At each nonleaf ask
  /// choose(degree) for a zero-based child, inspect *all* children in order,
  /// then descend to the chosen handle. Invalid choices throw out_of_range
  /// before inspection. Same lifetime and exception rules as
  /// traverse_inspected. Returns the number of visits handled by the index.
  template <class Visit, class Inspect, class Choose>
  size_t inspect_path(Visit&& visit, Inspect&& inspect, Choose&& choose);
  /// Cached navigation has the same semantics as the bound tree's facade.
  /// Invalid numerical coordinates throw out_of_range.
  size_t degree(Node h) { return view(h).degree(); }
  bool is_leaf(Node h) { return view(h).is_leaf(); }
  Node child(Node h, size_t k) { return view(h).child(k); }
  Node first_child(Node h) { return view(h).first_child(); }
  Node parent(Node h) { return view(h).parent(); }
  /// Same inclusive path and buffer contract as the tree facade. Parent
  /// steps share this reader's decoder slots, including across queries.
  void tree_path(Node u,
                 Node v,
                 std::vector<Node>& out,
                 std::vector<Node>& scratch) {
    if (!t_) {
      throw std::logic_error("unbound reader");
    }
    detail::collect_tree_path(
        u, v, out, scratch, [this](Node n) { return t_->path_key_impl(n); },
        [this](Node n) { return parent(n); });
  }
  bool is_last_child(Node h) { return view(h).is_last_child(); }
  Node next_sibling(Node h) { return view(h).next_sibling(); }
  template <class F>
  void for_each_child(Node h, F&& f) {
    view(h).for_each_child(f);
  }
};

/// Optional immutable breadth-first prefix for inspected paths. Keep both
/// this object and its original immutable tree alive and unmoved during use.
/// No tree-wide lookup map or reader state. DFS visits the whole tree and
/// does not use this prefix.
template <>
class RootNavigationIndex<ForestCoverTree> {
  friend class TreeTraversalContext<ForestCoverTree>;
  using Node = ForestCoverTree::Node;
  using U = uint32_t;
  static constexpr U missing = UINT32_MAX;
  struct Entry {
    Node node;
    U first = missing, degree = 0;
  };
  const ForestCoverTree& tree_;
  std::vector<Entry> entries_;

 public:
  /// Budget includes object and entry capacity when active; a disabled
  /// object still occupies sizeof(RootNavigationIndex). Build cache is
  /// temporary.
  explicit RootNavigationIndex(const ForestCoverTree& tree,
                               size_t budget = kDefaultTraversalBytes)
      : tree_(tree) {
    const size_t cap =
        std::min({size_t(32768), tree.size(),
                  (budget - std::min(budget, sizeof(*this))) / sizeof(Entry)});
    if (!cap) {
      return;
    }
    entries_.reserve(cap);
    if (entries_.capacity() > cap) {
      std::vector<Entry>().swap(entries_);
      return;
    }
    TreeTraversalContext<ForestCoverTree> r(tree);
    entries_.push_back({tree.root()});
    for (size_t i = 0; i < entries_.size(); ++i) {
      auto& e = entries_[i];
      e.degree = r.degree(e.node);
      if (e.degree && e.degree <= cap - entries_.size()) {
        e.first = static_cast<U>(entries_.size());
        r.for_each_child(e.node, [&](Node c) { entries_.push_back({c}); });
      }
    }
  }
  RootNavigationIndex(const ForestCoverTree&&,
                      size_t = kDefaultTraversalBytes) = delete;
  RootNavigationIndex(const RootNavigationIndex&) = delete;
  const ForestCoverTree& tree() const { return tree_; }
  /// Number of retained handles, including unexpanded frontier nodes.
  size_t size() const { return entries_.size(); }
  /// Active retained bytes. Excludes allocator bookkeeping and build scratch.
  size_t budget_usage_bytes() const {
    return size() ? memory_usage_bytes() : 0;
  }
  size_t memory_usage_bytes() const {
    return sizeof(*this) + entries_.capacity() * sizeof(Entry);
  }
  /// Copy cached BFS rows, leaving unexpanded child lists incomplete.
  TreePart inspect_part() const {
    std::vector<TreePart::Node> rows(entries_.size());
    for (size_t i = 0; i < entries_.size(); ++i) {
      const auto& e = entries_[i];
      auto& row = rows[i];
      row.source_id = tree_.number(e.node);
      row.children_complete = !e.degree || e.first != missing;
      if (e.first != missing) {
        for (U j = 0; j < e.degree; ++j) {
          row.children.push_back(e.first + j);
        }
      }
    }
    return {std::move(rows), entries_.empty() ? 0u : 1u};
  }

 private:
  // Enumerate cached children first; resume the exact walk at its frontier.
  template <class Visit, class Inspect, class Choose>
  size_t inspect_path(TreeTraversalContext<ForestCoverTree>& r,
                      Visit&& visit,
                      Inspect&& inspect,
                      Choose&& choose) const {
    size_t hits = 0;
    Node frontier = tree_.root();
    if (!entries_.empty()) {
      U i = 0;
      for (;;) {
        const auto e = entries_[i];
        if (e.first == missing && e.degree) {
          frontier = e.node;
          break;
        }
        ++hits;
        visit(Node(e.node));
        if (!e.degree) {
          return hits;
        }
        size_t k = choose(size_t(e.degree));
        if (k >= e.degree) {
          throw std::out_of_range("child index");
        }
        for (U j = 0; j < e.degree; ++j) {
          inspect(Node(entries_[e.first + j].node));
        }
        i = e.first + U(k);
      }
    }
    r.walk<true>(&tree_, frontier, visit, inspect, choose);
    return hits;
  }
};

inline TreeTraversalContext<ForestCoverTree>::TreeTraversalContext(
    const ForestCoverTree& tree,
    size_t cache_bytes,
    const Index* index)
    : t_(&tree), index_(index), cache_(tree.stream_, cache_bytes) {
  if (index && &index->tree() != &tree) {
    throw std::invalid_argument("forest index tree mismatch");
  }
}
inline void TreeTraversalContext<ForestCoverTree>::bind(
    const ForestCoverTree& tree,
    const Index* index) {
  if (index && &index->tree() != &tree) {
    throw std::invalid_argument("forest index tree mismatch");
  }
  cache_.bind(tree.stream_);
  t_ = &tree;
  index_ = index;
}
template <class Visit, class Inspect, class Choose>
size_t TreeTraversalContext<ForestCoverTree>::inspect_path(Visit&& visit,
                                                           Inspect&& inspect,
                                                           Choose&& choose) {
  if (index_) {
    return index_->inspect_path(*this, visit, inspect, choose);
  }
  walk<true>(t_, root(), visit, inspect, choose);
  return 0;
}
}  // namespace pixie::tree
