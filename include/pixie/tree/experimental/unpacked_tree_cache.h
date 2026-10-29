#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>

namespace pixie::detail {
// Unsupported encoding access: decoding stays with the representation, never
// in the benchmark. The caller retains the immutable source for the tree.
struct NavigationChunkAccess {
  template <class B, class C>
  static auto words(const BPTree<B, C>& tree) {
    return std::span(tree.words_, (tree.num_bits_ + 63) / 64);
  }
  template <class B, class C>
  static auto words(const DFUDSTree<B, C>& tree) {
    return std::span(tree.words_, (tree.num_bits_ + 63) / 64);
  }
  template <class B, class C>
  static auto words(const LoudsTree<B, C>& tree) {
    return std::span(tree.words_, (tree.bv.size() + 63) / 64);
  }
};
}  // namespace pixie::detail

namespace pixie::experimental::detail {
template <class Tree>
struct ChunkEncoding;
template <class B, class C>
struct ChunkEncoding<BPTree<B, C>> {
  static constexpr int kind = 0;
};
template <class B, class C>
struct ChunkEncoding<DFUDSTree<B, C>> {
  static constexpr int kind = 1;
};
template <class B, class C>
struct ChunkEncoding<LoudsTree<B, C>> {
  static constexpr int kind = 2;
};

struct ChunkHandle {
  size_t number = 0, pos = 0;
  TreeNode node() const { return {number, pos}; }
};
// Selectively retained node result. A complete degree and a partial prefix are
// sufficient to avoid repeated queries without imposing a maximum node degree.
struct DecodedNode {
  size_t degree = 0, count = 0;
  std::array<ChunkHandle, 8> children{};
};

template <class Tree, size_t Bits>
struct NavigationChunk {
  static constexpr int kind = ChunkEncoding<Tree>::kind;
  pixie::detail::tree_unpack::BoundedSubtreeCache<Bits> subtree;
  ChunkHandle anchor;
  size_t degree = 0, count = 0;
  // LOUDS stores a bounded window of child handles, not a preorder subtree.
  std::array<ChunkHandle, kind == 2 ? Bits : 8> children{};
  bool complete_subtree = false;

  bool contains(TreeNode node) const {
    return complete_subtree ? subtree.contains(node)
                            : node.number == anchor.number;
  }
  void load(const Tree& tree, TreeNode node) {
    anchor = {node.number, node.pos};
    count = 0;
    complete_subtree = false;
    const auto words = pixie::detail::NavigationChunkAccess::words(tree);
    if constexpr (kind != 2) {
      complete_subtree = subtree.template load<kind == 0>(
          tree, words, 2 * tree.size() - (kind == 0 ? 0 : 1), node);
      if (complete_subtree) {
        return;
      }
    }
    typename Tree::LocalNodeView local(tree, node);
    degree = local.degree();
    if (!degree) {
      return;
    }
    if constexpr (kind == 2) {
      // One select locates the window. All following handles come from SIMD
      // delimiter extraction over at most Bits bits; long runs end the prefix.
      const auto start = local.first_child();
      const size_t length = std::min(Bits, 2 * tree.size() - 1 - start.pos);
      alignas(64) uint8_t offsets[64];
      for (size_t base = 0; base < length && count < degree; base += 64) {
        const auto n = pixie::detail::tree_unpack::zero_offsets(
            ~pixie::detail::tree_unpack::word_at(words, start.pos + base),
            std::min<size_t>(64, length - base), offsets);
        for (unsigned i = 0; i < n && count < degree; ++i) {
          children[count] = {start.number + count,
                             start.pos + base + offsets[i]};
          ++count;
        }
      }
    } else {
      // A non-fitting subtree retains a bounded child prefix. Its degree may
      // still require ordinary BP navigation; the bit probe itself is bounded.
      auto child = local.first_child();
      while (count < std::min(degree, children.size())) {
        children[count++] = {child.number, child.pos};
        if (count < std::min(degree, children.size())) {
          child = tree.next_sibling(child);
        }
      }
    }
  }
  DecodedNode record(TreeNode node) const {
    DecodedNode out;
    out.degree = complete_subtree ? subtree.degree(node) : degree;
    out.count = std::min(out.degree, out.children.size());
    if (complete_subtree) {
      // At most eight links; no rank/select calls for decoded nodes.
      for (size_t i = 0; i < out.count; ++i) {
        auto child = subtree.child(node, i);
        out.children[i] = {child.number, child.pos};
      }
    } else {
      out.count = std::min(out.count, count);
      std::copy_n(children.begin(), out.count, out.children.begin());
    }
    return out;
  }
};

// Fixed-capacity CLOCK table. Only anchor keys are indexed for chunks; a
// current-chunk hint serves all its interior nodes without a directory lookup.
// No pointer into this table may survive a user callback.
template <class Value, size_t Budget>
class ChunkClock {
  struct Entry {
    size_t key = 0;
    Value value;
    bool recent = false;
  };

 public:
  static constexpr size_t capacity =
      std::bit_floor((Budget - 64) / (sizeof(Entry) + 2 * sizeof(uint32_t)));
  static_assert(capacity > 0);

 private:
  struct Storage {
    explicit Storage(size_t count)
        : entries(std::make_unique<Entry[]>(count)),
          directory(std::make_unique<uint32_t[]>(2 * count)) {}
    std::unique_ptr<Entry[]> entries;
    std::unique_ptr<uint32_t[]> directory;
  };
  std::unique_ptr<Storage> storage_;
  size_t used_ = 0, hand_ = 0;
  size_t capacity_ = 0, mask = 0;
  size_t hash(size_t x) const {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return (x ^ (x >> 31)) & mask;
  }
  void erase(size_t key) {
    auto& d = storage_->directory;
    size_t hole = hash(key);
    while (storage_->entries[d[hole] - 1].key != key) {
      hole = (hole + 1) & mask;
    }
    for (size_t pos = (hole + 1) & mask; d[pos]; pos = (pos + 1) & mask) {
      auto home = hash(storage_->entries[d[pos] - 1].key);
      if (((hole - home) & mask) < ((pos - home) & mask)) {
        d[hole] = d[pos];
        hole = pos;
      }
    }
    d[hole] = 0;
  }

 public:
  explicit ChunkClock(size_t budget = Budget) {
    budget = std::min(budget, Budget);
    constexpr auto overhead = sizeof(ChunkClock) + sizeof(Storage);
    capacity_ =
        budget > overhead
            ? std::min(capacity,
                       std::bit_floor((budget - overhead) /
                                      (sizeof(Entry) + 2 * sizeof(uint32_t))))
            : 0;
    mask = capacity_ ? 2 * capacity_ - 1 : 0;
  }
  bool enabled() const { return capacity_ != 0; }
  Value* find(size_t key) {
    if (!storage_) {
      return nullptr;
    }
    auto& d = storage_->directory;
    for (size_t p = hash(key); d[p]; p = (p + 1) & mask) {
      auto& e = storage_->entries[d[p] - 1];
      if (e.key == key) {
        e.recent = true;
        return &e.value;
      }
    }
    return nullptr;
  }
  template <class Fill>
  Value& insert(size_t key, Fill&& fill) {
    if (auto* existing = find(key)) {
      return *existing;
    }
    if (!storage_) {
      storage_ = std::make_unique<Storage>(capacity_);
    }
    size_t slot;
    if (used_ < capacity_) {
      slot = used_++;
    } else {
      while (storage_->entries[hand_].recent) {
        storage_->entries[hand_].recent = false;
        hand_ = (hand_ + 1) % capacity_;
      }
      slot = hand_;
      hand_ = (hand_ + 1) % capacity_;
      erase(storage_->entries[slot].key);
    }
    auto& e = storage_->entries[slot];
    e.key = key;
    e.recent = true;
    fill(e.value);
    auto p = hash(key);
    while (storage_->directory[p]) {
      p = (p + 1) & mask;
    }
    storage_->directory[p] = slot + 1;
    return e.value;
  }
  size_t bytes() const {
    return sizeof(ChunkClock) +
           (storage_ ? sizeof(Storage) +
                           capacity_ * (sizeof(Entry) + 2 * sizeof(uint32_t))
                     : 0);
  }
};
}  // namespace pixie::experimental::detail

namespace pixie::experimental {
/**
 * @brief Experimental combination of bounded decoding and CLOCK retention.
 * @tparam RetainChunks False retains only requested node records (degree and
 * eight children) from a temporary decoded chunk; true retains entire chunks.
 * @tparam ByteBudget Ceiling for context, retained payload and temporary query
 * buffers, including callback-safe copies. Excludes allocator bookkeeping and
 * caller views. The tree constructor may lower this ceiling; insufficient
 * budgets use local navigation without shared retention.
 * @tparam BitCapacity Maximum source bits probed per decoded subtree/window.
 * @details Implements tree/cache.h. BP/DFUDS decode complete fitting subtrees;
 * LOUDS decodes a child-handle window. Non-fitting subtrees use ordinary
 * navigation with a partial child prefix. CLOCK hits promote recent entries.
 * Chunk lookup indexes anchors, with a current-chunk hint for interior nodes.
 * Sequential navigation remains locally memoized. Degree, indexed children,
 * and enumeration use decoding and shared retention. Copies/moves start empty.
 * Views use stable handles; callbacks may query/evict or throw. The tree and
 * immutable source outlive views. Shared access needs external synchronization.
 * Plain streaming DFS keeps its existing encoding kernel; inspected BP/DFUDS
 * DFS can consume retained complete chunks through inspect_decoded().
 */
template <bool RetainChunks,
          size_t ByteBudget = 512 * 1024,
          size_t BitCapacity = 256>
struct UnpackedTreeCache {
  /** @brief Tree-owned decoding and retention context. */
  template <class Tree>
  class Cache {
    using Chunk = detail::NavigationChunk<Tree, BitCapacity>;
    using Record = detail::DecodedNode;
    using Value = std::conditional_t<RetainChunks, Chunk, Record>;
    static_assert(ByteBudget > sizeof(Chunk) + 256);
    using Store = detail::ChunkClock<Value, ByteBudget - sizeof(Chunk) - 128>;
    mutable Store store_;
    mutable std::unique_ptr<Chunk> scratch_;
    mutable Chunk* current_ = nullptr;
    size_t budget_;

    // Callback-safe children/subtree copies and decoder-local load arrays.
    static constexpr size_t kScratchBytes =
        sizeof(Chunk) + BitCapacity * sizeof(detail::ChunkHandle) +
        BitCapacity * 8 + 256;

    static size_t store_budget(size_t bytes) {
      const size_t overhead = sizeof(Cache) - sizeof(Store) + kScratchBytes +
                              (RetainChunks ? 0 : sizeof(Chunk));
      return bytes > overhead ? bytes - overhead : 0;
    }

    Chunk& chunk(const Tree& tree, TreeNode node) const {
      if constexpr (RetainChunks) {
        if (current_ && current_->contains(node)) {
          store_.find(current_->anchor.number);  // Promote the whole chunk.
          return *current_;
        }
        if (auto* found = store_.find(node.number)) {
          return *(current_ = found);
        }
        current_ =
            &store_.insert(node.number, [&](Chunk& c) { c.load(tree, node); });
        return *current_;
      } else {
        if (!scratch_) {
          scratch_ = std::make_unique<Chunk>();
        }
        if (!current_ || !current_->contains(node)) {
          scratch_->load(tree, node);
        }
        current_ = &*scratch_;
        return *current_;
      }
    }
    Record record(const Tree& tree, TreeNode node) const {
      if constexpr (RetainChunks) {
        return chunk(tree, node).record(node);
      } else {
        if (auto* found = store_.find(node.number)) {
          return *found;
        }
        const auto decoded = chunk(tree, node).record(node);
        return store_.insert(node.number, [&](Record& out) { out = decoded; });
      }
    }

   public:
    /** @brief Construct empty; persistent table allocation is lazy. */
    explicit Cache(size_t budget = ByteBudget)
        : store_(store_budget(std::min(budget, ByteBudget))),
          budget_(std::min(budget, ByteBudget)) {}
    /** @brief Copy/move policy without retaining references to another tree. */
    Cache(const Cache& other) : Cache(other.budget_) {}
    /** @brief Start empty when moved; do not borrow the source context. */
    Cache(Cache&& other) noexcept : Cache(other.budget_) {}
    /** @brief Assignment discards all retained data. */
    Cache& operator=(const Cache& other) {
      budget_ = other.budget_;
      store_ = Store{store_budget(budget_)};
      scratch_.reset();
      current_ = nullptr;
      return *this;
    }
    /** @brief Move assignment discards all retained data. */
    Cache& operator=(Cache&& other) noexcept {
      return operator=(static_cast<const Cache&>(other));
    }
    /** @brief Context, retained allocations and reserved temporary buffers;
     * zero when shared retention is disabled. */
    size_t memory_usage_bytes() const {
      return store_.enabled()
                 ? sizeof(Cache) - sizeof(Store) + kScratchBytes +
                       (RetainChunks ? 0 : sizeof(Chunk)) + store_.bytes()
                 : 0;
    }
    /** @brief Navigation view with the tree/cache.h semantics. */
    class NodeView {
      const Tree* tree_;
      const Cache* owner_;
      typename Tree::LocalNodeView local_;
      friend class Cache;
      NodeView(const Tree& tree, TreeNode node, const Cache& owner)
          : tree_(&tree), owner_(&owner), local_(tree, node) {}

     public:
      /** @brief Current stable node handle. */
      TreeNode node() const { return local_.node(); }
      /** @brief Reset to a valid handle of this tree. */
      void reset(TreeNode node) { local_.reset(node); }
      /** @brief Root/leaf predicates; no shared-cache admission. */
      bool is_root() const { return local_.is_root(); }
      /** @brief Test for a leaf without shared-cache admission. */
      bool is_leaf() const { return local_.is_leaf(); }
      /** @brief Exact child count, decoding/retaining on a miss. */
      size_t degree() const {
        if (!owner_->store_.enabled()) {
          return local_.degree();
        }
        if constexpr (RetainChunks) {
          const auto& c = owner_->chunk(*tree_, node());
          return c.complete_subtree ? c.subtree.degree(node()) : c.degree;
        } else {
          return owner_->record(*tree_, node()).degree;
        }
      }
      /** @brief Child at zero-based i; requires i < degree(). */
      TreeNode child(size_t i) const {
        if (!owner_->store_.enabled()) {
          return local_.child(i);
        }
        if constexpr (RetainChunks) {
          const auto& c = owner_->chunk(*tree_, node());
          if (c.complete_subtree) {
            return c.subtree.child(node(), i);
          }
          if (i < c.count) {
            return c.children[i].node();
          }
        } else {
          const auto r = owner_->record(*tree_, node());
          if (i < r.count) {
            return r.children[i].node();
          }
        }
        return local_.child(i);
      }
      /** @brief Local sequential navigation, with TreeBase preconditions. */
      TreeNode first_child() const { return local_.first_child(); }
      /** @brief Parent, or root itself for root; uses local memoization. */
      TreeNode parent() const { return local_.parent(); }
      /** @brief Next sibling; requires a non-root, non-last child. */
      TreeNode next_sibling() const { return local_.next_sibling(); }
      /** @brief Test for the last sibling; requires a non-root node. */
      bool is_last_child() const { return local_.is_last_child(); }
      /** @brief Enumerate all children in order; callbacks may evict entries.
       * Partial prefixes fall back to complete representation enumeration. */
      template <class Visitor>
      void for_each_child(Visitor&& visit) const {
        if (!owner_->store_.enabled()) {
          local_.for_each_child(std::forward<Visitor>(visit));
          return;
        }
        if constexpr (RetainChunks) {
          const auto& c = owner_->chunk(*tree_, node());
          const auto count =
              c.complete_subtree ? c.subtree.degree(node()) : c.degree;
          if (c.complete_subtree || count == c.count) {
            // Copy only child handles before callbacks; no borrowed entry lives
            // across arbitrary nested queries.
            std::array<detail::ChunkHandle, BitCapacity> children;
            if (c.complete_subtree) {
              size_t i = 0;
              c.subtree.children(node(), [&](TreeNode n) {
                children[i++] = {n.number, n.pos};
              });
            } else {
              std::copy_n(c.children.begin(), count, children.begin());
            }
            for (size_t i = 0; i < count; ++i) {
              visit(children[i].node());
            }
            return;
          }
        } else {
          const auto r = owner_->record(*tree_, node());
          if (r.count == r.degree) {
            for (size_t i = 0; i < r.count; ++i) {
              visit(r.children[i].node());
            }
            return;
          }
        }
        local_.for_each_child(std::forward<Visitor>(visit));
      }
      /** @brief Consume a complete decoded subtree, or return zero if it does
       * not fit. Returns its node count. Visitors follow unpack kernel events;
       * inspect receives each child once before descent. No callback is made
       * on failure. Callback-safe scratch copies are not persistent storage. */
      template <class Visitor, class Inspector>
      size_t inspect_decoded(Visitor& visit, Inspector& inspect) const {
        if (!owner_->store_.enabled()) {
          return 0;
        }
        if constexpr (detail::ChunkEncoding<Tree>::kind == 2) {
          return 0;
        } else {
          const auto& c = owner_->chunk(*tree_, node());
          // A containing chunk may start at an ancestor, not this node.
          if (!c.complete_subtree || c.anchor.number != node().number) {
            return 0;
          }
          const auto copy = c.subtree;
          if constexpr (!RetainChunks) {
            owner_->record(*tree_, node());
          }
          copy.inspect_depth_first(visit, inspect);
          return copy.size();
        }
      }
      /** @brief Requested persistent bytes, bounded by ByteBudget. */
      size_t cache_memory_usage_bytes() const {
        return owner_->memory_usage_bytes();
      }
    };
    /** @brief Bind a view; tree and context must outlive it. */
    NodeView navigate(const Tree& tree, TreeNode node) const {
      return {tree, node, *this};
    }
  };
};
}  // namespace pixie::experimental
