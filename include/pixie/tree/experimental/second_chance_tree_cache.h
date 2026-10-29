#pragma once

// Experiment snapshot, 2026-10-24. Default budget: 512 KiB. Operation-based
// admission keeps sequential queries in per-view memoization; degree/indexed
// child requests activate shared CLOCK caching. Tree defaults are unchanged.
// Release -O3 -DNDEBUG -march=native, CPU 0, five interleaved repetitions,
// 0.2s minimum, 0.1s warmup. Even powers 2^10..2^26 plus 10M random nodes.
// Worst median CLOCK/local CPU-time ratio across the ten sizes (<=1.2 target):
//
// | Tree  | DFS    | DFS + children | Allocated cache |
// | ----- | ------ | -------------- | --------------- |
// | BP    | 0.7166 | 0.8237         | 288.02 KiB      |
// | DFUDS | 1.0076 | 0.8988         | 304.02 KiB      |
// | LOUDS | 1.0081 | 1.0087         | 272.02 KiB      |
//
// 60 comparisons / 600 repetitions passed. Allocated memory measured in mixed
// path/DFS runs; pure sequential DFS does not allocate the shared table.
// This bound covers the existing sequential DFS workloads, not DFS written
// using degree/indexed-child calls at every node. Full evidence and commands:
// profiling/clock512-2026-10-24/README.md

#include <pixie/tree/cache.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

namespace pixie::experimental {

/**
 * @brief Experimental fixed-budget CLOCK (second-chance) navigation policy.
 * @tparam ByteBudget Ceiling for context, allocated payload and temporary query
 * buffers, excluding allocator bookkeeping and caller-owned node views.
 * The tree constructor may lower this ceiling. Storage is allocated on the
 * first shared query and never grows; insufficient budgets use local views.
 * @tparam ChildPrefix Maximum child handles retained per node. Larger lists
 * remain partially cached and use representation navigation for missing data.
 * @tparam TrackStatistics Enable hit/miss/eviction counters for diagnostics.
 * @tparam SnapshotChildren Copy complete child lists before invoking callbacks,
 * avoiding repeated lookups while allowing nested queries to evict entries.
 * @details Implements tree/cache.h. Hits set a recent bit. A persistent clock
 * hand clears recent bits and evicts the first unmarked entry. Views contain
 * stable node handles, never pointers into entries. Enumeration copies handles
 * before callbacks; nested queries may evict entries safely. No entries need
 * pinning because this API exposes no borrowed cache storage. The tree and its
 * source must outlive views and remain unmoved. Shared access, even through
 * different views, requires external synchronization. Copying or moving a
 * context starts with an empty cache, avoiding references to the old tree.
 * Direct tree queries bypass the policy. No change to default tree policies.
 * Sequential operations (leaf/parent/sibling/first-child) use per-view local
 * memoization without shared lookup or admission. Degree and indexed-child
 * queries activate shared caching for that view until reset. Enumeration uses
 * shared caching only after that activation; otherwise it decodes locally.
 * This operation-based admission avoids polluting shared storage with streams
 * of one-use nodes. It does not inspect workload names or detect traversals.
 */
template <std::size_t ByteBudget = 512 * 1024,
          std::size_t ChildPrefix = 8,
          bool TrackStatistics = false,
          bool SnapshotChildren = false>
struct SecondChanceTreeCache {
  /** @brief Query lookup counters; zero when instrumentation is disabled. */
  struct Statistics {
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t evictions = 0;
  };

  /** @brief Per-tree context with a fixed directory and fixed-size entries. */
  template <class Tree>
  class Cache {
    struct Handle {
      std::size_t number = 0;
      std::size_t pos = 0;
    };
    struct Entry {
      std::optional<typename Tree::LocalNodeView> local;
      std::optional<std::size_t> degree;
      std::array<Handle, ChildPrefix> children{};
      std::size_t children_size = 0;
      bool recent = false;
    };
    static constexpr size_t kScratchBytes =
        ChildPrefix * sizeof(Handle) + sizeof(typename Tree::LocalNodeView);
    // Reserve space for the context itself (pointer, clock hand, counters).
    static_assert(ByteBudget >= 128 + sizeof(Entry) + 2 * sizeof(uint32_t));

   public:
    /** @brief Number of resident nodes; rounded down to a power of two. */
    static constexpr std::size_t capacity = std::bit_floor(
        (ByteBudget - 128) / (sizeof(Entry) + 2 * sizeof(uint32_t)));

   private:
    static_assert(capacity < UINT32_MAX / 2);
    struct Storage {
      explicit Storage(size_t count)
          : entries(std::make_unique<Entry[]>(count)),
            directory(std::make_unique<uint32_t[]>(2 * count)) {}
      std::unique_ptr<Entry[]> entries;
      // Zero is empty; otherwise the entry index plus one.
      std::unique_ptr<uint32_t[]> directory;
    };
    struct Empty {};
    mutable std::unique_ptr<Storage> storage_;
    mutable std::size_t hand_ = 0;
    mutable std::size_t used_ = 0;
    std::size_t capacity_ = 0;
    std::size_t kMask = 0;
    [[no_unique_address]] mutable std::
        conditional_t<TrackStatistics, Statistics, Empty> statistics_;

    std::size_t hash(std::size_t key) const {
      uint64_t x = key;
      x ^= x >> 30;
      x *= 0xbf58476d1ce4e5b9ULL;
      x ^= x >> 27;
      x *= 0x94d049bb133111ebULL;
      return (x ^ (x >> 31)) & kMask;
    }

    void erase(std::size_t key) const {
      auto& directory = storage_->directory;
      auto hole = hash(key);
      while (storage_->entries[directory[hole] - 1].local->node().number !=
             key) {
        hole = (hole + 1) & kMask;
      }
      // Backward-shift deletion keeps probe chains valid without tombstones.
      for (auto pos = (hole + 1) & kMask; directory[pos];
           pos = (pos + 1) & kMask) {
        const auto home =
            hash(storage_->entries[directory[pos] - 1].local->node().number);
        if (((hole - home) & kMask) < ((pos - home) & kMask)) {
          directory[hole] = directory[pos];
          hole = pos;
        }
      }
      directory[hole] = 0;
    }

    Entry& get(const Tree& tree, TreeNode node) const {
      if (!storage_) {
        storage_ = std::make_unique<Storage>(capacity_);
      }
      auto& directory = storage_->directory;
      auto bucket = hash(node.number);
      while (directory[bucket]) {
        auto& entry = storage_->entries[directory[bucket] - 1];
        if (entry.local->node().number == node.number) {
          entry.recent = true;
          if constexpr (TrackStatistics) {
            ++statistics_.hits;
          }
          return entry;
        }
        bucket = (bucket + 1) & kMask;
      }
      if constexpr (TrackStatistics) {
        ++statistics_.misses;
      }
      std::size_t slot;
      if (used_ < capacity_) {
        slot = used_++;
      } else {
        while (storage_->entries[hand_].recent) {
          storage_->entries[hand_].recent = false;
          hand_ = (hand_ + 1) % capacity_;
        }
        slot = hand_;
        hand_ = (hand_ + 1) % capacity_;
        erase(storage_->entries[slot].local->node().number);
        if constexpr (TrackStatistics) {
          ++statistics_.evictions;
        }
      }
      auto& entry = storage_->entries[slot];
      entry.local.emplace(tree, node);
      entry.degree.reset();
      entry.children_size = 0;
      entry.recent = true;
      bucket = hash(node.number);
      while (directory[bucket]) {
        bucket = (bucket + 1) & kMask;
      }
      directory[bucket] = static_cast<uint32_t>(slot + 1);
      return entry;
    }

   public:
    /** @brief Construct an empty context; first access allocates its storage.
     */
    explicit Cache(std::size_t budget = ByteBudget) {
      budget = std::min(budget, ByteBudget);
      constexpr auto overhead = sizeof(Cache) + sizeof(Storage) + kScratchBytes;
      capacity_ =
          budget > overhead
              ? std::min(capacity,
                         std::bit_floor((budget - overhead) /
                                        (sizeof(Entry) + 2 * sizeof(uint32_t))))
              : 0;
      kMask = capacity_ ? 2 * capacity_ - 1 : 0;
    }
    /** @brief Copy the policy, starting with an empty context. */
    Cache(const Cache& other)
        : capacity_(other.capacity_), kMask(other.kMask) {}
    /** @brief Move the policy, starting with an empty context. */
    Cache(Cache&& other) noexcept : Cache(static_cast<const Cache&>(other)) {}
    /** @brief Discard cached data; existing borrowed ranges are not exposed. */
    Cache& operator=(const Cache& other) {
      storage_.reset();
      hand_ = used_ = 0;
      statistics_ = {};
      capacity_ = other.capacity_;
      kMask = other.kMask;
      return *this;
    }
    /** @brief Move assignment discards cached data. */
    Cache& operator=(Cache&& other) noexcept {
      return operator=(static_cast<const Cache&>(other));
    }
    /** @brief Context, requested storage and reserved query-buffer bytes;
     * excludes allocator overhead. Zero when shared retention is disabled. */
    std::size_t memory_usage_bytes() const {
      return capacity_ ? sizeof(Cache) + kScratchBytes +
                             (storage_ ? sizeof(Storage) +
                                             capacity_ * (sizeof(Entry) +
                                                          2 * sizeof(uint32_t))
                                       : 0)
                       : 0;
    }
    /** @brief Return counters, or zeros when TrackStatistics is false. */
    Statistics statistics() const {
      if constexpr (TrackStatistics) {
        return statistics_;
      } else {
        return {};
      }
    }

    /** @brief Stable-handle node view resolving evicted data on the next query.
     */
    class NodeView {
     public:
      /** @brief Return the current node handle by value. */
      TreeNode node() const { return node_; }
      /** @brief Move to a valid node of the same tree. */
      void reset(TreeNode node) {
        node_ = node;
        local_.reset(node);
        shared_ = false;
      }
      /** @brief Test whether the current node is root. */
      bool is_root() const { return tree_->is_root(node_); }
      /** @brief Test whether the current node has no children. */
      bool is_leaf() const { return local_.is_leaf(); }
      /** @brief Return the current node's child count. */
      std::size_t degree() const {
        if (!owner_->capacity_) {
          return local_.degree();
        }
        auto& entry = get();
        if (!entry.degree) {
          entry.degree = entry.local->degree();
        }
        return *entry.degree;
      }
      /** @brief Return child i, requiring i in [0, degree()). */
      TreeNode child(std::size_t i) const {
        if (!owner_->capacity_) {
          return local_.child(i);
        }
        auto& entry = get();
        if (i >= ChildPrefix) {
          return entry.local->child(i);
        }
        while (entry.children_size <= i) {
          const auto n = entry.children_size;
          const auto child =
              n == 0
                  ? entry.local->first_child()
                  : tree_->next_sibling(TreeNode(entry.children[n - 1].number,
                                                 entry.children[n - 1].pos));
          entry.children[n] = {child.number, child.pos};
          ++entry.children_size;
        }
        const auto result = entry.children[i];
        return TreeNode(result.number, result.pos);
      }
      /** @brief Return the first child; requires a non-leaf node. */
      TreeNode first_child() const { return local_.first_child(); }
      /** @brief Return the parent, or root itself for root. */
      TreeNode parent() const { return local_.parent(); }
      /** @brief Return next sibling; requires a non-root, non-last child. */
      TreeNode next_sibling() const { return local_.next_sibling(); }
      /** @brief Test for the last sibling; requires a non-root node. */
      bool is_last_child() const { return local_.is_last_child(); }
      /**
       * @brief Enumerate children with TreeBase::for_each_child semantics.
       * @details Nested queries and exceptions are supported. Each callback
       * receives a value; no cache reference survives a callback. The callback
       * must not reset this view or mutate its tree.
       */
      template <class Visitor>
      void for_each_child(Visitor&& visitor) const {
        if (!shared_) {
          local_.for_each_child(std::forward<Visitor>(visitor));
          return;
        }
        const auto& entry = get();
        if (entry.degree && *entry.degree == entry.children_size) {
          const auto count = *entry.degree;
          if constexpr (SnapshotChildren) {
            // Copy before callbacks: nested queries may evict this entry.
            std::array<Handle, ChildPrefix> children;
            std::copy_n(entry.children.begin(), count, children.begin());
            for (std::size_t i = 0; i < count; ++i) {
              visitor(TreeNode(children[i].number, children[i].pos));
            }
          } else {
            for (std::size_t i = 0; i < count; ++i) {
              visitor(child(i));
            }
          }
          return;
        }
        std::size_t count = 0;
        // Preserve per-node query results without borrowing an entry across
        // callbacks: a nested query may evict it while enumeration is active.
        auto local = *entry.local;
        local.for_each_child([&](TreeNode child) {
          if (count < ChildPrefix) {
            auto& current = get();
            // A nested callback may have evicted earlier prefix elements.
            if (current.children_size == count) {
              current.children[count] = {child.number, child.pos};
              ++current.children_size;
            }
          }
          ++count;
          visitor(child);
        });
        get().degree = count;
      }
      /** @brief Return the owning context's diagnostic counters. */
      Statistics cache_statistics() const { return owner_->statistics(); }
      /** @brief Return context plus allocated cache payload bytes. */
      std::size_t cache_memory_usage_bytes() const {
        return owner_->memory_usage_bytes();
      }

     private:
      friend class Cache;
      NodeView(const Tree& tree, TreeNode node, const Cache& owner)
          : tree_(&tree), node_(node), owner_(&owner), local_(tree, node) {}
      Entry& get() const {
        shared_ = true;
        return owner_->get(*tree_, node_);
      }
      const Tree* tree_;
      TreeNode node_;
      const Cache* owner_;
      typename Tree::LocalNodeView local_;
      mutable bool shared_ = false;
    };
    /** @brief Bind a view to a valid node; context and tree must outlive it. */
    NodeView navigate(const Tree& tree, TreeNode node) const {
      return NodeView(tree, node, *this);
    }
  };
};

}  // namespace pixie::experimental
