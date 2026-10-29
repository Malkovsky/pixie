#pragma once

#include <pixie/tree/cache.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

namespace pixie::experimental {

/**
 * @brief Bounded auxiliary navigation index for root-to-leaf walks.
 * @details The index borrows an immutable tree, which must outlive the index
 * and its cursors. It expands nodes in decreasing probability under uniform
 * child choices, using topology only (no sampled queries). Child lists are
 * stored contiguously. Unindexed descendants use the original tree.
 * Construction allocates at most Capacity entries within a runtime byte budget,
 * plus temporary build state. Zero disables indexing. Cursors require
 * this index to remain alive and unmoved for their lifetime.
 */
template <class Tree, uint32_t Capacity = 32768>
class TreeNavigationIndex {
  static_assert(Capacity > 0);
  static constexpr uint32_t kMissing = std::numeric_limits<uint32_t>::max();
  struct Entry {
    TreeNode node;
    size_t degree;
    uint32_t first = kMissing;
  };

 public:
  /** @brief Build from a nonempty tree that remains alive and unchanged. */
  explicit TreeNavigationIndex(const Tree& tree,
                               size_t byte_limit = kDefaultTreeCacheBytes)
      : tree_(tree) {
    if constexpr (requires { tree.cache_limit_bytes(); }) {
      byte_limit = std::min(byte_limit, tree.cache_limit_bytes());
    }
    const size_t capacity =
        byte_limit > sizeof(*this)
            ? std::min({size_t(Capacity), tree.size(),
                        (byte_limit - sizeof(*this)) / sizeof(Entry)})
            : 0;
    if (!capacity) {
      return;
    }
    entries_.reserve(capacity);
    // A conforming vector may overallocate reserve(); never exceed the cap.
    if (entries_.capacity() > capacity) {
      std::vector<Entry>().swap(entries_);
      return;
    }
    entries_.push_back({tree.root(), tree.degree(tree.root())});
    struct Pending {
      double probability;
      uint32_t index;
      bool operator<(const Pending& other) const {
        if (probability != other.probability) {
          return probability < other.probability;
        }
        return index > other.index;
      }
    };
    std::priority_queue<Pending> pending;
    pending.push({1.0, 0});
    while (!pending.empty()) {
      const auto current = pending.top();
      pending.pop();
      const auto node = entries_[current.index].node;
      const size_t degree = entries_[current.index].degree;
      if (!degree || degree > capacity - entries_.size()) {
        continue;
      }
      entries_[current.index].first = static_cast<uint32_t>(entries_.size());
      const double probability = current.probability / degree;
      tree.for_each_child(node, [&](TreeNode child) {
        const auto index = static_cast<uint32_t>(entries_.size());
        const size_t count = tree.is_leaf(child) ? 0 : tree.degree(child);
        entries_.push_back({child, count});
        if (count) {
          pending.push({probability, index});
        }
      });
    }
  }

  /** @brief Index object and reserved entry bytes, excluding the borrowed tree.
   */
  size_t memory_usage_bytes() const {
    return entries_.empty()
               ? 0
               : sizeof(*this) + entries_.capacity() * sizeof(Entry);
  }

  /** @brief Independent root-to-leaf cursor borrowing this index and its tree.
   */
  class Cursor {
   public:
    /** @brief Return the current valid node handle. */
    TreeNode node() const {
      return index_ == kMissing ? fallback_ : owner_->entries_[index_].node;
    }
    /** @brief Return the current node's number of immediate children. */
    size_t degree() const {
      if (index_ != kMissing) {
        return owner_->entries_[index_].degree;
      }
      return owner_->tree_.is_leaf(fallback_) ? 0
                                              : owner_->tree_.degree(fallback_);
    }
    /** @brief Descend to child i, requiring i in [0, degree()). */
    void move_to_child(size_t i) {
      if (index_ != kMissing) {
        const auto& entry = owner_->entries_[index_];
        if (entry.first != kMissing) {
          index_ = entry.first + static_cast<uint32_t>(i);
          return;
        }
      }
      fallback_ = uncached_child(i);
      index_ = kMissing;
    }
    /**
     * @brief Descend using a handle already obtained while inspecting children.
     * @param i Valid index in [0, degree()).
     * @param child Handle equal to child(node(), i); mismatches are undefined.
     */
    void move_to_child(size_t i, TreeNode child) {
      if (index_ != kMissing && owner_->entries_[index_].first != kMissing) {
        index_ = owner_->entries_[index_].first + static_cast<uint32_t>(i);
      } else {
        fallback_ = child;
        index_ = kMissing;
      }
    }
    /** @brief Visit immediate children with TreeBase::for_each_child semantics.
     */
    template <class Visitor>
    void for_each_child(Visitor&& visitor) const {
      if (index_ != kMissing) {
        const auto& entry = owner_->entries_[index_];
        if (entry.first != kMissing) {
          const auto* current = owner_->entries_.data() + entry.first;
          const auto* end = current + entry.degree;
          for (; current != end; ++current) {
            visitor(current->node);
          }
          return;
        }
        if (!entry.degree) {
          return;
        }
      }
      uncached_for_each_child(visitor);
    }

    /**
     * @brief Inspect all immediate children, then descend to child i.
     * @param i Child index in [0, degree()); invalid indices are undefined.
     * @param visitor Borrowed callback with for_each_child semantics.
     * @details Children are emitted in child-index order. Exceptions propagate
     * without moving the cursor. The callback must not modify this cursor or
     * its backing tree. Indexed lists need no per-child selection branch;
     * fallback enumeration retains the selected handle without another query.
     */
    template <class Visitor>
    void inspect_children_and_move_to_child(size_t i, Visitor&& visitor) {
      if (index_ != kMissing) {
        const auto& entry = owner_->entries_[index_];
        if (entry.first != kMissing) {
          const uint32_t selected = entry.first + static_cast<uint32_t>(i);
          const auto* current = owner_->entries_.data() + entry.first;
          const auto* end = current + entry.degree;
          for (; current != end; ++current) {
            visitor(current->node);
          }
          index_ = selected;
          return;
        }
      }
      auto selected = node();
      size_t offset = 0;
      auto inspect = [&](TreeNode child) {
        visitor(child);
        if (offset++ == i) {
          selected = child;
        }
      };
      uncached_for_each_child(inspect);
      fallback_ = selected;
      index_ = kMissing;
    }

   private:
    friend class TreeNavigationIndex;
    explicit Cursor(const TreeNavigationIndex& owner) : owner_(&owner) {
      if (owner.entries_.empty()) {
        index_ = kMissing;
        fallback_ = owner.tree_.root();
      }
    }
    // Keep uncommon succinct navigation out of the compact indexed hot path.
    [[gnu::noinline]] size_t uncached_degree() const {
      return owner_->tree_.is_leaf(fallback_) ? 0
                                              : owner_->tree_.degree(fallback_);
    }
    [[gnu::noinline]] TreeNode uncached_child(size_t i) const {
      return owner_->tree_.child(node(), i);
    }
    template <class Visitor>
    [[gnu::noinline]] void uncached_for_each_child(Visitor& visitor) const {
      owner_->tree_.for_each_child(node(), visitor);
    }
    const TreeNavigationIndex* owner_;
    uint32_t index_ = 0;
    TreeNode fallback_{0, 0};
  };

  /** @brief Return an independent cursor initially positioned at the root. */
  Cursor cursor() const { return Cursor(*this); }

 private:
  const Tree& tree_;
  std::vector<Entry> entries_;
};
}  // namespace pixie::experimental
