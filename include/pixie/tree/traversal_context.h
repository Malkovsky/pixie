#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/detail/indexed_tree_part.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::tree {
namespace detail::tree_traversal {
template <class Tree>
struct Encoding;
template <class B, detail::TreeEncodingStorage S>
struct Encoding<BPTree<B, S>> {
  static constexpr int value = 0;
};
template <class B, detail::TreeEncodingStorage S>
struct Encoding<DFUDSTree<B, S>> {
  static constexpr int value = 1;
};
template <class B, detail::TreeEncodingStorage S>
struct Encoding<LoudsTree<B, S>> {
  static constexpr int value = 2;
};
}  // namespace detail::tree_traversal

/**
 * @brief Single-operation traversal context borrowing a tree and optional
 * index.
 * @details The byte limit covers this context, reusable fetch scratch and one
 * output under construction; the immutable index is independently budgeted.
 * Fetched parts own their topology. Their source accessor borrows the
 * original tree, independently of this context.
 * Multiple results can coexist; operations on one context must not overlap.
 * Sequential thread handoff is allowed. No locks or shared cache lookup exist.
 */
template <class Tree>
class TreeTraversalContext {
  static constexpr int encoding = detail::tree_traversal::Encoding<Tree>::value;
  const Tree& tree_;
  std::span<const uint64_t> words_;
  const RootNavigationIndex<Tree>* index_ = nullptr;
  size_t byte_limit_;
  static constexpr size_t missing = std::numeric_limits<size_t>::max();
  struct Frame {
    TreeNode node;
    size_t parent, children = 0, degree = missing;
  };
  std::vector<Frame> scratch_;

 public:
  /** @brief Configure the byte limit for fetch().
   * @details Borrows tree and an optional independently budgeted index. They
   * must remain alive and unmoved while this context or its cursors use them.
   * Mutating operations on one context must not overlap.
   * @throws std::invalid_argument for an empty tree or an index for another
   * tree.
   */
  TreeTraversalContext(const Tree& tree,
                       size_t byte_limit = kDefaultTraversalBytes,
                       const RootNavigationIndex<Tree>* index = nullptr)
      : tree_(tree),
        words_(tree.encoding_words()),
        index_(index),
        byte_limit_(byte_limit) {
    if (tree.empty()) {
      throw std::invalid_argument("empty tree");
    }
    if (index && &index->tree_base() != &tree) {
      throw std::invalid_argument("root index belongs to another tree");
    }
  }

  /** @brief Fetch an owning TreePart rooted at a valid source node. Mutable
   * scratch is reused; operations on one context must not overlap. Handles
   * must be valid for tree(). max_bits is the sum of logical bits examined,
   * including repeated examinations across regions. Backend navigation used
   * to locate LOUDS child regions is separate from this scan charge.
   * Memory/bit exhaustion returns partial topology, never invalidates earlier
   * results, and may return an empty result. Output allocation is bounded by
   * byte_limit; caller-held output memory is accounted separately from context
   * scratch. No general cache lookup or deduplication is performed.
   */
  auto fetch(TreeNode root, size_t max_bits) {
    using Part = detail::IndexedTreePart<Tree, encoding == 2>;
    const size_t fixed = sizeof(*this) + sizeof(Part);
    const size_t available = byte_limit_ > fixed ? byte_limit_ - fixed : 0;
    const size_t capacity =
        std::min({max_bits, tree_.size(),
                  available / (Part::bytes_per_node + sizeof(Frame))});
    if (!capacity) {
      return Part(tree_, {});
    }
    scratch_.clear();
    scratch_.reserve(capacity);
    const size_t nbits =
        encoding == 0 ? 2 * tree_.size() : 2 * tree_.size() - 1;
    auto bit = [&](size_t pos) {
      --max_bits;
      return bool((words_[pos / 64] >> (pos % 64)) & 1);
    };
    auto add = [&](TreeNode node, size_t parent) {
      if (parent != missing) {
        ++scratch_[parent].children;
      }
      scratch_.push_back({node, parent});
      return scratch_.size() - 1;
    };
    if constexpr (encoding == 0) {
      // A contiguous BP prefix. Open nodes at the boundary stay incomplete.
      size_t parent = missing;
      for (size_t pos = root.pos; pos < nbits && max_bits; ++pos) {
        if (bit(pos)) {
          if (scratch_.size() == capacity) {
            break;
          }
          parent = add({root.number + scratch_.size(), pos}, parent);
        } else {
          auto& row = scratch_[parent];
          row.degree = row.children;
          parent = row.parent;
          if (parent == missing) {
            break;
          }
        }
      }
    } else {
      // A started unary record represents a node, even if its terminating
      // zero falls beyond the budget. Only the terminator proves its degree.
      auto decode = [&](size_t number, size_t& pos, size_t parent) {
        const size_t id = add({number, pos}, parent);
        const size_t start = pos;
        while (pos < nbits && max_bits) {
          if (!bit(pos++)) {
            scratch_[id].degree = pos - start - 1;
            break;
          }
        }
        return id;
      };
      size_t pos = root.pos;
      if constexpr (encoding == 1) {
        size_t parent = missing;
        while (max_bits && scratch_.size() < capacity) {
          parent = decode(root.number + scratch_.size(), pos, parent);
          while (parent != missing &&
                 scratch_[parent].children == scratch_[parent].degree) {
            parent = scratch_[parent].parent;
          }
          if (parent == missing) {
            break;
          }
        }
      } else {
        decode(root.number, pos, missing);
        for (size_t i = 0;
             i < scratch_.size() && max_bits && scratch_.size() < capacity;
             ++i) {
          const auto degree = scratch_[i].degree;
          if (!degree) {
            continue;
          }
          const auto child = tree_.first_child(scratch_[i].node);
          pos = child.pos;
          for (size_t k = 0;
               k < degree && max_bits && scratch_.size() < capacity; ++k) {
            decode(child.number + k, pos, i);
          }
        }
      }
    }
    std::vector<typename Part::Record> records;
    records.reserve(scratch_.size());
    // LOUDS emits sibling runs in consecutive local order starting at node 1.
    // For preorder scans the offsets index the explicit child buffer instead.
    size_t child_offset = encoding == 2 ? 1 : 0;
    for (auto& frame : scratch_) {
      const auto count = frame.children;
      records.emplace_back(frame.node, child_offset, count,
                           frame.degree == count);
      if constexpr (encoding != 2) {
        frame.children = child_offset;
      }
      child_offset += count;
    }
    typename Part::ChildStorage children;
    if constexpr (encoding != 2) {
      children.resize(child_offset);
      for (size_t i = 1; i < scratch_.size(); ++i) {
        children[scratch_[scratch_[i].parent].children++] = i;
      }
    }
    scratch_.clear();
    return Part(tree_, std::move(records), std::move(children));
  }
  /** @brief Bytes of allocated reusable scan scratch. */
  size_t scan_scratch_bytes() const {
    return scratch_.capacity() * sizeof(Frame);
  }

  /** @brief Active fetch budget charge. Zero before scan scratch is allocated.
   */
  size_t budget_usage_bytes() const {
    return scratch_.capacity() ? memory_usage_bytes() : 0;
  }
  /** @brief Context and allocated scan scratch, excluding index and outputs. */
  size_t memory_usage_bytes() const {
    return sizeof(*this) + scan_scratch_bytes();
  }
  /** @brief Reject a temporary tree because this context borrows it. */
  TreeTraversalContext(Tree&&,
                       size_t = kDefaultTraversalBytes,
                       const RootNavigationIndex<Tree>* = nullptr) = delete;
  /** @brief Reject a temporary const tree for the same reason. */
  TreeTraversalContext(const Tree&&,
                       size_t = kDefaultTraversalBytes,
                       const RootNavigationIndex<Tree>* = nullptr) = delete;
  /** @brief Independent root cursor borrowing this index and its tree. */
  auto cursor() const& {
    return index_ ? index_->cursor()
                  : typename RootNavigationIndex<Tree>::Cursor(tree_);
  }
  /** @brief Reject cursors borrowing a temporary context. */
  auto cursor() const&& = delete;

  /** @brief Borrow the immutable tree; this context does not extend its
   * lifetime. */
  const Tree& tree() const { return tree_; }
  /** @brief Complete encoding words borrowed by the context, including final
   * padding. */
  std::span<const uint64_t> source_words() const { return words_; }
  /** @brief Borrow the immutable root index; cursors must not outlive it. */
  const RootNavigationIndex<Tree>* root_index() const { return index_; }
};
}  // namespace pixie::tree
