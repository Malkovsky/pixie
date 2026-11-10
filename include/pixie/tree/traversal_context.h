#pragma once

#include <pixie/tree/bp.h>
#include <pixie/tree/decoded_tree.h>
#include <pixie/tree/detail/bp_micro.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/louds_inspection_iterative.h>
#include <pixie/tree/detail/traversal_profiles.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
#include <pixie/tree/root_navigation_index.h>

#include <memory>

namespace pixie {
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
struct SizingVisitor {
  void operator()(TreeNode) const {}
};
}  // namespace detail::tree_traversal

/**
 * @brief Single-operation traversal context borrowing a tree and optional
 * index.
 * @details The byte limit covers this context, reusable scan scratch and its
 * latest retained result; the immutable index is independently budgeted.
 * Legacy whole-walk decoder scratch is a reserved allowance, allocated by the
 * walker. scan/fetch use context-owned scratch. Returned results own their
 * records but borrow the tree; callers must preserve its lifetime and address.
 * Multiple results can coexist; operations on one context must not overlap.
 * Sequential thread handoff is allowed. No locks or shared cache lookup exist.
 */
template <class Tree>
class TreeTraversalContext {
  static constexpr int encoding = detail::tree_traversal::Encoding<Tree>::value;
  static constexpr size_t kStreamBytes = 4096;
  struct Plan {
    size_t scratch = 0;
    size_t bits = 0;
    size_t depth = 0;
  };
  template <size_t Bits, class V, class I>
  using ParenthesisWalk = std::conditional_t<
      encoding == 0,
      detail::bp_micro::BpMicroTraversal<Bits, Tree, V, I>,
      detail::dfuds_micro::DfudsMicroTraversal<Bits, Tree, V, I>>;
  using SizingVisitor = detail::tree_traversal::SizingVisitor;
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Plan plan_;
  const RootNavigationIndex<Tree>* index_ = nullptr;
  size_t byte_limit_;
  struct Frame {
    size_t node, remaining;
  };
  std::vector<Frame> scratch_;
  std::shared_ptr<const DecodedTree<Tree>> retained_;

  static constexpr size_t overhead() { return sizeof(TreeTraversalContext); }
  template <size_t Bits>
  static Plan parentheses_plan(size_t allowance) {
    using Walk = ParenthesisWalk<Bits, SizingVisitor, SizingVisitor>;
    if (Walk::scratch_bytes(8) > allowance) {
      return {};
    }
    const size_t depth =
        std::min(size_t{128},
                 (allowance - Walk::scratch_bytes(0)) / (64 * sizeof(size_t)));
    return {std::max(kStreamBytes, Walk::scratch_bytes(depth)), Bits, depth};
  }
  template <size_t Depth, size_t Children>
  static Plan louds_plan(size_t allowance) {
    using Walk = detail::tree_unpack::LoudsIterativeTraversal<
        Tree, SizingVisitor, SizingVisitor, true, Depth, Children>;
    constexpr size_t bytes = std::max(kStreamBytes, Walk::scratch_bytes());
    return bytes <= allowance ? Plan{bytes, 0, Depth} : Plan{};
  }
  static Plan choose_plan(size_t budget) {
    if (budget < sizeof(TreeTraversalContext) + kStreamBytes) {
      return {};
    }
    const size_t available = budget - sizeof(TreeTraversalContext);
    const size_t allowance = std::min(size_t{encoding == 2 ? 80 : 192} * 1024,
                                      std::max(kStreamBytes, available / 2));
    using Profiles =
        std::conditional_t<encoding == 2, detail::tree_traversal::LoudsProfiles,
                           detail::tree_traversal::ParenthesisProfiles>;
    Plan plan;
    if (Profiles::visit_until([&](auto profile) {
          using Profile = decltype(profile);
          if constexpr (encoding == 2) {
            plan = louds_plan<Profile::depth, Profile::children>(allowance);
          } else {
            plan = parentheses_plan<Profile::bits>(allowance);
          }
          return plan.depth != 0;
        })) {
      return plan;
    }
    return {kStreamBytes, 0, 0};
  }

 public:
  /** @brief Configure scan storage and whole-walk scratch allowance.
   * @details Borrows tree and an optional independently budgeted index. They
   * must remain alive and unmoved while this context or its cursors use them.
   * Mutating operations on one context must not overlap. Returned decoded
   * results borrow only the tree and may outlive the context.
   * @throws std::invalid_argument for an empty tree or an index for another
   * tree.
   */
  TreeTraversalContext(const Tree& tree,
                       size_t byte_limit = kDefaultTraversalBytes,
                       const RootNavigationIndex<Tree>* index = nullptr)
      : tree_(tree),
        words_(tree.encoding_words()),
        plan_(choose_plan(byte_limit)),
        index_(index),
        byte_limit_(byte_limit) {
    if (tree.empty()) {
      throw std::invalid_argument("empty tree");
    }
    if (index && &index->tree_base() != &tree) {
      throw std::invalid_argument("root index belongs to another tree");
    }
  }

  /** @brief Decode bounded topology without retaining the result. Mutable
   * scratch is reused; operations on one context must not overlap. Handles
   * must be valid for tree(). max_bits is the sum of logical bits examined,
   * including repeated examinations across regions. Backend navigation used
   * to locate LOUDS child regions is separate from this scan charge.
   * Memory/bit exhaustion returns partial topology, never invalidates earlier
   * results, and may return an empty result. Output allocation is bounded by
   * byte_limit; caller-held output memory is accounted separately from context
   * scratch. No general cache lookup or deduplication is performed.
   */
  std::shared_ptr<const DecodedTree<Tree>> scan(TreeNode root,
                                                size_t max_bits) {
    using Result = DecodedTree<Tree>;
    auto out = std::make_shared<Result>(tree_);
    const size_t fixed = sizeof(*this) + sizeof(Result) + plan_.scratch;
    const size_t available = byte_limit_ > fixed ? byte_limit_ - fixed : 0;
    const size_t capacity = std::min(
        {max_bits, tree_.size(),
         available / (sizeof(typename Result::Record) + sizeof(Frame))});
    if (!capacity) {
      return out;
    }
    out->records_.reserve(capacity);
    scratch_.clear();
    scratch_.reserve(capacity);
    const size_t nbits =
        encoding == 0 ? 2 * tree_.size() : 2 * tree_.size() - 1;
    auto bit = [&](size_t pos) {
      ++out->bits_;
      return bool((words_[pos / 64] >> (pos % 64)) & 1);
    };
    if constexpr (encoding == 0) {
      // A contiguous BP prefix. Open nodes at the boundary stay incomplete.
      size_t number = root.number;
      ++out->regions_;
      for (size_t pos = root.pos; pos < nbits && out->bits_ < max_bits; ++pos) {
        if (bit(pos)) {
          if (out->records_.size() == capacity) {
            break;
          }
          auto parent =
              scratch_.empty() ? Result::missing : scratch_.back().node;
          auto id = out->add({number++, pos}, parent);
          scratch_.push_back({id, 0});
        } else {
          out->records_[scratch_.back().node].complete = true;
          scratch_.pop_back();
          if (scratch_.empty()) {
            break;
          }
        }
      }
    } else if constexpr (encoding == 1) {
      // DFUDS is a contiguous preorder stream of unary degrees.
      size_t number = root.number, pos = root.pos;
      ++out->regions_;
      while (pos < nbits && out->bits_ < max_bits && out->size() < capacity) {
        auto parent = scratch_.empty() ? Result::missing : scratch_.back().node;
        auto id = out->add({number++, pos}, parent);
        if (!scratch_.empty()) {
          if (--scratch_.back().remaining == 0) {
            out->records_[parent].complete = true;
          }
        }
        size_t degree = 0;
        bool known = false;
        while (pos < nbits && out->bits_ < max_bits) {
          if (!bit(pos++)) {
            known = true;
            break;
          }
          ++degree;
        }
        if (!known) {
          break;
        }
        if (degree) {
          scratch_.push_back({id, degree});
        } else {
          out->records_[id].complete = true;
        }
        while (!scratch_.empty() && scratch_.back().remaining == 0) {
          scratch_.pop_back();
        }
        if (scratch_.empty()) {
          break;
        }
      }
    } else {
      // Each child run is a contiguous region; locate runs with native select.
      // Decode delimiters directly, rather than selecting each child handle.
      auto decode = [&](TreeNode node, size_t parent) {
        if (out->bits_ == max_bits || out->size() == capacity) {
          return Result::missing;
        }
        ++out->regions_;
        const size_t id = out->add(node, parent);
        size_t pos = node.pos, degree = 0;
        bool known = false;
        while (pos < nbits && out->bits_ < max_bits) {
          if (!bit(pos++)) {
            known = true;
            break;
          }
          ++degree;
        }
        scratch_.push_back({pos, known ? degree : Result::missing});
        if (known && degree == 0) {
          out->records_[id].complete = true;
        }
        return id;
      };
      decode(root, Result::missing);
      for (size_t parent = 0; parent < out->size(); ++parent) {
        size_t degree = scratch_[parent].remaining;
        if (!degree || degree == Result::missing) {
          continue;
        }
        auto child = tree_.first_child(out->records_[parent].node);
        size_t count = 0;
        for (; count < degree; ++count) {
          size_t id = decode(child, parent);
          if (id == Result::missing) {
            break;
          }
          child = TreeNode(child.number + 1, scratch_[id].node);
          if (scratch_[id].remaining == Result::missing) {
            ++count;
            break;
          }
        }
        out->records_[parent].complete = count == degree;
        if (out->bits_ == max_bits || out->size() == capacity) {
          break;
        }
      }
    }
    return out;
  }
  /** @brief Scan and retain the latest result. Older returned results survive.
   * No lookup is performed. See scan() for bounds and lifetime requirements. */
  std::shared_ptr<const DecodedTree<Tree>> fetch(TreeNode node,
                                                 size_t max_bits) {
    auto result = scan(node, max_bits);
    retained_ = result->empty() ? nullptr : result;
    return result;
  }
  /** @brief Release the context's result reference, retaining reusable scratch.
   */
  void clear() { retained_.reset(); }
  /** @brief Bytes of allocated reusable scan scratch. */
  size_t scan_scratch_bytes() const {
    return scratch_.capacity() * sizeof(Frame);
  }
  /** @brief Bytes of the latest fetched result; may also be held by callers. */
  size_t retained_bytes() const {
    return retained_ ? retained_->memory_usage_bytes() : 0;
  }

  /** @brief Peak explicit traversal-buffer allowance, including decoder loads.
   */
  size_t scratch_bytes() const { return plan_.scratch; }
  /** @brief Whether plain DFS uses its streaming decoder. */
  bool streaming_enabled() const { return plan_.scratch >= kStreamBytes; }
  /** @brief Whether inspected DFS can use decoded subtrees/level cursors. */
  bool optimized_inspection_enabled() const { return plan_.depth != 0; }
  /** @brief BP/DFUDS chunk size in bits; zero for LOUDS or direct inspection.
   */
  size_t decoded_chunk_bits() const { return plan_.bits; }
  /** @brief Cached traversal depth; BP/DFUDS keep decoding below this depth.
   */
  size_t cached_depth() const { return plan_.depth; }
  /** @brief Active index/context/scratch budget charge; zero when disabled. */
  size_t budget_usage_bytes() const {
    const size_t charged =
        plan_.scratch + scan_scratch_bytes() + retained_bytes();
    return charged ? sizeof(*this) + charged : size_t{0};
  }
  /** @brief Context bytes and reserved/allocated scratch plus retained result;
   * excludes the borrowed index and results held only by callers. */
  size_t memory_usage_bytes() const {
    return std::max(sizeof(*this), budget_usage_bytes());
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
}  // namespace pixie
