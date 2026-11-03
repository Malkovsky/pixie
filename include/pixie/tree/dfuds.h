#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/construction.h>
#include <pixie/tree/detail/dfuds_micro.h>
#include <pixie/tree/detail/dfuds_node_view.h>
#include <pixie/tree/detail/encoding_storage.h>
#include <pixie/tree/detail/serialization.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief A tree class based on the depth-first unary degree sequence (DFUDS)
 * representation
 */
template <class RMMTree,
          detail::TreeEncodingStorage EncodingStorage = AlignedStorage>
class DFUDSTree;

template <class B, detail::TreeEncodingStorage S>
struct detail::ConstructionTraversal<DFUDSTree<B, S>> {
  using Tree = DFUDSTree<B, S>;
  template <class Enter, class Leave>
  static void walk(const Tree& tree, Enter& enter, Leave& leave) {
    inspected_construction_dfs(
        tree, enter, leave, [&](auto& visit, auto& child) {
          dfuds_micro::DfudsMicroTraversal<
              8192, Tree, std::remove_reference_t<decltype(visit)>,
              std::remove_reference_t<decltype(child)>>
              traversal(tree, tree.encoding_words(), visit, child, 128);
          traversal.walk(tree.root());
        });
  }
  static constexpr size_t scratch_bytes() {
    using Callback = void (*)(TreeNode);
    return dfuds_micro::DfudsMicroTraversal<8192, Tree, Callback,
                                            Callback>::scratch_bytes(128);
  }
};

template <class B, detail::TreeEncodingStorage S>
struct TreeTraits<DFUDSTree<B, S>> {
  using TreeView = DFUDSTree<B, ReadOnlyStorageView>;
  using NodeView = detail::DfudsNodeView<DFUDSTree<B, S>>;
};

template <class RMMTree, detail::TreeEncodingStorage EncodingStorage>
class DFUDSTree : public TreeBase<DFUDSTree<RMMTree, EncodingStorage>> {
 private:
  DFUDSTree(EncodingStorage storage, size_t bits, RMMTree backend)
      : encoding_(std::move(storage)),
        num_bits_(bits),
        rmm_(std::move(backend)) {}

  EncodingStorage encoding_;
  size_t num_bits_;
  RMMTree rmm_;

  // Sibling edges match the parent's unary ones from right to left. The last
  // child's edge is therefore the first one of that run, preceded by zero
  // (or at bit zero for the root's degree run).
  bool starts_degree_run(size_t open) const {
    return open == 0 || rmm_.bit(open - 1) == 0;
  }

  // For a non-root first child, the parent's degree run ends in the matching
  // adjacent "10" pair immediately before the child's record.
  size_t incoming_open(const TreeNode& node) const {
    return rmm_.bit(node.pos - 2) ? node.pos - 2 : rmm_.open(node.pos - 1);
  }

  TreeNode parent_from_open(const TreeNode& node, size_t open) const {
    // [open, node.pos) is balanced: half its bits are zero delimiters.
    // Subtract those records from the child's preorder number.
    const auto number = node.number - (node.pos - open) / 2;
    // For a last child, the matching one is already the parent's record start.
    const auto pos = number == 0               ? 0
                     : starts_degree_run(open) ? open
                                               : rmm_.select0(number) + 1;
    return TreeNode(number, pos);
  }

 public:
  /**
   * @brief Construct an owning copy of another ordered tree's topology.
   * @details Preserves child order, but assigns destination node numbers and
   * bit positions. Source must satisfy detail::OrderedTreeSource and describe
   * a valid tree of size() nodes. Its handles must remain valid during this
   * call; no source storage is retained. Navigation metadata is rebuilt.
   * Succinct sources use decoded DFS; other sources use O(height) navigation
   * scratch. LOUDS output additionally stages O(size()) packed bits and
   * O(height) level records; vector capacity adds allocation overhead.
   * Construction scratch is independent of the query traversal-cache budget.
   * Allocation and unrepresentable-size exceptions propagate.
   */
  template <detail::OrderedTreeSource Source>
    requires(std::same_as<EncodingStorage, AlignedStorage> &&
             !std::same_as<std::remove_cvref_t<Source>, DFUDSTree>)
  explicit DFUDSTree(const Source& source)
      : DFUDSTree(
            detail::construct_tree_encoding<detail::TreeEncodingKind::kDFUDS>(
                source),
            source.size()) {}

  /** @brief Write one tree artifact: encoding, node count and nested backend.
   * @details External root indexes and active traversal state are not included.
   */
  void serialize_impl(BinaryWriter& writer) const
    requires Serializable<RMMTree>
  {
    detail::serialize_tree<detail::TreeEncodingKind::kDFUDS>(
        writer, this->size(), encoding_words_impl(), rmm_);
  }

  /** @brief Restore a complete tree without rebuilding its navigation index.
   * @param reader Cursor advanced only after successful restoration.
   * @param validation Quick checks framing, sizes and backend metadata; it
   * assumes the saved topology and metadata correspond. Full also scans and
   * validates topology and source-derived backend data.
   * @details Owning trees copy the encoding. TreeView borrows word-aligned
   * backing bytes, which must remain alive and immutable. Backend ownership
   * follows its deserializer. No hashes or traversal caches are restored.
   * @throws std::invalid_argument for malformed or incompatible artifacts;
   * allocation and size exceptions propagate. Failure leaves reader unchanged.
   */
  static DFUDSTree deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires Deserializable<RMMTree,
                            std::span<const uint64_t>,
                            DeserializationValidation>
  {
    return detail::deserialize_tree<detail::TreeEncodingKind::kDFUDS,
                                    EncodingStorage, RMMTree>(
        reader, validation,
        [](EncodingStorage&& encoding, size_t bits, RMMTree&& backend) {
          return DFUDSTree(std::move(encoding), bits, std::move(backend));
        });
  }

  using Node = TreeNode;

  using typename TreeBase<DFUDSTree>::NodeView;

 private:
  friend class detail::DfudsNodeView<DFUDSTree>;

 public:
  /** @brief Implement TreeBase::navigate without mutating the tree. */
  NodeView navigate_impl(TreeNode node) const { return NodeView(*this, node); }

  /** @brief Object, owned encoding and backend allocations; excludes borrowed
   * words. */
  size_t memory_usage_bytes_impl() const
    requires HasMemoryUsageBytes<RMMTree>
  {
    return sizeof(*this) + encoding_owned_bytes_impl() +
           nested_owned_memory_bytes(rmm_);
  }

  /**
   * @brief Construct from a valid encoding of tree_size nodes.
   * @details Owning trees copy complete source words; TreeView borrows them.
   * A view's source must stay alive, immutable and at a stable address.
   * Invalid topology is undefined behavior. A short source throws
   * std::invalid_argument; an unrepresentable bit count throws
   * std::length_error.
   */
  explicit DFUDSTree(std::span<const uint64_t> words, size_t tree_size)
      : DFUDSTree(detail::tree_encoding_from_words<EncodingStorage>(
                      words,
                      detail::tree_encoding_bits(tree_size, false)),
                  tree_size) {}

  /**
   * @brief Take storage by value: moves transfer ownership, lvalues copy.
   * @details ReadOnlyStorageView always borrows. AlignedStorage may have a
   * partial final logical byte/word; its allocation supplies word padding.
   * View storage must contain aligned complete words. Short/misaligned storage
   * throws std::invalid_argument. Topology and size must agree.
   */
  explicit DFUDSTree(EncodingStorage storage, size_t tree_size)
      : encoding_(std::move(storage)),
        num_bits_(detail::tree_encoding_bits(tree_size, false)),
        rmm_(encoding_words_impl(), num_bits_) {}

  /** @brief Deep-copy owning encodings and rebuild their backend references.
   * @details TreeView copies backend metadata and borrows the same encoding.
   */
  DFUDSTree(const DFUDSTree& other)
      : encoding_(other.encoding_), num_bits_(other.num_bits_), rmm_([&] {
          if constexpr (std::same_as<EncodingStorage, AlignedStorage>) {
            return RMMTree(encoding_words_impl(), num_bits_);
          } else {
            return other.rmm_;
          }
        }()) {}
  /** @brief Copy independently; invalidates existing borrowers of *this. */
  DFUDSTree& operator=(const DFUDSTree& other) {
    if (this != &other) {
      DFUDSTree copy(other);
      *this = std::move(copy);
    }
    return *this;
  }
  /** @brief Transfer storage/indexes; invalidates borrowers of the source.
   * @details The moved-from object may only be destroyed or assigned.
   */
  DFUDSTree(DFUDSTree&&) = default;
  /** @brief Transfer storage/indexes; invalidates borrowers of both trees. */
  DFUDSTree& operator=(DFUDSTree&&) = default;

  /** @brief Reject dangling views of temporary vectors. */
  DFUDSTree(std::vector<uint64_t>&&, size_t)
    requires std::same_as<EncodingStorage, ReadOnlyStorageView>
  = delete;
  /** @brief Reject dangling views of temporary const vectors. */
  DFUDSTree(const std::vector<uint64_t>&&, size_t)
    requires std::same_as<EncodingStorage, ReadOnlyStorageView>
  = delete;

  /** @brief Implement TreeBase::encoding_words. */
  std::span<const uint64_t> encoding_words_impl() const {
    return detail::tree_encoding_words(encoding_, num_bits_);
  }
  /** @brief Implement TreeBase::encoding_owned_bytes. */
  size_t encoding_owned_bytes_impl() const {
    return detail::tree_encoding_owned_bytes(encoding_);
  }

  /**
   * @brief Returns the root node
   */
  Node root_impl() const { return Node(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (num_bits_ + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const { return rmm_.bit(node.pos) == 0; }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const Node& node) const {
    return rmm_.select0(node.number + 1) - node.pos;
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = rmm_.select0(node.number + 1);
    size_t num = node.number + 1;
    return Node(num, pos + 1);
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    auto parent = NodeView(*this, node);
    const auto count = parent.degree();
    if (count == 0) {
      return;
    }
    auto child = NodeView(*this, parent.first_child());
    for (size_t i = 0; i < count; ++i) {
      visitor(child.node());
      if (i + 1 < count) {
        child.reset(child.next_sibling());
      }
    }
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    // A leaf occupies exactly one zero and the next record is its sibling.
    if (is_leaf_impl(node)) {
      return Node(node.number + 1, node.pos + 1);
    }
    size_t end = rmm_.fwdsearch(node.pos, -1);
    size_t pos = end + 1;
    // Preorder numbering advances by the subtree's node count.
    size_t num = node.number + (pos - node.pos + 1) / 2;
    return Node(num, pos);
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_open(node, incoming_open(node));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    return starts_degree_run(incoming_open(node));
  }
};

}  // namespace pixie
