#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rmm/btree.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/bp_node_view.h>
#include <pixie/tree/detail/encoding_storage.h>
#include <pixie/tree/detail/serialization.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief A tree class based on the balances parentheses (BP)
 * representation
 */
template <class RMMTree,
          detail::TreeEncodingStorage EncodingStorage = AlignedStorage>
class BPTree;

template <class B, detail::TreeEncodingStorage S>
struct TreeTraits<BPTree<B, S>> {
  using TreeView = BPTree<B, ReadOnlyStorageView>;
  using NodeView = detail::BpNodeView<BPTree<B, S>>;
};

template <class RMMTree, detail::TreeEncodingStorage EncodingStorage>
class BPTree : public TreeBase<BPTree<RMMTree, EncodingStorage>> {
 private:
  BPTree(EncodingStorage storage, size_t bits, RMMTree backend)
      : encoding_(std::move(storage)),
        num_bits_(bits),
        rmm_(std::move(backend)) {}

  EncodingStorage encoding_;
  size_t num_bits_;
  RMMTree rmm_;

 public:
  /** @brief Write one tree artifact: encoding, node count and nested backend.
   * @details External root indexes and active traversal state are not included.
   */
  void serialize_impl(BinaryWriter& writer) const
    requires Serializable<RMMTree>
  {
    const auto begin = detail::begin_frame(
        writer, detail::tree_serialization_tag(detail::TreeEncodingKind::kBP));
    writer.write_size(this->size());
    detail::write_tree_encoding(writer, encoding_words_impl(),
                                detail::tree_encoding_bits(this->size(), true));
    rmm_.serialize(writer);
    detail::end_frame(writer, begin);
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
  static BPTree deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires Deserializable<RMMTree,
                            std::span<const uint64_t>,
                            DeserializationValidation>
  {
    auto candidate = reader;
    auto payload = detail::read_frame(
        candidate,
        detail::tree_serialization_tag(detail::TreeEncodingKind::kBP));
    const auto nodes = payload.read_size();
    const auto bits = detail::tree_encoding_bits(nodes, true);
    auto encoding = detail::read_tree_encoding<EncodingStorage>(payload, bits);
    const auto words = detail::tree_encoding_words(encoding, bits);
    auto backend = RMMTree::deserialize(payload, words, validation);
    if (backend.size() != bits) {
      throw std::invalid_argument("Tree backend size mismatch");
    }
    if constexpr (requires {
                    backend.supports_select1();
                    backend.supports_select0();
                  }) {
      if (!backend.supports_select1() || !backend.supports_select0()) {
        throw std::invalid_argument(
            "Tree backend requires both select directions");
      }
    }
    const auto ones = backend.rank1(bits);
    if (ones != nodes) {
      throw std::invalid_argument("Tree backend population mismatch");
    }
    if (validation == DeserializationValidation::kFull) {
      detail::validate_tree_topology(words, nodes,
                                     detail::TreeEncodingKind::kBP);
    }
    payload.require_zero_padding(7);
    BPTree result(std::move(encoding), bits, std::move(backend));
    reader = candidate;
    return result;
  }

  using typename TreeBase<BPTree>::NodeView;

 private:
  friend class detail::BpNodeView<BPTree>;

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
   * @brief A node class of BP tree
   */
  using Node = TreeNode;

  /**
   * @brief Construct from a valid encoding of tree_size nodes.
   * @details Owning trees copy complete source words; TreeView borrows them.
   * A view's source must stay alive, immutable and at a stable address.
   * Invalid topology is undefined behavior. A short source throws
   * std::invalid_argument; an unrepresentable bit count throws
   * std::length_error.
   */
  explicit BPTree(std::span<const uint64_t> words, size_t tree_size)
      : BPTree(detail::tree_encoding_from_words<EncodingStorage>(
                   words,
                   detail::tree_encoding_bits(tree_size, true)),
               tree_size) {}

  /**
   * @brief Take storage by value: moves transfer ownership, lvalues copy.
   * @details ReadOnlyStorageView always borrows. AlignedStorage may have a
   * partial final logical byte/word; its allocation supplies word padding.
   * View storage must contain aligned complete words. Short/misaligned storage
   * throws std::invalid_argument. Topology and size must agree.
   */
  explicit BPTree(EncodingStorage storage, size_t tree_size)
      : encoding_(std::move(storage)),
        num_bits_(detail::tree_encoding_bits(tree_size, true)),
        rmm_(encoding_words_impl(), num_bits_) {}

  /** @brief Deep-copy owning encodings and rebuild their backend references.
   * @details TreeView copies backend metadata and borrows the same encoding.
   */
  BPTree(const BPTree& other)
      : encoding_(other.encoding_), num_bits_(other.num_bits_), rmm_([&] {
          if constexpr (std::same_as<EncodingStorage, AlignedStorage>) {
            return RMMTree(encoding_words_impl(), num_bits_);
          } else {
            return other.rmm_;
          }
        }()) {}
  /** @brief Copy independently; invalidates existing borrowers of *this. */
  BPTree& operator=(const BPTree& other) {
    if (this != &other) {
      BPTree copy(other);
      *this = std::move(copy);
    }
    return *this;
  }
  /** @brief Transfer storage/indexes; invalidates borrowers of the source.
   * @details The moved-from object may only be destroyed or assigned.
   */
  BPTree(BPTree&&) = default;
  /** @brief Transfer storage/indexes; invalidates borrowers of both trees. */
  BPTree& operator=(BPTree&&) = default;

  /** @brief Reject dangling views of temporary vectors. */
  BPTree(std::vector<uint64_t>&&, size_t)
    requires std::same_as<EncodingStorage, ReadOnlyStorageView>
  = delete;
  /** @brief Reject dangling views of temporary const vectors. */
  BPTree(const std::vector<uint64_t>&&, size_t)
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
  size_t size_impl() const { return num_bits_ / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const Node& node) const {
    return rmm_.bit(node.pos + 1) == 0;
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const Node& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   *    this method has O(d) time complexity!
   *
   *    TODO try make this faster
   */
  size_t degree_impl(const Node& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    Node child = first_child_impl(node);
    size_t child_count = 1;
    while (true) {
      if (is_last_child_impl(child)) {
        return child_count;
      }
      child = next_sibling_impl(child);
      child_count++;
    }
  }

  /**
   * @brief Returns first child of a @p node
   */
  Node first_child_impl(const Node& node) const {
    size_t pos = node.pos + 1;
    size_t num = node.number + 1;
    return Node(num, pos);
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    if (is_leaf_impl(node)) {
      return;
    }
    auto child = NodeView(*this, first_child_impl(node));
    while (true) {
      visitor(child.node());
      if (child.is_last_child()) {
        return;
      }
      child.reset(child.next_sibling());
    }
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   *    this method has O(i) time complexity!
   *
   *    TODO try make this faster
   */
  Node child_impl(const Node& node, size_t i) const {
    Node child = first_child_impl(node);
    while (i--) {
      child = next_sibling_impl(child);
    }
    return child;
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  Node parent_impl(const Node& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    // An opening immediately before this node identifies its parent. It is
    // also the preceding preorder node, so neither enclose nor rank is needed.
    if (rmm_.bit(node.pos - 1)) {
      return Node(node.number - 1, node.pos - 1);
    }
    size_t pos = rmm_.enclose(node.pos);
    // Between the parent's open and this child are complete sibling subtrees,
    // each occupying two bits per node, plus the parent's own opening bit.
    size_t num = node.number - (node.pos - pos + 1) / 2;
    return Node(num, pos);
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const Node& node) const {
    const size_t end = is_leaf_impl(node) ? node.pos + 1 : rmm_.close(node.pos);
    // A non-root node's close is followed by a sibling open or a parent close.
    return rmm_.bit(end + 1) == 0;
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  Node next_sibling_impl(const Node& node) const {
    // A leaf is exactly "10"; its sibling is the next preorder node.
    if (is_leaf_impl(node)) {
      return Node(node.number + 1, node.pos + 2);
    }
    size_t pos = rmm_.close(node.pos) + 1;
    size_t num = node.number + (pos - node.pos) / 2;
    return Node(num, pos);
  }
};

std::vector<uint64_t> adj_to_bp(size_t tree_size,
                                const std::vector<std::vector<size_t>>& adj) {
  size_t bp_size = tree_size * 2;
  std::vector<uint64_t> bp((bp_size + 63) / 64, 0);
  std::vector<std::pair<size_t, size_t>> stack;
  stack.push_back(std::make_pair(0, 0));
  size_t pos = 0;
  bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
  while (!stack.empty()) {
    auto& [v, p] = stack.back();
    p++;
    if (p >= adj[v].size()) {
      pos++;
      stack.pop_back();
      continue;
    }
    pos++;
    bp[pos >> 6] = bp[pos >> 6] | (1ULL << (pos & 63));
    stack.push_back(std::make_pair(adj[v][p], 0));
  }
  return bp;
}
}  // namespace pixie
