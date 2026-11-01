#pragma once

#include <pixie/memory_usage.h>
#include <pixie/rank_select/support.h>
#include <pixie/tree.h>
#include <pixie/tree/detail/encoding_storage.h>
#include <pixie/tree/detail/louds_children.h>
#include <pixie/tree/detail/louds_node_view.h>
#include <pixie/tree/detail/serialization.h>

#include <bit>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief A node class of LOUDS tree
 */
using LoudsNode = TreeNode;

/**
 * @brief A tree class based on the level order unary degree sequence (LOUDS)
 * representation
 * @tparam RankSelect Rank/select backend constructible from packed source words
 * and a bit count, with both select directions enabled by default. It must
 * provide size(), operator[], select() and select0() with RankSelectBase
 * semantics. EncodingStorage owns the source by default; the TreeView
 * specialization borrows immutable source words.
 * @details Each BFS-ordered node is encoded as one followed by one zero per
 * child. A handle stores its zero-based BFS number and the position of its one.
 * Siblings therefore occupy consecutive one delimiters. The parent delimiter
 * is the nearest preceding one before the node's incoming zero edge. Degree,
 * sibling and parent navigation use a single-word search for those delimiters,
 * falling back to indexed select when the delimiter lies outside that word.
 */
template <class RankSelect = RankSelectSupport<>,
          detail::TreeEncodingStorage EncodingStorage = AlignedStorage>
class LoudsTree;

template <class B, detail::TreeEncodingStorage S>
struct TreeTraits<LoudsTree<B, S>> {
  using TreeView = LoudsTree<B, ReadOnlyStorageView>;
  using NodeView = detail::LoudsNodeView<LoudsTree<B, S>>;
};

template <class RankSelect, detail::TreeEncodingStorage EncodingStorage>
class LoudsTree : public TreeBase<LoudsTree<RankSelect, EncodingStorage>> {
 private:
  LoudsTree(EncodingStorage storage, size_t bits, RankSelect backend)
      : encoding_(std::move(storage)),
        bv(std::move(backend)),
        words_(detail::tree_encoding_words(encoding_, bits).data()) {}

  EncodingStorage encoding_;
  RankSelect bv;
  const std::uint64_t* words_;

  // Node delimiters are ones. Search only the remainder of this word so a
  // high-degree node cannot turn navigation into an unbounded linear scan.
  size_t next_delimiter(size_t position, size_t rank) const {
    const size_t start = position + 1;
    const auto word = words_[start / 64] >> (start % 64);
    return word ? start + std::countr_zero(word) : bv.select(rank);
  }

  size_t parent_delimiter(size_t zero_position, size_t rank) const {
    // zero_position itself is zero, so retaining it avoids a shift by 64.
    const auto word = words_[zero_position / 64] &
                      (~std::uint64_t{0} >> (63 - zero_position % 64));
    return word ? (zero_position / 64) * 64 + 63 - std::countl_zero(word)
                : bv.select(rank);
  }

  LoudsNode parent_from_edge(const LoudsNode& node, size_t zero_pos) const {
    const size_t parent_number = zero_pos - node.number;
    return LoudsNode(parent_number,
                     parent_delimiter(zero_pos, parent_number + 1));
  }

 public:
  /** @brief Write one tree artifact: encoding, node count and nested backend.
   * @details External root indexes and active traversal state are not included.
   */
  void serialize_impl(BinaryWriter& writer) const
    requires Serializable<RankSelect>
  {
    const auto begin = detail::begin_frame(
        writer,
        detail::tree_serialization_tag(detail::TreeEncodingKind::kLOUDS));
    writer.write_size(this->size());
    detail::write_tree_encoding(
        writer, encoding_words_impl(),
        detail::tree_encoding_bits(this->size(), false));
    bv.serialize(writer);
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
  static LoudsTree deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires Deserializable<RankSelect,
                            std::span<const uint64_t>,
                            DeserializationValidation>
  {
    auto candidate = reader;
    auto payload = detail::read_frame(
        candidate,
        detail::tree_serialization_tag(detail::TreeEncodingKind::kLOUDS));
    const auto nodes = payload.read_size();
    const auto bits = detail::tree_encoding_bits(nodes, false);
    auto encoding = detail::read_tree_encoding<EncodingStorage>(payload, bits);
    const auto words = detail::tree_encoding_words(encoding, bits);
    auto backend = RankSelect::deserialize(payload, words, validation);
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
    const auto ones = backend.rank(bits);
    if (ones != nodes) {
      throw std::invalid_argument("Tree backend population mismatch");
    }
    if (validation == DeserializationValidation::kFull) {
      detail::validate_tree_topology(words, nodes,
                                     detail::TreeEncodingKind::kLOUDS);
    }
    payload.require_zero_padding(7);
    LoudsTree result(std::move(encoding), bits, std::move(backend));
    reader = candidate;
    return result;
  }

  using typename TreeBase<LoudsTree>::NodeView;

 private:
  friend class detail::LoudsNodeView<LoudsTree>;

 public:
  /** @brief Implement TreeBase::navigate without mutating the tree. */
  NodeView navigate_impl(TreeNode node) const { return NodeView(*this, node); }

  /** @brief Object, owned encoding and backend allocations; excludes borrowed
   * words. */
  size_t memory_usage_bytes_impl() const
    requires HasMemoryUsageBytes<RankSelect>
  {
    return sizeof(*this) + encoding_owned_bytes_impl() +
           nested_owned_memory_bytes(bv);
  }

  /**
   * @brief Construct from a valid encoding of tree_size nodes.
   * @details Owning trees copy complete source words; TreeView borrows them.
   * A view's source must stay alive, immutable and at a stable address.
   * Invalid topology is undefined behavior. A short source throws
   * std::invalid_argument; an unrepresentable bit count throws
   * std::length_error.
   */
  explicit LoudsTree(std::span<const uint64_t> words, size_t tree_size)
      : LoudsTree(detail::tree_encoding_from_words<EncodingStorage>(
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
  explicit LoudsTree(EncodingStorage storage, size_t tree_size)
      : encoding_(std::move(storage)),
        bv(encoding_words_for_size(tree_size),
           detail::tree_encoding_bits(tree_size, false)),
        words_(encoding_words_for_size(tree_size).data()) {}

  /** @brief Deep-copy owning encodings and rebuild their backend references.
   * @details TreeView copies backend metadata and borrows the same encoding.
   */
  LoudsTree(const LoudsTree& other)
      : encoding_(other.encoding_),
        bv([&] {
          if constexpr (std::same_as<EncodingStorage, AlignedStorage>) {
            return RankSelect(encoding_words_for_size(other.size()),
                              other.bv.size());
          } else {
            return other.bv;
          }
        }()),
        words_(encoding_words_for_size(other.size()).data()) {}
  /** @brief Copy independently; invalidates existing borrowers of *this. */
  LoudsTree& operator=(const LoudsTree& other) {
    if (this != &other) {
      LoudsTree copy(other);
      *this = std::move(copy);
    }
    return *this;
  }
  /** @brief Transfer storage/indexes; invalidates borrowers of the source.
   * @details The moved-from object may only be destroyed or assigned.
   */
  LoudsTree(LoudsTree&&) = default;
  /** @brief Transfer storage/indexes; invalidates borrowers of both trees. */
  LoudsTree& operator=(LoudsTree&&) = default;

  /** @brief Reject dangling views of temporary vectors. */
  LoudsTree(std::vector<uint64_t>&&, size_t)
    requires std::same_as<EncodingStorage, ReadOnlyStorageView>
  = delete;
  /** @brief Reject dangling views of temporary const vectors. */
  LoudsTree(const std::vector<uint64_t>&&, size_t)
    requires std::same_as<EncodingStorage, ReadOnlyStorageView>
  = delete;

  /** @brief Implement TreeBase::encoding_words. */
  std::span<const uint64_t> encoding_words_impl() const {
    return detail::tree_encoding_words(encoding_, bv.size());
  }
  /** @brief Implement TreeBase::encoding_owned_bytes. */
  size_t encoding_owned_bytes_impl() const {
    return detail::tree_encoding_owned_bytes(encoding_);
  }

 private:
  std::span<const uint64_t> encoding_words_for_size(size_t nodes) const {
    return detail::tree_encoding_words(
        encoding_, detail::tree_encoding_bits(nodes, false));
  }

 public:
  /**
   * @brief Returns the root node
   */
  LoudsNode root_impl() const { return LoudsNode(0, 0); }

  /**
   * @brief Returns the size of the tree
   */
  size_t size_impl() const { return (bv.size() + 1) / 2; }

  /**
   * @brief Indicates if @p node is a leaf
   */
  bool is_leaf_impl(const LoudsNode& node) const {
    return (node.pos + 1 == bv.size()) or bv[node.pos + 1];
  }

  /**
   * @brief Indicates if @p node is a root
   */
  bool is_root_impl(const LoudsNode& node) const { return node.number == 0; }

  /**
   * @brief Returns the number of children of a @p node
   */
  size_t degree_impl(const LoudsNode& node) const {
    if (is_leaf_impl(node)) {
      return 0;
    }
    return next_delimiter(node.pos, node.number + 2) - node.pos - 1;
  }

  /** @brief Implement TreeBase::for_each_child with local traversal state. */
  template <class Visitor>
  void for_each_child_impl(const TreeNode& node, Visitor&& visitor) const {
    detail::tree_unpack::louds_children(
        *this,
        std::span<const uint64_t>(words_,
                                  bv.size() / 64 + (bv.size() % 64 != 0)),
        bv.size(), node, 0, degree_impl(node), visitor);
  }

  /**
   * @brief Returns the i-th child of @p node
   * Indexing starts at 0
   */
  LoudsNode child_impl(const LoudsNode& node, size_t i) const {
    size_t zeros = node.pos + i + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns first child of a @p node
   */
  LoudsNode first_child_impl(const LoudsNode& node) const {
    size_t zeros = node.pos + 1 - node.number;
    return LoudsNode(zeros, bv.select(zeros + 1));
  }

  /**
   * @brief Returns the parent of a @p node if @p node is not root,
   * else returns root
   */
  LoudsNode parent_impl(const LoudsNode& node) const {
    if (node.number == 0) {
      return root_impl();
    }
    return parent_from_edge(node, bv.select0(node.number));
  }

  /**
   * @brief Indicates if @p node is last child
   */
  bool is_last_child_impl(const LoudsNode& node) const {
    size_t zero_pos = bv.select0(node.number);
    return bv[zero_pos + 1];
  }

  /**
   * @brief Returns next sibling of a @p node
   */
  LoudsNode next_sibling_impl(const LoudsNode& node) const {
    size_t sibling_number = node.number + 1;
    return LoudsNode(sibling_number,
                     next_delimiter(node.pos, sibling_number + 1));
  }
};

}  // namespace pixie
