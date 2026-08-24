#pragma once

#include <pixie/detail/serialization.h>
#include <pixie/packed_bit_builder.h>
#include <pixie/rank_select/support.h>
#include <pixie/wavelet_tree.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <span>
#include <vector>

namespace pixie {

template <StorageImplementation Storage = AlignedStorage>
class WaveletTreeIndex : public WaveletTreeBase<WaveletTreeIndex<Storage>> {
 private:
  using node_index_t = size_t;
  static constexpr node_index_t npos = std::numeric_limits<node_index_t>::max();
  static constexpr std::array<std::uint8_t, 8> kSerializationMagic = {
      'P', 'X', 'W', 'A', 'V', 'E', 'T', '\0'};
  static constexpr std::uint32_t kSerializationVersion = 4;
  static constexpr std::size_t kSerializationHeaderBytes = 24;

  struct PreWaveletNode {
    node_index_t parent = npos;
    node_index_t left_child = npos;
    node_index_t right_child = npos;
    uint64_t middle;
    PackedBitBuilder stream;
    explicit PreWaveletNode(uint64_t middle) : middle(middle) {}
  };

  /**
   * @brief Node of the wavelet tree
   * @details
   * Node with its bitvector representing division of characters corresponding
   * to this node: 0 - left, 1 - right, in their original order. It also has
   * indices of its children and parent for routing both up-down and bottom-up.
   *
   */
  struct WaveletNode {
    node_index_t parent, left_child, right_child;
    uint64_t middle;
    Storage bit_vector_data;
    RankSelectSupport<Storage> data;

    /** @brief Manually turns std::vector<uint64_t> into AlignedStorage */
    static AlignedStorage align(std::vector<uint64_t>&& data) {
      AlignedStorage result(data.size() * 64);
      auto view = result.writable_words64();
      std::copy(data.begin(), data.end(), view.begin());
      return result;
    }

    WaveletNode() = default;

    WaveletNode(PreWaveletNode&& node)
      requires(std::same_as<Storage, AlignedStorage>)
        : parent(node.parent),
          left_child(node.left_child),
          right_child(node.right_child),
          middle(node.middle) {
      const std::size_t bit_count = node.stream.size_bits();
      bit_vector_data = align(node.stream.take_words());
      data =
          RankSelectSupport<Storage>(bit_vector_data.as_words64(), bit_count);
    }

    /** @brief Write one node in canonical little-endian form. */
    void serialize(BinaryWriter& writer) const {
      writer.write_size(parent);
      writer.write_size(left_child);
      writer.write_size(right_child);
      writer.write_u64(middle);
      bit_vector_data.serialize(writer);
      data.serialize(writer);
    }

    /** @brief Construct one checked, non-owning node from @p reader. */
    static WaveletNode deserialize(BinaryReader& reader,
                                   DeserializationValidation validation)
      requires(std::same_as<Storage, ReadOnlyStorageView>)
    {
      WaveletNode result;
      result.parent = reader.read_size();
      result.left_child = reader.read_size();
      result.right_child = reader.read_size();
      result.middle = reader.read_u64();
      result.bit_vector_data = ReadOnlyStorageView::deserialize(reader);
      result.data = RankSelectSupport<Storage>::deserialize(
          result.bit_vector_data.as_words64(), reader, validation);
      return result;
    }
  };

  size_t alphabet_size_, data_size_;
  node_index_t root_;
  std::vector<WaveletNode> nodes_;
  std::vector<node_index_t> leaves_;
  std::vector<size_t> permutation_, inverse_permutation_;

  void validate_deserialized_topology(
      DeserializationValidation validation) const {
    if (root_ == npos) {
      if (validation == DeserializationValidation::kFull &&
          std::ranges::any_of(leaves_,
                              [](node_index_t leaf) { return leaf != npos; })) {
        throw std::invalid_argument(
            "Invalid serialized empty wavelet-tree leaves");
      }
      return;
    }
    if (nodes_[root_].parent != npos) {
      throw std::invalid_argument("Serialized wavelet-tree root has a parent");
    }

    std::vector<std::uint8_t> incoming_edges(nodes_.size());
    for (node_index_t parent = 0; parent < nodes_.size(); ++parent) {
      const WaveletNode& node = nodes_[parent];
      for (const node_index_t child : {node.left_child, node.right_child}) {
        if (child == npos) {
          continue;
        }
        if (nodes_[child].parent != parent) {
          throw std::invalid_argument(
              "Serialized wavelet-tree parent/child links disagree");
        }
        if (incoming_edges[child] != 0) {
          throw std::invalid_argument(
              "Serialized wavelet-tree node has multiple parents");
        }
        ++incoming_edges[child];
      }
    }

    for (node_index_t node = 0; node < nodes_.size(); ++node) {
      const std::size_t expected_edges = node == root_ ? 0 : 1;
      if (incoming_edges[node] != expected_edges) {
        throw std::invalid_argument(
            "Serialized wavelet-tree node is detached from its parent");
      }
    }

    struct PendingNode {
      node_index_t node;
      std::size_t symbol_begin;
      std::size_t symbol_end;
    };
    std::vector<bool> reached(nodes_.size());
    std::vector<PendingNode> pending = {{root_, 0, alphabet_size_}};
    while (!pending.empty()) {
      const PendingNode current = pending.back();
      pending.pop_back();
      const node_index_t node = current.node;
      reached[node] = true;
      const WaveletNode& metadata = nodes_[node];
      if (metadata.middle <= current.symbol_begin ||
          metadata.middle >= current.symbol_end) {
        throw std::invalid_argument(
            "Serialized wavelet-tree split is outside its symbol range");
      }

      const std::size_t one_count =
          validation == DeserializationValidation::kFull
              ? metadata.data.rank(metadata.data.size())
              : 0;
      const std::size_t zero_count =
          validation == DeserializationValidation::kFull
              ? metadata.data.size() - one_count
              : 0;
      const auto validate_branch = [&](node_index_t child,
                                       std::size_t symbol_begin,
                                       std::size_t symbol_end,
                                       std::size_t expected_size) {
        if (child != npos) {
          if (validation == DeserializationValidation::kFull &&
              nodes_[child].data.size() != expected_size) {
            throw std::invalid_argument(
                "Serialized wavelet-tree child has the wrong length");
          }
          pending.push_back({child, symbol_begin, symbol_end});
          return;
        }
        if (validation == DeserializationValidation::kFull) {
          for (std::size_t symbol = symbol_begin; symbol < symbol_end;
               ++symbol) {
            if (leaves_[symbol] != node) {
              throw std::invalid_argument(
                  "Serialized wavelet-tree leaf map disagrees with topology");
            }
          }
        }
      };
      validate_branch(metadata.left_child, current.symbol_begin,
                      metadata.middle, zero_count);
      validate_branch(metadata.right_child, metadata.middle, current.symbol_end,
                      one_count);
    }
    if (std::ranges::find(reached, false) != reached.end()) {
      throw std::invalid_argument(
          "Serialized wavelet-tree contains unreachable nodes");
    }
  }

  /**
   * @brief Recursive building of the nodes
   *
   * @tparam F get_middle function typename
   * @param begin Begin of the characters segment
   * @param end End of the characters segment
   * @param parent Node index of the parent node
   * @param get_middle Offset of separating cut from the segment beginning if
   * there is precomputed tree structure and npos otherwise
   * @param prefix_sum Span with i-th value equal to number of characters less
   * than i in the original data
   * @param nodes Resulting vector of PreWaveletNodes
   * @return Index of the node built
   *
   */
  template <typename F>
  node_index_t build_node(size_t begin,
                          size_t end,
                          node_index_t parent,
                          const F& get_middle,
                          std::span<const size_t> prefix_sum,
                          std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    if (end - begin == 1) {
      leaves_[begin] = parent;
      return npos;
    }
    if (prefix_sum[end] == prefix_sum[begin]) {
      for (size_t symbol = begin; symbol < end; symbol++) {
        leaves_[symbol] = parent;
      }
      return npos;
    }

    node_index_t result = nodes.size();
    size_t middle = get_middle(result);
    middle = begin + (middle == npos ? (end - begin) / 2 : middle);

    nodes.emplace_back(middle);
    nodes[result].stream.reserve_bits(prefix_sum[end] - prefix_sum[begin]);
    nodes[result].parent = parent;
    nodes[result].left_child =
        build_node(begin, middle, result, get_middle, prefix_sum, nodes);
    nodes[result].right_child =
        build_node(middle, end, result, get_middle, prefix_sum, nodes);

    return result;
  }

  /**
   * @brief Recursively copies segment of the original data corresponding to the
   * node
   *
   * @param node Index of the node
   * @param begin Begin of the segment
   * @param end End of the segment
   * @param dst Destination of the copy
   * @param tmp Temporary memory buffer
   *
   * @details
   * Copies left and right childs segments to the tmp and then builds up target
   * segment in the dst
   *
   */
  void copy_segment_content(node_index_t node,
                            size_t begin,
                            size_t end,
                            std::span<uint64_t> dst,
                            std::span<uint64_t> tmp) const {
    if (begin == end) {
      return;
    }
    const size_t rank = nodes_[node].data.rank(begin), rank0 = begin - rank;
    const size_t right = nodes_[node].data.rank(end) - rank,
                 left = (end - begin) - right;

    if (nodes_[node].left_child == npos) {
      std::fill_n(tmp.begin(), static_cast<long long>(left),
                  inverse_permutation_[nodes_[node].middle - 1]);
    } else {
      copy_segment_content(nodes_[node].left_child, rank0, rank0 + left,
                           tmp.subspan(0, left), dst.subspan(0, left));
    }
    if (nodes_[node].right_child == npos) {
      std::fill(tmp.begin() + static_cast<long long>(left), tmp.end(),
                inverse_permutation_[nodes_[node].middle]);
    } else {
      copy_segment_content(nodes_[node].right_child, rank, rank + right,
                           tmp.subspan(left, right), dst.subspan(left, right));
    }

    size_t j = 0, k = left;
    const auto& bit_vector = nodes_[node].bit_vector_data.as_words64();
    for (size_t i = begin; i < end; i++) {
      if ((bit_vector[i / 64] >> (i % 64)) & 1) {
        dst[i - begin] = tmp[k++];
      } else {
        dst[i - begin] = tmp[j++];
      }
    }
  }

  WaveletTreeIndex() = default;

 public:
  /**
   * @param alphabet_size Size of the alphabet
   * @param data Original text. Its characters are from the
   * range [0, alphabet_size)
   * @param build_type Either Standard or Huffman. This effects on how the
   * wavelet tree builds: like segment tree on trivially sorted characters or
   * like in Huffman algorithm
   *
   * @details
   * Standard: Just calls build_node
   * Huffman: Reorders characters with respect to Huffman algorithm and then
   * calls build_node with specific get_middle function
   *
   */
  WaveletTreeIndex(
      size_t alphabet_size,
      std::span<const uint64_t> data,
      const WaveletTreeBuildType build_type = WaveletTreeBuildType::Standard)
    requires(std::same_as<Storage, AlignedStorage>)
      : alphabet_size_(alphabet_size),
        data_size_(data.size()),
        leaves_(alphabet_size_, npos) {
    if (alphabet_size == 0) {
      root_ = npos;
      return;
    }
    std::vector<PreWaveletNode> nodes;
    nodes.reserve(alphabet_size_);
    std::vector<size_t> nodes_structure;
    nodes_structure.reserve(alphabet_size_);

    if (build_type == WaveletTreeBuildType::Standard) {
      permutation_.resize(alphabet_size);
      inverse_permutation_.resize(alphabet_size);
      std::iota(permutation_.begin(), permutation_.end(), 0);
      std::iota(inverse_permutation_.begin(), inverse_permutation_.end(), 0);
      nodes_structure.resize(alphabet_size_, npos);
    } else {
      struct Node {
        size_t size, left, right;
      };
      std::vector<Node> huffman_nodes(alphabet_size_, {0, 0, 0});
      for (auto symb : data) {
        huffman_nodes[symb].size++;
      }

      using elem_t = std::pair<size_t, size_t>;
      std::priority_queue<elem_t, std::vector<elem_t>, std::greater<>> queue;
      for (size_t i = 0; i < alphabet_size_; i++) {
        queue.emplace(huffman_nodes[i].size, i);
      }
      while (queue.size() >= 2) {
        auto right = queue.top().second;
        queue.pop();
        auto left = queue.top().second;
        queue.pop();
        huffman_nodes.push_back(
            {huffman_nodes[left].size + huffman_nodes[right].size, left + 1,
             right + 1});
        queue.emplace(huffman_nodes.back().size, huffman_nodes.size() - 1);
      }

      std::function<size_t(size_t)> enumerate = [&](size_t index) -> size_t {
        const auto& [size, left, right] = huffman_nodes[index];
        if (left == 0 || right == 0) {
          permutation_[index] = inverse_permutation_.size();
          inverse_permutation_.push_back(index);
          return 1;
        }
        size_t ind = nodes_structure.size(), subtree = 0;
        if (size > 0) {
          nodes_structure.push_back(0);
        }
        subtree += enumerate(left - 1);
        if (size > 0) {
          nodes_structure[ind] = subtree;
        }
        subtree += enumerate(right - 1);
        return subtree;
      };

      permutation_.resize(alphabet_size_);
      inverse_permutation_.reserve(alphabet_size_);
      enumerate(huffman_nodes.size() - 1);
    }

    std::vector<size_t> prefix_sum(alphabet_size + 1);
    for (auto symbol : data) {
      prefix_sum[permutation_[symbol] + 1]++;
    }
    for (size_t i = 0; i < alphabet_size_; i++) {
      prefix_sum[i + 1] += prefix_sum[i];
    }

    root_ = build_node(
        0, alphabet_size_, npos,
        [&](node_index_t node) { return nodes_structure[node]; }, prefix_sum,
        nodes);
    for (auto symbol : data) {
      auto index = permutation_[symbol];
      for (node_index_t current = root_; current != npos;) {
        auto& node = nodes[current];
        bool go_right = index >= node.middle;
        node.stream.write_bit(go_right);
        if (go_right) {
          current = node.right_child;
        } else {
          current = node.left_child;
        }
      }
    }
    nodes_.reserve(nodes.size());
    for (auto& node : nodes) {
      nodes_.emplace_back(std::move(node));
    }
  }

  /**
   * @brief Rank of specified symbol up to position pos (exclusive)
   *
   * @param symbol The character that the query is about
   * @param pos Character index in [0, size()]
   * @return Number of specified symbols in [0, pos)
   *
   */
  size_t rank_impl(uint64_t symbol, size_t pos) const {
    if (symbol >= alphabet_size_) [[unlikely]] {
      return 0;
    }
    symbol = permutation_[symbol];
    for (node_index_t current = root_; current != npos;) {
      const WaveletNode& node = nodes_[current];
      if (symbol < node.middle) {
        pos = node.data.rank0(pos);
        current = node.left_child;
      } else {
        pos = node.data.rank(pos);
        current = node.right_child;
      }
    }
    return pos;
  }

  /**
   * @brief Select the position of the rank-th specified symbol (1-indexed)
   *
   * @param symbol The character that the query is about
   * @param rank 1-indexed rank of specified symbol to select
   * @return Symbol index, or size() if rank is out of range
   *
   */
  size_t select_impl(uint64_t symbol, size_t rank) const {
    if (symbol >= alphabet_size_ || data_size_ == 0) [[unlikely]] {
      return data_size_;
    }
    symbol = permutation_[symbol];
    node_index_t current = leaves_[symbol];
    for (; current != npos; current = nodes_[current].parent) {
      const WaveletNode& node = nodes_[current];
      if (symbol < node.middle) {
        rank = node.data.select0(rank) + 1;
      } else {
        rank = node.data.select(rank) + 1;
      }
    }
    return rank - 1;
  }

  /**
   * @brief Accumulates the original data segment
   *
   * @param begin Begin of the segment
   * @param end End of the segment
   * @return Queried segment of data
   *
   * @details Queries packed bit vectors and rank/select metadata directly
   * through the node storage. A deserialized view does not consult or retain
   * its `BinaryReader`. The current implementation materializes the requested
   * output and an equally sized temporary buffer, for peak auxiliary and
   * result storage of two `uint64_t` values per returned symbol.
   *
   */
  std::vector<uint64_t> get_segment_impl(size_t begin, size_t end) const {
    if (alphabet_size_ == 0 || data_size_ == 0 || begin >= end) [[unlikely]] {
      return {};
    }
    auto length = static_cast<long long>(end - begin);
    std::vector<uint64_t> result(2 * length);
    copy_segment_content(root_, begin, end,
                         std::span{result.begin(), result.begin() + length},
                         std::span{result.begin() + length, result.end()});
    result.resize(length);
    return result;
  }

  /**
   * @return Returns the number of characters in data
   *
   */
  size_t size_impl() const { return data_size_; }

  /**
   * @brief Write a versioned canonical little-endian wavelet-tree artifact.
   *
   * @throws std::invalid_argument if the artifact would not begin at an
   * eight-byte-aligned writer offset required by zero-copy deserialization.
   */
  void serialize(BinaryWriter& writer) const {
    if (writer.size_bytes() % alignof(std::uint64_t) != 0) {
      throw std::invalid_argument(
          "Wavelet-tree serialization requires an aligned writer offset");
    }
    const std::size_t artifact_begin = writer.size_bytes();
    detail::write_magic(writer, kSerializationMagic);
    writer.write_u32(kSerializationVersion);
    writer.write_u32(0);
    const std::size_t artifact_size_position = writer.write_u64_placeholder();

    writer.write_size(alphabet_size_);
    writer.write_size(data_size_);
    writer.write_size(root_);
    writer.write_size(nodes_.size());
    for (const WaveletNode& node : nodes_) {
      node.serialize(writer);
    }
    for (const node_index_t leaf : leaves_) {
      writer.write_size(leaf);
    }
    for (const size_t idx : permutation_) {
      writer.write_size(idx);
    }

    const std::size_t unpadded_size = writer.size_bytes() - artifact_begin;
    writer.write_zeros(
        (sizeof(std::uint64_t) - unpadded_size % sizeof(std::uint64_t)) %
        sizeof(std::uint64_t));
    writer.patch_u64(
        artifact_size_position,
        static_cast<std::uint64_t>(writer.size_bytes() - artifact_begin));
  }

  /**
   * @brief Restore one checked, non-owning wavelet-tree artifact.
   *
   * @details The result retains views into the reader's backing bytes. Those
   * bytes must remain alive, immutable, and aligned for 64-bit access for the
   * result's lifetime. On success @p reader advances past exactly one framed
   * artifact; on failure it is unchanged. @p validation selects quick
   * structural checks or exact bitvector-derived metadata validation.
   *
   * @param reader Input cursor, advanced only after successful validation.
   * @param validation Quick structural or full bitvector-derived validation.
   *
   * @throws std::invalid_argument for malformed, truncated, incompatible, or
   * structurally inconsistent metadata.
   * @throws std::length_error when an encoded count is not representable.
   */
  static WaveletTreeIndex<ReadOnlyStorageView> deserialize(
      BinaryReader& reader,
      DeserializationValidation validation =
          DeserializationValidation::kQuick) {
    BinaryReader candidate = reader;
    if (reinterpret_cast<std::uintptr_t>(candidate.remaining_bytes().data()) %
            alignof(std::uint64_t) !=
        0) {
      throw std::invalid_argument(
          "Serialized wavelet-tree artifact is not word aligned");
    }
    const std::size_t available_size = candidate.remaining();
    detail::require_magic(candidate, kSerializationMagic);
    if (candidate.read_u32() != kSerializationVersion ||
        candidate.read_u32() != 0) {
      throw std::invalid_argument(
          "Incompatible serialized wavelet-tree artifact");
    }
    const std::size_t artifact_size = detail::checked_artifact_size(
        candidate.read_u64(), kSerializationHeaderBytes, available_size);
    BinaryReader payload =
        candidate.read_subreader(artifact_size - kSerializationHeaderBytes);

    WaveletTreeIndex<ReadOnlyStorageView> result;
    result.alphabet_size_ = payload.read_size();
    result.data_size_ = payload.read_size();
    result.root_ = payload.read_size();
    const std::size_t node_count = payload.read_size();
    const std::vector<WaveletNode> empty_nodes;
    if (node_count > empty_nodes.max_size()) {
      throw std::length_error(
          "Serialized wavelet-tree node count is too large");
    }
    constexpr std::size_t kMinimumNodeBytes =
        4 * sizeof(std::uint64_t) + sizeof(std::uint64_t);
    if (node_count > payload.remaining() / kMinimumNodeBytes) {
      throw SerializationError("Truncated serialized wavelet-tree nodes",
                               payload.byte_offset());
    }
    result.nodes_.resize(node_count);
    for (auto& node : result.nodes_) {
      node = WaveletNode::deserialize(payload, validation);
    }
    const std::vector<node_index_t> empty_indices;
    if (result.alphabet_size_ > empty_indices.max_size() ||
        result.alphabet_size_ >
            payload.remaining() / (2 * sizeof(std::uint64_t))) {
      throw std::length_error("Serialized wavelet-tree alphabet is too large");
    }
    result.leaves_.resize(result.alphabet_size_);
    for (node_index_t& leaf : result.leaves_) {
      leaf = payload.read_size();
    }
    result.permutation_.resize(result.alphabet_size_);
    for (size_t& index : result.permutation_) {
      index = payload.read_size();
    }
    result.inverse_permutation_.resize(result.alphabet_size_);
    std::vector<bool> seen(result.alphabet_size_);
    for (size_t i = 0; i < result.alphabet_size_; i++) {
      if (result.permutation_[i] >= result.alphabet_size_ ||
          seen[result.permutation_[i]]) {
        throw std::invalid_argument(
            "Invalid serialized wavelet-tree permutation");
      }
      seen[result.permutation_[i]] = true;
      result.inverse_permutation_[result.permutation_[i]] = i;
    }
    const auto valid_node_index = [&result](node_index_t index) {
      return index == npos || index < result.nodes_.size();
    };
    if (!valid_node_index(result.root_) ||
        (result.nodes_.empty() != (result.root_ == npos))) {
      throw std::invalid_argument("Invalid serialized wavelet-tree root");
    }
    for (const node_index_t leaf : result.leaves_) {
      if (!valid_node_index(leaf)) {
        throw std::invalid_argument("Invalid serialized wavelet-tree leaf");
      }
    }
    for (const WaveletNode& node : result.nodes_) {
      if (!valid_node_index(node.parent) ||
          !valid_node_index(node.left_child) ||
          !valid_node_index(node.right_child) ||
          node.data.size() > node.bit_vector_data.size_bits() ||
          node.middle == 0 || node.middle >= result.alphabet_size_) {
        throw std::invalid_argument("Invalid serialized wavelet-tree node");
      }
    }
    result.validate_deserialized_topology(validation);
    if (result.root_ != npos &&
        result.nodes_[result.root_].data.size() != result.data_size_) {
      throw std::invalid_argument(
          "Serialized wavelet-tree root has the wrong length");
    }
    payload.require_zero_padding(sizeof(std::uint64_t) - 1);
    reader = candidate;
    return result;
  }

  /**
   * @brief Restore one artifact from @p data and advance it on success.
   * @param data Input bytes, advanced only after successful validation.
   * @param validation Quick structural or full bitvector-derived validation.
   */
  static WaveletTreeIndex<ReadOnlyStorageView> deserialize(
      std::span<const std::byte>& data,
      DeserializationValidation validation =
          DeserializationValidation::kQuick) {
    BinaryReader reader(data);
    auto result = deserialize(reader, validation);
    data = data.subspan(reader.position());
    return result;
  }
};

using WaveletTree = WaveletTreeIndex<AlignedStorage>;
using WaveletTreeView = WaveletTreeIndex<ReadOnlyStorageView>;

}  // namespace pixie
