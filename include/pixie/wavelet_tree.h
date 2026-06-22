#pragma once

#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <span>
#include <vector>

#include "pixie/bitvector.h"

namespace pixie {

enum class WaveletTreeBuildType { Standard, Huffman };

template <typename Storage = AlignedStorage>
class WaveletTreeBase {
 private:
  using node_index_t = size_t;
  static constexpr node_index_t npos = std::numeric_limits<node_index_t>::max();

  struct PreWaveletNode {
    node_index_t parent = npos;
    node_index_t left_child = npos;
    node_index_t right_child = npos;
    uint64_t middle;
    OutputBitStream stream;
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
    BitVectorBase<Storage> data;

    /** @brief Manually turns std::vector<uint64_t> into AlignedStorage */
    static AlignedStorage align(std::vector<uint64_t>&& data) {
      AlignedStorage result(data.size() * 64);
      auto view = result.As64BitInts();
      std::copy(data.begin(), data.end(), view.begin());
      return std::move(result);
    }

    WaveletNode(PreWaveletNode&& node)
      requires(std::same_as<Storage, AlignedStorage>)
        : parent(node.parent),
          left_child(node.left_child),
          right_child(node.right_child),
          middle(node.middle),
          bit_vector_data(std::move(align(node.stream.extract()))),
          data(bit_vector_data.AsConst64BitInts(), node.stream.size()) {}

    /**
     * @brief Writes a node serialization to the bit stream
     *
     */
    void serialize(pixie::OutputBitStream& bs) const {
      bs << parent << left_child << right_child << middle;
      bit_vector_data.serialize(bs);
      data.serialize(bs);
    }

    static WaveletNode deserialize(std::span<const std::byte>& data)
      requires(std::same_as<Storage, MmapViewStorage>)
    {
      WaveletNode result;
      auto read = [&data](auto& value) {
        constexpr size_t length = sizeof(value);
        std::memcpy(&value, data.data(), length);
        data = data.subspan(length);
      };
      read(result.parent);
      read(result.left_child);
      read(result.right_child);
      read(result.middle);
      result.bit_vector_data = MmapViewStorage::deserialize(data);
      result.data = BitVectorBase<Storage>::deserialize(
          result.bit_vector_data.AsConst64BitInts(), data);
      return result;
    }
  };

  size_t alphabet_size_, data_size_;
  node_index_t root_;
  std::vector<WaveletNode> nodes_;
  std::vector<node_index_t> leaves_;
  std::vector<size_t> permutation_, inverse_permutation_;

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
    nodes[result].stream.reserve(prefix_sum[end] - prefix_sum[begin]);
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
    const auto& bit_vector = nodes_[node].bit_vector_data.AsConst64BitInts();
    for (size_t i = begin; i < end; i++) {
      if ((bit_vector[i / 64] >> (i % 64)) & 1) {
        dst[i - begin] = tmp[k++];
      } else {
        dst[i - begin] = tmp[j++];
      }
    }
  }

  WaveletTreeBase() = default;

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
  WaveletTreeBase(
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
        node.stream << go_right;
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
  size_t rank(uint64_t symbol, size_t pos) const {
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
  size_t select(uint64_t symbol, size_t rank) const {
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
   */
  std::vector<uint64_t> get_segment(size_t begin, size_t end) const {
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
  size_t size() { return data_size_; }

  /**
   * @brief Writes a wavelet tree serialization to the bit stream
   *
   */
  void serialize(pixie::OutputBitStream& bs) const {
    bs << alphabet_size_ << data_size_ << root_ << nodes_.size();
    for (const WaveletNode& node : nodes_) {
      node.serialize(bs);
    }
    for (const node_index_t leaf : leaves_) {
      bs << leaf;
    }
    for (const size_t idx : permutation_) {
      bs << idx;
    }
  }

  static WaveletTreeBase<MmapViewStorage> deserialize(
      std::span<const std::byte>& data) {
    WaveletTreeBase<MmapViewStorage> result;
    auto read = [&data](auto& value) {
      constexpr size_t length = sizeof(value);
      std::memcpy(&value, data.data(), length);
      data = data.subspan(length);
    };
    read(result.alphabet_size_);
    read(result.data_size_);
    read(result.root_);
    size_t size;
    read(size);
    result.nodes_.resize(size);
    for (auto& node: result.nodes_) {
      node = WaveletNode::deserialize(data);
    }
    result.leaves_.resize(result.alphabet_size_);
    for (node_index_t& leaf: result.leaves_) {
      read(leaf);
    }
    result.permutation_.resize(result.alphabet_size_);
    for (size_t& index: result.permutation_) {
      read(index);
    }
    result.inverse_permutation_.resize(result.alphabet_size_);
    for (size_t i = 0; i < result.alphabet_size_; i++) {
      result.inverse_permutation_[result.permutation_[i]] = i;
    }
    return result;
  }
};

typedef WaveletTreeBase<AlignedStorage> WaveletTree;

}  // namespace pixie
