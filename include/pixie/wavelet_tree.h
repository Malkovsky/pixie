#pragma once

#include <pixie/bitvector.h>

#include <limits>

namespace pixie {

enum WaveletTreeBuildType { Standard, Huffman };

class WaveletTree {
  using node_index_t = size_t;
  struct WaveletNode {
    static constexpr node_index_t kNil =
        std::numeric_limits<node_index_t>::max();
    node_index_t parent = kNil, left_child = kNil, right_child = kNil;
    uint64_t middle;
    std::vector<uint64_t> bit_vector_;
    BitVector data;
    WaveletNode(uint64_t middle,
                std::vector<uint64_t>&& bit_vector,
                size_t num_bits)
        : middle(middle),
          bit_vector_(std::move(bit_vector)),
          data(std::span{bit_vector_}, num_bits) {}
  };

  size_t alphabet_size_, data_size_;
  node_index_t root_;
  std::vector<WaveletNode> nodes_;
  std::vector<node_index_t> leaves_;
  std::vector<size_t> permutation_, inverse_permutation_;

  node_index_t BuildNode(
      size_t begin,
      size_t end,
      std::span<uint64_t> data,
      node_index_t parent,
      const std::function<size_t(node_index_t)>& get_middle = [](auto) {
        return -1ull;
      }) {
    if (end - begin == 1) {
      leaves_[begin] = parent;
      return WaveletNode::kNil;
    }
    if (data.empty()) {
      for (size_t symb = begin; symb < end; symb++) {
        leaves_[symb] = parent;
      }
      return WaveletNode::kNil;
    }

    node_index_t result = nodes_.size();
    size_t middle = get_middle(result);
    middle = begin + (middle == -1ull ? (end - begin) / 2 : middle);
    std::vector<uint64_t> bit_vector;
    bit_vector.resize((data.size() + 63) / 64);
    for (size_t i = 0; i < data.size(); i++) {
      if (permutation_[data[i]] >= middle) {
        bit_vector[i / 64] |= 1ull << (i % 64);
      }
    }

    nodes_.emplace_back(middle, std::move(bit_vector), data.size());
    auto cut = std::stable_partition(
        data.begin(), data.end(),
        [middle, this](uint64_t x) { return permutation_[x] < middle; });
    nodes_[result].parent = parent;
    nodes_[result].left_child = BuildNode(
        begin, middle, std::span{data.begin(), cut}, result, get_middle);
    nodes_[result].right_child =
        BuildNode(middle, end, std::span{cut, data.end()}, result, get_middle);

    return result;
  }

  void copySegmentContent(node_index_t node,
                          size_t begin,
                          size_t end,
                          std::span<uint64_t> dst,
                          std::span<uint64_t> tmp) const {
    if (begin == end) {
      return;
    }
    size_t rank = nodes_[node].data.rank(begin), rank0 = begin - rank;
    size_t right = nodes_[node].data.rank(end) - rank,
           left = (end - begin) - right;

    if (nodes_[node].left_child == WaveletNode::kNil) {
      std::fill(tmp.begin(), tmp.begin() + static_cast<long long>(left),
                inverse_permutation_[nodes_[node].middle - 1]);
    } else {
      copySegmentContent(nodes_[node].left_child, rank0, rank0 + left,
                         tmp.subspan(0, left), dst.subspan(0, left));
    }
    if (nodes_[node].right_child == WaveletNode::kNil) {
      std::fill(tmp.begin() + static_cast<long long>(left), tmp.end(),
                inverse_permutation_[nodes_[node].middle]);
    } else {
      copySegmentContent(nodes_[node].right_child, rank, rank + right,
                         tmp.subspan(left, right), dst.subspan(left, right));
    }

    size_t j = 0, k = left;
    const auto& bit_vector = nodes_[node].bit_vector_;
    for (size_t i = begin; i < end; i++) {
      if ((bit_vector[i / 64] >> (i % 64)) & 1) {
        dst[i - begin] = tmp[k++];
      } else {
        dst[i - begin] = tmp[j++];
      }
    }
  }

 public:
  WaveletTree(size_t alphabet_size,
              std::span<const uint64_t> data,
              WaveletTreeBuildType build_type = WaveletTreeBuildType::Standard)
      : alphabet_size_(alphabet_size),
        data_size_(data.size()),
        leaves_(alphabet_size_, WaveletNode::kNil) {
    if (alphabet_size == 0) {
      root_ = WaveletNode::kNil;
      return;
    }
    nodes_.reserve(alphabet_size_);
    if (build_type == WaveletTreeBuildType::Standard) {
      permutation_.resize(alphabet_size);
      inverse_permutation_.resize(alphabet_size);
      std::iota(permutation_.begin(), permutation_.end(), 0);
      std::iota(inverse_permutation_.begin(), inverse_permutation_.end(), 0);

      std::vector<uint64_t> data_copy(data.begin(), data.end());
      root_ = BuildNode(0, alphabet_size_, data_copy, WaveletNode::kNil);
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

      std::vector<size_t> nodes_structure;
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
      nodes_structure.reserve(alphabet_size_);
      inverse_permutation_.reserve(alphabet_size_);
      enumerate(huffman_nodes.size() - 1);

      std::vector<uint64_t> data_copy(data.begin(), data.end());
      root_ =
          BuildNode(0, alphabet_size_, data_copy, WaveletNode::kNil,
                    [&](node_index_t node) { return nodes_structure[node]; });
    }
  }

  size_t rank(uint64_t symbol, size_t pos) const {
    if (symbol >= alphabet_size_) [[unlikely]] {
      return 0;
    }
    symbol = permutation_[symbol];
    for (node_index_t current = root_; current != WaveletNode::kNil;) {
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

  size_t select(uint64_t symbol, size_t rank) const {
    if (symbol >= alphabet_size_ || data_size_ == 0) [[unlikely]] {
      return data_size_;
    }
    symbol = permutation_[symbol];
    node_index_t current = leaves_[symbol];
    for (; current != WaveletNode::kNil; current = nodes_[current].parent) {
      const WaveletNode& node = nodes_[current];
      if (symbol < node.middle) {
        rank = node.data.select0(rank) + 1;
      } else {
        rank = node.data.select(rank) + 1;
      }
    }
    return rank - 1;
  }

  std::vector<uint64_t> getSegment(size_t begin, size_t end) const {
    if (alphabet_size_ == 0 || data_size_ == 0) [[unlikely]] {
      return {};
    }
    auto length = static_cast<long long>(end - begin);
    std::vector<uint64_t> result(2 * length);
    copySegmentContent(root_, begin, end,
                       std::span{result.begin(), result.begin() + length},
                       std::span{result.begin() + length, result.end()});
    result.resize(length);
    return result;
  }
};

}  // namespace pixie
