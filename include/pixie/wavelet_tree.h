#pragma once

#include <pixie/bitvector.h>

#include <limits>

namespace pixie {

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

  node_index_t BuildNode(size_t begin,
                         size_t end,
                         std::span<uint64_t> data,
                         node_index_t parent) {
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

    size_t middle = begin + (end - begin) / 2;
    std::vector<uint64_t> bit_vector;
    bit_vector.resize((data.size() + 63) / 64);
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i] >= middle) {
        bit_vector[i / 64] |= 1ull << (i % 64);
      }
    }

    node_index_t result = nodes_.size();
    nodes_.emplace_back(middle, std::move(bit_vector), data.size());
    auto cut = std::stable_partition(
        data.begin(), data.end(), [middle](uint64_t x) { return x < middle; });
    nodes_[result].parent = parent;
    nodes_[result].left_child =
        BuildNode(begin, middle, std::span{data.begin(), cut}, result);
    nodes_[result].right_child =
        BuildNode(middle, end, std::span{cut, data.end()}, result);

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
                nodes_[node].middle - 1);
    } else {
      copySegmentContent(nodes_[node].left_child, rank0, rank0 + left,
                         tmp.subspan(0, left), dst.subspan(0, left));
    }
    if (nodes_[node].right_child == WaveletNode::kNil) {
      std::fill(tmp.begin() + static_cast<long long>(left), tmp.end(),
                nodes_[node].middle);
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
  WaveletTree(size_t alphabet_size, std::span<const uint64_t> data)
      : alphabet_size_(alphabet_size),
        data_size_(data.size()),
        leaves_(alphabet_size_, WaveletNode::kNil) {
    if (alphabet_size > 0) {
      std::vector<uint64_t> data_copy(data.begin(), data.end());
      nodes_.reserve(alphabet_size_);
      root_ = BuildNode(0, alphabet_size_, data_copy, WaveletNode::kNil);
    } else {
      root_ = WaveletNode::kNil;
    }
  }

  size_t rank(uint64_t symbol, size_t pos) const {
    if (symbol >= alphabet_size_) [[unlikely]] {
      return 0;
    }
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
