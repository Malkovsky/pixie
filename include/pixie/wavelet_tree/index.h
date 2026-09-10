#pragma once

#include <pixie/detail/byte_histogram.h>
#include <pixie/detail/huffman_build_table.h>
#include <pixie/detail/serialization.h>
#include <pixie/detail/wavelet_partition.h>
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
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pixie {

/**
 * @brief Storage-backed wavelet tree over an unsigned symbol type.
 * @tparam Symbol Unsigned symbol type. Its value range must cover the dense
 * alphabet `[0, alphabet_size)`.
 * @tparam Storage Owning aligned storage or a non-owning read-only view.
 */
template <WaveletTreeSymbol Symbol,
          StorageImplementation Storage = AlignedStorage>
class WaveletTreeIndex
    : public WaveletTreeBase<WaveletTreeIndex<Symbol, Storage>, Symbol>,
      public SerializationBase<WaveletTreeIndex<Symbol, Storage>> {
 private:
  using node_index_t = size_t;
  static constexpr node_index_t npos = std::numeric_limits<node_index_t>::max();
  static constexpr std::array<std::uint8_t, 8> kSerializationMagic = {
      'P', 'X', 'W', 'A', 'V', 'E', 'T', '\0'};
  static constexpr std::uint32_t kSerializationVersion = 5;
  static constexpr std::size_t kSerializationHeaderBytes = 24;
#if defined(__APPLE__) && defined(__aarch64__)
  static constexpr std::size_t kByteConstructionBlockSize = 16 * 1024;
#else
  static constexpr std::size_t kByteConstructionBlockSize = 32 * 1024;
#endif

  struct PreWaveletNode {
    node_index_t parent = npos;
    node_index_t left_child = npos;
    node_index_t right_child = npos;
    std::size_t middle;
    std::size_t left_size = 0;
    PackedBitBuilder stream;
    explicit PreWaveletNode(std::size_t middle) : middle(middle) {}
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
    std::size_t middle;
    Storage bit_vector_data;
    RankSelectSupport<Storage> data;

    WaveletNode() = default;

    WaveletNode(const WaveletNode& node)
        : parent(node.parent),
          left_child(node.left_child),
          right_child(node.right_child),
          middle(node.middle),
          bit_vector_data(node.bit_vector_data),
          data([&] {
            if constexpr (std::same_as<Storage, AlignedStorage>) {
              return RankSelectSupport<Storage>(bit_vector_data.as_words64(),
                                                node.data.size());
            } else {
              return node.data;
            }
          }()) {}

    WaveletNode& operator=(const WaveletNode& node) {
      if (this != &node) {
        WaveletNode copy(node);
        *this = std::move(copy);
      }
      return *this;
    }

    WaveletNode(WaveletNode&&) noexcept = default;
    WaveletNode& operator=(WaveletNode&&) noexcept = default;

    WaveletNode(PreWaveletNode&& node)
      requires(std::same_as<Storage, AlignedStorage>)
        : parent(node.parent),
          left_child(node.left_child),
          right_child(node.right_child),
          middle(node.middle) {
      const std::size_t bit_count = node.stream.size_bits();
      const std::vector<std::uint64_t> words = node.stream.take_words();
      bit_vector_data = AlignedStorage(std::span<const std::uint64_t>(words));
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

    /** @brief Construct one checked node from @p reader. */
    static WaveletNode deserialize(BinaryReader& reader,
                                   DeserializationValidation validation) {
      WaveletNode result;
      result.parent = reader.read_size();
      result.left_child = reader.read_size();
      result.right_child = reader.read_size();
      result.middle = reader.read_u64();
      result.bit_vector_data = Storage::deserialize(reader);
      result.data = RankSelectSupport<Storage>::deserialize(
          reader, result.bit_vector_data.as_words64(), validation);
      return result;
    }
  };

  size_t alphabet_size_ = 0;
  size_t data_size_ = 0;
  node_index_t root_ = npos;
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
          bool has_mapped_symbol = false;
          for (std::size_t symbol = symbol_begin; symbol < symbol_end;
               ++symbol) {
            if (leaves_[symbol] == node) {
              has_mapped_symbol = true;
            } else if (leaves_[symbol] != npos) {
              throw std::invalid_argument(
                  "Serialized wavelet-tree leaf map disagrees with topology");
            }
          }
          if (expected_size != 0 && !has_mapped_symbol) {
            throw std::invalid_argument(
                "Serialized wavelet-tree leaf is missing from its map");
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
    nodes[result].left_size = prefix_sum[middle] - prefix_sum[begin];
    nodes[result].stream.reserve_bits(prefix_sum[end] - prefix_sum[begin]);
    nodes[result].parent = parent;
    nodes[result].left_child =
        build_node(begin, middle, result, get_middle, prefix_sum, nodes);
    nodes[result].right_child =
        build_node(middle, end, result, get_middle, prefix_sum, nodes);

    return result;
  }

  static std::size_t count_huffman_subtree(
      const detail::HuffmanBuildTable& table,
      std::int16_t source_node,
      std::span<const std::size_t> symbol_counts,
      std::array<std::size_t, detail::kMaximumHuffmanNodes>& subtree_counts) {
    const detail::HuffmanBuildNode& node = table.tree[source_node];
    if (node.symbol >= 0) {
      return subtree_counts[source_node] = symbol_counts[node.symbol];
    }
    const std::size_t left =
        count_huffman_subtree(table, node.left, symbol_counts, subtree_counts);
    const std::size_t right =
        count_huffman_subtree(table, node.right, symbol_counts, subtree_counts);
    return subtree_counts[source_node] = left + right;
  }

  node_index_t materialize_huffman_node(
      const detail::HuffmanBuildTable& table,
      std::int16_t source_node,
      node_index_t parent,
      const std::array<std::size_t, detail::kMaximumHuffmanNodes>&
          subtree_counts,
      std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Symbol, std::uint8_t> &&
             std::same_as<Storage, AlignedStorage>)
  {
    const detail::HuffmanBuildNode& source = table.tree[source_node];
    if (source.symbol >= 0) {
      leaves_[table.symbol_to_rank[source.symbol]] = parent;
      return npos;
    }

    const node_index_t result = nodes.size();
    const std::size_t middle =
        static_cast<std::size_t>(table.split_rank[source_node]) + 1;
    nodes.emplace_back(middle);
    nodes[result].parent = parent;
    nodes[result].left_size = subtree_counts[source.left];
    nodes[result].stream.reserve_bits(subtree_counts[source_node]);
    nodes[result].left_child = materialize_huffman_node(
        table, source.left, result, subtree_counts, nodes);
    nodes[result].right_child = materialize_huffman_node(
        table, source.right, result, subtree_counts, nodes);
    return result;
  }

  void build_huffman_byte_topology(std::span<const std::size_t> symbol_counts,
                                   std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Symbol, std::uint8_t> &&
             std::same_as<Storage, AlignedStorage>)
  {
    const detail::HuffmanBuildTable table =
        detail::build_huffman_table(symbol_counts);
    permutation_.resize(alphabet_size_);
    inverse_permutation_.resize(alphabet_size_);

    std::size_t unused_rank = table.symbol_count;
    for (std::size_t symbol = 0; symbol < alphabet_size_; ++symbol) {
      const std::size_t rank = symbol_counts[symbol] == 0
                                   ? unused_rank++
                                   : table.symbol_to_rank[symbol];
      permutation_[symbol] = rank;
      inverse_permutation_[rank] = symbol;
    }
    if (table.symbol_count < 2) {
      root_ = npos;
      return;
    }

    std::array<std::size_t, detail::kMaximumHuffmanNodes> subtree_counts{};
    count_huffman_subtree(table, table.root, symbol_counts, subtree_counts);
    nodes.reserve(table.symbol_count - 1);
    root_ = materialize_huffman_node(table, table.root, npos, subtree_counts,
                                     nodes);
  }

  void build_node_streams(node_index_t node,
                          std::span<Symbol> input,
                          std::span<Symbol> output,
                          std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    PreWaveletNode& current = nodes[node];
    detail::partition_wavelet_ranks(
        std::span<const Symbol>(input), static_cast<Symbol>(current.middle),
        output, current.left_size, current.left_child != npos,
        current.right_child != npos, current.stream);

    if (current.left_child != npos) {
      build_node_streams(current.left_child, output.first(current.left_size),
                         input.first(current.left_size), nodes);
    }
    if (current.right_child != npos) {
      build_node_streams(current.right_child, output.subspan(current.left_size),
                         input.subspan(current.left_size), nodes);
    }
  }

  // PivCo-style block walk: keep the left subsequence in place, put the right
  // subsequence in scratch, and omit scatter for leaf children. Separate node
  // builders make every block append bit-contiguously without wire padding.
  void build_byte_block_streams(node_index_t node,
                                std::span<std::uint8_t> ranks,
                                std::span<std::uint8_t> scratch,
                                std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Symbol, std::uint8_t> &&
             std::same_as<Storage, AlignedStorage>)
  {
    PreWaveletNode& current = nodes[node];
    const bool write_left = current.left_child != npos;
    const bool write_right = current.right_child != npos;
    const std::size_t left_size = detail::partition_wavelet_byte_block(
        ranks, static_cast<std::uint8_t>(current.middle), scratch, write_left,
        write_right, current.stream);
    const std::size_t right_size = ranks.size() - left_size;

    if (write_left && left_size != 0) {
      build_byte_block_streams(current.left_child, ranks.first(left_size),
                               scratch.subspan(right_size, left_size), nodes);
    }
    if (write_right && right_size != 0) {
      build_byte_block_streams(current.right_child, scratch.first(right_size),
                               ranks.first(right_size), nodes);
    }
  }

  template <class ForEachSymbol>
  void build_byte_bit_streams(std::span<const std::size_t> symbol_counts,
                              ForEachSymbol& for_each_symbol,
                              std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Symbol, std::uint8_t> &&
             std::same_as<Storage, AlignedStorage>)
  {
    std::array<std::uint8_t, detail::kByteAlphabetSize> symbol_to_rank{};
    std::array<std::uint16_t, detail::kByteAlphabetSize> symbol_to_high_rank{};
    for (std::size_t symbol = 0; symbol < alphabet_size_; ++symbol) {
      const auto rank = static_cast<std::uint8_t>(permutation_[symbol]);
      symbol_to_rank[symbol] = rank;
      symbol_to_high_rank[symbol] = static_cast<std::uint16_t>(rank) << 8;
    }

    const std::size_t block_capacity = std::min(
        kByteConstructionBlockSize, std::max<std::size_t>(data_size_, 1));
    std::vector<std::uint8_t> block;
    std::vector<std::uint8_t> scratch(block_capacity);
    block.reserve(block_capacity);
    detail::ByteHistogram actual_histogram;
    const auto encode_block = [&] {
      if (block.empty()) {
        return;
      }
      if (root_ != npos) {
        build_byte_block_streams(root_, block,
                                 std::span(scratch).first(block.size()), nodes);
      }
      block.clear();
    };

    std::size_t emitted_size = 0;
    const auto append_symbols = [&](std::span<const Symbol> symbols) {
      if (emitted_size > data_size_ ||
          symbols.size() > data_size_ - emitted_size) {
        throw std::invalid_argument(
            "Wavelet-tree emitted symbols do not match their counts");
      }
      emitted_size += symbols.size();
      actual_histogram.add(std::as_bytes(symbols));
      while (!symbols.empty()) {
        const std::size_t copied =
            std::min(block_capacity - block.size(), symbols.size());
        const std::size_t destination = block.size();
        block.resize(destination + copied);
        detail::map_wavelet_byte_ranks(
            symbols.first(copied),
            std::span(block).subspan(destination, copied), symbol_to_rank,
            symbol_to_high_rank);
        symbols = symbols.subspan(copied);
        if (block.size() == block_capacity) {
          encode_block();
        }
      }
    };
    for_each_symbol([&](auto&& emitted) {
      using Emitted = std::remove_cvref_t<decltype(emitted)>;
      if constexpr (std::same_as<Emitted, Symbol>) {
        append_symbols(std::span<const Symbol>(&emitted, 1));
      } else {
        append_symbols(std::span<const Symbol>(emitted));
      }
    });
    encode_block();

    const auto actual_counts = actual_histogram.counts();
    if (emitted_size != data_size_ ||
        !std::ranges::equal(actual_counts.begin(),
                            actual_counts.begin() + alphabet_size_,
                            symbol_counts.begin(), symbol_counts.end()) ||
        std::ranges::any_of(actual_counts.begin() + alphabet_size_,
                            actual_counts.end(),
                            [](std::size_t count) { return count != 0; })) {
      throw std::invalid_argument(
          "Wavelet-tree emitted symbols do not match their counts");
    }
  }

  template <class ForEachSymbol>
  void build_generic_bit_streams(std::span<const std::size_t> symbol_counts,
                                 ForEachSymbol& for_each_symbol,
                                 std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    std::vector<std::size_t> actual_counts(alphabet_size_);
    std::vector<Symbol> ranks;
    if (root_ != npos) {
      ranks.reserve(data_size_);
    }
    for_each_symbol([&](Symbol symbol) {
      const std::size_t original = checked_symbol_index(symbol, alphabet_size_);
      if (actual_counts[original] == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("Wavelet-tree symbol count is too large");
      }
      ++actual_counts[original];
      if (root_ != npos) {
        ranks.push_back(static_cast<Symbol>(permutation_[original]));
      }
    });
    if (!std::ranges::equal(actual_counts, symbol_counts)) {
      throw std::invalid_argument(
          "Wavelet-tree emitted symbols do not match their counts");
    }
    if (root_ == npos) {
      return;
    }

    const PreWaveletNode& root = nodes[root_];
    const bool root_needs_output =
        root.left_child != npos || root.right_child != npos;
    std::vector<Symbol> scratch(root_needs_output ? data_size_ : 0);
    build_node_streams(root_, ranks, scratch, nodes);
  }

  template <class ForEachSymbol>
  void build_bit_streams(std::span<const std::size_t> symbol_counts,
                         ForEachSymbol& for_each_symbol,
                         std::vector<PreWaveletNode>& nodes)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    if constexpr (std::same_as<Symbol, std::uint8_t>) {
      build_byte_bit_streams(symbol_counts, for_each_symbol, nodes);
    } else {
      build_generic_bit_streams(symbol_counts, for_each_symbol, nodes);
    }
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
                            std::span<Symbol> dst,
                            std::span<Symbol> tmp) const {
    if (begin == end) {
      return;
    }
    const size_t rank = nodes_[node].data.rank(begin), rank0 = begin - rank;
    const size_t right = nodes_[node].data.rank(end) - rank,
                 left = (end - begin) - right;

    if (nodes_[node].left_child == npos) {
      std::fill_n(
          tmp.begin(), static_cast<long long>(left),
          static_cast<Symbol>(inverse_permutation_[nodes_[node].middle - 1]));
    } else {
      copy_segment_content(nodes_[node].left_child, rank0, rank0 + left,
                           tmp.subspan(0, left), dst.subspan(0, left));
    }
    if (nodes_[node].right_child == npos) {
      std::fill(tmp.begin() + static_cast<long long>(left), tmp.end(),
                static_cast<Symbol>(inverse_permutation_[nodes_[node].middle]));
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

  static void validate_alphabet_size(std::size_t alphabet_size) {
    if (alphabet_size != 0 &&
        alphabet_size - 1 >
            static_cast<std::size_t>(std::numeric_limits<Symbol>::max())) {
      throw std::invalid_argument(
          "Wavelet-tree alphabet does not fit its symbol type");
    }
  }

  static std::size_t checked_symbol_index(Symbol symbol,
                                          std::size_t alphabet_size) {
    const std::size_t index = static_cast<std::size_t>(symbol);
    if (index >= alphabet_size) {
      throw std::invalid_argument(
          "Wavelet-tree symbol is outside the alphabet");
    }
    return index;
  }

  template <class ForEachSymbol>
  void build_from_counts(std::size_t alphabet_size,
                         std::span<const std::size_t> symbol_counts,
                         ForEachSymbol&& for_each_symbol,
                         WaveletTreeBuildType build_type)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    validate_alphabet_size(alphabet_size);
    if (symbol_counts.size() != alphabet_size) {
      throw std::invalid_argument(
          "Wavelet-tree symbol counts must match the alphabet size");
    }
    alphabet_size_ = alphabet_size;
    for (const std::size_t count : symbol_counts) {
      if (count > std::numeric_limits<std::size_t>::max() - data_size_) {
        throw std::length_error("Wavelet-tree input is too large");
      }
      data_size_ += count;
    }
    leaves_.assign(alphabet_size_, npos);

    std::vector<PreWaveletNode> nodes;
    if constexpr (std::same_as<Symbol, std::uint8_t>) {
      if (build_type == WaveletTreeBuildType::Huffman) {
        if (alphabet_size_ != 0) {
          build_huffman_byte_topology(symbol_counts, nodes);
        }
        build_bit_streams(symbol_counts, for_each_symbol, nodes);
        nodes_.reserve(nodes.size());
        for (auto& node : nodes) {
          nodes_.emplace_back(std::move(node));
        }
        return;
      }
    }

    std::vector<std::size_t> nodes_structure;
    if (alphabet_size_ != 0) {
      nodes.reserve(alphabet_size_);
      nodes_structure.reserve(alphabet_size_);

      if (build_type == WaveletTreeBuildType::Standard) {
        permutation_.resize(alphabet_size_);
        inverse_permutation_.resize(alphabet_size_);
        std::iota(permutation_.begin(), permutation_.end(), 0);
        std::iota(inverse_permutation_.begin(), inverse_permutation_.end(), 0);
        nodes_structure.resize(alphabet_size_, npos);
      } else {
        struct HuffmanNode {
          std::size_t size;
          std::size_t left;
          std::size_t right;
        };
        std::vector<HuffmanNode> huffman_nodes(alphabet_size_, {0, 0, 0});
        for (std::size_t symbol = 0; symbol < alphabet_size_; ++symbol) {
          huffman_nodes[symbol].size = symbol_counts[symbol];
        }

        using QueueElement = std::pair<std::size_t, std::size_t>;
        std::priority_queue<QueueElement, std::vector<QueueElement>,
                            std::greater<>>
            queue;
        for (std::size_t symbol = 0; symbol < alphabet_size_; ++symbol) {
          queue.emplace(huffman_nodes[symbol].size, symbol);
        }
        while (queue.size() >= 2) {
          const std::size_t right = queue.top().second;
          queue.pop();
          const std::size_t left = queue.top().second;
          queue.pop();
          huffman_nodes.push_back(
              {huffman_nodes[left].size + huffman_nodes[right].size, left + 1,
               right + 1});
          queue.emplace(huffman_nodes.back().size, huffman_nodes.size() - 1);
        }

        std::function<std::size_t(std::size_t)> enumerate =
            [&](std::size_t index) -> std::size_t {
          const auto& [size, left, right] = huffman_nodes[index];
          if (left == 0 || right == 0) {
            permutation_[index] = inverse_permutation_.size();
            inverse_permutation_.push_back(index);
            return 1;
          }
          const std::size_t node = nodes_structure.size();
          std::size_t subtree = 0;
          if (size > 0) {
            nodes_structure.push_back(0);
          }
          subtree += enumerate(left - 1);
          if (size > 0) {
            nodes_structure[node] = subtree;
          }
          subtree += enumerate(right - 1);
          return subtree;
        };

        permutation_.resize(alphabet_size_);
        inverse_permutation_.reserve(alphabet_size_);
        enumerate(huffman_nodes.size() - 1);
      }

      std::vector<std::size_t> prefix_sum(alphabet_size_ + 1);
      for (std::size_t symbol = 0; symbol < alphabet_size_; ++symbol) {
        prefix_sum[permutation_[symbol] + 1] = symbol_counts[symbol];
      }
      std::partial_sum(prefix_sum.begin(), prefix_sum.end(),
                       prefix_sum.begin());
      root_ = build_node(
          0, alphabet_size_, npos,
          [&](node_index_t node) { return nodes_structure[node]; }, prefix_sum,
          nodes);
    }

    build_bit_streams(symbol_counts, for_each_symbol, nodes);

    nodes_.reserve(nodes.size());
    for (auto& node : nodes) {
      nodes_.emplace_back(std::move(node));
    }
  }

  WaveletTreeIndex() = default;

 public:
  using symbol_type = Symbol;

  /**
   * @brief Construct from a contiguous sequence of typed symbols.
   * @details Construction temporarily owns up to two buffers of one `Symbol`
   * per input symbol. The buffers are released before the permanent node
   * indexes are materialized.
   * @param alphabet_size Dense alphabet size; every symbol must be smaller.
   * @param data Input symbols retained only for the duration of construction.
   * @param build_type Standard or Huffman-shaped construction.
   * @throws std::invalid_argument if the alphabet or a symbol is invalid.
   */
  WaveletTreeIndex(
      std::size_t alphabet_size,
      std::span<const Symbol> data,
      const WaveletTreeBuildType build_type = WaveletTreeBuildType::Standard)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    validate_alphabet_size(alphabet_size);
    std::vector<std::size_t> counts;
    if constexpr (std::same_as<Symbol, std::uint8_t>) {
      detail::ByteHistogram histogram;
      histogram.add(std::as_bytes(data));
      const auto byte_counts = histogram.counts();
      if (std::ranges::any_of(byte_counts.begin() + alphabet_size,
                              byte_counts.end(),
                              [](std::size_t count) { return count != 0; })) {
        throw std::invalid_argument(
            "Wavelet-tree symbol is outside the alphabet");
      }
      counts.assign(byte_counts.begin(), byte_counts.begin() + alphabet_size);
    } else {
      counts.assign(alphabet_size, 0);
      for (const Symbol symbol : data) {
        ++counts[checked_symbol_index(symbol, alphabet_size)];
      }
    }
    build_from_counts(
        alphabet_size, counts,
        [&](auto&& emit) {
          if constexpr (requires { emit(data); }) {
            emit(data);
          } else {
            for (const Symbol symbol : data) {
              emit(symbol);
            }
          }
        },
        build_type);
  }

  /**
   * @brief Construct from counts and one streamed pass over the symbols.
   * @details @p for_each_symbol is invoked exactly once with a consumer that
   * accepts one `Symbol`. Emitted symbols must exactly match @p symbol_counts;
   * this permits callers to scan a replayable source once for counts and once
   * for construction without materializing the sequence themselves.
   * Construction temporarily owns up to two buffers of one `Symbol` per
   * emitted symbol. The buffers are released before the permanent node indexes
   * are materialized.
   * @param alphabet_size Dense alphabet size.
   * @param symbol_counts Count for every symbol in alphabet order.
   * @param for_each_symbol Callable accepting the construction consumer.
   * @param build_type Standard or Huffman-shaped construction.
   * @throws std::invalid_argument for invalid or inconsistent input.
   */
  template <class ForEachSymbol>
  WaveletTreeIndex(
      std::size_t alphabet_size,
      std::span<const std::size_t> symbol_counts,
      ForEachSymbol&& for_each_symbol,
      const WaveletTreeBuildType build_type = WaveletTreeBuildType::Standard)
    requires(std::same_as<Storage, AlignedStorage>)
  {
    build_from_counts(alphabet_size, symbol_counts,
                      std::forward<ForEachSymbol>(for_each_symbol), build_type);
  }

  /**
   * @brief Rank of specified symbol up to position pos (exclusive)
   *
   * @param symbol The character that the query is about
   * @param pos Character index in [0, size()]
   * @return Number of specified symbols in [0, pos)
   *
   */
  size_t rank_impl(Symbol symbol, size_t pos) const {
    std::size_t symbol_index = static_cast<std::size_t>(symbol);
    if (symbol_index >= alphabet_size_) [[unlikely]] {
      return 0;
    }
    symbol_index = permutation_[symbol_index];
    if (root_ == npos) [[unlikely]] {
      return data_size_ != 0 && symbol_index == 0 ? pos : 0;
    }
    if (leaves_[symbol_index] == npos) [[unlikely]] {
      return 0;
    }
    for (node_index_t current = root_; current != npos;) {
      const WaveletNode& node = nodes_[current];
      if (symbol_index < node.middle) {
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
  size_t select_impl(Symbol symbol, size_t rank) const {
    std::size_t symbol_index = static_cast<std::size_t>(symbol);
    if (symbol_index >= alphabet_size_ || data_size_ == 0 || rank == 0)
        [[unlikely]] {
      return data_size_;
    }
    symbol_index = permutation_[symbol_index];
    if (root_ == npos) [[unlikely]] {
      return symbol_index == 0 && rank <= data_size_ ? rank - 1 : data_size_;
    }
    if (leaves_[symbol_index] == npos) [[unlikely]] {
      return data_size_;
    }
    node_index_t current = leaves_[symbol_index];
    for (; current != npos; current = nodes_[current].parent) {
      const WaveletNode& node = nodes_[current];
      if (symbol_index < node.middle) {
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
   * result storage of two `Symbol` values per returned symbol.
   *
   */
  std::vector<Symbol> get_segment_impl(size_t begin, size_t end) const {
    if (alphabet_size_ == 0 || data_size_ == 0 || begin >= end) [[unlikely]] {
      return {};
    }
    const std::size_t length = end - begin;
    if (root_ == npos) [[unlikely]] {
      return std::vector<Symbol>(
          length, static_cast<Symbol>(inverse_permutation_.front()));
    }
    if (length > std::vector<Symbol>().max_size() / 2) {
      throw std::length_error("Wavelet-tree segment is too large");
    }
    std::vector<Symbol> result(2 * length);
    copy_segment_content(root_, begin, end, std::span(result).first(length),
                         std::span(result).subspan(length));
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
  void serialize_impl(BinaryWriter& writer) const {
    if (writer.size_bytes() % alignof(std::uint64_t) != 0) {
      throw std::invalid_argument(
          "Wavelet-tree serialization requires an aligned writer offset");
    }
    const std::size_t artifact_begin = writer.size_bytes();
    detail::write_magic(writer, kSerializationMagic);
    writer.write_u32(kSerializationVersion);
    writer.write_u32(std::numeric_limits<Symbol>::digits);
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
   * @brief Restore one checked wavelet-tree artifact.
   *
   * @details The aligned-storage specialization copies all restored metadata.
   * The read-only specialization retains views into the reader's backing
   * bytes, which must remain alive, immutable, and aligned for 64-bit access.
   * On success @p reader advances past exactly one framed artifact; on failure
   * it is unchanged. @p validation selects quick structural checks or exact
   * bitvector-derived metadata validation.
   *
   * @param reader Input cursor, advanced only after successful validation.
   * @param validation Quick structural or full bitvector-derived validation.
   *
   * @throws std::invalid_argument for malformed, truncated, incompatible, or
   * structurally inconsistent metadata.
   * @throws std::length_error when an encoded count is not representable.
   */
  static WaveletTreeIndex deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires(std::same_as<Storage, AlignedStorage> ||
             std::same_as<Storage, ReadOnlyStorageView>)
  {
    BinaryReader candidate = reader;
    if constexpr (std::same_as<Storage, ReadOnlyStorageView>) {
      if (reinterpret_cast<std::uintptr_t>(candidate.remaining_bytes().data()) %
              alignof(std::uint64_t) !=
          0) {
        throw std::invalid_argument(
            "Serialized wavelet-tree artifact is not word aligned");
      }
    }
    const std::size_t available_size = candidate.remaining();
    detail::require_magic(candidate, kSerializationMagic);
    if (candidate.read_u32() != kSerializationVersion ||
        candidate.read_u32() != std::numeric_limits<Symbol>::digits) {
      throw std::invalid_argument(
          "Incompatible serialized wavelet-tree artifact");
    }
    const std::size_t artifact_size = detail::checked_artifact_size(
        candidate.read_u64(), kSerializationHeaderBytes, available_size);
    BinaryReader payload =
        candidate.read_subreader(artifact_size - kSerializationHeaderBytes);

    WaveletTreeIndex result;
    result.alphabet_size_ = payload.read_size();
    result.validate_alphabet_size(result.alphabet_size_);
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
};

template <WaveletTreeSymbol Symbol>
using WaveletTree = WaveletTreeIndex<Symbol, AlignedStorage>;

template <WaveletTreeSymbol Symbol>
using WaveletTreeView = WaveletTreeIndex<Symbol, ReadOnlyStorageView>;

}  // namespace pixie
