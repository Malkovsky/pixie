#pragma once

/**
 * @file pivco_huffman.h
 * @brief Simple scalar PivCo-Huffman codec.
 *
 * Stores a Huffman-shaped tree of per-node routing bitmaps (the "tree of
 * bitmaps" layout shared with wavelet trees). Encoding walks each input symbol
 * from the root, appending one direction bit to every internal node on its
 * path. Decoding merges child symbol streams bottom-up from leaves to root,
 * using the node bitmap as a selector.
 *
 * This is a deliberately simple, unoptimized reference implementation: node
 * bitmaps are stored as packed `std::uint64_t` words, traversal is scalar,
 * and there are no flat-subtree, non-canonical, SIMD, or selective-ANS
 * optimizations. It is intended as a correct baseline for the PivCo-Huffman
 * layout and for future optimized variants.
 *
 * Range and ownership conventions follow `HuffmanBase`: symbols are bytes,
 * the compressed stream is a byte view, and the codec owns its serialized form.
 */

#include <pixie/huffman.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <queue>
#include <span>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Simple scalar PivCo-Huffman codec.
 * @details Implements `HuffmanBase` by building a Huffman tree over the
 *          byte alphabet, storing one routing bitmap per internal node as
 *          packed 64-bit words, and decoding by bottom-up merging of child
 *          symbol streams. The serialized form is a header followed by the
 *          tree structure and packed node bitmaps.
 */
class PivCoHuffman : public HuffmanBase<PivCoHuffman> {
 public:
  /** @brief Symbol type handled by the codec: one byte per symbol. */
  using symbol_type = std::uint8_t;

  /**
   * @brief Build a codec by encoding @p input.
   * @param input Symbol sequence to compress.
   */
  explicit PivCoHuffman(std::span<const symbol_type> input) {
    build(input);
    serialize();
  }

  /**
   * @brief Load a codec from a previously serialized compressed stream.
   * @param compressed Compressed byte view produced by another instance.
   */
  explicit PivCoHuffman(std::span<const std::byte> compressed) {
    deserialize(compressed);
  }

  // --- CRTP extension points ----------------------------------------------

  /** @brief Size of the original symbol stream. */
  std::size_t uncompressed_size_impl() const { return uncompressed_size_; }

  /** @brief Size of the compressed byte stream. */
  std::size_t compressed_size_impl() const { return compressed_.size(); }

  /** @brief Read-only view of the compressed byte stream. */
  std::span<const std::byte> compressed_data_impl() const {
    return compressed_;
  }

  /** @brief Reconstruct the original symbol sequence. */
  std::vector<symbol_type> decode_impl() const { return decode_from_tree(); }

 private:
  /// @brief Sentinel node index meaning "no node".
  static constexpr std::size_t kNpos = static_cast<std::size_t>(-1);

  /// @brief Byte-alphabet size (256 symbols).
  static constexpr std::size_t kAlphabet = 256;

  /// @brief Bits per word used by the packed bitmap.
  static constexpr std::size_t kWordBits = 64;

  /**
   * @brief Packed routing bitmap stored as 64-bit words.
   * @details Bit `i` lives in word `i / 64` at position `i % 64`. This avoids
   *          the slow bit-proxy access of `std::vector<bool>` and allows
   *          bulk serialization via direct word `memcpy`.
   */
  struct Bitmap {
    std::vector<std::uint64_t> words;
    std::size_t count = 0;

    /** @brief Number of 64-bit words backing the bitmap. */
    std::size_t word_count() const {
      return (count + kWordBits - 1) / kWordBits;
    }

    /** @brief Number of set (1) bits in the bitmap. */
    std::size_t popcount() const {
      std::size_t total = 0;
      for (std::uint64_t w : words) {
        total += std::popcount(w);
      }
      return total;
    }
  };

  /**
   * @brief A node in the Huffman tree of bitmaps.
   * @details Internal nodes carry a routing bitmap with one bit per symbol
   *          that passes through the node (0 = left, 1 = right), in input
   *          order. Leaves carry only their symbol.
   */
  struct Node {
    std::size_t left = kNpos;
    std::size_t right = kNpos;
    std::uint8_t symbol = 0;
    bool is_leaf = false;
    Bitmap bits;
  };

  std::size_t uncompressed_size_ = 0;
  std::size_t root_ = kNpos;
  std::vector<Node> nodes_;
  std::vector<std::byte> compressed_;

  // --- construction --------------------------------------------------------

  /** @brief Build the Huffman tree and per-node bitmaps from @p input. */
  void build(std::span<const symbol_type> input) {
    uncompressed_size_ = input.size();
    if (uncompressed_size_ == 0) {
      root_ = kNpos;
      return;
    }

    std::array<std::size_t, kAlphabet> freq{};
    for (auto s : input) {
      freq[s]++;
    }

    struct HeapItem {
      std::size_t weight;
      std::size_t idx;
    };
    auto cmp = [](const HeapItem& a, const HeapItem& b) {
      return a.weight > b.weight;
    };
    std::priority_queue<HeapItem, std::vector<HeapItem>, decltype(cmp)> heap(
        cmp);

    // At most 2*kAlphabet-1 nodes (leaves + internal). Pre-allocate to avoid
    // reallocation during tree construction.
    nodes_.reserve(2 * kAlphabet - 1);

    for (std::size_t s = 0; s < kAlphabet; s++) {
      if (freq[s] > 0) {
        Node n;
        n.is_leaf = true;
        n.symbol = static_cast<std::uint8_t>(s);
        nodes_.push_back(std::move(n));
        heap.push({freq[s], nodes_.size() - 1});
      }
    }

    while (heap.size() > 1) {
      HeapItem a = heap.top();
      heap.pop();
      HeapItem b = heap.top();
      heap.pop();
      Node n;
      n.is_leaf = false;
      n.left = a.idx;
      n.right = b.idx;
      std::size_t w = a.weight + b.weight;
      n.bits.count = w;
      n.bits.words.assign((w + kWordBits - 1) / kWordBits, 0);
      nodes_.push_back(std::move(n));
      heap.push({w, nodes_.size() - 1});
    }
    root_ = heap.top().idx;

    // Record, for each symbol, the path of (node, direction-bit) from root.
    std::array<std::vector<std::pair<std::size_t, bool>>, kAlphabet> paths;
    std::vector<std::pair<std::size_t, bool>> path;
    assign_paths(root_, path, paths);

    // Top-down partitioning: recursively split the symbol stream and fill
    // pre-allocated bitmaps. Leaves are skipped (zero overhead).
    encode_node(root_, 0, input, paths);
  }

  /** @brief Depth-first assignment of root-to-leaf paths. */
  void assign_paths(
      std::size_t idx,
      std::vector<std::pair<std::size_t, bool>>& path,
      std::array<std::vector<std::pair<std::size_t, bool>>, kAlphabet>& paths) {
    if (nodes_[idx].is_leaf) {
      paths[nodes_[idx].symbol] = path;
      return;
    }
    path.emplace_back(idx, false);
    assign_paths(nodes_[idx].left, path, paths);
    path.pop_back();
    path.emplace_back(idx, true);
    assign_paths(nodes_[idx].right, path, paths);
    path.pop_back();
  }

  /** @brief Top-down encoding: recursively partition symbols and fill bitmaps.
   *  @details At each internal node, fill the pre-allocated bitmap with the
   *           direction bit of each symbol's code at this depth, and split the
   *           symbol stream into left (0) and right (1) child vectors. Leaves
   *           return immediately — zero overhead.
   *  @param idx     Current node index.
   *  @param depth   Current tree depth (indexes into each symbol's code path).
   *  @param symbols Symbol stream arriving at this node, in order.
   *  @param paths   Per-symbol code paths from `assign_paths`. */
  void encode_node(std::size_t idx,
                   std::size_t depth,
                   std::span<const symbol_type> symbols,
                   const std::array<std::vector<std::pair<std::size_t, bool>>,
                                    kAlphabet>& paths) {
    if (nodes_[idx].is_leaf) {
      return;
    }

    auto& n = nodes_[idx];
    std::vector<symbol_type> left_symbols, right_symbols;
    std::size_t left_size = nodes_[n.left].bits.count;
    std::size_t right_size = nodes_[n.right].bits.count;
    left_symbols.reserve(left_size);
    right_symbols.reserve(right_size);

    // Fill bitmap via incremental word/bit tracking — avoids per-symbol
    // division and modulo.
    std::size_t word_idx = 0;
    std::size_t bit_pos = 0;
    for (std::size_t i = 0; i < symbols.size(); i++) {
      bool bit = paths[symbols[i]][depth].second;
      if (bit) {
        n.bits.words[word_idx] |= (1ull << bit_pos);
        right_symbols.push_back(symbols[i]);
      } else {
        left_symbols.push_back(symbols[i]);
      }
      if (++bit_pos == kWordBits) {
        bit_pos = 0;
        ++word_idx;
      }
    }

    encode_node(n.left, depth + 1, left_symbols, paths);
    encode_node(n.right, depth + 1, right_symbols, paths);
  }

  // --- decode --------------------------------------------------------------

  /** @brief Reconstruct the original sequence by bottom-up merging. */
  std::vector<symbol_type> decode_from_tree() const {
    if (root_ == kNpos) {
      return std::vector<symbol_type>();
    }
    return decode_node(root_, uncompressed_size_);
  }

  /** @brief Recursively decode a node's output stream bottom-up.
   *  @param weight Number of symbols this node must produce. For internal
   *         nodes this equals `bits.count`; for leaves it is passed down from
   *         the parent (no bitmap to read it from). */
  std::vector<symbol_type> decode_node(std::size_t idx,
                                       std::size_t weight) const {
    const auto& n = nodes_[idx];
    if (n.is_leaf) {
      return std::vector<symbol_type>(weight, n.symbol);
    }

    std::size_t right_weight = n.bits.popcount();
    std::size_t left_weight = n.bits.count - right_weight;
    std::vector<symbol_type> left_symbols = decode_node(n.left, left_weight);
    std::vector<symbol_type> right_symbols = decode_node(n.right, right_weight);

    std::vector<symbol_type> out(n.bits.count);
    std::size_t c_left = 0;
    std::size_t c_right = 0;
    // Word-wise merge: read each 64-bit word once and shift through its bits,
    // avoiding per-bit division/modulo.
    for (std::size_t w = 0; w < n.bits.words.size(); w++) {
      std::uint64_t word = n.bits.words[w];
      std::size_t base = w * kWordBits;
      std::size_t limit = base + kWordBits;
      if (limit > n.bits.count) {
        limit = n.bits.count;
      }
      for (std::size_t i = base; i < limit; i++) {
        if (word & 1ull) {
          out[i] = right_symbols[c_right++];
        } else {
          out[i] = left_symbols[c_left++];
        }
        word >>= 1;
      }
    }
    return out;
  }

  // --- serialization -------------------------------------------------------

  /** @brief Write the in-memory tree into the serialized byte buffer. */
  void serialize() {
    compressed_.clear();
    write(uncompressed_size_);
    if (uncompressed_size_ == 0) {
      return;
    }
    write(root_);
    write(nodes_.size());
    for (const auto& n : nodes_) {
      write(static_cast<std::uint8_t>(n.is_leaf ? 1 : 0));
      write(n.symbol);
      write(n.left);
      write(n.right);
      write(n.bits.count);
      write_words(n.bits.words);
    }
  }

  /** @brief Rebuild the in-memory tree from a serialized byte buffer. */
  void deserialize(std::span<const std::byte> data) {
    compressed_.assign(data.begin(), data.end());
    nodes_.clear();
    if (data.empty()) {
      uncompressed_size_ = 0;
      root_ = kNpos;
      return;
    }
    std::size_t pos = 0;
    uncompressed_size_ = read<std::size_t>(data, pos);
    if (uncompressed_size_ == 0) {
      root_ = kNpos;
      return;
    }
    root_ = read<std::size_t>(data, pos);
    const std::size_t count = read<std::size_t>(data, pos);
    nodes_.resize(count);
    for (std::size_t i = 0; i < count; i++) {
      const std::uint8_t flags = read<std::uint8_t>(data, pos);
      nodes_[i].is_leaf = (flags & 1) != 0;
      nodes_[i].symbol = read<std::uint8_t>(data, pos);
      nodes_[i].left = read<std::size_t>(data, pos);
      nodes_[i].right = read<std::size_t>(data, pos);
      nodes_[i].bits.count = read<std::size_t>(data, pos);
      nodes_[i].bits.words = read_words(data, pos, nodes_[i].bits.word_count());
    }
  }

  /** @brief Append a little-endian, fixed-width value to the byte buffer. */
  template <class T>
  void write(const T& value) {
    const std::size_t old = compressed_.size();
    compressed_.resize(old + sizeof(T));
    std::memcpy(compressed_.data() + old, &value, sizeof(T));
  }

  /** @brief Append a packed word array directly to the byte buffer. */
  void write_words(const std::vector<std::uint64_t>& words) {
    const std::size_t bytes = words.size() * sizeof(std::uint64_t);
    const std::size_t old = compressed_.size();
    compressed_.resize(old + bytes);
    if (bytes > 0) {
      std::memcpy(compressed_.data() + old, words.data(), bytes);
    }
  }

  /** @brief Read a fixed-width value from @p data at @p pos. */
  template <class T>
  static T read(std::span<const std::byte> data, std::size_t& pos) {
    T value;
    std::memcpy(&value, data.data() + pos, sizeof(T));
    pos += sizeof(T);
    return value;
  }

  /** @brief Read a packed word array of @p word_count words from @p data. */
  static std::vector<std::uint64_t> read_words(std::span<const std::byte> data,
                                               std::size_t& pos,
                                               std::size_t word_count) {
    std::vector<std::uint64_t> words(word_count);
    const std::size_t bytes = word_count * sizeof(std::uint64_t);
    if (bytes > 0) {
      std::memcpy(words.data(), data.data() + pos, bytes);
      pos += bytes;
    }
    return words;
  }
};

}  // namespace pixie
