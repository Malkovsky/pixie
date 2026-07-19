#pragma once

/**
 * @file pivco_huffman.h
 * @brief Simple scalar PivCo-Huffman codec.
 *
 * Stores a Huffman-shaped tree of per-node routing bitmaps (the "tree of
 * bitmaps" layout shared with wavelet trees). Encoding walks each input symbol
 * from the root, appending one direction bit to every internal node on its
 * path. Decoding walks each output position from the root down to a leaf,
 * consuming one bit per internal node.
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
 *          packed 64-bit words, and decoding by top-down per-position
 *          traversal. The serialized form is a header followed by the tree
 *          structure and packed node bitmaps.
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

    /** @brief Append one bit to the end of the bitmap. */
    void push_back(bool bit) {
      if (count % kWordBits == 0) {
        words.push_back(0);
      }
      if (bit) {
        words.back() |= (1ull << (count % kWordBits));
      }
      ++count;
    }

    /** @brief Read the bit at position @p i. */
    bool get(std::size_t i) const {
      return ((words[i / kWordBits] >> (i % kWordBits)) & 1ull) != 0;
    }

    /** @brief Number of stored bits. */
    std::size_t size() const { return count; }

    /** @brief Number of 64-bit words backing the bitmap. */
    std::size_t word_count() const {
      return (count + kWordBits - 1) / kWordBits;
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
      nodes_.push_back(std::move(n));
      heap.push({a.weight + b.weight, nodes_.size() - 1});
    }
    root_ = heap.top().idx;

    // Record, for each symbol, the path of (node, direction-bit) from root.
    std::array<std::vector<std::pair<std::size_t, bool>>, kAlphabet> paths;
    std::vector<std::pair<std::size_t, bool>> path;
    assign_paths(root_, path, paths);

    // Append one direction bit per visited node, in input order.
    for (auto s : input) {
      for (const auto& [idx, bit] : paths[s]) {
        nodes_[idx].bits.push_back(bit);
      }
    }
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

  // --- decode --------------------------------------------------------------

  /** @brief Reconstruct the original sequence by top-down traversal. */
  std::vector<symbol_type> decode_from_tree() const {
    std::vector<symbol_type> out(uncompressed_size_);
    if (root_ == kNpos) {
      return out;
    }
    if (nodes_[root_].is_leaf) {
      const std::uint8_t sym = nodes_[root_].symbol;
      for (auto& s : out) {
        s = sym;
      }
      return out;
    }
    std::vector<std::size_t> cursor(nodes_.size(), 0);
    for (std::size_t i = 0; i < uncompressed_size_; i++) {
      std::size_t idx = root_;
      while (!nodes_[idx].is_leaf) {
        const bool bit = nodes_[idx].bits.get(cursor[idx]++);
        idx = bit ? nodes_[idx].right : nodes_[idx].left;
      }
      out[i] = nodes_[idx].symbol;
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
