#pragma once

/**
 * @file pivco_huffman.h
 * @brief Simple scalar PivCo-Huffman codec with canonical Huffman codes.
 *
 * Stores a Huffman-shaped tree of per-node routing bitmaps (the "tree of
 * bitmaps" layout shared with wavelet trees). Encoding walks each input symbol
 * from the root, appending one direction bit to every internal node on its
 * path. Decoding merges child symbol streams bottom-up from leaves to root,
 * using the node bitmap as a selector.
 *
 * The tree uses canonical Huffman code lengths (<= 15 bits, enforced via the
 * package-merge length-limited algorithm): the wire format stores only the 256
 * per-symbol lengths (128 bytes as 4-bit nibbles), and the tree structure is
 * reconstructed deterministically at load time. This groups same-length codes
 * together, which is a prerequisite for future flat-subtree optimizations and
 * reduces serialized header size.
 *
 * This is a deliberately simple, unoptimized reference implementation: node
 * bitmaps are stored as packed `std::uint64_t` words, traversal is scalar,
 * and there are no flat-subtree, SIMD, or selective-ANS optimizations.
 *
 * Range and ownership conventions follow `HuffmanBase`: symbols are bytes,
 * the compressed stream is a byte view, and the codec owns its serialized form.
 */

#include <pixie/huffman.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Simple scalar PivCo-Huffman codec with canonical Huffman codes.
 * @details Implements `HuffmanBase` by building a Huffman tree over the
 *          byte alphabet, canonicalizing it (grouping same-length codes),
 *          storing one routing bitmap per internal node as packed 64-bit
 *          words, and decoding by bottom-up merging of child symbol streams.
 *          The serialized form is 256 code lengths followed by packed node
 *          bitmaps; the tree is reconstructed from the lengths at load time.
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

  /// @brief Maximum supported Huffman code length (matches common decoders).
  static constexpr std::size_t kMaxCodeLength = 15;

  /**
   * @brief Lightweight descriptor for a node's routing bitmap.
   * @details The packed 64-bit words live in a single shared arena
   *          (`arena_`); this struct only records where the node's words begin
   *          (`offset`, in words) and how many symbols they describe
   *          (`count`). Splitting descriptor from storage lets one contiguous
   *          allocation back every internal node, eliminating per-node heap
   *          traffic during decode and improving cache locality.
   */
  struct Bitmap {
    std::size_t offset = 0;
    std::size_t count = 0;

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
  /// @brief One contiguous allocation backing every internal node's bitmap
  ///        words. Node `i` owns `nodes_[i].bits.word_count()` words starting
  ///        at `arena_[nodes_[i].bits.offset]`.
  std::vector<std::uint64_t> arena_;
  /// @brief Reusable decode scratch: two ping-pong halves of size
  ///        `uncompressed_size_`. Allocated once and reused across decode
  ///        calls; `mutable` because `decode_impl()` is const.
  mutable std::vector<symbol_type> workspace_;
  std::array<std::uint8_t, kAlphabet> code_lengths_{};

  // --- construction --------------------------------------------------------

  /** @brief Build the canonical Huffman tree and per-node bitmaps from @p
   * input.
   */
  void build(std::span<const symbol_type> input) {
    uncompressed_size_ = input.size();
    if (uncompressed_size_ == 0) {
      root_ = kNpos;
      return;
    }

    // 1. Count symbol frequencies.
    std::array<std::size_t, kAlphabet> freq{};
    for (auto s : input) {
      freq[s]++;
    }

    // 2. Build a Huffman tree to determine code lengths, then extract lengths.
    compute_code_lengths(freq, code_lengths_);

    // 3. Rebuild a canonical tree from the lengths (groups same-length codes).
    //    Structure only: bitmap weights/offsets are assigned below, not here.
    build_canonical_tree(code_lengths_);

    // Single-symbol input: root is a leaf, no internal nodes, no arena.
    if (root_ == kNpos || nodes_[root_].is_leaf) {
      return;
    }

    // 4. Assign exact per-node weights (subtree frequency sums) and arena
    //    offsets in a single post-order walk, then allocate the arena once.
    std::size_t arena_words = 0;
    assign_weights_and_offsets(root_, freq, arena_words);
    arena_.assign(arena_words, 0);

    // 5. Assign per-symbol code paths, then fill node bitmaps via top-down
    //    in-place partition into the reusable workspace. This mirrors decode:
    //    the root reads its input from workspace half 0, writes its bitmap,
    //    and partitions symbols into half 1 for its children. Each level
    //    ping-pongs between the two halves — no per-node allocations.
    std::array<std::vector<std::pair<std::size_t, bool>>, kAlphabet> paths;
    std::vector<std::pair<std::size_t, bool>> path;
    assign_paths(root_, path, paths);
    workspace_.resize(uncompressed_size_ * 2);
    std::copy(input.begin(), input.end(), workspace_.data());
    encode_partition(root_, 0, 0, paths, freq);
  }

  /** @brief Compute length-limited Huffman code lengths via package-merge.
   *  @details Implements the Larmore-Hirschberg package-merge algorithm to
   *           produce optimal code lengths subject to `length <=
   * kMaxCodeLength`. This guarantees the 4-bit nibble header never truncates a
   * length, and keeps all downstream canonical-tree, flat-subtree, and
   *           reshaping optimizations valid.
   *
   *  Algorithm: maintain a sorted list of "items" per level (1..L). A level-1
   *  item is a single coin (one symbol). At each subsequent level, form
   *  "packages" by pairing consecutive items from the previous level (a
   *  package commits to selecting both its constituents), then merge with the
   *  original coins and keep the cheapest `2n-2` items. After L-1 iterations,
   *  summing each symbol's count across the final `2n-2` items yields its code
   *  length (each ≤ L, and the Kraft inequality holds exactly).
   *
   *  Each item carries a per-symbol count vector (`Counts`) so that summing at
   *  the end is a flat array reduction; with `n <= 256` and `L = 15` the total
   *  work is ~2M additions (< 1 ms), negligible relative to the bitmap fill. */
  void compute_code_lengths(const std::array<std::size_t, kAlphabet>& freq,
                            std::array<std::uint8_t, kAlphabet>& lengths) {
    lengths.fill(0);

    // Collect present symbols, sorted by weight ascending (the coin list S,
    // reused at every level of the package-merge).
    struct Coin {
      std::uint8_t symbol;
      std::size_t weight;
    };
    std::vector<Coin> coins;
    coins.reserve(kAlphabet);
    for (std::size_t s = 0; s < kAlphabet; s++) {
      if (freq[s] > 0) {
        coins.push_back({static_cast<std::uint8_t>(s), freq[s]});
      }
    }
    std::sort(coins.begin(), coins.end(),
              [](const Coin& a, const Coin& b) { return a.weight < b.weight; });

    const std::size_t n = coins.size();
    if (n == 0) {
      return;
    }
    // Single distinct symbol: must still get a 1-bit code to be decodable.
    if (n == 1) {
      lengths[coins[0].symbol] = 1;
      return;
    }

    // Per-item symbol-count vector. uint16 suffices: a symbol's count within a
    // single package can grow across levels, but final lengths are <= L = 15.
    using Counts = std::array<std::uint16_t, kAlphabet>;
    struct Item {
      std::size_t weight;
      Counts counts;
    };

    // Level-1 list: one coin per symbol.
    std::vector<Item> cur;
    cur.reserve(2 * n);
    for (const auto& c : coins) {
      Item it{};
      it.weight = c.weight;
      it.counts[c.symbol] = 1;
      cur.push_back(std::move(it));
    }

    // Keep the cheapest 2n-2 items at each level (the final selection size;
    // cheaper items at a level can never displace a more expensive one in an
    // optimal solution, so pruning to 2n-2 is safe).
    const std::size_t keep = 2 * n - 2;

    // Iterate levels 2..L. Coins are reused at every level via the merge.
    for (std::size_t level = 1; level < kMaxCodeLength; level++) {
      // Package: pair consecutive items of `cur` (already weight-sorted).
      std::vector<Item> packages;
      packages.reserve(cur.size() / 2);
      for (std::size_t i = 0; i + 1 < cur.size(); i += 2) {
        Item pkg{};
        pkg.weight = cur[i].weight + cur[i + 1].weight;
        for (std::size_t j = 0; j < kAlphabet; j++) {
          pkg.counts[j] = cur[i].counts[j] + cur[i + 1].counts[j];
        }
        packages.push_back(std::move(pkg));
      }
      // Merge coins (sorted) with packages (sorted — pairing preserves order),
      // keeping the cheapest `keep` items. Ties prefer coins (deterministic).
      // At early levels the combined list may hold fewer than `keep` items;
      // keep all of them in that case (the list grows toward `keep` over
      // levels as packages proliferate).
      const std::size_t available = n + packages.size();
      const std::size_t keep_here = std::min(keep, available);
      std::vector<Item> next;
      next.reserve(keep_here);
      std::size_t ci = 0;
      std::size_t pi = 0;
      while (next.size() < keep_here) {
        bool use_coin;
        if (ci >= n) {
          use_coin = false;
        } else if (pi >= packages.size()) {
          use_coin = true;
        } else {
          use_coin = coins[ci].weight <= packages[pi].weight;
        }
        if (use_coin) {
          Item it{};
          it.weight = coins[ci].weight;
          it.counts[coins[ci].symbol] = 1;
          next.push_back(std::move(it));
          ci++;
        } else {
          next.push_back(std::move(packages[pi]));
          pi++;
        }
      }
      cur = std::move(next);
    }

    // Sum symbol counts across the final 2n-2 selected items -> code lengths.
    Counts total{};
    for (const auto& it : cur) {
      for (std::size_t j = 0; j < kAlphabet; j++) {
        total[j] += it.counts[j];
      }
    }
    for (const auto& c : coins) {
      lengths[c.symbol] = static_cast<std::uint8_t>(total[c.symbol]);
    }
  }

  /** @brief Rebuild a canonical Huffman tree from code lengths.
   *  @details Assigns canonical codes (sorted by length, then symbol), then
   *           builds the tree top-down by inserting each (symbol, code) pair.
   *           Same-length codes are grouped together, which is a prerequisite
   *           for future flat-subtree optimization. */
  void build_canonical_tree(
      const std::array<std::uint8_t, kAlphabet>& lengths) {
    nodes_.clear();
    nodes_.reserve(2 * kAlphabet);

    // Collect symbols present, sorted by (length, symbol).
    struct SymbolEntry {
      std::uint8_t symbol;
      std::uint8_t length;
    };
    std::vector<SymbolEntry> entries;
    entries.reserve(kAlphabet);
    for (std::size_t s = 0; s < kAlphabet; s++) {
      if (lengths[s] > 0) {
        entries.push_back({static_cast<std::uint8_t>(s), lengths[s]});
      }
    }
    std::sort(entries.begin(), entries.end(),
              [](const SymbolEntry& a, const SymbolEntry& b) {
                if (a.length != b.length) {
                  return a.length < b.length;
                }
                return a.symbol < b.symbol;
              });

    // Special case: single distinct symbol → root is a leaf (no internal nodes,
    // no bitmaps). decode_node returns a constant-filled vector.
    if (entries.size() == 1) {
      nodes_.push_back(Node{});
      root_ = 0;
      nodes_[0].is_leaf = true;
      nodes_[0].symbol = entries[0].symbol;
      return;
    }

    // Assign canonical codes: first code at each length is 0, increment
    // within a length, shift left when length increases.
    std::vector<std::pair<std::uint64_t, std::uint8_t>>
        codes;  // (code, length)
    codes.reserve(entries.size());
    std::uint64_t code = 0;
    std::uint8_t prev_len = 0;
    for (const auto& e : entries) {
      if (prev_len == 0) {
        prev_len = e.length;
      }
      while (prev_len < e.length) {
        code <<= 1;
        prev_len++;
      }
      codes.push_back({code, e.length});
      code++;
    }

    // Build tree by inserting each (symbol, code) pair via bit-by-bit descent.
    // Create root.
    nodes_.push_back(Node{});
    root_ = 0;

    for (std::size_t i = 0; i < entries.size(); i++) {
      std::uint64_t c = codes[i].first;
      std::uint8_t len = codes[i].second;
      std::size_t cur = root_;

      // Descend len-1 bits, creating internal nodes as needed.
      for (std::size_t bit = 0; bit + 1 < len; bit++) {
        // Read MSB first (canonical codes are assigned MSB-first).
        bool go_right = (c >> (len - 1 - bit)) & 1;
        std::size_t& child = go_right ? nodes_[cur].right : nodes_[cur].left;
        if (child == kNpos) {
          nodes_.push_back(Node{});
          child = nodes_.size() - 1;
        }
        cur = child;
      }
      // Last bit → leaf.
      bool go_right = (c >> 0) & 1;
      std::size_t& child = go_right ? nodes_[cur].right : nodes_[cur].left;
      nodes_.push_back(Node{});
      child = nodes_.size() - 1;
      nodes_[child].is_leaf = true;
      nodes_[child].symbol = entries[i].symbol;
    }

    // Bitmap weights and arena offsets are assigned by the caller: `build()`
    // derives them from symbol frequencies (exact), `deserialize()` reads
    // them from the wire. This keeps the tree builder free of bitmap storage.
  }

  /** @brief Post-order assignment of exact per-node weights and arena offsets.
   *  @details The weight of an internal node is the sum of leaf frequencies in
   *           its subtree (exact, not a parent-derived upper bound). Assigns
   *           `bits.count` and `bits.offset` for each internal node and
   *           accumulates the total arena word count. Returns the subtree
   *           weight so each parent can sum its two children. */
  std::size_t assign_weights_and_offsets(
      std::size_t idx,
      const std::array<std::size_t, kAlphabet>& freq,
      std::size_t& arena_words) {
    if (idx == kNpos) {
      return 0;
    }
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return freq[n.symbol];
    }
    std::size_t left_weight =
        assign_weights_and_offsets(n.left, freq, arena_words);
    std::size_t right_weight =
        assign_weights_and_offsets(n.right, freq, arena_words);
    n.bits.count = left_weight + right_weight;
    n.bits.offset = arena_words;
    arena_words += n.bits.word_count();
    return n.bits.count;
  }

  /** @brief Exact weight of a child node, used to place partitions.
   *  @details Internal nodes carry their precomputed `bits.count`; leaves
   *           carry their symbol frequency (no bitmap stored). */
  std::size_t child_weight(
      std::size_t idx,
      const std::array<std::size_t, kAlphabet>& freq) const {
    const Node& n = nodes_[idx];
    return n.is_leaf ? freq[n.symbol] : n.bits.count;
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

  /** @brief Top-down in-place partition encoding into the reusable workspace.
   *  @details The dual of decode: at each internal node, read the incoming
   *           symbol stream from workspace half `depth`, fill the node's
   *           arena bitmap with direction bits, and partition symbols into
   *           half `depth+1` (left at `out_base`, right at
   *           `out_base + left_weight`). Leaves return immediately. No
   *           per-node allocations — the workspace is pre-sized and
   *           ping-ponged by depth, exactly mirroring decode.
   *  @param idx      Current node index.
   *  @param depth    Current tree depth (selects workspace half + path index).
   *  @param out_base Offset within the current half where this node's stream
   *                  begins.
   *  @param paths    Per-symbol code paths from `assign_paths`.
   *  @param freq     Symbol frequencies (sizing child partitions). */
  void encode_partition(
      std::size_t idx,
      std::size_t depth,
      std::size_t out_base,
      const std::array<std::vector<std::pair<std::size_t, bool>>, kAlphabet>&
          paths,
      const std::array<std::size_t, kAlphabet>& freq) {
    if (nodes_[idx].is_leaf) {
      return;
    }

    Node& n = nodes_[idx];
    const std::size_t left_weight = child_weight(n.left, freq);
    const std::size_t weight = n.bits.count;

    const symbol_type* src = workspace_half(depth) + out_base;
    symbol_type* left_dst = workspace_half(depth + 1) + out_base;
    symbol_type* right_dst = workspace_half(depth + 1) + out_base + left_weight;

    std::uint64_t* bits = arena_.data() + n.bits.offset;
    std::size_t word_idx = 0;
    std::size_t bit_pos = 0;
    std::size_t c_left = 0;
    std::size_t c_right = 0;
    for (std::size_t i = 0; i < weight; i++) {
      symbol_type s = src[i];
      if (paths[s][depth].second) {
        bits[word_idx] |= (1ull << bit_pos);
        right_dst[c_right++] = s;
      } else {
        left_dst[c_left++] = s;
      }
      if (++bit_pos == kWordBits) {
        bit_pos = 0;
        ++word_idx;
      }
    }

    encode_partition(n.left, depth + 1, out_base, paths, freq);
    encode_partition(n.right, depth + 1, out_base + left_weight, paths, freq);
  }

  // --- decode --------------------------------------------------------------

  /** @brief Reconstruct the original sequence by bottom-up merging.
   *  @details Uses a single reusable workspace buffer (sized to twice the
   *           uncompressed length, allocated once and reused across calls)
   *           and ping-pongs between its two halves by tree depth. Each node
   *           reads its children's outputs from one half and writes its merged
   *           output to the other, so no per-node allocation occurs after the
   *           workspace is sized. The final result is copied out of half 0. */
  std::vector<symbol_type> decode_from_tree() const {
    if (root_ == kNpos) {
      return std::vector<symbol_type>();
    }
    if (nodes_[root_].is_leaf) {
      return std::vector<symbol_type>(uncompressed_size_, nodes_[root_].symbol);
    }
    const std::size_t n = uncompressed_size_;
    workspace_.resize(n * 2);
    decode_node(root_, n, 0, 0);
    return std::vector<symbol_type>(workspace_.begin(), workspace_.begin() + n);
  }

  /** @brief Recursively decode a node's output stream bottom-up into the
   *         reusable workspace.
   *  @param weight   Number of symbols this node must produce.
   *  @param depth     Tree depth, selecting which workspace half receives this
   *                   node's output; children write to the opposite half.
   *  @param out_base  Offset within the selected half where this node's output
   *                   begins. Children are placed at `out_base` (left) and
   *                   `out_base + left_weight` (right) in the opposite half,
   *                   then merged into `[out_base, out_base + weight)` here. */
  void decode_node(std::size_t idx,
                   std::size_t weight,
                   std::size_t depth,
                   std::size_t out_base) const {
    const Node& n = nodes_[idx];
    symbol_type* dst = workspace_half(depth) + out_base;
    if (n.is_leaf) {
      std::fill(dst, dst + weight, n.symbol);
      return;
    }

    const std::size_t right_weight = bitmap_popcount(n.bits);
    const std::size_t left_weight = n.bits.count - right_weight;
    decode_node(n.left, left_weight, depth + 1, out_base);
    decode_node(n.right, right_weight, depth + 1, out_base + left_weight);

    // Children wrote to the opposite half: [out_base, out_base + left_weight)
    // holds the left stream, [out_base + left_weight, out_base + weight)
    // holds the right stream. Merge both into this node's half by the bitmap.
    const symbol_type* src = workspace_half(depth + 1) + out_base;
    const symbol_type* left_out = src;
    const symbol_type* right_out = src + left_weight;

    std::size_t out_i = 0;
    std::size_t c_left = 0;
    std::size_t c_right = 0;
    const std::uint64_t* bits = arena_.data() + n.bits.offset;
    const std::size_t words = n.bits.word_count();
    for (std::size_t w = 0; w < words; w++) {
      std::uint64_t word = bits[w];
      std::size_t limit = kWordBits;
      if (w + 1 == words) {
        limit = n.bits.count - w * kWordBits;
      }
      for (std::size_t i = 0; i < limit; i++) {
        dst[out_i++] =
            (word & 1ull) ? right_out[c_right++] : left_out[c_left++];
        word >>= 1;
      }
    }
  }

  /** @brief Pointer to the workspace half that holds depth-@p depth outputs. */
  symbol_type* workspace_half(std::size_t depth) const {
    return workspace_.data() + (depth % 2) * uncompressed_size_;
  }

  /** @brief Number of set bits in a node's arena-backed bitmap. */
  std::size_t bitmap_popcount(const Bitmap& b) const {
    const std::uint64_t* p = arena_.data() + b.offset;
    std::size_t total = 0;
    const std::size_t words = b.word_count();
    for (std::size_t i = 0; i < words; i++) {
      total += std::popcount(p[i]);
    }
    return total;
  }

  // --- serialization -------------------------------------------------------

  /** @brief Write the in-memory tree into the serialized byte buffer.
   *  @details Wire format:
   *           - uncompressed_size (8 bytes)
   *           - 256 code lengths as 128 bytes of 4-bit nibbles (0 = absent)
   *           - For each internal node (pre-order): bits.count (8 bytes) +
   *             packed words. Leaves are implied by the tree structure
   *             reconstructed from lengths. */
  void serialize() {
    compressed_.clear();
    write(uncompressed_size_);
    if (uncompressed_size_ == 0) {
      return;
    }

    // Pack stored code lengths as nibbles: 2 symbols per byte.
    for (std::size_t i = 0; i < kAlphabet; i += 2) {
      std::uint8_t byte = static_cast<std::uint8_t>(
          code_lengths_[i] | (code_lengths_[i + 1] << 4));
      write(byte);
    }

    // Write internal node bitmaps in pre-order (matches deserialize order).
    serialize_node(root_);
  }

  /** @brief Serialize a node's bitmap and recurse (pre-order). */
  void serialize_node(std::size_t idx) {
    if (idx == kNpos) {
      return;
    }
    const auto& n = nodes_[idx];
    if (!n.is_leaf) {
      write(n.bits.count);
      write_words(n.bits.offset, n.bits.word_count());
      serialize_node(n.left);
      serialize_node(n.right);
    }
  }

  /** @brief Rebuild the in-memory tree from a serialized byte buffer. */
  void deserialize(std::span<const std::byte> data) {
    compressed_.assign(data.begin(), data.end());
    nodes_.clear();
    arena_.clear();
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

    // Read 256 code lengths from nibbles.
    for (std::size_t i = 0; i < kAlphabet; i += 2) {
      std::uint8_t byte = read<std::uint8_t>(data, pos);
      code_lengths_[i] = byte & 0x0F;
      code_lengths_[i + 1] = (byte >> 4) & 0x0F;
    }

    // Rebuild the canonical tree from lengths.
    build_canonical_tree(code_lengths_);

    // Single-symbol input: root is a leaf, no internal nodes, no arena.
    if (root_ == kNpos || nodes_[root_].is_leaf) {
      return;
    }

    // Pass 1: read each internal node's `count` (and skip its word bytes on
    // the wire), assign contiguous arena offsets, and size the arena exactly.
    // `scan_pos` is discarded — pass 2 re-walks the same region to read words.
    const std::size_t bitmap_start = pos;
    std::size_t arena_words = 0;
    scan_node_counts(root_, data, pos, arena_words);
    arena_.assign(arena_words, 0);

    // Pass 2: rewind to the bitmap section and read each node's packed words
    // into its arena slot.
    pos = bitmap_start;
    read_node_words(root_, data, pos);
  }

  /** @brief Pass 1: read each internal node's `count` from the wire, assign its
   *         arena offset, and advance past its word bytes without storing
   *         them. Sums the total arena word count so the caller allocates once.
   *  @details Pre-order traversal matches the serialization order. */
  void scan_node_counts(std::size_t idx,
                        std::span<const std::byte> data,
                        std::size_t& pos,
                        std::size_t& arena_words) {
    if (idx == kNpos) {
      return;
    }
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }
    n.bits.count = read<std::size_t>(data, pos);
    n.bits.offset = arena_words;
    arena_words += n.bits.word_count();
    pos += n.bits.word_count() * sizeof(std::uint64_t);
    scan_node_counts(n.left, data, pos, arena_words);
    scan_node_counts(n.right, data, pos, arena_words);
  }

  /** @brief Pass 2: read each internal node's packed words from the wire into
   *         its arena slot (pre-order, matching serialization). */
  void read_node_words(std::size_t idx,
                       std::span<const std::byte> data,
                       std::size_t& pos) {
    if (idx == kNpos) {
      return;
    }
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }
    pos += sizeof(std::size_t);  // count, already read in pass 1
    const std::size_t words = n.bits.word_count();
    const std::size_t bytes = words * sizeof(std::uint64_t);
    if (bytes > 0) {
      std::memcpy(arena_.data() + n.bits.offset, data.data() + pos, bytes);
    }
    pos += bytes;
    read_node_words(n.left, data, pos);
    read_node_words(n.right, data, pos);
  }

  /** @brief Append a little-endian, fixed-width value to the byte buffer. */
  template <class T>
  void write(const T& value) {
    const std::size_t old = compressed_.size();
    compressed_.resize(old + sizeof(T));
    std::memcpy(compressed_.data() + old, &value, sizeof(T));
  }

  /** @brief Append a node's arena-backed words directly to the byte buffer. */
  void write_words(std::size_t offset, std::size_t word_count) {
    const std::size_t bytes = word_count * sizeof(std::uint64_t);
    const std::size_t old = compressed_.size();
    compressed_.resize(old + bytes);
    if (bytes > 0) {
      std::memcpy(compressed_.data() + old, arena_.data() + offset, bytes);
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
};

}  // namespace pixie
