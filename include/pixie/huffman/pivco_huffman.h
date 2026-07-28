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
 * The tree uses canonical Huffman code lengths (<= 15 bits, enforced with a
 * zlib-style overflow correction): the wire format stores one file-level set
 * of 256 per-symbol lengths (128 bytes as 4-bit nibbles), and the tree
 * structure is reconstructed once at load time and reused by every block.
 * Same-length symbols are reshaped into flat subtrees as described by the
 * PivCo-Huffman paper.
 *
 * This implementation uses packed `std::uint64_t` bitmaps, flat subtrees,
 * contiguous arenas, and bottom-up decoding. The merge/partition primitives are
 * vectorized in `pivco_simd.h` (AVX-512 VBMI2, AVX2, and scalar fallbacks);
 * selective ANS remains future work.
 *
 * Range and ownership conventions follow `HuffmanBase`: symbols are bytes,
 * the compressed stream is a byte view, and the codec owns its serialized form.
 */

#include <pixie/huffman.h>
#include <pixie/huffman/pivco_simd.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief Simple scalar PivCo-Huffman codec with canonical Huffman codes.
 * @details Implements `HuffmanBase` by building a Huffman tree over the
 *          byte alphabet, canonicalizing it (grouping same-length codes),
 *          storing one routing bitmap per internal node as packed 64-bit
 *          words, and decoding by bottom-up merging of child symbol streams.
 *          The serialized form is one file-level set of 256 code lengths
 *          followed by block-local packed node bitmaps; the shared tree is
 *          reconstructed from the lengths at load time.
 */
class PivCoHuffman : public HuffmanBase<PivCoHuffman> {
 public:
  /** @brief Symbol type handled by the codec: one byte per symbol. */
  using symbol_type = std::uint8_t;

  /**
   * @brief Build a codec by encoding @p input.
   * @param input Symbol sequence to compress.
   */
  explicit PivCoHuffman(std::span<const symbol_type> input) { build(input); }

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

  /** @brief Reconstruct the original symbol sequence by decoding each block. */
  std::vector<symbol_type> decode_impl() const {
    if (uncompressed_size_ == 0) {
      return std::vector<symbol_type>();
    }
    // SIMD merge kernels may touch one vector past their logical output. Keep
    // the same slack as the internal workspace; resize it away without copying
    // after the last block has decoded.
    std::vector<symbol_type> result(uncompressed_size_ + kWorkspaceSlack);
    std::size_t pos = blocks_offset_;
    std::size_t output_offset = 0;
    while (output_offset < uncompressed_size_) {
      deserialize_block_payload(pos);
      decode_block_into(std::span<symbol_type>(result).subspan(
          output_offset, block_uncompressed_size_));
      output_offset += block_uncompressed_size_;
    }
    result.resize(uncompressed_size_);
    return result;
  }

 private:
  /// @brief Sentinel node index meaning "no node".
  static constexpr std::size_t kNpos = static_cast<std::size_t>(-1);

  /// @brief Byte-alphabet size (256 symbols).
  static constexpr std::size_t kAlphabet = 256;

  /// @brief Bits per word used by the packed bitmap.
  static constexpr std::size_t kWordBits = 64;

  /// @brief Maximum supported Huffman code length (matches common decoders).
  static constexpr std::size_t kMaxCodeLength = 15;

  /// @brief Maximum depth of a flat subtree (table size = 2^depth ≤ 256).
  /// @details A flat subtree replaces `2^d - 1` internal nodes (each with a
  ///          1-bit-per-symbol bitmap and a merge pass) with a single node
  ///          storing `d` bits per symbol and a lookup table. This turns `d`
  ///          decode merge passes into a single table-lookup pass.
  static constexpr std::uint8_t kMaxFlatDepth = 8;

  /// @brief Maximum nodes in a full binary tree over the byte alphabet.
  static constexpr std::size_t kMaxTreeNodes = 2 * kAlphabet - 1;

  /// @brief Default block size for block-based processing (64 KiB).
  /// @details All blocks share one file-level Huffman tree but carry their own
  ///          routing bitmaps. This keeps the decode workspace (2 × block_size
  ///          = 128 KiB) within the per-core L2 while avoiding per-block tree
  ///          construction and serialization overhead.
  static constexpr std::size_t kBlockSize = 64 * 1024;

  /// @brief Extra bytes appended to the ping-pong workspace beyond `2 * block`.
  /// @details Vectorized merge/partition kernels read and write a full 8- or
  ///          64-byte vector on each iteration even when the final group is
  ///          shorter; this slack makes those over-reads/over-writes stay
  ///          within the allocated buffer rather than spilling past the end of
  ///          a child region that abuts the workspace boundary.
  static constexpr std::size_t kWorkspaceSlack = 64;

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
    std::size_t count = 0;  ///< Number of symbols (not bits).

    /** @brief Number of 64-bit words backing the bitmap.
     *  @param flat_depth 0 for normal nodes (1 bit/symbol); >0 for flat
     *         nodes (`flat_depth` bits/symbol). */
    std::size_t word_count(std::uint8_t flat_depth = 0) const {
      std::size_t bits_per = flat_depth > 0 ? flat_depth : 1;
      return (count * bits_per + kWordBits - 1) / kWordBits;
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
    /// @brief 0 = normal node (1 bit/symbol); >0 = flat subtree of this depth
    ///        (stores `flat_depth` bits/symbol + a lookup table).
    std::uint8_t flat_depth = 0;
    /// @brief Whether this flat table is the balanced bit-reversal mapping.
    bool flat_bitreverse = false;
    /// @brief Offset of this flat node's decode table in `flat_tables_`.
    std::size_t flat_table_offset = 0;
  };

  /// @brief A symbol with its frequency weight, used for Huffman construction.
  struct Leaf {
    std::size_t weight;
    std::uint8_t symbol;
  };

  /** @brief Compact root-to-leaf code used by the scalar encoder.
   *  @details Codes are stored MSB-first in the low `length` bits. The
   *           15-bit length limit keeps every code in one 16-bit value and
   *           avoids the per-symbol dynamic path vectors used previously. */
  struct SymbolCode {
    std::uint16_t bits = 0;
    std::uint8_t length = 0;
  };

  using DirectionTables =
      std::array<std::array<std::uint8_t, kAlphabet>, kMaxCodeLength>;
  using FrequencyTable = std::array<std::size_t, kAlphabet>;

  /// @brief Total uncompressed size across all blocks (set during
  ///        build/deserialize, read by `uncompressed_size_impl()`).
  std::size_t uncompressed_size_ = 0;
  /// @brief Per-block serialized size from the wire (fixed; last block may
  ///        be smaller). Used by decode to iterate blocks.
  std::size_t block_size_ = kBlockSize;
  /// @brief Serialized compressed stream (header + per-block data).
  std::vector<std::byte> compressed_;
  /// @brief Byte offset of the first block payload in `compressed_`.
  std::size_t blocks_offset_ = 0;

  // --- shared tree and current-block state --------------------------------
  //    The topology and flat tables are file-level. Bitmap counts, offsets,
  //    and words are replaced for each block during const decoding.
  mutable std::size_t root_ = kNpos;
  mutable std::vector<Node> nodes_;
  /// @brief One contiguous allocation backing every internal node's bitmap
  ///        words. Node `i` owns `nodes_[i].bits.word_count()` words starting
  ///        at `arena_[nodes_[i].bits.offset]`.
  mutable std::vector<std::uint64_t> arena_;
  /// @brief Reusable decode/encode scratch: two ping-pong halves of size
  ///        `block_uncompressed_size_`. Reused across blocks; `mutable`
  ///        because `decode_impl()` is const.
  mutable std::vector<symbol_type> workspace_;
  /// @brief Uncompressed size of the current block being processed.
  mutable std::size_t block_uncompressed_size_ = 0;
  mutable std::array<std::uint8_t, kAlphabet> code_lengths_{};
  /// @brief Contiguous storage for all flat-subtree decode tables.
  mutable std::vector<symbol_type> flat_tables_;

  // --- construction --------------------------------------------------------

  /** @brief Build and serialize a shared tree and all block payloads.
   *  @details Counts each block once, combines those counts to build one
   *           file-level Huffman tree, then reuses its topology and direction
   *           tables while encoding every block. The wire format is:
   *           - total_uncompressed_size (8 bytes)
   *           - block_size (8 bytes)
   *           - 256 file-level code lengths as 128 bytes of nibbles
   *           - Per block: block_uncompressed_size (8 bytes) followed by
   *             per-internal-node bitmaps. */
  void build(std::span<const symbol_type> input) {
    uncompressed_size_ = input.size();
    compressed_.clear();
    write(uncompressed_size_);
    write(block_size_);
    blocks_offset_ = compressed_.size();
    if (input.empty()) {
      return;
    }

    const std::size_t block_count =
        (input.size() + block_size_ - 1) / block_size_;
    std::vector<FrequencyTable> block_frequencies;
    block_frequencies.reserve(block_count);
    FrequencyTable file_frequencies{};

    std::size_t offset = 0;
    while (offset < input.size()) {
      const std::size_t this_block =
          std::min(input.size() - offset, block_size_);
      block_frequencies.push_back(
          count_frequencies(input.subspan(offset, this_block)));
      const FrequencyTable& frequencies = block_frequencies.back();
      for (std::size_t symbol = 0; symbol < kAlphabet; ++symbol) {
        file_frequencies[symbol] += frequencies[symbol];
      }
      offset += this_block;
    }

    compute_code_lengths(file_frequencies, code_lengths_);
    build_noncanonical_tree(code_lengths_);
    detect_flat_subtrees();
    serialize_code_lengths();
    blocks_offset_ = compressed_.size();

    std::array<SymbolCode, kAlphabet> codes{};
    DirectionTables directions{};
    if (root_ != kNpos && !nodes_[root_].is_leaf) {
      assign_codes(root_, 0, 0, codes);
      directions = build_direction_tables(codes);
    }

    offset = 0;
    std::size_t block_index = 0;
    while (offset < input.size()) {
      const std::size_t this_block =
          std::min(input.size() - offset, block_size_);
      build_block_payload(input.subspan(offset, this_block),
                          block_frequencies[block_index], directions);
      serialize_block_payload();
      offset += this_block;
      ++block_index;
    }
  }

  /** @brief Build one block's bitmaps using the file-level tree. */
  void build_block_payload(std::span<const symbol_type> input,
                           const FrequencyTable& freq,
                           const DirectionTables& directions) {
    block_uncompressed_size_ = input.size();
    if (block_uncompressed_size_ == 0) {
      return;
    }

    // Single-symbol input: root is a leaf, no internal nodes, no arena.
    if (root_ == kNpos || nodes_[root_].is_leaf) {
      return;
    }

    // Assign exact per-node weights (subtree frequency sums) and arena
    //    offsets in a single post-order walk, then allocate the arena once.
    std::size_t arena_words = 0;
    assign_weights_and_offsets(root_, freq, arena_words);
    arena_.assign(arena_words, 0);

    // Fill node bitmaps via top-down in-place partition into the reusable
    //    workspace. This mirrors decode:
    //    the root reads its input from workspace half 0, writes its bitmap,
    //    and partitions symbols into half 1 for its children. Each level
    //    ping-pongs between the two halves — no per-node allocations.
    workspace_.resize(block_uncompressed_size_ * 2 + kWorkspaceSlack);
    std::copy(input.begin(), input.end(), workspace_.data());
    encode_partition(root_, 0, 0, directions, freq);
  }

  /** @brief Count a block's byte frequencies with independent accumulators.
   *  @details Eight partial histograms break the dependency chain caused by
   *           repeatedly incrementing a hot symbol. The block is at most
   *           64 KiB, so 32-bit partial counters cannot overflow. */
  static FrequencyTable count_frequencies(std::span<const symbol_type> input) {
    constexpr std::size_t kLanes = 8;
    std::array<std::array<std::uint32_t, kAlphabet>, kLanes> partial{};
    std::size_t i = 0;
    for (; i + kLanes <= input.size(); i += kLanes) {
      partial[0][input[i]]++;
      partial[1][input[i + 1]]++;
      partial[2][input[i + 2]]++;
      partial[3][input[i + 3]]++;
      partial[4][input[i + 4]]++;
      partial[5][input[i + 5]]++;
      partial[6][input[i + 6]]++;
      partial[7][input[i + 7]]++;
    }
    for (; i < input.size(); ++i) {
      partial[0][input[i]]++;
    }

    FrequencyTable frequencies{};
    for (std::size_t symbol = 0; symbol < kAlphabet; ++symbol) {
      for (std::size_t lane = 0; lane < kLanes; ++lane) {
        frequencies[symbol] += partial[lane][symbol];
      }
    }
    return frequencies;
  }

  /** @brief Compute Huffman code lengths, limited to ≤ kMaxCodeLength.
   *  @details Two-phase approach (same as zlib/zstd):
   *           1. Sort present symbols by frequency and build a standard
   *              Huffman tree with the linear two-queue algorithm, then
   *              extract per-leaf depths by walking parent pointers.
   *           2. If any depth exceeds kMaxCodeLength, apply the zlib-style
   *              "overflow fixup": collapse all lengths > L into L, then
   *              rebalance the Kraft sum by borrowing from shorter codes.
   *
   *  This is O(n log n) for sorting and O(n + L) after that, with fixed-size
   *  storage for the byte alphabet. */
  void compute_code_lengths(const std::array<std::size_t, kAlphabet>& freq,
                            std::array<std::uint8_t, kAlphabet>& lengths) {
    lengths.fill(0);

    std::array<Leaf, kAlphabet> leaves{};
    std::size_t leaf_count = 0;
    for (std::size_t s = 0; s < kAlphabet; s++) {
      if (freq[s] > 0) {
        leaves[leaf_count++] = {freq[s], static_cast<std::uint8_t>(s)};
      }
    }
    if (leaf_count == 0) {
      return;
    }
    if (leaf_count == 1) {
      lengths[leaves[0].symbol] = 1;
      return;
    }

    std::sort(leaves.begin(), leaves.begin() + leaf_count,
              [](const Leaf& a, const Leaf& b) {
                if (a.weight != b.weight) {
                  return a.weight < b.weight;
                }
                return a.symbol < b.symbol;
              });

    // Build the Huffman tree with the linear two-queue algorithm. Sorted
    // leaves form one queue; newly-created internal nodes form the other.
    // Internal weights are produced in nondecreasing order, so selecting the
    // smaller queue front replaces all binary-heap maintenance.
    struct TmpNode {
      std::size_t weight = 0;
      std::size_t parent = kNpos;
    };
    std::array<TmpNode, kMaxTreeNodes> tree{};
    for (std::size_t i = 0; i < leaf_count; ++i) {
      tree[i].weight = leaves[i].weight;
    }

    std::size_t leaf_head = 0;
    std::size_t internal_head = leaf_count;
    std::size_t next_internal = leaf_count;
    auto take_smallest = [&]() {
      const bool take_leaf =
          leaf_head < leaf_count &&
          (internal_head == next_internal ||
           tree[leaf_head].weight <= tree[internal_head].weight);
      return take_leaf ? leaf_head++ : internal_head++;
    };

    const std::size_t node_count = 2 * leaf_count - 1;
    while (next_internal < node_count) {
      const std::size_t left = take_smallest();
      const std::size_t right = take_smallest();
      tree[next_internal].weight = tree[left].weight + tree[right].weight;
      tree[left].parent = next_internal;
      tree[right].parent = next_internal;
      ++next_internal;
    }

    // Extract code lengths by walking from each leaf to root.
    std::uint8_t max_len = 0;
    for (std::size_t i = 0; i < leaf_count; ++i) {
      std::uint8_t depth = 0;
      std::size_t cur = i;
      while (tree[cur].parent != kNpos) {
        cur = tree[cur].parent;
        depth++;
      }
      lengths[leaves[i].symbol] = depth;
      if (depth > max_len) {
        max_len = depth;
      }
    }

    if (max_len <= kMaxCodeLength) {
      return;
    }

    // Length limiting via zlib-style overflow fixup.
    limit_code_lengths(lengths,
                       std::span<const Leaf>(leaves.data(), leaf_count));
  }

  /** @brief Limit Huffman code lengths to ≤ kMaxCodeLength.
   *  @details Collapses all lengths > L into L, rebalances the length
   *           histogram, then assigns its longest codes to the least frequent
   *           symbols. */
  void limit_code_lengths(std::array<std::uint8_t, kAlphabet>& lengths,
                          std::span<const Leaf> leaves) {
    // Count codes at each length and track how many exceeded the limit. This
    // is the same overflow definition used by zlib's gen_bitlen().
    std::array<int, kMaxCodeLength + 2> bl_count{};
    int overflow = 0;
    for (const Leaf& leaf : leaves) {
      std::uint8_t len = lengths[leaf.symbol];
      if (len > kMaxCodeLength) {
        len = kMaxCodeLength;
        ++overflow;
      }
      ++bl_count[len];
    }

    // Move one leaf down from the deepest available shorter level, creating
    // two children there, and remove one overlong code from level L.
    while (overflow > 0) {
      int bits = static_cast<int>(kMaxCodeLength) - 1;
      while (bits > 0 && bl_count[bits] == 0) {
        --bits;
      }
      --bl_count[bits];
      bl_count[bits + 1] += 2;
      --bl_count[kMaxCodeLength];
      overflow -= 2;
    }

    // Leaves arrive sorted by ascending frequency. Assign the longest codes
    // first so the fixup preserves the minimum weighted path length for its
    // corrected length histogram.
    lengths.fill(0);
    std::size_t leaf_index = 0;
    for (int len = static_cast<int>(kMaxCodeLength); len >= 1; --len) {
      for (int count = 0; count < bl_count[len]; ++count) {
        lengths[leaves[leaf_index++].symbol] = static_cast<std::uint8_t>(len);
      }
    }
  }

  // --- bit helpers ---------------------------------------------------------
  //  Small functions for reading/writing individual bits and multi-bit
  //  values from/to packed `std::uint64_t` word arrays. Used by both the
  //  normal (1-bit-per-symbol) and flat (d-bits-per-symbol) code paths.

  /** @brief Reverse the low @p width bits of an at-most-eight-bit value. */
  static std::uint8_t reverse_low_bits(std::uint8_t value, std::uint8_t width) {
    value = static_cast<std::uint8_t>(((value & 0x55u) << 1) |
                                      ((value & 0xAAu) >> 1));
    value = static_cast<std::uint8_t>(((value & 0x33u) << 2) |
                                      ((value & 0xCCu) >> 2));
    value = static_cast<std::uint8_t>((value << 4) | (value >> 4));
    return value >> (8 - width);
  }

  // --- flat subtree detection ----------------------------------------------

  /** @brief Sentinel value meaning "subtree has non-uniform leaf depths". */
  static constexpr std::uint8_t kNotUniform = 0xFF;

  /** @brief Detect maximal flat subtrees in two linear tree passes.
   *  @details The post-order pass computes uniform leaf depths without
   *           allocating per node. A top-down pass then selects only maximal
   *           eligible subtrees and appends their tables to one contiguous
   *           arena. This follows the paper's flat-subtree optimization while
   *           avoiding temporary child tables that a bottom-up marker would
   *           immediately discard when their parent is also flat. */
  void detect_flat_subtrees() const {
    flat_tables_.clear();
    if (root_ != kNpos) {
      std::array<std::uint8_t, kMaxTreeNodes> uniform_depths{};
      compute_uniform_depths(root_, uniform_depths);
      mark_flat_topdown(root_, uniform_depths);
    }
  }

  /** @brief Record each subtree's uniform leaf depth in post-order. */
  std::uint8_t compute_uniform_depths(
      std::size_t idx,
      std::array<std::uint8_t, kMaxTreeNodes>& depths) const {
    const Node& n = nodes_[idx];
    if (n.is_leaf) {
      return depths[idx] = 0;
    }

    const std::uint8_t left_depth = compute_uniform_depths(n.left, depths);
    const std::uint8_t right_depth = compute_uniform_depths(n.right, depths);

    if (left_depth == kNotUniform || right_depth == kNotUniform ||
        left_depth != right_depth) {
      return depths[idx] = kNotUniform;
    }
    return depths[idx] = static_cast<std::uint8_t>(left_depth + 1);
  }

  /** @brief Mark maximal flat subtrees and append their lookup tables. */
  void mark_flat_topdown(
      std::size_t idx,
      const std::array<std::uint8_t, kMaxTreeNodes>& depths) const {
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }

    const std::uint8_t depth = depths[idx];
    if (depth >= 2 && depth <= kMaxFlatDepth) {
      n.flat_depth = depth;
      n.flat_table_offset = flat_tables_.size();
      flat_tables_.resize(flat_tables_.size() + (std::size_t{1} << depth));
      n.flat_bitreverse = build_flat_table(
          idx, depth, 0, flat_tables_.data() + n.flat_table_offset);
      return;
    }

    mark_flat_topdown(n.left, depths);
    mark_flat_topdown(n.right, depths);
  }

  /** @brief Populate the flat lookup table by walking the subtree.
   *  @param idx             Current node.
   *  @param flat_depth      Total depth of the flat subtree.
   *  @param code            Suffix accumulated so far (MSB-first: left=0,
   * right=1).
   *  @param table           Output table of size `2^total_depth`.
   *  @return True when the completed table is the balanced bit-reversal
   *          permutation used by the specialized SIMD kernels. */
  bool build_flat_table(std::size_t idx,
                        std::uint8_t flat_depth,
                        std::uint32_t code,
                        symbol_type* table) const {
    const Node& n = nodes_[idx];
    if (n.is_leaf) {
      table[reverse_low_bits(static_cast<std::uint8_t>(code), flat_depth)] =
          n.symbol;
      return n.symbol == code;
    }
    const bool left_bitreverse =
        build_flat_table(n.left, flat_depth, code << 1, table);
    const bool right_bitreverse =
        build_flat_table(n.right, flat_depth, (code << 1) | 1, table);
    return left_bitreverse && right_bitreverse;
  }

  /** @brief Total frequency weight of a subtree (sum of leaf frequencies). */
  std::size_t subtree_weight(
      std::size_t idx,
      const std::array<std::size_t, kAlphabet>& freq) const {
    if (idx == kNpos) {
      return 0;
    }
    const Node& n = nodes_[idx];
    if (n.is_leaf) {
      return freq[n.symbol];
    }
    return subtree_weight(n.left, freq) + subtree_weight(n.right, freq);
  }

  using SymbolsByLength =
      std::array<std::array<std::uint8_t, kAlphabet>, kMaxCodeLength + 1>;
  using LengthCounts = std::array<std::size_t, kMaxCodeLength + 1>;

  /** @brief Rebuild a Huffman tree from code lengths, reshaped to maximize
   *         flat-subtree opportunities.
   *  @details Instead of building the canonical tree bit-by-bit (which scatters
   *           same-length symbols across the tree), this builder groups
   *           same-length symbols into power-of-two flat subtrees, then
   *           assembles them into a tree whose root-to-leaf path lengths still
   *           equal the original code lengths. This is the paper's "non-
   *           canonical subtree" optimization (Section 2.2.4): same average
   *           code lengths, different (better) tree shape. */
  void build_noncanonical_tree(
      const std::array<std::uint8_t, kAlphabet>& lengths) const {
    nodes_.clear();
    nodes_.reserve(2 * kAlphabet);

    SymbolsByLength by_length{};
    LengthCounts length_counts{};
    std::size_t symbol_count = 0;
    std::uint8_t only_symbol = 0;
    for (std::size_t symbol = 0; symbol < kAlphabet; ++symbol) {
      const std::uint8_t length = lengths[symbol];
      if (length == 0) {
        continue;
      }
      by_length[length][length_counts[length]++] =
          static_cast<std::uint8_t>(symbol);
      only_symbol = static_cast<std::uint8_t>(symbol);
      ++symbol_count;
    }

    // Special case: single distinct symbol → root is a leaf.
    if (symbol_count == 1) {
      nodes_.push_back(Node{});
      root_ = 0;
      nodes_[0].is_leaf = true;
      nodes_[0].symbol = only_symbol;
      return;
    }

    // Build the reshaped tree using the power-of-two grouping strategy.
    build_reshaped_tree(by_length, length_counts);
  }

  /** @brief Build a non-canonical tree that groups same-length symbols into
   *         power-of-two flat subtrees.
   *  @details Uses a recursive slot-filling approach. At each depth, we have
   *           `slots` open positions. We fill them with leaves (grouped into
   *           power-of-2 balanced subtrees for flat-node detection) and
   *           internal nodes (which recurse deeper). The split is determined
   *           by the Kraft equality: `leaves_at_depth + 2*internal_at_depth =
   *           slots`, and `internal_at_depth` creates the slots for the next
   *           depth. */
  /** @brief Build a balanced binary subtree of depth @p d with @p count
   *         leaves (count must be 2^d).
   *  @param symbols Array of symbol values (count entries).
   *  @param count Number of symbols (must be 2^d).
   *  @param d Depth of the subtree (0 = single leaf).
   *  @returns Node index of the subtree root. */
  std::size_t build_balanced_subtree(const std::uint8_t* symbols,
                                     std::size_t count,
                                     std::uint8_t d) const {
    if (d == 0) {
      nodes_.push_back(Node{});
      std::size_t idx = nodes_.size() - 1;
      nodes_[idx].is_leaf = true;
      nodes_[idx].symbol = symbols[0];
      return idx;
    }
    std::size_t half = count / 2;
    std::size_t left = build_balanced_subtree(symbols, half, d - 1);
    std::size_t right = build_balanced_subtree(symbols + half, half, d - 1);
    nodes_.push_back(Node{});
    std::size_t idx = nodes_.size() - 1;
    nodes_[idx].left = left;
    nodes_[idx].right = right;
    return idx;
  }

  void build_reshaped_tree(const SymbolsByLength& by_length,
                           const LengthCounts& length_counts) const {
    // Per-depth consumption cursor.
    LengthCounts cursor{};
    root_ = build_slots(by_length, length_counts, cursor, 1, 2);
  }

  /** @brief Fill `slots` open positions at `depth` with leaves and internal
   *         nodes, returning the root of the resulting balanced subtree.
   *  @param by_length  Symbols grouped by code length.
   *  @param cursor     Per-depth consumption cursor (how many symbols used).
   *  @param depth      Current tree depth (1 = children of root).
   *  @param slots      Number of open positions (must be a power of 2). */
  std::size_t build_slots(const SymbolsByLength& by_length,
                          const LengthCounts& length_counts,
                          LengthCounts& cursor,
                          std::uint8_t depth,
                          std::size_t slots) const {
    if (slots == 0) {
      return kNpos;
    }

    // How many leaves to place at this depth.
    std::size_t available = length_counts[depth] - cursor[depth];
    std::size_t leaves_here = std::min(available, slots);
    std::size_t internal = slots - leaves_here;

    // Build children left-to-right: first leaves (grouped), then internal.
    std::array<std::size_t, 2> children{};
    std::size_t child_count = 0;

    // Group leaves into power-of-2 balanced subtrees.
    std::size_t leaf_off = cursor[depth];
    std::size_t remaining = leaves_here;
    while (remaining > 0) {
      std::uint8_t d = 0;
      std::size_t subset = 1;
      while (subset * 2 <= remaining && d < kMaxFlatDepth) {
        subset *= 2;
        d++;
      }
      children[child_count++] =
          build_balanced_subtree(by_length[depth].data() + leaf_off, subset, d);
      leaf_off += subset;
      remaining -= subset;
    }
    cursor[depth] += leaves_here;

    // Internal nodes: each recurses with 2 slots at depth+1.
    for (std::size_t i = 0; i < internal; i++) {
      children[child_count++] =
          build_slots(by_length, length_counts, cursor, depth + 1, 2);
    }

    // Assemble children into a balanced binary tree.
    return assemble_balanced(
        std::span<const std::size_t>(children.data(), child_count));
  }

  /** @brief Assemble a list of subtree roots into a balanced binary tree.
   *  @details Recursively pairs children left-right, creating internal nodes.
   */
  std::size_t assemble_balanced(std::span<const std::size_t> children) const {
    if (children.empty()) {
      return kNpos;
    }
    if (children.size() == 1) {
      return children[0];
    }
    // Pair up: left half vs right half.
    std::size_t mid = children.size() / 2;
    std::size_t left = assemble_balanced(children.first(mid));
    std::size_t right = assemble_balanced(children.subspan(mid));
    nodes_.push_back(Node{});
    std::size_t idx = nodes_.size() - 1;
    nodes_[idx].left = left;
    nodes_[idx].right = right;
    return idx;
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

    if (n.flat_depth > 0) {
      // Flat node: compute weight from subtree, assign arena for its
      // d-bits-per-symbol bitmap, but DON'T recurse into children (their
      // bitmaps aren't stored — the flat node replaces them).
      std::size_t weight = subtree_weight(idx, freq);
      n.bits.count = weight;
      n.bits.offset = arena_words;
      arena_words += n.bits.word_count(n.flat_depth);
      return weight;
    }

    std::size_t left_weight =
        assign_weights_and_offsets(n.left, freq, arena_words);
    std::size_t right_weight =
        assign_weights_and_offsets(n.right, freq, arena_words);
    n.bits.count = left_weight + right_weight;
    n.bits.offset = arena_words;
    arena_words += n.bits.word_count(0);
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

  /** @brief Assign compact MSB-first codes with a depth-first tree walk. */
  void assign_codes(std::size_t idx,
                    std::uint16_t bits,
                    std::uint8_t length,
                    std::array<SymbolCode, kAlphabet>& codes) const {
    if (nodes_[idx].is_leaf) {
      codes[nodes_[idx].symbol] = {bits, length};
      return;
    }
    assign_codes(nodes_[idx].left, static_cast<std::uint16_t>(bits << 1),
                 length + 1, codes);
    assign_codes(nodes_[idx].right, static_cast<std::uint16_t>((bits << 1) | 1),
                 length + 1, codes);
  }

  /** @brief Build the per-depth direction tables shared by all tree nodes.
   *  @details A symbol can occur in only one node at a given depth, so one
   *           256-byte table per depth replaces rebuilding the same direction
   *           information separately for every internal node. */
  static DirectionTables build_direction_tables(
      const std::array<SymbolCode, kAlphabet>& codes) {
    DirectionTables directions{};
    for (std::size_t symbol = 0; symbol < kAlphabet; ++symbol) {
      const SymbolCode code = codes[symbol];
      for (std::uint8_t depth = 0; depth < code.length; ++depth) {
        directions[depth][symbol] = static_cast<std::uint8_t>(
            (code.bits >> (code.length - depth - 1)) & 1u);
      }
    }
    return directions;
  }

  /** @brief Pack one flat node's fixed-width wire codes.
   *  @details Inverts the node's small decode table once, then packs directly
   *           from the symbol-to-wire-code lookup. Widths dividing a byte get
   *           dedicated loops; other widths track word position without a
   *           division per symbol. */
  static void encode_flat_values(std::uint64_t* bits,
                                 const symbol_type* src,
                                 std::size_t count,
                                 std::uint8_t flat_depth,
                                 const symbol_type* table,
                                 bool known_bitreverse) {
#if defined(PIXIE_PIVCO_NEON_SUPPORT)
    if (known_bitreverse) {
      pivco_flat_encode_neon(bits, src, nullptr, count, flat_depth, true);
      return;
    }
#endif

    // Three zero pad bytes make overlapping dword gathers at symbol 255 safe
    // in the AVX2 flat encoder.
    std::array<std::uint8_t, kAlphabet + 3> wire_code{};
    const std::size_t alphabet_size = std::size_t{1} << flat_depth;
    for (std::size_t code = 0; code < alphabet_size; ++code) {
      wire_code[table[code]] = static_cast<std::uint8_t>(code);
    }

#if defined(__AVX2__) && defined(__SSE4_1__) && !defined(__AVX512VBMI2__)
    if (flat_depth == 8 &&
        pivco_flat_encode_d8_bitreverse(bits, src, wire_code.data(), count,
                                        known_bitreverse)) {
      return;
    }
#if defined(__BMI2__)
    if (flat_depth >= 2 && flat_depth <= 7 &&
        pivco_flat_encode_bitreverse(bits, src, wire_code.data(), count,
                                     flat_depth, known_bitreverse)) {
      return;
    }
    if (flat_depth == 3) {
      pivco_flat_encode_d3(bits, src, wire_code.data(), count);
      return;
    }
#endif
#endif

#if defined(PIXIE_PIVCO_NEON_SUPPORT)
    pivco_flat_encode_neon(bits, src, wire_code.data(), count, flat_depth,
                           false);
    return;
#endif

    auto* bytes = reinterpret_cast<std::uint8_t*>(bits);
    if (flat_depth == 8) {
      for (std::size_t i = 0; i < count; ++i) {
        bytes[i] = wire_code[src[i]];
      }
      return;
    }
    if (flat_depth == 4) {
      std::size_t i = 0;
      for (; i + 2 <= count; i += 2) {
        bytes[i / 2] = static_cast<std::uint8_t>(wire_code[src[i]] |
                                                 (wire_code[src[i + 1]] << 4));
      }
      if (i < count) {
        bytes[i / 2] = wire_code[src[i]];
      }
      return;
    }
    if (flat_depth == 2) {
      std::size_t i = 0;
      for (; i + 4 <= count; i += 4) {
        bytes[i / 4] = static_cast<std::uint8_t>(
            wire_code[src[i]] | (wire_code[src[i + 1]] << 2) |
            (wire_code[src[i + 2]] << 4) | (wire_code[src[i + 3]] << 6));
      }
      if (i < count) {
        std::uint8_t packed = 0;
        for (std::size_t lane = 0; i + lane < count; ++lane) {
          packed |=
              static_cast<std::uint8_t>(wire_code[src[i + lane]] << (2 * lane));
        }
        bytes[i / 4] = packed;
      }
      return;
    }

    std::size_t word_index = 0;
    std::size_t bit_offset = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const std::uint64_t packed = wire_code[src[i]];
      bits[word_index] |= packed << bit_offset;
      if (bit_offset + flat_depth > kWordBits) {
        bits[word_index + 1] |= packed >> (kWordBits - bit_offset);
      }
      bit_offset += flat_depth;
      if (bit_offset >= kWordBits) {
        bit_offset -= kWordBits;
        ++word_index;
      }
    }
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
   *  @param directions Per-depth symbol direction lookup tables.
   *  @param freq     Symbol frequencies (sizing child partitions). */
  void encode_partition(std::size_t idx,
                        std::size_t depth,
                        std::size_t out_base,
                        const DirectionTables& directions,
                        const std::array<std::size_t, kAlphabet>& freq) {
    if (nodes_[idx].is_leaf) {
      return;
    }

    Node& n = nodes_[idx];
    if (n.bits.count == 0) {
      return;
    }

    // Flat node: write d-bit suffix per symbol, no partition, no recursion.
    if (n.flat_depth > 0) {
      std::uint64_t* bits = arena_.data() + n.bits.offset;
      const symbol_type* src = workspace_half(depth) + out_base;
      if (n.flat_depth == 1) {
        pivco_partition_encode_none(bits, src, directions[depth].data(),
                                    n.bits.count);
      } else {
        const symbol_type* table = flat_tables_.data() + n.flat_table_offset;
        encode_flat_values(bits, src, n.bits.count, n.flat_depth, table,
                           n.flat_bitreverse);
      }
      return;
    }

    // Normal node: 1-bit-per-symbol bitmap + partition into children.
    const std::size_t left_weight = child_weight(n.left, freq);
    const std::size_t weight = n.bits.count;

    const symbol_type* src = workspace_half(depth) + out_base;
    symbol_type* left_dst = workspace_half(depth + 1) + out_base;
    symbol_type* right_dst = workspace_half(depth + 1) + out_base + left_weight;

    std::uint64_t* bits = arena_.data() + n.bits.offset;
    const std::uint8_t* dir = directions[depth].data();
    const bool left_is_leaf = nodes_[n.left].is_leaf;
    const bool right_is_leaf = nodes_[n.right].is_leaf;
    if (left_is_leaf && right_is_leaf) {
      pivco_partition_encode_cst_cst(bits, src, nodes_[n.right].symbol, weight);
    } else if (right_is_leaf) {
#if ((defined(__AVX2__) && defined(__SSE4_1__)) || \
     defined(PIXIE_PIVCO_NEON_SUPPORT)) &&         \
    !defined(__AVX512VBMI2__)
      pivco_partition_encode_vec_cst(
          left_dst, bits, src, nodes_[n.right].symbol, weight, left_weight);
#else
      pivco_partition_encode_left(left_dst, bits, src, dir, weight,
                                  left_weight);
#endif
    } else if (left_is_leaf) {
#if ((defined(__AVX2__) && defined(__SSE4_1__)) || \
     defined(PIXIE_PIVCO_NEON_SUPPORT)) &&         \
    !defined(__AVX512VBMI2__)
      pivco_partition_encode_cst_vec(right_dst, bits, src,
                                     nodes_[n.left].symbol, weight,
                                     weight - left_weight);
#else
      pivco_partition_encode_right(right_dst, bits, src, dir, weight,
                                   weight - left_weight);
#endif
    } else {
      pivco_partition_encode(left_dst, right_dst, bits, src, dir, weight,
                             left_weight, weight - left_weight);
    }

    encode_partition(n.left, depth + 1, out_base, directions, freq);
    encode_partition(n.right, depth + 1, out_base + left_weight, directions,
                     freq);
  }

  // --- decode --------------------------------------------------------------

  /** @brief Decode the current block directly into its final output slice.
   *  @details Internal levels still ping-pong through the reusable two-block
   *           workspace, but the root writes into @p output. This avoids a
   *           temporary block vector and the subsequent append/copy. */
  void decode_block_into(std::span<symbol_type> output) const {
    if (root_ == kNpos || output.empty()) {
      return;
    }
    if (nodes_[root_].is_leaf) {
      std::fill(output.begin(), output.end(), nodes_[root_].symbol);
      return;
    }
    workspace_.resize(output.size() * 2 + kWorkspaceSlack);
    decode_node(root_, output.size(), 0, 0, output.data());
  }

  /** @brief Recursively decode a node's output stream bottom-up into the
   *         reusable workspace.
   *  @param weight   Number of symbols this node must produce.
   *  @param depth     Tree depth. The root writes directly to the caller's
   *                   output; deeper levels select a workspace half.
   *  @param out_base  Offset within the selected half where this node's output
   *                   begins. Children are placed at `out_base` (left) and
   *                   `out_base + left_weight` (right) in the opposite half,
   *                   then merged into `[out_base, out_base + weight)` here.
   *  @param block_output Final output address used by the root node. */
  void decode_node(std::size_t idx,
                   std::size_t weight,
                   std::size_t depth,
                   std::size_t out_base,
                   symbol_type* block_output) const {
    if (weight == 0) {
      return;
    }
    const Node& n = nodes_[idx];
    symbol_type* dst =
        depth == 0 ? block_output : workspace_half(depth) + out_base;
    if (n.is_leaf) {
      std::fill(dst, dst + weight, n.symbol);
      return;
    }

    // Flat node: read d-bit suffixes, lookup symbols, no merge.
    if (n.flat_depth > 0) {
      const std::uint64_t* bits = arena_.data() + n.bits.offset;
      const symbol_type* table = flat_tables_.data() + n.flat_table_offset;
      pivco_flat_decode(dst, bits, table, weight, n.flat_depth,
                        n.flat_bitreverse);
      return;
    }

    // Normal node: decode non-leaf children, then select the paper's matching
    // MVV/MCV/MVC/MCC merge primitive. Internal children already carry their
    // exact stream size, so no bitmap popcount pass is needed.
    const Node& left = nodes_[n.left];
    const Node& right = nodes_[n.right];
    const std::uint64_t* bits = arena_.data() + n.bits.offset;
    if (left.is_leaf && right.is_leaf) {
      pivco_merge_decode_cst_cst(dst, left.symbol, right.symbol, bits, weight);
      return;
    }

    const std::size_t left_weight =
        left.is_leaf ? weight - right.bits.count : left.bits.count;
    const std::size_t right_weight = weight - left_weight;
    if (!left.is_leaf) {
      decode_node(n.left, left_weight, depth + 1, out_base, block_output);
    }
    if (!right.is_leaf) {
      decode_node(n.right, right_weight, depth + 1, out_base + left_weight,
                  block_output);
    }

    const symbol_type* src = workspace_half(depth + 1) + out_base;
    const symbol_type* left_out = src;
    const symbol_type* right_out = src + left_weight;
    if (left.is_leaf) {
      pivco_merge_decode_cst_vec(dst, left.symbol, right_out, bits, weight);
    } else if (right.is_leaf) {
      pivco_merge_decode_vec_cst(dst, left_out, right.symbol, bits, weight);
    } else {
      pivco_merge_decode(dst, left_out, right_out, bits, weight);
    }
  }

  /** @brief Pointer to the workspace half that holds depth-@p depth outputs. */
  symbol_type* workspace_half(std::size_t depth) const {
    return workspace_.data() + (depth % 2) * block_uncompressed_size_;
  }

  // --- serialization -------------------------------------------------------

  /** @brief Append the shared file-level code lengths as packed nibbles. */
  void serialize_code_lengths() {
    for (std::size_t i = 0; i < kAlphabet; i += 2) {
      std::uint8_t byte = static_cast<std::uint8_t>(
          code_lengths_[i] | (code_lengths_[i + 1] << 4));
      write(byte);
    }
  }

  /** @brief Append one block's serialized payload to `compressed_`.
   *  @details Per-block wire format:
   *           - block_uncompressed_size (8 bytes)
   *           - For each internal node (pre-order): bits.count (8 bytes) +
   *             packed words. Leaves and topology are implied by the shared
   *             file-level code lengths. */
  void serialize_block_payload() {
    write(block_uncompressed_size_);
    if (block_uncompressed_size_ == 0) {
      return;
    }

    // Write internal node bitmaps in pre-order (matches deserialize order).
    serialize_node(root_);
  }

  /** @brief Serialize a node's bitmap and recurse (pre-order).
   *  @details Flat nodes are terminal: their children's bitmaps are not
   *           serialized (the flat node replaces them). */
  void serialize_node(std::size_t idx) {
    if (idx == kNpos) {
      return;
    }
    const auto& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }
    write(n.bits.count);
    write_words(n.bits.offset, n.bits.word_count(n.flat_depth));
    if (n.flat_depth > 0) {
      return;  // flat: children not serialized
    }
    serialize_node(n.left);
    serialize_node(n.right);
  }

  /** @brief Load the serialized stream and reconstruct its shared tree.
   *  @details Copies @p data into `compressed_` and reads the top-level
   *           header and file-level code lengths. Per-block bitmap payloads
   *           are loaded on demand during `decode_impl()`. */
  void deserialize(std::span<const std::byte> data) {
    compressed_.assign(data.begin(), data.end());
    nodes_.clear();
    arena_.clear();
    root_ = kNpos;
    if (data.empty()) {
      uncompressed_size_ = 0;
      blocks_offset_ = 0;
      return;
    }
    std::size_t pos = 0;
    uncompressed_size_ = read<std::size_t>(compressed_, pos);
    block_size_ = read<std::size_t>(compressed_, pos);
    if (uncompressed_size_ == 0) {
      blocks_offset_ = pos;
      return;
    }

    for (std::size_t i = 0; i < kAlphabet; i += 2) {
      const std::uint8_t byte = read<std::uint8_t>(compressed_, pos);
      code_lengths_[i] = byte & 0x0F;
      code_lengths_[i + 1] = (byte >> 4) & 0x0F;
    }
    build_noncanonical_tree(code_lengths_);
    detect_flat_subtrees();
    blocks_offset_ = pos;
  }

  /** @brief Load one block's arena from `compressed_` at @p pos.
   *  @details Reuses the file-level tree, reads the block's uncompressed size,
   *           and loads its arena words via a two-pass wire walk. Advances
   *           @p pos past the block payload. */
  void deserialize_block_payload(std::size_t& pos) const {
    block_uncompressed_size_ = read<std::size_t>(compressed_, pos);
    if (block_uncompressed_size_ == 0) {
      return;
    }

    // Single-symbol input: root is a leaf, no internal nodes, no arena.
    if (root_ == kNpos || nodes_[root_].is_leaf) {
      return;
    }

    // Pass 1: read each internal node's `count` (and skip its word bytes on
    // the wire), assign contiguous arena offsets, and size the arena exactly.
    const std::size_t bitmap_start = pos;
    std::size_t arena_words = 0;
    scan_node_counts(root_, compressed_, pos, arena_words);
    arena_.assign(arena_words, 0);

    // Pass 2: rewind to the bitmap section and read each node's packed words
    // into its arena slot.
    pos = bitmap_start;
    read_node_words(root_, compressed_, pos);
  }

  /** @brief Pass 1: read each internal node's `count` from the wire, assign its
   *         arena offset, and advance past its word bytes without storing
   *         them. Sums the total arena word count so the caller allocates once.
   *  @details Pre-order traversal matches the serialization order. Flat nodes
   *           are terminal (children not visited). */
  void scan_node_counts(std::size_t idx,
                        std::span<const std::byte> data,
                        std::size_t& pos,
                        std::size_t& arena_words) const {
    if (idx == kNpos) {
      return;
    }
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }
    n.bits.count = read<std::size_t>(data, pos);
    n.bits.offset = arena_words;
    arena_words += n.bits.word_count(n.flat_depth);
    pos += n.bits.word_count(n.flat_depth) * sizeof(std::uint64_t);
    if (n.flat_depth > 0) {
      return;  // flat: children not on the wire
    }
    scan_node_counts(n.left, data, pos, arena_words);
    scan_node_counts(n.right, data, pos, arena_words);
  }

  /** @brief Pass 2: read each internal node's packed words from the wire into
   *         its arena slot (pre-order, matching serialization). Flat nodes are
   *         terminal. */
  void read_node_words(std::size_t idx,
                       std::span<const std::byte> data,
                       std::size_t& pos) const {
    if (idx == kNpos) {
      return;
    }
    Node& n = nodes_[idx];
    if (n.is_leaf) {
      return;
    }
    pos += sizeof(std::size_t);  // count, already read in pass 1
    const std::size_t words = n.bits.word_count(n.flat_depth);
    const std::size_t bytes = words * sizeof(std::uint64_t);
    if (bytes > 0) {
      std::memcpy(arena_.data() + n.bits.offset, data.data() + pos, bytes);
    }
    pos += bytes;
    if (n.flat_depth > 0) {
      return;  // flat: children not on the wire
    }
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
