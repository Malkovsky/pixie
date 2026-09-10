#pragma once

/**
 * @file huffman_build_table.h
 * @brief Length-limited byte Huffman tree construction for wavelet indexes.
 *
 * This file is derived from PivCo's `huffman_table.c` and has been modified
 * for Pixie's header-only C++ representation. Wire-format, flat-subtree,
 * entropy-coding, and decoder tables were removed. The two-queue length
 * builder, length limiter, fused canonical tree shaping, and in-order rank
 * assignment are retained. Both projects are distributed under Apache-2.0.
 */

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace pixie::detail {

inline constexpr std::size_t kByteAlphabetSize = 256;
inline constexpr std::size_t kMaximumHuffmanNodes = 2 * kByteAlphabetSize - 1;
inline constexpr std::size_t kMaximumHuffmanCodeLength = 11;

struct HuffmanBuildNode {
  std::int16_t symbol = -1;
  std::int16_t left = -1;
  std::int16_t right = -1;
};

enum class HuffmanBuildNodeType : std::uint8_t {
  kInternal,
  kBothLeaves,
  kLeftLeaf,
  kLeaf,
};

struct HuffmanBuildTable {
  std::array<HuffmanBuildNode, kMaximumHuffmanNodes> tree{};
  std::array<std::uint8_t, kMaximumHuffmanNodes> split_rank{};
  std::array<HuffmanBuildNodeType, kMaximumHuffmanNodes> node_type{};
  std::array<std::uint8_t, kByteAlphabetSize> symbol_to_rank{};
  std::int16_t root = -1;
  std::size_t node_count = 0;
  std::size_t symbol_count = 0;
};

namespace huffman_build_detail {

struct FrequencyLeaf {
  std::size_t frequency;
  std::uint16_t symbol;
};

// Stable LSD radix sort. The input is seeded in symbol order, so equal
// frequencies retain PivCo's (frequency, symbol) tie discipline.
inline void sort_leaves_by_frequency(
    std::span<FrequencyLeaf> leaves,
    std::array<FrequencyLeaf, kByteAlphabetSize>& temporary) {
  std::size_t maximum = 0;
  for (const FrequencyLeaf leaf : leaves) {
    maximum = std::max(maximum, leaf.frequency);
  }
  std::size_t byte_count = 0;
  while (maximum != 0) {
    ++byte_count;
    maximum >>= 8;
  }

  FrequencyLeaf* source = leaves.data();
  FrequencyLeaf* destination = temporary.data();
  for (std::size_t byte = 0; byte < byte_count; ++byte) {
    const std::size_t shift = byte * 8;
    std::array<std::size_t, 256> counts{};
    for (std::size_t index = 0; index < leaves.size(); ++index) {
      ++counts[(source[index].frequency >> shift) & 0xff];
    }
    std::size_t prefix = 0;
    for (std::size_t& count : counts) {
      const std::size_t current = count;
      count = prefix;
      prefix += current;
    }
    for (std::size_t index = 0; index < leaves.size(); ++index) {
      const std::size_t bucket = (source[index].frequency >> shift) & 0xff;
      destination[counts[bucket]++] = source[index];
    }
    std::swap(source, destination);
  }
  if (source != leaves.data()) {
    std::copy_n(source, leaves.size(), leaves.data());
  }
}

inline std::array<std::uint8_t, kByteAlphabetSize> build_code_lengths(
    std::span<const std::size_t> frequencies,
    std::span<const std::uint16_t> used_symbols) {
  std::array<std::uint8_t, kByteAlphabetSize> lengths{};
  if (used_symbols.size() < 2) {
    if (!used_symbols.empty()) {
      lengths[used_symbols.front()] = 1;
    }
    return lengths;
  }

  std::array<FrequencyLeaf, kByteAlphabetSize> leaves{};
  for (std::size_t index = 0; index < used_symbols.size(); ++index) {
    const std::uint16_t symbol = used_symbols[index];
    leaves[index] = {frequencies[symbol], symbol};
  }
  std::array<FrequencyLeaf, kByteAlphabetSize> sort_temporary{};
  sort_leaves_by_frequency(std::span(leaves).first(used_symbols.size()),
                           sort_temporary);

  std::array<std::size_t, kMaximumHuffmanNodes> node_frequency{};
  std::array<std::uint16_t, kMaximumHuffmanNodes> parent{};
  const std::size_t leaf_count = used_symbols.size();
  for (std::size_t index = 0; index < leaf_count; ++index) {
    node_frequency[index] = leaves[index].frequency;
  }

  std::size_t next_leaf = 0;
  std::size_t internal_head = leaf_count;
  std::size_t next_internal = leaf_count;
  const auto take_minimum = [&]() {
    if (next_leaf < leaf_count &&
        (internal_head == next_internal ||
         node_frequency[next_leaf] <= node_frequency[internal_head])) {
      return next_leaf++;
    }
    return internal_head++;
  };
  for (std::size_t remaining = leaf_count; remaining > 1; --remaining) {
    const std::size_t left = take_minimum();
    const std::size_t right = take_minimum();
    node_frequency[next_internal] =
        node_frequency[left] + node_frequency[right];
    parent[left] = static_cast<std::uint16_t>(next_internal);
    parent[right] = static_cast<std::uint16_t>(next_internal);
    ++next_internal;
  }

  const std::size_t root = next_internal - 1;
  std::array<std::uint8_t, kMaximumHuffmanNodes> depth{};
  for (std::size_t index = root; index-- > 0;) {
    depth[index] = static_cast<std::uint8_t>(depth[parent[index]] + 1);
  }
  for (std::size_t index = 0; index < leaf_count; ++index) {
    lengths[leaves[index].symbol] = depth[index];
  }
  return lengths;
}

// PivCo's DEFLATE-style limiter caps tree height so byte-rank construction has
// a small, predictable recursion and scratch bound.
inline void limit_code_lengths(
    std::array<std::uint8_t, kByteAlphabetSize>& lengths) {
  std::array<std::size_t, kByteAlphabetSize> length_counts{};
  std::size_t maximum = 0;
  for (const std::uint8_t length : lengths) {
    if (length != 0) {
      ++length_counts[length];
      maximum = std::max(maximum, static_cast<std::size_t>(length));
    }
  }
  if (maximum <= kMaximumHuffmanCodeLength) {
    return;
  }

  for (std::size_t length = maximum; length > kMaximumHuffmanCodeLength;
       --length) {
    length_counts[kMaximumHuffmanCodeLength] += length_counts[length];
    length_counts[length] = 0;
  }

  std::size_t kraft = 0;
  for (std::size_t length = 1; length <= kMaximumHuffmanCodeLength; ++length) {
    kraft += length_counts[length] << (kMaximumHuffmanCodeLength - length);
  }
  constexpr std::size_t kKraftTarget = std::size_t{1}
                                       << kMaximumHuffmanCodeLength;
  while (kraft > kKraftTarget) {
    std::size_t best = kMaximumHuffmanCodeLength - 1;
    while (best > 0 && length_counts[best] == 0) {
      --best;
    }
    if (best == 0) {
      break;
    }
    --length_counts[best];
    ++length_counts[best + 1];
    kraft -= std::size_t{1} << (kMaximumHuffmanCodeLength - best - 1);
  }
  while (kraft < kKraftTarget &&
         length_counts[kMaximumHuffmanCodeLength] != 0) {
    bool shortened = false;
    for (std::size_t length = kMaximumHuffmanCodeLength - 1; length > 0;
         --length) {
      const std::size_t delta =
          (std::size_t{1} << (kMaximumHuffmanCodeLength - length)) - 1;
      if (kraft + delta <= kKraftTarget) {
        --length_counts[kMaximumHuffmanCodeLength];
        ++length_counts[length];
        kraft += delta;
        shortened = true;
        break;
      }
    }
    if (!shortened) {
      break;
    }
  }

  struct LengthSymbol {
    std::uint8_t length;
    std::uint8_t symbol;
  };
  std::array<LengthSymbol, kByteAlphabetSize> ordered{};
  std::size_t count = 0;
  for (std::size_t symbol = 0; symbol < lengths.size(); ++symbol) {
    if (lengths[symbol] != 0) {
      ordered[count++] = {
          std::min<std::uint8_t>(lengths[symbol], kMaximumHuffmanCodeLength),
          static_cast<std::uint8_t>(symbol)};
    }
  }
  for (std::size_t index = 1; index < count; ++index) {
    const LengthSymbol current = ordered[index];
    std::size_t position = index;
    while (position != 0 && ordered[position - 1].length > current.length) {
      ordered[position] = ordered[position - 1];
      --position;
    }
    ordered[position] = current;
  }
  std::size_t ordered_index = 0;
  for (std::size_t length = 1; length <= kMaximumHuffmanCodeLength; ++length) {
    for (std::size_t index = 0; index < length_counts[length]; ++index) {
      lengths[ordered[ordered_index++].symbol] =
          static_cast<std::uint8_t>(length);
    }
  }
}

struct FusedChunk {
  std::uint16_t suffix_bits;
  std::uint16_t depth;
  std::uint16_t symbol_count;
  std::uint16_t root_code;
  std::size_t symbol_begin;
};

inline std::uint16_t assign_inorder_ranks(HuffmanBuildTable& table,
                                          std::int16_t node_id,
                                          std::uint16_t rank) {
  const HuffmanBuildNode& node = table.tree[node_id];
  if (node.symbol >= 0) {
    table.symbol_to_rank[node.symbol] = static_cast<std::uint8_t>(rank);
    return static_cast<std::uint16_t>(rank + 1);
  }
  rank = assign_inorder_ranks(table, node.left, rank);
  table.split_rank[node_id] = static_cast<std::uint8_t>(rank - 1);
  return assign_inorder_ranks(table, node.right, rank);
}

inline void build_fused_tree(
    const std::array<std::uint8_t, kByteAlphabetSize>& lengths,
    HuffmanBuildTable& table) {
  std::array<std::uint16_t, kMaximumHuffmanCodeLength + 1> length_counts{};
  std::size_t maximum_length = 0;
  for (const std::uint8_t length : lengths) {
    if (length != 0) {
      ++length_counts[length];
      maximum_length =
          std::max(maximum_length, static_cast<std::size_t>(length));
    }
  }

  std::array<std::uint8_t, kByteAlphabetSize> ordered_symbols{};
  std::array<std::size_t, kMaximumHuffmanCodeLength + 2> length_begin{};
  std::array<std::size_t, kMaximumHuffmanCodeLength + 2> cursor{};
  std::size_t accumulated = 0;
  for (std::size_t length = 1; length <= maximum_length; ++length) {
    length_begin[length] = accumulated;
    cursor[length] = accumulated;
    accumulated += length_counts[length];
  }
  length_begin[maximum_length + 1] = accumulated;
  for (std::size_t symbol = 0; symbol < lengths.size(); ++symbol) {
    const std::uint8_t length = lengths[symbol];
    if (length != 0) {
      ordered_symbols[cursor[length]++] = static_cast<std::uint8_t>(symbol);
    }
  }

  std::array<FusedChunk, kByteAlphabetSize> chunks{};
  std::size_t chunk_count = 0;
  for (std::size_t length = 1; length <= maximum_length; ++length) {
    std::size_t symbol = length_begin[length];
    for (std::size_t pair = 0; pair < length_counts[length] / 2; ++pair) {
      chunks[chunk_count++] = {1, static_cast<std::uint16_t>(length - 1), 2, 0,
                               symbol};
      symbol += 2;
    }
    if ((length_counts[length] & 1U) != 0) {
      chunks[chunk_count++] = {0, static_cast<std::uint16_t>(length), 1, 0,
                               symbol};
    }
  }
  for (std::size_t index = 1; index < chunk_count; ++index) {
    const FusedChunk current = chunks[index];
    std::size_t position = index;
    while (position != 0 && chunks[position - 1].depth > current.depth) {
      chunks[position] = chunks[position - 1];
      --position;
    }
    chunks[position] = current;
  }

  std::uint32_t code = 0;
  std::size_t previous_depth = 0;
  for (std::size_t index = 0; index < chunk_count; ++index) {
    FusedChunk& chunk = chunks[index];
    code <<= chunk.depth - previous_depth;
    chunk.root_code = static_cast<std::uint16_t>(code);
    ++code;
    previous_depth = chunk.depth;
  }

  table.root = 0;
  table.node_count = 1;
  for (std::size_t index = 0; index < chunk_count; ++index) {
    const FusedChunk& chunk = chunks[index];
    std::int16_t node_id = table.root;
    for (std::size_t bit = chunk.depth; bit-- > 0;) {
      const bool right = ((chunk.root_code >> bit) & 1U) != 0;
      std::int16_t& child =
          right ? table.tree[node_id].right : table.tree[node_id].left;
      if (child < 0) {
        child = static_cast<std::int16_t>(table.node_count++);
      }
      node_id = child;
    }
    if (chunk.suffix_bits == 1) {
      HuffmanBuildNode& node = table.tree[node_id];
      node.left = static_cast<std::int16_t>(table.node_count++);
      table.tree[node.left].symbol = ordered_symbols[chunk.symbol_begin];
      node.right = static_cast<std::int16_t>(table.node_count++);
      table.tree[node.right].symbol = ordered_symbols[chunk.symbol_begin + 1];
    } else {
      table.tree[node_id].symbol = ordered_symbols[chunk.symbol_begin];
    }
  }

  assign_inorder_ranks(table, table.root, 0);
  for (std::size_t index = 0; index < table.node_count; ++index) {
    const HuffmanBuildNode& node = table.tree[index];
    if (node.symbol >= 0) {
      table.node_type[index] = HuffmanBuildNodeType::kLeaf;
      continue;
    }
    const bool left_leaf = table.tree[node.left].symbol >= 0;
    const bool right_leaf = table.tree[node.right].symbol >= 0;
    if (left_leaf && right_leaf) {
      table.node_type[index] = HuffmanBuildNodeType::kBothLeaves;
    } else if (left_leaf) {
      table.node_type[index] = HuffmanBuildNodeType::kLeftLeaf;
    } else {
      table.node_type[index] = HuffmanBuildNodeType::kInternal;
    }
  }
}

}  // namespace huffman_build_detail

/**
 * @brief Build PivCo's fused, non-flat byte Huffman tree.
 * @param frequencies Frequencies for a dense alphabet of at most 256 symbols.
 * @return Tree, node dispatch classes, symbol ranks, and split ranks.
 */
inline HuffmanBuildTable build_huffman_table(
    std::span<const std::size_t> frequencies) {
  HuffmanBuildTable table;
  std::array<std::uint16_t, kByteAlphabetSize> used_symbols{};
  for (std::size_t symbol = 0; symbol < frequencies.size(); ++symbol) {
    if (frequencies[symbol] != 0) {
      used_symbols[table.symbol_count++] = static_cast<std::uint16_t>(symbol);
    }
  }
  if (table.symbol_count == 0) {
    return table;
  }
  if (table.symbol_count == 1) {
    table.root = 0;
    table.node_count = 1;
    table.tree[0].symbol = used_symbols[0];
    table.node_type[0] = HuffmanBuildNodeType::kLeaf;
    table.symbol_to_rank[used_symbols[0]] = 0;
    return table;
  }

  auto lengths = huffman_build_detail::build_code_lengths(
      frequencies, std::span(used_symbols).first(table.symbol_count));
  huffman_build_detail::limit_code_lengths(lengths);
  huffman_build_detail::build_fused_tree(lengths, table);
  return table;
}

}  // namespace pixie::detail
