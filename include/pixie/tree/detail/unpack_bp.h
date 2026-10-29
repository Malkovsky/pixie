#pragma once
#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <span>

namespace pixie::detail::tree_unpack {
template <size_t Size, class Tree, class Visitor>
void bp_deep_word(const Tree& tree,
                  std::span<const uint64_t> words,
                  size_t base,
                  std::array<uint64_t, Size>& numbers,
                  std::array<uint64_t, Size>& positions,
                  size_t& depth,
                  size_t& number,
                  TreeNode& current,
                  Visitor& visit) {
  uint64_t word = word_at(words, base);
  bool lookahead = word_at(words, base + 64) & 1;
  unsigned offset = 0;
  while (offset < 64) {
    bool opening = word & 1;
    unsigned count = std::min<unsigned>(
        std::countr_zero(opening ? ~word : word), 64 - offset);
    const auto start = base + offset;
    offset += count;
    word = count == 64 ? 0 : word >> count;
    if (opening) {
      for (unsigned i = 0; i < count; ++i) {
        current = TreeNode(number++, start + i);
        if (depth < 128) {
          numbers[depth + 8] = current.number;
          positions[depth + 8] = current.pos;
        }
        ++depth;
        visit(current);
      }
    } else {
      for (unsigned i = 0; i < count; ++i) {
        --depth;
        if (depth &&
            (i + 1 < count || !(offset < 64 ? (word & 1) : lookahead))) {
          if (depth <= 128) {
            current = TreeNode(numbers[depth + 7], positions[depth + 7]);
          } else {
            current = tree.parent(current);
          }
          visit(current);
        }
      }
    }
  }
}
template <class Tree, class Visitor>
void bp_byte(const Tree& tree,
             std::span<const uint64_t> words,
             size_t nbits,
             Visitor& visit) {
  alignas(64) std::array<uint64_t, 136> numbers{}, positions{};
  size_t depth = 0, number = 0, base = 0;
  TreeNode current(0, 0);
  alignas(64) uint64_t event_numbers[8], event_positions[8];
  for (; base < nbits; base += std::min<size_t>(8, nbits - base)) {
    if (depth > 120 && nbits - base > 64) {
      bp_deep_word(tree, words, base, numbers, positions, depth, number,
                   current, visit);
      base += 56;
      continue;
    }
    unsigned pattern = unsigned(word_at(words, base) & 511);
    if ((pattern & 255) && nbits - base > 8 && depth <= 120) {
      const auto& p = byte_programs[pattern];
      decode_byte(p, base, number, depth, numbers, positions, event_numbers,
                  event_positions);
      for (unsigned i = 0; i < p.count; ++i) {
        current = TreeNode(event_numbers[i], event_positions[i]);
        visit(current);
      }
    } else {
      for (unsigned i = 0; i < 8 && base + i < nbits; ++i) {
        if ((pattern >> i) & 1) {
          current = TreeNode(number++, base + i);
          if (depth < 128) {
            numbers[depth + 8] = current.number;
            positions[depth + 8] = current.pos;
          }
          ++depth;
          visit(current);
        } else {
          --depth;
          if (depth && !((pattern >> (i + 1)) & 1)) {
            if (depth <= 128) {
              current = TreeNode(numbers[depth + 7], positions[depth + 7]);
            } else {
              current = tree.parent(current);
            }
            visit(current);
          }
        }
      }
    }
  }
}

}  // namespace pixie::detail::tree_unpack
