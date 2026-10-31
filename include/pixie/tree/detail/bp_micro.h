#pragma once

// Selected BP decoder: absolute handles, 24-byte records, reservoir loading,
// macro-node links and compiled seven-node topology programs.
#include <pixie/tree/detail/bp_micro_tables.h>
#include <pixie/tree/detail/bp_micro_vector.h>
#include <pixie/tree/detail/micro_traversal.h>
namespace pixie::detail::bp_micro {
template <size_t Bits>
class BpMicroSubtree {
  static_assert(Bits >= 128 && Bits % 128 == 0 && Bits <= 32768);
  static constexpr uint16_t missing = Bits / 2;
  using Position = size_t;
  struct Record {
    Position number, pos;
    uint16_t parent, next, program, returns;
  };
  std::array<Record, missing + 1> records_;
  size_t count_ = 0, nodes_ = 0;
  TreeNode node(size_t i) const {
    return {records_[i].number, records_[i].pos};
  }

 public:
  static constexpr size_t bit_size(size_t nodes) { return 2 * nodes; }
  static constexpr size_t load_scratch_bytes() {
    return (Bits / 64 + 1) * sizeof(uint64_t) +
           2 * (Bits / 2) * sizeof(uint16_t);
  }
  bool load(std::span<const uint64_t> words, size_t nbits, TreeNode root) {
    using namespace pixie::detail::tree_unpack;
    std::array<uint64_t, Bits / 64 + 1> local{};
    const size_t available = std::min<size_t>(Bits, nbits - root.pos);
    for (size_t offset = 0; offset < available; offset += 64) {
      local[offset / 64] = word_at(words, root.pos + offset);
    }
    size_t length = 0;
    int target = 0;
    for (size_t offset = 0; offset < available; offset += 128) {
      int delta;
      const size_t found =
          forward_search_128(local.data() + offset / 64, target, 0, &delta);
      if (found < std::min<size_t>(128, available - offset)) {
        length = offset + found + 1;
        break;
      }
      target -= delta;
    }
    if (!length) {
      return false;
    }
    nodes_ = length / 2;
    count_ = 0;
    std::array<uint16_t, missing> stack, last;
    size_t depth = 0, number = 0, bit = 0;
    auto peek = [&](size_t pos) {
      const size_t word = pos / 64;
      const unsigned shift = pos % 64;
      auto value = local[word] >> shift;
      if (shift) {
        value |= local[word + 1] << (64 - shift);
      }
      return value;
    };
    // Reservoir mode reuses the unconsumed source bits between records. Refill
    // before a 14-bit shape lookup; bound closing runs by the valid bit count.
    uint64_t pattern = 0;
    size_t valid = 0;
    while (bit < length) {
      if (valid < 14) {
        pattern = peek(bit);
        valid = 64;
      }
      const auto closes = std::min<size_t>(std::countr_zero(pattern),
                                           std::min(length - bit, valid));
      if (count_) {
        records_[count_ - 1].returns += closes;
      }
      depth -= closes;
      bit += closes;
      if (bit == length) {
        break;
      }
      pattern = closes == 64 ? 0 : pattern >> closes;
      valid -= closes;
      if (valid < 14) {
        pattern = peek(bit);
        valid = 64;
      }
      if (!(pattern & 1)) {
        continue;
      }
      const uint16_t key = (pattern >> 1) & 8191;
      const auto header = bp_micro_tables.headers[key];
      const auto length = header & 15;
      const size_t i = count_++;
      records_[i] = {Position(number + root.number),
                     Position(bit + root.pos),
                     0,
                     missing,
                     uint16_t(header >> 4),
                     0};

      if (depth) {
        auto parent = stack[depth - 1];
        records_[i].parent = parent;
        records_[last[parent]].next = i;
        last[parent] = i;
      }
      if (length) {
        bit += length;
        number += length / 2;
      } else {
        records_[i].program = 256;
        last[i] = missing;
        stack[depth++] = i;
        ++bit;
        ++number;
      }
      const size_t consumed = length ? length : 1;
      pattern >>= consumed;
      valid -= consumed;
    }
    return true;
  }
  size_t size() const { return nodes_; }
  template <class Visit, class Inspect>
  void walk(Visit& visit, Inspect& inspect) const {
    using pixie::detail::tree_unpack::emit;
    for (size_t i = 0; i < count_; ++i) {
      const auto& record = records_[i];
      if (record.program & 256) {
        emit(visit, node(i), true);
        for (size_t c = i + 1; c != missing; c = records_[c].next) {
          inspect(node(c));
        }
      } else {
        static constexpr auto dispatch = vector_micro::dispatch<Visit, Inspect>(
            std::make_index_sequence<197>{});
        dispatch[record.program](vector_micro::pack(node(i)), visit, inspect);
        size_t completed = i;
        for (size_t r = 0; r < record.returns; ++r) {
          completed = records_[completed].parent;
          emit(visit, node(completed), false);
        }
      }
    }
  }
};

// Shared orchestration with an encoding-specific subtree decoder.
template <size_t Bits, class Tree, class Visit, class Inspect>
using BpMicroTraversal =
    tree_unpack::MicroTraversal<BpMicroSubtree<Bits>, Tree, Visit, Inspect>;
}  // namespace pixie::detail::bp_micro
