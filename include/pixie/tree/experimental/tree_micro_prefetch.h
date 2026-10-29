#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_cache.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

namespace pixie::experimental::detail {
// Cache decoding by topology, not by node address. Trailing bits beyond the
// first complete subtree are irrelevant. The table holds only relative handles.
struct BpMicroProgram {
  std::array<uint8_t, 24> events{};
  uint32_t returns = 0;
  uint8_t bits = 0, nodes = 0, count = 0;
};
inline constexpr auto bp_micro_programs = [] {
  std::array<BpMicroProgram, 8192> table{};
  for (unsigned key = 0; key < table.size(); ++key) {
    const unsigned bits = (key << 1) | 1;
    int balance = 0;
    unsigned length = 0;
    for (unsigned b = 0; b < 14; ++b) {
      balance += ((bits >> b) & 1) ? 1 : -1;
      if (!balance) {
        length = b + 1;
        break;
      }
    }
    if (!length) {
      continue;
    }
    auto& p = table[key];
    p.bits = length;
    std::array<uint8_t, 7> pos{}, parent{}, stack{}, degree{}, next{};
    next.fill(7);
    unsigned depth = 0;
    for (unsigned b = 0; b < length; ++b) {
      if ((bits >> b) & 1) {
        const auto i = p.nodes++;
        pos[i] = b;
        if (depth) {
          parent[i] = stack[depth - 1];
          ++degree[parent[i]];
        }
        stack[depth++] = i;
      } else {
        --depth;
      }
    }
    for (unsigned i = 1; i < p.nodes; ++i) {
      for (unsigned j = i + 1; j < p.nodes; ++j) {
        if (parent[i] == parent[j]) {
          next[i] = j;
          break;
        }
      }
    }
    auto event = [&](unsigned i, bool inspect, bool returning) {
      if (returning) {
        p.returns |= uint32_t{1} << p.count;
      }
      p.events[p.count++] = i | (pos[i] << 3) | (inspect ? 128 : 0);
    };
    for (unsigned i = 0; i < p.nodes; ++i) {
      event(i, false, false);
      if (degree[i]) {
        for (unsigned c = i + 1; c != 7; c = next[c]) {
          event(c, true, false);
        }
      } else {
        auto completed = i;
        while (completed && next[completed] == 7) {
          completed = parent[completed];
          event(completed, false, true);
        }
      }
    }
  }
  return table;
}();

struct BpMicroTables {
  std::array<uint16_t, 8192> headers{};
  std::array<BpMicroProgram, 197> programs{};
  size_t count = 0;
};
inline constexpr auto bp_micro_tables = [] {
  BpMicroTables out;
  std::array<uint16_t, 8192> ids{};
  for (size_t key = 0; key < 8192; ++key) {
    const auto& p = bp_micro_programs[key];
    if (!p.bits) {
      continue;
    }
    const auto canonical = key & ((size_t{1} << (p.bits - 1)) - 1);
    if (canonical == key) {
      ids[key] = out.count;
      out.programs[out.count++] = p;
    }
    out.headers[key] = p.bits | (ids[canonical] << 4);
  }
  return out;
}();
static_assert(bp_micro_tables.count == 197);

template <size_t Id, class Visit, class Inspect, size_t... I>
[[gnu::always_inline]] inline void RunMicroEvents(TreeNode root,
                                                  Visit& visit,
                                                  Inspect& inspect,
                                                  std::index_sequence<I...>) {
  auto event = [&]<size_t J>() {
    constexpr auto code = bp_micro_tables.programs[Id].events[J];
    const auto node =
        TreeNode(root.number + (code & 7), root.pos + ((code >> 3) & 15));
    if constexpr (code & 128) {
      inspect(node);
    } else {
      pixie::detail::tree_unpack::emit(
          visit, node,
          !(bp_micro_tables.programs[Id].returns & (uint32_t{1} << J)));
    }
  };
  (event.template operator()<I>(), ...);
}
template <size_t Id, class Visit, class Inspect>
[[gnu::noinline]] void RunMicroProgram(TreeNode root,
                                       Visit& visit,
                                       Inspect& inspect) {
  RunMicroEvents<Id>(
      root, visit, inspect,
      std::make_index_sequence<bp_micro_tables.programs[Id].count>{});
}
template <class Visit, class Inspect, size_t... I>
constexpr auto MicroDispatch(std::index_sequence<I...>) {
  return std::array{&RunMicroProgram<I, Visit, Inspect>...};
}

template <size_t Bits,
          unsigned Compiled = 0,
          bool Reservoir = false,
          bool Absolute = false>
class BpMicroSubtree {
  static_assert(Compiled == 0 || Compiled == 1 || Compiled == 3 ||
                Compiled == 4);
  static_assert(Bits >= 128 && Bits % 128 == 0 && Bits <= 32768);
  static constexpr uint16_t missing = Bits / 2;
  using Position = std::conditional_t<Absolute, size_t, uint16_t>;
  struct Record {
    Position number, pos;
    uint16_t parent, next, degree, program, returns;
  };
  std::array<Record, missing + 1> records_;
  size_t count_ = 0, nodes_ = 0;
  TreeNode root_{0, 0};
  TreeNode node(size_t i) const {
    if constexpr (Absolute) {
      return {records_[i].number, records_[i].pos};
    } else {
      return {root_.number + records_[i].number, root_.pos + records_[i].pos};
    }
  }

 public:
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
    root_ = root;
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
      if constexpr (Reservoir) {
        if (valid < 14) {
          pattern = peek(bit);
          valid = 64;
        }
      } else {
        pattern = peek(bit);
      }
      const auto closes = std::min<size_t>(
          std::countr_zero(pattern),
          std::min(length - bit, Reservoir ? valid : size_t{64}));
      if (count_) {
        records_[count_ - 1].returns += closes;
      }
      depth -= closes;
      bit += closes;
      if (bit == length) {
        break;
      }
      if constexpr (Reservoir) {
        pattern = closes == 64 ? 0 : pattern >> closes;
        valid -= closes;
        if (valid < 14) {
          pattern = peek(bit);
          valid = 64;
        }
        if (!(pattern & 1)) {
          continue;
        }
      } else {
        if (closes == 64) {
          continue;
        }
        if (closes > 50) {
          pattern = peek(bit);
        } else {
          pattern >>= closes;
        }
      }
      const uint16_t key = (pattern >> 1) & 8191;
      const auto header = bp_micro_tables.headers[key];
      const auto length = header & 15;
      const size_t i = count_++;
      records_[i] = {Position(number + (Absolute ? root.number : 0)),
                     Position(bit + (Absolute ? root.pos : 0)),
                     0,
                     missing,
                     0,
                     uint16_t(header >> 4),
                     0};
      last[i] = missing;
      if (depth) {
        auto parent = stack[depth - 1];
        records_[i].parent = parent;
        records_[last[parent]].next = i;
        last[parent] = i;
        ++records_[parent].degree;
      }
      if (length) {
        bit += length;
        number += length / 2;
      } else {
        stack[depth++] = i;
        ++bit;
        ++number;
      }
      if constexpr (Reservoir) {
        const size_t consumed = length ? length : 1;
        pattern >>= consumed;
        valid -= consumed;
      }
    }
    return true;
  }
  size_t size() const { return nodes_; }
  template <class Visit, class Inspect>
  void walk(Visit& visit, Inspect& inspect) const {
    using pixie::detail::tree_unpack::emit;
    for (size_t i = 0; i < count_; ++i) {
      const auto& record = records_[i];
      if (record.degree) {
        emit(visit, node(i), true);
        for (size_t c = i + 1; c != missing; c = records_[c].next) {
          inspect(node(c));
        }
      } else {
        if constexpr (Compiled == 4) {
          if (record.program == 0) {
            RunMicroEvents<0>(
                node(i), visit, inspect,
                std::make_index_sequence<bp_micro_tables.programs[0].count>{});
          } else if (record.program == 1) {
            RunMicroEvents<1>(
                node(i), visit, inspect,
                std::make_index_sequence<bp_micro_tables.programs[1].count>{});
          } else {
            const auto& p = bp_micro_tables.programs[record.program];
            alignas(64) uint64_t handles[48];
            const auto root = node(i);
            pixie::detail::tree_unpack::decode_relative_events(
                p.events.data(), p.count, root.number, root.pos, handles);
            for (unsigned j = 0; j < p.count; ++j) {
              TreeNode current(handles[2 * j], handles[2 * j + 1]);
              if (p.events[j] & 128) {
                inspect(current);
              } else {
                emit(visit, current, !((p.returns >> j) & 1));
              }
            }
          }
        } else if constexpr (Compiled == 3) {
          switch (record.program) {
            case 0:
              RunMicroEvents<0>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    bp_micro_tables.programs[0].count>{});
              break;
            case 1:
              RunMicroEvents<1>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    bp_micro_tables.programs[1].count>{});
              break;
            case 2:
              RunMicroEvents<2>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    bp_micro_tables.programs[2].count>{});
              break;
            case 3:
              RunMicroEvents<3>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    bp_micro_tables.programs[3].count>{});
              break;
            default: {
              static constexpr auto dispatch = MicroDispatch<Visit, Inspect>(
                  std::make_index_sequence<197>{});
              dispatch[record.program](node(i), visit, inspect);
            }
          }
        } else if constexpr (Compiled == 1) {
          static constexpr auto dispatch =
              MicroDispatch<Visit, Inspect>(std::make_index_sequence<197>{});
          dispatch[record.program](node(i), visit, inspect);
        } else {
          const auto& p = bp_micro_tables.programs[record.program];
          const auto root = node(i);
          for (unsigned j = 0; j < p.count; ++j) {
            const unsigned code = p.events[j];
            TreeNode current(root.number + (code & 7),
                             root.pos + ((code >> 3) & 15));
            if (code & 128) {
              inspect(current);
            } else {
              emit(visit, current, !((p.returns >> j) & 1));
            }
          }
        }
        size_t completed = i;
        for (size_t r = 0; r < record.returns; ++r) {
          completed = records_[completed].parent;
          emit(visit, node(completed), false);
        }
      }
    }
  }
};

// Local traversal scratch; topology programs are immutable and shared. No
// cache admission/lookups inside a microtree, and no tree-specific table build.
template <size_t Bits,
          class Tree,
          class Visit,
          class Inspect,
          unsigned Compiled = 0>
class BpMicroTraversal {
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Visit& visit_;
  Inspect& inspect_;
  BpMicroSubtree<Bits,
                 Compiled >= 7 ? 1 : Compiled,
                 Compiled == 7 || Compiled == 8,
                 Compiled == 8>
      cache_;

 public:
  size_t nodes = 0, children = 0;
  BpMicroTraversal(const Tree& tree,
                   std::span<const uint64_t> words,
                   Visit& visit,
                   Inspect& inspect)
      : tree_(tree), words_(words), visit_(visit), inspect_(inspect) {}
  void walk(TreeNode node, size_t depth = 0) {
    using namespace pixie::detail::tree_unpack;
    if (tree_.is_leaf(node)) {
      ++nodes;
      emit(visit_, node, true);
      return;
    }
    if (cache_.load(words_, 2 * tree_.size(), node)) {
      cache_.walk(visit_, inspect_);
      nodes += cache_.size();
      children += cache_.size() - 1;
      return;
    }
    if (depth == 128) {
      auto events = [&](TreeNode n, bool entry) {
        emit(visit_, n, entry);
        if (entry) {
          ++nodes;
          tree_.for_each_child(n, [&](auto c) {
            ++children;
            inspect_(c);
          });
        }
      };
      subtree(tree_, node, events);
      return;
    }
    ++nodes;
    emit(visit_, node, true);
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    tree_.for_each_child(node, [&](auto c) {
      inspect_(c);
      if (count < 32) {
        numbers[count] = c.number;
        positions[count] = c.pos;
      }
      ++count;
    });
    children += count;
    if (count <= 32) {
      for (size_t i = 0; i < count; ++i) {
        walk({numbers[i], positions[i]}, depth + 1);
      }
    } else {
      tree_.for_each_child(node, [&](auto c) { walk(c, depth + 1); });
    }
    emit(visit_, node, false);
  }
};
}  // namespace pixie::experimental::detail
