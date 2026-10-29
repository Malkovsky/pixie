#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_cache.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

namespace pixie::experimental::dfuds_detail {
// Cache decoding by topology, not by node address. DFUDS degree runs encode
// complete subtrees. Trailing bits beyond the first complete subtree are
// irrelevant. The table holds only relative handles.
struct DfudsMicroProgram {
  std::array<uint8_t, 24> events{};
  uint32_t returns = 0;
  uint8_t bits = 0, nodes = 0, count = 0;
};
inline constexpr auto dfuds_micro_programs = [] {
  std::array<DfudsMicroProgram, 8192> table{};
  for (unsigned key = 0; key < table.size(); ++key) {
    const unsigned bits = key;
    int balance = 1;
    unsigned length = 0;
    for (unsigned b = 0; b < 13; ++b) {
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
    std::array<uint8_t, 7> pos{}, parent{}, stack{}, degree{}, next{},
        remaining{};
    next.fill(7);
    unsigned depth = 0;
    for (unsigned b = 0; b < length;) {
      const auto i = p.nodes++;
      pos[i] = b;
      while ((bits >> b) & 1) {
        ++degree[i];
        ++b;
      }
      ++b;
      if (depth) {
        parent[i] = stack[depth - 1];
        --remaining[depth - 1];
      }
      if (degree[i]) {
        stack[depth] = i;
        remaining[depth++] = degree[i];
      } else {
        while (depth && !remaining[depth - 1]) {
          --depth;
        }
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

struct DfudsMicroTables {
  std::array<uint16_t, 8192> headers{};
  std::array<DfudsMicroProgram, 197> programs{};
  size_t count = 0;
};
inline constexpr auto dfuds_micro_tables = [] {
  DfudsMicroTables out;
  std::array<uint16_t, 8192> ids{};
  for (size_t key = 0; key < 8192; ++key) {
    const auto& p = dfuds_micro_programs[key];
    if (!p.bits) {
      continue;
    }
    const auto canonical = key & ((size_t{1} << p.bits) - 1);
    if (canonical == key) {
      ids[key] = out.count;
      out.programs[out.count++] = p;
    }
    out.headers[key] = p.bits | (ids[canonical] << 4);
  }
  return out;
}();
static_assert(dfuds_micro_tables.count == 197);

template <size_t Id, class Visit, class Inspect, size_t... I>
[[gnu::always_inline]] inline void RunMicroEvents(TreeNode root,
                                                  Visit& visit,
                                                  Inspect& inspect,
                                                  std::index_sequence<I...>) {
  auto event = [&]<size_t J>() {
    constexpr auto code = dfuds_micro_tables.programs[Id].events[J];
    const auto node =
        TreeNode(root.number + (code & 7), root.pos + ((code >> 3) & 15));
    if constexpr (code & 128) {
      inspect(node);
    } else {
      pixie::detail::tree_unpack::emit(
          visit, node,
          !(dfuds_micro_tables.programs[Id].returns & (uint32_t{1} << J)));
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
      std::make_index_sequence<dfuds_micro_tables.programs[Id].count>{});
}
template <class Visit, class Inspect, size_t... I>
constexpr auto MicroDispatch(std::index_sequence<I...>) {
  return std::array{&RunMicroProgram<I, Visit, Inspect>...};
}

template <size_t Bits,
          unsigned Compiled = 0,
          bool Reservoir = false,
          bool Absolute = false,
          bool Excess = false,
          bool Compact = false>
class DfudsMicroSubtree {
  static_assert(Compiled == 0 || Compiled == 1 || Compiled == 3 ||
                Compiled == 4);
  static_assert(Bits >= 128 && Bits % 128 == 0 && Bits <= 32768);
  static constexpr uint16_t missing = Bits / 2;
  using Position = std::conditional_t<Absolute, size_t, uint16_t>;
  struct FullRecord {
    Position number, pos;
    uint16_t parent, next, degree, program, returns;
  };
  struct CompactRecord {
    Position number, pos;
    uint16_t parent, next, program, returns;
  };
  using Record = std::conditional_t<Compact, CompactRecord, FullRecord>;
  static_assert(!Compact || !Absolute || sizeof(Record) == 24);
  static bool is_macro(const Record& record) {
    if constexpr (Compact) {
      return record.program & 256;
    } else {
      return record.degree != 0;
    }
  }
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
    int target = -1;
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
    nodes_ = (length + 1) / 2;
    count_ = 0;
    std::array<uint16_t, missing> stack, last;
    std::array<std::conditional_t<Excess, int16_t, uint16_t>, missing>
        remaining;
    // A fitting chunk has at most Bits / 2 nodes. Completion targets range
    // from -1 to fewer than that many pending edges, fitting int16_t.
    [[maybe_unused]] int excess = 0;
    size_t depth = 0, number = 0, bit = 0;
    while (bit < length) {
      const auto pattern = word_at(local, bit);
      const auto header = dfuds_micro_tables.headers[pattern & 8191];
      const size_t micro_length = header & 15;
      size_t degree = 0;
      if (!micro_length) {
        degree = std::countr_one(pattern);
        while (degree % 64 == 0) {
          const size_t more = std::countr_one(word_at(local, bit + degree));
          degree += more;
          if (more < 64) {
            break;
          }
        }
      }
      const size_t i = count_++;
      if constexpr (Compact) {
        records_[i] = {Position(number + (Absolute ? root.number : 0)),
                       Position(bit + (Absolute ? root.pos : 0)),
                       0,
                       missing,
                       uint16_t((header >> 4) | (micro_length ? 0 : 256)),
                       0};
      } else {
        records_[i] = {Position(number + (Absolute ? root.number : 0)),
                       Position(bit + (Absolute ? root.pos : 0)),
                       0,
                       missing,
                       uint16_t(degree),
                       uint16_t(header >> 4),
                       0};
      }
      if (depth) {
        const size_t parent = stack[depth - 1];
        if constexpr (!Excess) {
          --remaining[depth - 1];
        }
        records_[i].parent = parent;
        records_[last[parent]].next = i;
        last[parent] = i;
      }
      if (micro_length) {
        bit += micro_length;
        number += (micro_length + 1) / 2;
        if constexpr (Excess) {
          --excess;
          while (depth && remaining[depth - 1] == excess) {
            --depth;
            ++records_[i].returns;
          }
        } else {
          while (depth && !remaining[depth - 1]) {
            --depth;
            ++records_[i].returns;
          }
        }
      } else {
        // Only macro records become parents in the decoded skeleton.
        last[i] = missing;
        stack[depth] = i;
        if constexpr (Excess) {
          remaining[depth++] = excess - 1;
          excess += int(degree) - 1;
        } else {
          remaining[depth++] = degree;
        }
        bit += degree + 1;
        ++number;
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
      if (is_macro(record)) {
        emit(visit, node(i), true);
        for (size_t c = i + 1; c != missing; c = records_[c].next) {
          inspect(node(c));
        }
      } else {
        if constexpr (Compiled == 4) {
          if (record.program == 0) {
            RunMicroEvents<0>(node(i), visit, inspect,
                              std::make_index_sequence<
                                  dfuds_micro_tables.programs[0].count>{});
          } else if (record.program == 1) {
            RunMicroEvents<1>(node(i), visit, inspect,
                              std::make_index_sequence<
                                  dfuds_micro_tables.programs[1].count>{});
          } else {
            const auto& p = dfuds_micro_tables.programs[record.program];
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
                                    dfuds_micro_tables.programs[0].count>{});
              break;
            case 1:
              RunMicroEvents<1>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    dfuds_micro_tables.programs[1].count>{});
              break;
            case 2:
              RunMicroEvents<2>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    dfuds_micro_tables.programs[2].count>{});
              break;
            case 3:
              RunMicroEvents<3>(node(i), visit, inspect,
                                std::make_index_sequence<
                                    dfuds_micro_tables.programs[3].count>{});
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
          const auto& p = dfuds_micro_tables.programs[record.program];
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
class DfudsMicroTraversal {
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Visit& visit_;
  Inspect& inspect_;
  size_t max_depth_;
  DfudsMicroSubtree<Bits,
                    Compiled >= 7 ? 1 : Compiled,
                    Compiled == 7 || Compiled == 8,
                    Compiled >= 8,
                    Compiled >= 9,
                    Compiled == 10>
      cache_;

 public:
  size_t nodes = 0, children = 0;
  // Decoded records, load-local arrays and recursive child buffers.
  static constexpr size_t scratch_bytes(size_t depth = 128) {
    return sizeof(DfudsMicroTraversal) +
           std::min(depth, size_t{128}) * 64 * sizeof(size_t) +
           (Bits / 64 + 1) * sizeof(uint64_t) +
           3 * (Bits / 2) * sizeof(uint16_t) + 1024;
  }
  DfudsMicroTraversal(const Tree& tree,
                      std::span<const uint64_t> words,
                      Visit& visit,
                      Inspect& inspect,
                      size_t max_depth = 128)
      : tree_(tree),
        words_(words),
        visit_(visit),
        inspect_(inspect),
        max_depth_(std::min(max_depth, size_t{128})) {}
  void walk(TreeNode node, size_t depth = 0) {
    using namespace pixie::detail::tree_unpack;
    if (tree_.is_leaf(node)) {
      ++nodes;
      emit(visit_, node, true);
      return;
    }
    if (cache_.load(words_, 2 * tree_.size() - 1, node)) {
      cache_.walk(visit_, inspect_);
      nodes += cache_.size();
      children += cache_.size() - 1;
      return;
    }
    if (depth == max_depth_) {
      walk_uncached(node);
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

 private:
  // Reuse the decoder at arbitrary depth without growing child-buffer state.
  void walk_uncached(TreeNode root) {
    using namespace pixie::detail::tree_unpack;
    auto view = tree_.navigate(root);
    bool descending = true;
    while (true) {
      if (descending && !view.is_leaf() &&
          cache_.load(words_, 2 * tree_.size() - 1, view.node())) {
        cache_.walk(visit_, inspect_);
        nodes += cache_.size();
        children += cache_.size() - 1;
      } else {
        emit(visit_, view.node(), descending);
        if (descending) {
          ++nodes;
          tree_.for_each_child(view.node(), [&](auto c) {
            ++children;
            inspect_(c);
          });
          if (!view.is_leaf()) {
            view.reset(view.first_child());
            continue;
          }
        }
      }
      descending = false;
      if (view.node().number == root.number) {
        return;
      }
      if (view.is_last_child()) {
        view.reset(view.parent());
      } else {
        view.reset(view.next_sibling());
        descending = true;
      }
    }
  }
};
}  // namespace pixie::experimental::dfuds_detail
