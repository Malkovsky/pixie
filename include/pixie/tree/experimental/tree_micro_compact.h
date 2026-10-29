#pragma once

// Experimental compact BP decoded records. Full absolute handles plus four
// 16-bit navigation/program fields occupy 24 bytes. MacroLast initializes
// last-child slots only for macro nodes. This is an opt-in traversal scratch
// layout, not a persistent cache or a change to BPTree defaults.
// 2026-10-26 registered 10M inspected DFS: 81.33 ms versus retained 84.66 ms
// (Release native, CPU 0, seven interleaved repetitions, 0.3 s minimum).
// BFS pointers remain faster at 71.07 ms in a separate reference run.
// See profiling/bp-compact-registered-2026-10-26/REPORT.md for scope and
// limits.
#include <pixie/tree/experimental/tree_micro_prefetch.h>
#include <pixie/tree/experimental/tree_micro_vector.h>
namespace pixie::experimental::compact_bp {
using namespace pixie::experimental::detail;
template <size_t Bits,
          unsigned Compiled = 0,
          bool Reservoir = false,
          bool Absolute = false,
          bool MacroLast = false,
          unsigned VectorPolicy = 0>
class BpMicroSubtree {
  static_assert(Compiled == 0 || Compiled == 1 || Compiled == 3 ||
                Compiled == 4);
  static_assert(Bits >= 128 && Bits % 128 == 0 && Bits <= 32768);
  static constexpr uint16_t missing = Bits / 2;
  using Position = std::conditional_t<Absolute, size_t, uint16_t>;
  struct Record {
    Position number, pos;
    uint16_t parent, next, program, returns;
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
                     uint16_t(header >> 4),
                     0};
      if constexpr (!MacroLast) {
        last[i] = missing;
      }
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
        if constexpr (MacroLast) {
          last[i] = missing;
        }
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
      if (record.program & 256) {
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
          if constexpr (VectorPolicy) {
            static constexpr auto dispatch =
                vector_micro::dispatch<VectorPolicy, Visit, Inspect>(
                    std::make_index_sequence<197>{});
            dispatch[record.program](vector_micro::pack(node(i)), visit,
                                     inspect);
          } else {
            static constexpr auto dispatch =
                MicroDispatch<Visit, Inspect>(std::make_index_sequence<197>{});
            dispatch[record.program](node(i), visit, inspect);
          }
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
          unsigned Compiled = 0,
          bool MacroLast = false,
          unsigned VectorPolicy = 0>
class BpMicroTraversal {
  const Tree& tree_;
  std::span<const uint64_t> words_;
  Visit& visit_;
  Inspect& inspect_;
  size_t max_depth_;
  BpMicroSubtree<Bits,
                 Compiled >= 7 ? 1 : Compiled,
                 Compiled == 7 || Compiled == 8,
                 Compiled == 8,
                 MacroLast,
                 VectorPolicy>
      cache_;

 public:
  size_t nodes = 0, children = 0;
  // Includes load-local arrays and all simultaneously live recursive child
  // buffers. Compiler call-frame overhead and user callback state are separate.
  static constexpr size_t scratch_bytes(size_t depth = 128) {
    return sizeof(BpMicroTraversal) +
           std::min(depth, size_t{128}) * 64 * sizeof(size_t) +
           (Bits / 64 + 1) * sizeof(uint64_t) +
           2 * (Bits / 2) * sizeof(uint16_t) + 1024;
  }
  BpMicroTraversal(const Tree& tree,
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
    if (cache_.load(words_, 2 * tree_.size(), node)) {
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
  // Keep decoding below the child-buffer depth limit. The navigation view
  // supplies ascent/sibling state without another recursive child array.
  void walk_uncached(TreeNode root) {
    using namespace pixie::detail::tree_unpack;
    auto view = tree_.navigate(root);
    bool descending = true;
    while (true) {
      if (descending && !view.is_leaf() &&
          cache_.load(words_, 2 * tree_.size(), view.node())) {
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
}  // namespace pixie::experimental::compact_bp
