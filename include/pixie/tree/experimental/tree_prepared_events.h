#pragma once

#include <pixie/tree/detail/unpack.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace pixie::experimental {
/**
 * @brief Bounded immutable DFS event blocks for a complete small tree.
 * @details Construction admits nonempty trees of at most MaxNodes nodes.
 * Events preserve entry, immediate-child inspection, and ancestor-return
 * order. Emitted handles keep the source tree's lifetime requirements.
 * Empty and oversize trees produce misses. Replay uses no mutable cache state.
 */
template <size_t MaxNodes = 4096>
class PreparedTreeEvents {
  static_assert(MaxNodes <= (std::numeric_limits<size_t>::max() - 5) / 3);
  struct Block {
    std::array<size_t, 16> handles{};
    uint8_t inspections = 0, returns = 0, count = 0;
  };
  std::vector<Block> blocks_, visits_;
  size_t nodes_ = 0;

  static void append(std::vector<Block>& blocks,
                     TreeNode node,
                     bool inspection,
                     bool returning) {
    if (blocks.empty() || blocks.back().count == 8) {
      blocks.emplace_back();
    }
    auto& b = blocks.back();
    const unsigned i = b.count++;
    b.handles[2 * i] = node.number;
    b.handles[2 * i + 1] = node.pos;
    b.inspections |= unsigned(inspection) << i;
    b.returns |= unsigned(returning) << i;
  }

  template <size_t Mask, class Visitor, class Inspector, size_t... I>
  static void run_events(const Block& b,
                         Visitor& visit,
                         Inspector& inspect,
                         std::index_sequence<I...>) {
    auto event = [&]<size_t J>() {
      const TreeNode node(b.handles[2 * J], b.handles[2 * J + 1]);
      if constexpr ((Mask >> J) & 1) {
        inspect(node);
      } else {
        pixie::detail::tree_unpack::emit(visit, node, !((b.returns >> J) & 1));
      }
    };
    (event.template operator()<I>(), ...);
  }
  template <size_t Mask, class Visitor, class Inspector>
  [[gnu::noinline]] static void run(const Block& b, Visitor& v, Inspector& i) {
    run_events<Mask>(b, v, i, std::make_index_sequence<8>{});
  }
  template <class Visitor, class Inspector, size_t... I>
  static constexpr auto dispatch(std::index_sequence<I...>) {
    return std::array{&run<I, Visitor, Inspector>...};
  }

 public:
  /** @brief Prepare complete root events; empty/oversize trees are not read. */
  template <class Tree>
  explicit PreparedTreeEvents(const Tree& tree) {
    if (tree.empty() || tree.size() > MaxNodes) {
      return;
    }
    nodes_ = tree.size();
    blocks_.reserve((3 * nodes_ - 2 + 7) / 8);
    visits_.reserve((2 * nodes_ - 1 + 7) / 8);
    auto event = [&](TreeNode node, bool entering) {
      append(blocks_, node, false, !entering);
      append(visits_, node, false, !entering);
      if (entering) {
        tree.for_each_child(
            node, [&](TreeNode child) { append(blocks_, child, true, false); });
      }
    };
    pixie::detail::tree_unpack::subtree(tree, tree.root(), event);
  }
  /** @brief Admitted node count, or zero for a miss. */
  size_t size() const { return nodes_; }
  /** @brief Object and reserved event storage; excludes shared dispatch tables.
   */
  size_t memory_usage_bytes() const {
    return sizeof(*this) +
           (blocks_.capacity() + visits_.capacity()) * sizeof(Block);
  }
  /**
   * @brief Replay events, or return false without callbacks for a miss.
   * @details Compiled selects unrolled eight-event blocks; the final partial
   * block is interpreted. InspectChildren=false selects a visit-only stream
   * and never invokes the inspector. Callback exceptions propagate. Calls may
   * nest.
   */
  template <bool Compiled = false,
            bool InspectChildren = true,
            class Visitor,
            class Inspector>
  bool try_walk(Visitor& visit, Inspector& inspect) const {
    if (!nodes_) {
      return false;
    }
    const auto& blocks = InspectChildren ? blocks_ : visits_;
    for (const auto& b : blocks) {
      if constexpr (Compiled) {
        if (b.count == 8) {
          static constexpr auto table =
              dispatch<Visitor, Inspector>(std::make_index_sequence<256>{});
          table[b.inspections](b, visit, inspect);
          continue;
        }
      }
      for (unsigned i = 0; i < b.count; ++i) {
        const TreeNode node(b.handles[2 * i], b.handles[2 * i + 1]);
        if (InspectChildren && ((b.inspections >> i) & 1)) {
          inspect(node);
        } else {
          pixie::detail::tree_unpack::emit(visit, node,
                                           !((b.returns >> i) & 1));
        }
      }
    }
    return true;
  }
};
}  // namespace pixie::experimental
