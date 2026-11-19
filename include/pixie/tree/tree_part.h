#pragma once

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace pixie::tree {

/** @brief Owning inspection snapshot of the nodes and edges in one tree part.
 * @details Roots, child links and event operands index nodes() from zero.
 * A local root has no parent in the part; its source parent may exist. Node
 * order is producer-specific. No record borrows a tree, index or decoder.
 * Inspection allocations are separate from navigation and cache budgets.
 */
class TreePart {
 public:
  /** @brief Retained index, decoded forest unit, or complete decoder hit. */
  enum class Kind { kRootIndex, kForestUnit, kBpEvents, kDfudsEvents };
  /** @brief Zero-based half-open interval [begin,end). */
  struct Range {
    std::size_t begin, end;
    bool operator==(const Range&) const = default;
  };
  /** @brief A source node and the ordered subset of its children in this part.
   */
  struct Node {
    /// Number in the source encoding, when available; not necessarily preorder.
    std::optional<std::size_t> source_id;
    /// Original encoding bit position, when meaningful.
    std::optional<std::size_t> bit_position;
    /// Total source degree, including children outside the part.
    std::size_t degree = 0;
    /// True exactly when children contains every source child. An empty
    /// complete list denotes a leaf; an empty incomplete list does not.
    bool children_complete = true;
    std::vector<std::size_t> children;
    /// For forest portals, the range of external units holding child roots.
    std::optional<Range> child_units;
  };
  /** @brief One inspected-DFS callback: enter/return visits name the visited
   * node, while inspect names the child observed before descending. */
  struct Event {
    enum class Type { kEnter, kInspect, kReturn };
    Type type;
    std::size_t node;
    bool operator==(const Event&) const = default;
  };

  /** @brief Take ownership of producer records. The producer supplies valid
   * local indices and lists in source child order; no validation is performed.
   */
  TreePart(Kind kind,
           std::vector<Node> nodes,
           std::vector<std::size_t> roots,
           std::vector<Event> events = {},
           std::optional<Range> bits = {})
      : kind_(kind),
        nodes_(std::move(nodes)),
        roots_(std::move(roots)),
        events_(std::move(events)),
        bits_(bits) {}

  /** @brief Producer and boundary convention. */
  Kind kind() const { return kind_; }
  /** @brief True when no nodes are represented, including a disabled index. */
  bool empty() const { return nodes_.empty(); }
  /** @brief Local node records, owned by this snapshot. */
  const std::vector<Node>& nodes() const { return nodes_; }
  /** @brief Local roots in source order. */
  const std::vector<std::size_t>& roots() const { return roots_; }
  /** @brief Decoder callback sequence; empty for non-event parts. */
  const std::vector<Event>& events() const { return events_; }
  /** @brief Exact source bit interval of an event part; absent otherwise. */
  std::optional<Range> source_bits() const { return bits_; }

 private:
  Kind kind_;
  std::vector<Node> nodes_;
  std::vector<std::size_t> roots_;
  std::vector<Event> events_;
  std::optional<Range> bits_;
};
}  // namespace pixie::tree
