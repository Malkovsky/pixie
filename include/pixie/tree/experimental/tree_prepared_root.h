#pragma once

#include <pixie/tree.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>

namespace pixie::experimental {
/**
 * @brief Retain a complete decoded root subtree when its bits fit a fixed cap.
 * @details Decoder supplies load(words, nbits, root), size(), and const
 * walk(visitor, inspector). Construction requires a valid nonempty immutable
 * tree root and its logical bit length, excluding padding. Handles emitted by
 * callbacks retain the underlying tree's lifetime requirements. Oversize roots
 * are not decoded. Independent calls use immutable records and may nest.
 */
template <class Decoder, size_t MaxBits>
class PreparedRootSubtree {
 public:
  PreparedRootSubtree(std::span<const uint64_t> words,
                      size_t nbits,
                      TreeNode root) {
    if (!nbits || nbits > MaxBits) {
      return;
    }
    auto decoded = std::make_unique<Decoder>();
    if (decoded->load(words, nbits, root)) {
      decoded_ = std::move(decoded);
    }
  }

  /** @brief Cached node count, or zero when the root was not admitted. */
  size_t size() const { return decoded_ ? decoded_->size() : 0; }

  /** @brief Object and decoded allocation bytes; excludes shared tables. */
  size_t memory_usage_bytes() const {
    return sizeof(*this) + (decoded_ ? sizeof(Decoder) : 0);
  }

  /**
   * @brief Replay exact DFS entry, child-inspection, and return events.
   * @return False on a cache miss without invoking either callback. On a hit,
   * returns true after the complete traversal. Callback exceptions propagate.
   */
  template <class Visitor, class Inspector>
  bool try_walk(Visitor& visit, Inspector& inspect) const {
    if (!decoded_) {
      return false;
    }
    decoded_->walk(visit, inspect);
    return true;
  }

 private:
  std::unique_ptr<Decoder> decoded_;
};
}  // namespace pixie::experimental
