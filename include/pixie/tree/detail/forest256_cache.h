#pragma once

#include <pixie/tree/detail/forest256_simd.h>

#include <cstring>
#include <memory>
#include <utility>

namespace pixie::tree::detail::forest256 {

/// Local BFS coordinates; [first, first + degree) remains usable after
/// eviction.
struct NodeInfo {
  uint16_t roots, first, degree;
  uint8_t parent;  ///< kNoParent for roots.
  bool last;       ///< Includes the last root in the root sibling group.
};

/// Borrows the stream; the reader does not extend its owner's lifetime. Keep
/// the words alive and unchanged during lookup. After an in-place edit, owner
/// move, or storage replacement, bind each affected reader to a live stream
/// before its next lookup, even if the address and size have not changed.
/// Copies borrow the same stream and allocate their own empty slots (O(K)).
/// Moves transfer the slots and binding, leaving an unbound zero-slot reader.
/// A reader is not thread-safe. Child enumeration needs only a copied NodeInfo.
class Forest256Cache {
  struct Slot {
    LevelCsr256 csr;
    size_t unit;
    bool valid = false;
  };

 public:
  /// All requested heap bytes, including tag, validity and alignment rounding,
  /// are charged. Allocator-internal metadata is outside this byte model.
  static constexpr size_t slot_bytes() noexcept { return sizeof(Slot); }
  explicit Forest256Cache(std::span<const uint64_t> stream, size_t budget)
      : stream_(stream),
        k_(budget / slot_bytes()),
        slots_(k_ ? new Slot[k_] : nullptr) {}

  Forest256Cache(const Forest256Cache& rhs)
      : Forest256Cache(rhs.stream_, rhs.allocated_bytes()) {}
  Forest256Cache(Forest256Cache&& rhs) noexcept
      : stream_(std::exchange(rhs.stream_, {})),
        k_(std::exchange(rhs.k_, 0)),
        slots_(std::move(rhs.slots_)) {}
  /// Adopts the source capacity, not the destination's original budget. Copy
  /// assignment allocates the replacement before releasing the old slots, so
  /// both allocations coexist briefly. The byte budget bounds retained slots,
  /// not assignment peak memory. Move assignment allocates no slots.
  Forest256Cache& operator=(Forest256Cache rhs) noexcept {
    std::swap(stream_, rhs.stream_);
    std::swap(k_, rhs.k_);
    slots_.swap(rhs.slots_);
    return *this;
  }

  /// Invalidates all tags even when the new backing address is unchanged.
  void bind(std::span<const uint64_t> stream) noexcept {
    stream_ = stream;
    for (size_t i = 0; i < k_; ++i) {
      slots_[i].valid = false;
    }
  }
  size_t slot_count() const noexcept { return k_; }
  /// Add sizeof(Forest256Cache) for fixed per-reader resident storage.
  size_t allocated_bytes() const noexcept { return k_ * slot_bytes(); }

  /// Optional caller-owned diagnostics; packed bits count the logical input.
  struct Stats {
    size_t hits{}, misses{}, evictions{}, packed_bits{};
  };

  /// Each unit must consistently identify one (codec, offset, count) within
  /// the current binding. local is in [0,count). No CSR reference escapes.
  NodeInfo lookup(size_t unit,
                  Codec codec,
                  size_t offset,
                  size_t count,
                  size_t local) {
    return get<false>(unit, codec, offset, count, local, nullptr);
  }
  NodeInfo lookup_counted(size_t unit,
                          Codec codec,
                          size_t offset,
                          size_t count,
                          size_t local,
                          Stats& stats) {
    return get<true>(unit, codec, offset, count, local, &stats);
  }

  /// Copy a decoded unit into caller-owned walk scratch. The same unit-key
  /// contract as lookup applies. The copy survives eviction and rebinding.
  void snapshot(size_t unit,
                Codec codec,
                size_t offset,
                size_t count,
                LevelCsr256& out) {
    if (k_) {
      const auto& decoded = cached<false>(unit, codec, offset, count, nullptr);
      std::memcpy(&out, &decoded, sizeof(out));
    } else {
      unpack_dense(codec, stream_, offset, count, out);
    }
  }

 private:
  static NodeInfo inspect(const LevelCsr256& csr, uint8_t u) noexcept {
    return {static_cast<uint16_t>(csr.root_count()),
            static_cast<uint16_t>(csr.offset_m1[u] + 1),
            static_cast<uint16_t>(csr.degree(u)), csr.parent_of(u),
            csr.is_last_child(u)};
  }
  template <bool Count>
  NodeInfo get(size_t unit,
               Codec codec,
               size_t offset,
               size_t count,
               size_t local,
               Stats* stats) {
    if (count > 256 || local >= count) {
      throw std::out_of_range("Forest256 local");
    }
    if (k_) {
      return inspect(cached<Count>(unit, codec, offset, count, stats), local);
    }
    LevelCsr256 transient;
    unpack_dense(codec, stream_, offset, count, transient);
    if constexpr (Count) {
      ++stats->misses;
      stats->packed_bits += 2 * count;
    }
    return inspect(transient, local);
  }
  template <bool Count>
  const LevelCsr256& cached(size_t unit,
                            Codec codec,
                            size_t offset,
                            size_t count,
                            Stats* stats) {
    auto& slot = slots_[unit % k_];
    if (slot.valid && slot.unit == unit) {
      if constexpr (Count) {
        ++stats->hits;
      }
    } else {
      // Extent/population checks throw before any output stores.
      unpack_dense(codec, stream_, offset, count, slot.csr);
      if constexpr (Count) {
        stats->evictions += slot.valid;
        ++stats->misses;
        stats->packed_bits += 2 * count;
      }
      slot.unit = unit;
      slot.valid = true;
    }
    return slot.csr;
  }

  std::span<const uint64_t> stream_;
  size_t k_;
  std::unique_ptr<Slot[]> slots_;
};
}  // namespace pixie::tree::detail::forest256
