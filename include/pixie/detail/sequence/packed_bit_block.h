#pragma once

#include <pixie/detail/sequence/bit_block.h>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

/**
 * @brief Bit block whose complete aligned object fits a physical-bit budget.
 * @details Two uint64_t fields hold valid length and circular origin; all
 * remaining whole words are payload. There are no hidden tree fields or dynamic
 * allocations. On an eight-bit-byte platform PackedBitBlock<2048> is exactly
 * 256 bytes and holds 1920 bits. Local contracts are inherited from BitBlock.
 * @tparam StorageBits Complete object budget, a positive CacheLine multiple.
 */
template <std::size_t StorageBits = 2048>
class PackedBitBlock : public BitBlock<StorageBits - 128> {
  static_assert(StorageBits >= kAlignedStorageLineBits &&
                StorageBits % kAlignedStorageLineBits == 0);
  static_assert(sizeof(BitBlock<StorageBits - 128>) == StorageBits / 8);

 public:
  /** @brief Construct empty without touching unused payload storage. */
  PackedBitBlock() noexcept {}
  using BitBlock<StorageBits - 128>::BitBlock;
};

}  // namespace pixie::detail::sequence
/// @endcond
