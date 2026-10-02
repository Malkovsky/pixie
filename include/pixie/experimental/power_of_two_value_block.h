#pragma once

#include <pixie/permutations/detail/packed_value_block.h>

#include <bit>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace pixie::experimental {

/**
 * @brief Packed unsigned leaf with a power-of-two element capacity.
 * @details Uses the shared owning packed block operations. Metadata and
 * cache-line padding are additional to the payload; for 256 uint64_t values
 * the complete block occupies 2112 bytes. Indexed reads return values, copies
 * own independent storage, and local mutations follow the internal block
 * preconditions. This is an experimental SequenceTree block, not a public
 * container.
 * @tparam T Unsigned integral value type, at most 64 bits.
 * @tparam Values Power-of-two capacity whose payload is a multiple of 64 bits.
 */
template <class T = std::uint64_t, std::size_t Values = 256>
class PowerOfTwoValueBlock
    : public permutations::detail::PackedValueBlock<
          T,
          sizeof(permutations::detail::BitBlock<
                 Values * std::numeric_limits<T>::digits>) *
              CHAR_BIT,
          std::numeric_limits<T>::digits,
          permutations::detail::BitBlock<Values *
                                         std::numeric_limits<T>::digits>> {
  static_assert(std::has_single_bit(Values));
  using Base = permutations::detail::PackedValueBlock<
      T,
      sizeof(permutations::detail::BitBlock<Values *
                                            std::numeric_limits<T>::digits>) *
          CHAR_BIT,
      std::numeric_limits<T>::digits,
      permutations::detail::BitBlock<Values * std::numeric_limits<T>::digits>>;

 public:
  using Base::Base;
  /** @brief Construct empty without initializing unused payload capacity. */
  PowerOfTwoValueBlock() noexcept {}
};

}  // namespace pixie::experimental
