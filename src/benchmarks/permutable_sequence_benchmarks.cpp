#include <pixie/permutations/sequence_implementations.h>

#include "sequence_facade_benchmarks.h"

namespace {
using namespace pixie;
template <typename T,
          ElementStorage Storage = ElementStorage::automatic,
          std::size_t Bits = 2048,
          std::size_t F = 8,
          LengthLayout Layout = LengthLayout::cumulative,
          std::size_t ChunkBytes = 4096>
struct ElementVariant : Configuration<T, Bits, F, Layout> {
  using sequence_type =
      PermutableSequence<T, Storage, Bits, F, Layout, ChunkBytes>;
  static constexpr bool rebased = false;
  static constexpr bool indirect =
      Storage == ElementStorage::indirect || !std::is_unsigned_v<T>;
  static constexpr std::size_t chunk_bytes = indirect ? ChunkBytes : 0;
  static auto Make(std::size_t n) {
    // Prvalue elements include move-only values. No full input/pointer vector.
    auto input =
        std::views::iota(std::size_t{0}, n) |
        std::views::transform([](std::size_t i) { return Value<T>(i); });
    return sequence_type::from_range(input);
  }
};

const bool registered = [] {
  Register<ElementVariant<std::uint64_t, ElementStorage::packed>, true>(
      "Packed64_B256_F8_Prefix");
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect>, true>(
      "Indirect64_B256_F8_Prefix_C4096");
  Register<ElementVariant<bool, ElementStorage::packed>>(
      "PackedBool_B256_F8_Prefix");
  Register<ElementVariant<bool, ElementStorage::indirect>, true>(
      "IndirectBool_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<8>>>("Payload8_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<64>>, true>("Payload64_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<256>>>("Payload256_B256_F8_Prefix_C4096");
  Register<ElementVariant<MoveOnlyPayload>, true>(
      "MoveOnly64_B256_F8_Prefix_C4096");

  Register<UntaggedVariant<std::uint16_t>>("Raw16_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint32_t>>("Raw32_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint64_t>, true>("Raw64_B256_F8_Prefix");
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect, 2048, 8,
                          LengthLayout::cumulative, 1024>,
           false, true>("Indirect64_B256_F8_Prefix_C1024");
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect, 2048, 8,
                          LengthLayout::cumulative, 16384>,
           false, true>("Indirect64_B256_F8_Prefix_C16384");

  return true;
}();
}  // namespace
