#include <pixie/permutations/permutation_implementations.h>

#include "sequence_facade_benchmarks.h"

namespace {
using namespace pixie;
template <typename T,
          std::size_t Bits = 2048,
          std::size_t F = 8,
          LengthLayout Layout = LengthLayout::cumulative>
struct PermutationVariant : Configuration<T, Bits, F, Layout> {
  using sequence_type = Permutation<T, Bits, F, Layout>;
  static constexpr bool rebased = true;
  static constexpr bool indirect = false;
  static constexpr std::size_t chunk_bytes = 0;
  static auto Make(std::size_t n) { return sequence_type::identity(n); }
};

const bool registered = [] {
  Register<PermutationVariant<std::uint16_t>, true>("Perm16_B256_F8_Prefix");
  Register<PermutationVariant<std::uint32_t>>("Perm32_B256_F8_Prefix");
  Register<PermutationVariant<std::uint64_t>, true>("Perm64_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint16_t>>("Raw16_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint32_t>>("Raw32_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint64_t>, true>("Raw64_B256_F8_Prefix");
  // One-axis controls; pair with Perm64_B256_F8_Prefix at the same N.
  Register<PermutationVariant<std::uint64_t, 2048, 4>, false, true>(
      "Perm64_B256_F4_Prefix");
  Register<PermutationVariant<std::uint64_t, 2048, 16>, false, true>(
      "Perm64_B256_F16_Prefix");
  Register<PermutationVariant<std::uint64_t, 2048, 8, LengthLayout::individual>,
           false, true>("Perm64_B256_F8_Raw");
  Register<PermutationVariant<std::uint64_t, 2048, 4, LengthLayout::individual>,
           false, true>("Perm64_B256_F4_Raw");
  Register<
      PermutationVariant<std::uint64_t, 2048, 16, LengthLayout::individual>,
      false, true>("Perm64_B256_F16_Raw");
  Register<PermutationVariant<std::uint64_t, 1024>, false, true>(
      "Perm64_B128_F8_Prefix");
  Register<PermutationVariant<std::uint64_t, 4096>, false, true>(
      "Perm64_B512_F8_Prefix");
  return true;
}();
}  // namespace
