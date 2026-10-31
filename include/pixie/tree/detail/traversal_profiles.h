#pragma once

#include <cstddef>
#include <tuple>

namespace pixie::detail::tree_traversal {
template <size_t Bits>
struct ParenthesisProfile {
  static constexpr size_t bits = Bits;
};
template <size_t Depth, size_t Children>
struct LoudsProfile {
  static constexpr size_t depth = Depth;
  static constexpr size_t children = Children;
};

// Ordered from most to least scratch. Both budget selection and execution
// dispatch use this list; stop at the first profile accepted by the callback.
template <class... Profiles>
struct ProfileList {
  static constexpr size_t size = sizeof...(Profiles);
  template <size_t Index>
  using At = std::tuple_element_t<Index, std::tuple<Profiles...>>;
  template <class Accept>
  static bool visit_until(Accept&& accept) {
    return (accept(Profiles{}) || ...);
  }
};
using ParenthesisProfiles = ProfileList<ParenthesisProfile<8192>,
                                        ParenthesisProfile<2048>,
                                        ParenthesisProfile<512>,
                                        ParenthesisProfile<128>>;
using LoudsProfiles = ProfileList<LoudsProfile<128, 32>,
                                  LoudsProfile<64, 16>,
                                  LoudsProfile<32, 16>,
                                  LoudsProfile<16, 8>,
                                  LoudsProfile<8, 8>>;
}  // namespace pixie::detail::tree_traversal
