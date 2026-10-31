#pragma once
// Compiled BP events carry both full-width root coordinates in a SIMD pair.
// Scalar fallback is local; callbacks always receive ordinary TreeNode handles.
#include <pixie/tree/detail/bp_micro_tables.h>
#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
#include <immintrin.h>
#endif

namespace pixie::detail::bp_micro::vector_micro {
struct Pair {
#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
  __m128i lanes;
#else
  uint64_t number, position;
#endif
};
inline Pair pack(TreeNode root) {
#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
  return {_mm_set_epi64x(root.pos, root.number)};
#else
  return {root.number, root.pos};
#endif
}
template <unsigned Number, unsigned Position>
[[gnu::always_inline]] inline TreeNode offset(Pair root) {
#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
  const auto value =
      _mm_add_epi64(root.lanes, _mm_set_epi64x(Position, Number));
  uint64_t fields[2];
  _mm_storeu_si128(reinterpret_cast<__m128i*>(fields), value);
  return {fields[0], fields[1]};
#else
  return {root.number + Number, root.position + Position};
#endif
}
template <size_t Id, class Visit, class Inspect, size_t... I>
[[gnu::always_inline]] inline void events(Pair root,
                                          Visit& visit,
                                          Inspect& inspect,
                                          std::index_sequence<I...>) {
  auto event = [&]<size_t J>() {
    constexpr auto code = bp_micro_tables.programs[Id].events[J];
    const auto node = offset<code & 7, (code >> 3) & 15>(root);
    if constexpr (code & 128) {
      inspect(node);
    } else {
      pixie::detail::tree_unpack::emit(
          visit, node,
          !(bp_micro_tables.programs[Id].returns & (uint32_t{1} << J)));
    }
  };
  (event.template operator()<I>(), ...);
}
template <size_t Id, class Visit, class Inspect>
[[gnu::noinline]] void program(Pair root, Visit& visit, Inspect& inspect) {
  {
    events<Id>(root, visit, inspect,
               std::make_index_sequence<bp_micro_tables.programs[Id].count>{});
  }
}
template <class Visit, class Inspect, size_t... I>
constexpr auto dispatch(std::index_sequence<I...>) {
  return std::array{&program<I, Visit, Inspect>...};
}
}  // namespace pixie::detail::bp_micro::vector_micro
