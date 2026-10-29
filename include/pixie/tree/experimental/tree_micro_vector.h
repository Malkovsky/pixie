#pragma once
// Internal compiled-event helper for experimental BP traversal. The vector
// argument keeps full root coordinates in SIMD registers on supported builds.
// Public callbacks continue to observe TreeNode; scalar fallback is local.
// Registered 10M inspected DFS, 2026-10-26: 76.26 ms versus 78.86 ms compact
// scalar-argument control (Release native, CPU 0, seven interleaved repeats).
// See profiling/bp-vector-registered-2026-10-26/REPORT.md for scope and limits.
#include <pixie/tree/experimental/tree_micro_prefetch.h>
#if defined(PIXIE_AVX2_SUPPORT) || defined(PIXIE_AVX512_SUPPORT)
#include <immintrin.h>
#endif

namespace pixie::experimental::vector_micro {
using namespace pixie::experimental::detail;
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
template <size_t Id, unsigned Policy, class Visit, class Inspect>
[[gnu::noinline]] void program(Pair root, Visit& visit, Inspect& inspect) {
  constexpr auto n = bp_micro_tables.programs[Id].nodes;
  if constexpr (Policy == 1 || (Policy == 2 && n <= 3) ||
                (Policy == 3 && n >= 4)) {
    events<Id>(root, visit, inspect,
               std::make_index_sequence<bp_micro_tables.programs[Id].count>{});
  } else {
    RunMicroEvents<Id>(
        offset<0, 0>(root), visit, inspect,
        std::make_index_sequence<bp_micro_tables.programs[Id].count>{});
  }
}
template <unsigned Policy, class Visit, class Inspect, size_t... I>
constexpr auto dispatch(std::index_sequence<I...>) {
  return std::array{&program<I, Policy, Visit, Inspect>...};
}
}  // namespace pixie::experimental::vector_micro
