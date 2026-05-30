#pragma once

#include <pixie/bits.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>

/**
 * Benchmark note:
 *
 * This header keeps experimental and historical excess_min/excess_positions
 * variants for comparison benchmarks. Production excess code lives in
 * pixie/bits.h; a benchmark win here should be treated as a candidate to port,
 * not evidence that callers use this variant already.
 *
 * Diagnostic run, 2026-05-30:
 *   taskset -c 0 build/benchmarks-diagnostic_local/excess_positions_benchmarks
 *     --benchmark_filter='BM_ExcessMin128(_|/|$)'
 *     --benchmark_repetitions=3
 *     --benchmark_perf_counters=CYCLES,INSTRUCTIONS,CACHE-MISSES
 *     --benchmark_counters_tabular=true
 *
 * Production excess_min_128:
 *   range      ns    cycles  instructions  cache misses
 *   0-128     8.1    35.8   117.9         0.001
 *   0-127    11.6    50.9   174.9         0.001
 *   0-16      3.7    16.4    88.5         0.000
 *   0-32      5.6    24.7   124.8         0.000
 *   0-48      9.1    40.5   122.1         0.000
 *   0-64      8.1    36.2   121.4         0.000
 *   0-31     12.1    54.2   177.6         0.000
 *   1-17     14.4    64.6   213.3         0.000
 *   3-35     13.5    60.1   212.4         0.000
 *   5-37     12.9    57.5   215.3         0.000
 *   32-64     7.4    33.1   153.0         0.000
 *   33-65    14.3    63.6   214.6         0.001
 *   64-96     7.8    34.2   148.9         0.000
 *   61-93    16.0    69.4   218.5         0.002
 *   96-128    6.5    28.1   151.9         0.001
 *   56-72     5.2    22.7   116.5         0.000
 *   60-68     8.5    38.1   148.7         0.000
 *   63-64     3.1    13.3    84.0         0.000
 *   17-17     2.1     9.1    47.0         0.000
 *   Random   11.3    49.7   197.0         0.001
 *
 * New non-aligned/random range timings, ns:
 *   range   Prod  Scalar  Nibble  Byte  Hybrid  Expand16  Lane64  Split64  Skip
 *   1-17    14.4    14.6    10.1  12.5    11.9      30.3    13.3     15.9  16.3
 *   3-35    13.5    27.1    12.5  14.2    16.1      46.1    13.1     12.7  14.1
 *   5-37    12.9    27.0    14.6  15.2    18.5      44.7    13.5     13.0  14.7
 *   33-65   14.3    26.6    13.3  17.8    17.4      42.1    14.5     13.3  13.6
 *   61-93   16.0    26.6    15.4  14.8    17.7      43.8    14.0     13.6  11.4
 *   Random  11.3    39.9    18.9  14.7    11.5      68.2    11.9     12.4  12.6
 *
 * Random range hardware counters:
 *   variant       ns    cycles  instructions
 *   Production   11.3    49.7   197.0
 *   ScalarBits   39.9   175.9   721.9
 *   NibbleLUT    18.9    83.2   291.2
 *   ByteLUT      14.7    63.9   261.2
 *   HybridLUT    11.5    50.8   209.8
 *   Expand16     68.2   300.8   879.8
 *   Lane64SSE    11.9    52.1   224.1
 *   Split64SSE   12.4    54.1   236.9
 *   ShortSkip    12.6    55.7   252.5
 *
 * Diagnostic run, 2026-05-30:
 *   taskset -c 0 build/benchmarks-diagnostic_local/excess_positions_benchmarks
 *     --benchmark_filter='BM_ExcessPositions512'
 *     --benchmark_repetitions=3
 *     --benchmark_perf_counters=CYCLES,INSTRUCTIONS,CACHE-MISSES
 *     --benchmark_counters_tabular=true
 *
 * excess_positions_512 random target in [-128, 128]:
 *   variant        ns    cycles  instructions  cache misses
 *   Production    11.4    50.8   188.9         0.001
 *   LUTAVX512     12.8    56.8   195.5         0.002
 *   BranchingLUT  16.7    73.4   261.4         0.003
 *   ExpandAVX512  21.0    93.6   266.6         0.003
 *   Expand8       24.6   109.9   449.7         0.002
 *   Expand        46.8   207.7   784.8         0.006
 *   ByteLUT       49.7   221.4   754.5         0.008
 *   Scalar       374.2  1656.0  7716.6         0.041
 *
 * excess_positions_512 fixed-target timings, ns:
 *   variant       -64   -8    0     8    64
 *   Production    11.6  18.0  18.3  19.1  12.3
 *   LUTAVX512     13.4  17.9  19.2  21.3  13.2
 *   BranchingLUT  19.1  28.9  28.7  28.3  16.8
 *   ExpandAVX512  22.7  36.6  36.3  36.2  22.5
 *   Expand8       17.6  52.7  47.5  46.9  17.7
 *   Expand        51.0  85.9  86.1  85.5  53.9
 *   ByteLUT       34.7  77.4  76.9  79.0  34.3
 *   Scalar       367.2 433.5 466.3 428.1 364.6
 */

namespace pixie::experimental {

namespace detail {

constexpr int8_t nibble_delta(uint8_t x) {
  return static_cast<int8_t>(2 * std::popcount(x) - 4);
}

constexpr int8_t byte_delta(uint8_t x) {
  return static_cast<int8_t>(2 * std::popcount(x) - 8);
}

constexpr int8_t min_prefix(uint8_t x, int bits) {
  int cur = 0;
  int best = 0;
  for (int bit = 0; bit < bits; ++bit) {
    cur += ((x >> bit) & 1u) != 0 ? 1 : -1;
    if (bit == 0 || cur < best) {
      best = cur;
    }
  }
  return static_cast<int8_t>(best);
}

constexpr int8_t min_prefix_offset(uint8_t x, int bits) {
  int cur = 0;
  int best = 0;
  int best_offset = 1;
  for (int bit = 0; bit < bits; ++bit) {
    cur += ((x >> bit) & 1u) != 0 ? 1 : -1;
    if (bit == 0 || cur < best) {
      best = cur;
      best_offset = bit + 1;
    }
  }
  return static_cast<int8_t>(best_offset);
}

template <size_t N, typename Fn>
constexpr std::array<int8_t, N> make_lut(Fn fn) {
  std::array<int8_t, N> out{};
  for (size_t i = 0; i < N; ++i) {
    out[i] = fn(static_cast<uint8_t>(i));
  }
  return out;
}

static inline constexpr std::array<int8_t, 16> kNibbleDelta =
    make_lut<16>([](uint8_t x) { return nibble_delta(x); });
static inline constexpr std::array<int8_t, 16> kNibbleMin =
    make_lut<16>([](uint8_t x) { return min_prefix(x, 4); });
static inline constexpr std::array<int8_t, 16> kNibbleMinOffset =
    make_lut<16>([](uint8_t x) { return min_prefix_offset(x, 4); });
static inline constexpr std::array<int8_t, 256> kByteDelta =
    make_lut<256>([](uint8_t x) { return byte_delta(x); });
static inline constexpr std::array<int8_t, 256> kByteMin =
    make_lut<256>([](uint8_t x) { return min_prefix(x, 8); });
static inline constexpr std::array<int8_t, 256> kByteMinOffset =
    make_lut<256>([](uint8_t x) { return min_prefix_offset(x, 8); });

static inline void scan_bit(const uint64_t* s,
                            size_t bit,
                            int& current,
                            int& best,
                            size_t& best_offset) noexcept {
  current += ((s[bit >> 6] >> (bit & 63)) & 1ull) != 0 ? 1 : -1;
  const size_t offset = bit + 1;
  if (current < best) {
    best = current;
    best_offset = offset;
  }
}

}  // namespace detail

/**
 * @brief Reference scalar excess_min_128 implementation.
 *
 * @details Workflow:
 *
 *   prefix(left) -> scan bits left..right-1 -> first strict minimum
 *
 * The value at offset left is included before scanning any bits, matching the
 * production inclusive prefix range [left, right]. Each scanned bit advances to
 * the next prefix offset. Ties are intentionally ignored, so the first minimum
 * offset is preserved.
 */
static inline ExcessResult excess_min_128_scalar_bits(const uint64_t* s,
                                                      size_t left,
                                                      size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  int current = best;
  for (size_t bit = left; bit < right; ++bit) {
    detail::scan_bit(s, bit, current, best, best_offset);
  }
  return {best, best_offset};
}

/**
 * @brief Scalar 4-bit LUT excess_min_128 experiment.
 *
 * @details Workflow:
 *
 *   unaligned bits -> full nibbles -> trailing bits
 *                      | delta
 *                      | local min
 *                      ` first local min offset
 *
 * Full nibbles use lookup tables for the nibble delta, the minimum prefix value
 * inside positions 1..4, and the first local bit offset that reaches that
 * minimum. Boundary bits are scanned scalar so the LUT never observes bits
 * outside the query range.
 */
static inline ExcessResult excess_min_128_nibble_lut(const uint64_t* s,
                                                     size_t left,
                                                     size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 3u) != 0; ++bit) {
    detail::scan_bit(s, bit, current, best, best_offset);
  }

  for (; bit + 4 <= right; bit += 4) {
    const uint8_t nibble =
        static_cast<uint8_t>((s[bit >> 6] >> (bit & 63)) & 0xFu);
    const int candidate = current + detail::kNibbleMin[nibble];
    if (candidate < best) {
      best = candidate;
      best_offset = bit + static_cast<size_t>(detail::kNibbleMinOffset[nibble]);
    }
    current += detail::kNibbleDelta[nibble];
  }

  for (; bit < right; ++bit) {
    detail::scan_bit(s, bit, current, best, best_offset);
  }
  return {best, best_offset};
}

/**
 * @brief Scalar 8-bit LUT excess_min_128 experiment.
 *
 * @details Workflow:
 *
 *   unaligned bits -> full bytes -> trailing bits
 *                     | delta
 *                     | local min
 *                     ` first local min offset
 *
 * Full bytes use lookup tables for byte delta, minimum prefix value inside
 * positions 1..8, and the first local bit offset that reaches that minimum.
 * This reduces loop iterations on byte-aligned ranges but pays scalar boundary
 * work on unaligned ranges.
 */
static inline ExcessResult excess_min_128_byte_lut(const uint64_t* s,
                                                   size_t left,
                                                   size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 7u) != 0; ++bit) {
    detail::scan_bit(s, bit, current, best, best_offset);
  }

  for (; bit + 8 <= right; bit += 8) {
    const uint8_t byte =
        static_cast<uint8_t>((s[bit >> 6] >> (bit & 63)) & 0xFFu);
    const int candidate = current + detail::kByteMin[byte];
    if (candidate < best) {
      best = candidate;
      best_offset = bit + static_cast<size_t>(detail::kByteMinOffset[byte]);
    }
    current += detail::kByteDelta[byte];
  }

  for (; bit < right; ++bit) {
    detail::scan_bit(s, bit, current, best, best_offset);
  }
  return {best, best_offset};
}

/**
 * @brief Hybrid dispatch over scalar, byte-LUT, nibble-LUT, and production.
 *
 * @details Workflow:
 *
 *   width <= 2                  -> scalar bits
 *   width <= 64 and byte aligned -> byte LUT
 *   width <= 32                 -> nibble LUT
 *   otherwise                   -> production excess_min_128
 *
 * This variant probes whether the fastest implementation depends primarily on
 * query width and boundary alignment. It keeps production behavior for wider
 * ranges where the AVX2 production path usually wins.
 */
static inline ExcessResult excess_min_128_hybrid_lut(const uint64_t* s,
                                                     size_t left,
                                                     size_t right) noexcept {
  if (left > right) {
    return {};
  }
  const size_t clamped_left = std::min<size_t>(left, 128);
  const size_t clamped_right = std::min<size_t>(right, 128);
  const size_t width = clamped_right - clamped_left;

  if (width <= 2) {
    return excess_min_128_scalar_bits(s, left, right);
  }
  if (width <= 64 && (clamped_left & 7u) == 0 && (clamped_right & 7u) == 0) {
    return excess_min_128_byte_lut(s, left, right);
  }
  if (width <= 32) {
    return excess_min_128_nibble_lut(s, left, right);
  }
  return excess_min_128(s, left, right);
}

#ifdef PIXIE_AVX2_SUPPORT
// clang-format off
static inline const __m128i excess_lut_delta_128 = _mm_setr_epi8(
    -4, -2, -2,  0,
    -2,  0,  0,  2,
    -2,  0,  0,  2,
     0,  2,  2,  4);
static inline const __m128i excess_lut_min_128 = _mm_setr_epi8(
    -4, -2, -2,  0,
    -2,  0, -1,  1,
    -3, -1, -1,  1,
    -2,  0, -1,  1);
static inline const __m128i excess_lut_nibble_index_128 = _mm_setr_epi8(
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12, 13, 14, 15);
static inline const __m128i excess_lut_nibble_mask_128 = _mm_set1_epi8(0x0F);
// clang-format on

namespace detail {

static inline __m128i excess_nibbles_64_sse(uint64_t word) noexcept {
  const __m128i word_vec = _mm_cvtsi64_si128(static_cast<int64_t>(word));
  const __m128i lo_nibbles =
      _mm_and_si128(word_vec, excess_lut_nibble_mask_128);
  const __m128i hi_nibbles =
      _mm_and_si128(_mm_srli_epi16(word_vec, 4), excess_lut_nibble_mask_128);
  return _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
}

static inline __m128i excess_prefix_sum_16x_i8(__m128i v) noexcept {
  __m128i x = v;
  __m128i t = _mm_slli_si128(x, 1);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 2);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 4);
  x = _mm_add_epi8(x, t);
  t = _mm_slli_si128(x, 8);
  return _mm_add_epi8(x, t);
}

static inline void scan_full_nibbles_64_sse(uint64_t word,
                                            int lane_base_excess,
                                            size_t lane_base_offset,
                                            size_t first_nibble,
                                            size_t last_nibble,
                                            int& best,
                                            size_t& best_offset) noexcept {
  if (first_nibble >= last_nibble) {
    return;
  }

  const __m128i nibbles = excess_nibbles_64_sse(word);
  __m128i ps =
      excess_prefix_sum_16x_i8(_mm_shuffle_epi8(excess_lut_delta_128, nibbles));
  const __m128i excl_ps = _mm_alignr_epi8(ps, _mm_setzero_si128(), 15);
  const __m128i candidates = _mm_add_epi8(
      _mm_add_epi8(_mm_set1_epi8(static_cast<int8_t>(lane_base_excess)),
                   excl_ps),
      _mm_shuffle_epi8(excess_lut_min_128, nibbles));

  const __m128i idx = excess_lut_nibble_index_128;
  const __m128i first_minus_one =
      _mm_set1_epi8(static_cast<int8_t>(static_cast<int>(first_nibble) - 1));
  const __m128i last = _mm_set1_epi8(static_cast<int8_t>(last_nibble));
  const __m128i active = _mm_and_si128(_mm_cmpgt_epi8(idx, first_minus_one),
                                       _mm_cmpgt_epi8(last, idx));
  const __m128i masked_candidates =
      _mm_blendv_epi8(_mm_set1_epi8(127), candidates, active);

  __m128i min128 = masked_candidates;
  min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 8));
  min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 4));
  min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 2));
  min128 = _mm_min_epi8(min128, _mm_alignr_epi8(min128, min128, 1));

  const int candidate_min =
      static_cast<int>(static_cast<int8_t>(_mm_extract_epi8(min128, 0)));
  if (candidate_min < best) {
    const __m128i equal_min = _mm_cmpeq_epi8(
        masked_candidates, _mm_set1_epi8(static_cast<int8_t>(candidate_min)));
    const uint32_t equal_mask =
        static_cast<uint32_t>(_mm_movemask_epi8(equal_min));
    const uint32_t nibble_index = std::countr_zero(equal_mask);
    const uint8_t nibble =
        static_cast<uint8_t>((word >> (nibble_index * 4u)) & 0xFu);
    best = candidate_min;
    best_offset = lane_base_offset + static_cast<size_t>(nibble_index) * 4u +
                  static_cast<size_t>(kNibbleMinOffset[nibble]);
  }
}

static inline ExcessResult excess_min_128_split64_sse_impl(
    const uint64_t* s,
    size_t left,
    size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  int current = best;
  size_t bit = left;
  for (; bit < right && (bit & 3u) != 0; ++bit) {
    scan_bit(s, bit, current, best, best_offset);
  }

  size_t first_full_nibble = bit >> 2;
  const size_t last_full_nibble = right >> 2;
  while (first_full_nibble < last_full_nibble) {
    const size_t word_index = first_full_nibble >> 4;
    const size_t lane_first = first_full_nibble & 15u;
    const size_t lane_last =
        std::min<size_t>(last_full_nibble - word_index * 16u, 16);
    const size_t lane_base_offset = word_index * 64u;
    scan_full_nibbles_64_sse(
        s[word_index], prefix_excess_128(s, lane_base_offset), lane_base_offset,
        lane_first, lane_last, best, best_offset);
    first_full_nibble = word_index * 16u + lane_last;
  }

  bit = std::max(bit, first_full_nibble * 4u);
  current = prefix_excess_128(s, bit);
  for (; bit < right; ++bit) {
    scan_bit(s, bit, current, best, best_offset);
  }

  return {best, best_offset};
}

}  // namespace detail

/**
 * @brief Single-64-bit-lane SSE excess_min_128 experiment.
 *
 * @details Workflow:
 *
 *   scalar boundary -> one 64-bit word as 16 nibbles -> scalar tail
 *                     fallback to production if full nibbles cross words
 *
 * This variant tests whether short ranges benefit from avoiding the 128-bit
 * cross-lane prefix work used by broader vector paths. It only handles the
 * single-word full-nibble case; multi-word ranges fall back to production.
 */
static inline ExcessResult excess_min_128_lane64_sse(const uint64_t* s,
                                                     size_t left,
                                                     size_t right) noexcept {
  if (left > right) {
    return {};
  }
  const size_t clamped_left = std::min<size_t>(left, 128);
  const size_t clamped_right = std::min<size_t>(right, 128);
  const size_t first_full_nibble = ((clamped_left + 3u) & ~size_t{3}) >> 2;
  const size_t last_full_nibble = clamped_right >> 2;
  if (first_full_nibble < last_full_nibble &&
      (first_full_nibble >> 4) != ((last_full_nibble - 1u) >> 4)) {
    return excess_min_128(s, left, right);
  }
  return detail::excess_min_128_split64_sse_impl(s, left, right);
}

/**
 * @brief Split-64-bit-lane SSE excess_min_128 experiment.
 *
 * @details Workflow:
 *
 *   scalar boundary -> word 0 full nibbles -> word 1 full nibbles -> tail
 *                       16-nibble SSE         16-nibble SSE
 *
 * Each 64-bit word is processed as an independent 16-nibble vector scan. The
 * base excess for a word is recomputed from prefix_excess_128, avoiding
 * vector-prefix carry propagation across the 64-bit boundary.
 */
static inline ExcessResult excess_min_128_split64_sse(const uint64_t* s,
                                                      size_t left,
                                                      size_t right) noexcept {
  return detail::excess_min_128_split64_sse_impl(s, left, right);
}

/**
 * @brief Short-range dispatch experiment.
 *
 * @details Workflow:
 *
 *   width <= 2              -> scalar bits
 *   full nibbles in 1 word  -> lane64 SSE
 *   width <= 80             -> split64 SSE
 *   otherwise               -> production excess_min_128
 *
 * This variant tests two ideas together: avoid 128-bit lane crossing when a
 * query is contained in one 64-bit word, and skip a few production iterations
 * for medium ranges where split-lane scans may be cheaper.
 */
static inline ExcessResult excess_min_128_short_skip(const uint64_t* s,
                                                     size_t left,
                                                     size_t right) noexcept {
  if (left > right) {
    return {};
  }
  const size_t clamped_left = std::min<size_t>(left, 128);
  const size_t clamped_right = std::min<size_t>(right, 128);
  const size_t width = clamped_right - clamped_left;
  if (width <= 2) {
    return excess_min_128_scalar_bits(s, left, right);
  }

  const size_t first_full_nibble = ((clamped_left + 3u) & ~size_t{3}) >> 2;
  const size_t last_full_nibble = clamped_right >> 2;
  if (first_full_nibble < last_full_nibble &&
      (first_full_nibble >> 4) == ((last_full_nibble - 1u) >> 4)) {
    return excess_min_128_lane64_sse(s, left, right);
  }
  if (width <= 80) {
    return excess_min_128_split64_sse(s, left, right);
  }
  return excess_min_128(s, left, right);
}

// clang-format off
static inline const __m256i excess_branch_lut_em4 = _mm256_setr_epi8(
    0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em3 = _mm256_setr_epi8(
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em2 = _mm256_setr_epi8(
    0x02, 0x08, 0x08, 0x00, 0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00,
    0x02, 0x08, 0x08, 0x00, 0x0A, 0x00, 0x00, 0x00,
    0x0A, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00);

static inline const __m256i excess_branch_lut_em1 = _mm256_setr_epi8(
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00,
    0x01, 0x04, 0x05, 0x00, 0x05, 0x00, 0x01, 0x00);

static inline const __m256i excess_branch_lut_e0 = _mm256_setr_epi8(
    0x00, 0x02, 0x02, 0x08, 0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00, 0x08, 0x02, 0x02, 0x00,
    0x00, 0x02, 0x02, 0x08, 0x00, 0x0A, 0x0A, 0x00,
    0x00, 0x0A, 0x0A, 0x00, 0x08, 0x02, 0x02, 0x00);

static inline const __m256i excess_branch_lut_e1 = _mm256_setr_epi8(
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01,
    0x00, 0x01, 0x00, 0x05, 0x00, 0x05, 0x04, 0x01);

static inline const __m256i excess_branch_lut_e2 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A, 0x00, 0x08, 0x08, 0x02,
    0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x0A,
    0x00, 0x00, 0x00, 0x0A, 0x00, 0x08, 0x08, 0x02);

static inline const __m256i excess_branch_lut_e3 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04);

static inline const __m256i excess_branch_lut_e4 = _mm256_setr_epi8(
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08);
// clang-format on

static inline __m256i excess_bit_masks_16x() noexcept {
  return _mm256_setr_epi16(0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020,
                           0x0040, 0x0080, 0x0100, 0x0200, 0x0400, 0x0800,
                           0x1000, 0x2000, 0x4000, (int16_t)0x8000);
}

static inline __m256i excess_prefix_sum_16x_i16(__m256i v) noexcept {
  __m256i x = v;
  __m256i t = _mm256_slli_si256(x, 2);
  x = _mm256_add_epi16(x, t);
  t = _mm256_slli_si256(x, 4);
  x = _mm256_add_epi16(x, t);
  t = _mm256_slli_si256(x, 8);
  x = _mm256_add_epi16(x, t);

  __m128i lo = _mm256_extracti128_si256(x, 0);
  __m128i hi = _mm256_extracti128_si256(x, 1);
  const int16_t carry = (int16_t)_mm_extract_epi16(lo, 7);
  hi = _mm_add_epi16(hi, _mm_set1_epi16(carry));

  __m256i out = _mm256_castsi128_si256(lo);
  out = _mm256_inserti128_si256(out, hi, 1);
  return out;
}

static inline int16_t excess_last_prefix_16x_i16(__m256i pref) noexcept {
  __m128i hi = _mm256_extracti128_si256(pref, 1);
  return (int16_t)_mm_extract_epi16(hi, 7);
}

static inline __m256i excess_bit_masks_32x8() noexcept {
  return _mm256_setr_epi8(0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80,
                          0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (char)0x80);
}

static inline __m256i excess_byte_selectors_32x8() noexcept {
  return _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                          2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
}

static inline __m256i excess_prefix_sum_32x_i8(__m256i v) noexcept {
  __m256i x = v;
  __m256i t = _mm256_slli_si256(x, 1);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 2);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 4);
  x = _mm256_add_epi8(x, t);
  t = _mm256_slli_si256(x, 8);
  x = _mm256_add_epi8(x, t);

  __m128i lo = _mm256_extracti128_si256(x, 0);
  __m128i hi = _mm256_extracti128_si256(x, 1);
  const int8_t carry = (int8_t)_mm_extract_epi8(lo, 15);
  hi = _mm_add_epi8(hi, _mm_set1_epi8(carry));

  __m256i out = _mm256_castsi128_si256(lo);
  out = _mm256_inserti128_si256(out, hi, 1);
  return out;
}

static inline int8_t excess_last_prefix_32x_i8(__m256i pref) noexcept {
  __m128i hi = _mm256_extracti128_si256(pref, 1);
  return (int8_t)_mm_extract_epi8(hi, 15);
}

/**
 * @brief AVX2 expand-to-i16 excess_min_128 experiment.
 *
 * @details Workflow:
 *
 *   16 input bits -> 16 x i16 +/-1 -> vector prefix sum -> store -> scalar min
 *
 * The implementation scans eight 16-bit chunks. For chunks overlapping the
 * query, bits are expanded to signed +/-1 i16 lanes and prefix summed with the
 * running carry. The vector result is stored to memory, then relevant lanes are
 * checked scalar for the first strict minimum.
 */
static inline ExcessResult excess_min_128_expand16_avx2(const uint64_t* s,
                                                        size_t left,
                                                        size_t right) noexcept {
  if (left > right) {
    return {};
  }
  left = std::min<size_t>(left, 128);
  right = std::min<size_t>(right, 128);

  int best = prefix_excess_128(s, left);
  size_t best_offset = left;
  if (left == right) {
    return {best, best_offset};
  }

  const __m256i masks = excess_bit_masks_16x();
  const __m256i zero = _mm256_setzero_si256();
  const __m256i pos = _mm256_set1_epi16(1);
  const __m256i neg = _mm256_set1_epi16(-1);

  int carry = 0;
  alignas(32) int16_t prefix_values[16];
  for (size_t chunk = 0; chunk < 8; ++chunk) {
    const size_t chunk_bit = chunk * 16;
    const uint16_t bits =
        chunk < 4
            ? static_cast<uint16_t>((s[0] >> (chunk * 16)) & 0xFFFFu)
            : static_cast<uint16_t>((s[1] >> ((chunk - 4) * 16)) & 0xFFFFu);
    const int delta = 2 * static_cast<int>(std::popcount(bits)) - 16;

    if (chunk_bit + 1 <= right && chunk_bit + 16 >= left) {
      const __m256i selected = _mm256_and_si256(
          _mm256_set1_epi16(static_cast<int16_t>(bits)), masks);
      const __m256i is_zero = _mm256_cmpeq_epi16(selected, zero);
      const __m256i steps = _mm256_blendv_epi8(pos, neg, is_zero);
      const __m256i pref =
          _mm256_add_epi16(excess_prefix_sum_16x_i16(steps),
                           _mm256_set1_epi16(static_cast<int16_t>(carry)));
      _mm256_store_si256(reinterpret_cast<__m256i*>(prefix_values), pref);

      for (size_t lane = 0; lane < 16; ++lane) {
        const size_t offset = chunk_bit + lane + 1;
        if (offset < left || offset > right) {
          continue;
        }
        const int value = prefix_values[lane];
        if (value < best) {
          best = value;
          best_offset = offset;
        }
      }
    }
    carry += delta;
  }

  return {best, best_offset};
}

/**
 * @brief Historical AVX2 branching-LUT excess_positions_512 variant.
 *
 * @details Workflow:
 *
 *   128-bit block -> 32 nibbles -> prefix sums -> branch by target-relative
 *                                 -> LUT masks -> packed output bits
 *
 * Each 128-bit block is converted to nibbles, prefix-summed, and filtered by
 * reachability. A family of per-target LUTs produces within-nibble match masks,
 * which are packed back to the output words.
 */
static inline void excess_positions_512_branching_lut(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  int cur = 0;
  const __m256i vdelta =
      _mm256_setr_epi8(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4, -4,
                       -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4);
  const __m256i vmult = _mm256_set1_epi16(0x1001);
  const __m128i vnibble_mask = _mm_set1_epi8(0x0F);

  for (int k = 0; k < 4; ++k) {
    __m128i word_vec = _mm_loadu_si128((const __m128i*)&s[2 * k]);
    __m128i lo_nibbles = _mm_and_si128(word_vec, vnibble_mask);
    __m128i hi_nibbles =
        _mm_and_si128(_mm_srli_epi16(word_vec, 4), vnibble_mask);

    __m128i unpack_lo = _mm_unpacklo_epi8(lo_nibbles, hi_nibbles);
    __m128i unpack_hi = _mm_unpackhi_epi8(lo_nibbles, hi_nibbles);
    __m256i nibbles = _mm256_inserti128_si256(_mm256_castsi128_si256(unpack_lo),
                                              unpack_hi, 1);

    __m256i ps = _mm256_shuffle_epi8(vdelta, nibbles);
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 1));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 2));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 4));
    ps = _mm256_add_epi8(ps, _mm256_slli_si256(ps, 8));

    __m128i ps_lo = _mm256_castsi256_si128(ps);
    __m128i ps_hi = _mm256_extracti128_si256(ps, 1);
    __m128i carry = _mm_set1_epi8((int8_t)_mm_extract_epi8(ps_lo, 15));
    ps_hi = _mm_add_epi8(ps_hi, carry);
    ps = _mm256_inserti128_si256(_mm256_castsi128_si256(ps_lo), ps_hi, 1);

    __m256i b = _mm256_permute2x128_si256(ps, ps, 0x08);
    __m256i excl_ps = _mm256_alignr_epi8(ps, b, 15);

    int target_rel = target_x - cur;
    int block_delta =
        2 * (std::popcount(s[2 * k]) + std::popcount(s[2 * k + 1])) - 128;

    const int d = 2 * target_rel - block_delta;
    if (d < -128 || d > 128) {
      cur += block_delta;
      continue;
    }

    if (target_rel == 128 || target_rel == -128) {
      out[2 * k + 1] |= (uint64_t{1} << 63);
      cur += block_delta;
      continue;
    }

    __m256i t = _mm256_sub_epi8(_mm256_set1_epi8((int8_t)target_rel), excl_ps);
    __m256i total_match = _mm256_setzero_si256();
    __m256i t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-4));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em4, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-3));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em3, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-2));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em2, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(-1));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_em1, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(0));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e0, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(1));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e1, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(2));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e2, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(3));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e3, nibbles)));
    t_eq = _mm256_cmpeq_epi8(t, _mm256_set1_epi8(4));
    total_match = _mm256_or_si256(
        total_match,
        _mm256_and_si256(t_eq,
                         _mm256_shuffle_epi8(excess_branch_lut_e4, nibbles)));

    __m256i res = _mm256_maddubs_epi16(total_match, vmult);
    __m128i packed = _mm_packus_epi16(_mm256_castsi256_si128(res),
                                      _mm256_extracti128_si256(res, 1));
    _mm_storeu_si128((__m128i*)&out[2 * k], packed);

    cur += block_delta;
  }
}
#else
/**
 * @brief Scalar fallback for the branching-LUT positions variant.
 *
 * @details Used when AVX2 is not enabled. Delegates to production
 * excess_positions_512 so callers can benchmark the same symbol across build
 * configurations.
 */
static inline void excess_positions_512_branching_lut(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  excess_positions_512(s, target_x, out);
}
#endif

#ifdef PIXIE_AVX512_SUPPORT
static inline __m512i excess_lut_delta_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-4, -2, -2, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 2, 2, 4));
}

static inline __m512i excess_lut_pos0_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1));
}

static inline __m512i excess_lut_pos1_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-2, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 2));
}

static inline __m512i excess_lut_pos2_64x() noexcept {
  return _mm512_broadcast_i32x4(
      _mm_setr_epi8(-3, -1, -1, 1, -1, 1, 1, 3, -3, -1, -1, 1, -1, 1, 1, 3));
}

static inline __m512i excess_bit_masks_64x8() noexcept {
  alignas(64) static constexpr int8_t masks[64] = {
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80,
      0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, (int8_t)0x80};
  return _mm512_load_si512((const void*)masks);
}

static inline __m512i excess_byte_selectors_64x8() noexcept {
  alignas(64) static constexpr int8_t selectors[64] = {
      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
      2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5,
      5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7};
  return _mm512_load_si512((const void*)selectors);
}

static inline __m512i excess_prefix_sum_64x_i8(__m512i v) noexcept {
  __m512i x = v;
  __m512i t = _mm512_bslli_epi128(x, 1);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 2);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 4);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 8);
  x = _mm512_add_epi8(x, t);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(x, last_byte);
  const __m512i shift1_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11);
  const __m512i shift2_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7);
  const __m512i shift3_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3);

  __m512i lane_base =
      _mm512_maskz_permutexvar_epi32(0xFFF0, shift1_idx, lane_carry);
  lane_base = _mm512_add_epi8(lane_base, _mm512_maskz_permutexvar_epi32(
                                             0xFF00, shift2_idx, lane_carry));
  lane_base = _mm512_add_epi8(lane_base, _mm512_maskz_permutexvar_epi32(
                                             0xF000, shift3_idx, lane_carry));
  return _mm512_add_epi8(x, lane_base);
}

static inline __m512i excess_prefix_sum_2x32_i8(__m512i v) noexcept {
  __m512i x = v;
  __m512i t = _mm512_bslli_epi128(x, 1);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 2);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 4);
  x = _mm512_add_epi8(x, t);
  t = _mm512_bslli_epi128(x, 8);
  x = _mm512_add_epi8(x, t);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(x, last_byte);
  const __m512i prev_lane_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8);
  const __mmask16 carry_to_second_lane_of_each_half = 0xF0F0;
  const __m512i lane_base = _mm512_maskz_permutexvar_epi32(
      carry_to_second_lane_of_each_half, prev_lane_idx, lane_carry);
  return _mm512_add_epi8(x, lane_base);
}

static inline __m512i excess_nibbles_64x_from_256(__m256i words) noexcept {
  const __m512i bytes16 = _mm512_cvtepu8_epi16(words);
  const __m512i low = _mm512_and_si512(bytes16, _mm512_set1_epi16(0x000F));
  const __m512i high = _mm512_and_si512(_mm512_srli_epi16(bytes16, 4),
                                        _mm512_set1_epi16(0x000F));
  return _mm512_or_si512(low, _mm512_slli_epi16(high, 8));
}

static inline __m512i excess_exclusive_prefix_2x32_i8(__m512i pref) noexcept {
  const __m512i zero = _mm512_setzero_si512();
  __m512i out = _mm512_alignr_epi8(pref, zero, 15);

  const __m512i last_byte = _mm512_set1_epi8(15);
  const __m512i lane_carry = _mm512_shuffle_epi8(pref, last_byte);
  const __m512i prev_lane_idx =
      _mm512_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8);
  const __m512i carry_dwords =
      _mm512_permutexvar_epi32(prev_lane_idx, lane_carry);
  const __mmask64 first_byte_of_second_lane_in_each_half =
      (uint64_t{1} << 16) | (uint64_t{1} << 48);
  return _mm512_or_si512(
      out, _mm512_maskz_mov_epi8(first_byte_of_second_lane_in_each_half,
                                 carry_dwords));
}

static inline uint64_t excess_repeat_byte(int value) noexcept {
  return uint64_t{0x0101010101010101} *
         static_cast<uint8_t>(static_cast<int8_t>(value));
}

/**
 * @brief AVX-512 nibble-LUT excess_positions_512 experiment.
 *
 * @details Workflow:
 *
 *   256 input bits -> 64 nibbles -> two 128-bit logical halves
 *                  -> nibble prefix/LUT matches -> packed output masks
 *
 * The implementation processes four words at a time. It computes reachability
 * for each 128-bit half, converts bytes to nibbles, builds exclusive prefix
 * sums, compares nibble-local positions against the target, and packs matches
 * back to four output words.
 */
static inline void excess_positions_512_lut_avx512(const uint64_t* s,
                                                   int target_x,
                                                   uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  static const __m512i vdelta = excess_lut_delta_64x();
  static const __m512i vpos0 = excess_lut_pos0_64x();
  static const __m512i vpos1 = excess_lut_pos1_64x();
  static const __m512i vpos2 = excess_lut_pos2_64x();
  static const __m512i vbit0 = _mm512_set1_epi8(1);
  static const __m512i vbit1 = _mm512_set1_epi8(2);
  static const __m512i vbit2 = _mm512_set1_epi8(4);
  static const __m512i vbit3 = _mm512_set1_epi8(8);
  static const __m512i vmult = _mm512_set1_epi16(0x1001);

  for (int k = 0; k < 2; ++k) {
    const int base_word = 4 * k;
    const int delta0 =
        2 * (std::popcount(s[base_word]) + std::popcount(s[base_word + 1])) -
        128;
    const int delta1 = 2 * (std::popcount(s[base_word + 2]) +
                            std::popcount(s[base_word + 3])) -
                       128;
    const int target0 = target_x;
    const int target1 = target_x - delta0;
    const bool reachable0 = [&] {
      const int d = 2 * target0 - delta0;
      return -128 <= d && d <= 128;
    }();
    const bool reachable1 = [&] {
      const int d = 2 * target1 - delta1;
      return -128 <= d && d <= 128;
    }();

    if (!reachable0 && !reachable1) {
      target_x -= delta0 + delta1;
      continue;
    }

    const __m256i words =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&s[base_word]));
    const __m512i nibbles = excess_nibbles_64x_from_256(words);
    const __m512i ps =
        excess_prefix_sum_2x32_i8(_mm512_shuffle_epi8(vdelta, nibbles));
    const __m512i excl_ps = excess_exclusive_prefix_2x32_i8(ps);
    const uint64_t repeated0 = excess_repeat_byte(target0);
    const uint64_t repeated1 = excess_repeat_byte(target1);
    const __m512i vtgt =
        _mm512_setr_epi64(repeated0, repeated0, repeated0, repeated0, repeated1,
                          repeated1, repeated1, repeated1);
    const __m512i t = _mm512_sub_epi8(vtgt, excl_ps);

    const __mmask64 cmp0 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos0, nibbles), t);
    const __mmask64 cmp1 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos1, nibbles), t);
    const __mmask64 cmp2 =
        _mm512_cmpeq_epi8_mask(_mm512_shuffle_epi8(vpos2, nibbles), t);
    const __mmask64 cmp3 = _mm512_cmpeq_epi8_mask(ps, vtgt);
    __m512i total_match = _mm512_maskz_mov_epi8(cmp0, vbit0);
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp1, vbit1));
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp2, vbit2));
    total_match =
        _mm512_or_si512(total_match, _mm512_maskz_mov_epi8(cmp3, vbit3));

    const __mmask64 active =
        (reachable0 ? __mmask64{0x00000000FFFFFFFFull} : __mmask64{0}) |
        (reachable1 ? __mmask64{0xFFFFFFFF00000000ull} : __mmask64{0});
    total_match = _mm512_maskz_mov_epi8(active, total_match);

    const __m512i res = _mm512_maddubs_epi16(total_match, vmult);
    const __m256i packed = _mm512_cvtepi16_epi8(res);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(&out[base_word]), packed);

    target_x -= delta0 + delta1;
  }
}
#else
/**
 * @brief Fallback for the AVX-512 LUT positions variant.
 *
 * @details Used when AVX-512 is not enabled. Delegates to production
 * excess_positions_512 so callers can benchmark the same symbol across build
 * configurations.
 */
static inline void excess_positions_512_lut_avx512(const uint64_t* s,
                                                   int target_x,
                                                   uint64_t* out) noexcept {
  excess_positions_512(s, target_x, out);
}
#endif

/**
 * @brief Expand-to-i16 excess_positions_512 experiment.
 *
 * @details Workflow:
 *
 *   16 input bits -> 16 x i16 +/-1 -> vector prefix sum -> compare target
 *                 -> pext mask -> output word
 *
 * With AVX2, this variant scans 16-bit chunks, expands them to i16 prefix
 * lanes, compares absolute prefix values against the target, and compresses the
 * comparison mask into output bits. Without AVX2 it uses a scalar scan.
 */
static inline void excess_positions_512_expand(const uint64_t* s,
                                               int target_x,
                                               uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX2_SUPPORT
  static const __m256i masks = excess_bit_masks_16x();
  static const __m256i vzero = _mm256_setzero_si256();
  static const __m256i vallones = _mm256_cmpeq_epi16(vzero, vzero);
  static const __m256i vminus1 = _mm256_set1_epi16(-1);
  static const __m256i vtwo = _mm256_set1_epi16(2);
  const __m256i vtarget = _mm256_set1_epi16((int16_t)target_x);

  int cur = 0;
  for (int block = 0; block < 4; ++block) {
    const int target_rel = target_x - cur;
    if (target_rel <= -64 || target_rel >= 64) {
      const int block_delta =
          2 * (std::popcount(s[2 * block]) + std::popcount(s[2 * block + 1])) -
          128;
      const int reachability = 2 * target_rel - block_delta;
      if (reachability < -128 || reachability > 128) {
        cur += block_delta;
        continue;
      }
    }

    for (int j = 0; j < 8; ++j) {
      const int k = 8 * block + j;
      const size_t word_idx = size_t(k) >> 2;
      const size_t shift = size_t(k & 3) * 16;
      const uint16_t bits16 =
          static_cast<uint16_t>((s[word_idx] >> shift) & 0xFFFFull);

      const __m256i vb = _mm256_set1_epi16((int16_t)bits16);
      const __m256i m = _mm256_and_si256(vb, masks);
      const __m256i is_zero = _mm256_cmpeq_epi16(m, vzero);
      const __m256i is_set = _mm256_andnot_si256(is_zero, vallones);
      const __m256i steps =
          _mm256_add_epi16(vminus1, _mm256_and_si256(is_set, vtwo));

      const __m256i pref_rel = excess_prefix_sum_16x_i16(steps);
      const __m256i base = _mm256_set1_epi16((int16_t)cur);
      const __m256i pref_abs = _mm256_add_epi16(pref_rel, base);
      const __m256i cmp = _mm256_cmpeq_epi16(pref_abs, vtarget);

      const uint32_t m32 = (uint32_t)_mm256_movemask_epi8(cmp);
      const uint16_t m16 = (uint16_t)_pext_u32(m32, 0xAAAAAAAAu);

      out[word_idx] |= uint64_t(m16) << shift;
      cur += (int)excess_last_prefix_16x_i16(pref_rel);
    }
  }
#else
  int cur = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
#endif
}

/**
 * @brief AVX2 expand-to-i8 excess_positions_512 experiment.
 *
 * @details Workflow:
 *
 *   32 input bits -> 32 x i8 +/-1 -> byte prefix sum -> compare target
 *                 -> movemask -> output word
 *
 * This variant handles 32-bit chunks. It first checks whether the target is
 * reachable within the chunk, then expands bits to byte lanes and uses the
 * vector comparison mask directly as output bits.
 */
static inline void excess_positions_512_expand8(const uint64_t* s,
                                                int target_x,
                                                uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX2_SUPPORT
  static const __m256i byte_selectors = excess_byte_selectors_32x8();
  static const __m256i masks = excess_bit_masks_32x8();
  static const __m256i vzero = _mm256_setzero_si256();
  static const __m256i vallones = _mm256_cmpeq_epi8(vzero, vzero);
  static const __m256i vminus1 = _mm256_set1_epi8(-1);
  static const __m256i vtwo = _mm256_set1_epi8(2);

  int cur = 0;
  for (int k = 0; k < 16; ++k) {
    const size_t word_idx = size_t(k) >> 1;
    const size_t shift = size_t(k & 1) * 32;
    const uint32_t bits32 =
        static_cast<uint32_t>((s[word_idx] >> shift) & 0xFFFFFFFFull);

    const int target_rel = target_x - cur;
    if (target_rel < -32 || target_rel > 32) {
      cur += 2 * static_cast<int>(std::popcount(bits32)) - 32;
      continue;
    }

    const __m256i src = _mm256_set1_epi32((int)bits32);
    const __m256i bytes = _mm256_shuffle_epi8(src, byte_selectors);
    const __m256i m = _mm256_and_si256(bytes, masks);
    const __m256i is_zero = _mm256_cmpeq_epi8(m, vzero);
    const __m256i is_set = _mm256_andnot_si256(is_zero, vallones);
    const __m256i steps =
        _mm256_add_epi8(vminus1, _mm256_and_si256(is_set, vtwo));

    const __m256i pref_rel = excess_prefix_sum_32x_i8(steps);
    const __m256i vtarget = _mm256_set1_epi8((int8_t)target_rel);
    const __m256i cmp = _mm256_cmpeq_epi8(pref_rel, vtarget);
    const uint32_t mask = static_cast<uint32_t>(_mm256_movemask_epi8(cmp));

    out[word_idx] |= uint64_t(mask) << shift;
    cur += static_cast<int>(excess_last_prefix_32x_i8(pref_rel));
  }
#else
  int cur = 0;
  for (size_t i = 0; i < 512; ++i) {
    const uint64_t w = s[i >> 6];
    const int bit = int((w >> (i & 63)) & 1ull);
    cur += bit ? +1 : -1;
    if (cur == target_x) {
      out[i >> 6] |= (uint64_t{1} << (i & 63));
    }
  }
#endif
}

/**
 * @brief AVX-512 expand-to-i8 excess_positions_512 experiment.
 *
 * @details Workflow:
 *
 *   64 input bits -> 64 x i8 +/-1 -> byte prefix sum -> k-mask output
 *
 * This variant processes one 64-bit word per vector. If the target is
 * unreachable in the word, it advances by popcount delta only. Otherwise it
 * expands the word to byte lanes, prefix-sums the lanes, compares against the
 * target, and stores the resulting AVX-512 mask as the output word.
 */
static inline void excess_positions_512_expand_avx512(const uint64_t* s,
                                                      int target_x,
                                                      uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

#ifdef PIXIE_AVX512_SUPPORT
  static const __m512i byte_selectors = excess_byte_selectors_64x8();
  static const __m512i masks = excess_bit_masks_64x8();
  static const __m512i vzero = _mm512_setzero_si512();
  static const __m512i vallones = _mm512_set1_epi8(-1);
  static const __m512i vminus1 = _mm512_set1_epi8(-1);
  static const __m512i vtwo = _mm512_set1_epi8(2);

  int cur = 0;
  for (int k = 0; k < 8; ++k) {
    const uint64_t bits64 = s[k];
    const int target_rel = target_x - cur;
    if (target_rel < -64 || target_rel > 64) {
      cur += 2 * static_cast<int>(std::popcount(bits64)) - 64;
      continue;
    }

    const __m512i src = _mm512_set1_epi64(static_cast<int64_t>(bits64));
    const __m512i bytes = _mm512_shuffle_epi8(src, byte_selectors);
    const __m512i m = _mm512_and_si512(bytes, masks);
    const __mmask64 is_zero = _mm512_cmpeq_epi8_mask(m, vzero);
    const __m512i is_set = _mm512_maskz_mov_epi8(~is_zero, vallones);
    const __m512i steps =
        _mm512_add_epi8(vminus1, _mm512_and_si512(is_set, vtwo));

    const __m512i pref_rel = excess_prefix_sum_64x_i8(steps);
    const __mmask64 match =
        _mm512_cmpeq_epi8_mask(pref_rel, _mm512_set1_epi8((int8_t)target_rel));
    out[k] = static_cast<uint64_t>(match);
    cur += 2 * static_cast<int>(std::popcount(bits64)) - 64;
  }
#else
  excess_positions_512_expand8(s, target_x, out);
#endif
}

struct ExcessByteLut {
  uint8_t masks[256][17];  // target index: T + 8
  int8_t deltas[256];

  constexpr ExcessByteLut() : masks{}, deltas{} {
    for (int b = 0; b < 256; ++b) {
      int pop = 0;
      for (int i = 0; i < 8; ++i) {
        if ((b >> i) & 1) {
          pop++;
        }
      }
      deltas[b] = static_cast<int8_t>(2 * pop - 8);

      for (int t = -8; t <= 8; ++t) {
        uint8_t mask = 0;
        int cur_pop = 0;
        for (int i = 0; i < 8; ++i) {
          if ((b >> i) & 1) {
            cur_pop++;
          }
          int excess = 2 * cur_pop - (i + 1);
          if (excess == t) {
            mask |= (1 << i);
          }
        }
        masks[b][t + 8] = mask;
      }
    }
  }
};

inline constexpr ExcessByteLut kExcessByteLut;

/**
 * @brief Scalar byte-LUT excess_positions_512 experiment.
 *
 * @details Workflow:
 *
 *   byte -> relative target in [-8, 8] -> LUT match mask
 *        -> byte delta -> next byte base excess
 *
 * The table stores, for each byte, all bit positions that reach each local
 * target and the byte delta. The scan walks 64 bytes, emits a mask when the
 * relative target is in range, and then advances the running excess by the
 * byte delta.
 */
static inline void excess_positions_512_byte_lut(const uint64_t* s,
                                                 int target_x,
                                                 uint64_t* out) noexcept {
  out[0] = out[1] = out[2] = out[3] = 0;
  out[4] = out[5] = out[6] = out[7] = 0;

  if (target_x < -512 || target_x > 512) {
    return;
  }

  const uint8_t* bytes = reinterpret_cast<const uint8_t*>(s);
  uint8_t* out_bytes = reinterpret_cast<uint8_t*>(out);

  int cur = 0;
  for (int i = 0; i < 64; ++i) {
    const uint8_t b = bytes[i];
    const int target_rel = target_x - cur;
    if (target_rel >= -8 && target_rel <= 8) {
      out_bytes[i] = kExcessByteLut.masks[b][target_rel + 8];
    }
    cur += kExcessByteLut.deltas[b];
  }
}

}  // namespace pixie::experimental
