#pragma once

// clang-format off
/**
 * @file permutation_implementations.h
 * @brief Permutation implementation catalog for benchmarks.
 * @details Consumers include pixie/permutations/permutation.h directly.
 *
 * Current snapshot, 2026-09-21: Ryzen 7 8845HS, WSL/Linux, GCC 13.3,
 * -O3 -DNDEBUG -march=native, Google Benchmark 1.9.4, CPU 0, 16 MiB L3.
 * uint64_t, B256/F8/cumulative. Ranges span two current-code CPU medians
 * (5 repetitions at 0.15s, 9 at 0.2s), not confidence intervals. Host drift
 * precludes close rankings or a claim of zero promotion overhead.
 *
 * | Representation | Read 1 MiB ns | Read 32 MiB ns | Read 128 MiB ns | Owned at 32 MiB |
 * | -------------- | ------------- | -------------- | --------------- | --------------- |
 * | Perm64         |            46 |        254-269 |         412-466 |       46.32 MiB |
 * | Raw64 control  |         43-44 |        197-207 |         335-351 |       36.57 MiB |
 *
 * | Representation | Global 32 MiB us        | Local65 32 MiB us       |
 * | -------------- | ----------------------- | ----------------------- |
 * | Perm64         |               37.8-38.4 |               25.8-26.4 |
 * | Raw64 control  |               32.0-32.5 |               23.7-24.3 |
 *
 * Useful MiB describes index payload, not total footprint. Reads include hash
 * generation, 64 reads/iteration. Rotations use 128 fixed-stream operations,
 * 8 batches/repetition; reconstruction/destruction excluded, reconstructed data
 * pre-touched. Local65 means elements and can cross leaves, not a hot kernel.
 * Equal rebased merge at combined N=4194304: 18.3-19.6 us/merge, two
 * independent merges/iteration, 32 iterations; rotated input setup/destruction
 * excluded. Short timed mutation totals remain sensitive to host/harness
 * variation.
 *
 * Fresh-tree requested storage includes owner, capacity, metadata, tags, and
 * padding, not allocator overhead, RSS, or temporary peaks. Default tagged
 * leaf/node allocations are 320/192 bytes versus untagged 256/128; owner is
 * 24 bytes. The 9.75 MiB tag/alignment cost at 32 MiB buys lazy donor rebasing.
 * Defaults remain conservative, configurable, and not claimed optimal.
 *
 * Reproduce each filter separately, using the local /bench wrapper:
 * @code{.sh}
 * flock /tmp/kilo/pixie-experiment-timing.lock \
 *   env CPP_BENCH_CPU=0 CPP_JOBS=4 \
 *   /home/user/.config/kilo/scripts/cpp-bench native \
 *   permutation_benchmarks '<filter>' 9 0.2s
 * @endcode
 * Read filter:
 * `^FacadeAccessIndependent/(Perm64|Raw64)_B256_F8_Prefix/N:(131072|4194304|16777216)$`
 * Rotation filter:
 * `^FacadeRotate(Global|Local65)Batch/(Perm64|Raw64)_B256_F8_Prefix/N:(131072|4194304)/iterations:8$`
 * Merge filter:
 * `^FacadeMergeEqual/Perm64_B256_F8_Prefix/N:4194304/iterations:32$`
 * Fixed-iteration rows ignore the minimum-time argument. Serialize entire
 * wrapper invocations; do not overlap builds or tests with timing.
 */
// clang-format on
#include <pixie/permutation.h>
#include <pixie/permutations/permutation.h>
