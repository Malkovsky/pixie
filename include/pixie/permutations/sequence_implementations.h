#pragma once

// clang-format off
/**
 * @file sequence_implementations.h
 * @brief Permutable-sequence implementation catalog for benchmarks.
 * @details Consumers include pixie/permutations/sequence.h directly.
 *
 * Current snapshot, 2026-09-21: Ryzen 7 8845HS, WSL/Linux, GCC 13.3,
 * -O3 -DNDEBUG -march=native, Google Benchmark 1.9.4, CPU 0, 16 MiB L3.
 * uint64_t, B256/F8/cumulative, indirect C4096. Ranges span two current-code
 * CPU medians (5 repetitions at 0.15s, 9 at 0.2s), not confidence intervals.
 * Repeated raw controls also drifted materially; no close ranking or universal
 * storage-policy winner is established.
 *
 * | Representation | Read 1 MiB ns | Read 32 MiB ns | Read 128 MiB ns | Owned at 32 MiB |
 * | -------------- | ------------- | -------------- | --------------- | --------------- |
 * | Packed64       |            43 |        179-204 |         320-356 |       36.57 MiB |
 * | Indirect64     |         41-43 |        169-191 |         326-339 |       68.82 MiB |
 * | Raw64 control  |            42 |        191-211 |         311-368 |       36.57 MiB |
 *
 * | Storage  | Global 32 MiB us        | Local65 32 MiB us       | Equal merge 32 MiB us |
 * | -------- | ----------------------- | ----------------------- | --------------------- |
 * | Packed   |               32.3-33.9 |               23.6-23.9 |             12.5-13.1 |
 * | Indirect |               32.9-34.0 |               23.6-24.1 |             12.4-13.2 |
 *
 * Useful MiB describes value payload, not total footprint. Reads include hash
 * generation, 64 reads/iteration. Rotation uses 128 fixed-stream operations,
 * 8 batches/repetition. Local65 means elements and can cross leaves. Equal
 * merge uses combined N=4194304, two independent merges/iteration and 32
 * iterations. Mutation setup/destruction is excluded; rebuilding pre-touches
 * storage. Short timed totals remain sensitive to host/harness variation.
 *
 * Fresh-tree requested bytes include alignment and capacity but exclude
 * allocator overhead, RSS, temporary peaks, and allocations inside T. Packed
 * and indirect owners are 24 and 40 bytes. Bulk indirect storage at 32 MiB has
 * 8192 chunks, 0.25 MiB chunk headers, and no vector slack. Constructing and
 * merging 4096 singleton uint64_t sequences retains 4096 chunks and 16.16 MiB,
 * including 15.97 MiB vector slack (2.53-2.57 us per singleton build+merge).
 * Stable references deliberately retain this slack; bulk construction or an
 * explicit smaller ChunkBytes avoids it without implicit compaction. Tests,
 * not these timings, establish that merge/rotation do not move payload objects.
 * Larger-object read rows inspect only the first word, not complete objects.
 *
 * Reproduce each filter separately, using the local /bench wrapper:
 * @code{.sh}
 * flock /tmp/kilo/pixie-experiment-timing.lock \
 *   env CPP_BENCH_CPU=0 CPP_JOBS=4 \
 *   /home/user/.config/kilo/scripts/cpp-bench native \
 *   permutable_sequence_benchmarks '<filter>' 9 0.2s
 * @endcode
 * Read filter:
 * `^FacadeAccessIndependent/(Packed64|Indirect64|Raw64)_B256_F8_Prefix(_C4096)?/N:(131072|4194304|16777216)$`
 * Rotation filter:
 * `^FacadeRotate(Global|Local65)Batch/(Packed64|Indirect64|Raw64)_B256_F8_Prefix(_C4096)?/N:(131072|4194304)/iterations:8$`
 * Merge filter:
 * `^FacadeMergeEqual/(Packed64|Indirect64)_B256_F8_Prefix(_C4096)?/N:4194304/iterations:32$`
 * Singleton filter:
 * `^FacadeSingletonBuildMerge/Indirect64_B256_F8_Prefix_C4096/N:4096/iterations:8$`
 * Fixed-iteration rows ignore the minimum-time argument. Serialize entire
 * wrapper invocations; do not overlap builds or tests with timing.
 */
// clang-format on
#include <pixie/permutable_sequence.h>
#include <pixie/permutations/sequence.h>
