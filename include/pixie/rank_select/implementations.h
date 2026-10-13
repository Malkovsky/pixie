#pragma once

/**
 * @file implementations.h
 * @brief All rank/select support implementations provided by Pixie.
 *
 * - `RankSelectSupport<MetadataStorage>`: non-owning source bits with
 *   storage-backed rank/select metadata.
 * - `RankSelectSupportView`: non-owning source bits with read-only metadata.
 * - `RankSelectSupportCompressed<>`: interleaved packed ranks and adaptive
 *   Elias--Fano select samples with validated source prefetch. Borrows source
 *   words and owns metadata; no serialization.
 * - `PastaRankSelectSupport`: optional Pasta backend owning a converted source
 *   copy and its rank/select metadata.
 */

// clang-format off
/*
 * Rank/select snapshot, 2026-10-06: final compressed select with source prefetch.
 * 2^34 bits, nine source families; Ryzen 7 8845HS / WSL2, CPU 0, GCC 13.3,
 * -O3 -g -DNDEBUG -march=native. Median WALL ns/query, three independent
 * launches per input/backend, rotated backend/operation order, 0.05 s warmup,
 * 0.5 s minimum timing. One backend per process. Construction, generation,
 * validation and source-cache I/O excluded; zero measured-query major faults.
 * 8192-entry shifted query pools, identical native generators/seeds, repetition
 * zero. Four pool passes per supported operation validated before timing.
 * Tables retain baseline, final compressed <4096,16> and Spider only; measured
 * before public-type relocation, experiment cleanup and metadata-storage
 * parameterization, not remeasured after. Secondary-directory allocations now
 * also use the metadata backend; recorded auxiliary bytes are historical.
 * Auxiliary MiB excludes borrowed source words (input occupies 2048 MiB).
 * Spider reports requested allocations minus N/8, including embedded rank bits.
 * Spider interleaved read_and_build_rank(...,4), revision
 * dbf6024cc1c7ac1b02b860ee4bdb6d2ff2117850, temporary wrapper only: no license
 * found, no source vendored. Wrapper reverses each byte's bits and translates
 * rank boundaries. Two select arrays include 128-byte sentinel padding each.
 * Spider has rank1/select1 only; absent operations are '-'. R1+S1 adds the two
 * separately measured medians. Overall=(R1+R0+S1+S0)/2 for Pixie only.
 * Sparse ~= 1/1024 ones; sparse0 complements it. Runs alternate 2^20 bits.
 * Mixed repeats a half-full 2^14-bit block then 64 sparse blocks; mixed0
 * complements it. FNBP uses the Cartesian RMQ monotone-stack generator; this
 * snapshot's temporary harness extended it to 2^34, now registered natively.
 * WSL variation is substantial; small differences do not establish a winner.
 *
 * baseline
 *
 * | Input   | Rank1 | Rank0 | Select1 | Select0 | R1+S1 | Overall | Aux MiB |
 * | :------ | ----: | ----: | ------: | ------: | ----: | ------: | ------: |
 * | 12p5    | 178.8 | 160.6 |   358.4 |   457.5 | 537.2 |   577.6 |  74.000 |
 * | 50      | 177.4 | 153.7 |   413.7 |   487.7 | 591.0 |   616.2 |  74.000 |
 * | 87p5    | 154.8 | 147.8 |   371.7 |   411.4 | 526.5 |   542.9 |  74.000 |
 * | sparse  | 137.6 | 164.6 |   544.1 |   393.0 | 681.7 |   619.7 |  74.000 |
 * | sparse0 | 154.6 | 161.3 |   368.6 |   745.2 | 523.2 |   714.9 |  74.000 |
 * | runs    | 147.1 | 151.6 |   350.4 |   372.3 | 497.6 |   510.7 |  74.000 |
 * | mixed   | 149.8 | 153.2 |   365.4 |   400.0 | 515.2 |   534.2 |  74.000 |
 * | mixed0  | 160.6 | 157.8 |   382.5 |   484.3 | 543.2 |   592.7 |  74.000 |
 * | fnbp    | 149.1 | 152.1 |   373.6 |   394.5 | 522.7 |   534.6 |  74.000 |
 *
 * compressed
 *
 * | Input   | Rank1 | Rank0 | Select1 | Select0 | R1+S1 | Overall | Aux MiB |
 * | :------ | ----: | ----: | ------: | ------: | ----: | ------: | ------: |
 * | 12p5    | 137.3 | 131.7 |   335.5 |   384.5 | 472.7 |   494.4 |  73.126 |
 * | 50      | 127.0 | 120.1 |   392.7 |   374.6 | 519.7 |   507.2 |  73.001 |
 * | 87p5    | 115.3 | 116.0 |   343.5 |   331.1 | 458.9 |   453.0 |  73.126 |
 * | sparse  | 114.6 | 114.7 |   469.4 |   240.5 | 583.9 |   469.5 |  69.074 |
 * | sparse0 | 124.4 | 117.1 |   267.6 |   518.7 | 392.0 |   513.9 |  69.074 |
 * | runs    | 136.5 | 117.9 |     9.4 |     8.9 | 145.9 |   136.3 |  67.001 |
 * | mixed   | 126.4 | 136.6 |   394.5 |   303.5 | 520.9 |   480.5 |  69.328 |
 * | mixed0  | 124.3 | 128.3 |   257.8 |   380.0 | 382.0 |   445.1 |  69.328 |
 * | fnbp    | 124.9 | 119.2 |   315.0 |   357.5 | 439.9 |   458.3 |  73.001 |
 *
 * Spider
 *
 * | Input   | Rank1 | Rank0 | Select1 | Select0 | R1+S1 | Overall | Aux MiB |
 * | :------ | ----: | ----: | ------: | ------: | ----: | ------: | ------: |
 * | 12p5    | 122.1 |     - |   340.6 |       - | 462.7 |       - |  78.129 |
 * | 50      | 148.1 |     - |   354.0 |       - | 502.1 |       - |  78.129 |
 * | 87p5    | 134.4 |     - |   342.8 |       - | 477.2 |       - |  76.879 |
 * | sparse  | 115.6 |     - |   363.3 |       - | 478.9 |       - |  78.129 |
 * | sparse0 | 142.5 |     - |   382.8 |       - | 525.3 |       - |  78.120 |
 * | runs    | 135.3 |     - |   354.6 |       - | 489.8 |       - |  78.129 |
 * | mixed   | 133.5 |     - |   397.3 |       - | 530.9 |       - |  73.668 |
 * | mixed0  | 124.5 |     - |   377.3 |       - | 501.8 |       - |  78.043 |
 * | fnbp    | 124.5 |     - |   367.2 |       - | 491.7 |       - |  78.129 |
 *
 * Equal-weight input geomean, excluding runs: compressed/Spider R1+S1 =
 * 0.942 (5.8% less query time across eight inputs).
 */
// clang-format on

#include <pixie/rank_select.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rank_select/support.h>

#ifdef PIXIE_PASTA_SUPPORT
#include <pixie/rank_select/pasta.h>
#endif
