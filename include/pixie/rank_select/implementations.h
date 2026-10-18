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
 * - `SdslRankSelectSupport`: optional SDSL rank-v and MCL select backend.
 * - `SuccinctFastRankSelectSupport`: optional FastBitVec backend; no
 * zero-select.
 * - `OrzoRankSelectSupport`: optional Orzo backend; no zero-select.
 * - `ThreeStarRankSelectSupport`: optional 3-star backend; no zero-select,
 *   process-global state, and only one live instance permitted.
 */

// clang-format off
/*
 * Rank/select snapshot, 2026-10-13--16: random inputs, 2^34 bits.
 * Ryzen 7 8845HS / WSL2, CPU 0, GCC 13.3, -O3 -DNDEBUG -march=native.
 * Median WALL ns/query over three repetitions, 20 million queries per
 * repetition; fixed input seeds and 8192-entry shifted query pools.
 * Build is median wall seconds over ten repetitions, two builds per
 * repetition. SDSL was measured separately on October 16. These measurements
 * precede the optional-adapter CRTP migration and header relocation.
 * Aux is logical index bytes / input bytes (%), excluding owned source copies;
 * build includes source copies. SDSL uses rank-v + MCL for both bit values and
 * serialized directory sizes. Missing operations are marked with a dash.
 * FastBitVec is omitted here: its benchmark size grid stops at 2^30.
 * Spider was not measured in this run; historical Spider results are omitted.
 * WSL variation is substantial; small differences do not establish a winner.
 * Full size curves, machine details and reproduction commands:
 * src/docs/benchmarks/rank_select.md (Doxygen: Benchmarks / Rank/select).
 *
 * 12.5% ones
 *
 * | Backend    | Rank1 | Rank0 | Select1 | Select0 | Build s |  Aux % |
 * | :--------- | ----: | ----: | ------: | ------: | ------: | -----: |
 * | Baseline   |  99.0 |  97.5 |   285.7 |   315.2 |    0.66 |  3.613 |
 * | Compressed |  93.8 |  97.5 |   248.6 |   292.4 |    0.67 |  3.571 |
 * | Pasta      | 119.1 | 122.2 |   329.8 |   413.8 |    2.60 |  3.516 |
 * | 3-star     | 267.6 | 264.3 |   447.7 |       - |    8.16 |  0.930 |
 * | Orzo       | 132.1 | 141.5 |   441.4 |       - |    2.90 |  2.346 |
 * | SDSL       |  92.4 |  98.1 |   453.5 |   460.6 |    6.89 | 48.221 |
 *
 * 50% ones
 *
 * | Backend    | Rank1 | Rank0 | Select1 | Select0 | Build s |  Aux % |
 * | :--------- | ----: | ----: | ------: | ------: | ------: | -----: |
 * | Baseline   |  95.8 | 122.1 |   317.4 |   300.8 |    0.65 |  3.613 |
 * | Compressed |  86.2 |  98.2 |   277.6 |   275.5 |    0.50 |  3.564 |
 * | Pasta      | 139.5 | 119.3 |   248.4 |   564.2 |    2.36 |  3.516 |
 * | 3-star     | 295.0 | 257.5 |   373.1 |       - |    7.64 |  0.785 |
 * | Orzo       | 139.4 | 131.2 |   385.8 |       - |    2.72 |  2.493 |
 * | SDSL       |  91.2 | 103.2 |   360.3 |   439.0 |   10.79 | 48.537 |
 *
 * 87.5% ones
 *
 * | Backend    | Rank1 | Rank0 | Select1 | Select0 | Build s |  Aux % |
 * | :--------- | ----: | ----: | ------: | ------: | ------: | -----: |
 * | Baseline   | 105.8 | 123.7 |   332.1 |   282.7 |    0.65 |  3.613 |
 * | Compressed |  89.5 |  94.9 |   271.5 |   245.0 |    0.70 |  3.571 |
 * | Pasta      | 122.4 | 127.1 |   307.0 |   773.7 |    2.67 |  3.516 |
 * | 3-star     | 287.2 | 278.5 |   465.5 |       - |    8.21 |  0.785 |
 * | Orzo       | 140.6 | 147.8 |   482.3 |       - |    2.94 |  2.639 |
 * | SDSL       |  95.5 |  94.9 |   435.4 |   448.7 |    7.70 | 48.244 |
 */
// clang-format on

#include <pixie/rank_select.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rank_select/support.h>

#ifdef PIXIE_PASTA_SUPPORT
#include <pixie/rank_select/pasta.h>
#endif
#ifdef SDSL_SUPPORT
#include <pixie/rank_select/sdsl.h>
#endif
#ifdef PIXIE_3STAR_SUPPORT
#include <pixie/rank_select/three_star.h>
#endif
#ifdef PIXIE_SUCCINCT_BITVEC_SUPPORT
#include <pixie/rank_select/succinct_bitvec.h>
#endif
#ifdef PIXIE_ORZO_SUPPORT
#include <pixie/rank_select/orzo.h>
#endif
