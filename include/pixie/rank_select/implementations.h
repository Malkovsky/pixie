#pragma once

/**
 * @file implementations.h
 * @brief All rank/select support implementations provided by Pixie.
 *
 * - `RankSelectSupport<MetadataStorage>`: non-owning source bits with
 *   storage-backed rank/select metadata.
 */

// clang-format off
/*
 * Rank/select benchmark snapshot, 2026-07-13.
 *
 * Tables report one pinned Release pass on CPU 0, with Google Benchmark's
 * 0.2 s warmup and 1.0 s minimum time. Query setup, support construction,
 * and the 512 KiB query pool are outside the timed query region. Every row
 * constructs RankSelectSupport with SelectSupport::kBoth, so auxiliary
 * memory includes both select directions and excludes the external source.
 *
 * Random sources use the deterministic benchmark generator (seed 42) with the
 * stated expected one-fill percentage. FNBP is the deterministic
 * Ferrada-Navarro balanced-parentheses source produced by the Cartesian RMQ
 * backward monotone-stack construction. FNBP rows reach 2^30; its source
 * construction is linear in encoded bits.
 *
 * Query CPU time, ns.
 *
 * | source       |    N |   rank1 |   rank0 | select1 | select0 |
 * | :----------- | ---: | ------: | ------: | ------: | ------: |
 * | random 12.5% | 2^10 |   2.731 |   2.864 |  19.286 |  15.712 |
 * | random 12.5% | 2^14 |   2.726 |   2.904 |  20.145 |  23.044 |
 * | random 12.5% | 2^18 |   2.826 |   2.937 |  19.492 |  16.752 |
 * | random 12.5% | 2^22 |   3.466 |   3.250 |  20.785 |  18.527 |
 * | random 12.5% | 2^26 |   7.378 |  12.393 |  60.703 |  59.822 |
 * | random 12.5% | 2^30 |  27.321 |  22.617 | 156.325 | 119.631 |
 * | random 12.5% | 2^34 |  53.798 |  64.239 | 275.460 | 300.452 |
 * | random 50%   | 2^10 |   2.774 |   2.891 |  20.117 |  22.126 |
 * | random 50%   | 2^14 |   2.750 |   2.877 |  20.415 |  19.468 |
 * | random 50%   | 2^18 |   2.840 |   3.027 |  19.088 |  19.561 |
 * | random 50%   | 2^22 |   3.172 |   3.352 |  21.740 |  20.402 |
 * | random 50%   | 2^26 |  10.587 |   7.433 |  55.366 |  59.901 |
 * | random 50%   | 2^30 |  22.155 |  22.636 | 154.189 | 122.074 |
 * | random 50%   | 2^34 |  55.595 |  59.473 | 288.319 | 299.892 |
 * | random 87.5% | 2^10 |   2.795 |   2.934 |  19.643 |  22.013 |
 * | random 87.5% | 2^14 |   2.769 |   2.931 |  20.337 |  18.226 |
 * | random 87.5% | 2^18 |   2.833 |   2.992 |  18.922 |  21.248 |
 * | random 87.5% | 2^22 |   3.119 |   3.285 |  21.576 |  24.707 |
 * | random 87.5% | 2^26 |   8.828 |  11.294 |  70.811 |  57.801 |
 * | random 87.5% | 2^30 |  22.653 |  22.509 | 147.442 | 129.210 |
 * | random 87.5% | 2^34 |  53.092 |  56.994 | 331.554 | 241.458 |
 * | FNBP         | 2^10 |   3.795 |   4.173 |  19.302 |  21.880 |
 * | FNBP         | 2^14 |   3.799 |   4.115 |  20.284 |  19.068 |
 * | FNBP         | 2^18 |   3.871 |   4.156 |  18.976 |  16.639 |
 * | FNBP         | 2^22 |   4.710 |   4.665 |  21.420 |  18.803 |
 * | FNBP         | 2^26 |   7.478 |   8.516 |  53.575 |  45.020 |
 * | FNBP         | 2^30 |  20.484 |  29.839 | 142.972 | 120.746 |
 *
 * Construction CPU time, ms, and owned auxiliary metadata. Random rows use
 * the 50% source; this kBoth benchmark reserves combined select-sample
 * capacity from the total bit count, independently of the fill ratio.
 *
 * | source     |    N |   build | aux MiB | aux bits ratio |
 * | :--------- | ---: | ------: | ------: | -------------: |
 * | random 50% | 2^10 |   0.009 |   0.002 |         19.500 |
 * | random 50% | 2^14 |   0.009 |   0.002 |          1.219 |
 * | random 50% | 2^18 |   0.022 |   0.005 |          0.145 |
 * | random 50% | 2^22 |   0.141 |   0.020 |          0.041 |
 * | random 50% | 2^26 |   2.172 |   0.291 |          0.036 |
 * | random 50% | 2^30 |  35.674 |   4.627 |          0.036 |
 * | random 50% | 2^34 | 639.913 |  74.002 |          0.036 |
 * | FNBP       | 2^10 |   0.009 |   0.002 |         19.500 |
 * | FNBP       | 2^14 |   0.009 |   0.002 |          1.219 |
 * | FNBP       | 2^18 |   0.021 |   0.005 |          0.145 |
 * | FNBP       | 2^22 |   0.139 |   0.020 |          0.041 |
 * | FNBP       | 2^26 |   2.289 |   0.291 |          0.036 |
 * | FNBP       | 2^30 |  36.093 |   4.627 |          0.036 |
 */
// clang-format on

#include <pixie/rank_select.h>
#include <pixie/rank_select/support.h>
