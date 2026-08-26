#pragma once

/**
 * @file implementations.h
 * @brief All integer-vector implementations provided by Pixie.
 *
 * - `PackedIntegerVector`: owning runtime-width packed integers.
 * - `PackedIntegerVectorView`: zero-copy packed-integer view.
 * - `PackedMonotoneIntegerVector`: owning validated monotone packed integers.
 * - `PackedMonotoneIntegerVectorView`: zero-copy monotone packed-integer view.
 */

// clang-format off
/*
 * Packed monotone integer-vector benchmark snapshot, 2026-08-25.
 *
 * The table reports one pinned Release pass on CPU 0 of an AMD Ryzen 7 8845HS,
 * with Google Benchmark's 0.1 s warmup and 0.5 s minimum time. Construction,
 * serialization, view restoration, and the deterministic 2^16-query pool are
 * outside the timed region. Times are CPU nanoseconds per lower-bound query,
 * rounded to the nearest nanosecond.
 *
 * Dense values are `i`; duplicate-heavy values are `i / 16`; sparse values
 * are deterministic cumulative increments in `[1, 1024]` (seed 42). Queries
 * alternate between sampled present values and their successor.
 *
 * | dataset         |    N | owner | view |
 * | :-------------- | ---: | ----: | ---: |
 * | dense           | 2^10 |    67 |   67 |
 * | dense           | 2^14 |    90 |   91 |
 * | dense           | 2^18 |   142 |  142 |
 * | dense           | 2^22 |   515 |  566 |
 * | dense           | 2^26 |  1546 | 1659 |
 * | duplicate-heavy | 2^10 |    45 |   43 |
 * | duplicate-heavy | 2^14 |    68 |   70 |
 * | duplicate-heavy | 2^18 |   118 |  122 |
 * | duplicate-heavy | 2^22 |   446 |  486 |
 * | duplicate-heavy | 2^26 |  1517 | 1546 |
 * | sparse          | 2^10 |    60 |   68 |
 * | sparse          | 2^14 |    94 |   93 |
 * | sparse          | 2^18 |   151 |  170 |
 * | sparse          | 2^22 |   746 |  884 |
 * | sparse          | 2^26 |  1956 | 1907 |
 */
// clang-format on

#include <pixie/integer_vector.h>
#include <pixie/integer_vector/monotone.h>
#include <pixie/integer_vector/packed.h>
#include <pixie/integer_vector/packed_monotone.h>
