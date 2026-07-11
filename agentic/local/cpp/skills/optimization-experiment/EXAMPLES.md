# Pixie Optimization Experiment Examples

Use these notes with `agentic/cpp/skills/optimization-experiment/SKILL.md` for
Pixie-specific optimization work.

## Example: `excess_min_128`

The `excess_min_128` experiment is a good template for future hot-path work:

1. **Lock semantics first**: preserve inclusive prefix range `[left, right]`,
   prefix offsets `[0, 128]`, first-min tie breaking, and invalid-range sentinel
   behavior.
2. **Add focused benchmarks before trusting results**: include full-range,
   RMQ-style `[0, 127]`, aligned short ranges, non-aligned short ranges,
   cross-word ranges, point ranges, and a reproducible random-range benchmark.
3. **Keep candidates same-API**: experimental variants should accept the same
   parameters and return the same result type as production, so benchmarks can
   swap implementations without adapter logic.
4. **Use experimental space for ideas**: keep exploratory variants in
   `include/pixie/experimental/` unless they are clearly production-ready.
5. **Promote narrowly**: when a variant wins only under a specific condition,
   promote only that condition. For `excess_min_128`, byte-LUT was useful only
   for byte-aligned short ranges, so the production fallback was narrowed instead
   of applied to every short range.
6. **Mark unlikely cold dispatches**: if a promoted optimization is for a narrow
   special case, mark the branch unlikely when that matches expected workload.
7. **Persist only accepted benchmark evidence in the repository**: add or update
   a top-of-file benchmark note when the measurement is the current production
   state, an accepted baseline, or the final comparison for a kept experimental
   variant. Include metric units, fixed decimal precision, and visible
   separators between logical result blocks. Do not include machine-local JSON
   paths such as `/tmp/...` in source comments, and do not accumulate every
   discarded probe in production headers.

## Benchmark Pattern

Prefer a diagnostic run that can explain timing changes, not just report them:

```bash
taskset -c 0 build/benchmarks-diagnostic_local/excess_positions_benchmarks \
  --benchmark_filter='BM_ExcessMin128(_|/|$)' \
  --benchmark_repetitions=3 \
  --benchmark_perf_counters=CYCLES,INSTRUCTIONS,CACHE-MISSES \
  --benchmark_counters_tabular=true \
  --benchmark_out=/tmp/pixie_excess_min_counters.json \
  --benchmark_out_format=json
```

For final comparison tables, include at least:

- CPU time
- cycles
- instructions
- cache misses when collected
- a random or mixed workload row
- the rows that justified any narrowed production dispatch condition

## Example: RMQ Optimization Report

The current RMQ work is the preferred template for larger optimization
experiments. It records the command shape, units, precision, and the main
tradeoff dimensions: query latency, build time, absolute memory, and normalized
memory.

Keep this kind of report at the top of the relevant header or experiment log:

```text
RMQ benchmark snapshot, 2026-06-11.

Query command:
  taskset -c 0 ./build/release/bench_rmq
    --benchmark_filter='^(rmq_sparse_table|rmq_segment_tree|rmq_hybrid_btree|rmq_cartesian_rmm|rmq_cartesian_hybrid_btree|rmq_sdsl_sct)/'
    --benchmark_min_time=0.25s

Build/memory command:
  taskset -c 0 ./build/release/bench_rmq
    --benchmark_filter='^(rmq_build_sparse_table|rmq_build_segment_tree|rmq_build_hybrid_btree|rmq_build_cartesian_rmm|rmq_build_cartesian_hybrid_btree|rmq_build_sdsl_sct)/'
    --benchmark_min_time=0.20s

Query CPU time. "max width" is the benchmark's maximum sampled query width.
Sparse-table rows are intentionally not registered above 2^22. "cartesian hybrid"
is the `rmq_cartesian_hybrid_btree` row.

| N      | max width   | sparse table (ns)   | segment tree (ns)   | hybrid btree (ns)   | cartesian hybrid (ns)   | sdsl sct (ns)   |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^10 |          64 |                10.1 |                35.9 |                    34.6 |                66.4 |           111.2 |
|   2^10 |        2^10 |                13.1 |                58.8 |                    41.0 |               118.6 |           231.3 |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^14 |          64 |                14.7 |                43.6 |                    33.5 |                67.3 |           125.4 |
|   2^14 |        4096 |                15.9 |                76.7 |                    38.1 |                66.3 |           430.8 |
|   2^14 |        2^14 |                15.8 |                93.5 |                    33.7 |                49.2 |           532.1 |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^18 |          64 |                63.0 |               103.7 |                    50.9 |                74.7 |           141.3 |
|   2^18 |        4096 |                52.2 |               183.0 |                    52.5 |                72.6 |           586.7 |
|   2^18 |        2^18 |                40.9 |               299.1 |                    46.0 |                27.6 |           791.0 |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^22 |          64 |                97.2 |               169.1 |                   137.4 |               256.9 |           263.2 |
|   2^22 |        4096 |                93.4 |               298.7 |                   101.9 |               251.9 |           726.7 |
|   2^22 |        2^18 |                68.0 |               437.3 |                    82.2 |                86.4 |          1053.0 |
|   2^22 |        2^22 |                60.7 |               485.3 |                    45.5 |                22.3 |          1129.7 |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^24 |          64 |                   - |               208.0 |                   166.6 |               480.4 |           549.6 |
|   2^24 |        4096 |                   - |               357.2 |                   222.2 |               535.4 |           864.8 |
|   2^24 |        2^18 |                   - |               472.3 |                   143.0 |                83.5 |          2824.2 |
|   2^24 |        2^22 |                   - |               510.7 |                    47.1 |                25.3 |          1169.8 |
|   2^24 |        2^24 |                   - |               546.5 |                    43.4 |                20.1 |          1239.6 |
| -----: | ----------: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^26 |          64 |                   - |               273.5 |                   342.3 |               567.5 |           710.4 |
|   2^26 |        4096 |                   - |               938.4 |                   308.9 |               688.5 |          1701.6 |
|   2^26 |        2^18 |                   - |               802.5 |                   175.6 |               195.1 |          2609.4 |
|   2^26 |        2^22 |                   - |              1690.7 |                    64.0 |                45.2 |          2387.0 |
|   2^26 |        2^26 |                   - |              1003.7 |                    51.0 |                24.1 |          2152.0 |

Build CPU time. Sparse-table build rows are intentionally not registered
above 2^22.

| N      | sparse table (ms)   | segment tree (ms)   | hybrid btree (ms)   | cartesian hybrid (ms)   | sdsl sct (ms)   |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^10 |               0.005 |               0.001 |                   0.001 |               0.009 |           0.007 |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^14 |               1.096 |               0.015 |                   0.018 |               0.058 |           0.162 |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^18 |              28.405 |               0.292 |                   0.301 |               1.985 |           2.484 |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^22 |             860.619 |              92.812 |                   6.077 |              31.552 |          40.262 |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^24 |                   - |             426.733 |                  21.802 |             124.868 |         158.091 |
| -----: | ------------------: | ------------------: | ----------------------: | ------------------: | --------------: |
|   2^26 |                   - |            1432.849 |                  86.641 |             518.439 |         647.360 |

Owned auxiliary memory. Benchmarks use 64-bit signed integer values and
64-bit indexes. The external input values are not owned by RMQ indexes and
are excluded.

| N      | sparse table (MiB)   | segment tree (MiB)   | hybrid btree (MiB)   | cartesian hybrid (MiB)   | sdsl sct (MiB)   |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^10 |                0.071 |                0.016 |                    0.009 |                0.016 |            0.001 |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^14 |                1.626 |                0.250 |                    0.011 |                0.022 |            0.005 |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^18 |               34.001 |                4.000 |                    0.065 |                0.186 |            0.079 |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^22 |              672.001 |               64.000 |                    0.801 |                2.632 |            1.262 |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^24 |                    - |              256.000 |                    3.155 |                7.525 |            5.051 |
| -----: | -------------------: | -------------------: | -----------------------: | -------------------: | ---------------: |
|   2^26 |                    - |             1024.000 |                    8.340 |               30.530 |           20.224 |

Owned auxiliary memory normalized by indexed value count (bits/value).

| N      | sparse table   | segment tree   | hybrid btree   | cartesian hybrid   | sdsl sct   |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^10 |        579.188 |        128.438 |             71.438 |        127.500 |      7.312 |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^14 |        832.262 |        128.027 |              5.434 |         11.102 |      2.770 |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^18 |       1088.020 |        128.002 |              2.089 |          5.957 |      2.534 |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^22 |       1344.002 |        128.000 |              1.603 |          5.264 |      2.524 |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^24 |              - |        128.000 |              1.577 |          3.763 |      2.526 |
| -----: | -------------: | -------------: | -----------------: | -------------: | ---------: |
|   2^26 |              - |        128.000 |              1.042 |          3.816 |      2.528 |
```

## Documentation Pattern

For experimental algorithms, use Doxygen block comments:

```cpp
/**
 * @brief Short algorithm name and purpose.
 *
 * @details Workflow:
 *
 *   input -> transform -> candidate result -> final reduction
 *
 * Explain the non-obvious tradeoff, such as avoiding lane crossing, reducing
 * loop iterations, or paying scalar boundary work.
 */
```

Keep long benchmark tables as a top-level `/** ... */` note when they describe
the whole experimental header rather than one symbol. For final or accepted
results, do not leave measurements only in chat or `/tmp`; persist a
top-of-file snapshot before finishing. For short-lived diagnostics, summarize
the finding without adding source-file churn.

## Non-Obvious Lessons From This Session

- A faster micro-kernel can regress neighboring ranges if dispatch is too broad;
  benchmark aligned and non-aligned cases separately.
- Random-range benchmarks are useful as a sanity check against overfitting fixed
  rows.
- Hardware counters helped distinguish lower-level cost: cycles and instructions
  were more actionable than cache misses for the `excess_min_128` candidates.
- Production and experimental code must stay visually distinct. A benchmark win
  in `pixie::experimental` is only a candidate, not a caller-visible change.
- When result names encode width unnecessarily, prefer the stable semantic name
  (`ExcessResult`) over width-specific detail (`ExcessResult128`) unless the API
  truly needs multiple widths.
