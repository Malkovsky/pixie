# Benchmark Examples

This file contains Pixie-specific benchmark command examples. Keep the general
benchmark workflow in `SKILL.md`; add concrete repository examples here.

## Standard Pixie Benchmark Builds

Pure timing benchmarks:

```bash
cmake -B build/benchmarks_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON
cmake --build build/benchmarks_${BUILD_SUFFIX} --config Release -j
```

Diagnostic benchmarks with hardware counter support:

```bash
cmake -B build/benchmarks-diagnostic_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=RelWithDebInfo -DPIXIE_BENCHMARKS=ON -DBENCHMARK_ENABLE_LIBPFM=ON -DPIXIE_DIAGNOSTICS=ON
cmake --build build/benchmarks-diagnostic_${BUILD_SUFFIX} --config RelWithDebInfo -j
```

## Pixie Benchmark Binaries

| Binary | What it covers |
|--------|---------------|
| `benchmarks` | BitVector rank/select |
| `bench_rmm` | RmM Tree operations |
| `bench_rmm_btree` | Experimental RmM btree operations |
| `bench_rmm_sdsl` | RmM vs sdsl-lite comparison |
| `louds_tree_benchmarks` | LOUDS Tree traversal |
| `alignment_comparison` | Memory alignment effects |

## RmM Btree vs SDSL Focused Comparison

Build the RmM comparison binaries with third-party backends enabled:

```bash
cmake -B build/benchmarks_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON -DPIXIE_THIRD_PARTY_BACKENDS=ON
cmake --build build/benchmarks_${BUILD_SUFFIX} --config Release --target bench_rmm_btree bench_rmm_sdsl -j
```

Use this configuration when comparing `bench_rmm_btree` against
`bench_rmm_sdsl`. If `bench_rmm_sdsl` does not list the expected operations,
reconfigure with `-DPIXIE_THIRD_PARTY_BACKENDS=ON` before rebuilding; stale
builds can silently omit comparison operations.

For RmM operation-level comparisons, prefer explicit operation and size filters
instead of broad benchmark regexes. First confirm both binaries expose the same
operations:

```bash
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/bench_rmm_btree \
  --ops=close,open,enclose,fwdsearch \
  --explicit_sizes=65536 \
  --benchmark_list_tests=true

${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/bench_rmm_sdsl \
  --ops=close,open,enclose,fwdsearch \
  --explicit_sizes=65536 \
  --benchmark_list_tests=true
```

Then run both binaries with identical query count, sizes, operations, and
repetitions. Use JSON output so the rows can be paired side by side:

```bash
OPS=close,open,enclose,fwdsearch
SIZES=65536,1048576,8388608,67108864
Q=32768
REPS=10

${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/bench_rmm_btree \
  --ops=${OPS} \
  --explicit_sizes=${SIZES} \
  --queries=${Q} \
  --benchmark_repetitions=${REPS} \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/rmm_btree_${BUILD_SUFFIX}.json

${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/bench_rmm_sdsl \
  --ops=${OPS} \
  --explicit_sizes=${SIZES} \
  --queries=${Q} \
  --benchmark_repetitions=${REPS} \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/rmm_sdsl_${BUILD_SUFFIX}.json
```

When reporting, pair rows by `(size, operation)` and include both absolute times
and ratio (`btree / sdsl`). Include coefficient of variation when repetitions
are used; high CV means the ratio should be treated as directional rather than
precise.
