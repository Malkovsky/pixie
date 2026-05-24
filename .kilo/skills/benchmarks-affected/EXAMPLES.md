# Benchmarks Affected Examples

This file contains Pixie-specific examples for the `benchmarks-affected` skill.
Keep general analyzer workflow in `SKILL.md`.

## Pixie Compile Database Setup

Configure a benchmark-enabled build with compile commands:

```bash
BUILD_SUFFIX=local
cmake -B build/benchmarks-all_${BUILD_SUFFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DPIXIE_BENCHMARKS=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build/benchmarks-all_${BUILD_SUFFIX} --config Release -j
```

Run affected benchmark analysis:

```bash
python3 .kilo/skills/benchmarks-affected/analyze_benchmarks_affected.py \
  --baseline main \
  --compile-commands build/benchmarks-all_local/compile_commands.json \
  --format json
```

Example filtered benchmark run from analyzer output:

```bash
FILTER='^(BM_RankNonInterleaved|BM_SelectNonInterleaved)(/|$)'
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" build/benchmarks-all_local/benchmarks --benchmark_filter="${FILTER}"
```
