---
name: benchmarks
description: Run Google Benchmark binaries for the Pixie project, including filtering, hardware counters, and perf profiling.
---

# Benchmarks Skill

You now have expertise in running and interpreting Pixie benchmarks. Follow these workflows:

## Build Directory Convention

Uses the same suffix convention as the cmake skill:

```bash
BUILD_SUFFIX=local
```

## CRITICAL: Never Run Benchmarks from a Debug Build

> **Always pass `--config Release` (or `--config RelWithDebInfo`) to `cmake --build`.**
> Multi-config generators (MSVC, Xcode) default to `Debug` if no `--config` is given.
> Google Benchmark will print `***WARNING*** Library was built as DEBUG` and timings will
> be 3-10x slower and meaningless. Always verify the binary path contains `Release/` or
> `RelWithDebInfo/`, never `Debug/`.

## Step 1 — Build

If benchmarks affected by the changes are easily tractable build only related targets.

**Pure timing (benchmarks-all, Release):**
```bash
cmake -B build/benchmarks-all_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON
cmake --build build/benchmarks-all_${BUILD_SUFFIX} --config Release -j
```

**Hardware counters / verbose report (benchmarks-diagnostic, RelWithDebInfo, Linux only):**
```bash
cmake -B build/benchmarks-diagnostic_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=RelWithDebInfo -DPIXIE_BENCHMARKS=ON -DBENCHMARK_ENABLE_LIBPFM=ON -DPIXIE_DIAGNOSTICS=ON
cmake --build build/benchmarks-diagnostic_${BUILD_SUFFIX} --config RelWithDebInfo -j
```

## Step 2 — Run

Prefer running benchmarks with filtering passing the benchmarks that should be affected.

### Available benchmark binaries

| Binary | What it covers |
|--------|---------------|
| `benchmarks` | BitVector rank/select |
| `bench_rmm` | RmM Tree operations |
| `bench_rmm_sdsl` | RmM vs sdsl-lite comparison |
| `louds_tree_benchmarks` | LOUDS Tree traversal |
| `alignment_comparison` | Memory alignment effects |

Binary paths vary by generator type:

| Generator | Path pattern |
|-----------|-------------|
| MSVC / Xcode (multi-config) | `build/<preset>_<suffix>/Release/<binary>` |
| Ninja / Make (single-config) | `build/<preset>_<suffix>/<binary>` |

### Run all benchmarks in a binary

```bash
# Multi-config (MSVC/Xcode)
build/benchmarks-all_${BUILD_SUFFIX}/Release/benchmarks

# Single-config (Ninja/Make)
build/benchmarks-all_${BUILD_SUFFIX}/benchmarks
```

### Filter benchmarks with a regex (FILTER parameter)

```bash
FILTER="BM_Rank"   # change to match benchmark names, e.g. "BM_Select", "BM_Louds", ""

# Multi-config
build/benchmarks-all_${BUILD_SUFFIX}/Release/benchmarks --benchmark_filter="${FILTER}"

# Single-config
build/benchmarks-all_${BUILD_SUFFIX}/benchmarks --benchmark_filter="${FILTER}"
```

Examples:
```bash
# Only rank benchmarks
... --benchmark_filter="BM_Rank"

# Only select on non-interleaved layouts
... --benchmark_filter="BM_Select.*NonInterleaved"

# List all available benchmark names without running
... --benchmark_list_tests=true
```

### Run with hardware counters (benchmarks-diagnostic build, Linux only)

```bash
build/benchmarks-diagnostic_${BUILD_SUFFIX}/RelWithDebInfo/benchmarks \
  --benchmark_filter="${FILTER}" \
  --benchmark_counters_tabular=true
```

### Save results to file

```bash
build/benchmarks-all_${BUILD_SUFFIX}/Release/benchmarks \
  --benchmark_filter="${FILTER}" \
  --benchmark_format=json \
  --benchmark_out=results.json
```

## Step 3 — Profile with perf (Linux only)

Use when hardware counters alone are not enough and you need a full call-graph profile for post-processing.

**Record:**
```bash
perf record -g -F 999 \
  build/benchmarks-diagnostic_${BUILD_SUFFIX}/benchmarks \
  --benchmark_filter="${FILTER}" \
  --benchmark_min_time=5s
```

**Quick report (terminal):**
```bash
perf report --stdio
```

**Flame graph (requires FlameGraph scripts):**
```bash
perf script | stackcollapse-perf.pl | flamegraph.pl > flamegraph.html
```

**Export for external tools (Hotspot, Firefox Profiler):**
```bash
perf script -F +pid > perf.data.txt
# or open with `hotspot perf.data`
```

## Useful Benchmark Flags

| Flag | Purpose |
|------|---------|
| `--benchmark_filter=<regex>` | Run only matching benchmarks |
| `--benchmark_list_tests=true` | List names without running |
| `--benchmark_repetitions=<n>` | Repeat each benchmark n times |
| `--benchmark_min_time=<Ns\|Xs>` | Minimum run time per benchmark |
| `--benchmark_format=json` | Machine-readable output |
| `--benchmark_out=<file>` | Save output to file |
| `--benchmark_counters_tabular=true` | Align hardware counter columns |
| `--benchmark_time_unit=ms` | Change display unit (ns/us/ms/s) |

## Best Practices

1. **Never run from a Debug binary**: always use `--config Release` at build time; check path contains `Release/`
2. **Use benchmarks-all for clean timing**: Release optimizations, no debug info, no libpfm overhead
3. **Use benchmarks-diagnostic for hardware counters**: RelWithDebInfo + libpfm; Linux only
4. **Use perf for deep profiling**: when counters point to a hotspot but don't explain it
5. **Pin CPU frequency** before timing runs: `sudo cpupower frequency-set -g performance`
6. **Filter to reduce noise**: narrow the filter regex to the benchmark under investigation
7. **Save JSON output** when comparing before/after changes: use `--benchmark_out` and diff the files
