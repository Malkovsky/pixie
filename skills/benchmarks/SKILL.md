---
name: benchmarks
description: Run Google Benchmark binaries, including filtering, hardware counters, and perf profiling.
---

# Benchmarks Skill

You now have expertise in running and interpreting Google Benchmark suites.
Follow these workflows:

## Build Directory Convention

Use a short commit hash suffix for committed revisions:

```bash
BUILD_SUFFIX=$(git rev-parse --short HEAD)
```

If the worktree has uncommitted changes, append a descriptive suffix so results
cannot be confused with a clean HEAD build:

```bash
BUILD_SUFFIX=$(git rev-parse --short HEAD)-dirty
```

If not a git repository, use

```bash
BUILD_SUFFIX=agent
```

## CRITICAL: Never Run Benchmarks from a Debug Build

> **Always pass `--config Release` (or `--config RelWithDebInfo`) to `cmake --build`.**
> Multi-config generators (MSVC, Xcode) default to `Debug` if no `--config` is given.
> Google Benchmark will print `***WARNING*** Library was built as DEBUG` and timings will
> be 3-10x slower and meaningless. Always verify the binary path contains `Release/` or
> `RelWithDebInfo/`, never `Debug/`.

## Step 1 — Build

If benchmarks affected by the changes are easily tractable build only related targets.

**Pure timing (benchmarks, Release):**
```bash
cmake -B build/benchmarks_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release
cmake --build build/benchmarks_${BUILD_SUFFIX} --config Release -j
```

**Hardware counters / verbose report (benchmarks-diagnostic, RelWithDebInfo, Linux only):**
```bash
cmake -B build/benchmarks-diagnostic_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBENCHMARK_ENABLE_LIBPFM=ON
cmake --build build/benchmarks-diagnostic_${BUILD_SUFFIX} --config RelWithDebInfo -j
```

For repository-specific benchmark examples, check
`agentic/local/cpp/skills/benchmarks/EXAMPLES.md` when present.

## Step 2 — Run

Prefer running benchmarks with filtering passing the benchmarks that should be affected.

Unless the user explicitly asks otherwise, pin benchmark execution to one CPU with
`taskset` to reduce scheduler noise. Use CPU 0 by default, or override with
`BENCH_CPU=<id>` when a better isolated/performance core is known:

```bash
BENCH_CPU=${BENCH_CPU:-0}
BENCH_RUN="taskset -c ${BENCH_CPU}"
```

If `taskset` is unavailable or fails on the host, report that benchmark results
are unpinned and more noisy.

Execution guardrails:
- Run benchmark commands sequentially in CI.
- Avoid background jobs (`nohup`, `&`) for benchmark collection.
- Always write machine-readable results with `--benchmark_out` when data is later parsed.

### Available benchmark binaries

Discover benchmark binary names from the repository's build system. Common
locations include `build/**/<binary>` for single-config generators and
`build/**/Release/<binary>` for multi-config generators. Repository-specific
binary lists belong in the repo-local benchmark examples overlay.

Binary paths vary by generator type:

| Generator | Path pattern |
|-----------|-------------|
| MSVC / Xcode (multi-config) | `build/<preset>_<suffix>/Release/<binary>` |
| Ninja / Make (single-config) | `build/<preset>_<suffix>/<binary>` |

### Run all benchmarks in a binary

```bash
# Multi-config (MSVC/Xcode)
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/Release/benchmarks

# Single-config (Ninja/Make)
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/benchmarks
```

### Filter benchmarks with a regex (FILTER parameter)

```bash
FILTER="BM_Foo"   # change to match benchmark names in the target binary

# Multi-config
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/Release/benchmarks --benchmark_filter="${FILTER}"

# Single-config
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/benchmarks --benchmark_filter="${FILTER}"
```

Examples:
```bash
# Only one benchmark family
... --benchmark_filter="BM_Foo"

# Only one layout/parameter family
... --benchmark_filter="BM_Foo.*Variant"

# List all available benchmark names without running
... --benchmark_list_tests=true
```

### Run with hardware counters (benchmarks-diagnostic build, Linux only)

The `--benchmark_perf_counters` flag requests hardware counter collection via libpfm. Counter names are platform-specific but common ones include `CYCLES`, `INSTRUCTIONS`, `CACHE-MISSES`, `CACHE-REFERENCES`, `BRANCH-MISSES`, `BRANCH-INSTRUCTIONS`.

```bash
${BENCH_RUN} build/benchmarks-diagnostic_${BUILD_SUFFIX}/RelWithDebInfo/benchmarks \
  --benchmark_filter="${FILTER}" \
  --benchmark_perf_counters=CYCLES,INSTRUCTIONS,CACHE-MISSES \
  --benchmark_counters_tabular=true
```

### Save results to file

```bash
${BENCH_RUN} build/benchmarks_${BUILD_SUFFIX}/Release/benchmarks \
  --benchmark_filter="${FILTER}" \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=results.json
```

Validate output before consuming:
```bash
python3 -m json.tool results.json > /dev/null
```

## Step 3 — Profile with perf (Linux only)

Use when hardware counters alone are not enough and you need a full call-graph profile for post-processing.

**Record:**
```bash
perf record -g -F 999 \
  -- ${BENCH_RUN} build/benchmarks-diagnostic_${BUILD_SUFFIX}/benchmarks \
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
| `--benchmark_perf_counters=CYCLES,INSTRUCTIONS,...` | Collect hardware perf counters (requires libpfm build) |
| `--benchmark_counters_tabular=true` | Align user/perf counter columns into a table |
| `--benchmark_time_unit=ms` | Change display unit (ns/us/ms/s) |

## Best Practices

1. **Never run from a Debug binary**: always use `--config Release` at build time; check path contains `Release/`
2. **Use benchmarks for clean timing**: Release optimizations, no debug info, no libpfm overhead
3. **Use benchmarks-diagnostic for hardware counters**: RelWithDebInfo + libpfm; Linux only
4. **Use perf for deep profiling**: when counters point to a hotspot but don't explain it
5. **Pin benchmark process** with `taskset -c ${BENCH_CPU:-0}` unless unavailable
6. **Pin CPU frequency** before timing runs: `sudo cpupower frequency-set -g performance`
7. **Filter to reduce noise**: narrow the filter regex to the benchmark under investigation
8. **Save JSON output** when comparing before/after changes: use `--benchmark_out` and diff the files
9. **Fail fast on environment issues**: precheck Python deps used by compare tooling (`numpy`, `scipy`)
10. **Use explicit retry limits**: on timeout, narrow scope and retry once; avoid repeated full-suite attempts
11. **Preflight perf counters**: run a tiny counter-enabled benchmark first; if counters unavailable, skip counter workflow
