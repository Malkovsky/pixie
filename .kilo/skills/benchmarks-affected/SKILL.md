---
name: benchmarks-affected
description: Analyze current branch versus a baseline and extract affected benchmark targets and benchmark functions using compile_commands and clang AST.
---

# Benchmarks Affected Skill

Use this skill to identify exactly which benchmark binaries and benchmark functions are affected by code changes on the current branch.

It implements a two-stage workflow:

1. `compile_commands.json` analysis to determine affected compile targets.
2. Clang AST analysis to determine affected benchmark functions.

## Goal

Given `HEAD` and a baseline branch (default `main`), produce:

- Changed files.
- Affected targets (with emphasis on benchmark targets).
- Exact benchmark functions impacted by the changes.
- A ready-to-use Google Benchmark filter regex.

## Prerequisites

1. Build tree with benchmarks enabled and compile database exported:

```bash
BUILD_SUFFIX=local
cmake -B build/benchmarks-all_${BUILD_SUFFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DPIXIE_BENCHMARKS=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build/benchmarks-all_${BUILD_SUFFIX} --config Release -j
```

2. `clang++` must be available on `PATH` (used for AST dump).

## Run

```bash
python3 .kilo/skills/benchmarks-affected/analyze_benchmarks_affected.py \
  --baseline main \
  --compile-commands build/benchmarks-all_local/compile_commands.json \
  --format json
```

If `--compile-commands` is omitted, the script auto-selects the most recently modified `build/**/compile_commands.json`.
Working tree changes are included by default. Use `--no-include-working-tree` to restrict analysis to `<baseline>...HEAD` only.

## Output

The analyzer reports:

- `affected_targets`: impacted CMake targets inferred from compile dependency analysis.
- `affected_benchmark_targets`: subset of benchmark binaries impacted.
- `affected_benchmarks`: precise benchmark function names from AST-level call analysis.
- `suggested_filter_regex`: regex to pass as `--benchmark_filter`.

## How to Use Findings

1. Build only impacted benchmark binaries where feasible.
2. Run benchmark binaries with the suggested filter:

```bash
FILTER='^(BM_RankNonInterleaved|BM_SelectNonInterleaved)(/|$)'
build/benchmarks-all_local/benchmarks --benchmark_filter="${FILTER}"
```

3. If impact mapping is broad/uncertain, run full binary for selected benchmark target(s).

## Guardrails

1. Keep baseline comparison at merge-base style diff: `<baseline>...HEAD`.
2. Use Release binaries for timing runs.
3. If AST parse fails for a TU, still trust compile target impact and mark benchmark-function scope as partial.
4. If benchmark infra (`CMakeLists.txt`, benchmark source layout) changed, fall back to conservative benchmark selection.
