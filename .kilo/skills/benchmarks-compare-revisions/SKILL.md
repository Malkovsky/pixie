---
name: benchmarks-compare-revisions
description: Compare benchmark performance between two git revisions via Google Benchmark compare.py, with optional hardware-counter comparison from diagnostic libpfm builds.
---

# Benchmarks Compare Revisions Skill

Use this skill to compare performance between two git revisions.

This workflow now depends on:

1. `.kilo/skills/benchmarks-affected/SKILL.md` to determine affected benchmark targets/functions and produce a benchmark filter.
2. `.kilo/skills/benchmarks/SKILL.md` for build/run operational details.

## Goal

Build two separate benchmark binaries using short commit hashes as build suffixes, compare timing results with Google Benchmark compare.py, and optionally compare hardware counters across the same revisions.

## Step 0 — Choose revisions, hashes, and options

Pick a baseline and a contender revision. Use short commit hashes to suffix build directories so builds do not collide.

Optional behavior flags:

- `COLLECT_COUNTERS=1` to enable hardware-counter collection and analysis in addition to timing comparison.
- `COLLECT_COUNTERS=0` to run timing-only comparison.

Counter collection is Linux-only and requires:

- diagnostic builds with `BENCHMARK_ENABLE_LIBPFM=ON`
- perf permissions on the host (for access to performance counters)

Example:
```bash
BASELINE=abc1234
CONTENDER=def5678
```

## Step 1 — Compute affected benchmark scope first

Run `benchmarks-affected` from the contender checkout to derive the compare scope.

Do not duplicate `benchmarks-affected` internals here (compile database selection, AST analysis, or fallback heuristics). Follow that skill directly and consume only its outputs.

Inputs to pass through:

- `--baseline ${BASELINE}`
- optional compile-commands path if auto-detection is not desired
- optional output format (`json` recommended for parsing)

Consume these outputs from `benchmarks-affected`:

- `suggested_filter_regex` -> set `FILTER`
- `affected_benchmark_targets` -> optionally constrain which benchmark binary/binaries to run
- `affected_benchmarks` -> function-level scope for validation/reporting

If `FILTER` is empty, fall back to full benchmark binary compare (conservative mode).

## Step 2 — Build both revisions

Use the existing benchmarks skill build steps, but set the build suffix to include the short hash for each revision.

Always build Release timing binaries.

If `COLLECT_COUNTERS=1`, also build diagnostic binaries (RelWithDebInfo + libpfm) for both revisions.

```bash
# Baseline
BUILD_SUFFIX=bench_${BASELINE}
git checkout ${BASELINE}
# Follow .kilo/skills/benchmarks/SKILL.md timing build instructions with this suffix
# If COLLECT_COUNTERS=1, also follow the diagnostic build instructions with this suffix

# Contender
BUILD_SUFFIX=bench_${CONTENDER}
git checkout ${CONTENDER}
# Follow .kilo/skills/benchmarks/SKILL.md timing build instructions with this suffix
# If COLLECT_COUNTERS=1, also follow the diagnostic build instructions with this suffix
```

Expected build trees:

- Timing: `build/benchmarks-all_bench_<short-hash>`
- Counters (optional): `build/benchmarks-diagnostic_bench_<short-hash>`

## Step 3 — Compare using compare.py

Use Google Benchmark compare tooling with a JSON-first flow to avoid long-running binary-vs-binary retries.

Locate compare.py from the Google Benchmark dependency (installed under the build tree):
```bash
COMPARE_PY=build/benchmarks-all_bench_${BASELINE}/_deps/googlebenchmark-src/tools/compare.py
```

Verify Python deps once (compare.py imports numpy/scipy):
```bash
python3 -c "import numpy, scipy"
```

Generate baseline/contender JSON sequentially with explicit file outputs:
```bash
BASE_JSON=/tmp/bench_${BASELINE}.json
CONT_JSON=/tmp/bench_${CONTENDER}.json

build/benchmarks-all_bench_${BASELINE}/benchmarks \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=${BASE_JSON} > /tmp/bench_${BASELINE}.log 2>&1

build/benchmarks-all_bench_${CONTENDER}/benchmarks \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=${CONT_JSON} > /tmp/bench_${CONTENDER}.log 2>&1
```

Validate JSON before comparing:
```bash
python3 -m json.tool ${BASE_JSON} > /dev/null
python3 -m json.tool ${CONT_JSON} > /dev/null
```

Run the comparison:
```bash
python3 ${COMPARE_PY} -a benchmarks ${BASE_JSON} ${CONT_JSON}
```

Use the affected filter from Step 1 when generating JSON files:
```bash
if [ -n "${FILTER}" ]; then
  FILTER_ARG="--benchmark_filter=${FILTER}"
else
  FILTER_ARG=""
fi

build/benchmarks-all_bench_${BASELINE}/benchmarks ${FILTER_ARG} --benchmark_report_aggregates_only=true --benchmark_display_aggregates_only=true ...
build/benchmarks-all_bench_${CONTENDER}/benchmarks ${FILTER_ARG} --benchmark_report_aggregates_only=true --benchmark_display_aggregates_only=true ...
```

## Step 3b — Compare hardware counters (optional, Linux only)

Run this step only when `COLLECT_COUNTERS=1`.

1. Preflight first with one tiny counter-enabled benchmark run from a diagnostic binary. If output includes warnings such as `Failed to get a file descriptor for performance counter`, mark counters unavailable and skip counter collection.
2. Run baseline and contender diagnostic binaries sequentially with explicit JSON outputs and the same filter scope:

```bash
BASE_COUNTERS_JSON=/tmp/bench_counters_${BASELINE}.json
CONT_COUNTERS_JSON=/tmp/bench_counters_${CONTENDER}.json

build/benchmarks-diagnostic_bench_${BASELINE}/benchmarks \
  ${FILTER_ARG} \
  --benchmark_counters_tabular=true \
  --benchmark_format=json \
  --benchmark_out=${BASE_COUNTERS_JSON} > /tmp/bench_counters_${BASELINE}.log 2>&1

build/benchmarks-diagnostic_bench_${CONTENDER}/benchmarks \
  ${FILTER_ARG} \
  --benchmark_counters_tabular=true \
  --benchmark_format=json \
  --benchmark_out=${CONT_COUNTERS_JSON} > /tmp/bench_counters_${CONTENDER}.log 2>&1
```

3. Validate JSON files before consuming:

```bash
python3 -m json.tool ${BASE_COUNTERS_JSON} > /dev/null
python3 -m json.tool ${CONT_COUNTERS_JSON} > /dev/null
```

4. Collect and compare these counter families when present:

- `instructions`, `cycles`
- `cache-misses`, `cache-references`
- `branch-misses`, `branches`
- `L1-dcache-load-misses`

5. Compute derived metrics when denominators are non-zero:

- IPC = `instructions / cycles`
- Cache miss rate = `cache-misses / cache-references`
- Branch mispredict rate = `branch-misses / branches`

6. Pair baseline and contender rows by benchmark name, compute deltas, and flag anomalies where timing direction conflicts with key counter direction.

7. Emit a canonical summary table for downstream consumers:

```markdown
| Benchmark | IPC (base -> new) | Cache Miss Rate (base -> new) | Branch Mispredict (base -> new) | Anomaly? |
|---|---:|---:|---:|---|
```

## Retry and Timeout Policy

1. Run benchmarks sequentially; do not background with `nohup`/`&`.
2. If a run times out, narrow filter and retry once.
3. Maximum retries per benchmark group: 1.
4. If still failing, emit blocked/partial findings instead of repeated attempts.

Apply this policy to both timing and counter runs.

## Step 4 — Record findings

Capture and return:

- compare.py output (terminal transcript or redirected file)
- effective filter used
- timing JSON artifacts for baseline and contender
- `counters_available` (`true`/`false`)
- if `counters_available=false`, a reason string (unsupported OS, missing libpfm, perf permission denied, preflight failure)
- if counters are available: counter JSON artifacts, derived metrics table, and anomaly list

## Best Practices / Guardrails

1. **Release only**: never compare Debug binaries.
2. **Short hash suffixes**: keep build dirs isolated per revision (example: `bench_<short-hash>`).
3. **Same host, same conditions**: do not compare across different machines or power profiles.
4. **Filter from analysis**: use `benchmarks-affected` output instead of hand-crafted filters whenever possible.
5. **Pin frequency**: for stable numbers, follow benchmark skill guidance on CPU governor.
6. **Counter collection is optional and Linux-only**: when unavailable, return timing-only outputs with `counters_available=false`.
7. **Always preflight counters**: do not run full counter collection if preflight fails.
8. **Keep build types separated**: timing uses `benchmarks-all_*` Release builds; counters use `benchmarks-diagnostic_*` RelWithDebInfo builds; never Debug.
