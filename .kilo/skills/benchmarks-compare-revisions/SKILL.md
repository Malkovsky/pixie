---
name: benchmarks-compare-revisions
description: Compare benchmark performance between two git revisions by building separate binaries with short-hash build suffixes and using Google Benchmark compare.py.
---

# Benchmarks Compare Revisions Skill

Use this skill to compare performance between two git revisions. It focuses on the compare workflow and relies on the existing benchmarks skill for build/run details (see .kilo/skills/benchmarks/SKILL.md).

## Goal

Build two separate benchmark binaries using short commit hashes as build suffixes, then compare results with Google Benchmark compare.py.

## Step 0 — Choose revisions and short hashes

Pick a baseline and a contender revision. Use short commit hashes to suffix build directories so builds do not collide.

Example:
```bash
BASELINE=abc1234
CONTENDER=def5678
```

## Step 1 — Build both revisions (Release only)

Use the existing benchmarks skill build steps, but set the build suffix to include the short hash for each revision:

```bash
# Baseline
BUILD_SUFFIX=bench_${BASELINE}
git checkout ${BASELINE}
# Follow .kilo/skills/benchmarks/SKILL.md build instructions with this suffix

# Contender
BUILD_SUFFIX=bench_${CONTENDER}
git checkout ${CONTENDER}
# Follow .kilo/skills/benchmarks/SKILL.md build instructions with this suffix
```

## Step 2 — Compare using compare.py

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

### Optional: filter to reduce noise and runtime

Pass filter when generating JSON files:
```bash
FILTER="BM_Rank"
build/benchmarks-all_bench_${BASELINE}/benchmarks --benchmark_filter="${FILTER}" --benchmark_report_aggregates_only=true --benchmark_display_aggregates_only=true ...
build/benchmarks-all_bench_${CONTENDER}/benchmarks --benchmark_filter="${FILTER}" --benchmark_report_aggregates_only=true --benchmark_display_aggregates_only=true ...
```

## Retry and Timeout Policy

1. Run benchmarks sequentially; do not background with `nohup`/`&`.
2. If a run times out, narrow filter and retry once.
3. Maximum retries per benchmark group: 1.
4. If still failing, emit blocked/partial findings instead of repeated attempts.

## Step 3 — Record findings

Capture the compare.py output (terminal transcript or redirected file) and note any statistically significant regressions or wins.

## Best Practices / Guardrails

1. **Release only**: never compare Debug binaries.
2. **Short hash suffixes**: keep build dirs isolated per revision (example: `bench_<short-hash>`).
3. **Same host, same conditions**: do not compare across different machines or power profiles.
4. **Filter for focus**: narrow `--benchmark_filter` to the changed area.
5. **Pin frequency**: for stable numbers, follow benchmark skill guidance on CPU governor.
