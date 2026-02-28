---
name: benchmarks-compare-revisions
description: Compare benchmark performance between two git revisions by building separate binaries with short-hash build suffixes and using Google Benchmark compare.py.
---

# Benchmarks Compare Revisions Skill

Use this skill to compare performance between two git revisions. It focuses on the compare workflow and relies on the existing benchmarks skill for build/run details (see .opencode/skills/benchmarks/SKILL.md).

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
# Follow .opencode/skills/benchmarks/SKILL.md build instructions with this suffix

# Contender
BUILD_SUFFIX=bench_${CONTENDER}
git checkout ${CONTENDER}
# Follow .opencode/skills/benchmarks/SKILL.md build instructions with this suffix
```

## Step 2 — Compare using compare.py

Use Google Benchmark’s compare.py to run both binaries and compute a statistical comparison.

Locate compare.py from the Google Benchmark dependency (installed under the build tree):
```bash
COMPARE_PY=build/benchmarks-all_bench_${BASELINE}/_deps/googlebenchmark-src/tools/compare.py
```

Run the comparison (benchmarks mode):
```bash
python3 ${COMPARE_PY} benchmarks \
  build/benchmarks-all_bench_${BASELINE}/benchmarks \
  build/benchmarks-all_bench_${CONTENDER}/benchmarks
```

### Optional: filter to reduce noise

Pass benchmark options after the binaries so compare.py forwards them:
```bash
python3 ${COMPARE_PY} benchmarks \
  build/benchmarks-all_bench_${BASELINE}/benchmarks \
  build/benchmarks-all_bench_${CONTENDER}/benchmarks \
  --benchmark_filter="BM_Rank"
```

## Step 3 — Record findings

Capture the compare.py output (terminal transcript or redirected file) and note any statistically significant regressions or wins.

## Best Practices / Guardrails

1. **Release only**: never compare Debug binaries.
2. **Short hash suffixes**: keep build dirs isolated per revision (example: `bench_<short-hash>`).
3. **Same host, same conditions**: do not compare across different machines or power profiles.
4. **Filter for focus**: narrow `--benchmark_filter` to the changed area.
5. **Pin frequency**: for stable numbers, follow benchmark skill guidance on CPU governor.
