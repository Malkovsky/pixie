---
name: optimization-experiment
description: Run iterative C++ optimization experiments for a target function or class by adding same-API experimental variants, validating correctness, benchmarking, comparing results, and deciding whether to promote a faster implementation.
---

# Optimization Experiment Skill

Use this skill when a user wants to improve performance of a specific C++
function, class, algorithm, or hot path through benchmark-driven experiments.

This workflow depends on:

1. `../benchmarks/SKILL.md` for Google Benchmark build/run commands, JSON output,
   hardware counters, pinning, and perf profiling.
2. `../benchmarks-affected/SKILL.md` when changes need an affected benchmark
   scope.
3. `../benchmarks-compare-revisions/SKILL.md` when comparing committed
   revisions.

## Goal

Iterate from a production implementation to one or more experimental
implementations, prove semantic equivalence, measure the impact, and decide
whether a candidate is worth promoting.

The standard loop is:

```text
target -> benchmark baseline -> experimental same-API variant
       -> correctness check -> benchmark compare -> keep / revise / discard
```

Stop when a candidate is clearly better on the intended workload without
correctness or maintenance regressions, or when the remaining ideas are too weak
to justify more iteration.

## Step 1 - Identify the Target and Contract

Start from the requested function/class and inspect the real implementation.
Record:

- public signature/API and call sites that must not change
- input domains, invalid-input behavior, boundary conditions, and tie-breaking
- compile-time feature gates such as SIMD flags or platform-specific paths
- existing tests and reference implementations
- existing benchmarks that should move if the optimization succeeds

Do not optimize before the contract is clear. If behavior is ambiguous, add or
find tests before changing implementation.

## Step 2 - Establish Benchmark Coverage

Find benchmark rows that directly exercise the target. Prefer narrow benchmark
filters over full-suite runs during iteration.

If coverage is missing or too broad, add focused benchmark cases before adding
the optimized implementation. Include cases for:

- the expected common path
- boundary and alignment-sensitive paths
- short, medium, and long ranges or sizes when width matters
- random or mixed workloads when real calls are not fixed-shape
- current production behavior and each experimental variant

Capture a baseline JSON before implementation changes:

```bash
BENCH_CPU=${BENCH_CPU:-0}
taskset -c "${BENCH_CPU}" <benchmark-binary> \
  --benchmark_filter="${FILTER}" \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_out=/tmp/<target>_baseline.json \
  --benchmark_out_format=json
```

Use `../benchmarks/SKILL.md` for exact build directories, Release versus
diagnostic builds, hardware-counter setup, and retry policy.

## Step 3 - Add Experimental Same-API Variants

Add candidate implementations beside production code in an experimental area,
namespace, header, or benchmark-local adapter that is already consistent with
the repository.

Rules:

- keep the callable signature/API identical to production where practical
- preserve public semantics exactly, including invalid inputs and tie-breaking
- keep production callers unchanged during experiments
- make variants benchmark-selectable by name
- avoid unrelated refactors while measuring
- keep losing variants only when they document useful evidence or support future
  comparison

For C++ libraries with feature-gated implementations, provide correct fallbacks
for unsupported targets or compile configurations.

## Step 4 - Validate Correctness Before Timing

Run relevant tests before trusting benchmark numbers. Add tests when the
experimental implementation introduces new risk.

Prefer:

- fixed edge cases for boundaries, empty/sentinel behavior, and exact ties
- randomized differential tests against a scalar or naive reference
- tests for feature-gated fallback builds when the code has SIMD or platform
  branches
- targeted regression tests for any bug found during benchmarking

Do not compare performance for a candidate that has not passed the correctness
checks for the same semantics as production.

## Step 5 - Benchmark and Compare

Run timing benchmarks from Release builds. Save JSON for every meaningful
baseline and candidate.

Use diagnostic builds with hardware counters when timing changes need
explanation:

- cycles and instructions for core execution cost
- cache counters for memory behavior
- branch counters when early exits or dispatch logic are involved

Compare both absolute timings and relative deltas. Watch for cases where a
candidate wins the cherry-picked row but regresses neighboring or realistic
workloads.

When results are noisy:

- pin to a CPU with `taskset` when available
- increase repetitions or minimum benchmark time
- rerun the narrow benchmark filter once
- avoid changing benchmark scope between baseline and candidate

## Step 6 - Iterate Deliberately

For each candidate, decide one of:

- **Promote**: repeatedly faster on intended rows, no important regressions,
  correct and maintainable.
- **Keep experimental**: interesting or workload-specific, but not production
  ready.
- **Discard**: slower, too complex, too narrow, or semantically risky.

Use benchmark data to choose the next idea. Examples:

- higher instruction count suggests fewer operations or simpler dispatch
- lower instructions but higher cycles suggests stalls, memory, or dependency
  chains
- short-range regressions suggest a narrower dispatch condition
- alignment-sensitive rows suggest splitting aligned and unaligned paths

When no idea wins convincingly, document the best result and stop rather than
overfitting.

## Step 7 - Finalize the Result

If promoting a candidate to production:

- keep the public API unchanged unless the user explicitly requested otherwise
- keep or update tests that protect the optimized behavior
- remove accidental benchmark-only scaffolding from production code
- preserve experimental variants only when useful for future research

If leaving work experimental:

- add a short note near the experimental code with benchmark date, command, and
  the relevant table or JSON artifact path
- clearly state that production callers do not use the experimental variant
- explain which workload the variant helps and where it loses

The final response should include:

- what changed
- correctness checks run
- benchmark command or JSON artifacts
- concise result table
- recommendation: promote, keep experimenting, or stop

## Guardrails

1. Benchmark before optimizing; otherwise there is no trustworthy baseline.
2. Never change semantics to win a benchmark.
3. Never compare Debug timings.
4. Keep production and experimental code paths distinguishable.
5. Prefer focused benchmark filters during iteration, then broaden before
   promotion.
6. Treat hardware counters as explanatory data, not a replacement for timing.
7. Record enough benchmark context that future agents do not confuse
   experimental wins with production behavior.
