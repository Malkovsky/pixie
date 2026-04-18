---
description: Benchmark-driven PR performance review versus target branch
---

# Perf Review Workflow

You are performing a performance review for the current PR branch.

Non-negotiable requirements:
1. Benchmark timing plus profiling data is the highest-priority judgment tool.
2. Compare source branch versus target branch and report relevant benchmark metric changes.
3. Provide analysis and a final verdict: does the PR improve performance or not.

## Inputs

- Optional argument `--target <branch>`: target branch override.
- Optional argument `--filter <regex>`: benchmark filter regex.

If arguments are omitted:
- Default target branch to PR base branch from `gh pr view --json baseRefName` when available.
- Fall back target branch to `main`.

Filter handling:
- If `--filter` is provided, pass it through.
- Else use the filter produced by `benchmarks-affected` through `benchmarks-compare-revisions`.
- If no filter can be derived, run conservative full-binary compare for impacted binaries.

## Step 1 - Resolve branches and hashes

1. Resolve contender from current checkout (`HEAD`) and compute short hash.
2. Resolve baseline branch using precedence: `--target` -> PR base from `gh pr view --json baseRefName` -> `main`.
3. Resolve baseline short hash.
4. Print branch/hash mapping before benchmark execution.

## Step 2 - Run timing comparison via skill (single source of truth)

Use `benchmarks-compare-revisions` as the single source of truth for revision builds, benchmark scope, compare.py flow, retry policy, and guardrails.

Do not duplicate or override its internal build/run steps in this command.

Pass-through inputs:
- Baseline ref/hash from Step 1.
- Contender ref/hash from Step 1.
- Optional `--filter` override.

Consume outputs from `benchmarks-compare-revisions`:
- Baseline and contender benchmark JSON artifacts.
- compare.py output per binary.
- Effective filter used.
- Scope metadata from `benchmarks-affected` (`affected_benchmark_targets`, `affected_benchmarks`) when available.

Execution guardrails:
- Run benchmarks sequentially.
- No background jobs (`nohup`, `&`).
- Use Release timing builds only.
- If timing comparison fails, return blocked verdict with exact failure points.

## Step 3 - Collect hardware counter profiles (Linux only, optional)

Run a preflight first to avoid wasted attempts:
1. Execute one tiny benchmark with perf counters (e.g. one benchmark case) and inspect output for counter availability.
2. If output includes warnings like `Failed to get a file descriptor for performance counter`, mark counters unavailable and skip counter collection.

If preflight passes and Linux profiling build is available, run both baseline and contender diagnostic binaries with counter output:

- `--benchmark_counters_tabular=true`
- `--benchmark_format=json`
- `--benchmark_out=<output.json>`
- Include `--benchmark_filter` when provided.

Collect and compare at least these counter families when present:
- `instructions`, `cycles` (for IPC)
- `cache-misses`, `cache-references` (cache miss rate)
- `branch-misses`, `branches` (branch mispredict rate)
- `L1-dcache-load-misses` (L1 data cache pressure)

Compute derived metrics when denominators are non-zero:
- IPC = instructions / cycles
- Cache miss rate = cache-misses / cache-references
- Branch mispredict rate = branch-misses / branches

If profiling is unavailable (non-Linux, libpfm missing, or perf permissions blocked), continue with timing-only review and explicitly mark profiling as unavailable in the report.

## Step 4 - Analyze timing and counter data

Timing classification per benchmark entry:
- Improvement: time delta < -5%
- Regression: time delta > +5%
- Neutral: between -5% and +5%

Aggregate per binary:
- Number of improvements/regressions/neutral
- Net average percentage change
- Largest regression and largest improvement

Counter correlation:
- Use hardware counters to explain major timing changes.
- Flag anomalies (timing improves while key counters degrade, or opposite).

Judgment priority:
- Base verdict primarily on benchmark timing comparison.
- Use counter data as explanatory evidence and confidence signal.

Noise-control expectations:
- Include at least one control benchmark family expected to be unaffected by the code change.
- Treat isolated swings without pattern as noise unless reproduced across related sizes/fill ratios.

## Step 5 - Produce final markdown report

Return a structured markdown report with this shape:

```markdown
## Performance Review: <contender_branch> vs <baseline_branch>

### Configuration
- Baseline: <branch> (<hash>)
- Contender: <branch> (<hash>)
- Platform: <os/arch>
- Benchmarks run: <binary list>
- Filter: <regex or none>
- Hardware counters: available / unavailable

### Timing Summary
| Binary | Improvements | Regressions | Neutral | Net Change |
|---|---:|---:|---:|---:|
| ... | N | N | N | +/-X% |

### Detailed Timing Results
<Annotated compare.py outputs by binary>

### Hardware Counter Profile (if available)
| Benchmark | IPC (base->new) | Cache Miss Rate (base->new) | Branch Mispredict (base->new) |
|---|---:|---:|---:|
| ... | X.XX -> Y.YY | A.A% -> B.B% | C.C% -> D.D% |

### Key Findings
- <Most important regressions/improvements>
- <Counter-based explanations for key timing shifts>

### Verdict
**[IMPROVES PERFORMANCE | REGRESSES PERFORMANCE | NO SIGNIFICANT CHANGE]**

<1-2 sentence justification grounded in benchmark metrics, with profiling context if available>
```

Verdict rules:
- `IMPROVES PERFORMANCE`: improvements outnumber regressions, no severe regression (>10%), and net average change is favorable.
- `REGRESSES PERFORMANCE`: any severe regression (>10%) or regressions dominate with net unfavorable average.
- `NO SIGNIFICANT CHANGE`: mostly neutral changes or mixed results that approximately cancel out.

## Failure Handling

- If required builds fail or timing comparison cannot run, output a blocked review with exact failure points and no misleading verdict.
- If only profiling fails, continue with timing-based verdict and explicitly list profiling limitation.
- If JSON output is invalid/truncated, discard it and rerun that benchmark command once with tighter filter and explicit output redirection.
