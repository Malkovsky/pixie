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
- Optional argument `--no-counters`: disable hardware-counter collection.

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

## Step 2 - Run timing and hardware-counter comparison via skill (single source of truth)

Use `benchmarks-compare-revisions` as the single source of truth for revision builds, benchmark scope, compare.py flow, retry policy, and guardrails.

Pass-through inputs:
- Baseline ref/hash from Step 1.
- Contender ref/hash from Step 1.
- Optional `--filter` override.
- Counter mode: default on (`COLLECT_COUNTERS=1`) on Linux, disabled when `--no-counters` is provided.

Consume outputs from `benchmarks-compare-revisions`:
- Baseline and contender benchmark JSON artifacts.
- compare.py output per binary.
- Effective filter used.
- Scope metadata from `benchmarks-affected` (`affected_benchmark_targets`, `affected_benchmarks`) when available.
- `counters_available` status and, when unavailable, explicit reason.
- Baseline and contender counter JSON artifacts (when available).
- Derived counter metrics per benchmark (IPC, cache miss rate, branch mispredict rate).
- Counter anomaly list and ready-to-embed counter summary table.

Execution guardrails:
- Run benchmarks sequentially.
- No background jobs (`nohup`, `&`).
- Use Release timing builds only.
- If timing comparison fails, return blocked verdict with exact failure points.

## Step 3 - Consume delegated hardware-counter outputs

Hardware-counter collection is delegated to `benchmarks-compare-revisions`.

Pass-through inputs:
- `COLLECT_COUNTERS=1` by default on Linux (unless `--no-counters` is provided).
- Same baseline/contender refs and effective filter used in Step 2.

Consume outputs:
- Counter preflight result.
- Counter JSON artifacts for both revisions.
- Derived metrics (IPC, cache miss rate, branch mispredict rate).
- Anomaly list and counter summary table for report embedding.

If counters are unavailable (`counters_available=false`), continue with timing-only review and explicitly mark profiling as unavailable in the report.

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
- Use skill-provided hardware counter summary and anomaly list to explain major timing changes.
- Do not recompute derived counter metrics in this command.

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
- If only profiling fails (`counters_available=false` from delegated skill output), continue with timing-based verdict and explicitly list profiling limitation.
- If JSON output is invalid/truncated, discard it and rerun that benchmark command once with tighter filter and explicit output redirection.
