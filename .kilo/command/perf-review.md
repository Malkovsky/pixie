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
- Default filter to empty (run full selected benchmark suites).

## Step 1 - Resolve Branches and Revisions

1. Identify contender branch and hash:
   - Contender branch: current checked-out branch (or `HEAD` if detached).
   - Contender hash: `git rev-parse --short HEAD`.
2. Identify baseline branch:
   - Use `--target` if provided.
   - Else use PR base branch from GitHub CLI when available.
   - Else use `main`.
3. Resolve baseline hash with `git rev-parse --short <baseline-ref>`.
4. Print branch and hash mapping before running benchmarks.

## Step 2 - Select Relevant Benchmark Binaries

Inspect changed files with:

`git diff --name-only <baseline-ref>...HEAD`

Map file paths to benchmark binaries:

| Changed path pattern | Benchmark binary | Coverage |
|---|---|---|
| `include/bit_vector*`, `include/interleaved*` | `benchmarks` | BitVector rank/select |
| `include/rmm*` | `bench_rmm` | RmM tree operations |
| `include/louds*` | `louds_tree_benchmarks` | LOUDS traversal |
| `include/simd*`, `include/aligned*` | `alignment_comparison` | SIMD and alignment |
| `include/misc/*` | all relevant | Differential helpers |
| `CMakeLists.txt`, benchmark infra, broad/unknown changes | all benchmarks | Conservative full run |

Available benchmark binaries:
- `benchmarks`
- `bench_rmm`
- `bench_rmm_sdsl`
- `louds_tree_benchmarks`
- `alignment_comparison`

If the mapping is ambiguous, run all benchmark binaries.
If `--filter` is provided, pass it through as `--benchmark_filter`.
Print selected binaries and why they were selected.

## Step 3 - Build Both Revisions (Timing and Profiling Builds)

Use isolated build directories per short hash.

1. Capture original ref (`git rev-parse --abbrev-ref HEAD` or detached `HEAD`).
2. If worktree is dirty, stash safely with untracked files:
   - `git stash push -u -m "perf-review-auto-stash"`
3. Build baseline revision:
   - `git checkout <baseline-hash-or-ref>`
   - Timing build (required):
     - `cmake -B build/benchmarks-all_bench_<baseline_hash> -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON`
     - `cmake --build build/benchmarks-all_bench_<baseline_hash> --config Release -j`
   - Profiling build (Linux only, recommended):
     - `cmake -B build/benchmarks-diagnostic_bench_<baseline_hash> -DCMAKE_BUILD_TYPE=RelWithDebInfo -DPIXIE_BENCHMARKS=ON -DBENCHMARK_ENABLE_LIBPFM=ON -DPIXIE_DIAGNOSTICS=ON`
     - `cmake --build build/benchmarks-diagnostic_bench_<baseline_hash> --config RelWithDebInfo -j`
4. Build contender revision:
   - `git checkout <contender-hash-or-original-ref>`
   - Repeat timing and profiling build with contender hash suffix.
5. Restore original ref and restore stashed state if a stash was created.

Critical guardrails:
- Never use Debug binaries for timing review.
- Timing comparisons must use `benchmarks-all` Release builds.
- Profiling counters should use `benchmarks-diagnostic` RelWithDebInfo builds.

## Step 4 - Resolve Binary Paths

Support both generator layouts:

- Multi-config: `build/<dir>/Release/<binary>` or `build/<dir>/RelWithDebInfo/<binary>`
- Single-config: `build/<dir>/<binary>`

For each needed binary, detect the existing executable path before running.
If a required binary is missing, report failure and stop with a blocked verdict.

## Step 5 - Run Timing Comparison (Primary Judgment)

Locate compare script from baseline timing build:

`build/benchmarks-all_bench_<baseline_hash>/_deps/googlebenchmark-src/tools/compare.py`

For each selected benchmark binary, run:

`python3 <compare.py> benchmarks <baseline_binary> <contender_binary> [--benchmark_filter="<filter>"]`

Capture full output for each binary and keep it for report details.

## Step 6 - Collect Hardware Counter Profiles (Linux Only)

If Linux profiling build is available, run both baseline and contender diagnostic binaries with counter output:

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

If profiling is unavailable (non-Linux or libpfm not available), continue with timing-only review and explicitly mark profiling as unavailable in the report.

## Step 7 - Analyze Timing and Counter Data

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

## Step 8 - Produce Final Markdown Report

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
