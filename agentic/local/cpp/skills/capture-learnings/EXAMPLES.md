# Pixie Capture Learnings Examples

Use these examples with
`agentic/cpp/skills/capture-learnings/SKILL.md`. They show how to route
learnings from Pixie's succinct-data-structure research without putting Pixie
paths or benchmark names into the shared skill.

## Before Requested Compaction

When the user requests compaction during an RMQ research session, invoke this
skill before writing the handoff.

Persist durable findings first. For example, a corrected benchmark-filter rule
belongs in `agentic/local/cpp/skills/benchmarks/EXAMPLES.md`, and a reusable
shape-forcing test strategy may belong in `AGENTS.md` when the evidence supports
it.

Then keep transient state only in the compaction checkpoint, such as:

- the current branch and commit
- staged, unstaged, and untracked work
- tests already run and tests still required
- unresolved PR feedback
- the next concrete action

Do not put those transient identifiers into the capture skill or benchmark
guidance.

## RMQ Optimization Session

Suppose an RMQ experiment produces:

- a promoted same-leaf query optimization
- pinned benchmark results across several range widths
- a discarded variant that wins only one synthetic row
- a corrected Google Benchmark filter
- new deterministic tests for a border traversal

Capture the results as follows:

| Learning | Evidence | Destination |
|---|---|---|
| The promoted query path is faster for its intended workload | Repeated pinned Release measurements plus passing differential tests | Current benchmark snapshot near the implementation or benchmark |
| The discarded variant regresses realistic neighboring widths | Saved comparison results with matching build and workload | Experiment log only when it prevents repeating the idea |
| Registered benchmark names append configuration suffixes | `--benchmark_list_tests=true` output and previously empty filtered runs | `agentic/local/cpp/skills/benchmarks/EXAMPLES.md` |
| Border correction needs a shape-forcing test | A reproduced failure or previously uncovered traversal branch | Test plus the testing guidance in `AGENTS.md` when broadly reusable |

Do not copy every intermediate timing table into `include/pixie/rmq.h`. That
header should contain the accepted current snapshot, while exploratory history
belongs in an experiment log or remains temporary.

## Coverage Learning

When a coverage investigation shows that random RMQ tests miss a specific
HybridBTree selector path:

1. Add a deterministic input that forces the path.
2. Verify the public result against a naive implementation.
3. Run `./scripts/coverage_report.sh` before trusting the aggregate.
4. Update `AGENTS.md` only when the shape-forcing strategy is useful for future
   RMQ or RmM work, not merely because one line changed coverage.

Do not preserve private defensive branches or impossible allocation failures as
coverage objectives.

## SIMD or Platform Result

For a SIMD experiment, keep the conclusion bounded by the tested configuration.
For example:

```text
Claim: Candidate X reduces query time for 128-bit excess minima on the tested
AVX2 path and input distribution.

Evidence: pinned Release benchmark comparison, differential tests against the
scalar implementation, and a successful DISABLE_AVX512 build.

Scope: Pixie's AVX2 implementation for the measured workload; not a general
claim about all CPUs or range distributions.
```

Route a reusable experimental procedure to the shared optimization skill. Keep
the concrete target, compiler flags, CPU, benchmark rows, and measured values in
Pixie's benchmark snapshot or experiment log.

## Proposal Format

Before updating durable guidance, present proposals in this compact form:

```text
Learning: Shape-forcing RMQ tests reveal traversal bugs that random tests can
miss.
Evidence: <test or coverage artifact>
Scope: Pixie RMQ/RmM implementations.
Destination: AGENTS.md testing guidance.
Action: Extend the existing rule; do not create a parallel testing document.
```

Reject a candidate learning when it is already obvious from the code, supported
only by one noisy timing run, tied to a temporary build directory, or unlikely
to change a future decision.
