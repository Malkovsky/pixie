---
description: Scan current branch and report impacted benchmark targets/functions.
---

# Benchmarks Affected

Identify which benchmark binaries and benchmark functions are affected by changes on the current branch.

Use the `benchmarks-affected` skill as the single source of truth for workflow details and guardrails.
Do not duplicate or override the skill instructions in this command.

## Inputs

- Optional `--baseline <ref>` (default: `main`)
- Optional `--compile-commands <path>`
- Optional `--no-include-working-tree`
- Optional `--format <text|json>` (default: `text`)

## Workflow

1. Execute the `benchmarks-affected` skill workflow.
2. Pass through command inputs to the analyzer invocation defined by the skill.
3. Report results with these sections:
   - Changed files
   - Affected benchmark targets
   - Affected benchmark functions
   - Suggested `--benchmark_filter` regex
   - Any warnings/failures

## Output rules

1. If `affected_benchmarks` is non-empty, prioritize those names.
2. If `affected_benchmarks` is empty but benchmark targets are affected, mark result as partial and include target-level impact.
3. Do not run full benchmark suites in this command; this command is for impact discovery only.
