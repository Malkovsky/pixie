# Pixie Optimization Experiment Examples

Use these notes with `agentic/cpp/skills/optimization-experiment/SKILL.md` for
Pixie-specific optimization work.

## Example: `excess_min_128`

The `excess_min_128` experiment is a good template for future hot-path work:

1. **Lock semantics first**: preserve inclusive prefix range `[left, right]`,
   prefix offsets `[0, 128]`, first-min tie breaking, and invalid-range sentinel
   behavior.
2. **Add focused benchmarks before trusting results**: include full-range,
   RMQ-style `[0, 127]`, aligned short ranges, non-aligned short ranges,
   cross-word ranges, point ranges, and a reproducible random-range benchmark.
3. **Keep candidates same-API**: experimental variants should accept the same
   parameters and return the same result type as production, so benchmarks can
   swap implementations without adapter logic.
4. **Use experimental space for ideas**: keep exploratory variants in
   `include/pixie/experimental/` unless they are clearly production-ready.
5. **Promote narrowly**: when a variant wins only under a specific condition,
   promote only that condition. For `excess_min_128`, byte-LUT was useful only
   for byte-aligned short ranges, so the production fallback was narrowed instead
   of applied to every short range.
6. **Mark unlikely cold dispatches**: if a promoted optimization is for a narrow
   special case, mark the branch unlikely when that matches expected workload.
7. **Persist benchmark evidence in the repository**: add a top-of-file
   benchmark note in `/** ... */` or `/* ... */` form to the experimental
   header, production header, benchmark file, or experiment log. Include metric
   units, fixed decimal precision, JSON paths when available, and visible
   separators between logical result blocks.

## Benchmark Pattern

Prefer a diagnostic run that can explain timing changes, not just report them:

```bash
taskset -c 0 build/benchmarks-diagnostic_local/excess_positions_benchmarks \
  --benchmark_filter='BM_ExcessMin128(_|/|$)' \
  --benchmark_repetitions=3 \
  --benchmark_perf_counters=CYCLES,INSTRUCTIONS,CACHE-MISSES \
  --benchmark_counters_tabular=true \
  --benchmark_out=/tmp/pixie_excess_min_counters.json \
  --benchmark_out_format=json
```

For final comparison tables, include at least:

- CPU time
- cycles
- instructions
- cache misses when collected
- a random or mixed workload row
- the rows that justified any narrowed production dispatch condition

## Documentation Pattern

For experimental algorithms, use Doxygen block comments:

```cpp
/**
 * @brief Short algorithm name and purpose.
 *
 * @details Workflow:
 *
 *   input -> transform -> candidate result -> final reduction
 *
 * Explain the non-obvious tradeoff, such as avoiding lane crossing, reducing
 * loop iterations, or paying scalar boundary work.
 */
```

Keep long benchmark tables as a top-level `/** ... */` note when they describe
the whole experimental header rather than one symbol. Do not leave measurements
only in chat or `/tmp`; persist a top-of-file snapshot before finishing.

## Non-Obvious Lessons From This Session

- A faster micro-kernel can regress neighboring ranges if dispatch is too broad;
  benchmark aligned and non-aligned cases separately.
- Random-range benchmarks are useful as a sanity check against overfitting fixed
  rows.
- Hardware counters helped distinguish lower-level cost: cycles and instructions
  were more actionable than cache misses for the `excess_min_128` candidates.
- Production and experimental code must stay visually distinct. A benchmark win
  in `pixie::experimental` is only a candidate, not a caller-visible change.
- When result names encode width unnecessarily, prefer the stable semantic name
  (`ExcessResult`) over width-specific detail (`ExcessResult128`) unless the API
  truly needs multiple widths.
