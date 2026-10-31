# Pixie profiling

Use `agentic/cpp/skills/profiling/SKILL.md` for the workflow and scripts.
Build with the `benchmarks-profile` preset for optimized code and debug info;
it also enables diagnostics and libpfm. Check the actual flags and select the
desired target rather than rebuilding every benchmark.

The rank/select harness supports a query window through
`src/benchmarks/perf_query_window.h`. The capture helper supplies its variables:

```sh
python3 agentic/cpp/skills/profiling/scripts/record_perf.py \
  --output-dir /tmp/rank-select-profile --cpu 0 \
  --control-prefix PIXIE_RANK_SELECT_PERF -- \
  ./build/benchmarks-profile/rank_select_benchmarks \
  --benchmark_filter='^rank_select_compressed_select1_50/N:1073741824(/|$)' \
  --benchmark_repetitions=1 \
  --benchmark_out=/tmp/rank-select-profile/workload.json \
  --benchmark_out_format=json
```

Check the row exists with `--benchmark_list_tests=true`. The control variables
disable query warmup, and `perf_queries`/`perf_phases` report cumulative captured
work. Calibration can still create multiple phases: inspect those counters and
use the established fixed-work runner when exact work is needed. A function
filter in the viewer alone does not establish a query-only capture.

On hosts where the PATH `perf` wrapper cannot find tools for the running kernel,
pass a verified native executable via `--perf-bin`; do not silently substitute
a software event or interpret unavailable PMU events as zero.

Generate the report only on user request:

```sh
python3 agentic/cpp/skills/profiling/scripts/export_html.py \
  --perf-data /tmp/rank-select-profile/perf.data \
  --binary ./build/benchmarks-profile/rank_select_benchmarks \
  --source-root . --metadata /tmp/rank-select-profile/capture.json \
  --output-dir /tmp/rank-select-profile/html --title 'Rank/select profile'
```

Keep captures and generated HTML outside tracked files. In select profiles,
cycle samples on a consuming popcount can reflect a preceding data load's stall;
inspect annotated assembly and counters before attributing the cost to popcount.

## Report presentation

The approved Pixie report format is the light, interactive viewer exemplified
by `build/profiles/select1-34-50/index.html`. The shared profiling skill's bundled
viewer already carries that layout; use it as the primary report rather than
introducing a separate dashboard. The old build artifact is a visual reference,
not a required input for future exports.

Use Pixie branding and workload-specific titles where appropriate, while keeping
the flame graph, frame/source panel, hotspots and controls in their standard
order. Additional analysis or heatmaps belong in companion views. Honor the
user's requested output directory and ensure the delivered entry-point link
opens the interactive viewer.
