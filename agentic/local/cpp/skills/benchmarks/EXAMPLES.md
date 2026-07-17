# Pixie Benchmark Examples

Use these notes with `agentic/cpp/skills/benchmarks/SKILL.md` for Pixie's local
benchmark binaries and reporting conventions.

## Rank/Select Benchmarking

The primary rank/select benchmark binary is:

```bash
./build/benchmarks/rank_select_benchmarks
```

It registers construction (`rank_select_build_both_*`) and query rows
(`rank1`, `rank0`, `select1`, `select0`) for independently generated packed bit
vectors with expected 12.5%, 50%, and 87.5% one-fill. Query rows build the
support and query pool before timing. The support owns only index metadata, so
the external bit sequence is reported separately as `input_bytes`. Every row
constructs `RankSelectSupport` with `SelectSupport::kBoth`, so its memory
counters include select1 and select0 metadata.

Use a pinned Release run over the documented table sizes:

```bash
taskset -c 0 ./build/benchmarks/rank_select_benchmarks \
  --benchmark_filter='^rank_select_(rank1|rank0|select1|select0|build_both)_(12p5|50|87p5)/N:(1024|67108864|17179869184)(/|$)' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true
```

Rows append Google Benchmark's `min_time` and warmup suffixes, so filters must
accept `(/|$)` after the size. Build rows report CPU milliseconds; query rows
report CPU nanoseconds. Both expose `aux_bytes`, `aux_mib`, and
`aux_bits_per_input_bit` for metadata owned by the support.

### FNBP Source Rows

`rank_select_fnbp_*` generates the Ferrada-Navarro balanced-parentheses source
used by the Cartesian RMQ backends. It covers `2^10` through `2^30` bits,
including a large-source case with 16,384 rank superblocks.

```bash
taskset -c 0 ./build/benchmarks/rank_select_benchmarks \
  --benchmark_filter='^rank_select_fnbp_(build_both|rank1|rank0|select1|select0)/N:(1024|16384|262144|4194304|67108864|1073741824)(/|$)' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true
```

Do not commit raw benchmark JSON or ad hoc result tables for stable
implementations. A user-requested current snapshot may live in the relevant
implementation catalog header; keep only rounded results and reproducible
methodology, align tables visually in source, and exclude local JSON, failed
probes, and before/after experiment history. Persist other results only for an
explicitly experimental implementation that has a registered benchmark.

## RMQ Benchmark Tables

The primary RMQ benchmark binary is usually:

```bash
./build/benchmark-all-backends/rmq_benchmarks
```

When the table in `include/pixie/rmq.h` needs updating, use pinned Release
runs and keep the header as the current snapshot only. Do not store local JSON
paths, failed probes, or "before/after" experiment history in that header.

For query rows:

```bash
taskset -c 0 ./build/benchmark-all-backends/rmq_benchmarks \
  --benchmark_filter='^rmq_(hybrid_btree|cartesian_hybrid_btree|cartesian_btree)/(1024|16384|262144|4194304|16777216|67108864)/' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/rmq_query.json
```

For build and memory rows:

```bash
taskset -c 0 ./build/benchmark-all-backends/rmq_benchmarks \
  --benchmark_filter='^rmq_build_(hybrid_btree|cartesian_hybrid_btree|cartesian_btree)/(1024|16384|262144|4194304|16777216|67108864)(/|$)' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/rmq_build_memory.json
```

Google Benchmark appends suffixes such as `/min_time:2.000/min_warmup_time:0.200`
to registered benchmark names. Exact filters that end at the size often miss all
rows; use a suffix-tolerant tail such as `(/|$)` or inspect names first:

```bash
./build/benchmark-all-backends/rmq_benchmarks --benchmark_list_tests=true | rg cartesian_rmm
```

## Memory Counters

RMQ benchmarks expose owned auxiliary memory through user counters when the
implementation provides `memory_usage_bytes()`:

- `aux_mib`: owned auxiliary memory in MiB.
- `aux_bits_per_value`: owned auxiliary memory normalized by indexed value
  count.
- `input_bytes`: external value array size; this is not counted as owned RMQ
  auxiliary memory.

Build benchmark rows are convenient for refreshing memory tables because they
construct the structure and emit the memory counters even when query timings are
not needed. For example, to measure `CartesianRmM` memory through `2^30`:

```bash
taskset -c 0 ./build/benchmark-all-backends/rmq_benchmarks \
  --benchmark_filter='^rmq_build_cartesian_rmm/(1024|16384|262144|4194304|16777216|67108864|268435456|1073741824)(/|$)' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/cartesian_rmm_memory.json
```

Parse the JSON counters directly and round to the table precision already used
in `include/pixie/rmq.h`; do not infer memory from `sizeof(*this)` alone,
because vector capacities and nested owned structures matter.
