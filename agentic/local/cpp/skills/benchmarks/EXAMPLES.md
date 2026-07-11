# Pixie Benchmark Examples

Use these notes with `agentic/cpp/skills/benchmarks/SKILL.md` for Pixie's local
benchmark binaries and reporting conventions.

## RMQ Benchmark Tables

The primary RMQ benchmark binary is usually:

```bash
./build/release-third-party/bench_rmq
```

When the table in `include/pixie/rmq.h` needs updating, use pinned Release
runs and keep the header as the current snapshot only. Do not store local JSON
paths, failed probes, or "before/after" experiment history in that header.

For query rows:

```bash
taskset -c 0 ./build/release-third-party/bench_rmq \
  --benchmark_filter='^rmq_(hybrid_btree|cartesian_hybrid_btree|cartesian_btree)/(1024|16384|262144|4194304|16777216|67108864)/' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/rmq_query.json
```

For build and memory rows:

```bash
taskset -c 0 ./build/release-third-party/bench_rmq \
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
./build/release-third-party/bench_rmq --benchmark_list_tests=true | rg cartesian_rmm
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
taskset -c 0 ./build/release-third-party/bench_rmq \
  --benchmark_filter='^rmq_build_cartesian_rmm/(1024|16384|262144|4194304|16777216|67108864|268435456|1073741824)(/|$)' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=/tmp/cartesian_rmm_memory.json
```

Parse the JSON counters directly and round to the table precision already used
in `include/pixie/rmq.h`; do not infer memory from `sizeof(*this)` alone,
because vector capacities and nested owned structures matter.
