# Benchmark Results

These results were generated on 2026-05-18 from `build/release` binaries on the local benchmark host. JSON inputs are kept under `src/docs/benchmarks`.

## Excess Positions

Command:

```sh
./build/release/excess_positions_benchmarks \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=src/docs/benchmarks/excess_positions.json
python3 scripts/excess_benchmark_table.py \
  src/docs/benchmarks/excess_positions.json \
  -o src/docs/excess_positions_benchmark_results.md
```

| Method | X=-64 | X=-8 | X=0 | X=8 | X=64 |
|---|---:|---:|---:|---:|---:|
| BranchingLUT | 15.71 ns | 26.56 ns | 26.55 ns | 26.43 ns | 15.37 ns |
| Current | 10.73 ns | 18.09 ns | 18.58 ns | 18.24 ns | 10.43 ns |
| Expand | 60.70 ns | 88.41 ns | 87.38 ns | 88.50 ns | 56.77 ns |
| Expand8 | 19.13 ns | 53.05 ns | 47.83 ns | 49.35 ns | 17.44 ns |
| ExpandAVX512 | 23.13 ns | 38.39 ns | 38.56 ns | 39.24 ns | 23.19 ns |
| LUTAVX512 | 12.33 ns | 18.34 ns | 18.06 ns | 18.21 ns | 12.75 ns |
| Scalar | 304.42 ns | 389.58 ns | 446.94 ns | 399.80 ns | 316.23 ns |

## BitVector Size Sweep

The BitVector plot uses the 50/50 fill variants for rank/select over the registered benchmark size grid. The benchmark definitions use fixed repeats, so the command-line repetition value is not used for this binary.

Command:

```sh
./build/release/benchmarks \
  --benchmark_filter='BM_(RankInterleaved|RankNonInterleaved|RankZeroNonInterleaved|SelectNonInterleaved|SelectZeroNonInterleaved)/' \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=src/docs/benchmarks/bitvector_size.json
python3 scripts/plot_size_benchmarks.py \
  src/docs/benchmarks/bitvector_size.json \
  -o src/docs/images/benchmarks/bitvector_size.png \
  --size-key n \
  --title 'BitVector benchmark time vs size'
```

![BitVector benchmark time vs size](images/benchmarks/bitvector_size.png)

## RmM Tree Size Sweep

The RmM comparison uses operations available in both Pixie and sdsl-lite over the same power-of-two tree sizes. Pixie's benchmark harness only constructs query pools needed by the selected operations.

Pixie command:

```sh
./build/release/bench_rmm \
  --ops=rank1,rank0,select1,excess,range_min_query_pos,range_min_query_val,close,open,enclose \
  --explicit_sizes=16384,32768,65536,131072,262144,524288,1048576,2097152,4194304 \
  --Q=32768 \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=src/docs/benchmarks/rmm_tree_size.json
```

sdsl-lite command:

```sh
./build/release/bench_rmm_sdsl \
  --ops=rank1,rank0,select1,excess,range_min_query_pos,range_min_query_val,close,open,enclose \
  --explicit_sizes=16384,32768,65536,131072,262144,524288,1048576,2097152,4194304 \
  --Q=32768 \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true \
  --benchmark_display_aggregates_only=true \
  --benchmark_format=json \
  --benchmark_out=src/docs/benchmarks/rmm_tree_sdsl_size.json
```

Plot command:

```sh
python3 scripts/plot_rmm.py \
  src/docs/benchmarks/rmm_tree_size.json \
  --sdsl-json src/docs/benchmarks/rmm_tree_sdsl_size.json \
  --save-dir src/docs/images/benchmarks/rmm_comparison \
  --logx
```

![rank1 comparison](images/benchmarks/rmm_comparison/rank1.png)

![rank0 comparison](images/benchmarks/rmm_comparison/rank0.png)

![select1 comparison](images/benchmarks/rmm_comparison/select1.png)

![excess comparison](images/benchmarks/rmm_comparison/excess.png)

![range_min_query_pos comparison](images/benchmarks/rmm_comparison/range_min_query_pos.png)

![range_min_query_val comparison](images/benchmarks/rmm_comparison/range_min_query_val.png)

![close comparison](images/benchmarks/rmm_comparison/close.png)

![open comparison](images/benchmarks/rmm_comparison/open.png)

![enclose comparison](images/benchmarks/rmm_comparison/enclose.png)
