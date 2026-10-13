# Reproducing rank/select benchmarks

Run from the repository root (CMake 3.25+). This measures all available backends
on 12.5%, 50%, and 87.5% random inputs, using the registered size grid.

```bash
cmake --preset benchmark-all-backends
cmake --build --preset benchmark-all-backends --target rank_select_benchmarks -j 2
export PIXIE_RANK_SELECT_FIXED_SEED=1
bench=build/benchmark-all-backends/rank_select_benchmarks

taskset -c 0 "$bench" \
  --benchmark_filter='^rank_select.*_(rank1|rank0|select1|select0)_(12p5|50|87p5)/' \
  --benchmark_repetitions=3 \
  --benchmark_out=build/rank-select-queries.json --benchmark_out_format=json

taskset -c 0 "$bench" \
  --benchmark_filter='^rank_select.*_build_both_(12p5|50|87p5)/' \
  --benchmark_repetitions=10 \
  --benchmark_out=build/rank-select-build.json --benchmark_out_format=json

# Requires Matplotlib; plots median wall times and auxiliary space percentages.
python3 scripts/plot_rank_select.py \
  build/rank-select-queries.json build/rank-select-build.json -o build/rank-select-plots
```

Build timings include owning adapters' source copies; auxiliary space excludes
those copies. FastBitVec stops at 2^30; 3-star and Orzo require BMI/BMI2.
