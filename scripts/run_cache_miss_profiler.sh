perf record -e cache-misses -g -- build/release-with-deb/benchmarks \
    --benchmark_filter=RankNonInterleaved/n:1717 \
    --benchmark_report_aggregates_only=true \
    --benchmark_perf_counters=INSTRUCTIONS,CYCLES,CACHE-MISSES