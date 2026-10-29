# Archived 128-node BP SIMD prototype

Extracted from `profiling/microtree-level-simd-2026-10-26` before cleanup.
This later prototype includes the earlier 256-bit matching experiments plus
level masks and AVX-512 VBMI2 byte compression. It is an experimental source
attachment to the notes, not part of the tree API or the selected traversal
adapter. The algorithms and historical results are described in
[Decoding and SIMD](../decoding-and-simd.md).

- `micro256.h`: scalar, SIMD matching, depth-mask and level-wise decoders;
  AVX2/scalar fallbacks and guarded VBMI2 variants.
- `probe.cpp`: exact event/child-handle comparisons over all 2,056 rooted
  ordered shapes through nine nodes, random trees, chains, stars, boundary
  sizes and dirty tail bits. Its optional timing mode compares full traversals.
- `kernels.cpp`: optional isolated decode/replay benchmarks with a rotating
  input pool and compiled microprogram controls.

Only include paths, formatting and standalone build plumbing were adjusted
during extraction. The original absolute dependency paths were replaced with
an optional installed Google Benchmark package. No historical binary or raw
timing output is bundled. The controls use the current Pixie headers, so a
new timing run is not an exact reconstruction of the historical comparison.

From the repository root, run correctness checks without external dependencies:

```sh
cmake -S include/pixie/tree/experimental/agentic_notes/microtree-simd \
  -B build/tree-microtree-notes -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=-march=native
cmake --build build/tree-microtree-notes --target probe -j2
ctest --test-dir build/tree-microtree-notes --output-on-failure --no-tests=error
```

For the optional timing probes, configure with `-DMICRO_BUILD_BENCHMARKS=ON`
and `-Dbenchmark_DIR=<directory-containing-benchmarkConfig.cmake>`, then build
both targets. `probe` without arguments checks correctness; with arguments it
accepts a node count, repetition-block count and optional mode. `kernels` uses
Google Benchmark flags. Inspect the source's mode mapping before comparing
results; not every mode runs by default. Native builds exercise only host
supported instructions; a separate baseline ISA build checks scalar fallbacks.

The generic SIMD decoders lost to compiled microprograms on large random
inspected DFS. A wide/shallow specialization remains a hypothesis. This source
is retained to avoid having to reconstruct the kernels, not as a recommendation
to enable them in the library.
