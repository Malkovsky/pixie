# Pixie focused validation

Read the shared `focused-validation` skill first. During an experiment, scope
the current header/operation, not all changes accumulated on this feature branch.
The normal `release` preset provides tests; `benchmarks` and
`benchmark-all-backends` provide benchmark targets. Analyze the appropriate
database for each. Reuse a compatible configured tree instead of configuring
and building a fresh integration tree for every edit.

## Native metadata and full-suite comparison

Checked 2026-09-28 on `build/immer-review-release`, GCC 13.3, Release,
Unix Makefiles, CMake/CTest 3.28.3, optional adapters enabled, `-j 2`.
Request `codemodel-v2` and `cmakeFiles-v1` in the existing tree using the shared
skill's File API setup, then refresh the CTest inventory after building.

```sh
ctest --test-dir build/immer-review-release -C Release \
  --show-only=json-v1 > /tmp/pixie-ctest-inventory.json
python3 -B agentic/cpp/skills/focused-validation/analyze_scope.py \
  --compile-commands build/immer-review-release/compile_commands.json \
  --changed-file include/pixie/experimental/slot_order16.h \
  --kind tests --ctest-json /tmp/pixie-ctest-inventory.json --format summary
```

No target names or runtime filters are needed to discover the affected family
slice. Use `--format json` redirected outside the repo for an executable plan.
The three automatically selected scopes were:

- `slot_order16.h`: `bit_sequence_tests`, `permutable_sequence_tests`.
- `immer_sequence.h`: `permutable_sequence_tests`.
- `sequence_tree.h`: `bit_sequence_tests`, `permutable_sequence_tests`,
  `permutation_tests`, `sorted_sequence_tests`.

GNU Make dry runs, treating each header as newly changed (`-n -W` against each
test target's generated `build.make`), independently selected exactly the same
translation-unit rebuild targets across all 19 configured test targets.

| Validation scope                 | Test targets | Test cases | Warm build (s) | CTest (s) |
| -------------------------------- | -----------: | ---------: | -------------: | --------: |
| Whole configured repository      |           19 |        876 |          14.58 |     41.40 |
| Slot-order affected families     |            2 |        449 |           1.83 |      6.84 |
| Immer-adapter affected family    |            1 |        228 |           1.02 |      2.83 |
| Shared sequence-tree families    |            4 |        507 |           3.28 |      6.89 |
| SlotOrder16 explicit case filter |            1 |          2 |           0.97 |      0.42 |

All runs passed; executed counts were checked against the inventory. These are
single-run wall times including process startup, not statistical benchmarks.
The table's warm builds performed **zero C++ compilations**. The separate first
whole-test build took 243.34 s and compiled 15 missing/unbuilt test translation
units; the four sequence/permutation targets were already built. That is not a
clean-build speedup comparison. A healthy incremental whole build would already
recompile only the affected sources; focused selection avoids unrelated target
traversal, missing-binary builds, and test execution. Runtime filtering still
does not reduce compilation within a selected family translation unit.

Final unrestricted test-only discovery took 7.61 s: 21 depfile hits, zero
compiler scans, zero AST invocations, and authoritative CMake ownership. This
cost is separate from the table. `--kind tests` avoids scanning unbuilt benchmark
translation units; use `all` or `benchmarks` when those are needed. During a
known experiment, an explicit target restriction can narrow discovery further,
but excludes other targets from the impact claim.

For large typed suites, use the generated verified-label selection rather than
an enormous name regex. CTest 3.28 has a regex-size limit and can otherwise
report no tests without failing. The analyzer now adds `--no-tests=error`,
uses exact label coverage when available, and falls back to inventory indices
for large partial selections. The numeric fallback was also checked against
this tree's CTest inventory. Regenerate index-based plans after builds or
registration changes. Small mounted-filesystem clock-skew warnings occurred
during some build checks; no repeated compilation occurred in the warm runs.

## Checked example: 16-slot child order

This example was exercised on the existing combined Release tree
`build/immer-review-release` (tests, benchmarks, and optional adapters enabled).
Substitute the actual configured tree; do not require that ad hoc directory in
another checkout.

```sh
ctest --test-dir build/immer-review-release -C Release \
  --show-only=json-v1 > /tmp/pixie-ctest-inventory.json
python3 agentic/cpp/skills/focused-validation/analyze_scope.py \
  --compile-commands build/immer-review-release/compile_commands.json \
  --changed-file include/pixie/experimental/slot_order16.h \
  --target bit_sequence_tests --target permutable_sequence_benchmarks \
  --ctest-json /tmp/pixie-ctest-inventory.json \
  --test-regex 'SlotOrder16[.]' --benchmark-regex '^SlotOrder16/' --format json
```

Review and execute the emitted build/CTest commands. The test inventory contains
target-prefixed names such as
`bit_sequence_tests.SlotOrder16.EveryRotationRangeAndDistance`; raw Google Test
names alone are not exact CTest names. This slice selects two primitive tests.

The slot-order benchmarks are in **permutable_sequence_benchmarks**, not
bit_sequence_benchmarks. Inventory the exact filter before timing:

```sh
build/immer-review-release/permutable_sequence_benchmarks \
  --benchmark_list_tests=true --benchmark_filter='^SlotOrder16/'
python3 scripts/benchmark_sequences.py \
  --build-dir build/immer-review-release \
  --output-dir /tmp/pixie-slot-order-screen \
  --binary permutable_sequence_benchmarks --filter '^SlotOrder16/' \
  --cpu 4 --repetitions 3 --min-time 0.01s --timeout 30
```

This selects four rows: packed order and pointer-array controls, each with
rotate/read and dependent-read probes. These timings screen the primitive;
they do not validate general tree mutations or compare the short 16-slot
dependency cycle with the full-period 32/64 cycles.

## Avoid rebuilding a large family translation unit in the edit loop

`permutable_sequence_benchmarks.cpp` instantiates many layouts and payload types.
A four-row runtime filter does not make its recompilation small. For the first
iterations, create a scratch CMake project under `/tmp` with:

- A correctness executable including only
  `<pixie/experimental/slot_order16.h>` and the standard library. Compare every
  range/distance for counts 0..16 with `std::rotate`, including a non-identity
  starting order. Use checks that remain active with `NDEBUG`, and register the
  executable with CTest.
- A small Google Benchmark executable including that concrete header, with
  packed-order and direct-pointer controls performing equivalent work. Reuse
  the already-built Release Google Benchmark library/includes from the selected
  tree; do not fetch/build another copy. Keep scratch results outside the repo.
- The real target's compiler and relevant options, not an arbitrary default.
  The checked setup used `/usr/bin/c++`, C++20, `-O3 -DNDEBUG`, and
  `-march=native`. Read the current compile database/link command when recreating
  it. Use the sanitizer/fallback configuration separately when that is the
  behavior under test.

Configure once, then rebuild just the scratch target after an edit. Never put
the full implementation catalog or unrelated payload types into the probe.
Once a candidate survives the screen, rebuild the real family targets and run
their shared specifications and relevant registered integrated benchmarks.

The 2026-09-27 workflow check used a fresh scratch build: the correctness probe
compiled in 1.59 s and passed 25,194 differential rotation cases through CTest;
the two-case benchmark probe compiled in 2.10 s and ran in 0.10 s. Its dependency
library was already built. No speedup claim about the data structure comes from
these workflow checks.

The four registered benchmark rows took 0.44 s in the primitive workflow check.
Keep benchmark screening separate from compilation and test costs.

## Broaden deliberately

Without the explicit target restriction, the checked dependency graph identified:

- `slot_order16.h`: bit-sequence tests, permutable-sequence tests, and
  permutable-sequence benchmarks. The two-target experiment above intentionally
  defers the integrated permutable-sequence tests.
- `immer_sequence.h`: permutable-sequence tests and benchmarks only.
- `sequence_tree.h`: all four sequence/permutation test targets and all four
  corresponding benchmark targets, including sorted wrappers.

Regenerate the graph when includes or CMake configuration change. Full
test-and-benchmark discovery scans more translation units than `--kind tests`;
that broader discovery is not required on every edit once the family is established.

Use family labels for selected-candidate validation; use the full registered
matrix only for explicit broad requests/finalization. Its largest fixtures and
setup/reset work dominate wall time even with a small `--benchmark_min_time`.
Retain the RAM guards, provide actual L3 size for beyond-L3 probes, and report
skips. Separate builds/tests from benchmark timing. Investigate clock-skew
warnings on the mounted workspace before attributing no-op rebuilds to C++ work.
