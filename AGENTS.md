# AGENTS.md - AI Coding Assistant Guidelines for Pixie

## Project Overview

Pixie is a **header-only succinct data structures library** written in C++20.
It targets practical performance for data sizes up to `2^64` bits while keeping
space usage close to the theoretical minimum.

Current library families are:

- rank/select support over packed bit sequences;
- range min-max (RmM) indexes;
- static range-minimum-query (RMQ) indexes;
- rooted-tree encodings (LOUDS, balanced parentheses, and DFUDS);
- aligned owning storage, read-only storage views, and mapped-file resources;
- storage-backed wavelet-tree indexes.

Planned additions include FM-indexes, compressed suffix arrays, and other
succinct data structures.

## Shared Agentic Guidance

`agentic/cpp` is the canonical source of reusable C++ agent guidance, skills,
and commands. Read `agentic/cpp/AGENTS.md` for shared workflow and skill policy
before relying on a shared skill or command. Do not duplicate that procedural
guidance in Pixie-specific instructions.

Pixie-specific examples and additions live in `agentic/local/cpp`:

- `agentic/local/cpp/skills/<skill>/EXAMPLES.md` supplies local examples for a
  matching shared skill.
- `agentic/local/cpp/commands/` supplies Pixie-specific command notes.

When a task matches a skill, read the shared skill first and then its local
example, if present:

1. `agentic/cpp/skills/<skill>/SKILL.md`
2. `agentic/local/cpp/skills/<skill>/EXAMPLES.md`

Use the shared skill or command as the canonical workflow. The local overlay
adds Pixie context; it does not replace the shared guidance.

## Architecture

### Project Layout Conventions

- **`include/pixie/`**: Header-only public API and implementations.
- **`include/pixie/<family>.h`**: Lightweight CRTP contract for one library
  family, such as `rank_select.h`, `rmq.h`, or `storage.h`.
- **`include/pixie/<family>/`**: Concrete implementations for that family.
- **`include/pixie/<family>/implementations.h`**: Catalog and umbrella include
  for the family's concrete implementations. It is the single family-level
  include for tests and benchmarks; a concrete implementation must not include
  its own catalog.
- **`include/pixie/experimental/`**: Isolated experimental primitives and
  implementations. Do not promote an experiment without tests and a registered
  benchmark where performance is relevant.
- **`include/pixie/io/`**: Resource-owning I/O helpers, such as `MappedFile`.
- **`include/references/`**: Naive reference implementations used for
  differential testing.
- **`src/tests/`**: Google Test suites. Family tests are consolidated into
  shared specification suites rather than one executable per implementation.
- **`src/benchmarks/`**: Google Benchmark harnesses and common benchmark
  helpers.
- **`src/docs/`**: Doxygen configuration and generated-documentation inputs.
- **`academic/`**: Durable scholarly and technical documentation: papers,
  presentations, reports, research and lecture notes, shared figures, and the
  canonical BibTeX bibliography. Follow `academic/AGENTS.md` within that tree.

### Family Interface Pattern

Public data-structure families use CRTP contracts. The current contracts are
`RankSelectBase`, `RmMBase`, `pixie::rmq::RmqBase`, `TreeBase`, `StorageBase`,
and `WaveletTreeBase`.

1. Define or extend the public contract in `include/pixie/<family>.h`.
   Public facade methods delegate to a clearly named `*_impl()` method on the
   concrete type.
2. Make each concrete implementation inherit its family contract and implement
   the required `*_impl()` methods. Do not add virtual dispatch for this API.
3. Add the concrete header to the corresponding `implementations.h` catalog.
   The catalog includes the contract and concrete headers, not the reverse.
4. Include the catalog in the family test and benchmark harness, then put every
   compatible implementation through the same typed specification suite.

The contract is the source of truth for observable semantics. Every public
facade operation and every extension-point requirement needs Doxygen
documentation that states its range convention, rank base, ties, sentinel
value, ownership/lifetime constraint, and invalid-input behavior where
applicable. There are no virtual public interfaces: document the CRTP facade
and its `*_impl()` contract instead of duplicating undocumented behavior in
each implementation.

### Key Design Decisions

1. **Header-only library**: implementations live under `include/`; there are
   no Pixie `.cpp` library sources.
2. **Non-owning indexes**: indexes use external immutable data where
   appropriate, commonly `std::span<const std::uint64_t>`. The caller keeps
   that data alive and stable for the index lifetime.
3. **SIMD conditional compilation**: use `PIXIE_AVX512_SUPPORT` and
   `PIXIE_AVX2_SUPPORT` guards with a scalar fallback.
4. **Target domain**: optimize for bit sequences and indexes up to `2^64`
   bits.
5. **Platform**: Linux/Unix is the current target platform. `MappedFile`
   requires POSIX `mmap`; there is no Windows implementation yet.
6. **Range conventions**: rank/select, RMQ, wavelet-tree, public tests, and
   benchmarks normally use zero-based half-open ranges `[left, right)`. Do not
   generalize this to every low-level primitive: RmM relative-excess range
   queries explicitly use their documented inclusive bit range. State the
   convention at every public boundary.
7. **Storage ownership**: keep mutable owners, read-only non-owning views, and
   OS-backed resource owners as separate types. Storage subranges are
   byte-oriented. A view never extends its backing owner's lifetime and can be
   invalidated by owner resize or destruction.
8. **Optional adapters**: third-party implementations stay behind their build
   option and must not become a library runtime dependency.

### Why Header-Only?

Pixie is header-only by design. Similar libraries, including sdsl-lite, have
made the same trade-off.

| Benefit | Explanation |
| --- | --- |
| SIMD flexibility | Users compile with their target `-march` flags and enable the best available code path. |
| Better inlining | The compiler sees hot rank/select and navigation code. |
| No ABI issues | The library works across compiler and standard-library versions. |
| Easy integration | Consumers include headers and do not link a Pixie binary. |
| Template-friendly | Templates need no explicit instantiations. |

Reconsider a compiled component only for substantial construction algorithms,
large static tables that should not be duplicated, global runtime state, or an
unavoidable external C library dependency.

## Technology Stack

- **Language**: C++20 (`std::span`, `std::popcount`, and `<bit>` are required).
- **Build**: CMake 3.18 or newer.
- **Testing**: Google Test 1.17.0, registered with CTest.
- **Benchmarking**: Google Benchmark 1.9.4.
- **SIMD**: AVX-512 when available, AVX2 fallback, and scalar fallback.
- **Style**: Chromium C++ style (`.clang-format`).

### Dependencies

Pixie has no runtime or default third-party dependency. CMake fetches build-time
dependencies only as needed. Direct CMake defaults and their effects are:

| Option | Default | Effect |
| --- | --- | --- |
| `PIXIE_TESTS` | `ON` standalone, `OFF` downstream | Builds tests and fetches Google Test. |
| `PIXIE_BENCHMARKS` | `OFF` | Builds native Google Benchmark targets. |
| `PIXIE_THIRD_PARTY_BACKENDS` | `OFF` | Enables optional SDSL adapters and their comparison targets. |
| `PIXIE_DIAGNOSTICS` | `OFF` | Enables diagnostic logging for profiling experiments. |
| `PIXIE_DOCS` | `OFF` | Enables the Doxygen `docs` target. |
| `PIXIE_COVERAGE` | `OFF` | Adds GCC coverage instrumentation. |

`MappedFile` uses native POSIX memory mapping on Linux/Unix. A default
FetchContent consumer therefore receives no third-party dependency, while a
standalone default build fetches Google Test. Enabling third-party backends also
fetches SDSL and pasta-toolbox dependencies, but only SDSL currently has a
registered Pixie adapter/comparison benchmark. Do not describe pasta-toolbox as
an available backend until Pixie adds and registers one.

## Build and Test Presets

Use CMake presets for supported builds. They keep build directories isolated
under `build/<preset>` and export `compile_commands.json`.

| Preset | Purpose |
| --- | --- |
| `debug` | Debug tests. |
| `release` | Optimized tests. |
| `asan` | Debug tests with AddressSanitizer and AVX-512 disabled. |
| `coverage` | Debug tests with GCC coverage instrumentation. |
| `benchmarks` | Native Pixie benchmarks. |
| `benchmark-all-backends` | Native benchmarks plus every currently registered optional backend target. |
| `benchmarks-profile` | RelWithDebInfo benchmarks with diagnostics and libpfm support. |
| `docs` | Doxygen documentation only. |

```bash
# Standard optimized test build.
cmake --preset release
cmake --build --preset release
ctest --preset release

# Native benchmark build.
cmake --preset benchmarks
cmake --build --preset benchmarks

# Comparison benchmarks, including the optional SDSL targets.
cmake --preset benchmark-all-backends
cmake --build --preset benchmark-all-backends

# Documentation.
cmake --preset docs
cmake --build --preset docs
```

For an ad hoc configuration, use a distinct build directory and pass explicit
options. Do not overwrite a supported preset build merely to test an unrelated
configuration.

```bash
cmake -S . -B build/local-avx2 -DCMAKE_BUILD_TYPE=Release \
  -DPIXIE_TESTS=ON -DDISABLE_AVX512=ON
cmake --build build/local-avx2 -j
ctest --test-dir build/local-avx2 --output-on-failure
```

## Testing

CTest is the canonical test entry point. Google Test discovery uses
`DISCOVERY_MODE PRE_TEST`, so CTest discovers and runs each test binary.

```bash
ctest --preset debug
ctest --preset release
ctest --preset asan
ctest --preset coverage

# Run only a registered test-target label.
ctest --preset release -L rank_select_tests
```

The registered test executables are `bit_algorithms_unittests`,
`rank_select_unittests`, `rank_select_tests`, `benchmark_tests`, `test_rmm`,
`tree_tests`, `wavelet_tree_tests`, `storage_tests`,
`excess_positions_tests`, `excess_record_lows_tests`, and `rmq_tests`. Run an
executable directly only when debugging a focused Google Test filter.

### Test Configuration via Environment Variables

- `EXCESS_POS_CASES`: randomized cases in `excess_positions_tests` (default
  `1000`).
- `EXCESS_POS_SEED`: random seed in `excess_positions_tests` (default `42`).
- `RECORD_LOWS_CASES`: randomized cases in `excess_record_lows_tests` (default
  `10000`).
- `RECORD_LOWS_SEED`: random seed in `excess_record_lows_tests` (default `42`).

### Testing Patterns

- **Shared specifications**: add every conforming concrete implementation to
  the family `TYPED_TEST` suite. A new implementation is not complete until it
  passes the same public contract as its peers.
- **Differential testing**: compare behavior with an appropriate naive type in
  `include/references/`.
- **Randomized testing**: use deterministic seeds and expose a seed only when
  it helps reproduce a failure.
- **Exhaustive short inputs**: enumerate small bit patterns and tree shapes
  when feasible.
- **Shape-forcing tests**: for RMQ and RmM, favor deterministic traversal-path
  cases over more undirected random input. Cover border correction, same-leaf
  paths, prefix/suffix selectors, sparse-overlay hit/miss paths, partial final
  blocks, and first-minimum ties.
- **Storage tests**: test owners and views through the same specification,
  including nested byte subranges, alignment constraints, serialization, and
  owner/view lifetime rules.

### Coverage and Codecov

Use the repository script for local coverage:

```bash
./scripts/coverage_report.sh
```

The script configures and builds the `coverage` preset, deletes stale
`.gcda`/`.gcov` files, runs `ctest --preset coverage`, and writes
`build/coverage/coverage.txt`.

- The report uses `gcov -pb`, so Codecov receives branch probabilities as well
  as line hits. A Codecov partial line executed but did not take every branch.
- Do not compare Codecov's branch-aware percentage with a line-only local
  summary.
- After recompiling instrumented code, stale `.gcda` files can cause an
  `overwriting an existing profile data with a different checksum` warning.
  Run the script before trusting the report.
- Header-only templates can appear in multiple gcov translation units. For a
  header such as `rmm/btree.h`, inspect its `.gcov` output or the Codecov
  aggregate instead of trusting the first block in `coverage.txt`.
- Prioritize public behavior paths over unreachable defensive or allocator
  failure paths.

## Code Style and Documentation

1. Run `clang-format` before committing; use the repository's Chromium style.
2. Put library code in the `pixie` namespace; RMQ contracts and implementations
   use `pixie::rmq`.
3. Use Doxygen comments for every public API. Base/contract headers must fully
   document each public facade operation and its implementation requirement.
4. Use `constexpr` for compile-time values.
5. Be aware of alignment. Prefer the 64-byte-aligned storage facilities where
   a hot data structure benefits from cache-line alignment rather than adding
   ad hoc aligned allocation code.
6. Keep public contracts lightweight. Do not include a family catalog from a
   concrete header, and do not add compatibility forwarding headers unless the
   user explicitly requests a compatibility layer.

## CI/CD Workflows

- **`build-test.yml`**: runs the AddressSanitizer fallback configuration
  through its presets.
- **`linter.yml`**: checks formatting for C/C++ sources.
- **`benchmarks-test.yml`**: runs benchmark checks.
- **`doxygen.yml`**: configures and builds the `docs` preset, then publishes
  `build/docs/docs/html`.

## Common Tasks for AI Agents

### Adding a Concrete Implementation

1. Select the existing family contract. Add an operation to that contract only
   when every implementation should expose it; document the semantics there.
2. Add the concrete header under `include/pixie/<family>/`, inherit the CRTP
   base, and implement the required `*_impl()` methods.
3. Add it to `<family>/implementations.h`; guard optional third-party adapters
   with the established build macro.
4. Add it to the existing shared family specification suite and compare it with
   a reference implementation where possible. Create a separate test target
   only for a genuinely new family.
5. Add comparable rows to the established family benchmark harness. Create a
   new benchmark target only for a genuinely new family.
6. Update `CMakeLists.txt` registration only when adding a new test or
   benchmark target, then run formatting and the relevant preset tests.

### Modifying SIMD Code

1. Keep an AVX-512 implementation, AVX2 implementation where useful, and a
   scalar fallback behind the existing feature guards.
2. Include `<immintrin.h>` only in translation units or headers that use SIMD
   intrinsics; do not make it a generic benchmark dependency.
3. Validate the fallback with the `asan` preset or an isolated
   `DISABLE_AVX512=ON` build.
4. Build the relevant benchmark preset before claiming a performance result.
   Use `benchmarks-profile` for hardware counters when supported by the host.

### Adding Tests

1. Test the public CRTP facade, not private helper implementation details.
2. Add a new implementation to the shared typed specification suite.
3. Add reference, boundary, random, and shape-forcing cases proportional to
   the changed behavior.
4. Use CTest preset runs for final validation.

## Performance Philosophy

- **Target domain**: data sizes up to `2^64` bits with corresponding queries.
- **Goal**: best practical performance across the target domain, not only
  asymptotic complexity.
- **Approach**: benchmark-driven optimization using Google Benchmark.
- **SIMD**: use vectorized operations only when benchmarked behavior justifies
  their complexity.
- **Cache efficiency**: respect 64-byte cache-line alignment for hot storage.

### Benchmark Result Retention

- Do not commit raw benchmark result files or ad hoc result tables for stable
  implementations.
- A user-requested current snapshot may live in the relevant implementation
  catalog header. Keep it rounded and reproducible, and visually align every
  Markdown table column in source, including headers, separators, and numeric
  values.
- Do not retain local JSON, failed probes, or before/after experiment history.
- Persist results for an experimental implementation only when that
  implementation is in the experimental area and has a registered benchmark.

## Academic Materials and References

`academic/` is Pixie's durable research and technical documentation area. Put
new papers, presentations, reports, research or lecture notes, figures, and
scholarly references there rather than scattering them across implementation
guidance or temporary root-level notes. Its canonical bibliography is
`academic/bibliography/references.bib`.

Academic material currently renders as independent Quarto documents. It is a
potential first-class part of Pixie's public documentation, but is not an input
to the CMake/Doxygen `docs` target today. Do not couple an academic document to
that pipeline incidentally; propose a dedicated publication and navigation
change when integrating it into the public documentation site.
