---
name: cmake
description: Compile and build CMake projects, including configuring build types, options, and running test binaries.
---

# CMake Build Skill

You now have expertise in building and configuring CMake projects. Follow these workflows:

## Build Directory Convention

Use a short commit hash suffix for committed revisions:

```bash
BUILD_SUFFIX=$(git rev-parse --short HEAD)
```

If the worktree has uncommitted changes, append a descriptive suffix so generated
artifacts cannot be confused with a clean HEAD build:

```bash
BUILD_SUFFIX=$(git rev-parse --short HEAD)-dirty
```

If not a git repository, use

```bash
BUILD_SUFFIX=agent
```

Build directories follow the pattern `build/<preset_name>_<suffix>`.

## Using Presets (Preferred When Available)

> **Important**: `cmake --preset` sets cache variables and generator but its `binaryDir` cannot be
> overridden from the command line. To use a preset's settings with a custom build dir, pass the
> relevant `-D` flags explicitly together with `-B`. Use `--preset` only to discover what flags a
> preset applies.

**List available presets:**
```bash
cmake --list-presets
```

**Replicate a preset's settings with a custom suffix build dir:**

Release:
```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release
cmake --build build/release_${BUILD_SUFFIX} -j
```

Debug:
```bash
cmake -B build/debug_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug_${BUILD_SUFFIX} -j
```

AddressSanitizer (mirrors `asan` preset):
```bash
cmake -B build/asan_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Debug -DPIXIE_BENCHMARKS=OFF -DENABLE_ADDRESS_SANITIZER=ON
cmake --build build/asan_${BUILD_SUFFIX} -j
```

Coverage (mirrors `coverage` preset):
```bash
cmake -B build/coverage_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Debug -DPIXIE_BENCHMARKS=OFF -DPIXIE_COVERAGE=ON
cmake --build build/coverage_${BUILD_SUFFIX} -j
```

Benchmarks (mirrors `benchmarks` preset):
```bash
cmake -B build/benchmarks_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON
cmake --build build/benchmarks_${BUILD_SUFFIX} -j
```

## Additional Feature Options

**Disable AVX-512 (use AVX2 fallback):**
```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DDISABLE_AVX512=ON
cmake --build build/release_${BUILD_SUFFIX} -j
```

**Custom march flag:**
```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DMARCH=icelake-client
cmake --build build/release_${BUILD_SUFFIX} -j
```

**Tests only (no benchmarks or third-party deps):**
```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=OFF
cmake --build build/release_${BUILD_SUFFIX} -j
```

## Best Practices

1. **Use out-of-source builds**: Keep build artifacts in `build/<preset_name>_<suffix>` directories
2. **Presets fix binaryDir**: `--preset` cannot be combined with `-B` to change the build dir; replicate `-D` flags manually with `-B` instead
3. **Reconfigure when options change**: Rerun the `cmake -B ...` step when toggling options
4. **Clean build directory when needed**: Delete the entire build folder for a fresh configuration
5. **Match build type to task**: Release for performance work, Debug/ASan for correctness
