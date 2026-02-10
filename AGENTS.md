# AGENTS.md - AI Coding Assistant Guidelines for Pixie

## Project Overview

Pixie is a **succinct data structures library** written in C++20. It provides space-efficient data structures that use close to the theoretical minimum space while supporting efficient queries. The library targets practical performance for data sizes up to 2^64 bits.

Current implementations include BitVector, RmM Tree, and LOUDS Tree. Planned additions include wavelet trees, FM-index, compressed suffix arrays, and other succinct data structures.

## Architecture

### Project Layout Conventions

- **`include/`**: Header-only library API (all implementations here, no `.cpp` files)
- **`include/misc/`**: Naive reference implementations for differential testing
- **`src/*_tests.cpp`**: Unit tests (Google Test)
- **`src/*_benchmarks.cpp`**: Performance benchmarks (Google Benchmark)
- **`src/docs/`**: Doxygen configuration

> **Note**: Future versions may reorganize headers under `include/pixie/` to support `#include <pixie/bitvector.h>` style imports.

### Key Design Decisions

1. **Header-only library** (see rationale below)
2. **Non-owning spans**: Data structures use `std::span<const uint64_t>` for external data where appropriate.
3. **SIMD conditional compilation**: Uses `#ifdef PIXIE_AVX512_SUPPORT` / `PIXIE_AVX2_SUPPORT` with scalar fallbacks.
4. **Target domain**: Optimized for data sizes up to 2^64 bits.
5. **Platform**: Linux/Unix is the primary target platform.

### Why Header-Only?

The library is header-only by design. This decision is based on analysis of similar libraries (notably sdsl-lite, which migrated from compiled to header-only in v3).

**Advantages for Pixie:**

| Benefit | Explanation |
|---------|-------------|
| **SIMD flexibility** | Users compile with their target `-march` flags, enabling optimal SIMD code paths |
| **Better inlining** | Compiler sees full implementation, enabling aggressive optimization for small hot functions (rank, select) |
| **No ABI issues** | Works across different compilers and standard library versions |
| **Easy integration** | Users just `#include` headers; no library linking or installation required |
| **Template-friendly** | Templates work naturally without explicit instantiation |

**When to reconsider:**

A compiled component may be warranted if Pixie adds:
- Heavy construction algorithms (e.g., suffix array construction)
- Large static lookup tables that shouldn't be duplicated
- Global runtime state (memory tracking, configuration)
- External C library dependencies

## Technology Stack

- **Language**: C++20 (required features: `std::span`, `std::popcount`, `<bit>`)
- **Build**: CMake >= 3.18
- **Testing**: Google Test v1.17.0
- **Benchmarking**: Google Benchmark v1.9.4
- **SIMD**: AVX-512 (primary), AVX2 (fallback), scalar fallbacks
- **Style**: Chromium C++ style (`.clang-format`)

### Dependencies

Pixie itself is header-only and has **no runtime dependencies**. Build-time dependencies are managed via CMake FetchContent and controlled by two options:

| Option | Default | What it enables |
|--------|---------|-----------------|
| `PIXIE_TESTS` | `ON` | Unit tests (fetches Google Test) |
| `PIXIE_BENCHMARKS` | `ON` | Benchmarks + comparison benchmarks (fetches Google Benchmark, pasta-toolbox, sdsl-lite v3) |

Third-party libraries (pasta-toolbox, sdsl-lite) are used **only** for comparison benchmarks, not by the library itself.

## Build Commands

```bash
# Standard build (Release)
cmake -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j

# Debug build
cmake -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug -j

# Without AVX-512 (AVX2 fallback)
cmake -B build/release -DDISABLE_AVX512=ON
cmake --build build/release -j

# With AddressSanitizer
cmake -B build/asan -DENABLE_ADDRESS_SANITIZER=ON
cmake --build build/asan -j

# Custom march flag
cmake -B build/release -DMARCH=icelake-client
cmake --build build/release -j

# Tests only (no benchmarks or third-party deps)
cmake -B build/release -DPIXIE_BENCHMARKS=OFF
cmake --build build/release -j

```

## Testing

### Running Tests

```bash
cmake-build/release/unittests           # BitVector tests
cmake-build/release/test_rmm            # RmM Tree tests
cmake-build/release/louds_tree_tests    # LOUDS Tree tests
```

### Test Configuration via Environment Variables

For `test_rmm.cpp`:
- `RMM_CASES`: Number of test cases
- `RMM_OPS`: Number of operations per case
- `RMM_MAX_N`: Maximum input size
- `RMM_SEED`: Random seed for reproducibility

### Testing Patterns

- **Differential testing**: Compare against naive reference implementations in `include/misc/`
- **Randomized testing**: Random bit vectors and balanced parentheses sequences
- **Exhaustive short inputs**: Test all patterns for small sizes

## Code Style Guidelines

1. **Formatting**: Run `clang-format` before committing (Chromium style)
2. **Namespace**: All library code in `pixie` namespace
3. **Documentation**: Use Doxygen-style comments for public API
4. **Constants**: Use `constexpr` for compile-time values
5. **Alignment**: Be ware of data alignment, in most cases it is preferable to use 64-byte aligned array allocations

## CI/CD Workflows

- **build-test.yml**: Builds with/without AVX-512, runs tests with ASan
- **linter.yml**: Clang-format checks on all C/C++ files
- **benchmarks-test.yml**: Performance benchmark runs
- **doxygen.yml**: Documentation generation

## Common Tasks for AI Agents

### Adding a New Data Structure

1. Create header in `pixie/include/` with Doxygen documentation
2. Add unit tests in `src/tests/<name>_tests.cpp`
3. Add benchmarks in `src/benchmarks/<name>_benchmarks.cpp`
4. Update `CMakeLists.txt` with new executables
5. Run `clang-format` on new files

### Modifying SIMD Code

1. Provide implementations for:
   - AVX-512 (`#ifdef PIXIE_AVX512_SUPPORT`)
   - AVX2 (`#ifdef PIXIE_AVX2_SUPPORT`)
   - Scalar fallback
2. Test with `-DDISABLE_AVX512=ON` to verify fallback works
3. Benchmark to ensure performance is maintained

### Adding Tests

1. Use Google Test framework
2. Include naive reference implementation for differential testing
3. Add edge cases: empty input, single element, boundary conditions
4. Use random testing with configurable seed for reproducibility

## Performance Philosophy

- **Target domain**: Data sizes up to 2^64 bits with corresponding queries
- **Goal**: Best practical performance across the target domain (not just asymptotic complexity)
- **Approach**: Benchmark-driven optimization using Google Benchmark
- **SIMD**: Leverage vectorized operations where beneficial
- **Cache efficiency**: Align data structures to cache line boundaries (64 bytes)

## References

- [SPIDER Paper](https://github.com/williams-cs/spider) - Reference for BitVector implementation
- [pasta-toolbox](https://github.com/pasta-toolbox/bit_vector) - Reference implementation for comparison
