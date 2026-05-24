# CMake Examples

This file contains Pixie-specific CMake command examples. Keep the reusable
CMake workflow in `SKILL.md`.

## Pixie Build Options

AddressSanitizer without benchmarks:

```bash
cmake -B build/asan_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Debug -DPIXIE_BENCHMARKS=OFF -DENABLE_ADDRESS_SANITIZER=ON
cmake --build build/asan_${BUILD_SUFFIX} -j
```

Coverage without benchmarks:

```bash
cmake -B build/coverage_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Debug -DPIXIE_BENCHMARKS=OFF -DPIXIE_COVERAGE=ON
cmake --build build/coverage_${BUILD_SUFFIX} -j
```

Benchmarks:

```bash
cmake -B build/benchmarks_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=ON
cmake --build build/benchmarks_${BUILD_SUFFIX} -j
```

Disable AVX-512 and use AVX2 fallback:

```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DDISABLE_AVX512=ON
cmake --build build/release_${BUILD_SUFFIX} -j
```

Custom march flag:

```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DMARCH=icelake-client
cmake --build build/release_${BUILD_SUFFIX} -j
```

Tests only, without benchmarks or third-party benchmark dependencies:

```bash
cmake -B build/release_${BUILD_SUFFIX} -DCMAKE_BUILD_TYPE=Release -DPIXIE_BENCHMARKS=OFF
cmake --build build/release_${BUILD_SUFFIX} -j
```
