# Pixie

<img src="https://raw.githubusercontent.com/Malkovsky/pixie/master/src/docs/images/logo.png" alt="Pixie logo" width="256" align="left" style="float: left; margin-right: 16px; margin-bottom: 8px;" />

[![Build & Test](https://github.com/Malkovsky/pixie/actions/workflows/build-test.yml/badge.svg?branch=main)](https://github.com/Malkovsky/pixie/actions/workflows/build-test.yml)
[![codecov](https://codecov.io/github/Malkovsky/pixie/graph/badge.svg?token=413VBX7M2U)](https://codecov.io/github/Malkovsky/pixie)
[![Documentation](https://img.shields.io/badge/docs-doxygen-green.svg)](https://malkovsky.github.io/pixie/)

`pixie` is a **succinct data structures library**.

<br clear="left" />

---

## Features

* **BitVector**
  * Data structure with 3.61% overhead supporting rank and select for 1 bits.
  * Supports:
    * `rank(i)`: number of set bits (`1`s) up to position `i`.
    * `select(k)`: position of the `k`-th set bit.
    * Similar operations `rank0/select0` for `0`.
  * Implementation mainly follows [1] with SIMD optimizations similar to [2]
  * Optimized via AVX-512/AVX-2, for large binary sequences performance is I/O bounded.
* **RmMTree**
  * Implementation of a range min-max tree, it supports `rank`, `select` and `excess`-related operations allowing for a fast navigation in DFUDS/BP trees.
  
---

## Requirements

* C++20
* [CMake](https://cmake.org/) ≥ 3.18.

---

## Build Instructions

```sh
git clone https://github.com/Malkovsky/pixie.git
cd pixie
cmake --preset release
cmake --build --preset release
```

Manual alternative:

```sh
mkdir -p build/release
cmake -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j
```

Tests are enabled by default (`PIXIE_TESTS=ON`). Benchmarks are opt-in; enable with `-DPIXIE_BENCHMARKS=ON` or configure with the `benchmarks-all` preset, you can use `benchmark-diagnostic` preset for performance diagnostics (Release with debug info + performance counters support). 

---

## Running Tests

After building with presets, binaries are located in `build/release`.

### BitVector

```sh
./build/release/unittests
```

### RmM Tree

```sh
./build/release/test_rmm
```

---

## Coverage

Configure a coverage build with GCC (benchmarks disabled):

```sh
cmake --preset coverage
cmake --build --preset coverage
```

Run tests and generate the gcov text report:

```sh
./scripts/coverage_report.sh
```

---

## Running Benchmarks

Before running benchmarks, configure with presets:

```sh
cmake --preset benchmarks-all
cmake --build --preset release
```

For a RelWithDebInfo diagnostic build, use:

```sh
cmake --preset benchmarks-diagnostic
cmake --build --preset release
```

### BitVector

Benchmarks are random 50/50 0-1 bitvectors up to $2^{34}$ bits.

```sh
./build/release/benchmarks
```

### RmM Tree

```sh
./build/release/bench_rmm
```

For comparison with range min-max tree implementation from [sdsl-lite](https://github.com/simongog/sdsl-lite) (Release build required; use the release preset or `-DCMAKE_BUILD_TYPE=Release`):

```bash
sudo cpupower frequency-set --governor performance
./build/release/bench_rmm_sdsl --benchmark_out=rmm_bench_sdsl.json
```

For visualization, write the JSON output to a file using `--benchmark_out=<file>` (e.g. `./build/release/bench_rmm --benchmark_out=rmm_bench.json`) and plot it with `scripts/plot_rmm.py` (add `--sdsl-json rmm_bench_sdsl.json` for comparison).

---

## Example Usage

```cpp
#include <pixie/bitvector.h>
#include <vector>
#include <iostream>

using namespace pixie;

int main() {
    std::vector<uint64_t> bits = {0b101101}; // 6 bits
    BitVector bv(bits, 6);

    std::cout << "bv: " << bv.to_string() << "\n";     // "101101"
    std::cout << "rank(4): " << bv.rank(4) << "\n";    // number of ones in first 4 bits
    std::cout << "select(2): " << bv.select(2) << "\n"; // position of 2nd one-bit
}
```

```cpp
#include <pixie/rmm_tree.h>
#include <string>
#include <iostream>

using namespace pixie;

int main() {
    // root
    // ├─ A
    // │  ├─ a1
    // │  └─ a2
    // ├─ B
    // └─ C
    //    └─ c1
    std::string bits = "11101001011000";
    RmMTree t(bits);

    std::cout << "close(1): " << t.close(1) << "\n";     // expected 6 (A)
    std::cout << "open(3): " << t.open(3) << "\n";       // expected 2 (a1)
    std::cout << "enclose(1): " << t.enclose(1) << "\n"; // expected 0 (root)
}
```

---

## References

  - [1] Laws et al., *SPIDER: Improved Succinct Rank and Select Performance* [SPIDER](https://github.com/williams-cs/spider)

  - [2] Kurpicz, *Engineering compact data structures for rank and select queries on bit vectors* [pasta-toolbox/bit\_vector](https://github.com/pasta-toolbox/bit_vector)

---
