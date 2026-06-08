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

Tests are enabled by default (`PIXIE_TESTS=ON`). Benchmarks are opt-in; enable with `-DPIXIE_BENCHMARKS=ON` or configure with the `benchmarks-all` preset. Use `benchmarks-third-party` for comparison backends such as sdsl-lite, and `benchmarks-diagnostic` for performance diagnostics (Release with debug info + performance counters support).

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

Write JSON and plot size-scaled benchmark curves with a log-scaled x-axis:

```sh
./build/release/benchmarks --benchmark_out=bitvector_bench.json --benchmark_out_format=json
python3 scripts/plot_size_benchmarks.py bitvector_bench.json -o graphs/bitvector_size.png --size-key n
```

### Excess Positions

```sh
./build/release/excess_positions_benchmarks --benchmark_out=excess_positions.json --benchmark_out_format=json
python3 scripts/excess_benchmark_table.py excess_positions.json -o src/docs/excess_positions_benchmark_results.md
```

Generated benchmark documentation can be written to `src/docs/benchmark_results.md`;
the documentation pipeline does not run benchmarks.

### Adding an RMQ Benchmark

Value RMQ implementations are benchmarked through the common CRTP interface in
`pixie::rmq::RmqBase`. To add a comparable backend, implement a non-owning index
that can be constructed from `std::span<const T>` and provides:

* `size_impl()`
* `arg_min_impl(left, right)` for half-open ranges `[left, right)`
* `value_at_impl(position)`

The public `size()`, `empty()`, `arg_min()`, and `range_min()` methods are then
provided by `RmqBase`. Ties should return the smaller original position.

Minimal example:

```cpp
#include <pixie/rmq/rmq_base.h>

#include <cstddef>
#include <functional>
#include <span>

namespace pixie::rmq {

template <class T, class Compare = std::less<T>>
class LinearRmq : public RmqBase<LinearRmq<T, Compare>, T> {
 public:
  using Self = LinearRmq<T, Compare>;
  static constexpr std::size_t npos = RmqBase<Self, T>::npos;

  explicit LinearRmq(std::span<const T> values, Compare compare = Compare())
      : values_(values), compare_(compare) {}

  std::size_t size_impl() const { return values_.size(); }

  std::size_t arg_min_impl(std::size_t left, std::size_t right) const {
    if (left >= right || right > values_.size()) {
      return npos;
    }
    std::size_t best = left;
    for (std::size_t i = left + 1; i < right; ++i) {
      if (compare_(values_[i], values_[best])) {
        best = i;
      }
    }
    return best;
  }

  T value_at_impl(std::size_t position) const { return values_[position]; }

 private:
  std::span<const T> values_;
  Compare compare_;
};

}  // namespace pixie::rmq
```

Then add rows to `src/benchmarks/bench_rmq.cpp` in `register_benchmarks()`.
Use `run_value_rmq_build` for construction cost and `run_queries` for query
cost:

```cpp
benchmark::RegisterBenchmark(
    "rmq_build_linear",
    &run_value_rmq_build<pixie::rmq::LinearRmq<
        std::int64_t, std::less<std::int64_t>>>)
    ->Arg(static_cast<std::int64_t>(size))
    ->Unit(benchmark::kMillisecond);

benchmark::RegisterBenchmark(
    "rmq_linear",
    &run_queries<pixie::rmq::LinearRmq<
        std::int64_t, std::less<std::int64_t>>>)
    ->Args({static_cast<std::int64_t>(size),
            static_cast<std::int64_t>(width)})
    ->Unit(benchmark::kNanosecond);
```

The RMQ benchmark harness rotates through several value arrays so results are
less dependent on the global-minimum position. The `width` argument is the
maximum query width, not an exact width.

To compare the new backend with `SegmentBTreeXL`, run both benchmark families
with a Google Benchmark filter. For example, after registering the new backend
as `rmq_linear` and `rmq_build_linear`:

```sh
./build/release/bench_rmq \
  --benchmark_filter='^(rmq_linear|rmq_segment_btree_xl)/(4194304|16777216)/(64|4096|262144|4194304|16777216)$'

./build/release/bench_rmq \
  --benchmark_filter='^(rmq_build_linear|rmq_build_segment_btree_xl)/(262144|4194304|16777216)$'
```

The first command compares query time for `2^22` and `2^24` input sizes across
the common RMQ widths. The second command compares construction time for the
same implementations.

### RmM Tree

```sh
./build/release/bench_rmm
```

For focused runs, `bench_rmm` accepts `--ops` with a comma-separated operation list. The benchmark harness only builds the query pools needed by the selected operations, so subset runs avoid much of the setup cost:

```sh
./build/release/bench_rmm --ops=rank1,select1 --benchmark_out=rmm_rank_select.json --benchmark_out_format=json
```

By default, RmM benchmarks step through sizes by powers of two. Use `--per_octave=<n>` for finer sampling between adjacent powers of two, or `--explicit_sizes=<csv>` for an exact size list.

Google Benchmark filters are also used to limit RmM setup when `--ops` is not provided:

```sh
./build/release/bench_rmm --benchmark_filter='^rank1$' --benchmark_out=rmm_rank1.json --benchmark_out_format=json
```

For comparison with range min-max tree implementation from [sdsl-lite](https://github.com/simongog/sdsl-lite), use the third-party benchmark preset. This defines `SDSL_SUPPORT` and builds `bench_rmm_sdsl`:

```bash
cmake --preset benchmarks-third-party
cmake --build --preset benchmarks-third-party
sudo cpupower frequency-set --governor performance
./build/release-third-party/bench_rmm_sdsl --benchmark_out=rmm_bench_sdsl.json
```

For visualization, write the JSON output to a file using `--benchmark_out=<file>` (e.g. `./build/release/bench_rmm --benchmark_out=rmm_bench.json`) and plot it with `scripts/plot_rmm.py` (add `--sdsl-json rmm_bench_sdsl.json` for per-operation sdsl-lite comparison plots). For size-scaled tree plots, use:

```sh
python3 scripts/plot_size_benchmarks.py rmm_bench.json -o graphs/rmm_size.png --size-key N
```

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
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

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
    std::vector<std::uint64_t> words((bits.size() + 63) / 64);
    for (std::size_t i = 0; i < bits.size(); ++i) {
        if (bits[i] == '1') {
            words[i / 64] |= std::uint64_t{1} << (i % 64);
        }
    }

    // RmMTree is non-owning: keep words alive and immutable while using t.
    RmMTree t(words, bits.size());

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
