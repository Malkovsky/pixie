# Pixie

<img src="https://raw.githubusercontent.com/Malkovsky/pixie/main/src/docs/images/logo.png" alt="Pixie logo" width="256" align="left" style="float: left; margin-right: 16px; margin-bottom: 8px;" />

[![Build & Test](https://github.com/Malkovsky/pixie/actions/workflows/build-test.yml/badge.svg?branch=main)](https://github.com/Malkovsky/pixie/actions/workflows/build-test.yml)
[![codecov](https://codecov.io/github/Malkovsky/pixie/graph/badge.svg?token=413VBX7M2U)](https://codecov.io/github/Malkovsky/pixie)
[![Documentation](https://img.shields.io/badge/docs-doxygen-green.svg)](https://malkovsky.github.io/pixie/)

`pixie` is a **succinct data structures library**. A general concept behind this kind of structures is to be very compact yet still processible without total decompression. One might want this kind of structures in a case when storage space is limited and CPU performance is not tight as typically succinct structures are behind in performance compared to traditional data structures.

<br clear="left" />

## Outline

- **Rank/select support**: fundamental primitive to efficiently map bits and their positions in a bit vector.
- **Range min-max tree**: fundamental primitive for efficient computation of excess queries mainly for navigation in balanced parenthesis sequences.
- **Succinct trees**: static variants of $2$-bit per entry trees, i.e. LOUDS, DFUDS, BP (Based on Euler tour and Ferrada-Navarro style).
- **Wavelet tree**, i.e. static structure that supposts rank/select on arbitrary finite alphabets, supports building a Huffman archieve with fast extraction of arbitrary segment.
- Succinct **cartesian tree** and a state of the art solution to static **RMQ** (array is immutable, queries are not known in advance).

---

## Getting started

Pixie requires C++20 and CMake 3.18 or newer. It is a header-only library and
is intended to be used via CMake FetchContent and the `pixie::pixie` target:

```cmake
include(FetchContent)

FetchContent_Declare(
  pixie
  GIT_REPOSITORY https://github.com/Malkovsky/pixie.git
  GIT_TAG <release-tag-or-commit>
)
FetchContent_MakeAvailable(pixie)

target_link_libraries(my_target PRIVATE pixie::pixie)
```

Pixie exposes lightweight CRTP contracts and catalogs of concrete
implementations. For example, check the `pixie/rmq/implementations.h` to see all the available implementations for RMQ and include what you need.

```cpp
#include <array>
#include <span>

#include <pixie/rmq/implementations.h>

const std::array<int, 6> values = {7, 3, 5, 1, 4, 1};
const pixie::rmq::HybridBTree<int> rmq{std::span<const int>(values)};

const auto position = rmq.arg_min(1, 6);  // 3: first minimum in [1, 6)
const auto minimum = rmq.range_min(1, 6); // 1
```

---

## Tests and benchmarks

Pixie uses CTest for its internal suite:

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

To benchmark another implementation against Pixie, implement the corresponding
interface and register it in the family benchmark harness.

Here's an example using the optional SDSL RmM available for  comparison through the
all-backends preset:

```sh
cmake --preset benchmark-all-backends
cmake --build --preset benchmark-all-backends
./build/benchmark-all-backends/rmm_sdsl_benchmarks
```

`SdslRmMTree` adapts `sdsl::bp_support_sada<>` to the same CRTP facade as
Pixie's native RmM implementations. The following abridged excerpt shows the
connection; the full adapter is in `include/pixie/rmm/sdsl.h`.

```cpp
#ifdef SDSL_SUPPORT
class SdslRmMTree : public RmMBase<SdslRmMTree> {
 public:
  using BpSupport = sdsl::bp_support_sada<>;
  static constexpr std::size_t npos = RmMBase<SdslRmMTree>::npos;

  std::size_t size_impl() const { return size_; }

  std::size_t rank1_impl(std::size_t end_position) const {
    if (size_ == 0 || end_position == 0) {
      return 0;
    }
    return tree_.rank(std::min(end_position, size_) - 1);
  }

  std::size_t close_impl(std::size_t open_position) const {
    if (size_ == 0) {
      return 0;
    }
    const std::size_t position = tree_.find_close(open_position);
    return position < size_ ? position : npos;
  }

 private:
  std::size_t size_{};
  sdsl::bit_vector bits_;
  BpSupport tree_;
};
#endif
```

The remaining `*_impl()` methods complete the `RmMBase` contract; benchmarks
call the inherited public facade, not SDSL-specific methods.

---

## Example Usage

```cpp
#include <array>
#include <cstdint>
#include <span>

#include <pixie/wavelet_tree/implementations.h>

int main() {
  const std::array<std::uint64_t, 6> text = {2, 0, 1, 2, 1, 0};
  pixie::WaveletTree tree(3, std::span<const std::uint64_t>(text));

  const auto ones_before_five = tree.rank(1, 5);  // 2
  const auto second_two = tree.select(2, 2);      // 3
  const auto segment = tree.get_segment(1, 4);    // {0, 1, 2}
}
```

## License

Copyright 2026 Pixie contributors.

Pixie is licensed under the [Apache License 2.0](LICENSE). Optional
third-party benchmark and backend integrations retain their own licenses.
