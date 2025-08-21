# pixie

`pixie` is a **succinct data structures library**.
---

## Features

* **BitVector**
  * Data structure with 3.61% overhead supporting rank and select for 1 bits. Select support for 0 bits require additional 0.39%, currently not implemented
  * Supports:
    * `rank(i)`: number of set bits (`1`s) up to position `i`.
    * `select(k)`: position of the `k`-th set bit.
  * Implementation mainly follows [1] with SIMD optimizations similar to [2] 
  * AVX-512 support is mandatory for now and thus will not compile without it.
---

## Requirements

* C++20
* Compiler with AVX-512 support recommended for best performance.
* [CMake](https://cmake.org/) ≥ 3.15.

---

## Build Instructions

```bash
git clone https://github.com/Malkovsky/pixie.git
cd pixie
mkdir build && cd build
cmake ..
make -j
```

This will build the library along with benchmarks and tests.

---

## Running Tests

After building:

```bash
/unittests
```

---

## Running Benchmarks

Bencharks are random 50/50 0-1 bitvectors up to $2^34$ bits.

```bash
./benchmarks
```

---

## Example Usage

```cpp
#include "bitvector.h"
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

---

## References

[1] Laws et al., *SPIDER: Improved Succinct Rank and Select Performance* [SPIDER](https://github.com/williams-cs/spider)
[2] Kurpicz, *Engineering compact data structures for rank and select queries on bit vectors* [pasta-toolbox/bit\_vector](https://github.com/pasta-toolbox/bit_vector)

---
