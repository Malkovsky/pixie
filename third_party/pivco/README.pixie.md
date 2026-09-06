# PivCo-Huffman backend

This directory contains the Apache-2.0 PivCo-Huffman C implementation used by
Pixie's Huffman codec. Its `LICENSE` file is retained alongside the sources.
Pixie builds the plain PivCo-Huffman file codec and the applicable scalar/SIMD
backend. The optional FSE integration, standalone tools, tests, paper sources,
and historical results are not included.

Pixie-specific build logic lives in `cmake/PivCo.cmake`. The public C++ adapter
is `include/pixie/huffman/pivco_huffman.h` and is available through the Huffman
implementation catalog. Changes to the vendored source files should carry the
notices required by the Apache-2.0 license.
