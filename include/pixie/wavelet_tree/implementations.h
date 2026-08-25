#pragma once

/**
 * @file implementations.h
 * @brief All wavelet-tree implementations provided by Pixie.
 *
 * - `WaveletTreeIndex<Symbol, Storage>`: typed, storage-parameterized tree.
 * - `WaveletTree<Symbol>`: owning aligned-storage alias.
 * - `WaveletTreeView<Symbol>`: non-owning read-only storage view alias.
 */

#include <pixie/wavelet_tree.h>
#include <pixie/wavelet_tree/index.h>
