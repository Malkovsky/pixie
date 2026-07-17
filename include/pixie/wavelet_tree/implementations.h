#pragma once

/**
 * @file implementations.h
 * @brief All wavelet-tree implementations provided by Pixie.
 *
 * - `WaveletTreeIndex<Storage>`: storage-parameterized wavelet tree.
 * - `WaveletTree`: owning aligned-storage alias.
 * - `WaveletTreeView`: non-owning read-only storage view alias.
 */

#include <pixie/wavelet_tree.h>
#include <pixie/wavelet_tree/index.h>
