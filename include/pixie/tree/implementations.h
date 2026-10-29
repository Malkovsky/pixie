#pragma once

/**
 * @file implementations.h
 * @brief All rooted ordered-tree encodings provided by Pixie.
 *
 * - `LoudsTree<RankSelect>`: level-order unary degree sequence with a
 * selectable rank/select backend (defaults to `RankSelectSupport<>`).
 * - `BPTree<RMMTree>`: balanced-parentheses encoding.
 * - `DFUDSTree<RMMTree>`: depth-first unary degree sequence.
 *
 * Include the individual concrete header in consumers. Traversal workloads
 * live in the benchmark harness; encoding-specific unpacking kernels remain
 * internal to the corresponding implementations.
 */

#include <pixie/tree.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>

// Bounded experimental auxiliary index for root-to-leaf navigation.
#include <pixie/tree/experimental/second_chance_tree_cache.h>
#include <pixie/tree/experimental/tree_navigation_index.h>
#include <pixie/tree/experimental/unpacked_tree_cache.h>
