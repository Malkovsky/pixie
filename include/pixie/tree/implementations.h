#pragma once

/**
 * @file implementations.h
 * @brief All rooted ordered-tree encodings provided by Pixie.
 *
 * - `LoudsTree`: level-order unary degree sequence.
 * - `BPTree<RMMTree>`: balanced-parentheses encoding.
 * - `DFUDSTree<RMMTree>`: depth-first unary degree sequence.
 */

#include <pixie/tree.h>
#include <pixie/tree/bp.h>
#include <pixie/tree/dfuds.h>
#include <pixie/tree/louds.h>
