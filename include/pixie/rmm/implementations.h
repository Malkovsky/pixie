#pragma once

/**
 * @file implementations.h
 * @brief All RmM implementations provided or adapted by Pixie.
 *
 * - `RmMBTree`: cache-oriented B-tree implementation.
 * - `SdslRmMTree`: optional SDSL adapter when `SDSL_SUPPORT` is enabled.
 */

#include <pixie/rmm.h>
#include <pixie/rmm/btree.h>

#ifdef SDSL_SUPPORT
#include <pixie/rmm/sdsl.h>
#endif
