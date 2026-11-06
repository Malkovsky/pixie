#pragma once

/**
 * @file implementations.h
 * @brief Catalog of Pixie storage implementations.
 *
 * - `AlignedStorage`: owning, fixed-size, 64-byte-aligned storage.
 * - `AlignedStorageBuilder`: growable aligned construction storage.
 * - `ReadOnlyStorageView`: non-owning read-only byte storage.
 * - `SlidingWindowStorage`: owning fixed-capacity sliding byte storage.
 */

#include <pixie/storage/aligned.h>
#include <pixie/storage/aligned_builder.h>
#include <pixie/storage/read_only_view.h>
#include <pixie/storage/sliding_window.h>
