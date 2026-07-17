#pragma once

/**
 * @file implementations.h
 * @brief Catalog of Pixie storage implementations.
 *
 * - `AlignedStorage`: owning, mutable, 64-byte-aligned storage.
 * - `ReadOnlyStorageView`: non-owning read-only byte storage.
 */

#include <pixie/storage/aligned.h>
#include <pixie/storage/read_only_view.h>
