#pragma once

/**
 * @file implementations.h
 * @brief Native Pixie storage-backed file-archive types.
 *
 * - `FileArchiveIndex<Storage>`: storage-parameterized archive.
 * - `FileArchive`: owning aligned-storage alias.
 * - `FileArchiveView`: non-owning read-only storage-view alias.
 */

#include <pixie/file_archive.h>
