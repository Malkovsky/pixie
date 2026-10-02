#pragma once

/**
 * @file sequence_implementations.h
 * @brief Production permutable-sequence catalog for benchmarks.
 * @details Consumers include pixie/permutations/sequence.h directly.
 * Experimental implementations and detailed comparison notes live under
 * pixie/experimental/ and are included explicitly by benchmark translation
 * units. The current snapshot is retained in
 * pixie/experimental/benchmarks/sequence_snapshot.md.
 */
#include <pixie/permutable_sequence.h>
#include <pixie/permutations/sequence.h>
#ifdef PIXIE_IMMER_SUPPORT
#include <pixie/permutations/immer_sequence.h>
#endif
