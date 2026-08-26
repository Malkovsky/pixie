#pragma once

/**
 * @file packed_monotone.h
 * @brief Monotone integer vectors backed by runtime-width packing.
 */

#include <pixie/integer_vector/monotone.h>
#include <pixie/integer_vector/packed.h>

#include <cstdint>

namespace pixie {

/** @brief Owning validated monotone adapter over a packed integer vector. */
template <IntegerVectorValue Value = std::uint64_t>
using PackedMonotoneIntegerVector =
    MonotoneIntegerVector<PackedIntegerVector<Value>>;

/**
 * @brief Validated monotone adapter over a zero-copy packed-vector view.
 * @details A deserialized instance borrows its aligned immutable artifact
 * bytes.
 */
template <IntegerVectorValue Value = std::uint64_t>
using PackedMonotoneIntegerVectorView =
    MonotoneIntegerVector<PackedIntegerVectorView<Value>>;

}  // namespace pixie
