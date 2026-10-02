#pragma once

/**
 * @file sequence_options.h
 * @brief Shared compile-time sequence tuning options.
 * @details This header contains no storage engine or concrete implementation.
 */
namespace pixie {
/**
 * @brief Representation of tree child lengths.
 * @details Cumulative stores prefix ends; individual stores separate lengths.
 * Both preserve identical public sequence semantics.
 */
enum class LengthLayout { cumulative, individual };
}  // namespace pixie
