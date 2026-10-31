#pragma once

#include <cstddef>

namespace pixie {
/** @brief Default allowance for an external traversal context and its scratch.
 */
inline constexpr std::size_t kDefaultTraversalBytes = 512 * 1024;
}  // namespace pixie
