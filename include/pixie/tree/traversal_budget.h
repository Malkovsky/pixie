#pragma once

#include <cstddef>

namespace pixie::tree {
/** @brief Default allowance for an external traversal context and its scratch.
 */
inline constexpr std::size_t kDefaultTraversalBytes = 512 * 1024;
}  // namespace pixie::tree
