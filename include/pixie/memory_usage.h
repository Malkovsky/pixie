#pragma once

#include <concepts>
#include <cstddef>
#include <optional>
#include <vector>

namespace pixie {

/**
 * @brief Return bytes reserved by a vector's heap buffer.
 *
 * @details This intentionally uses capacity rather than size because the goal
 * is to report practical owned memory after construction. The vector control
 * block itself is counted by the owning object's `sizeof(*this)`.
 */
template <class T, class Allocator>
std::size_t vector_capacity_bytes(const std::vector<T, Allocator>& values) {
  return values.capacity() * sizeof(T);
}

/**
 * @brief Concept for types that expose total owned memory usage in bytes.
 */
template <class T>
concept HasMemoryUsageBytes = requires(const T& value) {
  { value.memory_usage_bytes() } -> std::convertible_to<std::size_t>;
};

/**
 * @brief Return heap bytes owned below an inline nested object.
 *
 * @details `outer.sizeof(*this)` already counts the inline nested object
 * storage. This helper subtracts that inline storage from the nested object's
 * total memory usage and leaves only buffers owned below it.
 */
template <HasMemoryUsageBytes T>
std::size_t nested_owned_memory_bytes(const T& value) {
  const std::size_t total = value.memory_usage_bytes();
  return total > sizeof(T) ? total - sizeof(T) : 0;
}

/**
 * @brief Return heap bytes owned below an engaged optional nested object.
 */
template <HasMemoryUsageBytes T>
std::size_t optional_nested_owned_memory_bytes(const std::optional<T>& value) {
  return value.has_value() ? nested_owned_memory_bytes(*value) : 0;
}

}  // namespace pixie
