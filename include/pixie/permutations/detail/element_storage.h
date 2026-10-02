#pragma once

#include <pixie/permutable_sequence.h>

#include <concepts>
#include <limits>
#include <type_traits>

namespace pixie::permutations::detail {

// Resolve storage and read semantics once for the facade and its CRTP base.
template <class T, ElementStorage Requested>
struct ElementStorageTraits {
  static_assert(std::is_object_v<T> && !std::is_array_v<T> &&
                std::same_as<T, std::remove_cv_t<T>>);
  static constexpr bool packable = std::is_integral_v<T> &&
                                   std::is_unsigned_v<T> &&
                                   std::numeric_limits<T>::digits <= 64;
  static_assert(Requested == ElementStorage::automatic ||
                Requested == ElementStorage::packed ||
                Requested == ElementStorage::indirect);
  static_assert(Requested != ElementStorage::packed || packable,
                "Packed elements must be bool or unsigned integers <=64 bits");
  static constexpr ElementStorage storage =
      Requested == ElementStorage::automatic
          ? (packable ? ElementStorage::packed : ElementStorage::indirect)
          : Requested;
  static constexpr bool packed = storage == ElementStorage::packed;
  static_assert(
      packed || (std::is_nothrow_move_constructible_v<T> &&
                 std::is_nothrow_destructible_v<T>),
      "Indirect elements require noexcept move construction/destruction");
  using const_reference = std::conditional_t<packed, T, const T&>;
};

}  // namespace pixie::permutations::detail
