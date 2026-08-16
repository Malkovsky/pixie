#pragma once

#include <pixie/serialization.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pixie::detail {

inline constexpr std::uint8_t kLittleEndianMarker = 1;

template <std::size_t Size>
void write_magic(BinaryWriter& stream,
                 const std::array<std::uint8_t, Size>& magic) {
  for (const std::uint8_t byte : magic) {
    stream.write_u8(byte);
  }
}

template <std::size_t Size>
void require_magic(BinaryReader& reader,
                   const std::array<std::uint8_t, Size>& expected) {
  for (const std::uint8_t byte : expected) {
    const std::size_t offset = reader.byte_offset();
    if (reader.read_u8() != byte) {
      throw SerializationError("Serialized artifact has the wrong magic",
                               offset);
    }
  }
}

template <std::integral T>
void write_integral(BinaryWriter& stream, T value) {
  if constexpr (std::same_as<T, bool>) {
    stream.write_u8(static_cast<std::uint8_t>(value));
  } else if constexpr (std::same_as<T, std::size_t>) {
    stream.write_size(value);
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 1) {
    stream.write_u8(static_cast<std::uint8_t>(value));
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 2) {
    stream.write_u16(static_cast<std::uint16_t>(value));
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 4) {
    stream.write_u32(static_cast<std::uint32_t>(value));
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 8) {
    stream.write_u64(static_cast<std::uint64_t>(value));
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 1) {
    stream.write_i8(static_cast<std::int8_t>(value));
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 2) {
    stream.write_i16(static_cast<std::int16_t>(value));
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 4) {
    stream.write_i32(static_cast<std::int32_t>(value));
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 8) {
    stream.write_i64(static_cast<std::int64_t>(value));
  } else {
    static_assert(sizeof(T) == 0, "Unsupported serialized integer type");
  }
}

template <std::integral T>
T read_integral(BinaryReader& reader) {
  if constexpr (std::same_as<T, bool>) {
    const std::uint8_t value = reader.read_u8();
    if (value > 1) {
      throw SerializationError("Invalid serialized boolean value",
                               reader.byte_offset() - 1);
    }
    return value != 0;
  } else if constexpr (std::same_as<T, std::size_t>) {
    return reader.read_size();
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 1) {
    return static_cast<T>(reader.read_u8());
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 2) {
    return static_cast<T>(reader.read_u16());
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 4) {
    return static_cast<T>(reader.read_u32());
  } else if constexpr (std::is_unsigned_v<T> && sizeof(T) == 8) {
    return static_cast<T>(reader.read_u64());
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 1) {
    return static_cast<T>(reader.read_i8());
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 2) {
    return static_cast<T>(reader.read_i16());
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 4) {
    return static_cast<T>(reader.read_i32());
  } else if constexpr (std::is_signed_v<T> && sizeof(T) == 8) {
    return static_cast<T>(reader.read_i64());
  } else {
    static_assert(sizeof(T) == 0, "Unsupported serialized integer type");
  }
}

template <std::integral T>
inline constexpr std::size_t kSerializedIntegralBytes =
    std::same_as<T, std::size_t> ? sizeof(std::uint64_t) : sizeof(T);

template <std::integral T>
void write_vector(BinaryWriter& stream, std::span<const T> values) {
  stream.write_size(values.size());
  for (const T value : values) {
    write_integral(stream, value);
  }
}

template <std::integral T>
std::vector<T> read_vector(BinaryReader& reader) {
  const std::size_t count = reader.read_size();
  const std::vector<T> empty;
  if (count > empty.max_size()) {
    throw std::length_error("Serialized vector is too large");
  }
  if (count > reader.remaining() / kSerializedIntegralBytes<T>) {
    throw SerializationError("Truncated serialized vector",
                             reader.byte_offset());
  }
  std::vector<T> result;
  result.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    result.push_back(read_integral<T>(reader));
  }
  return result;
}

inline std::size_t checked_artifact_size(std::uint64_t encoded_size,
                                         std::size_t header_size,
                                         std::size_t available_size) {
  if (encoded_size > std::numeric_limits<std::size_t>::max()) {
    throw std::length_error("Serialized artifact size does not fit in size_t");
  }
  const std::size_t size = static_cast<std::size_t>(encoded_size);
  if (size < header_size || size > available_size) {
    throw std::invalid_argument("Truncated serialized artifact");
  }
  if (size % sizeof(std::uint64_t) != 0) {
    throw std::invalid_argument("Serialized artifact is not word padded");
  }
  return size;
}

}  // namespace pixie::detail
