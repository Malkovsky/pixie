#pragma once

#include <pixie/storage.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace pixie::detail::compressed {

template <class Storage>
concept MutableMetadataStorage =
    StorageImplementation<Storage> && std::default_initializable<Storage> &&
    requires(Storage& storage, const Storage& readonly) {
      storage.resize(std::size_t{});
      { storage.writable_bytes() } -> std::same_as<std::span<std::byte>>;
      { readonly.as_bytes() } -> std::same_as<std::span<const std::byte>>;
    };

// Fixed-size typed metadata; validate storage once, before initializing bytes.
template <class T, StorageImplementation Storage>
class MetadataArray {
  static_assert(std::is_trivially_copyable_v<T>);

 public:
  MetadataArray() = default;

  explicit MetadataArray(std::size_t count)
    requires MutableMetadataStorage<Storage>
  {
    std::ranges::fill(allocate(count), std::byte{0});
  }

  explicit MetadataArray(std::span<const T> values)
    requires MutableMetadataStorage<Storage>
  {
    auto bytes = allocate(values.size());
    if (!values.empty()) {
      std::memcpy(bytes.data(), values.data(), values.size_bytes());
    }
  }

  const T* data() const {
    return reinterpret_cast<const T*>(storage_.as_bytes().data());
  }

  T* data()
    requires MutableMetadataStorage<Storage>
  {
    return reinterpret_cast<T*>(storage_.writable_bytes().data());
  }

  const T& operator[](std::size_t index) const { return data()[index]; }

  T& operator[](std::size_t index)
    requires MutableMetadataStorage<Storage>
  {
    return data()[index];
  }

  std::size_t size() const { return storage_.size_bytes() / sizeof(T); }

  bool empty() const { return size() == 0; }

  std::size_t allocated_bytes() const
    requires requires(const Storage& storage) { storage.allocated_bytes(); }
  {
    return storage_.allocated_bytes();
  }

 private:
  Storage storage_;

  std::span<std::byte> allocate(std::size_t count)
    requires MutableMetadataStorage<Storage>
  {
    if (count > std::numeric_limits<std::size_t>::max() / (sizeof(T) * 8)) {
      throw std::length_error("Compressed metadata allocation too large");
    }
    if (count == 0) {
      return {};
    }
    const auto byte_count = count * sizeof(T);
    storage_.resize(byte_count * 8);
    auto bytes = storage_.writable_bytes();
    const auto readonly = storage_.as_bytes();
    if (bytes.size() != byte_count || readonly.size() != byte_count ||
        readonly.data() != bytes.data() ||
        storage_.size_bytes() != byte_count ||
        reinterpret_cast<std::uintptr_t>(bytes.data()) %
                std::max(alignof(T), alignof(std::uint64_t)) !=
            0) {
      throw std::invalid_argument(
          "Compressed metadata storage must be contiguous and word aligned");
    }
    return bytes;
  }
};

}  // namespace pixie::detail::compressed
