#pragma once

/**
 * @file storage.h
 * @brief Common interface for byte-addressable storage.
 *
 * Include `<pixie/storage/implementations.h>` to use Pixie's concrete storage
 * types.
 */

#include <pixie/bit_stream.h>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

namespace pixie {

/**
 * @brief CRTP facade for byte-addressable storage.
 *
 * @tparam Impl Concrete storage implementation.
 */
template <class Impl>
class StorageBase {
 public:
  /** @brief Return the exposed storage size in bytes. */
  std::size_t size_bytes() const { return impl().size_bytes_impl(); }

  /** @brief Return the exposed storage size in bits. */
  std::size_t size_bits() const { return size_bytes() * 8; }

  /** @brief Check whether the storage is empty. */
  bool empty() const { return size_bytes() == 0; }

  /** @brief Return a read-only view of all exposed bytes. */
  std::span<const std::byte> as_bytes() const { return impl().as_bytes_impl(); }

  /**
   * @brief Return a read-only view as 16-bit words.
   * @throws std::invalid_argument if the data is misaligned or its size is not
   * divisible by the word size.
   */
  std::span<const std::uint16_t> as_words16() const {
    return as_words<std::uint16_t>();
  }

  /**
   * @brief Return a read-only view as 64-bit words.
   * @throws std::invalid_argument if the data is misaligned or its size is not
   * divisible by the word size.
   */
  std::span<const std::uint64_t> as_words64() const {
    return as_words<std::uint64_t>();
  }

  /** @brief Return a non-owning read-only view of all exposed bytes. */
  auto view() const { return impl().view_impl(0, size_bytes()); }

  /**
   * @brief Return a non-owning read-only byte subrange.
   * @param offset_bytes First byte in the view.
   * @param count_bytes Number of bytes in the view.
   * @throws std::out_of_range if the subrange is outside this storage.
   */
  auto view(std::size_t offset_bytes, std::size_t count_bytes) const {
    return impl().view_impl(offset_bytes, count_bytes);
  }

  /** @brief Serialize the exposed byte sequence with a size prefix. */
  void serialize(OutputBitStream& stream) const {
    stream << size_bytes();
    for (const std::byte byte : as_bytes()) {
      stream << static_cast<std::uint8_t>(byte);
    }
  }

  /** @brief Resize mutable storage to hold at least @p size_bits bits. */
  void resize(std::size_t size_bits)
    requires requires(Impl& value) { value.resize_impl(size_bits); }
  {
    impl().resize_impl(size_bits);
  }

  /** @brief Return writable storage bytes. */
  auto writable_bytes()
    requires requires(Impl& value) { value.writable_bytes_impl(); }
  {
    return impl().writable_bytes_impl();
  }

  /** @brief Return writable storage as 16-bit words. */
  auto writable_words16()
    requires requires(Impl& value) { value.writable_words16_impl(); }
  {
    return impl().writable_words16_impl();
  }

  /** @brief Return writable storage as 64-bit words. */
  auto writable_words64()
    requires requires(Impl& value) { value.writable_words64_impl(); }
  {
    return impl().writable_words64_impl();
  }

  /** @brief Return bytes reserved by an owning storage implementation. */
  std::size_t allocated_bytes() const
    requires requires(const Impl& value) { value.allocated_bytes_impl(); }
  {
    return impl().allocated_bytes_impl();
  }

  /** @brief Request release of unused reserved storage. */
  void shrink_to_fit()
    requires requires(Impl& value) { value.shrink_to_fit_impl(); }
  {
    impl().shrink_to_fit_impl();
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }

  /** @brief Return this facade as its mutable CRTP implementation. */
  Impl& impl() { return static_cast<Impl&>(*this); }

  template <class Word>
  std::span<const Word> as_words() const {
    const auto bytes = as_bytes();
    if (bytes.size() % sizeof(Word) != 0 ||
        reinterpret_cast<std::uintptr_t>(bytes.data()) % alignof(Word) != 0) {
      throw std::invalid_argument("Storage is not aligned to the word type");
    }
    return {reinterpret_cast<const Word*>(bytes.data()),
            bytes.size() / sizeof(Word)};
  }
};

/** @brief A concrete CRTP implementation of `StorageBase`. */
template <class Storage>
concept StorageImplementation =
    std::derived_from<Storage, StorageBase<Storage>>;

}  // namespace pixie
