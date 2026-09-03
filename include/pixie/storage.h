#pragma once

/**
 * @file storage.h
 * @brief Common interface for byte-addressable storage.
 *
 * Include `<pixie/storage/implementations.h>` to use Pixie's concrete storage
 * types.
 */

#include <pixie/serialization.h>
#include <pixie/split_span.h>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

namespace pixie {

/**
 * @brief CRTP facade for byte-addressable storage.
 *
 * @details `Impl` must provide `size_bytes_impl()` and the following required
 * extension points.
 *
 * @par Required: `begin_position_impl() const`
 * Returns the absolute logical position of the first exposed byte as a
 * `position_type`. The result may advance when an implementation evicts bytes,
 * but positions of retained bytes do not change.
 *
 * @par Required: `end_position_impl() const`
 * Returns the absolute logical position one past the last exposed byte as a
 * `position_type`. It must be at least `begin_position_impl()`, and their
 * difference must be representable by `std::size_t` and equal
 * `size_bytes_impl()`.
 *
 * @par Required: `segments_impl(position, count_bytes) const`
 * Accepts a `position_type` and a `std::size_t` and returns
 * `SplitSpan<const std::byte>` containing exactly the bytes in `[position,
 * position + count_bytes)`. The physical spans are ordered by logical position
 * and their combined size equals `count_bytes`. Before calling this hook, the
 * facade verifies without overflow that `position` is in the closed range from
 * `begin_position_impl()` through `end_position_impl()` and that `count_bytes
 * <= end_position_impl() - position`; an empty range at the end position is
 * valid. The returned spans borrow the implementation's backing storage. The
 * implementation must document the backing storage's ownership requirements
 * and the operations that invalidate its spans.
 *
 * @par Optional: `segments_impl(position, count_bytes)`
 * A mutable implementation may provide this extension point to enable
 * writable segment access. It accepts the same parameter types as the const
 * overload and returns `SplitSpan<std::byte>`. It has the same range, ordering,
 * size, lifetime, and invalidation contract as the const overload, and writes
 * through the returned spans modify the corresponding logical bytes. If the
 * implementation also provides `prepare_segments_impl(position, count_bytes)`,
 * the facade calls that hook before validating the range. The preparation hook
 * must make the complete requested range available or throw without changing
 * the implementation.
 *
 * @tparam Impl Concrete storage implementation.
 */
template <class Impl>
class StorageBase : public SerializationBase<Impl> {
 public:
  /** @brief Monotonic logical byte position used to address storage ranges. */
  using position_type = std::uint64_t;

  /**
   * @brief Return the logical exposed storage size in bytes.
   * @details An owning implementation may reserve or pad more memory; use
   * `allocated_bytes()` when that physical allocation size is required.
   */
  std::size_t size_bytes() const { return impl().size_bytes_impl(); }

  /** @brief Return the exposed storage size in bits. */
  std::size_t size_bits() const { return size_bytes() * 8; }

  /** @brief Check whether the storage is empty. */
  bool empty() const { return size_bytes() == 0; }

  /** @brief Return the first logical byte position currently exposed. */
  position_type begin_position() const { return impl().begin_position_impl(); }

  /** @brief Return the position one past the last logical byte exposed. */
  position_type end_position() const { return impl().end_position_impl(); }

  /** @brief Return whether a complete logical range is currently exposed. */
  bool contains(position_type position, std::size_t count_bytes) const {
    const position_type begin = begin_position();
    const position_type end = end_position();
    return position >= begin && position <= end &&
           count_bytes <= end - position;
  }

  /** @brief Return a contiguous read-only view of all logical exposed bytes. */
  std::span<const std::byte> as_bytes() const
    requires requires(const Impl& value) { value.as_bytes_impl(); }
  {
    return impl().as_bytes_impl();
  }

  /**
   * @brief Return all logical bytes as one or two writable physical spans.
   * @details Available only for mutable storage implementations. Mutating or
   * resizing the storage may invalidate the returned descriptor.
   */
  SplitSpan<std::byte> segments()
    requires requires(Impl& value) {
      {
        value.segments_impl(position_type{}, std::size_t{})
      } -> std::same_as<SplitSpan<std::byte>>;
    }
  {
    return segments(begin_position(), size_bytes());
  }

  /**
   * @brief Return all logical bytes as one or two physical spans.
   * @details Contiguous storage returns one segment. A ring-backed storage can
   * return its tail followed by its head without allocation or copying.
   */
  SplitSpan<const std::byte> segments() const {
    return segments(begin_position(), size_bytes());
  }

  /**
   * @brief Return a checked writable logical byte range as one or two spans.
   * @details An implementation with a preparation hook may make a future range
   * available before checking it. Any newly exposed bytes remain the caller's
   * responsibility to initialize.
   * @param position First logical byte position in the range.
   * @param count_bytes Number of logical bytes in the range.
   * @throws std::out_of_range if the range cannot be made available.
   */
  SplitSpan<std::byte> segments(position_type position, std::size_t count_bytes)
    requires requires(Impl& value) {
      {
        value.segments_impl(position_type{}, std::size_t{})
      } -> std::same_as<SplitSpan<std::byte>>;
    }
  {
    if constexpr (requires(Impl& value) {
                    value.prepare_segments_impl(position_type{}, std::size_t{});
                  }) {
      impl().prepare_segments_impl(position, count_bytes);
    }
    validate_range(position, count_bytes);
    return impl().segments_impl(position, count_bytes);
  }

  /**
   * @brief Return a checked logical byte range as one or two physical spans.
   * @param position First logical byte position in the range.
   * @param count_bytes Number of logical bytes in the range.
   * @throws std::out_of_range if the range is outside this storage.
   */
  SplitSpan<const std::byte> segments(position_type position,
                                      std::size_t count_bytes) const {
    validate_range(position, count_bytes);
    return impl().segments_impl(position, count_bytes);
  }

  /**
   * @brief Return a read-only view as 16-bit words.
   * @throws std::invalid_argument if the data is misaligned or its size is not
   * divisible by the word size.
   */
  std::span<const std::uint16_t> as_words16() const
    requires requires(const Impl& value) { value.as_bytes_impl(); }
  {
    return as_words<std::uint16_t>();
  }

  /**
   * @brief Return a read-only view as 64-bit words.
   * @throws std::invalid_argument if the data is misaligned or its size is not
   * divisible by the word size.
   */
  std::span<const std::uint64_t> as_words64() const
    requires requires(const Impl& value) { value.as_bytes_impl(); }
  {
    return as_words<std::uint64_t>();
  }

  /** @brief Return a non-owning read-only view of all exposed bytes. */
  auto view() const
    requires requires(const Impl& value) {
      value.view_impl(std::size_t{}, std::size_t{});
    }
  {
    return impl().view_impl(0, size_bytes());
  }

  /**
   * @brief Return a non-owning read-only byte subrange.
   * @param offset_bytes First byte in the view.
   * @param count_bytes Number of bytes in the view.
   * @throws std::out_of_range if the subrange is outside this storage.
   */
  auto view(std::size_t offset_bytes, std::size_t count_bytes) const
    requires requires(const Impl& value) {
      value.view_impl(std::size_t{}, std::size_t{});
    }
  {
    return impl().view_impl(offset_bytes, count_bytes);
  }

  /**
   * @brief Serialize the exposed bytes with a 64-bit little-endian size prefix.
   */
  void serialize_impl(BinaryWriter& writer) const {
    writer.write_size(size_bytes());
    for (const std::span<const std::byte> segment : segments()) {
      writer.write_bytes(segment);
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
  /** @brief Validate a logical byte range without overflowing. */
  void validate_range(position_type position, std::size_t count_bytes) const {
    if (!contains(position, count_bytes)) {
      throw std::out_of_range("Storage range is outside the working window");
    }
  }

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
