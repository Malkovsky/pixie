#pragma once

#include <pixie/storage.h>
#include <pixie/storage/read_only_view.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>

namespace pixie {

inline constexpr std::size_t kAlignedStorageLineBytes = 64;
inline constexpr std::size_t kAlignedStorageLineBits =
    kAlignedStorageLineBytes * 8;
inline constexpr std::size_t kAlignedStorageLineWords64 =
    kAlignedStorageLineBytes / sizeof(std::uint64_t);
inline constexpr std::size_t kAlignedStorageLineWords16 =
    kAlignedStorageLineBytes / sizeof(std::uint16_t);

/** @brief A 64-byte aligned storage block. */
struct alignas(kAlignedStorageLineBytes) CacheLine {
  std::array<std::byte, kAlignedStorageLineBytes> data{};
};

static_assert(alignof(CacheLine) == kAlignedStorageLineBytes);
static_assert(sizeof(CacheLine) == kAlignedStorageLineBytes);

/**
 * @brief Owning storage with a logical byte size and 64-byte-aligned backing.
 *
 * @details Construction fixes the logical byte size; contents remain writable.
 * Exposed byte and word views cover `ceil(size_bits / 8)` logical bytes, while
 * the backing allocation remains rounded up to complete cache lines. Assignment
 * or destruction invalidates its views. Use AlignedStorageBuilder when the
 * allocation needs to grow. Writes are visible through existing views.
 */
class AlignedStorage : public StorageBase<AlignedStorage> {
 public:
  AlignedStorage() = default;

  /** @brief Deep-copy bytes into a distinct allocation. */
  AlignedStorage(const AlignedStorage& other)
      : AlignedStorage(other.as_bytes()) {}
  /** @brief Transfer allocation without invalidating views into that
   * allocation. */
  AlignedStorage(AlignedStorage&& other) noexcept
      : logical_size_bytes_(std::exchange(other.logical_size_bytes_, 0)),
        data_(std::move(other.data_)) {}
  /** @brief Replace bytes with an independent copy; failure leaves this object
   * unchanged. */
  AlignedStorage& operator=(const AlignedStorage& other) {
    if (this != &other) {
      AlignedStorage copy(other);
      swap(copy);
    }
    return *this;
  }
  /** @brief Replace the allocation by transfer; invalidate views of the
   * previous allocation. */
  AlignedStorage& operator=(AlignedStorage&& other) noexcept {
    if (this != &other) {
      data_ = std::move(other.data_);
      logical_size_bytes_ = std::exchange(other.logical_size_bytes_, 0);
    }
    return *this;
  }
  /** @brief Exchange allocations and logical sizes without moving bytes. */
  void swap(AlignedStorage& other) noexcept {
    std::swap(logical_size_bytes_, other.logical_size_bytes_);
    data_.swap(other.data_);
  }

  /** @brief Construct zero-filled storage for at least @p size_bits bits. */
  explicit AlignedStorage(std::size_t size_bits)
      : logical_size_bytes_(bytes_for_bits(size_bits)),
        data_(allocate(logical_size_bytes_)) {}

  /**
   * @brief Copy @p bytes into aligned owning storage.
   * @details The new storage does not retain a view of the source bytes.
   */
  explicit AlignedStorage(std::span<const std::byte> bytes)
      : logical_size_bytes_(bytes.size()),
        data_(allocate(logical_size_bytes_)) {
    std::ranges::copy(bytes, reinterpret_cast<std::byte*>(data_.get()));
  }

  /** @brief Copy complete 64-bit words into aligned owning storage. */
  explicit AlignedStorage(std::span<const std::uint64_t> words)
      : AlignedStorage(std::as_bytes(words)) {}

  /** @brief Return the logical number of exposed bytes. */
  std::size_t size_bytes_impl() const { return logical_size_bytes_; }

  /** @brief Return the first logical byte position. */
  position_type begin_position_impl() const { return 0; }

  /** @brief Return the position one past the final logical byte. */
  position_type end_position_impl() const { return logical_size_bytes_; }

  /** @brief Return the logical number of exposed bytes. */
  std::size_t logical_size_bytes() const { return logical_size_bytes_; }

  /** @brief Return the cache-line-rounded backing size in bytes. */
  std::size_t padded_size_bytes() const {
    return lines_for_bytes(logical_size_bytes_) * kAlignedStorageLineBytes;
  }

  /**
   * @brief Return a non-owning view of the complete cache-line backing.
   * @details This view includes allocation padding and is intended for
   * word-oriented indexes that require a complete final word. Serialization
   * and ordinary storage views continue to expose only logical bytes.
   */
  ReadOnlyStorageView padded_view() const {
    return ReadOnlyStorageView(std::as_bytes(std::span<const CacheLine>(
        data_.get(), lines_for_bytes(logical_size_bytes_))));
  }

  /** @brief Return the logical bytes as a read-only span. */
  std::span<const std::byte> as_bytes_impl() const {
    return std::as_bytes(std::span<const CacheLine>(
                             data_.get(), lines_for_bytes(logical_size_bytes_)))
        .first(logical_size_bytes_);
  }

  /** @brief Return a checked logical byte range as one physical segment. */
  SplitSpan<const std::byte> segments_impl(std::size_t offset_bytes,
                                           std::size_t count_bytes) const {
    return SplitSpan(as_bytes_impl().subspan(offset_bytes, count_bytes));
  }

  /** @brief Return a checked logical byte range as one writable segment. */
  SplitSpan<std::byte> segments_impl(std::size_t offset_bytes,
                                     std::size_t count_bytes) {
    return SplitSpan(writable_bytes_impl().subspan(offset_bytes, count_bytes));
  }

  /** @brief Return writable logical bytes, excluding allocation padding. */
  std::span<std::byte> writable_bytes_impl() {
    return std::as_writable_bytes(as_lines()).first(logical_size_bytes_);
  }

  /** @brief Return writable complete 16-bit words within the logical bytes. */
  std::span<std::uint16_t> writable_words16_impl() {
    return {reinterpret_cast<std::uint16_t*>(data_.get()),
            logical_size_bytes_ / sizeof(std::uint16_t)};
  }

  /** @brief Return writable complete 64-bit words within the logical bytes. */
  std::span<std::uint64_t> writable_words64_impl() {
    return {reinterpret_cast<std::uint64_t*>(data_.get()),
            logical_size_bytes_ / sizeof(std::uint64_t)};
  }

  /** @brief Return a checked read-only byte subrange. */
  ReadOnlyStorageView view_impl(std::size_t offset_bytes,
                                std::size_t count_bytes) const {
    if (offset_bytes > size_bytes_impl() ||
        count_bytes > size_bytes_impl() - offset_bytes) {
      throw std::out_of_range("Storage view is outside the allocation");
    }
    return ReadOnlyStorageView(
        as_bytes_impl().subspan(offset_bytes, count_bytes));
  }

  /** @brief Return exact cache-line-rounded allocation size in bytes. */
  std::size_t allocated_bytes_impl() const { return padded_size_bytes(); }

  /** @brief Return writable cache-line blocks, including allocation padding. */
  std::span<CacheLine> as_lines() {
    return {data_.get(), lines_for_bytes(logical_size_bytes_)};
  }

  /** @brief Return read-only cache-line blocks. */
  std::span<const CacheLine> as_lines() const {
    return {data_.get(), lines_for_bytes(logical_size_bytes_)};
  }

  /** @brief Restore an owning copy of one size-prefixed byte sequence. */
  static AlignedStorage deserialize_impl(BinaryReader& reader) {
    const std::size_t size = reader.read_size();
    if (size > std::numeric_limits<std::size_t>::max() / 8) {
      throw std::length_error("Serialized aligned storage is too large");
    }
    return AlignedStorage(reader.read_bytes(size));
  }

 private:
  static constexpr std::size_t bytes_for_bits(std::size_t size_bits) {
    return size_bits / 8 + (size_bits % 8 != 0);
  }

  static constexpr std::size_t lines_for_bytes(std::size_t size_bytes) {
    return size_bytes / kAlignedStorageLineBytes +
           (size_bytes % kAlignedStorageLineBytes != 0);
  }

  std::size_t logical_size_bytes_ = 0;
  std::unique_ptr<CacheLine[]> data_;

  static std::unique_ptr<CacheLine[]> allocate(std::size_t bytes) {
    const auto lines = lines_for_bytes(bytes);
    if (lines > std::numeric_limits<size_t>::max() / kAlignedStorageLineBytes) {
      throw std::length_error("Aligned storage allocation is too large");
    }
    return lines ? std::make_unique<CacheLine[]>(lines) : nullptr;
  }
};

static_assert(sizeof(AlignedStorage) == sizeof(void*) + sizeof(std::size_t));

}  // namespace pixie
