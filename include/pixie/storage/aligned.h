#pragma once

#include <pixie/storage.h>
#include <pixie/storage/read_only_view.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

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
 * @details Construction and resize accept a logical bit count. Exposed byte
 * and word views cover `ceil(size_bits / 8)` logical bytes, while the backing
 * allocation remains rounded up to complete cache lines. Resizing or
 * destroying this object invalidates its read-only views.
 */
class AlignedStorage : public StorageBase<AlignedStorage> {
 public:
  AlignedStorage() = default;

  /** @brief Construct storage for at least @p size_bits bits. */
  explicit AlignedStorage(std::size_t size_bits)
      : logical_size_bytes_(bytes_for_bits(size_bits)),
        data_(lines_for_bits(size_bits)) {}

  /** @brief Copy complete 64-bit words into aligned owning storage. */
  explicit AlignedStorage(std::span<const std::uint64_t> words)
      : AlignedStorage(bit_size_for_words(words.size())) {
    std::copy(words.begin(), words.end(), writable_words64_impl().begin());
  }

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
    return data_.size() * kAlignedStorageLineBytes;
  }

  /**
   * @brief Return a non-owning view of the complete cache-line backing.
   * @details This view includes allocation padding and is intended for
   * word-oriented indexes that require a complete final word. Serialization
   * and ordinary storage views continue to expose only logical bytes.
   */
  ReadOnlyStorageView padded_view() const {
    return ReadOnlyStorageView(
        std::as_bytes(std::span<const CacheLine>(data_)));
  }

  /** @brief Return the logical bytes as a read-only span. */
  std::span<const std::byte> as_bytes_impl() const {
    return std::as_bytes(std::span<const CacheLine>(data_))
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

  /** @brief Resize to hold at least @p size_bits bits. */
  void resize_impl(std::size_t size_bits) {
    data_.resize(lines_for_bits(size_bits));
    logical_size_bytes_ = bytes_for_bits(size_bits);
  }

  /** @brief Return writable logical bytes. */
  std::span<std::byte> writable_bytes_impl() {
    return std::as_writable_bytes(std::span<CacheLine>(data_))
        .first(logical_size_bytes_);
  }

  /** @brief Return writable logical storage as 16-bit words. */
  std::span<std::uint16_t> writable_words16_impl() {
    return {reinterpret_cast<std::uint16_t*>(data_.data()),
            logical_size_bytes_ / sizeof(std::uint16_t)};
  }

  /** @brief Return writable logical storage as 64-bit words. */
  std::span<std::uint64_t> writable_words64_impl() {
    return {reinterpret_cast<std::uint64_t*>(data_.data()),
            logical_size_bytes_ / sizeof(std::uint64_t)};
  }

  /** @brief Return bytes reserved by the underlying vector. */
  std::size_t allocated_bytes_impl() const {
    return data_.capacity() * kAlignedStorageLineBytes;
  }

  /** @brief Request release of unused vector capacity. */
  void shrink_to_fit_impl() { data_.shrink_to_fit(); }

  /** @brief Return mutable cache-line blocks. */
  std::span<CacheLine> as_lines() { return data_; }

  /** @brief Return read-only cache-line blocks. */
  std::span<const CacheLine> as_lines() const { return data_; }

  /** @brief Restore an owning copy of one size-prefixed byte sequence. */
  static AlignedStorage deserialize_impl(BinaryReader& reader) {
    const std::size_t size = reader.read_size();
    if (size > std::numeric_limits<std::size_t>::max() / 8) {
      throw std::length_error("Serialized aligned storage is too large");
    }
    AlignedStorage result(size * 8);
    const std::span<const std::byte> bytes = reader.read_bytes(size);
    std::ranges::copy(bytes, result.writable_bytes().begin());
    return result;
  }

 private:
  static std::size_t bit_size_for_words(std::size_t word_count) {
    constexpr std::size_t kWordBits =
        std::numeric_limits<std::uint64_t>::digits;
    if (word_count > std::numeric_limits<std::size_t>::max() / kWordBits) {
      throw std::length_error("Aligned storage word sequence is too large");
    }
    return word_count * kWordBits;
  }

  static constexpr std::size_t bytes_for_bits(std::size_t size_bits) {
    return size_bits / 8 + (size_bits % 8 != 0);
  }

  static constexpr std::size_t lines_for_bits(std::size_t size_bits) {
    return size_bits / kAlignedStorageLineBits +
           (size_bits % kAlignedStorageLineBits != 0);
  }

  std::size_t logical_size_bytes_ = 0;
  std::vector<CacheLine> data_;
};

}  // namespace pixie
