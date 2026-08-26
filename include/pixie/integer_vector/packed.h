#pragma once

/**
 * @file packed.h
 * @brief Runtime-width packed immutable integer vectors.
 */

#include <pixie/detail/serialization.h>
#include <pixie/integer_vector.h>
#include <pixie/serialization.h>
#include <pixie/storage.h>
#include <pixie/storage/aligned.h>
#include <pixie/storage/read_only_view.h>

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace pixie {

/**
 * @brief Immutable runtime-width packed integer vector over configurable
 * storage.
 *
 * @tparam Storage `AlignedStorage` for an owner or `ReadOnlyStorageView` for a
 * zero-copy view.
 * @tparam Value Unsigned element type other than `bool`, no wider than 64 bits.
 *
 * @details Values are packed LSB-first into 64-bit words. An owning vector
 * stores cache-line-padded allocation internally, but only logical packed
 * words are serialized. A view borrows its serialized payload; the backing
 * bytes must remain alive, immutable, and aligned for the view's lifetime.
 */
template <StorageImplementation Storage,
          IntegerVectorValue Value = std::uint64_t>
  requires(std::same_as<Storage, AlignedStorage> ||
           std::same_as<Storage, ReadOnlyStorageView>)
class BasicPackedIntegerVector
    : public IntegerVectorBase<BasicPackedIntegerVector<Storage, Value>, Value>,
      public SerializationBase<BasicPackedIntegerVector<Storage, Value>> {
 public:
  /** @brief Stored unsigned integer type. */
  using value_type = Value;

  /** @brief Type used for element counts, widths, and positions. */
  using size_type = std::size_t;

  /** @brief Construct an empty width-zero vector. */
  BasicPackedIntegerVector() = default;

  /**
   * @brief Pack @p values using their minimum required width.
   * @details Empty and all-zero inputs use width zero and no payload.
   * @throws std::length_error if the packed dimensions overflow `size_t`.
   */
  explicit BasicPackedIntegerVector(std::span<const value_type> values)
    requires std::same_as<Storage, AlignedStorage>
      : BasicPackedIntegerVector(values,
                                 inferred_width(values),
                                 InferredTag{}) {}

  /**
   * @brief Pack @p values after validating an explicit maximum @p width.
   *
   * @details Width must be in `[0, digits(value_type)]`. Every value must fit;
   * values are never truncated. Empty and all-zero inputs are normalized to
   * canonical width zero.
   *
   * @throws std::invalid_argument if the width is invalid or a value does not
   * fit.
   * @throws std::length_error if the packed dimensions overflow `size_t`.
   */
  BasicPackedIntegerVector(std::span<const value_type> values, size_type width)
    requires std::same_as<Storage, AlignedStorage>
      : BasicPackedIntegerVector(values, width, ExplicitTag{}) {}

  /** @brief Return the logical runtime width in bits, in `[0, digits(Value)]`.
   */
  size_type width() const { return width_; }

  /** @brief Return the logical number of elements. */
  size_type size_impl() const { return size_; }

  /** @brief Read one valid zero-based position without bounds checking. */
  value_type value_at_impl(size_type position) const {
    return read_field(logical_words(), position, width_);
  }

  /** @brief Copy a range already validated by the public facade. */
  void copy_to_impl(size_type begin, std::span<value_type> output) const {
    for (size_type i = 0; i < output.size(); ++i) {
      output[i] = value_at_impl(begin + i);
    }
  }

  /**
   * @brief Return inline bytes plus owned storage capacity.
   * @details A read-only view excludes all borrowed backing bytes.
   */
  size_type memory_usage_bytes_impl() const {
    if constexpr (std::same_as<Storage, AlignedStorage>) {
      return sizeof(*this) + storage_.allocated_bytes();
    } else {
      return sizeof(*this);
    }
  }

  /**
   * @brief Write a version-1 canonical little-endian packed-vector artifact.
   *
   * @details Exactly `ceil(size() * width() / 64)` logical words are written;
   * owner allocation padding is excluded. The artifact can later be restored
   * as either an independent owner or a zero-copy read-only view.
   *
   * @throws std::invalid_argument if the artifact would not begin at an
   * eight-byte-aligned writer offset.
   */
  void serialize_impl(BinaryWriter& writer) const {
    if (writer.size_bytes() % alignof(std::uint64_t) != 0) {
      throw std::invalid_argument(
          "Packed integer-vector serialization requires an aligned offset");
    }

    const Dimensions dimensions = checked_dimensions(size_, width_);
    const std::size_t artifact_begin = writer.size_bytes();
    detail::write_magic(writer, kSerializationMagic);
    writer.write_u32(kSerializationVersion);
    writer.write_u8(detail::kLittleEndianMarker);
    writer.write_u8(static_cast<std::uint8_t>(kValueDigits));
    writer.write_u8(static_cast<std::uint8_t>(width_));
    writer.write_u8(0);
    const std::size_t size_position = writer.write_u64_placeholder();
    writer.write_size(size_);

    const auto words = logical_words();
    for (std::size_t i = 0; i < dimensions.word_count; ++i) {
      std::uint64_t word = words[i];
      if (i + 1 == dimensions.word_count && dimensions.bit_count % 64 != 0) {
        word &= low_bits_mask(dimensions.bit_count % 64);
      }
      writer.write_u64(word);
    }
    writer.patch_u64(size_position, static_cast<std::uint64_t>(
                                        writer.size_bytes() - artifact_begin));
  }

  /**
   * @brief Restore one packed-vector artifact and advance @p reader on success.
   *
   * @details Owner restoration decodes canonical words into independent
   * aligned storage. View restoration retains an exact span into the reader's
   * backing bytes, which must remain alive and immutable. Views require native
   * little-endian word order and eight-byte-aligned payload bytes. Quick mode
   * validates complete framing and safe dimensions; full mode also requires
   * unused high bits in the final logical word to be zero. Failure leaves
   * @p reader unchanged.
   *
   * @throws std::invalid_argument for malformed, incompatible, unaligned, or
   * noncanonical input.
   * @throws std::length_error for unrepresentable dimensions.
   */
  static BasicPackedIntegerVector deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation =
          DeserializationValidation::kQuick) {
    BinaryReader candidate = reader;
    const std::size_t available_size = candidate.remaining();
    detail::require_magic(candidate, kSerializationMagic);
    if (candidate.read_u32() != kSerializationVersion ||
        candidate.read_u8() != detail::kLittleEndianMarker ||
        candidate.read_u8() != kValueDigits) {
      throw std::invalid_argument(
          "Incompatible serialized packed integer vector");
    }
    const std::size_t width = candidate.read_u8();
    if (candidate.read_u8() != 0 || width > kValueDigits) {
      throw std::invalid_argument(
          "Invalid serialized packed integer-vector width or reserved field");
    }
    const std::size_t artifact_size = detail::checked_artifact_size(
        candidate.read_u64(), kSerializationHeaderBytes, available_size);
    const std::size_t count = candidate.read_size();
    const Dimensions dimensions = checked_dimensions(count, width);
    if ((count == 0 && width != 0) ||
        dimensions.byte_count != artifact_size - kSerializationHeaderBytes) {
      throw std::invalid_argument(
          "Serialized packed integer vector has inconsistent dimensions");
    }
    BinaryReader payload = candidate.read_subreader(dimensions.byte_count);

    if (validation == DeserializationValidation::kFull) {
      BinaryReader validation_reader(payload.remaining_bytes());
      bool has_nonzero_word = false;
      for (std::size_t i = 0; i < dimensions.word_count; ++i) {
        const std::uint64_t word = validation_reader.read_u64();
        has_nonzero_word = has_nonzero_word || word != 0;
        if (i + 1 == dimensions.word_count && dimensions.bit_count % 64 != 0 &&
            (word & ~low_bits_mask(dimensions.bit_count % 64)) != 0) {
          throw std::invalid_argument(
              "Serialized packed integer vector has non-zero unused bits");
        }
      }
      if (count != 0 && width != 0 && !has_nonzero_word) {
        throw std::invalid_argument(
            "Serialized all-zero integer vector has noncanonical width");
      }
    }

    Storage storage;
    if constexpr (std::same_as<Storage, AlignedStorage>) {
      storage = AlignedStorage(word_aligned_bit_count(dimensions.word_count));
      auto words = storage.writable_words64();
      for (std::size_t i = 0; i < dimensions.word_count; ++i) {
        words[i] = payload.read_u64();
      }
    } else {
      if constexpr (std::endian::native != std::endian::little) {
        throw std::invalid_argument(
            "Packed integer-vector views require little-endian word order");
      }
      const auto bytes = payload.read_bytes(dimensions.byte_count);
      if (dimensions.byte_count != 0 &&
          reinterpret_cast<std::uintptr_t>(bytes.data()) %
                  alignof(std::uint64_t) !=
              0) {
        throw std::invalid_argument(
            "Serialized packed integer-vector payload is not word aligned");
      }
      storage = ReadOnlyStorageView(bytes);
    }
    if (!payload.empty()) {
      throw std::invalid_argument(
          "Serialized packed integer vector has trailing payload bytes");
    }

    BasicPackedIntegerVector result(std::move(storage), count, width,
                                    LoadTag{});
    reader = candidate;
    return result;
  }

 private:
  struct Dimensions {
    std::size_t bit_count;
    std::size_t word_count;
    std::size_t byte_count;
  };
  struct InferredTag {};
  struct ExplicitTag {};
  struct LoadTag {};

  inline static constexpr std::array<std::uint8_t, 8> kSerializationMagic = {
      'P', 'X', 'I', 'N', 'T', 'V', 'E', 'C'};
  static constexpr std::uint32_t kSerializationVersion = 1;
  static constexpr std::size_t kSerializationHeaderBytes = 32;
  static constexpr std::size_t kValueDigits =
      std::numeric_limits<value_type>::digits;

  BasicPackedIntegerVector(std::span<const value_type> values,
                           size_type width,
                           InferredTag)
    requires std::same_as<Storage, AlignedStorage>
      : BasicPackedIntegerVector(values, width, ExplicitTag{}) {}

  BasicPackedIntegerVector(std::span<const value_type> values,
                           size_type width,
                           ExplicitTag)
    requires std::same_as<Storage, AlignedStorage>
      : size_(values.size()) {
    if (width > kValueDigits) {
      throw std::invalid_argument("Packed integer-vector width is too large");
    }
    bool all_zero = true;
    for (const value_type value : values) {
      all_zero = all_zero && value == 0;
      if (required_width(value) > width) {
        throw std::invalid_argument(
            "Packed integer-vector value does not fit the width");
      }
    }
    width_ = all_zero ? 0 : width;
    const Dimensions dimensions = checked_dimensions(size_, width_);
    storage_ = AlignedStorage(word_aligned_bit_count(dimensions.word_count));
    auto words = storage_.writable_words64();
    for (std::size_t i = 0; i < values.size(); ++i) {
      write_field(words, i, width_, values[i]);
    }
  }

  BasicPackedIntegerVector(Storage storage,
                           size_type size,
                           size_type width,
                           LoadTag)
      : storage_(std::move(storage)), size_(size), width_(width) {}

  static size_type required_width(value_type value) {
    return static_cast<size_type>(std::bit_width(value));
  }

  static constexpr std::uint64_t low_bits_mask(size_type count) {
    return count >= 64 ? std::numeric_limits<std::uint64_t>::max()
                       : (std::uint64_t{1} << count) - 1;
  }

  static size_type inferred_width(std::span<const value_type> values) {
    size_type width = 0;
    for (const value_type value : values) {
      width = std::max(width, required_width(value));
    }
    return width;
  }

  static Dimensions checked_dimensions(size_type count, size_type width) {
    if (width > kValueDigits) {
      throw std::invalid_argument("Packed integer-vector width is too large");
    }
    if (width != 0 && count > std::numeric_limits<size_type>::max() / width) {
      throw std::length_error("Packed integer-vector bit count is too large");
    }
    const size_type bit_count = count * width;
    const size_type word_count = bit_count == 0 ? 0 : 1 + (bit_count - 1) / 64;
    if (word_count >
        std::numeric_limits<size_type>::max() / sizeof(std::uint64_t)) {
      throw std::length_error("Packed integer-vector payload is too large");
    }
    return {bit_count, word_count, word_count * sizeof(std::uint64_t)};
  }

  static size_type word_aligned_bit_count(size_type word_count) {
    constexpr size_type kWordBits = std::numeric_limits<std::uint64_t>::digits;
    if (word_count > std::numeric_limits<size_type>::max() / kWordBits) {
      throw std::length_error("Packed integer-vector storage is too large");
    }
    return word_count * kWordBits;
  }

  std::span<const std::uint64_t> logical_words() const {
    const Dimensions dimensions = checked_dimensions(size_, width_);
    return storage_.as_words64().first(dimensions.word_count);
  }

  static value_type read_field(std::span<const std::uint64_t> words,
                               size_type position,
                               size_type width) {
    if (width == 0) {
      return 0;
    }
    if (width == 64) {
      return static_cast<value_type>(words[position]);
    }
    const size_type bit_position = position * width;
    const size_type word = bit_position / 64;
    const size_type offset = bit_position % 64;
    std::uint64_t value = words[word] >> offset;
    if (offset + width > 64) {
      value |= words[word + 1] << (64 - offset);
    }
    return static_cast<value_type>(value & low_bits_mask(width));
  }

  static void write_field(std::span<std::uint64_t> words,
                          size_type position,
                          size_type width,
                          value_type value) {
    if (width == 0) {
      return;
    }
    if (width == 64) {
      words[position] = static_cast<std::uint64_t>(value);
      return;
    }
    const size_type bit_position = position * width;
    const size_type word = bit_position / 64;
    const size_type offset = bit_position % 64;
    words[word] |= static_cast<std::uint64_t>(value) << offset;
    if (offset + width > 64) {
      words[word + 1] |= static_cast<std::uint64_t>(value) >> (64 - offset);
    }
  }

  Storage storage_;
  size_type size_ = 0;
  size_type width_ = 0;
};

/** @brief Owning aligned packed integer vector. */
template <IntegerVectorValue Value = std::uint64_t>
using PackedIntegerVector = BasicPackedIntegerVector<AlignedStorage, Value>;

/**
 * @brief Zero-copy read-only packed integer-vector view.
 * @details A deserialized view borrows its aligned immutable artifact bytes.
 */
template <IntegerVectorValue Value = std::uint64_t>
using PackedIntegerVectorView =
    BasicPackedIntegerVector<ReadOnlyStorageView, Value>;

}  // namespace pixie
