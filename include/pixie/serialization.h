#pragma once

/**
 * @file serialization.h
 * @brief Checked byte-oriented binary serialization primitives.
 */

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Error raised while decoding malformed serialized data.
 */
class SerializationError : public std::invalid_argument {
 public:
  /**
   * @brief Construct an error for @p byte_offset in the input artifact.
   */
  SerializationError(std::string message, std::size_t byte_offset)
      : std::invalid_argument(std::move(message) + " at byte " +
                              std::to_string(byte_offset)),
        byte_offset_(byte_offset) {}

  /** @brief Return the byte offset associated with the error. */
  std::size_t byte_offset() const noexcept { return byte_offset_; }

 private:
  std::size_t byte_offset_;
};

/**
 * @brief Contract for an append-only output that supports checked overwrites.
 *
 * @details A sink starts empty. `write()` appends all bytes or throws;
 * `write_at()` replaces bytes within the already-written prefix without
 * changing the append position; and `finish()` completes the output. A sink
 * may become unusable after `finish()`.
 */
template <class Sink>
concept SeekableBinaryOutputSink = requires(Sink& sink,
                                            std::size_t position,
                                            std::span<const std::byte> bytes) {
  { sink.write(bytes) } -> std::same_as<void>;
  { sink.write_at(position, bytes) } -> std::same_as<void>;
  { sink.finish() } -> std::same_as<void>;
};

/**
 * @brief Explicit owning sink for binary data collected in memory.
 *
 * @details This sink intentionally grows with the complete output. Prefer a
 * fixed span or file sink when bounded additional memory is required.
 */
class VectorOutputSink {
 public:
  /** @brief Return whether no output bytes have been written. */
  bool empty() const noexcept { return data_.empty(); }

  /** @brief Return the number of output bytes written. */
  std::size_t size_bytes() const noexcept { return data_.size(); }

  /** @brief Return a view of all output bytes. */
  std::span<const std::byte> bytes() const noexcept { return data_; }

  /** @brief Reserve memory for at least @p size_bytes output bytes. */
  void reserve_bytes(std::size_t size_bytes) { data_.reserve(size_bytes); }

  /** @brief Append @p bytes to the output. */
  void write(std::span<const std::byte> bytes) {
    if (bytes.size() > data_.max_size() - data_.size()) {
      throw std::length_error("Binary output is too large");
    }
    data_.insert(data_.end(), bytes.begin(), bytes.end());
  }

  /**
   * @brief Replace output bytes beginning at @p position.
   * @throws std::out_of_range if the range is outside the existing output.
   */
  void write_at(std::size_t position, std::span<const std::byte> bytes) {
    if (position > data_.size() || bytes.size() > data_.size() - position) {
      throw std::out_of_range("Binary output patch is outside the output");
    }
    std::ranges::copy(bytes, data_.begin() + position);
  }

  /** @brief Complete this in-memory output. */
  void finish() noexcept {}

  /** @brief Transfer the bytes and reset this sink to an empty state. */
  std::vector<std::byte> take() { return std::exchange(data_, {}); }

 private:
  std::vector<std::byte> data_;
};

static_assert(SeekableBinaryOutputSink<VectorOutputSink>);

/**
 * @brief Non-owning fixed-capacity output sink.
 *
 * @details The caller keeps the writable span alive and stable until the
 * associated writer is finished.
 */
class SpanOutputSink {
 public:
  /** @brief Construct an empty output over caller-owned storage. */
  explicit SpanOutputSink(std::span<std::byte> storage) : storage_(storage) {}

  /** @brief Return whether no output bytes have been written. */
  bool empty() const noexcept { return size_bytes_ == 0; }

  /** @brief Return the number of output bytes written. */
  std::size_t size_bytes() const noexcept { return size_bytes_; }

  /** @brief Return the fixed output capacity in bytes. */
  std::size_t capacity_bytes() const noexcept { return storage_.size(); }

  /** @brief Return the prefix written so far. */
  std::span<const std::byte> bytes() const noexcept {
    return storage_.first(size_bytes_);
  }

  /**
   * @brief Append @p bytes to the output.
   * @throws std::length_error if the fixed output capacity is insufficient.
   */
  void write(std::span<const std::byte> bytes) {
    if (bytes.size() > storage_.size() - size_bytes_) {
      throw std::length_error("Fixed binary output capacity exceeded");
    }
    std::ranges::copy(bytes, storage_.begin() + size_bytes_);
    size_bytes_ += bytes.size();
  }

  /**
   * @brief Replace output bytes beginning at @p position.
   * @throws std::out_of_range if the range is outside the written prefix.
   */
  void write_at(std::size_t position, std::span<const std::byte> bytes) {
    if (position > size_bytes_ || bytes.size() > size_bytes_ - position) {
      throw std::out_of_range("Binary output patch is outside the output");
    }
    std::ranges::copy(bytes, storage_.begin() + position);
  }

  /** @brief Complete this fixed-span output. */
  void finish() noexcept {}

 private:
  std::span<std::byte> storage_;
  std::size_t size_bytes_ = 0;
};

static_assert(SeekableBinaryOutputSink<SpanOutputSink>);

/**
 * @brief Bounded-buffer writer for canonical little-endian binary data.
 *
 * @details The writer does not own its sink. Small fields are accumulated in
 * a fixed-size staging buffer; large caller-owned byte spans may be sent
 * directly to the sink. Integer methods encode a fixed number of bytes,
 * independent of the host ABI. Call `finish()` to report output errors and
 * deliver the final buffered bytes before inspecting or consuming the sink.
 * Destruction does not perform I/O.
 */
class BinaryWriter {
 public:
  /** @brief Default bounded staging-buffer size. */
  static constexpr std::size_t kDefaultBufferBytes = 64 * 1024;

  /**
   * @brief Construct a writer with an owned, fixed-size staging buffer.
   *
   * @param sink Empty sink that remains alive and unmoved until `finish()`.
   * @param buffer_bytes Maximum staging-buffer memory owned by the writer.
   * @throws std::invalid_argument if @p buffer_bytes is zero.
   */
  template <SeekableBinaryOutputSink Sink>
  explicit BinaryWriter(Sink& sink,
                        std::size_t buffer_bytes = kDefaultBufferBytes)
      : owned_buffer_(allocate_buffer(buffer_bytes)),
        buffer_(owned_buffer_.get(), buffer_bytes) {
    bind(sink);
  }

  /**
   * @brief Construct a writer over a caller-owned staging buffer.
   *
   * @param sink Empty sink that remains alive and unmoved until `finish()`.
   * @param buffer Writable staging memory that remains alive and unmoved until
   * `finish()`.
   * @throws std::invalid_argument if @p buffer is empty.
   */
  template <SeekableBinaryOutputSink Sink>
  BinaryWriter(Sink& sink, std::span<std::byte> buffer) : buffer_(buffer) {
    if (buffer.empty()) {
      throw std::invalid_argument(
          "BinaryWriter staging buffer must be non-empty");
    }
    bind(sink);
  }

  BinaryWriter(const BinaryWriter&) = delete;
  BinaryWriter& operator=(const BinaryWriter&) = delete;
  BinaryWriter(BinaryWriter&&) = delete;
  BinaryWriter& operator=(BinaryWriter&&) = delete;

  /** @brief Return whether no bytes have been written. */
  bool empty() const noexcept { return size_bytes() == 0; }

  /** @brief Return the logical number of bytes written. */
  std::size_t size_bytes() const noexcept {
    return flushed_bytes_ + buffered_bytes_;
  }

  /** @brief Return the maximum staging-buffer size. */
  std::size_t buffer_size_bytes() const noexcept { return buffer_.size(); }

  /** @brief Write an unsigned eight-bit integer. */
  void write_u8(std::uint8_t value) { write_unsigned(value); }

  /** @brief Write an unsigned 16-bit integer in little-endian order. */
  void write_u16(std::uint16_t value) { write_unsigned(value); }

  /** @brief Write an unsigned 32-bit integer in little-endian order. */
  void write_u32(std::uint32_t value) { write_unsigned(value); }

  /** @brief Write an unsigned 64-bit integer in little-endian order. */
  void write_u64(std::uint64_t value) { write_unsigned(value); }

  /** @brief Write a signed eight-bit integer. */
  void write_i8(std::int8_t value) {
    write_unsigned(std::bit_cast<std::uint8_t>(value));
  }

  /** @brief Write a signed 16-bit integer in little-endian order. */
  void write_i16(std::int16_t value) {
    write_unsigned(std::bit_cast<std::uint16_t>(value));
  }

  /** @brief Write a signed 32-bit integer in little-endian order. */
  void write_i32(std::int32_t value) {
    write_unsigned(std::bit_cast<std::uint32_t>(value));
  }

  /** @brief Write a signed 64-bit integer in little-endian order. */
  void write_i64(std::int64_t value) {
    write_unsigned(std::bit_cast<std::uint64_t>(value));
  }

  /**
   * @brief Write a platform size as an unsigned 64-bit integer.
   * @throws std::length_error if `size_t` is wider than the wire type.
   */
  void write_size(std::size_t value) {
    if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
      if (value > std::numeric_limits<std::uint64_t>::max()) {
        throw std::length_error("Size does not fit in the wire format");
      }
    }
    write_u64(static_cast<std::uint64_t>(value));
  }

  /** @brief Append @p bytes without interpretation. */
  void write_bytes(std::span<const std::byte> bytes) {
    require_open();
    require_output_size(bytes.size());
    while (!bytes.empty()) {
      if (buffered_bytes_ != 0) {
        const std::size_t count =
            std::min(bytes.size(), buffer_.size() - buffered_bytes_);
        std::ranges::copy(bytes.first(count),
                          buffer_.begin() + buffered_bytes_);
        buffered_bytes_ += count;
        bytes = bytes.subspan(count);
        if (buffered_bytes_ == buffer_.size()) {
          flush_buffer();
        }
        continue;
      }

      if (bytes.size() >= buffer_.size()) {
        const std::size_t count = bytes.size();
        write_to_sink(bytes);
        flushed_bytes_ += count;
        return;
      }

      std::ranges::copy(bytes, buffer_.begin());
      buffered_bytes_ = bytes.size();
      return;
    }
  }

  /** @brief Append @p count zero bytes. */
  void write_zeros(std::size_t count) {
    require_open();
    require_output_size(count);
    while (count != 0) {
      const std::size_t writable = buffer_.size() - buffered_bytes_;
      const std::size_t chunk = std::min(count, writable);
      std::fill_n(buffer_.begin() + buffered_bytes_, chunk, std::byte{0});
      buffered_bytes_ += chunk;
      count -= chunk;
      if (buffered_bytes_ == buffer_.size()) {
        flush_buffer();
      }
    }
  }

  /**
   * @brief Pad with zero bytes to the next multiple of @p alignment.
   * @throws std::invalid_argument if @p alignment is zero.
   */
  void align_to(std::size_t alignment) {
    require_open();
    if (alignment == 0) {
      throw std::invalid_argument("Serialization alignment must be non-zero");
    }
    const std::size_t remainder = size_bytes() % alignment;
    if (remainder != 0) {
      write_zeros(alignment - remainder);
    }
  }

  /**
   * @brief Write a zero 64-bit field and return its byte position.
   *
   * @details Use `patch_u64()` to fill framed lengths after writing a payload.
   */
  std::size_t write_u64_placeholder() {
    const std::size_t position = size_bytes();
    write_u64(0);
    return position;
  }

  /**
   * @brief Replace an existing 64-bit field with @p value.
   * @throws std::out_of_range if the field is outside the current output.
   */
  void patch_u64(std::size_t position, std::uint64_t value) {
    require_open();
    if (position > size_bytes() || size_bytes() - position < sizeof(value)) {
      throw std::out_of_range("Serialization patch is outside the output");
    }

    std::array<std::byte, sizeof(value)> encoded{};
    for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
      encoded[byte] = static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
    }

    if (position >= flushed_bytes_) {
      const std::size_t buffered_position = position - flushed_bytes_;
      std::ranges::copy(encoded, buffer_.begin() + buffered_position);
      return;
    }

    if (position + encoded.size() > flushed_bytes_) {
      flush_buffer();
    }
    write_at_sink(position, encoded);
  }

  /**
   * @brief Deliver currently buffered bytes while keeping the writer open.
   * @throws Any exception reported by the sink. A failed writer is unusable.
   */
  void flush() {
    require_open();
    flush_buffer();
  }

  /**
   * @brief Deliver all bytes and complete the sink.
   *
   * @details This operation is idempotent after success. A sink failure is
   * propagated and leaves the writer unusable; already-written output is not
   * rolled back.
   */
  void finish() {
    if (state_ == State::kFinished) {
      return;
    }
    require_open();
    flush_buffer();
    try {
      finish_sink_(sink_);
      state_ = State::kFinished;
    } catch (...) {
      state_ = State::kFailed;
      throw;
    }
  }

 private:
  enum class State { kOpen, kFinished, kFailed };

  using WriteSink = void (*)(void*, std::span<const std::byte>);
  using WriteAtSink = void (*)(void*, std::size_t, std::span<const std::byte>);
  using FinishSink = void (*)(void*);

  static std::unique_ptr<std::byte[]> allocate_buffer(
      std::size_t buffer_bytes) {
    if (buffer_bytes == 0) {
      throw std::invalid_argument(
          "BinaryWriter staging buffer must be non-empty");
    }
    return std::make_unique<std::byte[]>(buffer_bytes);
  }

  template <SeekableBinaryOutputSink Sink>
  void bind(Sink& sink) noexcept {
    sink_ = std::addressof(sink);
    write_sink_ = [](void* output, std::span<const std::byte> bytes) {
      static_cast<Sink*>(output)->write(bytes);
    };
    write_at_sink_ = [](void* output, std::size_t position,
                        std::span<const std::byte> bytes) {
      static_cast<Sink*>(output)->write_at(position, bytes);
    };
    finish_sink_ = [](void* output) { static_cast<Sink*>(output)->finish(); };
  }

  void require_open() const {
    if (state_ == State::kFinished) {
      throw std::logic_error("BinaryWriter is already finished");
    }
    if (state_ == State::kFailed) {
      throw std::logic_error("BinaryWriter is unusable after sink failure");
    }
  }

  void require_output_size(std::size_t count) const {
    if (count > std::numeric_limits<std::size_t>::max() - size_bytes()) {
      throw std::length_error("Serialized output is too large");
    }
  }

  void flush_buffer() {
    if (buffered_bytes_ == 0) {
      return;
    }
    const std::size_t count = buffered_bytes_;
    write_to_sink(buffer_.first(count));
    flushed_bytes_ += count;
    buffered_bytes_ = 0;
  }

  void write_to_sink(std::span<const std::byte> bytes) {
    try {
      write_sink_(sink_, bytes);
    } catch (...) {
      state_ = State::kFailed;
      throw;
    }
  }

  void write_at_sink(std::size_t position, std::span<const std::byte> bytes) {
    try {
      write_at_sink_(sink_, position, bytes);
    } catch (...) {
      state_ = State::kFailed;
      throw;
    }
  }

  template <std::unsigned_integral T>
  void write_unsigned(T value) {
    std::array<std::byte, sizeof(T)> encoded{};
    for (std::size_t byte = 0; byte < sizeof(T); ++byte) {
      encoded[byte] = static_cast<std::byte>((value >> (byte * 8)) & T{0xff});
    }
    write_bytes(encoded);
  }

  void* sink_ = nullptr;
  WriteSink write_sink_ = nullptr;
  WriteAtSink write_at_sink_ = nullptr;
  FinishSink finish_sink_ = nullptr;
  std::unique_ptr<std::byte[]> owned_buffer_;
  std::span<std::byte> buffer_;
  std::size_t flushed_bytes_ = 0;
  std::size_t buffered_bytes_ = 0;
  State state_ = State::kOpen;
};

/**
 * @brief Bounds-checked reader for canonical little-endian binary data.
 *
 * @details The reader does not own its input. Copies have independent cursors,
 * which allows a compound deserializer to parse into a temporary reader and
 * assign it back only after validation succeeds.
 */
class BinaryReader {
 public:
  /** @brief Construct a reader over caller-owned @p data. */
  explicit BinaryReader(std::span<const std::byte> data) : data_(data) {}

  /** @brief Return whether all input bytes have been consumed. */
  bool empty() const noexcept { return remaining() == 0; }

  /** @brief Return the total number of bytes in this reader. */
  std::size_t size_bytes() const noexcept { return data_.size(); }

  /** @brief Return the number of bytes consumed by this reader. */
  std::size_t position() const noexcept { return offset_; }

  /** @brief Return the current byte offset in the outermost input. */
  std::size_t byte_offset() const noexcept { return absolute_position(); }

  /** @brief Return the number of unconsumed bytes. */
  std::size_t remaining() const noexcept { return data_.size() - offset_; }

  /** @brief Return all currently unconsumed bytes. */
  std::span<const std::byte> remaining_bytes() const noexcept {
    return data_.subspan(offset_);
  }

  /** @brief Read an unsigned eight-bit integer. */
  std::uint8_t read_u8() { return read_unsigned<std::uint8_t>(); }

  /** @brief Read an unsigned little-endian 16-bit integer. */
  std::uint16_t read_u16() { return read_unsigned<std::uint16_t>(); }

  /** @brief Read an unsigned little-endian 32-bit integer. */
  std::uint32_t read_u32() { return read_unsigned<std::uint32_t>(); }

  /** @brief Read an unsigned little-endian 64-bit integer. */
  std::uint64_t read_u64() { return read_unsigned<std::uint64_t>(); }

  /** @brief Read a signed eight-bit integer. */
  std::int8_t read_i8() {
    return std::bit_cast<std::int8_t>(read_unsigned<std::uint8_t>());
  }

  /** @brief Read a signed little-endian 16-bit integer. */
  std::int16_t read_i16() {
    return std::bit_cast<std::int16_t>(read_unsigned<std::uint16_t>());
  }

  /** @brief Read a signed little-endian 32-bit integer. */
  std::int32_t read_i32() {
    return std::bit_cast<std::int32_t>(read_unsigned<std::uint32_t>());
  }

  /** @brief Read a signed little-endian 64-bit integer. */
  std::int64_t read_i64() {
    return std::bit_cast<std::int64_t>(read_unsigned<std::uint64_t>());
  }

  /**
   * @brief Read an unsigned 64-bit size and convert it to `size_t`.
   * @throws std::length_error if the value is not representable.
   */
  std::size_t read_size() {
    BinaryReader candidate = *this;
    const std::size_t field_offset = candidate.absolute_position();
    const std::uint64_t encoded = candidate.read_u64();
    if (encoded > std::numeric_limits<std::size_t>::max()) {
      throw std::length_error("Serialized size at byte " +
                              std::to_string(field_offset) +
                              " does not fit in size_t");
    }
    *this = candidate;
    return static_cast<std::size_t>(encoded);
  }

  /**
   * @brief Read exactly @p count uninterpreted bytes.
   * @throws SerializationError if the input is truncated.
   */
  std::span<const std::byte> read_bytes(std::size_t count) {
    require(count);
    const std::span<const std::byte> result = data_.subspan(offset_, count);
    offset_ += count;
    return result;
  }

  /**
   * @brief Read a bounded region as an independent child reader.
   * @throws SerializationError if the input is truncated.
   */
  BinaryReader read_subreader(std::size_t count) {
    const std::size_t origin = absolute_position();
    return BinaryReader(read_bytes(count), origin);
  }

  /**
   * @brief Advance by exactly @p count bytes.
   * @throws SerializationError if the input is truncated.
   */
  void skip(std::size_t count) {
    require(count);
    offset_ += count;
  }

  /**
   * @brief Consume at most @p maximum trailing zero-padding bytes.
   * @throws SerializationError for excessive or non-zero trailing bytes.
   */
  void require_zero_padding(std::size_t maximum) {
    if (remaining() > maximum) {
      throw SerializationError("Unexpected serialized payload bytes",
                               absolute_position());
    }
    BinaryReader candidate = *this;
    while (!candidate.empty()) {
      const std::size_t byte_offset = candidate.absolute_position();
      if (candidate.read_u8() != 0) {
        throw SerializationError("Non-zero serialized payload padding",
                                 byte_offset);
      }
    }
    *this = candidate;
  }

 private:
  BinaryReader(std::span<const std::byte> data, std::size_t origin)
      : data_(data), origin_(origin) {}

  std::size_t absolute_position() const noexcept { return origin_ + offset_; }

  void require(std::size_t count) const {
    if (count > remaining()) {
      throw SerializationError("Truncated serialized input",
                               absolute_position());
    }
  }

  template <std::unsigned_integral T>
  T read_unsigned() {
    require(sizeof(T));
    T value = 0;
    for (std::size_t byte = 0; byte < sizeof(T); ++byte) {
      value |=
          static_cast<T>(std::to_integer<std::uint8_t>(data_[offset_ + byte]))
          << (byte * 8);
    }
    offset_ += sizeof(T);
    return value;
  }

  std::span<const std::byte> data_;
  std::size_t offset_ = 0;
  std::size_t origin_ = 0;
};

}  // namespace pixie
