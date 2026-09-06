#pragma once

/**
 * @file pivco_huffman.h
 * @brief Pixie adapter for the PivCo-Huffman file codec.
 *
 * The codec implementation is Apache-2.0-licensed and lives in
 * `third_party/pivco`. Its compiled C implementation preserves the PivCo wire
 * format, block layout, and architecture-specific kernels while this header
 * exposes the standard Pixie Huffman CRTP facade.
 */

#include <pivcohuf_file.h>
#include <pixie/huffman.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace pixie {

/**
 * @brief Owning PivCo-Huffman codec.
 *
 * Encoding uses 32 KiB file blocks, plain PivCo-Huffman bitmaps (FSE disabled),
 * the hybrid vertical flat layout, and the SIMD backend selected at build
 * time. Construction from compressed bytes copies the serialized stream, so
 * returned views remain valid for this object's lifetime.
 *
 * @note Consumers linking `pixie::pixie` receive the required compiled backend
 *       transitively.
 */
class PivCoHuffman : public HuffmanBase<PivCoHuffman> {
 public:
  /** @brief Symbol type handled by the codec: one byte per symbol. */
  using symbol_type = std::uint8_t;

  /** @brief Number of input symbols encoded in each full block. */
  static constexpr std::size_t kBlockSize = 32 * 1024;

  /**
   * @brief Compress @p input into an owned PivCo stream.
   * @param input Byte sequence to copy into the compressed representation; an
   * empty span creates an empty codec.
   * @throws std::runtime_error if the backend cannot encode the input.
   */
  explicit PivCoHuffman(std::span<const symbol_type> input) { encode(input); }

  /**
   * @brief Load an owning copy of a serialized PivCo stream.
   * @param compressed Serialized bytes to copy; an empty span represents an
   * empty codec.
   * @throws std::runtime_error if a non-empty stream has an invalid header.
   */
  explicit PivCoHuffman(std::span<const std::byte> compressed)
      : compressed_(compressed.begin(), compressed.end()) {
    if (compressed_.empty()) {
      uncompressed_size_ = 0;
      return;
    }
    const int status = pivcohuf_peek_uncompressed_size(
        bytes(compressed_.data()), compressed_.size(), &uncompressed_size_);
    check(status, "inspect compressed stream");
  }

  /** @brief Return the logical uncompressed symbol count. */
  std::size_t uncompressed_size_impl() const { return uncompressed_size_; }

  /** @brief Return the serialized compressed size in bytes. */
  std::size_t compressed_size_impl() const { return compressed_.size(); }

  /** @brief Return a view over the owned PivCo stream. */
  std::span<const std::byte> compressed_data_impl() const {
    return compressed_;
  }

  /**
   * @brief Decode and return the complete input symbol sequence.
   * @throws std::runtime_error if the serialized stream is malformed or its
   * decoded size does not match its header.
   */
  std::vector<symbol_type> decode_impl() const {
    if (uncompressed_size_ == 0) {
      return {};
    }
    std::vector<symbol_type> output(uncompressed_size_);
    std::size_t output_size = output.size();
    const int status =
        pivcohuf_decompress(bytes(compressed_.data()), compressed_.size(),
                            output.data(), &output_size);
    check(status, "decode compressed stream");
    if (output_size != uncompressed_size_) {
      throw std::runtime_error(
          "PivCo decoder returned an unexpected output size");
    }
    return output;
  }

 private:
  static std::uint8_t* bytes(std::byte* ptr) {
    return reinterpret_cast<std::uint8_t*>(ptr);
  }

  static const std::uint8_t* bytes(const std::byte* ptr) {
    return reinterpret_cast<const std::uint8_t*>(ptr);
  }

  static void check(int status, const char* operation) {
    if (status != PIVCOHUF_OK) {
      throw std::runtime_error(std::string("failed to ") + operation +
                               " with PivCo status " + std::to_string(status));
    }
  }

  void encode(std::span<const symbol_type> input) {
    uncompressed_size_ = input.size();
    const std::size_t bound =
        pivcohuf_compress_bound_blk(input.size(), kBlockSize);
    compressed_.resize(bound);
    std::size_t compressed_size = compressed_.size();
    const int status = pivcohuf_compress_blk(
        input.data(), input.size(), bytes(compressed_.data()), &compressed_size,
        0, kBlockSize, nullptr);
    check(status, "encode input");
    compressed_.resize(compressed_size);
  }

  std::size_t uncompressed_size_ = 0;
  std::vector<std::byte> compressed_;
};

}  // namespace pixie
