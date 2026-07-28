#pragma once

/**
 * @file huffman.h
 * @brief Common CRTP contract for Huffman entropy codecs.
 *
 * PivCo-Huffman ("Pivot-Coded Huffman") reuses the wavelet-tree "tree of
 * bitmaps" layout to turn sequential, bit-by-bit Huffman-tree traversals into
 * vectorizable operations. A codec encodes a byte sequence into a compressed
 * stream by building a Huffman-shaped tree of per-node bitmaps, and decodes it
 * back by traversing that tree. Concrete implementations live under
 * `<pixie/huffman/implementations.h>`.
 *
 * @see Marcin Zukowski, "PivCo-Huffman", v1.0 (2026).
 */

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pixie {

/**
 * @brief CRTP facade for Huffman entropy codecs.
 *
 * The contract is the source of truth for observable codec semantics. Each
 * concrete implementation inherits `HuffmanBase<Impl>` and supplies the
 * required `*_impl()` extension points. There is no virtual dispatch: the
 * facade delegates statically through CRTP, mirroring the other Pixie families.
 *
 * Range and ownership conventions:
 *  - Symbol sequences are zero-based byte streams (`symbol_type`).
 *  - Compressed streams are byte-oriented views (`std::byte`).
 *  - The caller keeps any non-owning view alive for the codec lifetime.
 *
 * @tparam Impl Concrete codec type implementing the `*_impl()` contract.
 *
 * @see `<pixie/huffman/implementations.h>` for the available concrete
 *      implementations.
 */
template <class Impl>
class HuffmanBase {
 public:
  /** @brief Symbol type handled by the codec: one byte per symbol. */
  using symbol_type = std::uint8_t;

  /**
   * @brief Number of symbols in the original (uncompressed) stream.
   * @return Logical uncompressed symbol count.
   */
  std::size_t uncompressed_size() const {
    return impl().uncompressed_size_impl();
  }

  /**
   * @brief Number of bytes in the compressed representation.
   * @return Compressed stream size in bytes.
   */
  std::size_t compressed_size() const { return impl().compressed_size_impl(); }

  /**
   * @brief Check whether the codec holds no data.
   * @return `true` when `uncompressed_size() == 0`.
   */
  bool empty() const { return uncompressed_size() == 0; }

  /**
   * @brief Read-only view of the compressed byte stream.
   * @return Span over the serialized representation.
   *
   * @note The returned view is invalidated by codec destruction.
   */
  std::span<const std::byte> compressed_data() const {
    return impl().compressed_data_impl();
  }

  /**
   * @brief Reconstruct the original symbol sequence.
   * @return Decoded symbols of length `uncompressed_size()`.
   */
  std::vector<symbol_type> decode() const { return impl().decode_impl(); }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
