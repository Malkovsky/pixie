#pragma once

#include <pixie/bits.h>
#include <pixie/storage/aligned.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

/**
 * @brief Standalone, fixed-capacity, circular packed-bit block.
 * @details Owns typed words without allocation. The ring wraps at size(), never
 * at capacity. Metadata and CacheLine alignment padding are additional to the
 * payload budget. Copies own independent contents. Local mutations require
 * valid arguments and never allocate or throw; checked entry points are
 * supplied by SequenceTree. Unused bits are not logical elements.
 * @tparam PayloadBits Payload capacity, a positive multiple of 64.
 */
template <std::size_t PayloadBits = 2048>
class alignas(CacheLine) BitBlock {
  static_assert(PayloadBits >= 64 && PayloadBits % 64 == 0);
  static_assert(PayloadBits <= std::numeric_limits<std::size_t>::max() / 2);

 public:
  /** @brief Indexed reads return bits by value. */
  using value_type = bool;
  /** @brief Maximum valid element count. */
  static constexpr std::size_t capacity = PayloadBits;
  /** @brief Typed, LSB-first payload words. */
  using Payload = std::array<std::uint64_t, PayloadBits / 64>;
  /** @brief Byte offset of the payload; allocation alignment is alignof(Block).
   */
  static constexpr std::size_t payload_offset_bytes = 2 * sizeof(std::uint64_t);
  /** @brief Payload requires word alignment, not separate CacheLine alignment.
   */
  static constexpr std::size_t payload_alignment = alignof(std::uint64_t);

  /** @brief Construct an empty block without allocation. */
  BitBlock() noexcept {}
  /** @brief Copy the live physical words, preserving the circular origin.
   * @details Unused capacity is not read or initialized. Copies own their data.
   */
  BitBlock(const BitBlock& other) noexcept : n_(other.n_), h_(other.h_) {
    std::copy_n(other.bits_.begin(), n_ / 64 + (n_ % 64 != 0), bits_.begin());
  }
  /** @brief Replace contents with an independent copy of the live words.
   * @param other Source block; self-assignment is a no-op. @return *this.
   */
  BitBlock& operator=(const BitBlock& other) noexcept {
    if (this != &other) {
      n_ = other.n_;
      h_ = other.h_;
      std::copy_n(other.bits_.begin(), n_ / 64 + (n_ % 64 != 0), bits_.begin());
    }
    return *this;
  }
  /**
   * @brief Copy n LSB-first bits, ignoring extra words and dirty high padding.
   * @details The source need not outlive this block; initial origin is zero.
   * @param words Source word span.
   * @param n Logical count in [0, capacity].
   * @throws std::invalid_argument If n exceeds capacity or words is too short.
   */
  BitBlock(std::span<const std::uint64_t> words, std::size_t n) {
    if (n > capacity || n / 64 + (n % 64 != 0) > words.size()) {
      throw std::invalid_argument("BitBlock: input size");
    }
    n_ = n;
    std::copy_n(words.begin(), n / 64 + (n % 64 != 0), bits_.begin());
  }
  /** @brief Return the valid bit count. @return Count excluding padding. */
  std::size_t size() const noexcept { return n_; }
  /** @brief Test emptiness. @return Whether size() is zero. */
  bool empty() const noexcept { return n_ == 0; }
  /**
   * @brief Read a bit by value.
   * @param i Zero-based index in [0, size()).
   * @return Logical bit, not a writable reference.
   * @throws std::out_of_range If i >= size().
   */
  bool operator[](std::size_t i) const {
    if (i >= n_) {
      throw std::out_of_range("BitBlock: index");
    }
    const auto sum = h_ + i;
    const auto p = sum >= n_ ? sum - n_ : sum;
    return (bits_[p / 64] >> (p % 64)) & 1;
  }
  /**
   * @brief Read at most 64 consecutive logical bits without flattening.
   * @details Requires offset <= size() and width <= min(64, size()-offset).
   * A zero width returns zero, including on an empty block. The result is
   * LSB-first with zero high padding. Reads only words containing logical bits,
   * even across the circular boundary; width 64 never shifts by 64. Invalid
   * arguments violate the primitive precondition (asserted in debug builds).
   * @param offset First logical bit.
   * @param width Number of bits to decode.
   * @return Unsigned field value.
   */
  std::uint64_t read_bits(std::size_t offset,
                          std::size_t width) const noexcept {
    assert(offset <= n_ && width <= 64 && width <= n_ - offset);
    if (width == 0) {
      return 0;
    }
    // Both operands are below n_, and capacity <= SIZE_MAX / 2. The
    // circular sum cannot overflow and needs at most one subtraction.
    const auto sum = h_ + offset;
    const auto physical = sum >= n_ ? sum - n_ : sum;
    const auto first = std::min(width, n_ - physical);
    auto read = [&](std::size_t p, std::size_t count) {
      const auto shift = p % 64;
      auto value = bits_[p / 64] >> shift;
      if (count > 64 - shift) {
        value |= bits_[p / 64 + 1] << (64 - shift);
      }
      return count == 64 ? value : value & ((std::uint64_t{1} << count) - 1);
    };
    const auto value = read(physical, first);
    return first == width ? value : value | (read(0, width - first) << first);
  }
  /**
   * @brief Read a whole word from a word-aligned logical ring.
   * @details Requires size(), circular origin, and offset to be multiples of
   * 64, with offset < size(). Invalid inputs violate the primitive
   * precondition.
   * @param offset Zero-based logical bit offset.
   * @return The 64-bit word at offset, without bit-field decoding.
   */
  std::uint64_t read_aligned_word(std::size_t offset) const noexcept {
    assert(n_ % 64 == 0 && h_ % 64 == 0 && offset % 64 == 0 && offset < n_);
    // The common normalized case has no metadata-dependent value address.
    if (h_ == 0) {
      return bits_[offset / 64];
    }
    const auto sum = h_ + offset;
    const auto physical = sum >= n_ ? sum - n_ : sum;
    return bits_[physical / 64];
  }
  /**
   * @brief Replace a word in a word-aligned logical ring.
   * @pre size(), origin and offset are multiples of 64; offset < size().
   * @param offset Logical bit offset. @param value Replacement word.
   */
  void write_aligned_word(std::size_t offset, std::uint64_t value) noexcept {
    assert(n_ % 64 == 0 && h_ % 64 == 0 && offset % 64 == 0 && offset < n_);
    const auto sum = h_ + offset;
    bits_[(sum >= n_ ? sum - n_ : sum) / 64] = value;
  }
  /**
   * @brief Replace up to 64 logical bits, preserving all other bits.
   * @pre offset <= size() and width <= min(64, size()-offset).
   * @param offset First logical bit. @param width Number of bits to replace.
   * @param value LSB-first replacement; high bits outside width are ignored.
   */
  void write_bits(std::size_t offset,
                  std::size_t width,
                  std::uint64_t value) noexcept {
    assert(offset <= n_ && width <= 64 && width <= n_ - offset);
    if (width == 0) {
      return;
    }
    const auto sum = h_ + offset;
    const auto p = sum >= n_ ? sum - n_ : sum;
    const auto first = std::min(width, n_ - p);
    copy_packed_bits(&value, 0, bits_.data(), p, first);
    copy_packed_bits(&value, first, bits_.data(), 0, width - first);
  }
  /**
   * @brief Materialize logical words without changing representation.
   * @return LSB-first words with all padding zero.
   */
  Payload flatten() const noexcept {
    Payload out{};
    copy(0, out, 0, n_);
    return out;
  }
  /**
   * @brief Rotate a valid half-open local interval without allocation.
   * @details Requires left <= right <= size(). Empty ranges are no-ops; a whole
   * block rotation adjusts only the origin. Invalid arguments violate the
   * primitive precondition (asserted in debug builds).
   * @param left Inclusive start.
   * @param right Exclusive end.
   * @param distance Left distance, reduced modulo the nonempty range length.
   */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= n_);
    const auto length = right - left;
    if (length == 0 || (distance %= length) == 0) {
      return;
    }
    if (length == n_) {
      h_ = (h_ + distance) % n_;
      return;
    }
    Payload scratch{};
    copy(left + distance, scratch, 0, length - distance);
    copy(left, scratch, length - distance, distance);
    const auto p = (h_ + left) % n_;
    const auto prefix = std::min(length, n_ - p);
    copy_packed_bits(scratch.data(), 0, bits_.data(), p, prefix);
    copy_packed_bits(scratch.data(), prefix, bits_.data(), 0, length - prefix);
  }
  /**
   * @brief Repartition two blocks, preserving their concatenated logical bits.
   * @details Requires distinct blocks, left_size <= capacity, and
   * left_size <= size()+rhs.size() <= left_size+capacity. Nonallocating and
   * nonthrowing. Both circular origins may normalize. An empty rhs permits
   * splitting; left_size equal to the total permits combining.
   * @param rhs Right block; receives the remaining suffix.
   * @param left_size Desired valid count in this block.
   */
  void redistribute(BitBlock& rhs, std::size_t left_size) noexcept {
    assert(this != &rhs);
    const auto total = n_ + rhs.n_;
    assert(left_size <= capacity && left_size <= total &&
           total - left_size <= capacity);
    if (left_size == total && h_ == 0) {
      // A normalized prefix stays in place when absorbing the whole donor.
      // Copy only the incoming bits, including a possibly circular donor.
      // Partial-word writes preserve destination padding, so initialize newly
      // exposed words before the copy kernel can read those padding bits.
      std::fill(bits_.begin() + n_ / 64 + (n_ % 64 != 0),
                bits_.begin() + total / 64 + (total % 64 != 0), 0);
      rhs.copy(0, bits_, n_, rhs.n_);
      n_ = total;
      rhs.n_ = rhs.h_ = 0;
      return;
    }
    Payload a{}, b{};
    auto gather = [&](std::size_t start, Payload& out, std::size_t count) {
      const auto first = start < n_ ? std::min(count, n_ - start) : 0;
      copy(start, out, 0, first);
      rhs.copy(start + first - std::min(start + first, n_), out, first,
               count - first);
    };
    gather(0, a, left_size);
    gather(left_size, b, total - left_size);
    bits_ = a;
    rhs.bits_ = b;
    n_ = left_size;
    rhs.n_ = total - left_size;
    h_ = rhs.h_ = 0;
  }
  /** @brief Materialize logical order at physical origin zero, without
   * allocation.
   * @details Preserves size and contents; copies at most capacity bits.
   */
  void normalize() noexcept {
    if (h_ == 0) {
      return;
    }
    Payload scratch{};
    copy(0, scratch, 0, n_);
    bits_ = scratch;
    h_ = 0;
  }
  /** @brief Read one word from a normalized block without metadata loads.
   * @pre Origin is zero; offset is a multiple of 64 and offset+64 <= size().
   * @param offset Zero-based bit offset. @return The 64-bit word.
   */
  std::uint64_t read_normalized_word(std::size_t offset) const noexcept {
    assert(h_ == 0 && offset % 64 == 0 && offset <= n_ && n_ - offset >= 64);
    return bits_[offset / 64];
  }
  /** @brief Report fixed payload storage. @return Bytes including unused bits.
   */
  std::size_t payload_capacity_bytes() const noexcept { return sizeof(bits_); }
  /** @brief Report nonpayload storage. @return Metadata and alignment padding.
   */
  std::size_t metadata_bytes() const noexcept {
    return sizeof(*this) - sizeof(bits_);
  }

 private:
  void copy(std::size_t offset,
            Payload& out,
            std::size_t destination,
            std::size_t count) const noexcept {
    if (count == 0) {
      return;
    }
    const auto p = (h_ + offset) % n_;
    const auto prefix = std::min(count, n_ - p);
    copy_packed_bits(bits_.data(), p, out.data(), destination, prefix);
    copy_packed_bits(bits_.data(), 0, out.data(), destination + prefix,
                     count - prefix);
  }
  std::uint64_t n_ = 0;
  std::uint64_t h_ = 0;
  Payload bits_;
};

}  // namespace pixie::detail::sequence
/// @endcond
