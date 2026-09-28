#pragma once

#include <pixie/detail/sequence/packed_bit_block.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

/// @cond PIXIE_SEQUENCE_INTERNAL
namespace pixie::detail::sequence {

/**
 * @brief Fixed-width unsigned values over an owning circular packed bit block.
 * @details No adapter metadata or allocation: the underlying block records the
 * bit count and origin. Logical boundaries are multiples of Width. Copies own
 * independent data; reads return values. Local mutations require valid inputs
 * and never allocate or throw. Narrow Width is an unsupported internal
 * control, not a codec: construction rejects values which do not fit.
 * @tparam T Bool or an unsigned integer with at most 64 value bits.
 * @tparam StorageBits Complete local block budget, not a tree leaf budget.
 * @tparam Width Positive field width, at most the value width of T.
 * @tparam BitStorage Owning bit block providing the bounded payload storage.
 */
template <class T,
          std::size_t StorageBits = 2048,
          std::size_t Width = std::numeric_limits<T>::digits,
          class BitStorage = PackedBitBlock<StorageBits>>
class PackedValueBlock {
  static_assert(std::is_integral_v<T> && std::is_unsigned_v<T> &&
                std::numeric_limits<T>::digits <= 64);
  static_assert(Width > 0 && Width <= std::numeric_limits<T>::digits);
  using Bits = BitStorage;

 public:
  /** @brief Immutable indexed result. @details Includes bool by value. */
  using value_type = T;
  /** @brief Maximum element count. @details Unused trailing bits are padding.
   */
  static constexpr std::size_t capacity = Bits::capacity / Width;
  /** @brief Stored bits per value. @details Independent of positional counts.
   */
  static constexpr std::size_t width = Width;
  /** @brief Payload byte offset. @details Inherited from the sole member. */
  static constexpr std::size_t payload_offset_bytes =
      Bits::payload_offset_bytes;
  /** @brief Payload alignment. @details The block itself is cache-line aligned.
   */
  static constexpr std::size_t payload_alignment = Bits::payload_alignment;
  static_assert(capacity > 0);

  /** @brief Construct empty. @details Performs no allocation. */
  PackedValueBlock() noexcept {}
  /**
   * @brief Copy values into independent bounded storage.
   * @details The source is needed only during construction; initial origin is
   * zero.
   * @param values At most capacity representable fields.
   * @throws std::invalid_argument If count or a value exceeds its capacity.
   */
  explicit PackedValueBlock(std::span<const T> values) {
    if (values.size() > capacity) {
      throw std::invalid_argument("PackedValueBlock: input size");
    }
    typename Bits::Payload words;
    if constexpr (Width != 64) {
      words.fill(0);
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (static_cast<std::uint64_t>(values[i]) > field_max) {
        throw std::invalid_argument("PackedValueBlock: value width");
      }
      if constexpr (Width == 64) {
        words[i] = values[i];
      } else {
        encode(words, i, values[i]);
      }
    }
    bits_ = Bits(words, values.size() * Width);
  }
  /** @brief Return element count. @details Excludes unused trailing bits.
   * @return Number of logical values.
   */
  std::size_t size() const noexcept { return bits_.size() / Width; }
  /** @brief Test emptiness. @details Equivalent to size() == 0.
   * @return Whether the block contains no values.
   */
  bool empty() const noexcept { return bits_.empty(); }
  /**
   * @brief Read a logical value without flattening the payload.
   * @details Decodes at most 64 bits, including wrapped fields.
   * @param i Zero-based element position.
   * @return The decoded value at i.
   * @throws std::out_of_range If i >= size().
   */
  T operator[](std::size_t i) const {
    if (i >= size()) {
      throw std::out_of_range("PackedValueBlock: index");
    }
    return read_unchecked(i);
  }
  /** @brief Read a value whose position was validated by the owning tree.
   * @pre i < size(). @param i Zero-based position. @return The decoded value.
   * @details Full-width words stay normalized, allowing a direct payload read.
   */
  T read_unchecked(std::size_t i) const noexcept {
    assert(i < size());
    if constexpr (Width == 1) {
      return static_cast<T>(bits_[i]);
    } else if constexpr (Width == 64) {
      return static_cast<T>(bits_.read_normalized_word(i * 64));
    } else {
      return static_cast<T>(bits_.read_bits(i * Width, Width));
    }
  }
  /**
   * @brief Replace a logical field without changing its circular layout.
   * @pre i < size() and value is representable in Width bits.
   * @param i Zero-based position. @param value Replacement field.
   */
  void set_at(std::size_t i, T value) noexcept {
    assert(i < size() && static_cast<std::uint64_t>(value) <= field_max);
    if constexpr (Width == 64) {
      bits_.write_aligned_word(i * Width, value);
    } else {
      bits_.write_bits(i * Width, Width, value);
    }
  }
  /**
   * @brief Rotate a valid half-open element interval without allocation.
   * @details Requires left <= right <= size(). Reduces distance in element
   * units before converting to bits; empty intervals do nothing.
   * @param left Inclusive starting element position.
   * @param right Exclusive ending element position.
   * @param distance Leftward distance in elements, reduced modulo length.
   */
  void rotate_left(std::size_t left,
                   std::size_t right,
                   std::size_t distance) noexcept {
    assert(left <= right && right <= size());
    if constexpr (Width == 1) {
      bits_.rotate_left(left, right, distance);
    } else if (left != right) {
      bits_.rotate_left(left * Width, right * Width,
                        (distance % (right - left)) * Width);
      if constexpr (Width == 64) {
        bits_.normalize();
      }
    }
  }
  /**
   * @brief Insert one value at a local offset without allocation.
   * @details Requires offset <= size() < capacity. Shifts [offset,size()) one
   * position right using bounded stack scratch and places value at offset.
   * Does not allocate or throw. This primitive supports sequence-tree leaf
   * insertion, not a public arbitrary transformation on a container.
   * @param offset Local insertion position in [0,size()].
   * @param value Field to insert.
   */
  void insert_at(std::size_t offset, T value) noexcept {
    assert(offset <= size() && size() < capacity);
    typename Bits::Payload words{};
    const auto n = size();
    if constexpr (Width == 64) {
      for (std::size_t i = 0; i < n; ++i) {
        words[i] = bits_.read_aligned_word(i * 64);
      }
      std::memmove(words.data() + offset + 1, words.data() + offset,
                   (n - offset) * sizeof(std::uint64_t));
      words[offset] = static_cast<std::uint64_t>(value);
    } else {
      for (std::size_t i = 0; i < offset; ++i) {
        encode(words, i, bits_.read_bits(i * Width, Width));
      }
      encode(words, offset, value);
      for (std::size_t i = offset; i < n; ++i) {
        encode(words, i + 1, bits_.read_bits(i * Width, Width));
      }
    }
    bits_ = Bits(words, (n + 1) * Width);
  }
  /**
   * @brief Repartition concatenated values between distinct blocks.
   * @details Requires left_size <= capacity, left_size <= size()+rhs.size(),
   * and size()+rhs.size()-left_size <= capacity. Preserves order, normalizes
   * origins, and neither allocates nor throws.
   * @param rhs Distinct following block; receives the remaining suffix.
   * @param left_size Desired number of elements in this block afterward.
   */
  void redistribute(PackedValueBlock& rhs, std::size_t left_size) noexcept {
    assert(this != &rhs && left_size <= capacity &&
           left_size <= size() + rhs.size() &&
           size() + rhs.size() - left_size <= capacity);
    bits_.redistribute(rhs.bits_, left_size * Width);
  }
  /**
   * @brief Materialize a nonnegative index bias using bounded stack scratch.
   * @details Requires every decoded value plus bias to fit Width. Normalizes
   * the origin; does not allocate or throw. This primitive supports tagged
   * permutation leaves, not a public arbitrary transformation on a container.
   * @param bias Nonnegative value added to each logical field.
   */
  void add_bias(std::uint64_t bias) noexcept {
    typename Bits::Payload words{};
    const auto n = size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto value = bits_.read_bits(i * Width, Width);
      assert(bias <= field_max && value <= field_max - bias);
      encode(words, i, value + bias);
    }
    bits_ = Bits(words, n * Width);
  }
  /** @brief Report payload capacity. @details Includes field/word slack.
   * @return Bytes reserved for packed payload.
   */
  std::size_t payload_capacity_bytes() const noexcept {
    return bits_.payload_capacity_bytes();
  }
  /** @brief Report metadata bytes. @details Includes local alignment padding.
   * @return Non-payload bytes in the complete block object.
   */
  std::size_t metadata_bytes() const noexcept { return bits_.metadata_bytes(); }

 private:
  static constexpr std::uint64_t field_max =
      std::numeric_limits<std::uint64_t>::max() >> (64 - Width);
  static void encode(typename Bits::Payload& words,
                     std::size_t i,
                     std::uint64_t value) noexcept {
    const auto bit = i * Width;
    const auto shift = bit % 64;
    words[bit / 64] |= value << shift;
    if (Width > 64 - shift) {
      words[bit / 64 + 1] |= value >> (64 - shift);
    }
  }
  Bits bits_;
};

}  // namespace pixie::detail::sequence
/// @endcond
