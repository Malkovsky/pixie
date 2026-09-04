#pragma once

/**
 * @file monotone.h
 * @brief Validated monotone adapters for immutable integer vectors.
 */

#include <pixie/detail/serialization.h>
#include <pixie/integer_vector.h>
#include <pixie/memory_usage.h>
#include <pixie/serialization.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

namespace pixie {

/**
 * @brief Owning duplicate-preserving monotone adapter over an integer vector.
 *
 * @tparam Vector Immutable integer-vector implementation owned by value.
 *
 * @details Public construction scans once and rejects decreasing adjacent
 * values. The underlying vector remains private so the nondecreasing invariant
 * cannot be invalidated. Bounds use logarithmic index-based binary search for
 * this adapter.
 */
template <class Vector>
  requires requires { typename Vector::value_type; } &&
               IntegerVectorValue<typename Vector::value_type> &&
               std::derived_from<
                   Vector,
                   IntegerVectorBase<Vector, typename Vector::value_type>>
class MonotoneIntegerVector
    : public MonotoneIntegerVectorBase<MonotoneIntegerVector<Vector>,
                                       typename Vector::value_type>,
      public SerializationBase<MonotoneIntegerVector<Vector>> {
 public:
  /** @brief Stored unsigned integer type inherited from the nested vector. */
  using value_type = typename Vector::value_type;

  /** @brief Type used for element counts and positions. */
  using size_type = std::size_t;

  /**
   * @brief Take ownership of @p vector after validating nondecreasing order.
   * @throws std::invalid_argument if any value is smaller than its predecessor.
   */
  explicit MonotoneIntegerVector(Vector vector) : vector_(std::move(vector)) {
    validate_monotone();
  }

  /** @brief Return the number of retained values. */
  size_type size_impl() const { return vector_.size(); }

  /** @brief Read one valid zero-based position without bounds checking. */
  value_type value_at_impl(size_type position) const {
    return vector_[position];
  }

  /** @brief Copy a range already validated by the public facade. */
  void copy_to_impl(size_type begin, std::span<value_type> output) const {
    vector_.copy_to(begin, output);
  }

  /** @brief Return inline bytes plus memory owned below the nested vector. */
  size_type memory_usage_bytes_impl() const
    requires requires(const Vector& vector) {
      { vector.memory_usage_bytes() } -> std::convertible_to<size_type>;
    }
  {
    return sizeof(*this) + nested_owned_memory_bytes(vector_);
  }

  /** @brief Return the first index whose value is at least @p x. */
  size_type lower_bound_index_impl(value_type x) const {
    size_type first = 0;
    size_type count = vector_.size();
    while (count != 0) {
      const size_type step = count / 2;
      const size_type middle = first + step;
      if (vector_[middle] < x) {
        first = middle + 1;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    return first;
  }

  /** @brief Return the first index whose value is greater than @p x. */
  size_type upper_bound_index_impl(value_type x) const {
    size_type first = 0;
    size_type count = vector_.size();
    while (count != 0) {
      const size_type step = count / 2;
      const size_type middle = first + step;
      if (x < vector_[middle]) {
        count = step;
      } else {
        first = middle + 1;
        count -= step + 1;
      }
    }
    return first;
  }

  /**
   * @brief Write a version-1 monotone wrapper around one nested vector
   * artifact.
   * @throws std::invalid_argument if the wrapper would not start at an
   * eight-byte-aligned writer offset.
   */
  void serialize_impl(BinaryWriter& writer) const
    requires Serializable<Vector>
  {
    if (writer.size_bytes() % alignof(std::uint64_t) != 0) {
      throw std::invalid_argument(
          "Monotone integer-vector serialization requires an aligned offset");
    }
    const std::size_t artifact_begin = writer.size_bytes();
    detail::write_magic(writer, kSerializationMagic);
    writer.write_u32(kSerializationVersion);
    writer.write_u32(0);
    const std::size_t size_position = writer.write_u64_placeholder();
    vector_.serialize(writer);
    writer.patch_u64(size_position, static_cast<std::uint64_t>(
                                        writer.size_bytes() - artifact_begin));
  }

  /**
   * @brief Restore one framed monotone vector and advance @p reader on success.
   *
   * @details Quick validation checks both frames but trusts the encoded
   * ordering. Full validation also scans all decoded values and rejects a
   * decrease. The nested vector determines whether restored storage is owned
   * or borrowed. Failure leaves @p reader unchanged.
   */
  static MonotoneIntegerVector deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires Deserializable<Vector, DeserializationValidation>
  {
    BinaryReader candidate = reader;
    const std::size_t available_size = candidate.remaining();
    detail::require_magic(candidate, kSerializationMagic);
    if (candidate.read_u32() != kSerializationVersion ||
        candidate.read_u32() != 0) {
      throw std::invalid_argument(
          "Incompatible serialized monotone integer vector");
    }
    const std::size_t artifact_size = detail::checked_artifact_size(
        candidate.read_u64(), kSerializationHeaderBytes, available_size);
    BinaryReader payload =
        candidate.read_subreader(artifact_size - kSerializationHeaderBytes);
    Vector vector = Vector::deserialize(payload, validation);
    payload.require_zero_padding(0);

    MonotoneIntegerVector result(std::move(vector), TrustedTag{});
    if (validation == DeserializationValidation::kFull) {
      result.validate_monotone();
    }
    reader = candidate;
    return result;
  }

 private:
  struct TrustedTag {};

  inline static constexpr std::array<std::uint8_t, 8> kSerializationMagic = {
      'P', 'X', 'M', 'O', 'N', 'O', 'V', '1'};
  static constexpr std::uint32_t kSerializationVersion = 1;
  static constexpr std::size_t kSerializationHeaderBytes = 24;

  MonotoneIntegerVector(Vector vector, TrustedTag)
      : vector_(std::move(vector)) {}

  void validate_monotone() const {
    for (size_type i = 1; i < vector_.size(); ++i) {
      if (vector_[i] < vector_[i - 1]) {
        throw std::invalid_argument(
            "Monotone integer vector contains a decreasing pair");
      }
    }
  }

  Vector vector_;
};

}  // namespace pixie
