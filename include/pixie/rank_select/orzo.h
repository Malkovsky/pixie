#pragma once

#ifndef PIXIE_ORZO_SUPPORT
#error "OrzoRankSelectSupport requires the optional BMI/BMI2-enabled backend"
#endif

#include <orzo/orzo.h>
#include <pixie/rank_select.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace pixie {

/**
 * @brief Optional rank/select implementation over the default Orzo
 * configuration.
 * @details Owns a masked source copy padded to whole 5,632-bit lower blocks.
 * This keeps upstream construction and select within complete basic blocks,
 * including partial logical tails. Uses Orzo's half-open one rank and one-based
 * one select. Zero rank is derived because upstream rank0 adds one to its
 * half-open answer; zero select is unavailable. One-select occurrence zero
 * returns zero and absent one occurrences return size(). Requires BMI and BMI2.
 * The prepared upstream header fixes directory allocation rounding and its
 * terminal entry; query code is unchanged. Copy and move are disabled to avoid
 * copying upstream's raw owning pointers.
 */
class OrzoRankSelectSupport : public RankSelectBase<OrzoRankSelectSupport> {
 public:
  /**
   * @brief Build Orzo from least-significant-bit-first packed words.
   * @param source_words External complete 64-bit words, copied during
   * construction.
   * @param num_bits Logical bit count, clamped to available words.
   * @throws std::length_error If source padding is unrepresentable.
   */
  explicit OrzoRankSelectSupport(std::span<const std::uint64_t> source_words,
                                 std::size_t num_bits) {
    num_bits_ = num_bits;
    if (source_words.size() <= std::numeric_limits<std::size_t>::max() / 64) {
      num_bits_ = std::min(num_bits_, source_words.size() * 64);
    }
    if (num_bits_ == 0) {
      return;
    }
    if (num_bits_ >
        std::numeric_limits<std::size_t>::max() - (kLowerBits - 1)) {
      throw std::length_error("Orzo source padding is unrepresentable");
    }
    const std::size_t physical_bits =
        ((num_bits_ + kLowerBits - 1) / kLowerBits) * kLowerBits;
    words_.resize(physical_bits / 64, 0);
    const std::size_t logical_words = num_bits_ / 64 + (num_bits_ % 64 != 0);
    std::copy_n(source_words.begin(), logical_words, words_.begin());
    if (num_bits_ % 64 != 0) {
      words_[logical_words - 1] &= (std::uint64_t{1} << (num_bits_ % 64)) - 1;
    }

    std::size_t bucket_ones = 0;
    std::size_t sample_entries = 0;
    for (std::size_t word = 0; word < words_.size(); ++word) {
      bucket_ones += std::popcount(words_[word]);
      if ((word + 1) % (kSelectUpperBits / 64) == 0 ||
          word + 1 == words_.size()) {
        sample_entries +=
            std::max<std::size_t>(1, bucket_ones / kSampleFrequency +
                                         (bucket_ones % kSampleFrequency != 0));
        bucket_ones = 0;
      }
    }
    const std::size_t upper_count = ceiling(physical_bits, kUpperBits);
    const std::size_t lower_count = physical_bits / kLowerBits;
    const std::size_t select_upper_count =
        ceiling(physical_bits, kSelectUpperBits);
    metadata_bytes_ = (upper_count + 2) * sizeof(std::uint64_t) +
                      (lower_count + 1) * sizeof(__uint128_t) +
                      (select_upper_count + 2) * sizeof(std::uint64_t) +
                      select_upper_count * sizeof(std::uint32_t*) +
                      sample_entries * sizeof(std::uint32_t);
    support_.emplace(words_.data(), physical_bits);
    one_count_ = support_->get_one_count();
  }

  OrzoRankSelectSupport(const OrzoRankSelectSupport&) = delete;
  OrzoRankSelectSupport& operator=(const OrzoRankSelectSupport&) = delete;
  OrzoRankSelectSupport(OrzoRankSelectSupport&&) = delete;
  OrzoRankSelectSupport& operator=(OrzoRankSelectSupport&&) = delete;

  /** @brief Return the logical, unpadded source length. */
  std::size_t size_impl() const { return num_bits_; }
  /** @brief Return a bit; position must be below size(). */
  int bit_impl(std::size_t position) const {
    return (words_[position / 64] >> (position % 64)) & 1;
  }

  /** @brief Count ones in [0, end), clamping end to size(). */
  std::uint64_t rank_impl(std::size_t end) const {
    return end >= num_bits_ ? one_count_ : support_->rank1(words_.data(), end);
  }

  /** @brief Return a one-based one's position, zero for zero, size() if absent.
   */
  std::uint64_t select_impl(std::size_t occurrence) const {
    if (occurrence == 0) {
      return 0;
    }
    return occurrence > one_count_
               ? num_bits_
               : support_->select1(words_.data(), occurrence);
  }

  /** @brief Orzo always builds one-select support in this configuration. */
  bool supports_select1_impl() const { return true; }
  /** @brief The pinned upstream implementation does not provide zero select. */
  bool supports_select0_impl() const { return false; }
  /** @brief Return zero for occurrence zero, otherwise size(): disabled. */
  std::uint64_t select0_impl(std::size_t occurrence) const {
    return occurrence == 0 ? 0 : size();
  }
  /** @brief Return allocated source-copy bytes, including zero padding. */
  std::size_t source_copy_bytes() const {
    return words_.capacity() * sizeof(std::uint64_t);
  }
  /** @brief Return wrapper, upstream object, and requested directory bytes. */
  std::size_t index_logical_bytes() const {
    return sizeof(*this) + metadata_bytes_;
  }
  /** @brief Return object and allocated index/source bytes. */
  std::size_t memory_usage_bytes_impl() const {
    return index_logical_bytes() + source_copy_bytes();
  }

 private:
  using Index = Orzo<>;
  static constexpr std::size_t kLowerBits = 5632;
  static constexpr std::size_t kUpperBits = 259072;
  static constexpr std::size_t kSelectUpperBits = 4294895616;
  static constexpr std::size_t kSampleFrequency = 8192;

  static std::size_t ceiling(std::size_t size, std::size_t block) {
    return size / block + (size % block != 0);
  }

  // Upstream's queries lack const signatures, but read source and index only.
  mutable std::vector<std::uint64_t> words_;
  mutable std::optional<Index> support_;
  std::size_t num_bits_ = 0;
  std::size_t one_count_ = 0;
  std::size_t metadata_bytes_ = 0;
};

}  // namespace pixie
