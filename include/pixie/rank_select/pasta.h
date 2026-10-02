#pragma once

#ifndef PIXIE_PASTA_SUPPORT
#error "PastaRankSelectSupport requires PIXIE_THIRD_PARTY_BACKENDS"
#endif

#include <pixie/rank_select.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <pasta/bit_vector/support/l12_type.hpp>
#include <span>

namespace pixie {

/**
 * @brief Pasta flat rank/select adapter over a copied packed bit sequence.
 *
 * @details The constructor owns a normalized copy of the logical source words;
 * source mutations do not affect this index. Both select directions are always
 * available. Public positions are zero-based, prefix ranks use `[0, end)`, and
 * selects use one-based ranks with the Pixie rank-zero and out-of-range
 * sentinels.
 */
class PastaRankSelectSupport : public RankSelectBase<PastaRankSelectSupport> {
 private:
  using PastaSupport = pasta::FlatRankSelect<>;

  static uint64_t logical_word(std::span<const uint64_t> source,
                               std::size_t word_index,
                               std::size_t num_bits) {
    uint64_t word = source[word_index];
    const std::size_t remaining = num_bits - word_index * 64;
    if (remaining < 64) {
      word &= (uint64_t{1} << remaining) - 1;
    }
    return word;
  }

  pasta::BitVector bits_;
  PastaSupport support_;
  std::size_t num_bits_ = 0;
  std::size_t max_rank_ = 0;

 public:
  /**
   * @brief Construct an owning Pasta index from caller-provided packed words.
   * @param source_words Packed least-significant-bit-first source words.
   * @param num_bits Number of valid input bits, clamped to @p source_words.
   */
  explicit PastaRankSelectSupport(std::span<const uint64_t> source_words,
                                  std::size_t num_bits)
      : bits_(std::min(num_bits, source_words.size() * 64)),
        num_bits_(std::min(num_bits, source_words.size() * 64)) {
    std::fill(bits_.data().begin(), bits_.data().end(), 0);
    const std::size_t word_count = (num_bits_ + 63) / 64;
    for (std::size_t word_index = 0; word_index < word_count; ++word_index) {
      bits_.data()[word_index] =
          logical_word(source_words, word_index, num_bits_);
    }
    support_ = PastaSupport(bits_);
    max_rank_ = num_bits_ == 0 ? 0 : support_.rank1(num_bits_);
  }

  /** @brief Return the number of valid source bits. */
  std::size_t size_impl() const { return num_bits_; }

  /** @brief Return the bit at a zero-based valid position. */
  int bit_impl(std::size_t position) const {
    return static_cast<bool>(bits_[position]);
  }

  /** @brief Return the one rank in `[0, end_position)`. */
  uint64_t rank_impl(std::size_t end_position) const {
    return end_position >= num_bits_ ? max_rank_ : support_.rank1(end_position);
  }

  /** @brief Return the zero-based position of a one-based one rank. */
  uint64_t select_impl(std::size_t rank) const {
    if (rank == 0) {
      return 0;
    }
    return rank > max_rank_ ? num_bits_ : support_.select1(rank);
  }

  /** @brief Return the zero-based position of a one-based zero rank. */
  uint64_t select0_impl(std::size_t rank) const {
    if (rank == 0) {
      return 0;
    }
    return rank > num_bits_ - max_rank_ ? num_bits_ : support_.select0(rank);
  }

  /** @brief Pasta constructs support for one selects unconditionally. */
  bool supports_select1_impl() const { return true; }

  /** @brief Pasta constructs support for zero selects unconditionally. */
  bool supports_select0_impl() const { return true; }

  /**
   * @brief Return bytes owned by the converted source and Pasta metadata.
   * @details Pasta's `FlatRankSelect::space_usage()` omits its inherited L1/L2
   * table, so account for that table explicitly from its documented layout.
   */
  std::size_t memory_usage_bytes_impl() const {
    return source_copy_bytes() + index_logical_bytes() +
           (sizeof(*this) - sizeof(bits_) - sizeof(support_));
  }

  /** @brief Return bytes in Pasta's owned normalized source copy. */
  std::size_t source_copy_bytes() const {
    return bits_.space_usage() - sizeof(pasta::BitVector);
  }

  /**
   * @brief Return Pasta's logical index bytes, excluding its owned source.
   * @details Pasta reports logical vector sizes rather than allocator capacity.
   */
  std::size_t index_logical_bytes() const {
    const std::size_t l12_entries =
        bits_.data().size() / pasta::FlatRankSelectConfig::L1_WORD_SIZE + 1;
    return (support_.space_usage() - sizeof(PastaSupport)) +
           l12_entries * sizeof(pasta::BigL12Type);
  }
};

}  // namespace pixie
