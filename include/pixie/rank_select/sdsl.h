#pragma once

#ifndef SDSL_SUPPORT
#error "SdslRankSelectSupport requires PIXIE_THIRD_PARTY_BACKENDS"
#endif

#include <pixie/rank_select.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
// SDSL undefines ALIGNMENT; preserve the optional 3-star configuration.
#pragma push_macro("ALIGNMENT")
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#pragma pop_macro("ALIGNMENT")
#include <span>

namespace pixie {

/**
 * @brief Optional SDSL rank-v index with MCL select in both directions.
 * @details Owns a normalized source copy. Prefix ranks are half-open and
 * selects use one-based occurrences with zero-based positions and Pixie's
 * sentinels. Copy and move are disabled because SDSL support objects point into
 * bits_.
 */
class SdslRankSelectSupport : public RankSelectBase<SdslRankSelectSupport> {
 public:
  /**
   * @brief Copy packed source words and build rank1, select1 and select0.
   * @param source_words Complete least-significant-bit-first packed words.
   * @param num_bits Logical bit count, clamped to the supplied source capacity.
   * @details Ignores dirty trailing bits. Allocation failure propagates.
   */
  explicit SdslRankSelectSupport(std::span<const std::uint64_t> source_words,
                                 std::size_t num_bits)
      : bits_(copy_source(source_words, num_bits)),
        rank_(&bits_),
        select1_(&bits_),
        select0_(&bits_),
        one_count_(bits_.empty() ? 0 : rank_(bits_.size())) {}

  SdslRankSelectSupport(const SdslRankSelectSupport&) = delete;
  SdslRankSelectSupport& operator=(const SdslRankSelectSupport&) = delete;
  SdslRankSelectSupport(SdslRankSelectSupport&&) = delete;
  SdslRankSelectSupport& operator=(SdslRankSelectSupport&&) = delete;

  /** @brief Return the logical bit count. */
  std::size_t size_impl() const { return bits_.size(); }

  /** @brief Read a zero-based bit; position must be below size(). */
  int bit_impl(std::size_t position) const { return bits_[position]; }

  /** @brief Count ones in [0, end), clamping end to size(). */
  std::uint64_t rank_impl(std::size_t end) const {
    return end >= size() ? one_count_ : rank_(end);
  }

  /** @brief Locate a one-based one occurrence, or size() if absent. */
  std::uint64_t select_impl(std::size_t occurrence) const {
    return occurrence == 0           ? 0
           : occurrence > one_count_ ? size()
                                     : select1_(occurrence);
  }

  /** @brief Locate a one-based zero occurrence, or size() if absent. */
  std::uint64_t select0_impl(std::size_t occurrence) const {
    return occurrence == 0                    ? 0
           : occurrence > size() - one_count_ ? size()
                                              : select0_(occurrence);
  }

  /** @brief SDSL constructs one-select support unconditionally. */
  bool supports_select1_impl() const { return true; }

  /** @brief SDSL constructs zero-select support unconditionally. */
  bool supports_select0_impl() const { return true; }

  /** @brief Return allocated bytes in the owned source buffer. */
  std::size_t source_copy_bytes() const { return bits_.bit_capacity() / 8; }

  /**
   * @brief Return object bytes plus SDSL's serialized logical directory sizes.
   * @details Uses SDSL's size_in_bytes convention, excluding the owned source
   * payload and allocator bookkeeping. This inspects the directories; benchmark
   * callers must collect it outside measured construction/query work.
   */
  std::size_t index_logical_bytes() const {
    return sizeof(*this) + sdsl::size_in_bytes(rank_) +
           sdsl::size_in_bytes(select1_) + sdsl::size_in_bytes(select0_);
  }

  /** @brief Return the reported index footprint and owned source bytes. */
  std::size_t memory_usage_bytes_impl() const {
    return source_copy_bytes() + index_logical_bytes();
  }

 private:
  static sdsl::bit_vector copy_source(
      std::span<const std::uint64_t> source_words,
      std::size_t num_bits) {
    if (source_words.size() <= std::numeric_limits<std::size_t>::max() / 64) {
      num_bits = std::min(num_bits, source_words.size() * 64);
    }
    sdsl::bit_vector result(num_bits, 0);
    const std::size_t words = num_bits / 64 + (num_bits % 64 != 0);
    std::copy_n(source_words.begin(), words, result.data());
    if (num_bits % 64 != 0) {
      result.data()[words - 1] &= (std::uint64_t{1} << (num_bits % 64)) - 1;
    }
    return result;
  }

  sdsl::bit_vector bits_;
  sdsl::rank_support_v<> rank_;
  sdsl::select_support_mcl<1> select1_;
  sdsl::select_support_mcl<0> select0_;
  std::size_t one_count_;
};

}  // namespace pixie
