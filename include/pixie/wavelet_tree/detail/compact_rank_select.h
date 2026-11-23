#pragma once

#include <pixie/bits.h>
#include <pixie/storage.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

namespace pixie::detail {

/**
 * @brief Directly queryable rank/select view over compact wavelet-node
 * metadata.
 *
 * @details The serialized representation omits derivable first basic-block
 * ranks, single-superblock ranks, and initial select samples. Large blocks
 * retain the same 512-bit rank/select primitives as RankSelectSupport; only
 * the metadata search layout differs.
 */
template <StorageImplementation Storage>
class CompactWaveletRankSelect {
 public:
  static constexpr std::size_t kWordBits = 64;
  static constexpr std::size_t kBasicBlockBits = 512;
  static constexpr std::size_t kWordsPerBasicBlock = 8;
  static constexpr std::size_t kSuperBlockBits = 65536;
  static constexpr std::size_t kBasicBlocksPerSuperBlock = 128;
  static constexpr std::size_t kSelectSampleFrequency = 16384;

  CompactWaveletRankSelect() = default;

  CompactWaveletRankSelect(std::span<const std::uint64_t> source_words,
                           std::size_t bit_count,
                           std::size_t one_count,
                           Storage basic_ranks,
                           Storage super_ranks,
                           Storage select_samples)
      : bits_(source_words),
        basic_ranks_(std::move(basic_ranks)),
        super_ranks_(std::move(super_ranks)),
        select_samples_(std::move(select_samples)),
        bit_count_(bit_count),
        one_count_(one_count) {
    validate_layout();
  }

  CompactWaveletRankSelect(std::span<const std::uint64_t> source_words,
                           const CompactWaveletRankSelect& source)
      : bits_(source_words),
        basic_ranks_(source.basic_ranks_),
        super_ranks_(source.super_ranks_),
        select_samples_(source.select_samples_),
        bit_count_(source.bit_count_),
        one_count_(source.one_count_) {
    validate_layout();
  }

  CompactWaveletRankSelect(const CompactWaveletRankSelect&) = default;
  CompactWaveletRankSelect(CompactWaveletRankSelect&&) noexcept = default;
  CompactWaveletRankSelect& operator=(const CompactWaveletRankSelect&) =
      default;
  CompactWaveletRankSelect& operator=(CompactWaveletRankSelect&&) noexcept =
      default;

  [[nodiscard]] std::size_t size() const { return bit_count_; }
  [[nodiscard]] std::size_t one_count() const { return one_count_; }

  [[nodiscard]] int bit(std::size_t position) const {
    return static_cast<int>(
        (bits_[position / kWordBits] >> (position % kWordBits)) & 1U);
  }

  [[nodiscard]] std::size_t rank(std::size_t position) const {
    if (position >= bit_count_) [[unlikely]] {
      return one_count_;
    }
    const std::size_t basic_block = position / kBasicBlockBits;
    const std::size_t super_block = position / kSuperBlockBits;
    return super_rank(super_block) + basic_rank(basic_block) +
           rank_in_basic_block(basic_block, position % kBasicBlockBits);
  }

  [[nodiscard]] std::size_t rank0(std::size_t position) const {
    return position - rank(position);
  }

  [[nodiscard]] std::size_t select(std::size_t wanted_rank) const {
    if (wanted_rank == 0) [[unlikely]] {
      return 0;
    }
    if (wanted_rank > one_count_) [[unlikely]] {
      return bit_count_;
    }

    const std::size_t super_block = find_super_block(wanted_rank, false);
    std::size_t local_rank = wanted_rank - super_rank(super_block);
    const std::size_t basic_block =
        find_basic_block(local_rank, super_block, false);
    local_rank -= basic_rank(basic_block);
    return select_in_basic_block(basic_block, local_rank, false);
  }

  [[nodiscard]] std::size_t select0(std::size_t wanted_rank) const {
    if (wanted_rank == 0) [[unlikely]] {
      return 0;
    }
    if (wanted_rank > bit_count_ - one_count_) [[unlikely]] {
      return bit_count_;
    }

    const std::size_t super_block = find_super_block(wanted_rank, true);
    std::size_t local_rank =
        wanted_rank - (super_block * kSuperBlockBits - super_rank(super_block));
    const std::size_t basic_block =
        find_basic_block(local_rank, super_block, true);
    const std::size_t block_in_super = basic_block % kBasicBlocksPerSuperBlock;
    local_rank -= block_in_super * kBasicBlockBits - basic_rank(basic_block);
    return select_in_basic_block(basic_block, local_rank, true);
  }

  /** @brief Validate every compact prefix and sample against the source bits.
   */
  void validate_contents() const {
    const auto basic = basic_ranks_.as_words16();
    const auto super = super_ranks_.as_words64();
    const auto samples = select_samples_.as_words64();
    const std::size_t super_count = super_block_count(bit_count_);
    std::size_t rank = 0;
    std::size_t super_begin_rank = 0;
    std::size_t basic_index = 0;
    for (std::size_t block = 0; block < basic_block_count(bit_count_);
         ++block) {
      if (block % kBasicBlocksPerSuperBlock == 0) {
        const std::size_t super_block = block / kBasicBlocksPerSuperBlock;
        if (super_count > 1 && super[super_block] != rank) {
          throw std::invalid_argument(
              "Compact wavelet superblock rank disagrees with its bits");
        }
        super_begin_rank = rank;
      } else if (basic[basic_index++] != rank - super_begin_rank) {
        throw std::invalid_argument(
            "Compact wavelet basic-block rank disagrees with its bits");
      }

      const std::size_t first_word = block * kWordsPerBasicBlock;
      const std::size_t words =
          std::min(kWordsPerBasicBlock, word_count(bit_count_) - first_word);
      for (std::size_t word = 0; word < words; ++word) {
        rank += std::popcount(logical_word(first_word + word));
      }
    }
    if (rank != one_count_ || (super_count > 1 && super.back() != one_count_)) {
      throw std::invalid_argument(
          "Compact wavelet total rank disagrees with its bits");
    }

    const std::size_t one_sample_count = one_count_ / kSelectSampleFrequency;
    const std::size_t zero_sample_count =
        (bit_count_ - one_count_) / kSelectSampleFrequency;
    std::size_t next_one = kSelectSampleFrequency;
    std::size_t next_zero = kSelectSampleFrequency;
    std::size_t one_sample = 0;
    std::size_t zero_sample = one_sample_count;
    std::size_t ones_before = 0;
    std::size_t zeros_before = 0;
    for (std::size_t word_index = 0; word_index < word_count(bit_count_);
         ++word_index) {
      const std::uint64_t word = logical_word(word_index);
      const std::size_t word_bits = logical_word_bits(word_index);
      const std::size_t ones = std::popcount(word);
      const std::size_t zeros = word_bits - ones;
      while (one_sample < one_sample_count && next_one <= ones_before + ones) {
        const std::size_t position =
            word_index * kWordBits +
            select_64(word, next_one - ones_before - 1);
        if (samples[one_sample++] != position / kSuperBlockBits) {
          throw std::invalid_argument(
              "Compact wavelet one-select sample disagrees with its bits");
        }
        if (one_sample < one_sample_count) {
          next_one = (one_sample + 1) * kSelectSampleFrequency;
        }
      }
      const std::uint64_t zero_word =
          ~word & first_bits_mask(static_cast<std::uint32_t>(word_bits));
      while (zero_sample < one_sample_count + zero_sample_count &&
             next_zero <= zeros_before + zeros) {
        const std::size_t position =
            word_index * kWordBits +
            select_64(zero_word, next_zero - zeros_before - 1);
        if (samples[zero_sample++] != position / kSuperBlockBits) {
          throw std::invalid_argument(
              "Compact wavelet zero-select sample disagrees with its bits");
        }
        if (zero_sample < one_sample_count + zero_sample_count) {
          next_zero =
              (zero_sample - one_sample_count + 1) * kSelectSampleFrequency;
        }
      }
      ones_before += ones;
      zeros_before += zeros;
    }
    if (one_sample != one_sample_count || zero_sample != samples.size()) {
      throw std::invalid_argument(
          "Compact wavelet select sample count disagrees with its bits");
    }
  }

  [[nodiscard]] static std::size_t word_bytes(std::size_t bit_count) {
    return word_count(bit_count) * sizeof(std::uint64_t);
  }

  [[nodiscard]] static std::size_t basic_rank_bytes(std::size_t bit_count) {
    return stored_basic_rank_count(bit_count) * sizeof(std::uint16_t);
  }

  [[nodiscard]] static std::size_t super_rank_bytes(std::size_t bit_count) {
    const std::size_t count = super_block_count(bit_count);
    return (count <= 1 ? 0 : count + 1) * sizeof(std::uint64_t);
  }

  [[nodiscard]] static std::size_t select_sample_bytes(std::size_t bit_count,
                                                       std::size_t one_count) {
    return (one_count / kSelectSampleFrequency +
            (bit_count - one_count) / kSelectSampleFrequency) *
           sizeof(std::uint64_t);
  }

 private:
  inline static constexpr std::array<std::uint16_t, 32> kBasicDeltas = [] {
    std::array<std::uint16_t, 32> result{};
    for (std::size_t index = 0; index < result.size(); ++index) {
      result[index] = static_cast<std::uint16_t>(index * kBasicBlockBits);
    }
    return result;
  }();

  static std::size_t word_count(std::size_t bit_count) {
    return bit_count == 0 ? 0 : 1 + (bit_count - 1) / kWordBits;
  }

  static std::size_t basic_block_count(std::size_t bit_count) {
    return bit_count == 0 ? 0 : 1 + (bit_count - 1) / kBasicBlockBits;
  }

  static std::size_t super_block_count(std::size_t bit_count) {
    return bit_count == 0 ? 0 : 1 + (bit_count - 1) / kSuperBlockBits;
  }

  static std::size_t stored_basic_rank_count(std::size_t bit_count) {
    return basic_block_count(bit_count) - super_block_count(bit_count);
  }

  static std::size_t basic_rank_storage_index(std::size_t basic_block) {
    return basic_block - basic_block / kBasicBlocksPerSuperBlock - 1;
  }

  void validate_layout() const {
    if (one_count_ > bit_count_ || bits_.size() != word_count(bit_count_) ||
        basic_ranks_.size_bytes() != basic_rank_bytes(bit_count_) ||
        super_ranks_.size_bytes() != super_rank_bytes(bit_count_) ||
        select_samples_.size_bytes() !=
            select_sample_bytes(bit_count_, one_count_)) {
      throw std::invalid_argument("Invalid compact wavelet rank/select layout");
    }

    const std::size_t super_count = super_block_count(bit_count_);
    if (super_count > 1) {
      const auto ranks = super_ranks_.as_words64();
      if (ranks.front() != 0 || ranks.back() != one_count_) {
        throw std::invalid_argument(
            "Invalid compact wavelet superblock endpoints");
      }
      for (std::size_t super = 0; super < super_count; ++super) {
        const std::size_t bits =
            std::min(kSuperBlockBits, bit_count_ - super * kSuperBlockBits);
        if (ranks[super] > ranks[super + 1] ||
            ranks[super + 1] - ranks[super] > bits) {
          throw std::invalid_argument(
              "Invalid compact wavelet superblock ranks");
        }
      }
    }

    const auto basic = basic_ranks_.as_words16();
    std::size_t basic_index = 0;
    const std::size_t blocks = basic_block_count(bit_count_);
    for (std::size_t super = 0; super < super_count; ++super) {
      const std::size_t first = super * kBasicBlocksPerSuperBlock;
      const std::size_t count =
          std::min(kBasicBlocksPerSuperBlock, blocks - first);
      std::size_t previous = 0;
      for (std::size_t block = 1; block < count; ++block) {
        const std::size_t value = basic[basic_index++];
        if (value < previous || value - previous > kBasicBlockBits ||
            value > block * kBasicBlockBits) {
          throw std::invalid_argument(
              "Invalid compact wavelet basic-block ranks");
        }
        previous = value;
      }
    }

    const auto samples = select_samples_.as_words64();
    const std::size_t one_samples = one_count_ / kSelectSampleFrequency;
    for (std::size_t index = 0; index < samples.size(); ++index) {
      const std::uint64_t sample = samples[index];
      if (sample >= super_count) {
        throw std::invalid_argument("Invalid compact wavelet select sample");
      }

      const bool zeros = index >= one_samples;
      const std::size_t ordinal = zeros ? index - one_samples + 1 : index + 1;
      const std::size_t wanted_rank = ordinal * kSelectSampleFrequency;
      const std::size_t ones_before = super_rank(sample);
      const std::size_t ones_after = super_rank(sample + 1);
      const std::size_t bits_before = sample * kSuperBlockBits;
      const std::size_t bits_after = sample + 1 == super_count
                                         ? bit_count_
                                         : (sample + 1) * kSuperBlockBits;
      const std::size_t rank_before =
          zeros ? bits_before - ones_before : ones_before;
      const std::size_t rank_after =
          zeros ? bits_after - ones_after : ones_after;
      if (rank_before >= wanted_rank || rank_after < wanted_rank) {
        throw std::invalid_argument("Invalid compact wavelet select sample");
      }
    }
  }

  std::size_t super_rank(std::size_t super_block) const {
    const std::size_t count = super_block_count(bit_count_);
    if (count <= 1) {
      return super_block == 0 ? 0 : one_count_;
    }
    return super_ranks_.as_words64()[super_block];
  }

  std::size_t basic_rank(std::size_t basic_block) const {
    if (basic_block % kBasicBlocksPerSuperBlock == 0) {
      return 0;
    }
    return basic_ranks_.as_words16()[basic_rank_storage_index(basic_block)];
  }

  std::size_t logical_word_bits(std::size_t word_index) const {
    const std::size_t begin = word_index * kWordBits;
    return begin >= bit_count_ ? 0 : std::min(kWordBits, bit_count_ - begin);
  }

  std::uint64_t logical_word(std::size_t word_index) const {
    const std::size_t bits = logical_word_bits(word_index);
    if (bits == 0) {
      return 0;
    }
    if (bits == kWordBits) {
      return bits_[word_index];
    }
    return bits_[word_index] & first_bits_mask(bits);
  }

  std::size_t rank_in_basic_block(std::size_t basic_block,
                                  std::size_t offset) const {
    if (offset == 0) {
      return 0;
    }
    const std::size_t first_word = basic_block * kWordsPerBasicBlock;
    if (first_word + kWordsPerBasicBlock <= bits_.size()) {
      return rank_512(&bits_[first_word], offset);
    }

    std::size_t result = 0;
    std::size_t word_index = first_word;
    while (offset >= kWordBits) {
      result += std::popcount(logical_word(word_index++));
      offset -= kWordBits;
    }
    if (offset != 0) {
      result +=
          std::popcount(logical_word(word_index) &
                        first_bits_mask(static_cast<std::uint32_t>(offset)));
    }
    return result;
  }

  std::size_t select_sample(std::size_t sample_index, bool zeros) const {
    if (sample_index == 0) {
      return 0;
    }
    const auto samples = select_samples_.as_words64();
    const std::size_t one_samples = one_count_ / kSelectSampleFrequency;
    const std::size_t offset = zeros ? one_samples : 0;
    return samples[offset + sample_index - 1];
  }

  std::size_t find_super_block(std::size_t wanted_rank, bool zeros) const {
    const std::size_t count = super_block_count(bit_count_);
    if (count <= 1) {
      return 0;
    }

    std::size_t left =
        select_sample(wanted_rank / kSelectSampleFrequency, zeros);
    const auto ranks = super_ranks_.as_words64();
    if (!zeros) {
      while (left + 7 < ranks.size()) {
        const std::size_t length = lower_bound_8x64(&ranks[left], wanted_rank);
        if (length < 8) {
          return left + length - 1;
        }
        left += 8;
      }
      if (left + 3 < ranks.size()) {
        const std::size_t length = lower_bound_4x64(&ranks[left], wanted_rank);
        if (length < 4) {
          return left + length - 1;
        }
        left += 4;
      }
      while (left < ranks.size() && ranks[left] < wanted_rank) {
        ++left;
      }
      return left - 1;
    }

    while (left < ranks.size() &&
           left * kSuperBlockBits - ranks[left] < wanted_rank) {
      ++left;
    }
    return left - 1;
  }

  std::size_t find_basic_block(std::size_t local_rank,
                               std::size_t super_block,
                               bool zeros) const {
    const std::size_t first = super_block * kBasicBlocksPerSuperBlock;
    const std::size_t count = std::min(kBasicBlocksPerSuperBlock,
                                       basic_block_count(bit_count_) - first);
    if (count == 1) {
      return first;
    }
    if (local_rank == kSuperBlockBits) {
      return first + count - 1;
    }

    const auto ranks = basic_ranks_.as_words16();
    const std::size_t stored_begin = super_block * 127;
    const std::size_t stored_count = count - 1;
    std::size_t less = 1;
    std::size_t position = 0;
    while (position + 32 <= stored_count) {
      const std::size_t found =
          zeros ? lower_bound_delta_32x16(&ranks[stored_begin + position],
                                          local_rank, kBasicDeltas.data(),
                                          kBasicBlockBits * (position + 1))
                : lower_bound_32x16(&ranks[stored_begin + position],
                                    static_cast<std::uint16_t>(local_rank));
      less += found;
      if (found != 32) {
        return first + less - 1;
      }
      position += 32;
    }
    for (; position < stored_count; ++position) {
      const std::size_t block = position + 1;
      const std::size_t prefix =
          zeros ? block * kBasicBlockBits - ranks[stored_begin + position]
                : ranks[stored_begin + position];
      if (prefix >= local_rank) {
        break;
      }
      ++less;
    }
    return first + less - 1;
  }

  std::size_t select_in_basic_block(std::size_t basic_block,
                                    std::size_t wanted_rank,
                                    bool zeros) const {
    const std::size_t first_word = basic_block * kWordsPerBasicBlock;
    const std::size_t first_bit = first_word * kWordBits;
    if (first_bit + kBasicBlockBits <= bit_count_ &&
        first_word + kWordsPerBasicBlock <= bits_.size()) {
      return first_bit +
             (zeros ? select0_512(&bits_[first_word], wanted_rank - 1)
                    : select_512(&bits_[first_word], wanted_rank - 1));
    }

    for (std::size_t word_index = first_word;
         word_index < word_count(bit_count_); ++word_index) {
      const std::uint64_t word = logical_word(word_index);
      const std::uint64_t candidates =
          zeros ? (~word & first_bits_mask(logical_word_bits(word_index)))
                : word;
      const std::size_t count = std::popcount(candidates);
      if (wanted_rank > count) {
        wanted_rank -= count;
        continue;
      }
      return word_index * kWordBits + select_64(candidates, wanted_rank - 1);
    }
    return bit_count_;
  }

  std::span<const std::uint64_t> bits_;
  Storage basic_ranks_;
  Storage super_ranks_;
  Storage select_samples_;
  std::size_t bit_count_ = 0;
  std::size_t one_count_ = 0;
};

}  // namespace pixie::detail
