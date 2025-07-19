#pragma once

#include <bit>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace pixie {

/**
 * Non-interleaved non owning version of the bit vector class with rank and
 * select operations.
 *
 * Implementation follows the paper
 *
 * "SPIDER: Improved Succint Rank and Select Performance"
 * Matthew D. Laws, Jocelyn Bliven, Kit Conklin, Elyes Laalai, Samuel McCauley,
 * Zach S. Sturdevant
 *
 */
class BitVector {
private:
  /**
   *
   */
  const size_t kWordSize = 64;
  const size_t kSuperBlockRankIntSize = 64;
  const size_t kBasicBlockRankIntSize = 16;
  std::vector<uint64_t> super_block_rank;
  std::vector<uint16_t> basic_block_rank;
  const size_t kBasicBlockSize = 512;
  const size_t kSuperBlockSize = 65536;

  const size_t num_bits_;
  const size_t padded_size_;

  std::span<const uint64_t> bits;

  inline size_t word_index(size_t bit_pos) const { return bit_pos / kWordSize; }

  inline size_t bit_offset(size_t bit_pos) const { return bit_pos % kWordSize; }

  inline uint64_t first_bits_mask(size_t num) const {
    return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
  }

public:
  /**
   * Constructor from an external array of uint64_t.
   */
  explicit BitVector(std::span<const uint64_t> bit_vector, size_t num_bits)
      : num_bits_(std::min(num_bits, bit_vector.size() * kWordSize)),
        padded_size_(((num_bits_ + kWordSize - 1) / kWordSize) * kWordSize),
        bits(bit_vector) {
    build_rank();
  }

  /**
   * Returns the number of bits in the vector
   */
  size_t size() const { return num_bits_; }

  /**
   * Get the bit at the specified position
   */
  int operator[](size_t pos) const {
    size_t word_idx = word_index(pos);
    size_t bit_off = bit_offset(pos);

    return (bits[word_idx] & (1ULL << bit_off)) != 0;
  }

  /**
   * Precompute rank for fast queries
   */
  void build_rank() {
    size_t num_superblocks = 1 + (padded_size_ - 1) / kSuperBlockSize;
    super_block_rank.resize(num_superblocks);
    size_t num_basicblocks = 1 + (padded_size_ - 1) / kBasicBlockSize;
    basic_block_rank.resize(num_basicblocks);

    uint64_t super_block_sum = 0;
    uint16_t basic_block_sum = 0;

    for (size_t i = 0; i + kWordSize < padded_size_; i += kWordSize) {
      if (i % kSuperBlockSize == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank[i / kSuperBlockSize] = super_block_sum;
        basic_block_sum = 0;
      }
      if (i % kBasicBlockSize == 0) {
        basic_block_rank[i / kBasicBlockSize] = basic_block_sum;
      }
      basic_block_sum += std::popcount(bits[i / kWordSize]);
    }
  }

  /**
   * Rank operation: count the number of 1-bits up to position pos (exclusive)
   * rank_1(pos) = number of 1s in positions [0...pos-1]
   */
  uint64_t rank(size_t pos) const {
    // Super block rank
    uint64_t result = super_block_rank[pos / kSuperBlockSize];
    // Basic block rank
    result += basic_block_rank[pos / kBasicBlockSize];

    // Assuming compiler is smart enough to optimize it into
    // corresponding SIMD popcount for 128/256/512 bits.
    // Should probably do it manually. Manual implementation is
    // __m512i vec = _mm512_loadu_epi64(&ptr64[i]);
    // vec = _mm512_popcnt_epi64(vec);
    // cnt = _mm512_add_epi64(cnt, vec);
    for (size_t i = pos - (pos % kBasicBlockSize); i < pos; i += kWordSize) {
      result += std::popcount(bits[i / kWordSize] & first_bits_mask(pos - i));
    }
    return result;
  }

  /**
   * Convert the bit vector to a binary string, for debug purpose only
   */
  std::string to_string() const {
    std::string result;
    result.reserve(num_bits_);

    for (size_t i = 0; i < num_bits_; i++) {
      result.push_back(operator[](i) ? '1' : '0');
    }

    return result;
  }
};

} // namespace pixie
