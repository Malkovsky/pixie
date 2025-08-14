#pragma once

#include <algorithm>
#include <bit>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "bits.h"

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
  constexpr static size_t kWordSize = 64;
  constexpr static size_t kSuperBlockRankIntSize = 64;
  constexpr static size_t kBasicBlockRankIntSize = 16;
  constexpr static size_t kBasicBlockSize = 512;
  constexpr static size_t kSuperBlockSize = 65536;
  constexpr static size_t kWordsPerBlock = 8;


  std::vector<uint64_t> super_block_rank;
  std::vector<uint16_t> basic_block_rank;
  const size_t num_bits_;
  const size_t padded_size_;

  std::span<const uint64_t> bits;

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
    size_t word_idx = pos / kWordSize;
    size_t bit_off = pos % kWordSize;

    return (bits[word_idx] >> bit_off) & 1;
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
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = pos / kSuperBlockSize;
    // Precomputed rank
    uint64_t result = super_block_rank[s_block] + basic_block_rank[b_block];
    // Basic block tail
    result += popcount_512(&bits[b_block * kWordsPerBlock],
                           pos - (b_block * kBasicBlockSize));
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

/**
 * Iterleaved owning version of the bit vector class with rank and
 * select operations. Interleaving in SDS is the optimization
 * that stores auxiliary index and original in the same array which
 * reduces the average number of cache misses.
 *
 * If we can afford copying original data then interleaved version
 * is prefered over regular one.
 *
 * Implementation follows the paper
 *
 * "SPIDER: Improved Succint Rank and Select Performance"
 * Matthew D. Laws, Jocelyn Bliven, Kit Conklin, Elyes Laalai, Samuel McCauley,
 * Zach S. Sturdevant
 *
 */
class BitVectorInterleaved {
private:
  constexpr static size_t kWordSize = 64;
  constexpr static size_t kSuperBlockRankIntSize = 64;
  constexpr static size_t kBasicBlockRankIntSize = 16;
  /**
   * 496 bits data + 16 bit local rank
   */
  constexpr static size_t kBasicBlockSize = 496;
  /**
   * 63488 = 496 * 128, so position of superblock can be obtained
   * from the position of basic block by dividing on 128 or
   * right shift on 7 bits which is cheaper then performing another
   * division.
   */
  constexpr static size_t kSuperBlockSize = 63488;
  constexpr static size_t kBlocksPerSuperBlock = 128;
  constexpr static size_t kWordsPerBlock = 8;
  
  const size_t num_bits_;
  std::vector<uint64_t> bits_interleaved;
  std::vector<uint64_t> super_block_rank;
  
  class BitReader {
    size_t iterator_64_ = 0;
    size_t offset_size_ = 0;
    size_t offset_bits_ = 0;
    std::span<const uint64_t> bits_;

  public:
    BitReader(std::span<const uint64_t> bits) : bits_(bits) {}
    uint64_t ReadBits64(size_t num_bits) {
      if (num_bits > 64) {
        num_bits = 64;
      }
      uint64_t result = offset_bits_ & first_bits_mask(num_bits);
      if (offset_size_ >= num_bits) {
        offset_bits_ >>= num_bits;
        offset_size_ -= num_bits;
        return result;
      }
      uint64_t next = bits_[iterator_64_++];
      result ^= (next & first_bits_mask(num_bits - offset_size_))
                << offset_size_;
      offset_bits_ = (num_bits - offset_size_ == 64)
                         ? 0
                         : next >> (num_bits - offset_size_);
      offset_size_ = 64 - (num_bits - offset_size_);
      return result;
    }
  };

public:
  /**
   * Constructor from an external array of uint64_t.
   */
  explicit BitVectorInterleaved(std::span<const uint64_t> bit_vector,
                                size_t num_bits)
      : num_bits_(std::min(num_bits, bit_vector.size() * kWordSize)) {
    build_rank_interleaved(bit_vector, num_bits);
  }

  static inline uint64_t first_bits_mask(size_t num) {
    return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
  }

  /**
   * Returns the number of bits in the vector
   */
  size_t size() const { return num_bits_; }

  /**
   * Get the bit at the specified position
   */
  int operator[](size_t pos) const {
    size_t block_id = pos / kBasicBlockSize;
    size_t block_bit = pos - block_id * kBasicBlockSize;
    size_t word_id = block_id * kWordsPerBlock + block_bit / kWordSize;
    size_t word_bit = block_bit % kWordSize;
    kWordSize;

    return (bits_interleaved[word_id] >> word_bit) & 1;
  }

  /**
   * Precompute rank for fast queries
   */
  void build_rank_interleaved(std::span<const uint64_t> bits, size_t num_bits) {
    size_t num_superblocks = 1 + (num_bits_ - 1) / kSuperBlockSize;
    super_block_rank.resize(num_superblocks);
    size_t num_basicblocks = 1 + (num_bits_ - 1) / kBasicBlockSize;
    bits_interleaved.resize(num_basicblocks * (512 / kWordSize));

    uint64_t super_block_sum = 0;
    uint16_t basic_block_sum = 0;
    auto bit_reader = BitReader(bits);

    for (size_t i = 0; i * kBasicBlockSize < num_bits; ++i) {
      if (i % (kSuperBlockSize / kBasicBlockSize) == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank[i / (kSuperBlockSize / kBasicBlockSize)] =
            super_block_sum;
        basic_block_sum = 0;
      }
      bits_interleaved[i * (kWordsPerBlock) + 7] =
          static_cast<uint64_t>(basic_block_sum) << 48;

      for (size_t j = 0; j < 7 && kWordSize * (i + j) < num_bits; ++j) {
        bits_interleaved[i * (kWordsPerBlock) + j] = bit_reader.ReadBits64(
            std::min(64ull, num_bits - i * kBasicBlockSize + j * kWordSize));
        basic_block_sum +=
            std::popcount(bits_interleaved[i * (kWordsPerBlock) + j]);
      }
      if ((i + 7) * kWordSize < num_bits) {
        auto v = bit_reader.ReadBits64(
            std::min(48ull, num_bits - (i * kBasicBlockSize + 7 * kWordSize)));
        bits_interleaved[i * (kWordsPerBlock) + 7] ^= v;
        basic_block_sum += std::popcount(v);
      }
    }
  }

  /**
   * Rank operation: count the number of 1-bits up to position pos (exclusive)
   * rank_1(pos) = number of 1s in positions [0...pos-1]
   */
  uint64_t rank(size_t pos) const {
    // Super block rank
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = b_block / kBlocksPerSuperBlock;
    // Basic block rank. kSuperBlockSize / kBasicBlockSize is 128
    // and is known at compile time
    uint64_t result =
        super_block_rank[s_block];
    result += bits_interleaved[b_block * 8 + 7] >> 48;

    result += popcount_512(&bits_interleaved[b_block * kWordsPerBlock],
                           pos - (b_block * kBasicBlockSize));
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
