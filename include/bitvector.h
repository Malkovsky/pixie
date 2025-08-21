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
 * @brief Non-interleaved non owning version of the bit vector class with rank
 * and select operations.
 *
 * @details
 * Implementation follows ideas from
 *
 * {1} "SPIDER: Improved Succint Rank and Select Performance"
 * Matthew D. Laws, Jocelyn Bliven, Kit Conklin, Elyes Laalai, Samuel McCauley,
 * Zach S. Sturdevant
 * https://github.com/williams-cs/spider
 *
 * {2} "Engineering compact data structures for rank and select queries on
 * bit vectors" Kurpicz F.
 * https://github.com/pasta-toolbox/bit_vector
 *
 * Essentially it is a 2 level structure of
 * - Super blocks of 2^16 bits
 * - Basic blocks of 512 bits
 * - Super block ranks of 64 bits resulting in ~ 0.98% overhead
 * - Basic block ranks of 16 bits resulting in ~ 3.125% overhead
 * - 64 bit Select samples of each 16384 bit ~ 0.39% overhead
 *
 * Rank is 2 table lookups + SIMD popcount in a 512 block.
 *
 * Select is
 * - Initial super block guess by precomputed sampling
 * - SIMD supported linear scan to find superblock
 * --- For random data we expect the scan to quickly find the block
 * --- As 8 64-bit ranks fits into a cache line we probably expect a single
 * --- cache miss to find a superblock. SIMD allows to perform a single scan
 * --- can of 8 ranks in 3-4 ops.
 * - SIMD supported linear scan to find basicblock
 * --- Super block consits of 128 basic blocks 16 bits rank each
 * --- This fits into 4x512-bit cache lines. The best way seems to
 * --- be just to linearly scan linearly, AVX-512 allowes to scan
 * --- the whole cache line at once.
 *
 * Currently I didn't compare it against {1} or {2}, performance is compatible
 * and is bottlenecked by cache misses. Usage of AVX-512 makes it limited for
 * current generation processors.
 *
 * This implementation does not involve interleaving, it is notable that
 * interleaving makes linear scans harder. I'd suggest that if interleaving is
 * to be involved it's probably super-basic rank interleaving like in {2}. We
 * lose ability to effectively perform linear scan on superblocks as we already
 * expect to have good guesses from select samples. Interleanign original data
 * with basic block ranks like in {1} some special guessing technique should be
 * involved (like {1} does).
 */
class BitVector {
private:
  constexpr static size_t kWordSize = 64;
  constexpr static size_t kSuperBlockRankIntSize = 64;
  constexpr static size_t kBasicBlockRankIntSize = 16;
  constexpr static size_t kBasicBlockSize = 512;
  constexpr static size_t kSuperBlockSize = 65536;
  constexpr static size_t kWordsPerBlock = 8;
  constexpr static size_t kSelectSampleFrequency = 16384;
  constexpr static size_t kBlocksPerSuperBlock = 128;

  std::vector<uint64_t> super_block_rank;
  std::vector<uint16_t> basic_block_rank;
  std::vector<uint64_t> select_samples;
  const size_t num_bits_;
  const size_t padded_size_;
  size_t max_rank;

  std::span<const uint64_t> bits;

  inline uint64_t first_bits_mask(size_t num) const {
    return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
  }

  /**
   * Precompute rank for fast queries
   */
  void build_rank() {
    size_t num_superblocks = 8 + (padded_size_ - 1) / kSuperBlockSize;
    // Add more blocks to ease SIMD processing
    // num_basicblocks to fully cover superblock, i.e. 128
    // This reduces branching in select
    num_superblocks = ((num_superblocks + 7) / 8) * 8;
    size_t num_basicblocks = num_superblocks * kBlocksPerSuperBlock;
    super_block_rank.resize(num_superblocks);
    basic_block_rank.resize(num_basicblocks);

    uint64_t super_block_sum = 0;
    uint16_t basic_block_sum = 0;

    for (size_t i = 0; i / kBasicBlockSize < basic_block_rank.size();
         i += kWordSize) {
      if (i % kSuperBlockSize == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank[i / kSuperBlockSize] = super_block_sum;
        basic_block_sum = 0;
      }
      if (i % kBasicBlockSize == 0) {
        basic_block_rank[i / kBasicBlockSize] = basic_block_sum;
      }
      if (i / kWordSize < bits.size()) {
        basic_block_sum += std::popcount(bits[i / kWordSize]);
      }
    }
    max_rank = super_block_sum + basic_block_sum;
  }

  /**
   * @brief Calculate select samples
   */
  void build_select() {
    uint64_t milestone = kSelectSampleFrequency;
    uint64_t rank = 0;
    select_samples.emplace_back(0);
    for (size_t i = 0; i < bits.size(); ++i) {
      auto ones = std::popcount(bits[i]);
      if (rank + ones >= milestone) {
        auto pos = select_64(bits[i], milestone - rank - 1);
        select_samples.emplace_back((64 * i + pos) / kSuperBlockSize);
        milestone += kSelectSampleFrequency;
      }
      rank += ones;
    }
  }

  /**
   * @brief first step of the select operation
   */
  uint64_t find_superblock(uint64_t rank) const {
    uint64_t left = select_samples[rank / kSelectSampleFrequency];
    auto rank_8 = _mm512_set1_epi64(rank);

    while (left + 7 < super_block_rank.size()) {
      auto reg_512 = _mm512_loadu_epi64(&super_block_rank[left]);
      auto cmp = _mm512_cmpge_epu64_mask(reg_512, rank_8);
      if (cmp) {
        return left + _tzcnt_u16(cmp) - 1;
      }
      left += 8;
    }
    if (left + 3 < super_block_rank.size()) {
      auto reg_256 = _mm256_loadu_epi64(&super_block_rank[left]);
      auto cmp = _mm256_cmpge_epu64_mask(reg_256, *(__m256i *)(&rank_8));
      if (cmp) {
        return left + _tzcnt_u16(cmp) - 1;
      }
      left += 4;
    }
    while (left < super_block_rank.size() && super_block_rank[left] < rank)
      left++;
    return left - 1;
  }

  uint64_t find_basicblock(uint16_t local_rank, uint64_t s_block) const {
    auto rank_32 = _mm512_set1_epi16(local_rank);
    auto pos = 0;

    for (size_t i = 0; i < 4; ++i) {
      auto batch =
          _mm512_loadu_epi16(&basic_block_rank[128 * s_block + 32 * i]);
      auto cmp = _mm512_cmplt_epu16_mask(batch, rank_32);
      auto count = std::popcount(cmp);
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
      pos += 32;
    }
    return kBlocksPerSuperBlock * s_block + pos - 1;
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
    build_select();
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
   * Rank operation: count the number of 1-bits up to position pos (exclusive)
   * rank_1(pos) = number of 1s in positions [0...pos-1]
   */
  uint64_t rank(size_t pos) const {
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = pos / kSuperBlockSize;
    // Precomputed rank
    uint64_t result = super_block_rank[s_block] + basic_block_rank[b_block];
    // Basic block tail
    result += rank_512(&bits[b_block * kWordsPerBlock],
                       pos - (b_block * kBasicBlockSize));
    return result;
  }

  uint64_t select(size_t rank) const {
    if (rank > max_rank) [[unlikely]] {
      return num_bits_;
    }
    if (rank == 0) [[unlikely]] {
      return 0;
    }
    uint64_t s_block = find_superblock(rank);
    rank -= super_block_rank[s_block];
    auto pos = find_basicblock(rank, s_block);
    rank -= basic_block_rank[pos];
    pos *= kWordsPerBlock;

    // Final search
    if (pos + kWordsPerBlock - 1 < kWordsPerBlock) [[unlikely]] {
      size_t ones = std::popcount(bits[pos]);
      while (pos < bits.size() && ones < rank) {
        rank -= ones;
        ones = std::popcount(bits[++pos]);
      }
      return kWordSize * pos + select_64(bits[pos], rank - 1);
    }
    return kWordSize * pos + select_512(&bits[pos], rank - 1);
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
    // Multiplication/devisions
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = b_block / kBlocksPerSuperBlock;
    uint64_t b_block_pos = b_block * kWordsPerBlock;
    // Super block rank
    uint64_t result = super_block_rank[s_block];
    /**
     * Ok, so here's quite the important factor to load 512-bit region
     * at &bits_interleaved[b_block_pos], we store local rank as 16 last
     * bits of it. Prefetch should guarantee but seems like there is no
     * need for it.
     */
    // __builtin_prefetch(&bits_interleaved[b_block_pos]);
    result += rank_512(&bits_interleaved[b_block_pos],
                       pos - (b_block * kBasicBlockSize));
    result += bits_interleaved[b_block_pos + 7] >> 48;
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
