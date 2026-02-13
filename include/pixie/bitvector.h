#pragma once

#include <pixie/bits.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace pixie {

/**
 * @brief Non-interleaved, non-owning bit vector with rank and select.
 *
 *
 * @details
 * This is a two-level rank/select index for a bit vector stored
 * externally as
 * 64-bit words. The layout follows ideas from:
 *
 * {1}
 * "SPIDER: Improved Succinct Rank and Select Performance"
 * Matthew D. Laws,
 * Jocelyn Bliven, Kit Conklin, Elyes Laalai, Samuel McCauley,
 * Zach S.
 * Sturdevant
 * https://github.com/williams-cs/spider
 *
 * {2} "Engineering
 * compact data structures for rank and select queries on
 * bit vectors"
 * Kurpicz F.
 * https://github.com/pasta-toolbox/bit_vector
 *
 * Structure
 * overview:
 * - Super blocks of 2^16 bits with 64-bit ranks (~0.98%
 * overhead).
 * - Basic blocks of 512 bits with 16-bit ranks (~3.125%
 * overhead).
 * - Select samples every 16384 bits (~0.39% overhead).
 *
 *
 * Rank: 2 table lookups plus SIMD popcount in the 512-bit block.
 *
 * Select:

 * * - Start from a sampled super block.
 * - SIMD linear scan to find the super
 * block.
 * - SIMD linear scan to find the basic block.
 *
 * This variant does
 * not interleave data and index, favoring simpler scans.
 */
class BitVector {
 private:
  constexpr static size_t kWordSize = 64;
  constexpr static size_t kSuperBlockRankIntSize = 64;
  constexpr static size_t kBasicBlockRankIntSize = 16;
  constexpr static size_t kBasicBlockSize = 512;
  constexpr static size_t kWordsPerBlock = 8;
  constexpr static size_t kSuperBlockSize = 65536;
  constexpr static size_t kBlocksPerSuperBlock = 128;
  constexpr static size_t kSelectSampleFrequency = 16384;

  alignas(64) uint64_t delta_super[8];
  alignas(64) uint16_t delta_basic[32];

  std::vector<uint64_t> super_block_rank_;
  std::vector<uint16_t> basic_block_rank_;
  std::vector<uint64_t> select1_samples_;
  std::vector<uint64_t> select0_samples_;
  const size_t num_bits_;
  const size_t padded_size_;
  size_t max_rank_;

  std::span<const uint64_t> bits_;

  /**
   * @brief Precompute rank for fast queries.
   */
  void build_rank() {
    size_t num_superblocks = 8 + (padded_size_ - 1) / kSuperBlockSize;
    // Add more blocks to ease SIMD processing
    // num_basicblocks to fully cover superblock, i.e. 128
    // This reduces branching in select
    num_superblocks = ((num_superblocks + 7) / 8) * 8;
    size_t num_basicblocks = num_superblocks * kBlocksPerSuperBlock;
    super_block_rank_.resize(num_superblocks);
    basic_block_rank_.resize(num_basicblocks);

    uint64_t super_block_sum = 0;
    uint16_t basic_block_sum = 0;

    for (size_t i = 0; i / kBasicBlockSize < basic_block_rank_.size();
         i += kWordSize) {
      if (i % kSuperBlockSize == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank_[i / kSuperBlockSize] = super_block_sum;
        basic_block_sum = 0;
      }
      if (i % kBasicBlockSize == 0) {
        basic_block_rank_[i / kBasicBlockSize] = basic_block_sum;
      }
      if (i / kWordSize < bits_.size()) {
        basic_block_sum += std::popcount(bits_[i / kWordSize]);
      }
    }
    max_rank_ = super_block_sum + basic_block_sum;
  }

  /**
   * @brief Calculate select samples.
   */
  void build_select() {
    uint64_t milestone = kSelectSampleFrequency;
    uint64_t milestone0 = kSelectSampleFrequency;
    uint64_t rank = 0;
    uint64_t rank0 = 0;
    select1_samples_.emplace_back(0);
    select0_samples_.emplace_back(0);
    for (size_t i = 0; i < bits_.size(); ++i) {
      auto ones = std::popcount(bits_[i]);
      auto zeros = 64 - ones;
      if (rank + ones >= milestone) {
        auto pos = select_64(bits_[i], milestone - rank - 1);
        // TODO: try including global rank into select samples to save
        //       a cache miss on global rank scan
        select1_samples_.emplace_back((64 * i + pos) / kSuperBlockSize);
        milestone += kSelectSampleFrequency;
      }
      if (rank0 + zeros >= milestone0) {
        auto pos = select_64(~bits_[i], milestone0 - rank0 - 1);
        select0_samples_.emplace_back((64 * i + pos) / kSuperBlockSize);
        milestone0 += kSelectSampleFrequency;
      }
      rank += ones;
      rank0 += zeros;
    }

    for (size_t i = 0; i < 8; ++i) {
      delta_super[i] = i * kSuperBlockSize;
    }
    for (size_t i = 0; i < 32; ++i) {
      delta_basic[i] = i * kBasicBlockSize;
    }
  }

  /**
   * @brief First step of the select operation.
   * @param rank 1-based
   * rank of the 1-bit to locate.
   */
  uint64_t find_superblock(uint64_t rank) const {
    uint64_t left = select1_samples_[rank / kSelectSampleFrequency];

    while (left + 7 < super_block_rank_.size()) {
      auto len = lower_bound_8x64(&super_block_rank_[left], rank);
      if (len < 8) {
        return left + len - 1;
      }
      left += 8;
    }
    if (left + 3 < super_block_rank_.size()) {
      auto len = lower_bound_4x64(&super_block_rank_[left], rank);
      if (len < 4) {
        return left + len - 1;
      }
      left += 4;
    }
    while (left < super_block_rank_.size() && super_block_rank_[left] < rank) {
      left++;
    }
    return left - 1;
  }

  /**
   * @brief First step of the select0 operation.
   * @param rank0 1-based
   * rank of the 0-bit to locate.
   */
  uint64_t find_superblock_zeros(uint64_t rank0) const {
    uint64_t left = select0_samples_[rank0 / kSelectSampleFrequency];

    while (left + 7 < super_block_rank_.size()) {
      auto len = lower_bound_delta_8x64(&super_block_rank_[left], rank0,
                                        delta_super, kSuperBlockSize * left);
      if (len < 8) {
        return left + len - 1;
      }
      left += 8;
    }
    if (left + 3 < super_block_rank_.size()) {
      auto len = lower_bound_delta_4x64(&super_block_rank_[left], rank0,
                                        delta_super, kSuperBlockSize * left);
      if (len < 4) {
        return left + len - 1;
      }
      left += 4;
    }
    while (left < super_block_rank_.size() &&
           kSuperBlockSize * left - super_block_rank_[left] < rank0) {
      left++;
    }
    return left - 1;
  }

  /**
   * @brief SIMD-optimized linear scan.
   * @param local_rank Rank within
   * the super block.
   * @param s_block Super block index.
   * @details
   *
   * Processes 32 16-bit entries at once (full cache line), so there is at most

   * * 4 iterations.
   */
  uint64_t find_basicblock(uint16_t local_rank, uint64_t s_block) const {
    for (size_t pos = 0; pos < kBlocksPerSuperBlock; pos += 32) {
      auto count = lower_bound_32x16(
          &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank);
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
    }
    return kBlocksPerSuperBlock * s_block + kBlocksPerSuperBlock - 1;
  }

  /**
   * @brief SIMD-optimized linear scan.
   * @param local_rank0 Rank of
   * zeros within the super block.
   * @param s_block Super block index.
   *
   * @details
   * Processes 32 16-bit entries at once (full cache line), so
   * there is at most
   * 4 iterations.
   */
  uint64_t find_basicblock_zeros(uint16_t local_rank0, uint64_t s_block) const {
    for (size_t pos = 0; pos < kBlocksPerSuperBlock; pos += 32) {
      auto count = lower_bound_delta_32x16(
          &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank0,
          delta_basic, kBasicBlockSize * pos);
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
    }
    return kBlocksPerSuperBlock * s_block + kBlocksPerSuperBlock - 1;
  }

  /**
   * @brief Interpolation search with SIMD optimization.
   * @param
   * local_rank Rank within the super block.
   * @param s_block Super block
   * index.
   * @details
   * Similar to find_basicblock but initial guess is
   * based on linear
   * interpolation, for random data it should make initial
   * guess correct
   * most of the times, we start from the 32 wide block with
   * interpolation
   * guess at the center, if we see that select result lie in
   * lower blocks
   * we backoff to find_basicblock
   */
  uint64_t find_basicblock_is(uint16_t local_rank, uint64_t s_block) const {
    auto lower = super_block_rank_[s_block];
    auto upper = super_block_rank_[s_block + 1];

    uint64_t pos = kBlocksPerSuperBlock * local_rank / (upper - lower);
    pos = pos + 16 < 32 ? 0 : (pos - 16);
    pos = pos > 96 ? 96 : pos;
    while (pos < 96) {
      auto count = lower_bound_32x16(
          &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank);
      if (count == 0) {
        return find_basicblock(local_rank, s_block);
      }
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
      pos += 32;
    }
    pos = 96;
    auto count = lower_bound_32x16(
        &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank);
    if (count == 0) {
      return find_basicblock(local_rank, s_block);
    }
    return kBlocksPerSuperBlock * s_block + pos + count - 1;
  }

  /**
   * @brief Interpolation search with SIMD optimization.
   * @param
   * local_rank0 Rank of zeros within the super block.
   * @param s_block Super
   * block index.
   * @details
   * Similar to find_basicblock_zeros but
   * initial guess is based on linear
   * interpolation, for random data it
   * should make initial guess correct
   * most of the times, we start from the
   * 32 wide block with interpolation
   * guess at the center, if we see that
   * select result lie in lower blocks
   * we backoff to find_basicblock_zeros

   */
  uint64_t find_basicblock_is_zeros(uint16_t local_rank0,
                                    uint64_t s_block) const {
    auto lower = kSuperBlockSize * s_block - super_block_rank_[s_block];
    auto upper =
        kSuperBlockSize * (s_block + 1) - super_block_rank_[s_block + 1];

    uint64_t pos = kBlocksPerSuperBlock * local_rank0 / (upper - lower);
    pos = pos + 16 < 32 ? 0 : (pos - 16);
    pos = pos > 96 ? 96 : pos;
    while (pos < 96) {
      auto count = lower_bound_delta_32x16(
          &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank0,
          delta_basic, kBasicBlockSize * pos);
      if (count == 0) {
        return find_basicblock_zeros(local_rank0, s_block);
      }
      if (count < 32) {
        return kBlocksPerSuperBlock * s_block + pos + count - 1;
      }
      pos += 32;
    }
    pos = 96;
    auto count = lower_bound_delta_32x16(
        &basic_block_rank_[kBlocksPerSuperBlock * s_block + pos], local_rank0,
        delta_basic, kBasicBlockSize * pos);
    if (count == 0) {
      return find_basicblock_zeros(local_rank0, s_block);
    }
    return kBlocksPerSuperBlock * s_block + pos + count - 1;
  }

 public:
  /**
   * @brief Construct from an external array of 64-bit words.
   * @param
   * bit_vector Backing data, not owned.
   * @param num_bits Number of valid
   * bits in the vector.
   */
  explicit BitVector(std::span<const uint64_t> bit_vector, size_t num_bits)
      : num_bits_(std::min(num_bits, bit_vector.size() * kWordSize)),
        padded_size_(((num_bits_ + kWordSize - 1) / kWordSize) * kWordSize),
        bits_(bit_vector) {
    build_rank();
    build_select();
  }

  /**
   * @brief Returns the number of valid bits.
   */
  size_t size() const { return num_bits_; }

  /**
   * @brief Returns the bit at the given position.
   * @param pos Bit
   * index in [0, size()).
   */
  int operator[](size_t pos) const {
    size_t word_idx = pos / kWordSize;
    size_t bit_off = pos % kWordSize;

    return (bits_[word_idx] >> bit_off) & 1;
  }

  /**
   * @brief Rank of 1s up to position pos (exclusive).
   * @param pos Bit
   * index in [0, size()].
   * @return Number of 1s in [0, pos).
   */
  uint64_t rank(size_t pos) const {
    if (pos >= bits_.size() * kWordSize) [[unlikely]] {
      return max_rank_;
    }
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = pos / kSuperBlockSize;
    // Precomputed rank
    uint64_t result = super_block_rank_[s_block] + basic_block_rank_[b_block];
    // Basic block tail
    result += rank_512(&bits_[b_block * kWordsPerBlock],
                       pos - (b_block * kBasicBlockSize));
    return result;
  }

  /**
   * @brief Rank of 0s up to position pos (exclusive).
   * @param pos Bit
   * index in [0, size()].
   * @return Number of 0s in [0, pos).
   */
  uint64_t rank0(size_t pos) const {
    if (pos >= bits_.size() * kWordSize) [[unlikely]] {
      return bits_.size() * kWordSize - max_rank_;
    }
    return pos - rank(pos);
  }

  /**
   * @brief Select the position of the rank-th 1-bit (1-indexed).
   *
   * @param rank 1-based rank of the 1-bit to select.
   * @return Bit index, or
   * size() if rank is out of range.
   */
  uint64_t select(size_t rank) const {
    if (rank > max_rank_) [[unlikely]] {
      return num_bits_;
    }
    if (rank == 0) [[unlikely]] {
      return 0;
    }
    uint64_t s_block = find_superblock(rank);
    rank -= super_block_rank_[s_block];
    auto pos = find_basicblock_is(rank, s_block);
    rank -= basic_block_rank_[pos];
    pos *= kWordsPerBlock;

    // Final search
    if (pos + kWordsPerBlock - 1 < kWordsPerBlock) [[unlikely]] {
      size_t ones = std::popcount(bits_[pos]);
      while (pos < bits_.size() && ones < rank) {
        rank -= ones;
        ones = std::popcount(bits_[++pos]);
      }
      return kWordSize * pos + select_64(bits_[pos], rank - 1);
    }
    return kWordSize * pos + select_512(&bits_[pos], rank - 1);
  }

  /**
   * @brief Select the position of the rank0-th 0-bit (1-indexed).
   *
   * @param rank0 1-based rank of the 0-bit to select.
   * @return Bit index,
   * or size() if rank0 is out of range.
   */
  uint64_t select0(size_t rank0) const {
    if (rank0 > num_bits_ - max_rank_) [[unlikely]] {
      return num_bits_;
    }
    if (rank0 == 0) [[unlikely]] {
      return 0;
    }
    uint64_t s_block = find_superblock_zeros(rank0);
    rank0 -= kSuperBlockSize * s_block - super_block_rank_[s_block];
    auto pos = find_basicblock_is_zeros(rank0, s_block);
    auto pos_in_super_block = pos & (kBlocksPerSuperBlock - 1);
    rank0 -= kBasicBlockSize * pos_in_super_block - basic_block_rank_[pos];
    pos *= kWordsPerBlock;

    // Final search
    if (pos + kWordsPerBlock - 1 < kWordsPerBlock) [[unlikely]] {
      size_t zeros = std::popcount(~bits_[pos]);
      while (pos < bits_.size() && zeros < rank0) {
        rank0 -= zeros;
        zeros = std::popcount(~bits_[++pos]);
      }
      return kWordSize * pos + select_64(~bits_[pos], rank0 - 1);
    }
    return kWordSize * pos + select0_512(&bits_[pos], rank0 - 1);
  }

  /**
   * @brief Convert to a binary string (debug helper).
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
 * @brief Interleaved, owning bit vector with rank and select.
 *
 *
 * @details
 * This variant interleaves data with local rank metadata to reduce
 * cache
 * misses for rank queries. It copies input bits into an interleaved
 * layout.
 *
 * Based on:
 * "SPIDER: Improved Succinct Rank and Select
 * Performance"
 * Matthew D. Laws, Jocelyn Bliven, Kit Conklin, Elyes Laalai,
 * Samuel McCauley,
 * Zach S. Sturdevant
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
  std::vector<uint64_t> super_block_rank_;

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
      uint64_t next = iterator_64_ < bits_.size() ? bits_[iterator_64_++] : 0;
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
   * @brief Construct from an external array of 64-bit words.
   * @param
   * bit_vector Backing data to copy and interleave.
   * @param num_bits Number
   * of valid bits in the vector.
   */
  explicit BitVectorInterleaved(std::span<const uint64_t> bit_vector,
                                size_t num_bits)
      : num_bits_(std::min(num_bits, bit_vector.size() * kWordSize)) {
    build_rank_interleaved(bit_vector, num_bits);
  }

  /**
   * @brief Mask with the lowest num bits set.
   */
  static inline uint64_t first_bits_mask(size_t num) {
    return num >= 64 ? UINT64_MAX : ((1llu << num) - 1);
  }

  /**
   * @brief Returns the number of valid bits.
   */
  size_t size() const { return num_bits_; }

  /**
   * @brief Returns the bit at the given position.
   * @param pos Bit
   * index in [0, size()).
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
   * @brief Build the interleaved layout and rank index.
   * @param bits
   * Source bit vector as 64-bit words.
   * @param num_bits Number of valid
   * bits in the source.
   */
  void build_rank_interleaved(std::span<const uint64_t> bits, size_t num_bits) {
    size_t num_superblocks = 1 + (num_bits_ - 1) / kSuperBlockSize;
    super_block_rank_.resize(num_superblocks);
    size_t num_basicblocks = 1 + (num_bits_ - 1) / kBasicBlockSize;
    bits_interleaved.resize(num_basicblocks * (512 / kWordSize));

    uint64_t super_block_sum = 0;
    uint16_t basic_block_sum = 0;
    auto bit_reader = BitReader(bits);

    for (size_t i = 0; i * kBasicBlockSize < num_bits; ++i) {
      if (i % (kSuperBlockSize / kBasicBlockSize) == 0) {
        super_block_sum += basic_block_sum;
        super_block_rank_[i / (kSuperBlockSize / kBasicBlockSize)] =
            super_block_sum;
        basic_block_sum = 0;
      }
      bits_interleaved[i * (kWordsPerBlock) + 7] =
          static_cast<uint64_t>(basic_block_sum) << 48;

      for (size_t j = 0; j < 7 && kWordSize * (i + j) < num_bits; ++j) {
        bits_interleaved[i * (kWordsPerBlock) + j] =
            bit_reader.ReadBits64(std::min<uint64_t>(
                64ull, num_bits - i * kBasicBlockSize + j * kWordSize));
        basic_block_sum +=
            std::popcount(bits_interleaved[i * (kWordsPerBlock) + j]);
      }
      if ((i + 7) * kWordSize < num_bits) {
        auto v = bit_reader.ReadBits64(std::min<uint64_t>(
            48ull, num_bits - (i * kBasicBlockSize + 7 * kWordSize)));
        bits_interleaved[i * (kWordsPerBlock) + 7] ^= v;
        basic_block_sum += std::popcount(v);
      }
    }
  }

  /**
   * @brief Rank of 1s up to position pos (exclusive).
   * @param pos Bit
   * index in [0, size()].
   * @return Number of 1s in [0, pos).
   */
  uint64_t rank(size_t pos) const {
    // Multiplication/devisions
    uint64_t b_block = pos / kBasicBlockSize;
    uint64_t s_block = b_block / kBlocksPerSuperBlock;
    uint64_t b_block_pos = b_block * kWordsPerBlock;
    // Super block rank
    uint64_t result = super_block_rank_[s_block];
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
   * @brief Convert to a binary string (debug helper).
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

}  // namespace pixie
