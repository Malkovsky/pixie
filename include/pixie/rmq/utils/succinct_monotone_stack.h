#pragma once

/**
 * Monotone-stack construction experiment, 2026-06-13.
 *
 * SIMD investigation: for the RMQ benchmark's uniform-random value arrays, the
 * Cartesian construction pop-run distribution is short: roughly 50% of
 * iterations pop 0 entries, 25% pop 1, 12.5% pop 2, and only about 6.25% pop
 * 4 or more. A SIMD path that gathers several stack-top values would therefore
 * pay gather/setup overhead on mostly short runs. The retained improvement is
 * bit-parallel rather than SIMD: 64-bit key blocks plus compact non-empty block
 * summaries replace the old 63-bit data/pointer-word layout.
 *
 * Command shape:
 *   taskset -c 0 ./build/release/bench_rmq
 *     --benchmark_filter='^rmq_build_(cartesian_rmm|cartesian_hybrid_btree|sdsl_sct)/(4194304|67108864)$'
 *     --benchmark_repetitions=5
 *
 * CPU mean, milliseconds.
 *
 * | N    | row                      | before | after   |
 * | ---: | ------------------------ | -----: | ------: |
 * | 2^22 | CartesianRmM             | 44.731 |  40.274 |
 * | 2^22 | CartesianHybridBTree     | 47.943 |  45.019 |
 * | 2^22 | SdslSct control          | 39.911 |  39.628 |
 * | ---- | ------------------------ | ------ | ------- |
 * | 2^26 | CartesianRmM             | 750.316 | 678.925 |
 * | 2^26 | CartesianHybridBTree     | 772.157 | 732.582 |
 * | 2^26 | SdslSct control          | 642.540 | 644.627 |
 */

#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace pixie::rmq::utils {

/**
 * @brief Construction-only stack for monotonically increasing integer keys.
 *
 * @details The stack stores the currently live keys in 64-bit blocks and keeps
 * compact summary bitsets for non-empty blocks. It is intended for
 * Cartesian-tree construction passes where pushed keys are strictly increasing.
 * Under that precondition the logical stack top is the largest stored key.
 *
 * The helper keeps construction memory succinct in the worst case: an input of
 * `capacity` keys uses one bit per key plus about 1/64 extra summary bits,
 * instead of one full index per live stack entry. It does not own or inspect
 * input values and deliberately does not define the Cartesian-tree tie rule;
 * callers keep their existing pop predicate.
 */
class SuccinctIncreasingStack {
 public:
  explicit SuccinctIncreasingStack(std::size_t capacity)
      : words_(word_count(capacity)),
        nonempty_words_(word_count(words_.size())),
        nonempty_super_words_(word_count(nonempty_words_.size())),
        capacity_(capacity) {}

  /**
   * @brief Whether the logical stack is empty.
   */
  bool empty() const { return size_ == 0; }

  /**
   * @brief Return the current top key.
   *
   * @pre The stack is not empty.
   */
  std::size_t top() const {
    assert(!empty());
    return top_key_;
  }

  /**
   * @brief Push a key larger than the current top key.
   *
   * @details Keys are expected to be pushed in strictly increasing order over
   * the lifetime of this stack. The assertion checks the live top, which is the
   * only ordering constraint needed by this bitset representation.
   */
  void push(std::size_t key) {
    assert(key < capacity_);
    assert(empty() || top() < key);

    const std::size_t block = block_index(key);
    words_[block] |= std::uint64_t{1} << block_position(key);
    set_nonempty(block);
    top_key_ = key;
    ++size_;
  }

  /**
   * @brief Remove the current top key.
   *
   * @pre The stack is not empty.
   */
  void pop() {
    assert(!empty());

    const std::size_t block = block_index(top_key_);
    std::uint64_t word = words_[block];
    word &= ~(std::uint64_t{1} << block_position(top_key_));
    words_[block] = word;

    --size_;
    if (size_ == 0) {
      if (word == 0) {
        clear_nonempty(block);
      }
      top_key_ = 0;
      return;
    }

    if (word != 0) {
      top_key_ = block * kBlockBits + highest_bit(word);
      return;
    }

    clear_nonempty(block);
    const std::size_t previous_block = previous_nonempty_block(block);
    const std::uint64_t previous_word = words_[previous_block];
    assert(previous_word != 0);
    top_key_ = previous_block * kBlockBits + highest_bit(previous_word);
  }

 private:
  static constexpr std::size_t kBlockBits = 64;
  static constexpr std::size_t kBlockShift = 6;
  static constexpr std::size_t kBlockMask = kBlockBits - 1;

  static std::size_t word_count(std::size_t bit_count) {
    return (bit_count + kBlockBits - 1) / kBlockBits;
  }

  static std::size_t block_index(std::size_t key) { return key >> kBlockShift; }

  static std::size_t block_position(std::size_t key) {
    return key & kBlockMask;
  }

  static std::size_t highest_bit(std::uint64_t word) {
    assert(word != 0);
    return static_cast<std::size_t>(std::bit_width(word) - 1);
  }

  static std::uint64_t low_bits(std::size_t count) {
    return count == 0 ? 0 : ((std::uint64_t{1} << count) - 1);
  }

  void set_nonempty(std::size_t block) {
    const std::size_t summary = block_index(block);
    const std::uint64_t bit = std::uint64_t{1} << block_position(block);
    if ((nonempty_words_[summary] & bit) == 0) {
      nonempty_words_[summary] |= bit;
      nonempty_super_words_[block_index(summary)] |= std::uint64_t{1}
                                                     << block_position(summary);
    }
  }

  void clear_nonempty(std::size_t block) {
    const std::size_t summary = block_index(block);
    nonempty_words_[summary] &= ~(std::uint64_t{1} << block_position(block));
    if (nonempty_words_[summary] == 0) {
      nonempty_super_words_[block_index(summary)] &=
          ~(std::uint64_t{1} << block_position(summary));
    }
  }

  std::size_t previous_nonempty_block(std::size_t block) const {
    assert(block > 0);
    std::size_t summary = block_index(block);
    const std::size_t summary_position = block_position(block);
    if (summary_position != 0) {
      const std::uint64_t same_summary =
          nonempty_words_[summary] & low_bits(summary_position);
      if (same_summary != 0) {
        return summary * kBlockBits + highest_bit(same_summary);
      }
    }

    std::size_t super = block_index(summary);
    const std::size_t super_position = block_position(summary);
    std::uint64_t super_mask =
        super_position == 0
            ? 0
            : nonempty_super_words_[super] & low_bits(super_position);
    while (super_mask == 0) {
      assert(super > 0);
      --super;
      super_mask = nonempty_super_words_[super];
    }

    const std::size_t previous_summary =
        super * kBlockBits + highest_bit(super_mask);
    assert(previous_summary < nonempty_words_.size());
    const std::uint64_t previous_summary_word =
        nonempty_words_[previous_summary];
    assert(previous_summary_word != 0);
    return previous_summary * kBlockBits + highest_bit(previous_summary_word);
  }

  std::vector<std::uint64_t> words_;
  std::vector<std::uint64_t> nonempty_words_;
  std::vector<std::uint64_t> nonempty_super_words_;
  std::size_t capacity_ = 0;
  std::size_t size_ = 0;
  std::size_t top_key_ = 0;
};

}  // namespace pixie::rmq::utils
