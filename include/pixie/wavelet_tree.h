#pragma once

/**
 * @file wavelet_tree.h
 * @brief Common interface for wavelet-tree indexes.
 *
 * Include `<pixie/wavelet_tree/implementations.h>` to use Pixie's concrete
 * storage-backed wavelet trees.
 */

#include <cstddef>
#include <cstdint>
#include <vector>

namespace pixie {

/** @brief Construction strategy for a wavelet-tree implementation. */
enum class WaveletTreeBuildType { Standard, Huffman };

/**
 * @brief CRTP facade for wavelet-tree queries.
 *
 * @see `<pixie/wavelet_tree/implementations.h>` for the available concrete
 * implementations.
 */
template <class Impl>
class WaveletTreeBase {
 public:
  /**
   * @brief Return the number of symbols in the indexed sequence.
   * @return Logical sequence length.
   */
  std::size_t size() const { return impl().size_impl(); }

  /**
   * @brief Check whether the indexed sequence is empty.
   * @return `true` when `size() == 0`.
   */
  bool empty() const { return size() == 0; }

  /**
   * @brief Count occurrences of @p symbol in `[0, end_position)`.
   * @param symbol Symbol in the indexed alphabet.
   * @param end_position Prefix boundary in `[0, size()]`.
   * @return Number of occurrences, or zero for a symbol outside the alphabet.
   */
  std::size_t rank(std::uint64_t symbol, std::size_t end_position) const {
    return impl().rank_impl(symbol, end_position);
  }

  /**
   * @brief Return the position of the rank-th occurrence of @p symbol.
   * @param symbol Symbol in the indexed alphabet.
   * @param rank One-based occurrence rank.
   * @return Zero-based sequence position, or `size()` when absent.
   */
  std::size_t select(std::uint64_t symbol, std::size_t rank) const {
    return impl().select_impl(symbol, rank);
  }

  /**
   * @brief Reconstruct the sequence range `[@p begin, @p end)`.
   * @param begin First zero-based sequence position.
   * @param end One past the last sequence position; must not exceed `size()`.
   * @return Symbols in the requested half-open range.
   */
  std::vector<std::uint64_t> get_segment(std::size_t begin,
                                         std::size_t end) const {
    return impl().get_segment_impl(begin, end);
  }

 private:
  /** @brief Return this facade as its concrete CRTP implementation. */
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
