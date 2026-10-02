#pragma once

/**
 * @file byte_histogram.h
 * @brief Dependency-reduced histogram construction for byte streams.
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace pixie::detail {

/** @brief Incremental 256-bin byte histogram with four update lanes. */
class ByteHistogram {
 public:
  /** @brief Add one contiguous chunk to the histogram. */
  void add(std::span<const std::byte> bytes) {
    std::size_t index = 0;
    for (; index + kLaneCount <= bytes.size(); index += kLaneCount) {
      ++partial_[0][std::to_integer<std::uint8_t>(bytes[index])];
      ++partial_[1][std::to_integer<std::uint8_t>(bytes[index + 1])];
      ++partial_[2][std::to_integer<std::uint8_t>(bytes[index + 2])];
      ++partial_[3][std::to_integer<std::uint8_t>(bytes[index + 3])];
    }
    for (; index < bytes.size(); ++index) {
      ++partial_[index & (kLaneCount - 1)]
                [std::to_integer<std::uint8_t>(bytes[index])];
    }
  }

  /** @brief Return the combined symbol counts. */
  std::array<std::size_t, 256> counts() const {
    std::array<std::size_t, 256> result{};
    for (std::size_t symbol = 0; symbol < result.size(); ++symbol) {
      for (const auto& lane : partial_) {
        result[symbol] += lane[symbol];
      }
    }
    return result;
  }

 private:
  static constexpr std::size_t kLaneCount = 4;
  static_assert((kLaneCount & (kLaneCount - 1)) == 0);

  alignas(64) std::array<std::array<std::size_t, 256>, kLaneCount> partial_{};
};

}  // namespace pixie::detail
