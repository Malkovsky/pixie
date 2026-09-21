#include <gtest/gtest.h>
#include <pixie/detail/sequence/sequence_tree_kernels.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <span>

namespace {

using pixie::detail::sequence::node_select;
using pixie::detail::sequence::node_select_scalar;

constexpr auto kMax = std::numeric_limits<std::size_t>::max();
constexpr auto kTop = std::size_t{1}
                      << (std::numeric_limits<std::size_t>::digits - 1);

void CheckSelection(std::span<const std::size_t> ends, std::size_t index) {
  const auto expected =
      ends.empty() ? 0
                   : static_cast<std::size_t>(
                         std::upper_bound(ends.begin(), ends.end(), index) -
                         ends.begin());
  EXPECT_EQ(node_select_scalar(ends, index), expected)
      << "count=" << ends.size() << " index=" << index;
  EXPECT_EQ(node_select(ends, index), expected)
      << "count=" << ends.size() << " index=" << index;
}

TEST(SequenceTreeKernels, EmptyNullBacking) {
  const std::span<const std::size_t> ends(
      static_cast<const std::size_t*>(nullptr), std::size_t{0});
  for (const auto index :
       {std::size_t{0}, std::size_t{1}, kTop - 1, kTop, kMax}) {
    CheckSelection(ends, index);
  }
}

TEST(SequenceTreeKernels, ExactExtentCountsAndAlignments) {
  for (std::size_t count = 1; count <= 15; ++count) {
    for (std::size_t offset = 0; offset < 8; ++offset) {
      // Every tested span ends at the allocation boundary: ASan can detect
      // unmasked vector/tail overreads even when the starting address is
      // skewed.
      auto storage = std::make_unique<std::size_t[]>(offset + count);
      std::span<std::size_t> ends(storage.get() + offset, count);
      for (const auto base : {std::size_t{0}, kTop - 32, kMax - 64}) {
        for (std::size_t i = 0; i < count; ++i) {
          ends[i] = base + 4 * i;
        }
        for (std::size_t delta = 0; delta <= 64; ++delta) {
          CheckSelection(ends, base + delta);
        }
        CheckSelection(ends, 0);
        CheckSelection(ends, kMax);
      }
      for (const auto value : {std::size_t{0}, kTop, kMax}) {
        std::fill(ends.begin(), ends.end(), value);
        CheckSelection(ends, value);
        if (value != 0) {
          CheckSelection(ends, value - 1);
        }
        if (value != kMax) {
          CheckSelection(ends, value + 1);
        }
      }
    }
  }
}

TEST(SequenceTreeKernels, ExhaustiveShortMonotoneEnds) {
  for (unsigned mask = 0; mask < 256; ++mask) {
    const auto count = static_cast<std::size_t>(std::popcount(mask));
    auto storage = count ? std::make_unique<std::size_t[]>(count) : nullptr;
    std::span<std::size_t> ends(storage.get(), count);
    for (const auto base : {std::size_t{0}, kTop - 4, kMax - 7}) {
      std::size_t i = 0;
      for (unsigned bit = 0; bit < 8; ++bit) {
        if ((mask >> bit) & 1u) {
          ends[i++] = base + bit;
        }
      }
      for (std::size_t delta = 0; delta < 8; ++delta) {
        CheckSelection(ends, base + delta);
      }
      CheckSelection(ends, 0);
      CheckSelection(ends, kMax);
    }
  }
}

TEST(SequenceTreeKernels, RandomFullWidthAndRepeatedEnds) {
  std::mt19937_64 random(0x5e1ec7);
  for (std::size_t count = 1; count <= 15; ++count) {
    for (std::size_t trial = 0; trial < 128; ++trial) {
      const auto offset = trial % 8;
      auto storage = std::make_unique<std::size_t[]>(offset + count);
      std::span<std::size_t> ends(storage.get() + offset, count);
      for (auto& end : ends) {
        switch (trial % 4) {
          case 0:
            end = static_cast<std::size_t>(random());
            break;
          case 1:
            end = kTop - 16 + random() % 33;
            break;
          case 2:
            end = kMax - random() % 64;
            break;
          default:
            end = random() % 8;
            break;
        }
      }
      std::sort(ends.begin(), ends.end());
      for (const auto end : ends) {
        CheckSelection(ends, end);
        if (end != 0) {
          CheckSelection(ends, end - 1);
        }
        if (end != kMax) {
          CheckSelection(ends, end + 1);
        }
      }
      CheckSelection(ends, 0);
      CheckSelection(ends, kMax);
      for (unsigned i = 0; i < 16; ++i) {
        CheckSelection(ends, static_cast<std::size_t>(random()));
      }
    }
  }
}

}  // namespace
