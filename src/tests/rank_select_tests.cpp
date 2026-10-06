#include <gtest/gtest.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/rank_select/compressed.h>
#include <pixie/rank_select/support.h>
#ifdef PIXIE_PASTA_SUPPORT
#include <pixie/rank_select/pasta.h>
#endif

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace {

void overwrite_u64(std::vector<std::byte>& bytes,
                   std::size_t offset,
                   std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
  }
}

void overwrite_u32(std::vector<std::byte>& bytes,
                   std::size_t offset,
                   std::uint32_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    bytes[offset + byte] =
        static_cast<std::byte>((value >> (byte * 8)) & 0xffu);
  }
}

struct RankSelectSampleLayout {
  std::size_t select1_begin;
  std::size_t select1_count;
  std::size_t select0_begin;
  std::size_t select0_count;
  std::size_t super_bytes;
  std::size_t basic_bytes;
  std::size_t sample_bytes;
  std::size_t storage_offset;
};

RankSelectSampleLayout locate_rank_select_samples(
    std::span<const std::byte> artifact) {
  pixie::BinaryReader reader(artifact);
  reader.skip(3 * sizeof(std::uint64_t));
  RankSelectSampleLayout result;
  result.select1_begin = reader.read_size();
  result.select1_count = reader.read_size();
  result.select0_begin = reader.read_size();
  result.select0_count = reader.read_size();
  reader.skip(2 * sizeof(std::uint32_t));
  result.super_bytes = reader.read_size();
  reader.skip(result.super_bytes);
  result.basic_bytes = reader.read_size();
  reader.skip(result.basic_bytes);
  result.sample_bytes = reader.read_size();
  result.storage_offset = reader.position();
  reader.skip(result.sample_bytes);
  return result;
}

template <class Support>
class RankSelectSpecificationTest : public testing::Test {};

using RankSelectImplementations =
    testing::Types<pixie::RankSelectSupport<>,
                   pixie::RankSelectSupportCompressed<>,
                   pixie::RankSelectSupportCompressed<2048>
#ifdef PIXIE_PASTA_SUPPORT
                   ,
                   pixie::PastaRankSelectSupport
#endif
                   >;
TYPED_TEST_SUITE(RankSelectSpecificationTest, RankSelectImplementations);

TYPED_TEST(RankSelectSpecificationTest, EmptySequence) {
  const std::vector<std::uint64_t> words;
  const TypeParam bits(words, 0);

  EXPECT_TRUE(bits.empty());
  EXPECT_EQ(bits.size(), 0u);
  EXPECT_EQ(bits.rank(0), 0u);
  EXPECT_EQ(bits.rank0(0), 0u);
  EXPECT_EQ(bits.select(1), 0u);
  EXPECT_EQ(bits.select0(1), 0u);
  EXPECT_EQ(bits.to_string(), "");
}

TYPED_TEST(RankSelectSpecificationTest, RankSelectAndAccessMatchPackedBits) {
  constexpr std::size_t sizes[] = {1,    2,    63,    64,    65,    127,
                                   495,  496,  497,   511,   512,   513,
                                   1024, 4097, 65535, 65536, 65537, 131073};
  std::mt19937_64 rng(42);

  for (const std::size_t size : sizes) {
    std::vector<std::uint64_t> words((size + 63) / 64);
    for (auto& word : words) {
      word = rng();
    }
    const TypeParam bits(words, size);
    ASSERT_FALSE(bits.empty());
    ASSERT_EQ(bits.size(), size);
    ASSERT_TRUE(bits.supports_select1());
    ASSERT_TRUE(bits.supports_select0());

    std::size_t ones = 0;
    std::size_t zeros = 0;
    std::string expected;
    expected.reserve(size);
    for (std::size_t position = 0; position < size; ++position) {
      const bool value = ((words[position / 64] >> (position % 64)) & 1u) != 0;
      EXPECT_EQ(bits.rank(position), ones) << "size=" << size;
      EXPECT_EQ(bits.rank0(position), zeros) << "size=" << size;
      EXPECT_EQ(bits[position], value) << "size=" << size;
      expected.push_back(value ? '1' : '0');
      if (value) {
        ++ones;
        EXPECT_EQ(bits.select(ones), position) << "size=" << size;
      } else {
        ++zeros;
        EXPECT_EQ(bits.select0(zeros), position) << "size=" << size;
      }
    }
    EXPECT_EQ(bits.rank(size), ones);
    EXPECT_EQ(bits.rank0(size), zeros);
    EXPECT_EQ(bits.select(ones + 1), size);
    EXPECT_EQ(bits.select0(zeros + 1), size);
    EXPECT_EQ(bits.to_string(), expected);
  }
}

TEST(CompressedRankSelectTest, ExhaustiveShortPatterns) {
  using Support = pixie::RankSelectSupportCompressed<8, 4>;
  for (std::size_t size = 0; size <= 10; ++size) {
    for (std::uint64_t pattern = 0; pattern < (std::uint64_t{1} << size);
         ++pattern) {
      // Keep dirty padding to exercise the logical final-word mask.
      const std::array words{pattern | ~::first_bits_mask(size)};
      const Support support(words, size);
      std::size_t ones = 0;
      std::size_t zeros = 0;
      EXPECT_EQ(support.select(0), 0u);
      EXPECT_EQ(support.select0(0), 0u);
      for (std::size_t position = 0; position < size; ++position) {
        EXPECT_EQ(support.rank(position), ones);
        EXPECT_EQ(support.rank0(position), zeros);
        if ((pattern >> position) & 1) {
          EXPECT_EQ(support.select(++ones), position);
        } else {
          EXPECT_EQ(support.select0(++zeros), position);
        }
      }
      EXPECT_EQ(support.rank(size + 1), ones);
      EXPECT_EQ(support.select(ones + 1), size);
      EXPECT_EQ(support.select0(zeros + 1), size);
    }
  }
}

TYPED_TEST(RankSelectSpecificationTest,
           HomogeneousSuperblocksAndTerminalLocalRank) {
  using Support = TypeParam;
  for (const std::size_t size : {65536u, 65553u}) {
    for (const std::uint64_t value : {std::uint64_t{0}, ~std::uint64_t{0}}) {
      const std::vector<std::uint64_t> words((size + 63) / 64, value);
      const Support support(words, size);
      for (std::size_t position = 0; position <= size; ++position) {
        EXPECT_EQ(support.rank(position), value == 0 ? 0 : position);
        EXPECT_EQ(
            value == 0 ? support.select0(position) : support.select(position),
            position == 0 ? 0 : position - 1);
      }
      EXPECT_EQ(value == 0 ? support.select(1) : support.select0(1), size);
    }
  }
}

TEST(CompressedRankSelectTest,
     WideIntervalsSkipEmptySuperblocksInBothDirections) {
  using Support = pixie::RankSelectSupportCompressed<8, 4>;
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 45 * kSuper + 13;
  const std::array<std::size_t, 18> positions = {
      5 * kSuper + 1,  5 * kSuper + 2,  5 * kSuper + 3,  5 * kSuper + 4,
      9 * kSuper + 1,  9 * kSuper + 2,  20 * kSuper + 1, 20 * kSuper + 2,
      20 * kSuper + 3, 21 * kSuper + 1, 21 * kSuper + 2, 21 * kSuper + 3,
      30 * kSuper + 1, 30 * kSuper + 2, 30 * kSuper + 3, 40 * kSuper + 1,
      40 * kSuper + 2, 40 * kSuper + 3};
  for (const bool invert : {false, true}) {
    std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
    for (const auto position : positions) {
      words[position / 64] |= std::uint64_t{1} << (position % 64);
    }
    if (invert) {
      for (auto& word : words) {
        word = ~word;
      }
    }
    const Support support(words, kSize);
    const pixie::RankSelectSupport<> reference(words, kSize);
    for (std::size_t i = 0; i < positions.size(); ++i) {
      EXPECT_EQ(invert ? support.select0(i + 1) : support.select(i + 1),
                positions[i]);
    }
    EXPECT_EQ(invert ? support.select0(positions.size() + 1)
                     : support.select(positions.size() + 1),
              kSize);
    std::mt19937_64 rng(42);
    for (std::size_t i = 0; i < 10000; ++i) {
      const auto position = rng() % (kSize + 1);
      EXPECT_EQ(support.rank(position), reference.rank(position));
      const auto rank = 1 + rng() % (kSize - positions.size());
      EXPECT_EQ(invert ? support.select(rank) : support.select0(rank),
                invert ? reference.select(rank) : reference.select0(rank));
    }
  }
}

TEST(CompressedRankSelectTest, ExactThresholdAndSampleGroupBoundaries) {
  using Support = pixie::RankSelectSupportCompressed<8, 4>;
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 150 * kSuper + 37;
  std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
  std::vector<std::size_t> positions;
  // Eight occurrences every four superblocks: every coarse gap meets the
  // threshold exactly, and >32 checkpoints cross the packed-group boundary.
  for (std::size_t group = 0; group < 36; ++group) {
    for (std::size_t i = 0; i < 8; ++i) {
      const auto position = group * 4 * kSuper + i * 67;
      words[position / 64] |= std::uint64_t{1} << (position % 64);
      positions.push_back(position);
    }
  }
  const Support support(words, kSize);
  for (std::size_t i = 0; i < positions.size(); ++i) {
    EXPECT_EQ(support.select(i + 1), positions[i]);
    EXPECT_EQ(support.rank(positions[i]), i);
    EXPECT_EQ(support.rank(positions[i] + 1), i + 1);
  }
  EXPECT_EQ(support.select(positions.size() + 1), kSize);
}

TYPED_TEST(RankSelectSpecificationTest,
           LocalCorrectionHandlesConcentratedPopulations) {
  constexpr std::size_t kSize = 16384;
  for (const bool last : {false, true}) {
    for (const bool invert : {false, true}) {
      std::vector<std::uint64_t> words(kSize / 64, 0);
      const std::size_t first_word = last ? words.size() - 8 : 0;
      std::fill_n(words.begin() + first_word, 8, ~std::uint64_t{0});
      if (invert) {
        for (auto& word : words) {
          word = ~word;
        }
      }
      const TypeParam support(words, kSize);
      std::size_t ones = 0;
      std::size_t zeros = 0;
      for (std::size_t position = 0; position < kSize; ++position) {
        EXPECT_EQ(support.rank(position), ones);
        if ((words[position / 64] >> (position % 64)) & 1) {
          EXPECT_EQ(support.select(++ones), position);
        } else {
          EXPECT_EQ(support.select0(++zeros), position);
        }
      }
    }
  }
}

TEST(CompressedRankSelectTest, ModesClampingAndCopyMove) {
  using Support = pixie::RankSelectSupportCompressed<>;
  using Mode = Support::SelectSupport;
  const std::array<std::uint64_t, 1> words{0b101101};
  const Support empty;
  EXPECT_EQ(empty.rank(0), 0u);
  EXPECT_FALSE(empty.supports_select1());
  for (const auto mode :
       {Mode::kNone, Mode::kSelect1, Mode::kSelect0, Mode::kBoth}) {
    const Support support(words, 6, mode);
    EXPECT_EQ(support.rank(100), 4u);
    EXPECT_EQ(support.rank0(100), 2u);
    EXPECT_EQ(support.select(0), 0u);
    EXPECT_EQ(support.select0(0), 0u);
    EXPECT_EQ(support.select(4), support.supports_select1() ? 5u : 6u);
    EXPECT_EQ(support.select0(2), support.supports_select0() ? 4u : 6u);
    Support copy = support;
    const Support moved = std::move(copy);
    EXPECT_EQ(moved.select(4), support.select(4));
    EXPECT_EQ(moved.select0(2), support.select0(2));
    EXPECT_TRUE(copy.empty());
    EXPECT_FALSE(copy.supports_select1());
    EXPECT_FALSE(copy.supports_select0());
    EXPECT_EQ(copy.rank(1), 0u);
    EXPECT_EQ(copy.rank0(1), 0u);
    EXPECT_EQ(copy.select(1), 0u);
    EXPECT_EQ(copy.select0(1), 0u);
  }
  const Support clamped(words, 65);
  EXPECT_EQ(clamped.size(), 64u);
  EXPECT_THROW((Support(words, 6, static_cast<Mode>(4))),
               std::invalid_argument);
}

TEST(CompressedRankSelectTest,
     ModesClampingAndCopyMovePreservePositionSamples) {
  const auto check = []<class Support>() {
    using Mode = typename Support::SelectSupport;
    constexpr std::size_t kSize = 65536 + 19;
    std::vector<std::uint64_t> words((kSize + 63) / 64);
    std::mt19937_64 random(42);
    for (auto& word : words) {
      word = random();
    }
    const pixie::RankSelectSupport<> reference(words, kSize);
    for (const auto mode :
         {Mode::kNone, Mode::kSelect1, Mode::kSelect0, Mode::kBoth}) {
      const Support original(words, kSize, mode);
      Support copy(original);
      Support assigned;
      assigned = copy;
      Support moved(std::move(copy));
      Support move_assigned;
      move_assigned = std::move(assigned);
      for (const Support* support :
           std::array<const Support*, 3>{&original, &moved, &move_assigned}) {
        EXPECT_EQ(support->rank(kSize + 100), reference.rank(kSize));
        EXPECT_EQ(support->rank0(kSize + 100), reference.rank0(kSize));
        EXPECT_EQ(support->select(0), 0u);
        EXPECT_EQ(support->select0(0), 0u);
        for (std::size_t rank = 1; rank <= reference.rank(kSize); rank += 137) {
          EXPECT_EQ(support->select(rank), support->supports_select1()
                                               ? reference.select(rank)
                                               : kSize);
        }
        for (std::size_t rank = 1; rank <= reference.rank0(kSize);
             rank += 137) {
          EXPECT_EQ(support->select0(rank), support->supports_select0()
                                                ? reference.select0(rank)
                                                : kSize);
        }
      }
      for (const Support* empty : {&copy, &assigned}) {
        EXPECT_TRUE(empty->empty());
        EXPECT_FALSE(empty->supports_select1());
        EXPECT_FALSE(empty->supports_select0());
        EXPECT_EQ(empty->select(1), 0u);
        EXPECT_EQ(empty->select0(1), 0u);
      }
    }
    const Support clamped(words, words.size() * 64 + 1);
    EXPECT_EQ(clamped.size(), words.size() * 64);
    EXPECT_THROW((Support(words, kSize, static_cast<Mode>(4))),
                 std::invalid_argument);
  };
  check.template operator()<pixie::RankSelectSupportCompressed<>>();
}

TEST(CompressedRankSelectTest, ExactSamplesGroupEdgesGapsAndPartialTerminal) {
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 100 * kSuper + 37;
  std::mt19937_64 random(731);
  for (int shape = 0; shape < 3; ++shape) {
    std::vector<std::uint64_t> words((kSize + 63) / 64);
    std::array<std::vector<std::uint64_t>, 2> positions;
    for (std::size_t p = 0; p < kSize; ++p) {
      const auto super = p / kSuper;
      const bool bit = shape == 0 ? (random() & 1)
                       : shape == 1
                           ? ((super < 32 || super % 32 == 0) && (random() & 1))
                           : (super % 7 < 3);
      words[p / 64] |= std::uint64_t{bit} << (p % 64);
      positions[!bit].push_back(p);
    }
    // Ignore dirty source padding, including during final-sample construction.
    words.back() |= ~::first_bits_mask(kSize % 64);
    const auto check = [&]<std::size_t Frequency>() {
      const pixie::RankSelectSupportCompressed<Frequency> index(words, kSize);
      for (bool zero : {false, true}) {
        for (std::size_t k = 1; k <= positions[zero].size(); ++k) {
          const auto result = zero ? index.select0(k) : index.select(k);
          ASSERT_EQ(result, positions[zero][k - 1])
              << "shape=" << shape << " frequency=" << Frequency
              << " zero=" << zero << " rank=" << k;
        }
        EXPECT_EQ(zero ? index.select0(positions[zero].size() + 1)
                       : index.select(positions[zero].size() + 1),
                  kSize);
      }
    };
    check.template operator()<8>();
    check.template operator()<4096>();
    check.template operator()<65536>();
  }
}

TEST(CompressedRankSelectTest, CopyMoveAssignmentPreservesMixedDirectories) {
  using Support = pixie::RankSelectSupportCompressed<>;
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 257 * kSuper + 17;
  std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
  for (std::size_t super = 0; super < 257; ++super) {
    if (super % 65 == 0) {
      std::fill_n(words.begin() + super * kSuper / 64, 128, ~std::uint64_t{0});
    } else {
      words[super * kSuper / 64] = 0xffff;
    }
  }
  words.back() =
      ~std::uint64_t{0};  // Ignore dirty padding after 17 valid bits.
  const Support reference(words, kSize);
  Support copy;
  {
    const Support original(words, kSize);
    copy = original;
  }
  const auto check = [&](const Support& support) {
    EXPECT_EQ(support.size(), reference.size());
    EXPECT_EQ(support.supports_select1(), reference.supports_select1());
    EXPECT_EQ(support.supports_select0(), reference.supports_select0());
    for (const auto position : {std::size_t{0}, kSuper, 65 * kSuper, kSize,
                                std::numeric_limits<std::size_t>::max()}) {
      EXPECT_EQ(support.rank(position), reference.rank(position));
      EXPECT_EQ(support.rank0(position), reference.rank0(position));
    }
    for (const bool zero : {false, true}) {
      const auto total = zero ? reference.rank0(kSize) : reference.rank(kSize);
      for (const auto rank : {std::size_t{0}, std::size_t{1}, std::size_t{4096},
                              std::size_t{4097}, total, total + 1}) {
        EXPECT_EQ(zero ? support.select0(rank) : support.select(rank),
                  zero ? reference.select0(rank) : reference.select(rank));
      }
    }
  };
  check(copy);
  auto& self = copy;
  copy = self;
  check(copy);
  copy = std::move(self);
  check(copy);
  Support moved(words, 1024, Support::SelectSupport::kNone);
  moved = std::move(copy);
  check(moved);
  EXPECT_TRUE(copy.empty());
  EXPECT_FALSE(copy.supports_select1());
  EXPECT_FALSE(copy.supports_select0());
  EXPECT_EQ(copy.to_string(), "");
  EXPECT_EQ(copy.rank(1), 0u);
  EXPECT_EQ(copy.rank0(1), 0u);
  EXPECT_EQ(copy.select(1), 0u);
  EXPECT_EQ(copy.select0(1), 0u);
  copy = moved;
  check(copy);
}

TEST(CompressedRankSelectTest, ContiguousAndGappedWideIntervals) {
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 600 * kSuper + 17;
  for (const bool gaps : {false, true}) {
    for (const bool invert : {false, true}) {
      std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
      std::vector<std::size_t> positions;
      for (std::size_t super = 0; super < 600; ++super) {
        // The first wide interval is contiguous. Later ones may contain gaps.
        if (gaps && super >= 300 && super % 5 == 0) {
          continue;
        }
        words[super * kSuper / 64] = 0xffff;
        for (std::size_t bit = 0; bit < 16; ++bit) {
          positions.push_back(super * kSuper + bit);
        }
      }
      words.back() = ~std::uint64_t{0};
      for (std::size_t bit = 0; bit < 17; ++bit) {
        positions.push_back(600 * kSuper + bit);
      }
      if (invert) {
        for (auto& word : words) {
          word = ~word;
        }
      }
      const pixie::RankSelectSupportCompressed<> support(words, kSize);
      const pixie::RankSelectSupport<> reference(words, kSize);
      for (std::size_t i = 0; i < positions.size(); ++i) {
        EXPECT_EQ(invert ? support.select0(i + 1) : support.select(i + 1),
                  positions[i]);
        EXPECT_EQ(support.rank(positions[i]), invert ? positions[i] - i : i);
      }
      EXPECT_EQ(invert ? support.select0(positions.size() + 1)
                       : support.select(positions.size() + 1),
                kSize);
      const auto total =
          invert ? reference.rank(kSize) : reference.rank0(kSize);
      for (const auto rank : {std::size_t{1}, std::size_t{4096},
                              std::size_t{4097}, total, total + 1}) {
        EXPECT_EQ(invert ? support.select(rank) : support.select0(rank),
                  invert ? reference.select(rank) : reference.select0(rank));
      }
    }
  }
}

TYPED_TEST(RankSelectSpecificationTest, DefaultSamplingOnSparsePartialSource) {
  constexpr std::size_t kSize = (1 << 22) + 17;
  std::mt19937_64 rng(42);
  for (const bool invert : {false, true}) {
    std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
    std::vector<std::size_t> positions;
    for (std::size_t position = 0; position < kSize; ++position) {
      if (rng() % 1024 == 0 || position == kSize - 1) {
        words[position / 64] |= std::uint64_t{1} << (position % 64);
        positions.push_back(position);
      }
    }
    if (invert) {
      for (auto& word : words) {
        word = ~word;
      }
    }
    const TypeParam support(words, kSize);
    for (std::size_t i = 0; i < positions.size(); ++i) {
      EXPECT_EQ(invert ? support.select0(i + 1) : support.select(i + 1),
                positions[i]);
      const auto expected = invert ? positions[i] - i : i;
      EXPECT_EQ(support.rank(positions[i]), expected);
    }
    EXPECT_EQ(invert ? support.select0(positions.size() + 1)
                     : support.select(positions.size() + 1),
              kSize);
  }
}

TEST(CompressedRankSelectTest, BitmapLookupOnNoncontiguousMixedSource) {
  constexpr std::size_t kBits = 521 * 16384 + 19;
  std::vector<std::uint64_t> words((kBits + 63) / 64);
  std::mt19937_64 random(81);
  for (std::size_t i = 0; i < words.size(); ++i) {
    const auto super = i / 256;
    words[i] = super % 65 == 0 ? random()
               : super % 5 == 0
                   ? 0
                   : (i % 16 == 0 ? std::uint64_t{1} << (random() % 64) : 0);
  }
  for (bool complement : {false, true}) {
    if (complement) {
      for (auto& word : words) {
        word = ~word;
      }
    }
    const pixie::RankSelectSupport<> reference(words, kBits);
    const pixie::RankSelectSupportCompressed<> support(words, kBits);
    for (unsigned trial = 0; trial < 20000; ++trial) {
      const auto position = random() % (kBits + 1);
      const auto one = 1 + random() % reference.rank(kBits);
      const auto zero = 1 + random() % reference.rank0(kBits);
      EXPECT_EQ(support.rank(position), reference.rank(position));
      EXPECT_EQ(support.rank0(position), reference.rank0(position));
      EXPECT_EQ(support.select(one), reference.select(one));
      EXPECT_EQ(support.select0(zero), reference.select0(zero));
    }
  }
}

template <class Samples>
class MonotoneSampleEncodingTest : public testing::Test {};
using MonotoneSampleImplementations =
    testing::Types<pixie::detail::compressed::BlockedMonotoneSamples,
                   pixie::detail::compressed::PositionSamples>;
TYPED_TEST_SUITE(MonotoneSampleEncodingTest, MonotoneSampleImplementations);

TYPED_TEST(MonotoneSampleEncodingTest,
           RepeatedValuesWideOffsetsAndPartialGroups) {
  for (const auto step : {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{32},
                          std::uint64_t{1} << 40}) {
    std::vector<std::uint64_t> values;
    for (std::size_t i = 0; i < 99; ++i) {
      values.push_back((std::uint64_t{1} << 48) + (i / 3) * step);
    }
    const TypeParam samples(values);
    for (std::size_t i = 0; i < values.size(); ++i) {
      EXPECT_EQ(samples[i], values[i]);
      if constexpr (requires { samples.bounds(i); }) {
        if (i + 1 < values.size()) {
          EXPECT_EQ(samples.bounds(i),
                    std::make_pair(values[i], values[i + 1]));
        }
      }
    }
  }
  const std::array<std::uint64_t, 3> extremes{
      0, 0, std::numeric_limits<std::uint64_t>::max()};
  const TypeParam samples(extremes);
  for (std::size_t i = 0; i < extremes.size(); ++i) {
    EXPECT_EQ(samples[i], extremes[i]);
    if constexpr (requires { samples.bounds(i); }) {
      if (i + 1 < extremes.size()) {
        EXPECT_EQ(samples.bounds(i),
                  std::make_pair(extremes[i], extremes[i + 1]));
      }
    }
  }
  const std::array<std::uint64_t, 2> descending{1, 0};
  EXPECT_THROW((TypeParam(descending)), std::invalid_argument);
}

TYPED_TEST(MonotoneSampleEncodingTest, EveryLowWidthAndCrossWordOffset) {
  for (unsigned width = 0; width <= 59; ++width) {
    std::array<std::uint64_t, 32> values{};
    const auto maximum = width == 0 ? 31 : (std::uint64_t{32} << (width - 1));
    for (std::size_t i = 0; i < values.size(); ++i) {
      // Division before multiplication keeps full-width synthetic values valid.
      values[i] = (maximum / 31) * i;
    }
    values.back() = maximum;
    const TypeParam samples(values);
    for (std::size_t i = 0; i < values.size(); ++i) {
      EXPECT_EQ(samples[i], values[i]) << "low width=" << width;
      if constexpr (requires { samples.bounds(i); }) {
        if (i + 1 < values.size()) {
          EXPECT_EQ(samples.bounds(i), std::make_pair(values[i], values[i + 1]))
              << "low width=" << width;
        }
      }
    }
  }
}

template <class Support>
void check_full_sampling_endpoint() {
  // A clipped endpoint of 65,536 needs 17 bits even though local search keys
  // often fit in 16. Force cold localization in both select directions.
  constexpr std::size_t kSize = 65536 + 17;
  for (const bool invert : {false, true}) {
    const std::vector<std::uint64_t> words((kSize + 63) / 64,
                                           invert ? 0 : ~std::uint64_t{0});
    const Support support(words, kSize);
    for (const auto rank :
         {std::size_t{1}, std::size_t{16384}, std::size_t{32768},
          std::size_t{65535}, std::size_t{65536}, kSize}) {
      EXPECT_EQ(invert ? support.select0(rank) : support.select(rank),
                rank - 1);
    }
  }
}

TEST(CompressedRankSelectTest, WideEndpointAtFullSamplingInterval) {
  check_full_sampling_endpoint<pixie::RankSelectSupportCompressed<65536, 1>>();
}

TEST(CompressedRankSelectTest, UnsignedPackedEndpointGroups) {
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 65 * kSuper;
  using Support = pixie::RankSelectSupportCompressed<65536, 1>;
  for (const bool invert : {false, true}) {
    std::vector<std::uint64_t> words(kSize / 64, 0);
    for (std::size_t super = 0; super < 65; ++super) {
      std::fill_n(words.begin() + super * kSuper / 64, 16, ~std::uint64_t{0});
    }
    if (invert) {
      for (auto& word : words) {
        word = ~word;
      }
    }
    const Support support(words, kSize);
    // Cross both vector groups, the unsigned sign bit, and the coarse-sample
    // boundary. The last interval has just one occupied superblock.
    for (const std::size_t rank : {1u, 1024u, 16384u, 32767u, 32768u, 32769u,
                                   65535u, 65536u, 65537u, 66560u}) {
      const auto expected = (rank - 1) / 1024 * kSuper + (rank - 1) % 1024;
      EXPECT_EQ(invert ? support.select0(rank) : support.select(rank),
                expected);
    }
  }
}

TEST(CompressedRankSelectTest, HomogeneousSamplesSkipGapsAndHandlePartialTail) {
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 80 * kSuper + 37;
  std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
  for (const auto super : {0u, 40u, 79u}) {
    std::fill_n(words.begin() + super * kSuper / 64, kSuper / 64,
                ~std::uint64_t{0});
  }
  words.back() = ~std::uint64_t{0};
  const pixie::RankSelectSupportCompressed<> support(words, kSize);
  const pixie::RankSelectSupport<> reference(words, kSize);
  std::mt19937_64 rng(42);
  for (unsigned i = 0; i < 10000; ++i) {
    const auto one_rank = 1 + rng() % reference.rank(kSize);
    const auto zero_rank = 1 + rng() % reference.rank0(kSize);
    EXPECT_EQ(support.select(one_rank), reference.select(one_rank));
    EXPECT_EQ(support.select0(zero_rank), reference.select0(zero_rank));
  }
  EXPECT_EQ(support.select(reference.rank(kSize)), kSize - 1);
  EXPECT_EQ(support.select(reference.rank(kSize) + 1), kSize);
}

TYPED_TEST(RankSelectSpecificationTest, WideEndpointGroupsAndPartialTails) {
  constexpr std::size_t kSuper = 16384;
  for (const auto count : {31u, 32u, 33u, 63u, 64u, 65u}) {
    for (const bool gaps : {false, true}) {
      const std::size_t stride = gaps ? 2 : 1;
      const auto size = (count * stride + 5) * kSuper + 13;
      std::vector<std::uint64_t> words((size + 63) / 64, 0);
      std::vector<std::size_t> positions;
      for (std::size_t i = 0; i < count; ++i) {
        const auto position = (3 + i * stride) * kSuper + (i % 32) * 512 + 1;
        words[position / 64] |= std::uint64_t{1} << (position % 64);
        positions.push_back(position);
      }
      for (const bool invert : {false, true}) {
        const TypeParam support(words, size);
        const pixie::RankSelectSupport<> reference(words, size);
        for (std::size_t i = 0; i < count; ++i) {
          EXPECT_EQ(invert ? support.select0(i + 1) : support.select(i + 1),
                    positions[i]);
        }
        std::mt19937_64 rng(42);
        for (unsigned i = 0; i < 256; ++i) {
          const auto rank = 1 + rng() % (size - count);
          EXPECT_EQ(invert ? support.select(rank) : support.select0(rank),
                    invert ? reference.select(rank) : reference.select0(rank));
        }
        for (auto& word : words) {
          word = ~word;
        }
      }
    }
  }
}

TYPED_TEST(RankSelectSpecificationTest, MixedWideDirectoryWithNarrowTail) {
  constexpr std::size_t kSuper = 16384;
  constexpr std::size_t kSize = 80 * kSuper + 37;
  std::vector<std::uint64_t> words((kSize + 63) / 64, 0);
  for (const auto super : {0u, 40u, 79u}) {
    std::fill_n(words.begin() + super * kSuper / 64, kSuper / 64,
                ~std::uint64_t{0});
  }
  // Force mixed populations so construction cannot use the whole-index
  // homogeneous shortcut. Both wide samples cross an empty gap.
  words[3] &= ~(std::uint64_t{1} << 17);
  const TypeParam support(words, kSize);
  const pixie::RankSelectSupport<> reference(words, kSize);
  for (std::size_t rank = 1; rank <= 3 * kSuper; ++rank) {
    if (rank % 4096 <= 1 || rank % 4096 == 4095) {
      EXPECT_EQ(support.select(rank), reference.select(rank));
    }
  }
  std::mt19937_64 rng(42);
  for (unsigned i = 0; i < 1024; ++i) {
    const auto rank = 1 + rng() % (kSize - 3 * kSuper);
    EXPECT_EQ(support.select0(rank), reference.select0(rank));
  }
  EXPECT_EQ(support.select(3 * kSuper + 1), kSize);
}

TEST(CompressedRankSelectTest, LoadedBlockHandlesEmptyAndFullWordLanes) {
  std::array<std::uint64_t, 8> words{};
  constexpr std::array<unsigned, 8> counts{0, 1, 2, 7, 31, 32, 63, 64};
  for (std::size_t i = 0; i < words.size(); ++i) {
    words[i] = ::first_bits_mask(counts[i]);
  }
  for (unsigned round = 0; round < 2; ++round) {
    std::size_t ones = 0;
    std::size_t zeros = 0;
    for (std::size_t position = 0; position < 512; ++position) {
      if ((words[position / 64] >> (position % 64)) & 1) {
        EXPECT_EQ(::pixie::detail::compressed::select_loaded_block<false>(
                      words.data(), ones++),
                  position);
      } else {
        EXPECT_EQ(::pixie::detail::compressed::select_loaded_block<true>(
                      words.data(), zeros++),
                  position);
      }
    }
    for (auto& word : words) {
      word = ~word;
    }
  }
}

TEST(RankSelectSupportTest, AcceptsStorageSourceWithoutCopying) {
  pixie::AlignedStorage source(64);
  source.writable_words64()[0] = 0b101101;
  const pixie::RankSelectSupport<> support(source, 6);

  EXPECT_EQ(support.rank(6), 4u);
  EXPECT_EQ(support.select(4), 5u);
  source.writable_words64()[0] = 0b000001;
  EXPECT_EQ(support[2], 0);
}

TEST(RankSelectSupportTest, AcceptsPartialWordAlignedStorageSource) {
  pixie::AlignedStorage storage(1);
  storage.writable_bytes()[0] = std::byte{1};

  const pixie::RankSelectSupport support(storage, 1);

  EXPECT_EQ(support.size(), 1u);
  EXPECT_EQ(support.rank(1), 1u);
  EXPECT_EQ(support.select(1), 0u);

  const pixie::RankSelectSupport clamped(storage, 65);
  EXPECT_EQ(clamped.size(), 8u);
}

TEST(RankSelectSupportTest, OwningMetadataDeserializationRoundTrips) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260718);
  for (std::uint64_t& word : words) {
    word = rng();
  }

  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  std::vector<std::byte> artifact = output.take();
  pixie::BinaryReader reader(artifact);
  const pixie::RankSelectSupport<> restored =
      pixie::RankSelectSupport<>::deserialize(reader, words);
  EXPECT_TRUE(reader.empty());

  artifact.clear();
  artifact.shrink_to_fit();
  for (std::size_t position = 0; position <= kBitCount; ++position) {
    EXPECT_EQ(restored.rank(position), original.rank(position));
    EXPECT_EQ(restored.rank0(position), original.rank0(position));
  }
  const std::size_t ones = original.rank(kBitCount);
  const std::size_t zeros = original.rank0(kBitCount);
  for (std::size_t rank = 1; rank <= ones + 1; ++rank) {
    EXPECT_EQ(restored.select(rank), original.select(rank));
  }
  for (std::size_t rank = 1; rank <= zeros + 1; ++rank) {
    EXPECT_EQ(restored.select0(rank), original.select0(rank));
  }
}

TEST(RankSelectSupportTest, SerializesOnlySimdRoundedBasicBlockMetadata) {
  constexpr std::array<std::size_t, 6> bit_counts = {
      0, 1, 32 * 512, 32 * 512 + 1, 65536, 65537};
  for (const std::size_t bit_count : bit_counts) {
    SCOPED_TRACE(bit_count);
    const std::vector<std::uint64_t> words((bit_count + 63) / 64);
    const pixie::RankSelectSupport<> support(
        words, bit_count, pixie::RankSelectSupport<>::SelectSupport::kNone);
    pixie::VectorOutputSink output;
    pixie::BinaryWriter writer(output);
    support.serialize(writer);
    writer.finish();
    const std::vector<std::byte> artifact = output.take();
    const RankSelectSampleLayout layout = locate_rank_select_samples(artifact);
    const std::size_t data_superblocks =
        bit_count == 0 ? 0 : 1 + (bit_count - 1) / 65536;
    const std::size_t data_basicblocks =
        bit_count == 0 ? 0 : 1 + (bit_count - 1) / 512;
    const std::size_t stored_basicblocks = (data_basicblocks + 31) / 32 * 32;
    EXPECT_EQ(layout.super_bytes,
              (data_superblocks + 1) * sizeof(std::uint64_t));
    EXPECT_EQ(layout.basic_bytes, stored_basicblocks * sizeof(std::uint16_t));
    EXPECT_EQ(layout.sample_bytes, 0u);
  }
}

TEST(RankSelectSupportTest, QuickAndFullValidationAcceptValidMetadata) {
  constexpr std::size_t kBitCount = 65537;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260816);
  for (std::uint64_t& word : words) {
    word = rng();
  }
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> artifact = output.take();

  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader reader(artifact);
    const auto restored =
        pixie::RankSelectSupport<>::deserialize(reader, words, validation);
    EXPECT_TRUE(reader.empty());
    EXPECT_EQ(restored.rank(kBitCount), original.rank(kBitCount));
  }
}

TEST(RankSelectSupportTest,
     FullValidationRejectsSourceAndMetadataMismatchesTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();

  std::vector<std::uint64_t> different_words = words;
  different_words[0] ^= 1;
  pixie::BinaryReader quick_reader(valid);
  EXPECT_NO_THROW((void)pixie::RankSelectSupport<>::deserialize(
      quick_reader, different_words, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(quick_reader.empty());

  pixie::BinaryReader full_reader(valid);
  EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(
                   full_reader, different_words,
                   pixie::DeserializationValidation::kFull),
               std::invalid_argument);
  EXPECT_EQ(full_reader.position(), 0u);

  std::vector<std::byte> bad_rank = valid;
  overwrite_u64(bad_rank, 2 * sizeof(std::uint64_t),
                original.rank(kBitCount) + 1);
  pixie::BinaryReader metadata_quick_reader(bad_rank);
  EXPECT_NO_THROW((void)pixie::RankSelectSupport<>::deserialize(
      metadata_quick_reader, words, pixie::DeserializationValidation::kQuick));
  EXPECT_TRUE(metadata_quick_reader.empty());

  pixie::BinaryReader metadata_full_reader(bad_rank);
  EXPECT_THROW(
      (void)pixie::RankSelectSupport<>::deserialize(
          metadata_full_reader, words, pixie::DeserializationValidation::kFull),
      std::invalid_argument);
  EXPECT_EQ(metadata_full_reader.position(), 0u);
}

TEST(RankSelectSupportTest, SerializesDirectlyToMappedFile) {
  constexpr std::size_t kBitCount = 1025;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64);
  std::mt19937_64 rng(20260721);
  for (std::uint64_t& word : words) {
    word = rng();
  }
  const pixie::RankSelectSupport<> original(words, kBitCount);
  const auto path = std::filesystem::temp_directory_path() /
                    "pixie_rank_select_serialization_test.bin";
  std::filesystem::remove(path);
  {
    pixie::io::FileOutputSink output(path);
    std::array<std::byte, 13> staging{};
    pixie::BinaryWriter writer(output, staging);
    original.serialize(writer);
    writer.finish();
  }

  pixie::io::MappedFile file(path);
  for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                pixie::DeserializationValidation::kFull}) {
    pixie::BinaryReader reader(file.as_bytes());
    const auto restored =
        pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(
            reader, words, validation);
    EXPECT_TRUE(reader.empty());
    for (std::size_t position = 0; position <= kBitCount; position += 17) {
      EXPECT_EQ(restored.rank(position), original.rank(position));
      EXPECT_EQ(restored.rank0(position), original.rank0(position));
    }
  }
  std::filesystem::remove(path);
}

TEST(RankSelectSupportTest,
     DeserializationRejectsOutOfRangeSelectSamplesTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();
  const RankSelectSampleLayout layout = locate_rank_select_samples(valid);
  ASSERT_NE(layout.select1_count, 0u);
  ASSERT_NE(layout.select0_count, 0u);

  const auto expect_rejected = [&](std::size_t sample) {
    std::vector<std::byte> artifact = valid;
    overwrite_u64(artifact,
                  layout.storage_offset + sample * sizeof(std::uint64_t),
                  std::numeric_limits<std::uint64_t>::max());

    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader owning_reader(artifact);
      EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(
                       owning_reader, words, validation),
                   std::invalid_argument);
      EXPECT_EQ(owning_reader.position(), 0u);

      pixie::BinaryReader view_reader(artifact);
      EXPECT_THROW(
          (void)
              pixie::RankSelectSupport<pixie::ReadOnlyStorageView>::deserialize(
                  view_reader, words, validation),
          std::invalid_argument);
      EXPECT_EQ(view_reader.position(), 0u);
    }
  };

  expect_rejected(layout.select1_begin);
  expect_rejected(layout.select0_begin);
}

TEST(RankSelectSupportTest,
     DeserializationRejectsMalformedFixedMetadataTransactionally) {
  constexpr std::size_t kBitCount = 4097;
  std::vector<std::uint64_t> words((kBitCount + 63) / 64,
                                   0xaaaaaaaaaaaaaaaaULL);
  const pixie::RankSelectSupport<> original(words, kBitCount);
  pixie::VectorOutputSink output;
  pixie::BinaryWriter writer(output);
  original.serialize(writer);
  writer.finish();
  const std::vector<std::byte> valid = output.take();

  const auto expect_rejected = [&](std::vector<std::byte> artifact) {
    for (const auto validation : {pixie::DeserializationValidation::kQuick,
                                  pixie::DeserializationValidation::kFull}) {
      pixie::BinaryReader reader(artifact);
      EXPECT_THROW((void)pixie::RankSelectSupport<>::deserialize(reader, words,
                                                                 validation),
                   std::invalid_argument);
      EXPECT_EQ(reader.position(), 0u);
    }
  };

  auto bad_padded_size = valid;
  overwrite_u64(bad_padded_size, sizeof(std::uint64_t), 0);
  expect_rejected(std::move(bad_padded_size));

  auto excessive_rank = valid;
  overwrite_u64(excessive_rank, 2 * sizeof(std::uint64_t), kBitCount + 1);
  expect_rejected(std::move(excessive_rank));

  auto invalid_select1_begin = valid;
  overwrite_u64(invalid_select1_begin, 3 * sizeof(std::uint64_t),
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(invalid_select1_begin));

  auto missing_select1_samples = valid;
  overwrite_u64(missing_select1_samples, 4 * sizeof(std::uint64_t), 0);
  expect_rejected(std::move(missing_select1_samples));

  auto invalid_select0_begin = valid;
  overwrite_u64(invalid_select0_begin, 5 * sizeof(std::uint64_t),
                std::numeric_limits<std::uint64_t>::max());
  expect_rejected(std::move(invalid_select0_begin));

  auto missing_select0_samples = valid;
  overwrite_u64(missing_select0_samples, 6 * sizeof(std::uint64_t), 0);
  expect_rejected(std::move(missing_select0_samples));

  auto invalid_support = valid;
  overwrite_u32(invalid_support, 7 * sizeof(std::uint64_t), 4);
  expect_rejected(std::move(invalid_support));

  auto invalid_boolean = valid;
  overwrite_u32(invalid_boolean,
                7 * sizeof(std::uint64_t) + sizeof(std::uint32_t), 2);
  expect_rejected(std::move(invalid_boolean));
}

TEST(CompressedPositionSampleTest, GroupEdgesLowWordCrossingsAndExtremes) {
  for (unsigned width : {0u, 1u, 7u, 13u, 31u, 58u}) {
    std::vector<std::uint64_t> values;
    for (std::size_t i = 0; i < 67; ++i) {
      values.push_back((i / 2) << width);
    }
    values.back() = std::numeric_limits<std::uint64_t>::max();
    const pixie::detail::compressed::PositionSamples samples(values);
    for (std::size_t i = 0; i + 1 < values.size(); ++i) {
      EXPECT_EQ(samples.bounds(i), std::make_pair(values[i], values[i + 1]));
    }
  }
}

TEST(CompressedPositionSampleTest, InterleavedBitmapAcrossGroupBoundaries) {
  std::vector<std::uint64_t> values(100);
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = i * 17003;
  }
  const std::array<std::uint64_t, 7> wide{0, 17, 31, 32, 63, 64, 98};
  const pixie::detail::compressed::PositionSamples samples(values, wide);
  const pixie::detail::compressed::PositionSamples without_bitmap(values);
  std::size_t next = 0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    EXPECT_EQ(samples[i], values[i]);
    EXPECT_FALSE(without_bitmap.wide_entry(i));
    if (next < wide.size() && i == wide[next]) {
      EXPECT_EQ(samples.wide_entry(i), next++);
    } else {
      EXPECT_FALSE(samples.wide_entry(i));
    }
    if (i + 1 < values.size()) {
      EXPECT_EQ(samples.bounds(i), std::make_pair(values[i], values[i + 1]));
    }
  }
}

TEST(CompressedRankSelectTest, SkewedBlocksGapsAndTerminalSamples) {
  constexpr std::size_t kSuperBits = 16384;
  constexpr std::size_t kBits = 80 * kSuperBits + 137;
  constexpr std::array<std::size_t, 10> kPopulations{0,     1,  4095, 0, 0,
                                                     16384, 17, 1024, 0, 8192};
  std::vector<std::uint64_t> words((kBits + 63) / 64);
  for (std::size_t super = 0; super < 81; ++super) {
    const auto count = std::min(kPopulations[super % kPopulations.size()],
                                kBits - super * kSuperBits);
    for (std::size_t bit = 0; bit < count; ++bit) {
      const auto position = super * kSuperBits + bit;
      words[position / 64] |= std::uint64_t{1} << (position % 64);
    }
  }
  // Both trailing empty blocks and a partially occupied last block exercise
  // the terminal sample. Complementing also swaps the skew between directions.
  for (unsigned pattern = 0; pattern < 3; ++pattern) {
    const pixie::RankSelectSupport<> reference(words, kBits);
    const auto check = [&]<class Support>() {
      const Support support(words, kBits);
      for (std::size_t k = 1; k <= reference.rank(kBits); ++k) {
        ASSERT_EQ(support.select(k), reference.select(k)) << k;
      }
      for (std::size_t k = 1; k <= reference.rank0(kBits); ++k) {
        ASSERT_EQ(support.select0(k), reference.select0(k)) << k;
      }
    };
    check.operator()<pixie::RankSelectSupportCompressed<>>();
    check.operator()<pixie::RankSelectSupportCompressed<64>>();
    if (pattern == 0) {
      words.back() |= 1;
    } else {
      for (auto& word : words) {
        word = ~word;
      }
    }
  }
}

}  // namespace
