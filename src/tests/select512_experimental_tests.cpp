#include <gtest/gtest.h>
#include <pixie/experimental/select.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

using pixie::experimental::select0_512_avx2_binary_lut;
using pixie::experimental::select0_512_avx2_broadword_lut;
using pixie::experimental::select0_512_avx2_byte_lut;
using pixie::experimental::select0_512_avx2_no_pdep;
using pixie::experimental::select0_512_avx2_pdep;
using pixie::experimental::select0_512_avx2_vigna_broadword;
using pixie::experimental::select0_512_avx512_binary_lut;
using pixie::experimental::select0_512_avx512_broadword_lut;
using pixie::experimental::select0_512_avx512_byte_lut;
using pixie::experimental::select0_512_avx512_experimental_default;
using pixie::experimental::select0_512_avx512_hybrid1_pdep;
using pixie::experimental::select0_512_avx512_hybrid2_no_pdep;
using pixie::experimental::select0_512_avx512_hybrid2_pdep;
using pixie::experimental::select0_512_avx512_hybrid4_pdep;
using pixie::experimental::select0_512_avx512_no_pdep;
using pixie::experimental::select0_512_avx512_pdep;
using pixie::experimental::select0_512_avx512_vigna_broadword;
using pixie::experimental::select0_512_experimental_default;
using pixie::experimental::select0_512_scalar_binary_lut;
using pixie::experimental::select0_512_scalar_broadword_lut;
using pixie::experimental::select0_512_scalar_byte_lut;
using pixie::experimental::select0_512_scalar_experimental_default;
using pixie::experimental::select0_512_scalar_no_pdep;
using pixie::experimental::select0_512_scalar_pdep;
using pixie::experimental::select0_512_scalar_vigna_broadword;
using pixie::experimental::select_512_avx2_binary_lut;
using pixie::experimental::select_512_avx2_broadword_lut;
using pixie::experimental::select_512_avx2_byte_lut;
using pixie::experimental::select_512_avx2_no_pdep;
using pixie::experimental::select_512_avx2_pdep;
using pixie::experimental::select_512_avx2_vigna_broadword;
using pixie::experimental::select_512_avx512_binary_lut;
using pixie::experimental::select_512_avx512_broadword_lut;
using pixie::experimental::select_512_avx512_byte_lut;
using pixie::experimental::select_512_avx512_experimental_default;
using pixie::experimental::select_512_avx512_hybrid1_pdep;
using pixie::experimental::select_512_avx512_hybrid2_no_pdep;
using pixie::experimental::select_512_avx512_hybrid2_pdep;
using pixie::experimental::select_512_avx512_hybrid4_pdep;
using pixie::experimental::select_512_avx512_no_pdep;
using pixie::experimental::select_512_avx512_pdep;
using pixie::experimental::select_512_avx512_vigna_broadword;
using pixie::experimental::select_512_experimental_default;
using pixie::experimental::select_512_scalar_binary_lut;
using pixie::experimental::select_512_scalar_broadword_lut;
using pixie::experimental::select_512_scalar_byte_lut;
using pixie::experimental::select_512_scalar_experimental_default;
using pixie::experimental::select_512_scalar_no_pdep;
using pixie::experimental::select_512_scalar_pdep;
using pixie::experimental::select_512_scalar_vigna_broadword;
using pixie::experimental::select_64_binary_lut;
using pixie::experimental::select_64_broadword_lut;
using pixie::experimental::select_64_byte_lut;
using pixie::experimental::select_64_experimental_default;
using pixie::experimental::select_64_no_pdep;
using pixie::experimental::select_64_pdep;
using pixie::experimental::select_64_vigna_broadword;

namespace {

using Block = std::array<uint64_t, 8>;
using SelectFn = uint64_t (*)(const uint64_t*, uint64_t);

struct SelectVariant {
  const char* name;
  SelectFn select1;
  SelectFn select0;
};

static const SelectVariant kVariants[] = {
    {"scalar_pdep", select_512_scalar_pdep, select0_512_scalar_pdep},
    {"scalar_byte_lut", select_512_scalar_byte_lut,
     select0_512_scalar_byte_lut},
    {"scalar_binary_lut", select_512_scalar_binary_lut,
     select0_512_scalar_binary_lut},
    {"scalar_broadword_lut", select_512_scalar_broadword_lut,
     select0_512_scalar_broadword_lut},
    {"scalar_vigna_broadword", select_512_scalar_vigna_broadword,
     select0_512_scalar_vigna_broadword},
    {"scalar_no_pdep", select_512_scalar_no_pdep, select0_512_scalar_no_pdep},
    {"scalar_experimental_default", select_512_scalar_experimental_default,
     select0_512_scalar_experimental_default},
    {"avx2_pdep", select_512_avx2_pdep, select0_512_avx2_pdep},
    {"avx2_byte_lut", select_512_avx2_byte_lut, select0_512_avx2_byte_lut},
    {"avx2_binary_lut", select_512_avx2_binary_lut,
     select0_512_avx2_binary_lut},
    {"avx2_broadword_lut", select_512_avx2_broadword_lut,
     select0_512_avx2_broadword_lut},
    {"avx2_vigna_broadword", select_512_avx2_vigna_broadword,
     select0_512_avx2_vigna_broadword},
    {"avx2_no_pdep", select_512_avx2_no_pdep, select0_512_avx2_no_pdep},
    {"avx512_pdep", select_512_avx512_pdep, select0_512_avx512_pdep},
    {"avx512_byte_lut", select_512_avx512_byte_lut,
     select0_512_avx512_byte_lut},
    {"avx512_binary_lut", select_512_avx512_binary_lut,
     select0_512_avx512_binary_lut},
    {"avx512_broadword_lut", select_512_avx512_broadword_lut,
     select0_512_avx512_broadword_lut},
    {"avx512_vigna_broadword", select_512_avx512_vigna_broadword,
     select0_512_avx512_vigna_broadword},
    {"avx512_no_pdep", select_512_avx512_no_pdep, select0_512_avx512_no_pdep},
    {"avx512_hybrid1_pdep", select_512_avx512_hybrid1_pdep,
     select0_512_avx512_hybrid1_pdep},
    {"avx512_hybrid2_pdep", select_512_avx512_hybrid2_pdep,
     select0_512_avx512_hybrid2_pdep},
    {"avx512_hybrid4_pdep", select_512_avx512_hybrid4_pdep,
     select0_512_avx512_hybrid4_pdep},
    {"avx512_hybrid2_no_pdep", select_512_avx512_hybrid2_no_pdep,
     select0_512_avx512_hybrid2_no_pdep},
    {"avx512_experimental_default", select_512_avx512_experimental_default,
     select0_512_avx512_experimental_default},
    {"experimental_default", select_512_experimental_default,
     select0_512_experimental_default},
};

uint64_t naive_select_512(const uint64_t* s, uint64_t rank, bool value) {
  uint64_t seen = 0;
  for (uint64_t i = 0; i < 512; ++i) {
    const bool bit = ((s[i >> 6] >> (i & 63)) & 1ull) != 0;
    if (bit == value) {
      if (seen == rank) {
        return i;
      }
      ++seen;
    }
  }
  return 512;
}

uint64_t count_ones(const Block& block) {
  uint64_t count = 0;
  for (uint64_t word : block) {
    count += std::popcount(word);
  }
  return count;
}

std::vector<uint64_t> sample_ranks(uint64_t count) {
  std::vector<uint64_t> ranks;
  if (count == 0) {
    return ranks;
  }
  ranks.push_back(0);
  ranks.push_back(count / 2);
  ranks.push_back(count - 1);
  std::sort(ranks.begin(), ranks.end());
  ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
  return ranks;
}

void expect_select1(const Block& block,
                    uint64_t rank,
                    const char* label,
                    int case_id) {
  const uint64_t expected = naive_select_512(block.data(), rank, true);
  for (const SelectVariant& variant : kVariants) {
    ASSERT_EQ(variant.select1(block.data(), rank), expected)
        << variant.name << " select1 rank=" << rank << " case=" << case_id
        << " label=" << label;
  }
}

void expect_select0(const Block& block,
                    uint64_t rank,
                    const char* label,
                    int case_id) {
  const uint64_t expected = naive_select_512(block.data(), rank, false);
  for (const SelectVariant& variant : kVariants) {
    ASSERT_EQ(variant.select0(block.data(), rank), expected)
        << variant.name << " select0 rank=" << rank << " case=" << case_id
        << " label=" << label;
  }
}

void check_sampled_ranks(const Block& block, const char* label, int case_id) {
  const uint64_t ones = count_ones(block);
  const uint64_t zeros = 512 - ones;
  for (uint64_t rank : sample_ranks(ones)) {
    expect_select1(block, rank, label, case_id);
  }
  for (uint64_t rank : sample_ranks(zeros)) {
    expect_select0(block, rank, label, case_id);
  }
}

}  // namespace

TEST(Select64Experimental, ByteLutMatchesPdep) {
  std::mt19937_64 rng(42);
  for (int t = 0; t < 1000; ++t) {
    const uint64_t word = rng();
    const uint64_t ones = std::popcount(word);
    for (uint64_t rank = 0; rank < ones; ++rank) {
      ASSERT_EQ(select_64_byte_lut(word, rank), select_64_pdep(word, rank))
          << "byte_lut case=" << t << " rank=" << rank;
      ASSERT_EQ(select_64_binary_lut(word, rank), select_64_pdep(word, rank))
          << "binary_lut case=" << t << " rank=" << rank;
      ASSERT_EQ(select_64_broadword_lut(word, rank), select_64_pdep(word, rank))
          << "broadword_lut case=" << t << " rank=" << rank;
      ASSERT_EQ(select_64_vigna_broadword(word, rank),
                select_64_pdep(word, rank))
          << "vigna_broadword case=" << t << " rank=" << rank;
      ASSERT_EQ(select_64_no_pdep(word, rank), select_64_pdep(word, rank))
          << "no_pdep case=" << t << " rank=" << rank;
      ASSERT_EQ(select_64_experimental_default(word, rank),
                select_64_pdep(word, rank))
          << "experimental_default case=" << t << " rank=" << rank;
    }
  }
}

TEST(Select512Experimental, AllOnes) {
  const Block block = {std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max(),
                       std::numeric_limits<uint64_t>::max()};
  for (uint64_t rank = 0; rank < 512; ++rank) {
    expect_select1(block, rank, "all_ones", 0);
  }
}

TEST(Select512Experimental, AllZeros) {
  const Block block = {0, 0, 0, 0, 0, 0, 0, 0};
  for (uint64_t rank = 0; rank < 512; ++rank) {
    expect_select0(block, rank, "all_zeros", 0);
  }
}

TEST(Select512Experimental, Alternating) {
  const Block block = {0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                       0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                       0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull,
                       0xAAAAAAAAAAAAAAAAull, 0x5555555555555555ull};
  for (uint64_t rank = 0; rank < 256; ++rank) {
    expect_select1(block, rank, "alternating", 0);
    expect_select0(block, rank, "alternating", 0);
  }
}

TEST(Select512Experimental, SingleBitPositions) {
  for (int position = 0; position < 512; ++position) {
    Block block = {0, 0, 0, 0, 0, 0, 0, 0};
    block[position >> 6] = uint64_t{1} << (position & 63);
    expect_select1(block, 0, "single_bit", position);
    check_sampled_ranks(block, "single_bit", position);
  }
}

TEST(Select512Experimental, EmptyWordGaps) {
  const Block block = {0, 0x8000000000000000ull,
                       0, 0x0123456789ABCDEFull,
                       0, std::numeric_limits<uint64_t>::max(),
                       0, 0x0000000000000001ull};
  const uint64_t ones = count_ones(block);
  const uint64_t zeros = 512 - ones;
  for (uint64_t rank = 0; rank < ones; ++rank) {
    expect_select1(block, rank, "empty_word_gaps", 0);
  }
  for (uint64_t rank : sample_ranks(zeros)) {
    expect_select0(block, rank, "empty_word_gaps", 0);
  }
}

TEST(Select512Experimental, ExhaustiveLow16) {
  for (uint64_t pattern = 0; pattern < (1ull << 16); ++pattern) {
    Block block = {pattern, 0, 0, 0, 0, 0, 0, 0};
    const uint64_t ones = std::popcount(pattern);
    const uint64_t low_zeros = 16 - ones;
    const uint64_t zeros = 512 - ones;

    for (uint64_t rank = 0; rank < ones; ++rank) {
      expect_select1(block, rank, "exhaustive_low16",
                     static_cast<int>(pattern));
    }
    for (uint64_t rank = 0; rank < low_zeros; ++rank) {
      expect_select0(block, rank, "exhaustive_low16",
                     static_cast<int>(pattern));
    }
    expect_select0(block, zeros - 1, "exhaustive_low16",
                   static_cast<int>(pattern));
  }
}

TEST(Select512Experimental, Random) {
  const int cases = [] {
    const char* env = std::getenv("SELECT512_CASES");
    return env ? std::atoi(env) : 1000;
  }();
  const uint64_t seed = [] {
    const char* env = std::getenv("SELECT512_SEED");
    return env ? std::strtoull(env, nullptr, 10) : 42ull;
  }();

  std::mt19937_64 rng(seed);
  for (int t = 0; t < cases; ++t) {
    Block block;
    for (uint64_t& word : block) {
      word = rng();
    }

    check_sampled_ranks(block, "random", t);

    const uint64_t ones = count_ones(block);
    const uint64_t zeros = 512 - ones;
    for (int i = 0; i < 8; ++i) {
      if (ones != 0) {
        expect_select1(block, rng() % ones, "random", t);
      }
      if (zeros != 0) {
        expect_select0(block, rng() % zeros, "random", t);
      }
    }
  }
}
