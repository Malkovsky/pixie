#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/experimental/select.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
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

namespace {

using Block = std::array<uint64_t, 8>;

struct SelectCase {
  Block block;
  uint64_t rank1;
  uint64_t rank0;
};

uint64_t make_word(std::mt19937_64& rng, int fill_permille) {
  switch (fill_permille) {
    case 125:
      return rng() & rng() & rng();
    case 875:
      return rng() | rng() | rng();
    case 500:
    default:
      return rng();
  }
}

uint64_t count_ones(const Block& block) {
  uint64_t count = 0;
  for (uint64_t word : block) {
    count += std::popcount(word);
  }
  return count;
}

void ensure_has_ones_and_zeros(Block& block) {
  const uint64_t ones = count_ones(block);
  if (ones == 0) {
    block[0] = 1;
  } else if (ones == 512) {
    block[0] = ~uint64_t{1};
  }
}

uint64_t rank_for_count(uint64_t count, int rank_mode, std::mt19937_64& rng) {
  switch (rank_mode) {
    case 0:
      return 0;
    case 1:
      return count / 2;
    case 2:
      return count - 1;
    case 3:
    default:
      return rng() % count;
  }
}

std::vector<SelectCase> make_cases(int fill_permille,
                                   int rank_mode,
                                   size_t num_blocks = 4096) {
  std::mt19937_64 rng(42 + static_cast<uint64_t>(fill_permille) * 17 +
                      static_cast<uint64_t>(rank_mode) * 131);
  std::vector<SelectCase> cases(num_blocks);
  for (SelectCase& select_case : cases) {
    for (uint64_t& word : select_case.block) {
      word = make_word(rng, fill_permille);
    }
    ensure_has_ones_and_zeros(select_case.block);

    const uint64_t ones = count_ones(select_case.block);
    const uint64_t zeros = 512 - ones;
    select_case.rank1 = rank_for_count(ones, rank_mode, rng);
    select_case.rank0 = rank_for_count(zeros, rank_mode, rng);
  }
  return cases;
}

}  // namespace

#define PIXIE_SELECT512_ARGS               \
  ->ArgNames({"FillPermille", "RankMode"}) \
      ->ArgsProduct({{125, 500, 875}, {0, 1, 2, 3}})

#define PIXIE_DEFINE_SELECT512_BENCHMARK(Name, Function, RankMember)        \
  static void BM_##Name(benchmark::State& state) {                          \
    const auto cases = make_cases(static_cast<int>(state.range(0)),         \
                                  static_cast<int>(state.range(1)));        \
    const size_t num_cases = cases.size();                                  \
    size_t idx = 0;                                                         \
    uint64_t result = 0;                                                    \
    for (auto _ : state) {                                                  \
      const SelectCase& select_case = cases[idx % num_cases];               \
      result ^= Function(select_case.block.data(), select_case.RankMember); \
      benchmark::DoNotOptimize(result);                                     \
      ++idx;                                                                \
    }                                                                       \
    state.SetItemsProcessed(state.iterations());                            \
  }                                                                         \
  BENCHMARK(BM_##Name) PIXIE_SELECT512_ARGS

PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_Current, ::select_512, rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarPDEP,
                                 select_512_scalar_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarByteLUT,
                                 select_512_scalar_byte_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarBinaryLUT,
                                 select_512_scalar_binary_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarBroadwordLUT,
                                 select_512_scalar_broadword_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarVignaBroadword,
                                 select_512_scalar_vigna_broadword,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarNoPDEP,
                                 select_512_scalar_no_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ScalarExperimentalDefault,
                                 select_512_scalar_experimental_default,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2PDEP,
                                 select_512_avx2_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2ByteLUT,
                                 select_512_avx2_byte_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2BinaryLUT,
                                 select_512_avx2_binary_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2BroadwordLUT,
                                 select_512_avx2_broadword_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2VignaBroadword,
                                 select_512_avx2_vigna_broadword,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX2NoPDEP,
                                 select_512_avx2_no_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512PDEP,
                                 select_512_avx512_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512ByteLUT,
                                 select_512_avx512_byte_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512BinaryLUT,
                                 select_512_avx512_binary_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512BroadwordLUT,
                                 select_512_avx512_broadword_lut,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512VignaBroadword,
                                 select_512_avx512_vigna_broadword,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512NoPDEP,
                                 select_512_avx512_no_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512Hybrid1PDEP,
                                 select_512_avx512_hybrid1_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512Hybrid2PDEP,
                                 select_512_avx512_hybrid2_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512Hybrid4PDEP,
                                 select_512_avx512_hybrid4_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512Hybrid2NoPDEP,
                                 select_512_avx512_hybrid2_no_pdep,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_AVX512ExperimentalDefault,
                                 select_512_avx512_experimental_default,
                                 rank1);
PIXIE_DEFINE_SELECT512_BENCHMARK(Select512_ExperimentalDefault,
                                 select_512_experimental_default,
                                 rank1);

PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_Current, ::select0_512, rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarPDEP,
                                 select0_512_scalar_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarByteLUT,
                                 select0_512_scalar_byte_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarBinaryLUT,
                                 select0_512_scalar_binary_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarBroadwordLUT,
                                 select0_512_scalar_broadword_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarVignaBroadword,
                                 select0_512_scalar_vigna_broadword,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarNoPDEP,
                                 select0_512_scalar_no_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ScalarExperimentalDefault,
                                 select0_512_scalar_experimental_default,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2PDEP,
                                 select0_512_avx2_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2ByteLUT,
                                 select0_512_avx2_byte_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2BinaryLUT,
                                 select0_512_avx2_binary_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2BroadwordLUT,
                                 select0_512_avx2_broadword_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2VignaBroadword,
                                 select0_512_avx2_vigna_broadword,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX2NoPDEP,
                                 select0_512_avx2_no_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512PDEP,
                                 select0_512_avx512_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512ByteLUT,
                                 select0_512_avx512_byte_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512BinaryLUT,
                                 select0_512_avx512_binary_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512BroadwordLUT,
                                 select0_512_avx512_broadword_lut,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512VignaBroadword,
                                 select0_512_avx512_vigna_broadword,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512NoPDEP,
                                 select0_512_avx512_no_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512Hybrid1PDEP,
                                 select0_512_avx512_hybrid1_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512Hybrid2PDEP,
                                 select0_512_avx512_hybrid2_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512Hybrid4PDEP,
                                 select0_512_avx512_hybrid4_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512Hybrid2NoPDEP,
                                 select0_512_avx512_hybrid2_no_pdep,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_AVX512ExperimentalDefault,
                                 select0_512_avx512_experimental_default,
                                 rank0);
PIXIE_DEFINE_SELECT512_BENCHMARK(SelectZero512_ExperimentalDefault,
                                 select0_512_experimental_default,
                                 rank0);

#undef PIXIE_DEFINE_SELECT512_BENCHMARK
#undef PIXIE_SELECT512_ARGS
