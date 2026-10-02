#include "three_star_adapter.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kWordBits = 64;
constexpr std::size_t kWordsPerSuperblock = 1024;
constexpr std::size_t kSentinelSuperblocks = 5;

std::uint64_t logical_word(std::span<const std::uint64_t> source,
                           std::size_t word_index,
                           std::size_t num_bits) {
  std::uint64_t word = source[word_index];
  const std::size_t remaining = num_bits - word_index * kWordBits;
  if (remaining < kWordBits) {
    word &= (std::uint64_t{1} << remaining) - 1;
  }
  return word;
}

}  // namespace

uint64_t SAMPLEDIST_EXP = 0;
uint64_t SAMPLEDIST = 0;

bool show_overhead = false;
std::string BITGEN;
std::string BVNAME;
bool doBench = false;
bool doValid = false;
bool doValidw = false;
bool doMultibench = false;
uint64_t synth_bits = 0;
bool useSynth = false;
bool do_manual_bits = false;
bool do_fractals = false;
bool do_chunky = false;
bool do_smooth = false;
bool do_quicksmooth = false;
bool do_gidneysmooth = false;
bool do_alternate = false;
bool do_every_kth = false;
bool do_bimodal = false;
bool do_patches = false;
bool do_sinus = false;
int bimodal_factor = 0;
int every_kth_spacing = 0;
int synth_seqsize = 0;
int synth_01ratio = 0;
bool invert_all_bits = false;
bool do_random_queries = false;
uint64_t synth_accesses = 0;
uint64_t synth_ranks = 0;
uint64_t synth_selects = 0;
uint64_t synth_select1s = 0;
uint64_t synth_select0s = 0;
std::vector<uint64_t> query_type_counter;
uint64_t N_queries = 0;
uint64_t N_bits = 0;
uint64_t N_ones = 0;
uint64_t N_zeros = 0;
uint64_t N_words = 0;
uint64_t seed = 0;
double eff_01ratio = 0.0;
int parameter_cycles = 1;
int curr_parameter_cycle = 0;
int curr_instance = 0;
double sin_threshhold = 0.0;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>> bits;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>> queries;
std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGNMENT>>
    bits_64bit_reversed;
uint64_t total_used_space_in_bits = 0;
tp tCopyFinished;
double dRead = 0.0;
double dCopy = 0.0;
double dBuild = 0.0;
double dQueries = 0.0;
double dPerQuery = 0.0;

std::string get_unixtimestamp() {
  return {};
}

namespace pixie::benchmarks {

ThreeStarRankSelectSupport::ThreeStarRankSelectSupport(
    std::span<const std::uint64_t> source_words,
    std::size_t num_bits)
    : num_bits_(std::min(num_bits, source_words.size() * kWordBits)) {
  const std::size_t source_word_count = (num_bits_ + kWordBits - 1) / kWordBits;
  const std::size_t padded_word_count =
      ((source_word_count + kWordsPerSuperblock - 1) / kWordsPerSuperblock) *
      kWordsPerSuperblock;
  const std::size_t stored_words =
      padded_word_count + kSentinelSuperblocks * kWordsPerSuperblock;
  N_bits = num_bits_;
  // Upstream builds complete superblocks and copies sentinel blocks without
  // checking the input vector's end. Supply initialized padding for both.
  N_words = padded_word_count;
  bits.assign(stored_words, 0);
  one_count_ = 0;
  for (std::size_t word_index = 0; word_index < source_word_count;
       ++word_index) {
    const std::uint64_t word =
        logical_word(source_words, word_index, num_bits_);
    bits[word_index] = word;
    one_count_ += std::popcount(word);
  }
  N_ones = one_count_;
  N_zeros = num_bits_ - one_count_;

  // Match the paper's uncompressed robust configuration: L0=2048, a*,
  // alpha=16, and its two-level L0/L1 summary tree.
  ALPHA = 16;
  SUMMARY_LEVELS = 1;
  TREE3_STRATEGY = GET_LIFTED_CUBIC_THEORY_PARAMS;
  BV_COMPRESSION = 0;

  support_.build_auxiliaries();

  source_copy_bytes_ = stored_words * sizeof(std::uint64_t);
}

}  // namespace pixie::benchmarks
