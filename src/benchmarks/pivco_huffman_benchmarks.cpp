#include <benchmark/benchmark.h>
#include <pixie/huffman/implementations.h>

#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

using pixie::PivCoHuffman;

// ---------------------------------------------------------------------------
// Configuration:
//   PIVCO_BENCH_SEED - random seed (default 42)
//
// Input sizes are swept by Google Benchmark's range mechanism, not an env
// var: 1 KiB -> 1 MiB -> 1 GiB (multiplier 1024). This exercises the codec
// from cache-resident to gigabyte scale in a single run.
//
// To run a subset of datasets or sizes, use Google Benchmark's filter flag,
// e.g. `--benchmark_filter="BM_Encode/Uniform256"`.
// ---------------------------------------------------------------------------

namespace {

/// @brief RNG seed, overridable via `PIVCO_BENCH_SEED`.
std::uint64_t bench_seed() {
  if (const char* env = std::getenv("PIVCO_BENCH_SEED")) {
    try {
      return std::stoull(env);
    } catch (...) {
    }
  }
  return 42;
}

/// @brief Signature of a dataset generator: produces @p n bytes from @p rng.
using DatasetMaker = std::vector<std::uint8_t> (*)(std::size_t,
                                                   std::mt19937_64&);

// --- built-in datasets -----------------------------------------------------

/// @brief Uniformly distributed bytes over the full 256-symbol alphabet.
std::vector<std::uint8_t> make_uniform256(std::size_t n, std::mt19937_64& rng) {
  std::vector<std::uint8_t> data(n);
  for (auto& b : data) {
    b = static_cast<std::uint8_t>(rng());
  }
  return data;
}

/// @brief Binary alphabet {0, 1}: shortest possible codes (depth 1).
std::vector<std::uint8_t> make_binary(std::size_t n, std::mt19937_64& rng) {
  std::vector<std::uint8_t> data(n);
  for (auto& b : data) {
    b = static_cast<std::uint8_t>(rng() & 1u);
  }
  return data;
}

/// @brief Four-symbol alphabet {0,1,2,3}: code depth 2.
std::vector<std::uint8_t> make_low_entropy_4(std::size_t n,
                                             std::mt19937_64& rng) {
  std::vector<std::uint8_t> data(n);
  for (auto& b : data) {
    b = static_cast<std::uint8_t>(rng() & 3u);
  }
  return data;
}

/// @brief English-text-like distribution: ~17% space, ~80% lowercase letters,
///        ~2% digits, ~1% punctuation.
std::vector<std::uint8_t> make_text_like(std::size_t n, std::mt19937_64& rng) {
  std::vector<std::uint8_t> data(n);
  std::uniform_int_distribution<int> pick(0, 99);
  std::uniform_int_distribution<int> letter('a', 'z');
  std::uniform_int_distribution<int> digit('0', '9');
  std::uniform_int_distribution<int> punct('!', '/');
  for (auto& b : data) {
    const int r = pick(rng);
    if (r < 17) {
      b = ' ';
    } else if (r < 97) {
      b = static_cast<std::uint8_t>(letter(rng));
    } else if (r < 99) {
      b = static_cast<std::uint8_t>(digit(rng));
    } else {
      b = static_cast<std::uint8_t>(punct(rng));
    }
  }
  return data;
}

/// @brief Highly skewed: 99% one symbol, 1% random over 256. Exercises a deep,
///        lopsided Huffman tree.
std::vector<std::uint8_t> make_skewed99(std::size_t n, std::mt19937_64& rng) {
  std::vector<std::uint8_t> data(n, 0);
  std::uniform_int_distribution<int> coin(0, 99);
  for (auto& b : data) {
    if (coin(rng) < 1) {
      b = static_cast<std::uint8_t>(rng());
    }
  }
  return data;
}

/// @brief Degenerate single-symbol input (root is a leaf): best-case path.
std::vector<std::uint8_t> make_single_symbol(std::size_t n,
                                             std::mt19937_64& rng) {
  (void)rng;
  return std::vector<std::uint8_t>(n, 123);
}

/// @brief Compression efficiency in bits per input byte (8 = no compression).
void report_ratio(benchmark::State& state,
                  std::size_t compressed_bytes,
                  std::size_t input_bytes) {
  state.counters["bpb"] = static_cast<double>(compressed_bytes) * 8.0 /
                          static_cast<double>(input_bytes);
}

/// @brief Common range of input sizes: 1 KiB, 1 MiB, 1 GiB.
/// @details Applied to every registered benchmark. The 1 GiB tier exercises
///         gigabyte-scale inputs; at that size Google Benchmark runs a single
///         iteration per case.
constexpr std::int64_t kSizeLo = 1ull << 10;  // 1 KiB
constexpr std::int64_t kSizeHi = 1ull << 30;  // 1 GiB
constexpr std::int64_t kSizeMult = 1024;

// --- benchmarks ------------------------------------------------------------

/// @brief Measure full encoding: build Huffman tree, fill node bitmaps, and
///        serialize the compressed stream. The compressed size is read from
///        the first timed iteration (no separate probe build), which keeps
///        gigabyte-scale setup cost to a single encode.
void BM_Encode(benchmark::State& state, DatasetMaker make) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  std::mt19937_64 rng(bench_seed());
  const std::vector<std::uint8_t> data = make(size, rng);

  std::size_t compressed = 0;
  for (auto _ : state) {
    PivCoHuffman codec(data);
    compressed = codec.compressed_size();
    benchmark::DoNotOptimize(codec.compressed_data().data());
    benchmark::ClobberMemory();
  }
  if (compressed > 0) {
    report_ratio(state, compressed, size);
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations() * size));
}

/// @brief Measure decoding: reconstruct the original byte stream from the
///        in-memory tree. Encoding is done once, outside the timed loop.
void BM_Decode(benchmark::State& state, DatasetMaker make) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  std::mt19937_64 rng(bench_seed());
  const std::vector<std::uint8_t> data = make(size, rng);

  const PivCoHuffman codec(data);
  report_ratio(state, codec.compressed_size(), size);

  for (auto _ : state) {
    std::vector<std::uint8_t> out = codec.decode();
    benchmark::DoNotOptimize(out.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations() * size));
}

}  // namespace

// Register encode + decode for every built-in dataset, each swept over the
// full size range. Each registration is named `BM_<Op>/<Dataset>/<size>`, so
// `--benchmark_filter` selects subsets freely.
#define PIXIE_PIVCO_REGISTER(Op, Dataset, Maker) \
  BENCHMARK_CAPTURE(BM_##Op, Dataset, Maker)     \
      ->Unit(benchmark::kMillisecond)            \
      ->UseRealTime()                            \
      ->Repetitions(5)                           \
      ->ReportAggregatesOnly(true)               \
      ->RangeMultiplier(kSizeMult)               \
      ->Range(kSizeLo, kSizeHi);

PIXIE_PIVCO_REGISTER(Encode, Uniform256, make_uniform256)
PIXIE_PIVCO_REGISTER(Encode, Binary, make_binary)
PIXIE_PIVCO_REGISTER(Encode, LowEntropy4, make_low_entropy_4)
PIXIE_PIVCO_REGISTER(Encode, TextLike, make_text_like)
PIXIE_PIVCO_REGISTER(Encode, Skewed99, make_skewed99)
PIXIE_PIVCO_REGISTER(Encode, SingleSymbol, make_single_symbol)

PIXIE_PIVCO_REGISTER(Decode, Uniform256, make_uniform256)
PIXIE_PIVCO_REGISTER(Decode, Binary, make_binary)
PIXIE_PIVCO_REGISTER(Decode, LowEntropy4, make_low_entropy_4)
PIXIE_PIVCO_REGISTER(Decode, TextLike, make_text_like)
PIXIE_PIVCO_REGISTER(Decode, Skewed99, make_skewed99)
PIXIE_PIVCO_REGISTER(Decode, SingleSymbol, make_single_symbol)
