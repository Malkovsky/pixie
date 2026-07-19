#include <benchmark/benchmark.h>
#include <pixie/huffman/implementations.h>

#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using pixie::PivCoHuffman;

// ---------------------------------------------------------------------------
// Configuration via environment variables (follows the Pixie convention used
// by EXCESS_POS_CASES / RECORD_LOWS_CASES etc.):
//   PIVCO_BENCH_SIZE - input size in bytes (default 1<<16 = 64 KiB)
//   PIVCO_BENCH_SEED - random seed (default 42)
//
// To run a subset of datasets, use Google Benchmark's built-in filter flag,
// e.g. `--benchmark_filter=".*Skewed.*|.*Text.*"`. No code change needed.
// ---------------------------------------------------------------------------

namespace {

/// @brief Input size in bytes, overridable via `PIVCO_BENCH_SIZE`.
std::size_t bench_size() {
  if (const char* env = std::getenv("PIVCO_BENCH_SIZE")) {
    try {
      const unsigned long long v = std::stoull(env);
      if (v != 0) {
        return static_cast<std::size_t>(v);
      }
    } catch (...) {
    }
  }
  return 1ull << 16;
}

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

// --- benchmarks ------------------------------------------------------------

/// @brief Measure full encoding: build Huffman tree, fill node bitmaps, and
///        serialize the compressed stream.
void BM_Encode(benchmark::State& state, DatasetMaker make) {
  const std::size_t size = bench_size();
  std::mt19937_64 rng(bench_seed());
  const std::vector<std::uint8_t> data = make(size, rng);

  // Encode once outside timing to report compression ratio.
  const PivCoHuffman probe(data);
  report_ratio(state, probe.compressed_size(), size);

  for (auto _ : state) {
    PivCoHuffman codec(data);
    benchmark::DoNotOptimize(codec.compressed_data().data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(static_cast<int64_t>(state.iterations() * size));
}

/// @brief Measure decoding: reconstruct the original byte stream from the
///        in-memory tree. Encoding is done once, outside the timed loop.
void BM_Decode(benchmark::State& state, DatasetMaker make) {
  const std::size_t size = bench_size();
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

// Register encode + decode for every built-in dataset. Each registration is
// named `BM_<Op>/<Dataset>`, so `--benchmark_filter` selects subsets freely.
BENCHMARK_CAPTURE(BM_Encode, Uniform256, make_uniform256);
BENCHMARK_CAPTURE(BM_Encode, Binary, make_binary);
BENCHMARK_CAPTURE(BM_Encode, LowEntropy4, make_low_entropy_4);
BENCHMARK_CAPTURE(BM_Encode, TextLike, make_text_like);
BENCHMARK_CAPTURE(BM_Encode, Skewed99, make_skewed99);
BENCHMARK_CAPTURE(BM_Encode, SingleSymbol, make_single_symbol);

BENCHMARK_CAPTURE(BM_Decode, Uniform256, make_uniform256);
BENCHMARK_CAPTURE(BM_Decode, Binary, make_binary);
BENCHMARK_CAPTURE(BM_Decode, LowEntropy4, make_low_entropy_4);
BENCHMARK_CAPTURE(BM_Decode, TextLike, make_text_like);
BENCHMARK_CAPTURE(BM_Decode, Skewed99, make_skewed99);
BENCHMARK_CAPTURE(BM_Decode, SingleSymbol, make_single_symbol);
