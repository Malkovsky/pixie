#include <benchmark/benchmark.h>
#include <pixie/huffman/implementations.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

using pixie::PivCoHuffman;

// ---------------------------------------------------------------------------
// File-based PivCo-Huffman benchmark.
//
// Reads a real-world text file (default: prose_pride.txt, "Pride and
// Prejudice" from Project Gutenberg) and measures encode/decode throughput
// on it. This matches the article's evaluation methodology, which uses real
// datasets rather than synthetic random data.
//
// Configuration:
//   PIVCO_BENCH_FILE - path to the input file (default: prose_pride.txt)
// ---------------------------------------------------------------------------

namespace {

/// @brief Read the entire contents of @p path into a byte vector.
std::vector<std::uint8_t> read_file(const std::string& path) {
  std::ifstream file(path, std::ios::binary | std::ios::ate);
  if (!file) {
    std::fprintf(stderr, "Cannot open '%s'\n", path.c_str());
    std::exit(1);
  }
  const std::streamsize size = file.tellg();
  file.seekg(0);
  std::vector<std::uint8_t> data(static_cast<std::size_t>(size));
  if (size > 0) {
    file.read(reinterpret_cast<char*>(data.data()), size);
  }
  return data;
}

/// @brief Input file path, overridable via `PIVCO_BENCH_FILE`.
const std::string& bench_file() {
  static const std::string path = [] {
    if (const char* env = std::getenv("PIVCO_BENCH_FILE")) {
      return std::string(env);
    }
    return std::string("prose_pride.txt");
  }();
  return path;
}

/// @brief Compression efficiency in bits per input byte.
void report_ratio(benchmark::State& state,
                  std::size_t compressed_bytes,
                  std::size_t input_bytes) {
  state.counters["bpb"] = static_cast<double>(compressed_bytes) * 8.0 /
                          static_cast<double>(input_bytes);
}

/// @brief Compute Shannon entropy of the data in bits per byte, for reference.
double entropy(std::span<const std::uint8_t> data) {
  if (data.empty()) {
    return 0.0;
  }
  std::array<std::size_t, 256> freq{};
  for (auto b : data) {
    freq[b]++;
  }
  double h = 0.0;
  const double n = static_cast<double>(data.size());
  for (std::size_t i = 0; i < 256; i++) {
    if (freq[i] > 0) {
      const double p = static_cast<double>(freq[i]) / n;
      h -= p * std::log2(p);
    }
  }
  return h;
}

// --- benchmarks ------------------------------------------------------------

/// @brief Measure full encoding of the file.
void BM_EncodeFile(benchmark::State& state) {
  const std::vector<std::uint8_t> data = read_file(bench_file());
  const std::size_t size = data.size();

  // Report entropy once as a reference counter.
  state.counters["entropy_bpb"] = entropy(data);

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

/// @brief Measure decoding of the file. Encoding is done once, outside timing.
void BM_DecodeFile(benchmark::State& state) {
  const std::vector<std::uint8_t> data = read_file(bench_file());
  const std::size_t size = data.size();

  state.counters["entropy_bpb"] = entropy(data);

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

BENCHMARK(BM_EncodeFile)
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Repetitions(10)
    ->ReportAggregatesOnly(true);

BENCHMARK(BM_DecodeFile)
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Repetitions(10)
    ->ReportAggregatesOnly(true);
