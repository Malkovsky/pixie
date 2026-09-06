#include <benchmark/benchmark.h>
#include <pixie/huffman/implementations.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <span>
#include <string>
#include <vector>

namespace {

using pixie::PivCoHuffman;

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

const std::string& benchmark_file() {
  static const std::string path = [] {
    if (const char* value = std::getenv("PIVCO_BENCH_FILE")) {
      return std::string(value);
    }
    return std::string("prose_pride.txt");
  }();
  return path;
}

double entropy(std::span<const std::uint8_t> data) {
  if (data.empty()) {
    return 0.0;
  }
  std::array<std::size_t, 256> frequencies{};
  for (std::uint8_t symbol : data) {
    ++frequencies[symbol];
  }
  double result = 0.0;
  const double count = static_cast<double>(data.size());
  for (std::size_t frequency : frequencies) {
    if (frequency != 0) {
      const double probability = static_cast<double>(frequency) / count;
      result -= probability * std::log2(probability);
    }
  }
  return result;
}

double bits_per_symbol(std::size_t compressed_size,
                       std::size_t uncompressed_size) {
  if (uncompressed_size == 0) {
    return 0.0;
  }
  return static_cast<double>(compressed_size) * 8.0 /
         static_cast<double>(uncompressed_size);
}

template <class Codec>
void encode_file(benchmark::State& state) {
  const std::vector<std::uint8_t> data = read_file(benchmark_file());
  state.counters["entropy_bpb"] = entropy(data);
  std::size_t compressed_size = 0;
  for (auto _ : state) {
    Codec codec(data);
    compressed_size = codec.compressed_size();
    benchmark::DoNotOptimize(codec.compressed_data().data());
    benchmark::ClobberMemory();
  }
  state.counters["bpb"] = bits_per_symbol(compressed_size, data.size());
  state.SetBytesProcessed(
      static_cast<std::int64_t>(state.iterations() * data.size()));
}

template <class Codec>
void decode_file(benchmark::State& state) {
  const std::vector<std::uint8_t> data = read_file(benchmark_file());
  const Codec codec(data);
  state.counters["entropy_bpb"] = entropy(data);
  state.counters["bpb"] = bits_per_symbol(codec.compressed_size(), data.size());
  for (auto _ : state) {
    std::vector<std::uint8_t> decoded = codec.decode();
    benchmark::DoNotOptimize(decoded.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(
      static_cast<std::int64_t>(state.iterations() * data.size()));
}

BENCHMARK_TEMPLATE(encode_file, PivCoHuffman)
    ->Name("PivCoHuffman/EncodeFile")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Repetitions(10)
    ->ReportAggregatesOnly(true);
BENCHMARK_TEMPLATE(decode_file, PivCoHuffman)
    ->Name("PivCoHuffman/DecodeFile")
    ->Unit(benchmark::kMillisecond)
    ->UseRealTime()
    ->Repetitions(10)
    ->ReportAggregatesOnly(true);

}  // namespace

BENCHMARK_MAIN();
