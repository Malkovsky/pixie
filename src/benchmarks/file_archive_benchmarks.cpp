#include <benchmark/benchmark.h>
#include <pixie/file_archive/implementations.h>
#include <pixie/serialization.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

std::vector<std::byte> Bytes(std::string_view value) {
  const auto bytes = std::as_bytes(std::span(value.data(), value.size()));
  return {bytes.begin(), bytes.end()};
}

constexpr std::size_t kBuildFileBytes = 16 * 1024;

std::vector<std::byte> MakeContent(std::size_t size) {
  std::string content;
  content.reserve(size);
  while (content.size() < size) {
    content.append("The quick brown fox jumps over the lazy dog.\n");
  }
  content.resize(size, '\n');
  return Bytes(content);
}

std::vector<pixie::FileArchiveSource> MakeFiles(
    std::size_t count,
    std::size_t bytes_per_file = 256) {
  std::vector<pixie::FileArchiveSource> files;
  files.reserve(count);
  const std::vector<std::byte> content = MakeContent(bytes_per_file);
  for (std::size_t index = 0; index < count; ++index) {
    files.push_back(
        {.path = "src/file-" + std::to_string(index), .content = content});
  }
  return files;
}

std::vector<std::byte> Serialize(const pixie::FileArchive& archive) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  archive.serialize(writer);
  writer.finish();
  return sink.take();
}

pixie::FileArchiveView MakeView(std::vector<pixie::FileArchiveSource> sources,
                                std::vector<std::byte>& storage) {
  const pixie::FileArchive archive(std::move(sources));
  storage = Serialize(archive);
  pixie::BinaryReader reader(storage);
  return pixie::FileArchiveView::deserialize(reader);
}

pixie::WaveletTreeBuildType BuildType(const benchmark::State& state) {
  return state.range(1) == 0 ? pixie::WaveletTreeBuildType::Standard
                             : pixie::WaveletTreeBuildType::Huffman;
}

void SetArchiveSizeCounters(benchmark::State& state,
                            std::size_t input_bytes,
                            std::size_t serialized_bytes) {
  state.counters["serialized_bytes"] = static_cast<double>(serialized_bytes);
  state.counters["serialized_bytes_per_input_byte"] =
      static_cast<double>(serialized_bytes) / input_bytes;
}

struct StreamedFiles {
  std::vector<pixie::FileArchiveSourceMetadata> metadata;
  std::vector<std::vector<std::byte>> content;
};

StreamedFiles MakeStreamedFiles(std::size_t count) {
  StreamedFiles result;
  result.metadata.reserve(count);
  result.content.reserve(count);
  for (std::size_t index = 0; index < count; ++index) {
    result.metadata.push_back({.path = "src/file-" + std::to_string(index)});
    result.content.push_back(MakeContent(kBuildFileBytes));
  }
  return result;
}

template <class Consume>
void ReadStreamedFile(const StreamedFiles& files,
                      std::size_t source_index,
                      Consume&& consume) {
  constexpr std::size_t kChunkBytes = 4096;
  const auto& content = files.content[source_index];
  for (std::size_t offset = 0; offset < content.size(); offset += kChunkBytes) {
    consume(std::span(content).subspan(
        offset, std::min(kChunkBytes, content.size() - offset)));
  }
}

void BM_FileArchiveBuild(benchmark::State& state) {
  const std::size_t count = static_cast<std::size_t>(state.range(0));
  const auto build_type = BuildType(state);
  const std::vector<pixie::FileArchiveSource> source_files =
      MakeFiles(count, kBuildFileBytes);
  const std::size_t input_bytes = count * kBuildFileBytes;
  const std::size_t serialized_bytes = [&] {
    const pixie::FileArchive reference(source_files, build_type);
    return Serialize(reference).size();
  }();

  for (auto _ : state) {
    state.PauseTiming();
    std::vector<pixie::FileArchiveSource> files = source_files;
    state.ResumeTiming();
    std::optional<pixie::FileArchive> archive(std::in_place, std::move(files),
                                              build_type);
    benchmark::DoNotOptimize(archive->logical_size_bytes());
    state.PauseTiming();
    archive.reset();
    state.ResumeTiming();
  }
  state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(input_bytes));
  SetArchiveSizeCounters(state, input_bytes, serialized_bytes);
}

void BM_FileArchiveBuildStreaming(benchmark::State& state) {
  const std::size_t count = static_cast<std::size_t>(state.range(0));
  const auto build_type = BuildType(state);
  const StreamedFiles files = MakeStreamedFiles(count);
  const std::size_t input_bytes = count * kBuildFileBytes;
  const std::size_t serialized_bytes = [&] {
    std::size_t source_index = 0;
    const pixie::FileArchive reference(
        files.metadata,
        [&](const auto&, auto&& consume) {
          ReadStreamedFile(files, source_index++ % count,
                           std::forward<decltype(consume)>(consume));
        },
        build_type);
    return Serialize(reference).size();
  }();

  for (auto _ : state) {
    state.PauseTiming();
    std::vector<pixie::FileArchiveSourceMetadata> metadata = files.metadata;
    state.ResumeTiming();
    std::size_t source_index = 0;
    std::optional<pixie::FileArchive> archive(
        std::in_place, std::move(metadata),
        [&](const auto&, auto&& consume) {
          ReadStreamedFile(files, source_index++ % count,
                           std::forward<decltype(consume)>(consume));
        },
        build_type);
    benchmark::DoNotOptimize(archive->logical_size_bytes());
    state.PauseTiming();
    archive.reset();
    state.ResumeTiming();
  }
  state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(input_bytes));
  SetArchiveSizeCounters(state, input_bytes, serialized_bytes);
}

void BM_FileArchiveFind(benchmark::State& state) {
  const std::size_t count = static_cast<std::size_t>(state.range(0));
  std::vector<std::byte> storage;
  const pixie::FileArchiveView view = MakeView(MakeFiles(count), storage);
  std::vector<std::string> paths;
  paths.reserve(count);
  for (std::size_t index = 0; index < count; ++index) {
    paths.push_back("src/file-" + std::to_string(index));
  }

  std::size_t query = 0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(view.find(paths[query]));
    query = (query + 7919) % count;
  }
  state.SetItemsProcessed(state.iterations());
}

void BM_FileArchiveExtractLines(benchmark::State& state) {
  std::string content;
  for (std::size_t line = 0; line < 1U << 16U; ++line) {
    content.append("0123456789abcdef\n");
  }
  std::vector<std::byte> storage;
  const pixie::FileArchiveView view =
      MakeView({{.path = "source.cpp", .content = Bytes(content)}}, storage);

  std::size_t begin = 0;
  const std::size_t line_count = view.entry(0).line_count;
  const std::size_t extracted_lines = static_cast<std::size_t>(state.range(0));
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        view.extract_lines(0, begin, begin + extracted_lines));
    begin = (begin + 7919) % (line_count - extracted_lines + 1);
  }
  state.SetItemsProcessed(state.iterations() *
                          static_cast<std::int64_t>(extracted_lines));
}

void BM_FileArchiveExtract(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  std::vector<std::byte> storage;
  const pixie::FileArchiveView view =
      MakeView({{.path = "source.txt", .content = MakeContent(size)}}, storage);

  for (auto _ : state) {
    benchmark::DoNotOptimize(view.extract(0));
  }
  state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(size));
}

void BM_FileArchiveViewDeserialize(benchmark::State& state) {
  const std::size_t count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::byte> serialized = [&] {
    const pixie::FileArchive archive(MakeFiles(count, kBuildFileBytes));
    return Serialize(archive);
  }();

  for (auto _ : state) {
    pixie::BinaryReader reader(serialized);
    const pixie::FileArchiveView view = pixie::FileArchiveView::deserialize(
        reader, pixie::DeserializationValidation::kQuick);
    benchmark::DoNotOptimize(view.size());
  }
  state.SetItemsProcessed(state.iterations());
  state.counters["artifact_bytes"] = static_cast<double>(serialized.size());
}

BENCHMARK(BM_FileArchiveBuild)
    ->Args({16, 0})
    ->Args({16, 1})
    ->Args({256, 0})
    ->Args({256, 1})
    ->ArgNames({"files", "build_type"});
BENCHMARK(BM_FileArchiveBuildStreaming)
    ->Args({16, 0})
    ->Args({16, 1})
    ->Args({256, 0})
    ->Args({256, 1})
    ->ArgNames({"files", "build_type"});
BENCHMARK(BM_FileArchiveFind)->RangeMultiplier(8)->Range(64, 1U << 18U);
BENCHMARK(BM_FileArchiveExtractLines)->Arg(1)->Arg(8)->Arg(64)->Arg(512);
BENCHMARK(BM_FileArchiveExtract)
    ->Arg(4 * 1024)
    ->Arg(64 * 1024)
    ->Arg(1024 * 1024);
BENCHMARK(BM_FileArchiveViewDeserialize)->Arg(16)->Arg(256);

}  // namespace

BENCHMARK_MAIN();
