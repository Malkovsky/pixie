#include <benchmark/benchmark.h>
#include <pixie/file_archive/implementations.h>
#include <pixie/serialization.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::vector<std::byte> Bytes(std::string_view value) {
  const auto bytes = std::as_bytes(std::span(value.data(), value.size()));
  return {bytes.begin(), bytes.end()};
}

std::vector<pixie::FileArchiveSource> MakeFiles(std::size_t count) {
  std::vector<pixie::FileArchiveSource> files;
  files.reserve(count);
  const std::vector<std::byte> content = Bytes(std::string(255, 'x') + "\n");
  for (std::size_t index = 0; index < count; ++index) {
    files.push_back(
        {.path = "src/file-" + std::to_string(index), .content = content});
  }
  return files;
}

pixie::FileArchiveView MakeView(const pixie::FileArchive& archive,
                                std::vector<std::byte>& storage) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  archive.serialize(writer);
  writer.finish();
  storage = sink.take();
  pixie::BinaryReader reader(storage);
  return pixie::FileArchiveView::deserialize(reader);
}

void BM_FileArchiveFind(benchmark::State& state) {
  const std::size_t count = static_cast<std::size_t>(state.range(0));
  const pixie::FileArchive archive(MakeFiles(count));
  std::vector<std::byte> storage;
  const pixie::FileArchiveView view = MakeView(archive, storage);
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
  const pixie::FileArchive archive(
      {{.path = "source.cpp", .content = Bytes(content)}});
  std::vector<std::byte> storage;
  const pixie::FileArchiveView view = MakeView(archive, storage);

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

BENCHMARK(BM_FileArchiveFind)->RangeMultiplier(8)->Range(64, 1U << 18U);
BENCHMARK(BM_FileArchiveExtractLines)->Arg(1)->Arg(8)->Arg(64)->Arg(512);

}  // namespace

BENCHMARK_MAIN();
