#include <benchmark/benchmark.h>
#include <pixie/integer_vector/implementations.h>
#include <pixie/io/file_output_sink.h>
#include <pixie/io/mapped_file.h>
#include <pixie/rank_select/implementations.h>
#include <pixie/rmm/implementations.h>
#include <pixie/rmq/implementations.h>
#include <pixie/serialization.h>
#include <pixie/split_span.h>
#include <pixie/wavelet_tree/implementations.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kKiB = 1024;
constexpr std::size_t kMiB = 1024 * kKiB;
constexpr std::size_t kRecordBytes = 8 * sizeof(std::uint64_t);
constexpr std::size_t kDefaultStagingBytes = 64 * kKiB;
constexpr std::size_t kWaveletAlphabetSize = 256;
constexpr double kBenchmarkWarmupSeconds = 0.1;
constexpr double kBenchmarkMinSeconds = 0.5;

std::span<std::byte> as_writable_span(std::vector<std::byte>& bytes) {
  return {bytes.data(), bytes.size()};
}

std::span<const std::byte> as_const_span(const std::vector<std::byte>& bytes) {
  return {bytes.data(), bytes.size()};
}

std::uint64_t mix(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

std::vector<std::byte> make_payload(std::size_t size) {
  std::vector<std::byte> result(size);
  for (std::size_t i = 0; i < size; ++i) {
    result[i] = static_cast<std::byte>(mix(i) & 0xff);
  }
  return result;
}

std::vector<std::uint64_t> make_words(std::size_t bit_count) {
  const std::size_t word_count = bit_count == 0 ? 0 : 1 + (bit_count - 1) / 64;
  std::vector<std::uint64_t> result(word_count);
  for (std::size_t i = 0; i < word_count; ++i) {
    result[i] = mix(i);
  }
  return result;
}

std::vector<std::int64_t> make_values(std::size_t count) {
  std::vector<std::int64_t> result(count);
  for (std::size_t i = 0; i < count; ++i) {
    result[i] = static_cast<std::int64_t>(mix(i) & 0x7fffffffffffffffULL);
  }
  return result;
}

std::vector<std::uint64_t> make_symbols(std::size_t count) {
  std::vector<std::uint64_t> result(count);
  for (std::size_t i = 0; i < count; ++i) {
    result[i] = mix(i) % kWaveletAlphabetSize;
  }
  return result;
}

template <class Serializable>
std::vector<std::byte> serialize_to_vector(const Serializable& value) {
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  value.serialize(writer);
  writer.finish();
  return sink.take();
}

class AlignedArtifact {
 public:
  explicit AlignedArtifact(std::span<const std::byte> bytes)
      : words_((bytes.size() + sizeof(std::uint64_t) - 1) /
               sizeof(std::uint64_t)),
        size_bytes_(bytes.size()) {
    std::ranges::copy(bytes, writable_bytes().begin());
  }

  std::span<const std::byte> bytes() const {
    return std::as_bytes(std::span<const std::uint64_t>(words_))
        .first(size_bytes_);
  }

 private:
  std::span<std::byte> writable_bytes() {
    return std::as_writable_bytes(std::span<std::uint64_t>(words_));
  }

  std::vector<std::uint64_t> words_;
  std::size_t size_bytes_;
};

class ScopedTemporaryPath {
 public:
  explicit ScopedTemporaryPath(std::string filename)
      : path_(std::filesystem::temp_directory_path() / std::move(filename)) {
    remove();
  }

  ScopedTemporaryPath(const ScopedTemporaryPath&) = delete;
  ScopedTemporaryPath& operator=(const ScopedTemporaryPath&) = delete;

  ~ScopedTemporaryPath() { remove(); }

  const std::filesystem::path& path() const { return path_; }

 private:
  void remove() noexcept {
    std::error_code error;
    std::filesystem::remove(path_, error);
  }

  std::filesystem::path path_;
};

void set_throughput(benchmark::State& state, std::size_t bytes_per_iteration) {
  const auto iterations = static_cast<std::uint64_t>(state.iterations());
  const auto bytes = static_cast<std::uint64_t>(bytes_per_iteration);
  state.SetBytesProcessed(static_cast<std::int64_t>(iterations * bytes));
}

void set_artifact_counters(benchmark::State& state,
                           std::size_t item_count,
                           std::size_t artifact_bytes) {
  set_throughput(state, artifact_bytes);
  state.counters["artifact_bytes"] =
      benchmark::Counter(static_cast<double>(artifact_bytes));
  state.counters["items"] = benchmark::Counter(static_cast<double>(item_count));
}

template <class Serializable>
void serialize_iterations(benchmark::State& state,
                          const Serializable& value,
                          std::span<std::byte> destination,
                          std::span<std::byte> staging) {
  for (auto _ : state) {
    (void)_;
    pixie::SpanOutputSink sink(destination);
    pixie::BinaryWriter writer(sink, staging);
    value.serialize(writer);
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
    benchmark::ClobberMemory();
  }
}

template <class Deserialize>
void deserialize_iterations(benchmark::State& state,
                            std::span<const std::byte> artifact,
                            Deserialize deserialize) {
  for (auto _ : state) {
    (void)_;
    pixie::BinaryReader reader(artifact);
    auto value = deserialize(reader);
    benchmark::DoNotOptimize(value);
    benchmark::DoNotOptimize(reader.position());
  }
}

std::vector<std::byte> make_framed_artifact(std::size_t payload_bytes) {
  const std::size_t record_count = payload_bytes / kRecordBytes;
  pixie::VectorOutputSink sink;
  pixie::BinaryWriter writer(sink);
  for (std::size_t record = 0; record < record_count; ++record) {
    const std::size_t record_begin = writer.size_bytes();
    const std::size_t size_position = writer.write_u64_placeholder();
    for (std::size_t field = 0; field < 7; ++field) {
      writer.write_u64(mix(record * 7 + field));
    }
    writer.patch_u64(size_position, writer.size_bytes() - record_begin);
  }
  writer.finish();
  return sink.take();
}

void BM_BinaryWriterU64Span(benchmark::State& state) {
  const std::size_t requested_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t payload_bytes =
      requested_bytes - requested_bytes % sizeof(std::uint64_t);
  const std::size_t field_count = payload_bytes / sizeof(std::uint64_t);
  const std::size_t buffer_bytes = static_cast<std::size_t>(state.range(1));
  std::vector<std::uint64_t> fields(field_count);
  for (std::size_t field = 0; field < field_count; ++field) {
    fields[field] = mix(field);
  }
  std::vector<std::byte> destination(payload_bytes);
  std::vector<std::byte> staging(buffer_bytes);
  for (auto _ : state) {
    (void)_;
    pixie::SpanOutputSink sink(as_writable_span(destination));
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    for (std::size_t field = 0; field < field_count; ++field) {
      writer.write_u64(fields[field]);
    }
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
    benchmark::ClobberMemory();
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryReaderU64(benchmark::State& state) {
  const std::size_t requested_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t payload_bytes =
      requested_bytes - requested_bytes % sizeof(std::uint64_t);
  const std::size_t field_count = payload_bytes / sizeof(std::uint64_t);
  const std::vector<std::byte> source = make_framed_artifact(payload_bytes);
  for (auto _ : state) {
    (void)_;
    pixie::BinaryReader reader(as_const_span(source));
    std::uint64_t checksum = 0;
    for (std::size_t field = 0; field < field_count; ++field) {
      checksum ^= reader.read_u64();
    }
    benchmark::DoNotOptimize(checksum);
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryWriterBytesSpan(benchmark::State& state) {
  const std::size_t payload_bytes = static_cast<std::size_t>(state.range(0));
  const std::vector<std::byte> source = make_payload(payload_bytes);
  std::vector<std::byte> destination(payload_bytes);
  std::vector<std::byte> staging(kDefaultStagingBytes);
  for (auto _ : state) {
    (void)_;
    pixie::SpanOutputSink sink(as_writable_span(destination));
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    writer.write_bytes(as_const_span(source));
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
    benchmark::ClobberMemory();
  }
  set_throughput(state, payload_bytes);
}

template <bool Split>
void BM_BinaryWriterBytesSplitSpanImpl(benchmark::State& state) {
  const std::size_t payload_bytes = static_cast<std::size_t>(state.range(0));
  const std::vector<std::byte> source = make_payload(payload_bytes);
  std::vector<std::byte> destination(payload_bytes);
  std::vector<std::byte> staging(kDefaultStagingBytes);
  const std::size_t split_position = Split ? payload_bytes / 2 : payload_bytes;
  for (auto _ : state) {
    (void)_;
    pixie::SpanOutputSink sink(as_writable_span(destination));
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    const pixie::SplitSpan<const std::byte> segments(
        as_const_span(source).first(split_position),
        as_const_span(source).subspan(split_position));
    for (const std::span<const std::byte> segment : segments) {
      writer.write_bytes(segment);
    }
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
    benchmark::ClobberMemory();
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryWriterBytesSplitSpanOne(benchmark::State& state) {
  BM_BinaryWriterBytesSplitSpanImpl<false>(state);
}

void BM_BinaryWriterBytesSplitSpanTwo(benchmark::State& state) {
  BM_BinaryWriterBytesSplitSpanImpl<true>(state);
}

void BM_SpanDescriptor(benchmark::State& state) {
  const std::vector<std::byte> source = make_payload(64);
  for (auto _ : state) {
    (void)_;
    const std::span<const std::byte> segment = as_const_span(source);
    benchmark::DoNotOptimize(segment.data());
    benchmark::DoNotOptimize(segment.size());
  }
}

template <bool Split>
void BM_SplitSpanDescriptorImpl(benchmark::State& state) {
  const std::vector<std::byte> source = make_payload(64);
  constexpr std::size_t kSplitPosition = Split ? 32 : 64;
  for (auto _ : state) {
    (void)_;
    const pixie::SplitSpan<const std::byte> segments(
        as_const_span(source).first(kSplitPosition),
        as_const_span(source).subspan(kSplitPosition));
    std::size_t size = 0;
    for (const std::span<const std::byte> segment : segments) {
      benchmark::DoNotOptimize(segment.data());
      size += segment.size();
    }
    benchmark::DoNotOptimize(size);
  }
}

void BM_SplitSpanDescriptorOne(benchmark::State& state) {
  BM_SplitSpanDescriptorImpl<false>(state);
}

void BM_SplitSpanDescriptorTwo(benchmark::State& state) {
  BM_SplitSpanDescriptorImpl<true>(state);
}

void BM_BinaryWriterBytesVector(benchmark::State& state) {
  const std::size_t payload_bytes = static_cast<std::size_t>(state.range(0));
  const std::vector<std::byte> source = make_payload(payload_bytes);
  std::vector<std::byte> staging(kDefaultStagingBytes);
  for (auto _ : state) {
    (void)_;
    pixie::VectorOutputSink sink;
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    writer.write_bytes(as_const_span(source));
    writer.finish();
    benchmark::DoNotOptimize(sink.bytes().data());
    benchmark::DoNotOptimize(sink.size_bytes());
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryReaderByteSpans(benchmark::State& state) {
  const std::size_t payload_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t chunk_bytes = static_cast<std::size_t>(state.range(1));
  const std::vector<std::byte> source = make_payload(payload_bytes);
  for (auto _ : state) {
    (void)_;
    pixie::BinaryReader reader(as_const_span(source));
    std::uint64_t checksum = 0;
    while (!reader.empty()) {
      const std::size_t count = std::min(chunk_bytes, reader.remaining());
      const std::span<const std::byte> chunk = reader.read_bytes(count);
      checksum += std::to_integer<std::uint8_t>(chunk.front());
      checksum += std::to_integer<std::uint8_t>(chunk.back());
    }
    benchmark::DoNotOptimize(checksum);
  }
  set_throughput(state, payload_bytes);
  const std::size_t spans_per_iteration =
      (payload_bytes + chunk_bytes - 1) / chunk_bytes;
  state.SetItemsProcessed(static_cast<std::int64_t>(state.iterations()) *
                          static_cast<std::int64_t>(spans_per_iteration));
}

void BM_BinaryWriterFramedSpan(benchmark::State& state) {
  const std::size_t requested_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t payload_bytes =
      requested_bytes - requested_bytes % kRecordBytes;
  const std::size_t record_count = payload_bytes / kRecordBytes;
  const std::size_t buffer_bytes = static_cast<std::size_t>(state.range(1));
  std::vector<std::uint64_t> fields(record_count * 7);
  for (std::size_t field = 0; field < fields.size(); ++field) {
    fields[field] = mix(field);
  }
  std::vector<std::byte> destination(payload_bytes);
  std::vector<std::byte> staging(buffer_bytes);
  for (auto _ : state) {
    (void)_;
    pixie::SpanOutputSink sink(as_writable_span(destination));
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    for (std::size_t record = 0; record < record_count; ++record) {
      const std::size_t record_begin = writer.size_bytes();
      const std::size_t size_position = writer.write_u64_placeholder();
      for (std::size_t field = 0; field < 7; ++field) {
        writer.write_u64(fields[record * 7 + field]);
      }
      writer.patch_u64(size_position, writer.size_bytes() - record_begin);
    }
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
    benchmark::ClobberMemory();
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryReaderFramed(benchmark::State& state) {
  const std::size_t requested_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t payload_bytes =
      requested_bytes - requested_bytes % kRecordBytes;
  const std::vector<std::byte> source = make_framed_artifact(payload_bytes);
  for (auto _ : state) {
    (void)_;
    pixie::BinaryReader reader(as_const_span(source));
    std::uint64_t checksum = 0;
    while (!reader.empty()) {
      const std::size_t record_bytes = reader.read_size();
      pixie::BinaryReader record =
          reader.read_subreader(record_bytes - sizeof(std::uint64_t));
      while (!record.empty()) {
        checksum ^= record.read_u64();
      }
    }
    benchmark::DoNotOptimize(checksum);
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryWriterBytesFile(benchmark::State& state) {
  const std::size_t payload_bytes = static_cast<std::size_t>(state.range(0));
  const std::vector<std::byte> source = make_payload(payload_bytes);
  std::vector<std::byte> staging(kDefaultStagingBytes);
  ScopedTemporaryPath temporary("pixie_serialization_writer_benchmark.bin");
  for (auto _ : state) {
    (void)_;
    pixie::io::FileOutputSink sink(temporary.path());
    pixie::BinaryWriter writer(sink, as_writable_span(staging));
    writer.write_bytes(as_const_span(source));
    writer.finish();
    benchmark::DoNotOptimize(sink.size_bytes());
  }
  set_throughput(state, payload_bytes);
}

void BM_BinaryReaderU64MappedWarm(benchmark::State& state) {
  const std::size_t requested_bytes = static_cast<std::size_t>(state.range(0));
  const std::size_t payload_bytes =
      requested_bytes - requested_bytes % sizeof(std::uint64_t);
  const std::size_t field_count = payload_bytes / sizeof(std::uint64_t);
  const std::vector<std::byte> source = make_framed_artifact(payload_bytes);
  ScopedTemporaryPath temporary("pixie_serialization_reader_benchmark.bin");
  {
    pixie::io::FileOutputSink sink(temporary.path());
    sink.write(as_const_span(source));
    sink.finish();
  }
  pixie::io::MappedFile mapped(temporary.path());
  std::uint64_t warmup = 0;
  constexpr std::size_t kPageBytes = 4096;
  for (std::size_t offset = 0; offset < mapped.size_bytes();
       offset += kPageBytes) {
    warmup += std::to_integer<std::uint8_t>(mapped.as_bytes()[offset]);
  }
  benchmark::DoNotOptimize(warmup);
  for (auto _ : state) {
    (void)_;
    pixie::BinaryReader reader(mapped.as_bytes());
    std::uint64_t checksum = 0;
    for (std::size_t field = 0; field < field_count; ++field) {
      checksum ^= reader.read_u64();
    }
    benchmark::DoNotOptimize(checksum);
  }
  set_throughput(state, payload_bytes);
}

void BM_RankSelectSerialize(benchmark::State& state) {
  const std::size_t bit_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> words = make_words(bit_count);
  const pixie::RankSelectSupport<> index(words, bit_count);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  std::vector<std::byte> destination(artifact.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, index, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, bit_count, artifact.size());
}

template <pixie::DeserializationValidation Validation>
void BM_RankSelectDeserializeOwningImpl(benchmark::State& state) {
  const std::size_t bit_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> words = make_words(bit_count);
  const pixie::RankSelectSupport<> index(words, bit_count);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  deserialize_iterations(state, as_const_span(artifact),
                         [&](pixie::BinaryReader& reader) {
                           return pixie::RankSelectSupport<>::deserialize(
                               reader, words, Validation);
                         });
  set_artifact_counters(state, bit_count, artifact.size());
}

template <pixie::DeserializationValidation Validation>
void BM_RankSelectDeserializeViewImpl(benchmark::State& state) {
  const std::size_t bit_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> words = make_words(bit_count);
  const pixie::RankSelectSupport<> index(words, bit_count);
  const std::vector<std::byte> serialized = serialize_to_vector(index);
  const AlignedArtifact artifact(as_const_span(serialized));
  deserialize_iterations(state, artifact.bytes(),
                         [&](pixie::BinaryReader& reader) {
                           return pixie::RankSelectSupportView::deserialize(
                               reader, words, Validation);
                         });
  set_artifact_counters(state, bit_count, artifact.bytes().size());
}

void BM_RmMSerialize(benchmark::State& state) {
  const std::size_t bit_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> words = make_words(bit_count);
  const pixie::RmMTree index(words, bit_count);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  std::vector<std::byte> destination(artifact.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, index, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, bit_count, artifact.size());
}

template <pixie::DeserializationValidation Validation>
void BM_RmMDeserializeImpl(benchmark::State& state) {
  const std::size_t bit_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> words = make_words(bit_count);
  const pixie::RmMTree index(words, bit_count);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  deserialize_iterations(
      state, as_const_span(artifact), [&](pixie::BinaryReader& reader) {
        return pixie::RmMTree::deserialize(reader, words, Validation);
      });
  set_artifact_counters(state, bit_count, artifact.size());
}

using RmqIndex = pixie::rmq::CartesianHybridBTree<std::int64_t>;

void BM_RmqSerialize(benchmark::State& state) {
  const std::size_t value_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::int64_t> values = make_values(value_count);
  const RmqIndex index(values);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  std::vector<std::byte> destination(artifact.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, index, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, value_count, artifact.size());
}

template <pixie::DeserializationValidation Validation>
void BM_RmqDeserializeImpl(benchmark::State& state) {
  const std::size_t value_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::int64_t> values = make_values(value_count);
  const RmqIndex index(values);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  deserialize_iterations(
      state, as_const_span(artifact), [&](pixie::BinaryReader& reader) {
        return RmqIndex::deserialize(reader, values, Validation);
      });
  set_artifact_counters(state, value_count, artifact.size());
}

void BM_WaveletTreeSerialize(benchmark::State& state) {
  const std::size_t symbol_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> symbols = make_symbols(symbol_count);
  const pixie::WaveletTree<std::uint64_t> index(kWaveletAlphabetSize, symbols);
  const std::vector<std::byte> artifact = serialize_to_vector(index);
  std::vector<std::byte> destination(artifact.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, index, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, symbol_count, artifact.size());
}

std::vector<std::uint64_t> make_integer_vector_values(std::size_t size,
                                                      bool monotone) {
  std::vector<std::uint64_t> values(size);
  for (std::size_t i = 0; i < size; ++i) {
    values[i] = (i * 0x9e3779b97f4a7c15ULL) >> 56;
  }
  if (monotone) {
    std::ranges::sort(values);
  }
  return values;
}

template <class Owner>
Owner make_integer_vector_owner(std::size_t size);

template <>
pixie::PackedIntegerVector<> make_integer_vector_owner(std::size_t size) {
  const auto values = make_integer_vector_values(size, false);
  return pixie::PackedIntegerVector<>(values);
}

template <>
pixie::PackedMonotoneIntegerVector<> make_integer_vector_owner(
    std::size_t size) {
  const auto values = make_integer_vector_values(size, true);
  return pixie::PackedMonotoneIntegerVector<>{
      pixie::PackedIntegerVector<>(values)};
}

template <class Owner>
void integer_vector_serialize_owner(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const Owner owner = make_integer_vector_owner<Owner>(size);
  const std::vector<std::byte> artifact = serialize_to_vector(owner);
  std::vector<std::byte> destination(artifact.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, owner, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, size, artifact.size());
}

template <class Owner, class View>
void integer_vector_serialize_view(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const Owner owner = make_integer_vector_owner<Owner>(size);
  const std::vector<std::byte> serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(as_const_span(serialized));
  pixie::BinaryReader reader(artifact.bytes());
  const View view = View::deserialize(reader);
  std::vector<std::byte> destination(serialized.size());
  std::vector<std::byte> staging(kDefaultStagingBytes);
  serialize_iterations(state, view, as_writable_span(destination),
                       as_writable_span(staging));
  set_artifact_counters(state, size, serialized.size());
}

void BM_PackedIntegerVectorSerializeOwning(benchmark::State& state) {
  integer_vector_serialize_owner<pixie::PackedIntegerVector<>>(state);
}

void BM_PackedIntegerVectorSerializeView(benchmark::State& state) {
  integer_vector_serialize_view<pixie::PackedIntegerVector<>,
                                pixie::PackedIntegerVectorView<>>(state);
}

void BM_PackedMonotoneIntegerVectorSerializeOwning(benchmark::State& state) {
  integer_vector_serialize_owner<pixie::PackedMonotoneIntegerVector<>>(state);
}

void BM_PackedMonotoneIntegerVectorSerializeView(benchmark::State& state) {
  integer_vector_serialize_view<pixie::PackedMonotoneIntegerVector<>,
                                pixie::PackedMonotoneIntegerVectorView<>>(
      state);
}

template <class Owner, pixie::DeserializationValidation Validation>
void integer_vector_deserialize_owner(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const Owner owner = make_integer_vector_owner<Owner>(size);
  const std::vector<std::byte> artifact = serialize_to_vector(owner);
  deserialize_iterations(state, as_const_span(artifact),
                         [](pixie::BinaryReader& reader) {
                           return Owner::deserialize(reader, Validation);
                         });
  set_artifact_counters(state, size, artifact.size());
}

template <class Owner, class View, pixie::DeserializationValidation Validation>
void integer_vector_deserialize_view(benchmark::State& state) {
  const std::size_t size = static_cast<std::size_t>(state.range(0));
  const Owner owner = make_integer_vector_owner<Owner>(size);
  const std::vector<std::byte> serialized = serialize_to_vector(owner);
  const AlignedArtifact artifact(as_const_span(serialized));
  deserialize_iterations(state, artifact.bytes(),
                         [](pixie::BinaryReader& reader) {
                           return View::deserialize(reader, Validation);
                         });
  set_artifact_counters(state, size, artifact.bytes().size());
}

template <pixie::DeserializationValidation Validation>
void BM_PackedIntegerVectorDeserializeOwningImpl(benchmark::State& state) {
  integer_vector_deserialize_owner<pixie::PackedIntegerVector<>, Validation>(
      state);
}

template <pixie::DeserializationValidation Validation>
void BM_PackedIntegerVectorDeserializeViewImpl(benchmark::State& state) {
  integer_vector_deserialize_view<pixie::PackedIntegerVector<>,
                                  pixie::PackedIntegerVectorView<>, Validation>(
      state);
}

template <pixie::DeserializationValidation Validation>
void BM_PackedMonotoneIntegerVectorDeserializeOwningImpl(
    benchmark::State& state) {
  integer_vector_deserialize_owner<pixie::PackedMonotoneIntegerVector<>,
                                   Validation>(state);
}

template <pixie::DeserializationValidation Validation>
void BM_PackedMonotoneIntegerVectorDeserializeViewImpl(
    benchmark::State& state) {
  integer_vector_deserialize_view<pixie::PackedMonotoneIntegerVector<>,
                                  pixie::PackedMonotoneIntegerVectorView<>,
                                  Validation>(state);
}

template <pixie::DeserializationValidation Validation>
void BM_WaveletTreeDeserializeViewImpl(benchmark::State& state) {
  const std::size_t symbol_count = static_cast<std::size_t>(state.range(0));
  const std::vector<std::uint64_t> symbols = make_symbols(symbol_count);
  const pixie::WaveletTree<std::uint64_t> index(kWaveletAlphabetSize, symbols);
  const std::vector<std::byte> serialized = serialize_to_vector(index);
  const AlignedArtifact artifact(as_const_span(serialized));
  deserialize_iterations(
      state, artifact.bytes(), [](pixie::BinaryReader& reader) {
        return pixie::WaveletTreeView<std::uint64_t>::deserialize(reader,
                                                                  Validation);
      });
  set_artifact_counters(state, symbol_count, artifact.bytes().size());
}

#define PIXIE_DESERIALIZATION_WRAPPERS(name)                     \
  void name##Quick(benchmark::State& state) {                    \
    name##Impl<pixie::DeserializationValidation::kQuick>(state); \
  }                                                              \
  void name##Full(benchmark::State& state) {                     \
    name##Impl<pixie::DeserializationValidation::kFull>(state);  \
  }

PIXIE_DESERIALIZATION_WRAPPERS(BM_RankSelectDeserializeOwning)
PIXIE_DESERIALIZATION_WRAPPERS(BM_RankSelectDeserializeView)
PIXIE_DESERIALIZATION_WRAPPERS(BM_RmMDeserialize)
PIXIE_DESERIALIZATION_WRAPPERS(BM_RmqDeserialize)
PIXIE_DESERIALIZATION_WRAPPERS(BM_WaveletTreeDeserializeView)
PIXIE_DESERIALIZATION_WRAPPERS(BM_PackedIntegerVectorDeserializeOwning)
PIXIE_DESERIALIZATION_WRAPPERS(BM_PackedIntegerVectorDeserializeView)
PIXIE_DESERIALIZATION_WRAPPERS(BM_PackedMonotoneIntegerVectorDeserializeOwning)
PIXIE_DESERIALIZATION_WRAPPERS(BM_PackedMonotoneIntegerVectorDeserializeView)

#undef PIXIE_DESERIALIZATION_WRAPPERS

#define PIXIE_SERIALIZATION_TIMING() \
  MinWarmUpTime(kBenchmarkWarmupSeconds)->MinTime(kBenchmarkMinSeconds)

BENCHMARK(BM_BinaryWriterU64Span)
    ->ArgsProduct({{4 * kKiB, 1 * kMiB, 64 * kMiB},
                   {4 * kKiB, 64 * kKiB, 1 * kMiB}})
    ->ArgNames({"payload_bytes", "buffer_bytes"})
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryReaderU64)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterBytesSpan)
    ->Arg(1)
    ->Arg(2)
    ->Arg(64)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterBytesSplitSpanOne)
    ->Arg(1)
    ->Arg(2)
    ->Arg(64)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterBytesSplitSpanTwo)
    ->Arg(2)
    ->Arg(64)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_SpanDescriptor)->PIXIE_SERIALIZATION_TIMING();
BENCHMARK(BM_SplitSpanDescriptorOne)->PIXIE_SERIALIZATION_TIMING();
BENCHMARK(BM_SplitSpanDescriptorTwo)->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterBytesVector)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryReaderByteSpans)
    ->Args({1 * kMiB, 4 * kKiB})
    ->Args({1 * kMiB, 64 * kKiB})
    ->Args({64 * kMiB, 4 * kKiB})
    ->Args({64 * kMiB, 64 * kKiB})
    ->ArgNames({"payload_bytes", "chunk_bytes"})
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterFramedSpan)
    ->ArgsProduct({{4 * kKiB, 1 * kMiB, 64 * kMiB},
                   {4 * kKiB, 64 * kKiB, 1 * kMiB}})
    ->ArgNames({"payload_bytes", "buffer_bytes"})
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryReaderFramed)
    ->Arg(4 * kKiB)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryWriterBytesFile)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->UseRealTime()
    ->PIXIE_SERIALIZATION_TIMING();

BENCHMARK(BM_BinaryReaderU64MappedWarm)
    ->Arg(1 * kMiB)
    ->Arg(64 * kMiB)
    ->ArgName("payload_bytes")
    ->PIXIE_SERIALIZATION_TIMING();

#define PIXIE_STRUCTURE_BENCHMARK(function) \
  BENCHMARK(function)                       \
      ->Arg(1 << 10)                        \
      ->Arg(1 << 16)                        \
      ->Arg(1 << 20)                        \
      ->ArgName("N")                        \
      ->PIXIE_SERIALIZATION_TIMING()

PIXIE_STRUCTURE_BENCHMARK(BM_RankSelectSerialize);
PIXIE_STRUCTURE_BENCHMARK(BM_RankSelectDeserializeOwningQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_RankSelectDeserializeOwningFull);
PIXIE_STRUCTURE_BENCHMARK(BM_RankSelectDeserializeViewQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_RankSelectDeserializeViewFull);
PIXIE_STRUCTURE_BENCHMARK(BM_RmMSerialize);
PIXIE_STRUCTURE_BENCHMARK(BM_RmMDeserializeQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_RmMDeserializeFull);
PIXIE_STRUCTURE_BENCHMARK(BM_RmqSerialize);
PIXIE_STRUCTURE_BENCHMARK(BM_RmqDeserializeQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_RmqDeserializeFull);
PIXIE_STRUCTURE_BENCHMARK(BM_WaveletTreeSerialize);
PIXIE_STRUCTURE_BENCHMARK(BM_WaveletTreeDeserializeViewQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_WaveletTreeDeserializeViewFull);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorSerializeOwning);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorSerializeView);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorDeserializeOwningQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorDeserializeOwningFull);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorDeserializeViewQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedIntegerVectorDeserializeViewFull);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorSerializeOwning);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorSerializeView);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorDeserializeOwningQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorDeserializeOwningFull);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorDeserializeViewQuick);
PIXIE_STRUCTURE_BENCHMARK(BM_PackedMonotoneIntegerVectorDeserializeViewFull);

#undef PIXIE_STRUCTURE_BENCHMARK
#undef PIXIE_SERIALIZATION_TIMING

}  // namespace
