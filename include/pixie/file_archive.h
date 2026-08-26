#pragma once

/**
 * @file file_archive.h
 * @brief Self-contained byte-oriented file archives with line extraction.
 */

#include <pixie/detail/serialization.h>
#include <pixie/storage/read_only_view.h>
#include <pixie/wavelet_tree/index.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pixie {

/** @brief Kind of filesystem entry stored in a file archive. */
enum class FileArchiveEntryType : std::uint8_t {
  kRegular = 0,
  kSymlink = 1,
};

/** @brief Source entry used to construct an owning file archive. */
struct FileArchiveSource {
  std::string path;
  std::vector<std::byte> content;
  FileArchiveEntryType type = FileArchiveEntryType::kRegular;
  std::uint32_t mode = 0;
};

/** @brief Metadata for one replayable archive-construction source. */
struct FileArchiveSourceMetadata {
  std::string path;
  FileArchiveEntryType type = FileArchiveEntryType::kRegular;
  std::uint32_t mode = 0;
};

/** @brief Public metadata for one archive entry. */
struct FileArchiveEntry {
  std::string_view path;
  std::size_t content_offset = 0;
  std::size_t content_size = 0;
  std::size_t newline_rank_base = 0;
  std::size_t line_count = 0;
  std::uint32_t mode = 0;
  FileArchiveEntryType type = FileArchiveEntryType::kRegular;
  bool is_text = false;
};

namespace file_archive_detail {

inline constexpr std::array<std::uint8_t, 8> kMagic = {'P', 'I', 'X', 'A',
                                                       'R', 'C', 'H', '1'};
inline constexpr std::uint32_t kVersion = 5;
inline constexpr std::size_t kHeaderBytes = 3 * sizeof(std::uint64_t);
inline constexpr std::size_t kMetadataHeaderBytes =
    kHeaderBytes + 2 * sizeof(std::uint32_t) + 5 * sizeof(std::uint64_t);
inline constexpr std::size_t kRecordBytes =
    6 * sizeof(std::uint64_t) + sizeof(std::uint32_t) + 4;
inline constexpr std::size_t kByteAlphabetSize = 256;

struct FileRecord {
  std::size_t path_offset = 0;
  std::size_t path_size = 0;
  std::size_t content_offset = 0;
  std::size_t content_size = 0;
  std::size_t newline_rank_base = 0;
  std::size_t line_count = 0;
  std::uint32_t mode = 0;
  FileArchiveEntryType type = FileArchiveEntryType::kRegular;
  bool is_text = false;
};

/** @brief Incrementally validate a byte stream as UTF-8. */
class Utf8Validator {
 public:
  /** @brief Construct a validator, optionally rejecting embedded NUL bytes. */
  explicit Utf8Validator(bool reject_nul) : reject_nul_(reject_nul) {}

  /** @brief Consume the next contiguous chunk of the byte stream. */
  void Consume(std::span<const std::byte> bytes) {
    if (!valid_) {
      return;
    }
    for (const std::byte raw : bytes) {
      const std::uint8_t byte = std::to_integer<std::uint8_t>(raw);
      if (remaining_ == 0) {
        if (byte == 0 && reject_nul_) {
          valid_ = false;
          return;
        }
        if (byte <= 0x7fU) {
          continue;
        }
        if (byte >= 0xc2U && byte <= 0xdfU) {
          code_point_ = byte & 0x1fU;
          minimum_ = 0x80U;
          remaining_ = 1;
        } else if (byte >= 0xe0U && byte <= 0xefU) {
          code_point_ = byte & 0x0fU;
          minimum_ = 0x800U;
          remaining_ = 2;
        } else if (byte >= 0xf0U && byte <= 0xf4U) {
          code_point_ = byte & 0x07U;
          minimum_ = 0x10000U;
          remaining_ = 3;
        } else {
          valid_ = false;
          return;
        }
        continue;
      }
      if ((byte & 0xc0U) != 0x80U) {
        valid_ = false;
        return;
      }
      code_point_ = (code_point_ << 6U) | (byte & 0x3fU);
      --remaining_;
      if (remaining_ == 0 &&
          (code_point_ < minimum_ || code_point_ > 0x10ffffU ||
           (code_point_ >= 0xd800U && code_point_ <= 0xdfffU))) {
        valid_ = false;
        return;
      }
    }
  }

  /** @brief Return whether all consumed bytes form complete valid UTF-8. */
  bool valid() const { return valid_ && remaining_ == 0; }

 private:
  bool reject_nul_;
  bool valid_ = true;
  std::uint32_t code_point_ = 0;
  std::uint32_t minimum_ = 0;
  std::uint8_t remaining_ = 0;
};

inline bool IsUtf8(std::span<const std::byte> bytes, bool reject_nul) {
  Utf8Validator validator(reject_nul);
  validator.Consume(bytes);
  return validator.valid();
}

inline bool IsUtf8(std::string_view value, bool reject_nul) {
  return IsUtf8(std::as_bytes(std::span(value.data(), value.size())),
                reject_nul);
}

inline std::size_t CheckedSize(std::uint64_t value) {
  if (value > std::numeric_limits<std::size_t>::max()) {
    throw std::length_error("File-archive size does not fit in size_t");
  }
  return static_cast<std::size_t>(value);
}

inline void RequireZeroBytes(std::span<const std::byte> bytes,
                             std::size_t offset) {
  for (const std::byte byte : bytes) {
    if (byte != std::byte{0}) {
      throw SerializationError("Non-zero file-archive padding", offset);
    }
    ++offset;
  }
}

}  // namespace file_archive_detail

/**
 * @brief CRTP facade for file-archive lookup and extraction.
 * @details `Impl` must provide stable record and path views, a queryable
 * byte wavelet tree whose symbols cover all logical content, the logical
 * content size, construction type, and path-storage size. Views returned by an
 * implementation must remain valid for the facade object's lifetime. Owning
 * implementations own this state; read-only implementations may retain
 * non-owning views whose backing storage must outlive the facade.
 * @tparam Impl Owning or read-only concrete archive implementation.
 */
template <class Impl>
class FileArchiveBase {
 public:
  /** @brief Return the serialized file-archive format version. */
  static constexpr std::uint32_t format_version() {
    return file_archive_detail::kVersion;
  }

  /** @brief Return the number of archived entries. */
  std::size_t size() const { return impl().records_impl().size(); }

  /** @brief Return whether the archive has no entries. */
  bool empty() const { return size() == 0; }

  /** @brief Return the number of logical content bytes, excluding framing. */
  std::size_t logical_size_bytes() const { return impl().logical_size_impl(); }

  /** @brief Return the wavelet-tree construction strategy. */
  WaveletTreeBuildType build_type() const { return impl().build_type_impl(); }

  /** @brief Return serialized fixed-record bytes. */
  std::size_t file_table_bytes() const {
    return size() * file_archive_detail::kRecordBytes;
  }

  /** @brief Return serialized path-blob bytes. */
  std::size_t path_storage_bytes() const {
    return impl().path_storage_size_impl();
  }

  /** @brief Return framing, file-table, path, and alignment bytes. */
  std::size_t metadata_bytes() const {
    const std::size_t unaligned = file_archive_detail::kMetadataHeaderBytes +
                                  file_table_bytes() + path_storage_bytes();
    return unaligned +
           (alignof(std::uint64_t) - unaligned % alignof(std::uint64_t)) %
               alignof(std::uint64_t);
  }

  /**
   * @brief Find an exact archive-relative path.
   * @return Its zero-based entry index, or `std::nullopt` when absent.
   */
  std::optional<std::size_t> find(std::string_view path) const {
    const auto records = impl().records_impl();
    const auto position =
        std::lower_bound(records.begin(), records.end(), path,
                         [this](const file_archive_detail::FileRecord& record,
                                std::string_view wanted) {
                           return impl().path_impl(record) < wanted;
                         });
    if (position == records.end() || impl().path_impl(*position) != path) {
      return std::nullopt;
    }
    return static_cast<std::size_t>(position - records.begin());
  }

  /**
   * @brief Return metadata for an entry.
   * @throws std::out_of_range if @p index is not an archived entry.
   */
  FileArchiveEntry entry(std::size_t index) const {
    const auto records = impl().records_impl();
    if (index >= records.size()) {
      throw std::out_of_range("File-archive entry index is out of range");
    }
    const auto& record = records[index];
    return {.path = impl().path_impl(record),
            .content_offset = record.content_offset,
            .content_size = record.content_size,
            .newline_rank_base = record.newline_rank_base,
            .line_count = record.line_count,
            .mode = record.mode,
            .type = record.type,
            .is_text = record.is_text};
  }

  /**
   * @brief Reconstruct the complete byte content of one entry.
   * @throws std::out_of_range if @p index is not an archived entry.
   */
  std::vector<std::byte> extract(std::size_t index) const {
    const FileArchiveEntry metadata = entry(index);
    return extract_range(metadata.content_offset,
                         metadata.content_offset + metadata.content_size);
  }

  /**
   * @brief Reconstruct zero-based half-open line range `[left, right)`.
   * @details LF terminators are retained. Empty intervals are allowed. A
   * trailing LF does not create an additional line.
   * @throws std::out_of_range if @p index or `[left, right)` is out of range.
   * @throws std::invalid_argument if the entry is not a regular text file.
   */
  std::vector<std::byte> extract_lines(std::size_t index,
                                       std::size_t left,
                                       std::size_t right) const {
    const FileArchiveEntry metadata = entry(index);
    if (metadata.type != FileArchiveEntryType::kRegular) {
      throw std::invalid_argument("Line extraction requires a regular file");
    }
    if (!metadata.is_text) {
      throw std::invalid_argument("Line extraction requires a text file");
    }
    if (left > right || right > metadata.line_count) {
      throw std::out_of_range("File-archive line range is out of range");
    }
    const auto& tree = impl().tree_impl();
    const auto boundary = [&](std::size_t line) {
      if (line == 0) {
        return metadata.content_offset;
      }
      if (line == metadata.line_count) {
        return metadata.content_offset + metadata.content_size;
      }
      return tree.select('\n', metadata.newline_rank_base + line) + 1;
    };
    const std::size_t begin = boundary(left);
    const std::size_t end = left == right ? begin : boundary(right);
    return extract_range(begin, end);
  }

 private:
  const Impl& impl() const { return static_cast<const Impl&>(*this); }

  std::vector<std::byte> extract_range(std::size_t begin,
                                       std::size_t end) const {
    const std::vector<std::uint8_t> symbols =
        impl().tree_impl().get_segment(begin, end);
    std::vector<std::byte> bytes;
    bytes.reserve(symbols.size());
    for (const std::uint8_t symbol : symbols) {
      bytes.push_back(static_cast<std::byte>(symbol));
    }
    return bytes;
  }
};

/**
 * @brief Storage-backed byte-oriented file archive.
 * @tparam Storage Owning aligned storage or a non-owning read-only view.
 */
template <StorageImplementation Storage>
class FileArchiveIndex : public FileArchiveBase<FileArchiveIndex<Storage>>,
                         public SerializationBase<FileArchiveIndex<Storage>> {
 public:
  FileArchiveIndex()
    requires(std::same_as<Storage, ReadOnlyStorageView>)
  = default;

  /**
   * @brief Owning archives are move-only; archive views are copyable.
   * @details Internal rank indexes retain views into their wavelet-node bit
   * storage, so copying without rebuilding those indexes would leave views
   * bound to the source archive.
   */
  FileArchiveIndex(const FileArchiveIndex&)
    requires(std::same_as<Storage, ReadOnlyStorageView>)
  = default;
  FileArchiveIndex& operator=(const FileArchiveIndex&)
    requires(std::same_as<Storage, ReadOnlyStorageView>)
  = default;
  FileArchiveIndex(const FileArchiveIndex&)
    requires(std::same_as<Storage, AlignedStorage>)
  = delete;
  FileArchiveIndex& operator=(const FileArchiveIndex&)
    requires(std::same_as<Storage, AlignedStorage>)
  = delete;
  FileArchiveIndex(FileArchiveIndex&&) noexcept = default;
  FileArchiveIndex& operator=(FileArchiveIndex&&) noexcept = default;

  /**
   * @brief Construct an archive from file and symlink sources.
   * @details Sources are sorted by path. Paths must be non-empty, unique,
   * valid UTF-8 strings. Content is preserved byte-for-byte.
   */
  explicit FileArchiveIndex(
      std::vector<FileArchiveSource> sources,
      WaveletTreeBuildType build_type = WaveletTreeBuildType::Huffman)
    requires(std::same_as<Storage, AlignedStorage>)
      : build_type_(build_type) {
    build_sources(sources, [](const FileArchiveSource& source, auto&& consume) {
      consume(std::span<const std::byte>(source.content));
    });
  }

  /**
   * @brief Construct from metadata and a replayable chunk reader.
   * @details The reader is called twice per source in sorted path order. It
   * receives `(const FileArchiveSourceMetadata&, consumer)` and must pass the
   * same immutable content to `consumer` as byte spans on both calls. The
   * first pass derives metadata and symbol counts; the second constructs the
   * tree without retaining complete source contents.
   * @throws std::invalid_argument for invalid metadata or changed content.
   */
  template <class ReadSource>
  FileArchiveIndex(
      std::vector<FileArchiveSourceMetadata> sources,
      ReadSource&& read_source,
      WaveletTreeBuildType build_type = WaveletTreeBuildType::Huffman)
    requires(std::same_as<Storage, AlignedStorage>)
      : build_type_(build_type) {
    build_sources(sources, read_source);
  }

  /** @brief Serialize one native, framed Pixie file archive. */
  void serialize_impl(BinaryWriter& writer) const {
    if (writer.size_bytes() % alignof(std::uint64_t) != 0) {
      throw std::invalid_argument(
          "File-archive serialization requires an aligned writer offset");
    }
    const std::size_t begin = writer.size_bytes();
    detail::write_magic(writer, file_archive_detail::kMagic);
    writer.write_u32(file_archive_detail::kVersion);
    writer.write_u32(0);
    const std::size_t size_position = writer.write_u64_placeholder();
    writer.write_u32(static_cast<std::uint32_t>(build_type_));
    writer.write_u32(0);
    writer.write_size(records_.size());
    writer.write_size(logical_size_);
    writer.write_size(records_.size() * file_archive_detail::kRecordBytes);
    writer.write_size(paths_.size_bytes());
    writer.write_size(0);
    for (const auto& record : records_) {
      writer.write_size(record.path_offset);
      writer.write_size(record.path_size);
      writer.write_size(record.content_offset);
      writer.write_size(record.content_size);
      writer.write_size(record.newline_rank_base);
      writer.write_size(record.line_count);
      writer.write_u32(record.mode);
      writer.write_u8(static_cast<std::uint8_t>(record.type));
      writer.write_u8(static_cast<std::uint8_t>(record.is_text));
      writer.write_u16(0);
    }
    writer.write_bytes(paths_.as_bytes());
    writer.align_to(alignof(std::uint64_t));
    if (!tree_.has_value()) {
      throw std::logic_error("Cannot serialize an uninitialized file archive");
    }
    tree_->serialize(writer);
    writer.align_to(alignof(std::uint64_t));
    writer.patch_u64(size_position,
                     static_cast<std::uint64_t>(writer.size_bytes() - begin));
  }

  /**
   * @brief Deserialize one framed archive using this specialization's storage.
   * @details Owners copy paths and all wavelet-tree storage. Views retain
   * slices of the reader's backing bytes, which must outlive the result.
   */
  static FileArchiveIndex deserialize_impl(
      BinaryReader& reader,
      DeserializationValidation validation = DeserializationValidation::kQuick)
    requires(std::same_as<Storage, AlignedStorage> ||
             std::same_as<Storage, ReadOnlyStorageView>)
  {
    if constexpr (std::same_as<Storage, ReadOnlyStorageView>) {
      if (reinterpret_cast<std::uintptr_t>(reader.remaining_bytes().data()) %
              alignof(std::uint64_t) !=
          0) {
        throw std::invalid_argument(
            "Serialized file archive is not word aligned");
      }
    }
    const std::size_t available = reader.remaining();
    detail::require_magic(reader, file_archive_detail::kMagic);
    if (reader.read_u32() != file_archive_detail::kVersion ||
        reader.read_u32() != 0) {
      throw std::invalid_argument("Incompatible serialized file archive");
    }
    const std::size_t artifact_size = detail::checked_artifact_size(
        reader.read_u64(), file_archive_detail::kHeaderBytes, available);
    BinaryReader payload = reader.read_subreader(
        artifact_size - file_archive_detail::kHeaderBytes);

    FileArchiveIndex result(LoadTag{});
    const std::uint32_t build_type = payload.read_u32();
    if (build_type >
        static_cast<std::uint32_t>(WaveletTreeBuildType::Huffman)) {
      throw std::invalid_argument("Invalid file-archive build type");
    }
    result.build_type_ = static_cast<WaveletTreeBuildType>(build_type);
    if (payload.read_u32() != 0) {
      throw std::invalid_argument("Invalid file-archive reserved field");
    }
    const std::size_t count = payload.read_size();
    result.logical_size_ = payload.read_size();
    const std::size_t records_bytes = payload.read_size();
    const std::size_t paths_bytes = payload.read_size();
    if (payload.read_size() != 0 ||
        count > std::numeric_limits<std::size_t>::max() /
                    file_archive_detail::kRecordBytes ||
        records_bytes != count * file_archive_detail::kRecordBytes) {
      throw std::invalid_argument("Invalid file-archive section sizes");
    }

    BinaryReader records_reader = payload.read_subreader(records_bytes);
    result.records_.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
      file_archive_detail::FileRecord record;
      record.path_offset = records_reader.read_size();
      record.path_size = records_reader.read_size();
      record.content_offset = records_reader.read_size();
      record.content_size = records_reader.read_size();
      record.newline_rank_base = records_reader.read_size();
      record.line_count = records_reader.read_size();
      record.mode = records_reader.read_u32();
      const std::uint8_t type = records_reader.read_u8();
      const std::uint8_t is_text = records_reader.read_u8();
      if (type > static_cast<std::uint8_t>(FileArchiveEntryType::kSymlink) ||
          is_text > 1 || records_reader.read_u16() != 0) {
        throw std::invalid_argument("Invalid file-archive entry flags");
      }
      record.type = static_cast<FileArchiveEntryType>(type);
      record.is_text = is_text != 0;
      result.records_.push_back(record);
    }
    if (!records_reader.empty()) {
      throw std::invalid_argument("Invalid file-archive record section");
    }
    result.paths_ = MakeStorage(payload.read_bytes(paths_bytes));
    const std::size_t padding =
        (alignof(std::uint64_t) -
         payload.byte_offset() % alignof(std::uint64_t)) %
        alignof(std::uint64_t);
    const std::size_t padding_offset = payload.byte_offset();
    file_archive_detail::RequireZeroBytes(payload.read_bytes(padding),
                                          padding_offset);
    result.tree_.emplace(WaveletTreeIndex<std::uint8_t, Storage>::deserialize(
        payload, validation));
    payload.require_zero_padding(alignof(std::uint64_t) - 1);
    result.validate(validation == DeserializationValidation::kFull);
    return result;
  }

 private:
  struct LoadTag {};

  explicit FileArchiveIndex(LoadTag) {}

  static Storage MakeStorage(std::span<const std::byte> bytes) {
    if constexpr (std::same_as<Storage, ReadOnlyStorageView>) {
      return ReadOnlyStorageView(bytes);
    } else {
      if (bytes.size() > std::numeric_limits<std::size_t>::max() / 8) {
        throw std::length_error("File-archive path storage is too large");
      }
      AlignedStorage result(bytes.size() * 8);
      std::ranges::copy(bytes, result.writable_bytes().begin());
      return result;
    }
  }

  template <class Source, class ReadSource>
  void build_sources(std::vector<Source>& sources, ReadSource&& read_source) {
    std::sort(sources.begin(), sources.end(),
              [](const auto& left, const auto& right) {
                return left.path < right.path;
              });
    std::array<std::size_t, file_archive_detail::kByteAlphabetSize>
        symbol_counts{};
    std::vector<std::uint64_t> content_hashes;
    content_hashes.reserve(sources.size());
    std::string paths;
    std::size_t newline_rank = 0;
    for (std::size_t source_index = 0; source_index < sources.size();
         ++source_index) {
      const Source& source = sources[source_index];
      if (source.path.empty() ||
          !file_archive_detail::IsUtf8(source.path, true)) {
        throw std::invalid_argument(
            "File-archive paths must be non-empty UTF-8 strings");
      }
      if (source_index != 0 && sources[source_index - 1].path == source.path) {
        throw std::invalid_argument("File-archive paths must be unique");
      }
      if (source.type != FileArchiveEntryType::kRegular &&
          source.type != FileArchiveEntryType::kSymlink) {
        throw std::invalid_argument("Invalid file-archive entry type");
      }

      file_archive_detail::FileRecord record;
      record.path_offset = paths.size();
      record.path_size = source.path.size();
      record.content_offset = logical_size_;
      record.newline_rank_base = newline_rank;
      record.mode = source.mode;
      record.type = source.type;
      paths.append(source.path);

      file_archive_detail::Utf8Validator utf8(/*reject_nul=*/true);
      std::size_t newlines = 0;
      bool has_content = false;
      std::uint8_t last_byte = 0;
      std::uint64_t content_hash = 14695981039346656037ULL;
      read_source(source, [&](std::span<const std::byte> chunk) {
        if (chunk.size() >
            std::numeric_limits<std::size_t>::max() - logical_size_) {
          throw std::length_error("File-archive content is too large");
        }
        logical_size_ += chunk.size();
        utf8.Consume(chunk);
        for (const std::byte byte : chunk) {
          const std::uint8_t value = std::to_integer<std::uint8_t>(byte);
          ++symbol_counts[value];
          newlines += value == '\n' ? 1U : 0U;
          has_content = true;
          last_byte = value;
          content_hash ^= value;
          content_hash *= 1099511628211ULL;
        }
      });
      record.content_size = logical_size_ - record.content_offset;
      record.is_text = utf8.valid();
      if (source.type == FileArchiveEntryType::kRegular && has_content) {
        record.line_count =
            newlines + static_cast<std::size_t>(last_byte != '\n');
      }
      newline_rank += newlines;
      records_.push_back(record);
      content_hashes.push_back(content_hash);
    }

    paths_ = MakeStorage(
        std::as_bytes(std::span<const char>(paths.data(), paths.size())));

    tree_.emplace(
        file_archive_detail::kByteAlphabetSize, symbol_counts,
        [&](auto&& emit) {
          for (std::size_t index = 0; index < sources.size(); ++index) {
            std::size_t content_size = 0;
            std::uint64_t content_hash = 14695981039346656037ULL;
            read_source(sources[index], [&](std::span<const std::byte> chunk) {
              if (content_size > records_[index].content_size ||
                  chunk.size() > records_[index].content_size - content_size) {
                throw std::invalid_argument(
                    "File-archive source changed between build passes");
              }
              content_size += chunk.size();
              for (const std::byte byte : chunk) {
                const std::uint8_t value = std::to_integer<std::uint8_t>(byte);
                content_hash ^= value;
                content_hash *= 1099511628211ULL;
                emit(value);
              }
            });
            if (content_size != records_[index].content_size ||
                content_hash != content_hashes[index]) {
              throw std::invalid_argument(
                  "File-archive source changed between build passes");
            }
          }
        },
        build_type_);
  }

  /** @brief Permit the facade to use the documented extension points. */
  friend class FileArchiveBase<FileArchiveIndex>;

  void validate(bool full) const {
    if (tree_->size() != logical_size_) {
      throw std::invalid_argument("Invalid file-archive content size");
    }
    const std::size_t newline_count = tree_->rank('\n', logical_size_);
    std::size_t expected_content_offset = 0;
    std::size_t expected_newline_rank = 0;
    std::string_view previous_path;
    for (std::size_t index = 0; index < records_.size(); ++index) {
      const auto& record = records_[index];
      if (record.path_offset > paths_.size_bytes() ||
          record.path_size > paths_.size_bytes() - record.path_offset ||
          record.content_offset > logical_size_ ||
          record.content_size > logical_size_ - record.content_offset ||
          record.newline_rank_base > newline_count ||
          record.line_count > newline_count - record.newline_rank_base + 1 ||
          (record.type == FileArchiveEntryType::kSymlink &&
           record.line_count != 0)) {
        throw std::invalid_argument("Invalid file-archive entry bounds");
      }
      if (!full) {
        continue;
      }
      const std::string_view path = path_impl(record);
      if (path.empty() || !file_archive_detail::IsUtf8(path, true) ||
          (index != 0 && previous_path >= path) ||
          record.content_offset != expected_content_offset ||
          record.newline_rank_base != expected_newline_rank) {
        throw std::invalid_argument("Invalid file-archive entry metadata");
      }
      const std::vector<std::byte> content = this->extract(index);
      const bool is_text = file_archive_detail::IsUtf8(content, true);
      std::size_t newlines = 0;
      for (const std::byte byte : content) {
        newlines += byte == std::byte{'\n'} ? 1U : 0U;
      }
      std::size_t line_count = 0;
      if (record.type == FileArchiveEntryType::kRegular && !content.empty()) {
        line_count = newlines + static_cast<std::size_t>(content.back() !=
                                                         std::byte{'\n'});
      }
      if (record.is_text != is_text || record.line_count != line_count) {
        throw std::invalid_argument("Invalid file-archive derived metadata");
      }
      previous_path = path;
      expected_content_offset += content.size();
      expected_newline_rank += newlines;
    }
    if (full && expected_content_offset != logical_size_) {
      throw std::invalid_argument("Invalid file-archive content layout");
    }
  }

  /** @brief Return a stable view of validated, sorted records. */
  std::span<const file_archive_detail::FileRecord> records_impl() const {
    return records_;
  }
  /** @brief Return the path slice identified by a validated record. */
  std::string_view path_impl(
      const file_archive_detail::FileRecord& record) const {
    const std::span<const std::byte> bytes = paths_.as_bytes();
    const std::string_view paths =
        bytes.empty()
            ? std::string_view{}
            : std::string_view(reinterpret_cast<const char*>(bytes.data()),
                               bytes.size());
    return paths.substr(record.path_offset, record.path_size);
  }
  /** @brief Return the storage-backed byte tree containing archive content. */
  const WaveletTreeIndex<std::uint8_t, Storage>& tree_impl() const {
    return *tree_;
  }
  /** @brief Return the number of archived content bytes. */
  std::size_t logical_size_impl() const { return logical_size_; }
  /** @brief Return the serialized tree construction strategy. */
  WaveletTreeBuildType build_type_impl() const { return build_type_; }
  /** @brief Return bytes occupied by the retained path blob. */
  std::size_t path_storage_size_impl() const { return paths_.size_bytes(); }

  std::vector<file_archive_detail::FileRecord> records_;
  Storage paths_;
  std::optional<WaveletTreeIndex<std::uint8_t, Storage>> tree_;
  std::size_t logical_size_ = 0;
  WaveletTreeBuildType build_type_ = WaveletTreeBuildType::Huffman;
};

/** @brief Owning aligned-storage file archive. */
using FileArchive = FileArchiveIndex<AlignedStorage>;

/** @brief Read-only file archive retaining serialized storage views. */
using FileArchiveView = FileArchiveIndex<ReadOnlyStorageView>;

}  // namespace pixie
