#pragma once

#if !defined(__unix__) && !defined(__APPLE__)
#error "pixie/io/file_output_sink.h requires POSIX file support"
#endif

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <span>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace pixie::io {

/**
 * @brief Move-only seekable sink for a newly created or truncated POSIX file.
 *
 * @details Sequential writes use `write()`, while backpatches use `pwrite()`
 * without changing the append position. `finish()` closes the descriptor and
 * reports close errors, but does not provide durable `fsync()` semantics.
 */
class FileOutputSink {
 public:
  /**
   * @brief Open @p path for binary output, creating or truncating it.
   * @throws std::system_error if the file cannot be opened.
   */
  explicit FileOutputSink(const std::filesystem::path& path) {
    descriptor_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (descriptor_ == -1) {
      throw std::system_error(errno, std::generic_category(),
                              "Failed to open binary output file");
    }
  }

  FileOutputSink(const FileOutputSink&) = delete;
  FileOutputSink& operator=(const FileOutputSink&) = delete;

  /** @brief Transfer ownership of an open output file. */
  FileOutputSink(FileOutputSink&& other) noexcept
      : descriptor_(std::exchange(other.descriptor_, -1)),
        size_bytes_(std::exchange(other.size_bytes_, 0)) {}

  /** @brief Close the current file and transfer ownership from @p other. */
  FileOutputSink& operator=(FileOutputSink&& other) noexcept {
    if (this != &other) {
      close_noexcept();
      descriptor_ = std::exchange(other.descriptor_, -1);
      size_bytes_ = std::exchange(other.size_bytes_, 0);
    }
    return *this;
  }

  /** @brief Close the output file without reporting errors. */
  ~FileOutputSink() { close_noexcept(); }

  /** @brief Return the number of bytes appended to the output. */
  std::size_t size_bytes() const noexcept { return size_bytes_; }

  /**
   * @brief Append all @p bytes to the file.
   * @throws std::system_error for an operating-system write failure.
   * @throws std::length_error if the resulting file offset is unsupported.
   * @throws std::logic_error if the sink is already finished.
   */
  void write(std::span<const std::byte> bytes) {
    require_open();
    require_file_range(size_bytes_, bytes.size());
    while (!bytes.empty()) {
      const std::size_t chunk = std::min(
          bytes.size(),
          static_cast<std::size_t>(std::numeric_limits<ssize_t>::max()));
      const ssize_t written = ::write(descriptor_, bytes.data(), chunk);
      if (written == -1) {
        if (errno == EINTR) {
          continue;
        }
        throw std::system_error(errno, std::generic_category(),
                                "Failed to write binary output file");
      }
      if (written == 0) {
        throw std::system_error(EIO, std::generic_category(),
                                "Binary output file write made no progress");
      }
      const std::size_t count = static_cast<std::size_t>(written);
      size_bytes_ += count;
      bytes = bytes.subspan(count);
    }
  }

  /**
   * @brief Replace already-written bytes beginning at @p position.
   * @throws std::out_of_range if the range is outside the written file.
   * @throws std::system_error for an operating-system write failure.
   * @throws std::length_error if the file offset is unsupported.
   * @throws std::logic_error if the sink is already finished.
   */
  void write_at(std::size_t position, std::span<const std::byte> bytes) {
    require_open();
    if (position > size_bytes_ || bytes.size() > size_bytes_ - position) {
      throw std::out_of_range("Binary output patch is outside the file");
    }
    require_file_range(position, bytes.size());
    while (!bytes.empty()) {
      const std::size_t chunk = std::min(
          bytes.size(),
          static_cast<std::size_t>(std::numeric_limits<ssize_t>::max()));
      const off_t offset = static_cast<off_t>(position);
      const ssize_t written =
          ::pwrite(descriptor_, bytes.data(), chunk, offset);
      if (written == -1) {
        if (errno == EINTR) {
          continue;
        }
        throw std::system_error(errno, std::generic_category(),
                                "Failed to patch binary output file");
      }
      if (written == 0) {
        throw std::system_error(EIO, std::generic_category(),
                                "Binary output file patch made no progress");
      }
      const std::size_t count = static_cast<std::size_t>(written);
      position += count;
      bytes = bytes.subspan(count);
    }
  }

  /**
   * @brief Close the file and report any close error.
   * @details This operation is idempotent after success.
   * @throws std::system_error if closing the file fails.
   */
  void finish() {
    if (descriptor_ == -1) {
      return;
    }
    const int descriptor = std::exchange(descriptor_, -1);
    if (::close(descriptor) == -1) {
      throw std::system_error(errno, std::generic_category(),
                              "Failed to close binary output file");
    }
  }

 private:
  static void require_file_range(std::size_t position, std::size_t count) {
    constexpr auto kMaximumOffset = std::numeric_limits<off_t>::max();
    const auto maximum = static_cast<std::uintmax_t>(kMaximumOffset);
    if (position > maximum || count > maximum - position) {
      throw std::length_error("Binary output file offset is too large");
    }
  }

  void require_open() const {
    if (descriptor_ == -1) {
      throw std::logic_error("Binary output file is already finished");
    }
  }

  void close_noexcept() noexcept {
    if (descriptor_ != -1) {
      const int descriptor = std::exchange(descriptor_, -1);
      (void)::close(descriptor);
    }
  }

  int descriptor_ = -1;
  std::size_t size_bytes_ = 0;
};

}  // namespace pixie::io
