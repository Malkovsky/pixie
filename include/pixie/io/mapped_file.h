#pragma once

#if !defined(__unix__) && !defined(__APPLE__)
#error "pixie/io/mapped_file.h requires POSIX mmap support"
#endif

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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
 * @brief Move-only owner of a read-only POSIX memory-mapped file.
 * @details The mapping remains valid after the file is unlinked. Empty files
 * expose an empty byte view without creating a zero-length mapping.
 */
class MappedFile {
 public:
  /**
   * @brief Map @p path for reading.
   * @throws std::system_error if the file cannot be opened, inspected, or
   * mapped.
   * @throws std::length_error if the file does not fit in `std::size_t`.
   */
  explicit MappedFile(const std::filesystem::path& path) {
    const int descriptor = ::open(path.c_str(), O_RDONLY);
    if (descriptor == -1) {
      throw std::system_error(errno, std::generic_category(),
                              "Failed to open mapped file");
    }

    struct stat status {};
    if (::fstat(descriptor, &status) == -1) {
      const int error = errno;
      close_descriptor(descriptor);
      throw std::system_error(error, std::generic_category(),
                              "Failed to inspect mapped file");
    }

    const auto file_size = static_cast<std::uintmax_t>(status.st_size);
    if (status.st_size < 0 ||
        file_size > std::numeric_limits<std::size_t>::max()) {
      close_descriptor(descriptor);
      throw std::length_error("Mapped file is too large");
    }

    const auto size = static_cast<std::size_t>(file_size);
    if (size == 0) {
      close_descriptor(descriptor);
      return;
    }

    void* const mapping =
        ::mmap(nullptr, size, PROT_READ, MAP_SHARED, descriptor, 0);
    const int error = mapping == MAP_FAILED ? errno : 0;
    close_descriptor(descriptor);
    if (mapping == MAP_FAILED) {
      throw std::system_error(error, std::generic_category(),
                              "Failed to map file");
    }

    data_ = mapping;
    size_bytes_ = size;
  }

  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  /** @brief Transfer ownership of a mapping. */
  MappedFile(MappedFile&& other) noexcept
      : data_(std::exchange(other.data_, nullptr)),
        size_bytes_(std::exchange(other.size_bytes_, 0)) {}

  /** @brief Release this mapping and transfer ownership from @p other. */
  MappedFile& operator=(MappedFile&& other) noexcept {
    if (this != &other) {
      reset();
      data_ = std::exchange(other.data_, nullptr);
      size_bytes_ = std::exchange(other.size_bytes_, 0);
    }
    return *this;
  }

  /** @brief Unmap the file contents, if any. */
  ~MappedFile() { reset(); }

  /** @brief Return the mapped file contents. */
  std::span<const std::byte> as_bytes() const noexcept {
    if (data_ == nullptr) {
      return {};
    }
    return {static_cast<const std::byte*>(data_), size_bytes_};
  }

  /** @brief Return the mapped file size in bytes. */
  std::size_t size_bytes() const noexcept { return size_bytes_; }

 private:
  static void close_descriptor(int descriptor) noexcept {
    if (descriptor != -1) {
      ::close(descriptor);
    }
  }

  void reset() noexcept {
    if (data_ != nullptr) {
      ::munmap(data_, size_bytes_);
      data_ = nullptr;
      size_bytes_ = 0;
    }
  }

  void* data_ = nullptr;
  std::size_t size_bytes_ = 0;
};

}  // namespace pixie::io
