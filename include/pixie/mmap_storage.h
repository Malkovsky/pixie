#pragma once

#include <cstring>
#include <filesystem>
#include <mio/mmap.hpp>
#include <span>
#include <stdexcept>

#include "bit_stream.h"

namespace pixie {

class MmapFile {
 private:
  mio::mmap_source mmap_;

 public:
  explicit MmapFile(const std::filesystem::path& path) : mmap_(path.string()) {
    if (!mmap_.is_open()) {
      throw std::runtime_error("Failed to mmap file");
    }
  }

  MmapFile(const MmapFile&) = delete;
  MmapFile& operator=(const MmapFile&) = delete;

  MmapFile(MmapFile&&) noexcept = default;
  MmapFile& operator=(MmapFile&&) noexcept = default;

  std::span<const std::byte> AsBytes() const {
    return {reinterpret_cast<const std::byte*>(mmap_.data()), mmap_.size()};
  }

  size_t Size() const noexcept { return mmap_.size(); }
};

class MmapViewStorage {
 private:
  std::span<const std::byte> data_;

 public:
  MmapViewStorage() = default;
  explicit MmapViewStorage(std::span<const std::byte> data) : data_(data) {}

  void resize(size_t) {
    throw std::logic_error("Cannot resize read-only mmap storage");
  }

  std::span<const uint64_t> AsConst64BitInts() const {
    return {reinterpret_cast<const uint64_t*>(data_.data()), data_.size() / 8};
  }

  std::span<const std::uint16_t> AsConst16BitInts() const {
    return {reinterpret_cast<const uint16_t*>(data_.data()), data_.size() / 2};
  }

  std::span<const std::byte> AsConstBytes() const { return data_; }

  std::span<uint64_t> As64BitInts() {
    throw std::logic_error("Read-only storage");
  }
  std::span<std::uint16_t> As16BitInts() {
    throw std::logic_error("Read-only storage");
  }

  /** @brief Writes bit representation of data to OutputBitStream */
  void serialize(pixie::OutputBitStream& bs) const {
    bs << data_.size();
    for (const std::byte byte : data_) {
      bs << static_cast<uint8_t>(byte);
    }
  }

  static MmapViewStorage deserialize(std::span<const std::byte>& data) {
    MmapViewStorage result;
    constexpr size_t length = sizeof(size_t);
    size_t size;
    std::memcpy(&size, data.data(), length);
    result.data_ = data.subspan(length, size);
    data = data.subspan(length + size);
    return result;
  };
};

}  // namespace pixie
