#pragma once

#include <pixie/storage/aligned.h>
#include <pixie/storage/read_only_view.h>

#include <limits>

namespace pixie::detail {
enum class TreeEncodingKind : uint64_t { kBP = 1, kDFUDS = 2, kLOUDS = 3 };

template <class Storage>
concept TreeEncodingStorage = std::same_as<Storage, AlignedStorage> ||
                              std::same_as<Storage, ReadOnlyStorageView>;

inline size_t tree_encoding_bits(size_t nodes, bool bp) {
  if (nodes > std::numeric_limits<size_t>::max() / 2) {
    throw std::length_error("Tree encoding bit count is too large");
  }
  return nodes ? 2 * nodes - (bp ? 0 : 1) : 0;
}

template <TreeEncodingStorage Storage>
std::span<const uint64_t> tree_encoding_words(const Storage& storage,
                                              size_t bits) {
  const size_t bytes = bits / 8 + (bits % 8 != 0);
  if (storage.size_bytes() < bytes) {
    throw std::invalid_argument("Tree encoding is too short");
  }
  const auto words = [&] {
    if constexpr (std::same_as<Storage, AlignedStorage>) {
      return storage.padded_view().as_words64();
    } else {
      return storage.as_words64();
    }
  }();
  const size_t count = bits / 64 + (bits % 64 != 0);
  if (words.size() < count) {
    throw std::invalid_argument("Tree encoding needs complete words");
  }
  return words.first(count);
}

template <TreeEncodingStorage Storage>
Storage tree_encoding_from_words(std::span<const uint64_t> words, size_t bits) {
  const size_t count = bits / 64 + (bits % 64 != 0);
  if (words.size() < count) {
    throw std::invalid_argument("Tree encoding is too short");
  }
  return Storage(std::as_bytes(words.first(count)));
}

template <TreeEncodingStorage Storage>
size_t tree_encoding_owned_bytes(const Storage& storage) {
  if constexpr (std::same_as<Storage, AlignedStorage>) {
    return storage.allocated_bytes();
  } else {
    return 0;
  }
}
}  // namespace pixie::detail
