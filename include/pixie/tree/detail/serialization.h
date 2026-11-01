#pragma once

#include <pixie/detail/serialization_frame.h>
#include <pixie/tree/detail/encoding_storage.h>

namespace pixie::detail {

enum class TreeEncodingKind : uint64_t { kBP = 1, kDFUDS = 2, kLOUDS = 3 };

inline uint64_t tree_serialization_tag(TreeEncodingKind kind) {
  return 0x3152540058495050ULL + (static_cast<uint64_t>(kind) << 32);
}

inline void write_tree_encoding(BinaryWriter& writer,
                                std::span<const uint64_t> words,
                                size_t bits) {
  writer.write_size(words.size_bytes());
  for (size_t i = 0; i < words.size(); ++i) {
    auto value = words[i];
    if (i + 1 == words.size() && bits % 64) {
      value &= (uint64_t{1} << (bits % 64)) - 1;
    }
    writer.write_u64(value);
  }
}

template <TreeEncodingStorage Storage>
Storage read_tree_encoding(BinaryReader& reader, size_t bits) {
  const auto bytes = reader.read_size();
  const auto count = bits / 64 + (bits % 64 != 0);
  if (bytes != count * 8) {
    throw std::invalid_argument("Invalid serialized tree encoding size");
  }
  const auto payload = reader.read_bytes(bytes);
  Storage storage(payload);
  if constexpr (std::endian::native != std::endian::little) {
    if constexpr (std::same_as<Storage, ReadOnlyStorageView>) {
      throw std::invalid_argument(
          "Mapped tree encodings require a little-endian host");
    } else {
      BinaryReader words(payload);
      for (auto& word : storage.as_words64()) {
        word = words.read_u64();
      }
    }
  }
  const auto words = tree_encoding_words(storage, bits);
  if (bits % 64 && (words.back() >> (bits % 64))) {
    throw std::invalid_argument("Nonzero tree encoding padding");
  }
  return storage;
}

inline void validate_tree_topology(std::span<const uint64_t> words,
                                   size_t nodes,
                                   TreeEncodingKind kind) {
  const auto bits = tree_encoding_bits(nodes, kind == TreeEncodingKind::kBP);
  int64_t pending = kind == TreeEncodingKind::kBP ? 0 : 1;
  size_t ones = 0;
  for (size_t i = 0; i < bits; ++i) {
    const bool bit = (words[i / 64] >> (i % 64)) & 1;
    ones += bit;
    if (kind == TreeEncodingKind::kLOUDS) {
      if ((!i && !bit) || (bit && pending == 0)) {
        throw std::invalid_argument("Invalid LOUDS topology");
      }
      pending += bit ? -1 : 1;
    } else {
      pending += bit ? 1 : -1;
      if (pending < 0 || (!pending && i + 1 < bits)) {
        throw std::invalid_argument("Invalid parenthesis topology");
      }
    }
  }
  const auto expected_ones =
      kind == TreeEncodingKind::kDFUDS && nodes ? nodes - 1 : nodes;
  if ((nodes && pending != 0) || ones != expected_ones) {
    throw std::invalid_argument("Invalid tree topology counts");
  }
}

}  // namespace pixie::detail
