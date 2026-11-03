#pragma once

#include <pixie/detail/serialization_frame.h>
#include <pixie/tree/detail/encoding_storage.h>

namespace pixie::detail {

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

template <TreeEncodingKind Kind, class Backend>
void serialize_tree(BinaryWriter& writer,
                    size_t nodes,
                    std::span<const uint64_t> words,
                    const Backend& backend) {
  const auto begin = begin_frame(writer, tree_serialization_tag(Kind));
  writer.write_size(nodes);
  write_tree_encoding(writer, words,
                      tree_encoding_bits(nodes, Kind == TreeEncodingKind::kBP));
  backend.serialize(writer);
  end_frame(writer, begin);
}

// Construct the final tree before committing the reader. The factory transfers
// the encoding allocation and its borrowing backend together, without exposing
// an intermediate state whose copy could leave the backend bound to old words.
template <TreeEncodingKind Kind,
          TreeEncodingStorage Storage,
          class Backend,
          class Construct>
auto deserialize_tree(BinaryReader& reader,
                      DeserializationValidation validation,
                      Construct construct) {
  auto candidate = reader;
  auto payload = read_frame(candidate, tree_serialization_tag(Kind));
  const auto nodes = payload.read_size();
  const auto bits = tree_encoding_bits(nodes, Kind == TreeEncodingKind::kBP);
  auto encoding = read_tree_encoding<Storage>(payload, bits);
  const auto words = tree_encoding_words(encoding, bits);
  auto backend = Backend::deserialize(payload, words, validation);
  if (backend.size() != bits) {
    throw std::invalid_argument("Tree backend size mismatch");
  }
  if constexpr (requires {
                  backend.supports_select1();
                  backend.supports_select0();
                }) {
    if (!backend.supports_select1() || !backend.supports_select0()) {
      throw std::invalid_argument(
          "Tree backend requires both select directions");
    }
  }
  const auto ones = [&] {
    if constexpr (Kind == TreeEncodingKind::kLOUDS) {
      return backend.rank(bits);
    } else {
      return backend.rank1(bits);
    }
  }();
  const auto expected_ones =
      Kind == TreeEncodingKind::kDFUDS && nodes ? nodes - 1 : nodes;
  if (ones != expected_ones) {
    throw std::invalid_argument("Tree backend population mismatch");
  }
  if (validation == DeserializationValidation::kFull) {
    validate_tree_topology(words, nodes, Kind);
  }
  payload.require_zero_padding(7);
  auto result = construct(std::move(encoding), bits, std::move(backend));
  reader = candidate;
  return result;
}

}  // namespace pixie::detail
