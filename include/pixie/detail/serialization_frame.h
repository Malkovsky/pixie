#pragma once

#include <pixie/detail/serialization.h>

namespace pixie::detail {

// Each independently readable section has a type tag, version and bounded
// payload. Lengths are relative to the section, so nesting needs no relocation.
inline size_t begin_frame(BinaryWriter& writer, uint64_t tag) {
  const auto begin = writer.size_bytes();
  writer.write_u64(tag);
  writer.write_u64(1);
  writer.write_u64_placeholder();
  return begin;
}

inline void end_frame(BinaryWriter& writer, size_t begin) {
  const auto padding = (8 - (writer.size_bytes() - begin) % 8) % 8;
  writer.write_zeros(padding);
  writer.patch_u64(begin + 16, writer.size_bytes() - begin);
}

inline BinaryReader read_frame(BinaryReader& reader, uint64_t tag) {
  const auto available = reader.remaining();
  if (reader.read_u64() != tag || reader.read_u64() != 1) {
    throw std::invalid_argument("Incompatible serialized section");
  }
  const auto bytes = checked_artifact_size(reader.read_u64(), 24, available);
  return reader.read_subreader(bytes - 24);
}

}  // namespace pixie::detail
