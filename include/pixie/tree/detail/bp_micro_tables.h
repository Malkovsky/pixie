#pragma once

#include <pixie/tree.h>
#include <pixie/tree/detail/unpack_cache.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

namespace pixie::detail::bp_micro {
// Cache decoding by topology, not by node address. Trailing bits beyond the
// first complete subtree are irrelevant. The table holds only relative handles.
struct BpMicroProgram {
  std::array<uint8_t, 24> events{};
  uint32_t returns = 0;
  uint8_t bits = 0, nodes = 0, count = 0;
};
inline constexpr auto bp_micro_programs = [] {
  std::array<BpMicroProgram, 8192> table{};
  for (unsigned key = 0; key < table.size(); ++key) {
    const unsigned bits = (key << 1) | 1;
    int balance = 0;
    unsigned length = 0;
    for (unsigned b = 0; b < 14; ++b) {
      balance += ((bits >> b) & 1) ? 1 : -1;
      if (!balance) {
        length = b + 1;
        break;
      }
    }
    if (!length) {
      continue;
    }
    auto& p = table[key];
    p.bits = length;
    std::array<uint8_t, 7> pos{}, parent{}, stack{}, degree{}, next{};
    next.fill(7);
    unsigned depth = 0;
    for (unsigned b = 0; b < length; ++b) {
      if ((bits >> b) & 1) {
        const auto i = p.nodes++;
        pos[i] = b;
        if (depth) {
          parent[i] = stack[depth - 1];
          ++degree[parent[i]];
        }
        stack[depth++] = i;
      } else {
        --depth;
      }
    }
    for (unsigned i = 1; i < p.nodes; ++i) {
      for (unsigned j = i + 1; j < p.nodes; ++j) {
        if (parent[i] == parent[j]) {
          next[i] = j;
          break;
        }
      }
    }
    auto event = [&](unsigned i, bool inspect, bool returning) {
      if (returning) {
        p.returns |= uint32_t{1} << p.count;
      }
      p.events[p.count++] = i | (pos[i] << 3) | (inspect ? 128 : 0);
    };
    for (unsigned i = 0; i < p.nodes; ++i) {
      event(i, false, false);
      if (degree[i]) {
        for (unsigned c = i + 1; c != 7; c = next[c]) {
          event(c, true, false);
        }
      } else {
        auto completed = i;
        while (completed && next[completed] == 7) {
          completed = parent[completed];
          event(completed, false, true);
        }
      }
    }
  }
  return table;
}();

struct BpMicroTables {
  std::array<uint16_t, 8192> headers{};
  std::array<BpMicroProgram, 197> programs{};
  size_t count = 0;
};
inline constexpr auto bp_micro_tables = [] {
  BpMicroTables out;
  std::array<uint16_t, 8192> ids{};
  for (size_t key = 0; key < 8192; ++key) {
    const auto& p = bp_micro_programs[key];
    if (!p.bits) {
      continue;
    }
    const auto canonical = key & ((size_t{1} << (p.bits - 1)) - 1);
    if (canonical == key) {
      ids[key] = out.count;
      out.programs[out.count++] = p;
    }
    out.headers[key] = p.bits | (ids[canonical] << 4);
  }
  return out;
}();
static_assert(bp_micro_tables.count == 197);

}  // namespace pixie::detail::bp_micro
