#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pixie::tree::experimental::decomposition_detail {

inline constexpr size_t kAbsent = std::numeric_limits<size_t>::max();

inline size_t checked_add(size_t a, size_t b) {
  if (b > std::numeric_limits<size_t>::max() - a) {
    throw std::length_error("Decomposed tree size overflow");
  }
  return a + b;
}

template <class T>
size_t allocated_bytes(const std::vector<T>& values) {
  return values.capacity() * sizeof(T);
}

// All links use the complete byte domain, including 255. Validity is separate.
// Forest roots have no local parent; their common external parent is stored by
// the cover. No fictitious root is inserted in the packed BP stream.
struct DecodedForest {
  std::array<uint8_t, 256> parent{};
  std::array<uint8_t, 256> first_child{};
  std::array<uint8_t, 256> next_sibling{};
  std::bitset<256> has_parent;
  std::bitset<256> has_first_child;
  std::bitset<256> has_next_sibling;
  uint16_t size = 0;

  // Decode exactly count nodes from an arbitrary bit offset. No read crosses
  // the supplied span, even when the final BP bit ends in the final word.
  void decode(std::span<const uint64_t> words, size_t offset, size_t count) {
    if (count > 256 ||
        offset > std::numeric_limits<size_t>::max() - 2 * count) {
      throw std::invalid_argument("Invalid local BP forest size");
    }
    const size_t end = offset + 2 * count;
    if (end / 64 + (end % 64 != 0) > words.size()) {
      throw std::invalid_argument("Truncated local BP forest");
    }
    has_parent.reset();
    has_first_child.reset();
    has_next_sibling.reset();
    std::array<uint8_t, 256> stack{};
    std::array<uint8_t, 257> previous{};
    std::bitset<257> has_previous;
    size_t depth = 0;
    size_t nodes = 0;
    for (size_t bit = offset; bit < end; ++bit) {
      if ((words[bit / 64] >> (bit % 64)) & 1) {
        if (nodes == count) {
          throw std::invalid_argument("Too many local BP opens");
        }
        const auto local = static_cast<uint8_t>(nodes++);
        if (depth) {
          parent[local] = stack[depth - 1];
          has_parent.set(local);
          if (!has_first_child[parent[local]]) {
            first_child[parent[local]] = local;
            has_first_child.set(parent[local]);
          }
        }
        if (has_previous[depth]) {
          next_sibling[previous[depth]] = local;
          has_next_sibling.set(previous[depth]);
        }
        previous[depth] = local;
        has_previous.set(depth);
        stack[depth++] = local;
        has_previous.reset(depth);
      } else {
        if (!depth) {
          throw std::invalid_argument("Unbalanced local BP forest");
        }
        --depth;
      }
    }
    if (depth || nodes != count) {
      throw std::invalid_argument("Incomplete local BP forest");
    }
    size = static_cast<uint16_t>(nodes);
  }
};

struct PreorderRun {
  size_t component;
  uint16_t length;
  uint8_t local_start;
};

// Every group has at most 256 runs of at most 256 nodes. Exclusive starts
// fit in uint16_t, including a final start of 65,535. The terminal total is
// deliberately NOT stored in uint16_t and query keys are never narrowed.
struct RunDirectory {
  std::vector<uint16_t> starts;
  std::vector<uint64_t> samples;

  void build(std::span<const PreorderRun> runs) {
    size_t total = 0;
    size_t base = 0;
    starts.reserve(runs.size());
    samples.reserve(runs.size() / 256 + (runs.size() % 256 != 0));
    for (size_t i = 0; i < runs.size(); ++i) {
      if (!runs[i].length || runs[i].length > 256) {
        throw std::logic_error("Invalid preorder run length");
      }
      if (i % 256 == 0) {
        samples.push_back(total);
        base = total;
      }
      if (total - base > std::numeric_limits<uint16_t>::max()) {
        throw std::logic_error("Preorder bottom-directory overflow");
      }
      starts.push_back(static_cast<uint16_t>(total - base));
      total = checked_add(total, runs[i].length);
    }
  }

  size_t start(size_t run) const { return samples[run / 256] + starts[run]; }

  size_t in_group(size_t group, size_t global) const {
    const size_t begin = group * 256;
    const size_t end = begin + std::min(size_t{256}, starts.size() - begin);
    const size_t relative = global - samples[group];
    return std::upper_bound(starts.begin() + begin, starts.begin() + end,
                            relative) -
           starts.begin() - 1;
  }

  size_t find(size_t global) const {
    const size_t group =
        std::upper_bound(samples.begin(), samples.end(), global) -
        samples.begin() - 1;
    return in_group(group, global);
  }

  size_t bytes() const {
    return allocated_bytes(starts) + allocated_bytes(samples);
  }
};

// A balanced 256-way directory, independent of the source tree height. Bottom
// exclusive counters are narrow; upper global bases are full-width. Each level
// stores only starts, with no padded 256-entry allocation per routing node.
struct BalancedRouting {
  RunDirectory bottom;
  std::vector<std::vector<uint64_t>> levels;
  std::vector<std::array<size_t, 2>> destinations;

  void build(std::span<const PreorderRun> runs,
             std::span<const uint8_t>,
             std::vector<std::array<size_t, 2>> links) {
    bottom.build(runs);
    destinations = std::move(links);
    auto previous = std::span<const uint64_t>(bottom.samples);
    while (previous.size() > 256) {
      std::vector<uint64_t> next;
      next.reserve(previous.size() / 256 + (previous.size() % 256 != 0));
      for (size_t i = 0; i < previous.size(); i += 256) {
        next.push_back(previous[i]);
      }
      levels.push_back(std::move(next));
      previous = levels.back();
    }
  }

  size_t find(size_t global) const {
    size_t group = 0;
    auto descend = [&](const auto& level) {
      const size_t begin = group * 256;
      const size_t end = begin + std::min(size_t{256}, level.size() - begin);
      group =
          std::upper_bound(level.begin() + begin, level.begin() + end, global) -
          level.begin() - 1;
    };
    for (auto level = levels.rbegin(); level != levels.rend(); ++level) {
      descend(*level);
    }
    descend(bottom.samples);
    return bottom.in_group(group, global);
  }

  size_t destination(size_t component, size_t ordinal) const {
    return destinations[component][ordinal];
  }
  size_t start(size_t run) const { return bottom.start(run); }
  size_t height() const {
    return bottom.starts.empty() ? 0 : levels.size() + 2;
  }
  size_t boundary_bytes() const { return allocated_bytes(destinations); }
  size_t macro_bytes() const { return 0; }
  size_t preorder_bytes() const {
    size_t result = bottom.bytes() + allocated_bytes(levels);
    for (const auto& level : levels) {
      result += allocated_bytes(level);
    }
    return result;
  }
};

// Conventional LOUDS: degree one-bits followed by one zero per BFS component.
// Zero-select uses a prefix count every 512 bits, binary search, then at most
// eight word reads. The index owns no pointer into the words: normal tree
// copies and moves remain safe. This modest rank directory is counted, not
// called free.
struct MacroRouting {
  RunDirectory lookup;
  std::vector<uint64_t> words;
  std::vector<uint64_t> zero_samples;
  size_t bits = 0;

  void build(std::span<const PreorderRun> runs,
             std::span<const uint8_t> degrees,
             std::vector<std::array<size_t, 2>>) {
    lookup.build(runs);
    if (degrees.empty()) {
      return;
    }
    bits = checked_add(degrees.size(), degrees.size() - 1);
    words.resize(bits / 64 + (bits % 64 != 0));
    size_t position = 0;
    for (auto degree : degrees) {
      for (size_t i = 0; i < degree; ++i, ++position) {
        words[position / 64] |= uint64_t{1} << (position % 64);
      }
      ++position;
    }
    if (position != bits) {
      throw std::logic_error("Invalid component LOUDS degrees");
    }
    size_t zeros = 0;
    for (size_t word = 0; word < words.size(); ++word) {
      if (word % 8 == 0) {
        zero_samples.push_back(zeros);
      }
      zeros +=
          std::min(size_t{64}, bits - 64 * word) - std::popcount(words[word]);
    }
  }

  size_t select_zero(size_t rank) const {
    const size_t sample =
        std::upper_bound(zero_samples.begin(), zero_samples.end(), rank) -
        zero_samples.begin() - 1;
    size_t remaining = rank - zero_samples[sample];
    const size_t begin = sample * 8;
    const size_t end = begin + std::min(size_t{8}, words.size() - begin);
    for (size_t word = begin; word < end; ++word) {
      uint64_t zeros = ~words[word];
      const size_t valid = std::min(size_t{64}, bits - 64 * word);
      if (valid < 64) {
        zeros &= (uint64_t{1} << valid) - 1;
      }
      const size_t count = std::popcount(zeros);
      if (remaining < count) {
        while (remaining--) {
          zeros &= zeros - 1;
        }
        return 64 * word + std::countr_zero(zeros);
      }
      remaining -= count;
    }
    throw std::logic_error("Invalid macro zero-select");
  }

  size_t destination(size_t component, size_t ordinal) const {
    const size_t start = component ? select_zero(component - 1) + 1 : 0;
    return start - component + 1 + ordinal;
  }
  size_t find(size_t global) const { return lookup.find(global); }
  size_t start(size_t run) const { return lookup.start(run); }
  size_t height() const { return lookup.starts.empty() ? 0 : 2; }
  size_t boundary_bytes() const { return 0; }
  size_t macro_bytes() const {
    return allocated_bytes(words) + allocated_bytes(zero_samples);
  }
  size_t preorder_bytes() const { return lookup.bytes(); }
};

}  // namespace pixie::tree::experimental::decomposition_detail
