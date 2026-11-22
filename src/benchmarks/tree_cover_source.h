#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace pixie::benchmarks {
enum class CoverShape {
  kRandomRecursive,
  kStar,
  kChain,
  kComb,
  kBalanced,
  kCatalan,
  kBinaryCatalan
};

// Source IDs are insertion IDs, not preorder ranks. Registered sizes fit in
// uint32_t. The source is released before any query timing begins.
class CoverSource {
 public:
  CoverSource(size_t n, CoverShape shape, uint64_t seed)
      : first_(n, kNone), next_(n, kNone), degree_(n, 0) {
    std::mt19937_64 random(seed);
    if (shape == CoverShape::kBinaryCatalan && n) {
      // Uniform full-binary preorder with n internal nodes; discard nulls.
      std::vector<int8_t> steps(2 * n + 1, -1);
      std::fill_n(steps.begin(), n, 1);
      std::shuffle(steps.begin(), steps.end(), random);
      int level = 0, minimum = 0;
      size_t start = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        level += steps[i];
        if (level < minimum) {
          minimum = level;
          start = i + 1;
        }
      }
      std::vector<uint32_t> slots{kNone};
      uint32_t id = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        auto p = slots.back();
        slots.pop_back();
        if (steps[(start + i) % steps.size()] > 0) {
          next_[id] = p;
          if (p != kNone) {
            ++degree_[p];
          }
          slots.push_back(id);
          slots.push_back(id++);
        }
      }
    } else if (shape == CoverShape::kCatalan && n) {
      // Cycle lemma: rotate n-1 up and n down steps after the first
      // global minimum. Dropping the last down step gives a uniform Dyck path.
      std::vector<int8_t> steps(2 * n - 1, -1);
      std::fill_n(steps.begin(), n - 1, 1);
      std::shuffle(steps.begin(), steps.end(), random);
      int level = 0, minimum = 0;
      size_t start = 0;
      for (size_t i = 0; i < steps.size(); ++i) {
        level += steps[i];
        if (level < minimum) {
          minimum = level;
          start = i + 1;
        }
      }
      uint32_t node = 0, id = 0;
      for (size_t i = 0; i + 1 < steps.size(); ++i) {
        if (steps[(start + i) % steps.size()] == 1) {
          next_[++id] = node;
          ++degree_[node];
          node = id;
        } else {
          node = next_[node];
        }
      }
    }
    for (size_t i = 1; i < n; ++i) {
      size_t parent = 0;
      switch (shape) {
        case CoverShape::kCatalan:
        case CoverShape::kBinaryCatalan:
          continue;
        case CoverShape::kRandomRecursive:
          parent = random() % i;
          break;
        case CoverShape::kStar:
          break;
        case CoverShape::kChain:
          parent = i - 1;
          break;
        case CoverShape::kComb:
          parent = (i - 1) & ~size_t{1};
          break;
        case CoverShape::kBalanced:
          parent = (i - 1) / 2;
          break;
      }
      next_[i] = static_cast<uint32_t>(parent);
      ++degree_[parent];
    }
    for (size_t i = n; i-- > 1;) {
      const auto parent = next_[i];
      next_[i] = first_[parent];
      first_[parent] = static_cast<uint32_t>(i);
    }
  }
  size_t size() const { return first_.size(); }
  uint32_t root() const { return 0; }
  size_t degree(uint32_t node) const { return degree_[node]; }
  uint32_t first_child(uint32_t node) const { return first_[node]; }
  uint32_t next_sibling(uint32_t node) const { return next_[node]; }

 private:
  static constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();
  std::vector<uint32_t> first_, next_, degree_;
};

}  // namespace pixie::benchmarks
