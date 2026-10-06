#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace pixie::benchmarks {

// Deterministic inclusive-range queries. Like the RMQ generator, combine a
// precomputed pool with a random shift per pass and a per-query offset. RNG
// work stays outside the hot loop except once at each pool wrap. Shifting
// avoids replaying the pool verbatim; it does not guarantee unique queries.
class ShiftedQueryPool {
 public:
  static constexpr std::size_t kPoolBytes = 64 * 1024;
  static constexpr std::size_t kQueryCount = kPoolBytes / sizeof(std::size_t);

  ShiftedQueryPool(std::size_t first, std::size_t last, std::uint64_t seed)
      : first_(first), range_(last - first + 1), rng_state_(seed) {
    if (first > last) {
      throw std::invalid_argument("Empty query range");
    }
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::size_t> distribution(0, last - first);
    pool_.resize(kQueryCount);
    for (auto& query : pool_) {
      query = distribution(rng);
    }
    shift_ = next_shift();
  }

  std::size_t next() {
    if (index_ == pool_.size()) {
      index_ = 0;
      shift_ = next_shift();
    }
    const auto offset = add_mod(shift_, phase_);
    const auto query = first_ + add_mod(pool_[index_++], offset);
    if (++phase_ == range_) {
      phase_ = 0;
    }
    return query;
  }

 private:
  std::size_t next_shift() {
    std::uint64_t z = (rng_state_ += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    const auto value = static_cast<std::size_t>(z ^ (z >> 31));
    return range_ == 0 ? value : value % range_;
  }

  // Inputs are already reduced. Avoid overflow even for ranges near SIZE_MAX.
  // A zero range denotes the complete size_t domain, whose addition wraps.
  std::size_t add_mod(std::size_t a, std::size_t b) const {
    if (range_ == 0) {
      return a + b;
    }
    return a >= range_ - b ? a - (range_ - b) : a + b;
  }

  std::size_t first_;
  std::size_t range_;
  std::uint64_t rng_state_;
  std::vector<std::size_t> pool_;
  std::size_t index_ = 0;
  std::size_t phase_ = 0;
  std::size_t shift_ = 0;
};

}  // namespace pixie::benchmarks
