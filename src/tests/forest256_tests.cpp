#include <gtest/gtest.h>
#include <pixie/tree/detail/forest256_cache.h>

#include <cstdlib>
#include <limits>
#include <vector>

namespace {
using namespace pixie::tree::detail::forest256;

// Construct the two dense streams directly from BFS adjacency, without the
// production padded encoders. Leave every neighboring bit dirty.
std::vector<uint64_t> pack(Codec codec,
                           size_t offset,
                           size_t roots,
                           const std::vector<uint16_t>& degrees) {
  const size_t n = degrees.size();
  std::vector<uint64_t> words((offset + 2 * n + 63) / 64, ~uint64_t{});
  auto clear = [&](size_t i) {
    words[(offset + i) / 64] &= ~(uint64_t{1} << ((offset + i) % 64));
  };
  if (codec == Codec::louds) {
    size_t bit = roots;
    clear(bit++);
    for (auto d : degrees) {
      bit += d;
      if (bit < 2 * n) {
        clear(bit++);
      }
    }
  } else {
    size_t end = roots;
    for (size_t i = 0; i + 1 < roots; ++i) {
      clear(n + i);
    }
    for (size_t i = 0; i < n; ++i) {
      if (!degrees[i]) {
        clear(i);
      }
      for (size_t j = 1; j < degrees[i]; ++j) {
        clear(n + end + j - 1);
      }
      end += degrees[i];
    }
  }
  return words;
}

void check(size_t roots, const std::vector<uint16_t>& degrees, size_t offset) {
  const size_t n = degrees.size();
  std::vector<uint8_t> parents(roots, kNoParent);
  for (size_t i = 0; i < n; ++i) {
    parents.insert(parents.end(), degrees[i], i);
  }
  for (auto codec : {Codec::louds, Codec::twin}) {
    auto words = pack(codec, offset, roots, degrees);
    for (auto decode : {unpack_dense, unpack_dense_scalar}) {
      LevelCsr256 csr;
      decode(codec, words, offset, n, csr);
      ASSERT_EQ(csr.root_count(), roots);
      size_t next = roots;
      for (size_t i = 0; i < n; ++i) {
        if (csr.parent[i] != parents[i] || csr.offset_m1[i] + 1u != next ||
            csr.degree(i) != degrees[i] ||
            csr.is_last_child(i) !=
                (i + 1 == n || parents[i + 1] != parents[i])) {
          FAIL() << "codec=" << int(codec) << " n=" << n << " r=" << roots
                 << " offset=" << offset << " i=" << i;
        }
        next += degrees[i];
      }
      ASSERT_EQ(csr.offset_m1[n] + 1u, n);
    }
  }
}

TEST(Forest256Dense, CountsRootsOffsetsAndSentinels) {
  // The full displacement cross product is useful under ASan; the normal
  // suite keeps every (count, roots) and all offsets on sentinel boundaries.
  const bool all_offsets = std::getenv("FOREST256_ALL_OFFSETS") != nullptr;
  for (size_t n = 1; n <= 256; ++n) {
    for (size_t r = 1; r <= n; ++r) {
      std::vector<uint16_t> degrees(n);
      // Each of the first n-r nodes has a child; in particular 256/1
      // exercises parent 254, while n isolated roots exercise parent 255.
      std::fill_n(degrees.begin(), n - r, 1);
      if (all_offsets || r == 1 || r == n || r + 1 == n) {
        for (size_t shift = 0; shift < 64; ++shift) {
          check(r, degrees, 64 + shift);
        }
      } else {
        check(r, degrees, 64 + (n + r) % 64);
      }
    }
    std::vector<uint16_t> star(n);
    star[0] = n - 1;
    check(1, star, n % 64);
  }
}

TEST(Forest256Dense, ExhaustiveSmallForests) {
  for (size_t n = 1; n <= 7; ++n) {
    std::vector<uint16_t> d(n);
    for (size_t r = 1; r <= n; ++r) {
      auto visit = [&](auto&& self, size_t i, size_t end) -> void {
        if (i == n) {
          check(r, d, 13);
          return;
        }
        if (i == end) {
          return;
        }
        for (size_t k = 0; k <= n - end; ++k) {
          d[i] = k;
          self(self, i + 1, end + k);
        }
      };
      visit(visit, 0, r);
    }
  }
}

TEST(Forest256Dense, RejectsInvalidExtents) {
  LevelCsr256 csr;
  uint64_t words[1]{};
  for (auto c : {Codec::louds, Codec::twin}) {
    EXPECT_THROW(unpack_dense(c, words, 63, 1, csr), std::out_of_range);
    EXPECT_THROW(
        unpack_dense(c, words, std::numeric_limits<size_t>::max(), 256, csr),
        std::out_of_range);
    EXPECT_THROW(unpack_dense(c, words, 0, 0, csr), std::invalid_argument);
    EXPECT_THROW(unpack_dense(c, words, 0, 257, csr), std::invalid_argument);
  }
  EXPECT_THROW(unpack_dense(Codec::louds, words, 0, 1, csr),
               std::invalid_argument);
}

TEST(Forest256CacheTest, BudgetEvictionLifetimeAndRebinding) {
  constexpr auto c = Codec::louds;
  auto a = pack(c, 0, 1, {3, 0, 0, 0});
  auto b = pack(c, 0, 2, {0, 1, 0});
  a.push_back(b.front());
  const size_t s = Forest256Cache::slot_bytes();
  for (size_t budget : {size_t{0}, s - 1, s, 2 * s - 1, 2 * s}) {
    Forest256Cache cache(a, budget);
    Forest256Cache::Stats stats;
    EXPECT_EQ(cache.allocated_bytes(), budget / s * s);
    auto saved = cache.lookup_counted(0, c, 0, 4, 0, stats);
    EXPECT_EQ(saved.first, 1);
    auto nested = [&](size_t child) {
      auto other = cache.lookup_counted(1, c, 64, 3, 2, stats);
      EXPECT_EQ(other.parent, 1);
      EXPECT_EQ(cache.lookup_counted(0, c, 0, 4, child, stats).parent, 0);
    };
    for (size_t i = 0; i < saved.degree; ++i) {
      nested(saved.first + i);
    }
    EXPECT_EQ(stats.misses, budget < 2 * s ? 7 : 2);
    EXPECT_EQ(stats.hits, budget < 2 * s ? 0 : 5);
    EXPECT_EQ(stats.evictions, budget >= s && budget < 2 * s ? 6 : 0);
    EXPECT_EQ(stats.packed_bits, budget < 2 * s ? 50 : 14);
    Forest256Cache copy = cache;
    copy.bind(b);
    EXPECT_EQ(copy.lookup(0, c, 0, 3, 0).roots, 2);
    auto moved = std::move(copy);
    EXPECT_EQ(copy.slot_count(), 0);
    EXPECT_EQ(moved.lookup(0, c, 0, 3, 2).parent, 1);
    EXPECT_THROW(cache.lookup(2, c, 999, 3, 0), std::out_of_range);
    EXPECT_EQ(cache.lookup(0, c, 0, 4, 0).roots, 1);
    moved = cache;
    EXPECT_EQ(moved.lookup(0, c, 0, 4, 0).degree, 3);
  }
}

TEST(Forest256CacheTest, CollisionReplacesOnlyItsSlot) {
  constexpr auto c = Codec::louds;
  auto words = pack(c, 0, 1, {2, 0, 0});
  words.push_back(pack(c, 0, 2, {0, 1, 0}).front());
  words.push_back(pack(c, 0, 3, {0, 0, 0}).front());
  Forest256Cache cache(words, 2 * Forest256Cache::slot_bytes());
  Forest256Cache::Stats stats;
  auto roots = [&](size_t unit) {
    return cache.lookup_counted(unit, c, unit * 64, 3, 0, stats).roots;
  };
  EXPECT_EQ(roots(0), 1);
  EXPECT_EQ(roots(1), 2);
  EXPECT_EQ(roots(2), 3);  // Collides with 0, not 1.
  EXPECT_EQ(stats.evictions, 1u);
  EXPECT_EQ(roots(1), 2);
  EXPECT_EQ(stats.hits, 1u);
  EXPECT_EQ(roots(0), 1);
  EXPECT_EQ(stats.misses, 4u);
  EXPECT_EQ(stats.evictions, 2u);
}

TEST(Forest256CacheTest, RejectedReplacementLeavesWarmEntryUsable) {
  constexpr auto c = Codec::louds;
  auto words = pack(c, 0, 1, {2, 0, 0});
  words.push_back(0);  // Fits in the stream, but has invalid population.
  Forest256Cache cache(words, Forest256Cache::slot_bytes());
  ASSERT_EQ(cache.lookup(0, c, 0, 3, 0).degree, 2);
  EXPECT_THROW(cache.lookup(1, c, 64, 3, 0), std::invalid_argument);
  Forest256Cache::Stats probe;
  EXPECT_EQ(cache.lookup_counted(0, c, 0, 3, 0, probe).degree, 2);
  EXPECT_EQ(probe.hits, 1u);
}

TEST(Forest256CacheTest, SnapshotSurvivesEvictionAndRebinding) {
  for (auto c : {Codec::louds, Codec::twin}) {
    auto words = pack(c, 0, 1, {2, 0, 0});
    words.push_back(pack(c, 0, 2, {0, 1, 0}).front());
    auto rebound = pack(c, 0, 3, {0, 0, 0});
    for (size_t budget : {size_t{0}, Forest256Cache::slot_bytes()}) {
      Forest256Cache cache(words, budget);
      LevelCsr256 saved;
      cache.snapshot(0, c, 0, 3, saved);
      ASSERT_EQ(cache.lookup(1, c, 64, 3, 2).parent, 1);
      cache.bind(rebound);
      ASSERT_EQ(cache.lookup(0, c, 0, 3, 0).roots, 3);
      EXPECT_EQ(saved.root_count(), 1u);
      ASSERT_EQ(saved.degree(0), 2u);
      EXPECT_EQ(saved.child(0, 1), 2);
      EXPECT_EQ(saved.parent_of(2), 0);
      EXPECT_FALSE(saved.is_last_child(1));
    }
  }
}

TEST(Forest256CacheTest, LastOf256RootsHasEmptyChildRangeAtEnd) {
  for (auto c : {Codec::louds, Codec::twin}) {
    auto words = pack(c, 0, 256, std::vector<uint16_t>(256));
    Forest256Cache cache(words, Forest256Cache::slot_bytes());
    const auto last = cache.lookup(0, c, 0, 256, 255);
    EXPECT_EQ(last.roots, 256);
    EXPECT_EQ(last.parent, kNoParent);
    EXPECT_TRUE(last.last);
    EXPECT_EQ(last.first, 256);  // [256, 256), not a child to dereference.
    EXPECT_EQ(last.degree, 0);
    EXPECT_THROW(cache.lookup(0, c, 0, 256, 256), std::out_of_range);
  }
}

TEST(Forest256CacheTest, SameAddressRebindInvalidatesWarmSlot) {
  for (auto c : {Codec::louds, Codec::twin}) {
    auto stream = pack(c, 0, 1, {2, 0, 0});
    auto replacement = pack(c, 0, 3, {0, 0, 0});
    ASSERT_EQ(stream.size(), 1u);
    ASSERT_EQ(replacement.size(), 1u);
    Forest256Cache reader(stream, Forest256Cache::slot_bytes());
    Forest256Cache::Stats stats;
    auto root = [&] { return reader.lookup_counted(0, c, 0, 3, 0, stats); };
    ASSERT_EQ(root().degree, 2);
    ASSERT_EQ(root().roots, 1);
    ASSERT_EQ(stats.hits, 1u);

    stream.front() = replacement.front();
    reader.bind(stream);
    const auto fresh = root();
    EXPECT_EQ(fresh.degree, 0);
    EXPECT_EQ(fresh.roots, 3);
    EXPECT_EQ(stats.misses, 2u);
    EXPECT_EQ(stats.evictions, 0u);
  }
}

TEST(Forest256CacheTest, AssignmentAdoptsCapacityAndMoveEmptiesSource) {
  constexpr auto c = Codec::louds;
  const size_t s = Forest256Cache::slot_bytes();
  auto a = pack(c, 0, 1, {2, 0, 0});
  auto b = pack(c, 0, 3, {0, 0, 0});
  Forest256Cache source(a, 3 * s + 1), dest(b, s), zero(b, s - 1);
  ASSERT_EQ(source.lookup(0, c, 0, 3, 0).degree, 2);
  ASSERT_EQ(dest.lookup(0, c, 0, 3, 0).degree, 0);

  dest = source;
  EXPECT_EQ(dest.slot_count(), 3u);
  EXPECT_EQ(dest.allocated_bytes(), 3 * s);
  Forest256Cache::Stats stats;
  EXPECT_EQ(dest.lookup_counted(0, c, 0, 3, 0, stats).degree, 2);
  EXPECT_EQ(stats.misses, 1u);  // A copy starts cold.
  source = zero;
  EXPECT_EQ(source.allocated_bytes(), 0u);
  EXPECT_EQ(source.lookup(0, c, 0, 3, 0).degree, 0);

  source = std::move(dest);
  EXPECT_EQ(source.allocated_bytes(), 3 * s);
  EXPECT_EQ(source.lookup_counted(0, c, 0, 3, 0, stats).degree, 2);
  EXPECT_EQ(stats.hits, 1u);  // A move transfers warm slots.
  EXPECT_EQ(dest.slot_count(), 0u);
  EXPECT_EQ(dest.allocated_bytes(), 0u);
  EXPECT_THROW(dest.lookup(0, c, 0, 3, 0), std::out_of_range);
  dest.bind(b);
  EXPECT_EQ(dest.lookup(0, c, 0, 3, 0).roots, 3);
}
}  // namespace
