#include <gtest/gtest.h>
#include <pixie/rmm_tree.h>
#include <reference_implementations/naive_rmm_tree.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;

static std::string bits_to_parens(const std::string& bits) {
  std::string s;
  s.reserve(bits.size());
  for (char c : bits) {
    s.push_back(c == '1' ? '(' : ')');
  }
  return s;
}

static std::string vecbits_to_string(const std::vector<uint8_t>& v) {
  std::string s;
  s.resize(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    s[i] = v[i] ? '1' : '0';
  }
  return s;
}

static std::vector<std::uint64_t> pack_words_lsb_first(
    const std::string& bits) {
  const size_t n = bits.size();
  std::vector<std::uint64_t> w((n + 63) / 64, 0);
  for (size_t i = 0; i < n; ++i) {
    if (bits[i] == '1') {
      w[i >> 6] |= (std::uint64_t(1) << (i & 63));
    }
  }
  return w;
}

static std::string random_bits(std::mt19937_64& rng, size_t n) {
  std::uniform_int_distribution<int> b01(0, 1);
  std::string s;
  s.resize(n);
  for (size_t i = 0; i < n; i++) {
    s[i] = char('0' + b01(rng));
  }
  return s;
}

static std::string random_dyck_bits(std::mt19937_64& rng, size_t m) {
  std::size_t L = m << 1;
  std::string s;
  s.resize(L);
  if (L == 0) {
    return s;
  }
  std::size_t opens_left = m;
  std::size_t closes_left = m;
  int balance = 0;
  std::uniform_real_distribution<double> U(0., 1.);
  for (std::size_t i = 0; i < L; ++i) {
    if (opens_left == 0) {
      s[i] = '0';
      --closes_left;
      --balance;
      continue;
    }

    if (closes_left == 0) {
      s[i] = '1';
      --opens_left;
      ++balance;
      continue;
    }

    if (balance == 0) {
      s[i] = '1';
      --opens_left;
      ++balance;
      continue;
    }
    if (U(rng) < 0.5) {
      s[i] = '1';
      --opens_left;
      ++balance;
    } else {
      s[i] = '0';
      --closes_left;
      --balance;
    }
  }
  return s;
}

static constexpr uint64_t kSeed = 42;
static constexpr size_t kRandomCases = 20;
static constexpr size_t kOpsPerCase = 600;
static constexpr size_t kMaxBits = 20000;
static constexpr size_t kLongOpsPerCase = 2000;
static constexpr size_t kLongMaxBits = 65536;

static void run_case_and_compare(const std::string& bits,
                                 std::mt19937_64& rng,
                                 size_t ops_per_case,
                                 uint64_t seed) {
  const bool use_words = std::uniform_int_distribution<int>(0, 1)(rng);
  pixie::RmMTree rm;
  NaiveRmM nv;
  auto words = pack_words_lsb_first(bits);
  if (use_words) {
    rm = pixie::RmMTree(std::span<const std::uint64_t>(words), bits.size());
    nv = NaiveRmM(words, bits.size());
  } else {
    rm = pixie::RmMTree(std::span<const std::uint64_t>(words), bits.size());
    nv = NaiveRmM(bits);
  }

  const size_t N = bits.size();
  const size_t ones = nv.rank1(N);
  const size_t zeros = N - ones;
  const size_t pairs10 = (N >= 2 ? nv.rank10(N) : 0);

  std::uniform_int_distribution<size_t> pos_i(0, N);
  std::uniform_int_distribution<size_t> pos_i_nz(0, N ? N - 1 : 0);
  std::uniform_int_distribution<int> d_dist(-(int)std::min<size_t>(N, 200),
                                            (int)std::min<size_t>(N, 200));

  for (size_t q = 0; q < ops_per_case; ++q) {
    int which = std::uniform_int_distribution<int>(0, 17)(rng);

    size_t i = 0, j = 0;
    if (N > 0) {
      i = pos_i_nz(rng);
      j = pos_i_nz(rng);
      if (i > j) {
        std::swap(i, j);
      }
    }

    auto trace_header = ::testing::Message()
                        << "\nseed=" << seed << " | N=" << N
                        << " | use_words=" << (use_words ? 1 : 0)
                        << "\nbits:   " << bits
                        << "\nparens: " << bits_to_parens(bits)
                        << "\nwhich=" << which << " (0..17)\n";

    switch (which) {
      case 0: {  // rank1
        size_t x = pos_i(rng);
        auto a = nv.rank1(x), b = rm.rank1(x);
        EXPECT_EQ(a, b) << trace_header << "op=rank1(" << x << ")";
      } break;

      case 1: {  // rank0
        size_t x = pos_i(rng);
        auto a = nv.rank0(x), b = rm.rank0(x);
        EXPECT_EQ(a, b) << trace_header << "op=rank0(" << x << ")";
      } break;

      case 2: {  // select1
        size_t k = std::uniform_int_distribution<size_t>(0, ones + 3)(rng);
        auto a = nv.select1(k), b = rm.select1(k);
        EXPECT_EQ(a, b) << trace_header << "op=select1(" << k << ")";
      } break;

      case 3: {  // select0
        size_t k = std::uniform_int_distribution<size_t>(0, zeros + 3)(rng);
        auto a = nv.select0(k), b = rm.select0(k);
        EXPECT_EQ(a, b) << trace_header << "op=select0(" << k << ")";
      } break;

      case 4: {  // rank10
        size_t x =
            (N >= 2 ? std::uniform_int_distribution<size_t>(0, N)(rng) : 0);
        auto a = nv.rank10(x), b = rm.rank10(x);
        EXPECT_EQ(a, b) << trace_header << "op=rank10(" << x << ")";
      } break;

      case 5: {  // select10
        size_t k = std::uniform_int_distribution<size_t>(0, pairs10 + 3)(rng);
        auto a = nv.select10(k), b = rm.select10(k);
        EXPECT_EQ(a, b) << trace_header << "op=select10(" << k << ")";
      } break;

      case 6: {  // excess
        size_t x = pos_i(rng);
        auto a = nv.excess(x), b = rm.excess(x);
        EXPECT_EQ(a, b) << trace_header << "op=excess(" << x << ")";
      } break;

      case 7: {  // fwdsearch
        if (N == 0) {
          break;
        }
        size_t start = pos_i_nz(rng);
        int d = d_dist(rng);
        auto a = nv.fwdsearch(start, d), b = rm.fwdsearch(start, d);
        EXPECT_EQ(a, b) << trace_header << "op=fwdsearch(" << start << "," << d
                        << ")";
      } break;

      case 8: {  // bwdsearch
        if (N == 0) {
          break;
        }
        size_t start = pos_i_nz(rng);
        int d = d_dist(rng);
        auto a = nv.bwdsearch(start, d), b = rm.bwdsearch(start, d);
        EXPECT_EQ(a, b) << trace_header << "op=bwdsearch(" << start << "," << d
                        << ")";
      } break;

      case 9: {  // range_min_query_pos
        if (N == 0) {
          break;
        }
        auto a = nv.range_min_query_pos(i, j), b = rm.range_min_query_pos(i, j);
        EXPECT_EQ(a, b) << trace_header << "op=range_min_query_pos(" << i << ","
                        << j << ")";
      } break;

      case 10: {  // range_min_query_val
        if (N == 0) {
          break;
        }
        auto a = nv.range_min_query_val(i, j), b = rm.range_min_query_val(i, j);
        EXPECT_EQ(a, b) << trace_header << "op=range_min_query_val(" << i << ","
                        << j << ")";
      } break;

      case 11: {  // mincount
        if (N == 0) {
          break;
        }
        auto a = nv.mincount(i, j), b = rm.mincount(i, j);
        EXPECT_EQ(a, b) << trace_header << "op=mincount(" << i << "," << j
                        << ")";
      } break;

      case 12: {  // minselect
        if (N == 0) {
          break;
        }
        size_t cnt = nv.mincount(i, j);
        size_t k = cnt == 0
                       ? 1
                       : std::uniform_int_distribution<size_t>(1, cnt + 1)(rng);
        auto a = nv.minselect(i, j, k), b = rm.minselect(i, j, k);
        EXPECT_EQ(a, b) << trace_header << "op=minselect(" << i << "," << j
                        << "," << k << ")";
      } break;

      case 13: {  // range_max_query_pos
        if (N == 0) {
          break;
        }
        auto a = nv.range_max_query_pos(i, j), b = rm.range_max_query_pos(i, j);
        EXPECT_EQ(a, b) << trace_header << "op=range_max_query_pos(" << i << ","
                        << j << ")";
      } break;

      case 14: {  // range_max_query_val
        if (N == 0) {
          break;
        }
        auto a = nv.range_max_query_val(i, j), b = rm.range_max_query_val(i, j);
        EXPECT_EQ(a, b) << trace_header << "op=range_max_query_val(" << i << ","
                        << j << ")";
      } break;

      case 15: {  // close
        if (N == 0) {
          break;
        }
        size_t x = pos_i_nz(rng);
        auto a = nv.close(x), b = rm.close(x);
        EXPECT_EQ(a, b) << trace_header << "op=close(" << x << ")";
      } break;

      case 16: {  // open
        if (N == 0) {
          break;
        }
        size_t x = pos_i(rng);
        auto a = nv.open(x), b = rm.open(x);
        EXPECT_EQ(a, b) << trace_header << "op=open(" << x << ")";
      } break;

      case 17: {  // enclose
        if (N == 0) {
          break;
        }
        size_t x = pos_i(rng);
        auto a = nv.enclose(x), b = rm.enclose(x);
        EXPECT_EQ(a, b) << trace_header << "op=enclose(" << x << ")";
      } break;
    }
  }
}

class RmMRandomTest : public ::testing::Test {
 protected:
  std::mt19937_64 rng{kSeed};
};

TEST_F(RmMRandomTest, RandomMixedBits) {
  std::uniform_int_distribution<int> len_u(1, (int)kMaxBits);
  for (size_t t = 0; t < kRandomCases; ++t) {
    const size_t n = (size_t)len_u(rng);
    const std::string bits = random_bits(rng, n);
    run_case_and_compare(bits, rng, kOpsPerCase, kSeed);
  }
}

TEST_F(RmMRandomTest, RandomDyckBits) {
  std::uniform_int_distribution<int> len_even(0, (int)(kMaxBits / 2));
  for (size_t t = 0; t < kRandomCases; ++t) {
    size_t m = 1 + (size_t)len_even(rng);
    std::string bits = random_dyck_bits(rng, m);
    if (bits.empty()) {
      bits = "10";
    }
    run_case_and_compare(bits, rng, kOpsPerCase, kSeed);
  }
}

TEST_F(RmMRandomTest, ShortInputs) {
  for (size_t n = 1; n <= 8; ++n) {
    const size_t total = (n ? (1ull << n) : 1ull);
    for (size_t mask = 0; mask < total; ++mask) {
      std::string bits;
      bits.resize(n);
      for (size_t i = 0; i < n; ++i) {
        bits[i] = ((mask >> i) & 1ull) ? '1' : '0';
      }
      {
        SCOPED_TRACE(::testing::Message()
                     << "short-inputs:string bits=" << bits);
        run_case_and_compare(bits, rng, /*ops_per_case=*/200, kSeed);
      }
      {
        SCOPED_TRACE(::testing::Message()
                     << "short-inputs:words  bits=" << bits);
        run_case_and_compare(bits, rng, /*ops_per_case=*/200, kSeed);
      }
    }
  }
}

static void run_case_repeated_span_construction(const std::string& bits,
                                                std::mt19937_64& rng) {
  const size_t N = bits.size();
  auto span_words = pack_words_lsb_first(bits);
  pixie::RmMTree rm_s(std::span<const std::uint64_t>(span_words), N);
  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm_w(std::span<const std::uint64_t>(words), N);

  std::uniform_int_distribution<size_t> pos_i(0, N);
  std::uniform_int_distribution<size_t> pos_i_nz(0, N ? N - 1 : 0);
  std::uniform_int_distribution<int> d_dist(-(int)std::min<size_t>(N, 200),
                                            (int)std::min<size_t>(N, 200));

  for (int t = 0; t < 128; ++t) {
    int which = std::uniform_int_distribution<int>(0, 10)(rng);
    size_t i = 0, j = 0;
    if (N > 0) {
      i = pos_i_nz(rng);
      j = pos_i_nz(rng);
      if (i > j) {
        std::swap(i, j);
      }
    }
    switch (which) {
      case 0: {
        size_t x = pos_i(rng);
        EXPECT_EQ(rm_s.rank1(x), rm_w.rank1(x));
        break;
      }
      case 1: {
        size_t x = pos_i(rng);
        EXPECT_EQ(rm_s.rank0(x), rm_w.rank0(x));
        break;
      }
      case 2: {
        size_t k = std::uniform_int_distribution<size_t>(0, N + 3)(rng);
        EXPECT_EQ(rm_s.select1(k), rm_w.select1(k));
        break;
      }
      case 3: {
        size_t k = std::uniform_int_distribution<size_t>(0, N + 3)(rng);
        EXPECT_EQ(rm_s.select0(k), rm_w.select0(k));
        break;
      }
      case 4: {
        size_t x =
            (N >= 2 ? std::uniform_int_distribution<size_t>(0, N)(rng) : 0);
        EXPECT_EQ(rm_s.rank10(x), rm_w.rank10(x));
        break;
      }
      case 5: {
        size_t k = std::uniform_int_distribution<size_t>(0, N)(rng);
        EXPECT_EQ(rm_s.select10(k), rm_w.select10(k));
        break;
      }
      case 6: {
        size_t x = pos_i(rng);
        EXPECT_EQ(rm_s.excess(x), rm_w.excess(x));
        break;
      }
      case 7: {
        if (N == 0) {
          break;
        }
        size_t start = pos_i_nz(rng);
        int d = d_dist(rng);
        EXPECT_EQ(rm_s.fwdsearch(start, d), rm_w.fwdsearch(start, d));
        break;
      }
      case 8: {
        if (N == 0) {
          break;
        }
        size_t start = pos_i_nz(rng);
        int d = d_dist(rng);
        EXPECT_EQ(rm_s.bwdsearch(start, d), rm_w.bwdsearch(start, d));
        break;
      }
      case 9: {
        if (N == 0) {
          break;
        }
        EXPECT_EQ(rm_s.range_min_query_pos(i, j),
                  rm_w.range_min_query_pos(i, j));
        break;
      }
      case 10: {
        if (N == 0) {
          break;
        }
        EXPECT_EQ(rm_s.range_max_query_pos(i, j),
                  rm_w.range_max_query_pos(i, j));
        break;
      }
    }
  }
}

TEST_F(RmMRandomTest, RepeatedSpanConstruction) {
  std::uniform_int_distribution<int> len_u(0, (int)kMaxBits);
  for (size_t t = 0; t < kRandomCases; ++t) {
    const size_t n = (size_t)len_u(rng);
    const std::string bits = random_bits(rng, n);
    run_case_repeated_span_construction(bits, rng);
  }
}

TEST(RmMEdgeCases, EmptyInput) {
  std::vector<std::uint64_t> words;
  pixie::RmMTree rm(std::span<const std::uint64_t>(words), 0);
  NaiveRmM nv(std::string{});
  EXPECT_EQ(rm.rank1(0), nv.rank1(0));
  EXPECT_EQ(rm.rank0(0), nv.rank0(0));
  EXPECT_EQ(rm.rank10(0), nv.rank10(0));
  EXPECT_EQ(rm.select1(1), nv.select1(1));
  EXPECT_EQ(rm.select0(1), nv.select0(1));
  EXPECT_EQ(rm.fwdsearch(0, 0), nv.fwdsearch(0, 0));
  EXPECT_EQ(rm.bwdsearch(0, 0), nv.bwdsearch(0, 0));
  EXPECT_EQ(rm.range_min_query_pos(0, 0), nv.range_min_query_pos(0, 0));
  EXPECT_EQ(rm.range_max_query_pos(0, 0), nv.range_max_query_pos(0, 0));
}

static void expect_rank_select_equal(const pixie::RmMTree& rm,
                                     const NaiveRmM& nv,
                                     size_t n) {
  for (size_t x = 0; x <= n; ++x) {
    EXPECT_EQ(rm.rank1(x), nv.rank1(x)) << "rank1 x=" << x;
    EXPECT_EQ(rm.rank0(x), nv.rank0(x)) << "rank0 x=" << x;
    EXPECT_EQ(rm.rank10(x), nv.rank10(x)) << "rank10 x=" << x;
  }

  const size_t ones = nv.rank1(n);
  const size_t zeros = n - ones;
  const size_t pairs10 = (n >= 2 ? nv.rank10(n) : 0);

  for (size_t k = 1; k <= ones + 1; ++k) {
    EXPECT_EQ(rm.select1(k), nv.select1(k)) << "select1 k=" << k;
  }
  for (size_t k = 1; k <= zeros + 1; ++k) {
    EXPECT_EQ(rm.select0(k), nv.select0(k)) << "select0 k=" << k;
  }
  for (size_t k = 1; k <= pairs10 + 1; ++k) {
    EXPECT_EQ(rm.select10(k), nv.select10(k)) << "select10 k=" << k;
  }
}

static void expect_range_ops_equal(const pixie::RmMTree& rm,
                                   const NaiveRmM& nv,
                                   size_t n) {
  if (n == 0) {
    return;
  }
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<size_t> pos(0, n - 1);
  std::uniform_int_distribution<size_t> k_dist;
  for (int t = 0; t < 512; ++t) {
    size_t i = pos(rng);
    size_t j = pos(rng);
    if (i > j) {
      std::swap(i, j);
    }

    EXPECT_EQ(rm.range_min_query_pos(i, j), nv.range_min_query_pos(i, j));
    EXPECT_EQ(rm.range_min_query_val(i, j), nv.range_min_query_val(i, j));
    EXPECT_EQ(rm.range_max_query_pos(i, j), nv.range_max_query_pos(i, j));
    EXPECT_EQ(rm.range_max_query_val(i, j), nv.range_max_query_val(i, j));

    size_t cnt = nv.mincount(i, j);
    EXPECT_EQ(rm.mincount(i, j), cnt);
    k_dist.param(std::uniform_int_distribution<size_t>::param_type(1, cnt + 1));
    size_t k = k_dist(rng);
    EXPECT_EQ(rm.minselect(i, j, k), nv.minselect(i, j, k));
  }
}

TEST(RmMEdgeCases, MultiwordPattern10AcrossWordBoundaries) {
  const size_t n = 640;
  std::string bits(n, '1');

  for (size_t i = 0; i + 1 < n; i += 3) {
    bits[i] = '1';
    bits[i + 1] = '0';
  }
  for (size_t boundary = 63; boundary + 1 < n; boundary += 64) {
    bits[boundary] = '1';
    bits[boundary + 1] = '0';
  }

  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm(std::span<const std::uint64_t>(words), bits.size(),
                    /*leaf_block_bits=*/256);
  NaiveRmM nv(bits);

  expect_rank_select_equal(rm, nv, n);
  expect_range_ops_equal(rm, nv, n);
}

TEST(RmMEdgeCases, PartialLastLeafSelects) {
  const size_t n = 600;

  std::string mostly_zero(n, '0');
  for (size_t i = 576; i < n; ++i) {
    mostly_zero[i] = '1';
  }
  auto mostly_zero_words = pack_words_lsb_first(mostly_zero);
  pixie::RmMTree rm_select1(std::span<const std::uint64_t>(mostly_zero_words),
                            mostly_zero.size(),
                            /*leaf_block_bits=*/256);
  NaiveRmM nv_select1(mostly_zero);
  expect_rank_select_equal(rm_select1, nv_select1, n);

  std::string mostly_one(n, '1');
  for (size_t i = 576; i < n; ++i) {
    mostly_one[i] = '0';
  }
  auto mostly_one_words = pack_words_lsb_first(mostly_one);
  pixie::RmMTree rm_select0(std::span<const std::uint64_t>(mostly_one_words),
                            mostly_one.size(),
                            /*leaf_block_bits=*/256);
  NaiveRmM nv_select0(mostly_one);
  expect_rank_select_equal(rm_select0, nv_select0, n);
}

TEST(RmMEdgeCases, Select10OnIncompleteInternalNode) {
  constexpr size_t leaf_block_bits = 256;
  const size_t n = (leaf_block_bits * 2) + 32;  // exactly 3 leaves
  std::string bits(n, '1');

  // Put all "10" patterns into the last (partial) leaf.
  for (size_t i = leaf_block_bits * 2; i + 1 < n; i += 4) {
    bits[i] = '1';
    bits[i + 1] = '0';
  }

  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm(std::span<const std::uint64_t>(words), bits.size(),
                    leaf_block_bits);
  NaiveRmM nv(bits);

  const size_t pairs10 = nv.rank10(n);
  ASSERT_GT(pairs10, 0u);
  for (size_t k = 1; k <= pairs10 + 1; ++k) {
    EXPECT_EQ(rm.select10(k), nv.select10(k)) << "select10 k=" << k;
  }
}

/**
 * Invalid arguments should fail fast and return npos/0 as specified.
 * Covers bad ranks, bad ranges and out-of-bounds BP navigation calls.
 */
TEST(RmMEdgeCases, InvalidArgumentsGuards) {
  const size_t n = 600;
  std::string bits(n, '1');
  for (size_t i = 0; i < n; i += 5) {
    bits[i] = '0';
  }

  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm(std::span<const std::uint64_t>(words), bits.size(),
                    /*leaf_block_bits=*/256);

  EXPECT_EQ(rm.select1(0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.select0(0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.select10(0), pixie::RmMTree::npos);

  EXPECT_EQ(rm.fwdsearch(n, 0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.bwdsearch(0, 0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.bwdsearch(n + 1, 0), pixie::RmMTree::npos);

  EXPECT_EQ(rm.range_min_query_pos(10, 9), pixie::RmMTree::npos);
  EXPECT_EQ(rm.range_min_query_pos(0, n), pixie::RmMTree::npos);
  EXPECT_EQ(rm.range_max_query_pos(10, 9), pixie::RmMTree::npos);
  EXPECT_EQ(rm.range_max_query_pos(0, n), pixie::RmMTree::npos);
  EXPECT_EQ(rm.range_min_query_val(10, 9), 0);
  EXPECT_EQ(rm.range_max_query_val(10, 9), 0);
  EXPECT_EQ(rm.mincount(10, 9), 0);
  EXPECT_EQ(rm.minselect(10, 9, 1), pixie::RmMTree::npos);
  EXPECT_EQ(rm.minselect(0, n - 1, 0), pixie::RmMTree::npos);

  EXPECT_EQ(rm.close(n), pixie::RmMTree::npos);
  EXPECT_EQ(rm.open(0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.open(n + 1), pixie::RmMTree::npos);
  EXPECT_EQ(rm.enclose(0), pixie::RmMTree::npos);
  EXPECT_EQ(rm.enclose(n + 1), pixie::RmMTree::npos);
}

/**
 * bit_count is larger than the provided words buffer.
 * Verifies that non-owning construction rejects spans that are too short.
 */
TEST(RmMEdgeCases, SpanConstructorRejectsShortInputStorage) {
  std::vector<std::uint64_t> words = {0xAAAAAAAAAAAAAAAAull};
  const size_t bit_count = 300;

  EXPECT_THROW(
      (void)pixie::RmMTree(std::span<const std::uint64_t>(words), bit_count),
      std::invalid_argument);
}

/**
 * Same bitvector built through different configuration paths (auto vs explicit
 * leaf size and different overhead caps). Query
 * results must be identical.
 */
TEST(RmMEdgeCases, ExplicitBuildParametersAndOverheadCap) {
  std::mt19937_64 rng(42);
  const size_t n = 128;
  const std::string bits = random_bits(rng, n);
  NaiveRmM nv(bits);

  auto words = pack_words_lsb_first(bits);
  const std::span<const std::uint64_t> word_span(words);
  pixie::RmMTree rm_auto(word_span, n, 0, 1.f);
  pixie::RmMTree rm_explicit(word_span, n, 512, 2.f);
  pixie::RmMTree rm_span(word_span, n, 256, 1.f);

  expect_rank_select_equal(rm_auto, nv, n);
  expect_range_ops_equal(rm_auto, nv, n);
  expect_rank_select_equal(rm_explicit, nv, n);
  expect_range_ops_equal(rm_explicit, nv, n);
  expect_rank_select_equal(rm_span, nv, n);
  expect_range_ops_equal(rm_span, nv, n);
}

TEST(RmMEdgeCases, BoundaryHeavyQueries) {
  constexpr size_t leaf_block_bits = 128;
  const size_t n = 777;
  std::string bits(n, '0');
  for (size_t i = 0; i < n; ++i) {
    bits[i] = (((i * 37 + i / 7) % 11) < 5) ? '1' : '0';
  }
  for (size_t position : {63u, 64u, 65u, 127u, 128u, 129u, 255u, 256u}) {
    bits[position] = (position & 1u) ? '1' : '0';
  }

  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm(std::span<const std::uint64_t>(words), n, leaf_block_bits);
  NaiveRmM nv(bits);

  for (size_t end_position : std::array<size_t, 12>{0, 1, 63, 64, 65, 127, 128,
                                                    129, 255, 256, 512, n}) {
    EXPECT_EQ(rm.rank1(end_position), nv.rank1(end_position));
    EXPECT_EQ(rm.rank0(end_position), nv.rank0(end_position));
    EXPECT_EQ(rm.rank10(end_position), nv.rank10(end_position));
  }

  const size_t ones = nv.rank1(n);
  const size_t zeros = n - ones;
  const size_t pairs10 = nv.rank10(n);
  for (size_t rank : std::array<size_t, 7>{1, 2, 7, 31, 64, ones, ones + 1}) {
    EXPECT_EQ(rm.select1(rank), nv.select1(rank)) << "select1 rank=" << rank;
  }
  for (size_t rank : std::array<size_t, 7>{1, 2, 7, 31, 64, zeros, zeros + 1}) {
    EXPECT_EQ(rm.select0(rank), nv.select0(rank)) << "select0 rank=" << rank;
  }
  for (size_t rank :
       std::array<size_t, 7>{1, 2, 7, 31, 64, pairs10, pairs10 + 1}) {
    EXPECT_EQ(rm.select10(rank), nv.select10(rank)) << "select10 rank=" << rank;
  }

  for (auto [left, right] :
       std::array<std::pair<size_t, size_t>, 10>{{{0, 0},
                                                  {57, 75},
                                                  {63, 130},
                                                  {120, 260},
                                                  {128, 511},
                                                  {129, 511},
                                                  {250, 520},
                                                  {512, n - 1},
                                                  {0, n - 1},
                                                  {n - 9, n - 1}}}) {
    EXPECT_EQ(rm.range_min_query_pos(left, right),
              nv.range_min_query_pos(left, right));
    EXPECT_EQ(rm.range_min_query_val(left, right),
              nv.range_min_query_val(left, right));
    EXPECT_EQ(rm.range_max_query_pos(left, right),
              nv.range_max_query_pos(left, right));
    EXPECT_EQ(rm.range_max_query_val(left, right),
              nv.range_max_query_val(left, right));

    const size_t count = nv.mincount(left, right);
    EXPECT_EQ(rm.mincount(left, right), count);
    EXPECT_EQ(rm.minselect(left, right, 1), nv.minselect(left, right, 1));
    EXPECT_EQ(rm.minselect(left, right, count),
              nv.minselect(left, right, count));
    EXPECT_EQ(rm.minselect(left, right, count + 1),
              nv.minselect(left, right, count + 1));
  }

  for (size_t position : std::array<size_t, 12>{0, 1, 63, 64, 65, 127, 128, 129,
                                                255, 256, 512, n - 1}) {
    for (int delta : {-40, -9, -8, -2, -1, 0, 1, 2, 8, 9, 40}) {
      EXPECT_EQ(rm.fwdsearch(position, delta), nv.fwdsearch(position, delta))
          << "fwdsearch position=" << position << " delta=" << delta;
      EXPECT_EQ(rm.bwdsearch(position + 1, delta),
                nv.bwdsearch(position + 1, delta))
          << "bwdsearch position=" << (position + 1) << " delta=" << delta;
    }
  }
}

TEST_F(RmMRandomTest, LongRandom) {
  std::uniform_int_distribution<int> coin(0, 1);
  std::uniform_int_distribution<int> len_u(1, (int)kLongMaxBits);
  std::uniform_int_distribution<int> len_even(0, (int)(kLongMaxBits / 2));

  for (size_t iter = 0; iter < kRandomCases; ++iter) {
    std::string bits;
    if (coin(rng) == 0) {
      bits = random_bits(rng, (size_t)len_u(rng));
    } else {
      size_t m = 1 + (size_t)len_even(rng);
      bits = random_dyck_bits(rng, m);
      if (bits.empty()) {
        bits = "10";
      }
    }

    run_case_and_compare(bits, rng, kLongOpsPerCase, kSeed);
  }
}

TEST(RmMTest, RankBasic) {
  std::vector<std::uint64_t> bits = {0b10110};
  pixie::RmMTree rm(std::span<const std::uint64_t>(bits), 5);

  EXPECT_EQ(rm.rank1(0), 0);  // No bits
  EXPECT_EQ(rm.rank1(1), 0);  // 0
  EXPECT_EQ(rm.rank1(2), 1);  // 10
  EXPECT_EQ(rm.rank1(3), 2);  // 110
  EXPECT_EQ(rm.rank1(4), 2);  // 0110
  EXPECT_EQ(rm.rank1(5), 3);  // 10110
}

TEST(RmMTest, RankWithZeros) {
  std::vector<std::uint64_t> bits = {0};
  pixie::RmMTree rm(std::span<const std::uint64_t>(bits), 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(rm.rank1(i), 0);
  }
}

TEST(RmMTest, SelectBasic) {
  std::vector<std::uint64_t> bits = {0b1100010110010110};
  pixie::RmMTree rm(std::span<const std::uint64_t>(bits), 16);

  EXPECT_EQ(rm.select1(1), 1);
  EXPECT_EQ(rm.select1(2), 2);
  EXPECT_EQ(rm.select1(3), 4);
  EXPECT_EQ(rm.select1(4), 7);
  EXPECT_EQ(rm.select1(5), 8);
  EXPECT_EQ(rm.select1(6), 10);
  EXPECT_EQ(rm.select1(7), 14);
  EXPECT_EQ(rm.select1(8), 15);
}

TEST(RmMTest, MainRankTest) {
  std::mt19937_64 rng(42);
  std::vector<std::uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  size_t rm_size = 65536 * 32 * 64;
  pixie::RmMTree rm(std::span<const std::uint64_t>(bits), rm_size);
  size_t rank = 0;
  for (size_t i = 0; i < rm_size; ++i) {
    ASSERT_EQ(rank, rm.rank1(i));
    rank += (bits[i >> 6] >> (i & 63)) & 1u;
  }
}

TEST(RmMTest, MainSelectTest) {
  std::mt19937_64 rng(42);
  std::vector<std::uint64_t> bits(65536 * 32);
  for (size_t i = 0; i < 65536 * 32; i++) {
    bits[i] = rng();
  }

  size_t rm_size = 65536 * 32 * 64;
  pixie::RmMTree rm(std::span<const std::uint64_t>(bits), rm_size);
  size_t rank = 0;

  for (size_t i = 0; i < rm_size; ++i) {
    if ((bits[i >> 6] >> (i & 63)) & 1u) {
      ASSERT_EQ(rm.select1(++rank), i);
      ASSERT_EQ(rm.rank1(i), rank - 1);
      ASSERT_EQ(rm.rank1(i + 1), rank);
    }
  }
}
