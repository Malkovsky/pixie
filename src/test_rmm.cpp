#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "naive_rmm_tree.h"
#include "rmm_tree.h"

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
  std::string s(2 * m, '0');
  if (m == 0) {
    return s;
  }
  size_t opens_left = m;
  size_t closes_left = m;
  int h = 0;
  std::bernoulli_distribution coin(0.5);
  for (size_t pos = 0; pos < 2 * m; ++pos) {
    if (opens_left == 0) {
      s[pos] = '0';
      --closes_left;
      --h;
      continue;
    }
    if (closes_left == 0) {
      s[pos] = '1';
      --opens_left;
      ++h;
      continue;
    }
    if (h == 0) {
      s[pos] = '1';
      --opens_left;
      ++h;
      continue;
    }
    if ((size_t)h == closes_left) {
      s[pos] = '0';
      --closes_left;
      --h;
      continue;
    }
    if (coin(rng)) {
      s[pos] = '1';
      --opens_left;
      ++h;
    } else {
      s[pos] = '0';
      --closes_left;
      --h;
    }
  }
  return s;
}

struct Limits {
  size_t CASES = 150;
  size_t OPS_PER_CASE = 600;
  size_t MAX_N = 20000;
};

static Limits load_limits_from_env() {
  Limits L;
  if (const char* s = std::getenv("RMM_CASES")) {
    L.CASES = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }
  if (const char* s = std::getenv("RMM_OPS")) {
    L.OPS_PER_CASE = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }
  if (const char* s = std::getenv("RMM_MAX_N")) {
    L.MAX_N = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }
  return L;
}

static uint64_t choose_seed() {
  if (const char* s = std::getenv("RMM_SEED")) {
    return std::strtoull(s, nullptr, 10);
  }
  std::random_device rd;
  return ((uint64_t)rd() << 32) ^ (uint64_t)rd() ^
         (uint64_t)std::chrono::high_resolution_clock::now()
             .time_since_epoch()
             .count();
}

static void run_case_and_compare(const std::string& bits,
                                 std::mt19937_64& rng,
                                 size_t ops_per_case,
                                 uint64_t seed) {
  const bool use_words = std::uniform_int_distribution<int>(0, 1)(rng);
  pixie::RmMTree rm;
  NaiveRmM nv;
  if (use_words) {
    auto words = pack_words_lsb_first(bits);
    rm = pixie::RmMTree(words, bits.size());
    nv = NaiveRmM(words, bits.size());
  } else {
    rm = pixie::RmMTree(bits);
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
  uint64_t seed{};
  std::mt19937_64 rng;
  Limits L;
  void SetUp() override {
    L = load_limits_from_env();
    seed = choose_seed();
    rng.seed(seed);
    std::cerr << "[ RmMRandomTest ] seed=" << seed << " CASES=" << L.CASES
              << " OPS=" << L.OPS_PER_CASE << " MAX_N=" << L.MAX_N << "\n";
  }
};

TEST_F(RmMRandomTest, RandomMixedBits) {
  std::uniform_int_distribution<int> len_u(1, (int)L.MAX_N);
  for (size_t t = 0; t < L.CASES; ++t) {
    const size_t n = (size_t)len_u(rng);
    const std::string bits = random_bits(rng, n);
    run_case_and_compare(bits, rng, L.OPS_PER_CASE, seed);
  }
}

TEST_F(RmMRandomTest, RandomDyckBits) {
  std::uniform_int_distribution<int> len_even(0, (int)(L.MAX_N / 2));
  for (size_t t = 0; t < L.CASES; ++t) {
    size_t m = 1 + (size_t)len_even(rng);
    std::string bits = random_dyck_bits(rng, m);
    if (bits.empty()) {
      bits = "10";
    }
    run_case_and_compare(bits, rng, L.OPS_PER_CASE, seed);
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
        run_case_and_compare(bits, rng, /*ops_per_case=*/200,
                             /*seed=*/0xC0FFEEull);
      }
      {
        SCOPED_TRACE(::testing::Message()
                     << "short-inputs:words  bits=" << bits);
        run_case_and_compare(bits, rng, /*ops_per_case=*/200,
                             /*seed=*/0xC0FFEEull);
      }
    }
  }
}

static void run_case_string_vs_words(const std::string& bits,
                                     std::mt19937_64& rng) {
  const size_t N = bits.size();
  pixie::RmMTree rm_s(bits);
  auto words = pack_words_lsb_first(bits);
  pixie::RmMTree rm_w(words, N);

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

TEST_F(RmMRandomTest, WordsVsString) {
  std::uniform_int_distribution<int> len_u(0, (int)L.MAX_N);
  for (size_t t = 0; t < L.CASES; ++t) {
    const size_t n = (size_t)len_u(rng);
    const std::string bits = random_bits(rng, n);
    run_case_string_vs_words(bits, rng);
  }
}

TEST(RmMEdgeCases, EmptyInput) {
  pixie::RmMTree rm(std::string{});
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

TEST(RmMTreeStress, LongRandom) {
  Limits L;
  L.OPS_PER_CASE = 2000;
  L.MAX_N = 65536;

  if (const char* s = std::getenv("RMM_CASES")) {
    L.CASES = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }
  if (const char* s = std::getenv("RMM_OPS")) {
    L.OPS_PER_CASE = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }
  if (const char* s = std::getenv("RMM_MAX_N")) {
    L.MAX_N = std::max<size_t>(1, std::strtoull(s, nullptr, 10));
  }

  const uint64_t seed = choose_seed();
  std::mt19937_64 rng(seed);
  size_t LOG_EVERY = 10;
  if (const char* s = std::getenv("RMM_LOG_EVERY")) {
    size_t v = std::strtoull(s, nullptr, 10);
    if (v) {
      LOG_EVERY = v;
    }
  }

  std::cerr << "[ LongRandom ] seed=" << seed << " CASES=" << L.CASES
            << " OPS=" << L.OPS_PER_CASE << " MAX_N=" << L.MAX_N
            << " LOG_EVERY=" << LOG_EVERY << "\n";

  std::uniform_int_distribution<int> coin(0, 1);
  std::uniform_int_distribution<int> len_u(1, (int)L.MAX_N);
  std::uniform_int_distribution<int> len_even(0, (int)(L.MAX_N / 2));

  size_t total_ops = 0;

  for (size_t iter = 1; iter <= L.CASES; ++iter) {
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

    run_case_and_compare(bits, rng, L.OPS_PER_CASE, seed);
    total_ops += L.OPS_PER_CASE;

    if (iter % LOG_EVERY == 0) {
      std::cerr << "[ LongRandom ] iter=" << iter << " total_ops=" << total_ops
                << " last_N=" << bits.size() << " ok\n";
    }
  }
}

TEST(RmMTest, RankBasic) {
  std::vector<std::uint64_t> bits = {0b10110};
  pixie::RmMTree rm(bits, 5);

  EXPECT_EQ(rm.rank1(0), 0);  // No bits
  EXPECT_EQ(rm.rank1(1), 0);  // 0
  EXPECT_EQ(rm.rank1(2), 1);  // 10
  EXPECT_EQ(rm.rank1(3), 2);  // 110
  EXPECT_EQ(rm.rank1(4), 2);  // 0110
  EXPECT_EQ(rm.rank1(5), 3);  // 10110
}

TEST(RmMTest, RankWithZeros) {
  std::vector<std::uint64_t> bits = {0};
  pixie::RmMTree rm(bits, 5);

  for (size_t i = 0; i <= 5; i++) {
    EXPECT_EQ(rm.rank1(i), 0);
  }
}

TEST(RmMTest, SelectBasic) {
  std::vector<std::uint64_t> bits = {0b1100010110010110};
  pixie::RmMTree rm(bits, 16);

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
  pixie::RmMTree rm(bits, rm_size);
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
  pixie::RmMTree rm(bits, rm_size);
  size_t rank = 0;

  for (size_t i = 0; i < rm_size; ++i) {
    if ((bits[i >> 6] >> (i & 63)) & 1u) {
      ASSERT_EQ(rm.select1(++rank), i);
      ASSERT_EQ(rm.rank1(i), rank - 1);
      ASSERT_EQ(rm.rank1(i + 1), rank);
    }
  }
}
