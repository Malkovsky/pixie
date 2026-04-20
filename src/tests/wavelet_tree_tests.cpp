#include <gtest/gtest.h>
#include <pixie/utils.h>
#include <pixie/wavelet_tree.h>

#include <random>

using pixie::WaveletTree;

TEST(WaveletTreeTest, BasicSelect) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  std::vector<std::vector<size_t>> rank(alphabet_size);
  for (size_t i = 0; i < data_size; i++) {
    rank[data[i]].push_back(i);
  }

  WaveletTree wavelet_tree(alphabet_size, data);

  for (uint64_t symb = 0; symb < alphabet_size; symb++) {
    for (size_t i = 0; i <= rank[symb].size(); i++) {
      uint64_t exp = i == rank[symb].size() ? data_size : rank[symb][i];
      uint64_t act = wavelet_tree.select(symb, i + 1);
      EXPECT_EQ(act, exp);
    }
  }
}

TEST(WaveletTreeTest, BasicRank) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  std::vector<size_t> count(alphabet_size);

  WaveletTree wavelet_tree(alphabet_size, data);
  for (size_t i = 0; i <= data_size; i++) {
    for (uint64_t symb = 0; symb < alphabet_size; symb++) {
      uint64_t exp = count[symb];
      uint64_t act = wavelet_tree.rank(symb, i);
      EXPECT_EQ(act, exp);
    }

    if (i == data_size) {
      break;
    }
    count[data[i]]++;
  }
}

TEST(WaveletTreeTest, BasicSegment) {
  const std::vector<uint64_t> data = {3, 2, 0, 3, 1, 1, 2};
  size_t data_size = 7, alphabet_size = 4;

  WaveletTree wavelet_tree(alphabet_size, data);

  for (size_t begin = 0; begin <= data_size; begin++) {
    for (size_t end = begin; end <= data_size; end++) {
      auto segment = wavelet_tree.getSegment(begin, end);
      EXPECT_EQ(segment.size(), end - begin);
      for(size_t i = 0; i < end - begin; i++){
        EXPECT_EQ(segment[i], data[begin + i]);
      }
    }
  }
}

TEST(WaveletTreeTest, SmokeSelect) {
  std::vector<std::vector<size_t>> rank;
  for (size_t data_size = 8; data_size < (1 << 22); data_size <<= 1) {
    size_t alphabet_size = 1024;
    std::mt19937_64 rng(239);
    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);

    rank.assign(alphabet_size, {});
    for (size_t i = 0; i < data_size; i++) {
      rank[data[i]].push_back(i);
    }

    WaveletTree wavelet_tree(alphabet_size, data);

    for (uint64_t symb = 0; symb < alphabet_size; symb++) {
      for (size_t i = 0; i <= rank[symb].size(); i++) {
        uint64_t exp = i == rank[symb].size() ? data_size : rank[symb][i];
        uint64_t act = wavelet_tree.select(symb, i + 1);
        EXPECT_EQ(act, exp);
      }
    }
  }
}

TEST(WaveletTreeTest, SmokeRank) {
  std::vector<size_t> count;
  for (size_t data_size = 8; data_size < (1 << 22); data_size <<= 1) {
    size_t alphabet_size = 1024;
    std::mt19937_64 rng(239);
    std::vector<uint64_t> data =
        generate_random_data(data_size, alphabet_size, rng);
    std::vector<uint64_t> query =
        generate_random_data(data_size + 1, alphabet_size, rng);

     count.assign(alphabet_size, 0);

    WaveletTree wavelet_tree(alphabet_size, data);

    for (size_t i = 0; i <= data_size; i++) {
      uint64_t symb = query[i];
      uint64_t exp = count[symb];
      uint64_t act = wavelet_tree.rank(symb, i);
      EXPECT_EQ(act, exp);

      if (i == data_size) {
        break;
      }
      count[data[i]]++;
    }
  }
}


TEST(WaveletTreeTest, SmokeSegment) {
  size_t data_size = 256, alphabet_size = 100;

  std::mt19937_64 rng(239);
  std::vector<uint64_t> data =
      generate_random_data(data_size, alphabet_size, rng);

  WaveletTree wavelet_tree(alphabet_size, data);

  for (size_t begin = 0; begin <= data_size; begin++) {
    for (size_t end = begin; end <= data_size; end++) {
      auto segment = wavelet_tree.getSegment(begin, end);
      EXPECT_EQ(segment.size(), end - begin);
      for(size_t i = 0; i < end - begin; i++){
        EXPECT_EQ(segment[i], data[begin + i]);
      }
    }
  }
}