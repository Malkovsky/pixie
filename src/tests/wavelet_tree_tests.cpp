#include <gtest/gtest.h>
#include <pixie/utils.h>
#include <pixie/wavelet_tree/implementations.h>

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
      auto segment = wavelet_tree.get_segment(begin, end);
      EXPECT_EQ(segment.size(), end - begin);
      for (size_t i = 0; i < end - begin; i++) {
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

    for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                            pixie::WaveletTreeBuildType::Huffman}) {
      WaveletTree wavelet_tree(alphabet_size, data, build_type);

      for (uint64_t symb = 0; symb < alphabet_size; symb++) {
        for (size_t i = 0; i <= rank[symb].size(); i++) {
          uint64_t exp = i == rank[symb].size() ? data_size : rank[symb][i];
          uint64_t act = wavelet_tree.select(symb, i + 1);
          EXPECT_EQ(act, exp);
        }
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

    for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                            pixie::WaveletTreeBuildType::Huffman}) {
      count.assign(alphabet_size, 0);
      WaveletTree wavelet_tree(alphabet_size, data, build_type);

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
}

TEST(WaveletTreeTest, SmokeSegment) {
  size_t data_size = 256, alphabet_size = 100;

  std::mt19937_64 rng(239);
  std::vector<uint64_t> data =
      generate_random_data(data_size, alphabet_size, rng);

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree wavelet_tree(alphabet_size, data, build_type);

    for (size_t begin = 0; begin <= data_size; begin++) {
      for (size_t end = begin; end <= data_size; end++) {
        auto segment = wavelet_tree.get_segment(begin, end);
        EXPECT_EQ(segment.size(), end - begin);
        for (size_t i = 0; i < end - begin; i++) {
          EXPECT_EQ(segment[i], data[begin + i]);
        }
      }
    }
  }
}

TEST(WaveletTreeTest, AlphabetSizeOne) {
  // alphabet_size == 1: tree has no internal nodes (root_ == npos).
  // This previously crashed in get_segment_impl.
  const std::vector<uint64_t> data = {0, 0, 0, 0, 0};
  size_t data_size = 5, alphabet_size = 1;

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree wavelet_tree(alphabet_size, data, build_type);

    // rank: all positions have symbol 0
    for (size_t pos = 0; pos <= data_size; pos++) {
      EXPECT_EQ(wavelet_tree.rank(0, pos), pos);
    }
    // rank of out-of-range symbol returns 0
    EXPECT_EQ(wavelet_tree.rank(1, 3), 0);

    // select: the k-th occurrence of symbol 0 is at position k-1
    for (size_t rank = 1; rank <= data_size; rank++) {
      EXPECT_EQ(wavelet_tree.select(0, rank), rank - 1);
    }
    // out-of-range rank returns data_size
    EXPECT_EQ(wavelet_tree.select(0, data_size + 1), data_size);

    // get_segment: all elements are 0
    for (size_t begin = 0; begin <= data_size; begin++) {
      for (size_t end = begin; end <= data_size; end++) {
        auto segment = wavelet_tree.get_segment(begin, end);
        EXPECT_EQ(segment.size(), end - begin);
        for (size_t i = 0; i < end - begin; i++) {
          EXPECT_EQ(segment[i], 0);
        }
      }
    }
  }
}

TEST(WaveletTreeTest, AlphabetSizeZero) {
  size_t alphabet_size = 0;
  const std::vector<uint64_t> data = {};
  WaveletTree wavelet_tree(alphabet_size, data);

  EXPECT_EQ(wavelet_tree.size(), 0);
  EXPECT_TRUE(wavelet_tree.empty());
  EXPECT_EQ(wavelet_tree.rank(0, 0), 0);
  EXPECT_EQ(wavelet_tree.select(0, 1), 0);
  EXPECT_EQ(wavelet_tree.get_segment(0, 0).size(), 0);
}

TEST(WaveletTreeTest, EmptyData) {
  size_t alphabet_size = 4;
  const std::vector<uint64_t> data = {};
  WaveletTree wavelet_tree(alphabet_size, data);

  EXPECT_EQ(wavelet_tree.size(), 0);
  EXPECT_TRUE(wavelet_tree.empty());
  for (uint64_t symb = 0; symb < alphabet_size; symb++) {
    EXPECT_EQ(wavelet_tree.rank(symb, 0), 0);
    EXPECT_EQ(wavelet_tree.select(symb, 1), 0);
  }
}

TEST(WaveletTreeTest, SymbolWithZeroOccurrences) {
  // alphabet_size = 5, but symbol 4 never appears in the data.
  const std::vector<uint64_t> data = {0, 1, 2, 3, 0, 1, 2, 3};
  size_t data_size = 8, alphabet_size = 5;

  WaveletTree wavelet_tree(alphabet_size, data);

  // rank of absent symbol is 0 at every position
  for (size_t pos = 0; pos <= data_size; pos++) {
    EXPECT_EQ(wavelet_tree.rank(4, pos), 0);
  }

  // select of absent symbol returns data_size
  EXPECT_EQ(wavelet_tree.select(4, 1), data_size);
  EXPECT_EQ(wavelet_tree.select(4, 2), data_size);

  // get_segment still works correctly
  for (size_t begin = 0; begin <= data_size; begin++) {
    for (size_t end = begin; end <= data_size; end++) {
      auto segment = wavelet_tree.get_segment(begin, end);
      EXPECT_EQ(segment.size(), end - begin);
      for (size_t i = 0; i < end - begin; i++) {
        EXPECT_EQ(segment[i], data[begin + i]);
      }
    }
  }
}

TEST(WaveletTreeTest, AllSameSymbolLargeAlphabet) {
  // All data is symbol 0, but alphabet_size is large.
  // Most symbols have 0 occurrences.
  const std::vector<uint64_t> data = {0, 0, 0, 0, 0, 0, 0, 0};
  size_t data_size = 8, alphabet_size = 256;

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree wavelet_tree(alphabet_size, data, build_type);

    EXPECT_EQ(wavelet_tree.rank(0, data_size), data_size);
    for (uint64_t symb = 1; symb < alphabet_size; symb++) {
      EXPECT_EQ(wavelet_tree.rank(symb, data_size), 0);
      EXPECT_EQ(wavelet_tree.select(symb, 1), data_size);
    }

    auto segment = wavelet_tree.get_segment(0, data_size);
    for (auto s : segment) {
      EXPECT_EQ(s, 0);
    }
  }
}

TEST(WaveletTreeTest, SerializationSmoke) {
  size_t data_size = 4096, alphabet_size = 100;

  std::mt19937_64 rng(239);
  std::vector<uint64_t> data =
      generate_random_data(data_size, alphabet_size, rng);

  for (auto build_type : {pixie::WaveletTreeBuildType::Standard,
                          pixie::WaveletTreeBuildType::Huffman}) {
    WaveletTree orig_tree(alphabet_size, data, build_type);

    pixie::OutputBitStream bs;
    orig_tree.serialize(bs);
    std::vector<uint64_t> serialized_data = bs.extract();

    std::span<const std::byte> byte_span(
        reinterpret_cast<const std::byte*>(serialized_data.data()),
        serialized_data.size() * sizeof(uint64_t));

    auto view_tree = pixie::WaveletTreeView::deserialize(byte_span);

    for (size_t i = 0; i <= data_size; i += 16) {
      uint64_t symb = data[i == data_size ? 0 : i];
      EXPECT_EQ(orig_tree.rank(symb, i), view_tree.rank(symb, i));
    }

    std::vector<size_t> count(alphabet_size, 0);
    for (auto symb : data) {
      count[symb]++;
    }

    for (uint64_t symb = 0; symb < alphabet_size; symb++) {
      for (uint64_t rank = 1; rank <= count[symb]; rank++) {
        EXPECT_EQ(orig_tree.select(symb, rank), view_tree.select(symb, rank));
      }
    }

    auto orig_segment = orig_tree.get_segment(0, data_size);
    auto view_segment = view_tree.get_segment(0, data_size);
    EXPECT_EQ(orig_segment, view_segment);
  }
}
