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

    auto mmap_tree =
        pixie::WaveletTreeBase<pixie::MmapViewStorage>::deserialize(byte_span);

    for (size_t i = 0; i <= data_size; i += 16) {
      uint64_t symb = data[i == data_size ? 0 : i];
      EXPECT_EQ(orig_tree.rank(symb, i), mmap_tree.rank(symb, i));
    }

    std::vector<size_t> count(alphabet_size, 0);
    for (auto symb : data) {
      count[symb]++;
    }

    for (uint64_t symb = 0; symb < alphabet_size; symb++) {
      for (uint64_t rank = 1; rank <= count[symb]; rank++) {
        EXPECT_EQ(orig_tree.select(symb, rank), mmap_tree.select(symb, rank));
      }
    }

    auto orig_segment = orig_tree.get_segment(0, data_size);
    auto mmap_segment = mmap_tree.get_segment(0, data_size);
    EXPECT_EQ(orig_segment, mmap_segment);
  }
}
