cd #include<gtest / gtest.h>
#include <pixie/detail/sequence/bit_block.h>
#include <pixie/detail/sequence/packed_bit_block.h>
#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/experimental/grouped_slot_order256.h>
#include <pixie/experimental/permuted_bit_block.h>
#include <pixie/experimental/power_of_two_value_block.h>
#include <pixie/experimental/slot_order.h>
#include <pixie/experimental/slot_order16.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <ranges>
#include <span>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

    namespace {
  using namespace pixie;
  using namespace pixie::detail::sequence;
  using namespace pixie::experimental;

  template <class T>
  class WideSlotOrderSpec : public ::testing::Test {};
  using WideSlotOrders =
      ::testing::Types<SplitSlotOrder<32>, SplitSlotOrder<64>,
                       ByteSlotOrder<32>, ByteSlotOrder<64>,
                       ByteSlotOrder<32, false>, ByteSlotOrder<64, false>>;
  TYPED_TEST_SUITE(WideSlotOrderSpec, WideSlotOrders);

  TYPED_TEST(WideSlotOrderSpec, EveryFullNodeRangeAndDistance) {
    using Order = TypeParam;
    constexpr auto n = Order::capacity;
    Order initial(n);
    initial.rotate_left(0, n, 19);
    for (std::size_t left = 0; left <= n; ++left) {
      for (std::size_t right = left; right <= n; ++right) {
        const auto length = right - left;
        for (std::size_t distance = 0; distance <= 2 * length + 1; ++distance) {
          auto order = initial;
          order.rotate_left(left, right, distance);
          ASSERT_TRUE(order.valid());
          for (std::size_t i = 0; i < n; ++i) {
            const auto source = i >= left && i < right
                                    ? left + (i - left + distance) % length
                                    : i;
            ASSERT_EQ(order[i], initial[source])
                << left << ':' << right << ':' << distance << ':' << i;
          }
        }
      }
    }
  }

  TYPED_TEST(WideSlotOrderSpec,
             PartialNodesRepeatedRotationsKeepPhysicalSlots) {
    using Order = TypeParam;
    constexpr auto capacity = Order::capacity;
    std::array<unsigned, capacity> values{};
    std::mt19937 random(997);
    for (std::size_t n = 0; n <= capacity; ++n) {
      Order order(n);
      std::array<void*, capacity> pointers{};
      std::vector<void*> expected;
      for (std::size_t i = 0; i < n; ++i) {
        pointers[i] = &values[i];
        expected.push_back(pointers[i]);
      }
      const auto original = pointers;
      for (unsigned step = 0; step < 256; ++step) {
        const auto left = random() % (n + 1);
        const auto right = left + random() % (n - left + 1);
        const auto distance = random();
        ASSERT_EQ(order.rotate_children(pointers, left, right, distance), 0u);
        if (left != right) {
          std::rotate(expected.begin() + left,
                      expected.begin() + left + distance % (right - left),
                      expected.begin() + right);
        }
        ASSERT_TRUE(order.valid());
        ASSERT_EQ(pointers, original);
        for (std::size_t i = 0; i < n; ++i) {
          ASSERT_EQ(pointers[order[i]], expected[i]);
        }
        for (std::size_t i = 0; i < capacity; ++i) {
          ASSERT_EQ(order.occupied(i), i < n);
        }
      }
    }
  }

  TEST(GroupedSlotOrder256, FullGroupSpillsAndWholeGroupReordering) {
    std::array<unsigned, 256> payload{};
    for (std::size_t n : {16u, 32u, 48u, 240u, 256u}) {
      for (std::size_t distance = 0; distance <= 2 * n; ++distance) {
        GroupedSlotOrder256 order(n);
        std::array<void*, 256> pointers{};
        for (std::size_t i = 0; i < n; ++i) {
          pointers[i] = &payload[i];
        }
        const auto writes = order.rotate_children(pointers, 0, n, distance);
        const auto remainder = (distance % n) % 16;
        EXPECT_EQ(writes, n == 16
                              ? 0u
                              : (n / 16) * std::min(remainder, 16 - remainder));
        ASSERT_TRUE(order.valid());
        for (std::size_t i = 0; i < n; ++i) {
          ASSERT_EQ(pointers[order[i]], &payload[(i + distance) % n]);
        }
      }
    }
  }

  TEST(GroupedSlotOrder256, PartialGroupsAndRepeatedArbitraryRotations) {
    std::array<unsigned, 256> payload{};
    std::mt19937 random(211);
    for (std::size_t n = 0; n <= 256; ++n) {
      GroupedSlotOrder256 order(n);
      std::array<void*, 256> pointers{};
      std::vector<void*> expected;
      for (std::size_t i = 0; i < n; ++i) {
        pointers[i] = &payload[i];
        expected.push_back(pointers[i]);
      }
      for (unsigned step = 0; step < 256; ++step) {
        const auto left = random() % (n + 1);
        const auto right = left + random() % (n - left + 1);
        const auto distance = random();
        order.rotate_children(pointers, left, right, distance);
        if (left != right) {
          std::rotate(expected.begin() + left,
                      expected.begin() + left + distance % (right - left),
                      expected.begin() + right);
        }
        ASSERT_TRUE(order.valid());
        for (std::size_t i = 0; i < n; ++i) {
          ASSERT_EQ(pointers[order[i]], expected[i]);
        }
        for (std::size_t i = 0; i < 256; ++i) {
          ASSERT_EQ(order.occupied(i), i < n);
          if (!order.occupied(i)) {
            ASSERT_EQ(pointers[i], nullptr);
          }
        }
      }
    }
  }

  TEST(SlotOrder16, EveryRotationRangeAndDistance) {
    for (std::size_t n = 0; n <= 16; ++n) {
      SlotOrder16 initial(n);
      // A non-identity starting order also detects accidental slot renumbering.
      initial.rotate_left(0, n, 7);
      for (std::size_t left = 0; left <= n; ++left) {
        for (std::size_t right = left; right <= n; ++right) {
          for (std::size_t distance = 0; distance <= 2 * n + 1; ++distance) {
            auto actual = initial;
            std::vector<unsigned> expected;
            for (std::size_t i = 0; i < n; ++i) {
              expected.push_back(initial[i]);
            }
            if (left != right) {
              std::rotate(expected.begin() + left,
                          expected.begin() + left + distance % (right - left),
                          expected.begin() + right);
            }
            actual.rotate_left(left, right, distance);
            ASSERT_EQ(actual.size(), n);
            ASSERT_EQ(actual.occupied(), initial.occupied());
            for (std::size_t i = 0; i < n; ++i) {
              ASSERT_EQ(actual[i], expected[i]);
            }
          }
        }
      }
    }
  }

  TEST(SlotOrder16, ReusesHolesWithoutMovingSurvivingSlots) {
    SlotOrder16 actual;
    std::vector<unsigned> expected;
    std::mt19937 random(917);
    for (unsigned step = 0; step < 20000; ++step) {
      const auto n = expected.size();
      const auto operation = random() % 3;
      if (n == 0 || (n < 16 && operation == 0)) {
        const auto position = random() % (n + 1);
        const auto old_mask = actual.occupied();
        const auto slot = actual.insert(position);
        ASSERT_LT(slot, 16u);
        ASSERT_EQ(old_mask & (1u << slot), 0u);
        expected.insert(expected.begin() + position, slot);
      } else if (operation == 1) {
        const auto position = random() % n;
        ASSERT_EQ(actual.erase(position), expected[position]);
        expected.erase(expected.begin() + position);
      } else {
        const auto left = random() % (n + 1);
        const auto right = left + random() % (n - left + 1);
        const auto distance = random();
        actual.rotate_left(left, right, distance);
        if (left != right) {
          std::rotate(expected.begin() + left,
                      expected.begin() + left + distance % (right - left),
                      expected.begin() + right);
        }
      }
      ASSERT_EQ(actual.size(), expected.size());
      unsigned mask = 0;
      for (std::size_t i = 0; i < expected.size(); ++i) {
        ASSERT_EQ(actual[i], expected[i]);
        ASSERT_EQ(mask & (1u << actual[i]), 0u);
        mask |= 1u << actual[i];
      }
      ASSERT_EQ(actual.occupied(), mask);
    }
  }

  // No bit operations or exposed writable references: verifies the generic
  // contract independently of packed payload assumptions.
  struct IntegerBlock {
    using value_type = std::uint32_t;
    static constexpr std::size_t capacity = 8;
    std::array<value_type, capacity> values{};
    std::size_t n = 0;
    inline static std::size_t copied_elements = 0;
    IntegerBlock() noexcept = default;
    explicit IntegerBlock(std::span<const value_type> input) : n(input.size()) {
      assert(n <= capacity);
      std::copy(input.begin(), input.end(), values.begin());
    }
    std::size_t size() const noexcept { return n; }
    value_type operator[](std::size_t i) const {
      assert(i < n);
      return values[i];
    }
    void rotate_left(std::size_t l, std::size_t r, std::size_t d) noexcept {
      assert(l <= r && r <= n);
      if (l != r) {
        std::rotate(values.begin() + l, values.begin() + l + d % (r - l),
                    values.begin() + r);
      }
    }
    void redistribute(IntegerBlock& rhs, std::size_t left) noexcept {
      assert(this != &rhs && left <= capacity && left <= n + rhs.n &&
             n + rhs.n - left <= capacity);
      std::array<value_type, 2 * capacity> scratch;
      const auto total = n + rhs.n;
      std::copy_n(values.begin(), n, scratch.begin());
      std::copy_n(rhs.values.begin(), rhs.n, scratch.begin() + n);
      std::copy_n(scratch.begin(), left, values.begin());
      std::copy_n(scratch.begin() + left, total - left, rhs.values.begin());
      n = left;
      rhs.n = total - left;
      copied_elements += total;
    }
  };
  static_assert(SequenceBlock<IntegerBlock>);
  static_assert(SequenceBlock<BitBlock<>>);
  static_assert(SequenceBlock<PackedBitBlock<>>);
  static_assert(SequenceBlock<PermutedBitBlock<>>);

  struct MoveOnlyIntegerBlock : IntegerBlock {
    using IntegerBlock::IntegerBlock;
    MoveOnlyIntegerBlock() noexcept = default;
    MoveOnlyIntegerBlock(const MoveOnlyIntegerBlock&) = delete;
    MoveOnlyIntegerBlock& operator=(const MoveOnlyIntegerBlock&) = delete;
    MoveOnlyIntegerBlock(MoveOnlyIntegerBlock&&) noexcept = default;
    MoveOnlyIntegerBlock& operator=(MoveOnlyIntegerBlock&&) noexcept = default;
  };
  static_assert(SequenceBlock<MoveOnlyIntegerBlock>);

  struct OddIntegerBlock : IntegerBlock {
    using IntegerBlock::IntegerBlock;
    static constexpr std::size_t capacity = 7;
  };
  static_assert(SequenceBlock<OddIntegerBlock>);

  struct alignas(128) ConstOnlyAlignedBlock : IntegerBlock {
    using IntegerBlock::IntegerBlock;
    using IntegerBlock::operator[];
    value_type operator[](std::size_t) = delete;
  };
  static_assert(SequenceBlock<ConstOnlyAlignedBlock>);

  template <class Value>
  void rotate(std::vector<Value> & values, std::size_t l, std::size_t r,
              std::size_t d) {
    if (l != r) {
      std::rotate(values.begin() + l, values.begin() + l + d % (r - l),
                  values.begin() + r);
    }
  }

  std::vector<std::size_t> representative_origins(std::size_t size) {
    assert(size != 0);
    std::vector<std::size_t> origins = {
        0,
        1,
        2,
        63,
        64,
        65,
        127,
        128,
        129,
        size / 2,
        size > 1 ? size - 2 : 0,
        size - 1,
    };
    for (const auto numerator : {3U, 5U, 7U, 9U, 11U, 13U, 15U}) {
      origins.push_back(size * numerator / 16);
    }
    std::erase_if(origins,
                  [size](std::size_t origin) { return origin >= size; });
    std::ranges::sort(origins);
    origins.erase(std::unique(origins.begin(), origins.end()), origins.end());
    return origins;
  }

  std::vector<std::uint64_t> pack(const std::vector<bool>& bits) {
    std::vector<std::uint64_t> words((bits.size() + 63) / 64);
    for (std::size_t i = 0; i < bits.size(); ++i) {
      words[i / 64] |= std::uint64_t{bits[i]} << (i % 64);
    }
    if (bits.size() % 64) {
      words.back() |= ~std::uint64_t{0} << (bits.size() % 64);
    }
    return words;
  }
  template <class Tree>
  std::vector<typename Tree::value_type> data(std::size_t n,
                                              std::size_t seed = 42) {
    std::mt19937_64 random(seed);
    std::vector<typename Tree::value_type> result(n);
    for (std::size_t i = 0; i < n; ++i) {
      if constexpr (std::same_as<typename Tree::value_type, bool>) {
        result[i] = random() & 1;
      } else {
        result[i] = static_cast<std::uint32_t>(random());
      }
    }
    return result;
  }
  template <class Tree>
  Tree build(const std::vector<typename Tree::value_type>& values,
             std::size_t chunk = Tree::block_capacity) {
    using Block = typename Tree::block_type;
    // A lazy, single-pass-compatible producer avoids an auxiliary block buffer.
    auto blocks =
        std::views::iota(std::size_t{0}, (values.size() + chunk - 1) / chunk) |
        std::views::transform([&](std::size_t index) {
          const auto begin = index * chunk;
          const auto n = std::min(chunk, values.size() - begin);
          if constexpr (std::same_as<typename Tree::value_type, bool>) {
            std::vector<bool> local(values.begin() + begin,
                                    values.begin() + begin + n);
            return Block(pack(local), n);
          } else {
            return Block(std::span(values).subspan(begin, n));
          }
        });
    return Tree::from_blocks(blocks);
  }
  template <class Tree>
  void check(const Tree& tree,
             const std::vector<typename Tree::value_type>& values) {
    ASSERT_TRUE(tree.test_validate());
    ASSERT_EQ(tree.size(), values.size());
    ASSERT_EQ(tree.empty(), values.empty());
    for (std::size_t i = 0; i < values.size(); ++i) {
      ASSERT_EQ(tree[i], values[i]) << "index=" << i;
    }
    std::size_t offset = 0;
    tree.for_each_block([&](const auto& block) {
      for (std::size_t i = 0; i < block.size(); ++i) {
        ASSERT_EQ(block[i], values[offset++]);
      }
    });
    EXPECT_EQ(offset, values.size());
    const auto memory = tree.memory_usage();
    EXPECT_EQ(memory.total_bytes,
              sizeof(tree) + memory.block_bytes + memory.node_bytes);
    EXPECT_EQ(memory.block_bytes, memory.blocks * Tree::block_storage_bytes);
    EXPECT_EQ(memory.node_bytes, memory.nodes * Tree::node_storage_bytes);
    EXPECT_EQ(tree.node_count(), memory.nodes);
    EXPECT_EQ(tree.internal_memory_bytes(), memory.node_bytes);
    EXPECT_LE(memory.blocks, tree.size() / (Tree::block_capacity / 2) + 2);
    if (tree.empty()) {
      EXPECT_EQ(memory.total_bytes, sizeof(tree));
    }
  }

  template <class T>
  class SequenceTreeSpec : public ::testing::Test {};
  using TreeTypes = ::testing::Types<
      SequenceTree<BitBlock<128>, 4>,
      SequenceTree<BitBlock<192>, 8, LengthLayout::individual>,
      SequenceTree<BitBlock<>, 16>,
      SequenceTree<PackedBitBlock<512>, 4, LengthLayout::individual>,
      SequenceTree<PackedBitBlock<1024>, 8>,
      SequenceTree<PackedBitBlock<>, 16, LengthLayout::individual>,
      SequenceTree<PermutedBitBlock<>, 4>, SequenceTree<IntegerBlock, 4>,
      SequenceTree<IntegerBlock, 8, LengthLayout::individual>,
      SequenceTree<IntegerBlock, 16>,
      SequenceTree<PowerOfTwoValueBlock<std::uint64_t, 32>, 8>,
      SequenceTree<PowerOfTwoValueBlock<std::uint64_t, 256>, 32>,
      SequenceTree<IntegerBlock, 16, LengthLayout::individual>,
      SequenceTree<MoveOnlyIntegerBlock, 8>,
      SequenceTree<OddIntegerBlock, 4, LengthLayout::individual>,
      SequenceTree<ConstOnlyAlignedBlock, 4>>;
  TYPED_TEST_SUITE(SequenceTreeSpec, TreeTypes);

  TEST(SequenceTree, FullSubtreeReadsAtEveryCommonDepthAndBeyond) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    std::size_t n = Tree::block_capacity;
    for (std::size_t height = 1; height <= 6; ++height) {
      n *= 4;
      auto expected = data<Tree>(n, height);
      auto tree = build<Tree>(expected);
      ASSERT_EQ(tree.height(), height);
      check(tree, expected);
      auto suffix = tree.split_off(n / 2);
      suffix.merge(tree);
      std::rotate(expected.begin(), expected.begin() + n / 2, expected.end());
      check(suffix, expected);
      ASSERT_TRUE(suffix.test_validate());
    }
  }

  TYPED_TEST(SequenceTreeSpec, EmptyInvalidOwnershipAndAccounting) {
    using Tree = TypeParam;
    static_assert(!std::is_copy_constructible_v<Tree>);
    static_assert(std::is_nothrow_move_constructible_v<Tree>);
    static_assert(std::is_nothrow_move_assignable_v<Tree>);
    Tree empty;
    check(empty, {});
    empty.rotate_left(0, 0, -1);
    empty.merge(empty);
    EXPECT_THROW(empty[0], std::out_of_range);
    EXPECT_THROW(empty.split_off(1), std::out_of_range);
    EXPECT_THROW(empty.rotate_left(0, 1, 0), std::out_of_range);
    EXPECT_THROW(empty.rotate_left(1, 0, 0), std::out_of_range);
    auto expected = data<Tree>(Tree::block_capacity * 3 + 1);
    auto tree = build<Tree>(expected, 3);
    check(tree, expected);
    EXPECT_THROW(tree[tree.size()], std::out_of_range);
    EXPECT_THROW(tree.split_off(tree.size() + 1), std::out_of_range);
    EXPECT_THROW(tree.rotate_left(0, tree.size() + 1, 0), std::out_of_range);
    tree.merge(tree);
    tree.merge(empty);
    auto suffix = tree.split_off(tree.size());
    check(suffix, {});
    auto all = tree.split_off(0);
    check(tree, {});
    empty.merge(all);
    check(all, {});
    Tree moved(std::move(empty));
    check(empty, {});
    tree = std::move(moved);
    check(moved, {});
    auto& alias = tree;
    tree = std::move(alias);
    check(tree, expected);
  }

  TYPED_TEST(SequenceTreeSpec, EveryShortCutAndRotation) {
    using Tree = TypeParam;
    for (std::size_t n = 0; n <= 9; ++n) {
      const auto original = data<Tree>(n, n);
      for (std::size_t p = 0; p <= n; ++p) {
        auto tree = build<Tree>(original);
        auto suffix = tree.split_off(p);
        check(tree, {original.begin(), original.begin() + p});
        check(suffix, {original.begin() + p, original.end()});
        tree.merge(suffix);
        check(tree, original);
      }
      for (std::size_t l = 0; l <= n; ++l) {
        for (std::size_t r = l; r <= n; ++r) {
          for (std::size_t d = 0; d <= r - l + 1; ++d) {
            auto tree = build<Tree>(original);
            auto expected = original;
            tree.rotate_left(l, r, d);
            rotate(expected, l, r, d);
            check(tree, expected);
          }
        }
      }
    }
  }

  TYPED_TEST(SequenceTreeSpec, RandomCutsJoinsRotationsAndUnderfullSeams) {
    using Tree = TypeParam;
    std::mt19937_64 random(781293);
    auto expected = data<Tree>(Tree::block_capacity * 35 + 3);
    auto tree = build<Tree>(expected);
    for (std::size_t step = 0; step < 250; ++step) {
      SCOPED_TRACE(step);
      if (step % 7 == 0) {
        const auto p = random() % (tree.size() + 1);
        auto suffix = tree.split_off(p);
        check(tree, {expected.begin(), expected.begin() + p});
        check(suffix, {expected.begin() + p, expected.end()});
        if (step % 14 == 0) {
          suffix.merge(tree);
          tree = std::move(suffix);
          rotate(expected, 0, expected.size(), p);
        } else {
          tree.merge(suffix);
        }
      } else if (step % 11 == 0) {
        auto tail = data<Tree>(1 + random() % 17, step);
        auto donor = build<Tree>(tail);
        tree.merge(donor);
        expected.insert(expected.end(), tail.begin(), tail.end());
        check(donor, {});
      } else {
        auto l = random() % (tree.size() + 1);
        auto r = random() % (tree.size() + 1);
        if (l > r) {
          std::swap(l, r);
        }
        if (step % 5 == 0) {
          l = 1;
          r = tree.size() - 1;
        }
        const auto d =
            step % 3 == 0 ? std::numeric_limits<std::size_t>::max() : random();
        tree.rotate_left(l, r, d);
        rotate(expected, l, r, d);
      }
      check(tree, expected);
    }
  }

  TYPED_TEST(SequenceTreeSpec, UnequalHeightsTinyMergesRootBoundaries) {
    using Tree = TypeParam;
    auto expected = data<Tree>(Tree::block_capacity * 70 + 1);
    auto tree = build<Tree>(expected);
    const auto initial_height = tree.height();
    ASSERT_GT(initial_height, 0);
    for (const auto p :
         {std::size_t{1}, Tree::block_capacity - 1, Tree::block_capacity * 4,
          Tree::block_capacity * 16 + 1, tree.size() - 1}) {
      auto suffix = tree.split_off(p);
      check(tree, {expected.begin(), expected.begin() + p});
      check(suffix, {expected.begin() + p, expected.end()});
      tree.merge(suffix);
      check(tree, expected);
    }
    auto suffix = tree.split_off(1);
    EXPECT_EQ(tree.height(), 0);
    tree.merge(suffix);
    check(tree, expected);
    Tree tiny;
    std::vector<typename Tree::value_type> short_expected;
    for (std::size_t i = 0; i < 1100; ++i) {
      const auto one = data<Tree>(1, i);
      auto donor = build<Tree>(one);
      if (i % 2) {
        tiny.merge(donor);
        short_expected.push_back(one[0]);
      } else {
        donor.merge(tiny);
        tiny = std::move(donor);
        short_expected.insert(short_expected.begin(), one[0]);
      }
      ASSERT_TRUE(tiny.test_validate()) << i;
    }
    check(tiny, short_expected);
  }

  template <class B>
  class BitBlockSpec : public ::testing::Test {};
  using BlockTypes =
      ::testing::Types<BitBlock<128>, BitBlock<192>, BitBlock<>,
                       PackedBitBlock<512>, PackedBitBlock<>,
                       PermutedBitBlock<>, PermutedBitBlock<false>>;
  TYPED_TEST_SUITE(BitBlockSpec, BlockTypes);

  template <class Block>
  void check_block(const Block& b, const std::vector<bool>& expected) {
    ASSERT_EQ(b.size(), expected.size());
    const auto flat = b.flatten();
    for (std::size_t i = 0; i < Block::capacity; ++i) {
      ASSERT_EQ(bool((flat[i / 64] >> (i % 64)) & 1),
                i < expected.size() ? expected[i] : false)
          << i;
      if (i < expected.size()) {
        ASSERT_EQ(b[i], expected[i]) << i;
      }
    }
  }

  TYPED_TEST(BitBlockSpec, ExhaustiveShortPatterns) {
    using Block = TypeParam;
    EXPECT_THROW(Block({}, Block::capacity + 1), std::invalid_argument);
    EXPECT_THROW(Block({}, 1), std::invalid_argument);
    Block empty;
    empty.rotate_left(0, 0, -1);
    EXPECT_THROW(empty[0], std::out_of_range);
    for (std::size_t n = 0; n <= 6; ++n) {
      for (std::uint64_t pattern = 0; pattern < (1ULL << n); ++pattern) {
        const auto dirty = pattern | (~std::uint64_t{0} << n);
        for (std::size_t l = 0; l <= n; ++l) {
          for (std::size_t r = l; r <= n; ++r) {
            for (std::size_t d = 0; d <= r - l + 1; ++d) {
              Block b(std::span(&dirty, 1), n);
              std::vector<bool> expected(n);
              for (std::size_t i = 0; i < n; ++i) {
                expected[i] = (pattern >> i) & 1;
              }
              b.rotate_left(l, r, d);
              rotate(expected, l, r, d);
              auto words = pack(expected);
              const auto mask = n ? (1ULL << n) - 1 : 0;
              EXPECT_EQ(b.flatten()[0], words.empty() ? 0 : words[0] & mask);
              for (std::size_t i = 0; i < n; ++i) {
                ASSERT_EQ(b[i], expected[i]);
              }
            }
          }
        }
      }
    }
  }

  TEST(BitBlock, FieldReadsAcrossWordAndCircularBoundaries) {
    using Block = BitBlock<192>;
    const Block::Payload words{0x9a275b18ed430fc6ULL, 0xf037ac62819de54bULL,
                               0x76543210fedcba98ULL};
    EXPECT_EQ(Block{}.read_bits(0, 0), 0);
    for (std::size_t n : {1, 63, 64, 65, 127, 128, 129, 192}) {
      for (std::size_t origin = 0; origin < n; ++origin) {
        Block block(words, n);
        block.rotate_left(0, n, origin);
        for (std::size_t offset = 0; offset <= n; ++offset) {
          std::uint64_t expected = 0;
          EXPECT_EQ(block.read_bits(offset, 0), 0);
          for (std::size_t width = 1;
               width <= std::min<std::size_t>(64, n - offset); ++width) {
            const auto p = (origin + offset + width - 1) % n;
            expected |= ((words[p / 64] >> (p % 64)) & 1) << (width - 1);
            ASSERT_EQ(block.read_bits(offset, width), expected)
                << n << ':' << origin << ':' << offset << ':' << width;
          }
        }
      }
    }
  }

  TEST(BitBlock, NormalizationPreservesEveryCircularOrigin) {
    using Block = BitBlock<192>;
    const Block::Payload words{0x9a275b18ed430fc6ULL, 0xf037ac62819de54bULL,
                               0x76543210fedcba98ULL};
    Block empty;
    empty.normalize();
    EXPECT_TRUE(empty.empty());
    for (const std::size_t size : {1u, 63u, 64u, 65u, 127u, 128u, 192u}) {
      for (std::size_t origin = 0; origin < size; ++origin) {
        Block block(words, size);
        block.rotate_left(0, size, origin);
        block.normalize();
        ASSERT_EQ(block.size(), size);
        for (std::size_t i = 0; i < size; ++i) {
          const auto source = (origin + i) % size;
          ASSERT_EQ(block[i], (words[source / 64] >> (source % 64)) & 1);
        }
        for (std::size_t offset = 0; offset + 64 <= size; offset += 64) {
          EXPECT_EQ(block.read_normalized_word(offset),
                    block.read_bits(offset, 64));
        }
        const auto first_word =
            block.read_bits(0, std::min<std::size_t>(size, 64));
        block.normalize();
        EXPECT_EQ(block.read_bits(0, std::min<std::size_t>(size, 64)),
                  first_word);
      }
    }
  }

  TEST(BitBlock, FieldWritesPreserveOutsideBitsAcrossCircularBoundaries) {
    using Block = BitBlock<192>;
    const Block::Payload words{0x9a275b18ed430fc6ULL, 0xf037ac62819de54bULL,
                               0x76543210fedcba98ULL};
    constexpr std::uint64_t replacement = 0x123456789abcdef0ULL;
    Block{}.write_bits(0, 0, replacement);
    for (std::size_t n : {1, 63, 64, 65, 127, 128, 129, 192}) {
      for (std::size_t origin :
           {std::size_t{0}, std::size_t{1} % n, n / 2, n - 1}) {
        for (std::size_t offset = 0; offset <= n; ++offset) {
          for (std::size_t width : {0, 1, 7, 32, 63, 64}) {
            if (width > n - offset) {
              continue;
            }
            Block block(words, n);
            block.rotate_left(0, n, origin);
            block.write_bits(offset, width, replacement);
            for (std::size_t i = 0; i < n; ++i) {
              const auto p = (origin + i) % n;
              const bool expected = i >= offset && i - offset < width
                                        ? (replacement >> (i - offset)) & 1
                                        : (words[p / 64] >> (p % 64)) & 1;
              ASSERT_EQ(block[i], expected);
            }
          }
        }
      }
    }
  }

  TYPED_TEST(BitBlockSpec, RepresentativeOriginsRedistributionAndDirtyPadding) {
    using Block = TypeParam;
    std::mt19937_64 random(865124);
    for (const std::size_t n :
         {std::size_t{65}, Block::capacity - 3, Block::capacity}) {
      auto initial = data<SequenceTree<Block>>(n);
      for (const auto origin : representative_origins(n)) {
        Block a(pack(initial), n);
        auto expected = initial;
        a.rotate_left(0, n, origin);
        rotate(expected, 0, n, origin);
        const auto l = origin % 64;
        const auto r = std::min(n, l + 65);
        const auto d = random();
        a.rotate_left(l, r, d);
        rotate(expected, l, r, d);
        Block b;
        const auto cut = origin % (n + 1);
        a.redistribute(b, cut);
        check_block(a, {expected.begin(), expected.begin() + cut});
        check_block(b, {expected.begin() + cut, expected.end()});
        a.redistribute(b, n);
        check_block(a, expected);
        check_block(b, {});
      }
    }
    // Representative origins cover word and 128-bit chunk boundaries;
    // randomized cases cover arbitrary block sizes and redistribution cuts
    // without making the coverage configuration exhaust every physical origin.
    for (std::size_t step = 0; step < 128; ++step) {
      const auto n = random() % (Block::capacity + 1);
      const auto m = random() % (Block::capacity + 1);
      auto avalues = data<SequenceTree<Block>>(n, step);
      auto bvalues = data<SequenceTree<Block>>(m, step + 1);
      Block a(pack(avalues), n), b(pack(bvalues), m);
      a.rotate_left(0, n, 65);
      b.rotate_left(0, m, 127);
      rotate(avalues, 0, n, 65);
      rotate(bvalues, 0, m, 127);
      avalues.insert(avalues.end(), bvalues.begin(), bvalues.end());
      const auto low = n + m > Block::capacity ? n + m - Block::capacity : 0;
      const auto high = std::min(Block::capacity, n + m);
      const auto cut = low + random() % (high - low + 1);
      a.redistribute(b, cut);
      check_block(a, {avalues.begin(), avalues.begin() + cut});
      check_block(b, {avalues.begin() + cut, avalues.end()});
    }
  }

  TEST(PermutedBitBlock, MapKernelsPreserveRepresentationAndContents) {
    using Block = PermutedBitBlock<>;
    static_assert(std::is_trivially_copyable_v<Block>);
    std::mt19937_64 random(1282048);
    for (const std::size_t n : {257, 385, 513, 1921, 2045, 2048}) {
      auto expected = data<SequenceTree<Block>>(n);
      Block b(pack(expected), n);
      auto apply = [&](std::size_t l, std::size_t r, std::size_t d) {
        b.rotate_left(l, r, d);
        rotate(expected, l, r, d);
        check_block(b, expected);
      };
      for (std::size_t l = 0; l <= n / 128; ++l) {
        for (std::size_t r = l; r <= n / 128; ++r) {
          for (std::size_t d = 0; d <= r - l; ++d) {
            apply(l * 128, r * 128, d * 128);
          }
        }
      }
      // Whole rotations use varied deterministic distances so the map state
      // sees origins across the block without repeating the same branch
      // coverage n times. The aligned triples above retain exhaustive
      // map-kernel controls.
      for (std::size_t iteration = 0; iteration < 64; ++iteration) {
        apply(0, n, 1 + random() % n);
        apply(1, n - 1, random());
        apply(63, 257, 65);
      }
      Block::Payload words;
      words.fill(~std::uint64_t{0});
      Block uniform(words, n);
      uniform.rotate_left(0, 256, 128);
      uniform.rotate_left(0, n, 127);
      std::array<unsigned char, sizeof(Block)> before, after;
      std::memcpy(before.data(), &uniform, sizeof(uniform));
      uniform.rotate_left(1, n - 1, -1);
      uniform.rotate_left(63, 257, 65);
      std::memcpy(after.data(), &uniform, sizeof(uniform));
      EXPECT_EQ(before, after);
    }
  }

  TEST(PermutedBitBlock, WrappedMappedRedistributionSmoke) {
    using Block = PermutedBitBlock<>;
    auto expected = data<SequenceTree<Block>>(1921);
    Block a(pack(expected), expected.size());
    for (const auto h : {std::size_t{0}, std::size_t{1}, std::size_t{128},
                         std::size_t{1900}}) {
      a.rotate_left(0, 512, 128);
      rotate(expected, 0, 512, 128);
      a.rotate_left(0, a.size(), h);
      rotate(expected, 0, expected.size(), h);
      a.rotate_left(63, 257, 65);
      rotate(expected, 63, 257, 65);
      Block b;
      a.redistribute(b, 127);
      check_block(a, {expected.begin(), expected.begin() + 127});
      check_block(b, {expected.begin() + 127, expected.end()});
      a.redistribute(b, expected.size());
      check_block(a, expected);
    }
  }

  TEST(SequenceTreeSeam, TwoTinyLeavesRequireTheirNeighbors) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    for (std::size_t left_count : {1, 9, 33, 129}) {
      for (std::size_t right_count : {8, 16, 40, 136}) {
        auto avalues = data<Tree>(left_count, left_count);
        auto bvalues = data<Tree>(right_count, right_count);
        auto a = build<Tree>(avalues);
        auto b = build<Tree>(bvalues);
        auto tail = b.split_off(7);
        b = std::move(tail);
        bvalues.erase(bvalues.begin(), bvalues.begin() + 7);
        a.merge(b);
        avalues.insert(avalues.end(), bvalues.begin(), bvalues.end());
        check(a, avalues);
        check(b, {});
      }
    }
  }

  TEST(SequenceTreeLayout, TypedAlignedPhysicalBudgets) {
    static_assert(sizeof(PackedBitBlock<512>) == 64);
    static_assert(sizeof(PackedBitBlock<1024>) == 128);
    static_assert(sizeof(PackedBitBlock<2048>) == 256);
    static_assert(sizeof(PackedBitBlock<4096>) == 512);
    static_assert(PackedBitBlock<2048>::capacity == 1920);
    static_assert(alignof(PackedBitBlock<>) == alignof(pixie::CacheLine));
    static_assert(BitBlock<>::capacity == 2048);
    static_assert(sizeof(BitBlock<>) == 320);
    static_assert(BitBlock<>::payload_offset_bytes == 16);
    static_assert(BitBlock<>::payload_alignment == alignof(std::uint64_t));
    EXPECT_EQ((SequenceTree<IntegerBlock, 4>::node_storage_bytes), 64);
    EXPECT_EQ((SequenceTree<IntegerBlock, 8>::node_storage_bytes), 128);
    EXPECT_EQ((SequenceTree<IntegerBlock, 16>::node_storage_bytes), 256);
    using Tree = SequenceTree<PackedBitBlock<>>;
    auto tree = build<Tree>(data<Tree>(20000));
    EXPECT_TRUE(tree.test_validate());
    EXPECT_EQ(tree.payload_capacity_bytes(), tree.block_count() * 240);
    EXPECT_EQ(tree.metadata_bytes() + tree.payload_capacity_bytes(),
              tree.memory_usage_bytes());
  }

  TEST(SequenceTreeLayout, ConstOnlyReadsAndOverAlignedLeaves) {
    using Tree = SequenceTree<ConstOnlyAlignedBlock, 4>;
    static_assert(alignof(ConstOnlyAlignedBlock) == 128);
    static_assert(Tree::block_alignment == alignof(ConstOnlyAlignedBlock));
    auto expected = data<Tree>(257);
    auto tree = build<Tree>(expected);
    auto verify = [&] {
      check(std::as_const(tree), expected);
      std::size_t visited = 0;
      tree.for_each_block([&](const ConstOnlyAlignedBlock& block) {
        EXPECT_EQ(reinterpret_cast<std::uintptr_t>(&block) % 128, 0);
        ++visited;
      });
      EXPECT_EQ(visited, tree.block_count());
    };
    verify();
    auto tail = tree.split_off(17);
    check(tail, {expected.begin() + 17, expected.end()});
    tree.merge(tail);
    verify();
    tree.rotate_left(1, tree.size() - 1, 73);
    rotate(expected, 1, expected.size() - 1, 73);
    verify();
  }

  TEST(SequenceTreeConstruction, SinglePassAndThrowingProducer) {
    using Tree = SequenceTree<IntegerBlock>;
    std::istringstream stream("1 2 3 4 5 6 7 8 9 10 11 12 13");
    auto blocks = std::ranges::istream_view<std::uint32_t>(stream) |
                  std::views::transform([](std::uint32_t n) {
                    return IntegerBlock(std::span(&n, 1));
                  });
    auto tree = Tree::from_blocks(blocks);
    check(tree, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13});
    const auto live = Tree::test_counters.live_allocations;
    auto throwing = std::views::iota(0, 100) | std::views::transform([](int n) {
                      if (n == 47) {
                        throw std::runtime_error("producer");
                      }
                      const auto value = static_cast<std::uint32_t>(n);
                      return IntegerBlock(std::span(&value, 1));
                    });
    EXPECT_THROW(Tree::from_blocks(throwing), std::runtime_error);
    EXPECT_EQ(Tree::test_counters.live_allocations, live);
  }

  TEST(SequenceTreeConstruction, StreamingBuildHasLinearAllocationCount) {
    using Tree = SequenceTree<IntegerBlock, 8>;
    for (const std::size_t leaves : {1, 64, 4097, 65536}) {
      auto blocks = std::views::iota(std::size_t{0}, leaves) |
                    std::views::transform([](std::size_t index) {
                      std::array<std::uint32_t, IntegerBlock::capacity> values;
                      for (std::size_t i = 0; i < values.size(); ++i) {
                        values[i] = index * IntegerBlock::capacity + i;
                      }
                      return IntegerBlock(values);
                    });
      Tree::test_reset_counters();
      IntegerBlock::copied_elements = 0;
      auto tree = Tree::from_blocks(blocks);
      const auto counters = Tree::test_counters;
      std::size_t carry_nodes = 0;
      std::size_t levels = 1;
      for (auto remaining = leaves / 8; remaining; remaining /= 8) {
        carry_nodes += remaining;
        ++levels;
      }
      // One leaf allocation per block, geometric bottom-up carries, and just
      // one height-bounded spare pool for the final fringe, not per-leaf
      // merges.
      EXPECT_LE(counters.allocations,
                leaves + carry_nodes + (leaves > 1 ? 3 * levels : 0));
      EXPECT_LE(counters.node_visits, 4 * levels);
      EXPECT_LE(counters.child_transfers, 4 * leaves + 100 * levels);
      EXPECT_LE(IntegerBlock::copied_elements,
                2 * leaves * IntegerBlock::capacity);
      EXPECT_EQ(tree.size(), leaves * IntegerBlock::capacity);
      EXPECT_EQ(tree.block_count(), leaves);
      ASSERT_TRUE(tree.test_validate());
      std::size_t expected = 0;
      tree.for_each_block([&](const IntegerBlock& block) {
        for (std::size_t i = 0; i < block.size(); ++i) {
          ASSERT_EQ(block[i], expected++);
        }
      });
      EXPECT_EQ(expected, tree.size());
    }
  }

  template <class Tree>
  void every_preflight_allocation_is_strong_and_leak_free() {
    const auto original = data<Tree>(153);
    const auto donor_values = data<Tree>(35, 21);
    for (int operation = 0; operation < 4; ++operation) {
      std::size_t allocations = 0;
      {
        auto tree = build<Tree>(original);
        auto donor = build<Tree>(donor_values);
        // Create underfull leaves at both joining ends, including a third
        // neighbor.
        if (operation == 1) {
          auto discarded = tree.split_off(145);
          auto tail = donor.split_off(7);
          donor = std::move(tail);
        }
        Tree::test_reset_counters();
        if (operation == 0) {
          auto tail = tree.split_off(73);
        } else if (operation == 1) {
          tree.merge(donor);
        } else if (operation == 2) {
          tree.rotate_left(1, 151, 71);
        } else {
          tree.rotate_left(0, tree.size(), 71);
        }
        allocations = Tree::test_counters.allocations;
        ASSERT_GT(allocations, 0);
      }
      for (std::size_t failure = 0; failure <= allocations; ++failure) {
        SCOPED_TRACE(::testing::Message() << operation << ":" << failure);
        auto tree = build<Tree>(original);
        auto donor = build<Tree>(donor_values);
        auto expected = original;
        auto donor_expected = donor_values;
        if (operation == 1) {
          auto discarded = tree.split_off(145);
          expected.resize(145);
          auto tail = donor.split_off(7);
          donor = std::move(tail);
          donor_expected.erase(donor_expected.begin(),
                               donor_expected.begin() + 7);
        }
        const auto identities = tree.test_leaf_identities();
        const auto donor_identities = donor.test_leaf_identities();
        const auto node_identities = tree.test_internal_node_identities();
        const auto donor_node_identities =
            donor.test_internal_node_identities();
        const auto live = Tree::test_counters.live_allocations;
        const auto live_bytes = Tree::test_counters.live_bytes;
        Tree tail;
        auto apply = [&] {
          if (operation == 0) {
            tail = tree.split_off(73);
          } else if (operation == 1) {
            tree.merge(donor);
          } else if (operation == 2) {
            tree.rotate_left(1, 151, 71);
          } else {
            tree.rotate_left(0, tree.size(), 71);
          }
        };
        Tree::test_reset_counters();
        Tree::test_fail_after(static_cast<std::ptrdiff_t>(failure));
        if (failure < allocations) {
          EXPECT_THROW(apply(), std::bad_alloc);
        } else {
          // All preflight allocations succeed, but any allocation during commit
          // would fail. This exercises the nonallocating commit with recycling.
          EXPECT_NO_THROW(apply());
        }
        Tree::test_fail_after(-1);
        if (failure < allocations) {
          EXPECT_EQ(Tree::test_counters.live_allocations, live);
          EXPECT_EQ(Tree::test_counters.live_bytes, live_bytes);
          EXPECT_EQ(tree.test_leaf_identities(), identities);
          EXPECT_EQ(donor.test_leaf_identities(), donor_identities);
          EXPECT_EQ(tree.test_internal_node_identities(), node_identities);
          EXPECT_EQ(donor.test_internal_node_identities(),
                    donor_node_identities);
        } else {
          EXPECT_EQ(Tree::test_counters.allocations, allocations);
          if (operation == 0) {
            check(tail, {expected.begin() + 73, expected.end()});
            expected.resize(73);
          } else if (operation == 1) {
            expected.insert(expected.end(), donor_expected.begin(),
                            donor_expected.end());
            donor_expected.clear();
          } else if (operation == 2) {
            rotate(expected, 1, 151, 71);
          } else {
            rotate(expected, 0, expected.size(), 71);
          }
        }
        check(tree, expected);
        check(donor, donor_expected);
      }
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  TEST(SequenceTreeFailure, ConstructionFailureReclaimsPartialLevels) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    const auto expected = data<Tree>(511);
    Tree::test_reset_counters();
    { auto tree = build<Tree>(expected); }
    const auto allocations = Tree::test_counters.allocations;
    for (std::size_t i = 0; i < allocations; ++i) {
      Tree::test_fail_after(i);
      EXPECT_THROW(build<Tree>(expected), std::bad_alloc);
      Tree::test_fail_after(-1);
      EXPECT_EQ(Tree::test_counters.live_allocations, 0);
    }
  }

  TEST(SequenceTreeFailure, NoRepairMergeUsesHeightGapBudgetAndRollsBack) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    const std::array<std::pair<std::size_t, std::size_t>, 6> shapes{
        {{36, 36}, {512, 8}, {8, 512}, {512, 32}, {32, 512}, {512, 512}}};
    for (const auto& [left_size, right_size] : shapes) {
      SCOPED_TRACE(::testing::Message() << left_size << ":" << right_size);
      const auto left_values = data<Tree>(left_size, 13);
      const auto right_values = data<Tree>(right_size, 17);
      std::size_t budget = 0;
      {
        auto tree = build<Tree>(left_values);
        auto donor = build<Tree>(right_values);
        const auto taller = std::max(tree.height(), donor.height());
        const auto gap = taller - std::min(tree.height(), donor.height());
        budget = gap + 1;
        Tree::test_reset_counters();
        IntegerBlock::copied_elements = 0;
        tree.merge(donor);
        EXPECT_EQ(Tree::test_counters.allocations, budget);
        EXPECT_LT(Tree::test_counters.allocations, 64 * (taller + 2));
        EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
        EXPECT_EQ(IntegerBlock::copied_elements, 0);
        auto expected = left_values;
        expected.insert(expected.end(), right_values.begin(),
                        right_values.end());
        check(tree, expected);
        check(donor, {});
      }
      for (std::size_t failure = 0; failure < budget; ++failure) {
        SCOPED_TRACE(failure);
        auto tree = build<Tree>(left_values);
        auto donor = build<Tree>(right_values);
        const auto leaves = tree.test_leaf_identities();
        const auto donor_leaves = donor.test_leaf_identities();
        const auto nodes = tree.test_internal_node_identities();
        const auto donor_nodes = donor.test_internal_node_identities();
        const auto live = Tree::test_counters.live_allocations;
        const auto live_bytes = Tree::test_counters.live_bytes;
        Tree::test_fail_after(static_cast<std::ptrdiff_t>(failure));
        EXPECT_THROW(tree.merge(donor), std::bad_alloc);
        Tree::test_fail_after(-1);
        EXPECT_EQ(Tree::test_counters.live_allocations, live);
        EXPECT_EQ(Tree::test_counters.live_bytes, live_bytes);
        EXPECT_EQ(tree.test_leaf_identities(), leaves);
        EXPECT_EQ(donor.test_leaf_identities(), donor_leaves);
        EXPECT_EQ(tree.test_internal_node_identities(), nodes);
        EXPECT_EQ(donor.test_internal_node_identities(), donor_nodes);
        check(tree, left_values);
        check(donor, right_values);
      }
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  TEST(SequenceTree, BorrowedShortSeamBecomesExterior) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    constexpr auto capacity = Tree::block_capacity;
    const auto left_values = data<Tree>(2 * capacity, 31);
    auto prefix = build<Tree>(left_values);
    auto left = prefix.split_off(capacity - 1);
    auto discarded = left.split_off(2);
    ASSERT_EQ(left.block_count(), 2);
    const auto right_values = data<Tree>(65 * capacity, 37);
    auto right_prefix = build<Tree>(right_values);
    auto right = right_prefix.split_off(capacity - 1);
    // Both seam leaves and the one remaining left neighbor contain one value.
    // Borrowing the neighbor exhausts the left tree, so the resulting three
    // values may legally remain in an underfull exterior leaf.
    std::vector<Tree::value_type> expected(left_values.begin() + capacity - 1,
                                           left_values.begin() + capacity + 1);
    expected.insert(expected.end(), right_values.begin() + capacity - 1,
                    right_values.end());
    left.merge(right);
    check(left, expected);
    check(right, {});
  }

  TEST(SequenceTreeFailure, FittingLeafDonorNeedsNoAncestorAllocations) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    for (const std::size_t size : {17u, 65u, 513u}) {
      auto expected = data<Tree>(size, 31);
      auto tree = build<Tree>(expected);
      ASSERT_GT(tree.height(), 0);
      const auto nodes = tree.test_internal_node_identities();
      const auto leaves = tree.test_leaf_identities();
      const auto incoming = data<Tree>(1, 37);
      auto donor = build<Tree>(incoming);
      Tree::test_fail_after(0);
      EXPECT_NO_THROW(tree.merge(donor));
      Tree::test_fail_after(-1);
      expected.insert(expected.end(), incoming.begin(), incoming.end());
      check(tree, expected);
      check(donor, {});
      EXPECT_EQ(tree.test_internal_node_identities(), nodes);
      EXPECT_EQ(tree.test_leaf_identities(), leaves);
    }
  }

  TEST(SequenceTreeFailure, LeafRootMergeAllocatesOnlyRequiredParent) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    const std::array<std::pair<std::size_t, std::size_t>, 8> shapes{
        {{1, 1}, {1, 7}, {4, 4}, {3, 5}, {8, 8}, {8, 1}, {1, 8}, {4, 5}}};
    for (const auto& [left_size, right_size] : shapes) {
      SCOPED_TRACE(::testing::Message() << left_size << ":" << right_size);
      const auto left_values = data<Tree>(left_size, 31);
      const auto right_values = data<Tree>(right_size, 37);
      const std::size_t budget = left_size + right_size > Tree::block_capacity;
      for (std::size_t failure = 0; failure <= budget; ++failure) {
        auto tree = build<Tree>(left_values);
        auto donor = build<Tree>(right_values);
        ASSERT_EQ(tree.height(), 0);
        ASSERT_EQ(donor.height(), 0);
        const auto leaves = tree.test_leaf_identities();
        const auto donor_leaves = donor.test_leaf_identities();
        Tree::test_reset_counters();
        IntegerBlock::copied_elements = 0;
        const auto live = Tree::test_counters.live_allocations;
        const auto live_bytes = Tree::test_counters.live_bytes;
        Tree::test_fail_after(static_cast<std::ptrdiff_t>(failure));
        if (failure < budget) {
          EXPECT_THROW(tree.merge(donor), std::bad_alloc);
        } else {
          // With budget zero, failure injection rejects any attempted
          // allocation.
          EXPECT_NO_THROW(tree.merge(donor));
        }
        Tree::test_fail_after(-1);
        if (failure < budget) {
          EXPECT_EQ(Tree::test_counters.live_allocations, live);
          EXPECT_EQ(Tree::test_counters.live_bytes, live_bytes);
          EXPECT_EQ(tree.test_leaf_identities(), leaves);
          EXPECT_EQ(donor.test_leaf_identities(), donor_leaves);
          check(tree, left_values);
          check(donor, right_values);
        } else {
          EXPECT_EQ(Tree::test_counters.allocations, budget);
          if (left_size >= Tree::block_capacity / 2 &&
              right_size >= Tree::block_capacity / 2 && budget != 0) {
            EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
            EXPECT_EQ(IntegerBlock::copied_elements, 0);
          }
          EXPECT_EQ(tree.height(), budget);
          EXPECT_EQ(tree.node_count(), budget);
          EXPECT_EQ(tree.block_count(), budget + 1);
          auto expected = left_values;
          expected.insert(expected.end(), right_values.begin(),
                          right_values.end());
          check(tree, expected);
          check(donor, {});
        }
      }
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  std::size_t common_identity_count(std::vector<const void*> before,
                                    std::vector<const void*> after) {
    std::sort(before.begin(), before.end(), std::less<const void*>{});
    std::sort(after.begin(), after.end(), std::less<const void*>{});
    EXPECT_EQ(std::adjacent_find(before.begin(), before.end()), before.end());
    EXPECT_EQ(std::adjacent_find(after.begin(), after.end()), after.end());
    std::vector<const void*> common;
    std::set_intersection(before.begin(), before.end(), after.begin(),
                          after.end(), std::back_inserter(common),
                          std::less<const void*>{});
    return common.size();
  }

  TEST(SequenceTreeRecycling, FullSpineJoinRetainsDismantledAllocations) {
    using Tree = SequenceTree<IntegerBlock, 4>;
    for (const std::size_t donor_leaves : {1, 64}) {
      auto expected = data<Tree>(64 * Tree::block_capacity);
      const auto donor_values =
          data<Tree>(donor_leaves * Tree::block_capacity, 19);
      auto tree = build<Tree>(expected);
      auto donor = build<Tree>(donor_values);
      auto nodes = tree.test_internal_node_identities();
      const auto donor_nodes = donor.test_internal_node_identities();
      nodes.insert(nodes.end(), donor_nodes.begin(), donor_nodes.end());
      auto leaves = tree.test_leaf_identities();
      const auto donor_ids = donor.test_leaf_identities();
      leaves.insert(leaves.end(), donor_ids.begin(), donor_ids.end());
      const auto budget = tree.height() - donor.height() + 1;
      const auto bytes = Tree::test_counters.live_bytes;
      Tree::test_reset_counters();
      Tree::test_fail_after(budget);
      tree.merge(donor);
      Tree::test_fail_after(-1);
      EXPECT_EQ(Tree::test_counters.allocations, budget);
      EXPECT_EQ(Tree::test_counters.peak_bytes - bytes,
                budget * Tree::node_storage_bytes);
      EXPECT_EQ(
          common_identity_count(nodes, tree.test_internal_node_identities()),
          nodes.size());
      EXPECT_EQ(tree.test_leaf_identities(), leaves);
      expected.insert(expected.end(), donor_values.begin(), donor_values.end());
      check(tree, expected);
      check(donor, {});
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  TEST(SequenceTreeRecycling, BalancedRootsJoinWithoutTransferringChildren) {
    using Tree = SequenceTree<IntegerBlock, 8>;
    for (std::size_t left = 4; left <= 8; ++left) {
      for (std::size_t right = 4; right <= 8; ++right) {
        if (left + right <= 8) {
          continue;
        }
        SCOPED_TRACE(::testing::Message() << left << ':' << right);
        auto expected = data<Tree>(left * Tree::block_capacity);
        const auto incoming = data<Tree>(right * Tree::block_capacity, 19);
        auto tree = build<Tree>(expected);
        auto donor = build<Tree>(incoming);
        auto nodes = tree.test_internal_node_identities();
        const auto donor_nodes = donor.test_internal_node_identities();
        nodes.insert(nodes.end(), donor_nodes.begin(), donor_nodes.end());
        Tree::test_reset_counters();
        tree.merge(donor);
        EXPECT_EQ(Tree::test_counters.child_transfers, 2);
        EXPECT_EQ(
            common_identity_count(nodes, tree.test_internal_node_identities()),
            nodes.size());
        expected.insert(expected.end(), incoming.begin(), incoming.end());
        check(tree, expected);
        check(donor, {});
        // Exercise subsequent edits through the retained, possibly unequal
        // roots.
        tree.rotate_left(1, tree.size() - 1, Tree::block_capacity + 1);
        rotate(expected, 1, expected.size() - 1, Tree::block_capacity + 1);
        check(tree, expected);
      }
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  template <class T>
  class SequenceTreeRecyclingSpec : public ::testing::Test {};
  using RecyclingTreeTypes = ::testing::Types<
      SequenceTree<IntegerBlock, 4>,
      SequenceTree<IntegerBlock, 4, LengthLayout::individual>,
      SequenceTree<IntegerBlock, 8>,
      SequenceTree<IntegerBlock, 8, LengthLayout::individual>,
      SequenceTree<IntegerBlock, 16>,
      SequenceTree<IntegerBlock, 16, LengthLayout::individual>>;
  TYPED_TEST_SUITE(SequenceTreeRecyclingSpec, RecyclingTreeTypes);

  TYPED_TEST(SequenceTreeRecyclingSpec,
             EveryPreflightAllocationIsStrongAndLeakFree) {
    every_preflight_allocation_is_strong_and_leak_free<TypeParam>();
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, RotationDeficitAndExactLeafSpares) {
    using Tree = TypeParam;
    const auto original = data<Tree>(4096 * Tree::block_capacity);
    // Whole/partial, aligned/unaligned, and two cuts inside the same leaf.
    const std::array<std::array<std::size_t, 3>, 7> ranges{{
        {0, original.size(), 8},
        {0, original.size(), 9},
        {8, original.size() - 8, 16},
        {1, original.size() - 1, 17},
        {1, 10, 1},
        {0, original.size() - 1, 8},
        {1, original.size(), 7},
    }};
    for (const auto& [left, right, distance] : ranges) {
      SCOPED_TRACE(::testing::Message()
                   << left << ":" << right << ":" << distance);
      auto expected = original;
      auto tree = build<Tree>(original);
      const bool whole = left == 0 && right == original.size();
      const auto nodes =
          whole ? 5 * tree.height() + 4 : 15 * tree.height() + 24;
      const auto leaves = std::size_t(left % Tree::block_capacity != 0) +
                          ((left + distance) % Tree::block_capacity != 0) +
                          (right % Tree::block_capacity != 0);
      const auto bytes = Tree::test_counters.live_bytes;
      const auto ids = tree.test_leaf_identities();
      Tree::test_reset_counters();
      Tree::test_fail_after(nodes + leaves);
      tree.rotate_left(left, right, distance);
      Tree::test_fail_after(-1);
      EXPECT_EQ(Tree::test_counters.allocations, nodes + leaves);
      EXPECT_EQ(Tree::test_counters.peak_bytes - bytes,
                nodes * Tree::node_storage_bytes +
                    leaves * Tree::block_storage_bytes);
      rotate(expected, left, right, distance);
      check(tree, expected);
      if (leaves == 0) {
        EXPECT_EQ(common_identity_count(ids, tree.test_leaf_identities()),
                  ids.size());
        EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
      }
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, CompleteChildRotationsNeverAllocate) {
    using Tree = TypeParam;
    /**
     * @brief Force partial internal nodes and unequal child totals.
     * @details The last leaf has five elements, so every rotated position is
     * legal.
     */
    auto expected = data<Tree>(257 * Tree::block_capacity - 3);
    auto tree = build<Tree>(expected);
    const auto nodes = tree.test_internal_node_identities();
    const auto leaves = tree.test_leaf_identities();
    const auto shapes = tree.test_child_boundaries();
    for (auto index : {std::size_t{0}, shapes.size() / 2, shapes.size() - 1}) {
      const auto& boundaries = shapes[index];
      for (std::size_t a = 0; a + 2 < boundaries.size(); ++a) {
        for (std::size_t b = a + 1; b + 1 < boundaries.size(); ++b) {
          for (std::size_t c = b + 1; c < boundaries.size(); ++c) {
            const auto l = boundaries[a], m = boundaries[b], r = boundaries[c];
            SCOPED_TRACE(::testing::Message()
                         << index << ':' << l << ':' << m << ':' << r);
            Tree::test_reset_counters();
            const auto bytes = Tree::test_counters.live_bytes;
            Tree::test_fail_after(0);
            EXPECT_NO_THROW(tree.rotate_left(l, r, m - l));
            Tree::test_fail_after(-1);
            EXPECT_EQ(Tree::test_counters.allocations, 0);
            EXPECT_EQ(Tree::test_counters.peak_bytes, bytes);
            EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
            if (index == 0) {
              // One covering-root visit; only included global endpoints require
              // child-spine occupancy reads, without revisiting the root.
              EXPECT_EQ(
                  Tree::test_counters.node_visits,
                  1 + (tree.height() - 1) * ((l == 0) + (r == tree.size())));
            }
            rotate(expected, l, r, m - l);
            check(tree, expected);
            EXPECT_EQ(common_identity_count(
                          nodes, tree.test_internal_node_identities()),
                      nodes.size());
            EXPECT_EQ(
                common_identity_count(leaves, tree.test_leaf_identities()),
                leaves.size());
            /** @brief The inverse stays child aligned despite unequal totals.
             */
            Tree::test_fail_after(0);
            EXPECT_NO_THROW(tree.rotate_left(l, r, r - m));
            Tree::test_fail_after(-1);
            rotate(expected, l, r, r - m);
            EXPECT_EQ(tree.test_internal_node_identities(), nodes);
            EXPECT_EQ(tree.test_leaf_identities(), leaves);
            check(tree, expected);
          }
        }
      }
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, ChildBoundaryReadsOnlyNecessarySpines) {
    using Tree = TypeParam;
    const auto original = data<Tree>(4096 * Tree::block_capacity);
    auto tree = build<Tree>(original);
    const auto shapes = tree.test_child_boundaries();
    const auto leaves = tree.test_leaf_identities();
    const auto nodes = tree.test_internal_node_identities();
    // Root-first preorder starts with the first-child chain down to leaf
    // parents.
    for (std::size_t depth = 0; depth < tree.height(); ++depth) {
      const auto& boundaries = shapes[depth];
      ASSERT_GE(boundaries.size(), 5);
      for (const bool prefix : {false, true}) {
        for (const bool suffix : {false, true}) {
          const auto left = boundaries[prefix ? 0 : 1];
          const auto right = boundaries[boundaries.size() - (suffix ? 1 : 2)];
          const auto middle = boundaries[2];
          SCOPED_TRACE(::testing::Message()
                       << depth << ':' << prefix << ':' << suffix);
          Tree::test_reset_counters();
          const auto bytes = Tree::test_counters.live_bytes;
          Tree::test_fail_after(0);
          EXPECT_NO_THROW(tree.rotate_left(left, right, middle - left));
          Tree::test_fail_after(-1);
          EXPECT_EQ(Tree::test_counters.node_visits,
                    depth + 1 +
                        (tree.height() - depth - 1) *
                            ((left == 0) + (right == tree.size())));
          EXPECT_EQ(Tree::test_counters.allocations, 0);
          EXPECT_EQ(Tree::test_counters.peak_bytes, bytes);
          EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
          auto expected = original;
          rotate(expected, left, right, middle - left);
          check(tree, expected);
          Tree::test_fail_after(0);
          EXPECT_NO_THROW(tree.rotate_left(left, right, right - middle));
          Tree::test_fail_after(-1);
          EXPECT_EQ(tree.test_leaf_identities(), leaves);
          EXPECT_EQ(tree.test_internal_node_identities(), nodes);
          check(tree, original);
        }
      }
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, UnderfullExteriorChildFallsBack) {
    using Tree = TypeParam;
    for (const bool first : {false, true}) {
      for (const bool whole : {false, true}) {
        SCOPED_TRACE(::testing::Message() << first << ':' << whole);
        auto expected = data<Tree>(4096 * Tree::block_capacity);
        auto tree = build<Tree>(expected);
        if (first) {
          tree = tree.split_off(Tree::block_capacity - 1);
          expected.erase(expected.begin(),
                         expected.begin() + Tree::block_capacity - 1);
        } else {
          const auto keep = tree.size() - Tree::block_capacity + 1;
          auto discarded = tree.split_off(keep);
          expected.resize(keep);
        }
        std::vector<std::size_t> sizes;
        tree.for_each_block(
            [&](const auto& block) { sizes.push_back(block.size()); });
        ASSERT_EQ(first ? sizes.front() : sizes.back(), 1);
        const auto boundaries = tree.test_child_boundaries().front();
        ASSERT_GE(boundaries.size(), 4);
        const auto left = whole || first ? 0 : boundaries[1];
        const auto right =
            whole || !first ? tree.size() : boundaries[boundaries.size() - 2];
        const auto distance = boundaries[whole || first ? 1 : 2] - left;
        const auto leaves = tree.test_leaf_identities();
        const auto nodes = tree.test_internal_node_identities();
        Tree::test_reset_counters();
        const auto bytes = Tree::test_counters.live_bytes;
        Tree::test_fail_after(0);
        EXPECT_THROW(tree.rotate_left(left, right, distance), std::bad_alloc);
        Tree::test_fail_after(-1);
        EXPECT_EQ(Tree::test_counters.payload_mutations, 0);
        EXPECT_EQ(Tree::test_counters.child_transfers, 0);
        EXPECT_EQ(Tree::test_counters.live_bytes, bytes);
        EXPECT_EQ(tree.test_leaf_identities(), leaves);
        EXPECT_EQ(tree.test_internal_node_identities(), nodes);
        check(tree, expected);
        tree.rotate_left(left, right, distance);
        rotate(expected, left, right, distance);
        check(tree, expected);
      }
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, SingleLeafRotationUsesOneDescent) {
    using Tree = TypeParam;
    for (const std::size_t leaf_count : {1, 512}) {
      auto expected = data<Tree>(leaf_count * Tree::block_capacity);
      auto tree = build<Tree>(expected);
      const auto leaves = tree.test_leaf_identities();
      const auto nodes = tree.test_internal_node_identities();
      const auto start = (leaf_count / 2) * Tree::block_capacity;
      for (std::size_t l = 0; l < Tree::block_capacity; ++l) {
        for (auto r = l + 2; r <= Tree::block_capacity; ++r) {
          Tree::test_reset_counters();
          Tree::test_fail_after(0);
          EXPECT_NO_THROW(tree.rotate_left(start + l, start + r, 1));
          Tree::test_fail_after(-1);
          EXPECT_EQ(Tree::test_counters.node_visits, tree.height());
          EXPECT_EQ(Tree::test_counters.allocations, 0);
          EXPECT_EQ(Tree::test_counters.payload_mutations, 1);
          EXPECT_EQ(Tree::test_counters.child_transfers, 0);
          rotate(expected, start + l, start + r, 1);
          check(tree, expected);
          EXPECT_EQ(tree.test_leaf_identities(), leaves);
          EXPECT_EQ(tree.test_internal_node_identities(), nodes);
        }
      }
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec,
             SplitPrefixBudgetAcrossFragmentedShapes) {
    using Tree = TypeParam;
    /**
     * @brief Exhaust cuts in full and repeatedly fragmented trees.
     * @details Force empty, singleton and grouped siblings at different
     * heights, including F4's minimally occupied ancestors. Sanitizer builds
     * assert pool availability at every pop; fail-after forbids extra commit
     * allocation.
     */
    for (const std::size_t n : {63, 127, 257, 1025}) {
      for (const bool fragmented : {false, true}) {
        auto expected = data<Tree>(n);
        if (fragmented) {
          for (std::size_t i = 1; i <= 3; ++i) {
            rotate(expected, i, n - i, n / (i + 2));
          }
        }
        for (std::size_t p = 0; p <= n; ++p) {
          auto tree = build<Tree>(data<Tree>(n));
          if (fragmented) {
            for (std::size_t i = 1; i <= 3; ++i) {
              tree.rotate_left(i, n - i, n / (i + 2));
            }
          }
          Tree::test_reset_counters();
          const auto budget = 3 * tree.height() + 1;
          Tree::test_fail_after(budget);
          Tree right;
          EXPECT_NO_THROW(right = tree.split_off(p));
          Tree::test_fail_after(-1);
          EXPECT_LE(Tree::test_counters.allocations, budget);
          check(tree, {expected.begin(), expected.begin() + p});
          check(right, {expected.begin() + p, expected.end()});
        }
      }
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, LeafParentSplitPacksEachSideOnce) {
    using Tree = TypeParam;
    constexpr auto capacity = Tree::block_capacity;
    const auto expected = data<Tree>(4 * capacity);
    for (const auto p : {2 * capacity, 2 * capacity + 1}) {
      auto tree = build<Tree>(expected);
      ASSERT_EQ(tree.height(), 1u);
      Tree::test_reset_counters();
      auto right = tree.split_off(p);
      // Four incoming child links and at most five outgoing links after the
      // leaf cut: no intermediate sibling group is packed then unpacked.
      EXPECT_LE(Tree::test_counters.child_transfers, 9u);
      check(tree, {expected.begin(), expected.begin() + p});
      check(right, {expected.begin() + p, expected.end()});
      tree.merge(right);
      check(tree, expected);
      check(right, {});
    }
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, GroupedSeamWorkAndOffPathIdentity) {
    using Tree = TypeParam;
    constexpr auto capacity = Tree::block_capacity;
    const std::array<std::array<std::size_t, 2>, 4> shapes{
        {{4096, 4096}, {4096, 1}, {1, 4096}, {4096, 4096}}};
    for (std::size_t shape = 0; shape < shapes.size(); ++shape) {
      SCOPED_TRACE(shape);
      auto left_values = data<Tree>(shapes[shape][0] * capacity, 37);
      auto right_values = data<Tree>(shapes[shape][1] * capacity, 41);
      auto tree = build<Tree>(left_values);
      auto donor = build<Tree>(right_values);
      {
        const auto keep =
            tree.size() - capacity + (shape == 3 ? capacity / 2 - 1 : 1);
        auto discarded = tree.split_off(keep);
        left_values.resize(keep);
        if (shape != 3) {
          auto rest = donor.split_off(capacity - 1);
          donor = std::move(rest);
          right_values.erase(right_values.begin(),
                             right_values.begin() + capacity - 1);
        }
      }
      // Shape 1 appends into the exterior leaf without rebuilding its spine.
      // Shape 0 borrows one neighbor and compacts to two leaves; shape 2
      // compacts to one exterior leaf. Shape 3 repairs just the adjacent pair.
      // A seam root below F/2 occupancy must not escape as a nonroot node.
      const auto h = std::max(tree.height(), donor.height());
      const auto budget = shape == 1 ? 0 : 2 * h + 4;
      auto leaves = tree.test_leaf_identities();
      const auto donor_leaves = donor.test_leaf_identities();
      leaves.insert(leaves.end(), donor_leaves.begin(), donor_leaves.end());
      auto nodes = tree.test_internal_node_identities();
      const auto donor_nodes = donor.test_internal_node_identities();
      nodes.insert(nodes.end(), donor_nodes.begin(), donor_nodes.end());
      Tree::test_reset_counters();
      IntegerBlock::copied_elements = 0;
      const auto bytes = Tree::test_counters.live_bytes;
      Tree::test_fail_after(budget);
      tree.merge(donor);
      Tree::test_fail_after(-1);
      const auto counters = Tree::test_counters;
      EXPECT_EQ(counters.allocations, budget);
      EXPECT_EQ(counters.peak_bytes - bytes, budget * Tree::node_storage_bytes);
      // Work bounds for these full-spine fixtures, not empirical reserve
      // bounds. Already balanced joins may preserve nodes without visiting
      // their children. Reintroducing one join per seam leaf adds repeated
      // taller-spine visits.
      if (shape == 0) {
        EXPECT_LE(counters.node_visits, 18 * h - 3);
      } else if (shape == 1) {
        EXPECT_EQ(counters.node_visits, h);
      } else if (shape == 2) {
        EXPECT_LE(counters.node_visits, 6 * h - 3);
      }
      EXPECT_LE(counters.node_visits, 24 * h);
      EXPECT_LE(counters.payload_mutations, 4);
      EXPECT_LE(IntegerBlock::copied_elements, 8 * capacity);
      EXPECT_GE(common_identity_count(leaves, tree.test_leaf_identities()) + 4,
                leaves.size());
      EXPECT_GE(
          common_identity_count(nodes, tree.test_internal_node_identities()) +
              counters.node_visits,
          nodes.size());
      left_values.insert(left_values.end(), right_values.begin(),
                         right_values.end());
      check(tree, left_values);
      check(donor, {});
    }
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }

  TYPED_TEST(SequenceTreeRecyclingSpec, DeepFragmentedRandomOwnership) {
    using Tree = TypeParam;
    auto expected = data<Tree>(4096 * Tree::block_capacity + 3);
    auto tree = build<Tree>(expected);
    std::mt19937_64 random(716381);
    for (std::size_t step = 0; step < 80; ++step) {
      SCOPED_TRACE(step);
      const auto p = random() % (tree.size() + 1);
      auto tail = tree.split_off(p);
      check(tree, {expected.begin(), expected.begin() + p});
      check(tail, {expected.begin() + p, expected.end()});
      // Reverse the pieces and repeatedly relocate underfull exterior leaves.
      tail.merge(tree);
      EXPECT_TRUE(tree.empty());
      tree = std::move(tail);
      rotate(expected, 0, expected.size(), p);
      const auto left = random() % tree.size();
      const auto right = left + random() % (tree.size() - left + 1);
      const auto distance = random();
      tree.rotate_left(left, right, distance);
      rotate(expected, left, right, distance);
      check(tree, expected);
      const auto memory = tree.memory_usage();
      EXPECT_EQ(Tree::test_counters.live_allocations,
                memory.blocks + memory.nodes);
      EXPECT_EQ(Tree::test_counters.live_bytes,
                memory.block_bytes + memory.node_bytes);
    }
  }

  template <LengthLayout Layout>
  void locality() {
    using Tree = SequenceTree<IntegerBlock, 4, Layout>;
    for (const std::size_t leaves : {64, 1024, 16384}) {
      auto expected = data<Tree>(leaves * IntegerBlock::capacity);
      auto tree = build<Tree>(expected);
      const auto identities = tree.test_leaf_identities();
      const auto original_nodes = tree.test_internal_node_identities();
      EXPECT_EQ(original_nodes.size(), tree.node_count());
      const auto h = tree.height();
      auto bound = [&](std::size_t mutations) {
        EXPECT_LE(Tree::test_counters.node_visits, 160 * (h + 4));
        EXPECT_LE(Tree::test_counters.child_transfers, 160 * 4 * (h + 4));
        EXPECT_LE(Tree::test_counters.payload_mutations, mutations);
        EXPECT_LE(IntegerBlock::copied_elements,
                  mutations * 2 * IntegerBlock::capacity);
      };
      Tree::test_reset_counters();
      IntegerBlock::copied_elements = 0;
      auto tail = tree.split_off(tree.size() / 2 + 1);
      bound(1);
      EXPECT_LE(Tree::test_counters.allocations, 3 * h + 1);
      EXPECT_TRUE(tree.test_validate());
      EXPECT_TRUE(tail.test_validate());
      auto split_nodes = tree.test_internal_node_identities();
      const auto tail_nodes = tail.test_internal_node_identities();
      split_nodes.insert(split_nodes.end(), tail_nodes.begin(),
                         tail_nodes.end());
      EXPECT_EQ(split_nodes.size(), tree.node_count() + tail.node_count());
      const auto split_survivors =
          common_identity_count(original_nodes, split_nodes);
      EXPECT_GE(split_survivors + 32 * (h + 1), original_nodes.size());
      if (leaves == 16384) {
        EXPECT_GT(split_survivors, original_nodes.size() / 2);
      }
      Tree::test_reset_counters();
      IntegerBlock::copied_elements = 0;
      tree.merge(tail);
      bound(5);
      EXPECT_LE(Tree::test_counters.allocations, 2 * h + 4);
      check(tree, expected);
      const auto before_rotation_nodes = tree.test_internal_node_identities();
      Tree::test_reset_counters();
      IntegerBlock::copied_elements = 0;
      const auto rotation_height = tree.height();
      tree.rotate_left(1, tree.size() - 1, tree.size() / 3);
      bound(18);
      EXPECT_LE(Tree::test_counters.allocations, 15 * rotation_height + 24 + 3);
      rotate(expected, 1, expected.size() - 1, expected.size() / 3);
      check(tree, expected);
      const auto after_rotation_nodes = tree.test_internal_node_identities();
      EXPECT_EQ(after_rotation_nodes.size(), tree.node_count());
      const auto rotation_survivors =
          common_identity_count(before_rotation_nodes, after_rotation_nodes);
      EXPECT_GE(rotation_survivors + 128 * (rotation_height + 4),
                before_rotation_nodes.size());
      if (leaves == 16384) {
        EXPECT_GT(rotation_survivors, before_rotation_nodes.size() / 2);
      }
      EXPECT_GE(common_identity_count(identities, tree.test_leaf_identities()),
                leaves - 16);
    }
  }
  TEST(SequenceTreeLocality, CumulativeSpinesNotLeafScans) {
    locality<LengthLayout::cumulative>();
  }
  TEST(SequenceTreeLocality, IndividualSpinesNotLeafScans) {
    locality<LengthLayout::individual>();
  }

  // Uniform conceptual elements allow full-width counts without huge
  // allocations.
  struct WeightedBlock {
    using value_type = std::uint32_t;
    static constexpr std::size_t capacity =
        (std::numeric_limits<std::size_t>::max() / 4) & ~std::size_t{1};
    std::size_t n = 0;
    std::size_t size() const noexcept { return n; }
    value_type operator[](std::size_t i) const {
      assert(i < n);
      (void)i;
      return 7;
    }
    void rotate_left(std::size_t l, std::size_t r, std::size_t) noexcept {
      assert(l <= r && r <= n);
      (void)l;
      (void)r;
    }
    void redistribute(WeightedBlock& rhs, std::size_t left) noexcept {
      assert(this != &rhs && left <= capacity && left <= n + rhs.n &&
             n + rhs.n - left <= capacity);
      rhs.n = n + rhs.n - left;
      n = left;
    }
  };
  template <LengthLayout Layout>
  void full_width_counts() {
    using Tree = SequenceTree<WeightedBlock, 4, Layout>;
    constexpr auto maximum = std::numeric_limits<std::size_t>::max();
    {
      std::array<WeightedBlock, 4> full{{{WeightedBlock::capacity},
                                         {WeightedBlock::capacity},
                                         {WeightedBlock::capacity},
                                         {WeightedBlock::capacity - 1}}};
      auto aligned = Tree::from_blocks(full);
      const auto ids = aligned.test_leaf_identities();
      Tree::test_reset_counters();
      Tree::test_fail_after(0);
      EXPECT_NO_THROW(
          aligned.rotate_left(0, aligned.size(), 2 * WeightedBlock::capacity));
      Tree::test_fail_after(-1);
      EXPECT_EQ(Tree::test_counters.allocations, 0);
      EXPECT_TRUE(aligned.test_validate());
      auto expected = ids;
      std::rotate(expected.begin(), expected.begin() + 2, expected.end());
      EXPECT_EQ(aligned.test_leaf_identities(), expected);
    }
    std::array<WeightedBlock, 5> blocks{
        {{WeightedBlock::capacity},
         {WeightedBlock::capacity},
         {WeightedBlock::capacity},
         {WeightedBlock::capacity},
         {maximum - 4 * WeightedBlock::capacity}}};
    auto tree = Tree::from_blocks(blocks);
    ASSERT_EQ(tree.size(), maximum);
    ASSERT_TRUE(tree.test_validate());
    for (const auto p : {std::size_t{1}, WeightedBlock::capacity, maximum / 2,
                         maximum / 2 + 1, maximum - 1}) {
      EXPECT_EQ(tree[p], 7);
      auto suffix = tree.split_off(p);
      EXPECT_EQ(tree.size(), p);
      EXPECT_EQ(suffix.size(), maximum - p);
      EXPECT_TRUE(tree.test_validate());
      EXPECT_TRUE(suffix.test_validate());
      tree.merge(suffix);
      EXPECT_EQ(tree.size(), maximum);
      EXPECT_TRUE(tree.test_validate());
    }
    std::array<WeightedBlock, 1> extra{{{1}}};
    auto donor = Tree::from_blocks(extra);
    Tree::test_fail_after(0);
    EXPECT_THROW(tree.merge(donor), std::length_error);
    EXPECT_THROW(tree.rotate_left(1, 0, 0), std::out_of_range);
    Tree::test_fail_after(-1);
    EXPECT_EQ(tree.size(), maximum);
    EXPECT_EQ(donor.size(), 1);
    tree.rotate_left(1, maximum - 1, maximum / 2);
    EXPECT_TRUE(tree.test_validate());
    EXPECT_EQ(tree.size(), maximum);
    std::array<WeightedBlock, 6> overflow{{{WeightedBlock::capacity},
                                           {WeightedBlock::capacity},
                                           {WeightedBlock::capacity},
                                           {WeightedBlock::capacity},
                                           {7},
                                           {1}}};
    EXPECT_THROW(Tree::from_blocks(overflow), std::length_error);
  }
  TEST(SequenceTreeCounts, CumulativeUnsignedHighBitAndOverflow) {
    full_width_counts<LengthLayout::cumulative>();
  }
  TEST(SequenceTreeCounts, IndividualUnsignedHighBitAndOverflow) {
    full_width_counts<LengthLayout::individual>();
  }

}  // namespace
