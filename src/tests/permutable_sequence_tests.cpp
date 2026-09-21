#include <gtest/gtest.h>
#include <pixie/permutations/sequence.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using namespace pixie;
using namespace pixie::detail::sequence;

template <class Sequence>
struct ResetFailures {
  ~ResetFailures() {
    Sequence::test_fail_payload_after(-1);
    Sequence::test_tree_type::test_fail_after(-1);
  }
};

template <class T>
std::vector<T> values(std::size_t size, std::uint64_t seed = 42) {
  std::mt19937_64 random(seed);
  std::vector<T> result;
  result.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    if constexpr (std::same_as<T, std::string>) {
      result.push_back(std::to_string(random()) + std::string(40, 'x'));
    } else {
      result.push_back(static_cast<T>(random()));
    }
  }
  return result;
}

template <class T>
void rotate(std::vector<T>& input,
            std::size_t left,
            std::size_t right,
            std::size_t distance) {
  if (left != right) {
    std::rotate(input.begin() + left,
                input.begin() + left + distance % (right - left),
                input.begin() + right);
  }
}

template <class Sequence>
void check(const Sequence& sequence,
           const std::vector<typename Sequence::value_type>& expected) {
  ASSERT_TRUE(sequence.test_tree().test_validate());
  ASSERT_EQ(sequence.size(), expected.size());
  ASSERT_EQ(sequence.empty(), expected.empty());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    ASSERT_EQ(sequence[i], expected[i]) << "index=" << i;
  }
  const auto memory = sequence.memory_usage();
  EXPECT_EQ(memory.facade_bytes, sizeof(sequence));
  EXPECT_EQ(memory.total_bytes, sizeof(sequence) + memory.tree_block_bytes +
                                    memory.tree_node_bytes +
                                    memory.chunk_header_bytes +
                                    memory.vector_capacity_bytes);
  EXPECT_EQ(memory.vector_capacity_bytes,
            memory.vector_live_bytes + memory.vector_slack_bytes);
  EXPECT_EQ(sequence.memory_usage_bytes(), memory.total_bytes);
  if constexpr (Sequence::storage == ElementStorage::packed) {
    EXPECT_EQ(memory.order_bytes, 0);
    EXPECT_EQ(memory.chunks, 0);
    EXPECT_EQ(memory.vector_capacity_bytes, 0);
    EXPECT_GE(memory.packed_capacity_bits,
              sequence.size() *
                  std::numeric_limits<typename Sequence::value_type>::digits);
  } else {
    EXPECT_EQ(memory.packed_capacity_bits, 0);
    EXPECT_EQ(memory.order_bytes,
              memory.tree_block_bytes + memory.tree_node_bytes);
    EXPECT_EQ(memory.vector_live_bytes,
              sequence.size() * sizeof(typename Sequence::value_type));
  }
  if (sequence.empty()) {
    EXPECT_EQ(memory.total_bytes, sizeof(sequence));
    EXPECT_EQ(memory.chunks, 0);
  }
}

template <class Sequence>
class PermutableSequenceSpec : public ::testing::Test {};
using Sequences = ::testing::Types<
    PermutableSequence<bool>,
    PermutableSequence<bool,
                       ElementStorage::indirect,
                       512,
                       4,
                       LengthLayout::individual,
                       31>,
    PermutableSequence<std::uint8_t, ElementStorage::automatic, 512, 4>,
    PermutableSequence<std::uint16_t,
                       ElementStorage::packed,
                       1024,
                       8,
                       LengthLayout::individual>,
    PermutableSequence<std::uint32_t, ElementStorage::automatic, 2048, 16>,
    PermutableSequence<std::uint64_t,
                       ElementStorage::packed,
                       512,
                       16,
                       LengthLayout::individual>,
    PermutableSequence<std::uint32_t, ElementStorage::indirect, 512, 8>,
    PermutableSequence<std::int64_t,
                       ElementStorage::automatic,
                       512,
                       16,
                       LengthLayout::individual,
                       73>,
    PermutableSequence<std::string,
                       ElementStorage::automatic,
                       1024,
                       4,
                       LengthLayout::cumulative,
                       127>>;
TYPED_TEST_SUITE(PermutableSequenceSpec, Sequences);

template <class S, class Range>
concept CanConsume =
    requires(Range&& range) { S::from_range(std::forward<Range>(range)); };

TYPED_TEST(PermutableSequenceSpec, PublicContractAliasesNoexceptAndLayout) {
  using S = TypeParam;
  using T = typename S::value_type;
  using R = typename S::const_reference;
  using Base = PermutableSequenceBase<S, T, R>;
  static_assert(std::is_base_of_v<Base, S> && std::is_empty_v<Base>);
  static_assert(std::same_as<typename S::size_type, std::size_t>);
  static_assert(std::same_as<decltype(std::declval<const Base&>()[0]), R>);
  static_assert(noexcept(std::declval<const Base&>().size()));
  static_assert(noexcept(std::declval<const Base&>().empty()));
  static_assert(noexcept(std::declval<const Base&>().memory_usage_bytes()));
  static_assert(!noexcept(std::declval<const Base&>()[0]));
  static_assert(!CanConsume<S, int>);
  static_assert(!CanConsume<S, std::array<std::array<int, 2>, 2>&>);
  static_assert(
      sizeof(S) ==
      sizeof(typename S::test_tree_type) +
          (S::storage == ElementStorage::packed ? 0 : 2 * sizeof(void*)));
  auto input = values<T>(3);
  static_assert(std::same_as<decltype(Base::from_range(input)), S>);
  auto owner = Base::from_range(input);
  Base& contract = owner;
  EXPECT_EQ(contract.size(), 3);
  EXPECT_THROW(contract[3], std::out_of_range);
  EXPECT_THROW(contract.rotate_left(0, 4, 0), std::out_of_range);
  EXPECT_EQ(contract.memory_usage_bytes(), owner.memory_usage().total_bytes);
  if constexpr (std::is_reference_v<R>) {
    const T* first = &contract[0];
    S receiver;
    receiver.merge(owner);
    EXPECT_EQ(&receiver[0], first);
  }
}

TYPED_TEST(PermutableSequenceSpec, EmptyCheckedRangesAndExclusiveMoves) {
  using Sequence = TypeParam;
  using T = typename Sequence::value_type;
  static_assert(!std::is_copy_constructible_v<Sequence>);
  static_assert(!std::is_copy_assignable_v<Sequence>);
  static_assert(std::is_nothrow_move_constructible_v<Sequence>);
  static_assert(std::is_nothrow_move_assignable_v<Sequence>);
  if constexpr (Sequence::storage == ElementStorage::packed) {
    static_assert(
        std::same_as<decltype(std::declval<const Sequence&>()[0]), T>);
  } else {
    static_assert(
        std::same_as<decltype(std::declval<const Sequence&>()[0]), const T&>);
  }
  Sequence sequence;
  check(sequence, {});
  EXPECT_THROW(sequence[0], std::out_of_range);
  EXPECT_THROW(sequence.rotate_left(1, 0, 0), std::out_of_range);
  EXPECT_THROW(sequence.rotate_left(0, 1, 0), std::out_of_range);
  EXPECT_NO_THROW(sequence.rotate_left(0, 0, SIZE_MAX));
  auto input = values<T>(117);
  const auto expected = input;
  sequence = Sequence::from_range(input);
  check(sequence, expected);
  EXPECT_THROW(sequence[sequence.size()], std::out_of_range);
  EXPECT_THROW(sequence.rotate_left(10, 9, 0), std::out_of_range);
  EXPECT_THROW(sequence.rotate_left(0, sequence.size() + 1, 0),
               std::out_of_range);
  sequence.merge(sequence);
  auto* self = &sequence;
  sequence = std::move(*self);
  check(sequence, expected);
  Sequence moved(std::move(sequence));
  check(sequence, {});
  sequence.merge(moved);
  check(moved, {});
  check(sequence, expected);
  sequence.merge(moved);
  check(sequence, expected);
  auto replacement_input = values<T>(17, 71);
  const auto replacement_expected = replacement_input;
  moved = Sequence::from_range(replacement_input);
  sequence = std::move(moved);
  check(moved, {});
  check(sequence, replacement_expected);
}

TYPED_TEST(PermutableSequenceSpec, DifferentialRotationsAndUnrebasedMerges) {
  using Sequence = TypeParam;
  using T = typename Sequence::value_type;
  auto input = values<T>(2309);
  auto expected = input;
  auto sequence = Sequence::from_range(input);
  std::mt19937_64 random(174);
  for (std::size_t step = 0; step < 80; ++step) {
    if (step % 13 == 0) {
      auto addition = values<T>(step + 1, step);
      expected.insert(expected.end(), addition.begin(), addition.end());
      auto donor = Sequence::from_range(addition);
      sequence.merge(donor);
      check(donor, {});
    }
    auto left = random() % (expected.size() + 1);
    auto right = random() % (expected.size() + 1);
    if (left > right) {
      std::swap(left, right);
    }
    if (step % 7 == 0) {
      left = 0;
      right = expected.size();
    }
    const auto distance = random();
    sequence.rotate_left(left, right, distance);
    rotate(expected, left, right, distance);
    check(sequence, expected);
  }
}

TYPED_TEST(PermutableSequenceSpec, EveryShortRotationAndBoundaryConstruction) {
  using Sequence = TypeParam;
  using T = typename Sequence::value_type;
  for (std::size_t size = 0; size <= 7; ++size) {
    const auto original = values<T>(size);
    for (std::size_t left = 0; left <= size; ++left) {
      for (std::size_t right = left; right <= size; ++right) {
        for (std::size_t distance = 0; distance <= size + 1; ++distance) {
          auto input = original;
          auto expected = original;
          auto sequence = Sequence::from_range(input);
          sequence.rotate_left(left, right, distance);
          rotate(expected, left, right, distance);
          check(sequence, expected);
        }
      }
    }
  }
  constexpr auto capacity = Sequence::test_tree_type::block_capacity;
  for (auto size : {capacity - 1, capacity, capacity + 1, 2 * capacity + 1}) {
    auto input = values<T>(size);
    auto expected = input;
    auto sequence = Sequence::from_range(input);
    check(sequence, expected);
  }
}

TEST(PermutableSequence, CompileTimeSelectionAndExactPointerBudget) {
  static_assert(PermutableSequence<bool>::storage == ElementStorage::packed);
  static_assert(PermutableSequence<unsigned>::storage ==
                ElementStorage::packed);
  static_assert(PermutableSequence<int>::storage == ElementStorage::indirect);
  static_assert(PermutableSequence<double>::storage ==
                ElementStorage::indirect);
  static_assert(PermutableSequence<std::string>::storage ==
                ElementStorage::indirect);
  static_assert(SequenceBlock<pixie::detail::sequence::PointerBlock<int, 512>>);
  static_assert(sizeof(pixie::detail::sequence::PointerBlock<int, 512>) == 64);
  static_assert(sizeof(pixie::detail::sequence::PointerBlock<bool, 2048>) ==
                256);
  static_assert(pixie::detail::sequence::PointerBlock<int, 512>::capacity ==
                (64 - 2 * sizeof(std::size_t)) / sizeof(const int*));
  EXPECT_EQ(sizeof(PermutableSequence<unsigned>),
            sizeof(PermutableSequence<unsigned>::test_tree_type));
}

TEST(PermutableSequence, NativePointerBlockWrappedRedistribution) {
  using Block = pixie::detail::sequence::PointerBlock<int, 512>;
  constexpr auto capacity = Block::capacity;
  std::array<int, 2 * capacity + 1> objects{};
  std::array<const int*, 2 * capacity + 1> pointers{};
  for (std::size_t i = 0; i < pointers.size(); ++i) {
    pointers[i] = &objects[i];
  }
  EXPECT_THROW(
      (Block(std::span<const int* const>(pointers.data(), capacity + 1))),
      std::length_error);
  for (std::size_t a = 0; a <= capacity; ++a) {
    for (std::size_t b = 0; b <= capacity; ++b) {
      for (std::size_t split = 0; split <= a + b; ++split) {
        if (split > capacity || a + b - split > capacity) {
          continue;
        }
        Block left(std::span<const int* const>(pointers.data(), a));
        Block right(std::span<const int* const>(pointers.data() + a, b));
        std::vector<const int*> expected(pointers.begin(),
                                         pointers.begin() + a + b);
        left.rotate_left(0, a, 5);
        right.rotate_left(0, b, 3);
        rotate(expected, 0, a, 5);
        rotate(expected, a, a + b, 3);
        if (a > 2) {
          left.rotate_left(1, a, 1);
          rotate(expected, 1, a, 1);
        }
        left.redistribute(right, split);
        ASSERT_EQ(left.size(), split);
        ASSERT_EQ(right.size(), a + b - split);
        for (std::size_t i = 0; i < split; ++i) {
          EXPECT_EQ(left[i], expected[i]);
        }
        for (std::size_t i = split; i < a + b; ++i) {
          EXPECT_EQ(right[i - split], expected[i]);
        }
      }
    }
  }
}

TEST(PermutableSequence, SinglePassStreamInput) {
  std::istringstream stream("1 2 3 4 5 6 7 8 9");
  auto input = std::ranges::istream_view<unsigned>(stream);
  auto sequence = PermutableSequence<unsigned>::from_range(input);
  check(sequence, {1, 2, 3, 4, 5, 6, 7, 8, 9});
}

struct MovingInput {
  unsigned length = 23;
  unsigned position = 0;
  unsigned moves = 0;
  unsigned throw_increment = std::numeric_limits<unsigned>::max();
  struct Iterator {
    using value_type = unsigned;
    using difference_type = std::ptrdiff_t;
    using iterator_concept = std::input_iterator_tag;
    MovingInput* input;
    unsigned operator*() const { return input->position; }
    Iterator& operator++() {
      if (input->position == input->throw_increment) {
        throw std::runtime_error("input increment");
      }
      ++input->position;
      return *this;
    }
    void operator++(int) { ++*this; }
    bool operator==(std::default_sentinel_t) const {
      return input->position == input->length;
    }
    friend unsigned iter_move(const Iterator& it) {
      ++it.input->moves;
      return it.input->position + 100;
    }
  };
  Iterator begin() { return {this}; }
  std::default_sentinel_t end() const { return {}; }
};
static_assert(std::ranges::input_range<MovingInput>);
static_assert(!std::ranges::forward_range<MovingInput>);

TEST(PermutableSequence, HonorsCustomIteratorMoveAndCleansUpThrowingIncrement) {
  using Packed = PermutableSequence<std::uint64_t, ElementStorage::packed, 512>;
  using Indirect = PermutableSequence<unsigned, ElementStorage::indirect, 512,
                                      8, LengthLayout::cumulative, 16>;
  auto exercise = []<class Sequence>() {
    {
      MovingInput input;
      auto sequence = Sequence::from_range(input);
      ASSERT_EQ(sequence.size(), input.length);
      EXPECT_EQ(input.moves, input.length);
      for (unsigned i = 0; i < input.length; ++i) {
        EXPECT_EQ(sequence[i], i + 100);
      }
    }
    for (unsigned fail = 0; fail < 23; ++fail) {
      MovingInput input;
      input.throw_increment = fail;
      EXPECT_THROW(Sequence::from_range(input), std::runtime_error);
      EXPECT_EQ(input.moves, fail + 1);
      EXPECT_EQ(Sequence::test_tree_type::test_counters.live_allocations, 0);
    }
  };
  exercise.template operator()<Packed>();
  exercise.template operator()<Indirect>();
}

// No default construction, copying, or move assignment is available. Slot
// ownership therefore cannot silently depend on vector resizing or T swaps.
struct alignas(128) Tracked {
  inline static int live = 0;
  inline static int moves = 0;
  int value;
  explicit Tracked(int n) : value(n) { ++live; }
  Tracked(const Tracked&) = delete;
  Tracked& operator=(const Tracked&) = delete;
  Tracked(Tracked&& other) noexcept : value(std::exchange(other.value, -1)) {
    ++live;
    ++moves;
  }
  Tracked& operator=(Tracked&&) = delete;
  ~Tracked() noexcept { --live; }
};
static_assert(!std::is_default_constructible_v<Tracked>);
static_assert(!std::is_move_assignable_v<Tracked>);
using Indirect = PermutableSequence<Tracked,
                                    ElementStorage::automatic,
                                    512,
                                    4,
                                    LengthLayout::individual,
                                    512>;

Indirect tracked_sequence(int begin, int end) {
  auto input = std::views::iota(begin, end) |
               std::views::transform([](int i) { return Tracked(i); });
  return Indirect::from_range(input);
}

struct ThrowingTrackedInput {
  int length = 41;
  int position = 0;
  int comparisons = 0;
  int throw_increment = -1;
  int throw_comparison = -1;
  int live_at_failure = -1;
  std::size_t tree_allocations_at_failure = 0;

  [[noreturn]] void fail() {
    live_at_failure = Tracked::live;
    tree_allocations_at_failure =
        Indirect::test_tree_type::test_counters.live_allocations;
    throw std::runtime_error("tracked input failure");
  }
  struct Iterator {
    using value_type = Tracked;
    using difference_type = std::ptrdiff_t;
    using iterator_concept = std::input_iterator_tag;
    ThrowingTrackedInput* input;
    Tracked operator*() const { return Tracked(input->position); }
    Iterator& operator++() {
      if (input->position == input->throw_increment) {
        input->fail();
      }
      ++input->position;
      return *this;
    }
    void operator++(int) { ++*this; }
    bool operator==(std::default_sentinel_t) const {
      if (input->comparisons++ == input->throw_comparison) {
        input->fail();
      }
      return input->position == input->length;
    }
  };
  Iterator begin() { return {this}; }
  std::default_sentinel_t end() const { return {}; }
};
static_assert(std::ranges::input_range<ThrowingTrackedInput>);
static_assert(!std::ranges::forward_range<ThrowingTrackedInput>);

TEST(PermutableSequence,
     TrackedIncrementAndSentinelFailuresReclaimPartialChunks) {
  using Tree = Indirect::test_tree_type;
  ASSERT_EQ(Tracked::live, 0);
  const auto initial_allocations = Tree::test_counters.live_allocations;
  {
    auto existing = tracked_sequence(100, 107);
    const auto baseline_live = Tracked::live;
    const auto baseline_allocations = Tree::test_counters.live_allocations;
    const auto baseline_bytes = Tree::test_counters.live_bytes;
    ThrowingTrackedInput complete;
    {
      auto sequence = Indirect::from_range(complete);
      ASSERT_EQ(sequence.size(), complete.length);
      for (int i = 0; i < complete.length; ++i) {
        EXPECT_EQ(sequence[i].value, i);
      }
    }
    auto check_cleanup = [&] {
      EXPECT_EQ(Tracked::live, baseline_live);
      EXPECT_EQ(Tree::test_counters.live_allocations, baseline_allocations);
      EXPECT_EQ(Tree::test_counters.live_bytes, baseline_bytes);
      for (std::size_t i = 0; i < existing.size(); ++i) {
        EXPECT_EQ(existing[i].value, 100 + static_cast<int>(i));
      }
    };
    check_cleanup();
    bool increment_after_published_leaf = false;
    for (int failure = 0; failure < complete.length; ++failure) {
      SCOPED_TRACE(failure);
      ThrowingTrackedInput input;
      input.throw_increment = failure;
      EXPECT_THROW(Indirect::from_range(input), std::runtime_error);
      // ++ runs after emplacement and destruction of the iterator temporary,
      // so this includes the newly constructed, still unpublished slot.
      EXPECT_EQ(input.live_at_failure, baseline_live + failure + 1);
      increment_after_published_leaf |=
          input.tree_allocations_at_failure > baseline_allocations;
      check_cleanup();
    }
    EXPECT_TRUE(increment_after_published_leaf);
    bool comparison_with_partial_first_chunk = false;
    bool comparison_after_published_leaf = false;
    for (int failure = 0; failure < complete.comparisons; ++failure) {
      SCOPED_TRACE(failure);
      ThrowingTrackedInput input;
      input.throw_comparison = failure;
      EXPECT_THROW(Indirect::from_range(input), std::runtime_error);
      EXPECT_GE(input.live_at_failure, baseline_live);
      comparison_with_partial_first_chunk |=
          input.live_at_failure > baseline_live &&
          input.live_at_failure < baseline_live + 4;
      comparison_after_published_leaf |=
          input.tree_allocations_at_failure > baseline_allocations;
      check_cleanup();
    }
    EXPECT_TRUE(comparison_with_partial_first_chunk);
    EXPECT_TRUE(comparison_after_published_leaf);
  }
  EXPECT_EQ(Tracked::live, 0);
  EXPECT_EQ(Tree::test_counters.live_allocations, initial_allocations);
}

TEST(PermutableSequence, ConsumesExistingMoveOnlyObjectsExactlyOnce) {
  ASSERT_EQ(Tracked::live, 0);
  {
    std::vector<Tracked> input;
    input.reserve(17);
    for (int i = 0; i < 17; ++i) {
      input.emplace_back(i);
    }
    Tracked::moves = 0;
    auto sequence = Indirect::from_range(input);
    EXPECT_EQ(Tracked::moves, 17);
    EXPECT_EQ(Tracked::live, 34);
    for (int i = 0; i < 17; ++i) {
      EXPECT_EQ(input[i].value, -1);
      EXPECT_EQ(sequence[i].value, i);
      EXPECT_NE(&sequence[i], &input[i]);
    }
  }
  EXPECT_EQ(Tracked::live, 0);
}

template <class Sequence>
auto addresses(const Sequence& sequence) {
  std::vector<const typename Sequence::value_type*> result;
  for (std::size_t i = 0; i < sequence.size(); ++i) {
    result.push_back(std::addressof(sequence[i]));
  }
  return result;
}

template <class Sequence>
void check_addresses(
    const Sequence& sequence,
    const std::vector<const typename Sequence::value_type*>& expected) {
  ASSERT_EQ(sequence.size(), expected.size());
  EXPECT_TRUE(sequence.test_tree().test_validate());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    EXPECT_EQ(std::addressof(sequence[i]), expected[i]);
  }
}

TEST(PermutableSequence, NonemptyMoveAssignmentDestroysOnlyReplacedPayload) {
  using Tree = Indirect::test_tree_type;
  ASSERT_EQ(Tracked::live, 0);
  const auto initial_allocations = Tree::test_counters.live_allocations;
  {
    auto receiver = tracked_sequence(0, 17);
    auto source = tracked_sequence(100, 129);
    source.rotate_left(0, source.size(), 7);
    const auto incoming_addresses = addresses(source);
    std::vector<int> incoming_values;
    for (auto* value : incoming_addresses) {
      incoming_values.push_back(value->value);
    }
    const auto old_count = receiver.size();
    const auto old_memory = receiver.test_tree().memory_usage();
    const auto before_live = Tracked::live;
    const auto before_allocations = Tree::test_counters.live_allocations;
    Tracked::moves = 0;
    ResetFailures<Indirect> reset;
    Tree::test_fail_after(0);
    Indirect::test_fail_payload_after(0);
    EXPECT_EQ(&(receiver = std::move(source)), &receiver);
    check_addresses(receiver, incoming_addresses);
    for (std::size_t i = 0; i < incoming_values.size(); ++i) {
      EXPECT_EQ(receiver[i].value, incoming_values[i]);
    }
    EXPECT_EQ(Tracked::moves, 0);
    EXPECT_EQ(Tracked::live, before_live - static_cast<int>(old_count));
    EXPECT_EQ(Tree::test_counters.live_allocations,
              before_allocations - old_memory.blocks - old_memory.nodes);
    EXPECT_TRUE(source.empty());
    EXPECT_EQ(source.size(), 0);
    EXPECT_EQ(source.memory_usage_bytes(), sizeof(source));
    EXPECT_EQ(source.memory_usage().chunks, 0);
    EXPECT_TRUE(source.test_tree().test_validate());
  }
  EXPECT_EQ(Tracked::live, 0);
  EXPECT_EQ(Tree::test_counters.live_allocations, initial_allocations);
}

TEST(PermutableSequence,
     StableAlignedMoveOnlyPayloadNeverMovesAfterConstruction) {
  ASSERT_EQ(Tracked::live, 0);
  {
    auto sequence = tracked_sequence(0, 173);
    auto donor = tracked_sequence(173, 254);
    auto expected = addresses(sequence);
    auto donor_addresses = addresses(donor);
    Tracked::moves = 0;
    for (auto* address : expected) {
      EXPECT_EQ(reinterpret_cast<std::uintptr_t>(address) % 128, 0);
    }
    sequence.rotate_left(3, 171, 37);
    rotate(expected, 3, 171, 37);
    donor.rotate_left(0, donor.size(), 23);
    rotate(donor_addresses, 0, donor_addresses.size(), 23);
    sequence.merge(donor);
    expected.insert(expected.end(), donor_addresses.begin(),
                    donor_addresses.end());
    EXPECT_TRUE(donor.empty());
    EXPECT_EQ(donor.memory_usage_bytes(), sizeof(donor));
    auto moved = std::move(sequence);
    EXPECT_TRUE(sequence.empty());
    sequence = std::move(moved);
    EXPECT_TRUE(moved.empty());
    check_addresses(sequence, expected);
    EXPECT_EQ(Tracked::moves, 0);
    EXPECT_EQ(Tracked::live, 254);
    auto memory = sequence.memory_usage();
    EXPECT_EQ(memory.chunks, (173 + 3) / 4 + (81 + 3) / 4);
    EXPECT_EQ(memory.vector_live_bytes, 254 * sizeof(Tracked));
    EXPECT_EQ(memory.vector_capacity_bytes,
              memory.vector_live_bytes + memory.vector_slack_bytes);
    sequence = std::move(donor);
    EXPECT_TRUE(sequence.empty());
    EXPECT_EQ(Tracked::live, 0);
    EXPECT_EQ(Tracked::moves, 0);
  }
  EXPECT_EQ(Tracked::live, 0);
}

TEST(PermutableSequence, ForcedIndirectBoolHasRealStableObjectReferences) {
  using Sequence = PermutableSequence<bool, ElementStorage::indirect, 512, 8,
                                      LengthLayout::cumulative, 7>;
  std::array<bool, 23> input{};
  for (std::size_t i = 0; i < input.size(); ++i) {
    input[i] = i % 3 == 0;
  }
  auto sequence = Sequence::from_range(input);
  auto donor = Sequence::from_range(input);
  auto expected = addresses(sequence);
  const auto donor_addresses = addresses(donor);
  expected.insert(expected.end(), donor_addresses.begin(),
                  donor_addresses.end());
  Sequence::test_fail_payload_after(0);
  ResetFailures<Sequence> reset;
  sequence.merge(donor);
  sequence.rotate_left(1, 45, 13);
  rotate(expected, 1, 45, 13);
  check_addresses(sequence, expected);
  for (std::size_t i = 0; i < sequence.size(); ++i) {
    const bool& reference = sequence[i];
    EXPECT_EQ(reference, *expected[i]);
    for (std::size_t j = 0; j < i; ++j) {
      EXPECT_NE(expected[i], expected[j]);
    }
  }
  EXPECT_EQ(sequence.memory_usage().chunks, 8);
}

TEST(PermutableSequence,
     ThousandsOfSingletonMergesRetainChunksAndDestroyIteratively) {
  using Sequence = PermutableSequence<Tracked, ElementStorage::indirect, 512, 8,
                                      LengthLayout::cumulative, 1>;
  ASSERT_EQ(Tracked::live, 0);
  {
    Sequence sequence;
    std::vector<const Tracked*> expected;
    for (int i = 0; i < 10000; ++i) {
      auto input = std::views::iota(i, i + 1) |
                   std::views::transform([](int n) { return Tracked(n); });
      auto donor = Sequence::from_range(input);
      expected.push_back(&donor[0]);
      Tracked::moves = 0;
      Sequence::test_fail_payload_after(0);
      ResetFailures<Sequence> reset;
      sequence.merge(donor);
      EXPECT_EQ(Tracked::moves, 0);
      EXPECT_TRUE(donor.empty());
    }
    check_addresses(sequence, expected);
    const auto memory = sequence.memory_usage();
    EXPECT_EQ(memory.chunks, 10000);
    EXPECT_EQ(memory.vector_live_bytes, 10000 * sizeof(Tracked));
    EXPECT_EQ(memory.vector_slack_bytes, 0);
  }
  EXPECT_EQ(Tracked::live, 0);
}

TEST(PermutableSequence,
     EveryFailedIndirectMergeAllocationPreservesBothOwners) {
  using Tree = Indirect::test_tree_type;
  ResetFailures<Indirect> reset;
  ASSERT_EQ(Tracked::live, 0);
  bool succeeded = false;
  for (std::ptrdiff_t failure = 0; failure < 128 && !succeeded; ++failure) {
    Tree::test_fail_after(-1);
    auto sequence = tracked_sequence(0, 111);
    auto donor = tracked_sequence(111, 194);
    sequence.rotate_left(2, 109, 17);
    donor.rotate_left(0, donor.size(), 23);
    const auto original = addresses(sequence);
    const auto addition = addresses(donor);
    const auto live_allocations = Tree::test_counters.live_allocations;
    const auto before = sequence.memory_usage();
    const auto donor_before = donor.memory_usage();
    Tracked::moves = 0;
    Tree::test_fail_after(failure);
    Indirect::test_fail_payload_after(0);
    try {
      sequence.merge(donor);
      succeeded = true;
      auto expected = original;
      expected.insert(expected.end(), addition.begin(), addition.end());
      check_addresses(sequence, expected);
      EXPECT_TRUE(donor.empty());
      EXPECT_EQ(donor.memory_usage_bytes(), sizeof(donor));
    } catch (const std::bad_alloc&) {
      check_addresses(sequence, original);
      check_addresses(donor, addition);
      EXPECT_EQ(sequence.memory_usage_bytes(), before.total_bytes);
      EXPECT_EQ(donor.memory_usage_bytes(), donor_before.total_bytes);
      EXPECT_EQ(Tree::test_counters.live_allocations, live_allocations);
    }
    EXPECT_EQ(Tracked::moves, 0);
    EXPECT_EQ(Tracked::live, 194);
    Tree::test_fail_after(-1);
    Indirect::test_fail_payload_after(-1);
  }
  EXPECT_TRUE(succeeded);
  EXPECT_EQ(Tracked::live, 0);
  EXPECT_EQ(Tree::test_counters.live_allocations, 0);
}

TEST(PermutableSequence,
     EveryFailedIndirectRotationAllocationPreservesAddresses) {
  using Tree = Indirect::test_tree_type;
  ResetFailures<Indirect> reset;
  bool succeeded = false;
  for (std::ptrdiff_t failure = 0; failure < 256 && !succeeded; ++failure) {
    Tree::test_fail_after(-1);
    auto sequence = tracked_sequence(0, 113);
    const auto original = addresses(sequence);
    const auto live_allocations = Tree::test_counters.live_allocations;
    Tracked::moves = 0;
    Tree::test_fail_after(failure);
    try {
      sequence.rotate_left(1, 111, 41);
      succeeded = true;
      auto expected = original;
      rotate(expected, 1, 111, 41);
      check_addresses(sequence, expected);
    } catch (const std::bad_alloc&) {
      check_addresses(sequence, original);
      EXPECT_EQ(Tree::test_counters.live_allocations, live_allocations);
    }
    EXPECT_EQ(Tracked::moves, 0);
    EXPECT_EQ(Tracked::live, 113);
  }
  EXPECT_TRUE(succeeded);
  EXPECT_EQ(Tracked::live, 0);
  EXPECT_EQ(Tree::test_counters.live_allocations, 0);
}

TEST(PermutableSequence, EveryFailedPackedMergeAllocationPreservesValues) {
  using Sequence =
      PermutableSequence<std::uint64_t, ElementStorage::packed, 512, 4>;
  using Tree = Sequence::test_tree_type;
  ResetFailures<Sequence> reset;
  bool succeeded = false;
  for (std::ptrdiff_t failure = 0; failure < 128 && !succeeded; ++failure) {
    Tree::test_fail_after(-1);
    auto input = values<std::uint64_t>(111);
    auto addition = values<std::uint64_t>(83, 73);
    auto sequence = Sequence::from_range(input);
    auto donor = Sequence::from_range(addition);
    Tree::test_fail_after(failure);
    const auto live = Tree::test_counters.live_allocations;
    try {
      sequence.merge(donor);
      succeeded = true;
      input.insert(input.end(), addition.begin(), addition.end());
      check(sequence, input);
      check(donor, {});
    } catch (const std::bad_alloc&) {
      check(sequence, input);
      check(donor, addition);
      EXPECT_EQ(Tree::test_counters.live_allocations, live);
    }
  }
  EXPECT_TRUE(succeeded);
  EXPECT_EQ(Tree::test_counters.live_allocations, 0);
}

TEST(PermutableSequence, FittingMergeSucceedsWithEveryNextAllocationDisabled) {
  using Tree = Indirect::test_tree_type;
  ResetFailures<Indirect> reset;
  auto sequence = tracked_sequence(0, 1);
  auto donor = tracked_sequence(1, 2);
  const auto* first = &sequence[0];
  const auto* second = &donor[0];
  Tracked::moves = 0;
  Tree::test_fail_after(0);
  Indirect::test_fail_payload_after(0);
  EXPECT_NO_THROW(sequence.merge(donor));
  EXPECT_EQ(&sequence[0], first);
  EXPECT_EQ(&sequence[1], second);
  EXPECT_EQ(Tracked::moves, 0);
  EXPECT_TRUE(donor.empty());
}

TEST(PermutableSequence, ConstructorChunkAndVectorFailuresReclaimAllPayload) {
  ResetFailures<Indirect> reset;
  ASSERT_EQ(Tracked::live, 0);
  // 19 values at four values per chunk: exactly five new/reserve pairs.
  for (std::ptrdiff_t failure = 0; failure < 10; ++failure) {
    Indirect::test_fail_payload_after(failure);
    EXPECT_THROW(tracked_sequence(0, 19), std::bad_alloc);
    EXPECT_EQ(Tracked::live, 0);
    EXPECT_EQ(Indirect::test_tree_type::test_counters.live_allocations, 0);
  }
  Indirect::test_fail_payload_after(10);
  EXPECT_NO_THROW(tracked_sequence(0, 19));
  EXPECT_EQ(Tracked::live, 0);
}

TEST(PermutableSequence, EveryConstructorOrderAllocationFailureReclaimsChunks) {
  using Tree = Indirect::test_tree_type;
  ResetFailures<Indirect> reset;
  bool succeeded = false;
  for (std::ptrdiff_t failure = 0; failure < 256 && !succeeded; ++failure) {
    Tree::test_fail_after(failure);
    try {
      auto sequence = tracked_sequence(0, 113);
      EXPECT_EQ(sequence.size(), 113);
      succeeded = true;
    } catch (const std::bad_alloc&) {
    }
    EXPECT_EQ(Tracked::live, 0);
    EXPECT_EQ(Tree::test_counters.live_allocations, 0);
  }
  EXPECT_TRUE(succeeded);
}

struct CopyThrows {
  inline static int live = 0;
  inline static int copies_before_throw = -1;
  int value;
  explicit CopyThrows(int n) : value(n) { ++live; }
  CopyThrows(const CopyThrows& other) : value(other.value) {
    if (copies_before_throw == 0) {
      throw std::runtime_error("input copy");
    }
    if (copies_before_throw > 0) {
      --copies_before_throw;
    }
    ++live;
  }
  CopyThrows(CopyThrows&& other) noexcept : value(other.value) { ++live; }
  ~CopyThrows() noexcept { --live; }
};

TEST(PermutableSequence,
     ThrowingInputAndElementCopiesCleanUpPartialConstruction) {
  for (int failure = 0; failure < 23; ++failure) {
    auto input = std::views::iota(0, 23) | std::views::transform([&](int n) {
                   if (n == failure) {
                     throw std::runtime_error("input iterator");
                   }
                   return Tracked(n);
                 });
    EXPECT_THROW(Indirect::from_range(input), std::runtime_error);
    EXPECT_EQ(Tracked::live, 0);
    EXPECT_EQ(Indirect::test_tree_type::test_counters.live_allocations, 0);
  }
  using Sequence = PermutableSequence<CopyThrows, ElementStorage::automatic,
                                      512, 4, LengthLayout::cumulative, 16>;
  {
    std::vector<CopyThrows> input;
    input.reserve(23);
    for (int i = 0; i < 23; ++i) {
      input.emplace_back(i);
    }
    for (int failure = 0; failure < 23; ++failure) {
      CopyThrows::copies_before_throw = failure;
      EXPECT_THROW(Sequence::from_range(std::as_const(input)),
                   std::runtime_error);
      EXPECT_EQ(CopyThrows::live, 23);
      EXPECT_EQ(Sequence::test_tree_type::test_counters.live_allocations, 0);
    }
    CopyThrows::copies_before_throw = -1;
    auto sequence = Sequence::from_range(std::as_const(input));
    EXPECT_EQ(sequence.size(), 23);
  }
  EXPECT_EQ(CopyThrows::live, 0);
}

}  // namespace
