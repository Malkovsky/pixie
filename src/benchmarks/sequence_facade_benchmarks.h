#pragma once

#include <benchmark/benchmark.h>
#include <pixie/detail/sequence/packed_value_block.h>
#include <pixie/detail/sequence/sequence_tree.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>

// Shared benchmark-local workloads. Do NOT enable PIXIE_SEQUENCE_TREE_TESTING.
// Allocation events, peak requested bytes and payload move counts belong in a
// separate untimed diagnostic, not instrumented versions of these timing
// instantiations. memory_usage_bytes() reports final live requested bytes, not
// allocator/RSS or construction high-water usage. Category counters are
// descriptive, including overlapping subsets such as vector slack and order
// bytes; do not sum them all.
//
// Use pinned, serialized /bench runs, at least five repetitions, under
// /tmp/kilo/pixie-experiment-timing.lock. Select rows rather than running the
// entire registry. CPU ns/iteration is per BATCH for reads/rotations; divide by
// operations_per_iteration for ns/read or ns/rotation. items_per_second already
// uses operation units, except construction, which uses constructed elements.
// Payload access reads only the first uint64_t, NOT the whole payload object.
//
// N always counts elements, not bits or block budgets. Main footprint rows use
// 1/32/128 MiB of useful data (bool: one useful bit), plus a matched N=65536
// access row for width comparisons. Index16 never exceeds 65536 elements.
// PIXIE_FACADE_LARGE_MIB replaces 128; PIXIE_FACADE_NAVIGATION_MIB replaces
// 512. PIXIE_SEQUENCE_RAM_MIB caps the preflight budget (default cap: 4096
// MiB); half remains headroom. Linux MemAvailable is only a hint, not a cgroup
// limit. PIXIE_SEQUENCE_L3_MIB is required by the navigation-beyond-L3 row,
// which checks actual UNTAGGED internal-node bytes rather than useful payload
// bytes.

namespace {

using pixie::LengthLayout;
using pixie::detail::sequence::PackedValueBlock;
using pixie::detail::sequence::SequenceTree;

constexpr std::size_t kMiB = 1024 * 1024;
constexpr std::size_t kReads = 64;
constexpr std::size_t kRotations = 128;

std::uint64_t Mix(std::uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

std::size_t EnvironmentMiB(const char* name, std::size_t fallback = 0) {
  const char* text = std::getenv(name);
  if (text == nullptr || *text < '0' || *text > '9') {
    return fallback * kMiB;
  }
  errno = 0;
  char* end = nullptr;
  const auto value = std::strtoull(text, &end, 10);
  // Bound later conversion to element counts and signed benchmark arguments.
  constexpr auto limit =
      std::min<std::uint64_t>(std::numeric_limits<std::size_t>::max(),
                              std::numeric_limits<std::int64_t>::max()) /
      (8 * kMiB);
  if (errno != 0 || *end != '\0' || value == 0 || value > limit) {
    return fallback * kMiB;
  }
  return static_cast<std::size_t>(value) * kMiB;
}

template <std::size_t Bytes>
struct Payload {
  static_assert(Bytes >= 8 && Bytes % 8 == 0);
  std::array<std::uint64_t, Bytes / 8> words;
  explicit Payload(std::size_t i) noexcept {
    for (std::size_t j = 0; j < words.size(); ++j) {
      words[j] = Mix(i + 42 + j * 0x9e3779b97f4a7c15ULL);
    }
  }
};

struct MoveOnlyPayload : Payload<64> {
  using Payload<64>::Payload;
  MoveOnlyPayload(const MoveOnlyPayload&) = delete;
  MoveOnlyPayload& operator=(const MoveOnlyPayload&) = delete;
  MoveOnlyPayload(MoveOnlyPayload&&) noexcept = default;
  MoveOnlyPayload& operator=(MoveOnlyPayload&&) noexcept = default;
};
static_assert(sizeof(Payload<8>) == 8 && sizeof(Payload<64>) == 64 &&
              sizeof(Payload<256>) == 256 && sizeof(MoveOnlyPayload) == 64);

template <typename T>
T Value(std::size_t i) {
  if constexpr (std::is_same_v<T, bool>) {
    return (Mix(i + 42) & 1) != 0;
  } else {
    return T(i);
  }
}

template <typename T>
std::uint64_t Token(const T& value) {
  if constexpr (std::is_integral_v<T>) {
    return static_cast<std::uint64_t>(value);
  } else {
    return value.words[0];
  }
}

template <typename T, std::size_t Bits, std::size_t F, LengthLayout Layout>
struct Configuration {
  using value_type = T;
  static constexpr auto storage_bits = Bits;
  static constexpr auto fanout = F;
  static constexpr auto layout = Layout;
  static constexpr auto useful_bits =
      std::is_same_v<T, bool> ? 1 : sizeof(T) * 8;
};

template <typename Block>
Block MakeBlock(std::size_t first, std::size_t n) {
  using T = typename Block::value_type;
  std::array<T, Block::capacity> scratch{};
  for (std::size_t i = 0; i < n; ++i) {
    scratch[i] = Value<T>(first + i);
  }
  return Block(std::span<const T>(scratch.data(), n));
}

template <typename T,
          std::size_t Bits = 2048,
          std::size_t F = 8,
          LengthLayout Layout = LengthLayout::cumulative>
struct UntaggedVariant : Configuration<T, Bits, F, Layout> {
  using block_type = PackedValueBlock<T, Bits>;
  using sequence_type = SequenceTree<block_type, F, Layout>;
  static constexpr bool rebased = false;
  static constexpr bool indirect = false;
  static constexpr std::size_t chunk_bytes = 0;
  static auto Make(std::size_t n) {
    const auto blocks =
        n / block_type::capacity + (n % block_type::capacity != 0);
    auto input = std::views::iota(std::size_t{0}, blocks) |
                 std::views::transform([=](std::size_t i) {
                   const auto first = i * block_type::capacity;
                   return MakeBlock<block_type>(
                       first, std::min(block_type::capacity, n - first));
                 });
    return sequence_type::from_blocks(input);
  }
};

template <typename V>
bool Preflight(benchmark::State& state,
               std::size_t n,
               bool singletons = false,
               std::size_t live_fixtures = 1) {
  using T = typename V::value_type;
  if constexpr (V::rebased) {
    if (n != 0 && n - 1 > std::numeric_limits<T>::max()) {
      state.SkipWithError("Permutation index domain exceeded");
      return false;
    }
  }
  std::ifstream input("/proc/meminfo");
  std::string key, rest;
  std::size_t kib = 0, available = 0;
  while (input >> key >> kib) {
    if (key == "MemAvailable:") {
      available = kib * 1024;
      break;
    }
    std::getline(input, rest);
  }
  const auto explicit_cap = EnvironmentMiB("PIXIE_SEQUENCE_RAM_MIB");
  const auto cap = explicit_cap == 0 ? 4096 * kMiB : explicit_cap;
  available = available == 0 ? explicit_cap : std::min(available, cap);
  // Planning only: current blocks reserve 128 bits for local bookkeeping.
  // Budget an extra cache line per leaf/node for tags/padding, the minimum
  // branching factor, 2x occupancy/construction slack, and small-tree preflight
  // spares. Floating arithmetic avoids overflow for user-supplied sizes.
  // Singleton chunks may each reserve the entire chunk budget.
  constexpr auto field_bits = V::indirect ? sizeof(void*) * 8 : V::useful_bits;
  constexpr auto capacity =
      std::max<std::size_t>(1, (V::storage_bits - 128) / field_bits);
  const long double leaves = singletons ? n : 1 + n / capacity;
  const long double nodes = 1 + leaves / (V::fanout / 2 - 1);
  long double estimate =
      128 * 1024 + 2 * (leaves * (V::storage_bits / 8 + 64) +
                        nodes * (2 * V::fanout * sizeof(void*) + 64));
  if constexpr (V::indirect) {
    const auto slots = std::max<std::size_t>(1, V::chunk_bytes / sizeof(T));
    const long double chunks = singletons ? n : 1 + n / slots;
    estimate += 2 * chunks * (std::max(V::chunk_bytes, sizeof(T)) + 128);
  }
  estimate *= live_fixtures;
  state.counters["available_budget_bytes"] = available;
  state.counters["planning_bytes"] = static_cast<double>(estimate);
  if (available == 0 || estimate > available / 2) {
    state.SkipWithError(
        "RAM preflight refused row; review PIXIE_SEQUENCE_RAM_MIB");
    return false;
  }
  return true;
}

template <typename V>
void Memory(benchmark::State& state,
            const typename V::sequence_type& sequence) {
  const auto memory = sequence.memory_usage();
  const auto bytes = memory.total_bytes;
  state.counters["blocks"] = memory.blocks;
  state.counters["nodes"] = memory.nodes;
  if constexpr (requires { memory.tag_padding_bytes; }) {
    using Tree =
        SequenceTree<PackedValueBlock<typename V::value_type, V::storage_bits>,
                     V::fanout, V::layout, true>;
    state.counters["tree_block_bytes"] =
        memory.blocks * Tree::block_storage_bytes;
    state.counters["tree_node_bytes"] = memory.nodes * Tree::node_storage_bytes;
    state.counters["payload_capacity_bytes"] = memory.payload_capacity_bytes;
    state.counters["block_metadata_bytes"] = memory.block_metadata_bytes;
    state.counters["ordering_tree_bytes"] = memory.ordering_tree_bytes;
    state.counters["tag_padding_bytes"] = memory.tag_padding_bytes;
    state.counters["facade_bytes"] = memory.facade_bytes;
  } else if constexpr (requires { memory.chunk_header_bytes; }) {
    state.counters["tree_block_bytes"] = memory.tree_block_bytes;
    state.counters["tree_node_bytes"] = memory.tree_node_bytes;
    state.counters["order_bytes"] = memory.order_bytes;
    state.counters["packed_capacity_bits"] = memory.packed_capacity_bits;
    state.counters["chunks"] = memory.chunks;
    state.counters["chunk_header_bytes"] = memory.chunk_header_bytes;
    state.counters["vector_capacity_bytes"] = memory.vector_capacity_bytes;
    state.counters["vector_live_bytes"] = memory.vector_live_bytes;
    state.counters["vector_slack_bytes"] = memory.vector_slack_bytes;
    state.counters["facade_bytes"] = memory.facade_bytes;
  } else {
    state.counters["tree_block_bytes"] = memory.block_bytes;
    state.counters["tree_node_bytes"] = memory.node_bytes;
  }
  state.counters["elements"] = sequence.size();
  state.counters["useful_bytes"] =
      static_cast<double>(sequence.size()) * V::useful_bits / 8;
  state.counters["owned_bytes"] = bytes;
  state.counters["bytes_per_element"] =
      sequence.empty() ? 0 : static_cast<double>(bytes) / sequence.size();
  // sizeof is a separate descriptive counter, NOT an extra summand of bytes.
  state.counters["sizeof_facade_or_tree"] = sizeof(sequence);
  if constexpr (requires { V::immer; }) {
    state.counters["leaf_capacity_elements"] = V::sequence_type::leaf_capacity;
  } else {
    state.counters["block_budget_bytes"] = V::storage_bits / 8;
    state.counters["cumulative_lengths"] =
        V::layout == LengthLayout::cumulative;
  }
  state.counters["chunk_budget_bytes"] = V::chunk_bytes;
  state.counters["sizeof_value"] = sizeof(typename V::value_type);
  state.counters["useful_bits_per_element"] = V::useful_bits;
  state.counters["fanout"] = V::fanout;
  state.counters["rebased_merge"] = V::rebased;
  state.counters["indirect_storage"] = V::indirect;
}

template <typename V,
          bool Dependent,
          bool BeyondL3 = false,
          bool Edited = false,
          bool Aligned = false>
void Access(benchmark::State& state) {
  static_assert(!(Edited && Aligned));
  const std::size_t n = state.range(0);
  if (!Preflight<V>(state, n)) {
    return;
  }
  const auto l3 = EnvironmentMiB("PIXIE_SEQUENCE_L3_MIB");
  if constexpr (BeyondL3) {
    if (l3 == 0) {
      state.SkipWithError(
          "Declare measured shared L3 via PIXIE_SEQUENCE_L3_MIB");
      return;
    }
  }
  const auto sequence = [&] {
    auto result = V::Make(n);
    if constexpr (Edited) {
      for (std::size_t i = 0; i < kRotations; ++i) {
        const auto left = Mix(i + 731) % (n / 2);
        result.rotate_left(left, left + n / 2, 1 + Mix(i + 917) % (n / 2 - 1));
      }
    }
    if constexpr (Aligned) {
      // Force ordering metadata into use, rather than mostly repacking nodes
      // through arbitrary cuts. Every backend receives the same disjoint edits.
      for (std::size_t left = 0; left < n; left += 4096) {
        const auto length = std::min<std::size_t>(4096, n - left);
        result.rotate_left(left, left + length,
                           std::min<std::size_t>(256, length - 1));
      }
    }
    return result;
  }();
  if constexpr (Edited) {
    for (std::size_t i = 0; i < kReads; ++i) {
      const auto position = Mix(i) % n;
      auto original = position;
      for (std::size_t j = kRotations; j != 0; --j) {
        const auto left = Mix(j - 1 + 731) % (n / 2);
        const auto distance = 1 + Mix(j - 1 + 917) % (n / 2 - 1);
        if (original >= left && original < left + n / 2) {
          original = left + (original - left + distance) % (n / 2);
        }
      }
      if (Token(sequence[position]) !=
          Token(Value<typename V::value_type>(original))) {
        state.SkipWithError(
            "Edited read fixture failed inverse-position oracle");
        return;
      }
    }
    state.counters["setup_rotations"] = kRotations;
  }
  if constexpr (Aligned) {
    for (std::size_t i = 0; i < kReads; ++i) {
      const auto position = Mix(i) % n;
      const auto left = position / 4096 * 4096;
      const auto length = std::min<std::size_t>(4096, n - left);
      const auto original =
          left +
          (position - left + std::min<std::size_t>(256, length - 1)) % length;
      if (Token(sequence[position]) !=
          Token(Value<typename V::value_type>(original))) {
        state.SkipWithError(
            "Aligned read fixture failed inverse-position oracle");
        return;
      }
    }
    state.counters["setup_rotations"] = (n + 4095) / 4096;
    state.counters["setup_aligned_group_values"] = 4096;
  }
  Memory<V>(state, sequence);
  if constexpr (BeyondL3) {
    // Only instantiated for the existing untagged tree introspection API.
    const auto internal = sequence.internal_memory_bytes();
    state.counters["internal_bytes"] = internal;
    state.counters["declared_l3_bytes"] = l3;
    state.counters["internal_over_l3"] = static_cast<double>(internal) / l3;
    if (internal <= l3) {
      state.SkipWithError(
          "Navigation does not exceed L3; increase navigation MiB");
      return;
    }
  }
  state.SetLabel(
      Dependent
          ? "64 dependent reads; scalar/first word only; hash included"
          : "64 independent reads; scalar/first word only; hash included");
  // Fresh full-domain positions, not a small repeated query pool. A payload
  // reference is never copied merely to read its first word.
  std::uint64_t counter = 123, previous = 0, checksum = 0;
  for (auto _ : state) {
    for (std::size_t i = 0; i < kReads; ++i) {
      auto key = counter++;
      if constexpr (Dependent) {
        key ^= previous * 0x9e3779b97f4a7c15ULL;
      }
      auto value = Token(sequence[Mix(key) % n]);
      checksum += value;
      if constexpr (Dependent) {
        previous = value;
      }
    }
    // Consume arithmetic on the payload, not only a potentially folded load.
    // A batch sink also preserves the dependent-read chain through value.
    benchmark::DoNotOptimize(checksum);
  }
  state.counters["operations_per_iteration"] = kReads;
  state.counters["observed_bits_per_read"] =
      std::min<std::size_t>(64, V::useful_bits);
  state.SetItemsProcessed(state.iterations() * kReads);
}

enum class Rotation {
  Whole,
  Global,
  Local,
  ChildAligned16,
  ChildAligned32,
  ChildAligned64,
  GroupedSpill256
};
struct Query {
  std::size_t left, right, distance;
};

template <Rotation Mode>
auto RotationQueries(std::size_t n) {
  std::array<Query, kRotations> queries;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    const auto key = Mix(i + 128);
    if constexpr (Mode == Rotation::GroupedSpill256) {
      // One 256-child leaf-parent node, using the same element range for
      // every backend. The distance can require pointers to cross groups.
      const auto left = (key % (n / 65536)) * 65536;
      queries[i] = {left, left + 65536, (1 + Mix(key) % 255) * 256};
      continue;
    }
    if constexpr (Mode == Rotation::ChildAligned16 ||
                  Mode == Rotation::ChildAligned32 ||
                  Mode == Rotation::ChildAligned64) {
      // Identical element ranges for every backend. All three boundaries
      // lie on 256-value leaf boundaries inside the selected group size.
      constexpr std::size_t children = Mode == Rotation::ChildAligned16   ? 16
                                       : Mode == Rotation::ChildAligned32 ? 32
                                                                          : 64;
      const auto group = key % (n / (children * 256));
      const auto first = Mix(key) % (children / 2);
      const auto count = 2 + Mix(key + 1) % (children - 1 - first);
      const auto left = group * (children * 256) + first * 256;
      queries[i] = {left, left + count * 256,
                    (1 + Mix(key + 2) % (count - 1)) * 256};
      continue;
    }
    const auto length = Mode == Rotation::Whole   ? n
                        : Mode == Rotation::Local ? std::min<std::size_t>(65, n)
                                                  : n / 2 + key % (n / 4);
    const auto left = Mix(key) % (n - length + 1);
    queries[i] = {left, left + length, 1 + Mix(key + 1) % (length - 1)};
  }
  return queries;
}

template <typename V, Rotation Mode>
void RotateBatch(benchmark::State& state) {
  using Sequence = typename V::sequence_type;
  const std::size_t n = state.range(0);
  if (!Preflight<V>(state, n)) {
    return;
  }
  const auto queries = RotationQueries<Mode>(n);
  Sequence sequence;
  state.SetLabel(
      "128 rotations/batch; fixed stream; rebuild/destruction excluded; "
      "rebuild pre-touches data; allocation+repair included");
  for (auto _ : state) {
    state.PauseTiming();
    sequence = Sequence{};
    sequence = V::Make(n);
    state.ResumeTiming();
    for (const auto& q : queries) {
      sequence.rotate_left(q.left, q.right, q.distance);
      benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(sequence);
  }
  if (sequence.size() != n) {
    state.SkipWithError("Rotation changed size");
    return;
  }
  // Bounded inverse-position oracle; no complete reference input allocation.
  const auto matches = [&](std::size_t position) {
    auto original = position;
    for (auto q = queries.rbegin(); q != queries.rend(); ++q) {
      if (original >= q->left && original < q->right) {
        original =
            q->left + (original - q->left + q->distance) % (q->right - q->left);
      }
    }
    return Token(sequence[position]) ==
           Token(Value<typename V::value_type>(original));
  };
  for (std::size_t probe = 0; probe < 8; ++probe) {
    if (!matches(Mix(probe + 19) % n)) {
      state.SkipWithError(
          "Rotation disagrees with sampled inverse-position oracle");
      return;
    }
  }
  // Follow one affected output position from every operation through later
  // rotations. Uniform probes alone almost always miss the short local ranges.
  for (std::size_t i = 0; i < queries.size(); ++i) {
    auto position = queries[i].left;
    for (std::size_t j = i + 1; j < queries.size(); ++j) {
      const auto& q = queries[j];
      if (position >= q.left && position < q.right) {
        const auto offset = position - q.left;
        const auto length = q.right - q.left;
        position =
            q.left + (offset >= q.distance ? offset - q.distance
                                           : length - (q.distance - offset));
      }
    }
    if (!matches(position)) {
      state.SkipWithError("Rotation disagrees at an affected position");
      return;
    }
  }
  Memory<V>(state, sequence);
  state.counters["operations_per_iteration"] = kRotations;
  state.counters["query_pool_bytes"] = sizeof(queries);
  state.SetItemsProcessed(state.iterations() * kRotations);
}

template <typename V, bool Equal>
void Merge(benchmark::State& state) {
  using Sequence = typename V::sequence_type;
  const std::size_t n = state.range(0);
  // Keep the same batch count for equal logical representations, with at least
  // two independent merges per timing interval and bounded fixture storage.
  constexpr std::size_t max_batch = 8;
  const auto useful_bytes = static_cast<long double>(n) * V::useful_bits / 8;
  const std::size_t batch = static_cast<std::size_t>(std::clamp(
      64.0L * kMiB / useful_bytes, 2.0L, static_cast<long double>(max_batch)));
  if (!Preflight<V>(state, n, false, batch)) {
    return;
  }
  const auto donor_n = Equal ? n / 2 : n / 64;
  const auto receiver_n = n - donor_n;
  std::array<Sequence, max_batch> receivers, donors;
  state.SetLabel(
      V::rebased
          ? "batched rebased merges; independent rotated identities built "
            "untimed"
          : "batched unchanged-value merges; independent rotated ranges built "
            "untimed");
  for (auto _ : state) {
    state.PauseTiming();
    for (std::size_t i = 0; i < batch; ++i) {
      receivers[i] = Sequence{};
      donors[i] = Sequence{};
      receivers[i] = V::Make(receiver_n);
      donors[i] = V::Make(donor_n);
      receivers[i].rotate_left(0, receiver_n, receiver_n / 3);
      donors[i].rotate_left(0, donor_n, donor_n / 3);
    }
    state.ResumeTiming();
    for (std::size_t i = 0; i < batch; ++i) {
      receivers[i].merge(donors[i]);
      benchmark::DoNotOptimize(receivers[i]);
      benchmark::ClobberMemory();
    }
  }
  auto donor_first = Token(Value<typename V::value_type>(donor_n / 3));
  if constexpr (V::rebased) {
    donor_first += receiver_n;
  }
  std::size_t fixture_bytes = 0;
  for (std::size_t i = 0; i < batch; ++i) {
    const auto& receiver = receivers[i];
    if (receiver.size() != n || !donors[i].empty() ||
        Token(receiver[0]) !=
            Token(Value<typename V::value_type>(receiver_n / 3)) ||
        Token(receiver[receiver_n]) != donor_first) {
      state.SkipWithError("Merge size/donor/value invariant failed");
      return;
    }
    fixture_bytes +=
        receiver.memory_usage_bytes() + donors[i].memory_usage_bytes();
  }
  Memory<V>(state, receivers[0]);
  state.counters["active_fixture_bytes"] = fixture_bytes;
  state.counters["receiver_elements"] = receiver_n;
  state.counters["donor_elements"] = donor_n;
  state.counters["operations_per_iteration"] = batch;
  state.SetItemsProcessed(state.iterations() * batch);
}

template <typename V, bool Singletons = false>
void Build(benchmark::State& state) {
  using Sequence = typename V::sequence_type;
  const std::size_t n = state.range(0);
  if (!Preflight<V>(state, n, Singletons)) {
    return;
  }
  state.SetLabel(
      Singletons
          ? "singleton construction+merge; unchanged repeated values unless "
            "rebased; final destruction excluded"
          : "streaming value generation+build; final destruction excluded");
  bool reported = false;
  for (auto _ : state) {
    auto sequence = [&] {
      if constexpr (Singletons) {
        Sequence result;
        for (std::size_t i = 0; i < n; ++i) {
          auto donor = V::Make(1);
          result.merge(donor);
        }
        return result;
      } else {
        return V::Make(n);
      }
    }();
    benchmark::DoNotOptimize(sequence);
    benchmark::ClobberMemory();
    state.PauseTiming();
    if (!reported) {
      Memory<V>(state, sequence);
      reported = true;
    }
    const auto expected_last = Singletons && !V::rebased ? 0 : n - 1;
    const bool valid = sequence.size() == n &&
                       Token(sequence[n - 1]) ==
                           Token(Value<typename V::value_type>(expected_last));
    sequence = Sequence{};
    state.ResumeTiming();
    if (!valid) {
      state.SkipWithError("Construction size/value invariant failed");
      break;
    }
  }
  state.counters["operations_per_iteration"] = n;
  state.SetItemsProcessed(state.iterations() * n);
  state.SetBytesProcessed(state.iterations() * (n * V::useful_bits / 8));
}

// Position generation, fixture construction and verification are outside
// timing. The same growing-size position stream is used for every backend.
template <typename V>
void InsertBatch(benchmark::State& state) {
  using Sequence = typename V::sequence_type;
  constexpr std::size_t count = 128;
  const auto n = static_cast<std::size_t>(state.range(0));
  if (!Preflight<V>(state, n + count)) {
    return;
  }
  std::array<std::size_t, count> positions;
  for (std::size_t i = 0; i < count; ++i) {
    positions[i] = Mix(i + 731) % (n + i + 1);
  }
  Sequence sequence;
  state.SetLabel(
      "128 positional inserts; search/build/destruction excluded; "
      "allocation and leaf repair included");
  for (auto _ : state) {
    state.PauseTiming();
    sequence = Sequence{};
    sequence = V::Make(n);
    state.ResumeTiming();
    for (std::size_t i = 0; i < count; ++i) {
      sequence.insert_at(positions[i], Value<typename V::value_type>(n + i));
      benchmark::ClobberMemory();
    }
    benchmark::DoNotOptimize(sequence);
  }
  if (sequence.size() != n + count) {
    state.SkipWithError("Insertion size invariant failed");
    return;
  }
  for (std::size_t i = 0; i < count; ++i) {
    auto final_position = positions[i];
    for (std::size_t j = i + 1; j < count; ++j) {
      final_position += positions[j] <= final_position;
    }
    if (Token(sequence[final_position]) !=
        Token(Value<typename V::value_type>(n + i))) {
      state.SkipWithError("Inserted value disagrees at final position");
      return;
    }
  }
  Memory<V>(state, sequence);
  state.counters["operations_per_iteration"] = count;
  state.SetItemsProcessed(state.iterations() * count);
}

template <typename Block, bool Rotate>
void HotLeaf(benchmark::State& state) {
  const auto n = Block::capacity - 3;
  auto block = MakeBlock<Block>(0, n);
  const auto queries = RotationQueries<Rotation::Global>(n);
  state.SetLabel(Rotate
                     ? "128 hot-leaf rotations; rebuild excluded; no tree/tags"
                     : "64 hot-leaf reads; hash included; no tree/tags");
  std::uint64_t counter = 123;
  for (auto _ : state) {
    if constexpr (Rotate) {
      state.PauseTiming();
      block = MakeBlock<Block>(0, n);
      state.ResumeTiming();
      for (const auto& q : queries) {
        block.rotate_left(q.left, q.right, q.distance);
        benchmark::ClobberMemory();
      }
      benchmark::DoNotOptimize(block);
    } else {
      for (std::size_t i = 0; i < kReads; ++i) {
        auto value = block[Mix(counter++) % n];
        benchmark::DoNotOptimize(value);
      }
    }
  }
  state.counters["elements"] = n;
  state.counters["block_capacity_elements"] = Block::capacity;
  state.counters["sizeof_block"] = sizeof(Block);
  state.counters["operations_per_iteration"] = Rotate ? kRotations : kReads;
  state.SetItemsProcessed(state.iterations() * (Rotate ? kRotations : kReads));
}

// Only the default variants get the whole workload set. Secondary value/storage
// axes get reads, global rotation and construction, with merge comparisons for
// every payload size and index width. Layout knobs get reads and global
// rotation at 32 MiB only. This is a screen, not a Cartesian product.
template <typename V>
void Sizes(benchmark::internal::Benchmark* row, bool matched = false) {
  row->ArgName("N");
  using T = typename V::value_type;
  if constexpr (std::is_same_v<T, std::uint16_t>) {
    row->Arg(4096)->Arg(65536);
  } else {
    if (matched) {
      row->Arg(65536);
    }
    row->Arg(kMiB * 8 / V::useful_bits);
    row->Arg(32 * kMiB * 8 / V::useful_bits);
    const auto large = EnvironmentMiB("PIXIE_FACADE_LARGE_MIB", 128);
    if (large != kMiB && large != 32 * kMiB) {
      row->Arg(large * 8 / V::useful_bits);
    }
  }
}

template <typename V, bool Full = false, bool Control = false>
void Register(const char* variant) {
  const auto add = [=](const char* operation, auto function, bool fixed = false,
                       bool matched = false, std::int64_t iterations = 8) {
    const auto name = std::string(operation) + "/" + variant;
    auto* row = benchmark::RegisterBenchmark(name.c_str(), function);
    if constexpr (Control) {
      row->ArgName("N")->Arg(32 * kMiB * 8 / V::useful_bits);
    } else {
      Sizes<V>(row, matched);
    }
    // Avoid adaptive millions of untimed full reconstructions for a fast merge.
    // Each repetition starts with the same operations/topology, including all
    // rotation batches. Increase repetitions for marginal comparisons.
    if (fixed) {
      row->Iterations(iterations);
    }
  };
  add("FacadeAccessIndependent", Access<V, false>, false, true);
  constexpr bool is_bool = std::is_same_v<typename V::value_type, bool>;
  add("FacadeRotateGlobalBatch", RotateBatch<V, Rotation::Global>, true,
      is_bool);
  if constexpr (!Control) {
    add("FacadeBuild", Build<V>, false, is_bool);
  }
  if constexpr (Full) {
    add("FacadeAccessDependent", Access<V, true>);
    add("FacadeRotateWholeBatch", RotateBatch<V, Rotation::Whole>, true);
    add("FacadeRotateLocal65Batch", RotateBatch<V, Rotation::Local>, true);
    const auto name = std::string("FacadeSingletonBuildMerge/") + variant;
    benchmark::RegisterBenchmark(name.c_str(), Build<V, true>)
        ->ArgName("N")
        ->Arg(256)
        ->Arg(4096)
        ->Iterations(8);
  }
  if constexpr (!Control && (Full || V::indirect || V::rebased)) {
    add("FacadeMergeEqual", Merge<V, true>, true, true, 32);
    add("FacadeMergeUnequal", Merge<V, false>, true, true, 32);
  }
}

template <typename T, std::size_t Bits = 2048>
void RegisterLeaf(const char* variant) {
  using Block = PackedValueBlock<T, Bits>;
  const auto read = std::string("FacadeLeafRead/") + variant;
  const auto rotate = std::string("FacadeLeafRotateBatch/") + variant;
  benchmark::RegisterBenchmark(read.c_str(), HotLeaf<Block, false>);
  benchmark::RegisterBenchmark(rotate.c_str(), HotLeaf<Block, true>);
}

}  // namespace
