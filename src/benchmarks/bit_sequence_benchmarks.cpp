#include <benchmark/benchmark.h>
#include <pixie/bits.h>
#include <pixie/detail/sequence/bit_block.h>
#include <pixie/detail/sequence/packed_bit_block.h>
#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/experimental/permuted_bit_block.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <random>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include "packed_copy_kernels.h"

// MiB rows use fixed useful payload sizes; budget_MiB access rows instead bound
// requested live tree storage, excluding allocator bookkeeping and RSS.
// Select exact rows through /bench; do not run the entire size/workload matrix
// for an initial screen. Serialize timing with
// flock /tmp/kilo/pixie-experiment-timing.lock <benchmark wrapper command>.
// No results from the retired flat-directory implementation apply here.

namespace {

using CopyFunction =
    void (*)(const uint64_t*, size_t, uint64_t*, size_t, size_t);

// Fixed addresses and query stream across variants. Layout 0: word-aligned;
// 1: source shifted by 17; 2: mixed bit and cache-line offsets, variable
// lengths. Warm reusable buffers, not a streaming-memory benchmark. Timing
// includes the query lookup, indirect call, boundary work and store barrier,
// not allocation.
void PackedCopy(benchmark::State& state, CopyFunction copy) {
  const size_t bytes = state.range(0);
  const auto layout = state.range(1);
  alignas(64) static std::array<uint64_t, 16416> input, output;
  std::mt19937_64 random(9182);
  for (auto& word : input) {
    word = random();
  }
  struct Query {
    size_t source, destination, count;
  };
  std::array<Query, 256> queries;
  for (auto& q : queries) {
    q = layout == 0 ? Query{0, 0, bytes * 8}
        : layout == 1
            ? Query{17, 0, bytes * 8}
            : Query{random() % 512, random() % 512, bytes * 8 - random() % 128};
  }
  // Validate controls against a bitwise oracle outside timing, including bits
  // outside the copied range. Exhaustive exact-extent checks live in tests.
  for (const auto& q : queries) {
    output.fill(0xabcdef0123456789ULL);
    auto expected = output;
    for (size_t i = 0; i < q.count; ++i) {
      const auto s = q.source + i;
      const auto d = q.destination + i;
      const auto mask = uint64_t{1} << (d % 64);
      expected[d / 64] = (expected[d / 64] & ~mask) |
                         (((input[s / 64] >> (s % 64)) & 1) << (d % 64));
    }
    copy(input.data(), q.source, output.data(), q.destination, q.count);
    if (output != expected) {
      state.SkipWithError("Packed copy disagrees with bitwise oracle");
      return;
    }
  }
  size_t i = 0;
  int64_t bits = 0;
  for (auto _ : state) {
    const auto& q = queries[i++ % queries.size()];
    copy(input.data(), q.source, output.data(), q.destination, q.count);
    bits += q.count;
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(bits / 8);
}

const bool packed_copy_registered = [] {
  const auto add = [](const char* name, CopyFunction copy) {
    auto* benchmark = benchmark::RegisterBenchmark(name, PackedCopy, copy);
    benchmark->ArgNames({"bytes", "layout"});
    for (int bytes : {64, 128, 256, 512, 4096, 131072}) {
      for (int layout : {0, 1, 2}) {
        benchmark->Args({bytes, layout});
      }
    }
  };
  using namespace packed_copy_benchmark;
  add("PackedCopy/production", Production);
#if defined(__AVX512F__)
  add("PackedCopy/shift1", Copy<false, 1>);
  add("PackedCopy/shift2", Copy<false, 2>);
  add("PackedCopy/shift4", Copy<false, 4>);
#if defined(__AVX512VBMI2__)
  add("PackedCopy/funnel1", Copy<true, 1>);
  add("PackedCopy/funnel2", Copy<true, 2>);
  add("PackedCopy/funnel4", Copy<true, 4>);
#endif
#endif
  return true;
}();

using BitBlock = pixie::detail::sequence::BitBlock<2048>;
using Packed128 = pixie::detail::sequence::PackedBitBlock<128 * 8>;
using Packed256 = pixie::detail::sequence::PackedBitBlock<256 * 8>;
using Packed512 = pixie::detail::sequence::PackedBitBlock<512 * 8>;
using MatchedBitBlock = pixie::detail::sequence::BitBlock<Packed256::capacity>;
using DirectBlock = pixie::experimental::PermutedBitBlock<false>;
using MappedBlock = pixie::experimental::PermutedBitBlock<true>;
using pixie::LengthLayout;

// Keep benchmark traits here rather than requiring introspection aliases in
// the experimental public API. Screen axes, not the full template product.
template <typename Block, std::size_t Fanout, LengthLayout Layout>
struct Variant {
  using block_type = Block;
  using tree_type =
      pixie::detail::sequence::SequenceTree<Block, Fanout, Layout>;
  static constexpr auto fanout = Fanout;
  static constexpr auto layout = Layout;
};
using P256F8Raw = Variant<Packed256, 8, LengthLayout::individual>;
using P256F8Prefix = Variant<Packed256, 8, LengthLayout::cumulative>;
using P256F4Prefix = Variant<Packed256, 4, LengthLayout::cumulative>;
using P256F4Raw = Variant<Packed256, 4, LengthLayout::individual>;
using P256F16Prefix = Variant<Packed256, 16, LengthLayout::cumulative>;
using P256F16Raw = Variant<Packed256, 16, LengthLayout::individual>;
using P128F8Prefix = Variant<Packed128, 8, LengthLayout::cumulative>;
using P512F8Prefix = Variant<Packed512, 8, LengthLayout::cumulative>;
using DirectMatchedF8Prefix =
    Variant<MatchedBitBlock, 8, LengthLayout::cumulative>;

constexpr std::size_t kMiB = 1024 * 1024;
constexpr std::int64_t kQueriesPerIteration = 64;

// Stateless SplitMix64 finalizer: input words are reproducible by global word
// index, independent of block capacity, without a second sequence-sized buffer.
std::uint64_t Mix(std::uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

template <typename Block>
Block MakeBlock(std::size_t first_bit, std::size_t n) {
  std::array<std::uint64_t, (Block::capacity + 63) / 64> words{};
  for (std::size_t i = 0; i < (n + 63) / 64; ++i) {
    words[i] = Mix(first_bit / 64 + i + 42);
  }
  return Block(std::span<const std::uint64_t>(words), n);
}

template <typename V>
auto MakeTree(std::size_t n,
              std::size_t component_bits = V::block_type::capacity) {
  using Block = typename V::block_type;
  using Tree = typename V::tree_type;
  static_assert(Block::capacity % 64 == 0);
  const auto count = n / component_bits + (n % component_bits != 0);
  // The view yields owning blocks by value; no input/leaf directory is kept.
  // All registered component sizes are whole words, even for short components.
  auto blocks =
      std::views::iota(std::size_t{0}, count) |
      std::views::transform([=](std::size_t i) {
        const auto start = i * component_bits;
        return MakeBlock<Block>(start, std::min(component_bits, n - start));
      });
  return Tree::from_blocks(blocks);
}

std::size_t EnvironmentMiB(const char* name) {
  const char* text = std::getenv(name);
  if (text == nullptr || *text < '0' || *text > '9') {
    return 0;
  }
  char* end = nullptr;
  const auto value = std::strtoull(text, &end, 10);
  if (*end != '\0' || value > std::numeric_limits<std::size_t>::max() / kMiB) {
    return 0;
  }
  return static_cast<std::size_t>(value) * kMiB;
}

template <typename V>
bool CheckMemoryBudget(benchmark::State& state,
                       std::size_t n,
                       std::size_t extra_bytes = 0) {
  // Linux MemAvailable is a preflight hint, not a cgroup limit, reservation or
  // measured peak. For a tighter container/process budget set
  // PIXIE_SEQUENCE_RAM_MIB. Leave half that available budget as headroom.
  std::ifstream input("/proc/meminfo");
  std::string key, rest;
  std::size_t kib = 0;
  std::size_t available = 0;
  while (input >> key >> kib) {
    if (key == "MemAvailable:") {
      available = kib * 1024;
      break;
    }
    std::getline(input, rest);
  }
  const auto cap = EnvironmentMiB("PIXIE_SEQUENCE_RAM_MIB");
  if (cap != 0) {
    available = available == 0 ? cap : std::min(available, cap);
  }
  using Block = typename V::block_type;
  const auto leaves = n / Block::capacity + (n % Block::capacity != 0);
  // Planning estimate only: allow a 256-byte node per leaf and 2x slack for
  // repair/construction. Actual requested live storage is reported separately.
  const auto estimate = 2 * leaves * (sizeof(Block) + 256) + extra_bytes;
  state.counters["available_budget_bytes"] = available;
  state.counters["planning_bytes"] = estimate;
  if (available == 0 || estimate > available / 2) {
    state.SkipWithError(
        "Insufficient/unknown available RAM; check "
        "PIXIE_SEQUENCE_RAM_MIB before selecting a large row");
    return false;
  }
  return true;
}

template <typename V>
void TreeMemory(benchmark::State& state, const typename V::tree_type& tree) {
  using Block = typename V::block_type;
  using Tree = typename V::tree_type;
  // Final live requested bytes only: no allocator bookkeeping, RSS, allocation
  // event counts, or temporary high-water claims. Do not sum overlapping
  // fields. node_count() is expected to count internal nodes only; each node
  // and each block has one separate allocation under the agreed tree ownership
  // model.
  const auto memory = tree.memory_usage();
  const auto leaves = memory.blocks;
  const auto nodes = memory.nodes;
  state.counters["logical_bits"] = tree.size();
  state.counters["useful_bytes"] = static_cast<double>(tree.size()) / 8;
  state.counters["payload_capacity_bytes"] =
      leaves * ((Block::capacity + 7) / 8);
  state.counters["internal_bytes"] = memory.node_bytes;
  state.counters["owned_bytes"] = memory.total_bytes;
  state.counters["bytes_per_bit"] =
      static_cast<double>(memory.total_bytes) / tree.size();
  state.counters["leaf_object_bytes"] = leaves * sizeof(Block);
  state.counters["leaf_allocation_bytes"] = memory.block_bytes;
  state.counters["tree_object_bytes"] = sizeof(tree);
  state.counters["blocks"] = leaves;
  state.counters["internal_nodes"] = nodes;
  state.counters["live_allocations"] = leaves + nodes;
  state.counters["height"] = tree.height();
  state.counters["block_capacity_bits"] = Block::capacity;
  state.counters["sizeof_block"] = sizeof(Block);
  state.counters["sizeof_leaf_allocation"] = Tree::block_storage_bytes;
  state.counters["sizeof_node_allocation"] = Tree::node_storage_bytes;
  state.counters["alignof_block"] = alignof(Block);
  state.counters["leaf_occupancy"] =
      static_cast<double>(tree.size()) / (leaves * Block::capacity);
  state.counters["fanout"] = V::fanout;
  state.counters["cumulative_lengths"] = V::layout == LengthLayout::cumulative;
  const auto l3 = EnvironmentMiB("PIXIE_SEQUENCE_L3_MIB");
  if (l3 != 0) {
    state.counters["declared_l3_bytes"] = l3;
    state.counters["internal_over_l3"] =
        static_cast<double>(memory.node_bytes) / l3;
  }
}

template <typename V>
auto MakeBudgetTree(std::size_t budget) {
  using Block = typename V::block_type;
  // Start above the target: leaf objects alone use almost the entire budget.
  // Trim an owning suffix, at block boundaries, to include actual internal-node
  // storage in the budget. No simultaneous replacement tree is constructed.
  auto tree = MakeTree<V>((budget / sizeof(Block)) * Block::capacity);
  for (auto owned = tree.memory_usage_bytes(); owned > budget;
       owned = tree.memory_usage_bytes()) {
    const auto leaves = tree.size() / Block::capacity;
    const auto retained = std::min(leaves - 1, leaves * budget / owned);
    auto suffix = tree.split_off(retained * Block::capacity);
  }
  return tree;
}

template <typename V,
          bool Dependent,
          bool RequireColdNavigation = false,
          bool MemoryBudget = false>
void TreeAccess(benchmark::State& state) {
  const std::size_t bytes = state.range(0) * kMiB;
  const auto initial_bits =
      MemoryBudget
          ? (bytes / sizeof(typename V::block_type)) * V::block_type::capacity
          : bytes * 8;
  if (!CheckMemoryBudget<V>(state, initial_bits)) {
    return;
  }
  const auto l3 = EnvironmentMiB("PIXIE_SEQUENCE_L3_MIB");
  if constexpr (RequireColdNavigation) {
    if (l3 == 0) {
      state.SkipWithError(
          "Declare measured shared L3 via PIXIE_SEQUENCE_L3_MIB");
      return;
    }
  }
  auto tree = [&] {
    if constexpr (MemoryBudget) {
      return MakeBudgetTree<V>(bytes);
    } else {
      return MakeTree<V>(initial_bits);
    }
  }();
  const auto n = tree.size();
  TreeMemory<V>(state, tree);
  if constexpr (MemoryBudget) {
    state.counters["requested_budget_bytes"] = bytes;
    state.counters["budget_utilization"] =
        static_cast<double>(tree.memory_usage_bytes()) / bytes;
    if (tree.memory_usage_bytes() < bytes * 99 / 100) {
      state.SkipWithError("Requested-memory comparison undershot by over 1%");
      return;
    }
  }
  if constexpr (RequireColdNavigation) {
    if (tree.internal_memory_bytes() <= l3) {
      state.SkipWithError(
          "Internal navigation storage does not exceed declared L3");
      return;
    }
  }
  // Construction has written/touched every live allocation, outside timing.
  // A full-width counter supplies fresh positions over all N, not a query pool.
  // Both streams include hash/modulo generation. The dependent stream feeds
  // the returned random input bit into the next hash, not just a fixed branch.
  std::uint64_t counter = 123;
  std::uint64_t previous = 0;
  state.SetLabel(Dependent
                     ? "64 dependent reads/iteration; generator included"
                     : "64 independent reads/iteration; generator included");
  for (auto _ : state) {
    for (std::int64_t i = 0; i < kQueriesPerIteration; ++i) {
      auto key = counter++;
      if constexpr (Dependent) {
        key ^= previous * 0x9e3779b97f4a7c15ULL;
      }
      auto bit = tree[Mix(key) % n];
      benchmark::DoNotOptimize(bit);
      if constexpr (Dependent) {
        previous = static_cast<std::uint64_t>(bit);
      }
    }
  }
  state.counters["query_state_bytes"] =
      sizeof(counter) + (Dependent ? sizeof(previous) : 0);
  state.counters["queries_per_iteration"] = kQueriesPerIteration;
  state.SetItemsProcessed(state.iterations() * kQueriesPerIteration);
}

enum class Mutation {
  SplitRejoinRandom,
  SplitRejoinBoundary,
  Local,
  Whole,
  Aligned,
  Unaligned
};

template <typename V, Mutation Operation, bool OneLeaf = false>
void TreeMutation(benchmark::State& state) {
  using Block = typename V::block_type;
  const std::size_t n = OneLeaf ? state.range(0) : state.range(0) * kMiB * 8;
  if (!CheckMemoryBudget<V>(state, n)) {
    return;
  }
  auto tree = MakeTree<V>(n);
  std::uint64_t counter = 128;
  state.SetLabel(
      "steady state; generator+allocation+repair included; no reset");
  for (auto _ : state) {
    const auto random = Mix(counter++);
    if constexpr (Operation == Mutation::SplitRejoinRandom ||
                  Operation == Mutation::SplitRejoinBoundary) {
      // Grid boundaries match the initial build; repair can change leaf seams.
      const auto cut =
          Operation == Mutation::SplitRejoinRandom
              ? random % (n + 1)
              : (random % (n / Block::capacity + 1)) * Block::capacity;
      auto suffix = tree.split_off(cut);
      tree.merge(suffix);
    } else if constexpr (Operation == Mutation::Local) {
      const auto left = random % (n - 65 + 1);
      tree.rotate_left(left, left + 65, 1);
    } else if constexpr (Operation == Mutation::Whole) {
      tree.rotate_left(0, n, 1 + random % (n - 1));
    } else if constexpr (Operation == Mutation::Aligned) {
      // Fanout*capacity grid and distance align initial bottom-level subtrees.
      // Subsequent repair need not preserve that geometry; no reset is hidden.
      constexpr auto unit = V::fanout * Block::capacity;
      const auto units = n / unit;
      const auto length = units / 2;
      const auto left = (random % (units - length + 1)) * unit;
      tree.rotate_left(left, left + length * unit,
                       (1 + random % (length - 1)) * unit);
    } else {
      // Large, changing, generally unaligned intervals cover the full sequence.
      const auto length = n / 2 + random % (n / 4);
      const auto left = Mix(random) % (n - length + 1);
      tree.rotate_left(left, left + length, 1 + Mix(random + 1) % (length - 1));
    }
    benchmark::DoNotOptimize(tree);
    benchmark::ClobberMemory();
  }
  if (tree.size() != n) {
    state.SkipWithError("Mutation changed the logical size");
    return;
  }
  TreeMemory<V>(state, tree);
  state.counters["max_external_roots"] =
      Operation == Mutation::SplitRejoinRandom ||
              Operation == Mutation::SplitRejoinBoundary
          ? 2
          : 1;
  state.counters["query_state_bytes"] = sizeof(counter);
  state.SetItemsProcessed(state.iterations());
}

/**
 * @brief Rotate exact child boundaries in freshly built, identical trees.
 * @details Each timed batch mutates 64 independent full F^3-leaf trees once.
 * Level 3 covers the root; levels 1 and 2 cover its first descendant chain.
 * Full selects all children; partial selects children [1,F-1). Construction
 * and destruction are excluded; no evolving geometry or testing hooks occur.
 */
void TreeRotateChildBoundary(benchmark::State& state) {
  using V = P256F8Prefix;
  constexpr auto f = V::fanout;
  constexpr auto capacity = V::block_type::capacity;
  constexpr auto n = f * f * f * capacity;
  constexpr std::size_t batch = 64;
  std::array<V::tree_type, batch> trees;
  std::size_t unit = capacity;
  for (std::int64_t level = 1; level < state.range(0); ++level) {
    unit *= f;
  }
  const auto left = state.range(1) ? 0 : unit;
  const auto right = state.range(1) ? f * unit : (f - 1) * unit;
  state.SetLabel("64 identical trees; one rotation/tree; rebuild excluded");
  for (auto _ : state) {
    state.PauseTiming();
    for (auto& tree : trees) {
      tree = {};
      tree = MakeTree<V>(n);
    }
    state.ResumeTiming();
    for (auto& tree : trees) {
      tree.rotate_left(left, right, 2 * unit);
      benchmark::DoNotOptimize(tree);
    }
    benchmark::ClobberMemory();
  }
  state.counters["operations_per_iteration"] = batch;
  state.SetItemsProcessed(state.iterations() * batch);
}

template <typename V, bool Equal>
void TreeMergeResetSplit(benchmark::State& state) {
  const std::size_t n = state.range(0) * kMiB * 8;
  if (!CheckMemoryBudget<V>(state, n)) {
    return;
  }
  auto tree = MakeTree<V>(n);
  const auto cut = Equal ? n / 2 : n - n / 64;
  state.SetLabel(
      "merge only; reset split excluded, warms paths; final footprint");
  state.counters["donor_bits"] = n - cut;
  // Record actual input heights, not an assumption based on logical lengths.
  auto donor = tree.split_off(cut);
  state.counters["initial_receiver_height"] = tree.height();
  state.counters["initial_donor_height"] = donor.height();
  tree.merge(donor);
  for (auto _ : state) {
    state.PauseTiming();
    donor = tree.split_off(cut);
    state.ResumeTiming();
    tree.merge(donor);
    benchmark::DoNotOptimize(tree);
    benchmark::ClobberMemory();
  }
  if (tree.size() != n || !donor.empty()) {
    state.SkipWithError("Merge failed size/consuming-donor invariant");
    return;
  }
  TreeMemory<V>(state, tree);
  state.counters["max_external_roots"] = 2;
  state.SetItemsProcessed(state.iterations());
}

template <typename V>
void TreeSplitResetMerge(benchmark::State& state) {
  const std::size_t n = state.range(0) * kMiB * 8;
  if (!CheckMemoryBudget<V>(state, n)) {
    return;
  }
  auto tree = MakeTree<V>(n);
  std::uint64_t counter = 128;
  state.SetLabel("split+generator; reset merge excluded, warms paths");
  for (auto _ : state) {
    auto suffix = tree.split_off(Mix(counter++) % (n + 1));
    benchmark::DoNotOptimize(tree);
    benchmark::DoNotOptimize(suffix);
    state.PauseTiming();
    tree.merge(suffix);
    state.ResumeTiming();
  }
  TreeMemory<V>(state, tree);
  state.SetItemsProcessed(state.iterations());
}

template <typename V>
void TreeMergeSmallComponents(benchmark::State& state) {
  using Block = typename V::block_type;
  using Tree = typename V::tree_type;
  const std::size_t n = state.range(0) * 1024 * 8;
  if (!CheckMemoryBudget<V>(state, n)) {
    return;
  }
  state.SetLabel(
      "64-bit component generation+construction+consuming merges; "
      "final destruction excluded");
  bool reported = false;
  for (auto _ : state) {
    Tree tree;
    for (std::size_t first = 0; first < n; first += 64) {
      std::array<Block, 1> input{MakeBlock<Block>(first, 64)};
      auto component = Tree::from_blocks(input);
      tree.merge(component);
    }
    benchmark::DoNotOptimize(tree);
    benchmark::ClobberMemory();
    state.PauseTiming();
    if (!reported) {
      TreeMemory<V>(state, tree);
      reported = true;
    }
    tree = Tree{};
    state.ResumeTiming();
  }
  state.counters["components_per_iteration"] = n / 64;
  state.SetItemsProcessed(state.iterations() * (n / 64));
  state.SetBytesProcessed(state.iterations() * (n / 8));
}

template <typename V, bool SmallComponents>
void TreeBuild(benchmark::State& state) {
  using Tree = typename V::tree_type;
  const std::size_t n = state.range(0) * kMiB * 8;
  if (!CheckMemoryBudget<V>(state, n)) {
    return;
  }
  const auto component_bits = SmallComponents ? 64 : V::block_type::capacity;
  state.SetLabel("streaming input generation+build; destruction excluded");
  state.counters["input_component_bits"] = component_bits;
  state.counters["input_components"] =
      n / component_bits + (n % component_bits != 0);
  bool reported = false;
  for (auto _ : state) {
    auto tree = MakeTree<V>(n, component_bits);
    benchmark::DoNotOptimize(tree);
    benchmark::ClobberMemory();
    state.PauseTiming();
    if (!reported) {
      TreeMemory<V>(state, tree);
      reported = true;
    }
    tree = Tree{};
    state.ResumeTiming();
  }
  state.SetItemsProcessed(state.iterations() * n);
  state.SetBytesProcessed(state.iterations() * (n / 8));
}

template <typename V>
void TreeMaterializeTraversal(benchmark::State& state) {
  using Block = typename V::block_type;
  const std::size_t n = state.range(0) * kMiB * 8;
  const auto words = n / 64 + (n % 64 != 0);
  if (!CheckMemoryBudget<V>(state, n, words * sizeof(std::uint64_t))) {
    return;
  }
  const auto tree = MakeTree<V>(n);
  // This is a second full output buffer, not in-place finalization. Allocate,
  // initialize and pre-touch it outside timing; every valid bit is overwritten
  // during each traversal. No transient allocation/high-water claim is made.
  std::vector<std::uint64_t> output(words, 0);
  benchmark::DoNotOptimize(output.data());
  benchmark::ClobberMemory();
  state.SetLabel(
      "ordered traversal+block flatten+copy; second output buffer; "
      "allocation+initialization+destruction excluded");
  std::size_t offset = 0;
  for (auto _ : state) {
    offset = 0;
    tree.for_each_block([&](const Block& block) {
      const auto flat = block.flatten();
      pixie::copy_packed_bits(flat.data(), 0, output.data(), offset,
                              block.size());
      offset += block.size();
    });
    // Never copy a block's physical padding. Clear final output padding even
    // for a future non-word-sized case, without a redundant full-buffer clear.
    if (n % 64 != 0) {
      output.back() &= (std::uint64_t{1} << (n % 64)) - 1;
    }
    benchmark::DoNotOptimize(output.data());
    benchmark::ClobberMemory();
  }
  if (offset != n) {
    state.SkipWithError("Block traversal did not materialize the logical size");
    return;
  }
  for (std::size_t i = 0; i < words; ++i) {
    auto expected = Mix(i + 42);
    if (i + 1 == words && n % 64 != 0) {
      expected &= (std::uint64_t{1} << (n % 64)) - 1;
    }
    if (output[i] != expected) {
      state.SkipWithError("Materialized contents differ from generated input");
      return;
    }
  }
  TreeMemory<V>(state, tree);
  const auto output_bytes = output.capacity() * sizeof(std::uint64_t);
  state.counters["output_useful_bytes"] = (n + 7) / 8;
  state.counters["output_capacity_bytes"] = output_bytes;
  state.counters["live_output_allocations"] = 1;
  state.counters["tree_plus_output_bytes"] =
      tree.memory_usage_bytes() + output_bytes;
  state.SetItemsProcessed(state.iterations() * n);
  state.SetBytesProcessed(state.iterations() * ((n + 7) / 8));
}

void TreeSizes(benchmark::internal::Benchmark* benchmark) {
  benchmark->ArgName("MiB");
  for (const auto mib : {1, 8, 32, 128, 512}) {
    benchmark->Arg(mib);
  }
}

template <typename V>
void RegisterTreeVariant(const char* variant) {
  const auto add = [=](const char* operation, auto function) {
    const auto name = std::string(operation) + "/" + variant;
    benchmark::RegisterBenchmark(name.c_str(), function)->Apply(TreeSizes);
  };
  add("TreeAccessIndependent", TreeAccess<V, false>);
  add("TreeAccessDependent", TreeAccess<V, true>);
  const auto budget_name =
      std::string("TreeAccessBudgetIndependent/") + variant;
  auto* budget = benchmark::RegisterBenchmark(
      budget_name.c_str(), TreeAccess<V, false, false, true>);
  budget->ArgName("budget_MiB");
  for (const auto mib : {4, 16, 64, 256}) {
    budget->Arg(mib);
  }
  add("TreeSplitResetMerge", TreeSplitResetMerge<V>);
  add("TreeSplitRejoinRandom", TreeMutation<V, Mutation::SplitRejoinRandom>);
  add("TreeSplitRejoinCapacityGrid",
      TreeMutation<V, Mutation::SplitRejoinBoundary>);
  add("TreeRotateLocal65", TreeMutation<V, Mutation::Local>);
  add("TreeRotateWhole", TreeMutation<V, Mutation::Whole>);
  add("TreeRotateGlobalSubtreeGrid", TreeMutation<V, Mutation::Aligned>);
  add("TreeRotateGlobalUnaligned", TreeMutation<V, Mutation::Unaligned>);
  add("TreeMergeEqualResetSplit", TreeMergeResetSplit<V, true>);
  add("TreeMergeUnequalResetSplit", TreeMergeResetSplit<V, false>);
  add("TreeBuild", TreeBuild<V, false>);
}

const bool tree_registered = [] {
  benchmark::RegisterBenchmark(
      "TreeRotateLocal65/P256F8Prefix",
      TreeMutation<P256F8Prefix, Mutation::Local, true>)
      ->ArgName("Bits")
      ->Arg(65)
      ->Arg(128)
      ->Arg(P256F8Prefix::block_type::capacity);
  benchmark::RegisterBenchmark("TreeRotateChildBoundary/P256F8Prefix",
                               TreeRotateChildBoundary)
      ->ArgNames({"Level", "Full"})
      ->Args({1, 0})
      ->Args({1, 1})
      ->Args({2, 0})
      ->Args({2, 1})
      ->Args({3, 0})
      ->Args({3, 1})
      ->Iterations(8);
  RegisterTreeVariant<P256F8Raw>("P256F8Raw");
  RegisterTreeVariant<P256F8Prefix>("P256F8Prefix");
  RegisterTreeVariant<P256F4Prefix>("P256F4Prefix");
  RegisterTreeVariant<P256F4Raw>("P256F4Raw");
  RegisterTreeVariant<P256F16Prefix>("P256F16Prefix");
  RegisterTreeVariant<P256F16Raw>("P256F16Raw");
  RegisterTreeVariant<P128F8Prefix>("P128F8Prefix");
  RegisterTreeVariant<P512F8Prefix>("P512F8Prefix");
  RegisterTreeVariant<DirectMatchedF8Prefix>("DirectMatchedF8Prefix");
  // These expensive/conditional probes are baseline-only for the first screen.
  benchmark::RegisterBenchmark("TreeBuildSmallComponents/P256F8Prefix",
                               TreeBuild<P256F8Prefix, true>)
      ->Apply(TreeSizes);
  benchmark::RegisterBenchmark("TreeMergeSmallComponents/P256F8Prefix",
                               TreeMergeSmallComponents<P256F8Prefix>)
      ->ArgName("KiB")
      ->Arg(4)
      ->Arg(64)
      ->Arg(1024);
  benchmark::RegisterBenchmark("TreeMaterializeTraversal/P256F8Prefix",
                               TreeMaterializeTraversal<P256F8Prefix>)
      ->Apply(TreeSizes);
  benchmark::RegisterBenchmark("TreeNavigationBeyondL3/P256F8Prefix",
                               TreeAccess<P256F8Prefix, false, true>)
      ->ArgName("MiB")
      ->Arg(512);
  return true;
}();

template <typename Block>
void BlockMemory(benchmark::State& state, std::size_t n) {
  state.counters["logical_bits"] = n;
  state.counters["block_capacity_bits"] = Block::capacity;
  state.counters["sizeof_block"] = sizeof(Block);
  state.counters["alignof_block"] = alignof(Block);
  state.counters["payload_capacity_bytes"] = (Block::capacity + 7) / 8;
  state.counters["live_heap_allocations"] = 0;
}

struct BlockRotation {
  std::size_t left, right, distance;
};

// 256 precomputed operations, seed 128. Modes: aligned proper subranges,
// unaligned proper subranges, alternating 50/50, whole block, wrapped partial.
// Random generation and construction excluded; repeated mutation is timed.
// N is capacity (or capacity-3), not a fixed logical size across all types.
// Packed256 and MatchedBitBlock match capacity; the other controls need not.
template <typename Block>
void BlockRotate(benchmark::State& state) {
  const auto mode = state.range(0);
  constexpr std::size_t capacity = Block::capacity;
  const std::size_t n = mode == 4 ? capacity - 3 : capacity;
  auto block = MakeBlock<Block>(0, n);
  std::mt19937_64 random(128);
  std::array<BlockRotation, 256> queries;
  for (std::size_t i = 0; i < queries.size(); ++i) {
    auto& q = queries[i];
    const auto l = random() % 8;
    const auto r = l + 2 + random() % (capacity / 128 - 2 - l);
    q = {l * 128, r * 128, 128 * (1 + random() % (r - l - 1))};
    if (mode == 1 || mode == 4 || (mode == 2 && i % 2)) {
      q = {q.left + 1, q.right - 3, q.distance + 1};
    }
    if (mode == 3) {
      q = {0, n, 65};
    }
  }
  if (mode == 4) {
    block.rotate_left(0, n, n - 7);
  }
  std::size_t i = 0;
  for (auto _ : state) {
    const auto& q = queries[i++ % queries.size()];
    block.rotate_left(q.left, q.right, q.distance);
    benchmark::DoNotOptimize(block);
    benchmark::ClobberMemory();
  }
  BlockMemory<Block>(state, n);
  state.counters["query_pool_bytes"] = sizeof(queries);
  state.SetItemsProcessed(state.iterations());
}

template <typename Block, bool Flatten>
void BlockRead(benchmark::State& state) {
  const std::size_t n = Block::capacity - state.range(0);
  auto block = MakeBlock<Block>(0, n);
  const auto aligned_end = (Block::capacity / 128 - 1) * 128;
  block.rotate_left(128, aligned_end, 384);
  block.rotate_left(256, aligned_end - 128, 640);
  block.rotate_left(0, n, 65);
  // Keep the map/origin as runtime state rather than constant-folding setup.
  benchmark::DoNotOptimize(block);
  std::array<std::size_t, Flatten ? 0 : 1024> queries;
  if constexpr (!Flatten) {
    std::mt19937_64 random(123);
    for (auto& q : queries) {
      q = random() % n;
    }
  }
  std::size_t i = 0;
  for (auto _ : state) {
    if constexpr (Flatten) {
      auto flat = block.flatten();
      benchmark::DoNotOptimize(flat);
      benchmark::ClobberMemory();
    } else {
      benchmark::DoNotOptimize(block[queries[i++ % queries.size()]]);
    }
  }
  BlockMemory<Block>(state, n);
  state.counters["query_pool_bytes"] = Flatten ? 0 : sizeof(queries);
  state.SetItemsProcessed(state.iterations() * (Flatten ? n : 1));
  if constexpr (Flatten) {
    state.SetBytesProcessed(state.iterations() * ((n + 7) / 8));
  }
}

BENCHMARK_TEMPLATE(BlockRotate, DirectBlock)->DenseRange(0, 4);
BENCHMARK_TEMPLATE(BlockRotate, MappedBlock)->DenseRange(0, 4);
BENCHMARK_TEMPLATE(BlockRotate, BitBlock)->DenseRange(0, 4);
BENCHMARK_TEMPLATE(BlockRotate, Packed256)->DenseRange(0, 4);
BENCHMARK_TEMPLATE(BlockRotate, MatchedBitBlock)->DenseRange(0, 4);
BENCHMARK_TEMPLATE(BlockRead, DirectBlock, false)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, MappedBlock, false)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, BitBlock, false)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, Packed256, false)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, MatchedBitBlock, false)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, DirectBlock, true)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);
BENCHMARK_TEMPLATE(BlockRead, MappedBlock, true)
    ->ArgName("unused_bits")
    ->Arg(0)
    ->Arg(3);

}  // namespace
