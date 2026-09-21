#include <benchmark/benchmark.h>
#include <pixie/detail/sequence/sequence_tree_kernels.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>

namespace {

using SelectFunction = std::size_t (*)(std::span<const std::size_t>,
                                       std::size_t) noexcept;

void NodeSelect(benchmark::State& state, SelectFunction select) {
  const auto count = static_cast<std::size_t>(state.range(1));
  struct alignas(64) Node {
    std::array<std::size_t, 15> ends;
  };
  struct Query {
    std::size_t node;
    std::size_t index;
  };
  // Shared by both function-pointer controls, not separate template statics:
  // matching rows use identical input addresses and deterministic query pools.
  static std::array<Node, 64> nodes;
  static std::array<Query, 512> queries;
  constexpr auto max = std::numeric_limits<std::size_t>::max();
  constexpr auto top = std::size_t{1}
                       << (std::numeric_limits<std::size_t>::digits - 1);
  std::mt19937_64 random(0x5e1ec7);
  for (std::size_t n = 0; n < nodes.size(); ++n) {
    auto end = n % 3 == 0 ? 0 : n % 3 == 1 ? top - 256 : max - 4096;
    for (std::size_t i = 0; i < count; ++i) {
      end += 1 + random() % 128;
      nodes[n].ends[i] = end;
    }
  }
  for (std::size_t i = 0; i < queries.size(); ++i) {
    auto& q = queries[i];
    q.node = random() % nodes.size();
    const auto& ends = nodes[q.node].ends;
    const auto child = (i / 3) % (count + 1);
    // Cycle all child positions, including the implicit last child, with
    // before/equal/after endpoints in low, sign-crossing and near-SIZE_MAX
    // data.
    q.index = child == count ? ends[count - 1] + 1 + random() % 64
                             : ends[child] - 1 + i % 3;
    if (i % 64 == 62) {
      q.index = 0;
    } else if (i % 64 == 63) {
      q.index = max;
    }
  }
  for (const auto& q : queries) {
    const std::span<const std::size_t> ends(nodes[q.node].ends.data(), count);
    const auto expected = static_cast<std::size_t>(
        std::upper_bound(ends.begin(), ends.end(), q.index) - ends.begin());
    if (select(ends, q.index) != expected) {
      state.SkipWithError("Node selection disagrees with upper_bound");
      return;
    }
  }

  constexpr std::int64_t batch = 64;
  std::size_t cursor = 0;
  for (auto _ : state) {
    for (std::int64_t i = 0; i < batch; ++i) {
      const auto& q = queries[cursor++ % queries.size()];
      auto result = select({nodes[q.node].ends.data(), count}, q.index);
      benchmark::DoNotOptimize(result);
    }
  }
  state.SetItemsProcessed(state.iterations() * batch);
  state.SetLabel(
      "hot kernel; 64 queries/iteration; timed: query/node lookup, "
      "indirect call, result barrier; excludes setup/validation");
}

const bool node_select_registered = [] {
  const auto add = [](const char* name, SelectFunction select) {
    auto* benchmark = benchmark::RegisterBenchmark(name, NodeSelect, select);
    benchmark->ArgNames({"fanout", "ends"});
    for (const int fanout : {4, 8, 16}) {
      for (const int count : {fanout / 2 - 1, fanout - 2, fanout - 1}) {
        benchmark->Args({fanout, count});
      }
    }
  };
  add("NodeSelect/Scalar", pixie::detail::sequence::node_select_scalar);
  add("NodeSelect/Dispatched", pixie::detail::sequence::node_select);
  return true;
}();

}  // namespace
