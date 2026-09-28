#include <pixie/experimental/power_of_two_sequence.h>
#include <pixie/experimental/slot_order.h>
#include <pixie/experimental/slot_order16.h>
#include <pixie/permutations/sequence_implementations.h>

#include <random>

#include "sequence_facade_benchmarks.h"

namespace {
using namespace pixie;

// Primitive screen only: excludes tree traversal and subtree-size maintenance.
template <bool Permuted, bool Rotate, class Order = experimental::SlotOrder16>
void SlotOrderProbe(benchmark::State& state) {
  constexpr auto capacity = Order::capacity;
  std::array<std::uint64_t, capacity> payload;
  std::array<const std::uint64_t*, capacity> pointers;
  for (std::size_t i = 0; i < capacity; ++i) {
    // Preserve the original sixteen-slot screen. For wider probes, 5*x+1
    // has a full-period cycle modulo the power-of-two capacity.
    payload[i] =
        capacity == 16 ? (i * 7 + 3) % capacity : (i * 5 + 1) % capacity;
    pointers[i] = &payload[i];
  }
  Order order(capacity);
  struct Query {
    std::size_t left, right, distance;
  };
  std::array<Query, 256> queries;
  std::mt19937 random(917);
  for (auto& query : queries) {
    query.left = random() % capacity;
    query.right = query.left + 1 + random() % (capacity - query.left);
    query.distance = random() % (query.right - query.left);
  }
  std::uint64_t checksum = 0;
  for (auto _ : state) {
    for (const auto& query : queries) {
      if constexpr (Rotate) {
        if constexpr (Permuted) {
          order.rotate_left(query.left, query.right, query.distance);
        } else {
          std::rotate(pointers.begin() + query.left,
                      pointers.begin() + query.left + query.distance,
                      pointers.begin() + query.right);
        }
      }
      const auto index = Rotate ? query.left : checksum % capacity;
      const auto value = [&] {
        if constexpr (Permuted) {
          return *pointers[order[index]];
        } else {
          return *pointers[index];
        }
      }();
      if constexpr (!Rotate && capacity != 16) {
        checksum = value;
      } else {
        checksum += value;
      }
    }
    benchmark::DoNotOptimize(checksum);
  }
  state.SetItemsProcessed(state.iterations() * queries.size());
  state.counters["operations_per_iteration"] = queries.size();
  state.counters["order_object_bytes"] = Permuted ? sizeof(Order) : 0;
}

template <class Order, bool Permuted = true>
void RegisterSlotOrder(const char* prefix) {
  benchmark::RegisterBenchmark((std::string(prefix) + "/RotateRead").c_str(),
                               SlotOrderProbe<Permuted, true, Order>);
  benchmark::RegisterBenchmark((std::string(prefix) + "/DependentRead").c_str(),
                               SlotOrderProbe<Permuted, false, Order>);
}

const bool slot_order_registered = [] {
  benchmark::RegisterBenchmark("SlotOrder16/PointerArray/RotateRead",
                               SlotOrderProbe<false, true>);
  benchmark::RegisterBenchmark("SlotOrder16/Packed/RotateRead",
                               SlotOrderProbe<true, true>);
  benchmark::RegisterBenchmark("SlotOrder16/PointerArray/DependentRead",
                               SlotOrderProbe<false, false>);
  benchmark::RegisterBenchmark("SlotOrder16/Packed/DependentRead",
                               SlotOrderProbe<true, false>);
  RegisterSlotOrder<experimental::SplitSlotOrder<32>, false>(
      "SlotOrder32/PointerArray");
  RegisterSlotOrder<experimental::SplitSlotOrder<32>>(
      "SlotOrder32/SplitScalar");
  RegisterSlotOrder<experimental::ByteSlotOrder<32, false>>(
      "SlotOrder32/ByteScalar");
  RegisterSlotOrder<experimental::ByteSlotOrder<32>>("SlotOrder32/ByteSimd");
  RegisterSlotOrder<experimental::SplitSlotOrder<64>, false>(
      "SlotOrder64/PointerArray");
  RegisterSlotOrder<experimental::SplitSlotOrder<64>>(
      "SlotOrder64/SplitScalar");
  RegisterSlotOrder<experimental::ByteSlotOrder<64, false>>(
      "SlotOrder64/ByteScalar");
  RegisterSlotOrder<experimental::ByteSlotOrder<64>>("SlotOrder64/ByteSimd");
  return true;
}();

template <typename T,
          ElementStorage Storage = ElementStorage::automatic,
          std::size_t Bits = 2048,
          std::size_t F = 8,
          LengthLayout Layout = LengthLayout::cumulative,
          std::size_t ChunkBytes = 4096>
struct ElementVariant : Configuration<T, Bits, F, Layout> {
  using sequence_type =
      PermutableSequence<T, Storage, Bits, F, Layout, ChunkBytes>;
  static constexpr bool rebased = false;
  static constexpr bool indirect =
      Storage == ElementStorage::indirect || !std::is_unsigned_v<T>;
  static constexpr std::size_t chunk_bytes = indirect ? ChunkBytes : 0;
  static auto Make(std::size_t n) {
    // Prvalue elements include move-only values. No full input/pointer vector.
    auto input =
        std::views::iota(std::size_t{0}, n) |
        std::views::transform([](std::size_t i) { return Value<T>(i); });
    return sequence_type::from_range(input);
  }
};

template <std::size_t Values,
          std::size_t F,
          std::size_t ReservedNodes = 0,
          bool CacheRegularPrefix = false,
          bool PermutedLeaves = false>
struct PowerOfTwoVariant
    : Configuration<
          std::uint64_t,
          sizeof(experimental::PowerOfTwoValueBlock<std::uint64_t, Values>) * 8,
          F,
          LengthLayout::cumulative> {
  using block_type = experimental::PowerOfTwoValueBlock<std::uint64_t, Values>;
  using sequence_type =
      experimental::PowerOfTwoPermutableSequence<std::uint64_t,
                                                 Values,
                                                 F,
                                                 ReservedNodes,
                                                 CacheRegularPrefix,
                                                 PermutedLeaves>;
  static constexpr bool rebased = false;
  static constexpr bool indirect = false;
  static constexpr std::size_t chunk_bytes = 0;
  static auto Make(std::size_t n) {
    return sequence_type::from_range(
        std::views::iota(std::uint64_t{0}, std::uint64_t(n)));
  }
};

#ifdef PIXIE_IMMER_SUPPORT
struct ImmerVariant : Configuration<std::uint64_t,
                                    2048,
                                    ImmerPermutableSequence<>::fanout,
                                    LengthLayout::cumulative> {
  using sequence_type = ImmerPermutableSequence<>;
  static constexpr bool rebased = false;
  static constexpr bool indirect = false;
  static constexpr bool immer = true;
  static constexpr std::size_t chunk_bytes = 0;
  static auto Make(std::size_t n) {
    auto input = std::views::iota(std::uint64_t{0}, std::uint64_t(n));
    return sequence_type::from_range(input);
  }
};
#endif

template <typename V>
void RegisterInsertion(const char* variant) {
  const auto name = std::string("FacadeInsertRandomBatch/") + variant;
  auto* row = benchmark::RegisterBenchmark(name.c_str(), InsertBatch<V>);
  Sizes<V>(row);
  row->Iterations(8);
}

template <typename V>
void RegisterEditedReads(const char* variant) {
  const auto independent =
      std::string("FacadeAccessEditedIndependent/") + variant;
  Sizes<V>(benchmark::RegisterBenchmark(independent.c_str(),
                                        Access<V, false, false, true>));
  const auto dependent = std::string("FacadeAccessEditedDependent/") + variant;
  Sizes<V>(benchmark::RegisterBenchmark(dependent.c_str(),
                                        Access<V, true, false, true>));
}

template <class V>
void RegisterChildAligned(const char* variant) {
  const auto name = std::string("FacadeRotateChildAlignedBatch/") + variant;
  auto* row = benchmark::RegisterBenchmark(
      name.c_str(), RotateBatch<V, Rotation::ChildAligned16>);
  Sizes<V>(row);
  row->Iterations(8);
}

template <class V>
void RegisterGroupedSpill(const char* variant) {
  const auto name = std::string("FacadeRotateGroupedSpillBatch/") + variant;
  auto* row = benchmark::RegisterBenchmark(
      name.c_str(), RotateBatch<V, Rotation::GroupedSpill256>);
  Sizes<V>(row);
  row->Iterations(8);
}

template <std::size_t Bits, std::size_t F>
void RegisterPackedConfiguration(const char* name) {
  using V = ElementVariant<std::uint64_t, ElementStorage::packed, Bits, F>;
  Register<V, true>(name);
  RegisterInsertion<V>(name);
  RegisterEditedReads<V>(name);
}

template <class V>
void RegisterWideAligned(const char* variant) {
  const auto add = [&](const char* operation, auto function) {
    const auto name = std::string(operation) + "/" + variant;
    auto* row = benchmark::RegisterBenchmark(name.c_str(), function);
    Sizes<V>(row);
    row->Iterations(8);
  };
  add("FacadeRotateChildAligned32Batch",
      RotateBatch<V, Rotation::ChildAligned32>);
  add("FacadeRotateChildAligned64Batch",
      RotateBatch<V, Rotation::ChildAligned64>);
  const auto name = std::string("FacadeAccessAlignedIndependent/") + variant;
  Sizes<V>(benchmark::RegisterBenchmark(name.c_str(),
                                        Access<V, false, false, false, true>));
}

template <class V>
void RegisterWideConfiguration(const char* variant) {
  Register<V, true>(variant);
  RegisterInsertion<V>(variant);
  RegisterEditedReads<V>(variant);
  RegisterChildAligned<V>(variant);
  RegisterWideAligned<V>(variant);
}

const bool registered = [] {
  RegisterWideConfiguration<PowerOfTwoVariant<256, 32, 96, false, true>>(
      "PowerOfTwo64_L256_F32_R96_Bytes");
  RegisterWideConfiguration<PowerOfTwoVariant<256, 64, 48>>(
      "PowerOfTwo64_L256_F64_R48");
  RegisterWideConfiguration<PowerOfTwoVariant<256, 64, 48, false, true>>(
      "PowerOfTwo64_L256_F64_R48_Bytes");
  RegisterWideAligned<PowerOfTwoVariant<256, 32, 96>>(
      "PowerOfTwo64_L256_F32_R96");
#ifdef PIXIE_IMMER_SUPPORT
  RegisterWideAligned<ImmerVariant>("Immer64_F32_L32");
#endif
  Register<ElementVariant<std::uint64_t, ElementStorage::packed>, true>(
      "Packed64_B256_F8_Prefix");
  RegisterInsertion<ElementVariant<std::uint64_t, ElementStorage::packed>>(
      "Packed64_B256_F8_Prefix");
  RegisterEditedReads<ElementVariant<std::uint64_t, ElementStorage::packed>>(
      "Packed64_B256_F8_Prefix");
  Register<ElementVariant<std::uint64_t, ElementStorage::packed, 2048, 16>,
           true>("Packed64_B256_F16_Prefix");
  RegisterInsertion<
      ElementVariant<std::uint64_t, ElementStorage::packed, 2048, 16>>(
      "Packed64_B256_F16_Prefix");
  Register<ElementVariant<std::uint64_t, ElementStorage::packed, 2048, 32>,
           true>("Packed64_B256_F32_Prefix");
  RegisterInsertion<
      ElementVariant<std::uint64_t, ElementStorage::packed, 2048, 32>>(
      "Packed64_B256_F32_Prefix");
  RegisterPackedConfiguration<16384, 16>("Packed64_B2048_F16_Prefix");
  RegisterPackedConfiguration<16384, 32>("Packed64_B2048_F32_Prefix");
  Register<PowerOfTwoVariant<256, 32>, true>("PowerOfTwo64_L256_F32");
  Register<PowerOfTwoVariant<256, 256, 12>, true>("PowerOfTwo64_L256_F256_R12");
  RegisterInsertion<PowerOfTwoVariant<256, 256, 12>>(
      "PowerOfTwo64_L256_F256_R12");
  RegisterEditedReads<PowerOfTwoVariant<256, 256, 12>>(
      "PowerOfTwo64_L256_F256_R12");
  RegisterChildAligned<PowerOfTwoVariant<256, 256, 12>>(
      "PowerOfTwo64_L256_F256_R12");
  RegisterGroupedSpill<PowerOfTwoVariant<256, 256, 12>>(
      "PowerOfTwo64_L256_F256_R12");
  Register<PowerOfTwoVariant<256, 256, 12, false, true>, true>(
      "PowerOfTwo64_L256_F256_R12_Groups");
  RegisterInsertion<PowerOfTwoVariant<256, 256, 12, false, true>>(
      "PowerOfTwo64_L256_F256_R12_Groups");
  RegisterEditedReads<PowerOfTwoVariant<256, 256, 12, false, true>>(
      "PowerOfTwo64_L256_F256_R12_Groups");
  RegisterChildAligned<PowerOfTwoVariant<256, 256, 12, false, true>>(
      "PowerOfTwo64_L256_F256_R12_Groups");
  RegisterGroupedSpill<PowerOfTwoVariant<256, 256, 12, false, true>>(
      "PowerOfTwo64_L256_F256_R12_Groups");
  RegisterGroupedSpill<PowerOfTwoVariant<256, 16, 96, false, true>>(
      "PowerOfTwo64_L256_F16_R96_Slots");
  RegisterGroupedSpill<PowerOfTwoVariant<256, 32, 96>>(
      "PowerOfTwo64_L256_F32_R96");
#ifdef PIXIE_IMMER_SUPPORT
  RegisterGroupedSpill<ImmerVariant>("Immer64_F32_L32");
#endif
  Register<PowerOfTwoVariant<256, 16, 96>, true>("PowerOfTwo64_L256_F16_R96");
  RegisterInsertion<PowerOfTwoVariant<256, 16, 96>>(
      "PowerOfTwo64_L256_F16_R96");
  RegisterEditedReads<PowerOfTwoVariant<256, 16, 96>>(
      "PowerOfTwo64_L256_F16_R96");
  Register<PowerOfTwoVariant<256, 16, 96, false, true>, true>(
      "PowerOfTwo64_L256_F16_R96_Slots");
  RegisterInsertion<PowerOfTwoVariant<256, 16, 96, false, true>>(
      "PowerOfTwo64_L256_F16_R96_Slots");
  RegisterEditedReads<PowerOfTwoVariant<256, 16, 96, false, true>>(
      "PowerOfTwo64_L256_F16_R96_Slots");
  RegisterChildAligned<PowerOfTwoVariant<256, 16, 96>>(
      "PowerOfTwo64_L256_F16_R96");
  RegisterChildAligned<PowerOfTwoVariant<256, 16, 96, false, true>>(
      "PowerOfTwo64_L256_F16_R96_Slots");
  RegisterChildAligned<PowerOfTwoVariant<256, 32, 96>>(
      "PowerOfTwo64_L256_F32_R96");
#ifdef PIXIE_IMMER_SUPPORT
  RegisterChildAligned<ImmerVariant>("Immer64_F32_L32");
#endif
  RegisterInsertion<PowerOfTwoVariant<256, 32>>("PowerOfTwo64_L256_F32");
  RegisterEditedReads<PowerOfTwoVariant<256, 32>>("PowerOfTwo64_L256_F32");
  Register<PowerOfTwoVariant<256, 32, 0, true>, true>(
      "PowerOfTwo64_L256_F32_Cached");
  RegisterInsertion<PowerOfTwoVariant<256, 32, 0, true>>(
      "PowerOfTwo64_L256_F32_Cached");
  RegisterEditedReads<PowerOfTwoVariant<256, 32, 0, true>>(
      "PowerOfTwo64_L256_F32_Cached");
  Register<PowerOfTwoVariant<256, 32, 96, true>, true>(
      "PowerOfTwo64_L256_F32_R96_Cached");
  RegisterInsertion<PowerOfTwoVariant<256, 32, 96, true>>(
      "PowerOfTwo64_L256_F32_R96_Cached");
  RegisterEditedReads<PowerOfTwoVariant<256, 32, 96, true>>(
      "PowerOfTwo64_L256_F32_R96_Cached");
  Register<PowerOfTwoVariant<256, 32, 96>, true>("PowerOfTwo64_L256_F32_R96");
  RegisterInsertion<PowerOfTwoVariant<256, 32, 96>>(
      "PowerOfTwo64_L256_F32_R96");
  RegisterEditedReads<PowerOfTwoVariant<256, 32, 96>>(
      "PowerOfTwo64_L256_F32_R96");
  // Below the reserve threshold: keep allocation-free retention policy from
  // hiding a regression in ordinary small-tree rotation work.
  const auto small_rotation = []<class V>(const char* name) {
    benchmark::RegisterBenchmark(name, RotateBatch<V, Rotation::Global>)
        ->ArgName("N")
        ->Arg(4096)
        ->Arg(65536)
        ->Iterations(8);
  };
  small_rotation.template operator()<PowerOfTwoVariant<256, 32>>(
      "FacadeRotateGlobalSmallBatch/PowerOfTwo64_L256_F32");
  small_rotation.template operator()<PowerOfTwoVariant<256, 32, 96>>(
      "FacadeRotateGlobalSmallBatch/PowerOfTwo64_L256_F32_R96");
  Register<PowerOfTwoVariant<32, 32>, true>("PowerOfTwo64_L32_F32");
  RegisterInsertion<PowerOfTwoVariant<32, 32>>("PowerOfTwo64_L32_F32");
  RegisterEditedReads<PowerOfTwoVariant<32, 32>>("PowerOfTwo64_L32_F32");
  Register<PowerOfTwoVariant<64, 32>, true>("PowerOfTwo64_L64_F32");
  RegisterInsertion<PowerOfTwoVariant<64, 32>>("PowerOfTwo64_L64_F32");
  RegisterEditedReads<PowerOfTwoVariant<64, 32>>("PowerOfTwo64_L64_F32");
#ifdef PIXIE_IMMER_SUPPORT
  Register<ImmerVariant, true>("Immer64_F32_L32");
  RegisterInsertion<ImmerVariant>("Immer64_F32_L32");
  RegisterEditedReads<ImmerVariant>("Immer64_F32_L32");
  small_rotation.template operator()<ImmerVariant>(
      "FacadeRotateGlobalSmallBatch/Immer64_F32_L32");
#endif
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect>, true>(
      "Indirect64_B256_F8_Prefix_C4096");
  Register<ElementVariant<bool, ElementStorage::packed>>(
      "PackedBool_B256_F8_Prefix");
  Register<ElementVariant<bool, ElementStorage::indirect>, true>(
      "IndirectBool_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<8>>>("Payload8_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<64>>, true>("Payload64_B256_F8_Prefix_C4096");
  Register<ElementVariant<Payload<256>>>("Payload256_B256_F8_Prefix_C4096");
  Register<ElementVariant<MoveOnlyPayload>, true>(
      "MoveOnly64_B256_F8_Prefix_C4096");

  Register<UntaggedVariant<std::uint16_t>>("Raw16_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint32_t>>("Raw32_B256_F8_Prefix");
  Register<UntaggedVariant<std::uint64_t>, true>("Raw64_B256_F8_Prefix");
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect, 2048, 8,
                          LengthLayout::cumulative, 1024>,
           false, true>("Indirect64_B256_F8_Prefix_C1024");
  Register<ElementVariant<std::uint64_t, ElementStorage::indirect, 2048, 8,
                          LengthLayout::cumulative, 16384>,
           false, true>("Indirect64_B256_F8_Prefix_C16384");

  return true;
}();
}  // namespace
