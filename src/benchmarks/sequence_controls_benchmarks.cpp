#include "sequence_facade_benchmarks.h"

namespace {
const bool registered = [] {
  // At 512 MiB useful uint64_t data the default untagged navigation is
  // plausibly larger than a 16 MiB shared L3, with total requested storage
  // below 1 GiB. The actual node-byte check, not this expectation, decides
  // whether to run.
  benchmark::RegisterBenchmark(
      "FacadeNavigationBeyondL3/Raw64_B256_F8_Prefix",
      Access<UntaggedVariant<std::uint64_t>, false, true>)
      ->ArgName("N")
      ->Arg(EnvironmentMiB("PIXIE_FACADE_NAVIGATION_MIB", 512) /
            sizeof(std::uint64_t));

  RegisterLeaf<bool>("Bool_B256");
  RegisterLeaf<std::uint16_t>("U16_B256");
  RegisterLeaf<std::uint32_t>("U32_B256");
  RegisterLeaf<std::uint64_t>("U64_B256");
  RegisterLeaf<std::uint64_t, 1024>("U64_B128");
  RegisterLeaf<std::uint64_t, 4096>("U64_B512");
  return true;
}();
}  // namespace
