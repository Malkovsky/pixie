#include "tree_benchmark_registration.h"

namespace {
const bool registered = [] {
  using pixie::benchmarks::RegisterTreeWorkload;
  using pixie::benchmarks::TreeWorkload;
  RegisterTreeWorkload<TreeWorkload::kRandomPath>("random_path");
  RegisterTreeWorkload<TreeWorkload::kRandomPathChildren>(
      "random_path_children");
  return true;
}();
}  // namespace
