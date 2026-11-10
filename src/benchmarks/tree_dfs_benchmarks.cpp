#include "tree_benchmark_registration.h"

namespace {
const bool registered = [] {
  using pixie::benchmarks::RegisterTreeWorkload;
  using pixie::benchmarks::TreeWorkload;
  RegisterTreeWorkload<TreeWorkload::kDfs>("dfs");
  RegisterTreeWorkload<TreeWorkload::kDfsChildren>("dfs_children");
  return true;
}();
}  // namespace
