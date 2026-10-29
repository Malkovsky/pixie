#ifndef MICRO_CHECK_ONLY
#include <benchmark/benchmark.h>
#endif
#include <pixie/rmm/btree.h>
#include <pixie/tree/bp.h>
#include <time.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <tuple>

#include "micro256.h"
#include "src/benchmarks/tree_input.h"

using Tree = pixie::BPTree<pixie::experimental::RmMBTree<>>;
template <unsigned Mode, class V, class I>
void run(const Tree& tree,
         std::span<const uint64_t> words,
         V& visit,
         I& inspect) {
  if constexpr (Mode == 0 || Mode == 4) {
    pixie::experimental::detail::BpMicroTraversal<Mode == 0 ? 8192 : 256, Tree,
                                                  V, I, 8>
        traversal(tree, words, visit, inspect);
    traversal.walk(tree.root());
    if (traversal.nodes != tree.size() ||
        traversal.children != tree.size() - 1) {
      std::abort();
    }
  } else {
    micro256::Traversal<Mode >= 5 ? Mode - 2 : Mode - 1, Tree, V, I> traversal(
        tree, words, visit, inspect);
    traversal.walk(tree.root());
    if (traversal.nodes != tree.size() ||
        traversal.children != tree.size() - 1) {
      std::abort();
    }
  }
}

void check(std::vector<uint64_t> words, size_t n) {
  Tree tree(words, n);
  using Event = std::tuple<size_t, size_t, unsigned>;
  std::vector<Event> expected, actual;
  auto visit = [&](auto node, bool entry) {
    actual.emplace_back(node.number, node.pos, entry ? 0 : 1);
  };
  auto inspect = [&](auto node) {
    actual.emplace_back(node.number, node.pos, 2);
  };
  auto current = tree.root();
  bool entry = true;
  while (true) {
    visit(current, entry);
    if (entry) {
      tree.for_each_child(current, inspect);
      if (!tree.is_leaf(current)) {
        current = tree.first_child(current);
        continue;
      }
    }
    entry = false;
    if (tree.is_root(current)) {
      break;
    }
    if (tree.is_last_child(current)) {
      current = tree.parent(current);
    } else {
      current = tree.next_sibling(current);
      entry = true;
    }
  }
  expected.swap(actual);
  auto candidate = [&]<unsigned M>() {
    actual.clear();
    run<M>(tree, words, visit, inspect);
    if (actual != expected) {
      std::fprintf(stderr, "trace mismatch n=%zu mode=%u\n", n, M);
      std::abort();
    }
  };
#ifndef MICRO_CHECK_ONLY
  candidate.template operator()<0>();
  candidate.template operator()<4>();
#endif
  candidate.template operator()<1>();
  candidate.template operator()<2>();
  candidate.template operator()<3>();
  candidate.template operator()<5>();
  candidate.template operator()<6>();
  candidate.template operator()<7>();
  candidate.template operator()<8>();
  candidate.template operator()<9>();
}

double now() {
  timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return t.tv_sec * 1e9 + t.tv_nsec;
}
int main(int argc, char** argv) {
  if (argc == 1) {
    unsigned shapes = 0;
    for (unsigned n = 1; n <= 9; ++n) {
      for (uint64_t bits = 1; bits < (uint64_t{1} << (2 * n)); bits += 2) {
        int balance = 0;
        bool valid = true;
        for (unsigned p = 0; p < 2 * n; ++p) {
          balance += (bits >> p) & 1 ? 1 : -1;
          if (balance <= 0 && p + 1 != 2 * n) {
            valid = false;
            break;
          }
        }
        if (valid && balance == 0) {
          check({bits}, n);
          ++shapes;
        }
      }
    }
    for (size_t n : {31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257, 4095,
                     4096, 4097, 16385}) {
      for (unsigned seed = 0; seed < (n <= 129 ? 32u : 2u); ++seed) {
        auto w = pixie::benchmarks::GenerateRandomBp(n, seed + 42);
        check(w, n);
        if (2 * n % 64) {
          w.back() |= ~uint64_t{0} << (2 * n % 64);
        }
        check(w, n);
      }
      for (unsigned shape = 0; shape < 2; ++shape) {
        std::vector<uint64_t> w((2 * n + 63) / 64, 0);
        for (size_t k = 0; k < n; ++k) {
          size_t p = shape ? k : (k ? 2 * k - 1 : 0);
          w[p / 64] |= uint64_t{1} << (p % 64);
        }
        check(w, n);
        if (2 * n % 64) {
          w.back() |= ~uint64_t{0} << (2 * n % 64);
        }
        check(w, n);
      }
    }
    std::printf("%u shapes and boundary/random exact traces passed\n", shapes);
    return 0;
  }
#ifndef MICRO_CHECK_ONLY
  size_t n = std::strtoull(argv[1], nullptr, 10);
  unsigned blocks = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 12;
  auto words = pixie::benchmarks::GenerateRandomBp(n);
  Tree tree(words, n);
  auto observe = [](auto node) { benchmark::DoNotOptimize(node); };
  auto call = [&](unsigned mode) {
    switch (mode) {
      case 0:
        run<0>(tree, words, observe, observe);
        break;
      case 1:
        run<1>(tree, words, observe, observe);
        break;
      case 2:
        run<2>(tree, words, observe, observe);
        break;
      case 3:
        run<3>(tree, words, observe, observe);
        break;
      case 4:
        run<4>(tree, words, observe, observe);
        break;
      case 5:
        run<5>(tree, words, observe, observe);
        break;
      case 6:
        run<6>(tree, words, observe, observe);
        break;
      case 7:
        run<7>(tree, words, observe, observe);
        break;
      case 8:
        run<8>(tree, words, observe, observe);
        break;
      case 9:
        run<9>(tree, words, observe, observe);
        break;
    }
  };
  if (argc > 3) {
    unsigned mode = std::strtoul(argv[3], nullptr, 10);
    call(mode);
    std::ofstream control;
    std::ifstream ack;
    if (auto path = std::getenv("MICRO_PERF_CONTROL")) {
      control.open(path);
      ack.open(std::getenv("MICRO_PERF_ACK"));
      control << "enable\n" << std::flush;
      std::string line;
      std::getline(ack, line);
    }
    for (unsigned b = 0; b < blocks; ++b) {
      call(mode);
    }
    if (control.is_open()) {
      control << "disable\n" << std::flush;
      std::string line;
      std::getline(ack, line);
    }
    return 0;
  }
  for (unsigned m = 0; m < 10; ++m) {
    call(m);
  }
  std::puts("block,mode,nodes,cpu_ns");
  for (unsigned b = 0; b < blocks; ++b) {
    for (unsigned s = 0; s < 10; ++s) {
      unsigned m = (b + (b % 2 ? 9 - s : s)) % 10;
      auto start = now();
      call(m);
      std::printf("%u,%u,%zu,%.0f\n", b, m, n, now() - start);
      std::fflush(stdout);
    }
  }
#endif
}
