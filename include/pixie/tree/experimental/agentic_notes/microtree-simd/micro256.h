#pragma once

#include <pixie/tree/experimental/tree_micro_prefetch.h>

namespace micro256 {

// Inclusive excess modulo 256. Valid microtree depths are 0..128.
inline void excess(const uint64_t* bits, unsigned length, uint8_t* out) {
  int carry = 0;
#if defined(PIXIE_AVX512_SUPPORT)
  for (unsigned i = 0; i < length; i += 64) {
    const auto w = bits[i / 64];
    auto v =
        _mm512_mask_blend_epi8(w, _mm512_set1_epi8(-1), _mm512_set1_epi8(1));
    v = _mm512_add_epi8(v, _mm512_bslli_epi128(v, 1));
    v = _mm512_add_epi8(v, _mm512_bslli_epi128(v, 2));
    v = _mm512_add_epi8(v, _mm512_bslli_epi128(v, 4));
    v = _mm512_add_epi8(v, _mm512_bslli_epi128(v, 8));
    auto repeat = [](int x) {
      return uint64_t(uint8_t(x)) * 0x0101010101010101ULL;
    };
    auto a = repeat(carry);
    auto b = repeat(carry + 2 * std::popcount(w & 0xffffULL) - 16);
    auto c = repeat(carry + 2 * std::popcount(w & 0xffffffffULL) - 32);
    auto d = repeat(carry + 2 * std::popcount(w & 0xffffffffffffULL) - 48);
    _mm512_storeu_si512(
        out + i, _mm512_add_epi8(v, _mm512_setr_epi64(a, a, b, b, c, c, d, d)));
    carry += 2 * std::popcount(w) - 64;
  }
#elif defined(PIXIE_AVX2_SUPPORT)
  auto shuffle =
      _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                       2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
  auto masks = _mm256_set1_epi64x(0x8040201008040201ULL);
  for (unsigned i = 0; i < length; i += 32) {
    uint32_t w = bits[i / 64] >> (i % 64);
    auto v = _mm256_shuffle_epi8(_mm256_set1_epi32(w), shuffle);
    v = _mm256_cmpeq_epi8(_mm256_and_si256(v, masks), masks);
    v = _mm256_sub_epi8(
        _mm256_or_si256(_mm256_and_si256(v, _mm256_set1_epi8(2)),
                        _mm256_set1_epi8(1)),
        _mm256_set1_epi8(2));
    v = _mm256_add_epi8(v, _mm256_slli_si256(v, 1));
    v = _mm256_add_epi8(v, _mm256_slli_si256(v, 2));
    v = _mm256_add_epi8(v, _mm256_slli_si256(v, 4));
    v = _mm256_add_epi8(v, _mm256_slli_si256(v, 8));
    auto a = uint64_t(uint8_t(carry)) * 0x0101010101010101ULL;
    auto b = uint64_t(uint8_t(carry + 2 * std::popcount(w & 0xffffu) - 16)) *
             0x0101010101010101ULL;
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i),
                        _mm256_add_epi8(v, _mm256_setr_epi64x(a, a, b, b)));
    carry += 2 * std::popcount(w) - 32;
  }
#else
  for (unsigned i = 0; i < length; ++i) {
    carry += (bits[i / 64] >> (i % 64)) & 1 ? 1 : -1;
    out[i] = uint8_t(carry);
  }
#endif
}

inline uint64_t equal64(const uint8_t* p, uint8_t target) {
#if defined(PIXIE_AVX512_SUPPORT)
  return _mm512_cmpeq_epi8_mask(_mm512_loadu_si512(p),
                                _mm512_set1_epi8(target));
#elif defined(PIXIE_AVX2_SUPPORT)
  auto t = _mm256_set1_epi8(target);
  auto low = _mm256_cmpeq_epi8(
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p)), t);
  auto high = _mm256_cmpeq_epi8(
      _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + 32)), t);
  return uint32_t(_mm256_movemask_epi8(low)) |
         (uint64_t(uint32_t(_mm256_movemask_epi8(high))) << 32);
#else
  uint64_t mask = 0;
  for (unsigned i = 0; i < 64; ++i) {
    if (p[i] == target) {
      mask |= uint64_t{1} << i;
    }
  }
  return mask;
#endif
}

inline unsigned matching(const uint8_t* depths, unsigned p, unsigned length) {
  for (unsigned w = p / 64; w < (length + 63) / 64; ++w) {
    auto mask = equal64(depths + w * 64, depths[p] - 1);
    if (w == p / 64) {
      mask &= ~uint64_t{0} << (p % 64);
    }
    if (mask) {
      return w * 64 + std::countr_zero(mask);
    }
  }
  std::abort();
}

inline unsigned open_positions(const uint64_t* bits,
                               unsigned length,
                               uint8_t* positions) {
  unsigned n = 0;
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
  alignas(64) static constexpr auto indices = [] {
    std::array<uint8_t, 64> a{};
    for (unsigned i = 0; i < 64; ++i) {
      a[i] = i;
    }
    return a;
  }();
#endif
  for (unsigned w = 0; w < (length + 63) / 64; ++w) {
    auto mask = bits[w];
    if (w == length / 64 && length % 64) {
      mask &= (uint64_t{1} << (length % 64)) - 1;
    }
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
    auto v = _mm512_add_epi8(_mm512_load_si512(indices.data()),
                             _mm512_set1_epi8(w * 64));
    _mm512_mask_compressstoreu_epi8(positions + n, mask, v);
    n += std::popcount(mask);
#else
    while (mask) {
      positions[n++] = w * 64 + std::countr_zero(mask);
      mask &= mask - 1;
    }
#endif
  }
  return n;
}

template <bool RegisterStore>
inline unsigned append_positions(uint64_t mask, unsigned base, uint8_t* out) {
#if defined(PIXIE_AVX512_SUPPORT) && defined(__AVX512VBMI2__)
  alignas(64) static constexpr auto indices = [] {
    std::array<uint8_t, 64> a{};
    for (unsigned i = 0; i < 64; ++i) {
      a[i] = i;
    }
    return a;
  }();
  auto v = _mm512_add_epi8(_mm512_load_si512(indices.data()),
                           _mm512_set1_epi8(base));
  if constexpr (RegisterStore) {
    _mm512_storeu_si512(out, _mm512_maskz_compress_epi8(mask, v));
  } else {
    _mm512_mask_compressstoreu_epi8(out, mask, v);
  }
  return std::popcount(mask);
#else
  unsigned count = 0;
  while (mask) {
    out[count++] = base + std::countr_zero(mask);
    mask &= mask - 1;
  }
  return count;
#endif
}

// Mode 0: scalar stack. 1: SIMD matching, materialized links.
// Mode 2: SIMD depth masks, direct replay of compressed opening positions.
template <unsigned Mode>
class Block {
  std::array<uint64_t, 4> bits_{};
  alignas(64) std::array<uint8_t, 256> depths_{};
  std::array<uint8_t, 256> close_{}, number_{}, opener_{};
  unsigned length_ = 0;
  pixie::TreeNode root_{0, 0};
  bool bit(unsigned p) const { return (bits_[p / 64] >> (p % 64)) & 1; }
  pixie::TreeNode node(unsigned p) const {
    return {root_.number + number_[p], root_.pos + p};
  }

 public:
  bool load(std::span<const uint64_t> words,
            size_t nbits,
            pixie::TreeNode root) {
    using namespace pixie::detail::tree_unpack;
    root_ = root;
    unsigned available = std::min<size_t>(256, nbits - root.pos);
    for (unsigned i = 0; i < 4; ++i) {
      bits_[i] = word_at(words, root.pos + i * 64);
    }
    length_ = 0;
    int target = 0;
    for (unsigned i = 0; i < available; i += 128) {
      int delta;
      auto found = forward_search_128(bits_.data() + i / 64, target, 0, &delta);
      if (found < std::min(128u, available - i)) {
        length_ = i + found + 1;
        break;
      }
      target -= delta;
    }
    if (!length_) {
      return false;
    }
    if constexpr (Mode == 0) {
      unsigned depth = 0, rank = 0;
      uint8_t stack[128];
      for (unsigned p = 0; p < length_; ++p) {
        if (bit(p)) {
          number_[p] = rank++;
          stack[depth++] = p;
        } else {
          auto open = stack[--depth];
          close_[open] = p;
          opener_[p] = open;
        }
      }
    } else {
      excess(bits_.data(), length_, depths_.data());
      if constexpr (Mode >= 5) {
        // At a fixed depth, opening and closing parentheses have identical
        // left-to-right order. Pair whole-level masks, rather than searching
        // the depth vector separately for every opening.
        unsigned height = 0;
        for (unsigned p = 0; p < length_; ++p) {
          height = std::max(height, unsigned(depths_[p]));
        }
        if (Mode >= 6 && height > 16) {
          unsigned depth = 0, rank = 0;
          uint8_t stack[128];
          for (unsigned p = 0; p < length_; ++p) {
            if (bit(p)) {
              number_[p] = rank++;
              stack[depth++] = p;
            } else {
              auto open = stack[--depth];
              close_[open] = p;
              opener_[p] = open;
            }
          }
        } else {
          for (unsigned d = 1; d <= height; ++d) {
            // A full register store may write up to 63 unused trailing bytes.
            uint8_t opens[192], closes[192];
            unsigned no = 0, nc = 0;
            for (unsigned w = 0; w < (length_ + 63) / 64; ++w) {
              uint64_t valid = ~uint64_t{0};
              if (w == length_ / 64 && length_ % 64) {
                valid = (uint64_t{1} << (length_ % 64)) - 1;
              }
              auto o = equal64(depths_.data() + w * 64, d) & bits_[w] & valid;
              auto c =
                  equal64(depths_.data() + w * 64, d - 1) & ~bits_[w] & valid;
              if constexpr (Mode >= 6) {
                no += append_positions<Mode == 7>(o, w * 64, opens + no);
                nc += append_positions<Mode == 7>(c, w * 64, closes + nc);
              } else {
                while (o) {
                  opens[no++] = w * 64 + std::countr_zero(o);
                  o &= o - 1;
                }
                while (c) {
                  closes[nc++] = w * 64 + std::countr_zero(c);
                  c &= c - 1;
                }
              }
            }
            assert(no == nc);
            for (unsigned k = 0; k < no; ++k) {
              auto p = opens[k], q = closes[k];
              close_[p] = q;
              opener_[q] = p;
              number_[p] = (unsigned(p) + d - 1) / 2;
            }
          }
        }
      } else if constexpr (Mode == 1) {
        uint8_t positions[128];
        unsigned n = open_positions(bits_.data(), length_, positions);
        for (unsigned i = 0; i < n; ++i) {
          auto p = positions[i];
          auto q = matching(depths_.data(), p, length_);
          close_[p] = q;
          opener_[q] = p;
          number_[p] = i;
        }
      }
    }
    return true;
  }
  unsigned size() const { return length_ / 2; }
  template <class V, class I>
  void walk(V& visit, I& inspect) const {
    using pixie::detail::tree_unpack::emit;
    if constexpr (Mode >= 2 && Mode <= 4) {
      uint8_t stack[128], positions[128];
      uint64_t depth_masks[Mode == 4 ? 129 : 1][4];
      uint64_t valid_depths[3]{};
      auto masks_for = [&](unsigned depth) {
        const uint64_t flag = uint64_t{1} << (depth % 64);
        if (!(valid_depths[depth / 64] & flag)) {
          for (unsigned w = 0; w < (length_ + 63) / 64; ++w) {
            depth_masks[depth][w] = equal64(depths_.data() + w * 64, depth);
          }
          valid_depths[depth / 64] |= flag;
        }
        return depth_masks[depth];
      };
      unsigned n = open_positions(bits_.data(), length_, positions);
      auto handle = [&](unsigned p) {
        return pixie::TreeNode(
            root_.number + (p + unsigned(depths_[p]) - 1) / 2, root_.pos + p);
      };
      for (unsigned i = 0; i < n; ++i) {
        unsigned p = positions[i], d = depths_[p];
        stack[d - 1] = p;
        emit(visit, pixie::TreeNode(root_.number + i, root_.pos + p), true);
        if (bit(p + 1)) {
          if constexpr (Mode == 3 || Mode == 4) {
            const uint64_t* closing_masks = nullptr;
            const uint64_t* child_masks = nullptr;
            if constexpr (Mode == 4) {
              closing_masks = masks_for(d - 1);
              child_masks = masks_for(d + 1);
            }
            for (unsigned w = p / 64; w < (length_ + 63) / 64; ++w) {
              auto closes = Mode == 4 ? closing_masks[w]
                                      : equal64(depths_.data() + w * 64, d - 1);
              auto children =
                  (Mode == 4 ? child_masks[w]
                             : equal64(depths_.data() + w * 64, d + 1)) &
                  bits_[w];
              if (w == p / 64) {
                auto keep = ~uint64_t{0} << (p % 64);
                closes &= keep;
                children &= keep;
              }
              if (closes) {
                children &= (uint64_t{1} << std::countr_zero(closes)) - 1;
              }
              while (children) {
                inspect(handle(w * 64 + std::countr_zero(children)));
                children &= children - 1;
              }
              if (closes) {
                break;
              }
            }
          } else {
            unsigned q = matching(depths_.data(), p, length_);
            for (unsigned w = p / 64; w <= q / 64; ++w) {
              auto children =
                  equal64(depths_.data() + w * 64, d + 1) & bits_[w];
              if (w == p / 64) {
                children &= ~uint64_t{0} << (p % 64);
              }
              if (w == q / 64) {
                children &= (uint64_t{1} << (q % 64)) - 1;
              }
              while (children) {
                inspect(handle(w * 64 + std::countr_zero(children)));
                children &= children - 1;
              }
            }
          }
        } else {
          unsigned next_depth = i + 1 < n ? depths_[positions[i + 1]] : 1;
          while (d > next_depth) {
            --d;
            emit(visit, handle(stack[d - 1]), false);
          }
        }
      }
    } else {
      for (unsigned p = 0; p < length_; ++p) {
        if (bit(p)) {
          emit(visit, node(p), true);
          for (unsigned c = p + 1; c < close_[p]; c = unsigned(close_[c]) + 1) {
            inspect(node(c));
          }
        } else if (p != unsigned(opener_[p]) + 1) {
          emit(visit, node(opener_[p]), false);
        }
      }
    }
  }
};

template <unsigned Mode, class Tree, class V, class I>
class Traversal {
  const Tree& tree_;
  std::span<const uint64_t> words_;
  V& visit_;
  I& inspect_;
  Block<Mode> block_;

 public:
  size_t nodes = 0, children = 0;
  Traversal(const Tree& t, std::span<const uint64_t> w, V& v, I& i)
      : tree_(t), words_(w), visit_(v), inspect_(i) {}
  void walk(pixie::TreeNode node, size_t depth = 0) {
    using namespace pixie::detail::tree_unpack;
    if (tree_.is_leaf(node)) {
      ++nodes;
      emit(visit_, node, true);
      return;
    }
    if (block_.load(words_, 2 * tree_.size(), node)) {
      block_.walk(visit_, inspect_);
      nodes += block_.size();
      children += block_.size() - 1;
      return;
    }
    if (depth == 128) {
      auto events = [&](pixie::TreeNode n, bool entry) {
        emit(visit_, n, entry);
        if (entry) {
          ++nodes;
          tree_.for_each_child(n, [&](auto c) {
            ++children;
            inspect_(c);
          });
        }
      };
      subtree(tree_, node, events);
      return;
    }
    ++nodes;
    emit(visit_, node, true);
    std::array<size_t, 32> numbers, positions;
    size_t count = 0;
    tree_.for_each_child(node, [&](auto c) {
      inspect_(c);
      if (count < 32) {
        numbers[count] = c.number;
        positions[count] = c.pos;
      }
      ++count;
    });
    children += count;
    if (count <= 32) {
      for (size_t i = 0; i < count; ++i) {
        walk({numbers[i], positions[i]}, depth + 1);
      }
    } else {
      tree_.for_each_child(node, [&](auto c) { walk(c, depth + 1); });
    }
    emit(visit_, node, false);
  }
};
}  // namespace micro256
