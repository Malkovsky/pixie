#pragma once

#include <pixie/detail/sequence/packed_value_block.h>
#include <pixie/detail/sequence/pointer_block.h>
#include <pixie/detail/sequence/sequence_tree.h>
#include <pixie/permutable_sequence.h>

#include <algorithm>
#include <array>
#include <climits>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <ranges>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie {

/**
 * @brief Owning, move-only tree sequence with immutable element reads.
 * @details Packed storage returns T by value. Indirect storage returns const T&
 * and requires nonthrowing move construction and destruction, but neither
 * default construction nor move assignment. Signed integers initially use
 * indirect storage. Indirect rotations and concatenating merges never move
 * payload objects after construction. References to indirect elements survive
 * both, including references obtained from a consumed donor. They expire upon
 * destruction of the eventual owner or replacement by move assignment.
 *
 * Indirect payload lives in independent frozen vectors of actual object slots
 * (including actual bool objects). A singly linked owning chain is spliced in
 * O(1) after a successful order-tree merge and destroyed iteratively. Repeated
 * small merges deliberately retain small chunks and their unused capacity;
 * there is no compaction, contiguous-data promise, or reference rebasing.
 * Access, rotation and merge use O(Fanout*height) tree work plus bounded local
 * block operations. Construction, destruction and explicit accounting may be
 * linear. Single-leaf and eligible complete-child rotations do not allocate;
 * general rotations use transactional split/join with bounded leaf work.
 * Default order blocks occupy 256 bytes with 1920 payload bits; untagged F8
 * nodes occupy 128 bytes on 64-bit hosts. Dynamic metadata is not promised
 * strictly succinct. Internal detail types are unsupported.
 * No public split, writable access, or serialization is provided.
 * @tparam T Unqualified object type; see storage requirements above.
 * @tparam Storage Compile-time selection, never a runtime switch.
 * @tparam StorageBits Complete local block budget, not whole-tree memory.
 * @tparam Fanout Even tree fanout, at least four.
 * @tparam Layout Representation of tree child lengths.
 * @tparam ChunkBytes Indirect target bytes per vector, at least one T per
 * chunk.
 */
template <class T,
          ElementStorage Storage = ElementStorage::automatic,
          std::size_t StorageBits = 2048,
          std::size_t Fanout = 8,
          LengthLayout Layout = LengthLayout::cumulative,
          std::size_t ChunkBytes = 4096>
class PermutableSequence
    : public PermutableSequenceBase<
          PermutableSequence<T,
                             Storage,
                             StorageBits,
                             Fanout,
                             Layout,
                             ChunkBytes>,
          T,
          std::conditional_t<Storage == ElementStorage::packed ||
                                 (Storage == ElementStorage::automatic &&
                                  std::is_integral_v<T> &&
                                  std::is_unsigned_v<T> &&
                                  std::numeric_limits<T>::digits <= 64),
                             T,
                             const T&>> {
  static_assert(std::is_object_v<T> && !std::is_array_v<T> &&
                std::same_as<T, std::remove_cv_t<T>>);
  static constexpr bool packable = std::is_integral_v<T> &&
                                   std::is_unsigned_v<T> &&
                                   std::numeric_limits<T>::digits <= 64;
  static_assert(Storage == ElementStorage::automatic ||
                Storage == ElementStorage::packed ||
                Storage == ElementStorage::indirect);
  static_assert(Storage != ElementStorage::packed || packable,
                "Packed elements must be bool or unsigned integers <=64 bits");

 public:
  /**
   * @brief Stored element type.
   * @details Identical to T, irrespective of physical representation.
   */
  using value_type = T;
  /**
   * @brief Resolved compile-time storage strategy.
   * @details Always packed or indirect; automatic is resolved from T.
   */
  static constexpr ElementStorage storage =
      Storage == ElementStorage::automatic
          ? (packable ? ElementStorage::packed : ElementStorage::indirect)
          : Storage;
  /**
   * @brief Immutable indexed result type.
   * @details T by value when packed; const T& with stable indirect ownership.
   */
  using const_reference =
      std::conditional_t<storage == ElementStorage::packed, T, const T&>;

 private:
  friend class PermutableSequenceBase<PermutableSequence, T, const_reference>;
  static constexpr bool packed = storage == ElementStorage::packed;
  static_assert(
      packed || (std::is_nothrow_move_constructible_v<T> &&
                 std::is_nothrow_destructible_v<T>),
      "Indirect elements require noexcept move construction/destruction");
  using Block =
      std::conditional_t<packed,
                         detail::sequence::PackedValueBlock<T, StorageBits>,
                         detail::sequence::PointerBlock<T, StorageBits>>;
  static_assert(sizeof(Block) * CHAR_BIT == StorageBits);
  using Tree = detail::sequence::SequenceTree<Block, Fanout, Layout>;
  static constexpr std::size_t chunk_capacity =
      std::max(std::size_t{1}, ChunkBytes / sizeof(T));

  struct BoolSlot {
    bool value;
  };
  static_assert(std::is_trivial_v<BoolSlot>);
  using Slot = std::conditional_t<std::same_as<T, bool>, BoolSlot, T>;
  static_assert(sizeof(Slot) == sizeof(T));
  struct Chunk {
    std::vector<Slot> values;
    Chunk* next = nullptr;
  };
  struct Chunks {
    Chunk* head = nullptr;
    Chunk* tail = nullptr;
    Chunks() noexcept = default;
    Chunks(const Chunks&) = delete;
    Chunks& operator=(const Chunks&) = delete;
    Chunks(Chunks&& other) noexcept
        : head(std::exchange(other.head, nullptr)),
          tail(std::exchange(other.tail, nullptr)) {}
    Chunks& operator=(Chunks&& other) noexcept {
      if (this != &other) {
        clear();
        head = std::exchange(other.head, nullptr);
        tail = std::exchange(other.tail, nullptr);
      }
      return *this;
    }
    ~Chunks() { clear(); }
    void clear() noexcept {
      while (head) {
        auto* next = head->next;
        delete head;
        head = next;
      }
      tail = nullptr;
    }
    void append(Chunk* chunk) noexcept {
      if (tail) {
        tail->next = chunk;
      } else {
        head = chunk;
      }
      tail = chunk;
    }
    void splice(Chunks& donor) noexcept {
      if (!donor.head) {
        return;
      }
      if (tail) {
        tail->next = donor.head;
      } else {
        head = donor.head;
      }
      tail = donor.tail;
      donor.head = donor.tail = nullptr;
    }
  };
  struct NoChunks {};
  // Declaration order destroys the pointer tree before the pointed-to objects.
  [[no_unique_address]] std::conditional_t<packed, NoChunks, Chunks> chunks_;
  Tree tree_;

#ifdef PIXIE_SEQUENCE_TREE_TESTING
  inline static thread_local std::ptrdiff_t payload_failure_ = -1;
#endif
  static void payload_allocation() {
#ifdef PIXIE_SEQUENCE_TREE_TESTING
    if (payload_failure_ == 0) {
      throw std::bad_alloc();
    }
    if (payload_failure_ > 0) {
      --payload_failure_;
    }
#endif
  }
  static std::size_t add(std::size_t a, std::size_t b) noexcept {
    const auto max = std::numeric_limits<std::size_t>::max();
    return b > max - a ? max : a + b;
  }
  static std::size_t multiply(std::size_t a, std::size_t b) noexcept {
    const auto max = std::numeric_limits<std::size_t>::max();
    return b && a > max / b ? max : a * b;
  }

 public:
  /**
   * @brief Construct canonical empty.
   * @details Performs no allocation and owns no payload or tree nodes.
   */
  PermutableSequence() noexcept = default;
  /**
   * @brief Exclusive ownership disallows copying.
   * @details This deleted operation cannot duplicate payload ownership.
   */
  PermutableSequence(const PermutableSequence&) = delete;
  /**
   * @brief Exclusive ownership disallows copy assignment.
   * @details Transfer ownership with move assignment instead.
   */
  PermutableSequence& operator=(const PermutableSequence&) = delete;
  /**
   * @brief Transfer ownership without moving elements.
   * @details The source becomes canonical empty, without allocation. Existing
   * indirect references remain valid under the new owner's lifetime.
   * @param other Source whose ownership is consumed.
   */
  PermutableSequence(PermutableSequence&& other) noexcept = default;
  /**
   * @brief Replace ownership without moving elements; source becomes empty.
   * @details Self-move is a no-op. References to replaced elements expire.
   * References to transferred indirect elements remain valid. Reclaims the old
   * payload without recursive chunk-chain destruction; performs no allocation.
   * @param other Source whose ownership is consumed unless it is this object.
   * @return This sequence after ownership replacement.
   */
  PermutableSequence& operator=(PermutableSequence&& other) noexcept {
    if (this != &other) {
      tree_ = std::move(other.tree_);
      chunks_ = std::move(other.chunks_);
    }
    return *this;
  }
  /**
   * @brief Destroy order nodes and iteratively reclaim all owned chunks.
   * @details Invalidates references to owned elements; does not throw.
   */
  ~PermutableSequence() = default;

  /**
   * @brief Consume a possibly single-pass range through ranges::iter_move.
   * @details Builds the order tree incrementally, with block-sized stack
   * scratch and at most one unpublished chunk; there is no full-input temporary
   * or pointer directory. Each chunk reserves before constructing its elements
   * and is never resized after publishing addresses. Input need not be sized,
   * common, contiguous, default-constructible, or move-assignable. On
   * allocation, input, or element-construction failure all acquired ownership
   * is reclaimed, but already visited input may have been consumed. No
   * caller-owned references survive construction. Exceptions from input
   * iteration or element construction propagate without leaking ownership.
   * @tparam Range Input range whose iterator-move results can construct T.
   * @param range Source consumed in iteration order through ranges::iter_move.
   * @return Exclusively owned sequence preserving the input element order.
   * @throws std::length_error If the element count exceeds SIZE_MAX or the
   * configured chunk capacity exceeds the backing vector's maximum size.
   * @throws std::bad_alloc If chunk, vector, or tree allocation fails.
   */
 private:
  template <std::ranges::input_range Range>
    requires std::
        constructible_from<T, std::ranges::range_rvalue_reference_t<Range>>
      static PermutableSequence from_range_impl(Range&& range) {
    PermutableSequence result;
    auto it = std::ranges::begin(range);
    const auto end = std::ranges::end(range);
    Chunk* current = nullptr;
    std::size_t offset = 0;
    auto has_input = [&] {
      if constexpr (packed) {
        return it != end;
      } else {
        return (current && offset < current->values.size()) || it != end;
      }
    };
    auto next_block = [&] {
      std::array<typename Block::value_type, Block::capacity> scratch;
      std::size_t count = 0;
      while (count < Block::capacity && has_input()) {
        if constexpr (packed) {
          scratch[count++] = T(std::ranges::iter_move(it));
          ++it;
        } else {
          if (!current || offset == current->values.size()) {
            payload_allocation();
            auto chunk = std::make_unique<Chunk>();
            payload_allocation();
            chunk->values.reserve(chunk_capacity);
            while (chunk->values.size() < chunk_capacity && it != end) {
              if constexpr (std::same_as<T, bool>) {
                chunk->values.emplace_back(bool(std::ranges::iter_move(it)));
              } else {
                chunk->values.emplace_back(std::ranges::iter_move(it));
              }
              ++it;
            }
            current = chunk.get();
            result.chunks_.append(chunk.release());
            offset = 0;
          }
          if constexpr (std::same_as<T, bool>) {
            scratch[count++] = std::addressof(current->values[offset++].value);
          } else {
            scratch[count++] = std::addressof(current->values[offset++]);
          }
        }
      }
      return Block(
          std::span<const typename Block::value_type>(scratch.data(), count));
    };
    // A genuine single-pass cursor: dereferencing never advances the input.
    // A transform_view would incorrectly promise an equality-preserving map.
    struct Blocks {
      decltype(next_block)& next;
      Block current;
      struct Iterator {
        using value_type [[maybe_unused]] = Block;
        using difference_type [[maybe_unused]] = std::ptrdiff_t;
        using iterator_concept [[maybe_unused]] = std::input_iterator_tag;
        Blocks* source;
        Block& operator*() const noexcept { return source->current; }
        Iterator& operator++() {
          source->current = source->next();
          return *this;
        }
        void operator++(int) { ++*this; }
        bool operator==(std::default_sentinel_t) const noexcept {
          return source->current.size() == 0;
        }
      };
      Iterator begin() {
        current = next();
        return {this};
      }
      std::default_sentinel_t end() const noexcept { return {}; }
    } blocks{next_block, {}};
    result.tree_ = Tree::from_blocks(blocks);
    return result;
  }
  /**
   * @brief Return the logical element count.
   * @details Constant time; counts elements, not bytes or bits.
   * @return Count in [0,SIZE_MAX].
   */
  std::size_t size_impl() const noexcept { return tree_.size(); }
  /**
   * @brief Read immutable element i, using zero-based indexing.
   * @details Indirect references remain stable through rotation, merge and
   * ownership moves, until their eventual owner destroys or replaces them.
   * Access uses O(Fanout*height) tree work; no mutable reference is exposed.
   * @param i Zero-based element position in [0,size()).
   * @return Packed T by value, or a const T& to the owned indirect element.
   * @throws std::out_of_range If i >= size().
   */
  const_reference value_at_impl(std::size_t i) const {
    if constexpr (packed) {
      return tree_[i];
    } else {
      return *tree_[i];
    }
  }
  /**
   * @brief Insert one value at a zero-based position.
   * @details Packed storage calls tree insert directly. Indirect storage
   * constructs the value in a payload chunk (allocating a new chunk when the
   * tail is full), then inserts the stable pointer into the tree. On tree
   * failure the chunk entry is popped and the exception propagates. For
   * move-only T the value is consumed even on tree failure.
   * @param position Insertion position in [0,size()].
   * @param value Element to own and insert.
   * @throws std::out_of_range If position > size().
   * @throws std::bad_alloc If chunk, vector, or tree allocation fails.
   * @throws Any Exception from element construction.
   */
  void insert_at_impl(std::size_t position, T value) {
    if (position > tree_.size()) {
      throw std::out_of_range("PermutableSequence: insert position");
    }
    if constexpr (packed) {
      tree_.insert_at(position, std::move(value));
    } else {
      std::unique_ptr<Chunk> staged;
      if (!chunks_.tail || chunks_.tail->values.size() >= chunk_capacity) {
        payload_allocation();
        staged = std::make_unique<Chunk>();
        payload_allocation();
        staged->values.reserve(chunk_capacity);
      }
      auto& slots = staged ? staged->values : chunks_.tail->values;
      if constexpr (std::same_as<T, bool>) {
        slots.emplace_back(bool(std::move(value)));
      } else {
        slots.emplace_back(std::move(value));
      }
      const T* ptr;
      if constexpr (std::same_as<T, bool>) {
        ptr = std::addressof(slots.back().value);
      } else {
        ptr = std::addressof(slots.back());
      }
      try {
        tree_.insert_at(position, ptr);
      } catch (...) {
        slots.pop_back();
        throw;
      }
      if (staged) {
        chunks_.append(staged.release());
      }
    }
  }
  /**
   * @brief Rotate half-open [left,right) left, preserving all outside values.
   * @details Validates the range even for distance zero. Reduces distance
   * modulo nonempty range length; an empty valid range is a no-op. Allocation
   * failure leaves values and indirect addresses unchanged.
   * @param left Inclusive start of the half-open range.
   * @param right Exclusive end of the half-open range.
   * @param distance Left-rotation distance in elements, reduced when nonempty.
   * @throws std::out_of_range If left > right or right > size().
   * @throws std::bad_alloc If tree preflight allocation fails.
   */
  void rotate_left_impl(std::size_t left,
                        std::size_t right,
                        std::size_t distance) {
    tree_.rotate_left(left, right, distance);
  }
  /**
   * @brief Append unchanged donor values and consume its exclusive ownership.
   * @details Only identical configurations can merge. Self-merge and empty
   * donor are no-ops. On success donor is canonical empty. Allocation or size
   * failure leaves both participants unchanged. Only after tree commit are
   * payload chains spliced in
   * O(1), with no chunk traversal, payload move, chunk/vector allocation or
   * reallocation. References from either participant retain their addresses.
   * @param donor Sequence whose unchanged values follow this sequence's values.
   * @throws std::length_error If the combined element count exceeds SIZE_MAX.
   * @throws std::bad_alloc If tree preflight allocation fails.
   */
  void merge_impl(PermutableSequence& donor) {
    if (this == &donor) {
      return;
    }
    tree_.merge(donor.tree_);
    if constexpr (!packed) {
      chunks_.splice(donor.chunks_);
    }
  }

 public:
  /**
   * @brief Find the first position where pred(value) is true.
   * @details Uses the tree's guided descent with boundary-value pruning. The
   * predicate receives the logical element (T by value when packed, const T&
   * when indirect). Caller must guarantee monotonicity: all false results
   * precede all true results.
   * @param pred Predicate on T returning true at or past the search target.
   * @return Position in [0, size()]; size() means pred is false everywhere.
   */
  template <class Pred>
  std::size_t lower_bound(Pred pred) const {
    if constexpr (packed) {
      return tree_.lower_bound(pred);
    } else {
      return tree_.lower_bound([&](const T* p) { return pred(*p); });
    }
  }
  /**
   * @brief Explicit live requested-memory breakdown, excluding allocator/RSS.
   * @details Total is facade_bytes + tree_block_bytes + tree_node_bytes +
   * chunk_header_bytes + vector_capacity_bytes. Other byte/bit fields describe
   * subsets and must not be added again. Facade bytes include the embedded tree
   * and chain handles; chunk headers include each vector object exactly once.
   * Vector capacity includes live slots and slack, with native T alignment.
   * Allocations owned internally by T (for example string buffers) are
   * excluded. Arithmetic saturates at SIZE_MAX; temporary
   * construction/preflight peaks are excluded. This is not an implicit
   * mutation-time statistic.
   */
  struct MemoryUsage {
    /** @brief Inline owner bytes. @details Includes tree and chain handles. */
    std::size_t facade_bytes = sizeof(PermutableSequence);
    /** @brief Leaf count. @details Counts live order-tree leaves. */
    std::size_t blocks = 0;
    /** @brief Node count. @details Counts internal order-tree allocations. */
    std::size_t nodes = 0;
    /** @brief Leaf bytes. @details Includes capacity, metadata, and padding. */
    std::size_t tree_block_bytes = 0;
    /** @brief Node bytes. @details Includes internal allocation padding. */
    std::size_t tree_node_bytes = 0;
    /** @brief Indirect ordering bytes. @details Subset of tree bytes; zero
     * packed. */
    std::size_t order_bytes = 0;
    /** @brief Packed payload capacity in bits. @details Subset; zero indirect.
     */
    std::size_t packed_capacity_bits = 0;
    /** @brief Payload chunk count. @details Zero for packed storage. */
    std::size_t chunks = 0;
    /** @brief Chunk header bytes. @details Includes each vector object once. */
    std::size_t chunk_header_bytes = 0;
    /** @brief Payload allocation bytes. @details Includes live slots and slack.
     */
    std::size_t vector_capacity_bytes = 0;
    /** @brief Live payload slot bytes. @details Subset of vector capacity. */
    std::size_t vector_live_bytes = 0;
    /** @brief Unused payload slot bytes. @details Subset of vector capacity. */
    std::size_t vector_slack_bytes = 0;
    /** @brief Total live bytes. @details Saturated sum of disjoint categories.
     */
    std::size_t total_bytes = sizeof(PermutableSequence);
  };
  /**
   * @brief Enumerate tree allocations and chunks in linear time, without
   * writes.
   * @details Called only explicitly; merge and rotation never call accounting.
   * Performs no allocation and does not throw. Byte arithmetic saturates at
   * SIZE_MAX; see MemoryUsage for included and excluded allocation domains.
   * @return Requested live-memory breakdown, including this facade exactly
   * once.
   */
  MemoryUsage memory_usage() const noexcept {
    MemoryUsage result;
    const auto tree = tree_.memory_usage();
    result.blocks = tree.blocks;
    result.nodes = tree.nodes;
    result.tree_block_bytes = tree.block_bytes;
    result.tree_node_bytes = tree.node_bytes;
    const auto tree_bytes = add(tree.block_bytes, tree.node_bytes);
    if constexpr (packed) {
      result.packed_capacity_bits =
          multiply(multiply(tree.blocks, Block::capacity),
                   std::numeric_limits<T>::digits);
    } else {
      result.order_bytes = tree_bytes;
      for (auto* chunk = chunks_.head; chunk; chunk = chunk->next) {
        ++result.chunks;
        result.chunk_header_bytes =
            add(result.chunk_header_bytes, sizeof(Chunk));
        result.vector_capacity_bytes =
            add(result.vector_capacity_bytes,
                multiply(chunk->values.capacity(), sizeof(Slot)));
        result.vector_live_bytes =
            add(result.vector_live_bytes,
                multiply(chunk->values.size(), sizeof(Slot)));
        result.vector_slack_bytes =
            add(result.vector_slack_bytes,
                multiply(chunk->values.capacity() - chunk->values.size(),
                         sizeof(Slot)));
      }
    }
    result.total_bytes =
        add(sizeof(*this), add(tree_bytes, add(result.chunk_header_bytes,
                                               result.vector_capacity_bytes)));
    return result;
  }
  /**
   * @brief Return live requested bytes, including this object.
   * @details Explicit linear-time accounting, equivalent to
   * memory_usage().total_bytes; does not allocate or throw.
   * @return Requested live bytes, saturated at SIZE_MAX on diagnostic overflow.
   */
 private:
  std::size_t memory_usage_bytes_impl() const noexcept {
    return memory_usage().total_bytes;
  }

#ifdef PIXIE_SEQUENCE_TREE_TESTING
 public:
  /**
   * @brief Test-only order-tree type for existing allocation/structure probes.
   * @details All translation units must agree on PIXIE_SEQUENCE_TREE_TESTING.
   */
  using test_tree_type = Tree;
  /**
   * @brief Return a test-only read-only tree view.
   * @details No mutable facade escape hatch; the reference borrows this
   * object's lifetime and the tree contents may change upon facade mutation.
   * @return Const reference to the embedded order tree.
   */
  const Tree& test_tree() const noexcept { return tree_; }
  /**
   * @brief Fail the next payload allocation after n successful allocation
   * sites.
   * @details -1 disables injection. Sites are chunk new and vector reserve,
   * thread-local per facade instantiation; tree allocation has its own hook.
   * No global allocation interception or production counter is installed.
   * @param n Number of successful allocation sites allowed before bad_alloc;
   * pass -1 to disable injection.
   */
  static void test_fail_payload_after(std::ptrdiff_t n) noexcept {
    payload_failure_ = n;
  }
#endif
};

}  // namespace pixie
