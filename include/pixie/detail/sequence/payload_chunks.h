#pragma once

#include <concepts>
#include <type_traits>
#include <utility>
#include <vector>

namespace pixie::detail::sequence {

// Own frozen payload vectors independently of the order tree. Publishing a
// chunk transfers ownership; its reserved capacity must not grow afterwards.
// Chain splicing preserves addresses, and destruction uses constant stack.
template <class T>
struct PayloadChunks {
  // Avoid vector<bool>'s proxy storage: indirect reads need actual bool
  // objects.
  struct BoolSlot {
    bool value;
  };
  using Slot = std::conditional_t<std::same_as<T, bool>, BoolSlot, T>;
  static_assert(sizeof(Slot) == sizeof(T));
  struct Chunk {
    std::vector<Slot> values;
    Chunk* next = nullptr;
  };

  Chunk* head = nullptr;
  Chunk* tail = nullptr;

  PayloadChunks() noexcept = default;
  PayloadChunks(const PayloadChunks&) = delete;
  PayloadChunks& operator=(const PayloadChunks&) = delete;
  PayloadChunks(PayloadChunks&& other) noexcept
      : head(std::exchange(other.head, nullptr)),
        tail(std::exchange(other.tail, nullptr)) {}
  PayloadChunks& operator=(PayloadChunks&& other) noexcept {
    if (this != &other) {
      clear();
      head = std::exchange(other.head, nullptr);
      tail = std::exchange(other.tail, nullptr);
    }
    return *this;
  }
  ~PayloadChunks() { clear(); }

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
  void splice(PayloadChunks& donor) noexcept {
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

}  // namespace pixie::detail::sequence
