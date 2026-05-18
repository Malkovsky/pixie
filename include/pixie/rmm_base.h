#pragma once

#include <cstddef>
#include <limits>

namespace pixie {

template <class Impl>
class RmMBase {
 public:
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

  std::size_t size() const { return impl().size_impl(); }

  std::size_t rank1(std::size_t end_position) const {
    return impl().rank1_impl(end_position);
  }

  std::size_t rank0(std::size_t end_position) const {
    return impl().rank0_impl(end_position);
  }

  std::size_t select1(std::size_t rank) const {
    return impl().select1_impl(rank);
  }

  std::size_t select0(std::size_t rank) const {
    return impl().select0_impl(rank);
  }

  std::size_t rank10(std::size_t end_position) const {
    return impl().rank10_impl(end_position);
  }

  std::size_t select10(std::size_t rank) const {
    return impl().select10_impl(rank);
  }

  int excess(std::size_t end_position) const {
    return impl().excess_impl(end_position);
  }

  std::size_t fwdsearch(std::size_t start_position, int delta) const {
    return impl().fwdsearch_impl(start_position, delta);
  }

  std::size_t bwdsearch(std::size_t start_position, int delta) const {
    return impl().bwdsearch_impl(start_position, delta);
  }

  std::size_t range_min_query_pos(std::size_t range_begin,
                                  std::size_t range_end) const {
    return impl().range_min_query_pos_impl(range_begin, range_end);
  }

  int range_min_query_val(std::size_t range_begin,
                          std::size_t range_end) const {
    return impl().range_min_query_val_impl(range_begin, range_end);
  }

  std::size_t range_max_query_pos(std::size_t range_begin,
                                  std::size_t range_end) const {
    return impl().range_max_query_pos_impl(range_begin, range_end);
  }

  int range_max_query_val(std::size_t range_begin,
                          std::size_t range_end) const {
    return impl().range_max_query_val_impl(range_begin, range_end);
  }

  std::size_t mincount(std::size_t range_begin, std::size_t range_end) const {
    return impl().mincount_impl(range_begin, range_end);
  }

  std::size_t minselect(std::size_t range_begin,
                        std::size_t range_end,
                        std::size_t rank) const {
    return impl().minselect_impl(range_begin, range_end, rank);
  }

  std::size_t close(std::size_t open_position) const {
    return impl().close_impl(open_position);
  }

  std::size_t open(std::size_t close_position) const {
    return impl().open_impl(close_position);
  }

  std::size_t enclose(std::size_t open_position) const {
    return impl().enclose_impl(open_position);
  }

 private:
  const Impl& impl() const { return static_cast<const Impl&>(*this); }
};

}  // namespace pixie
