// Dijkstra with an indexed 4-ary heap and decrease-key: every vertex enters the heap at most
// once. O((n + m) log_4 n) time, O(n) heap space.

#include <cstdint>

#include "baseline/dijkstra_common.hpp"
#include "baseline/heaps/dary_heap.hpp"
#include "sssp/baseline.hpp"

namespace sssp::baseline {

template <class W>
SsspResult<W> dijkstra_4ary_heap(const CsrGraph<W>& g, Vertex source) {
  return detail::dijkstra_indexed<W, heaps::IndexedDaryHeap<W, 4>>(g, source);
}

template SsspResult<std::uint64_t> dijkstra_4ary_heap<std::uint64_t>(const CsrGraph<std::uint64_t>&,
                                                                     Vertex);
template SsspResult<double> dijkstra_4ary_heap<double>(const CsrGraph<double>&, Vertex);

}  // namespace sssp::baseline
