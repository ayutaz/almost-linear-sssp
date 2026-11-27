// Dijkstra with an indexed Fibonacci heap (Fredman–Tarjan): O(m + n log n) amortised time.

#include <cstdint>

#include "baseline/dijkstra_common.hpp"
#include "baseline/heaps/fibonacci_heap.hpp"
#include "sssp/baseline.hpp"

namespace sssp::baseline {

template <class W>
SsspResult<W> dijkstra_fibonacci_heap(const CsrGraph<W>& g, Vertex source) {
  return detail::dijkstra_indexed<W, heaps::FibonacciHeap<W>>(g, source);
}

template SsspResult<std::uint64_t> dijkstra_fibonacci_heap<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, Vertex);
template SsspResult<double> dijkstra_fibonacci_heap<double>(const CsrGraph<double>&, Vertex);

}  // namespace sssp::baseline
