// Dijkstra with an indexed radix heap (Ahuja–Mehlhorn–Orlin–Tarjan). Dijkstra pops keys in
// non-decreasing order and only inserts keys >= the popped key (weights are non-negative), so
// the radix heap's monotonicity precondition always holds. O(m + 64 n) time for 64-bit keys.

#include <cstdint>

#include "baseline/dijkstra_common.hpp"
#include "baseline/heaps/radix_heap.hpp"
#include "sssp/baseline.hpp"

namespace sssp::baseline {

SsspResult<std::uint64_t> dijkstra_radix_heap(const CsrGraph<std::uint64_t>& g, Vertex source) {
  return detail::dijkstra_indexed<std::uint64_t, heaps::RadixHeap>(g, source);
}

}  // namespace sssp::baseline
