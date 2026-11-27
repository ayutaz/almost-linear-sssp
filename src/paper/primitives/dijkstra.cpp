// Sec 3.1 / Lemma 3.2: explicit instantiations for the two supported length types.
// The templates live in the header so that callers may also use other search-identifier widths.
#include "paper/primitives/dijkstra.hpp"

namespace sssp::paper {

template class IndexedBinaryHeap<std::uint64_t>;
template class IndexedBinaryHeap<double>;
template class SearchArena<std::uint64_t>;
template class SearchArena<double>;
template SsspResult<std::uint64_t> dijkstra<std::uint64_t>(const CsrGraph<std::uint64_t>&, Vertex,
                                                           std::uint64_t);
template SsspResult<double> dijkstra<double>(const CsrGraph<double>&, Vertex, double);
template SsspResult<std::uint64_t> dijkstra_scan<std::uint64_t>(const CsrGraph<std::uint64_t>&,
                                                                Vertex, std::uint64_t);
template SsspResult<double> dijkstra_scan<double>(const CsrGraph<double>&, Vertex, double);

}  // namespace sssp::paper
