#pragma once
// Baseline SSSP algorithms. Each is instantiated for std::uint64_t and double
// (radix heap: unsigned integers only). All throw std::out_of_range for an invalid source.

#include <cstdint>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::baseline {

/// O(nm) independent oracle for small graphs (queue-free, round-based relaxation).
template <class W>
SsspResult<W> bellman_ford(const CsrGraph<W>& g, Vertex source);

/// Lazy-deletion binary heap Dijkstra (std::priority_queue semantics, own implementation).
template <class W>
SsspResult<W> dijkstra_binary_heap(const CsrGraph<W>& g, Vertex source);

/// Indexed 4-ary heap with decrease-key.
template <class W>
SsspResult<W> dijkstra_4ary_heap(const CsrGraph<W>& g, Vertex source);

/// Fibonacci heap with decrease-key (Fredman–Tarjan).
template <class W>
SsspResult<W> dijkstra_fibonacci_heap(const CsrGraph<W>& g, Vertex source);

/// Radix heap (Ahuja–Mehlhorn–Orlin–Tarjan), monotone integer keys only.
SsspResult<std::uint64_t> dijkstra_radix_heap(const CsrGraph<std::uint64_t>& g, Vertex source);

}  // namespace sssp::baseline
