#pragma once
// Edge-list utilities on the interchange format sssp::EdgeList<W> (graph preparation for
// datasets, generators and experiments; not part of any SSSP algorithm and not paper code).
//
//   keep_lightest_per_pair  one edge per unordered vertex pair: the lightest, ties by input index
//   relabel_vertices        apply an injective vertex map; unmapped vertices and their edges go
//   connected_components    canonical component labels (numbered by smallest member)
//   extract_component       the induced subgraph of one component, renumbered
//   component_of            the component that contains a given vertex
//   largest_component       the component with the most vertices (ties: smallest label)
//
// Common contract:
// - Every function validates its input like sssp::validate_edges (std::out_of_range for bad
//   endpoints, std::invalid_argument for negative / NaN / infinite weights, std::length_error
//   for too many vertices or edges).
// - Edges that survive keep their relative input order, their orientation (u, v) up to the
//   vertex renumbering, and their weight bit for bit (they are copied, never recomputed).
// - Each result records, for every output edge, the index of the input edge it came from
//   (`original_edge`), and for renumbered vertex sets the input vertex of every output vertex
//   (`original_vertex`).
// - Weights are only compared, never added (so uint64 weights near the maximum are fine).
// - Running time O(n + m); no hashing, no sorting by weight.
//
// Relation to existing code: sssp::io::DuplicateEdgePolicy::kKeepMinWeight (DIMACS reader)
// keeps the FIRST edge of a pair and overwrites its weight with the pair minimum, so the
// surviving edge is not an input edge; keep_lightest_per_pair keeps the lightest input edge
// itself. sssp::paper::lightest_edge_pairs (src/paper/primitives/edge_pairs.hpp) computes the
// same survivors on a span of edges for the Mode A algorithm and always drops loops; it is
// internal to the paper code and returns indices only.
//
// Instantiated for W = std::uint64_t and W = double.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

namespace sssp {

/// What keep_lightest_per_pair does with self-loops {u, u}.
enum class LoopPolicy {
  kKeep,  ///< a loop is its own pair: the lightest loop of each vertex is kept
  kDrop,  ///< every loop is removed (loops never lie on a shortest path)
};

/// An edge list derived from another one by removing (and possibly renumbering) edges.
template <class W>
struct FilteredEdgeList {
  EdgeList<W> graph;
  std::vector<EdgeId> original_edge;  ///< input edge index of each output edge (increasing)
};

/// The subgraph induced by a vertex set, with the vertices renumbered 0..k-1 in increasing
/// order of their input ids.
template <class W>
struct InducedSubgraph {
  EdgeList<W> graph;
  std::vector<Vertex> original_vertex;  ///< input vertex of each output vertex (increasing)
  std::vector<EdgeId> original_edge;    ///< input edge index of each output edge (increasing)
};

/// Connected components of an undirected edge list (self-loops do not join anything).
struct ComponentLabels {
  /// label[v] in [0, count). Canonical: components are numbered in increasing order of their
  /// smallest vertex, so two edge lists with the same components get the same labels.
  std::vector<std::uint32_t> label;
  std::size_t count = 0;          ///< number of components (isolated vertices included)
  std::vector<std::size_t> size;  ///< vertices per component, indexed by label
};

/// For every unordered pair {u, v} that has at least one edge, keeps exactly the least edge
/// under the (weight, input index) order: the lightest edge, ties broken by the smaller index
/// (weights compare with <, so -0.0 and +0.0 tie). num_vertices is unchanged.
template <class W>
FilteredEdgeList<W> keep_lightest_per_pair(const EdgeList<W>& g,
                                           LoopPolicy loops = LoopPolicy::kKeep);

/// Maps vertex v to new_id[v] in a graph of new_num_vertices vertices. new_id[v] == kNoVertex
/// removes v and every edge incident to it. Throws std::invalid_argument if new_id.size() !=
/// g.num_vertices or two vertices share a new id (that would be a contraction, not a
/// relabelling), std::out_of_range if a new id is >= new_num_vertices, std::length_error if
/// new_num_vertices does not fit a Vertex. Output vertices that no input vertex maps to are
/// isolated.
template <class W>
FilteredEdgeList<W> relabel_vertices(const EdgeList<W>& g, std::span<const Vertex> new_id,
                                     std::size_t new_num_vertices);

/// Canonical connected-component labels (see ComponentLabels).
template <class W>
ComponentLabels connected_components(const EdgeList<W>& g);

/// The induced subgraph of the component with the given label. `components` must be
/// connected_components(g): std::invalid_argument if its size does not match g, std::out_of_range
/// if label >= components.count.
template <class W>
InducedSubgraph<W> extract_component(const EdgeList<W>& g, const ComponentLabels& components,
                                     std::uint32_t label);

/// The induced subgraph of the component that contains v (std::out_of_range if v is not a
/// vertex of g).
template <class W>
InducedSubgraph<W> component_of(const EdgeList<W>& g, Vertex v);

/// The induced subgraph of a component with the most vertices; ties go to the smallest label,
/// i.e. the component whose smallest vertex is smallest. The empty graph gives an empty result.
template <class W>
InducedSubgraph<W> largest_component(const EdgeList<W>& g);

}  // namespace sssp
