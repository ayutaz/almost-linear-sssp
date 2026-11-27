#pragma once
// [Supplement, not in the paper] Threshold-SCC reference hierarchy (research notes 02, sec 2.11
// item 1). Test support: it lets Sec 4 (normalisation, child order, Evaluate) be exercised on
// arbitrary directed state graphs, independently of the Sec 5 constructions.
//
// Construction, for the distinct positive arc lengths W_1 < W_2 < ... in increasing order:
// the strongly connected components of the arcs of length <= W_k are unions of the components
// at W_{k-1}; every component that unites two or more previous components becomes an internal
// node of height W_k whose children are those previous components. Hence
//   * heights are exact and strictly increase towards the root (H2);
//   * the light arcs of a node of height W (length < W, between distinct children) are arcs of
//     the threshold-(< W) graph between distinct SCCs, so their quotient is acyclic (H3);
//   * the roots are the SCCs of the whole graph, so the root quotient is acyclic (H4);
//   * links are a spanning tree of the children taken from arcs of length <= W inside the new
//     SCC (it is strongly, hence weakly, connected by such arcs), so c0 = 1 (H5); every link
//     records its backing arc (H6);
//   * on physical_graph_of_arcs(D, pi, n) every arc is a physical edge, so H1 holds.
// Cost O(K (N + M) alpha) for K distinct lengths: fine for tests, not a paper bound.

#include <cstddef>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace sssp::paper {

/// Builds the threshold-SCC hierarchy (c0 = 1) of `graph` with state images `projection`.
/// Throws std::invalid_argument if projection.size() != graph.num_states or if some arc length
/// is not positive and finite (zero-length cycles would violate Def 4.1, research notes 02
/// sec 0 item 4; zero components must be contracted by the caller, Sec 5.1 l.962), and
/// std::out_of_range for an arc endpoint that is not a state. Self-loops are kept in D but
/// never affect the forest.
template <class W>
OrderedHierarchy<W> build_threshold_scc_hierarchy(StateGraph<W> graph,
                                                  std::vector<Vertex> projection);

/// The undirected physical multigraph on `num_physical` vertices with one edge
/// {pi(x), pi(y)} of length l per state arc x -> y of length l (self-loops included). On it,
/// every state arc projects to a physical walk of length <= its length (H1).
template <class W>
CsrGraph<W> physical_graph_of_arcs(const StateGraph<W>& graph,
                                   const std::vector<Vertex>& projection, std::size_t num_physical);

}  // namespace sssp::paper
