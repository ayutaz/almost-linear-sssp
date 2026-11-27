// The outer state graph, Sec 8.1 (Lemmas 8.1-8.3). See outer_states.hpp.

#include "paper/outer/outer_states.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "paper/owner/calendars.hpp"

namespace sssp::paper::outer {

namespace {

__extension__ typedef unsigned __int128 Wide;  // exact sums of 71-bit values (Q-06)

constexpr unsigned kScaleLog2 = owner::kCalendarScaleLog2;  // J' = 2^7 J (Lemma 7.3)
static_assert(kScaleLog2 == 7);

template <class W>
bool valid_positive(W x) {
  if constexpr (std::is_floating_point_v<W>) {
    return std::isfinite(x) && W{0} < x;
  } else {
    return x != WeightTraits<W>::infinity() && W{0} < x;
  }
}

/// x 2^7 by seven doublings (double only); std::overflow_error if not finite.
double scaled_double(double x) {
  for (unsigned i = 0; i < kScaleLog2; ++i) x += x;
  if (!std::isfinite(x)) {
    throw std::overflow_error("outer_radial_deleted: 2^7 w is not finite");
  }
  return x;
}

}  // namespace

template <class W>
bool outer_radial_deleted(W w, W detour, W scaled_jx, W scaled_jy) {
  if constexpr (std::is_floating_point_v<W>) {
    // l.1967 in calendar units: 2^7 w > J'_x + J'_y + 2^7 d_e (D-21).
    const W sw = scaled_double(w);
    const W sd = scaled_double(detour);
    return scaled_jx + scaled_jy + sd < sw;
  } else {
    static_assert(sizeof(W) <= sizeof(std::uint64_t));
    const Wide sw = static_cast<Wide>(w) << kScaleLog2;
    const Wide sd = static_cast<Wide>(detour) << kScaleLog2;
    return static_cast<Wide>(scaled_jx) + static_cast<Wide>(scaled_jy) + sd < sw;
  }
}

template <class W>
OuterStates<W> build_outer_states(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                  const owner::OuterEdgeClassification<W>& cls,
                                  const owner::AuxGraph<W>& aux, const std::vector<W>& scaled_j,
                                  Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  const std::size_t m = g.num_edges();
  const std::size_t na = aux.num_vertices();
  if (od.num_vertices() != n || od.reduced.size() != m) {
    throw std::invalid_argument("build_outer_states: decomposition of a different graph");
  }
  if (cls.edge_class.size() != m || cls.detour.size() != m || cls.name.size() != m) {
    throw std::invalid_argument("build_outer_states: classification of a different graph");
  }
  if (aux.aux_vertex.size() != n || na == 0 || aux.tree_parent.size() != na ||
      aux.tree_parent_edge.size() != na || aux.input_edge.size() != aux.graph.num_edges()) {
    throw std::invalid_argument("build_outer_states: auxiliary graph of a different graph");
  }
  if (scaled_j.size() != na) {
    throw std::invalid_argument("build_outer_states: one scaled J per auxiliary vertex required");
  }
  if (scaled_j[owner::kAuxRoot] != W{0}) {
    throw std::invalid_argument("build_outer_states: J of the root must be 0 (Lemma 7.2)");
  }
  for (Vertex a = 1; a < na; ++a) {
    if (!valid_positive(scaled_j[a])) {
      throw std::invalid_argument("build_outer_states: J of auxiliary vertex " + std::to_string(a) +
                                  " is not finite and positive (Lemma 7.2)");
    }
  }
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v) != (aux.aux_vertex[v] == owner::kAuxRoot)) {
      throw std::invalid_argument("build_outer_states: pivots and the auxiliary root disagree");
    }
  }

  OuterStates<W> out;
  out.hop_bound = od.hop_bound;

  // The coordinate J' per input vertex (0 at every pivot, identified with the root a).
  std::vector<W> coordinate(n);
  for (Vertex v = 0; v < n; ++v) coordinate[v] = scaled_j[aux.aux_vertex[v]];

  // Step 1 (l.1967-1968): delete radial edges with w > J_x + J_y + d_e, in calendar units.
  const auto edges = g.edges();
  EdgeClassification ec{std::vector<EdgeClass>(m, EdgeClass::kOmitted),
                        std::vector<std::uint32_t>(m, kNoCorridor)};
  out.deleted.assign(m, 0);
  for (EdgeId e = 0; e < m; ++e) {
    const auto& ed = edges[e];
    switch (cls.edge_class[e]) {
      case owner::OuterEdgeClass::kOmitted:
        break;
      case owner::OuterEdgeClass::kMacro:
        ec.edge_class[e] = EdgeClass::kMacro;
        ec.macro_key[e] = cls.name[e];  // the route name, replaced in Sec 8.2 (l.2077)
        break;
      case owner::OuterEdgeClass::kRadial:
        if (outer_radial_deleted(ed.weight, cls.detour[e], coordinate[ed.u], coordinate[ed.v])) {
          out.deleted[e] = 1;
          ++out.num_deleted;
        } else {
          ec.edge_class[e] = EdgeClass::kRadial;
        }
        break;
    }
  }

  // Step 2 (l.1973-1983): Def 5.4 with coordinate J' and K = 16.
  out.states = build_radial_states(g, coordinate, ec, kOuterShellRatioLog2, diag);
  auto& st = out.states;
  out.num_ordinary_states = st.num_states();
  out.num_ordinary_arcs = st.graph.arcs.size();
  out.num_ordinary_blocks = st.num_blocks;
  out.arc_kind.assign(out.num_ordinary_arcs, OuterArcKind::kOrdinary);

  // Step 3 (l.1984-1997): the physical owner tree from T' and the special states.
  out.owner_parent.assign(n, kNoVertex);
  out.owner_parent_edge.assign(n, kNoEdge);
  for (Vertex a = 1; a < na; ++a) {
    const Vertex v = aux.input_vertex[a];
    const EdgeId e = aux.input_edge[aux.tree_parent_edge[a]];
    const auto& ed = edges[e];
    if (ed.u != v && ed.v != v) {
      throw std::invalid_argument("build_outer_states: T' edge does not meet its vertex");
    }
    out.owner_parent[v] = ed.u == v ? ed.v : ed.u;
    out.owner_parent_edge[v] = e;
  }
  out.special_inward.assign(n, kNoState);
  out.special_outward.assign(n, kNoState);
  const auto add_state = [&](Vertex v, Phase ph) {
    const auto id = static_cast<StateId>(st.projection.size());
    st.projection.push_back(v);
    st.phase.push_back(ph);
    st.window.push_back(0);
    st.block.push_back(static_cast<std::uint32_t>(st.num_blocks++));  // singleton (l.1998)
    return id;
  };
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v)) continue;
    out.special_inward[v] = add_state(v, Phase::kInward);
    out.special_outward[v] = add_state(v, Phase::kOutward);
  }
  st.graph.num_states = st.projection.size();
  const auto add_arc = [&](StateId t, StateId h, W len, EdgeId e, OuterArcKind kind) {
    st.graph.arcs.push_back({t, h, len});
    st.arc_kind.push_back(RadialArcKind::kOwnerTree);
    st.arc_edge.push_back(e);
    out.arc_kind.push_back(kind);
  };
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v)) continue;
    const Vertex p = out.owner_parent[v];
    const EdgeId e = out.owner_parent_edge[v];
    const W w = edges[e].weight;
    const bool parent_is_pivot = od.is_pivot(p);
    // l.1985: an ordinary inward copy of v may enter v_in at zero.
    for (const StateId c : {st.native_inward[v], st.lower_inward[v]}) {
      if (c != kNoState) add_arc(c, out.special_inward[v], W{0}, kNoEdge, OuterArcKind::kTreeEntry);
    }
    // l.1986-1987: inward tree states move to their parents, ending at the ordinary inward
    // pivot copy.
    add_arc(out.special_inward[v], parent_is_pivot ? st.native_inward[p] : out.special_inward[p], w,
            e, OuterArcKind::kInwardTree);
    // l.1987-1989: from the ordinary outward pivot copy down the tree of outward tree states;
    // entered only from the parent tree state or the owner pivot.
    add_arc(parent_is_pivot ? st.native_outward[p] : out.special_outward[p], out.special_outward[v],
            w, e, OuterArcKind::kOutwardTree);
    // l.1988: at v the outward tree state may enter the native outward copy at zero.
    add_arc(out.special_outward[v], st.native_outward[v], W{0}, kNoEdge, OuterArcKind::kTreeExit);
  }

  // Step 4 (l.1998-2000): Sec 5.3 on the augmented graph.
  out.joins = compute_joins(od, st, diag);

  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["outer_states.calls"]);
    SSSP_COUNT(diag->counters["outer_states.edges_scanned"] += m);
    SSSP_COUNT(diag->counters["outer_states.deleted_radial_edges"] += out.num_deleted);
    SSSP_COUNT(diag->counters["outer_states.special_states"] +=
               st.num_states() - out.num_ordinary_states);
    SSSP_COUNT(diag->counters["outer_states.special_arcs"] +=
               st.graph.arcs.size() - out.num_ordinary_arcs);
    SSSP_COUNT(diag->counters["outer_states.states"] += st.num_states());
    SSSP_COUNT(diag->counters["outer_states.arcs"] += st.graph.arcs.size());
  }
  return out;
}

#define SSSP_INSTANTIATE_OUTER_STATES(W)                                                          \
  template bool outer_radial_deleted<W>(W, W, W, W);                                              \
  template OuterStates<W> build_outer_states<W>(                                                  \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, const owner::OuterEdgeClassification<W>&, \
      const owner::AuxGraph<W>&, const std::vector<W>&, Diagnostics*);
SSSP_INSTANTIATE_OUTER_STATES(std::uint64_t)
SSSP_INSTANTIATE_OUTER_STATES(double)
#undef SSSP_INSTANTIATE_OUTER_STATES

}  // namespace sssp::paper::outer
