#include "paper/outer/shortest_paths.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/model/params.hpp"

namespace sssp::paper::outer {

namespace {

void add_ops(OpCounters& into, const OpCounters& from) noexcept {
  into.vertex_scans += from.vertex_scans;
  into.edge_scans += from.edge_scans;
  into.relaxations += from.relaxations;
  into.comparisons += from.comparisons;
  into.additions += from.additions;
  into.queue_pushes += from.queue_pushes;
  into.queue_pops += from.queue_pops;
  into.decrease_keys += from.decrease_keys;
}

}  // namespace

void absorb_diagnostics(Diagnostics& into, Diagnostics&& from) {
  add_ops(into.ops, from.ops);
  for (const auto& [name, value] : from.counters) into.counters[name] += value;
  for (auto& f : from.fallbacks) into.fallbacks.push_back(std::move(f));
}

OuterSeparator paper_outer_separator(std::uint64_t n) {
  // Sec 7.1 l.1644: the owner decomposition of Sec 5.1 with b = L_N^64.
  return OuterSeparator{model::L(n), model::kOuterExponent};
}

std::vector<FallbackShare> fallback_shares(const Diagnostics& diag) {
  std::vector<FallbackShare> shares;
  for (const auto& fb : diag.fallbacks) {
    auto it = std::ranges::find_if(shares, [&](const auto& s) { return s.reason == fb.reason; });
    if (it == shares.end()) {
      shares.push_back(FallbackShare{fb.reason, 0, 0, 0});
      it = shares.end() - 1;
    }
    ++it->records;
    it->vertices += fb.vertices;
    it->edges += fb.edges;
  }
  return shares;
}

template <class W>
OuterConstruction<W> run_outer_construction(const CsrGraph<W>& g, Vertex source,
                                            model::WordRng& rng,
                                            const OuterConstructionParams& params,
                                            Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  if (source >= n) throw std::out_of_range("outer_construction: source is not a vertex");

  OuterConstruction<W> r;
  r.separator = params.separator ? *params.separator : paper_outer_separator(n);
  if (r.separator.base == 0 || r.separator.exponent == 0) {
    throw std::invalid_argument(
        "outer_construction: separator b = base^exponent needs base >= 1 and exponent >= 1");
  }
  // The word the decomposition works with; a b beyond the word range saturates, and any b >= N
  // marks nothing but the root's residual (one pivot).
  const auto b =
      static_cast<std::size_t>(model::saturating_pow(r.separator.base, r.separator.exponent));

  // Alg 8 l.8: MST (Lemma 3.5), owner cells, deletion of the strictly tree-dominated nontree
  // edges. Parallel edges stay (the lightest-per-pair reduction belongs to the base case, l.5).
  r.od = owner_decomposition(g, b, diag);
  // Alg 8 l.9-10: NameRoutes on the projected interowner edges, then radial / macro.
  r.cls = owner::classify_interowner_edges(g, r.od, rng, params.name_routes, diag);
  // Alg 8 l.10: G' and T'.
  r.aux = owner::build_aux_graph(g, r.od, r.cls, diag);
  // Alg 8 l.11: OwnerDistances (Alg 7) gives J.
  owner::OwnerDistancesOptions od_options;
  od_options.local = params.local;
  r.jr = owner::owner_distances(r.aux, diag, od_options);
  // Alg 8 l.12-13: delete radial e = xy with w_e > J_x + J_y + d_e; the outer states with their
  // joins (Sec 8.1).
  r.os = build_outer_states(g, r.od, r.cls, r.aux, r.jr.scaled_distance, diag);
  // Sec 8.2 l.2080: F from P and the exact h = b + 1 (Q-71), K = 16.
  r.factor = outer_scale_factor_power_separator(r.cls.p_log2, r.separator.base,
                                                r.separator.exponent, kOuterShellRatioLog2);
  // Alg 8 l.13: scales, grids, atoms and the partition hierarchy (Sec 8.2, Lemma 8.4).
  r.sc = build_outer_scales(g, r.od, r.cls, r.aux, r.jr.calendars, r.os, r.factor, diag);
  r.part = build_outer_partition(r.os, r.sc, diag);
  // Alg 8 l.13: capped local hierarchies in place of the atoms (Lemma 5.15).
  r.sub = substitute_atoms(r.part.hierarchy, diag, params.local);

  // Alg 8 l.14: seed the native inward copy of s; if it lies in an atom the offer moves to its
  // native local subcopy (replacement procedure, l.1318).
  const std::vector<StateOffer<W>> seed{{r.os.states.native_inward[source], W{0}}};
  r.offers = redirect_offers<W>(r.sub, seed);

  // Alg 8 l.15: Evaluate (Alg 3). The hierarchy is lent to the call and taken back.
  SeededOrderedHierarchy<W> seeded{std::move(r.sub.hierarchy), {{source, W{0}}}, r.offers};
  EvalResult<W> ev = evaluate(seeded);
  r.sub.hierarchy = std::move(seeded.hierarchy);
  if (diag != nullptr) absorb_diagnostics(*diag, std::move(ev.diag));
  r.state_value = std::move(ev.value);

  // Alg 8 l.15: the minimum emitted value per vertex over every state projecting to it (Q-73).
  const W inf = WeightTraits<W>::infinity();
  r.dist.assign(n, inf);
  const std::size_t new_states = r.sub.hierarchy.num_states();
  for (StateId x = 0; x < new_states; ++x) {
    const Vertex v = r.os.states.projection[r.sub.old_state[x]];
    if (r.state_value[x] < r.dist[v]) r.dist[v] = r.state_value[x];
  }
  for (Vertex v = 0; v < n; ++v) {
    if (r.dist[v] == inf) {
      // The input is connected, so every vertex has a finite distance and a protected copy
      // (Lemma 4.6). Never answered by infinity or by another search.
      throw std::logic_error("outer_construction: a vertex received no finite value");
    }
  }

  if (diag != nullptr) {
    SSSP_COUNT(diag->ops.comparisons += new_states);
    SSSP_COUNT(diag->counters["outer.constructions"] += 1);
    SSSP_COUNT(diag->counters["outer.vertices"] += n);
    SSSP_COUNT(diag->counters["outer.edges"] += g.num_edges());
    SSSP_COUNT(diag->counters["outer.pivots"] += r.od.pivots.size());
    SSSP_COUNT(diag->counters["outer.macro_edges"] += static_cast<std::uint64_t>(
                   std::ranges::count(r.cls.edge_class, owner::OuterEdgeClass::kMacro)));
    SSSP_COUNT(diag->counters["outer.radial_edges_deleted"] += r.os.num_deleted);
    SSSP_COUNT(diag->counters["outer.states"] += r.os.states.num_states());
    SSSP_COUNT(diag->counters["outer.new_states"] += new_states);
    SSSP_COUNT(diag->counters["outer.atoms"] += r.sc.atoms.size());
    SSSP_COUNT(diag->counters["outer.scale_factor_log2"] += r.factor.log2);
  }
  return r;
}

template <class W>
std::vector<W> outer_construction(const CsrGraph<W>& g, Vertex source, model::WordRng& rng,
                                  const OuterConstructionParams& params, Diagnostics* diag) {
  return std::move(run_outer_construction(g, source, rng, params, diag).dist);
}

template OuterConstruction<std::uint64_t> run_outer_construction(const CsrGraph<std::uint64_t>&,
                                                                 Vertex, model::WordRng&,
                                                                 const OuterConstructionParams&,
                                                                 Diagnostics*);
template OuterConstruction<double> run_outer_construction(const CsrGraph<double>&, Vertex,
                                                          model::WordRng&,
                                                          const OuterConstructionParams&,
                                                          Diagnostics*);
template std::vector<std::uint64_t> outer_construction(const CsrGraph<std::uint64_t>&, Vertex,
                                                       model::WordRng&,
                                                       const OuterConstructionParams&,
                                                       Diagnostics*);
template std::vector<double> outer_construction(const CsrGraph<double>&, Vertex, model::WordRng&,
                                                const OuterConstructionParams&, Diagnostics*);

}  // namespace sssp::paper::outer
