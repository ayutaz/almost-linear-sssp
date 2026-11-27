// Mode C (practical-optimized) SSSP through one local hierarchy, task P5-PRACTICAL-LOCAL.
// See local_sssp.hpp for the steps, what is and is not in the paper, and the option log.

#include "experimental/practical/local_sssp.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "experimental/paper_experimental/local_assembly.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/model/params.hpp"
#include "paper/outer/contract_restore.hpp"
#include "sssp/baseline.hpp"

namespace sssp::practical {

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

/// Appends the diagnostics of one stage: counts and counters are summed, fallbacks kept in order.
void absorb(Diagnostics& into, Diagnostics&& from) {
  add_ops(into.ops, from.ops);
  for (const auto& [name, value] : from.counters) into.counters[name] += value;
  for (auto& f : from.fallbacks) into.fallbacks.push_back(std::move(f));
}

/// Where the stages write their counters: nowhere when the counters are compiled out, so a
/// timed build does not pay for the bookkeeping behind a non-null pointer.
Diagnostics* stage_diag([[maybe_unused]] Diagnostics& diag) noexcept {
#if SSSP_DIAGNOSTICS
  return &diag;
#else
  return nullptr;
#endif
}

/// Step 3: Alg 4 l.6-10 on the connected positive-weight graph `g` (N >= 2 vertices). Each
/// stage output lives only until the next stage has consumed it.
template <class W>
paper::LocalHierarchyResult<W> build_hierarchy(const CsrGraph<W>& g,
                                               const PracticalOptions& options, Diagnostics* diag) {
  const std::uint64_t n = g.num_vertices();
  // Separator b of Def 5.1 (paper: L_k^4, Sec 5.1 l.967; here L_N or a constant, D-7).
  const auto b =
      static_cast<std::size_t>(options.separator == 0 ? paper::model::L(n) : options.separator);
  // Alg 4 l.6: MST, tree-dominated edges removed, owner cells (Sec 5.1). Rejects a
  // disconnected input and weights that are not positive and finite.
  const auto od = paper::owner_decomposition(g, b, diag);
  // Alg 4 l.7-8: classification, windows and macro resets (Sec 5.2).
  const auto rs = paper::local_radial_states(g, od, diag);
  // Alg 4 l.9-10: bottleneck minima and join records (Sec 5.3).
  const auto joins = paper::compute_joins(od, rs, diag);
  // Alg 4 l.10-12 with the Mode B assembly: D-11 (alignment unused) or D-8 with the grid.
  const auto assembly = options.assembly == Assembly::kExactHeight
                            ? paper_experimental::LocalAssembly::kExactHeight
                            : paper_experimental::LocalAssembly::kGrid;
  return std::move(paper_experimental::assemble_local_hierarchy(
                       od, rs, joins, paper_experimental::LocalAlignment::kLambda, assembly, diag)
                       .result);
}

/// Steps 3-4 on the contracted source component: the exact distance of every vertex from
/// `source`.
template <class W>
std::vector<W> hierarchy_distances(const CsrGraph<W>& g, Vertex source,
                                   const PracticalOptions& options, Diagnostics& diag) {
  Diagnostics* const sd = stage_diag(diag);
  paper::StateId entry = paper::kNoState;
  // Sec 4.1 (Lemma 4.3): the normalised forest carries its own copy of the state graph and the
  // projection, so the hierarchy of step 3 is released at the end of this block.
  paper::NormalizedHierarchy<W> nh = [&] {
    const paper::LocalHierarchyResult<W> lh = build_hierarchy(g, options, sd);
    entry = lh.native_entry[source];
    SSSP_COUNT(diag.counters["practical_local_states"] += lh.hierarchy.num_states());
    SSSP_COUNT(diag.counters["practical_local_arcs"] += lh.hierarchy.graph.arcs.size());
    SSSP_COUNT(diag.counters["practical_local_nodes"] += lh.hierarchy.nodes.size());
    return paper::normalize(lh.hierarchy, sd);
  }();
  // Sec 4.2 (Lemma 4.5): child orders, representatives and gaps.
  paper::order_children(nh, sd);
  // Alg 3 on the normalisation above, one initial offer: the native entry of the source at 0.
  const std::array<paper::StateOffer<W>, 1> offers{{{entry, W{0}}}};
  paper::EvalOptions<W> eval_options;
  eval_options.supplemental_offers = options.supplemental_offers;
  const paper::EvalResult<W> ev = paper::evaluate_normalized<W>(nh, offers, eval_options, sd);

  // Lemma 5.10 with Lemma 4.6: the least emission over the states of a vertex is its distance.
  const std::size_t n = g.num_vertices();
  const W inf = WeightTraits<W>::infinity();
  std::vector<W> dist(n, inf);
  const auto& projection = nh.forest.projection;
  for (std::size_t x = 0; x < projection.size(); ++x) {
    const W value = ev.value[x];
    if (value < dist[projection[x]]) dist[projection[x]] = value;
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (dist[v] == inf) {
      // The component is connected, so this contradicts Lemma 5.10 / 4.6 (or the arithmetic
      // left the exact model). A failure: no other search is substituted.
      throw std::logic_error("practical::local_sssp: Evaluate emitted no state of vertex " +
                             std::to_string(v) + " of the contracted source component");
    }
  }
  return dist;
}

}  // namespace

std::vector<OptionEntry> option_log(const PracticalOptions& o) {
  return {
      {"assembly", o.assembly == Assembly::kExactHeight ? "exact-height(D-11)" : "lambda-grid(D-8)",
       "assembly of the local hierarchy; the paper has the Sec 5.4 grid with c = 10 and N^16"},
      {"separator", o.separator == 0 ? std::string("L_N") : std::to_string(o.separator),
       "separator size b of the owner cells; the paper uses L_k^4 (D-7)"},
      {"direct_search_max_vertices", std::to_string(o.direct_search_max_vertices),
       "Mode C: not in the paper; contracted components up to this size go to binary-heap "
       "Dijkstra and are recorded as a fallback"},
      {"supplemental_offers", o.supplemental_offers ? "on" : "off(D-14)",
       "Alg 2 l.6 supplemental offers; off loses the visit bound of Lemma 4.7"},
      {"thresholds", "none(N>=2)",
       "Mode C: not in the paper; the size tests of Sec 2 (2 L_k^4, L_N^64) are not applied"},
      {"normalization", "once-reused",
       "Sec 4.1-4.2 run once and Alg 3 runs on that result (evaluate_normalized); no copy into "
       "a seeded hierarchy"},
      {"validation", "skipped", "Def 4.1 / 4.2 are not re-validated at run time"},
      {"connectivity_check", "skipped",
       "Mode C: not in the paper; the contracted source component is connected by construction"},
      {"intermediates", "released-per-stage",
       "Mode C: not in the paper; stage outputs are freed once the next stage has consumed them"},
      {"tree", "lemma-8.5-restore",
       "Alg 8 l.6: first-discovery tree of tight original arcs from the contracted labels"},
  };
}

void record_options(const PracticalOptions& options, Diagnostics& diag) {
  for (const auto& e : option_log(options)) {
    diag.counters["practical." + e.key + "=" + e.value] = 1;
  }
}

std::string direct_search_reason(std::uint64_t direct_search_max_vertices) {
  return "practical-local: N <= " + std::to_string(direct_search_max_vertices) +
         " (direct search, binary-heap Dijkstra)";
}

template <class W>
SsspResult<W> local_sssp(const CsrGraph<W>& g, Vertex s, const PracticalOptions& options) {
  Diagnostics diag;
  record_options(options, diag);  // the choices travel with the result, in every build
  SSSP_COUNT(++diag.counters["practical_local_calls"]);

  // Step 1 (Alg 8 l.2). Throws std::out_of_range for an invalid source.
  const paper::outer::ZeroContraction<W> c =
      paper::outer::contract_zero_edges(g, s, stage_diag(diag));
  const std::uint64_t big_n = c.graph.num_vertices();
  const std::uint64_t big_m = c.graph.num_edges();

  // Step 2 (Mode C: not in the paper): size dispatch on N.
  std::vector<W> contracted_dist;
  if (big_n == 1) {
    SSSP_COUNT(++diag.counters["practical_local_trivial"]);
    diag.fallbacks.push_back({std::string(kTrivialComponentReason), big_n, big_m});
    contracted_dist.assign(1, W{0});
  } else if (big_n <= options.direct_search_max_vertices) {
    SSSP_COUNT(++diag.counters["practical_local_direct_searches"]);
    diag.fallbacks.push_back(
        {direct_search_reason(options.direct_search_max_vertices), big_n, big_m});
    SsspResult<W> direct = baseline::dijkstra_binary_heap(c.graph, c.source);
    add_ops(diag.ops, direct.diag.ops);
    contracted_dist = std::move(direct.dist);
  } else {
    // Steps 3-4: the construction handles every vertex; no FallbackRecord.
    SSSP_COUNT(++diag.counters["practical_local_hierarchies"]);
    contracted_dist = hierarchy_distances(c.graph, c.source, options, diag);
  }

  // Step 5 (Alg 8 l.6, Lemma 8.5).
  SsspResult<W> r = paper::outer::restore_predecessor_tree(g, s, c, contracted_dist);
  absorb(diag, std::move(r.diag));
  r.diag = std::move(diag);
  return r;
}

template SsspResult<std::uint64_t> local_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                              const PracticalOptions&);
template SsspResult<double> local_sssp(const CsrGraph<double>&, Vertex, const PracticalOptions&);

}  // namespace sssp::practical
