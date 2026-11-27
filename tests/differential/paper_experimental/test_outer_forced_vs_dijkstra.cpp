// The outer construction of Algorithm 8 (arXiv 2611.34911v1, Sec 8.3, Alg 8 l.8-17) forced
// through Mode B parameters (D-7: outer threshold 0 * L^0 and a constant separator b), against
// independent oracles. Acceptance tests of docs/research/04_outer_scales_pipeline.md Sec 6.2:
//   PL-4   every generator family, b in {1, 2, 3, 4, 8, 16, 64}, seeds fixed and varied
//   PL-5   the stage checkers on the intermediate objects
//   PL-6   fallback records: only Alg 4's base cases (and D-10 when a cap is set)
//   PL-7   zero edges, unreachable vertices, loops, parallel edges
//   PL-8   seed independence where NameRoutes samples (D-17 parameters)
//   PL-10  double weights: dyadic exact; non-dyadic within the tolerance or an exception
//   PL-11  unaligned blocks end to end (SC-7 exact); weights above 2^62 (Q-79): every run of
//          SC-8, the 2^62 path, extreme_ratio (PL-7) and the minimal D-20 case must give the
//          outcome predict_range derives from the input (exact distances, the D-20 rejection,
//          or Evaluate's rejection of a state copy beyond the range), never "either"; each
//          throw is also judged against the overflow rule of include/sssp/paper.hpp
//   PL-12  the cases without [slow] are the reduced-size inputs for the sanitizer build
// PL-1 is in tests/differential/paper/test_paper_sssp_differential.cpp, PL-2, PL-3 and PL-9 in
// tests/unit/paper/outer/test_shortest_paths.cpp and tests/unit/paper_experimental/.
// Oracles: Bellman-Ford and baseline binary-heap Dijkstra on the ORIGINAL graph, and
// verify_sssp for the tree (root, feasibility of every arc, tight tree arcs, parents reach s,
// unreachable = infinity). Mode B failure rule: a run either returns a verified result or
// throws; no test accepts a returned result that the oracle rejects.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <initializer_list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "experimental/paper_experimental/sssp.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/outer/base_case.hpp"
#include "paper/outer/contract_restore.hpp"
#include "paper/outer/outer_states.hpp"
#include "paper/outer/scales.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "paper/owner/owner_distances.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"
#include "unit/paper/outer/outer_states_test_support.hpp"
#include "unit/paper/outer/shortest_paths_test_support.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::SsspResult;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::StateId;
namespace pe = sssp::paper_experimental;
namespace po = sssp::paper::outer;
namespace tp = sssp::test::pipeline;
namespace to = sssp::test::outer;

namespace {

constexpr std::uint64_t kSeparators[] = {1, 2, 3, 4, 8, 16, 64};

/// The minimal forcing set (note 04 Sec 3.4); `tune` may change further time-only constants.
template <class Tune>
pe::PaperParams forced(std::uint64_t b, Tune tune) {
  pe::OuterParams o;
  o.threshold = pe::SizeThreshold{0, 0};
  o.separator = pe::SeparatorRule::constant(b);
  tune(o);
  return {pe::LocalParams{}, o};
}

pe::PaperParams forced(std::uint64_t b) {
  return forced(b, [](pe::OuterParams&) {});
}

/// D-17 parameters under which NameRoutes samples already for t >= 6 pivots: |S| = t, loop
/// while |U| > 2t; P is the paper's, so trials are accepted and the names depend on the seed.
pe::PaperParams sampling(std::uint64_t b) {
  return forced(b, [](pe::OuterParams& o) {
    o.sample_factor = 1;
    o.residual_factor = 2;
  });
}

/// The same with P = 1, which lets interowner edges be radial (2 d_e <= w) and so reach the
/// deletion of Alg 8 l.12. With P below the paper's the acceptance bound of Lemma 6.8(2) is
/// lost, so the loop is capped (D-10); certificates stay valid either way (Lemma 6.8(1)).
pe::PaperParams sampling_p1(std::uint64_t b) {
  return forced(b, [](pe::OuterParams& o) {
    o.sample_factor = 1;
    o.residual_factor = 2;
    o.p_log2 = 0;
    o.retry_cap = 8;
  });
}

/// A dense graph with weights in 1..4: many ties, so many nontree edges survive the reduction
/// of Alg 8 l.8 and the pivot-routing instance has more than 2t distinct pairs.
template <class W>
CsrGraph<W> tied_dense_graph(std::uint64_t seed) {
  sssp::gen::GraphSpec spec;
  spec.type = "er";
  spec.n = 120;
  spec.avg_degree = 10.0;
  spec.options.seed = seed;
  spec.options.weights.dist = sssp::gen::WeightDist::kUniformInt;
  spec.options.weights.lo = 1.0;
  spec.options.weights.hi = 4.0;
  spec.options.connected = true;
  const auto base = sssp::to_csr(sssp::gen::generate<std::uint64_t>(spec));
  std::vector<Edge<W>> edges;
  for (const auto& e : base.edges()) edges.push_back({e.u, e.v, static_cast<W>(e.weight)});
  return CsrGraph<W>::from_edges(base.num_vertices(), edges);
}

template <class W>
CsrGraph<W> family(std::string_view type, std::size_t n, std::uint64_t seed,
                   const sssp::gen::WeightSpec& weights) {
  sssp::gen::GraphSpec spec;
  spec.type = std::string(type);
  spec.n = n;
  spec.avg_degree = 3.0;
  spec.degree = 3;
  spec.multi = sssp::gen::MultiEdgePolicy::kKeep;  // parallel edges and loops for regular
  spec.options.seed = seed;
  spec.options.weights = weights;
  return sssp::to_csr(sssp::gen::generate<W>(spec));
}

/// The same family with the integer weights of the std::uint64_t generator in both weight
/// types: for double they are integral values (dyadic, every sum below 2^53 exact; Q-55 plan
/// A). The double generators draw 53-bit reals for most distributions, which is PL-10's input.
template <class W>
CsrGraph<W> integral_family(std::string_view type, std::size_t n, std::uint64_t seed,
                            const sssp::gen::WeightSpec& weights) {
  const auto base = family<std::uint64_t>(type, n, seed, weights);
  std::vector<Edge<W>> edges;
  for (const auto& e : base.edges()) edges.push_back({e.u, e.v, static_cast<W>(e.weight)});
  return CsrGraph<W>::from_edges(base.num_vertices(), edges);
}

sssp::gen::WeightSpec weights_of(sssp::gen::WeightDist dist, double lo, double hi,
                                 double zero_probability = 0.0, int wide_log2 = 40) {
  sssp::gen::WeightSpec w;
  w.dist = dist;
  w.lo = lo;
  w.hi = hi;
  w.zero_probability = zero_probability;
  w.wide_log2 = wide_log2;
  return w;
}

/// The two independent oracles on (g, s).
template <class W>
struct Oracles {
  SsspResult<W> bf;
  SsspResult<W> dij;
};

template <class W>
Oracles<W> oracles(const CsrGraph<W>& g, Vertex s) {
  return {sssp::baseline::bellman_ford(g, s), sssp::baseline::dijkstra_binary_heap(g, s)};
}

/// Exact comparison with both oracles and the tree verification.
template <class W>
void require_exact(const CsrGraph<W>& g, Vertex s, const SsspResult<W>& r, const Oracles<W>& o) {
  REQUIRE(r.dist == o.bf.dist);
  REQUIRE(r.dist == o.dij.dist);
  const auto report = sssp::verify_sssp(g, s, r, &o.bf.dist);
  INFO(report.summary());
  REQUIRE(report.ok);
}

template <class W>
void require_exact(const CsrGraph<W>& g, Vertex s, const SsspResult<W>& r) {
  require_exact(g, s, r, oracles(g, s));
}

/// PL-6: nothing but Alg 4's base case (paper local parameters) is recorded, and no outer base
/// case.
template <class W>
void require_local_records_only(const SsspResult<W>& r) {
  const std::string local = sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams);
  for (const auto& fb : r.diag.fallbacks) REQUIRE(fb.reason == local);
  const auto shares = po::fallback_shares(r.diag);
  REQUIRE(shares.size() <= 1);  // one reason: its share is the representative value (Q-69)
}

/// Outcomes of runs outside the exactness guarantee: verified results, and exceptions counted
/// by type and by the component that raised them (the message up to the first ':').
struct Outcomes {
  std::size_t returned = 0;
  std::map<std::string, std::size_t> thrown;

  [[nodiscard]] std::size_t threw() const {
    std::size_t k = 0;
    for (const auto& [key, count] : thrown) k += count;
    return k;
  }
  [[nodiscard]] std::string summary() const {
    std::ostringstream os;
    os << "returned (verified) " << returned << ", threw " << threw();
    for (const auto& [key, count] : thrown) os << "\n  " << count << " x " << key;
    return os.str();
  }
};

/// Type name and raising component of an exception (std::out_of_range is not accepted: an index
/// outside a container is a defect, not a rejected input).
inline std::string exception_key(const std::exception& e) {
  std::string type = "std::exception";
  if (dynamic_cast<const std::overflow_error*>(&e) != nullptr) {
    type = "overflow_error";
  } else if (dynamic_cast<const std::domain_error*>(&e) != nullptr) {
    type = "domain_error";
  } else if (dynamic_cast<const std::invalid_argument*>(&e) != nullptr) {
    type = "invalid_argument";
  } else if (dynamic_cast<const std::out_of_range*>(&e) != nullptr) {
    type = "out_of_range";
  } else if (dynamic_cast<const std::logic_error*>(&e) != nullptr) {
    type = "logic_error";
  }
  const std::string what = e.what();
  return type + " (" + what.substr(0, what.find(':')) + ")";
}

/// For inputs outside the exactness guarantee (non-dyadic doubles; Q-06, Q-55 plan A, Q-63).
/// Mode B failure rule: the run may throw; a returned result must pass the tree verification
/// and agree with Bellman-Ford within gamma_{n-1} (docs/design/numerics.md Sec 3.1). Nothing in
/// between is accepted.
inline void verified_or_throws(const CsrGraph<double>& g, Vertex s, const pe::PaperParams& params,
                               std::uint64_t seed, Outcomes& out) {
  const auto bf = sssp::baseline::bellman_ford(g, s);
  std::string key;
  try {
    const auto r = pe::paper_sssp(g, s, params, seed);
    const auto report = sssp::verify_sssp(g, s, r, &bf.dist, sssp::Tolerance::for_graph(g));
    INFO(report.summary());
    REQUIRE(report.ok);
    ++out.returned;
    return;
  } catch (const std::exception& e) {
    key = exception_key(e);
  }
  INFO(key);
  REQUIRE(key.rfind("out_of_range", 0) != 0);
  ++out.thrown[key];
}

/// Runs the forced outer branch on `g` from s for every separator and checks it exactly. The
/// oracles are computed once per (g, s): Bellman-Ford alone takes about a minute on
/// adversarial_decrease_key at N = 2000 from s = N - 1 in the debug build.
template <class W>
std::size_t run_all_separators(const CsrGraph<W>& g, Vertex s, std::uint64_t seed) {
  const Oracles<W> o = oracles(g, s);
  std::size_t runs = 0;
  for (const std::uint64_t b : kSeparators) {
    INFO("b=" << b << " s=" << s << " seed=" << seed);
    const auto r = pe::paper_sssp(g, s, forced(b), seed);
    require_exact(g, s, r, o);
    require_local_records_only(r);
    // PL-6: the parameter set is part of the result (PaperParams::record), also with
    // SSSP_DIAGNOSTICS = 0: two deviations, both D-7.
    REQUIRE(r.diag.counters.at("params.deviations") == 2);
    REQUIRE(r.diag.counters.at("params.deviation.D-7.outer.threshold") == 1);
    REQUIRE(r.diag.counters.at("params.outer.separator=" + std::to_string(b)) == 1);
#if SSSP_DIAGNOSTICS
    REQUIRE(r.diag.counters.at("outer.constructions") == 1);
    REQUIRE(r.diag.counters.at("outer.vertices") ==
            r.diag.counters.at("contract_retained_vertices"));
    REQUIRE(r.diag.counters.count("outer_base_cases") == 0);
#endif
    ++runs;
  }
  return runs;
}

/// The stages of one forced run on the contracted source component of (g, s).
template <class W>
struct Staged {
  po::ZeroContraction<W> c;
  po::OuterConstruction<W> r;
  sssp::Diagnostics diag;
};

template <class W>
Staged<W> stages(const CsrGraph<W>& g, Vertex s, const pe::PaperParams& params,
                 std::uint64_t seed) {
  Staged<W> st;
  st.c = po::contract_zero_edges(g, s);
  sssp::paper::model::WordRng rng(seed);
  st.r = po::run_outer_construction(
      st.c.graph, st.c.source, rng,
      pe::outer_construction_params(params, st.c.graph.num_vertices()), &st.diag);
  return st;
}

/// Sec 8.1 on the intermediate objects, against references computed here in exact integers
/// (to::exact: integral or quarter-integer weights only).
///   * J: the coordinate of Def 5.4 is 2^7 J with J the G' distance from the root (O(n^2)
///     Dijkstra on G', independent of Alg 7);
///   * l.1967: a radial edge is deleted iff w > J_x + J_y + d_e; nothing else is deleted;
///   * l.1969-1970: no uncut owner-tree edge is deleted; Lemma 8.1: survivors w <= 2 (J_x + J_y);
///   * Lemma 8.1 (the ordinary graph realizes every geodesic) with the seed of Alg 8 l.14: the
///     least state distance from the native inward copy of s over the states projecting to v
///     is d(s, v) (naive Dijkstra on the state arcs, naive Dijkstra on the graph).
template <class W>
void require_outer_states(const Staged<W>& st) {
  const auto& g = st.c.graph;
  const auto& r = st.r;
  const std::size_t n = g.num_vertices();
  to::OuterCase<W> oc;
  oc.n = n;
  oc.aux = r.aux;
  const std::vector<to::Exact> j = to::reference_j(oc);
  const auto& states = r.os.states;
  for (Vertex v = 0; v < n; ++v) {
    INFO("vertex " << v);
    REQUIRE(to::exact(states.coordinate[v]) == to::kScale * j[v]);
  }
  REQUIRE(r.os.deleted.size() == g.num_edges());
  std::size_t deleted = 0;
  for (EdgeId e = 0; e < g.num_edges(); ++e) {
    const auto& ed = g.edge(e);
    const bool radial = r.cls.edge_class[e] == sssp::paper::owner::OuterEdgeClass::kRadial;
    const to::Exact w = to::exact(ed.weight);
    const bool expect = radial && w > j[ed.u] + j[ed.v] + to::exact(r.cls.detour[e]);
    INFO("edge " << e);
    REQUIRE((r.os.deleted[e] != 0) == expect);
    if (radial && !expect) REQUIRE(w <= 2 * (j[ed.u] + j[ed.v]));
    deleted += expect ? 1 : 0;
  }
  REQUIRE(r.os.num_deleted == deleted);
  for (const EdgeId e : r.od.tree_edges) {
    if (r.od.cut[e] == 0) REQUIRE(r.os.deleted[e] == 0);
  }
  std::vector<to::XArc> arcs;
  arcs.reserve(states.graph.arcs.size());
  for (const auto& a : states.graph.arcs) arcs.push_back({a.tail, a.head, to::exact(a.length)});
  const auto sd = to::state_distances(
      states.num_states(), arcs, static_cast<std::uint32_t>(states.native_inward[st.c.source]));
  std::vector<Edge<W>> edges(g.edges().begin(), g.edges().end());
  const auto phys = to::physical_distances(n, edges, st.c.source);
  std::vector<to::Exact> least(n, to::kExactInf);
  for (StateId x = 0; x < states.num_states(); ++x) {
    const Vertex v = states.projection[x];
    REQUIRE(sd[x] >= phys[v]);  // a state path projects to a physical walk
    least[v] = std::min(least[v], sd[x]);
  }
  for (Vertex v = 0; v < n; ++v) {
    INFO("vertex " << v);
    REQUIRE(least[v] == phys[v]);
  }
}

template <class W>
void require_empty(const std::vector<std::string>& messages, const char* what) {
  INFO(what);
  for (const auto& m : messages) INFO(m);
  REQUIRE(messages.empty());
}

/// PL-5: every stage's invariant checker on the intermediate objects.
template <class W>
void require_stage_invariants(const Staged<W>& st, bool route_names) {
  const auto& g = st.c.graph;
  const auto& r = st.r;
  // Sec 6.1, Lemma 6.8(1): the certificates of NameRoutes against all-pairs distances in G_P.
  if (route_names) {
    const auto pe_edges = sssp::paper::owner::projected_interowner_edges(g, r.od);
    sssp::paper::pivots::RouteNames<W> names;
    names.p_log2 = r.cls.p_log2;
    names.dictionary = r.cls.dictionary;
    for (const EdgeId e : pe_edges.input_edge) names.name.push_back(r.cls.name[e]);
    const auto dist = sssp::test::floyd_warshall<W>(pe_edges.instance.t, pe_edges.instance.edges);
    require_empty<W>(sssp::paper::pivots::check_route_names(pe_edges.instance, names, dist),
                     "check_route_names");
  }
  // Sec 7.1-7.2, Lemmas 7.1-7.2 with Alg 7's J.
  require_empty<W>(sssp::paper::owner::check_aux_metric(g, r.od, r.cls, r.aux, r.jr.distance),
                   "check_aux_metric");
  // Sec 7.3.
  require_empty<W>(sssp::paper::owner::check_calendars(r.aux, r.jr.calendars), "check_calendars");
  require_empty<W>(sssp::paper::owner::check_atoms(r.aux, r.jr.calendars, r.jr.atoms),
                   "check_atoms");
  // Sec 8.1 (P3-OUTER-STATES has no library checker; these are the test-side checks of
  // tests/property/paper/test_outer_states_property.cpp, against exact references).
  require_outer_states(st);
  // Sec 8.2, Lemma 8.4.
  require_empty<W>(po::check_outer_scales(g, r.od, r.cls, r.aux, r.jr.calendars, r.os, r.sc),
                   "check_outer_scales");
  // Def 4.1 / 4.2 on the substituted hierarchy with the seed of Alg 8 l.14.
  REQUIRE(r.sub.hierarchy.c0 <= 2);
  const sssp::paper::SeededOrderedHierarchy<W> seeded{
      r.sub.hierarchy, {{st.c.source, W{0}}}, r.offers};
  const auto report = sssp::paper::validate_seeded_ordered_hierarchy(seeded, g);
  INFO(report.summary());
  REQUIRE(report.ok());
}

}  // namespace

// ------------------------------------------------------- range of std::uint64_t (Q-79, D-20)

using Wide = unsigned __int128;

/// What a forced run on a std::uint64_t input must do, derived here from the input and from the
/// stages before Evaluate (each checked by its own tests), never from the run's result:
///   * kCalendars: an owner distance R of the Sec 5.1 decomposition is >= 2^55, so the top of
///     the Sec 7.3 grid, 2^9 R, does not fit (D-20); the run throws std::overflow_error
///     "calendars:". For 2^54 <= R < 2^55 the outcome also depends on the doublings from the
///     grid origin (measured: one edge of 2^54 runs, a path 2^54, 1 is rejected), so inputs in
///     that band are refused here (REQUIRE) instead of being predicted. Below 2^54 the
///     calendars accept the input.
///   * kEvaluate: some copy of a vertex in the outer state graph lies only on state walks of
///     length >= 2^64 - 1 (exact 128-bit Dijkstra from the native inward copy of s); Evaluate
///     throws "evaluate: state X received only offers beyond the range" for such a copy X.
///   * kExact: neither; the run returns the exact distances.
/// `tentative` is the overflow rule of include/sssp/paper.hpp: some reachable u has
/// d[u] + w >= 2^64 - 1 (exact distances in 128 bits). Only then MAY a run throw.
struct RangeCase {
  enum class Kind : std::uint8_t { kExact, kCalendars, kEvaluate };
  Kind kind = Kind::kExact;
  bool tentative = false;
  Wide max_owner_distance = 0;
  std::vector<to::Exact> dist;       ///< exact distances from s (to::kExactInf: unreachable)
  std::vector<std::uint8_t> beyond;  ///< per state of the substituted hierarchy (kEvaluate)
};

inline constexpr Wide kD20Accepted = Wide{1} << 54;  // R below this: accepted
inline constexpr Wide kD20Rejected = Wide{1} << 55;  // R at or above this: rejected
inline constexpr Wide kSentinel = Wide{WeightTraits<std::uint64_t>::infinity()};

/// `h` is the graph the outer construction sees (connected, positive), `b` the constant
/// separator of `params`, `seed` the NameRoutes seed of the run.
inline RangeCase predict_range(const CsrGraph<std::uint64_t>& h, Vertex s, std::uint64_t b,
                               const pe::PaperParams& params, std::uint64_t seed) {
  using W = std::uint64_t;
  RangeCase rc;
  const std::vector<Edge<W>> edges(h.edges().begin(), h.edges().end());
  rc.dist = to::physical_distances(h.num_vertices(), edges, s);
  for (const auto& e : edges) {
    for (const Vertex u : {e.u, e.v}) {
      if (rc.dist[u] != to::kExactInf && rc.dist[u] + Wide{e.weight} >= kSentinel) {
        rc.tentative = true;
      }
    }
  }
  const auto od = sssp::paper::owner_decomposition(h, static_cast<std::size_t>(b));
  for (const W r : od.owner_distance) {
    rc.max_owner_distance =
        std::max(rc.max_owner_distance, r == WeightTraits<W>::infinity() ? to::kExactInf : Wide{r});
  }
  INFO("max owner distance >= 2^54: " << (rc.max_owner_distance >= kD20Accepted));
  REQUIRE((rc.max_owner_distance < kD20Accepted || rc.max_owner_distance >= kD20Rejected));
  if (rc.max_owner_distance >= kD20Rejected) {
    rc.kind = RangeCase::Kind::kCalendars;
    return rc;
  }
  // Alg 8 l.8-13 as run_outer_construction runs them.
  const auto p = pe::outer_construction_params(params, h.num_vertices());
  sssp::paper::model::WordRng rng(seed);
  const auto cls = sssp::paper::owner::classify_interowner_edges(h, od, rng, p.name_routes);
  const auto aux = sssp::paper::owner::build_aux_graph(h, od, cls);
  sssp::paper::owner::OwnerDistancesOptions odo;
  odo.local = p.local;
  const auto jr = sssp::paper::owner::owner_distances(aux, nullptr, odo);
  const auto os = po::build_outer_states(h, od, cls, aux, jr.scaled_distance);
  const auto factor =
      po::outer_scale_factor_power_separator(cls.p_log2, b, 1, po::kOuterShellRatioLog2);
  const auto sc = po::build_outer_scales(h, od, cls, aux, jr.calendars, os, factor);
  const auto part = po::build_outer_partition(os, sc);
  const auto sub = sssp::paper::substitute_atoms(part.hierarchy, nullptr, p.local);
  std::vector<to::XArc> arcs;
  arcs.reserve(os.states.graph.arcs.size());
  for (const auto& a : os.states.graph.arcs) arcs.push_back({a.tail, a.head, Wide{a.length}});
  const auto sd = to::state_distances(os.states.num_states(), arcs,
                                      static_cast<std::uint32_t>(os.states.native_inward[s]));
  // Lemma 8.1 with the seed of l.14: the least copy of v realizes d(s, v).
  std::vector<Wide> least(h.num_vertices(), to::kExactInf);
  for (StateId x = 0; x < os.states.num_states(); ++x) {
    const Vertex v = os.states.projection[x];
    least[v] = std::min(least[v], sd[x]);
  }
  REQUIRE(least == rc.dist);
  rc.beyond.assign(sub.hierarchy.num_states(), 0);
  for (StateId x = 0; x < sub.hierarchy.num_states(); ++x) {
    const Wide d = sd[sub.old_state[x]];
    if (d != to::kExactInf && d >= kSentinel) {
      rc.beyond[x] = 1;
      rc.kind = RangeCase::Kind::kEvaluate;
    }
  }
  return rc;
}

/// Runs `run` (which returns distances indexed like `oracle`) and requires the outcome that
/// `rc` predicts: the exact distances, or exactly the predicted std::overflow_error. A throw is
/// counted in `violations` when the overflow rule does not allow it (rc.tentative false).
template <class Run>
void require_predicted(const RangeCase& rc, const std::vector<to::Exact>& oracle, Run run,
                       std::size_t& violations) {
  std::string message;
  try {
    const std::vector<std::uint64_t> dist = run();
    REQUIRE(dist.size() == oracle.size());
    for (Vertex v = 0; v < dist.size(); ++v) {
      INFO("vertex " << v);
      const Wide d =
          dist[v] == WeightTraits<std::uint64_t>::infinity() ? to::kExactInf : Wide{dist[v]};
      REQUIRE(d == oracle[v]);
    }
  } catch (const std::overflow_error& e) {
    message = e.what();
  }
  INFO("outcome: " << (message.empty() ? std::string("returned") : message));
  switch (rc.kind) {
    case RangeCase::Kind::kExact:
      REQUIRE(message.empty());
      return;
    case RangeCase::Kind::kCalendars:
      REQUIRE(message.rfind("calendars:", 0) == 0);
      break;
    case RangeCase::Kind::kEvaluate: {
      const std::string prefix = "evaluate: state ";
      REQUIRE(message.rfind(prefix, 0) == 0);
      REQUIRE(message.find("received only offers beyond the range") != std::string::npos);
      const auto x = static_cast<std::size_t>(std::stoul(message.substr(prefix.size())));
      REQUIRE(x < rc.beyond.size());
      REQUIRE(rc.beyond[x] == 1);
      break;
    }
  }
  if (!rc.tentative) ++violations;
}

// ------------------------------------------------------------------------------------ PL-4

TEMPLATE_TEST_CASE("PL-4: forced outer branch equals the oracles on every generator family",
                   "[differential][paper_experimental][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  // Integral weights in both types (dyadic doubles, exact sums): a narrow and a wide range.
  const sssp::gen::WeightSpec weight_cases[] = {
      weights_of(sssp::gen::WeightDist::kUniformInt, 1.0, 40.0),
      weights_of(sssp::gen::WeightDist::kWideRange, 1.0, 1.0, 0.0, 30),
  };
  std::size_t runs = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (const auto& wc : weight_cases) {
      for (const std::uint64_t seed : {std::uint64_t{1}, std::uint64_t{2}}) {
        const auto g = integral_family<W>(info.name, 96, seed, wc);
        const std::size_t n = g.num_vertices();
        for (const Vertex s : {Vertex{0}, static_cast<Vertex>(n / 2)}) {
          INFO(info.name << " " << sssp::gen::to_string(wc.dist) << " graph seed=" << seed);
          runs += run_all_separators(g, s, 1000 + seed);
        }
      }
    }
  }
  REQUIRE(runs == sssp::gen::generator_catalog().size() * 2 * 2 * 2 * std::size(kSeparators));
}

namespace {

/// PL-4 at N = 2000: every family with the graph seed `seed` and its own NameRoutes seed. One
/// test case per (weight type, seed) keeps each well under the ctest timeout (measured in the
/// debug build: about 35 s for seed 7, about 95 s for seed 8, most of it Bellman-Ford).
template <class W>
void forced_at_2000(std::uint64_t seed, Vertex (*source)(std::size_t n)) {
  const auto wc = weights_of(sssp::gen::WeightDist::kUniformInt, 1.0, 1000.0);
  std::size_t runs = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    const auto g = integral_family<W>(info.name, 2000, seed, wc);
    INFO(info.name << " graph seed=" << seed);
    runs += run_all_separators(g, source(g.num_vertices()), 70 + seed);
  }
  REQUIRE(runs == sssp::gen::generator_catalog().size() * std::size(kSeparators));
}

}  // namespace

TEMPLATE_TEST_CASE("PL-4: forced outer branch at N = 2000 on every generator family (seed 7)",
                   "[differential][paper_experimental][outer][pipeline][slow]", std::uint64_t,
                   double) {
  forced_at_2000<TestType>(7, [](std::size_t) { return Vertex{0}; });
}

TEMPLATE_TEST_CASE("PL-4: forced outer branch at N = 2000 on every generator family (seed 8)",
                   "[differential][paper_experimental][outer][pipeline][slow]", std::uint64_t,
                   double) {
  forced_at_2000<TestType>(8, [](std::size_t n) { return static_cast<Vertex>(n - 1); });
}

// ------------------------------------------------------------------------------------ PL-5

TEMPLATE_TEST_CASE("PL-5: the stage checkers pass on the intermediate objects",
                   "[differential][paper_experimental][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  const auto wc = weights_of(sssp::gen::WeightDist::kUniformInt, 1.0, 40.0);
  std::size_t checked = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    const auto g = integral_family<W>(info.name, 60, 3, wc);
    for (const std::uint64_t b : kSeparators) {
      INFO(info.name << " b=" << b);
      const auto st = stages(g, Vertex{0}, forced(b), 5);
      tp::require_oracle(st.c.graph, st.c.source, st.r.dist);
      require_stage_invariants(st, true);
      ++checked;
    }
  }
  // The same where NameRoutes samples (D-17) and where interowner edges are radial and some
  // are deleted at Alg 8 l.12 (P = 1).
  std::uint64_t deleted = 0;
  const auto dense = tied_dense_graph<W>(21);
  for (const std::uint64_t b : {std::uint64_t{2}, std::uint64_t{3}}) {
    for (const bool p1 : {false, true}) {
      INFO("tied dense b=" << b << " p1=" << p1);
      const auto st = stages(dense, Vertex{0}, p1 ? sampling_p1(b) : sampling(b), 5);
      tp::require_oracle(st.c.graph, st.c.source, st.r.dist);
      require_stage_invariants(st, true);
      deleted += st.r.os.num_deleted;
      ++checked;
    }
  }
  CHECK(deleted > 0);
  REQUIRE(checked == sssp::gen::generator_catalog().size() * std::size(kSeparators) + 4);
}

// ------------------------------------------------------------------------------------ PL-6

TEMPLATE_TEST_CASE("PL-6: one fallback record per Alg 4 base case and nothing else",
                   "[differential][paper_experimental][outer][pipeline][fallback]", std::uint64_t,
                   double) {
  using W = TestType;
  const auto wc = weights_of(sssp::gen::WeightDist::kUniformInt, 1.0, 40.0);
  const std::string local = sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams);
  std::uint64_t over_cover = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    const auto g = integral_family<W>(info.name, 80, 4, wc);
    for (const std::uint64_t b : kSeparators) {
      INFO(info.name << " b=" << b);
      const auto st = stages(g, Vertex{0}, forced(b), 9);
      const auto& r = st.r;
      // With the paper's local parameters every atom (fewer than b <= 64 vertices) takes Alg
      // 4's base case: one record per atom of Alg 7 (Sec 7.3) and per atom of Sec 8.2.
      const std::size_t alg7_atoms = r.jr.atoms.num_atoms();
      const std::size_t sec82_atoms = r.sc.atoms.size();
      REQUIRE(st.diag.fallbacks.size() == alg7_atoms + sec82_atoms);
      std::uint64_t vertices = 0;
      for (std::size_t i = 0; i < alg7_atoms; ++i) vertices += r.jr.atoms.members(i).size();
      for (const auto& atom : r.sc.atoms) vertices += atom.states.size();
      const auto shares = po::fallback_shares(st.diag);
      if (st.diag.fallbacks.empty()) {
        REQUIRE(shares.empty());
        continue;
      }
      REQUIRE(shares.size() == 1);
      CHECK(shares[0].reason == local);
      CHECK(shares[0].records == alg7_atoms + sec82_atoms);
      CHECK(shares[0].vertices == vertices);
      // Q-69: the vertex sum of a reason may exceed N (a vertex in an Alg 7 atom and in
      // several Sec 8.2 atoms); it is reported against N, never added across reasons.
      if (shares[0].vertices > st.c.graph.num_vertices()) ++over_cover;
    }
  }
  CHECK(over_cover > 0);
}

TEST_CASE("PL-6: a NameRoutes retry cap (D-10) is recorded under its own reason",
          "[differential][paper_experimental][outer][pipeline][fallback]") {
  using W = std::uint64_t;
  const std::string local = sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams);
  std::size_t capped = 0;
  for (const std::uint64_t graph_seed : {1U, 2U, 3U, 4U, 5U, 6U}) {
    const auto g = tied_dense_graph<W>(graph_seed);
    const auto params = forced(1, [](pe::OuterParams& o) {
      o.sample_factor = 1;
      o.residual_factor = 1;
      o.p_log2 = 0;
      o.retry_cap = 1;
    });
    const auto r = pe::paper_sssp(g, Vertex{0}, params, graph_seed);
    require_exact(g, Vertex{0}, r);  // certificates stay valid under the cap (Lemma 6.8(1))
    CHECK(r.diag.counters.at("params.deviation.D-10.outer.retry_cap") == 1);
    for (const auto& share : po::fallback_shares(r.diag)) {
      CHECK(share.reason != po::kOuterBaseCaseScanReason);
      if (share.reason != local) {
        ++capped;
        CHECK(share.reason.find("name_routes") != std::string::npos);
      }
    }
  }
  // P = 1 never passes the acceptance test of Alg 5 l.13, so the cap is reached every time.
  CHECK(capped == 6);
}

// ------------------------------------------------------------------------------------ PL-7

TEMPLATE_TEST_CASE("PL-7: zero edges, unreachable vertices, loops and parallel edges",
                   "[differential][paper_experimental][outer][pipeline]", std::uint64_t, double) {
  using W = TestType;
  std::size_t runs = 0;
  std::size_t range_runs = 0;
  std::size_t range_violations = 0;
  // The mandatory catalogue (single vertex, no edges, zero cycles, loops, multi-edges,
  // disconnected parts, ...), every source.
  for (const auto& gc : sssp::test::required_graphs<W>()) {
    const auto g = gc.graph();
    const std::size_t n = g.num_vertices();
    // Every source on the small cases; the given source and both ends on the large ones.
    std::vector<Vertex> sources;
    if (n <= 20) {
      for (Vertex s = 0; s < n; ++s) sources.push_back(s);
    } else {
      sources = {gc.source, Vertex{0}, static_cast<Vertex>(n - 1)};
    }
    for (const Vertex s : sources) {
      INFO(gc.name);
      if (gc.name == "extreme_ratio") {
        if constexpr (WeightTraits<W>::is_exact) {
          // Weights of 2^60: each run must give the outcome predict_range derives (exact, or
          // the D-20 rejection of the Sec 7.3 calendars).
          const auto c = po::contract_zero_edges(g, s);
          const std::vector<Edge<W>> edges(g.edges().begin(), g.edges().end());
          const auto oracle = to::physical_distances(n, edges, s);
          for (const std::uint64_t b : kSeparators) {
            INFO("b=" << b << " s=" << s);
            const auto rc = predict_range(c.graph, c.source, b, forced(b), 11);
            require_predicted(
                rc, oracle, [&] { return pe::paper_sssp(g, s, forced(b), 11).dist; },
                range_violations);
            ++range_runs;
          }
        } else {
          // 2^60 + 1 is not a double: sums round, outside the guarantee (Q-55 plan A).
          Outcomes out;
          for (const std::uint64_t b : kSeparators) verified_or_throws(g, s, forced(b), 11, out);
        }
        continue;
      }
      runs += run_all_separators(g, s, 11);
    }
  }
  // Generator families with zero weights mixed in (contraction changes N) and parallel edges.
  std::uint64_t contracted_away = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (const double p0 : {0.3, 0.8}) {
      // The integer draws in both weight types (kZeroMixed draws reals for double, outside the
      // exactness guarantee; Q-55 plan A, PL-10).
      const auto g = integral_family<W>(
          info.name, 90, 5, weights_of(sssp::gen::WeightDist::kZeroMixed, 1.0, 40.0, p0));
      INFO(info.name << " p0=" << p0);
      const Vertex s = static_cast<Vertex>(g.num_vertices() / 3);
      runs += run_all_separators(g, s, 12);
      const auto c = po::contract_zero_edges(g, s);
      contracted_away += g.num_vertices() - c.graph.num_vertices();
    }
  }
  CHECK(runs > 0);
  CHECK(contracted_away > 0);
  if constexpr (WeightTraits<W>::is_exact) {
    CHECK(range_runs == 5 * std::size(kSeparators));
    // No tentative length reaches 2^64 - 1 on extreme_ratio, so every D-20 rejection there is
    // outside the overflow rule of include/sssp/paper.hpp (reported; covered by D-20 only).
    WARN("PL-7 extreme_ratio: " << range_violations << " of " << range_runs
                                << " runs rejected by D-20 where the overflow rule does not "
                                   "allow a throw");
  }
}

// ------------------------------------------------------------------------------------ PL-8

namespace {

/// PL-8 on the tied dense graphs of `graph_seeds`, b in {2, 3}, with and without P = 1, for the
/// NameRoutes seeds 0 .. seeds-1.
void seed_independence(std::initializer_list<std::uint64_t> graph_seeds, std::uint64_t seeds) {
  using W = std::uint64_t;
  std::uint64_t radial_interowner = 0;
  std::uint64_t deleted = 0;
  std::uint64_t varied_configs = 0;
  for (const std::uint64_t graph_seed : graph_seeds) {
    const auto g = tied_dense_graph<W>(graph_seed);
    const auto bf = sssp::baseline::bellman_ford(g, Vertex{0});
    for (const std::uint64_t b : {std::uint64_t{2}, std::uint64_t{3}}) {
      for (const bool p1 : {false, true}) {
        const auto params = p1 ? sampling_p1(b) : sampling(b);
        std::vector<std::vector<W>> dictionaries;
        for (std::uint64_t seed = 0; seed < seeds; ++seed) {
          INFO("graph " << graph_seed << " b=" << b << " p1=" << p1 << " seed=" << seed);
          const auto r = pe::paper_sssp(g, Vertex{0}, params, seed);
          REQUIRE(r.dist == bf.dist);
          REQUIRE(sssp::verify_sssp(g, Vertex{0}, r, &bf.dist).ok);
          // What the run exercised, read from the stages (independent of SSSP_DIAGNOSTICS).
          const auto st = stages(g, Vertex{0}, params, seed);
          const std::size_t t = st.r.od.pivots.size();
          REQUIRE(t >= 6);
          // NameRoutes samples: more than 2t distinct pivot pairs (Alg 5 l.5 under D-17).
          const auto pe_edges = sssp::paper::owner::projected_interowner_edges(st.c.graph, st.r.od);
          REQUIRE(pe_edges.instance.edges.size() > 4 * t);
          const auto edges = st.c.graph.edges();
          for (EdgeId e = 0; e < edges.size(); ++e) {
            if (st.r.cls.interowner[e] != 0 &&
                st.r.cls.edge_class[e] == sssp::paper::owner::OuterEdgeClass::kRadial) {
              ++radial_interowner;
            }
          }
          deleted += st.r.os.num_deleted;
          dictionaries.push_back(st.r.cls.dictionary);
#if SSSP_DIAGNOSTICS
          REQUIRE(st.diag.counters.at("name_routes.trials") > 0);
#endif
        }
        std::ranges::sort(dictionaries);
        const auto distinct = std::ranges::unique(dictionaries).begin() - dictionaries.begin();
        varied_configs += distinct >= 2 ? 1 : 0;
      }
    }
  }
  // The seed changed the route names in some configuration, never the distances.
  CHECK(varied_configs > 0);
  CHECK(radial_interowner > 0);
  CHECK(deleted > 0);
}

}  // namespace

TEST_CASE("PL-8: the distances do not depend on the seed (reduced: one graph, 4 seeds)",
          "[differential][paper_experimental][outer][pipeline]") {
  seed_independence({21}, 4);
}

TEST_CASE("PL-8: the distances do not depend on the seed where NameRoutes samples (16 seeds)",
          "[differential][paper_experimental][outer][pipeline][slow]") {
  seed_independence({21, 22}, 16);
}

// ------------------------------------------------------------------------------------ PL-10

TEST_CASE("PL-10: dyadic double weights are exact",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = double;
  // Quarter-integers k / 4, 1 <= k <= 160: every sum below 2^53 is exact.
  std::size_t runs = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    auto base = family<std::uint64_t>(info.name, 90, 8,
                                      weights_of(sssp::gen::WeightDist::kUniformInt, 1.0, 160.0));
    std::vector<Edge<W>> edges;
    for (const auto& e : base.edges()) {
      edges.push_back({e.u, e.v, static_cast<W>(e.weight) * 0.25});
    }
    const auto g = CsrGraph<W>::from_edges(base.num_vertices(), edges);
    INFO(info.name);
    runs += run_all_separators(g, Vertex{1}, 31);
  }
  REQUIRE(runs == sssp::gen::generator_catalog().size() * std::size(kSeparators));
}

TEST_CASE("PL-10: non-dyadic double weights give a verified result or an exception",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = double;
  // Outside the exactness guarantee (Q-06, Q-55 plan A, Q-63). Mode B failure rule: the run
  // may throw; a returned result must pass the tree verification and agree with Bellman-Ford
  // within gamma_{n-1} (docs/design/numerics.md Sec 3.1). Nothing in between is accepted.
  Outcomes out;
  for (const auto& info : sssp::gen::generator_catalog()) {
    for (const auto dist :
         {sssp::gen::WeightDist::kUniformReal, sssp::gen::WeightDist::kLogUniform}) {
      const auto g = family<W>(info.name, 90, 13, weights_of(dist, 0.1, 97.3));
      for (const std::uint64_t b : kSeparators) {
        INFO(info.name << " " << sssp::gen::to_string(dist) << " b=" << b);
        verified_or_throws(g, Vertex{0}, forced(b), 3, out);
      }
    }
  }
  INFO(out.summary());
  CHECK(out.returned + out.threw() ==
        sssp::gen::generator_catalog().size() * 2 * std::size(kSeparators));
  // Printed so that the share of failures, and who raised them, is on record with the run.
  WARN("PL-10 non-dyadic doubles: " << out.summary());
}

// ------------------------------------------------------------------------------------ PL-11

namespace {

/// P = 1 (D-17), so that F is small enough for an unaligned block on std::uint64_t (note 04
/// Sec 0 item 2).
pe::PaperParams p1_only(std::uint64_t b) {
  return forced(b, [](pe::OuterParams& o) { o.p_log2 = 0; });
}

/// The outer construction on `g` itself (no contraction, so the vertex ids, the MST root and
/// the pivots are those of the hand analysis).
po::OuterConstruction<std::uint64_t> direct(const CsrGraph<std::uint64_t>& g, Vertex s,
                                            const pe::PaperParams& params) {
  sssp::paper::model::WordRng rng(3);
  return po::run_outer_construction(g, s, rng,
                                    pe::outer_construction_params(params, g.num_vertices()));
}

}  // namespace

TEST_CASE("PL-11: unaligned block end to end (SC-7 input, corridor 2^60)",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = std::uint64_t;
  // L3 with the corridor edge 0 - 1 at 2^60, P = 1, b = 2: f = 29 and F^2 s_U = 3 * 2^58 < 2^60,
  // so the blocks {1, 3} of the inward window and of outward window 0 are unaligned.
  const W big = W{1} << 60;
  const auto g = CsrGraph<W>::from_edges(4, {{0, 1, big}, {1, 2, 3}, {2, 3, 3}, {1, 3, 5}});
  for (Vertex s = 0; s < 4; ++s) {
    INFO("s=" << s);
    const auto r = direct(g, s, p1_only(2));
    CHECK(r.factor.log2 == 29);
    CHECK(tp::count_blocks(r.sc, po::OuterBlockKind::kUnaligned) == 2);
    tp::require_oracle(g, s, r.dist);
  }
}

TEST_CASE("PL-11: unaligned block with an atom end to end (SC-8 input, weights above 2^63)",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = std::uint64_t;
  // SC-8 as corrected by P3-OUTER-SCALES (test_outer_scales.cpp): the L4 tree with 0 - 3 at
  // 2^63 + 1, 3 - 4 at 2, 4 - 5 at 1, plus the nontree edge 2 - 5 at 2^63 + 2 (a macro with a
  // nonpivot head in each cell; without it no nonsingleton block is attached). P = 1, b = 3:
  // f = 31 and F^2 * 2 = 2^63 < a_U, so two blocks {4, 5} are unaligned; each is one atom.
  const W corridor = (W{1} << 63) + 1;
  const std::vector<Edge<W>> edges{{0, 1, 4}, {1, 2, 4}, {0, 3, corridor},
                                   {3, 4, 2}, {4, 5, 1}, {2, 5, corridor + 1}};
  const auto g = CsrGraph<W>::from_edges(6, edges);
  // Exact oracle in 128 bits (Floyd-Warshall written here): every true distance fits 64 bits,
  // but tentative sums such as d(0, 5) + w(2 - 5) do not, so the baselines may throw (the
  // overflow rule of include/sssp/paper.hpp) and cannot serve as the oracle.
  const Wide wide_inf = ~Wide{0} / 4;
  std::vector<std::vector<Wide>> d(6, std::vector<Wide>(6, wide_inf));
  for (std::size_t v = 0; v < 6; ++v) d[v][v] = 0;
  for (const auto& e : edges) {
    d[e.u][e.v] = std::min<Wide>(d[e.u][e.v], e.weight);
    d[e.v][e.u] = std::min<Wide>(d[e.v][e.u], e.weight);
  }
  for (std::size_t k = 0; k < 6; ++k) {
    for (std::size_t i = 0; i < 6; ++i) {
      for (std::size_t j = 0; j < 6; ++j) d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
    }
  }
  REQUIRE(d[0][5] == Wide{corridor} + 3);  // by hand: 0 - 3 - 4 - 5
  REQUIRE(d[2][5] == Wide{corridor} + 1);  // by hand: the edge 2 - 5

  // Weights above 2^62 (Q-79), with the Sec 4.2 (P3-CHILDORDER-OVERFLOW), Sec 5.1
  // (P3-OWNERDECOMP-OVERFLOW) and Evaluate (D-23) fixes on main. Every vertex distance fits,
  // but from s = 0 and s = 3 one copy of a vertex lies only on state walks longer than
  // 2^64 - 1 (its vertex is reached through another copy). Evaluate rejects such a copy
  // instead of leaving it at infinity. The overflow rule allows that throw there: from s = 0,
  // d(5) + w(5 - 2) = 2^64 + 6; from s = 3, d(0) + w(0 - 3) = 2^64 + 2. predict_range derives
  // the outcome per source; it must be exact from s = 1, 2, 4, 5.
  std::size_t exact_sources = 0;
  std::size_t rejected_sources = 0;
  std::size_t violations = 0;
  for (Vertex s = 0; s < 6; ++s) {
    INFO("s=" << s);
    const auto rc = predict_range(g, s, 3, p1_only(3), 3);
    for (Vertex v = 0; v < 6; ++v) REQUIRE(rc.dist[v] == d[s][v]);
    REQUIRE(rc.kind != RangeCase::Kind::kCalendars);
    require_predicted(
        rc, rc.dist,
        [&] {
          const auto r = direct(g, s, p1_only(3));
          CHECK(r.factor.log2 == 31);
          CHECK(tp::count_blocks(r.sc, po::OuterBlockKind::kUnaligned) == 2);
          CHECK(r.sc.atoms.size() == 3);
          return r.dist;
        },
        violations);
    (rc.kind == RangeCase::Kind::kExact ? exact_sources : rejected_sources) += 1;
    if (rc.kind == RangeCase::Kind::kEvaluate) CHECK(rc.tentative);
  }
  CHECK(exact_sources == 4);
  CHECK(rejected_sources == 2);
  CHECK(violations == 0);  // both rejections are allowed by the overflow rule
}

TEST_CASE("PL-11 / Q-79: a path of four edges of 2^62 (Sec 5.1 root distances)",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = std::uint64_t;
  // The reproduction of P3-OWNERDECOMP-OVERFLOW (docs/TASKS.md) through the outer
  // construction. From s = 1, 2, 3 every distance fits below 2^64 - 1. By hand:
  // d(s, v) = |s - v| 2^62.
  const W q = W{1} << 62;
  const std::vector<Edge<W>> edges{{0, 1, q}, {1, 2, q}, {2, 3, q}, {3, 4, q}};
  const auto g = CsrGraph<W>::from_edges(5, edges);
  std::size_t exact_runs = 0;
  std::size_t d20_runs = 0;
  std::size_t violations = 0;
  for (const std::uint64_t b : kSeparators) {
    for (const Vertex s : {Vertex{1}, Vertex{2}, Vertex{3}}) {
      INFO("b=" << b << " s=" << s);
      std::vector<to::Exact> by_hand(5);
      for (Vertex v = 0; v < 5; ++v) by_hand[v] = Wide{q} * (s > v ? s - v : v - s);
      // On g itself the MST root is vertex 0, so D_T(0, 4) = 2^64 is formed in Sec 5.1 (an
      // auxiliary quantity, no longer a reason to throw).
      const auto rd = predict_range(g, s, b, forced(b), 3);
      REQUIRE(rd.dist == by_hand);
      require_predicted(rd, by_hand, [&] { return direct(g, s, forced(b)).dist; }, violations);
      // Through Alg 8 l.2 the contracted vertices are numbered from s.
      const auto c = po::contract_zero_edges(g, s);
      const auto re = predict_range(c.graph, c.source, b, forced(b), 5);
      require_predicted(
          re, by_hand, [&] { return pe::paper_sssp(g, s, forced(b), 5).dist; }, violations);
      // b = 1: every vertex is a pivot, R = 0, exact. b >= 2: a cell spans an edge of 2^62,
      // R >= 2^62, rejected by D-20 in both entries.
      const auto want = b == 1 ? RangeCase::Kind::kExact : RangeCase::Kind::kCalendars;
      CHECK(rd.kind == want);
      CHECK(re.kind == want);
      (b == 1 ? exact_runs : d20_runs) += 2;
    }
  }
  CHECK(exact_runs == 6);
  CHECK(d20_runs == 36);
  // From s = 2 no tentative length reaches 2^64 - 1 (the largest d[u] + w is 3 * 2^62), so the
  // D-20 rejections from s = 2 are outside the overflow rule (2 entries x 6 separators).
  CHECK(violations == 12);
  // The Mode A entry takes the recorded base case (N < 2^64) and is exact on the same input.
  const auto a = sssp::paper::paper_sssp(g, Vertex{2});
  REQUIRE(a.dist == std::vector<W>{2 * q, q, 0, q, 2 * q});
}

TEST_CASE("D-20 versus the overflow rule: one edge of 2^55 (minimal)",
          "[differential][paper_experimental][outer][pipeline]") {
  using W = std::uint64_t;
  // The smallest input found on which D-20 rejects a run that the overflow rule of
  // include/sssp/paper.hpp does not allow to throw: two vertices, one edge of 2^55, b = 3 (one
  // pivot, R = 2^55). d(0, 1) = 2^55 and the largest d[u] + w is 2^56. One edge of 2^53 runs
  // to the end. Reported against D-20 (adopted), not fixed here.
  for (const int k : {53, 55}) {
    INFO("k=" << k);
    const W w = W{1} << k;
    const auto g = CsrGraph<W>::from_edges(2, {{0, 1, w}});
    for (const Vertex s : {Vertex{0}, Vertex{1}}) {
      INFO("s=" << s);
      const auto rc = predict_range(g, s, 3, forced(3), 5);
      REQUIRE_FALSE(rc.tentative);
      std::size_t violations = 0;
      const std::vector<to::Exact> by_hand =
          s == 0 ? std::vector<to::Exact>{0, Wide{w}} : std::vector<to::Exact>{Wide{w}, 0};
      require_predicted(
          rc, by_hand, [&] { return pe::paper_sssp(g, s, forced(3), 5).dist; }, violations);
      CHECK(rc.kind == (k == 53 ? RangeCase::Kind::kExact : RangeCase::Kind::kCalendars));
      CHECK(violations == (k == 53 ? 0U : 1U));
    }
  }
}
