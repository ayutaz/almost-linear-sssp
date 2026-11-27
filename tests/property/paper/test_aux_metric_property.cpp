// Property tests for Sec 7.1-7.2 (Lemma 7.1, Lemma 7.2) of arXiv 2611.34911v1 on seeded graph
// families, against oracles written in this file: brute-force tree paths (projected lengths),
// the classification recomputed from the names, G' rebuilt from the edge list, an O(n^2) Dijkstra
// on G' and on G, Floyd-Warshall on G_P, and Kruskal on G'.
//
// Separator sizes: the paper's outer b = L_N^64 saturates (one pivot on every feasible input), so
// b in {1, 2, 3, 4, ceil(log2 N), N, N + 1} is passed directly, as Mode B does (Q-03); the paper's
// b is exercised separately. Route names come from four configurations:
//   * exact names: for each projected edge the least instance weight lambda <= a_e with
//     P lambda >= d_{G_P}(x, y) (P = 1 and P = 2), the strongest valid names of Definition 6.1;
//     they make many interowner edges radial;
//   * NameRoutes with sample_factor = residual_factor = 1, P = 1 and retry cap 64 (Mode B,
//     D-10/D-17) so that sampling rounds run on small instances;
//   * NameRoutes with the paper's parameters (Mode A).
// Lemmas 7.1-7.2 hold for every valid assignment, so each source must satisfy every check.
// Weights are integers (std::uint64_t) or quarter-integers below 2^42 (double), so every sum
// formed here is exact in double as well (Q-06).
//
// Counters are compiled out with SSSP_DIAGNOSTICS=0 (ci-release, release): assertions on them run
// only when they are counted. The coverage checks (Coverage below) are counted from the returned
// structures, so they hold in every build.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::OwnerDecomposition;
using sssp::paper::owner::AuxGraph;
using sssp::paper::owner::OuterEdgeClass;
using sssp::paper::owner::OuterEdgeClassification;
using sssp::paper::pivots::PivotRoutingInstance;
using sssp::paper::pivots::RouteNames;
namespace lt = sssp::test::local;
namespace ow = sssp::paper::owner;

namespace {

enum class NameSource : std::uint8_t { kExactP1, kExactP2, kNameRoutesSmall, kNameRoutesPaper };

const char* name(NameSource s) {
  switch (s) {
    case NameSource::kExactP1:
      return "exact-P1";
    case NameSource::kExactP2:
      return "exact-P2";
    case NameSource::kNameRoutesSmall:
      return "name-routes-small";
    case NameSource::kNameRoutesPaper:
      return "name-routes-paper";
  }
  return "?";
}

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

/// lambda * 2^p, or infinity() if the true value does not fit below the sentinel.
template <class W>
W scaled(W lambda, unsigned p) {
  if constexpr (std::is_floating_point_v<W>) {
    W x = lambda;
    for (unsigned i = 0; i < p; ++i) x = x + x;
    return x;
  } else {
    const unsigned __int128 x = static_cast<unsigned __int128>(lambda) << p;
    if (x >= static_cast<unsigned __int128>(WeightTraits<W>::infinity())) {
      return WeightTraits<W>::infinity();
    }
    return static_cast<W>(x);
  }
}

/// 2 d <= w, exact (d may be the saturated infinity()).
template <class W>
bool radial_test(W d, W w) {
  if constexpr (std::is_floating_point_v<W>) {
    return d + d <= w;
  } else {
    if (d == WeightTraits<W>::infinity()) return false;
    return 2 * static_cast<unsigned __int128>(d) <= static_cast<unsigned __int128>(w);
  }
}

/// lambda <= (2h + 1) w, exact for integers; for double the small h are exact and the saturated
/// h only enlarge the right side.
template <class W>
bool lemma71_lower(W lambda, W w, std::size_t h) {
  if constexpr (std::is_floating_point_v<W>) {
    return lambda <= (2.0 * static_cast<double>(h) + 1.0) * w;
  } else {
    const long double rhs =
        (2.0L * static_cast<long double>(h) + 1.0L) * static_cast<long double>(w);
    if (rhs >= 0x1p100L) return true;
    const unsigned __int128 r = (2 * static_cast<unsigned __int128>(h) + 1) * w;
    return static_cast<unsigned __int128>(lambda) <= r;
  }
}

/// h r with h possibly saturated; R <= h r.
template <class W>
bool within_hops(W R, W r, std::size_t h) {
  if constexpr (std::is_floating_point_v<W>) {
    return R <= static_cast<double>(h) * r;
  } else {
    return static_cast<long double>(R) <= static_cast<long double>(h) * static_cast<long double>(r);
  }
}

/// Exact names (see file comment): the least instance weight lambda <= a_e with
/// 2^p lambda >= d_{G_P}(x, y).
template <class W>
RouteNames<W> exact_names(const PivotRoutingInstance<W>& inst, unsigned p) {
  const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
  std::vector<W> weights;
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::ranges::sort(weights);
  std::vector<W> chosen;
  for (const auto& e : inst.edges) {
    W best = e.weight;
    for (const W a : weights) {
      if (a > e.weight) break;
      if (scaled(a, p) >= dist[e.u][e.v]) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  RouteNames<W> names;
  names.p_log2 = p;
  names.dictionary = chosen;
  std::ranges::sort(names.dictionary);
  names.dictionary.erase(std::unique(names.dictionary.begin(), names.dictionary.end()),
                         names.dictionary.end());
  for (const W c : chosen) {
    const auto it = std::ranges::lower_bound(names.dictionary, c);
    names.name.push_back(static_cast<std::uint32_t>(it - names.dictionary.begin()));
  }
  return names;
}

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t multi_pivot = 0;
  std::uint64_t radial_interowner = 0;
  std::uint64_t macro = 0;
  std::uint64_t root_loops = 0;
  std::uint64_t j_below_owner_distance = 0;  ///< vertices with J < R
  std::uint64_t saturated = 0;
};

/// Runs every independent check on one instance (REQUIRE stops at the first violation).
template <class W>
void check_instance(const std::vector<Edge<W>>& edges, std::size_t n, std::size_t b,
                    NameSource source, std::uint64_t seed, Coverage& cov) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  const auto od = sssp::paper::owner_decomposition(g, b);
  const std::size_t h = od.hop_bound;
  const std::size_t t = od.pivots.size();
  const auto tree = lt::root_tree(n, edges, od.tree_edges);
  const auto pivot_of = [&](Vertex v) { return od.pivots[od.owner[v]]; };
  std::vector<W> R(n);
  std::vector<W> r(n);
  for (Vertex v = 0; v < n; ++v) {
    const auto path = lt::tree_path(tree, v, pivot_of(v));
    R[v] = lt::path_sum(edges, path);
    r[v] = lt::path_max(edges, path);
  }

  // Sec 7.1 l.1646: the projected interowner instance of the reduced graph.
  PivotRoutingInstance<W> expect_inst{t, {}};
  std::vector<EdgeId> expect_input;
  std::vector<std::uint8_t> reduced(edges.size(), 0);
  for (EdgeId id = 0; id < edges.size(); ++id) {
    const auto& e = edges[id];
    if (e.u == e.v) continue;
    reduced[id] = e.weight <= lt::path_sum(edges, lt::tree_path(tree, e.u, e.v)) ? 1 : 0;
    if (reduced[id] == 0 || od.owner[e.u] == od.owner[e.v]) continue;
    expect_inst.edges.push_back({od.owner[e.u], od.owner[e.v], R[e.u] + e.weight + R[e.v]});
    expect_input.push_back(id);
  }
  REQUIRE(od.reduced == reduced);
  const auto pe = ow::projected_interowner_edges(g, od);
  REQUIRE(pe.instance.t == t);
  REQUIRE(pe.input_edge == expect_input);
  REQUIRE(pe.instance.edges.size() == expect_inst.edges.size());
  for (std::size_t i = 0; i < expect_inst.edges.size(); ++i) {
    REQUIRE(pe.instance.edges[i].u == expect_inst.edges[i].u);
    REQUIRE(pe.instance.edges[i].v == expect_inst.edges[i].v);
    REQUIRE(pe.instance.edges[i].weight == expect_inst.edges[i].weight);
    REQUIRE(pe.instance_edge[expect_input[i]] == i);
  }

  // Names.
  Diagnostics diag;
  OuterEdgeClassification<W> cls;
  RouteNames<W> names;
  if (source == NameSource::kExactP1 || source == NameSource::kExactP2) {
    names = exact_names(expect_inst, source == NameSource::kExactP1 ? 0U : 1U);
    cls = ow::classify_outer_edges(g, od, pe, names, &diag);
  } else {
    sssp::paper::pivots::NameRoutesParams params;
    if (source == NameSource::kNameRoutesSmall) {
      params.sample_factor = 1;
      params.residual_factor = 1;
      params.p_log2 = 0;
      params.retry_cap = 64;
    }
    // The same seed reproduces the names that classify_interowner_edges draws.
    sssp::paper::model::WordRng rng_names(seed);
    names = sssp::paper::pivots::name_routes(pe.instance, rng_names, params);
    sssp::paper::model::WordRng rng(seed);
    cls = ow::classify_interowner_edges(g, od, rng, params, &diag);
  }
  const auto gp_dist = sssp::test::floyd_warshall<W>(t, expect_inst.edges);
  REQUIRE(sssp::paper::pivots::check_route_names(pe.instance, names, gp_dist).empty());
  REQUIRE(cls.dictionary == names.dictionary);
  REQUIRE(cls.p_log2 == names.p_log2);
  REQUIRE(cls.hop_bound == h);

  // Sec 7.1 l.1650-1651: classification recomputed from the names.
  for (EdgeId id = 0; id < edges.size(); ++id) {
    const auto& e = edges[id];
    INFO("edge " << id << " (" << e.u << ", " << e.v << ") w " << e.weight);
    if (reduced[id] == 0) {
      REQUIRE(cls.edge_class[id] == OuterEdgeClass::kOmitted);
      continue;
    }
    if (od.owner[e.u] == od.owner[e.v]) {
      REQUIRE(cls.interowner[id] == 0);
      REQUIRE(cls.edge_class[id] == OuterEdgeClass::kRadial);
      REQUIRE(cls.detour[id] == W{0});
      continue;
    }
    REQUIRE(cls.interowner[id] == 1);
    const W lambda = names.dictionary[names.name[pe.instance_edge[id]]];
    REQUIRE(cls.dictionary[cls.name[id]] == lambda);
    const W d = scaled(lambda, names.p_log2);
    REQUIRE(cls.detour[id] == d);
    if (d == WeightTraits<W>::infinity()) ++cov.saturated;
    const bool radial = radial_test(d, e.weight);
    REQUIRE(cls.edge_class[id] == (radial ? OuterEdgeClass::kRadial : OuterEdgeClass::kMacro));
    if (radial) {
      ++cov.radial_interowner;
    } else {
      ++cov.macro;
      // Lemma 7.1: lambda / (2h + 1) <= w < 2 P lambda.
      REQUIRE(lemma71_lower(lambda, e.weight, h));
      REQUIRE_FALSE(radial_test(d, e.weight));
    }
    // Lemma 7.1: every cut corridor edge is macro.
    if (od.cut[id] != 0) REQUIRE(cls.edge_class[id] == OuterEdgeClass::kMacro);
  }

  // Sec 7.2 l.1663-1667: G' rebuilt from the edge list.
  const AuxGraph<W> aux = ow::build_aux_graph(g, od, cls, &diag);
  std::vector<Vertex> aux_of(n, kNoVertex);
  std::size_t na = 1;
  for (Vertex v = 0; v < n; ++v) aux_of[v] = od.is_pivot(v) ? 0 : static_cast<Vertex>(na++);
  REQUIRE(aux.num_vertices() == na);
  REQUIRE(aux.aux_vertex == aux_of);
  std::vector<Edge<W>> expect_aux;
  std::vector<EdgeId> expect_aux_input;
  std::uint64_t loops = 0;
  for (EdgeId id = 0; id < edges.size(); ++id) {
    if (cls.edge_class[id] != OuterEdgeClass::kRadial) continue;
    const Vertex x = aux_of[edges[id].u];
    const Vertex y = aux_of[edges[id].v];
    if (x == y) {
      REQUIRE(x == 0);  // only two pivots can collapse
      ++loops;
      continue;
    }
    expect_aux.push_back({x, y, edges[id].weight + cls.detour[id]});
    expect_aux_input.push_back(id);
  }
  cov.root_loops += loops;
  REQUIRE(aux.input_edge == expect_aux_input);
  const auto got = aux.graph.edges();
  REQUIRE(got.size() == expect_aux.size());
  for (std::size_t i = 0; i < got.size(); ++i) {
    REQUIRE(got[i].u == expect_aux[i].u);
    REQUIRE(got[i].v == expect_aux[i].v);
    REQUIRE(got[i].weight == expect_aux[i].weight);
    REQUIRE(aux.aux_edge[expect_aux_input[i]] == i);
    // Lemma 7.2: w' <= 3w/2, i.e. 2w' <= 3w (exact here).
    const W w = edges[expect_aux_input[i]].weight;
    REQUIRE(got[i].weight + got[i].weight <= w + w + w);
  }

  // T' = owner spokes: the uncut tree edges, rooted at a, depth <= h.
  std::vector<EdgeId> expect_tree;
  for (const EdgeId id : od.tree_edges) {
    if (od.cut[id] == 0) expect_tree.push_back(aux.aux_edge[id]);
  }
  std::ranges::sort(expect_tree);
  REQUIRE(aux.tree_edges == expect_tree);
  REQUIRE(aux.tree_edges.size() + 1 == na);
  REQUIRE(aux.tree_parent[0] == kNoVertex);
  REQUIRE(aux.tree_order.size() == na);
  REQUIRE(aux.tree_order.front() == ow::kAuxRoot);
  std::vector<std::size_t> pos(na);
  for (std::size_t i = 0; i < na; ++i) pos[aux.tree_order[i]] = i;
  std::vector<std::uint32_t> depth(na, 0);
  for (const Vertex x : aux.tree_order) {
    if (x == 0) continue;
    const Vertex p = aux.tree_parent[x];
    const EdgeId pe_id = aux.tree_parent_edge[x];
    REQUIRE(p != kNoVertex);
    REQUIRE(pos[p] < pos[x]);
    REQUIRE(std::ranges::binary_search(aux.tree_edges, pe_id));
    const auto& te = got[pe_id];
    REQUIRE(((te.u == x && te.v == p) || (te.u == p && te.v == x)));
    depth[x] = depth[p] + 1;
    REQUIRE(depth[x] <= h);
    // Sec 7.3 l.1693: in T' the root bottleneck and root distance are r_v and R_v.
    const W rp = aux.radius[p];
    REQUIRE(aux.radius[x] == std::max(rp, te.weight));
    REQUIRE(aux.owner_distance[x] == aux.owner_distance[p] + te.weight);
    REQUIRE(aux.radius[x] == r[aux.input_vertex[x]]);
    REQUIRE(aux.owner_distance[x] == R[aux.input_vertex[x]]);
  }

  // J: reference Dijkstra on G' against the naive one written here.
  const auto J = ow::owner_coordinates_ref(aux);
  const auto J_naive = lt::dijkstra_naive<W>(na, expect_aux, {{0, W{0}}});
  REQUIRE(J == J_naive);

  // Lemma 7.2: r <= J <= R <= h r; |J_x - J_y| <= w' <= 3w/2 on radial edges.
  for (Vertex x = 1; x < na; ++x) {
    INFO("aux vertex " << x << " (input " << aux.input_vertex[x] << ")");
    REQUIRE(aux.radius[x] <= J[x]);
    REQUIRE(J[x] <= aux.owner_distance[x]);
    REQUIRE(within_hops(aux.owner_distance[x], aux.radius[x], h));
    if (J[x] < aux.owner_distance[x]) ++cov.j_below_owner_distance;
  }
  for (const auto& e : got) {
    REQUIRE(J[e.u] <= J[e.v] + e.weight);
    REQUIRE(J[e.v] <= J[e.u] + e.weight);
  }

  // Lemma 7.2: J_v bounds a physical route from v to its own pivot: d_G(v, p(v)) <= J_v.
  std::vector<std::vector<W>> from_pivot(t);
  for (std::size_t p = 0; p < t; ++p) {
    from_pivot[p] = lt::dijkstra_naive<W>(n, edges, {{od.pivots[p], W{0}}});
  }
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v)) continue;
    REQUIRE(from_pivot[od.owner[v]][v] <= J[aux_of[v]]);
  }

  // Lemma 7.2: T' is an MST of G': same total weight as Kruskal, and the cycle property.
  {
    std::vector<std::size_t> order(got.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(
        order, [&](std::size_t a, std::size_t c) { return got[a].weight < got[c].weight; });
    std::vector<Vertex> comp(na);
    std::iota(comp.begin(), comp.end(), Vertex{0});
    W kruskal{0};
    for (const std::size_t i : order) {
      const Vertex cu = comp[got[i].u];
      const Vertex cv = comp[got[i].v];
      if (cu == cv) continue;
      for (auto& c : comp) {
        if (c == cv) c = cu;
      }
      kruskal += got[i].weight;
    }
    W tree_total{0};
    for (const EdgeId id : aux.tree_edges) tree_total += got[id].weight;
    REQUIRE(tree_total == kruskal);
    for (EdgeId i = 0; i < got.size(); ++i) {
      if (std::ranges::binary_search(aux.tree_edges, i)) continue;
      Vertex x = got[i].u;
      Vertex y = got[i].v;
      W m{0};
      while (x != y) {
        if (depth[x] >= depth[y]) {
          m = std::max(m, got[aux.tree_parent_edge[x]].weight);
          x = aux.tree_parent[x];
        } else {
          m = std::max(m, got[aux.tree_parent_edge[y]].weight);
          y = aux.tree_parent[y];
        }
      }
      REQUIRE(m <= got[i].weight);
    }
  }

  REQUIRE(ow::check_aux_metric(g, od, cls, aux, J).empty());
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.counters["aux_metric.root_loops"] == loops);
#endif
  ++cov.cases;
  if (t > 1) ++cov.multi_pivot;
}

}  // namespace

TEMPLATE_TEST_CASE("Lemmas 7.1 and 7.2 on seeded graph families", "[property]", std::uint64_t,
                   double) {
  using W = TestType;
  Coverage cov;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights weights[] = {lt::Weights::kTies, lt::Weights::kUniform,
                                 lt::Weights::kExponential};
  const NameSource sources[] = {NameSource::kExactP1, NameSource::kExactP2,
                                NameSource::kNameRoutesSmall, NameSource::kNameRoutesPaper};
  const std::size_t sizes[] = {1, 2, 4, 7, 12, 20, 33};
  for (const auto f : families) {
    for (const auto wd : weights) {
      for (const std::size_t n : sizes) {
        for (std::uint64_t seed = 1; seed <= 2; ++seed) {
          const std::uint64_t gseed = seed * 1000003 + n * 31 + static_cast<std::uint64_t>(f) * 7 +
                                      static_cast<std::uint64_t>(wd);
          const auto edges = lt::make_graph<W>(f, wd, n, gseed);
          for (const std::size_t b : separator_sizes(n)) {
            for (const auto src : sources) {
              INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " seed=" << gseed
                               << " b=" << b << " names=" << name(src));
              check_instance(edges, n, b, src, gseed + b, cov);
              if constexpr (std::is_integral_v<W>) {
                // Heavy integer weights (2^14 times the exponential ones, below 2^55; every
                // tree sum stays below 2^61): the paper's P lambda no longer fits 64 bits, so
                // the saturated detours of the integer model are exercised.
                if (wd == lt::Weights::kExponential) {
                  auto heavy = edges;
                  for (auto& e : heavy) e.weight <<= 14;
                  check_instance(heavy, n, b, src, gseed + b, cov);
                }
              }
            }
          }
        }
      }
    }
  }
  INFO("cases " << cov.cases << ", radial interowner " << cov.radial_interowner << ", macro "
                << cov.macro << ", root loops " << cov.root_loops << ", J < R "
                << cov.j_below_owner_distance << ", saturated " << cov.saturated);
  CHECK(cov.multi_pivot > 0);
  CHECK(cov.radial_interowner > 0);
  CHECK(cov.macro > 0);
  CHECK(cov.root_loops > 0);
  CHECK(cov.j_below_owner_distance > 0);
  if constexpr (std::is_integral_v<W>) CHECK(cov.saturated > 0);
}

TEMPLATE_TEST_CASE("paper outer separator: one owner and J equals the distance from the pivot",
                   "[property]", std::uint64_t, double) {
  using W = TestType;
  const lt::Family families[] = {lt::Family::kSparse, lt::Family::kGrid, lt::Family::kComplete,
                                 lt::Family::kRandomTree};
  for (const auto f : families) {
    for (const std::size_t n : {1, 5, 17, 40}) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        INFO(lt::name(f) << " n=" << n << " seed=" << seed);
        const auto edges = lt::make_graph<W>(f, lt::Weights::kUniform, n, seed);
        const auto g = CsrGraph<W>::from_edges(n, edges);
        const auto od = sssp::paper::owner_decomposition(g, ow::outer_separator_size(n));
        REQUIRE(od.pivots.size() == 1);
        sssp::paper::model::WordRng rng(seed);
        Diagnostics diag;
        const auto cls = ow::classify_interowner_edges(g, od, rng, {}, &diag);
        CHECK(cls.interowner == std::vector<std::uint8_t>(edges.size(), 0));
#if SSSP_DIAGNOSTICS
        CHECK(diag.counters["aux_metric.interowner_edges"] == 0);
#endif
        CHECK(diag.fallbacks.empty());
        const auto aux = ow::build_aux_graph(g, od, cls);
        const auto J = ow::owner_coordinates_ref(aux);
        const auto ref = sssp::baseline::dijkstra_binary_heap(g, od.pivots[0]);
        for (Vertex v = 0; v < n; ++v) CHECK(J[aux.aux_vertex[v]] == ref.dist[v]);
        CHECK(ow::check_aux_metric(g, od, cls, aux, J).empty());
      }
    }
  }
}
