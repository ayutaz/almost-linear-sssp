// Sec 7.1-7.2 (arXiv 2611.34911v1, txt l.1642-1690). See aux_metric.hpp for the statement, the
// numbering and the numeric conventions.

#include "paper/owner/aux_metric.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/local/saturating.hpp"
#include "paper/model/params.hpp"
#include "paper/primitives/dijkstra.hpp"

namespace sssp::paper::owner {

namespace {

/// b = L_N^64 (Sec 7.1 l.1644).
constexpr unsigned kOuterSeparatorExponent = 64;

template <class W>
void require_matching(const CsrGraph<W>& g, const OwnerDecomposition<W>& od) {
  if (od.num_vertices() != g.num_vertices() || od.reduced.size() != g.num_edges() ||
      od.cut.size() != g.num_edges()) {
    throw std::invalid_argument("aux_metric: owner decomposition does not belong to the graph");
  }
}

template <class W>
void require_matching(const CsrGraph<W>& g, const OuterEdgeClassification<W>& cls) {
  const std::size_t m = g.num_edges();
  if (cls.edge_class.size() != m || cls.interowner.size() != m || cls.name.size() != m ||
      cls.detour.size() != m || cls.dictionary_detour.size() != cls.dictionary.size()) {
    throw std::invalid_argument("aux_metric: classification does not belong to the graph");
  }
}

/// Reduced edge (Alg 8 l.8) that is not a loop.
template <class W>
bool kept(const Edge<W>& e, const OwnerDecomposition<W>& od, EdgeId id) {
  return e.u != e.v && od.reduced[id] != 0;
}

void add_diagnostics(Diagnostics& into, const Diagnostics& from) {
  for (const auto& [k, v] : from.counters) {
    if (k.find(".max_") != std::string::npos) {
      into.counters[k] = std::max(into.counters[k], v);
    } else {
      into.counters[k] += v;
    }
  }
  into.ops.vertex_scans += from.ops.vertex_scans;
  into.ops.edge_scans += from.ops.edge_scans;
  into.ops.relaxations += from.ops.relaxations;
  into.ops.comparisons += from.ops.comparisons;
  into.ops.additions += from.ops.additions;
  into.ops.queue_pushes += from.ops.queue_pushes;
  into.ops.queue_pops += from.ops.queue_pops;
  into.ops.decrease_keys += from.ops.decrease_keys;
  into.fallbacks.insert(into.fallbacks.end(), from.fallbacks.begin(), from.fallbacks.end());
}

/// 2 x, saturating.
template <class W>
W twice(W x) {
  return sat::add(x, x);
}

/// The integer 2h + 1, saturating.
std::uint64_t two_h_plus_one(std::size_t h) {
  constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();
  const auto hh = static_cast<std::uint64_t>(h);
  return hh > (kMax - 1) / 2 ? kMax : 2 * hh + 1;
}

}  // namespace

std::size_t outer_separator_size(std::size_t n) {
  return static_cast<std::size_t>(model::saturating_pow(model::L(n), kOuterSeparatorExponent));
}

template <class W>
ProjectedEdges<W> projected_interowner_edges(const CsrGraph<W>& g,
                                             const OwnerDecomposition<W>& od) {
  require_matching(g, od);
  const auto edges = g.edges();
  ProjectedEdges<W> out;
  out.instance.t = od.pivots.size();
  out.instance_edge.assign(edges.size(), kNoInstanceEdge);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    const auto id = static_cast<EdgeId>(i);
    if (!kept(e, od, id) || od.owner[e.u] == od.owner[e.v]) continue;
    // Sec 7.1 l.1646: projected length A_e = R_x + w + R_y (a physical walk p(x) -> p(y)).
    const W a = WeightTraits<W>::checked_add(
        WeightTraits<W>::checked_add(od.owner_distance[e.u], e.weight), od.owner_distance[e.v]);
    out.instance_edge[i] = static_cast<std::uint32_t>(out.instance.edges.size());
    out.instance.edges.push_back({od.owner[e.u], od.owner[e.v], a});
    out.input_edge.push_back(id);
  }
  return out;
}

template <class W>
OuterEdgeClassification<W> classify_outer_edges(const CsrGraph<W>& g,
                                                const OwnerDecomposition<W>& od,
                                                const ProjectedEdges<W>& projected,
                                                const pivots::RouteNames<W>& names,
                                                Diagnostics* diag) {
  require_matching(g, od);
  const auto edges = g.edges();
  const std::size_t m = edges.size();
  if (projected.instance_edge.size() != m ||
      projected.input_edge.size() != projected.instance.edges.size()) {
    throw std::invalid_argument("classify_outer_edges: projected edges do not match the graph");
  }
  if (names.name.size() != projected.instance.edges.size()) {
    throw std::invalid_argument("classify_outer_edges: one name per projected edge is required");
  }
  for (std::size_t i = 0; i < names.name.size(); ++i) {
    if (names.name[i] >= names.dictionary.size()) {
      throw std::invalid_argument("classify_outer_edges: name outside the dictionary");
    }
    // Definition 6.1: a certificate is an input weight lambda_e <= a_e.
    if (projected.instance.edges[i].weight < names.dictionary[names.name[i]]) {
      throw std::invalid_argument("classify_outer_edges: name exceeds the projected length");
    }
  }

  OuterEdgeClassification<W> cls;
  cls.hop_bound = od.hop_bound;
  cls.p_log2 = names.p_log2;
  cls.dictionary = names.dictionary;
  // l.1648-1649: d = P lambda with a power of two P, computed once per name by doubling.
  std::uint64_t saturated = 0;
  cls.dictionary_detour.reserve(cls.dictionary.size());
  for (const W lambda : cls.dictionary) {
    const W d = sat::mul_pow2(lambda, names.p_log2);
    if (d == WeightTraits<W>::infinity()) ++saturated;
    cls.dictionary_detour.push_back(d);
  }

  cls.edge_class.assign(m, OuterEdgeClass::kOmitted);
  cls.interowner.assign(m, 0);
  cls.name.assign(m, kNoName);
  cls.detour.assign(m, W{0});
  std::uint64_t interowner = 0;
  std::uint64_t radial_inter = 0;
  std::uint64_t macro = 0;
  std::uint64_t intraowner = 0;
  std::uint64_t omitted = 0;
  for (std::size_t i = 0; i < m; ++i) {
    const auto& e = edges[i];
    if (!kept(e, od, static_cast<EdgeId>(i))) {
      ++omitted;
      continue;
    }
    if (od.owner[e.u] == od.owner[e.v]) {
      // l.1651-1652: intraowner edges are radial with d_e = 0.
      cls.edge_class[i] = OuterEdgeClass::kRadial;
      ++intraowner;
      continue;
    }
    const std::uint32_t k = projected.instance_edge[i];
    if (k == kNoInstanceEdge) {
      throw std::invalid_argument("classify_outer_edges: interowner edge without a projection");
    }
    ++interowner;
    cls.interowner[i] = 1;
    cls.name[i] = names.name[k];
    cls.detour[i] = cls.dictionary_detour[cls.name[i]];
    // l.1650-1651: radial iff 2 d_e <= w (a saturated 2 d_e exceeds every finite w).
    if (!(e.weight < twice(cls.detour[i]))) {
      cls.edge_class[i] = OuterEdgeClass::kRadial;
      ++radial_inter;
    } else {
      cls.edge_class[i] = OuterEdgeClass::kMacro;
      ++macro;
    }
  }
  if (diag != nullptr) {
    SSSP_COUNT(diag->counters["aux_metric.interowner_edges"] += interowner);
    SSSP_COUNT(diag->counters["aux_metric.radial_interowner_edges"] += radial_inter);
    SSSP_COUNT(diag->counters["aux_metric.macro_edges"] += macro);
    SSSP_COUNT(diag->counters["aux_metric.intraowner_edges"] += intraowner);
    SSSP_COUNT(diag->counters["aux_metric.omitted_edges"] += omitted);
    SSSP_COUNT(diag->counters["aux_metric.saturated_detours"] += saturated);
  }
  static_cast<void>(interowner);
  static_cast<void>(radial_inter);
  static_cast<void>(macro);
  static_cast<void>(intraowner);
  static_cast<void>(omitted);
  static_cast<void>(saturated);
  return cls;
}

template <class W>
OuterEdgeClassification<W> classify_interowner_edges(const CsrGraph<W>& g,
                                                     const OwnerDecomposition<W>& od,
                                                     model::WordRng& rng,
                                                     const pivots::NameRoutesParams& params,
                                                     Diagnostics* diag) {
  // Alg 8 l.9: NameRoutes on the projected interowner edges.
  const ProjectedEdges<W> projected = projected_interowner_edges(g, od);
  const pivots::RouteNames<W> names = pivots::name_routes(projected.instance, rng, params);
  if (diag != nullptr) add_diagnostics(*diag, names.diag);
  // Alg 8 l.10 (first half): classify macro and radial edges.
  return classify_outer_edges(g, od, projected, names, diag);
}

template <class W>
AuxGraph<W> build_aux_graph(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                            const OuterEdgeClassification<W>& cls, Diagnostics* diag) {
  require_matching(g, od);
  require_matching(g, cls);
  const auto edges = g.edges();
  const std::size_t n = g.num_vertices();
  const std::size_t m = edges.size();

  AuxGraph<W> aux;
  aux.hop_bound = od.hop_bound;
  // l.1663: identify all pivots into the common vertex a.
  aux.aux_vertex.assign(n, kAuxRoot);
  aux.input_vertex.push_back(kNoVertex);
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v)) continue;
    aux.aux_vertex[v] = static_cast<Vertex>(aux.input_vertex.size());
    aux.input_vertex.push_back(v);
  }
  const std::size_t na = aux.input_vertex.size();

  // l.1663-1666: keep the radial edges with w' = w + d_e; discard the positive loops at a.
  std::vector<Edge<W>> aux_edges;
  aux.aux_edge.assign(m, kNoEdge);
  std::uint64_t root_loops = 0;
  for (std::size_t i = 0; i < m; ++i) {
    if (cls.edge_class[i] != OuterEdgeClass::kRadial) continue;
    const auto& e = edges[i];
    const Vertex x = aux.aux_vertex[e.u];
    const Vertex y = aux.aux_vertex[e.v];
    if (x == y) {
      if (x != kAuxRoot)
        throw std::invalid_argument("build_aux_graph: a loop is classified radial");
      ++root_loops;
      continue;
    }
    aux.aux_edge[i] = static_cast<EdgeId>(aux_edges.size());
    aux_edges.push_back({x, y, WeightTraits<W>::checked_add(e.weight, cls.detour[i])});
    aux.input_edge.push_back(static_cast<EdgeId>(i));
  }

  // l.1664-1665: the owner spokes (the uncut MST edges, all intraowner and radial) form T'.
  std::vector<std::vector<std::pair<Vertex, EdgeId>>> adj(na);
  for (const EdgeId id : od.tree_edges) {
    if (od.cut[id] != 0) continue;
    const EdgeId a = aux.aux_edge[id];
    if (a == kNoEdge) throw std::logic_error("build_aux_graph: an owner spoke is not radial");
    aux.tree_edges.push_back(a);
    adj[aux_edges[a].u].emplace_back(aux_edges[a].v, a);
    adj[aux_edges[a].v].emplace_back(aux_edges[a].u, a);
  }
  std::ranges::sort(aux.tree_edges);
  aux.tree_parent.assign(na, kNoVertex);
  aux.tree_parent_edge.assign(na, kNoEdge);
  aux.radius.assign(na, W{0});
  aux.owner_distance.assign(na, W{0});
  std::vector<std::uint8_t> seen(na, 0);
  seen[kAuxRoot] = 1;
  aux.tree_order.reserve(na);
  aux.tree_order.push_back(kAuxRoot);
  for (std::size_t k = 0; k < aux.tree_order.size(); ++k) {
    const Vertex x = aux.tree_order[k];
    for (const auto& [y, a] : adj[x]) {
      if (seen[y] != 0) continue;
      seen[y] = 1;
      aux.tree_parent[y] = x;
      aux.tree_parent_edge[y] = a;
      aux.tree_order.push_back(y);
    }
  }
  if (aux.tree_order.size() != na || aux.tree_edges.size() + 1 != na) {
    throw std::logic_error("build_aux_graph: the owner spokes do not form a spanning tree");
  }
  // Sec 7.3 l.1693: in T' the root bottleneck and root distance at v are r_v and R_v.
  for (Vertex x = 1; x < na; ++x) {
    aux.radius[x] = od.radius[aux.input_vertex[x]];
    aux.owner_distance[x] = od.owner_distance[aux.input_vertex[x]];
  }

  aux.graph = CsrGraph<W>::from_edges(na, std::move(aux_edges));
  if (diag != nullptr) {
    SSSP_COUNT(diag->counters["aux_metric.aux_vertices"] += na);
    SSSP_COUNT(diag->counters["aux_metric.aux_edges"] += aux.graph.num_edges());
    SSSP_COUNT(diag->counters["aux_metric.root_loops"] += root_loops);
  }
  static_cast<void>(root_loops);
  return aux;
}

template <class W>
std::vector<W> owner_coordinates_ref(const AuxGraph<W>& aux) {
  return dijkstra(aux.graph, kAuxRoot).dist;
}

template <class W>
std::vector<std::string> check_aux_metric(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                          const OuterEdgeClassification<W>& cls,
                                          const AuxGraph<W>& aux, const std::vector<W>& J) {
  std::vector<std::string> out;
  const auto fail = [&out](const std::string& what) { out.push_back(what); };
  if (od.num_vertices() != g.num_vertices() || od.reduced.size() != g.num_edges() ||
      cls.edge_class.size() != g.num_edges() || cls.interowner.size() != g.num_edges() ||
      cls.name.size() != g.num_edges() || cls.detour.size() != g.num_edges() ||
      cls.dictionary_detour.size() != cls.dictionary.size() ||
      aux.aux_edge.size() != g.num_edges() || aux.aux_vertex.size() != g.num_vertices() ||
      J.size() != aux.num_vertices() || aux.graph.num_vertices() != aux.num_vertices()) {
    fail("sizes do not match");
    return out;
  }
  const auto edges = g.edges();
  const std::size_t h = od.hop_bound;
  const std::uint64_t k_lambda = two_h_plus_one(h);
  constexpr W kInf = WeightTraits<W>::infinity();

  // Sec 7.1 and Lemma 7.1.
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    const std::string id = "edge " + std::to_string(i) + ": ";
    const auto c = cls.edge_class[i];
    if (!kept(e, od, static_cast<EdgeId>(i))) {
      if (c != OuterEdgeClass::kOmitted) fail(id + "outside the reduced graph but classified");
      continue;
    }
    if (od.owner[e.u] == od.owner[e.v]) {
      if (c != OuterEdgeClass::kRadial || cls.interowner[i] != 0 || cls.detour[i] != W{0}) {
        fail(id + "intraowner edge is not radial with d = 0");
      }
      continue;
    }
    if (cls.interowner[i] == 0 || cls.name[i] >= cls.dictionary.size()) {
      fail(id + "interowner edge without a name");
      continue;
    }
    const W lambda = cls.dictionary[cls.name[i]];
    const W d = cls.detour[i];
    if (d != cls.dictionary_detour[cls.name[i]] ||
        cls.dictionary_detour[cls.name[i]] != sat::mul_pow2(lambda, cls.p_log2)) {
      fail(id + "d_e is not P lambda_e");
    }
    const bool radial = !(e.weight < twice(d));
    if (c != (radial ? OuterEdgeClass::kRadial : OuterEdgeClass::kMacro)) {
      fail(id + "class disagrees with the test 2 d_e <= w");
    }
    if (od.cut[i] != 0 && c != OuterEdgeClass::kMacro) fail(id + "cut corridor edge is not macro");
    if (c == OuterEdgeClass::kMacro) {
      // Lemma 7.1: lambda / (2h + 1) <= w < 2 P lambda.
      if (sat::mul_int(e.weight, k_lambda) < lambda) fail(id + "macro with lambda > (2h+1) w");
      if (!(e.weight < twice(d))) fail(id + "macro with w >= 2 P lambda");
    }
  }

  // G' and T' (Sec 7.2).
  const std::size_t na = aux.num_vertices();
  const auto aux_edges = aux.graph.edges();
  if (aux.input_edge.size() != aux_edges.size()) fail("G': edge maps do not match");
  for (std::size_t a = 0; a < aux_edges.size() && a < aux.input_edge.size(); ++a) {
    const EdgeId i = aux.input_edge[a];
    const std::string id = "aux edge " + std::to_string(a) + ": ";
    if (i >= edges.size() || cls.edge_class[i] != OuterEdgeClass::kRadial || aux.aux_edge[i] != a) {
      fail(id + "not a radial input edge");
      continue;
    }
    const auto& e = edges[i];
    const auto& f = aux_edges[a];
    if (f.u != aux.aux_vertex[e.u] || f.v != aux.aux_vertex[e.v]) fail(id + "wrong endpoints");
    if (f.weight != sat::add(e.weight, cls.detour[i])) fail(id + "w' is not w + d_e");
    // Lemma 7.2: |J_x - J_y| <= w' <= 3w/2.
    if (J[f.v] < kInf && sat::add(J[f.v], f.weight) < J[f.u]) fail(id + "J_x > J_y + w'");
    if (J[f.u] < kInf && sat::add(J[f.u], f.weight) < J[f.v]) fail(id + "J_y > J_x + w'");
    if (sat::add(sat::add(e.weight, e.weight), e.weight) < twice(f.weight)) {
      fail(id + "w' > 3w/2");
    }
  }
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (cls.edge_class[i] != OuterEdgeClass::kRadial || aux.aux_edge[i] != kNoEdge) continue;
    if (!(aux.aux_vertex[edges[i].u] == kAuxRoot && aux.aux_vertex[edges[i].v] == kAuxRoot)) {
      fail("edge " + std::to_string(i) + ": radial edge missing from G'");
    }
  }

  if (J[kAuxRoot] != W{0}) fail("J_a != 0");
  std::vector<std::uint8_t> in_tree(aux_edges.size(), 0);
  for (const EdgeId a : aux.tree_edges) {
    if (a < in_tree.size()) in_tree[a] = 1;
  }
  std::size_t spokes = 0;
  for (const EdgeId id : od.tree_edges) {
    if (od.cut[id] != 0) continue;
    ++spokes;
    const EdgeId a = aux.aux_edge[id];
    if (a == kNoEdge || in_tree[a] == 0) fail("T': owner spoke " + std::to_string(id) + " missing");
  }
  if (aux.tree_edges.size() != spokes || spokes + 1 != na || aux.tree_order.size() != na ||
      aux.tree_parent.size() != na || aux.tree_parent_edge.size() != na ||
      aux.radius.size() != na || aux.owner_distance.size() != na) {
    fail("T': not a spanning tree of the owner spokes");
    return out;
  }
  std::vector<std::size_t> depth(na, 0);
  std::vector<std::uint8_t> placed(na, 0);
  if (aux.tree_order.front() != kAuxRoot || aux.tree_parent[kAuxRoot] != kNoVertex) {
    fail("T': not rooted at a");
    return out;
  }
  placed[kAuxRoot] = 1;
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex x = aux.tree_order[k];
    const std::string id = "aux vertex " + std::to_string(x) + ": ";
    if (x >= na || placed[x] != 0) {
      fail("T': tree order is not a permutation");
      return out;
    }
    const Vertex p = aux.tree_parent[x];
    const EdgeId a = aux.tree_parent_edge[x];
    if (p >= na || placed[p] == 0 || a >= aux_edges.size() || in_tree[a] == 0 ||
        !((aux_edges[a].u == x && aux_edges[a].v == p) ||
          (aux_edges[a].u == p && aux_edges[a].v == x))) {
      fail(id + "bad parent in T'");
      return out;
    }
    placed[x] = 1;
    depth[x] = depth[p] + 1;
    if (depth[x] > h) fail(id + "deeper than h in T'");
    const W w = aux_edges[a].weight;
    if (aux.radius[x] != std::max(aux.radius[p], w)) fail(id + "r is not the T' root bottleneck");
    if (aux.owner_distance[x] != sat::add(aux.owner_distance[p], w)) {
      fail(id + "R is not the T' root distance");
    }
    // Lemma 7.2: r <= J <= R <= h r.
    if (J[x] < aux.radius[x]) fail(id + "J < r");
    if (aux.owner_distance[x] < J[x]) fail(id + "J > R");
    if (sat::mul_int(aux.radius[x], static_cast<std::uint64_t>(h)) < aux.owner_distance[x]) {
      fail(id + "R > h r");
    }
  }
  // Lemma 7.2: T' is an MST of G' (cycle property; tree paths have <= 2h edges).
  for (std::size_t a = 0; a < aux_edges.size(); ++a) {
    if (in_tree[a] != 0) continue;
    Vertex x = aux_edges[a].u;
    Vertex y = aux_edges[a].v;
    W heaviest{0};
    while (x != y) {
      Vertex& deeper = depth[x] >= depth[y] ? x : y;
      heaviest = std::max(heaviest, aux_edges[aux.tree_parent_edge[deeper]].weight);
      deeper = aux.tree_parent[deeper];
    }
    if (aux_edges[a].weight < heaviest) {
      fail("aux edge " + std::to_string(a) + ": violates the MST cycle property of T'");
    }
  }
  return out;
}

#define SSSP_AUX_METRIC_INSTANTIATE(W)                                                      \
  template ProjectedEdges<W> projected_interowner_edges<W>(const CsrGraph<W>&,              \
                                                           const OwnerDecomposition<W>&);   \
  template OuterEdgeClassification<W> classify_outer_edges<W>(                              \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, const ProjectedEdges<W>&,           \
      const pivots::RouteNames<W>&, Diagnostics*);                                          \
  template OuterEdgeClassification<W> classify_interowner_edges<W>(                         \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, model::WordRng&,                    \
      const pivots::NameRoutesParams&, Diagnostics*);                                       \
  template AuxGraph<W> build_aux_graph<W>(const CsrGraph<W>&, const OwnerDecomposition<W>&, \
                                          const OuterEdgeClassification<W>&, Diagnostics*); \
  template std::vector<W> owner_coordinates_ref<W>(const AuxGraph<W>&);                     \
  template std::vector<std::string> check_aux_metric<W>(                                    \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, const OuterEdgeClassification<W>&,  \
      const AuxGraph<W>&, const std::vector<W>&);

SSSP_AUX_METRIC_INSTANTIATE(std::uint64_t)
SSSP_AUX_METRIC_INSTANTIATE(double)

#undef SSSP_AUX_METRIC_INSTANTIATE

}  // namespace sssp::paper::owner
