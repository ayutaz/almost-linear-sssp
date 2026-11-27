// Property and differential tests for Sec 7.3 (Lemmas 7.3-7.6) of arXiv 2611.34911v1 on seeded
// graph families, against a brute-force reference written in this file in exact rational
// arithmetic (unsigned __int128).
//
// Inputs: the auxiliary graph G' of Sec 7.1-7.2 (owner_decomposition -> names ->
// classify_outer_edges -> build_aux_graph) for separator sizes b in {1, 2, 3, 4, ceil(log2 N), N,
// N + 1} (Mode B style, Q-03; the paper's b = L_N^64 gives one pivot, which b >= N covers) and
// two name sources: exact names with P = 1 (the strongest valid names of Definition 6.1, which
// keep many interowner edges radial) and NameRoutes with the paper's parameters.
//
// The reference recomputes every Sec 7.3 object from its definition with linear scans: overlap
// flags, tree pieces and graph pieces by label propagation, grids, 32 q <= r < 64 q and
// R < U+ <= 2R by scanning the grid, the calendars of the occurring (C, q), each chain as the
// sorted set of all values i q, the interpiece arcs U+_u -> 16 q_v, and the atoms by BFS over the
// defining edges. Lengths are converted exactly: std::uint64_t as is, double weights (quarter-
// integers, local_test_support.hpp) times 4, so the reference works on integers either way.
// Integer outputs must agree exactly. For double the grids, widths and upper markers must agree
// exactly (doublings); the chain values are sums and must agree within the recursive-summation
// bound (2M + 2) eps x of calendars.hpp (the number of inexact values is reported).
//
// The size bounds are the ones derived in the proofs: Lemma 7.4 gives
// cells < 136 #calendars + 256 #retained tree edges (U+ <= 2R, R_x <= 17 r_x / 16 < 68 q at a
// piece root, eligible edges < 64 q, sum over q > w/64 of w/q < 128), and Lemma 7.3 gives
// grid points <= 7 #C + 4.09 N (the union of [r/4, 4R] over a tree piece spans at most
// 4.0875 |P| octaves since R grows by at most 17x per retained edge).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::EdgeId;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::owner::AuxGraph;
using sssp::paper::pivots::PivotRoutingInstance;
using sssp::paper::pivots::RouteNames;
namespace lt = sssp::test::local;
namespace ow = sssp::paper::owner;

namespace {

/// REQUIRE that a checker returned no messages, showing the first ones otherwise.
#define REQUIRE_NO_MESSAGES(expr)                                 \
  do {                                                            \
    const auto msgs_ = (expr);                                    \
    for (std::size_t i_ = 0; i_ < msgs_.size() && i_ < 8; ++i_) { \
      UNSCOPED_INFO(msgs_[i_]);                                   \
    }                                                             \
    REQUIRE(msgs_.empty());                                       \
  } while (false)

using Exact = unsigned __int128;
constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

/// Exact integer image of a length: x (integers) or 4 x (quarter-integer doubles).
template <class W>
Exact exact(W x) {
  if constexpr (std::is_floating_point_v<W>) {
    const W y = x * W{4};
    REQUIRE(y == std::floor(y));
    REQUIRE(y < W{0x1p100});
    return static_cast<Exact>(y);
  } else {
    return static_cast<Exact>(x);
  }
}

double to_double(Exact x) {
  return static_cast<double>(x);
}

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

/// Exact names with P = 1 (see the file comment of test_aux_metric_property.cpp): the least
/// instance weight lambda <= a_e with lambda >= d_{G_P}(x, y).
template <class W>
RouteNames<W> exact_names(const PivotRoutingInstance<W>& inst) {
  const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
  std::vector<W> weights;
  for (const auto& e : inst.edges) weights.push_back(e.weight);
  std::ranges::sort(weights);
  std::vector<W> chosen;
  for (const auto& e : inst.edges) {
    W best = e.weight;
    for (const W a : weights) {
      if (a > e.weight) break;
      if (a >= dist[e.u][e.v]) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  RouteNames<W> names;
  names.p_log2 = 0;
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

/// Labels 0.. by first appearance in `order` of the classes of `label`.
std::vector<std::uint32_t> renumber(const std::vector<std::uint32_t>& label,
                                    const std::vector<Vertex>& order) {
  std::map<std::uint32_t, std::uint32_t> id;
  std::vector<std::uint32_t> out(label.size(), kNone);
  for (const Vertex v : order) {
    if (label[v] == kNone) continue;
    const auto [it, fresh] = id.try_emplace(label[v], static_cast<std::uint32_t>(id.size()));
    out[v] = it->second;
  }
  return out;
}

/// Components by label propagation to a fixpoint (O(n m), brute force).
std::vector<std::uint32_t> components(std::size_t n,
                                      const std::vector<std::pair<Vertex, Vertex>>& e,
                                      const std::vector<std::uint8_t>& active) {
  std::vector<std::uint32_t> label(n);
  for (std::size_t v = 0; v < n; ++v) label[v] = static_cast<std::uint32_t>(v);
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& [u, v] : e) {
      if (label[u] == label[v]) continue;
      const std::uint32_t m = std::min(label[u], label[v]);
      label[u] = label[v] = m;
      changed = true;
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (active[v] == 0) label[v] = kNone;
  }
  return label;
}

/// Light clusters under heavy spokes: vertex 0 joins c clusters by spokes of weight
/// 1000 * 2^U[0, 10]; each cluster is a random tree plus extra edges of weight 1..3 (quarter-
/// integers for double). With one pivot every cluster edge is far below 2 q, so the atoms are
/// whole clusters and their extra edges are defining chords.
template <class W>
std::vector<Edge<W>> heavy_spokes(std::size_t n, std::uint64_t seed) {
  sssp::test::SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  if (n <= 1) return e;
  const std::size_t clusters = 1 + rng.below(3);
  std::vector<std::vector<Vertex>> members(clusters);
  for (Vertex v = 1; v < n; ++v) members[rng.below(clusters)].push_back(v);
  for (const auto& mem : members) {
    if (mem.empty()) continue;
    e.push_back({0, mem[0], lt::to_weight<W>(1000ULL << rng.between(0, 10), rng)});
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], lt::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   lt::to_weight<W>(rng.between(1, 3), rng)});
    }
  }
  return e;
}

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t vertices = 0;
  std::uint64_t interpiece_arcs = 0;
  std::uint64_t pieces_with_two_widths = 0;
  std::uint64_t tree_pieces_sharing_c = 0;
  std::uint64_t shared_markers = 0;
  std::uint64_t multi_atoms = 0;
  std::uint64_t defining_chords = 0;
  std::uint64_t inexact_markers = 0;
  double max_cells_per_n = 0;
  double max_grid_per_n = 0;
};

template <class W>
void check_instance(const std::vector<Edge<W>>& edges, std::size_t n, std::size_t b, bool exact_p1,
                    std::uint64_t seed, Coverage& cov) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  const auto od = sssp::paper::owner_decomposition(g, b);
  const auto pe = ow::projected_interowner_edges(g, od);
  Diagnostics diag;
  ow::OuterEdgeClassification<W> cls;
  if (exact_p1) {
    cls = ow::classify_outer_edges(g, od, pe, exact_names(pe.instance));
  } else {
    sssp::paper::model::WordRng rng(seed);
    cls = ow::classify_interowner_edges(g, od, rng);
  }
  const AuxGraph<W> aux = ow::build_aux_graph(g, od, cls);
  const auto cal = ow::build_calendars(aux, &diag);
  const auto atoms = ow::find_atoms(aux, cal, &diag);
  REQUIRE_NO_MESSAGES(ow::check_calendars(aux, cal));
  REQUIRE_NO_MESSAGES(ow::check_atoms(aux, cal, atoms));

  const std::size_t na = aux.num_vertices();
  const auto ae = aux.graph.edges();
  const std::vector<Vertex>& order = aux.tree_order;
  std::vector<Exact> r(na, 0);
  std::vector<Exact> R(na, 0);
  for (Vertex v = 1; v < na; ++v) {
    r[v] = exact(aux.radius[v]);
    R[v] = exact(aux.owner_distance[v]);
  }
  std::vector<Exact> w(ae.size());
  for (std::size_t i = 0; i < ae.size(); ++i) w[i] = exact(ae[i].weight);
  std::vector<std::uint8_t> nonroot(na, 1);
  nonroot[0] = 0;

  // Overlap of the closed intervals [r/4, 4R] off the root.
  std::vector<std::uint8_t> overlap(ae.size(), 0);
  for (std::size_t i = 0; i < ae.size(); ++i) {
    const Vertex u = ae[i].u;
    const Vertex v = ae[i].v;
    if (u == 0 || v == 0 || u == v) continue;
    overlap[i] = (r[u] <= 16 * R[v] && r[v] <= 16 * R[u]) ? 1 : 0;
  }
  REQUIRE(cal.overlap == overlap);

  // Tree pieces: retained parent-child edges of T' - a.
  std::vector<std::pair<Vertex, Vertex>> kept_tree;
  std::size_t retained = 0;
  for (Vertex v = 1; v < na; ++v) {
    const Vertex p = aux.tree_parent[v];
    if (p != 0 && overlap[aux.tree_parent_edge[v]] != 0) {
      kept_tree.emplace_back(p, v);
      ++retained;
    }
  }
  const auto tp = renumber(components(na, kept_tree, nonroot), order);
  REQUIRE(cal.tree_piece == tp);
  std::vector<Vertex> roots;
  for (const Vertex v : order) {
    if (v != 0 && (aux.tree_parent[v] == 0 || tp[aux.tree_parent[v]] != tp[v])) roots.push_back(v);
  }
  REQUIRE(cal.tree_piece_root == roots);

  // Graph pieces: all overlapping edges.
  std::vector<std::pair<Vertex, Vertex>> kept_graph;
  for (std::size_t i = 0; i < ae.size(); ++i) {
    if (overlap[i] != 0) kept_graph.emplace_back(ae[i].u, ae[i].v);
  }
  const auto gp = renumber(components(na, kept_graph, nonroot), order);
  REQUIRE(cal.graph_piece == gp);
  const std::size_t num_c = cal.graph_pieces.size();
  std::size_t num_c_ref = 0;
  for (Vertex v = 1; v < na; ++v) num_c_ref = std::max<std::size_t>(num_c_ref, gp[v] + 1);
  REQUIRE(num_c == num_c_ref);
  for (std::size_t P = 0; P < roots.size(); ++P) {
    REQUIRE(cal.tree_piece_graph_piece[P] == gp[roots[P]]);
  }
  if (roots.size() > num_c) ++cov.tree_pieces_sharing_c;

  // Grids, widths and upper markers by linear scans (scaled exact units: 2^7 x).
  std::vector<Exact> origin(num_c, 0);
  std::vector<Exact> maxR(num_c, 0);
  std::vector<std::vector<Exact>> grid(num_c);
  for (Vertex v = 1; v < na; ++v) {
    const std::uint32_t c = gp[v];
    if (origin[c] == 0 || r[v] < origin[c]) origin[c] = r[v];
    maxR[c] = std::max(maxR[c], R[v]);
  }
  for (std::size_t c = 0; c < num_c; ++c) {
    for (Exact x = origin[c];; x *= 2) {
      grid[c].push_back(x);
      if (x >= 512 * maxR[c]) break;  // 4 max R, scaled by 2^7
    }
    const auto& got = cal.graph_pieces[c].grid;
    REQUIRE(got.size() == grid[c].size());
    for (std::size_t j = 0; j < got.size(); ++j) REQUIRE(exact(got[j]) == grid[c][j]);
    REQUIRE(exact(cal.graph_pieces[c].origin) == origin[c]);
  }
  std::vector<std::uint32_t> wj(na, kNone);
  std::vector<std::uint32_t> uk(na, kNone);
  for (Vertex v = 1; v < na; ++v) {
    INFO("aux vertex " << v);
    const auto& gr = grid[gp[v]];
    std::uint32_t i = 0;
    while (i + 1 < gr.size() && gr[i + 1] <= 128 * r[v]) ++i;
    REQUIRE(i >= 5);
    wj[v] = i - 5;
    std::uint32_t k = 0;
    while (gr[k] <= 128 * R[v]) ++k;
    uk[v] = k;
    // 32 q <= r < 64 q and R < U+ <= 2R, exactly.
    REQUIRE(32 * gr[wj[v]] <= 128 * r[v]);
    REQUIRE(128 * r[v] < 64 * gr[wj[v]]);
    REQUIRE(128 * R[v] < gr[k]);
    REQUIRE(gr[k] <= 256 * R[v]);
  }
  for (Vertex v = 1; v < na; ++v) {
    REQUIRE(cal.width_index[v] == wj[v]);
    REQUIRE(cal.upper_index[v] == uk[v]);
    REQUIRE(exact(cal.width[v]) == grid[gp[v]][wj[v]]);
    REQUIRE(exact(cal.upper[v]) == grid[gp[v]][uk[v]]);
  }

  // Calendars of the occurring (C, q), coarse to fine; chains = sorted sets of all i q.
  std::vector<std::map<std::uint32_t, std::uint32_t, std::greater<>>> cal_ref(num_c);
  for (Vertex v = 1; v < na; ++v) {
    auto [it, fresh] = cal_ref[gp[v]].try_emplace(wj[v], uk[v]);
    it->second = std::max(it->second, uk[v]);
  }
  std::size_t next_cal = 0;
  std::size_t marker_base = 0;
  std::size_t calendar_markers = 0;
  std::vector<std::vector<Exact>> chain(num_c);
  for (std::size_t c = 0; c < num_c; ++c) {
    if (cal_ref[c].size() > 1) ++cov.pieces_with_two_widths;
    std::set<Exact> values;
    for (const auto& [j, k] : cal_ref[c]) {
      const Exact q = grid[c][j];
      for (Exact x = 0; x <= grid[c][k]; x += q) values.insert(x);
    }
    chain[c].assign(values.begin(), values.end());
    const auto& gpc = cal.graph_pieces[c];
    REQUIRE(gpc.first_marker == marker_base);
    REQUIRE(gpc.end_marker == marker_base + chain[c].size());
    const std::size_t M = chain[c].size();
    for (std::size_t m = 0; m < M; ++m) {
      const W got = cal.marker_value[marker_base + m];
      const double want = to_double(chain[c][m]);
      if constexpr (std::is_integral_v<W>) {
        REQUIRE(exact(got) == chain[c][m]);
      } else {
        const double f = got * 4.0;
        if (f != want) ++cov.inexact_markers;
        REQUIRE(std::abs(f - want) <=
                static_cast<double>(2 * M + 2) * std::numeric_limits<double>::epsilon() * want);
      }
      REQUIRE(cal.marker_graph_piece[marker_base + m] == c);
    }
    REQUIRE(gpc.calendars.size() == cal_ref[c].size());
    std::size_t k_idx = 0;
    for (const auto& [j, k] : cal_ref[c]) {
      const auto& wc = cal.calendars[next_cal];
      REQUIRE(gpc.calendars[k_idx++] == next_cal);
      REQUIRE(wc.graph_piece == c);
      REQUIRE(wc.grid_index == j);
      const Exact q = grid[c][j];
      REQUIRE(Exact{wc.num_cells()} * q == grid[c][k]);
      for (std::size_t i = 0; i < wc.markers.size(); ++i) {
        const auto it = std::ranges::lower_bound(chain[c], Exact{i} * q);
        REQUIRE(wc.markers[i] == marker_base + static_cast<std::size_t>(it - chain[c].begin()));
      }
      calendar_markers += wc.markers.size();
      ++next_cal;
    }
    marker_base += M;
  }
  REQUIRE(cal.calendars.size() == next_cal);
  REQUIRE(cal.marker_value.size() == marker_base);
  cov.shared_markers += calendar_markers - marker_base;
  for (Vertex v = 1; v < na; ++v) {
    const auto& wc = cal.calendars[cal.calendar[v]];
    REQUIRE(wc.graph_piece == gp[v]);
    REQUIRE(wc.grid_index == wj[v]);
  }

  // Interpiece arcs U+_u -> 16 q_v with r_v > 16 R_u.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> arcs;
  for (const auto& c : cal.graph_pieces) {
    for (std::uint32_t m = c.first_marker; m + 1 < c.end_marker; ++m) arcs.emplace_back(m, m + 1);
  }
  std::size_t next_arc = 0;
  const auto marker_of = [&](std::uint32_t c, Exact x) {
    const auto it = std::ranges::lower_bound(chain[c], x);
    REQUIRE(it != chain[c].end());
    REQUIRE(*it == x);
    return static_cast<std::uint32_t>(cal.graph_pieces[c].first_marker + (it - chain[c].begin()));
  };
  for (std::size_t i = 0; i < ae.size(); ++i) {
    Vertex u = ae[i].u;
    Vertex v = ae[i].v;
    if (u == 0 || v == 0 || u == v || gp[u] == gp[v]) continue;
    if (!(r[v] > 16 * R[u])) std::swap(u, v);
    REQUIRE(r[v] > 16 * R[u]);
    REQUIRE(next_arc < cal.interpiece_arcs.size());
    const auto& a = cal.interpiece_arcs[next_arc++];
    REQUIRE(a.aux_edge == i);
    REQUIRE(a.tail == u);
    REQUIRE(a.head == v);
    REQUIRE(a.from == marker_of(gp[u], grid[gp[u]][uk[u]]));
    REQUIRE(a.to == marker_of(gp[v], 16 * grid[gp[v]][wj[v]]));
    // Lemma 7.6: U+_u <= 2 R_u < r_v / 8 < 16 q_v.
    REQUIRE(grid[gp[u]][uk[u]] < 16 * grid[gp[v]][wj[v]]);
    arcs.emplace_back(a.from, a.to);
  }
  REQUIRE(next_arc == cal.interpiece_arcs.size());
  cov.interpiece_arcs += next_arc;
  // The marker DAG: marker_order is topological for the reference arcs.
  {
    std::vector<std::size_t> pos(marker_base, marker_base);
    REQUIRE(cal.marker_order.size() == marker_base);
    for (std::size_t k = 0; k < marker_base; ++k) {
      REQUIRE(pos[cal.marker_order[k]] == marker_base);
      pos[cal.marker_order[k]] = k;
    }
    for (const auto& [x, y] : arcs) REQUIRE(pos[x] < pos[y]);
  }

  // Lemma 7.4 and Lemma 7.3 size bounds (see the file comment).
  const std::size_t N = na - 1;
  REQUIRE(cal.num_cells() < 136 * cal.calendars.size() + 256 * retained + 1);
  REQUIRE(100 * cal.num_grid_points() <= 700 * num_c + 409 * N);
  if (N > 0) {
    cov.max_cells_per_n = std::max(cov.max_cells_per_n,
                                   static_cast<double>(cal.num_cells()) / static_cast<double>(N));
    cov.max_grid_per_n = std::max(
        cov.max_grid_per_n, static_cast<double>(cal.num_grid_points()) / static_cast<double>(N));
  }
#if SSSP_DIAGNOSTICS
  REQUIRE(diag.counters["calendars.cells"] == cal.num_cells());
  REQUIRE(diag.counters["calendars.grid_points"] == cal.num_grid_points());
  // Lemma 7.3: O(1 + log log(b + 4)) comparisons per vertex (two binary searches in windows of
  // at most 2 (4 + ceil log2 h) + 1 and 1 + ceil log2 h indices).
  {
    std::uint64_t c_h = 0;
    while (c_h < 64 && (std::uint64_t{1} << c_h) < aux.hop_bound) ++c_h;
    std::uint64_t lg = 0;
    while ((std::uint64_t{1} << lg) < 2 * (4 + c_h) + 2) ++lg;
    REQUIRE(diag.counters["calendars.grid_search_steps"] <= 2 * lg * N);
  }
#endif

  // Atoms: BFS components of the defining edges 2^7 w' < 2 max(q_u, q_v).
  std::vector<std::pair<Vertex, Vertex>> defining;
  for (std::size_t i = 0; i < ae.size(); ++i) {
    const Vertex u = ae[i].u;
    const Vertex v = ae[i].v;
    if (u == 0 || v == 0 || u == v) continue;
    const Exact q = std::max(grid[gp[u]][wj[u]], grid[gp[v]][wj[v]]);
    if (128 * w[i] < 2 * q) {
      defining.emplace_back(u, v);
      REQUIRE(r[u] == r[v]);  // Lemma 7.5
    }
  }
  const auto atom_ref = renumber(components(na, defining, nonroot), order);
  REQUIRE(atoms.atom == atom_ref);
  REQUIRE(atoms.tree_edges.size() == N - atoms.num_atoms());
  cov.defining_chords += defining.size() - atoms.tree_edges.size();
  for (std::size_t a = 0; a < atoms.num_atoms(); ++a) {
    const auto mem = atoms.members(a);
    if (mem.size() > 1) ++cov.multi_atoms;
    REQUIRE(mem.size() < b);  // Lemma 7.5 with Lemma 5.2
    Exact len = 0;
    for (std::uint32_t k = atoms.tree_offsets[a]; k < atoms.tree_offsets[a + 1]; ++k) {
      len += w[atoms.tree_edges[k]];
    }
    REQUIRE(exact(atoms.length[a]) == len);
    const Vertex rep = atoms.representative[a];
    REQUIRE(rep == mem.front());
    // L_Y <= 2 q |Y| and R_v <= R_rep + L_Y.
    REQUIRE(128 * len <= 2 * grid[gp[rep]][wj[rep]] * mem.size());
    for (const Vertex v : mem) REQUIRE(R[v] <= R[rep] + len);
  }
  ++cov.cases;
  cov.vertices += N;
}

}  // namespace

TEMPLATE_TEST_CASE("Sec 7.3 calendars and atoms against an exact brute-force reference",
                   "[property][calendars]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights weights[] = {lt::Weights::kTies, lt::Weights::kUniform,
                                 lt::Weights::kExponential};
  const std::size_t sizes[] = {1, 2, 4, 7, 12, 20, 33};
  for (const auto f : families) {
    for (const auto wd : weights) {
      for (const std::size_t n : sizes) {
        for (std::uint64_t seed = 1; seed <= 2; ++seed) {
          const std::uint64_t gseed = seed * 1000003 + n * 31 + static_cast<std::uint64_t>(f) * 7 +
                                      static_cast<std::uint64_t>(wd);
          const auto edges = lt::make_graph<W>(f, wd, n, gseed);
          for (const std::size_t b : separator_sizes(n)) {
            for (const bool exact_p1 : {true, false}) {
              INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " seed=" << gseed
                               << " b=" << b << (exact_p1 ? " exact-P1" : " name-routes-paper"));
              check_instance(edges, n, b, exact_p1, gseed + b, cov);
            }
          }
        }
      }
    }
  }
  for (const std::size_t n : sizes) {
    for (std::uint64_t seed = 1; seed <= 4; ++seed) {
      const std::uint64_t gseed = seed * 7919 + n;
      const auto edges = heavy_spokes<W>(n, gseed);
      for (const std::size_t b : separator_sizes(n)) {
        for (const bool exact_p1 : {true, false}) {
          INFO("heavy-spokes n=" << n << " seed=" << gseed << " b=" << b
                                 << (exact_p1 ? " exact-P1" : " name-routes-paper"));
          check_instance(edges, n, b, exact_p1, gseed + b, cov);
        }
      }
    }
  }
  INFO("cases " << cov.cases << ", vertices " << cov.vertices << ", interpiece arcs "
                << cov.interpiece_arcs << ", graph pieces with >= 2 widths "
                << cov.pieces_with_two_widths << ", shared markers " << cov.shared_markers
                << ", multi-vertex atoms " << cov.multi_atoms << ", defining chords "
                << cov.defining_chords << ", inexact double markers " << cov.inexact_markers
                << ", max cells/N " << cov.max_cells_per_n << ", max grid points/N "
                << cov.max_grid_per_n);
  CHECK(cov.interpiece_arcs > 0);
  CHECK(cov.pieces_with_two_widths > 0);
  CHECK(cov.tree_pieces_sharing_c > 0);
  CHECK(cov.shared_markers > 0);
  CHECK(cov.multi_atoms > 0);
  CHECK(cov.defining_chords > 0);
  // Report of the measured ratios (WARN does not fail the test).
  WARN("calendars: cases " << cov.cases << ", max cells/N " << cov.max_cells_per_n
                           << ", max grid points/N " << cov.max_grid_per_n
                           << ", inexact double markers " << cov.inexact_markers);
}

TEST_CASE("double calendars with 53-bit weights satisfy the checkers within the rounding bound",
          "[property][calendars]") {
  // Weights with full 53-bit mantissas: the cell sums i q round, so check_calendars compares the
  // markers with the recursive-summation bound of calendars.hpp (Q-06). Every other statement
  // (pieces, grids, 32 q <= r < 64 q, R < U+ <= 2R, arcs, atoms) is exact in double.
  std::uint64_t rounded = 0;
  std::uint64_t cases = 0;
  const lt::Family families[] = {lt::Family::kRandomTree, lt::Family::kSparse, lt::Family::kGrid};
  for (const auto f : families) {
    for (const std::size_t n : {7, 20, 33}) {
      for (std::uint64_t seed = 1; seed <= 3; ++seed) {
        auto edges = lt::make_graph<double>(f, lt::Weights::kExponential, n, seed);
        sssp::test::SplitMix64 rng(seed * 104729 + n);
        for (auto& e : edges) e.weight = e.weight * (1.0 + rng.unit()) / 3.0;
        const auto g = CsrGraph<double>::from_edges(n, edges);
        for (const std::size_t b : separator_sizes(n)) {
          INFO(lt::name(f) << " n=" << n << " seed=" << seed << " b=" << b);
          const auto od = sssp::paper::owner_decomposition(g, b);
          const auto pe = ow::projected_interowner_edges(g, od);
          const auto cls = ow::classify_outer_edges(g, od, pe, exact_names(pe.instance));
          const auto aux = ow::build_aux_graph(g, od, cls);
          const auto cal = ow::build_calendars(aux);
          const auto atoms = ow::find_atoms(aux, cal);
          REQUIRE_NO_MESSAGES(ow::check_calendars(aux, cal));
          REQUIRE_NO_MESSAGES(ow::check_atoms(aux, cal, atoms));
          for (const auto& wc : cal.calendars) {
            for (std::size_t i = 0; i < wc.markers.size(); ++i) {
              if (cal.marker_value[wc.markers[i]] != static_cast<double>(i) * wc.width) ++rounded;
            }
          }
          ++cases;
        }
      }
    }
  }
  INFO("cases " << cases << ", markers differing from the single-rounded i q " << rounded);
  CHECK(cases > 0);
  CHECK(rounded > 0);  // the tolerance branch of check_calendars is exercised
}
