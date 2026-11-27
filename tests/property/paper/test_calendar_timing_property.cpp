// Sec 7.3 calendar and atom invariants that involve the owner coordinate J (txt l.1692-1779,
// Lemmas 7.5 and 7.6) of arXiv 2611.34911v1 [P4-PROPERTY-LEMMAS].
//
// Already covered elsewhere (not repeated here): test_calendars_property.cpp recomputes every
// Sec 7.3 object (overlap, pieces, grids, q, U+, calendars, chains, interpiece arcs, atoms) from
// its definition and compares it with build_calendars / find_atoms, and checks the size bounds of
// Lemmas 7.3 and 7.4; test_owner_distances_property.cpp checks Alg 7's use of them.
// New here, with J_v = d_{w'}(v, a) from an O(n^2) Dijkstra on G' written in this file (exact
// integers; quarter-integer doubles times 4):
//   * every nonroot v has its cell: i = floor(2^7 J_v / q_v) satisfies 32 <= i and
//     (i + 1) q_v <= U+_v, the calendar of (C_v, q_v) has cell i, and its two markers carry the
//     values i q_v and (i + 1) q_v (l.1699-1700: cells from zero up to the largest U+);
//   * Lemma 7.6 proof (l.1768-1773): for every interpiece arc U+_u -> 16 q_v, J_u < J_v, the cell
//     of J_u ends by the first marker, the cell of J_v starts at or beyond the second marker,
//     and the arc increases the marker value; the same holds for every J allowed when the arcs
//     are built (r <= J <= R): the cell of R_u ends by the first marker and the cell of r_v
//     starts at or beyond the second;
//   * l.1709: every edge of G' between two different atoms (both ends off the root) has
//     w' >= 2 max(q_u, q_v), i.e. at least twice either endpoint's cell width;
//   * Lemma 7.5 proof: a defining edge (w' < 2 max(q_u, q_v)) has r_u = r_v, and within an
//     atom Y every pair satisfies |J_x - J_y| <= L_Y (the discovery tree is a path in G');
//   * inputs: every connected graph on N <= 4 vertices with a two-scale weight set, every b, and
//     both name sources (exact names with P = 1, and NameRoutes with the paper's parameters);
//     seeded random clustered graphs through the shrinking driver.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/calendars.hpp"
#include "paper/pivots/name_routes.hpp"
#include "property/paper/shrink_support.hpp"
#include "sssp/graph.hpp"
#include "support/enumerate.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
namespace ow = sssp::paper::owner;
namespace pv = sssp::paper::pivots;
namespace sh = sssp::test::shrink;

namespace {

using Exact = unsigned __int128;
constexpr Exact kExactInf = ~Exact{0};

template <class W>
Exact exact(W x) {
  if constexpr (std::is_floating_point_v<W>) {
    const W y = x * W{4};
    if (!(y == std::floor(y)) || !(y < W{0x1p100})) throw std::domain_error("not exact");
    return static_cast<Exact>(y);
  } else {
    return static_cast<Exact>(x);
  }
}

std::string str(Exact x) {
  std::string s;
  do {
    s.insert(s.begin(), static_cast<char>('0' + static_cast<int>(x % 10)));
    x /= 10;
  } while (x != 0);
  return s;
}

/// Exact names with P = 1: the least instance weight lambda <= A_e with lambda >= d_{G_P}(x, y)
/// (Definition 6.1 with the strongest valid certificate).
template <class W>
pv::RouteNames<W> exact_names(const pv::PivotRoutingInstance<W>& inst) {
  const auto dist = sssp::test::floyd_warshall<W>(inst.t, inst.edges);
  std::vector<W> all;
  for (const auto& e : inst.edges) all.push_back(e.weight);
  std::sort(all.begin(), all.end());
  pv::RouteNames<W> names;
  std::vector<W> chosen;
  for (const auto& e : inst.edges) {
    W best = e.weight;
    for (const W a : all) {
      if (e.weight < a) break;
      if (!(a < dist[e.u][e.v])) {
        best = a;
        break;
      }
    }
    chosen.push_back(best);
  }
  names.dictionary = chosen;
  std::sort(names.dictionary.begin(), names.dictionary.end());
  names.dictionary.erase(std::unique(names.dictionary.begin(), names.dictionary.end()),
                         names.dictionary.end());
  for (const W c : chosen) {
    names.name.push_back(static_cast<std::uint32_t>(
        std::lower_bound(names.dictionary.begin(), names.dictionary.end(), c) -
        names.dictionary.begin()));
  }
  names.p_log2 = 0;
  return names;
}

struct Coverage {
  std::uint64_t cases = 0;
  std::uint64_t cells = 0;
  std::uint64_t interpiece_arcs = 0;
  std::uint64_t interatom_edges = 0;
  std::uint64_t defining_edges = 0;
  std::uint64_t multi_atoms = 0;
  std::uint64_t atom_spread = 0;  // atoms whose members have different J
};

template <class W>
std::optional<std::string> calendar_violation(const sh::Case<W>& c, bool exact_p1,
                                              std::uint64_t seed, Coverage* cov) {
  if (c.b == 0 || !sh::connected(c.n, c.edges)) return std::nullopt;
  for (const auto& e : c.edges) {
    if (!(W{0} < e.weight)) return std::nullopt;
  }
  const auto g = CsrGraph<W>::from_edges(c.n, c.edges);
  const auto od = sssp::paper::owner_decomposition(g, c.b);
  const auto pe = ow::projected_interowner_edges(g, od);
  ow::OuterEdgeClassification<W> cls;
  if (exact_p1) {
    cls = ow::classify_outer_edges(g, od, pe, exact_names(pe.instance));
  } else {
    sssp::paper::model::WordRng rng(seed);
    cls = ow::classify_interowner_edges(g, od, rng);
  }
  const auto aux = ow::build_aux_graph(g, od, cls);
  const auto cal = ow::build_calendars(aux);
  const auto atoms = ow::find_atoms(aux, cal);
  const std::size_t na = aux.num_vertices();
  std::ostringstream err;
  const auto fail = [&](const std::string& what) {
    err << what;
    return std::optional<std::string>(err.str());
  };
  if (na == 0 || aux.tree_parent[0] != sssp::kNoVertex) return fail("aux vertex 0 is not the root");

  // J by Dijkstra on G' from the root (O(n^2), exact).
  const auto ae = aux.graph.edges();
  std::vector<Exact> J(na, kExactInf);
  std::vector<bool> done(na, false);
  J[0] = 0;
  for (;;) {
    std::size_t best = na;
    for (std::size_t v = 0; v < na; ++v) {
      if (!done[v] && J[v] != kExactInf && (best == na || J[v] < J[best])) best = v;
    }
    if (best == na) break;
    done[best] = true;
    for (const auto& e : ae) {
      const Exact w = exact(e.weight);
      if (e.u == best && J[best] + w < J[e.v]) J[e.v] = J[best] + w;
      if (e.v == best && J[best] + w < J[e.u]) J[e.u] = J[best] + w;
    }
  }
  constexpr unsigned kScale = ow::kCalendarScaleLog2;  // calendar values are 2^7 lengths
  const auto marker = [&](std::uint32_t id) { return exact(cal.marker_value.at(id)); };

  for (Vertex v = 1; v < na; ++v) {
    err.str("");
    err << "vertex " << v << ": ";
    if (J[v] == kExactInf) return fail("unreachable in G'");
    const Exact r = exact(aux.radius[v]);
    const Exact R = exact(aux.owner_distance[v]);
    if (!(r <= J[v] && J[v] <= R)) return fail("r <= J <= R fails (Lemma 7.2)");
    const Exact q = exact(cal.width[v]);
    const Exact up = exact(cal.upper[v]);
    const Exact sj = J[v] << kScale;
    if (q == 0) return fail("zero width");
    const Exact i = sj / q;
    if (i < 32) return fail("J below 32 q");
    if ((i + 1) * q > up) return fail("the cell of J ends above U+");
    const auto& wc = cal.calendars.at(cal.calendar[v]);
    if (exact(wc.width) != q) return fail("calendar width differs from q_v");
    if (i >= wc.num_cells()) return fail("the calendar has no cell for J");
    if (marker(wc.markers[static_cast<std::size_t>(i)]) != i * q ||
        marker(wc.markers[static_cast<std::size_t>(i) + 1]) != (i + 1) * q) {
      return fail("the markers of the cell of J are not i q and (i + 1) q");
    }
    if (cov != nullptr) ++cov->cells;
  }
  for (const auto& a : cal.interpiece_arcs) {
    err.str("");
    err << "interpiece arc " << a.tail << " -> " << a.head << ": ";
    const Vertex u = a.tail;
    const Vertex v = a.head;
    if (u == 0 || v == 0 || u >= na || v >= na) return fail("endpoint is the root or invalid");
    if (cal.graph_piece[u] == cal.graph_piece[v]) return fail("both ends in one graph piece");
    if (!(J[u] < J[v])) return fail("J_u < J_v fails");
    const Exact qu = exact(cal.width[u]);
    const Exact qv = exact(cal.width[v]);
    const Exact end_u = ((J[u] << kScale) / qu + 1) * qu;
    const Exact start_v = ((J[v] << kScale) / qv) * qv;
    if (end_u > marker(a.from)) return fail("the cell of J_u ends after the first marker");
    if (marker(a.to) > start_v) return fail("the cell of J_v starts before the second marker");
    if (!(marker(a.from) < marker(a.to))) return fail("the arc does not increase time");
    // The same for every J allowed when the arcs are built (r <= J <= R, Lemma 7.2): the cell of
    // R_u ends by the first marker and the cell of r_v starts at or beyond the second.
    const Exact end_ru = ((exact(aux.owner_distance[u]) << kScale) / qu + 1) * qu;
    const Exact start_rv = ((exact(aux.radius[v]) << kScale) / qv) * qv;
    if (end_ru > marker(a.from)) return fail("the cell of R_u ends after the first marker");
    if (marker(a.to) > start_rv) return fail("the cell of r_v starts before the second marker");
    if (cov != nullptr) ++cov->interpiece_arcs;
  }
  for (const auto& e : ae) {
    if (e.u == 0 || e.v == 0 || e.u == e.v) continue;
    const Exact w = exact(e.weight) << kScale;
    const Exact qmax = std::max(exact(cal.width[e.u]), exact(cal.width[e.v]));
    err.str("");
    err << "aux edge " << e.u << "-" << e.v << ": ";
    if (w < 2 * qmax) {
      if (aux.radius[e.u] != aux.radius[e.v]) return fail("defining edge with unequal r");
      if (atoms.atom[e.u] != atoms.atom[e.v]) return fail("defining edge between two atoms");
      if (cov != nullptr) ++cov->defining_edges;
    } else if (atoms.atom[e.u] != atoms.atom[e.v] && cov != nullptr) {
      ++cov->interatom_edges;
    }
  }
  for (std::size_t y = 0; y < atoms.num_atoms(); ++y) {
    const auto mem = atoms.members(y);
    const Exact len = exact(atoms.length[y]);
    Exact lo = kExactInf;
    Exact hi = 0;
    for (const Vertex x : mem) {
      lo = std::min(lo, J[x]);
      hi = std::max(hi, J[x]);
    }
    if (hi - lo > len) {
      err.str("");
      err << "atom " << y << ": J spread " << str(hi - lo) << " > L_Y " << str(len);
      return err.str();
    }
    if (cov != nullptr) {
      cov->multi_atoms += mem.size() > 1 ? 1 : 0;
      cov->atom_spread += hi > lo ? 1 : 0;
    }
  }
  if (cov != nullptr) ++cov->cases;
  return std::nullopt;
}

template <class W>
std::vector<W> scaled(const std::vector<std::uint64_t>& xs) {
  std::vector<W> out;
  for (const auto x : xs) {
    if constexpr (std::is_floating_point_v<W>) {
      out.push_back(static_cast<W>(x) / W{4});
    } else {
      out.push_back(static_cast<W>(x));
    }
  }
  return out;
}

/// Clusters of light edges (weights 1..3) hung below vertex 0 by heavy spokes of very different
/// weights, plus random multi-scale chords: atoms, several graph pieces and interpiece arcs.
template <class W>
sh::Case<W> clustered(std::uint64_t seed) {
  sssp::test::SplitMix64 g(seed);
  sh::Case<W> c;
  c.n = 5 + g.below(26);
  c.b = 1 + g.below(4);
  const auto w = [&](std::uint64_t x) { return scaled<W>({x}).front(); };
  for (Vertex v = 1; v < c.n; ++v) {
    const bool spoke = g.below(4) == 0;
    const auto parent = static_cast<Vertex>(spoke ? 0 : g.below(v));
    c.edges.push_back({parent, v, w(spoke ? (100ULL << g.below(14)) : 1 + g.below(3))});
  }
  const std::size_t extra = g.below(c.n);
  for (std::size_t i = 0; i < extra; ++i) {
    const auto u = static_cast<Vertex>(g.below(c.n));
    const auto v = static_cast<Vertex>(g.below(c.n));
    if (u != v) c.edges.push_back({u, v, w((1 + g.below(3)) << (3 * g.below(5)))});
  }
  return c;
}

void require_coverage(const Coverage& cov) {
  WARN("cases " << cov.cases << ", cells of J checked " << cov.cells << ", interpiece arcs "
                << cov.interpiece_arcs << ", interatom edges " << cov.interatom_edges
                << ", defining edges " << cov.defining_edges << ", multi-vertex atoms "
                << cov.multi_atoms << ", atoms with J spread " << cov.atom_spread);
  CHECK(cov.interpiece_arcs > 0);
  CHECK(cov.interatom_edges > 0);
  CHECK(cov.defining_edges > 0);
  CHECK(cov.multi_atoms > 0);
  CHECK(cov.atom_spread > 0);
}

}  // namespace

TEMPLATE_TEST_CASE("Sec 7.3 cells of J, Lemma 7.6 arc timing and atoms on every graph with N <= 4",
                   "[property][calendars][lemma76]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  std::uint64_t graph = 0;
  for (std::size_t n = 2; n <= 4; ++n) {
    sssp::test::for_each_small_graph<W>(n, scaled<W>({1, 3, 700}), [&](const auto& edges) {
      ++graph;
      for (std::size_t b = 1; b <= n + 1; ++b) {
        for (const bool exact_p1 : {true, false}) {
          const sh::Case<W> c{n, edges, b};
          if (const auto v = calendar_violation(c, exact_p1, graph, &cov)) {
            FAIL(sh::describe(c) << (exact_p1 ? " exact-P1" : " name-routes") << " seed " << graph
                                 << ": " << *v);
          }
        }
      }
    });
  }
  CHECK(cov.cells > 0);
  CHECK(cov.defining_edges > 0);
  WARN("exhaustive: cases " << cov.cases << ", interpiece arcs " << cov.interpiece_arcs
                            << ", multi-vertex atoms " << cov.multi_atoms);
}

TEMPLATE_TEST_CASE("Sec 7.3 cells of J, Lemma 7.6 arc timing and atoms on seeded clustered graphs",
                   "[property][calendars][lemma76]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  for (const bool exact_p1 : {true, false}) {
    const auto out = sh::run_seeded<W>(1, 250, clustered<W>, [&](const sh::Case<W>& c) {
      return calendar_violation(c, exact_p1, 99, &cov);
    });
    INFO((exact_p1 ? "exact-P1: " : "name-routes: ") << out.report());
    CHECK_FALSE(out.failing_seed.has_value());
  }
  require_coverage(cov);
}
