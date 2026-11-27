// Property tests of the entry point of Algorithm 4 (LocalHierarchy, Sec 5.4 l.1244-1261,
// Lemma 5.10) with the paper's Mode A constants (kPaperLocalParams: k > 2 L_k^4, b = L_k^4),
// at sizes where the radial branch applies. The least such size is 209953 (L = 18,
// b = 18^4 = 104976); 262142 is the last size with L = 18 and 262143 the first with L = 19
// (b = 19^4 = 130321).
//
// Checked directly on the output of local_hierarchy(), without Evaluate:
//   * Alg 4 l.2 / fallback policy: no FallbackRecord on the radial branch (and, with diagnostic
//     counters, the radial branch was taken and the base case was not);
//   * Lemma 5.10: c0 <= 2; a native entry at every vertex (projecting to it); the hierarchy
//     satisfies Def 4.1, judged both by the project's validator and by the independent checker
//     in local_hierarchy_large_support.hpp;
//   * Lemma 5.10 "Every geodesic has an equal-length state realization from its native entry":
//     a Dijkstra on the state graph D from the native entry of s, minimised over the copies of
//     each vertex, equals the physical Dijkstra distance (both written in the test support);
//     every state distance is at least the physical distance of its image (H1);
//   * Sec 5.2 "always with the original physical length when an edge is used": a state arc that
//     names an original edge joins copies of its endpoints and has exactly its weight; an arc
//     that names no edge joins two copies of one vertex with length 0 (the phase switch);
//   * Lemma 5.7 / 5.10 size bounds: at most four copies per vertex (two inward and two outward
//     two-shell windows, Sec 5.2 l.1035-1041), fewer internal nodes than states (H2: at least two
//     children each), and at most 8m + 2N arcs (a constant number of copies per edge and
//     direction, Sec 5.2 "Each vertex and edge have a constant number of copies").
//
// These cases are tagged [slow]: each builds a hierarchy on >= 209953 vertices. They are not
// reduced for faster runs, because a smaller input falls below the Mode A threshold and takes
// the base case instead.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "property/paper/local_hierarchy_large_support.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::Edge;
using sssp::kNoEdge;
using sssp::Vertex;
namespace lt = sssp::test::local;
namespace ll = sssp::test::local_large;

namespace {

template <class W>
void check_radial_properties(std::size_t n, const std::vector<Edge<W>>& edges) {
  const auto g = CsrGraph<W>::from_edges(n, edges);
  const std::size_t m = g.num_edges();
  REQUIRE(sssp::paper::local_construction_applies(n, sssp::paper::kPaperLocalParams));

  Diagnostics diag;
  const auto lh = sssp::paper::local_hierarchy(g, &diag);
  const auto& h = lh.hierarchy;

  // Alg 4 l.6-12: the radial branch records no fallback.
  CHECK(diag.fallbacks.empty());
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters["local_hierarchy_calls"] == 1);
  CHECK(diag.counters["local_hierarchy_radial"] == 1);
  CHECK(diag.counters["local_hierarchy_base_case"] == 0);
  CHECK(diag.counters["local_base_case_calls"] == 0);
  CHECK(diag.counters["local_assemble_calls"] == 1);
  CHECK(diag.counters["local_hierarchy_states"] == h.num_states());
  CHECK(diag.counters["local_hierarchy_arcs"] == h.graph.arcs.size());
  CHECK(diag.counters["local_hierarchy_nodes"] == h.nodes.size());
#endif

  // Lemma 5.10: connection constant at most two; native entries.
  CHECK(h.c0 >= 1);
  CHECK(h.c0 <= 2);
  REQUIRE(lh.native_entry.size() == n);
  REQUIRE(h.projection.size() == h.num_states());
  std::vector<char> is_entry(h.num_states(), 0);
  for (Vertex v = 0; v < n; ++v) {
    const auto x = lh.native_entry[v];
    REQUIRE(x < h.num_states());
    REQUIRE(h.projection[x] == v);
    REQUIRE(is_entry[x] == 0);
    is_entry[x] = 1;
  }

  // Sizes (Sec 5.2, Lemma 5.7).
  CHECK(h.num_states() >= n);
  CHECK(h.num_states() <= 4 * n);
  CHECK(h.nodes.size() < 2 * h.num_states());
  CHECK(h.graph.arcs.size() <= 8 * m + 2 * n);

  // State arcs carry original lengths (Sec 5.2) or are zero-cost switches between copies.
  REQUIRE(lh.arc_edge.size() == h.graph.arcs.size());
  std::size_t bad_arcs = 0;
  for (std::size_t a = 0; a < h.graph.arcs.size(); ++a) {
    const auto& arc = h.graph.arcs[a];
    const Vertex pu = h.projection[arc.tail];
    const Vertex pv = h.projection[arc.head];
    const auto id = lh.arc_edge[a];
    bool ok = false;
    if (id == kNoEdge) {
      ok = pu == pv && arc.length == W{0};
    } else if (id < m) {
      const auto& e = edges[id];
      ok = ((e.u == pu && e.v == pv) || (e.u == pv && e.v == pu)) && arc.length == e.weight;
    }
    if (!ok) ++bad_arcs;
  }
  CHECK(bad_arcs == 0);

  // Def 4.1: the project's validator and the independent checker.
  const auto valid = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(valid.summary());
  CHECK(valid.ok());
  const auto rep = ll::check_def41(h, n, edges, 2);
  INFO(rep.summary());
  INFO("internal nodes " << rep.internal_nodes << ", depth " << rep.max_depth << ", links "
                         << rep.links << " (certificate-only " << rep.certificate_only_links
                         << "), light quotient arcs " << rep.light_quotient_arcs
                         << ", interroot arcs " << rep.interroot_arcs << ", walk searches "
                         << rep.walk_searches);
  CHECK(rep.ok());
  // Lemma 5.10 does not promise a single root (several SCCs of the final state graph are
  // processed in condensation order, Sec 5.4), so only a nonempty forest is required.
  CHECK(!h.roots.empty());

  // Lemma 5.10: path realization from the native entry, independent of Evaluate.
  const auto phys = ll::undirected(n, edges);
  const auto states = ll::state_digraph(h);
  for (const Vertex src : {Vertex{0}, static_cast<Vertex>(n / 2), static_cast<Vertex>(n - 1)}) {
    INFO("source " << src);
    const auto want = ll::dijkstra(phys, {{src, W{0}}});
    const auto ds = ll::dijkstra(states, {{lh.native_entry[src], W{0}}});
    std::vector<W> best(n, ll::infinity<W>());
    std::size_t below_physical = 0;
    for (std::size_t x = 0; x < ds.size(); ++x) {
      const Vertex v = h.projection[x];
      if (ds[x] < want[v]) ++below_physical;
      if (ds[x] < best[v]) best[v] = ds[x];
    }
    CHECK(below_physical == 0);
    std::size_t mismatches = 0;
    Vertex first_bad = sssp::kNoVertex;
    for (Vertex v = 0; v < n; ++v) {
      REQUIRE(want[v] < ll::infinity<W>());
      if (best[v] != want[v]) {
        if (mismatches == 0) first_bad = v;
        ++mismatches;
      }
    }
    INFO("first mismatching vertex " << first_bad);
    CHECK(mismatches == 0);
  }
}

}  // namespace

TEST_CASE("Alg 4 Mode A at k = 209953: path, uniform weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kPath, lt::Weights::kUniform, n, 11));
}

TEST_CASE("Alg 4 Mode A at k = 209953: random tree, tied weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kRandomTree, lt::Weights::kTies, n, 12));
}

TEST_CASE("Alg 4 Mode A at k = 209953: grid, exponential weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kGrid, lt::Weights::kExponential, n, 13));
}

TEST_CASE("Alg 4 Mode A at k = 209953: spider, uniform weights", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSpider, lt::Weights::kUniform, n, 14));
}

TEST_CASE("Alg 4 Mode A at k = 209953: clustered with heavy bridges", "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<std::uint64_t>(n, lt::clustered(n, 15));
}

TEST_CASE("Alg 4 Mode A at k = 209953: sparse multigraph, quarter-integer doubles",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstRadial;
  check_radial_properties<double>(
      n, lt::make_graph<double>(lt::Family::kSparse, lt::Weights::kUniform, n, 16));
}

TEST_CASE("Alg 4 Mode A at k = 262142 (last L = 18): sparse multigraph, exponential weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kLastL18;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSparse, lt::Weights::kExponential, n, 17));
}

TEST_CASE("Alg 4 Mode A at k = 262143 (first L = 19): sparse multigraph, tied weights",
          "[slow][local-large]") {
  const std::size_t n = ll::kFirstL19;
  check_radial_properties<std::uint64_t>(
      n, lt::make_graph<std::uint64_t>(lt::Family::kSparse, lt::Weights::kTies, n, 18));
}

namespace {

bool reports(const ll::Def41Report& rep, const std::string& prefix) {
  for (const auto& v : rep.violations) {
    if (v.rfind(prefix, 0) == 0) return true;
  }
  return false;
}

/// Some state below node c (the first leaf in child order).
template <class W>
sssp::paper::StateId some_state_below(const sssp::paper::OrderedHierarchy<W>& h,
                                      sssp::paper::NodeId c) {
  while (!h.nodes[c].children.empty()) c = h.nodes[c].children.front();
  return h.nodes[c].state;
}

}  // namespace

TEST_CASE("independent Def 4.1 checker detects each violated clause") {
  // The checker used above must not be vacuous: every clause is broken once on a small radial
  // hierarchy (Mode B separator b = 3 on a 6 x 6 grid, so it stays fast) and must be reported.
  using W = std::uint64_t;
  using sssp::paper::kNoNode;
  using sssp::paper::NodeId;
  std::vector<Edge<W>> e;
  for (Vertex v = 0; v < 36; ++v) {
    if (v % 6 != 5) e.push_back({v, v + 1, 1 + (v * 7) % 13});
    if (v + 6 < 36) e.push_back({v, v + 6, 1 + (v * 5) % 11});
  }
  const auto g = CsrGraph<W>::from_edges(36, e);
  const auto lh =
      sssp::paper::local_hierarchy(g, nullptr, sssp::paper::LocalHierarchyParams{1, 1, 3});
  const auto& h0 = lh.hierarchy;
  REQUIRE(ll::check_def41(h0, 36, e, 2).ok());

  // An internal node with at least two children and a parent.
  NodeId c = kNoNode;
  for (NodeId i = 0; i < h0.nodes.size(); ++i) {
    if (h0.nodes[i].children.size() >= 2 && h0.nodes[i].parent != kNoNode) {
      c = i;
      break;
    }
  }
  REQUIRE(c != kNoNode);
  const auto x = some_state_below(h0, h0.nodes[c].children[0]);
  const auto y = some_state_below(h0, h0.nodes[c].children[1]);

  SECTION("H1: an arc shorter than the physical distance") {
    auto h = h0;
    std::size_t a = 0;
    while (h.projection[h.graph.arcs[a].tail] == h.projection[h.graph.arcs[a].head]) ++a;
    h.graph.arcs[a].length = 0;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H1:"));
  }
  SECTION("H2: a non-positive height") {
    auto h = h0;
    h.nodes[c].height = 0;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H2:"));
  }
  SECTION("H2: a child higher than its parent") {
    auto h = h0;
    h.nodes[c].height = h.nodes[h.nodes[c].parent].height + 1;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H2:"));
  }
  SECTION("H3: a light two-cycle between distinct children") {
    auto h = h0;
    h.graph.arcs.push_back({x, y, 0});
    h.graph.arcs.push_back({y, x, 0});
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H3:"));
  }
  SECTION("H4: a cycle between two roots") {
    // The base case on two components has two roots.
    const std::vector<Edge<W>> two{{0, 1, 1}, {2, 3, 1}, {3, 4, 2}};
    const auto b = sssp::paper::local_hierarchy(CsrGraph<W>::from_edges(5, two));
    REQUIRE(b.hierarchy.roots.size() == 2);
    auto h = b.hierarchy;
    REQUIRE(ll::check_def41(h, 5, two, 1).ok());
    h.graph.arcs.push_back({b.native_entry[0], b.native_entry[2], 5});
    h.graph.arcs.push_back({b.native_entry[2], b.native_entry[0], 5});
    CHECK(reports(ll::check_def41(h, 5, two, 1), "H4:"));
  }
  SECTION("H5: a missing link") {
    auto h = h0;
    REQUIRE(!h.nodes[c].links.empty());
    h.nodes[c].links.pop_back();
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H5:"));
  }
  SECTION("H5: a link longer than c0 h(C)") {
    auto h = h0;
    h.nodes[c].links[0].length = 2 * h.nodes[c].height + 1;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H5: link longer"));
  }
  SECTION("H5: a link shorter than the physical distance") {
    auto h = h0;
    auto& link = h.nodes[c].links[0];
    if (h.projection[link.x] == h.projection[link.y]) {
      link.x = x;
      link.y = y;
    }
    REQUIRE(h.projection[link.x] != h.projection[link.y]);
    link.length = 0;
    link.arc = sssp::paper::kNoArc;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H5: link shorter"));
  }
  SECTION("H6: a backing arc that does not join the link endpoints") {
    auto h = h0;
    auto& link = h.nodes[c].links[0];
    sssp::paper::ArcId a = 0;
    while (h.graph.arcs[a].tail == link.x || h.graph.arcs[a].tail == link.y) ++a;
    link.arc = a;
    CHECK(reports(ll::check_def41(h, 36, e, 2), "H6:"));
  }
  SECTION("structure: a leaf map that points elsewhere") {
    auto h = h0;
    std::swap(h.leaf_of_state[0], h.leaf_of_state[1]);
    CHECK_FALSE(ll::check_def41(h, 36, e, 2).ok());
  }
}
