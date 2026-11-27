// Tests for the Def 4.1 / Def 4.2 checker (Sec 4, l.646-672).
// Every condition H1-H6 is violated one at a time on a hand-built valid example, and the
// checker must name the violated condition and its location. Physical distances that a test
// relies on are confirmed with an independent Floyd-Warshall written here.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::paper::ArcId;
using sssp::paper::Condition;
using sssp::paper::kNoArc;
using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::PhysicalLink;
using sssp::paper::PhysicalSource;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateArc;
using sssp::paper::StateGraph;
using sssp::paper::StateOffer;
using sssp::paper::ValidationReport;
using sssp::paper::ViolationKind;

namespace {

// Naive all-pairs distances on the undirected physical graph (test oracle only).
template <class W>
std::vector<std::vector<W>> floyd_warshall(std::size_t n, const std::vector<Edge<W>>& edges) {
  const W inf = sssp::WeightTraits<W>::infinity();
  std::vector<std::vector<W>> d(n, std::vector<W>(n, inf));
  for (std::size_t i = 0; i < n; ++i) d[i][i] = W{0};
  for (const auto& e : edges) {
    if (e.weight < d[e.u][e.v]) d[e.u][e.v] = d[e.v][e.u] = e.weight;
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        if (d[i][k] != inf && d[k][j] != inf && d[i][k] + d[k][j] < d[i][j])
          d[i][j] = d[i][k] + d[k][j];
  return d;
}

// Physical path 0 -1- 1 -2- 2 -4- 3; one state per vertex.
// Arc ids: 0: 0->1 (1), 1: 1->0 (6), 2: 1->2 (2), 3: 2->1 (5), 4: 2->3 (4), 5: 3->2 (4).
// Forest: A (h=1) = {L0, L1} with link (0,1,1); R (h=4) = {A, L2, L3} with links (1,2,2),
// (2,3,4). Light arcs at R (< 4 between distinct children): only 1->2 (A -> L2). Valid, c0 = 1.
template <class W>
struct Example {
  std::vector<Edge<W>> edges{{0, 1, W{1}}, {1, 2, W{2}}, {2, 3, W{4}}};
  CsrGraph<W> phys = CsrGraph<W>::from_edges(4, edges);
  OrderedHierarchy<W> h;
  NodeId a = 0;
  NodeId r = 0;

  Example() {
    StateGraph<W> d;
    d.num_states = 4;
    d.arcs = {{0, 1, W{1}}, {1, 0, W{6}}, {1, 2, W{2}}, {2, 1, W{5}}, {2, 3, W{4}}, {3, 2, W{4}}};
    h = sssp::paper::make_leaf_forest<W>(std::move(d), {0, 1, 2, 3});
    a = sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{1}, 0}});
    r = sssp::paper::add_internal_node<W>(h, W{4}, {a, 2, 3}, {{1, 2, W{2}, 2}, {2, 3, W{4}, 4}});
    sssp::paper::collect_roots(h);
  }

  [[nodiscard]] ValidationReport check() const {
    return sssp::paper::validate_ordered_hierarchy(h, phys);
  }
};

// Only conditions in `allowed` may appear, and `expected` must appear.
bool only(const ValidationReport& rep, Condition expected) {
  if (!rep.has(expected)) return false;
  for (const auto& v : rep.violations) {
    if (v.condition() != expected) return false;
  }
  return true;
}

}  // namespace

TEMPLATE_TEST_CASE("valid hand-built ordered hierarchy passes the Def 4.1 checker",
                   "[paper][hierarchy][def41]", std::uint64_t, double) {
  Example<TestType> ex;
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(rep.ok());
  REQUIRE(ex.h.roots == std::vector<NodeId>{ex.r});
  REQUIRE(ex.h.nodes[ex.a].parent == ex.r);
  REQUIRE(ex.h.leaf_of_state.size() == 4);
  for (sssp::paper::StateId x = 0; x < 4; ++x) {
    REQUIRE(ex.h.is_leaf(ex.h.leaf_of_state[x]));
    REQUIRE(ex.h.nodes[ex.h.leaf_of_state[x]].state == x);
  }
}

TEMPLATE_TEST_CASE("H1: an arc shorter than every physical walk is reported at that arc",
                   "[paper][hierarchy][def41][h1]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  const auto fw = floyd_warshall<W>(4, ex.edges);
  REQUIRE(fw[0][3] == W{7});
  // 0->3 of length 3 < d(0,3) = 7. It is light at R (3 < 4) but keeps R's quotient acyclic.
  ex.h.graph.arcs.push_back({0, 3, W{3}});
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH1));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kArcWalkTooLong);
  REQUIRE(rep.violations[0].arc == ArcId{6});
}

TEMPLATE_TEST_CASE("H1: an arc realised by a multi-edge physical walk is accepted",
                   "[paper][hierarchy][def41][h1]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.graph.arcs.push_back({0, 3, W{7}});  // exactly d(0,3)
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(rep.ok());
}

TEST_CASE("H1: NaN, negative and infinite arc lengths are rejected", "[paper][hierarchy][h1]") {
  for (const double bad :
       {std::numeric_limits<double>::quiet_NaN(), -1.0, std::numeric_limits<double>::infinity()}) {
    Example<double> ex;
    ex.h.graph.arcs[3].length = bad;
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(rep.has(ViolationKind::kArcLengthInvalid));
    REQUIRE(rep.has(Condition::kH1));
  }
}

TEMPLATE_TEST_CASE("H2: an internal node with a single child is reported",
                   "[paper][hierarchy][def41][h2]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  // Insert U between R and L3: U has one child.
  const auto u = static_cast<NodeId>(ex.h.nodes.size());
  sssp::paper::HierarchyNode<W> node;
  node.parent = ex.r;
  node.children = {3};
  node.height = W{4};
  ex.h.nodes.push_back(node);
  ex.h.nodes[3].parent = u;
  ex.h.nodes[ex.r].children = {ex.a, 2, u};
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH2));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kTooFewChildren);
  REQUIRE(rep.violations[0].node == u);
}

TEMPLATE_TEST_CASE("H2: a zero height is reported", "[paper][hierarchy][def41][h2]", std::uint64_t,
                   double) {
  using W = TestType;
  Example<W> ex;
  ex.h.nodes[ex.a].height = W{0};
  const auto rep = ex.check();
  INFO(rep.summary());
  bool found = false;
  for (const auto& v : rep.violations) {
    if (v.kind == ViolationKind::kNonPositiveHeight && v.node == ex.a) found = true;
  }
  REQUIRE(found);
}

TEST_CASE("H2: a NaN height is reported", "[paper][hierarchy][h2]") {
  Example<double> ex;
  ex.h.nodes[ex.r].height = std::numeric_limits<double>::quiet_NaN();
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(rep.has(ViolationKind::kNonPositiveHeight));
}

TEMPLATE_TEST_CASE("H2: a child higher than its parent is reported at the child",
                   "[paper][hierarchy][def41][h2]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.nodes[ex.a].height = W{5};  // parent R has height 4
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH2));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kHeightDecreases);
  REQUIRE(rep.violations[0].node == ex.a);
}

TEMPLATE_TEST_CASE("H3: a cycle of light arcs between children is reported at the parent",
                   "[paper][hierarchy][def41][h3]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  const auto fw = floyd_warshall<W>(4, ex.edges);
  REQUIRE(fw[2][0] == W{3});                // the new arc is physically feasible (H1 holds)
  ex.h.graph.arcs.push_back({2, 0, W{3}});  // L2 -> A light; with 1->2 (A -> L2) a cycle
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH3));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kLightArcCycle);
  REQUIRE(rep.violations[0].node == ex.r);
}

TEMPLATE_TEST_CASE("H3: arcs of length equal to the height are not light",
                   "[paper][hierarchy][def41][h3]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.graph.arcs.push_back({2, 0, W{4}});  // not < h(R) = 4, so the quotient stays acyclic
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(rep.ok());
}

TEMPLATE_TEST_CASE("H4: a cycle among roots is reported regardless of arc length",
                   "[paper][hierarchy][def41][h4]", std::uint64_t, double) {
  using W = TestType;
  const auto phys = CsrGraph<W>::from_edges(2, {{0, 1, W{1}}});
  StateGraph<W> d;
  d.num_states = 2;
  d.arcs = {{0, 1, W{10}}};
  auto h = sssp::paper::make_leaf_forest<W>(d, {0, 1});
  sssp::paper::collect_roots(h);
  REQUIRE(h.roots.size() == 2);
  {
    const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
    INFO(rep.summary());
    REQUIRE(rep.ok());
  }
  h.graph.arcs.push_back({1, 0, W{10}});
  const auto rep = sssp::paper::validate_ordered_hierarchy(h, phys);
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH4));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kRootArcCycle);
}

TEMPLATE_TEST_CASE("H5: a missing link leaves the children unconnected",
                   "[paper][hierarchy][def41][h5]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.nodes[ex.r].links.pop_back();
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH5));
  REQUIRE(rep.has(ViolationKind::kLinkCount));
  REQUIRE(rep.has(ViolationKind::kLinksNotSpanning));
  for (const auto& v : rep.violations) REQUIRE(v.node == ex.r);
}

TEMPLATE_TEST_CASE("H5: a redundant link closes a cycle in the link tree",
                   "[paper][hierarchy][def41][h5]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.nodes[ex.r].links.push_back({1, 2, W{2}, 2});
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH5));
  REQUIRE(rep.has(ViolationKind::kLinkCount));
  REQUIRE(rep.has(ViolationKind::kLinkCycle));
}

TEMPLATE_TEST_CASE("H5: a link longer than c0 times the height is reported",
                   "[paper][hierarchy][def41][h5]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  ex.h.nodes[ex.r].links[0] = {1, 2, W{5}, 2};  // 5 > 1 * h(R) = 4
  auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH5));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kLinkTooLong);
  REQUIRE(rep.violations[0].node == ex.r);
  REQUIRE(rep.violations[0].link == 0);
  ex.h.c0 = 2;  // 5 <= 2 * 4
  rep = ex.check();
  INFO(rep.summary());
  REQUIRE(rep.ok());
}

TEMPLATE_TEST_CASE("H5: a link shorter than every physical walk is reported",
                   "[paper][hierarchy][def41][h5]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  const auto fw = floyd_warshall<W>(4, ex.edges);
  REQUIRE(fw[1][2] == W{2});
  ex.h.nodes[ex.r].links[0] = {1, 2, W{1}, kNoArc};  // certificate-only, but 1 < d(1,2) = 2
  const auto rep = ex.check();
  INFO(rep.summary());
  REQUIRE(only(rep, Condition::kH5));
  REQUIRE(rep.violations.size() == 1);
  REQUIRE(rep.violations[0].kind == ViolationKind::kLinkWalkTooLong);
  REQUIRE(rep.violations[0].link == 0);
}

TEMPLATE_TEST_CASE("H5: link endpoints must lie in two distinct children of the node",
                   "[paper][hierarchy][def41][h5]", std::uint64_t, double) {
  using W = TestType;
  {
    Example<W> ex;
    ex.h.nodes[ex.r].links[0] = {0, 1, W{1}, 0};  // both endpoints inside child A
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(only(rep, Condition::kH5));
    REQUIRE(rep.has(ViolationKind::kLinkSameChild));
  }
  {
    Example<W> ex;
    ex.h.nodes[ex.a].links[0] = {0, 2, W{1}, kNoArc};  // state 2 is not below A
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(rep.has(ViolationKind::kLinkEndpointOutside));
    for (const auto& v : rep.violations) {
      if (v.condition() == Condition::kH5) REQUIRE(v.node == ex.a);
    }
  }
}

TEMPLATE_TEST_CASE("H6: a link that claims a backing state arc must match that arc",
                   "[paper][hierarchy][def41][h6]", std::uint64_t, double) {
  using W = TestType;
  {
    Example<W> ex;
    ex.h.nodes[ex.r].links[0].arc = 4;  // arc 4 is 2->3, not {1,2}
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(only(rep, Condition::kH6));
    REQUIRE(rep.violations.size() == 1);
    REQUIRE(rep.violations[0].kind == ViolationKind::kLinkArcMismatch);
    REQUIRE(rep.violations[0].node == ex.r);
    REQUIRE(rep.violations[0].link == 0);
  }
  {
    Example<W> ex;
    ex.h.nodes[ex.r].links[0].arc = 99;
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(only(rep, Condition::kH6));
    REQUIRE(rep.has(ViolationKind::kLinkArcOutOfRange));
  }
  {
    Example<W> ex;
    ex.h.nodes[ex.r].links[0].arc = kNoArc;       // certificate-only link is legal
    ex.h.nodes[ex.a].links[0] = {1, 0, W{1}, 0};  // reversed orientation of arc 0 is legal
    const auto rep = ex.check();
    INFO(rep.summary());
    REQUIRE(rep.ok());
  }
}

TEMPLATE_TEST_CASE("malformed forests are reported as structural errors without crashing",
                   "[paper][hierarchy][structure]", std::uint64_t, double) {
  using W = TestType;
  SECTION("projection has the wrong size") {
    Example<W> ex;
    ex.h.projection.pop_back();
    REQUIRE(ex.check().has(Condition::kStructure));
  }
  SECTION("projection points outside the physical graph") {
    Example<W> ex;
    ex.h.projection[3] = 9;
    REQUIRE(ex.check().has(ViolationKind::kProjectionOutOfRange));
  }
  SECTION("arc endpoint out of range") {
    Example<W> ex;
    ex.h.graph.arcs.push_back({0, 7, W{1}});
    const auto rep = ex.check();
    REQUIRE(rep.has(ViolationKind::kArcEndpointOutOfRange));
    REQUIRE(rep.violations[0].arc == ArcId{6});
  }
  SECTION("leaf_of_state disagrees with the leaves") {
    Example<W> ex;
    ex.h.leaf_of_state[1] = 0;
    REQUIRE(ex.check().has(ViolationKind::kLeafMismatch));
  }
  SECTION("parent and child lists disagree") {
    Example<W> ex;
    ex.h.nodes[2].parent = ex.a;
    REQUIRE(ex.check().has(ViolationKind::kParentChildMismatch));
  }
  SECTION("parent pointers form a cycle") {
    Example<W> ex;
    ex.h.nodes[ex.r].parent = ex.a;
    ex.h.nodes[ex.a].children.push_back(ex.r);
    REQUIRE(ex.check().has(Condition::kStructure));
  }
  SECTION("roots list is incomplete") {
    Example<W> ex;
    ex.h.roots.clear();
    REQUIRE(ex.check().has(ViolationKind::kRootsMismatch));
  }
}

TEMPLATE_TEST_CASE("Def 4.2: initial offers must be finite and feasible",
                   "[paper][hierarchy][def42]", std::uint64_t, double) {
  using W = TestType;
  Example<W> ex;
  SeededOrderedHierarchy<W> seeded;
  seeded.hierarchy = ex.h;
  // delta(v) = min(5 + d(0,v), 0 + d(3,v)) = (5, 6, 4, 0)
  seeded.sources = {PhysicalSource<W>{0, W{5}}, PhysicalSource<W>{3, W{0}}};
  seeded.offers = {StateOffer<W>{2, W{4}}, StateOffer<W>{0, W{5}}, StateOffer<W>{1, W{100}}};
  {
    const auto rep = sssp::paper::validate_seeded_ordered_hierarchy(seeded, ex.phys);
    INFO(rep.summary());
    REQUIRE(rep.ok());
  }
  seeded.offers.push_back(StateOffer<W>{2, W{3}});  // 3 < delta(2) = 4
  {
    const auto rep = sssp::paper::validate_seeded_ordered_hierarchy(seeded, ex.phys);
    INFO(rep.summary());
    REQUIRE(only(rep, Condition::kSeed));
    REQUIRE(rep.violations.size() == 1);
    REQUIRE(rep.violations[0].kind == ViolationKind::kOfferInfeasible);
    REQUIRE(rep.violations[0].state == 2);
  }
  seeded.offers.pop_back();
  seeded.offers.push_back(StateOffer<W>{9, W{3}});
  REQUIRE(sssp::paper::validate_seeded_ordered_hierarchy(seeded, ex.phys)
              .has(ViolationKind::kOfferStateOutOfRange));
  seeded.offers.pop_back();
  seeded.offers.push_back(StateOffer<W>{1, sssp::WeightTraits<W>::infinity()});
  REQUIRE(sssp::paper::validate_seeded_ordered_hierarchy(seeded, ex.phys)
              .has(ViolationKind::kOfferNotFinite));
  seeded.offers.pop_back();
  seeded.sources.push_back(PhysicalSource<W>{8, W{0}});
  REQUIRE(sssp::paper::validate_seeded_ordered_hierarchy(seeded, ex.phys)
              .has(ViolationKind::kSourceOutOfRange));
}

TEST_CASE("every violation kind maps to its Def 4.1 condition", "[paper][hierarchy]") {
  using sssp::paper::condition_of;
  REQUIRE(condition_of(ViolationKind::kArcWalkTooLong) == Condition::kH1);
  REQUIRE(condition_of(ViolationKind::kHeightDecreases) == Condition::kH2);
  REQUIRE(condition_of(ViolationKind::kLightArcCycle) == Condition::kH3);
  REQUIRE(condition_of(ViolationKind::kRootArcCycle) == Condition::kH4);
  REQUIRE(condition_of(ViolationKind::kLinkTooLong) == Condition::kH5);
  REQUIRE(condition_of(ViolationKind::kLinkArcMismatch) == Condition::kH6);
  REQUIRE(condition_of(ViolationKind::kOfferInfeasible) == Condition::kSeed);
  REQUIRE(condition_of(ViolationKind::kRootsMismatch) == Condition::kStructure);
}
