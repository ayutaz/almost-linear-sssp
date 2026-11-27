// verify_sssp must accept every correct result and reject each deliberately broken one with
// the specific failure kind. Results here are written by hand, not produced by an algorithm.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "support/floyd_warshall.hpp"
#include "support/graphs.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::FailureKind;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::SsspResult;
using sssp::Tolerance;
using sssp::verify_sssp;
using sssp::Vertex;

using U64 = std::uint64_t;

namespace {

template <class W>
SsspResult<W> make_result(std::vector<W> dist, std::vector<Vertex> parent,
                          std::vector<EdgeId> parent_edge) {
  SsspResult<W> r;
  r.dist = std::move(dist);
  r.parent = std::move(parent);
  r.parent_edge = std::move(parent_edge);
  return r;
}

//   0 --1-- 1 --2-- 2      e0 = {0,1,1}, e1 = {1,2,2}, e2 = {0,2,5}, e3 = {1,3,1}
//   |        \             vertex 4 is isolated
//   +---5----(2)   1 --1-- 3
template <class W>
CsrGraph<W> sample_graph() {
  return CsrGraph<W>::from_edges(5, {{0, 1, W{1}}, {1, 2, W{2}}, {0, 2, W{5}}, {1, 3, W{1}}});
}

template <class W>
SsspResult<W> sample_result() {
  const W inf = sssp::WeightTraits<W>::infinity();
  return make_result<W>({W{0}, W{1}, W{3}, W{2}, inf}, {kNoVertex, 0, 1, 1, kNoVertex},
                        {kNoEdge, 0, 1, 3, kNoEdge});
}

}  // namespace

TEMPLATE_TEST_CASE("verify accepts a correct shortest-path tree", "[verify]", U64, double) {
  const auto rep = verify_sssp(sample_graph<TestType>(), 0, sample_result<TestType>());
  INFO(rep.summary());
  REQUIRE(rep.ok);
  REQUIRE(rep.failures.empty());
  REQUIRE(rep.total_failures == 0);
}

TEMPLATE_TEST_CASE("verify detects broken results", "[verify]", U64, double) {
  using W = TestType;
  const W inf = sssp::WeightTraits<W>::infinity();
  const auto g = sample_graph<W>();
  auto r = sample_result<W>();

  SECTION("distance inconsistent with its tree edge") {
    r.dist[2] = W{4};
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE_FALSE(rep.ok);
    REQUIRE(rep.has(FailureKind::kTreeDistance));
    REQUIRE(rep.has(FailureKind::kNotOptimal));
  }
  SECTION("tree-consistent but not shortest (certificate violated)") {
    r.dist[2] = W{5};
    r.parent[2] = 0;
    r.parent_edge[2] = 2;
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE_FALSE(rep.ok);
    REQUIRE_FALSE(rep.has(FailureKind::kTreeDistance));
    REQUIRE(rep.has(FailureKind::kNotOptimal));
  }
  SECTION("distance smaller than any path") {
    r.dist[3] = W{1};
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kTreeDistance));
  }
  SECTION("parent references form a cycle") {
    r.parent[1] = 3;
    r.parent_edge[1] = 3;
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE_FALSE(rep.ok);
    REQUIRE(rep.has(FailureKind::kParentCycle));
  }
  SECTION("self parent is a cycle") {
    r.parent[3] = 3;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kParentCycle));
  }
  SECTION("parent edge id does not exist") {
    r.parent_edge[2] = 99;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kParentEdgeMissing));
  }
  SECTION("parent edge does not join parent and child") {
    r.parent_edge[2] = 0;  // e0 = {0,1}
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kParentEdgeMismatch));
  }
  SECTION("parent vertex disagrees with parent edge") {
    r.parent[2] = 0;  // e1 = {1,2} is still recorded
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kParentEdgeMismatch));
  }
  SECTION("parent vertex out of range") {
    r.parent[2] = 77;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kInvalidParent));
  }
  SECTION("finite distance without a parent") {
    r.parent[3] = kNoVertex;
    r.parent_edge[3] = kNoEdge;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kMissingParent));
  }
  SECTION("reachable vertex reported as infinite") {
    r.dist[3] = inf;
    r.parent[3] = kNoVertex;
    r.parent_edge[3] = kNoEdge;
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE_FALSE(rep.ok);
    REQUIRE(rep.has(FailureKind::kReachableButInfinite));
  }
  SECTION("unreachable vertex reported as finite") {
    r.dist[4] = W{7};
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE(rep.has(FailureKind::kUnreachableButFinite));
  }
  SECTION("unreachable vertex with a parent") {
    r.parent[4] = 0;
    r.parent_edge[4] = 0;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kUnreachableHasParent));
  }
  SECTION("source distance is not zero") {
    r.dist[0] = W{1};
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kSourceDistance));
  }
  SECTION("source has a parent") {
    r.parent[0] = 1;
    r.parent_edge[0] = 0;
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kSourceParent));
  }
  SECTION("result vectors have the wrong size") {
    r.parent_edge.pop_back();
    const auto rep = verify_sssp(g, 0, r);
    REQUIRE_FALSE(rep.ok);
    REQUIRE(rep.has(FailureKind::kSizeMismatch));
  }
  SECTION("source out of range") {
    REQUIRE(verify_sssp(g, 5, r).has(FailureKind::kInvalidSource));
  }
}

TEMPLATE_TEST_CASE("verify detects a zero-weight parent cycle that every distance check misses",
                   "[verify]", U64, double) {
  using W = TestType;
  // 0 -0- 1 -0- 2, all distances 0; parent[1] = 2 and parent[2] = 1 is locally consistent.
  const auto g = CsrGraph<W>::from_edges(3, {{0, 1, W{0}}, {1, 2, W{0}}});
  const auto r = make_result<W>({W{0}, W{0}, W{0}}, {kNoVertex, 2, 1}, {kNoEdge, 1, 1});
  const auto rep = verify_sssp(g, 0, r);
  REQUIRE_FALSE(rep.ok);
  REQUIRE(rep.has(FailureKind::kParentCycle));
  REQUIRE_FALSE(rep.has(FailureKind::kTreeDistance));
  REQUIRE_FALSE(rep.has(FailureKind::kNotOptimal));
}

TEMPLATE_TEST_CASE("verify does not require a particular tree among equal-length paths", "[verify]",
                   U64, double) {
  using W = TestType;
  // Two zero-weight routes 0-1-3 and 0-2-3.
  const auto g =
      CsrGraph<W>::from_edges(4, {{0, 1, W{0}}, {0, 2, W{0}}, {1, 3, W{0}}, {2, 3, W{0}}});
  const auto via1 =
      make_result<W>({W{0}, W{0}, W{0}, W{0}}, {kNoVertex, 0, 0, 1}, {kNoEdge, 0, 1, 2});
  const auto via2 =
      make_result<W>({W{0}, W{0}, W{0}, W{0}}, {kNoVertex, 0, 0, 2}, {kNoEdge, 0, 1, 3});
  REQUIRE(verify_sssp(g, 0, via1, &via2.dist).ok);
  REQUIRE(verify_sssp(g, 0, via2, &via1.dist).ok);
}

TEMPLATE_TEST_CASE("verify compares against a reference when given", "[verify]", U64, double) {
  using W = TestType;
  const W inf = sssp::WeightTraits<W>::infinity();
  const auto g = sample_graph<W>();
  const auto r = sample_result<W>();
  std::vector<W> ref{W{0}, W{1}, W{3}, W{2}, inf};
  REQUIRE(verify_sssp(g, 0, r, &ref).ok);

  ref[2] = W{4};
  REQUIRE(verify_sssp(g, 0, r, &ref).has(FailureKind::kReferenceMismatch));
  ref[2] = W{3};
  ref[4] = W{9};
  REQUIRE(verify_sssp(g, 0, r, &ref).has(FailureKind::kReferenceMismatch));
  ref.pop_back();
  REQUIRE(verify_sssp(g, 0, r, &ref).has(FailureKind::kReferenceSizeMismatch));
}

TEST_CASE("verify rejects NaN and negative double distances", "[verify][numerics]") {
  const auto g = sample_graph<double>();
  auto r = sample_result<double>();
  r.dist[3] = std::nan("");
  REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kInvalidDistance));
  r.dist[3] = -2.0;
  REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kInvalidDistance));
}

TEST_CASE("Tolerance follows the gamma_k bound for k-hop rounded sums", "[verify][numerics]") {
  constexpr double u = 0x1.0p-53;
  REQUIRE(Tolerance::exact().rel == 0.0);
  REQUIRE(Tolerance::exact().abs == 0.0);
  REQUIRE(Tolerance::for_path_hops(0).rel == 0.0);
  REQUIRE(Tolerance::for_path_hops(1).rel == u / (1 - u));
  REQUIRE(Tolerance::for_path_hops(1000).rel == 1000 * u / (1 - 1000 * u));
  REQUIRE(Tolerance::for_path_hops(1000).abs == 0.0);
  REQUIRE_THROWS_AS(Tolerance::for_path_hops(std::uint64_t{1} << 53), std::invalid_argument);
  const auto g = CsrGraph<double>::from_edges(7, {{0, 1, 1.0}});
  REQUIRE(Tolerance::for_graph(g).rel == Tolerance::for_path_hops(6).rel);
}

TEST_CASE("verify applies the tolerance to double distances only", "[verify][numerics]") {
  // 0 -0.1- 1 -0.2- 2. fl(0.1 + 0.2) = 0.30000000000000004 != 0.3.
  const auto g = CsrGraph<double>::from_edges(3, {{0, 1, 0.1}, {1, 2, 0.2}});
  const double fl_sum = 0.1 + 0.2;
  REQUIRE(fl_sum != 0.3);
  auto r = make_result<double>({0.0, 0.1, fl_sum}, {kNoVertex, 0, 1}, {kNoEdge, 0, 1});
  REQUIRE(verify_sssp(g, 0, r).ok);  // relaxation-style results are exact in fl arithmetic

  const std::vector<double> exact_ref{0.0, 0.1, 0.3};  // nearest doubles to the real distances
  REQUIRE(verify_sssp(g, 0, r, &exact_ref).has(FailureKind::kReferenceMismatch));
  REQUIRE(verify_sssp(g, 0, r, &exact_ref, Tolerance::for_graph(g)).ok);

  SECTION("one ulp off is rejected exactly but accepted within gamma_k") {
    r.dist[2] = std::nextafter(fl_sum, 1.0);
    REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kTreeDistance));
    REQUIRE(verify_sssp(g, 0, r, nullptr, Tolerance::for_graph(g)).ok);
  }
  SECTION("an error far beyond gamma_k is rejected") {
    r.dist[2] = fl_sum * (1 + 1e-9);
    const auto rep = verify_sssp(g, 0, r, &exact_ref, Tolerance::for_graph(g));
    REQUIRE(rep.has(FailureKind::kTreeDistance));
    REQUIRE(rep.has(FailureKind::kReferenceMismatch));
  }
  SECTION("an absolute tolerance widens the bound") {
    r.dist[2] = fl_sum + 1e-6;
    REQUIRE_FALSE(verify_sssp(g, 0, r).ok);
    REQUIRE(verify_sssp(g, 0, r, nullptr, Tolerance{0.0, 1e-5}).ok);
  }
}

TEST_CASE("verify is exact for integer weights regardless of tolerance", "[verify]") {
  const auto g = sample_graph<U64>();
  auto r = sample_result<U64>();
  r.dist[2] = 4;
  REQUIRE_FALSE(verify_sssp(g, 0, r, nullptr, Tolerance{0.5, 10.0}).ok);
}

TEST_CASE("verify reports integer overflow on a tree edge", "[verify]") {
  const U64 max = std::numeric_limits<U64>::max();
  const auto g = CsrGraph<U64>::from_edges(3, {{0, 1, max - 1}, {1, 2, 5}});
  const auto r = make_result<U64>({0, max - 1, 3}, {kNoVertex, 0, 1}, {kNoEdge, 0, 1});
  REQUIRE(verify_sssp(g, 0, r).has(FailureKind::kArithmeticOverflow));
}

TEST_CASE("verify caps recorded failures but counts all of them", "[verify]") {
  constexpr std::size_t kN = 2000;
  const auto g = CsrGraph<U64>::from_edges(kN, sssp::test::uniform_path_edges<U64>(kN, 1));
  SsspResult<U64> r;
  r.dist.assign(kN, 0);  // every vertex but the source violates the tree check
  r.parent.assign(kN, kNoVertex);
  r.parent_edge.assign(kN, kNoEdge);
  for (Vertex v = 1; v < kN; ++v) {
    r.parent[v] = v - 1;
    r.parent_edge[v] = v - 1;
  }
  const auto rep = verify_sssp(g, 0, r);
  REQUIRE_FALSE(rep.ok);
  REQUIRE(rep.failures.size() == sssp::VerificationReport::kMaxRecordedFailures);
  REQUIRE(rep.count(FailureKind::kTreeDistance) == kN - 1);
  REQUIRE(rep.total_failures >= kN - 1);
  REQUIRE(rep.summary().find("tree_distance") != std::string::npos);
}

TEST_CASE("verify accepts Bellman-Ford on the required catalogue", "[verify][bellman_ford]") {
  for (const auto& c : sssp::test::required_graphs<U64>()) {
    CAPTURE(c.name);
    const auto g = c.graph();
    const auto fw = sssp::test::floyd_warshall<U64>(c.n, c.edges);
    const auto rep =
        verify_sssp(g, c.source, sssp::baseline::bellman_ford(g, c.source), &fw[c.source]);
    INFO(rep.summary());
    REQUIRE(rep.ok);
  }
  for (const auto& c : sssp::test::required_graphs<double>()) {
    CAPTURE(c.name);
    const auto g = c.graph();
    const auto fw = sssp::test::floyd_warshall<double>(c.n, c.edges);
    const auto r = sssp::baseline::bellman_ford(g, c.source);
    INFO(verify_sssp(g, c.source, r).summary());
    REQUIRE(verify_sssp(g, c.source, r).ok);  // exact in fl arithmetic
    REQUIRE(verify_sssp(g, c.source, r, &fw[c.source], Tolerance::for_graph(g)).ok);
  }
}

TEST_CASE("write_result_json dumps round-trippable distances", "[verify][dump]") {
  SECTION("double") {
    const auto r = make_result<double>({0.0, 0.1, std::numeric_limits<double>::infinity()},
                                       {kNoVertex, 0, kNoVertex}, {kNoEdge, 0, kNoEdge});
    std::ostringstream os;
    sssp::write_result_json(os, 0, r);
    REQUIRE(os.str() ==
            R"({"n":3,"source":0,"dist":[0,0.10000000000000001,null],"parent":[-1,0,-1],)"
            R"("parent_edge":[-1,0,-1]})"
            "\n");
  }
  SECTION("uint64") {
    const auto r =
        make_result<U64>({0, 18446744073709551614ULL, sssp::WeightTraits<U64>::infinity()},
                         {kNoVertex, 0, kNoVertex}, {kNoEdge, 0, kNoEdge});
    std::ostringstream os;
    sssp::write_result_json(os, 0, r);
    REQUIRE(os.str() ==
            R"({"n":3,"source":0,"dist":[0,18446744073709551614,null],"parent":[-1,0,-1],)"
            R"("parent_edge":[-1,0,-1]})"
            "\n");
  }
}
