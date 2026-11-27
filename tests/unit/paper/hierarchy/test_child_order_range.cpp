// Range tests for the child order of Sec 4.2 (l.725-779, Lemma 4.4, Lemma 4.5): lengths close to
// the top of the weight type (task P3-CHILDORDER-OVERFLOW).
//
// What is specified (overflow rule of include/sssp/paper.hpp): the root distances of the link
// forest F and the running charge of the traversal are auxiliary quantities. They must not
// raise; a gap g(i, j) that fits the type is returned exactly, and a gap beyond the range is
// WeightTraits<W>::infinity() (an upper bound that is still valid).
//
// Expectations:
//   * hand examples: traced by hand through Sec 4.2 at unit scale (the trace of
//     test_child_order.cpp), then multiplied by the scale;
//   * homogeneity: Sec 4.1 and Sec 4.2 only compare and add lengths, so multiplying every length
//     by c > 0 keeps every decision (ranks, representatives) and multiplies every gap by c. The
//     scaled gaps are compared with c * g formed in unsigned __int128;
//   * hierarchies without light arcs (Alg 4 base case): every gap is the F distance between
//     consecutive representatives (l.747), recomputed here by a direct walk over the links in
//     unsigned __int128.
//
// Scaling after normalize(): on a hierarchy whose mass S(C) does not fit, normalize() of this
// branch still raises (the saturation D-23 belongs to P3-EVAL-ARITH-FIX). order_children() does
// not read the mass, so the tests normalise at unit scale, multiply every length of the result
// and set the masses to their saturated products.
//
// double: scaling by a power of two is exact, so the same expectations hold bit for bit as long
// as the unit-scale sums are exact (quarter-integer lengths). Root distances beyond the range of
// double are +infinity and make the gaps that use them +infinity.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"

using sssp::Vertex;
using sssp::paper::NodeId;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::OrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;

namespace {

__extension__ typedef unsigned __int128 U128;

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

/// a * c in the weight type, saturated at infinity.
template <class W>
W scaled(W a, W c) {
  if constexpr (std::is_integral_v<W>) {
    const U128 p = U128{a} * c;
    return p >= U128{kInf<W>} ? kInf<W> : static_cast<W>(p);
  } else {
    return a * c;  // c is a power of two: exact, or +inf
  }
}

/// Multiplies every length of a normalised hierarchy by c (see the file comment).
template <class W>
void scale_lengths(NormalizedHierarchy<W>& nh, W c) {
  for (auto& a : nh.forest.graph.arcs) a.length = scaled(a.length, c);
  for (auto& nd : nh.forest.nodes) {
    nd.height = scaled(nd.height, c);
    for (auto& lk : nd.links) lk.length = scaled(lk.length, c);
  }
  for (auto& v : nh.unit) v = scaled(v, c);
  for (auto& v : nh.original_height) v = scaled(v, c);
  for (auto& v : nh.mass) v = scaled(v, c);
}

/// The largest length order_children() reads.
template <class W>
W max_length(const NormalizedHierarchy<W>& nh) {
  W m{0};
  for (const auto& a : nh.forest.graph.arcs) m = a.length > m ? a.length : m;
  for (const auto& nd : nh.forest.nodes) {
    m = nd.height > m ? nd.height : m;
    for (const auto& lk : nd.links) m = lk.length > m ? lk.length : m;
  }
  return m;
}

/// States 0..5 on the link path 0 - 1 - 2 - 3 - 4 - 5; the hand example of
/// test_child_order.cpp ("cyclic start inside a child"), whose trace gives the rank order
/// Z, A, 5, 4 at the root with g(Z, A) = 4, g(A, 5) = 6, g(5, 4) = 10 and, inside A, 3, 2, 1
/// with gaps 1 and 1.
template <class W>
struct HandExample {
  NormalizedHierarchy<W> nh;
  NodeId x, a, z, l4, l5, l1, l2, l3;
};

template <class W>
HandExample<W> hand_example() {
  StateGraph<W> g;
  g.num_states = 6;
  g.arcs = {{0, 2, W{4}}, {3, 4, W{3}}, {1, 5, W{5}}};
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2, 3, 4, 5});
  const NodeId a =
      sssp::paper::add_internal_node<W>(h, W{1}, {1, 2, 3}, {{1, 2, W{1}}, {2, 3, W{1}}});
  sssp::paper::add_internal_node<W>(h, W{8}, {0, a, 4, 5},
                                    {{0, 1, W{2}}, {3, 4, W{2}}, {4, 5, W{2}}});
  sssp::paper::collect_roots(h);
  HandExample<W> e{sssp::paper::normalize(h), 0, 0, 0, 0, 0, 0, 0, 0};
  e.x = e.nh.forest.roots[0];
  e.a = e.nh.node_of_input[a];
  e.z = e.nh.node_of_input[0];
  e.l1 = e.nh.node_of_input[1];
  e.l2 = e.nh.node_of_input[2];
  e.l3 = e.nh.node_of_input[3];
  e.l4 = e.nh.node_of_input[4];
  e.l5 = e.nh.node_of_input[5];
  return e;
}

template <class W>
void check_hand_example(const HandExample<W>& e, W s) {
  const auto& nh = e.nh;
  REQUIRE(nh.forest.nodes[e.x].children == std::vector<NodeId>{e.z, e.a, e.l5, e.l4});
  CHECK(nh.rep[e.z] == StateId{0});
  CHECK(nh.rep[e.a] == StateId{2});
  CHECK(nh.rep[e.l5] == StateId{5});
  CHECK(nh.rep[e.l4] == StateId{4});
  CHECK(nh.gap_next[e.z] == scaled(W{4}, s));
  CHECK(nh.gap_next[e.a] == scaled(W{6}, s));
  CHECK(nh.gap_next[e.l5] == scaled(W{10}, s));
  CHECK(nh.gap_next[e.l4] == W{0});
  REQUIRE(nh.forest.nodes[e.a].children == std::vector<NodeId>{e.l3, e.l2, e.l1});
  CHECK(nh.gap_next[e.l3] == s);
  CHECK(nh.gap_next[e.l2] == s);
  CHECK(nh.gap_next[e.l1] == W{0});
}

/// One root of height 2 over the leaves 0, 1, 2 with the links (0, 1) and (0, 2) of length 2 and
/// no arcs. F is the star centred at state 0 with preorder 0, 1, 2; the three leaves are DFS
/// roots in that order and finish at once, so the rank order is 2, 1, 0 with
/// g(2, 1) = F(1, 2) = 4 and g(1, 0) = F(0, 1) = 2 (l.747).
template <class W>
NormalizedHierarchy<W> star_example() {
  StateGraph<W> g;
  g.num_states = 3;
  auto h = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
  sssp::paper::add_internal_node<W>(h, W{2}, {0, 1, 2}, {{0, 1, W{2}}, {0, 2, W{2}}});
  sssp::paper::collect_roots(h);
  return sssp::paper::normalize(h);
}

template <class W>
W length_from(std::mt19937_64& rng) {
  if constexpr (std::is_integral_v<W>) {
    return static_cast<W>(1 + rng() % 12);
  } else {
    return 0.25 * static_cast<double>(1 + rng() % 48);  // quarter-integers: exact sums
  }
}

/// Threshold-SCC reference hierarchy of a random directed state graph (has light arcs).
template <class W>
OrderedHierarchy<W> reference_instance(std::uint64_t seed, std::size_t n, std::size_t m) {
  std::mt19937_64 rng(seed);
  StateGraph<W> g;
  g.num_states = n;
  for (std::size_t i = 0; i < m; ++i) {
    const auto u = static_cast<StateId>(rng() % n);
    const auto v = static_cast<StateId>(rng() % n);
    g.arcs.push_back({u, v, length_from<W>(rng)});
  }
  std::vector<Vertex> pi(n);
  for (auto& p : pi) p = static_cast<Vertex>(rng() % n);
  return sssp::paper::build_threshold_scc_hierarchy<W>(std::move(g), std::move(pi));
}

/// Alg 4 base case of a random connected graph (no light arcs, Lemma 3.5).
template <class W>
OrderedHierarchy<W> base_case_instance(std::uint64_t seed) {
  sssp::gen::GenOptions opt;
  opt.seed = seed;
  opt.connected = true;
  opt.weights = {sssp::gen::WeightDist::kUniformInt, 1.0, 48.0};
  auto el = (seed % 2) == 0 ? sssp::gen::erdos_renyi<W>(10 + (seed % 25), 3.0, opt)
                            : sssp::gen::random_tree<W>(5 + (seed % 30), opt);
  if constexpr (!std::is_integral_v<W>) {
    for (auto& e : el.edges) e.weight *= 0.25;  // quarter-integers
  }
  const auto g = sssp::to_csr(std::move(el));
  return sssp::paper::local_hierarchy_base_case(g).hierarchy;
}

/// The scale that brings the largest length of `nh` to the top of the type: uint64, the largest
/// factor that keeps every single length below the sentinel; double, 2^900 (exact, and the sums
/// of these small instances stay finite).
template <class W>
W top_scale(const NormalizedHierarchy<W>& nh) {
  if constexpr (std::is_integral_v<W>) {
    return (kInf<W> - 1) / max_length(nh);
  } else {
    return std::ldexp(1.0, 900);
  }
}

/// F distances from every state by a direct search over the links, in unsigned __int128 for
/// integers (never overflows) and in double for double (exact for the scaled quarter-integers).
template <class W>
using Exact = std::conditional_t<std::is_integral_v<W>, U128, double>;

template <class W>
std::vector<std::vector<Exact<W>>> forest_distances(const OrderedHierarchy<W>& f) {
  const std::size_t n = f.num_states();
  std::vector<std::vector<std::pair<StateId, W>>> adj(n);
  for (const auto& nd : f.nodes) {
    for (const auto& lk : nd.links) {
      adj[lk.x].push_back({lk.y, lk.length});
      adj[lk.y].push_back({lk.x, lk.length});
    }
  }
  std::vector<std::vector<Exact<W>>> d(n, std::vector<Exact<W>>(n, Exact<W>{0}));
  for (StateId s = 0; s < n; ++s) {
    std::vector<std::uint8_t> seen(n, 0);
    std::vector<StateId> stack{s};
    seen[s] = 1;
    while (!stack.empty()) {
      const StateId u = stack.back();
      stack.pop_back();
      for (const auto& [v, w] : adj[u]) {
        if (seen[v] != 0) continue;
        seen[v] = 1;
        d[s][v] = d[s][u] + Exact<W>{w};
        stack.push_back(v);
      }
    }
  }
  return d;
}

template <class W>
W to_weight(Exact<W> v) {
  if constexpr (std::is_integral_v<W>) {
    return v >= U128{kInf<W>} ? kInf<W> : static_cast<W>(v);
  } else {
    return v;
  }
}

/// Orders `nh` at unit scale and at scale c and checks homogeneity. Returns the scaled result.
template <class W>
NormalizedHierarchy<W> check_homogeneous(const NormalizedHierarchy<W>& nh, W c) {
  auto small = nh;
  sssp::paper::order_children(small);
  auto big = nh;
  scale_lengths(big, c);
  REQUIRE_NOTHROW(sssp::paper::order_children(big));
  REQUIRE(big.rank == small.rank);
  REQUIRE(big.rep == small.rep);
  for (NodeId v = 0; v < nh.num_nodes(); ++v) {
    REQUIRE(big.forest.nodes[v].children == small.forest.nodes[v].children);
    INFO("node " << v);
    REQUIRE(big.gap_next[v] == scaled(small.gap_next[v], c));
  }
  return big;
}

}  // namespace

TEST_CASE("child order range: uint64 hand example at scale 5 * 2^58 keeps exact gaps",
          "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  // nu(X) = 8 s and every gap (at most 10 s) fit, but the segment between the marks of leaf 4
  // and leaf 5 (return over the arc 3 -> 4, move F(3, 1), arc 1 -> 5) is 14 s - 4 s when it is
  // kept as a pair of sums of root distances, and 14 s = 17.5 * 2^60 >= 2^64.
  constexpr W s = W{5} << 58;
  auto e = hand_example<W>();
  scale_lengths(e.nh, s);
  REQUIRE(e.nh.unit[e.x] == W{5} << 61);
  REQUIRE_NOTHROW(sssp::paper::order_children(e.nh));
  check_hand_example(e, s);
}

TEST_CASE("child order range: uint64 root distances beyond 2^64 are exact",
          "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  // Links 2, 1, 1, 2, 2 on the path 0 .. 5 at scale 2^62: rd = 0, 2, 3, 4, 6, 8 times 2^62.
  constexpr W s = W{1} << 62;
  auto e = hand_example<W>();
  scale_lengths(e.nh, s);
  sssp::paper::LinkForest<W> lf;
  REQUIRE_NOTHROW(lf = sssp::paper::build_link_forest(e.nh));
  const std::vector<unsigned> expect{0, 2, 3, 4, 6, 8};
  for (StateId x = 0; x < 6; ++x) {
    INFO("state " << x);
    CHECK(lf.root_distance[x] == U128{expect[x]} * s);
    CHECK(lf.preorder[x] == x);
  }
  // 4 s = 2^64 and 10 s do not fit: those gaps saturate; the others are exact.
  REQUIRE_NOTHROW(sssp::paper::order_children(e.nh));
  REQUIRE(e.nh.forest.nodes[e.x].children == std::vector<NodeId>{e.z, e.a, e.l5, e.l4});
  CHECK(e.nh.gap_next[e.z] == kInf<W>);   // 4 * 2^62
  CHECK(e.nh.gap_next[e.a] == kInf<W>);   // 6 * 2^62
  CHECK(e.nh.gap_next[e.l5] == kInf<W>);  // 10 * 2^62
  CHECK(e.nh.gap_next[e.l3] == s);
  CHECK(e.nh.gap_next[e.l2] == s);
}

TEST_CASE("child order range: a uint64 gap that does not fit saturates at infinity",
          "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  auto nh = star_example<W>();
  const NodeId x = nh.forest.roots[0];
  const NodeId l0 = nh.node_of_input[0];
  const NodeId l1 = nh.node_of_input[1];
  const NodeId l2 = nh.node_of_input[2];
  SECTION("links of 2^63: g(2, 1) = 2^64 saturates, g(1, 0) = 2^63 is exact") {
    scale_lengths(nh, W{1} << 62);
    REQUIRE_NOTHROW(sssp::paper::order_children(nh));
    REQUIRE(nh.forest.nodes[x].children == std::vector<NodeId>{l2, l1, l0});
    CHECK(nh.gap_next[l2] == kInf<W>);
    CHECK(nh.gap_next[l1] == W{1} << 63);
    CHECK(nh.gap_next[l0] == W{0});
  }
  SECTION("links of 2^63 - 1: g(2, 1) = 2^64 - 2 is the largest finite value and is exact") {
    for (auto& lk : nh.forest.nodes[x].links) lk.length = (W{1} << 63) - 1;
    nh.forest.nodes[x].height = nh.unit[x] = nh.original_height[x] = (W{1} << 63) - 1;
    nh.mass[x] = kInf<W>;
    REQUIRE_NOTHROW(sssp::paper::order_children(nh));
    CHECK(nh.gap_next[l2] == kInf<W> - 1);
    CHECK(nh.gap_next[l1] == (W{1} << 63) - 1);
  }
}

TEST_CASE("child order range: uint64 gaps are homogeneous up to the top of the type",
          "[paper][hierarchy][child_order]") {
  using W = std::uint64_t;
  std::size_t saturated = 0;
  std::size_t finite = 0;
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    INFO("seed " << seed);
    const auto nh =
        sssp::paper::normalize(reference_instance<W>(seed, 8 + seed % 30, 20 + 3 * seed));
    const auto big = check_homogeneous(nh, top_scale(nh));
    for (const W g : big.gap_next) (g == kInf<W> ? saturated : finite) += 1;
  }
  // Both outcomes occur: the family reaches past the range and stays within it.
  CHECK(saturated > 0);
  CHECK(finite > 0);
}

TEMPLATE_TEST_CASE("child order range: base-case gaps equal exact F distances at the top scale",
                   "[paper][hierarchy][child_order]", std::uint64_t, double) {
  using W = TestType;
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    INFO("seed " << seed);
    const auto nh = sssp::paper::normalize(base_case_instance<W>(seed));
    REQUIRE(nh.light_arcs.empty());
    const auto big = check_homogeneous(nh, top_scale(nh));
    const auto fdist = forest_distances(big.forest);
    const auto lf = sssp::paper::build_link_forest(big);
    for (StateId s = 0; s < big.forest.num_states(); ++s) {
      // F is rooted at the smallest state of each tree; the hierarchy is connected.
      REQUIRE(lf.root_distance[s] == fdist[0][s]);
    }
    for (NodeId x = 0; x < big.num_nodes(); ++x) {
      const auto& ch = big.forest.nodes[x].children;
      for (std::size_t k = 0; k + 1 < ch.size(); ++k) {
        INFO("node " << x << " child " << k);
        REQUIRE(big.gap_next[ch[k]] == to_weight<W>(fdist[big.rep[ch[k]]][big.rep[ch[k + 1]]]));
      }
    }
  }
}

TEST_CASE("child order range: double gaps of quarter-integer lengths scale exactly",
          "[paper][hierarchy][child_order]") {
  using W = double;
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    INFO("seed " << seed);
    const auto nh =
        sssp::paper::normalize(reference_instance<W>(seed, 8 + seed % 30, 20 + 3 * seed));
    const auto big = check_homogeneous(nh, top_scale(nh));
    for (const W g : big.gap_next) REQUIRE(std::isfinite(g));
  }
}

TEST_CASE("child order range: double hand example at scale 5 * 2^1018 keeps exact gaps",
          "[paper][hierarchy][child_order]") {
  using W = double;
  // Every root distance (at most 8 s) and every gap (at most 10 s) is finite, but 14 s on the
  // positive side of the pair (see the uint64 case) is 17.5 * 2^1020 >= 2^1024.
  const W s = std::ldexp(5.0, 1018);
  auto e = hand_example<W>();
  scale_lengths(e.nh, s);
  REQUIRE_NOTHROW(sssp::paper::order_children(e.nh));
  check_hand_example(e, s);
}

TEST_CASE("child order range: a double gap beyond the range is +infinity, never NaN",
          "[paper][hierarchy][child_order]") {
  using W = double;
  SECTION("star with links of 2^1023: g(2, 1) = 2^1024") {
    auto nh = star_example<W>();
    const NodeId x = nh.forest.roots[0];
    scale_lengths(nh, std::ldexp(1.0, 1022));
    REQUIRE_NOTHROW(sssp::paper::order_children(nh));
    const auto& ch = nh.forest.nodes[x].children;
    REQUIRE(ch.size() == 3);
    CHECK(nh.gap_next[ch[0]] == kInf<W>);
    CHECK(nh.gap_next[ch[1]] == std::ldexp(1.0, 1023));
    CHECK(nh.gap_next[ch[2]] == 0.0);
  }
  SECTION("root distances beyond the range: finite gaps below them stay exact") {
    // Path 0 .. 5 at scale 2^1021: rd = 0, 2, 3, 4, 6, 8 times 2^1021, so rd(5) = 2^1024 is
    // +infinity. Inside A (states 1, 2, 3) every root distance is finite and the gaps are
    // exact. At the root, every segment moves to or from state 5 or sums past the range, and
    // its gap is +infinity (documented limit of the double representation: g(Z, A) = 4 s would
    // fit).
    const W s = std::ldexp(1.0, 1021);
    auto e = hand_example<W>();
    scale_lengths(e.nh, s);
    sssp::paper::LinkForest<W> lf;
    REQUIRE_NOTHROW(lf = sssp::paper::build_link_forest(e.nh));
    CHECK(lf.root_distance[4] == 6.0 * s);
    CHECK(lf.root_distance[5] == kInf<W>);
    REQUIRE_NOTHROW(sssp::paper::order_children(e.nh));
    CHECK(e.nh.gap_next[e.l3] == s);
    CHECK(e.nh.gap_next[e.l2] == s);
    for (const NodeId v : {e.z, e.a, e.l5}) {
      INFO("node " << v);
      CHECK_FALSE(std::isnan(e.nh.gap_next[v]));
      CHECK(e.nh.gap_next[v] >= 4.0 * s);
    }
  }
}
