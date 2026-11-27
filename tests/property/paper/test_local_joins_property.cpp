// Property tests for the threshold components and joins of Sec 5.3 on seeded graph families:
// I0, O0, I*, O* against brute force (search / bottleneck Dijkstra per definition, l.1110-1117),
// the bounds of Lemma 5.8, and Lemma 5.9 exhaustively: for every distinct W, the union-find
// components of the joins of weight <= W (and < W) restricted to states equal the strongly
// connected components (Kosaraju, written here) of the state arcs of length <= W (< W).
//
// Separator sizes b in {1, 2, 3, 4, ceil(log2 N), N + 1} are passed directly (Mode B style,
// Q-03); K and the macro test are the paper's. Lemmas 5.8-5.9 use only N >= 2 and K >= 8h
// (research note 02 Sec 4.2).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::JoinKind;
using sssp::paper::RadialArcKind;
using sssp::paper::RadialStateGraph;
using sssp::paper::StateId;
using sssp::paper::ThresholdJoins;
namespace lt = sssp::test::local;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n + 1};
  return {s.begin(), s.end()};
}

// Bottleneck (minimax) search over the eligible arcs (all but macro), forward or backward.
template <class W>
std::vector<W> minimax(const RadialStateGraph<W>& rs, std::vector<W> d, bool reverse) {
  const std::size_t s = rs.num_states();
  std::vector<bool> done(s, false);
  for (;;) {
    std::size_t best = s;
    for (std::size_t x = 0; x < s; ++x) {
      if (!done[x] && d[x] != kInf<W> && (best == s || d[x] < d[best])) best = x;
    }
    if (best == s) break;
    done[best] = true;
    for (std::size_t a = 0; a < rs.graph.arcs.size(); ++a) {
      if (rs.arc_kind[a] == RadialArcKind::kMacro) continue;
      const auto& arc = rs.graph.arcs[a];
      const StateId from = reverse ? arc.head : arc.tail;
      const StateId to = reverse ? arc.tail : arc.head;
      if (from != best) continue;
      d[to] = std::min(d[to], std::max(d[best], arc.length));
    }
  }
  return d;
}

// Reachability over eligible arcs from one state (forward or backward).
template <class W>
std::vector<bool> reach(const RadialStateGraph<W>& rs, StateId src, bool reverse) {
  std::vector<bool> seen(rs.num_states(), false);
  std::vector<StateId> st{src};
  seen[src] = true;
  while (!st.empty()) {
    const StateId x = st.back();
    st.pop_back();
    for (std::size_t a = 0; a < rs.graph.arcs.size(); ++a) {
      if (rs.arc_kind[a] == RadialArcKind::kMacro) continue;
      const auto& arc = rs.graph.arcs[a];
      const StateId from = reverse ? arc.head : arc.tail;
      const StateId to = reverse ? arc.tail : arc.head;
      if (from == x && !seen[to]) {
        seen[to] = true;
        st.push_back(to);
      }
    }
  }
  return seen;
}

// SCC labels (Kosaraju) of the arcs accepted by `keep`.
template <class W, class Keep>
std::vector<int> scc(const RadialStateGraph<W>& rs, Keep keep) {
  const std::size_t s = rs.num_states();
  std::vector<std::vector<StateId>> fw(s), bw(s);
  for (const auto& a : rs.graph.arcs) {
    if (!keep(a.length)) continue;
    fw[a.tail].push_back(a.head);
    bw[a.head].push_back(a.tail);
  }
  std::vector<bool> seen(s, false);
  std::vector<StateId> finish;
  std::function<void(StateId)> dfs1 = [&](StateId x) {
    seen[x] = true;
    for (const StateId y : fw[x]) {
      if (!seen[y]) dfs1(y);
    }
    finish.push_back(x);
  };
  for (StateId x = 0; x < s; ++x) {
    if (!seen[x]) dfs1(x);
  }
  std::vector<int> label(s, -1);
  int next = 0;
  std::function<void(StateId)> dfs2 = [&](StateId x) {
    label[x] = next;
    for (const StateId y : bw[x]) {
      if (label[y] < 0) dfs2(y);
    }
  };
  for (std::size_t i = s; i-- > 0;) {
    if (label[finish[i]] < 0) {
      dfs2(finish[i]);
      ++next;
    }
  }
  return label;
}

// Union-find labels of the join vertices over the joins accepted by `keep`.
template <class W, class Keep>
std::vector<std::size_t> join_components(const ThresholdJoins<W>& j, Keep keep) {
  std::vector<std::size_t> p(j.num_join_vertices());
  std::iota(p.begin(), p.end(), std::size_t{0});
  const std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
    return p[x] == x ? x : p[x] = find(p[x]);
  };
  for (const auto& r : j.joins) {
    if (keep(r.weight)) p[find(r.a)] = find(r.b);
  }
  std::vector<std::size_t> out(j.num_states);
  for (std::size_t x = 0; x < j.num_states; ++x) out[x] = find(x);
  return out;
}

// Two labellings of the states describe the same partition.
template <class A, class B>
bool same_partition(const std::vector<A>& a, const std::vector<B>& b) {
  std::map<A, B> ab;
  std::map<B, A> ba;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (ab.emplace(a[i], b[i]).first->second != b[i]) return false;
    if (ba.emplace(b[i], a[i]).first->second != a[i]) return false;
  }
  return true;
}

template <class W>
void check_joins(std::size_t n, const sssp::paper::OwnerDecomposition<W>& od,
                 const RadialStateGraph<W>& rs, const ThresholdJoins<W>& j,
                 std::uint64_t& nontrivial_cores) {
  const std::size_t s = rs.num_states();
  const auto& arcs = rs.graph.arcs;
  REQUIRE(j.num_states == s);
  REQUIRE(j.num_pivots == od.pivots.size());

  // Brute force from the definitions.
  std::vector<W> i0(s, kInf<W>), o0(s, kInf<W>), seed_in(s, kInf<W>), seed_out(s, kInf<W>);
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (rs.arc_kind[a] != RadialArcKind::kMacro) continue;
    const W w = arcs[a].length;
    const auto fwd = reach(rs, arcs[a].head, false);
    const auto bwd = reach(rs, arcs[a].tail, true);
    for (StateId z = 0; z < s; ++z) {
      if (fwd[z]) i0[z] = std::min(i0[z], w);
      if (bwd[z]) o0[z] = std::min(o0[z], w);
    }
    seed_in[arcs[a].head] = std::min(seed_in[arcs[a].head], w);
    seed_out[arcs[a].tail] = std::min(seed_out[arcs[a].tail], w);
  }
  const auto istar = minimax(rs, seed_in, false);
  const auto ostar = minimax(rs, seed_out, true);
  CHECK(j.i0 == i0);
  CHECK(j.o0 == o0);
  CHECK(j.i_star == istar);
  CHECK(j.o_star == ostar);
  for (StateId z = 0; z < s; ++z) {
    CHECK(j.theta[z] == std::max(istar[z], ostar[z]));
    CHECK(j.b0[z] == std::max(i0[z], o0[z]));
    // Attaining macro edges carry the value; B0's named key is its corridor maximum.
    if (j.i0[z] != kInf<W>) CHECK(od.cut.size() > j.i0_edge[z]);
    if (j.b0[z] == kInf<W>) {
      CHECK(j.b0_edge[z] == sssp::kNoEdge);
      continue;
    }
    REQUIRE(j.b0_edge[z] != sssp::kNoEdge);
    CHECK(rs.classification.edge_class[j.b0_edge[z]] == sssp::paper::EdgeClass::kMacro);
    CHECK(j.b0_key[z] == rs.classification.macro_key[j.b0_edge[z]]);
  }

  // Lemma 5.8 (exact in double: see local_test_support; products here are far from tight).
  const double h = static_cast<double>(od.hop_bound);
  const double kk = static_cast<double>(std::uint64_t{1} << rs.shell_ratio_log2);
  const double n3 = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n);
  for (StateId z = 0; z < s; ++z) {
    CHECK((j.b0[z] == kInf<W>) == (j.theta[z] == kInf<W>));
    if (j.b0[z] == kInf<W>) continue;
    const double b0 = static_cast<double>(j.b0[z]);
    const double th = static_cast<double>(j.theta[z]);
    const double sg = static_cast<double>(od.corridors[j.b0_key[z]].weight);
    CHECK_FALSE(j.theta[z] < od.radius[rs.projection[z]]);
    CHECK(b0 <= th);
    CHECK(th <= 3.0 * h * kk * kk * b0);
    CHECK(sg <= th);
    CHECK(th <= 3.0 * h * kk * kk * n3 * sg);
  }

  // Join records are well formed.
  std::size_t owners = 0;
  for (const auto& r : j.joins) {
    REQUIRE(r.a < j.num_join_vertices());
    REQUIRE(r.b < j.num_join_vertices());
    if (r.kind == JoinKind::kOwner) {
      ++owners;
      CHECK(r.weight == j.theta[r.a]);
      CHECK(r.b == j.aux_vertex(od.owner[rs.projection[r.a]]));
    }
  }
  std::size_t finite = 0;
  for (StateId z = 0; z < s; ++z) finite += j.theta[z] != kInf<W> ? 1 : 0;
  CHECK(owners == finite);

  // Lemma 5.9 for every distinct positive W among the arc lengths and join weights, both the
  // non-strict and the strict version.
  std::set<W> values;
  for (const auto& a : arcs) {
    if (W{0} < a.length) values.insert(a.length);
  }
  for (const auto& r : j.joins) values.insert(r.weight);
  for (const W w : values) {
    const auto le = [w](W x) { return !(w < x); };
    const auto lt_ = [w](W x) { return x < w; };
    const auto s_le = scc(rs, le);
    const auto s_lt = scc(rs, lt_);
    INFO("W = " << w);
    CHECK(same_partition(s_le, join_components(j, le)));
    CHECK(same_partition(s_lt, join_components(j, lt_)));
    std::map<int, int> sizes;
    for (const int c : s_le) ++sizes[c];
    for (const auto& [c, k] : sizes) nontrivial_cores += k > 1 ? 1 : 0;
  }
}

}  // namespace

TEMPLATE_TEST_CASE("threshold joins match brute force and satisfy Lemmas 5.8-5.9 on families",
                   "[property]", std::uint64_t, double) {
  using W = TestType;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights dists[] = {lt::Weights::kTies, lt::Weights::kUniform,
                               lt::Weights::kExponential};
  std::uint64_t runs = 0;
  std::uint64_t cores = 0;
  std::uint64_t finite_positive_theta = 0;
  for (const auto f : families) {
    for (const auto wd : dists) {
      for (std::uint64_t seed = 1; seed <= 4; ++seed) {
        const std::size_t n = f == lt::Family::kComplete ? 3 + seed : 4 + 5 * seed;
        const auto edges = lt::make_graph<W>(f, wd, n, seed * 15485863 + static_cast<unsigned>(f));
        const auto g = CsrGraph<W>::from_edges(n, edges);
        for (const std::size_t b : separator_sizes(n)) {
          INFO(lt::name(f) << " / " << lt::name(wd) << " n=" << n << " seed=" << seed
                           << " b=" << b);
          const auto od = sssp::paper::owner_decomposition(g, b);
          const auto rs = sssp::paper::local_radial_states(g, od);
          const auto j = sssp::paper::compute_joins(od, rs);
          check_joins(n, od, rs, j, cores);
          for (StateId z = 0; z < rs.num_states(); ++z) {
            if (j.theta[z] != kInf<W> && !od.is_pivot(rs.projection[z])) ++finite_positive_theta;
          }
          ++runs;
        }
      }
    }
  }
  CHECK(runs > 300);
  CHECK(cores > 0);
  CHECK(finite_positive_theta > 0);
}
