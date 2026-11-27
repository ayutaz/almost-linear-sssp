// Property tests for the local radial state representation of Sec 5.2 (Def 5.4, Lemmas
// 5.5-5.7) built from the owner decomposition of Sec 5.1, on seeded graph families with
// brute-force oracles written in the test (naive Dijkstra on the physical and on the state
// graph, Kahn's algorithm on the block quotient).
//
// Separator sizes b in {1, 2, 3, 4, ceil(log2 N), N + 1} are passed directly (Mode B style,
// Q-03): the paper's b = L_N^4 > N on every testable size. K and the macro exponent are the
// paper's (K = least power of two >= 8h, macro iff w <= N^3 sigma). The proofs of Lemmas 5.5-5.7
// use only N >= 2 and K >= 8h (research note 02 Sec 4.2), so the lemmas are asserted as stated.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "sssp/graph.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeId;
using sssp::Vertex;
using sssp::paper::EdgeClass;
using sssp::paper::kNoState;
using sssp::paper::Phase;
using sssp::paper::RadialArcKind;
using sssp::paper::RadialStateGraph;
using sssp::paper::StateId;
namespace lt = sssp::test::local;

namespace {

std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n + 1};
  return {s.begin(), s.end()};
}

// Naive Dijkstra on the directed state graph from state offers.
template <class W>
std::vector<W> state_dijkstra(const RadialStateGraph<W>& rs,
                              const std::vector<std::pair<StateId, W>>& offers) {
  constexpr W kInf = sssp::WeightTraits<W>::infinity();
  const std::size_t s = rs.num_states();
  std::vector<std::vector<std::pair<StateId, W>>> out(s);
  for (const auto& a : rs.graph.arcs) out[a.tail].emplace_back(a.head, a.length);
  std::vector<W> d(s, kInf);
  for (const auto& [x, c] : offers) d[x] = std::min(d[x], c);
  std::vector<bool> done(s, false);
  for (;;) {
    std::size_t best = s;
    for (std::size_t x = 0; x < s; ++x) {
      if (!done[x] && d[x] != kInf && (best == s || d[x] < d[best])) best = x;
    }
    if (best == s) break;
    done[best] = true;
    for (const auto& [y, w] : out[best]) d[y] = std::min(d[y], d[best] + w);
  }
  return d;
}

// Shortest-path tree (parent edge per vertex) of the physical graph from s.
template <class W>
std::vector<EdgeId> sp_tree(std::size_t n, const std::vector<Edge<W>>& edges, Vertex s) {
  constexpr W kInf = sssp::WeightTraits<W>::infinity();
  std::vector<W> d(n, kInf);
  std::vector<EdgeId> pe(n, sssp::kNoEdge);
  std::vector<bool> done(n, false);
  d[s] = W{0};
  for (;;) {
    std::size_t best = n;
    for (std::size_t v = 0; v < n; ++v) {
      if (!done[v] && d[v] != kInf && (best == n || d[v] < d[best])) best = v;
    }
    if (best == n) break;
    done[best] = true;
    for (EdgeId id = 0; id < edges.size(); ++id) {
      const auto& e = edges[id];
      if (e.u != best && e.v != best) continue;
      const Vertex y = e.u == best ? e.v : e.u;
      if (d[best] + e.weight < d[y]) {
        d[y] = d[best] + e.weight;
        pe[y] = id;
      }
    }
  }
  return pe;
}

template <class W>
void check_radial(const std::vector<Edge<W>>& edges, std::size_t n, std::uint64_t seed,
                  const sssp::paper::OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs) {
  const std::size_t m = edges.size();
  const auto& cls = rs.classification;
  const W K = static_cast<W>(std::uint64_t{1} << rs.shell_ratio_log2);
  const W h = static_cast<W>(od.hop_bound);
  const W nn = static_cast<W>(n);
  const auto& r = od.radius;

  // Local classification (l.1021-1026): macro edges are intercell with sigma <= w <= N^3 sigma;
  // corridor edges are macro; radial edges obey w <= 3h max(r_x, r_y).
  for (EdgeId i = 0; i < m; ++i) {
    const auto& e = edges[i];
    if (od.reduced[i] == 0 || e.u == e.v) {
      CHECK(cls.edge_class[i] == EdgeClass::kOmitted);
      continue;
    }
    if (od.cut[i] != 0) CHECK(cls.edge_class[i] == EdgeClass::kMacro);
    if (cls.edge_class[i] == EdgeClass::kMacro) {
      REQUIRE(od.owner[e.u] != od.owner[e.v]);
      const W sg = od.sigma(od.owner[e.u], od.owner[e.v]);
      CHECK(od.corridors[cls.macro_key[i]].weight == sg);
      CHECK_FALSE(e.weight < sg);
      CHECK_FALSE(nn * nn * nn * sg < e.weight);
    } else {
      REQUIRE(cls.edge_class[i] == EdgeClass::kRadial);
      CHECK_FALSE(W{3} * h * std::max(r[e.u], r[e.v]) < e.weight);
      CHECK_FALSE((od.is_pivot(e.u) && od.is_pivot(e.v)));
    }
  }

  // Shells: a K^i <= r < a K^(i+1) within each near component (computed in double: powers of
  // two times quarter-integers below 2^42 are exact).
  for (Vertex v = 0; v < n; ++v) {
    if (od.is_pivot(v)) {
      CHECK(rs.near_component[v] == sssp::paper::kNoNearComponent);
      continue;
    }
    const auto c = rs.near_component[v];
    const double a = static_cast<double>(rs.component_anchor[c]);
    double lo = a;
    for (std::uint32_t i = 0; i < rs.shell[v]; ++i) lo *= static_cast<double>(K);
    CHECK(lo <= static_cast<double>(r[v]));
    CHECK(static_cast<double>(r[v]) < lo * static_cast<double>(K));
    CHECK_FALSE(r[v] < rs.component_anchor[c]);
  }
  for (EdgeId i = 0; i < m; ++i) {
    const auto& e = edges[i];
    if (cls.edge_class[i] != EdgeClass::kRadial || od.is_pivot(e.u) || od.is_pivot(e.v)) continue;
    const W lo = std::min(r[e.u], r[e.v]);
    const W hi = std::max(r[e.u], r[e.v]);
    if (!(K * K * lo < hi)) CHECK(rs.near_component[e.u] == rs.near_component[e.v]);
  }

  // Lemma 5.7 sizes: <= 4 copies per vertex; arcs <= 8 per edge plus 2 switches per vertex.
  CHECK(rs.num_states() <= 4 * n);
  CHECK(rs.graph.arcs.size() <= 8 * m + 2 * n);

  // Every arc projects to a physical walk of its length (Lemma 5.7, converse direction).
  for (std::size_t a = 0; a < rs.graph.arcs.size(); ++a) {
    const auto& arc = rs.graph.arcs[a];
    const Vertex pt = rs.projection[arc.tail];
    const Vertex ph = rs.projection[arc.head];
    if (rs.arc_kind[a] == RadialArcKind::kSwitch) {
      CHECK(pt == ph);
      CHECK(arc.length == W{0});
      CHECK(rs.phase[arc.tail] == Phase::kInward);
      CHECK(arc.head == rs.native_outward[pt]);
    } else {
      const auto& e = edges[rs.arc_edge[a]];
      CHECK(((e.u == pt && e.v == ph) || (e.u == ph && e.v == pt)));
      CHECK(arc.length == e.weight);
    }
  }

  // Blocks lie in one window; within-window arcs stay in a block; the block quotient without
  // macro arcs is a DAG (Lemma 5.7).
  std::vector<StateId> first(rs.num_blocks, kNoState);
  for (StateId x = 0; x < rs.num_states(); ++x) {
    const auto b = rs.block[x];
    if (first[b] == kNoState) {
      first[b] = x;
      continue;
    }
    const StateId y = first[b];
    CHECK(rs.phase[x] == rs.phase[y]);
    CHECK(rs.window[x] == rs.window[y]);
    CHECK(rs.near_component[rs.projection[x]] == rs.near_component[rs.projection[y]]);
    CHECK(rs.projection[x] != rs.projection[y]);
  }
  std::vector<std::vector<std::uint32_t>> qadj(rs.num_blocks);
  std::vector<std::size_t> indeg(rs.num_blocks, 0);
  for (std::size_t a = 0; a < rs.graph.arcs.size(); ++a) {
    const auto bt = rs.block[rs.graph.arcs[a].tail];
    const auto bh = rs.block[rs.graph.arcs[a].head];
    if (rs.arc_kind[a] == RadialArcKind::kWithinWindow) {
      CHECK(bt == bh);
      continue;
    }
    if (rs.arc_kind[a] == RadialArcKind::kMacro) continue;
    CHECK(bt != bh);
    qadj[bt].push_back(bh);
    ++indeg[bh];
  }
  std::vector<std::uint32_t> ready;
  for (std::uint32_t b = 0; b < rs.num_blocks; ++b) {
    if (indeg[b] == 0) ready.push_back(b);
  }
  std::size_t popped = 0;
  while (!ready.empty()) {
    const auto b = ready.back();
    ready.pop_back();
    ++popped;
    for (const auto c : qadj[b]) {
      if (--indeg[c] == 0) ready.push_back(c);
    }
  }
  CHECK(popped == rs.num_blocks);

  // Lemma 5.7: min over copies = physical distance from arbitrary weighted sources at the native
  // inward copies; no copy is below the physical distance.
  sssp::test::SplitMix64 rng(seed ^ 0x5eedULL);
  constexpr W kInf = sssp::WeightTraits<W>::infinity();
  for (int trial = 0; trial < 4; ++trial) {
    std::vector<std::pair<Vertex, W>> src;
    std::vector<std::pair<StateId, W>> offers;
    const std::size_t k = trial == 0 ? 1 : 1 + rng.below(3);
    for (std::size_t i = 0; i < k; ++i) {
      const auto v = static_cast<Vertex>(rng.below(n));
      const W c = trial == 0 ? W{0} : static_cast<W>(rng.below(50));
      src.emplace_back(v, c);
      offers.emplace_back(rs.native_inward[v], c);
    }
    const auto phys = lt::dijkstra_naive<W>(n, edges, src);
    const auto ds = state_dijkstra(rs, offers);
    std::vector<W> best(n, kInf);
    for (StateId x = 0; x < rs.num_states(); ++x) {
      const Vertex v = rs.projection[x];
      CHECK_FALSE(ds[x] < phys[v]);
      best[v] = std::min(best[v], ds[x]);
    }
    CHECK(best == phys);
  }

  // Lemmas 5.5 and 5.6 on every maximal radial run of a shortest path between any two vertices.
  for (Vertex s = 0; s < n; ++s) {
    const auto pe = sp_tree<W>(n, edges, s);
    for (Vertex t = 0; t < n; ++t) {
      std::vector<Vertex> path{t};
      std::vector<EdgeId> pedges;
      for (Vertex x = t; x != s;) {
        const EdgeId id = pe[x];
        REQUIRE(id != sssp::kNoEdge);
        pedges.push_back(id);
        x = edges[id].u == x ? edges[id].v : edges[id].u;
        path.push_back(x);
      }
      for (std::size_t i = 0; i < pedges.size();) {
        REQUIRE(cls.edge_class[pedges[i]] != EdgeClass::kOmitted);
        if (cls.edge_class[pedges[i]] != EdgeClass::kRadial) {
          ++i;
          continue;
        }
        std::size_t j = i;
        while (j < pedges.size() && cls.edge_class[pedges[j]] == EdgeClass::kRadial) ++j;
        // Run: vertices path[i..j], edges pedges[i..j-1].
        for (std::size_t a = i; a < j; ++a) {
          W mx = r[path[a]];
          for (std::size_t bnd = a + 1; bnd <= j; ++bnd) {
            mx = std::max({mx, r[path[bnd]], edges[pedges[bnd - 1]].weight});
            // Lemma 5.5: not both endpoint radii <= M / K.
            CHECK_FALSE((!(mx < K * r[path[a]]) && !(mx < K * r[path[bnd]])));
          }
        }
        std::size_t split = i;
        for (std::size_t a = i; a <= j; ++a) {
          if (r[path[a]] < r[path[split]]) split = a;
        }
        W run_min = r[path[i]];
        for (std::size_t a = i; a <= split; ++a) {
          run_min = std::min(run_min, r[path[a]]);
          CHECK_FALSE(K * run_min < r[path[a]]);  // Lemma 5.6, prefix
        }
        W run_max = r[path[split]];
        for (std::size_t a = split; a <= j; ++a) {
          run_max = std::max(run_max, r[path[a]]);
          CHECK_FALSE(K * r[path[a]] < run_max);  // Lemma 5.6, suffix
        }
        i = j;
      }
    }
  }
}

}  // namespace

TEMPLATE_TEST_CASE("local radial states satisfy Def 5.4 and Lemmas 5.5-5.7 on graph families",
                   "[property]", std::uint64_t, double) {
  using W = TestType;
  const lt::Family families[] = {lt::Family::kPath,   lt::Family::kSpider, lt::Family::kRandomTree,
                                 lt::Family::kSparse, lt::Family::kGrid,   lt::Family::kComplete};
  const lt::Weights dists[] = {lt::Weights::kTies, lt::Weights::kUniform,
                               lt::Weights::kExponential};
  std::uint64_t runs = 0;
  std::uint64_t far_arcs = 0;
  std::uint64_t near_transfers = 0;
  std::uint64_t macro_arcs = 0;
  for (const auto f : families) {
    for (const auto wd : dists) {
      for (std::uint64_t seed = 1; seed <= 5; ++seed) {
        const std::size_t n = f == lt::Family::kComplete ? 3 + seed : 4 + 5 * seed;
        const auto edges = lt::make_graph<W>(f, wd, n, seed * 104729 + static_cast<unsigned>(f));
        const auto g = CsrGraph<W>::from_edges(n, edges);
        for (const std::size_t b : separator_sizes(n)) {
          INFO(lt::name(f) << " / " << lt::name(wd) << " n=" << n << " seed=" << seed
                           << " b=" << b);
          const auto od = sssp::paper::owner_decomposition(g, b);
          const auto rs = sssp::paper::local_radial_states(g, od);
          check_radial(edges, n, seed, od, rs);
          for (const auto k : rs.arc_kind) {
            far_arcs += k == RadialArcKind::kFarTransfer ? 1 : 0;
            near_transfers += k == RadialArcKind::kNearTransfer ? 1 : 0;
            macro_arcs += k == RadialArcKind::kMacro ? 1 : 0;
          }
          ++runs;
        }
      }
    }
  }
  CHECK(runs > 400);
  // The families exercise every transfer kind.
  CHECK(far_arcs > 0);
  CHECK(near_transfers > 0);
  CHECK(macro_arcs > 0);
}
