// Property tests for the Sec 5.4 assembly (txt l.1182-1243, Lemma 5.10 l.1265-1300) on seeded
// connected graph families, with an oracle written independently of src/paper:
//   * Def 4.1 by validate_ordered_hierarchy (c0 = 2);
//   * the threshold-SCC characterisation used in the proof of Lemma 5.10 (l.1271-1275): the
//     states below every internal node C form one strongly connected component of the state
//     arcs shorter than its class upper bound (at most that bound for an exact merge), and the
//     states below each child of C form an SCC of the arcs shorter than h(C); the roots are the
//     SCCs of the whole state graph. SCCs come from a Kosaraju search written here;
//   * block scales (l.1183-1189, Q-34, Q-36): s_U <= r_x < K^2 s_U, intrastate edges at most
//     3hK^2 s_U, t_U in [s_U/2, s_U], t_U = s_U for separate-grid blocks, aligned t_U a grid
//     point; unaligned blocks have total MST length < sigma_U / 2 (l.1208);
//   * Y components (l.1189-1197): one owner, one radius r_Y >= s_U, fewer than b vertices, and
//     every local height below t_U;
//   * O(N) nodes: at most S - 1 internal nodes, children - 1 links each.
//
// Separator sizes are passed directly (Q-03); N >= 32 and h <= N as Sec 5.4 l.1206 requires.
// A clustered family (light clusters joined by edges of weight about 2^54) produces unaligned
// blocks, which need sigma_U > N^10 s_U.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/saturating.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::BlockScaleKind;
using sssp::paper::NodeId;
using sssp::paper::NodeOrigin;
using sssp::paper::RadialArcKind;
using sssp::paper::StateId;
namespace lt = sssp::test::local;

namespace {

/// Kosaraju SCC ids of the state arcs accepted by `keep`.
template <class W, class Keep>
std::vector<std::size_t> scc(const sssp::paper::StateGraph<W>& d, Keep keep) {
  const std::size_t s = d.num_states;
  std::vector<std::vector<StateId>> fwd(s), rev(s);
  for (const auto& a : d.arcs) {
    if (!keep(a.length)) continue;
    fwd[a.tail].push_back(a.head);
    rev[a.head].push_back(a.tail);
  }
  std::vector<std::size_t> order;
  std::vector<bool> seen(s, false);
  for (std::size_t r = 0; r < s; ++r) {
    if (seen[r]) continue;
    std::vector<std::pair<StateId, std::size_t>> st{{static_cast<StateId>(r), 0}};
    seen[r] = true;
    while (!st.empty()) {
      auto& [x, i] = st.back();
      if (i < fwd[x].size()) {
        const StateId y = fwd[x][i++];
        if (!seen[y]) {
          seen[y] = true;
          st.emplace_back(y, 0);
        }
      } else {
        order.push_back(x);
        st.pop_back();
      }
    }
  }
  std::vector<std::size_t> comp(s, s);
  std::size_t nc = 0;
  for (std::size_t i = order.size(); i-- > 0;) {
    const std::size_t r = order[i];
    if (comp[r] != s) continue;
    std::vector<StateId> st{static_cast<StateId>(r)};
    comp[r] = nc;
    while (!st.empty()) {
      const StateId x = st.back();
      st.pop_back();
      for (const StateId y : rev[x]) {
        if (comp[y] == s) {
          comp[y] = nc;
          st.push_back(y);
        }
      }
    }
    ++nc;
  }
  return comp;
}

template <class W>
std::vector<StateId> states_below(const sssp::paper::OrderedHierarchy<W>& h, NodeId v) {
  std::vector<StateId> out;
  std::vector<NodeId> st{v};
  while (!st.empty()) {
    const NodeId x = st.back();
    st.pop_back();
    if (h.is_leaf(x)) out.push_back(h.nodes[x].state);
    for (const NodeId c : h.nodes[x].children) st.push_back(c);
  }
  std::sort(out.begin(), out.end());
  return out;
}

/// The SCC of `x` as a sorted state list.
std::vector<StateId> members(const std::vector<std::size_t>& comp, StateId x) {
  std::vector<StateId> out;
  for (std::size_t y = 0; y < comp.size(); ++y) {
    if (comp[y] == comp[x]) out.push_back(static_cast<StateId>(y));
  }
  return out;
}

struct Coverage {
  std::size_t aligned = 0, unaligned = 0, unattached = 0, y_nontrivial = 0, exact_nodes = 0,
              separate_nodes = 0, global_nodes = 0, multi_root = 0;
};

template <class W>
void check_one(const CsrGraph<W>& g, std::size_t b, Coverage& cov) {
  const auto od = sssp::paper::owner_decomposition(g, b);
  const auto rs = sssp::paper::local_radial_states(g, od);
  const auto j = sssp::paper::compute_joins(od, rs);
  const auto a = sssp::paper::assemble_local_hierarchy(od, rs, j);
  const auto& h = a.result.hierarchy;
  const std::size_t s = rs.num_states();
  const std::size_t n = g.num_vertices();

  // Def 4.1 with c0 = 2 (Lemma 5.10).
  REQUIRE(h.c0 == 2);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(report.summary());
  REQUIRE(report.ok());

  // Output maps (Alg 4 l.12).
  REQUIRE(h.num_states() == s);
  REQUIRE(a.result.native_entry.size() == n);
  REQUIRE(a.result.arc_edge.size() == h.graph.arcs.size());
  for (Vertex v = 0; v < n; ++v) {
    REQUIRE(a.result.native_entry[v] == rs.native_inward[v]);
    REQUIRE(h.projection[a.result.native_entry[v]] == v);
  }

  // O(N) nodes; children - 1 links per internal node.
  std::size_t internal = 0;
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.is_leaf(static_cast<NodeId>(v))) continue;
    ++internal;
    REQUIRE(h.nodes[v].links.size() + 1 == h.nodes[v].children.size());
  }
  REQUIRE(internal + 1 <= std::max<std::size_t>(s, 1));
  REQUIRE(a.node_upper.size() == h.nodes.size());
  REQUIRE(a.node_origin.size() == h.nodes.size());

  // Threshold SCCs (proof of Lemma 5.10 l.1271-1275). Cache Kosaraju per threshold.
  std::map<std::pair<W, int>, std::vector<std::size_t>> cache;
  const auto comp_of = [&](W bound, int mode) -> const std::vector<std::size_t>& {
    auto it = cache.find({bound, mode});
    if (it != cache.end()) return it->second;
    auto c = mode == 0 ? scc(h.graph, [&](W l) { return l < bound; })
                       : scc(h.graph, [&](W l) { return l <= bound; });
    return cache.emplace(std::make_pair(bound, mode), std::move(c)).first->second;
  };
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    const auto node = static_cast<NodeId>(v);
    if (h.is_leaf(node)) continue;
    switch (a.node_origin[v]) {
      case NodeOrigin::kExactLocal:
        ++cov.exact_nodes;
        REQUIRE(a.node_upper[v] == h.nodes[v].height);
        REQUIRE(a.node_upper_inclusive[v] == 1);
        break;
      case NodeOrigin::kSeparateGrid:
        ++cov.separate_nodes;
        break;
      case NodeOrigin::kGlobalGrid:
        ++cov.global_nodes;
        break;
    }
    const auto below = states_below(h, node);
    const auto& up = comp_of(a.node_upper[v], a.node_upper_inclusive[v]);
    INFO("node " << v << " height " << h.nodes[v].height << " upper " << a.node_upper[v]);
    REQUIRE(members(up, below.front()) == below);
    const auto& low = comp_of(h.nodes[v].height, 0);
    for (const NodeId c : h.nodes[v].children) {
      const auto cb = states_below(h, c);
      REQUIRE(members(low, cb.front()) == cb);
    }
  }
  const auto all = scc(h.graph, [](W) { return true; });
  for (const NodeId r : h.roots) {
    const auto below = states_below(h, r);
    REQUIRE(members(all, below.front()) == below);
  }
  if (h.roots.size() > 1) ++cov.multi_root;

  // Block scales.
  std::vector<std::vector<StateId>> block_states(rs.num_blocks);
  for (StateId x = 0; x < s; ++x) block_states[rs.block[x]].push_back(x);
  REQUIRE(a.blocks.size() == rs.num_blocks);
  REQUIRE(a.y_component.size() == s);
  for (std::size_t u = 0; u < rs.num_blocks; ++u) {
    const auto& bs = a.blocks[u];
    const auto& xs = block_states[u];
    REQUIRE(bs.size == xs.size());
    if (xs.size() == 1) {
      REQUIRE(bs.kind == BlockScaleKind::kSingleton);
      continue;
    }
    REQUIRE(bs.kind != BlockScaleKind::kSingleton);
    INFO("block " << u << " s " << bs.s << " t " << bs.t);
    for (const StateId x : xs) {
      const W r = od.radius[rs.projection[x]];
      REQUIRE(bs.s <= r);  // Q-34 (local)
      REQUIRE(r < sssp::paper::sat::mul_pow2(bs.s, 2 * rs.shell_ratio_log2));
    }
    // Every intrastate edge has weight at most 3hK^2 s_U (l.1186).
    const std::uint64_t factor = (3 * std::uint64_t{od.hop_bound}) << (2 * rs.shell_ratio_log2);
    const W bound = sssp::paper::sat::mul_int(bs.s, factor);
    for (std::size_t arc = 0; arc < rs.graph.arcs.size(); ++arc) {
      if (rs.arc_kind[arc] != RadialArcKind::kWithinWindow) continue;
      if (rs.block[rs.graph.arcs[arc].tail] != u) continue;
      REQUIRE(rs.graph.arcs[arc].length <= bound);
    }
    // Q-36 (local): t_U in [s_U / 2, s_U].
    REQUIRE(bs.t <= bs.s);
    REQUIRE(bs.s <= bs.t + bs.t);
    W mst_total{0};
    for (const auto& jr : j.joins) {
      if (jr.kind == sssp::paper::JoinKind::kBlockMst && rs.block[jr.a] == u)
        mst_total += jr.weight;
    }
    switch (bs.kind) {
      case BlockScaleKind::kAligned: {
        ++cov.aligned;
        REQUIRE(bs.z != sssp::paper::kNoState);
        REQUIRE(std::find(a.grid.begin(), a.grid.end(), bs.t) != a.grid.end());
        break;
      }
      case BlockScaleKind::kUnaligned: {
        ++cov.unaligned;
        REQUIRE(bs.t == bs.s);
        const W sigma = od.corridors[bs.sigma_key].weight;
        REQUIRE(mst_total + mst_total < sigma);  // l.1208
        break;
      }
      case BlockScaleKind::kUnattached:
        ++cov.unattached;
        REQUIRE(bs.t == bs.s);
        REQUIRE(bs.z == sssp::paper::kNoState);
        break;
      case BlockScaleKind::kSingleton:
        break;
    }
    if (bs.z != sssp::paper::kNoState) {
      REQUIRE(rs.block[bs.z] == u);
      for (const StateId x : xs) REQUIRE(!(j.theta[x] < j.theta[bs.z]));
      REQUIRE(j.b0_key[bs.z] == bs.sigma_key);
    }
  }

  // Y components (l.1189-1197).
  std::map<std::uint32_t, std::vector<StateId>> ys;
  for (StateId x = 0; x < s; ++x) {
    if (a.y_component[x] != sssp::paper::kNoComponent) ys[a.y_component[x]].push_back(x);
  }
  for (const auto& [y, xs] : ys) {
    const auto u = rs.block[xs.front()];
    const auto& bs = a.blocks[u];
    std::set<std::uint32_t> owners;
    std::set<W> radii;
    for (const StateId x : xs) {
      REQUIRE(rs.block[x] == u);
      owners.insert(od.owner[rs.projection[x]]);
      radii.insert(od.radius[rs.projection[x]]);
    }
    if (xs.size() >= 2) {
      ++cov.y_nontrivial;
      REQUIRE(owners.size() == 1);
      REQUIRE(radii.size() == 1);
      REQUIRE(*radii.begin() >= bs.s);
      REQUIRE(xs.size() < b);
    }
  }
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.is_leaf(static_cast<NodeId>(v)) || a.node_origin[v] != NodeOrigin::kExactLocal) continue;
    const auto x = states_below(h, static_cast<NodeId>(v)).front();
    REQUIRE(h.nodes[v].height < a.blocks[rs.block[x]].t);
  }
}

template <class W>
void sweep(Coverage& cov, bool with_clusters) {
  for (const auto f : {lt::Family::kPath, lt::Family::kSpider, lt::Family::kRandomTree,
                       lt::Family::kSparse, lt::Family::kGrid}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t n : {32U, 45U}) {
        for (const std::size_t b : {1U, 2U, 3U, 6U, 31U}) {
          for (std::uint64_t seed = 1; seed <= 2; ++seed) {
            const auto g = CsrGraph<W>::from_edges(n, lt::make_graph<W>(f, wd, n, seed * 97 + b));
            INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " b=" << b
                             << " seed=" << seed);
            check_one(g, b, cov);
          }
        }
      }
    }
  }
  if (with_clusters) {
    if constexpr (std::is_same_v<W, std::uint64_t>) {
      for (const std::size_t n : {32U, 40U, 64U}) {
        for (const std::size_t b : {2U, 3U, 5U, 8U}) {
          for (std::uint64_t seed = 1; seed <= 6; ++seed) {
            const auto g = CsrGraph<W>::from_edges(n, lt::clustered(n, seed * 31 + b));
            INFO("clustered n=" << n << " b=" << b << " seed=" << seed);
            check_one(g, b, cov);
          }
        }
      }
    }
  }
}

}  // namespace

TEMPLATE_TEST_CASE("Sec 5.4 assembly: Def 4.1, threshold SCCs, block scales and Y components", "",
                   std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  sweep<W>(cov, true);
  // Every branch of Sec 5.4 is exercised.
  CHECK(cov.aligned > 0);
  CHECK(cov.unattached > 0);
  CHECK(cov.y_nontrivial > 0);
  CHECK(cov.exact_nodes > 0);
  CHECK(cov.separate_nodes > 0);
  CHECK(cov.global_nodes > 0);
  CHECK(cov.multi_root > 0);
  if constexpr (std::is_same_v<W, std::uint64_t>) CHECK(cov.unaligned > 0);
}
