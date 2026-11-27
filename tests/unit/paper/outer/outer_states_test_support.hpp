#pragma once
// Test support for the outer state graph (Sec 8.1 of arXiv 2611.34911v1, P3-OUTER-STATES).
//
// Builds the inputs of build_outer_states by the steps of Alg 8 l.8-11 (owner decomposition,
// route names, classification and G', OwnerDistances) and provides oracles written here,
// independently of src/paper: exact integer images of lengths (std::uint64_t as is, quarter-
// integer doubles times 4, held in unsigned __int128), an O(n^2) Dijkstra on G' for J, a naive
// Dijkstra on state graphs, Tarjan's SCCs, Kahn's acyclicity test and a union-find.
//
// double inputs are quarter-integers (Q-55 plan A: exactness is established for integer and
// dyadic weights only), so every sum below 2^53 is exact and comparisons are exact.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/local/local_hierarchy.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/model/random.hpp"
#include "paper/owner/aux_metric.hpp"
#include "paper/owner/owner_distances.hpp"
#include "paper/pivots/name_routes.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/local_test_support.hpp"

namespace sssp::test::outer {

using Exact = unsigned __int128;
inline constexpr Exact kExactInf = ~Exact{0};
inline constexpr Exact kScale = 128;  // 2^7, the calendar unit of Alg 7 (Lemma 7.3)

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

/// Exact image of a length, or kExactInf for WeightTraits<W>::infinity().
template <class W>
Exact exact_or_inf(W x) {
  return x == WeightTraits<W>::infinity() ? kExactInf : exact(x);
}

/// Separator sizes used by the Mode B style tests (D-7): the paper's outer b = L_N^64 saturates
/// and gives a single pivot, which is covered by Mode A cases separately.
inline std::vector<std::size_t> separator_sizes(std::size_t n) {
  std::size_t lg = 0;
  while ((std::size_t{1} << lg) < n) ++lg;
  std::set<std::size_t> s{1, 2, 3, 4, std::max<std::size_t>(lg, 1), n, n + 1};
  return {s.begin(), s.end()};
}

/// Exact route names with P = 1 (lambda_e = the least instance weight >= the G_P distance),
/// as in tests/property/paper/test_owner_distances_property.cpp.
template <class W>
paper::pivots::RouteNames<W> exact_names(const paper::pivots::PivotRoutingInstance<W>& inst) {
  const auto dist = floyd_warshall<W>(inst.t, inst.edges);
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
  paper::pivots::RouteNames<W> names;
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

enum class Names : std::uint8_t { kExactP1, kNameRoutesPaper };

/// The inputs of Sec 8.1: G, its outer owner decomposition, the classification of Sec 7.1, G'
/// with T' (Sec 7.2) and Alg 7's J (Sec 7.4).
template <class W>
struct OuterCase {
  std::size_t n = 0;
  std::vector<Edge<W>> edges;
  CsrGraph<W> g;
  paper::OwnerDecomposition<W> od;
  paper::owner::OuterEdgeClassification<W> cls;
  paper::owner::AuxGraph<W> aux;
  paper::owner::OwnerDistancesResult<W> jr;
  Diagnostics diag;  ///< diagnostics of the steps before Sec 8.1
};

template <class W>
OuterCase<W> make_outer_case(std::size_t n, const std::vector<Edge<W>>& edges, std::size_t b,
                             Names names, std::uint64_t seed,
                             const paper::LocalHierarchyParams& local = paper::kPaperLocalParams) {
  OuterCase<W> c;
  c.n = n;
  c.edges = edges;
  c.g = CsrGraph<W>::from_edges(n, edges);
  c.od = paper::owner_decomposition(c.g, b);
  if (names == Names::kExactP1) {
    const auto pe = paper::owner::projected_interowner_edges(c.g, c.od);
    c.cls = paper::owner::classify_outer_edges(c.g, c.od, pe, exact_names(pe.instance));
  } else {
    paper::model::WordRng rng(seed);
    c.cls = paper::owner::classify_interowner_edges(c.g, c.od, rng, {}, &c.diag);
  }
  c.aux = paper::owner::build_aux_graph(c.g, c.od, c.cls);
  paper::owner::OwnerDistancesOptions opt;
  opt.local = local;
  c.jr = paper::owner::owner_distances(c.aux, &c.diag, opt);
  return c;
}

/// Reference J (exact images, input units) per input vertex: O(n^2) Dijkstra on G' from the
/// root, mapped back through aux_vertex (every pivot gets 0).
template <class W>
std::vector<Exact> reference_j(const OuterCase<W>& c) {
  const auto& aux = c.aux;
  const std::size_t na = aux.num_vertices();
  std::vector<Exact> d(na, kExactInf);
  std::vector<std::uint8_t> done(na, 0);
  d[paper::owner::kAuxRoot] = 0;
  for (std::size_t it = 0; it < na; ++it) {
    std::size_t u = na;
    for (std::size_t v = 0; v < na; ++v) {
      if (done[v] == 0 && d[v] != kExactInf && (u == na || d[v] < d[u])) u = v;
    }
    if (u == na) break;
    done[u] = 1;
    for (const auto& arc : aux.graph.arcs(static_cast<Vertex>(u))) {
      const Exact x = d[u] + exact(arc.weight);
      if (x < d[arc.to]) d[arc.to] = x;
    }
  }
  std::vector<Exact> j(c.n);
  for (Vertex v = 0; v < c.n; ++v) j[v] = d[aux.aux_vertex[v]];
  return j;
}

/// Exact physical distances from `s` on the edge list (naive O(n^2) Dijkstra).
template <class W>
std::vector<Exact> physical_distances(std::size_t n, const std::vector<Edge<W>>& edges, Vertex s) {
  std::vector<std::vector<std::pair<Vertex, Exact>>> adj(n);
  for (const auto& e : edges) {
    adj[e.u].emplace_back(e.v, exact(e.weight));
    adj[e.v].emplace_back(e.u, exact(e.weight));
  }
  std::vector<Exact> d(n, kExactInf);
  std::vector<std::uint8_t> done(n, 0);
  d[s] = 0;
  for (std::size_t it = 0; it < n; ++it) {
    std::size_t u = n;
    for (std::size_t v = 0; v < n; ++v) {
      if (done[v] == 0 && d[v] != kExactInf && (u == n || d[v] < d[u])) u = v;
    }
    if (u == n) break;
    done[u] = 1;
    for (const auto& [y, w] : adj[u]) d[y] = std::min(d[y], d[u] + w);
  }
  return d;
}

/// A directed arc with an exact length (test-side copy of a state arc).
struct XArc {
  std::uint32_t tail;
  std::uint32_t head;
  Exact len;
};

/// Exact distances from `s` in a directed graph on `ns` states (naive O(ns^2 + arcs)).
inline std::vector<Exact> state_distances(std::size_t ns, const std::vector<XArc>& arcs,
                                          std::uint32_t s) {
  std::vector<std::vector<std::pair<std::uint32_t, Exact>>> adj(ns);
  for (const auto& a : arcs) adj[a.tail].emplace_back(a.head, a.len);
  std::vector<Exact> d(ns, kExactInf);
  std::vector<std::uint8_t> done(ns, 0);
  d[s] = 0;
  for (std::size_t it = 0; it < ns; ++it) {
    std::size_t u = ns;
    for (std::size_t v = 0; v < ns; ++v) {
      if (done[v] == 0 && d[v] != kExactInf && (u == ns || d[v] < d[u])) u = v;
    }
    if (u == ns) break;
    done[u] = 1;
    for (const auto& [y, w] : adj[u]) d[y] = std::min(d[y], d[u] + w);
  }
  return d;
}

/// Strongly connected components (Tarjan, iterative); returns a component id per vertex.
inline std::vector<std::uint32_t> scc_ids(
    std::size_t n, const std::vector<std::pair<std::uint32_t, std::uint32_t>>& arcs) {
  std::vector<std::vector<std::uint32_t>> adj(n);
  for (const auto& [t, h] : arcs) adj[t].push_back(h);
  constexpr std::uint32_t kUnset = 0xffffffffu;
  std::vector<std::uint32_t> index(n, kUnset), low(n, 0), comp(n, kUnset);
  std::vector<std::uint8_t> on_stack(n, 0);
  std::vector<std::uint32_t> stack;
  std::uint32_t next_index = 0;
  std::uint32_t next_comp = 0;
  std::vector<std::pair<std::uint32_t, std::size_t>> call;
  for (std::uint32_t root = 0; root < n; ++root) {
    if (index[root] != kUnset) continue;
    call.emplace_back(root, 0);
    index[root] = low[root] = next_index++;
    stack.push_back(root);
    on_stack[root] = 1;
    while (!call.empty()) {
      auto& [v, i] = call.back();
      if (i < adj[v].size()) {
        const std::uint32_t w = adj[v][i++];
        if (index[w] == kUnset) {
          index[w] = low[w] = next_index++;
          stack.push_back(w);
          on_stack[w] = 1;
          call.emplace_back(w, 0);
        } else if (on_stack[w] != 0) {
          low[v] = std::min(low[v], index[w]);
        }
        continue;
      }
      const std::uint32_t done = v;
      call.pop_back();
      if (!call.empty()) low[call.back().first] = std::min(low[call.back().first], low[done]);
      if (low[done] == index[done]) {
        for (;;) {
          const std::uint32_t x = stack.back();
          stack.pop_back();
          on_stack[x] = 0;
          comp[x] = next_comp;
          if (x == done) break;
        }
        ++next_comp;
      }
    }
  }
  return comp;
}

/// Kahn's algorithm: true iff the digraph on n vertices is acyclic (self-loops are cycles).
inline bool is_acyclic(std::size_t n,
                       const std::vector<std::pair<std::uint32_t, std::uint32_t>>& arcs) {
  std::vector<std::size_t> indeg(n, 0);
  std::vector<std::vector<std::uint32_t>> adj(n);
  for (const auto& [t, h] : arcs) {
    adj[t].push_back(h);
    ++indeg[h];
  }
  std::vector<std::uint32_t> queue;
  for (std::uint32_t v = 0; v < n; ++v) {
    if (indeg[v] == 0) queue.push_back(v);
  }
  for (std::size_t i = 0; i < queue.size(); ++i) {
    for (const auto h : adj[queue[i]]) {
      if (--indeg[h] == 0) queue.push_back(h);
    }
  }
  return queue.size() == n;
}

/// Union-find with path halving.
struct UnionFind {
  std::vector<std::uint32_t> parent;
  explicit UnionFind(std::size_t n) : parent(n) {
    std::iota(parent.begin(), parent.end(), std::uint32_t{0});
  }
  std::uint32_t find(std::uint32_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  }
  void unite(std::uint32_t a, std::uint32_t b) { parent[find(a)] = find(b); }
};

/// True iff the two labelings induce the same partition of 0..n-1.
inline bool same_partition(const std::vector<std::uint32_t>& a,
                           const std::vector<std::uint32_t>& b) {
  if (a.size() != b.size()) return false;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
  std::vector<std::uint32_t> map_ab(a.size() + 1, 0xffffffffu), map_ba(b.size() + 1, 0xffffffffu);
  // Relabel to dense ids first so that the maps fit.
  const auto dense = [](const std::vector<std::uint32_t>& x) {
    std::vector<std::uint32_t> ids(x.begin(), x.end());
    std::ranges::sort(ids);
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    std::vector<std::uint32_t> out(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
      out[i] = static_cast<std::uint32_t>(std::ranges::lower_bound(ids, x[i]) - ids.begin());
    }
    return out;
  };
  const auto da = dense(a);
  const auto db = dense(b);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (map_ab[da[i]] == 0xffffffffu) map_ab[da[i]] = db[i];
    if (map_ba[db[i]] == 0xffffffffu) map_ba[db[i]] = da[i];
    if (map_ab[da[i]] != db[i] || map_ba[db[i]] != da[i]) return false;
  }
  return true;
}

/// Light clusters below heavy edges to vertex 0 (as in the Alg 7 tests), with chords of
/// weights spread over 2^0..2^12: J ratios far above K^2 = 256, several shells, far edges.
template <class W>
std::vector<Edge<W>> heavy_spokes(std::size_t n, std::uint64_t seed) {
  SplitMix64 rng(seed);
  std::vector<Edge<W>> e;
  if (n <= 1) return e;
  const std::size_t clusters = 1 + rng.below(3);
  std::vector<std::vector<Vertex>> members(clusters);
  for (Vertex v = 1; v < n; ++v) members[rng.below(clusters)].push_back(v);
  for (const auto& mem : members) {
    if (mem.empty()) continue;
    e.push_back({0, mem[0], local::to_weight<W>(1000ULL << rng.between(0, 10), rng)});
    for (std::size_t i = 1; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(i)], mem[i], local::to_weight<W>(rng.between(1, 3), rng)});
    }
    for (std::size_t i = 0; i < mem.size(); ++i) {
      e.push_back({mem[rng.below(mem.size())], mem[rng.below(mem.size())],
                   local::to_weight<W>(rng.between(1, 3) << rng.between(0, 12), rng)});
    }
  }
  return e;
}

/// Vertex roles of three_pivot_corridor.
struct ThreePivotCorridor {
  std::size_t n = 0;
  std::size_t b = 0;  ///< separator size that marks exactly p2, q, p (Sec 5.1 postorder)
  Vertex root = 0;    ///< vertex 0, a leaf of p (T is rooted at vertex 0, Q-31)
  Vertex p = 0, q = 0, p2 = 0, y = 0, t = 0;
  std::vector<Vertex> path;  ///< a_1 .. a_k, a_1 adjacent to q, a_k adjacent to y
  EdgeId corridor_yp = 0;    ///< the cut corridor y - p (w = 256, macro)
  EdgeId corridor_qp2 = 0;   ///< the cut edge q - p2 between adjacent pivots (w = 200, macro)
  EdgeId shortcut_yq = 0;    ///< y - q (w = 260, nontree, radial, kept at equality 260 = J_y)
};

/// A hand-built input on which the owner-tree states of Sec 8.1 (l.1984-1997) are necessary
/// for Lemma 8.3. T (rooted at the leaf 0 of p):
///   0 - p (1000), leaves of p (1000), p - y (256),
///   y - a_k - ... - a_1 - q (255 each), q - t (17), leaves of q (1000), q - p2 (200),
///   leaves of p2 (1000);
/// nontree edges t - a_1 (256) and y - q (260). With b = k + 4 the postorder marks p2, q and p
/// (residuals b each) and nothing else, so A = {p, q, p2}; the corridor q..p is cut at y - p.
/// J (Lemma 7.2): J_t = 17, J_{a_i} = min(255 i, 260 + 255 (k + 1 - i)), J_y = 260. The near
/// component {t, a_1, .., a_k, y} has anchor 17 and shells [17, 272), [272, 4352), [4352, ...):
/// y and a_1 are in shell 0, the middle of the path (k >= 34) in shell 2. At the threshold
/// W = 256 the shortcut y - q is excluded, an ordinary inward walk from y cannot climb from
/// shell 0 to shell 2, and an ordinary outward walk from q cannot descend from shell 2 to y, so
/// only v_in / v_out connect y with q: the macro cycles through p (via y - p) and through p2
/// (via q - p2) form one SCC with the special states and two without them.
template <class W>
std::vector<Edge<W>> three_pivot_corridor(std::size_t k, ThreePivotCorridor& r) {
  REQUIRE(k >= 34);
  const std::size_t f3 = 2;      // leaves of p
  const std::size_t f1 = k + 2;  // leaves of q besides t
  const std::size_t f4 = k + 3;  // leaves of p2
  r.b = k + 4;
  std::vector<Edge<W>> e;
  Vertex next = 0;
  r.root = next++;
  r.p = next++;
  e.push_back({r.root, r.p, W{1000}});
  for (std::size_t i = 0; i < f3; ++i) e.push_back({r.p, next++, W{1000}});
  r.y = next++;
  r.corridor_yp = static_cast<EdgeId>(e.size());
  e.push_back({r.p, r.y, W{256}});
  r.path.assign(k, 0);
  Vertex prev = r.y;
  for (std::size_t i = k; i >= 1; --i) {  // a_k first, adjacent to y
    r.path[i - 1] = next++;
    e.push_back({prev, r.path[i - 1], W{255}});
    prev = r.path[i - 1];
  }
  r.q = next++;
  e.push_back({r.path[0], r.q, W{255}});
  r.t = next++;
  e.push_back({r.q, r.t, W{17}});
  for (std::size_t i = 0; i < f1; ++i) e.push_back({r.q, next++, W{1000}});
  r.p2 = next++;
  r.corridor_qp2 = static_cast<EdgeId>(e.size());
  e.push_back({r.q, r.p2, W{200}});
  for (std::size_t i = 0; i < f4; ++i) e.push_back({r.p2, next++, W{1000}});
  e.push_back({r.t, r.path[0], W{256}});
  r.shortcut_yq = static_cast<EdgeId>(e.size());
  e.push_back({r.y, r.q, W{260}});
  r.n = next;
  return e;
}

/// Calls f(OuterCase&) for every generated case: the local families x weight distributions x
/// sizes x seeds plus heavy spokes, for every separator size of separator_sizes(n) (Mode B
/// style, D-7) and the paper's outer b = L_N^64 (Mode A), each with exact P = 1 names and with
/// NameRoutes (paper parameters). The case is described by INFO.
template <class W, class F>
void for_each_case(const std::vector<std::size_t>& sizes, std::uint64_t seeds, F&& f) {
  const local::Family families[] = {local::Family::kPath,       local::Family::kSpider,
                                    local::Family::kRandomTree, local::Family::kSparse,
                                    local::Family::kGrid,       local::Family::kComplete};
  const local::Weights weights[] = {local::Weights::kTies, local::Weights::kUniform,
                                    local::Weights::kExponential};
  const auto run = [&](std::size_t n, const std::vector<Edge<W>>& edges, std::uint64_t gseed,
                       const std::string& what) {
    auto bs = separator_sizes(n);
    bs.push_back(paper::owner::outer_separator_size(n));
    for (const std::size_t b : bs) {
      for (const Names names : {Names::kExactP1, Names::kNameRoutesPaper}) {
        INFO(what << " n=" << n << " seed=" << gseed << " b=" << b
                  << (names == Names::kExactP1 ? " exact-P1" : " name-routes"));
        auto c = make_outer_case<W>(n, edges, b, names, gseed + b);
        f(c);
      }
    }
  };
  for (const std::size_t n : sizes) {
    for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
      for (const auto fam : families) {
        for (const auto wd : weights) {
          const std::uint64_t gseed = seed * 1000003 + n * 31 +
                                      static_cast<std::uint64_t>(fam) * 7 +
                                      static_cast<std::uint64_t>(wd);
          run(n, local::make_graph<W>(fam, wd, n, gseed), gseed,
              std::string(local::name(fam)) + " " + local::name(wd));
        }
      }
      const std::uint64_t hseed = seed * 7919 + n;
      run(n, heavy_spokes<W>(n, hseed), hseed, "heavy-spokes");
    }
  }
}

}  // namespace sssp::test::outer
