// Property tests for the edge-list utilities (P2-GRAPH-EDGELIST): every function is compared on
// random small multigraphs against a brute-force reference written from the contract in
// include/sssp/edge_list.hpp (quadratic scans and a boolean transitive closure; no shared code),
// and shortest-path distances are checked to be preserved with the Floyd-Warshall oracle.
//
// Exactness for double: the "small" weights are multiples of 1/4 below 4 and m <= 24, so every
// path sum, and every sum of two such sums formed by Floyd-Warshall, is a multiple of 1/4 below
// 2^8 (at most 10 significant bits). The sums are exact whatever the vertex order, so distances
// are compared with ==.
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/edge_list.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "support/floyd_warshall.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::EdgeId;
using sssp::EdgeList;
using sssp::kNoVertex;
using sssp::LoopPolicy;
using sssp::Vertex;
using sssp::test::SplitMix64;

namespace {

constexpr int kSeeds = 400;
constexpr std::size_t kMaxN = 9;
constexpr std::size_t kMaxM = 24;

enum class Weights { kSmall, kHuge };

template <class W>
W draw_weight(SplitMix64& rng, Weights kind) {
  if constexpr (std::is_integral_v<W>) {
    if (kind == Weights::kHuge) {
      return std::numeric_limits<W>::max() - rng.below(4);  // ties and near-max comparisons
    }
    return rng.below(4);  // 0..3: many ties, zero weights
  } else {
    if (kind == Weights::kHuge) {
      constexpr double kBig = 0x1.fffffffffffffp+1023;
      return rng.below(2) == 0 ? kBig : std::nextafter(kBig, 0.0);
    }
    const auto q = rng.below(16);  // 0, 1/4, ..., 15/4; a zero is sometimes -0.0
    if (q == 0) return rng.below(2) == 0 ? 0.0 : -0.0;
    return static_cast<double>(q) * 0.25;
  }
}

template <class W>
EdgeList<W> random_multigraph(SplitMix64& rng, Weights kind) {
  EdgeList<W> g;
  g.num_vertices = rng.below(kMaxN + 1);
  if (g.num_vertices == 0) return g;
  const auto m = rng.below(kMaxM + 1);
  for (std::uint64_t i = 0; i < m; ++i) {
    Vertex u = static_cast<Vertex>(rng.below(g.num_vertices));
    Vertex v = static_cast<Vertex>(rng.below(g.num_vertices));
    // Bias towards repeated pairs (either orientation) and loops so the dedupe has work to do.
    if (!g.edges.empty() && rng.below(3) == 0) {
      const Edge<W> prev = g.edges[rng.below(g.edges.size())];
      u = rng.below(2) == 0 ? prev.u : prev.v;
      v = rng.below(2) == 0 ? prev.u : prev.v;
    }
    g.edges.push_back({u, v, draw_weight<W>(rng, kind)});
  }
  return g;
}

template <class W>
bool same_edge(const Edge<W>& a, const Edge<W>& b) {
  if (a.u != b.u || a.v != b.v || !(a.weight == b.weight)) return false;
  if constexpr (std::is_floating_point_v<W>) {
    if (std::signbit(a.weight) != std::signbit(b.weight)) return false;
  }
  return true;
}

bool same_pair(Vertex a, Vertex b, Vertex c, Vertex d) {
  return (a == c && b == d) || (a == d && b == c);
}

/// Reference dedupe: edge i survives iff no edge j of the same unordered pair precedes it in
/// the (weight, index) order; loops are dropped under kDrop.
template <class W>
std::vector<EdgeId> brute_lightest(const EdgeList<W>& g, LoopPolicy loops) {
  std::vector<EdgeId> kept;
  for (std::size_t i = 0; i < g.edges.size(); ++i) {
    const auto& ei = g.edges[i];
    if (loops == LoopPolicy::kDrop && ei.u == ei.v) continue;
    bool beaten = false;
    for (std::size_t j = 0; j < g.edges.size() && !beaten; ++j) {
      const auto& ej = g.edges[j];
      if (j == i || !same_pair(ei.u, ei.v, ej.u, ej.v)) continue;
      beaten = ej.weight < ei.weight || (!(ei.weight < ej.weight) && j < i);
    }
    if (!beaten) kept.push_back(static_cast<EdgeId>(i));
  }
  return kept;
}

/// Reference reachability: boolean Warshall closure over the undirected edges.
template <class W>
std::vector<std::vector<bool>> reachability(const EdgeList<W>& g) {
  const std::size_t n = g.num_vertices;
  std::vector<std::vector<bool>> r(n, std::vector<bool>(n, false));
  for (std::size_t v = 0; v < n; ++v) r[v][v] = true;
  for (const auto& e : g.edges) {
    r[e.u][e.v] = true;
    r[e.v][e.u] = true;
  }
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (!r[i][k]) continue;
      for (std::size_t j = 0; j < n; ++j) {
        if (r[k][j]) r[i][j] = true;
      }
    }
  }
  return r;
}

/// Reference canonical labels: scanning vertices in increasing order, an unlabelled vertex opens
/// the next label and every vertex reachable from it receives that label.
std::vector<std::uint32_t> brute_labels(const std::vector<std::vector<bool>>& r) {
  const std::size_t n = r.size();
  std::vector<std::uint32_t> label(n, std::numeric_limits<std::uint32_t>::max());
  std::uint32_t next = 0;
  for (std::size_t v = 0; v < n; ++v) {
    if (label[v] != std::numeric_limits<std::uint32_t>::max()) continue;
    for (std::size_t x = 0; x < n; ++x) {
      if (r[v][x]) label[x] = next;
    }
    ++next;
  }
  return label;
}

/// Reference induced subgraph on the vertices with in[v], renumbered in increasing order.
template <class W>
sssp::InducedSubgraph<W> brute_induced(const EdgeList<W>& g, const std::vector<bool>& in) {
  sssp::InducedSubgraph<W> s;
  std::vector<Vertex> id(g.num_vertices, kNoVertex);
  for (std::size_t v = 0; v < g.num_vertices; ++v) {
    if (!in[v]) continue;
    id[v] = static_cast<Vertex>(s.original_vertex.size());
    s.original_vertex.push_back(static_cast<Vertex>(v));
  }
  s.graph.num_vertices = s.original_vertex.size();
  for (std::size_t i = 0; i < g.edges.size(); ++i) {
    const auto& e = g.edges[i];
    if (!in[e.u] || !in[e.v]) continue;
    s.graph.edges.push_back({id[e.u], id[e.v], e.weight});
    s.original_edge.push_back(static_cast<EdgeId>(i));
  }
  return s;
}

template <class W>
void require_same_subgraph(const sssp::InducedSubgraph<W>& got,
                           const sssp::InducedSubgraph<W>& want) {
  REQUIRE(got.graph.num_vertices == want.graph.num_vertices);
  REQUIRE(got.original_vertex == want.original_vertex);
  REQUIRE(got.original_edge == want.original_edge);
  REQUIRE(got.graph.edges.size() == want.graph.edges.size());
  for (std::size_t i = 0; i < got.graph.edges.size(); ++i) {
    REQUIRE(same_edge(got.graph.edges[i], want.graph.edges[i]));
  }
}

}  // namespace

TEMPLATE_TEST_CASE("keep_lightest_per_pair equals the brute-force (weight, index) minimum",
                   "[graph][edge_list][property]", std::uint64_t, double) {
  using W = TestType;
  for (const Weights kind : {Weights::kSmall, Weights::kHuge}) {
    for (int seed = 0; seed < kSeeds; ++seed) {
      SplitMix64 rng(static_cast<std::uint64_t>(seed) * 2 + (kind == Weights::kHuge ? 1 : 0));
      const auto g = random_multigraph<W>(rng, kind);
      for (const LoopPolicy loops : {LoopPolicy::kKeep, LoopPolicy::kDrop}) {
        INFO("seed " << seed << " huge " << (kind == Weights::kHuge) << " drop "
                     << (loops == LoopPolicy::kDrop));
        const auto r = sssp::keep_lightest_per_pair(g, loops);
        const auto want = brute_lightest(g, loops);
        REQUIRE(r.graph.num_vertices == g.num_vertices);
        REQUIRE(r.original_edge == want);
        REQUIRE(r.graph.edges.size() == want.size());
        for (std::size_t i = 0; i < want.size(); ++i) {
          REQUIRE(same_edge(r.graph.edges[i], g.edges[want[i]]));
        }
        // At most one edge per unordered pair remains.
        for (std::size_t i = 0; i < r.graph.edges.size(); ++i) {
          for (std::size_t j = i + 1; j < r.graph.edges.size(); ++j) {
            const auto& a = r.graph.edges[i];
            const auto& b = r.graph.edges[j];
            REQUIRE_FALSE(same_pair(a.u, a.v, b.u, b.v));
          }
        }
      }
    }
  }
}

TEMPLATE_TEST_CASE("keep_lightest_per_pair preserves all-pairs shortest-path distances",
                   "[graph][edge_list][property]", std::uint64_t, double) {
  using W = TestType;
  for (int seed = 0; seed < kSeeds; ++seed) {
    SplitMix64 rng(0x5EED0000ULL + static_cast<std::uint64_t>(seed));
    const auto g = random_multigraph<W>(rng, Weights::kSmall);
    const auto before = sssp::test::floyd_warshall<W>(g.num_vertices, g.edges);
    for (const LoopPolicy loops : {LoopPolicy::kKeep, LoopPolicy::kDrop}) {
      INFO("seed " << seed);
      const auto r = sssp::keep_lightest_per_pair(g, loops);
      const auto after = sssp::test::floyd_warshall<W>(r.graph.num_vertices, r.graph.edges);
      REQUIRE(after == before);
    }
  }
}

TEMPLATE_TEST_CASE("connected_components equals the brute-force closure labelling",
                   "[graph][edge_list][property]", std::uint64_t, double) {
  using W = TestType;
  for (int seed = 0; seed < kSeeds; ++seed) {
    SplitMix64 rng(0xC0C0ULL + static_cast<std::uint64_t>(seed));
    const auto g = random_multigraph<W>(rng, Weights::kSmall);
    INFO("seed " << seed);
    const auto reach = reachability(g);
    const auto want = brute_labels(reach);
    const auto c = sssp::connected_components(g);
    REQUIRE(c.label == want);
    std::size_t count = 0;
    for (const auto l : want) count = std::max<std::size_t>(count, std::size_t{l} + 1);
    REQUIRE(c.count == count);
    std::vector<std::size_t> size(count, 0);
    for (const auto l : want) ++size[l];
    REQUIRE(c.size == size);
  }
}

TEMPLATE_TEST_CASE("relabel_vertices equals the brute-force mapped edge list",
                   "[graph][edge_list][property]", std::uint64_t, double) {
  using W = TestType;
  for (int seed = 0; seed < kSeeds; ++seed) {
    SplitMix64 rng(0xAB1EULL + static_cast<std::uint64_t>(seed));
    const auto g = random_multigraph<W>(rng, Weights::kSmall);
    const std::size_t n = g.num_vertices;
    // Random injective map into [0, n + extra): a shuffled id pool, some vertices dropped.
    const std::size_t new_n = n + rng.below(3);
    std::vector<Vertex> pool(new_n);
    for (std::size_t i = 0; i < new_n; ++i) pool[i] = static_cast<Vertex>(i);
    for (std::size_t i = new_n; i > 1; --i) std::swap(pool[i - 1], pool[rng.below(i)]);
    const bool drops = rng.below(2) == 0;
    std::vector<Vertex> new_id(n);
    for (std::size_t v = 0; v < n; ++v) {
      new_id[v] = (drops && rng.below(4) == 0) ? kNoVertex : pool[v];
    }
    INFO("seed " << seed << " drops " << drops);
    const auto r = sssp::relabel_vertices(g, new_id, new_n);
    REQUIRE(r.graph.num_vertices == new_n);
    std::vector<EdgeId> want_ids;
    std::vector<Edge<W>> want_edges;
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
      const auto& e = g.edges[i];
      if (new_id[e.u] == kNoVertex || new_id[e.v] == kNoVertex) continue;
      want_ids.push_back(static_cast<EdgeId>(i));
      want_edges.push_back({new_id[e.u], new_id[e.v], e.weight});
    }
    REQUIRE(r.original_edge == want_ids);
    REQUIRE(r.graph.edges.size() == want_edges.size());
    for (std::size_t i = 0; i < want_edges.size(); ++i) {
      REQUIRE(same_edge(r.graph.edges[i], want_edges[i]));
    }
    if (!drops) {
      // A pure relabelling preserves every distance under the map.
      const auto d = sssp::test::floyd_warshall<W>(n, g.edges);
      const auto d2 = sssp::test::floyd_warshall<W>(new_n, r.graph.edges);
      for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = 0; b < n; ++b) REQUIRE(d2[new_id[a]][new_id[b]] == d[a][b]);
      }
    }
  }
}

TEMPLATE_TEST_CASE("component subgraphs equal the brute-force induced subgraphs",
                   "[graph][edge_list][property]", std::uint64_t, double) {
  using W = TestType;
  for (int seed = 0; seed < kSeeds; ++seed) {
    SplitMix64 rng(0xC0DEULL + static_cast<std::uint64_t>(seed));
    const auto g = random_multigraph<W>(rng, Weights::kSmall);
    const std::size_t n = g.num_vertices;
    INFO("seed " << seed);
    const auto reach = reachability(g);
    const auto labels = brute_labels(reach);
    const auto c = sssp::connected_components(g);
    const auto d = sssp::test::floyd_warshall<W>(n, g.edges);

    // Every component, through extract_component and through component_of of each member.
    for (std::uint32_t l = 0; l < c.count; ++l) {
      std::vector<bool> in(n);
      for (std::size_t v = 0; v < n; ++v) in[v] = labels[v] == l;
      const auto want = brute_induced(g, in);
      const auto got = sssp::extract_component(g, c, l);
      require_same_subgraph(got, want);
      for (const Vertex v : want.original_vertex)
        require_same_subgraph(sssp::component_of(g, v), want);
      // Distances inside a component are those of the whole graph.
      const auto ds = sssp::test::floyd_warshall<W>(got.graph.num_vertices, got.graph.edges);
      for (std::size_t a = 0; a < got.original_vertex.size(); ++a) {
        for (std::size_t b = 0; b < got.original_vertex.size(); ++b) {
          REQUIRE(ds[a][b] == d[got.original_vertex[a]][got.original_vertex[b]]);
        }
      }
    }

    // Largest: most vertices, ties by the smallest label.
    std::vector<std::size_t> size(c.count, 0);
    for (const auto l : labels) ++size[l];
    std::uint32_t best = 0;
    for (std::uint32_t l = 1; l < c.count; ++l) {
      if (size[l] > size[best]) best = l;
    }
    std::vector<bool> in(n);
    for (std::size_t v = 0; v < n; ++v) in[v] = c.count > 0 && labels[v] == best;
    require_same_subgraph(sssp::largest_component(g), brute_induced(g, in));
  }
}
