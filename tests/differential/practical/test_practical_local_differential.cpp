// Differential test of the Mode C (practical-optimized) local-hierarchy SSSP, task
// P5-PRACTICAL-LOCAL: practical::local_sssp against Bellman-Ford and binary-heap Dijkstra on
// the ORIGINAL graph, plus verify_sssp (shortest-path-tree certificate) and an independent
// union-find count of the contracted source component that fixes the expected FallbackRecord.
//
// Inputs: every generator family of generator_catalog() x kSeeds (20) seeds x 7 weight
// distributions x 7 option sets, for std::uint64_t and for dyadic double weights, on about 32
// vertices; the same on 150 vertices and 5 seeds on 400 vertices under [slow]; plus the
// boundary suite (n = 1, isolated source, disconnected, zero weights, multi-edges, self-loops)
// from every source and n = 0.
//
// Comparison rule (documented tolerance):
//   * std::uint64_t: exact equality.
//   * double: the weights are dyadic rationals k / 2^j with k < 2^31 and j in {0, 2, 10}, taken
//     from the integer instance of the same topology. Every sum the oracles and the algorithm
//     form is then an integer multiple of 2^-j below 2^53 * 2^-j, hence exact, so the tolerance
//     is ZERO: Tolerance::exact() and operator== on the distance vectors. Doubles whose sums
//     round are outside this test on purpose: Evaluate's supplemental offers and gaps are only
//     certified in exact arithmetic (evaluate.hpp, D-19 / D-15), see local_sssp.hpp "Numerics".
//
// Inputs left out on purpose (P3-EVAL-ARITH-FIX, in progress): weights near 2^63 for uint64
// (the mass of a node overflows inside Evaluate and raises std::overflow_error) and doubles at
// 2^53 + 1 (a slot becomes empty by rounding and raises std::logic_error). The weights here
// stay below 2^31 so that no mass, scaled threshold or path sum comes near either limit; the
// generator's wide_range(2^50/2^62) and the "heavy" weights of the baseline test are not used.
//
// Nothing is timed here.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "experimental/practical/local_sssp.hpp"
#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

using sssp::CsrGraph;
using sssp::EdgeList;
using sssp::Vertex;
using sssp::WeightTraits;
namespace pr = sssp::practical;
using U64 = std::uint64_t;

namespace {

constexpr std::uint64_t kSeeds = 20;

// ------------------------------------------------------------------ independent oracle parts

/// Size of the contracted source component, computed independently of the implementation:
/// union-find over zero edges restricted to the vertices Bellman-Ford reaches; edges are those
/// joining two different reached classes (loops discarded, parallel edges kept).
struct ContractedSize {
  std::uint64_t vertices = 0;
  std::uint64_t edges = 0;
};

template <class W>
ContractedSize contracted_size(const CsrGraph<W>& g, const std::vector<W>& bf_dist) {
  const std::size_t n = g.num_vertices();
  std::vector<std::size_t> root(n);
  std::iota(root.begin(), root.end(), std::size_t{0});
  const auto find = [&root](std::size_t x) {
    while (root[x] != x) x = root[x] = root[root[x]];
    return x;
  };
  const auto reached = [&](Vertex v) { return bf_dist[v] != WeightTraits<W>::infinity(); };
  for (const auto& e : g.edges()) {
    if (reached(e.u) && e.weight == W{0}) root[find(e.u)] = find(e.v);
  }
  ContractedSize sz;
  for (std::size_t v = 0; v < n; ++v) {
    if (reached(static_cast<Vertex>(v)) && find(v) == v) ++sz.vertices;
  }
  for (const auto& e : g.edges()) {
    if (reached(e.u) && find(e.u) != find(e.v)) ++sz.edges;
  }
  return sz;
}

/// The FallbackRecords the header promises for a contracted component of this size.
void require_fallbacks(const sssp::Diagnostics& d, const pr::PracticalOptions& o,
                       const ContractedSize& sz) {
  if (sz.vertices == 1) {
    REQUIRE(d.fallbacks.size() == 1);
    REQUIRE(d.fallbacks[0].reason == pr::kTrivialComponentReason);
    REQUIRE(d.fallbacks[0].vertices == 1);
    REQUIRE(d.fallbacks[0].edges == 0);
  } else if (sz.vertices <= o.direct_search_max_vertices) {
    REQUIRE(d.fallbacks.size() == 1);
    REQUIRE(d.fallbacks[0].reason == pr::direct_search_reason(o.direct_search_max_vertices));
    REQUIRE(d.fallbacks[0].vertices == sz.vertices);
    REQUIRE(d.fallbacks[0].edges == sz.edges);
  } else {
    REQUIRE(d.fallbacks.empty());
  }
}

// ------------------------------------------------------------------ option sets

struct NamedOptions {
  std::string name;
  pr::PracticalOptions options;
};

/// The default plus sets that force the hierarchy on every N >= 2 with each assembly, several
/// separators b (1, 2, 3, L_N and one above N) and the supplemental offers off.
std::vector<NamedOptions> option_sets(std::size_t n) {
  std::vector<NamedOptions> v;
  v.push_back({"default", {}});
  const auto add = [&v](std::string name, pr::Assembly a, std::uint64_t b, bool supp) {
    pr::PracticalOptions o;
    o.assembly = a;
    o.separator = b;
    o.direct_search_max_vertices = 0;
    o.supplemental_offers = supp;
    v.push_back({std::move(name), o});
  };
  add("exact b=L_N", pr::Assembly::kExactHeight, 0, true);
  add("exact b=1 no-supplemental", pr::Assembly::kExactHeight, 1, false);
  add("exact b=3", pr::Assembly::kExactHeight, 3, true);
  add("exact b=N+1", pr::Assembly::kExactHeight, n + 1, true);
  add("lambda-grid b=L_N", pr::Assembly::kLambdaGrid, 0, true);
  add("lambda-grid b=2", pr::Assembly::kLambdaGrid, 2, true);
  return v;
}

// ------------------------------------------------------------------ the check

struct Stats {
  std::uint64_t runs = 0;
  std::uint64_t hierarchy_runs = 0;  ///< no FallbackRecord: the construction handled it
  std::uint64_t direct_runs = 0;
  std::uint64_t trivial_runs = 0;
  std::uint64_t contracted_or_unreached = 0;  ///< inputs where N < n
};

/// Runs every option set on (g, s); each result must equal both oracles exactly, pass
/// verify_sssp with the reference under Tolerance::exact(), and carry the promised records.
template <class W>
void check_instance(const std::string& label, const CsrGraph<W>& g, Vertex s, Stats& st) {
  const std::size_t n = g.num_vertices();
  const auto bf = sssp::baseline::bellman_ford(g, s);
  const auto dij = sssp::baseline::dijkstra_binary_heap(g, s);
  REQUIRE(bf.dist == dij.dist);  // the inputs are exact: the oracles agree bit for bit
  const ContractedSize sz = contracted_size(g, bf.dist);
  if (sz.vertices < n) ++st.contracted_or_unreached;
  for (const auto& [name, o] : option_sets(n)) {
    INFO(label << " n=" << n << " m=" << g.num_edges() << " s=" << s << " options: " << name
               << " contracted N=" << sz.vertices << " m'=" << sz.edges);
    const auto r = pr::local_sssp(g, s, o);
    REQUIRE(r.dist == bf.dist);
    REQUIRE(r.dist == dij.dist);
    const auto rep = sssp::verify_sssp(g, s, r, &bf.dist, sssp::Tolerance::exact());
    INFO(rep.summary());
    REQUIRE(rep.ok);
    require_fallbacks(r.diag, o, sz);
    for (const auto& e : pr::option_log(o)) {
      REQUIRE(r.diag.counters.count("practical." + e.key + "=" + e.value) == 1);
    }
    ++st.runs;
    if (r.diag.fallbacks.empty()) {
      ++st.hierarchy_runs;
    } else if (sz.vertices == 1) {
      ++st.trivial_runs;
    } else {
      ++st.direct_runs;
    }
  }
}

// ------------------------------------------------------------------ inputs

/// The integer instance itself (uint64) or its dyadic image w / 2^shift (double, exact).
template <class W>
EdgeList<W> with_weight_type(const EdgeList<U64>& el, unsigned shift) {
  EdgeList<W> out{el.num_vertices, {}};
  out.edges.reserve(el.edges.size());
  for (const auto& e : el.edges) {
    if constexpr (std::is_same_v<W, double>) {
      out.edges.push_back(
          {e.u, e.v, static_cast<double>(e.weight) / static_cast<double>(U64{1} << shift)});
    } else {
      static_cast<void>(shift);
      out.edges.push_back({e.u, e.v, e.weight});
    }
  }
  return out;
}

struct NamedWeights {
  std::string name;
  sssp::gen::WeightSpec spec;
};

/// Integer weights below 2^31 (see the file comment): ties, a plain range, zero-heavy mixes
/// (contraction changes N), one value, and a wide range of magnitudes.
std::vector<NamedWeights> weight_specs() {
  using sssp::gen::WeightDist;
  const auto mk = [](WeightDist d, double lo, double hi) {
    sssp::gen::WeightSpec w;
    w.dist = d;
    w.lo = lo;
    w.hi = hi;
    return w;
  };
  std::vector<NamedWeights> v;
  v.push_back({"uniform_int[1,3]", mk(WeightDist::kUniformInt, 1, 3)});
  v.push_back({"uniform_int[0,2]", mk(WeightDist::kUniformInt, 0, 2)});
  v.push_back({"uniform_int[1,1000]", mk(WeightDist::kUniformInt, 1, 1000)});
  auto z3 = mk(WeightDist::kZeroMixed, 1, 5);
  z3.zero_probability = 0.3;
  v.push_back({"zero_mixed(0.3,[1,5])", z3});
  auto z9 = mk(WeightDist::kZeroMixed, 1, 5);
  z9.zero_probability = 0.9;
  v.push_back({"zero_mixed(0.9,[1,5])", z9});
  v.push_back({"constant(7)", mk(WeightDist::kConstant, 7, 7)});
  auto wide = mk(WeightDist::kWideRange, 1, 1);
  wide.wide_log2 = 30;
  v.push_back({"wide_range(2^30)", wide});
  return v;
}

/// Parameters for `family` on about n vertices. `variant` cycles the topology knobs; for
/// "regular" it also cycles the multi-edge policy, so parallel edges and loops are generated.
sssp::gen::GraphSpec make_spec(const std::string& family, std::size_t n, std::uint64_t variant) {
  sssp::gen::GraphSpec s;
  s.type = family;
  s.n = n;
  const double degrees[] = {1.5, 2.5, 4.0};
  s.avg_degree = degrees[variant % 3];
  s.degree = 3 + (variant % 2);
  if ((s.n * s.degree) % 2 != 0) ++s.n;
  const sssp::gen::MultiEdgePolicy policies[] = {sssp::gen::MultiEdgePolicy::kKeep,
                                                 sssp::gen::MultiEdgePolicy::kDrop,
                                                 sssp::gen::MultiEdgePolicy::kSwitchToSimple};
  s.multi = policies[variant % 3];
  s.attach = 1 + (variant % 3);
  s.exponent = 2.2 + 0.2 * static_cast<double>(variant % 4);
  return s;
}

/// family x weight distribution x seeds 1..seeds on about n vertices; the source cycles through
/// the first, the last and the middle vertex with the seed.
template <class W>
Stats run_family(const std::string& family, std::size_t n, std::uint64_t seeds = kSeeds) {
  Stats st;
  const auto specs = weight_specs();
  for (std::size_t wi = 0; wi < specs.size(); ++wi) {
    for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
      auto spec = make_spec(family, n, seed + wi);
      spec.options.weights = specs[wi].spec;
      spec.options.seed = seed;
      spec.options.connected = (seed % 2) == 0;
      const auto el = sssp::gen::generate<U64>(spec);
      const unsigned shifts[] = {0, 2, 10};
      const auto g = sssp::to_csr(with_weight_type<W>(el, shifts[seed % 3]));
      const std::size_t nv = g.num_vertices();
      REQUIRE(nv > 0);
      const std::string label = family + " weights=" + specs[wi].name +
                                " seed=" + std::to_string(seed) +
                                " shift=" + std::to_string(shifts[seed % 3]);
      const std::size_t sources[] = {0, nv - 1, nv / 2};
      check_instance<W>(label, g, static_cast<Vertex>(sources[(seed + wi) % 3]), st);
    }
  }
  return st;
}

std::vector<std::string> families() {
  std::vector<std::string> f;
  for (const auto& info : sssp::gen::generator_catalog()) f.emplace_back(info.name);
  return f;
}

}  // namespace

// ------------------------------------------------------------------ test cases

TEMPLATE_TEST_CASE("practical-local equals Bellman-Ford and Dijkstra on the boundary suite",
                   "[differential][practical][local]", std::uint64_t, double) {
  using W = TestType;
  Stats st;
  std::size_t graphs = 0;
  for (const auto& ng : sssp::gen::small_graph_suite<W>()) {
    if constexpr (std::is_same_v<W, double>) {
      // extreme_ratio uses 1e-12 and 1e12 for double: not dyadic, its sums round. It is covered
      // for uint64 (1 and 2^50) and by the dyadic wide_range family below.
      if (ng.name == "extreme_ratio") continue;
    }
    const auto g = sssp::to_csr(ng.graph);
    for (Vertex s = 0; s < g.num_vertices(); ++s) {
      check_instance<W>("small_graph_suite:" + ng.name, g, s, st);
    }
    ++graphs;
  }
  CHECK(graphs >= 16);
  CHECK(st.hierarchy_runs > 0);
  CHECK(st.direct_runs > 0);
  CHECK(st.trivial_runs > 0);
  CHECK(st.contracted_or_unreached > 0);
}

TEMPLATE_TEST_CASE("practical-local: hand-made boundary graphs", "[differential][practical][local]",
                   std::uint64_t, double) {
  using W = TestType;
  using E = std::vector<sssp::Edge<W>>;
  Stats st;
  const auto run = [&st](const std::string& name, std::size_t n, const E& edges) {
    const auto g = CsrGraph<W>::from_edges(n, edges);
    for (Vertex s = 0; s < n; ++s) check_instance<W>(name, g, s, st);
  };
  run("n=1", 1, {});
  run("n=1 with loops", 1, {{0, 0, W{0}}, {0, 0, W{3}}});
  run("n=2 no edge", 2, {});
  run("n=2 zero edge", 2, {{0, 1, W{0}}});
  run("n=2 triple edge", 2, {{0, 1, W{4}}, {1, 0, W{2}}, {0, 1, W{2}}});
  run("zero and positive parallel edges", 3, {{0, 1, W{0}}, {0, 1, W{5}}, {1, 2, W{1}}});
  run("all zero", 5, {{0, 1, W{0}}, {1, 2, W{0}}, {2, 3, W{0}}, {3, 4, W{0}}, {4, 0, W{0}}});
  run("zero classes joined by parallel edges and loops", 7,
      {{0, 1, W{0}},
       {2, 3, W{0}},
       {1, 2, W{6}},
       {0, 3, W{4}},
       {0, 2, W{4}},
       {3, 3, W{1}},
       {3, 4, W{2}},
       {4, 4, W{0}},
       {4, 5, W{2}},
       {5, 0, W{9}},
       {1, 5, W{8}}});
  run("two components and an isolated vertex", 9,
      {{0, 1, W{2}},
       {1, 2, W{2}},
       {2, 0, W{5}},
       {3, 4, W{1}},
       {4, 5, W{0}},
       {5, 6, W{1}},
       {6, 3, W{1}},
       {6, 7, W{3}}});
  // Clusters of light edges joined by heavy ones (unaligned blocks under D-8), 2^30 apart.
  {
    E e;
    const std::size_t clusters = 6, size = 5;
    for (std::size_t c = 0; c < clusters; ++c) {
      for (std::size_t i = 0; i + 1 < size; ++i) {
        e.push_back({static_cast<Vertex>(c * size + i), static_cast<Vertex>(c * size + i + 1),
                     static_cast<W>(1 + (c + i) % 3)});
      }
      e.push_back({static_cast<Vertex>(c * size), static_cast<Vertex>(c * size + size - 1), W{2}});
      if (c + 1 < clusters) {
        e.push_back({static_cast<Vertex>(c * size + 2), static_cast<Vertex>((c + 1) * size),
                     static_cast<W>((U64{1} << 30) + c)});
      }
    }
    run("light clusters joined by heavy edges", clusters * size, e);
  }
  CHECK(st.hierarchy_runs > 0);
  CHECK(st.trivial_runs > 0);

  const auto empty = CsrGraph<W>::from_edges(0, {});
  CHECK_THROWS_AS(pr::local_sssp(empty, 0, pr::PracticalOptions{}), std::out_of_range);
}

TEMPLATE_TEST_CASE("practical-local equals Bellman-Ford and Dijkstra on every generator family",
                   "[differential][practical][local]", std::uint64_t, double) {
  using W = TestType;
  Stats total;
  for (const auto& family : families()) {
    const Stats st = run_family<W>(family, 32);
    INFO(family << ": runs=" << st.runs << " hierarchy=" << st.hierarchy_runs
                << " direct=" << st.direct_runs << " trivial=" << st.trivial_runs);
    // 7 weight distributions x 20 seeds x 7 option sets.
    CHECK(st.runs == weight_specs().size() * kSeeds * option_sets(0).size());
    CHECK(st.hierarchy_runs > 0);
    total.runs += st.runs;
    total.hierarchy_runs += st.hierarchy_runs;
    total.direct_runs += st.direct_runs;
    total.trivial_runs += st.trivial_runs;
    total.contracted_or_unreached += st.contracted_or_unreached;
  }
  CHECK(families().size() >= 10);
  CHECK(total.direct_runs > 0);
  CHECK(total.contracted_or_unreached > 0);
  // The construction, not the direct search, handles most runs of this test.
  CHECK(2 * total.hierarchy_runs > total.runs);
}

TEMPLATE_TEST_CASE("practical-local equals Bellman-Ford and Dijkstra on larger generated graphs",
                   "[differential][practical][local][slow]", std::uint64_t, double) {
  using W = TestType;
  Stats total;
  for (const auto& family : families()) {
    // 20 seeds at n = 150 and 5 seeds at n = 400 (Bellman-Ford and 7 option sets per instance).
    for (const auto& [n, seeds] : {std::pair<std::size_t, std::uint64_t>{150, kSeeds}, {400, 5}}) {
      const Stats st = run_family<W>(family, n, seeds);
      INFO(family << " n=" << n << ": runs=" << st.runs << " hierarchy=" << st.hierarchy_runs);
      CHECK(st.runs == weight_specs().size() * seeds * option_sets(0).size());
      CHECK(st.hierarchy_runs > 0);
      total.runs += st.runs;
      total.hierarchy_runs += st.hierarchy_runs;
    }
  }
  CHECK(2 * total.hierarchy_runs > total.runs);
}
