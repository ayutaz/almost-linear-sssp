// Differential test of every baseline SSSP algorithm.
//
// For each input (the boundary-case suite and every generator family x several weight
// distributions x several seeds, from several sources) every baseline must
//   (a) pass verify_sssp WITHOUT a reference (shortest-path-tree certificate), and
//   (b) agree with Bellman-Ford on dist: exactly for uint64; bit for bit for double as well,
//       because every relaxation algorithm computes the minimum over paths of the left-to-right
//       rounded path sum (docs/design/numerics.md Sec 2), so no tolerance is needed.
// For uint64 an independent saturating oracle (in this file, sharing no code with the library)
// also fixes the exact distances and which inputs must throw std::overflow_error. The two
// overflow behaviours differ ON PURPOSE (docs/design/baselines.md, numerics.md Sec 6):
//   * Bellman-Ford throws iff some reachable vertex has a true distance >= 2^64 - 1;
//   * every Dijkstra throws iff some reachable u (with representable d*(u)) has an arc of
//     weight w with d*(u) + w >= 2^64 - 1, even if that arc is not on a shortest path.
//
// Seeds per (family, weight distribution): 2 by default (CI); set SSSP_DIFF_SEEDS=<k> to run k.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/verify.hpp"
#include "support/prng.hpp"

namespace {

using sssp::CsrGraph;
using sssp::EdgeList;
using sssp::SsspResult;
using sssp::Vertex;
using U64 = std::uint64_t;
constexpr U64 kMax = std::numeric_limits<U64>::max();

// ------------------------------------------------------------------ algorithms under test

template <class W>
struct Algorithm {
  std::string name;
  std::function<SsspResult<W>(const CsrGraph<W>&, Vertex)> run;
};

template <class W>
std::vector<Algorithm<W>> dijkstras() {
  std::vector<Algorithm<W>> a{
      {"dijkstra_binary_heap",
       [](const CsrGraph<W>& g, Vertex s) { return sssp::baseline::dijkstra_binary_heap(g, s); }},
      {"dijkstra_4ary_heap",
       [](const CsrGraph<W>& g, Vertex s) { return sssp::baseline::dijkstra_4ary_heap(g, s); }},
      {"dijkstra_fibonacci_heap",
       [](const CsrGraph<W>& g, Vertex s) {
         return sssp::baseline::dijkstra_fibonacci_heap(g, s);
       }},
  };
  if constexpr (std::is_same_v<W, U64>) {
    a.push_back({"dijkstra_radix_heap", [](const CsrGraph<W>& g, Vertex s) {
                   return sssp::baseline::dijkstra_radix_heap(g, s);
                 }});
  }
  return a;
}

// ------------------------------------------------------------------ independent uint64 oracle

/// Exact distances saturated at 2^64 - 1 ("not representable"), plus reachability. Round-based
/// relaxation over the edge list with saturating addition; min and saturating + commute with
/// saturation, so sat[v] = min(d*(v), 2^64 - 1) exactly.
struct SaturatedOracle {
  std::vector<U64> sat;
  std::vector<char> reachable;
};

U64 sat_add(U64 a, U64 b) {
  return a > kMax - b ? kMax : a + b;
}

SaturatedOracle saturated_oracle(const EdgeList<U64>& el, Vertex s) {
  const std::size_t n = el.num_vertices;
  SaturatedOracle o{std::vector<U64>(n, kMax), std::vector<char>(n, 0)};
  o.sat[s] = 0;
  o.reachable[s] = 1;
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& e : el.edges) {
      for (const auto [a, b] : {std::pair{e.u, e.v}, std::pair{e.v, e.u}}) {
        if (o.reachable[a] == 0) continue;
        const U64 c = sat_add(o.sat[a], e.weight);
        if (o.reachable[b] == 0 || c < o.sat[b]) {
          o.reachable[b] = 1;
          o.sat[b] = c;
          changed = true;
        }
      }
    }
  }
  return o;
}

struct OverflowExpectation {
  bool bellman_ford_throws = false;  ///< some reachable true distance >= 2^64 - 1
  bool dijkstra_throws = false;      ///< some scanned arc sum d*(u) + w >= 2^64 - 1
  std::vector<U64> dist;             ///< exact distances (infinity = unreachable) if representable
};

OverflowExpectation expect_u64(const EdgeList<U64>& el, Vertex s) {
  const auto o = saturated_oracle(el, s);
  OverflowExpectation x;
  x.dist.assign(el.num_vertices, sssp::WeightTraits<U64>::infinity());
  for (std::size_t v = 0; v < el.num_vertices; ++v) {
    if (o.reachable[v] == 0) continue;
    if (o.sat[v] == kMax) x.bellman_ford_throws = true;
    x.dist[v] = o.sat[v];
  }
  for (const auto& e : el.edges) {
    for (const Vertex u : {e.u, e.v}) {
      if (o.reachable[u] != 0 && o.sat[u] != kMax && sat_add(o.sat[u], e.weight) == kMax) {
        x.dijkstra_throws = true;
      }
    }
  }
  return x;
}

// ------------------------------------------------------------------ the differential check

template <class W>
std::string describe(const std::string& label, const EdgeList<W>& el, Vertex s) {
  std::string out = label + " n=" + std::to_string(el.num_vertices) +
                    " m=" + std::to_string(el.edges.size()) + " source=" + std::to_string(s);
  if (el.edges.size() <= 40) {
    out += " edges:";
    for (const auto& e : el.edges) {
      out += " {" + std::to_string(e.u) + "," + std::to_string(e.v) + "," +
             std::to_string(e.weight) + "}";
    }
  }
  return out;
}

struct Stats {
  std::uint64_t instances = 0;
  std::uint64_t runs = 0;
  std::uint64_t bf_overflows = 0;
  std::uint64_t dijkstra_overflows = 0;
};

/// Runs every baseline on (el, s) and checks (a) and (b). Returns false on the first mismatch
/// (after recording it with FAIL_CHECK).
template <class W>
bool check_instance(const std::string& label, const EdgeList<W>& el, Vertex s, Stats& st) {
  const auto g = sssp::to_csr(el);
  ++st.instances;
  const std::string where = describe(label, el, s);

  OverflowExpectation ex;
  if constexpr (std::is_same_v<W, U64>) {
    ex = expect_u64(el, s);
    // A non-representable distance is reached through a finite vertex whose arc overflows.
    if (ex.bellman_ford_throws && !ex.dijkstra_throws) {
      FAIL_CHECK("oracle inconsistency: " << where);
      return false;
    }
  }

  SsspResult<W> bf;
  bool bf_threw = false;
  try {
    bf = sssp::baseline::bellman_ford(g, s);
  } catch (const std::overflow_error&) {
    bf_threw = true;
  }
  ++st.runs;
  if (bf_threw != ex.bellman_ford_throws) {
    FAIL_CHECK("bellman_ford overflow behaviour: threw=" << bf_threw << " expected="
                                                         << ex.bellman_ford_throws << "\n"
                                                         << where);
    return false;
  }
  if (bf_threw) {
    ++st.bf_overflows;
  } else {
    const auto rep = sssp::verify_sssp(g, s, bf);
    if (!rep.ok) {
      FAIL_CHECK("bellman_ford fails verify_sssp: " << where << "\n" << rep.summary());
      return false;
    }
    if constexpr (std::is_same_v<W, U64>) {
      if (bf.dist != ex.dist) {
        FAIL_CHECK("bellman_ford differs from the saturating oracle: " << where);
        return false;
      }
    }
  }

  for (const auto& alg : dijkstras<W>()) {
    SsspResult<W> r;
    bool threw = false;
    try {
      r = alg.run(g, s);
    } catch (const std::overflow_error&) {
      threw = true;
    }
    ++st.runs;
    if (threw != ex.dijkstra_throws) {
      FAIL_CHECK(alg.name << " overflow behaviour: threw=" << threw
                          << " expected=" << ex.dijkstra_throws << "\n"
                          << where);
      return false;
    }
    if (threw) {
      ++st.dijkstra_overflows;
      continue;
    }
    const auto rep = sssp::verify_sssp(g, s, r);
    if (!rep.ok) {
      FAIL_CHECK(alg.name << " fails verify_sssp: " << where << "\n" << rep.summary());
      return false;
    }
    // Dijkstra did not throw => no arc sum overflowed => all distances representable.
    if (bf_threw) {
      FAIL_CHECK(alg.name << " finished but Bellman-Ford overflowed: " << where);
      return false;
    }
    if (r.dist != bf.dist) {
      for (std::size_t v = 0; v < r.dist.size(); ++v) {
        if (r.dist[v] != bf.dist[v]) {
          FAIL_CHECK(alg.name << " dist[" << v << "] = " << r.dist[v]
                              << " but Bellman-Ford = " << bf.dist[v] << "\n"
                              << where);
          break;
        }
      }
      return false;
    }
  }
  return true;
}

// ------------------------------------------------------------------ inputs

std::uint64_t seeds_per_combination() {
  const char* env = std::getenv("SSSP_DIFF_SEEDS");
  if (env == nullptr || *env == '\0') return 2;
  char* end = nullptr;
  const unsigned long long k = std::strtoull(env, &end, 10);
  if (end == env || *end != '\0' || k == 0) {
    throw std::invalid_argument("SSSP_DIFF_SEEDS must be a positive integer");
  }
  return k;
}

struct NamedWeights {
  std::string name;
  sssp::gen::WeightSpec spec;
  /// If true, every weight w is replaced by min(w, 3) * 2^62 after generation (exact for both
  /// weight types). uint64 path sums then overflow after a few heavy hops, which the generator's
  /// own distributions (at most 2^62 per edge, and mostly far smaller) almost never reach.
  bool heavy = false;
};

template <class W>
void make_heavy(EdgeList<W>& el) {
  for (auto& e : el.edges) {
    const W w = e.weight < W{3} ? e.weight : W{3};
    if constexpr (std::is_integral_v<W>) {
      e.weight = static_cast<W>(static_cast<U64>(w) << 62);
    } else {
      e.weight = w * 0x1p62;
    }
  }
}

std::vector<NamedWeights> weight_specs() {
  using sssp::gen::WeightDist;
  auto mk = [](WeightDist d) {
    sssp::gen::WeightSpec w;
    w.dist = d;
    return w;
  };
  std::vector<NamedWeights> v;
  auto ties = mk(WeightDist::kUniformInt);
  ties.lo = 0;
  ties.hi = 3;
  v.push_back({"uniform_int[0,3]", ties});
  auto big = mk(WeightDist::kUniformInt);
  big.lo = 1;
  big.hi = 0x1p53;
  v.push_back({"uniform_int[1,2^53]", big});
  auto real = mk(WeightDist::kUniformReal);
  real.lo = 0;
  real.hi = 10;
  v.push_back({"uniform_real[0,10]", real});
  auto logu = mk(WeightDist::kLogUniform);
  logu.lo = 1e-6;
  logu.hi = 1e6;
  v.push_back({"log_uniform[1e-6,1e6]", logu});
  auto expo = mk(WeightDist::kExponential);
  expo.mean = 10;
  v.push_back({"exponential(10)", expo});
  auto cst = mk(WeightDist::kConstant);
  cst.lo = 7;
  v.push_back({"constant(7)", cst});
  auto zm = mk(WeightDist::kZeroMixed);
  zm.lo = 1;
  zm.hi = 5;
  zm.zero_probability = 0.5;
  v.push_back({"zero_mixed(0.5,[1,5])", zm});
  auto wide = mk(WeightDist::kWideRange);
  wide.wide_log2 = 50;
  v.push_back({"wide_range(2^50)", wide});
  auto huge = mk(WeightDist::kWideRange);
  huge.wide_log2 = 62;  // uint64: path sums overflow 2^64 after a few hops
  v.push_back({"wide_range(2^62)", huge});
  v.push_back({"heavy(min(uniform_int[0,3],3)*2^62)", ties, true});
  return v;
}

/// Random but valid parameters for `type`, small enough for O(nm) Bellman-Ford.
sssp::gen::GraphSpec random_spec(const std::string& type, sssp::test::SplitMix64& rng) {
  sssp::gen::GraphSpec s;
  s.type = type;
  s.n = 1 + rng.below(60);
  if (type == "er") {
    const double max_deg = static_cast<double>(s.n - 1);
    const double choices[] = {0.5, 1.5, 3.0, 6.0};
    s.avg_degree = std::min(choices[rng.below(4)], max_deg);
  } else if (type.rfind("regular", 0) == 0) {
    s.type = "regular";
    s.n = 2 + rng.below(40);
    s.degree = 1 + rng.below(std::min<std::uint64_t>(4, s.n - 1));
    if ((s.n * s.degree) % 2 != 0) ++s.n;
    s.multi = sssp::gen::parse_multi_edge_policy(type.substr(std::string("regular_").size()));
  } else if (type == "ba") {
    s.n = 2 + rng.below(60);
    s.attach = 1 + rng.below(std::min<std::uint64_t>(3, s.n - 1));
  } else if (type == "grid" || type == "road") {
    s.rows = 1 + rng.below(8);
    s.cols = 1 + rng.below(8);
  } else if (type == "chung_lu") {
    // chung_lu requires m = round(n * avg_degree / 2) <= q = floor(n(n-1)/2) / 2 (integer
    // division, src/graph/generators.cpp); avg_degree <= (2q - 1) / n keeps n*avg/2 <= q - 1/2.
    s.n = 4 + rng.below(60);
    const auto q = static_cast<double>((s.n * (s.n - 1) / 2) / 2);
    s.avg_degree = std::min(1.0 + static_cast<double>(rng.below(3)),
                            (2.0 * q - 1.0) / static_cast<double>(s.n));
    s.exponent = 2.1 + rng.unit();
  } else if (type == "adversarial_decrease_key") {
    s.n = 1 + rng.below(30);
  } else if (type.rfind("geometric", 0) == 0) {
    // radius 0 (no edges), sparse, dense, and >= sqrt(2) (complete graph)
    s.type = "geometric";
    const double radii[] = {0.0, 0.1, 0.25, 0.6, 1.5};
    s.radius = radii[rng.below(5)];
    if (type == "geometric_euclidean") {
      s.geometric.weights = sssp::gen::GeometricWeights::kEuclidean;
    }
  }
  return s;
}

const std::vector<std::string> kFamilies{"er",
                                         "regular_keep",
                                         "regular_drop",
                                         "regular_switch",
                                         "ba",
                                         "grid",
                                         "road",
                                         "tree",
                                         "chung_lu",
                                         "adversarial_decrease_key",
                                         "adversarial_long_path",
                                         "geometric_random",
                                         "geometric_euclidean"};

/// Stable per-family seed component (std::hash is implementation-defined).
std::uint64_t family_id(const std::string& family) {
  for (std::size_t i = 0; i < kFamilies.size(); ++i) {
    if (kFamilies[i] == family) return i + 1;
  }
  throw std::invalid_argument("unknown family " + family);
}

template <class W>
void run_family(const std::string& family) {
  const std::uint64_t seeds = seeds_per_combination();
  Stats st;
  const auto specs = weight_specs();
  for (std::size_t wi = 0; wi < specs.size(); ++wi) {
    for (std::uint64_t seed = 1; seed <= seeds; ++seed) {
      sssp::test::SplitMix64 rng(family_id(family) * 0x10000 + wi * 0x100 + seed);
      auto spec = random_spec(family, rng);
      spec.options.weights = specs[wi].spec;
      spec.options.seed = seed;
      spec.options.connected = (seed % 2) == 0;
      auto el = sssp::gen::generate<W>(spec);
      if (specs[wi].heavy) make_heavy(el);
      const std::string label = family + " weights=" + specs[wi].name +
                                " seed=" + std::to_string(seed) +
                                (spec.options.connected ? " connected" : "");
      const std::size_t n = el.num_vertices;
      if (n == 0) continue;
      const std::vector<Vertex> sources{0, static_cast<Vertex>(n - 1),
                                        static_cast<Vertex>(rng.below(n))};
      for (const Vertex s : sources) {
        if (!check_instance<W>(label, el, s, st)) {
          REQUIRE(false);
        }
      }
    }
  }
  INFO("instances=" << st.instances << " runs=" << st.runs << " bellman_ford overflows="
                    << st.bf_overflows << " dijkstra overflows=" << st.dijkstra_overflows);
  CHECK(st.instances >= specs.size() * seeds);
  if constexpr (std::is_same_v<W, U64>) {
    // The heavy weights must really exercise the overflow paths (not only the happy path).
    CHECK(st.dijkstra_overflows > 0);
  }
}

}  // namespace

// ------------------------------------------------------------------ test cases

TEMPLATE_TEST_CASE("baselines agree on the small graph suite from every source",
                   "[differential][baselines]", std::uint64_t, double) {
  using W = TestType;
  Stats st;
  for (const auto& ng : sssp::gen::small_graph_suite<W>()) {
    for (Vertex s = 0; s < ng.graph.num_vertices; ++s) {
      REQUIRE(check_instance<W>("small_graph_suite:" + ng.name, ng.graph, s, st));
    }
  }
  CHECK(st.instances > 17);
}

TEMPLATE_TEST_CASE("baselines agree on generated graphs: er", "[differential][baselines]",
                   std::uint64_t, double) {
  run_family<TestType>("er");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: regular (keep/drop/switch)",
                   "[differential][baselines]", std::uint64_t, double) {
  run_family<TestType>("regular_keep");
  run_family<TestType>("regular_drop");
  run_family<TestType>("regular_switch");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: ba", "[differential][baselines]",
                   std::uint64_t, double) {
  run_family<TestType>("ba");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: grid and road",
                   "[differential][baselines]", std::uint64_t, double) {
  run_family<TestType>("grid");
  run_family<TestType>("road");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: tree", "[differential][baselines]",
                   std::uint64_t, double) {
  run_family<TestType>("tree");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: chung_lu", "[differential][baselines]",
                   std::uint64_t, double) {
  run_family<TestType>("chung_lu");
}
TEMPLATE_TEST_CASE("baselines agree on generated graphs: adversarial families",
                   "[differential][baselines]", std::uint64_t, double) {
  run_family<TestType>("adversarial_decrease_key");
  run_family<TestType>("adversarial_long_path");
}

TEMPLATE_TEST_CASE("baselines agree on generated graphs: geometric (random/euclidean weights)",
                   "[differential][baselines]", std::uint64_t, double) {
  run_family<TestType>("geometric_random");
  run_family<TestType>("geometric_euclidean");
}

TEST_CASE("every generator family is covered by the differential test",
          "[differential][baselines]") {
  for (const auto& info : sssp::gen::generator_catalog()) {
    bool covered = false;
    for (const auto& f : kFamilies) covered = covered || f.rfind(info.name, 0) == 0;
    INFO("generator " << info.name);
    CHECK(covered);
  }
}

TEST_CASE("overflow: Bellman-Ford and Dijkstra differ by design (uint64)",
          "[differential][baselines][overflow]") {
  Stats st;
  const U64 big = U64{1} << 63;
  // 1. true distance 2^64 is not representable: every algorithm throws.
  const EdgeList<U64> unrep{3, {{0, 1, big}, {1, 2, big}}};
  // 2. single edge of weight 2^64 - 1 = infinity(): d*(1) is not representable.
  const EdgeList<U64> max_edge{2, {{0, 1, kMax}}};
  // 3. a heavy non-shortest arc back to the source: d* = (0, 1, 2) is representable, Bellman-Ford
  //    returns it; Dijkstra throws when it scans vertex 2 (2 + (2^64 - 2) reaches infinity()).
  const EdgeList<U64> heavy_back{3, {{0, 1, 1}, {1, 2, 1}, {0, 2, kMax - 1}}};
  // 4. d*(1) = 2^64 - 2, the largest finite distance: Bellman-Ford returns it, Dijkstra throws
  //    on the arc 1 -> 0 (the same edge scanned from the other end).
  const EdgeList<U64> largest_finite{2, {{0, 1, kMax - 1}}};
  // 5. a non-representable vertex behind a representable one.
  const EdgeList<U64> chain{4, {{0, 1, 5}, {1, 2, kMax - 6}, {2, 3, 1}}};

  struct Case {
    const char* name;
    const EdgeList<U64>* g;
    bool bf_throws;
    bool dijkstra_throws;
  };
  const Case cases[] = {
      {"unrepresentable", &unrep, true, true},
      {"max_edge", &max_edge, true, true},
      {"heavy_back", &heavy_back, false, true},
      {"largest_finite", &largest_finite, false, true},
      {"chain", &chain, true, true},
  };
  for (const auto& c : cases) {
    INFO(c.name);
    // The hand-derived expectation must match the independent oracle ...
    const auto ex = expect_u64(*c.g, 0);
    CHECK(ex.bellman_ford_throws == c.bf_throws);
    CHECK(ex.dijkstra_throws == c.dijkstra_throws);
    // ... and every implementation must follow it.
    CHECK(check_instance<U64>(c.name, *c.g, 0, st));
  }
  // The representable cases: Bellman-Ford returns the exact distances.
  CHECK(sssp::baseline::bellman_ford(sssp::to_csr(heavy_back), 0).dist ==
        std::vector<U64>{0, 1, 2});
  CHECK(sssp::baseline::bellman_ford(sssp::to_csr(largest_finite), 0).dist ==
        std::vector<U64>{0, kMax - 1});
}

TEMPLATE_TEST_CASE("every baseline rejects an invalid source with std::out_of_range",
                   "[differential][baselines]", std::uint64_t, double) {
  using W = TestType;
  const auto empty = CsrGraph<W>::from_edges(0, {});
  const auto two = CsrGraph<W>::from_edges(2, {{0, 1, W{1}}});
  CHECK_THROWS_AS(sssp::baseline::bellman_ford(empty, 0), std::out_of_range);
  CHECK_THROWS_AS(sssp::baseline::bellman_ford(two, 2), std::out_of_range);
  for (const auto& alg : dijkstras<W>()) {
    INFO(alg.name);
    CHECK_THROWS_AS(alg.run(empty, 0), std::out_of_range);
    CHECK_THROWS_AS(alg.run(two, 2), std::out_of_range);
    CHECK_THROWS_AS(alg.run(two, sssp::kNoVertex), std::out_of_range);
  }
}

TEST_CASE("verify_sssp rejects every single-vertex distance perturbation of a baseline result",
          "[differential][baselines][verify]") {
  // Independent review of the certificate: if dist[v] is changed by +-1 for one reachable v != s,
  // either v's tree edge or a neighbour's optimality constraint must break.
  std::uint64_t perturbations = 0;
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    sssp::gen::GraphSpec spec;
    spec.type = "er";
    spec.n = 30;
    spec.avg_degree = 3;
    spec.options.seed = seed;
    spec.options.weights.lo = 0;
    spec.options.weights.hi = 4;
    const auto g = sssp::to_csr(sssp::gen::generate<U64>(spec));
    const auto r = sssp::baseline::dijkstra_binary_heap(g, 0);
    REQUIRE(sssp::verify_sssp(g, 0, r).ok);
    for (Vertex v = 1; v < g.num_vertices(); ++v) {
      if (r.dist[v] == sssp::WeightTraits<U64>::infinity()) continue;
      for (const int delta : {-1, +1}) {
        if (delta < 0 && r.dist[v] == 0) continue;
        auto bad = r;
        bad.dist[v] = delta < 0 ? bad.dist[v] - 1 : bad.dist[v] + 1;
        INFO("seed " << seed << " vertex " << v << " delta " << delta);
        CHECK_FALSE(sssp::verify_sssp(g, 0, bad).ok);
        ++perturbations;
      }
    }
  }
  CHECK(perturbations > 100);
}
