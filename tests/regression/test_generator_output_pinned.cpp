// Pinned generator output (P2-UNIONFIND-DEDUP).
//
// Not a minimised counterexample (so there is no cases/*.json record): this test freezes the
// exact output of every generator family, of connect_components / count_components and of
// connected_components (src/graph/edge_list.cpp) for fixed seeds. The generators promise that
// a seed reproduces the same graph (include/sssp/generators.hpp), and existing tests, benchmarks
// and recorded results depend on that. Any refactoring of src/graph (for example replacing the
// union-find used by connect_components, the road_like connectivity repair and
// connected_components) must keep every digest below.
//
// The expected digests were computed with the code as of commit 44482b2, before the union-find
// in src/graph was shared. A digest is FNV-1a (64-bit) over num_vertices, the edge count and,
// for each edge in order, u, v and the weight bits (std::uint64_t value, or the IEEE-754 bit
// pattern for double). The generators avoid platform-dependent floating-point (det_log /
// det_exp, no FMA contraction), so the digests are the same on every supported toolchain.

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/edge_list.hpp"
#include "sssp/generators.hpp"

namespace {

using sssp::EdgeList;
using sssp::gen::GraphSpec;

class Fnv1a {
 public:
  void add(std::uint64_t x) {
    for (int i = 0; i < 8; ++i) {
      h_ ^= (x >> (8 * i)) & 0xffU;
      h_ *= 0x100000001b3ULL;
    }
  }
  [[nodiscard]] std::uint64_t value() const { return h_; }

 private:
  std::uint64_t h_ = 0xcbf29ce484222325ULL;
};

template <class W>
std::uint64_t weight_bits(W w) {
  if constexpr (std::is_same_v<W, double>) {
    return std::bit_cast<std::uint64_t>(w);
  } else {
    return static_cast<std::uint64_t>(w);
  }
}

template <class W>
std::uint64_t digest(const EdgeList<W>& g) {
  Fnv1a h;
  h.add(g.num_vertices);
  h.add(g.edges.size());
  for (const auto& e : g.edges) {
    h.add(e.u);
    h.add(e.v);
    h.add(weight_bits(e.weight));
  }
  return h.value();
}

/// Digest of connected_components(g): count, then every label, then every size.
template <class W>
std::uint64_t component_digest(const EdgeList<W>& g) {
  const sssp::ComponentLabels c = sssp::connected_components(g);
  Fnv1a h;
  h.add(c.count);
  for (const auto l : c.label) h.add(l);
  for (const auto s : c.size) h.add(s);
  h.add(sssp::gen::count_components(g));
  return h.value();
}

std::string hex(std::uint64_t x) {
  std::ostringstream os;
  os << "0x" << std::hex << std::setw(16) << std::setfill('0') << x;
  return os.str();
}

GraphSpec spec(std::string type, std::size_t n, std::uint64_t seed, bool connected = false) {
  GraphSpec s;
  s.type = std::move(type);
  s.n = n;
  s.options.seed = seed;
  s.options.connected = connected;
  return s;
}

/// The pinned specifications. Names are stable; each is generated for both weight types.
std::vector<std::pair<std::string, GraphSpec>> pinned_specs() {
  std::vector<std::pair<std::string, GraphSpec>> out;
  auto put = [&](std::string name, GraphSpec s) {
    out.emplace_back(std::move(name), std::move(s));
  };

  for (const std::size_t n : {0U, 1U, 2U, 17U, 200U, 1000U}) {
    for (const bool conn : {false, true}) {
      for (const double deg : {0.0, 1.0, 4.0}) {
        // G(n, m) needs m = round(n * deg / 2) <= n(n-1)/2, i.e. deg <= n - 1.
        if (deg > static_cast<double>(n > 0 ? n - 1 : 0)) continue;
        GraphSpec s = spec("er", n, 7 + n, conn);
        s.avg_degree = deg;
        put("er n=" + std::to_string(n) + " deg=" + std::to_string(deg) +
                " conn=" + std::to_string(conn),
            s);
      }
    }
  }
  {
    GraphSpec s = spec("er", 500, 11, true);
    s.avg_degree = 2.0;
    s.options.weights.dist = sssp::gen::WeightDist::kExponential;
    put("er n=500 exponential conn=1", s);
    s.options.weights.dist = sssp::gen::WeightDist::kLogUniform;
    put("er n=500 log_uniform conn=1", s);
    s.options.weights.dist = sssp::gen::WeightDist::kZeroMixed;
    put("er n=500 zero_mixed conn=1", s);
    s.options.weights.dist = sssp::gen::WeightDist::kWideRange;
    put("er n=500 wide_range conn=1", s);
    s.options.weights.dist = sssp::gen::WeightDist::kUniformReal;
    put("er n=500 uniform_real conn=1", s);
  }
  for (const std::size_t n : {4U, 10U, 200U}) {
    for (const auto multi : {sssp::gen::MultiEdgePolicy::kKeep, sssp::gen::MultiEdgePolicy::kDrop,
                             sssp::gen::MultiEdgePolicy::kSwitchToSimple}) {
      for (const bool conn : {false, true}) {
        GraphSpec s = spec("regular", n, 3 * n + 1, conn);
        s.degree = 3;
        s.multi = multi;
        put("regular n=" + std::to_string(n) + " d=3 multi=" +
                std::to_string(static_cast<int>(multi)) + " conn=" + std::to_string(conn),
            s);
      }
    }
  }
  {
    GraphSpec s = spec("regular", 300, 5, true);
    s.degree = 1;
    s.multi = sssp::gen::MultiEdgePolicy::kDrop;
    put("regular n=300 d=1 drop conn=1", s);
  }
  for (const std::size_t n : {3U, 50U, 500U}) {
    GraphSpec s = spec("ba", n, 13 + n);
    s.attach = 2;
    put("ba n=" + std::to_string(n) + " attach=2", s);
  }
  for (const auto [r, c] : {std::pair<std::size_t, std::size_t>{1, 1}, {3, 7}, {20, 20}}) {
    GraphSpec s = spec("grid", r * c, 17);
    s.rows = r;
    s.cols = c;
    put("grid " + std::to_string(r) + "x" + std::to_string(c), s);
  }
  // road_like restores removed streets with a union-find (Kruskal order); low keep
  // probabilities force many restorations.
  for (const auto [r, c] :
       {std::pair<std::size_t, std::size_t>{1, 1}, {2, 2}, {5, 9}, {30, 30}, {1, 40}}) {
    for (const double keep : {0.9, 0.5, 0.1, 0.0}) {
      GraphSpec s = spec("road", r * c, 19 + r * c);
      s.rows = r;
      s.cols = c;
      s.road.keep_probability = keep;
      put("road " + std::to_string(r) + "x" + std::to_string(c) + " keep=" + std::to_string(keep),
          s);
    }
  }
  for (const std::size_t n : {1U, 2U, 50U, 1000U})
    put("tree n=" + std::to_string(n), spec("tree", n, 23 + n));
  for (const std::size_t n : {100U, 1000U}) {
    for (const bool conn : {false, true}) {
      GraphSpec s = spec("chung_lu", n, 29 + n, conn);
      s.avg_degree = 3.0;
      s.exponent = 2.3;
      put("chung_lu n=" + std::to_string(n) + " conn=" + std::to_string(conn), s);
    }
  }
  for (const std::size_t n : {1U, 2U, 10U, 64U}) {
    put("adversarial_decrease_key n=" + std::to_string(n), spec("adversarial_decrease_key", n, 1));
    put("adversarial_long_path n=" + std::to_string(n), spec("adversarial_long_path", n, 1));
  }
  for (const std::size_t n : {1U, 50U, 500U}) {
    for (const auto gw :
         {sssp::gen::GeometricWeights::kFromSpec, sssp::gen::GeometricWeights::kEuclidean}) {
      for (const bool conn : {false, true}) {
        GraphSpec s = spec("geometric", n, 31 + n, conn);
        s.avg_degree = 2.0;
        s.geometric.weights = gw;
        put("geometric n=" + std::to_string(n) + " weights=" +
                std::to_string(static_cast<int>(gw)) + " conn=" + std::to_string(conn),
            s);
      }
    }
  }
  return out;
}

/// (name, digest) for every pinned output, in a fixed order.
std::vector<std::pair<std::string, std::uint64_t>> actual_digests() {
  std::vector<std::pair<std::string, std::uint64_t>> out;
  for (const auto& [name, s] : pinned_specs()) {
    const auto gu = sssp::gen::generate<std::uint64_t>(s);
    const auto gd = sssp::gen::generate<double>(s);
    out.emplace_back(name + " /u64", digest(gu));
    out.emplace_back(name + " /f64", digest(gd));
    out.emplace_back(name + " /u64 components", component_digest(gu));
    out.emplace_back(name + " /f64 components", component_digest(gd));
  }
  // erdos_renyi_gnm is not reachable through generate().
  for (const std::size_t n : {0U, 5U, 300U}) {
    sssp::gen::GenOptions o;
    o.seed = 37 + n;
    o.connected = true;
    const std::size_t m = n / 2;  // <= n(n-1)/2 for every n; leaves many components
    out.emplace_back("gnm n=" + std::to_string(n) + " /u64",
                     digest(sssp::gen::erdos_renyi_gnm<std::uint64_t>(n, m, o)));
    out.emplace_back("gnm n=" + std::to_string(n) + " /f64",
                     digest(sssp::gen::erdos_renyi_gnm<double>(n, m, o)));
  }
  // connect_components called directly on edge lists with many components, with the
  // components' smallest vertices out of edge order (so a union-find whose roots are not the
  // smallest vertex would number them differently if the labelling depended on the roots).
  {
    EdgeList<std::uint64_t> g{60, {}};
    for (sssp::Vertex v = 59; v >= 3; v -= 3)
      g.edges.push_back({v, static_cast<sssp::Vertex>(v - 3), 1});
    for (sssp::Vertex v = 40; v < 59; v += 2)
      g.edges.push_back({v, static_cast<sssp::Vertex>(v + 1), 2});
    out.emplace_back("direct components before /u64", component_digest(g));
    sssp::gen::connect_components(g, sssp::gen::WeightSpec{}, 41);
    out.emplace_back("direct connect_components /u64", digest(g));
    out.emplace_back("direct components after /u64", component_digest(g));
  }
  {
    EdgeList<double> g{1, {}};
    sssp::gen::connect_components(g, sssp::gen::WeightSpec{}, 43);
    out.emplace_back("direct connect_components n=1 /f64", digest(g));
    EdgeList<double> e{0, {}};
    sssp::gen::connect_components(e, sssp::gen::WeightSpec{}, 43);
    out.emplace_back("direct connect_components n=0 /f64", digest(e));
    out.emplace_back("direct components n=0 /f64", component_digest(e));
  }
  return out;
}

// Computed with the code before P2-UNIONFIND-DEDUP (commit 44482b2); never edit to match output.
const std::vector<std::pair<std::string, std::uint64_t>> kExpected = {
    {"er n=0 deg=0.000000 conn=0 /u64", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=0 /f64", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=0 /u64 components", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=0 /f64 components", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=1 /u64", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=1 /f64", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=1 /u64 components", 0x88201fb960ff6465ULL},
    {"er n=0 deg=0.000000 conn=1 /f64 components", 0x88201fb960ff6465ULL},
    {"er n=1 deg=0.000000 conn=0 /u64", 0x392209f14dea4c24ULL},
    {"er n=1 deg=0.000000 conn=0 /f64", 0x392209f14dea4c24ULL},
    {"er n=1 deg=0.000000 conn=0 /u64 components", 0x372cac50b2ba8344ULL},
    {"er n=1 deg=0.000000 conn=0 /f64 components", 0x372cac50b2ba8344ULL},
    {"er n=1 deg=0.000000 conn=1 /u64", 0x392209f14dea4c24ULL},
    {"er n=1 deg=0.000000 conn=1 /f64", 0x392209f14dea4c24ULL},
    {"er n=1 deg=0.000000 conn=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"er n=1 deg=0.000000 conn=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"er n=2 deg=0.000000 conn=0 /u64", 0x261c4b49872994e7ULL},
    {"er n=2 deg=0.000000 conn=0 /f64", 0x261c4b49872994e7ULL},
    {"er n=2 deg=0.000000 conn=0 /u64 components", 0xc6ec3c0bb6a7c5c4ULL},
    {"er n=2 deg=0.000000 conn=0 /f64 components", 0xc6ec3c0bb6a7c5c4ULL},
    {"er n=2 deg=1.000000 conn=0 /u64", 0x64464d86148e4dabULL},
    {"er n=2 deg=1.000000 conn=0 /f64", 0xc6c66e8345f0b8d1ULL},
    {"er n=2 deg=1.000000 conn=0 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=2 deg=1.000000 conn=0 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=2 deg=0.000000 conn=1 /u64", 0x4bc8f65f1afd40c3ULL},
    {"er n=2 deg=0.000000 conn=1 /f64", 0xc7ea328346e82f17ULL},
    {"er n=2 deg=0.000000 conn=1 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=2 deg=0.000000 conn=1 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=2 deg=1.000000 conn=1 /u64", 0x64464d86148e4dabULL},
    {"er n=2 deg=1.000000 conn=1 /f64", 0xc6c66e8345f0b8d1ULL},
    {"er n=2 deg=1.000000 conn=1 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=2 deg=1.000000 conn=1 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"er n=17 deg=0.000000 conn=0 /u64", 0x290366727f3bd034ULL},
    {"er n=17 deg=0.000000 conn=0 /f64", 0x290366727f3bd034ULL},
    {"er n=17 deg=0.000000 conn=0 /u64 components", 0xb790c2b1e271dc54ULL},
    {"er n=17 deg=0.000000 conn=0 /f64 components", 0xb790c2b1e271dc54ULL},
    {"er n=17 deg=1.000000 conn=0 /u64", 0xb1001a97222ccae8ULL},
    {"er n=17 deg=1.000000 conn=0 /f64", 0x367c4ec79df9211fULL},
    {"er n=17 deg=1.000000 conn=0 /u64 components", 0xc61a1d524ff18527ULL},
    {"er n=17 deg=1.000000 conn=0 /f64 components", 0xc61a1d524ff18527ULL},
    {"er n=17 deg=4.000000 conn=0 /u64", 0xd57e8f91cf49df7eULL},
    {"er n=17 deg=4.000000 conn=0 /f64", 0x16a81a807cfdafa7ULL},
    {"er n=17 deg=4.000000 conn=0 /u64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=4.000000 conn=0 /f64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=0.000000 conn=1 /u64", 0x39b120bf923d7817ULL},
    {"er n=17 deg=0.000000 conn=1 /f64", 0x6be38d5a056a13a9ULL},
    {"er n=17 deg=0.000000 conn=1 /u64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=0.000000 conn=1 /f64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=1.000000 conn=1 /u64", 0xf0bdbad5acc2b5ddULL},
    {"er n=17 deg=1.000000 conn=1 /f64", 0xd186b2f88d39b67eULL},
    {"er n=17 deg=1.000000 conn=1 /u64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=1.000000 conn=1 /f64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=4.000000 conn=1 /u64", 0xd57e8f91cf49df7eULL},
    {"er n=17 deg=4.000000 conn=1 /f64", 0x16a81a807cfdafa7ULL},
    {"er n=17 deg=4.000000 conn=1 /u64 components", 0xb3f3259eccb04f54ULL},
    {"er n=17 deg=4.000000 conn=1 /f64 components", 0xb3f3259eccb04f54ULL},
    {"er n=200 deg=0.000000 conn=0 /u64", 0x3ea12408497a572dULL},
    {"er n=200 deg=0.000000 conn=0 /f64", 0x3ea12408497a572dULL},
    {"er n=200 deg=0.000000 conn=0 /u64 components", 0x2ad999fe2dbd8665ULL},
    {"er n=200 deg=0.000000 conn=0 /f64 components", 0x2ad999fe2dbd8665ULL},
    {"er n=200 deg=1.000000 conn=0 /u64", 0x4a926c88ae1c129fULL},
    {"er n=200 deg=1.000000 conn=0 /f64", 0x339466710818013eULL},
    {"er n=200 deg=1.000000 conn=0 /u64 components", 0x0f74a6d63089c243ULL},
    {"er n=200 deg=1.000000 conn=0 /f64 components", 0x0f74a6d63089c243ULL},
    {"er n=200 deg=4.000000 conn=0 /u64", 0x6073556ef3c6efe9ULL},
    {"er n=200 deg=4.000000 conn=0 /f64", 0xc3b5ff3afaec8847ULL},
    {"er n=200 deg=4.000000 conn=0 /u64 components", 0x5b815b46f9ea1861ULL},
    {"er n=200 deg=4.000000 conn=0 /f64 components", 0x5b815b46f9ea1861ULL},
    {"er n=200 deg=0.000000 conn=1 /u64", 0x81ab4b4cf681b7e2ULL},
    {"er n=200 deg=0.000000 conn=1 /f64", 0x6e760cb7c60aa65bULL},
    {"er n=200 deg=0.000000 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"er n=200 deg=0.000000 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"er n=200 deg=1.000000 conn=1 /u64", 0xedcb336472f798c8ULL},
    {"er n=200 deg=1.000000 conn=1 /f64", 0x3cee08a0dbce4ee7ULL},
    {"er n=200 deg=1.000000 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"er n=200 deg=1.000000 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"er n=200 deg=4.000000 conn=1 /u64", 0xbfabbba33267fbebULL},
    {"er n=200 deg=4.000000 conn=1 /f64", 0xde81322ab7f4dc4cULL},
    {"er n=200 deg=4.000000 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"er n=200 deg=4.000000 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"er n=1000 deg=0.000000 conn=0 /u64", 0xf92ebc7e57b3da5cULL},
    {"er n=1000 deg=0.000000 conn=0 /f64", 0xf92ebc7e57b3da5cULL},
    {"er n=1000 deg=0.000000 conn=0 /u64 components", 0xe54d5ab0da4004a5ULL},
    {"er n=1000 deg=0.000000 conn=0 /f64 components", 0xe54d5ab0da4004a5ULL},
    {"er n=1000 deg=1.000000 conn=0 /u64", 0x723cc7f8a3f47e2fULL},
    {"er n=1000 deg=1.000000 conn=0 /f64", 0x1582088ba1e8f855ULL},
    {"er n=1000 deg=1.000000 conn=0 /u64 components", 0x29765dc1f7049a97ULL},
    {"er n=1000 deg=1.000000 conn=0 /f64 components", 0x29765dc1f7049a97ULL},
    {"er n=1000 deg=4.000000 conn=0 /u64", 0x68b4fbd3eaf7cba0ULL},
    {"er n=1000 deg=4.000000 conn=0 /f64", 0xc73398b7e5fce481ULL},
    {"er n=1000 deg=4.000000 conn=0 /u64 components", 0x98cc4afaffb2295dULL},
    {"er n=1000 deg=4.000000 conn=0 /f64 components", 0x98cc4afaffb2295dULL},
    {"er n=1000 deg=0.000000 conn=1 /u64", 0xa19f43edd4d73077ULL},
    {"er n=1000 deg=0.000000 conn=1 /f64", 0x0a99b9549f2648f8ULL},
    {"er n=1000 deg=0.000000 conn=1 /u64 components", 0x3440d3453fc520dcULL},
    {"er n=1000 deg=0.000000 conn=1 /f64 components", 0x3440d3453fc520dcULL},
    {"er n=1000 deg=1.000000 conn=1 /u64", 0x839e3c602acca120ULL},
    {"er n=1000 deg=1.000000 conn=1 /f64", 0x67c33b65d9ee38eeULL},
    {"er n=1000 deg=1.000000 conn=1 /u64 components", 0x3440d3453fc520dcULL},
    {"er n=1000 deg=1.000000 conn=1 /f64 components", 0x3440d3453fc520dcULL},
    {"er n=1000 deg=4.000000 conn=1 /u64", 0xb63d2ab1d54b3903ULL},
    {"er n=1000 deg=4.000000 conn=1 /f64", 0x05352247032f707bULL},
    {"er n=1000 deg=4.000000 conn=1 /u64 components", 0x3440d3453fc520dcULL},
    {"er n=1000 deg=4.000000 conn=1 /f64 components", 0x3440d3453fc520dcULL},
    {"er n=500 exponential conn=1 /u64", 0x890207fa2390ad07ULL},
    {"er n=500 exponential conn=1 /f64", 0xa415dbcdc1a721f6ULL},
    {"er n=500 exponential conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 exponential conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 log_uniform conn=1 /u64", 0xf541d2e94de34967ULL},
    {"er n=500 log_uniform conn=1 /f64", 0xba2024607474a708ULL},
    {"er n=500 log_uniform conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 log_uniform conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 zero_mixed conn=1 /u64", 0xf3e47b72a6dbc4feULL},
    {"er n=500 zero_mixed conn=1 /f64", 0x546ed78928cabbf9ULL},
    {"er n=500 zero_mixed conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 zero_mixed conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 wide_range conn=1 /u64", 0x470f8c4f1194796dULL},
    {"er n=500 wide_range conn=1 /f64", 0x259b231c2deb40b7ULL},
    {"er n=500 wide_range conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 wide_range conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 uniform_real conn=1 /u64", 0x8ac5ef051028009bULL},
    {"er n=500 uniform_real conn=1 /f64", 0x0c72e31bf30a8f25ULL},
    {"er n=500 uniform_real conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"er n=500 uniform_real conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"regular n=4 d=3 multi=0 conn=0 /u64", 0xf02a66dc3c051807ULL},
    {"regular n=4 d=3 multi=0 conn=0 /f64", 0xc23f7a74efef262cULL},
    {"regular n=4 d=3 multi=0 conn=0 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=0 conn=0 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=0 conn=1 /u64", 0xf02a66dc3c051807ULL},
    {"regular n=4 d=3 multi=0 conn=1 /f64", 0xc23f7a74efef262cULL},
    {"regular n=4 d=3 multi=0 conn=1 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=0 conn=1 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=1 conn=0 /u64", 0xa9621928c49acf11ULL},
    {"regular n=4 d=3 multi=1 conn=0 /f64", 0xfa4ba77cc653c28fULL},
    {"regular n=4 d=3 multi=1 conn=0 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=1 conn=0 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=1 conn=1 /u64", 0xa9621928c49acf11ULL},
    {"regular n=4 d=3 multi=1 conn=1 /f64", 0xfa4ba77cc653c28fULL},
    {"regular n=4 d=3 multi=1 conn=1 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=1 conn=1 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=2 conn=0 /u64", 0x2427a7618d356247ULL},
    {"regular n=4 d=3 multi=2 conn=0 /f64", 0x4a1dbfbcdcb0b45cULL},
    {"regular n=4 d=3 multi=2 conn=0 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=2 conn=0 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=2 conn=1 /u64", 0x2427a7618d356247ULL},
    {"regular n=4 d=3 multi=2 conn=1 /f64", 0x4a1dbfbcdcb0b45cULL},
    {"regular n=4 d=3 multi=2 conn=1 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=4 d=3 multi=2 conn=1 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"regular n=10 d=3 multi=0 conn=0 /u64", 0x8a549109cfff6a8fULL},
    {"regular n=10 d=3 multi=0 conn=0 /f64", 0x5a2d76721048d53bULL},
    {"regular n=10 d=3 multi=0 conn=0 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=0 conn=0 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=0 conn=1 /u64", 0x8a549109cfff6a8fULL},
    {"regular n=10 d=3 multi=0 conn=1 /f64", 0x5a2d76721048d53bULL},
    {"regular n=10 d=3 multi=0 conn=1 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=0 conn=1 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=1 conn=0 /u64", 0x170ec258a7e0d7f7ULL},
    {"regular n=10 d=3 multi=1 conn=0 /f64", 0x7892ea859e2a9657ULL},
    {"regular n=10 d=3 multi=1 conn=0 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=1 conn=0 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=1 conn=1 /u64", 0x170ec258a7e0d7f7ULL},
    {"regular n=10 d=3 multi=1 conn=1 /f64", 0x7892ea859e2a9657ULL},
    {"regular n=10 d=3 multi=1 conn=1 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=1 conn=1 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=2 conn=0 /u64", 0xfe43b339b9d0c94fULL},
    {"regular n=10 d=3 multi=2 conn=0 /f64", 0x6bb4c35b9b4ca29bULL},
    {"regular n=10 d=3 multi=2 conn=0 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=2 conn=0 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=2 conn=1 /u64", 0xfe43b339b9d0c94fULL},
    {"regular n=10 d=3 multi=2 conn=1 /f64", 0x6bb4c35b9b4ca29bULL},
    {"regular n=10 d=3 multi=2 conn=1 /u64 components", 0x1ad12c733512a24fULL},
    {"regular n=10 d=3 multi=2 conn=1 /f64 components", 0x1ad12c733512a24fULL},
    {"regular n=200 d=3 multi=0 conn=0 /u64", 0x0d66964b59302caeULL},
    {"regular n=200 d=3 multi=0 conn=0 /f64", 0x6de721ba5566fbcdULL},
    {"regular n=200 d=3 multi=0 conn=0 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=0 conn=0 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=0 conn=1 /u64", 0x0d66964b59302caeULL},
    {"regular n=200 d=3 multi=0 conn=1 /f64", 0x6de721ba5566fbcdULL},
    {"regular n=200 d=3 multi=0 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=0 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=1 conn=0 /u64", 0x1177d4337c685e33ULL},
    {"regular n=200 d=3 multi=1 conn=0 /f64", 0x67cd5fc760f2e59cULL},
    {"regular n=200 d=3 multi=1 conn=0 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=1 conn=0 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=1 conn=1 /u64", 0x1177d4337c685e33ULL},
    {"regular n=200 d=3 multi=1 conn=1 /f64", 0x67cd5fc760f2e59cULL},
    {"regular n=200 d=3 multi=1 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=1 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=2 conn=0 /u64", 0x1a441889dbabe0ceULL},
    {"regular n=200 d=3 multi=2 conn=0 /f64", 0x4a1ed7a802410441ULL},
    {"regular n=200 d=3 multi=2 conn=0 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=2 conn=0 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=2 conn=1 /u64", 0x1a441889dbabe0ceULL},
    {"regular n=200 d=3 multi=2 conn=1 /f64", 0x4a1ed7a802410441ULL},
    {"regular n=200 d=3 multi=2 conn=1 /u64 components", 0xa92f48f0778f388dULL},
    {"regular n=200 d=3 multi=2 conn=1 /f64 components", 0xa92f48f0778f388dULL},
    {"regular n=300 d=1 drop conn=1 /u64", 0x17c5ffe7fe7b9d3fULL},
    {"regular n=300 d=1 drop conn=1 /f64", 0x3042b685dc6b1c98ULL},
    {"regular n=300 d=1 drop conn=1 /u64 components", 0x3b8ae904d210c9a2ULL},
    {"regular n=300 d=1 drop conn=1 /f64 components", 0x3b8ae904d210c9a2ULL},
    {"ba n=3 attach=2 /u64", 0xa63342d17150c5ddULL},
    {"ba n=3 attach=2 /f64", 0x246e5409761800c8ULL},
    {"ba n=3 attach=2 /u64 components", 0x25a3640df7015cc6ULL},
    {"ba n=3 attach=2 /f64 components", 0x25a3640df7015cc6ULL},
    {"ba n=50 attach=2 /u64", 0xfbc54cbfa9596b0dULL},
    {"ba n=50 attach=2 /f64", 0x30252f8213c6e599ULL},
    {"ba n=50 attach=2 /u64 components", 0x0e7d1eb77f939037ULL},
    {"ba n=50 attach=2 /f64 components", 0x0e7d1eb77f939037ULL},
    {"ba n=500 attach=2 /u64", 0x2f1102367016f00fULL},
    {"ba n=500 attach=2 /f64", 0x60e393552fa6630eULL},
    {"ba n=500 attach=2 /u64 components", 0x7d01a19212f47f6aULL},
    {"ba n=500 attach=2 /f64 components", 0x7d01a19212f47f6aULL},
    {"grid 1x1 /u64", 0x392209f14dea4c24ULL},
    {"grid 1x1 /f64", 0x392209f14dea4c24ULL},
    {"grid 1x1 /u64 components", 0x372cac50b2ba8344ULL},
    {"grid 1x1 /f64 components", 0x372cac50b2ba8344ULL},
    {"grid 3x7 /u64", 0x0b5867f432fc3605ULL},
    {"grid 3x7 /f64", 0x6042162a513630a4ULL},
    {"grid 3x7 /u64 components", 0x2b71a1e264bb8050ULL},
    {"grid 3x7 /f64 components", 0x2b71a1e264bb8050ULL},
    {"grid 20x20 /u64", 0x9be3407d34cec98fULL},
    {"grid 20x20 /f64", 0xc73427e3692e8724ULL},
    {"grid 20x20 /u64 components", 0x46421c849998d64eULL},
    {"grid 20x20 /f64 components", 0x46421c849998d64eULL},
    {"road 1x1 keep=0.900000 /u64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.900000 /f64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.900000 /u64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.900000 /f64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.500000 /u64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.500000 /f64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.500000 /u64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.500000 /f64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.100000 /u64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.100000 /f64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.100000 /u64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.100000 /f64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.000000 /u64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.000000 /f64", 0x392209f14dea4c24ULL},
    {"road 1x1 keep=0.000000 /u64 components", 0x372cac50b2ba8344ULL},
    {"road 1x1 keep=0.000000 /f64 components", 0x372cac50b2ba8344ULL},
    {"road 2x2 keep=0.900000 /u64", 0x92f0bee710c3e59aULL},
    {"road 2x2 keep=0.900000 /f64", 0xacb5936ca00361b7ULL},
    {"road 2x2 keep=0.900000 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.900000 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.500000 /u64", 0x97ea3ec7bf355dbbULL},
    {"road 2x2 keep=0.500000 /f64", 0x25c8487fd2d9a2e1ULL},
    {"road 2x2 keep=0.500000 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.500000 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.100000 /u64", 0x97ea3ec7bf355dbbULL},
    {"road 2x2 keep=0.100000 /f64", 0x25c8487fd2d9a2e1ULL},
    {"road 2x2 keep=0.100000 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.100000 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.000000 /u64", 0x6c45ac4db87d2cebULL},
    {"road 2x2 keep=0.000000 /f64", 0xbc36c62933736eadULL},
    {"road 2x2 keep=0.000000 /u64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 2x2 keep=0.000000 /f64 components", 0x533b4e0c2e4a73c1ULL},
    {"road 5x9 keep=0.900000 /u64", 0x71c3f794349d25eeULL},
    {"road 5x9 keep=0.900000 /f64", 0xf526e0a6a128be8bULL},
    {"road 5x9 keep=0.900000 /u64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.900000 /f64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.500000 /u64", 0x440441f2cea9a198ULL},
    {"road 5x9 keep=0.500000 /f64", 0x7a2cfa26f1a4a169ULL},
    {"road 5x9 keep=0.500000 /u64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.500000 /f64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.100000 /u64", 0x14d8800d77f0b63aULL},
    {"road 5x9 keep=0.100000 /f64", 0xd8273e080e1c8a98ULL},
    {"road 5x9 keep=0.100000 /u64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.100000 /f64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.000000 /u64", 0x8e110e07837603e0ULL},
    {"road 5x9 keep=0.000000 /f64", 0x8c484ae969a58f69ULL},
    {"road 5x9 keep=0.000000 /u64 components", 0x69b94cb7272a2228ULL},
    {"road 5x9 keep=0.000000 /f64 components", 0x69b94cb7272a2228ULL},
    {"road 30x30 keep=0.900000 /u64", 0x358c69c9aac155eeULL},
    {"road 30x30 keep=0.900000 /f64", 0x4ecdc0186c7987f3ULL},
    {"road 30x30 keep=0.900000 /u64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.900000 /f64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.500000 /u64", 0x6da8cd525b8c5d4cULL},
    {"road 30x30 keep=0.500000 /f64", 0xc9ab37591d4249a9ULL},
    {"road 30x30 keep=0.500000 /u64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.500000 /f64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.100000 /u64", 0x1977722ebe7e330bULL},
    {"road 30x30 keep=0.100000 /f64", 0x59182537750f7592ULL},
    {"road 30x30 keep=0.100000 /u64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.100000 /f64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.000000 /u64", 0xdcde537144f95e5bULL},
    {"road 30x30 keep=0.000000 /f64", 0x793a11e9b2e97377ULL},
    {"road 30x30 keep=0.000000 /u64 components", 0xa8719afc0a9af3b0ULL},
    {"road 30x30 keep=0.000000 /f64 components", 0xa8719afc0a9af3b0ULL},
    {"road 1x40 keep=0.900000 /u64", 0x4faedc14f87c1dc6ULL},
    {"road 1x40 keep=0.900000 /f64", 0x04f5fcc85f5cf6baULL},
    {"road 1x40 keep=0.900000 /u64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.900000 /f64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.500000 /u64", 0x096c48e7e3a3043fULL},
    {"road 1x40 keep=0.500000 /f64", 0xd41e553ae1a46455ULL},
    {"road 1x40 keep=0.500000 /u64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.500000 /f64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.100000 /u64", 0xb3d0a3d595e86e8bULL},
    {"road 1x40 keep=0.100000 /f64", 0x272b96dc0d1ce4b9ULL},
    {"road 1x40 keep=0.100000 /u64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.100000 /f64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.000000 /u64", 0x879ca8782513de37ULL},
    {"road 1x40 keep=0.000000 /f64", 0xed7c23eacfbad429ULL},
    {"road 1x40 keep=0.000000 /u64 components", 0x5932ca678aa170adULL},
    {"road 1x40 keep=0.000000 /f64 components", 0x5932ca678aa170adULL},
    {"tree n=1 /u64", 0x392209f14dea4c24ULL},
    {"tree n=1 /f64", 0x392209f14dea4c24ULL},
    {"tree n=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"tree n=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"tree n=2 /u64", 0x4bc8f65f1afd40c3ULL},
    {"tree n=2 /f64", 0xc7ea328346e82f17ULL},
    {"tree n=2 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"tree n=2 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"tree n=50 /u64", 0x495139abe3c503a8ULL},
    {"tree n=50 /f64", 0x8e3c68c4d1e30a15ULL},
    {"tree n=50 /u64 components", 0x0e7d1eb77f939037ULL},
    {"tree n=50 /f64 components", 0x0e7d1eb77f939037ULL},
    {"tree n=1000 /u64", 0xd80cc7b11b85026aULL},
    {"tree n=1000 /f64", 0xfec9cb0a1bf19a72ULL},
    {"tree n=1000 /u64 components", 0x3440d3453fc520dcULL},
    {"tree n=1000 /f64 components", 0x3440d3453fc520dcULL},
    {"chung_lu n=100 conn=0 /u64", 0xb10a796ab6ca701dULL},
    {"chung_lu n=100 conn=0 /f64", 0xbfcb452677039a4aULL},
    {"chung_lu n=100 conn=0 /u64 components", 0x3305079bc8b58dddULL},
    {"chung_lu n=100 conn=0 /f64 components", 0x3305079bc8b58dddULL},
    {"chung_lu n=100 conn=1 /u64", 0xca2a42250c4f289cULL},
    {"chung_lu n=100 conn=1 /f64", 0xddadf89269bedc4fULL},
    {"chung_lu n=100 conn=1 /u64 components", 0x8c2f5efa37b92be1ULL},
    {"chung_lu n=100 conn=1 /f64 components", 0x8c2f5efa37b92be1ULL},
    {"chung_lu n=1000 conn=0 /u64", 0x70c15195ad7dc3d8ULL},
    {"chung_lu n=1000 conn=0 /f64", 0x5b3392233f56ca97ULL},
    {"chung_lu n=1000 conn=0 /u64 components", 0xf78f5bd47e7183b3ULL},
    {"chung_lu n=1000 conn=0 /f64 components", 0xf78f5bd47e7183b3ULL},
    {"chung_lu n=1000 conn=1 /u64", 0x358025a4eae2d8a9ULL},
    {"chung_lu n=1000 conn=1 /f64", 0xc63afee7a8267371ULL},
    {"chung_lu n=1000 conn=1 /u64 components", 0x3440d3453fc520dcULL},
    {"chung_lu n=1000 conn=1 /f64 components", 0x3440d3453fc520dcULL},
    {"adversarial_decrease_key n=1 /u64", 0x392209f14dea4c24ULL},
    {"adversarial_decrease_key n=1 /f64", 0x392209f14dea4c24ULL},
    {"adversarial_decrease_key n=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"adversarial_decrease_key n=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"adversarial_long_path n=1 /u64", 0x392209f14dea4c24ULL},
    {"adversarial_long_path n=1 /f64", 0x392209f14dea4c24ULL},
    {"adversarial_long_path n=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"adversarial_long_path n=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"adversarial_decrease_key n=2 /u64", 0xa8b94b7a3bcb1f26ULL},
    {"adversarial_decrease_key n=2 /f64", 0xc4f117834461aa16ULL},
    {"adversarial_decrease_key n=2 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"adversarial_decrease_key n=2 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"adversarial_long_path n=2 /u64", 0xa8b94b7a3bcb1f26ULL},
    {"adversarial_long_path n=2 /f64", 0xc4f117834461aa16ULL},
    {"adversarial_long_path n=2 /u64 components", 0x5f3187dbaf3d3c47ULL},
    {"adversarial_long_path n=2 /f64 components", 0x5f3187dbaf3d3c47ULL},
    {"adversarial_decrease_key n=10 /u64", 0x4f5be87df3702082ULL},
    {"adversarial_decrease_key n=10 /f64", 0x50bf49b5cdb240e6ULL},
    {"adversarial_decrease_key n=10 /u64 components", 0x1ad12c733512a24fULL},
    {"adversarial_decrease_key n=10 /f64 components", 0x1ad12c733512a24fULL},
    {"adversarial_long_path n=10 /u64", 0x93b559b429576d36ULL},
    {"adversarial_long_path n=10 /f64", 0x418bf40bda6e4a76ULL},
    {"adversarial_long_path n=10 /u64 components", 0x1ad12c733512a24fULL},
    {"adversarial_long_path n=10 /f64 components", 0x1ad12c733512a24fULL},
    {"adversarial_decrease_key n=64 /u64", 0x454def80dd0f9f96ULL},
    {"adversarial_decrease_key n=64 /f64", 0xb80483679ebbb512ULL},
    {"adversarial_decrease_key n=64 /u64 components", 0xaf46597ed414b305ULL},
    {"adversarial_decrease_key n=64 /f64 components", 0xaf46597ed414b305ULL},
    {"adversarial_long_path n=64 /u64", 0x94cf21b73f953726ULL},
    {"adversarial_long_path n=64 /f64", 0x8f84fda2b563d166ULL},
    {"adversarial_long_path n=64 /u64 components", 0xaf46597ed414b305ULL},
    {"adversarial_long_path n=64 /f64 components", 0xaf46597ed414b305ULL},
    {"geometric n=1 weights=0 conn=0 /u64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=0 conn=0 /f64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=0 conn=0 /u64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=0 conn=0 /f64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=0 conn=1 /u64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=0 conn=1 /f64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=0 conn=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=0 conn=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=1 conn=0 /u64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=1 conn=0 /f64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=1 conn=0 /u64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=1 conn=0 /f64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=1 conn=1 /u64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=1 conn=1 /f64", 0x392209f14dea4c24ULL},
    {"geometric n=1 weights=1 conn=1 /u64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=1 weights=1 conn=1 /f64 components", 0x372cac50b2ba8344ULL},
    {"geometric n=50 weights=0 conn=0 /u64", 0x204478d6cdbed91eULL},
    {"geometric n=50 weights=0 conn=0 /f64", 0x4299832f5a3d1813ULL},
    {"geometric n=50 weights=0 conn=0 /u64 components", 0xf07d351ac9066597ULL},
    {"geometric n=50 weights=0 conn=0 /f64 components", 0xf07d351ac9066597ULL},
    {"geometric n=50 weights=0 conn=1 /u64", 0xc662cf851461bc91ULL},
    {"geometric n=50 weights=0 conn=1 /f64", 0xb9a927690e193be6ULL},
    {"geometric n=50 weights=0 conn=1 /u64 components", 0x0e7d1eb77f939037ULL},
    {"geometric n=50 weights=0 conn=1 /f64 components", 0x0e7d1eb77f939037ULL},
    {"geometric n=50 weights=1 conn=0 /u64", 0xea9dcd2e51da0df2ULL},
    {"geometric n=50 weights=1 conn=0 /f64", 0xeaa277002907c52eULL},
    {"geometric n=50 weights=1 conn=0 /u64 components", 0xf07d351ac9066597ULL},
    {"geometric n=50 weights=1 conn=0 /f64 components", 0xf07d351ac9066597ULL},
    {"geometric n=50 weights=1 conn=1 /u64", 0x29bbb378ef9b3a0bULL},
    {"geometric n=50 weights=1 conn=1 /f64", 0xc436566b562a2b43ULL},
    {"geometric n=50 weights=1 conn=1 /u64 components", 0x0e7d1eb77f939037ULL},
    {"geometric n=50 weights=1 conn=1 /f64 components", 0x0e7d1eb77f939037ULL},
    {"geometric n=500 weights=0 conn=0 /u64", 0x264d6d088e48f21bULL},
    {"geometric n=500 weights=0 conn=0 /f64", 0xe86b63731d82b59eULL},
    {"geometric n=500 weights=0 conn=0 /u64 components", 0xf904d1c619906993ULL},
    {"geometric n=500 weights=0 conn=0 /f64 components", 0xf904d1c619906993ULL},
    {"geometric n=500 weights=0 conn=1 /u64", 0xa4a36de5a4f77b15ULL},
    {"geometric n=500 weights=0 conn=1 /f64", 0x00cedc70b6dfae40ULL},
    {"geometric n=500 weights=0 conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"geometric n=500 weights=0 conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"geometric n=500 weights=1 conn=0 /u64", 0xe9016bf5f42c7eafULL},
    {"geometric n=500 weights=1 conn=0 /f64", 0x8ab3ac5209accfd9ULL},
    {"geometric n=500 weights=1 conn=0 /u64 components", 0xf904d1c619906993ULL},
    {"geometric n=500 weights=1 conn=0 /f64 components", 0xf904d1c619906993ULL},
    {"geometric n=500 weights=1 conn=1 /u64", 0x21ddd9199053e8adULL},
    {"geometric n=500 weights=1 conn=1 /f64", 0xedffb5c3a4b4cb01ULL},
    {"geometric n=500 weights=1 conn=1 /u64 components", 0x7d01a19212f47f6aULL},
    {"geometric n=500 weights=1 conn=1 /f64 components", 0x7d01a19212f47f6aULL},
    {"gnm n=0 /u64", 0x88201fb960ff6465ULL},
    {"gnm n=0 /f64", 0x88201fb960ff6465ULL},
    {"gnm n=5 /u64", 0x9d2023ce48a6c5bdULL},
    {"gnm n=5 /f64", 0x4b4d0556be452261ULL},
    {"gnm n=300 /u64", 0xbce7e335b6645589ULL},
    {"gnm n=300 /f64", 0xfda5a598b120e663ULL},
    {"direct components before /u64", 0x2af78679b053103dULL},
    {"direct connect_components /u64", 0xc50f78a163590d61ULL},
    {"direct components after /u64", 0xf3290988b10d2db9ULL},
    {"direct connect_components n=1 /f64", 0x392209f14dea4c24ULL},
    {"direct connect_components n=0 /f64", 0x88201fb960ff6465ULL},
    {"direct components n=0 /f64", 0x88201fb960ff6465ULL},
};

}  // namespace

TEST_CASE("generator, connect_components and connected_components output is pinned for fixed seeds",
          "[regression][generators]") {
  const auto actual = actual_digests();
  for (std::size_t i = 0; i < actual.size(); ++i) {
    // The message is the table row for this output, so a failure shows what changed.
    INFO("{\"" << actual[i].first << "\", " << hex(actual[i].second) << "ULL},");
    CHECK(i < kExpected.size());
    if (i >= kExpected.size()) continue;
    CHECK(actual[i].first == kExpected[i].first);
    CHECK(actual[i].second == kExpected[i].second);
  }
  CHECK(actual.size() == kExpected.size());
}
