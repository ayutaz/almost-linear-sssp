// Reproducible graph generators (see include/sssp/generators.hpp for the contract).

// Floating-point contraction (FMA) changes rounding and therefore the generated graphs.
// GCC in ISO mode (-std=c++20) already uses -ffp-contract=off; clang needs this pragma.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#include "sssp/generators.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sssp::gen {

namespace {
constexpr double kLn2 = 0.693147180559945309417232121458176568;
constexpr double kLn2Hi = 0x1.62e42feep-1;        // fdlibm ln2_hi
constexpr double kLn2Lo = 0x1.a39ef35793c76p-33;  // fdlibm ln2_lo
constexpr double kSqrtHalf = 0.707106781186547524400844362104849039;
}  // namespace

double det_log(double x) noexcept {
  // x = m * 2^e with m in [sqrt(1/2), sqrt(2)); frexp is exact.
  int e = 0;
  double m = std::frexp(x, &e);  // m in [0.5, 1)
  if (m < kSqrtHalf) {
    m *= 2.0;
    --e;
  }
  // log m = 2 atanh(s), s = (m-1)/(m+1), |s| <= 0.1716; the series converges geometrically
  // with ratio s^2 <= 0.0295, so 14 terms reach well below 2^-53.
  const double s = (m - 1.0) / (m + 1.0);
  const double s2 = s * s;
  double term = s;
  double sum = 0.0;
  for (int k = 0; k < 14; ++k) {
    sum += term / static_cast<double>(2 * k + 1);
    term *= s2;
  }
  return static_cast<double>(e) * kLn2 + 2.0 * sum;
}

double det_exp(double x) noexcept {
  // x = k ln2 + r with |r| <= ln2/2; exp(r) by Taylor series; ldexp is exact.
  const double kd = std::floor(x / kLn2 + 0.5);
  // Cody–Waite reduction: ln2 = kLn2Hi + kLn2Lo where kLn2Hi has 32 trailing zero bits, so
  // kd * kLn2Hi is exact for |kd| < 2^20 and the reduction loses no precision.
  const double r = (x - kd * kLn2Hi) - kd * kLn2Lo;
  double term = 1.0;
  double sum = 1.0;
  for (int i = 1; i < 24; ++i) {
    term *= r / static_cast<double>(i);
    sum += term;
  }
  return std::ldexp(sum, static_cast<int>(kd));
}

// ------------------------------------------------------------------ weights

namespace {

constexpr double kMaxExactInt = 9007199254740992.0;  // 2^53

void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("invalid weight spec: ") + what);
}

/// Rounds a non-negative real draw to the weight type (nearest integer, ties up, for integers).
template <class W>
W from_real(double x) {
  if constexpr (std::is_integral_v<W>) {
    return static_cast<W>(std::floor(x + 0.5));
  } else {
    return static_cast<W>(x);
  }
}

std::uint64_t int_lo(const WeightSpec& s) {
  return static_cast<std::uint64_t>(std::ceil(s.lo));
}
std::uint64_t int_hi(const WeightSpec& s) {
  return static_cast<std::uint64_t>(std::floor(s.hi));
}

}  // namespace

void validate_weight_spec(const WeightSpec& s) {
  require(std::isfinite(s.lo) && std::isfinite(s.hi), "lo and hi must be finite");
  require(s.lo >= 0.0, "lo must be >= 0");
  require(s.hi <= kMaxExactInt, "hi must be <= 2^53");
  switch (s.dist) {
    case WeightDist::kUniformInt:
      require(s.lo <= s.hi && int_lo(s) <= int_hi(s), "[lo, hi] must contain an integer");
      break;
    case WeightDist::kUniformReal:
      require(s.lo <= s.hi, "lo must be <= hi");
      break;
    case WeightDist::kLogUniform:
      require(s.lo > 0.0 && s.lo <= s.hi, "log_uniform needs 0 < lo <= hi");
      break;
    case WeightDist::kExponential:
      require(std::isfinite(s.mean) && s.mean > 0.0 && s.mean <= 0x1p47,
              "exponential needs 0 < mean <= 2^47");
      break;
    case WeightDist::kConstant:
      break;
    case WeightDist::kZeroMixed:
      require(s.zero_probability >= 0.0 && s.zero_probability <= 1.0,
              "zero_probability must be in [0, 1]");
      require(s.lo <= s.hi && int_lo(s) <= int_hi(s), "[lo, hi] must contain an integer");
      break;
    case WeightDist::kWideRange:
      require(s.wide_log2 >= 0 && s.wide_log2 <= 62, "wide_log2 must be in [0, 62]");
      break;
  }
}

template <class W>
W sample_weight(const WeightSpec& s, Xoshiro256ss& rng) {
  switch (s.dist) {
    case WeightDist::kUniformInt:
      return static_cast<W>(rng.uniform_range(int_lo(s), int_hi(s)));
    case WeightDist::kUniformReal:
      return from_real<W>(s.lo + (s.hi - s.lo) * rng.uniform01());
    case WeightDist::kLogUniform: {
      const double a = det_log(s.lo);
      const double b = det_log(s.hi);
      const double x = det_exp(a + (b - a) * rng.uniform01());
      return from_real<W>(std::clamp(x, s.lo, s.hi));
    }
    case WeightDist::kExponential:
      // 1 - U is in (0, 1], so the logarithm is finite.
      return from_real<W>(-s.mean * det_log(1.0 - rng.uniform01()));
    case WeightDist::kConstant:
      return from_real<W>(s.lo);
    case WeightDist::kZeroMixed:
      if (rng.bernoulli(s.zero_probability)) return W{0};
      if constexpr (std::is_integral_v<W>) {
        return static_cast<W>(rng.uniform_range(int_lo(s), int_hi(s)));
      } else {
        return s.lo + (s.hi - s.lo) * rng.uniform01();
      }
    case WeightDist::kWideRange: {
      const auto k = static_cast<std::uint64_t>(s.wide_log2);
      const auto e = rng.uniform_range(0, k);
      if (e == k) return static_cast<W>(std::uint64_t{1} << k);
      if constexpr (std::is_integral_v<W>) {
        const std::uint64_t base = std::uint64_t{1} << e;
        return static_cast<W>(rng.uniform_range(base, 2 * base - 1));
      } else {
        return std::ldexp(1.0 + rng.uniform01(), static_cast<int>(e));
      }
    }
  }
  return W{0};  // unreachable
}

template std::uint64_t sample_weight<std::uint64_t>(const WeightSpec&, Xoshiro256ss&);
template double sample_weight<double>(const WeightSpec&, Xoshiro256ss&);

namespace {
constexpr std::array<std::pair<WeightDist, std::string_view>, 7> kWeightNames{{
    {WeightDist::kUniformInt, "uniform_int"},
    {WeightDist::kUniformReal, "uniform_real"},
    {WeightDist::kLogUniform, "log_uniform"},
    {WeightDist::kExponential, "exponential"},
    {WeightDist::kConstant, "constant"},
    {WeightDist::kZeroMixed, "zero_mixed"},
    {WeightDist::kWideRange, "wide_range"},
}};
}  // namespace

std::string_view to_string(WeightDist d) noexcept {
  for (const auto& [k, name] : kWeightNames) {
    if (k == d) return name;
  }
  return "unknown";
}

WeightDist parse_weight_dist(std::string_view name) {
  for (const auto& [k, n] : kWeightNames) {
    if (n == name) return k;
  }
  throw std::invalid_argument("unknown weight distribution: " + std::string(name));
}

// ------------------------------------------------------------------ graph generators

namespace {

using Pair = std::pair<Vertex, Vertex>;

enum Stream : std::uint64_t { kTopology = 0, kWeights = 1, kConnect = 2 };

void check(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument(what);
}

void check_n(std::size_t n) {
  check(n < static_cast<std::size_t>(kNoVertex), "too many vertices");
}

/// Order-independent key of an unordered vertex pair (n < 2^32 so it fits in 64 bits).
std::uint64_t pair_key(Vertex a, Vertex b) {
  if (a > b) std::swap(a, b);
  return (static_cast<std::uint64_t>(a) << 32) | b;
}

/// Attaches weights drawn from the weight sub-stream, in edge order.
template <class W>
EdgeList<W> with_weights(std::size_t n, const std::vector<Pair>& pairs, const GenOptions& opt) {
  Xoshiro256ss rng(substream_seed(opt.seed, kWeights));
  EdgeList<W> g{n, {}};
  g.edges.reserve(pairs.size());
  for (const auto& [u, v] : pairs) g.edges.push_back({u, v, sample_weight<W>(opt.weights, rng)});
  return g;
}

template <class W>
EdgeList<W> finish(std::size_t n, const std::vector<Pair>& pairs, const GenOptions& opt) {
  auto g = with_weights<W>(n, pairs, opt);
  if (opt.connected) connect_components(g, opt.weights, substream_seed(opt.seed, kConnect));
  return g;
}

class UnionFind {
 public:
  explicit UnionFind(std::size_t n) : parent_(n) {
    for (std::size_t i = 0; i < n; ++i) parent_[i] = i;
  }
  std::size_t find(std::size_t x) {
    while (parent_[x] != x) x = parent_[x] = parent_[parent_[x]];
    return x;
  }
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    if (a < b) std::swap(a, b);
    parent_[a] = b;  // smaller index becomes the root: deterministic and order-independent
    return true;
  }

 private:
  std::vector<std::size_t> parent_;
};

std::size_t rounded_edge_count(std::size_t n, double avg_degree) {
  check(std::isfinite(avg_degree) && avg_degree >= 0.0, "avg_degree must be finite and >= 0");
  return static_cast<std::size_t>(std::floor(static_cast<double>(n) * avg_degree / 2.0 + 0.5));
}

/// Pair index t in [0, n(n-1)/2) -> (b, a) with b < a, t = a(a-1)/2 + b. Exact integer result
/// (the double estimate is only a starting point and is corrected).
Pair decode_pair(std::uint64_t t) {
  auto a = static_cast<std::uint64_t>((1.0 + std::sqrt(1.0 + 8.0 * static_cast<double>(t))) / 2.0);
  while (a > 1 && a * (a - 1) / 2 > t) --a;
  while ((a + 1) * a / 2 <= t) ++a;
  const std::uint64_t b = t - a * (a - 1) / 2;
  return {static_cast<Vertex>(b), static_cast<Vertex>(a)};
}

}  // namespace

template <class W>
std::size_t count_components(const EdgeList<W>& g) {
  UnionFind uf(g.num_vertices);
  std::size_t c = g.num_vertices;
  for (const auto& e : g.edges) {
    if (uf.unite(e.u, e.v)) --c;
  }
  return c;
}

template <class W>
void connect_components(EdgeList<W>& g, const WeightSpec& weights, std::uint64_t seed) {
  validate_weight_spec(weights);
  const std::size_t n = g.num_vertices;
  UnionFind uf(n);
  for (const auto& e : g.edges) uf.unite(e.u, e.v);
  // Components in order of their smallest vertex; members in increasing vertex order.
  std::vector<std::size_t> label(n, SIZE_MAX);
  std::vector<std::vector<Vertex>> members;
  for (std::size_t v = 0; v < n; ++v) {
    const auto r = uf.find(v);
    if (label[r] == SIZE_MAX) {
      label[r] = members.size();
      members.emplace_back();
    }
    members[label[r]].push_back(static_cast<Vertex>(v));
  }
  Xoshiro256ss rng(seed);
  Xoshiro256ss wrng(substream_seed(seed, kWeights));
  for (std::size_t i = 1; i < members.size(); ++i) {
    const auto j = static_cast<std::size_t>(rng.uniform_below(i));
    const Vertex a = members[i][rng.uniform_below(members[i].size())];
    const Vertex b = members[j][rng.uniform_below(members[j].size())];
    g.edges.push_back({b, a, sample_weight<W>(weights, wrng)});
  }
}

template <class W>
EdgeList<W> erdos_renyi_gnm(std::size_t n, std::size_t m, const GenOptions& opt) {
  check_n(n);
  validate_weight_spec(opt.weights);
  const std::uint64_t max_pairs = static_cast<std::uint64_t>(n) * (n > 0 ? n - 1 : 0) / 2;
  check(m <= max_pairs, "G(n,m): m exceeds n(n-1)/2");
  // Floyd's algorithm: a uniformly random m-subset of the pair indices in O(m) expected time.
  Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
  std::unordered_set<std::uint64_t> chosen;
  chosen.reserve(m * 2);
  std::vector<Pair> pairs;
  pairs.reserve(m);
  for (std::uint64_t j = max_pairs - m; j < max_pairs; ++j) {
    const std::uint64_t t = rng.uniform_range(0, j);
    const std::uint64_t pick = chosen.insert(t).second ? t : j;
    if (pick == j) chosen.insert(j);
    pairs.push_back(decode_pair(pick));
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> erdos_renyi(std::size_t n, double avg_degree, const GenOptions& opt) {
  return erdos_renyi_gnm<W>(n, rounded_edge_count(n, avg_degree), opt);
}

template <class W>
EdgeList<W> random_regular(std::size_t n, std::size_t d, MultiEdgePolicy policy,
                           const GenOptions& opt) {
  check_n(n);
  validate_weight_spec(opt.weights);
  check((n * d) % 2 == 0, "random_regular: n * d must be even");
  if (policy == MultiEdgePolicy::kSwitchToSimple) {
    check(d == 0 || d < n, "random_regular: a simple d-regular graph needs d < n");
  }
  Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
  std::vector<Vertex> stubs;
  stubs.reserve(n * d);
  for (std::size_t v = 0; v < n; ++v) stubs.insert(stubs.end(), d, static_cast<Vertex>(v));
  for (std::size_t i = stubs.size(); i > 1; --i) {  // Fisher–Yates
    std::swap(stubs[i - 1], stubs[rng.uniform_below(i)]);
  }
  std::vector<Pair> pairs;
  pairs.reserve(stubs.size() / 2);
  for (std::size_t i = 0; i + 1 < stubs.size(); i += 2) pairs.emplace_back(stubs[i], stubs[i + 1]);

  if (policy == MultiEdgePolicy::kDrop) {
    std::unordered_set<std::uint64_t> seen;
    std::vector<Pair> kept;
    for (const auto& [a, b] : pairs) {
      if (a != b && seen.insert(pair_key(a, b)).second) kept.emplace_back(a, b);
    }
    pairs = std::move(kept);
  } else if (policy == MultiEdgePolicy::kSwitchToSimple && !pairs.empty()) {
    // Degree-preserving double-edge swaps until no self-loop or multi-edge remains.
    std::unordered_map<std::uint64_t, std::size_t> mult;
    for (const auto& [a, b] : pairs) ++mult[pair_key(a, b)];
    auto bad = [&](const Pair& p) {
      return p.first == p.second || mult[pair_key(p.first, p.second)] > 1;
    };
    const std::size_t max_rounds = 1000 + 100 * pairs.size();
    std::size_t rounds = 0;
    for (bool dirty = true; dirty;) {
      dirty = false;
      for (std::size_t i = 0; i < pairs.size(); ++i) {
        if (!bad(pairs[i])) continue;
        dirty = true;
        check(++rounds < max_rounds, "random_regular: switching did not converge");
        const auto j = static_cast<std::size_t>(rng.uniform_below(pairs.size()));
        if (j == i) continue;
        auto [a, b] = pairs[i];
        auto [c, e] = pairs[j];
        if (rng.bernoulli(0.5)) std::swap(c, e);
        // proposal: {a,b},{c,e} -> {a,c},{b,e}
        if (a == c || b == e) continue;
        const auto k1 = pair_key(a, c), k2 = pair_key(b, e);
        if (k1 == k2 || mult[k1] > 0 || mult[k2] > 0) continue;
        --mult[pair_key(a, b)];
        --mult[pair_key(c, e)];
        ++mult[k1];
        ++mult[k2];
        pairs[i] = {a, c};
        pairs[j] = {b, e};
      }
    }
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> barabasi_albert(std::size_t n, std::size_t m, const GenOptions& opt) {
  check_n(n);
  validate_weight_spec(opt.weights);
  check(m >= 1 && m < n, "barabasi_albert: need 1 <= m < n");
  Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
  const std::size_t m0 = m + 1;
  std::vector<Pair> pairs;
  std::vector<Vertex> endpoints;  // each vertex appears deg(v) times: sampling ∝ degree
  for (std::size_t i = 0; i < m0; ++i) {
    for (std::size_t j = i + 1; j < m0; ++j) {
      pairs.emplace_back(static_cast<Vertex>(i), static_cast<Vertex>(j));
      endpoints.push_back(static_cast<Vertex>(i));
      endpoints.push_back(static_cast<Vertex>(j));
    }
  }
  std::vector<Vertex> targets;
  for (std::size_t v = m0; v < n; ++v) {
    targets.clear();
    while (targets.size() < m) {
      const Vertex t = endpoints[rng.uniform_below(endpoints.size())];
      if (std::find(targets.begin(), targets.end(), t) == targets.end()) targets.push_back(t);
    }
    for (Vertex t : targets) {
      pairs.emplace_back(t, static_cast<Vertex>(v));
      endpoints.push_back(t);
      endpoints.push_back(static_cast<Vertex>(v));
    }
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> grid_2d(std::size_t rows, std::size_t cols, const GenOptions& opt) {
  check(rows == 0 || cols <= static_cast<std::size_t>(kNoVertex) / rows, "grid too large");
  const std::size_t n = rows * cols;
  check_n(n);
  validate_weight_spec(opt.weights);
  std::vector<Pair> pairs;
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      const auto v = static_cast<Vertex>(r * cols + c);
      if (c + 1 < cols) pairs.emplace_back(v, v + 1);
      if (r + 1 < rows) pairs.emplace_back(v, static_cast<Vertex>(v + cols));
    }
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> road_like(std::size_t rows, std::size_t cols, const RoadLikeOptions& road,
                      const GenOptions& opt) {
  check(rows == 0 || cols <= static_cast<std::size_t>(kNoVertex) / rows, "grid too large");
  check(road.jitter >= 0.0 && road.jitter < 0.5, "road_like: jitter must be in [0, 0.5)");
  check(road.keep_probability >= 0.0 && road.keep_probability <= 1.0,
        "road_like: keep_probability must be in [0, 1]");
  check(road.shortcut_probability >= 0.0 && road.shortcut_probability <= 1.0,
        "road_like: shortcut_probability must be in [0, 1]");
  check(road.shortcut_radius >= 2, "road_like: shortcut_radius must be >= 2");
  check(std::isfinite(road.integer_scale) && road.integer_scale > 0.0,
        "road_like: integer_scale must be > 0");
  const std::size_t n = rows * cols;
  check_n(n);
  Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
  std::vector<double> x(n), y(n);
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      const std::size_t v = r * cols + c;
      y[v] = static_cast<double>(r) + road.jitter * (2.0 * rng.uniform01() - 1.0);
      x[v] = static_cast<double>(c) + road.jitter * (2.0 * rng.uniform01() - 1.0);
    }
  }
  std::vector<Pair> pairs, removed;
  std::unordered_set<std::uint64_t> present;
  auto add = [&](Vertex a, Vertex b) {
    pairs.emplace_back(a, b);
    present.insert(pair_key(a, b));
  };
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      const auto v = static_cast<Vertex>(r * cols + c);
      auto street = [&](Vertex w) {
        if (rng.bernoulli(road.keep_probability)) {
          add(v, w);
        } else {
          removed.emplace_back(v, w);
        }
      };
      if (c + 1 < cols) street(v + 1);
      if (r + 1 < rows) street(static_cast<Vertex>(v + cols));
    }
  }
  // Local shortcuts: a random target within Chebyshev radius R that is not a grid neighbour.
  const auto R = static_cast<std::uint64_t>(road.shortcut_radius);
  for (std::size_t v = 0; v < n; ++v) {
    if (!rng.bernoulli(road.shortcut_probability)) continue;
    const auto dr = static_cast<long long>(rng.uniform_range(0, 2 * R)) - static_cast<long long>(R);
    const auto dc = static_cast<long long>(rng.uniform_range(0, 2 * R)) - static_cast<long long>(R);
    if (std::llabs(dr) + std::llabs(dc) < 2) continue;
    const auto r2 = static_cast<long long>(v / cols) + dr;
    const auto c2 = static_cast<long long>(v % cols) + dc;
    if (r2 < 0 || c2 < 0 || r2 >= static_cast<long long>(rows) ||
        c2 >= static_cast<long long>(cols)) {
      continue;
    }
    const auto w =
        static_cast<Vertex>(static_cast<std::size_t>(r2) * cols + static_cast<std::size_t>(c2));
    if (!present.count(pair_key(static_cast<Vertex>(v), w))) add(static_cast<Vertex>(v), w);
  }
  // Connectivity: restore removed streets that join different components (Kruskal order).
  UnionFind uf(n);
  for (const auto& [a, b] : pairs) uf.unite(a, b);
  for (const auto& [a, b] : removed) {
    if (uf.unite(a, b)) add(a, b);
  }
  EdgeList<W> g{n, {}};
  g.edges.reserve(pairs.size());
  for (const auto& [a, b] : pairs) {
    const double dx = x[a] - x[b];
    const double dy = y[a] - y[b];
    const double len = std::sqrt(dx * dx + dy * dy);
    if constexpr (std::is_integral_v<W>) {
      g.edges.push_back(
          {a, b, std::max<W>(1, static_cast<W>(std::floor(len * road.integer_scale + 0.5)))});
    } else {
      g.edges.push_back({a, b, len});
    }
  }
  return g;
}

template <class W>
EdgeList<W> random_tree(std::size_t n, const GenOptions& opt) {
  check_n(n);
  validate_weight_spec(opt.weights);
  std::vector<Pair> pairs;
  if (n == 2) pairs.emplace_back(0, 1);
  if (n > 2) {
    Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
    std::vector<Vertex> code(n - 2);
    std::vector<std::size_t> degree(n, 1);
    for (auto& c : code) {
      c = static_cast<Vertex>(rng.uniform_below(n));
      ++degree[c];
    }
    // Linear-time Prüfer decoding.
    std::size_t ptr = 0;
    while (degree[ptr] != 1) ++ptr;
    std::size_t leaf = ptr;
    for (Vertex v : code) {
      pairs.emplace_back(static_cast<Vertex>(leaf), v);
      if (--degree[v] == 1 && v < ptr) {
        leaf = v;
      } else {
        ++ptr;
        while (degree[ptr] != 1) ++ptr;
        leaf = ptr;
      }
    }
    pairs.emplace_back(static_cast<Vertex>(leaf), static_cast<Vertex>(n - 1));
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> chung_lu(std::size_t n, double avg_degree, double exponent, const GenOptions& opt) {
  check_n(n);
  validate_weight_spec(opt.weights);
  check(std::isfinite(exponent) && exponent > 2.0, "chung_lu: exponent must be > 2");
  const std::size_t m = rounded_edge_count(n, avg_degree);
  const std::uint64_t max_pairs = static_cast<std::uint64_t>(n) * (n > 0 ? n - 1 : 0) / 2;
  check(m <= max_pairs / 2, "chung_lu: avg_degree too large for a simple graph");
  const double alpha = 1.0 / (exponent - 1.0);
  std::vector<double> prefix(n);
  double total = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    total += det_exp(-alpha * det_log(static_cast<double>(i + 1)));
    prefix[i] = total;
  }
  Xoshiro256ss rng(substream_seed(opt.seed, kTopology));
  auto draw = [&] {
    const double t = rng.uniform01() * total;
    const auto it = std::upper_bound(prefix.begin(), prefix.end(), t);
    return static_cast<Vertex>(std::min<std::size_t>(it - prefix.begin(), n - 1));
  };
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(2 * m);
  std::vector<Pair> pairs;
  pairs.reserve(m);
  const std::size_t max_attempts = 1000 + 100 * m;
  for (std::size_t attempt = 0; pairs.size() < m; ++attempt) {
    if (attempt >= max_attempts) throw std::runtime_error("chung_lu: too many rejected pairs");
    const Vertex a = draw();
    const Vertex b = draw();
    if (a != b && seen.insert(pair_key(a, b)).second) pairs.emplace_back(a, b);
  }
  return finish<W>(n, pairs, opt);
}

template <class W>
EdgeList<W> adversarial_decrease_key(std::size_t k) {
  check_n(k + 1);
  EdgeList<W> g{k + 1, {}};
  if (k >= 1) g.edges.push_back({0, 1, W{1}});
  for (std::size_t i = 2; i <= k; ++i) {
    g.edges.push_back({0, static_cast<Vertex>(i), static_cast<W>(2 * k)});
  }
  for (std::size_t i = 1; i <= k; ++i) {
    for (std::size_t j = i + 1; j <= k; ++j) {
      g.edges.push_back(
          {static_cast<Vertex>(i), static_cast<Vertex>(j), static_cast<W>(2 * (j - i) - 1)});
    }
  }
  return g;
}

template <class W>
EdgeList<W> adversarial_long_path(std::size_t n) {
  check_n(n);
  EdgeList<W> g{n, {}};
  for (std::size_t i = 0; i + 1 < n; ++i) {
    g.edges.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 1), W{1}});
    if (i + 2 < n) g.edges.push_back({static_cast<Vertex>(i), static_cast<Vertex>(i + 2), W{3}});
  }
  return g;
}

// ------------------------------------------------------------------ registry

namespace {
constexpr std::array<GeneratorInfo, 9> kCatalog{{
    {"er", "Erdos-Renyi G(n,m), m = round(n*avg_degree/2), simple", "n, avg_degree"},
    {"regular", "random d-regular graph (configuration model)", "n, degree, multi"},
    {"ba", "Barabasi-Albert preferential attachment, connected", "n, attach"},
    {"grid", "2D 4-neighbour grid, connected", "rows, cols (default floor(sqrt(n)))"},
    {"road", "road-like perturbed grid, Euclidean weights, connected (ignores --weights)",
     "rows, cols, road options"},
    {"tree", "uniformly random labelled tree (Pruefer)", "n"},
    {"chung_lu", "Chung-Lu scale-free graph with power-law exponent", "n, avg_degree, exponent"},
    {"adversarial_decrease_key", "Dijkstra decrease-key adversary, k = n-1 (fixed weights)", "n"},
    {"adversarial_long_path", "path with decoy skip edges, n-1 hop shortest path (fixed weights)",
     "n"},
}};

std::size_t side(std::size_t n) {
  std::size_t s = static_cast<std::size_t>(std::sqrt(static_cast<double>(n)));
  while (s * s > n) --s;
  while ((s + 1) * (s + 1) <= n) ++s;
  return s;
}
}  // namespace

std::span<const GeneratorInfo> generator_catalog() noexcept {
  return kCatalog;
}

MultiEdgePolicy parse_multi_edge_policy(std::string_view name) {
  if (name == "keep") return MultiEdgePolicy::kKeep;
  if (name == "drop") return MultiEdgePolicy::kDrop;
  if (name == "switch") return MultiEdgePolicy::kSwitchToSimple;
  throw std::invalid_argument("unknown multi-edge policy: " + std::string(name) +
                              " (keep|drop|switch)");
}

template <class W>
EdgeList<W> generate(const GraphSpec& s) {
  const auto& o = s.options;
  const std::size_t rows = s.rows ? s.rows : side(s.n);
  const std::size_t cols = s.cols ? s.cols : side(s.n);
  if (s.type == "er") return erdos_renyi<W>(s.n, s.avg_degree, o);
  if (s.type == "regular") return random_regular<W>(s.n, s.degree, s.multi, o);
  if (s.type == "ba") return barabasi_albert<W>(s.n, s.attach, o);
  if (s.type == "grid") return grid_2d<W>(rows, cols, o);
  if (s.type == "road") return road_like<W>(rows, cols, s.road, o);
  if (s.type == "tree") return random_tree<W>(s.n, o);
  if (s.type == "chung_lu") return chung_lu<W>(s.n, s.avg_degree, s.exponent, o);
  if (s.type == "adversarial_decrease_key") {
    check(s.n >= 1, "adversarial_decrease_key: n must be >= 1");
    return adversarial_decrease_key<W>(s.n - 1);
  }
  if (s.type == "adversarial_long_path") return adversarial_long_path<W>(s.n);
  throw std::invalid_argument("unknown generator type: " + s.type + " (see --list)");
}

#define SSSP_INSTANTIATE_GENERATORS(W)                                                  \
  template std::size_t count_components<W>(const EdgeList<W>&);                         \
  template void connect_components<W>(EdgeList<W>&, const WeightSpec&, std::uint64_t);  \
  template EdgeList<W> erdos_renyi<W>(std::size_t, double, const GenOptions&);          \
  template EdgeList<W> erdos_renyi_gnm<W>(std::size_t, std::size_t, const GenOptions&); \
  template EdgeList<W> random_regular<W>(std::size_t, std::size_t, MultiEdgePolicy,     \
                                         const GenOptions&);                            \
  template EdgeList<W> barabasi_albert<W>(std::size_t, std::size_t, const GenOptions&); \
  template EdgeList<W> grid_2d<W>(std::size_t, std::size_t, const GenOptions&);         \
  template EdgeList<W> road_like<W>(std::size_t, std::size_t, const RoadLikeOptions&,   \
                                    const GenOptions&);                                 \
  template EdgeList<W> random_tree<W>(std::size_t, const GenOptions&);                  \
  template EdgeList<W> chung_lu<W>(std::size_t, double, double, const GenOptions&);     \
  template EdgeList<W> adversarial_decrease_key<W>(std::size_t);                        \
  template EdgeList<W> adversarial_long_path<W>(std::size_t);                           \
  template EdgeList<W> generate<W>(const GraphSpec&);

SSSP_INSTANTIATE_GENERATORS(std::uint64_t)
SSSP_INSTANTIATE_GENERATORS(double)
#undef SSSP_INSTANTIATE_GENERATORS

}  // namespace sssp::gen
