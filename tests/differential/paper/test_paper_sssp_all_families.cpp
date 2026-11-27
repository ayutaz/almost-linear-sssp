// P4-DIFF-ALL: Mode A (Algorithm 8, arXiv 2611.34911v1 Sec 8.3, Theorem 1) on every generator
// family against three oracles.
//
// Under test: sssp::paper::paper_sssp with both variants, "paper-faithful" (scan Dijkstra base
// case, Alg 8 l.5 / Lemma 3.2) and "paper-faithful-d1" (binary-heap base case, deviation D-1,
// ADR 0008), for std::uint64_t and double.
//
// Oracles, none of which shares code with Mode A:
//   (1) sssp::baseline::bellman_ford,
//   (2) sssp::baseline::dijkstra_binary_heap,
//   (3) an exact rational reference written in this file: every weight is a dyadic rational
//       m * 2^e, all values are scaled to integers by a common power of two and summed in an
//       arbitrary-precision natural number (BigNat), and an O(n^2) scan Dijkstra runs on them.
//       It never rounds and never overflows, so it gives the true distances d*(v).
//
// Expectations (derived from the spec, not from an implementation run):
//   * uint64: dist == d* exactly, == Bellman-Ford, == Dijkstra (docs/design/numerics.md Sec 1).
//     Overflow (paper.hpp): if some reachable d*(v) >= 2^64 - 1 (the infinity sentinel) Mode A
//     must throw std::overflow_error. It MAY also throw when some reachable u with d*(u) < 2^64-1
//     has an arc whose sum d*(u) + w reaches the sentinel (a tentative value would reach it, as
//     in every Dijkstra baseline); in every other case it must return.
//   * double: Mode A is zero-edge contraction (adding 0 is exact), a lightest edge per pair
//     (monotone: a heavier parallel edge never gives a smaller rounded sum), Dijkstra (Lemma 3.2)
//     and the tight-arc restoration of Lemma 8.5, i.e. a relaxation algorithm in the rounded
//     algebra (double, min, fl(a + w)) of numerics.md Sec 2. Hence dist must equal Bellman-Ford
//     BIT FOR BIT and verify_sssp must pass with Tolerance::exact(). Against the exact reference,
//     |dist(v) - d*(v)| <= gamma_k * d*(v), gamma_k = k u / (1 - k u), u = 2^-53, with
//     k = max(hops of the reported tree path, hops of an exact shortest path) (numerics.md
//     Sec 3.1: lower bound from the tree path, upper bound from the exact path). The check is
//     done in exact integer arithmetic: |D - D*| * (2^53 - k) <= k * D*.
//     Overflow: Bellman-Ford throwing (some reachable rounded distance is +inf) implies Mode A
//     must throw; Mode A may also throw if some arc gives fl(d(u) + w) = +inf.
//   * Fallback policy (CLAUDE.md, ADR 0008): every returned run carries exactly one
//     FallbackRecord, reason as written in ADR 0008, vertices = N = number of zero-edge classes of
//     the source component, edges = m' = number of edges of the source component whose ends lie
//     in different classes (Alg 8 l.2-3; parallel edges kept, loops dropped). N and m' are
//     computed here by a BFS + union-find. No SSSP_DIAGNOSTICS guard: the record must also be
//     present in the release presets.
//
// Known, excluded on purpose: none of the inputs below reproduces the two Evaluate arithmetic
// counterexamples of P3-EVAL-ARITH-FIX (uint64 weight 2^63, double 2^53 + 1). Mode A for
// N < 2^64 never calls Evaluate (Alg 8 l.4-5 base case), so inputs of the same shape are kept
// here; if one of them ever fails, it is attributed to P3-EVAL-ARITH-FIX in the report rather
// than counted as a new counterexample.
//
// Set SSSP_P4_REPORT=1 to print the per-test statistics (runs, vertices, fallback distribution,
// maximal relative error of double results).

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "sssp/baseline.hpp"
#include "sssp/generators.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "support/enumerate.hpp"
#include "support/graphs.hpp"
#include "support/prng.hpp"

namespace {

using sssp::CsrGraph;
using sssp::Edge;
using sssp::EdgeList;
using sssp::kNoEdge;
using sssp::kNoVertex;
using sssp::SsspResult;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::PaperOptions;
using sssp::paper::PaperVariant;
using U64 = std::uint64_t;

// ------------------------------------------------------------------ arbitrary-precision naturals

/// Natural number in base 2^32, little-endian, no leading zero limbs (zero = no limbs).
class BigNat {
 public:
  BigNat() = default;
  explicit BigNat(U64 x) {
    while (x != 0) {
      limbs_.push_back(static_cast<std::uint32_t>(x));
      x >>= 32U;
    }
  }

  /// m * 2^shift.
  static BigNat shifted(U64 m, std::size_t shift) {
    BigNat r(m);
    if (r.limbs_.empty()) return r;
    r.limbs_.insert(r.limbs_.begin(), shift / 32, 0U);
    const unsigned bits = static_cast<unsigned>(shift % 32);
    if (bits != 0) {
      std::uint32_t carry = 0;
      for (auto& l : r.limbs_) {
        const std::uint32_t next = l >> (32U - bits);
        l = (l << bits) | carry;
        carry = next;
      }
      if (carry != 0) r.limbs_.push_back(carry);
    }
    return r;
  }

  [[nodiscard]] bool is_zero() const noexcept { return limbs_.empty(); }

  friend int compare(const BigNat& a, const BigNat& b) noexcept {
    if (a.limbs_.size() != b.limbs_.size()) return a.limbs_.size() < b.limbs_.size() ? -1 : 1;
    for (std::size_t i = a.limbs_.size(); i-- > 0;) {
      if (a.limbs_[i] != b.limbs_[i]) return a.limbs_[i] < b.limbs_[i] ? -1 : 1;
    }
    return 0;
  }
  friend bool operator==(const BigNat& a, const BigNat& b) noexcept { return compare(a, b) == 0; }
  friend bool operator<(const BigNat& a, const BigNat& b) noexcept { return compare(a, b) < 0; }
  friend bool operator<=(const BigNat& a, const BigNat& b) noexcept { return compare(a, b) <= 0; }

  friend BigNat operator+(const BigNat& a, const BigNat& b) {
    BigNat r;
    const std::size_t len = std::max(a.limbs_.size(), b.limbs_.size());
    r.limbs_.resize(len);
    U64 carry = 0;
    for (std::size_t i = 0; i < len; ++i) {
      const U64 s = carry + (i < a.limbs_.size() ? a.limbs_[i] : 0U) +
                    (i < b.limbs_.size() ? b.limbs_[i] : 0U);
      r.limbs_[i] = static_cast<std::uint32_t>(s);
      carry = s >> 32U;
    }
    if (carry != 0) r.limbs_.push_back(static_cast<std::uint32_t>(carry));
    return r;
  }

  /// |a - b|.
  friend BigNat abs_diff(const BigNat& a, const BigNat& b) {
    const bool swap = a < b;
    const BigNat& hi = swap ? b : a;
    const BigNat& lo = swap ? a : b;
    BigNat r = hi;
    std::int64_t borrow = 0;
    for (std::size_t i = 0; i < r.limbs_.size(); ++i) {
      std::int64_t d = static_cast<std::int64_t>(r.limbs_[i]) - borrow -
                       static_cast<std::int64_t>(i < lo.limbs_.size() ? lo.limbs_[i] : 0U);
      borrow = d < 0 ? 1 : 0;
      if (d < 0) d += std::int64_t{1} << 32;
      r.limbs_[i] = static_cast<std::uint32_t>(d);
    }
    r.trim();
    return r;
  }

  /// a * k.
  [[nodiscard]] BigNat times(U64 k) const {
    const BigNat lo = times32(static_cast<std::uint32_t>(k));
    BigNat hi = times32(static_cast<std::uint32_t>(k >> 32U));
    if (!hi.limbs_.empty()) hi.limbs_.insert(hi.limbs_.begin(), 0U);
    return lo + hi;
  }

  /// value ~= first * 2^second (for reporting only).
  [[nodiscard]] std::pair<double, int> approx() const {
    double m = 0.0;
    const std::size_t k = limbs_.size();
    const std::size_t take = std::min<std::size_t>(k, 3);
    for (std::size_t i = 0; i < take; ++i) m = m * 4294967296.0 + limbs_[k - 1 - i];
    return {m, static_cast<int>(32 * (k - take))};
  }

  /// Exact conversion to uint64 if it fits.
  [[nodiscard]] bool fits_u64(U64& out) const {
    if (limbs_.size() > 2) return false;
    out = 0;
    for (std::size_t i = limbs_.size(); i-- > 0;) out = (out << 32U) | limbs_[i];
    return true;
  }

 private:
  [[nodiscard]] BigNat times32(std::uint32_t k) const {
    BigNat r;
    if (k == 0 || limbs_.empty()) return r;
    r.limbs_.resize(limbs_.size());
    U64 carry = 0;
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      const U64 p = static_cast<U64>(limbs_[i]) * k + carry;
      r.limbs_[i] = static_cast<std::uint32_t>(p);
      carry = p >> 32U;
    }
    if (carry != 0) r.limbs_.push_back(static_cast<std::uint32_t>(carry));
    return r;
  }
  void trim() {
    while (!limbs_.empty() && limbs_.back() == 0) limbs_.pop_back();
  }
  std::vector<std::uint32_t> limbs_;
};

double ratio(const BigNat& a, const BigNat& b) {
  if (a.is_zero()) return 0.0;
  const auto [ma, ea] = a.approx();
  const auto [mb, eb] = b.approx();
  return std::ldexp(ma / mb, ea - eb);
}

// ------------------------------------------------------------------ dyadic decomposition

/// x = m * 2^e with m odd (or m = 0). Exact: frexp and ldexp by a power of two never round.
struct Dyadic {
  U64 m = 0;
  int e = 0;
};

Dyadic dyadic(U64 w) {
  Dyadic d{w, 0};
  while (d.m != 0 && (d.m & 1U) == 0) {
    d.m >>= 1U;
    ++d.e;
  }
  return d;
}

Dyadic dyadic(double w) {
  if (w == 0.0) return {};
  int k = 0;
  const double f = std::frexp(w, &k);  // w = f * 2^k, f in [0.5, 1)
  Dyadic d{static_cast<U64>(std::ldexp(f, 53)), k - 53};
  while ((d.m & 1U) == 0) {
    d.m >>= 1U;
    ++d.e;
  }
  return d;
}

/// Value of d times 2^-scale_exp as a natural number (requires d.e >= scale_exp).
BigNat scaled(const Dyadic& d, int scale_exp) {
  return BigNat::shifted(d.m, static_cast<std::size_t>(d.e - scale_exp));
}

// ------------------------------------------------------------------ exact reference

struct ExactReference {
  std::vector<char> reachable;
  std::vector<BigNat> dist;         ///< d*(v) * 2^-scale_exp
  std::vector<std::uint64_t> hops;  ///< hop count of one exact shortest path (fewest hops)
  std::vector<BigNat> weight;       ///< scaled weight per edge id
};

/// Binary-heap Dijkstra (lazy deletion) over exact naturals, on the undirected edge list; ties in
/// distance are broken towards fewer hops, so hops[v] is the least hop count of an exact
/// shortest path.
template <class W>
ExactReference exact_reference(const EdgeList<W>& el, Vertex s, int scale_exp) {
  const std::size_t n = el.num_vertices;
  ExactReference x;
  x.reachable.assign(n, 0);
  x.dist.assign(n, BigNat{});
  x.hops.assign(n, 0);
  x.weight.reserve(el.edges.size());
  std::vector<std::vector<std::pair<Vertex, std::size_t>>> adj(n);
  for (std::size_t i = 0; i < el.edges.size(); ++i) {
    const auto& e = el.edges[i];
    x.weight.push_back(scaled(dyadic(e.weight), scale_exp));
    adj[e.u].emplace_back(e.v, i);
    if (e.u != e.v) adj[e.v].emplace_back(e.u, i);
  }
  std::vector<char> done(n, 0);
  x.reachable[s] = 1;
  struct Item {
    BigNat d;
    std::uint64_t h;
    Vertex v;
  };
  const auto later = [](const Item& a, const Item& b) {  // min-heap on (d, h)
    const int c = compare(a.d, b.d);
    return c != 0 ? c > 0 : a.h > b.h;
  };
  std::priority_queue<Item, std::vector<Item>, decltype(later)> heap(later);
  heap.push(Item{BigNat{}, 0, s});
  while (!heap.empty()) {
    const Item top = heap.top();
    heap.pop();
    const Vertex best = top.v;
    if (done[best] != 0 || !(top.d == x.dist[best]) || top.h != x.hops[best]) continue;
    done[best] = 1;
    for (const auto& [v, id] : adj[best]) {
      if (done[v] != 0) continue;
      const BigNat c = x.dist[best] + x.weight[id];
      const std::uint64_t h = x.hops[best] + 1;
      if (x.reachable[v] == 0 || c < x.dist[v] || (c == x.dist[v] && h < x.hops[v])) {
        x.reachable[v] = 1;
        x.dist[v] = c;
        x.hops[v] = h;
        heap.push(Item{c, h, v});
      }
    }
  }
  return x;
}

// ------------------------------------------------------------------ expected fallback record

std::string expected_reason(PaperVariant v) {  // ADR 0008, Decision
  return v == PaperVariant::kFaithful
             ? "outer: N <= L_N^64 (paper base case, scan Dijkstra)"
             : "outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)";
}

struct Contracted {
  std::uint64_t vertices = 0;  ///< N
  std::uint64_t edges = 0;     ///< m'
};

/// Alg 8 l.2-3 sizes: source component by BFS over all edges, zero classes by union-find.
template <class W>
Contracted contracted_size(const EdgeList<W>& el, Vertex s) {
  const std::size_t n = el.num_vertices;
  std::vector<std::vector<Vertex>> adj(n);
  for (const auto& e : el.edges) {
    adj[e.u].push_back(e.v);
    adj[e.v].push_back(e.u);
  }
  std::vector<char> in(n, 0);
  std::vector<Vertex> q{s};
  in[s] = 1;
  for (std::size_t i = 0; i < q.size(); ++i) {
    for (const Vertex v : adj[q[i]]) {
      if (in[v] == 0) {
        in[v] = 1;
        q.push_back(v);
      }
    }
  }
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), std::size_t{0});
  const auto find = [&p](std::size_t a) {
    while (p[a] != a) a = p[a] = p[p[a]];
    return a;
  };
  for (const auto& e : el.edges) {
    if (e.weight == W{0}) p[find(e.u)] = find(e.v);
  }
  Contracted c;
  for (std::size_t v = 0; v < n; ++v) {
    if (in[v] != 0 && find(v) == v) ++c.vertices;
  }
  for (const auto& e : el.edges) {
    if (in[e.u] != 0 && find(e.u) != find(e.v)) ++c.edges;
  }
  return c;
}

// ------------------------------------------------------------------ statistics

struct Stats {
  std::uint64_t instances = 0;   ///< (graph, source) pairs
  std::uint64_t runs = 0;        ///< Mode A calls (instances x variants)
  std::uint64_t vertices = 0;    ///< sum of n over runs
  std::uint64_t edges = 0;       ///< sum of m over runs
  std::uint64_t mismatches = 0;  ///< runs with any failed check
  std::uint64_t required_throws = 0;
  std::uint64_t allowed_throws = 0;     ///< thrown where only a tentative sum overflows
  std::uint64_t fb_vertices = 0;        ///< sum of FallbackRecord::vertices
  std::uint64_t fb_edges = 0;           ///< sum of FallbackRecord::edges
  std::uint64_t contracted_runs = 0;    ///< N < n (zero classes merged or vertices unreached)
  std::uint64_t single_class_runs = 0;  ///< N == 1
  std::uint64_t full_runs = 0;          ///< N == n (nothing merged, everything reached)
  std::map<std::string, std::uint64_t> reasons;
  double max_rel_err = 0.0;
  std::uint64_t inexact_double_vertices = 0;  ///< dist != d* exactly (rounding)

  void print(const std::string& label) const {
    if (std::getenv("SSSP_P4_REPORT") == nullptr) return;
    std::ostringstream os;
    os << "[P4-DIFF-ALL] " << label << ": instances=" << instances << " runs=" << runs
       << " vertices=" << vertices << " edges=" << edges << " mismatches=" << mismatches
       << " required_throws=" << required_throws << " allowed_throws=" << allowed_throws
       << " fallback_vertices=" << fb_vertices << " fallback_edges=" << fb_edges
       << " runs_with_N<n=" << contracted_runs << " runs_with_N=1=" << single_class_runs
       << " runs_with_N=n=" << full_runs << " max_rel_err=" << max_rel_err
       << " inexact_double_vertices=" << inexact_double_vertices;
    for (const auto& [r, c] : reasons) os << " | reason \"" << r << "\" x" << c;
    std::cout << os.str() << '\n';
  }
};

template <class W>
std::string describe(const std::string& label, const EdgeList<W>& el, Vertex s) {
  std::ostringstream os;
  os.precision(17);
  os << label << " n=" << el.num_vertices << " m=" << el.edges.size() << " source=" << s;
  if (el.edges.size() <= 40) {
    os << " edges:";
    for (const auto& e : el.edges) os << " {" << e.u << "," << e.v << "," << e.weight << "}";
  }
  return os.str();
}

template <class F>
bool throws_overflow(F&& f) {
  try {
    static_cast<void>(f());
  } catch (const std::overflow_error&) {
    return true;
  }
  return false;
}

// ------------------------------------------------------------------ the differential check

/// Runs both Mode A variants on (el, s) and checks them against the three oracles, the tree
/// verifier and the fallback expectation. Every failure is reported with FAIL_CHECK and counted.
template <class W>
void check_instance(const std::string& label, const EdgeList<W>& el, Vertex s, Stats& st) {
  constexpr bool kExact = WeightTraits<W>::is_exact;
  const auto g = sssp::to_csr(el);
  const std::size_t n = el.num_vertices;
  const W inf = WeightTraits<W>::infinity();
  ++st.instances;

  // Oracles (1) and (2); their overflow behaviour is part of the expectation for double.
  SsspResult<W> bf;
  const bool bf_threw = throws_overflow([&] { return bf = sssp::baseline::bellman_ford(g, s); });
  SsspResult<W> dij;
  const bool dij_threw =
      throws_overflow([&] { return dij = sssp::baseline::dijkstra_binary_heap(g, s); });

  // Mode A runs.
  struct Run {
    PaperVariant v;
    bool threw = false;
    SsspResult<W> r;
  };
  std::vector<Run> runs;
  for (const PaperVariant v : sssp::paper::kPaperVariants) {
    Run run{v, false, {}};
    run.threw =
        throws_overflow([&] { return run.r = sssp::paper::paper_sssp(g, s, PaperOptions{v}); });
    runs.push_back(std::move(run));
  }

  // Oracle (3): common scale 2^scale_exp over every weight and every reported double distance.
  int scale_exp = 0;
  if constexpr (!kExact) {
    scale_exp = std::numeric_limits<int>::max();
    const auto see = [&scale_exp](double x) {
      if (x != 0.0 && std::isfinite(x)) scale_exp = std::min(scale_exp, dyadic(x).e);
    };
    for (const auto& e : el.edges) see(e.weight);
    for (const auto& run : runs) {
      if (!run.threw) {
        for (const W d : run.r.dist) see(d);
      }
    }
    if (scale_exp == std::numeric_limits<int>::max()) scale_exp = 0;
  }
  const ExactReference ex = exact_reference(el, s, scale_exp);
  const Contracted cs = contracted_size(el, s);

  // Overflow expectation.
  bool must_throw = false;
  bool may_throw = false;
  if constexpr (kExact) {
    const BigNat sentinel(inf);
    for (std::size_t v = 0; v < n; ++v) {
      if (ex.reachable[v] != 0 && sentinel <= ex.dist[v]) must_throw = true;
    }
    may_throw = must_throw;
    for (std::size_t i = 0; i < el.edges.size(); ++i) {
      const auto& e = el.edges[i];
      for (const Vertex u : {e.u, e.v}) {
        if (ex.reachable[u] != 0 && ex.dist[u] < sentinel &&
            sentinel <= ex.dist[u] + ex.weight[i]) {
          may_throw = true;
        }
      }
    }
    if (bf_threw != must_throw) {  // oracle consistency: (1) and (3) must agree on overflow
      ++st.mismatches;
      FAIL_CHECK("oracle disagreement on overflow: bellman_ford threw=" << bf_threw << " exact="
                                                                        << must_throw << "\n"
                                                                        << describe(label, el, s));
      return;
    }
  } else {
    must_throw = bf_threw;
    may_throw = must_throw;
    if (!bf_threw) {
      for (const auto& e : el.edges) {
        for (const auto [a, b] : {std::pair{e.u, e.v}, std::pair{e.v, e.u}}) {
          static_cast<void>(b);
          if (bf.dist[a] != inf && !std::isfinite(bf.dist[a] + e.weight)) may_throw = true;
        }
      }
    }
  }

  for (auto& run : runs) {
    ++st.runs;
    st.vertices += n;
    st.edges += el.edges.size();
    const std::string where =
        std::string(sssp::paper::algorithm_name(run.v)) + " " + describe(label, el, s);
    std::vector<std::string> fails;
    if (run.threw) {
      if (must_throw) {
        ++st.required_throws;
      } else if (may_throw) {
        ++st.allowed_throws;
      } else {
        fails.emplace_back("threw std::overflow_error although no sum reaches the sentinel");
      }
    } else if (must_throw) {
      fails.emplace_back("returned although a reachable distance reaches the sentinel");
    } else {
      const auto& r = run.r;
      // Tree verifier, alone and against oracles (1) and (2); exact also for double (Sec 2).
      const auto rep = sssp::verify_sssp(g, s, r);
      if (!rep.ok) fails.push_back("verify_sssp: " + rep.summary());
      if (!bf_threw) {
        const auto rb = sssp::verify_sssp(g, s, r, &bf.dist);
        if (!rb.ok) fails.push_back("verify_sssp vs bellman_ford: " + rb.summary());
        if (r.dist != bf.dist) fails.emplace_back("dist differs from bellman_ford (bitwise)");
      }
      if (!dij_threw) {
        const auto rd = sssp::verify_sssp(g, s, r, &dij.dist);
        if (!rd.ok) fails.push_back("verify_sssp vs dijkstra: " + rd.summary());
        if (r.dist != dij.dist) fails.emplace_back("dist differs from dijkstra (bitwise)");
      }
      // Oracle (3), and the hop count of the reported tree path.
      if (r.dist.size() == n && r.parent.size() == n && r.parent_edge.size() == n) {
        for (std::size_t v = 0; v < n && fails.size() < 8; ++v) {
          const bool finite = r.dist[v] != inf;
          if (finite != (ex.reachable[v] != 0)) {
            fails.push_back("vertex " + std::to_string(v) + ": reachability differs from exact");
            continue;
          }
          if (!finite) {
            if (r.parent[v] != kNoVertex || r.parent_edge[v] != kNoEdge) {
              fails.push_back("vertex " + std::to_string(v) + ": unreachable with a parent");
            }
            continue;
          }
          std::uint64_t tree_hops = 0;
          for (Vertex x = static_cast<Vertex>(v); x != s; x = r.parent[x]) {
            if (r.parent[x] >= n || ++tree_hops > n) {
              tree_hops = n;
              fails.push_back("vertex " + std::to_string(v) + ": parent chain misses source");
              break;
            }
          }
          if constexpr (kExact) {
            U64 d = 0;
            if (!ex.dist[v].fits_u64(d) || d != r.dist[v]) {
              fails.push_back("vertex " + std::to_string(v) + ": dist " +
                              std::to_string(r.dist[v]) + " != exact");
            }
          } else {
            const BigNat got = scaled(dyadic(r.dist[v]), scale_exp);
            const std::uint64_t k = std::max(tree_hops, ex.hops[v]);
            const BigNat err = abs_diff(got, ex.dist[v]);
            if (!err.is_zero()) ++st.inexact_double_vertices;
            const BigNat lhs = err.times((U64{1} << 53U) - k);
            const BigNat rhs = ex.dist[v].times(k);
            if (rhs < lhs) {
              std::ostringstream os;
              os.precision(17);
              os << "vertex " << v << ": dist " << r.dist[v] << " outside gamma_" << k
                 << " of the exact distance (rel err " << ratio(err, ex.dist[v]) << ")";
              fails.push_back(os.str());
            }
            if (!ex.dist[v].is_zero())
              st.max_rel_err = std::max(st.max_rel_err, ratio(err, ex.dist[v]));
          }
        }
      } else {
        fails.emplace_back("result vectors do not have n entries");
      }
      // Fallback record (CLAUDE.md fallback policy, ADR 0008).
      const auto& fb = r.diag.fallbacks;
      if (fb.size() != 1) {
        fails.push_back("expected exactly one FallbackRecord, got " + std::to_string(fb.size()));
      } else {
        ++st.reasons[fb[0].reason];
        st.fb_vertices += fb[0].vertices;
        st.fb_edges += fb[0].edges;
        if (fb[0].vertices < n) ++st.contracted_runs;
        if (fb[0].vertices == 1) ++st.single_class_runs;
        if (fb[0].vertices == n) ++st.full_runs;
        if (fb[0].reason != expected_reason(run.v)) {
          fails.push_back("fallback reason \"" + fb[0].reason + "\"");
        }
        if (fb[0].vertices != cs.vertices || fb[0].edges != cs.edges) {
          fails.push_back("fallback size (" + std::to_string(fb[0].vertices) + ", " +
                          std::to_string(fb[0].edges) + ") != contracted (" +
                          std::to_string(cs.vertices) + ", " + std::to_string(cs.edges) + ")");
        }
      }
    }
    if (!fails.empty()) {
      ++st.mismatches;
      std::string msg = where;
      for (const auto& f : fails) msg += "\n  " + f;
      FAIL_CHECK(msg);
    }
  }
}

// ------------------------------------------------------------------ generator families

struct Family {
  std::string label;  ///< catalog name plus the policy used
  std::string type;   ///< generator_catalog name
  sssp::gen::MultiEdgePolicy multi = sssp::gen::MultiEdgePolicy::kSwitchToSimple;
  sssp::gen::GeometricWeights geo = sssp::gen::GeometricWeights::kFromSpec;
};

/// Every catalog entry once; regular and geometric additionally with their second policy.
std::vector<Family> families() {
  std::vector<Family> f;
  for (const auto& info : sssp::gen::generator_catalog()) {
    f.push_back({std::string(info.name), std::string(info.name)});
  }
  f.push_back({"regular/keep", "regular", sssp::gen::MultiEdgePolicy::kKeep});
  f.push_back({"geometric/euclidean", "geometric", sssp::gen::MultiEdgePolicy::kSwitchToSimple,
               sssp::gen::GeometricWeights::kEuclidean});
  return f;
}

struct WeightCase {
  std::string label;
  sssp::gen::WeightSpec spec;
};

template <class W>
std::vector<WeightCase> weight_cases() {
  using sssp::gen::WeightDist;
  const auto mk = [](std::string label, WeightDist d, double lo, double hi) {
    sssp::gen::WeightSpec w;
    w.dist = d;
    w.lo = lo;
    w.hi = hi;
    return WeightCase{std::move(label), w};
  };
  std::vector<WeightCase> c;
  c.push_back(mk("int[1,4]", WeightDist::kUniformInt, 1, 4));  // many ties
  c.push_back(mk("int[1,2^40]", WeightDist::kUniformInt, 1, 0x1p40));
  {
    auto z = mk("zero_mixed p=0.3", WeightDist::kZeroMixed, 1, 9);
    z.spec.zero_probability = 0.3;
    c.push_back(z);
  }
  {
    auto z = mk("zero_mixed p=0.9", WeightDist::kZeroMixed, 1, 9);
    z.spec.zero_probability = 0.9;
    c.push_back(z);
  }
  {
    auto w = mk("wide_range 2^50", WeightDist::kWideRange, 1, 1);
    w.spec.wide_log2 = 50;
    c.push_back(w);
  }
  if constexpr (!WeightTraits<W>::is_exact) {
    c.push_back(mk("real[1,40)", WeightDist::kUniformReal, 1, 40));
    c.push_back(mk("log_uniform[1e-3,1e3]", WeightDist::kLogUniform, 1e-3, 1e3));
    {
      auto e = mk("exponential mean 10", WeightDist::kExponential, 0, 0);
      e.spec.mean = 10.0;
      c.push_back(e);
    }
  }
  return c;
}

template <class W>
EdgeList<W> generate(const Family& fam, const WeightCase& wc, std::size_t n, std::uint64_t seed) {
  sssp::gen::GraphSpec spec;
  spec.type = fam.type;
  spec.n = n;
  spec.avg_degree = fam.type == "geometric" ? 4.0 : 2.5;
  spec.degree = 3;
  spec.attach = 2;
  spec.exponent = 2.5;
  spec.multi = fam.multi;
  spec.geometric.weights = fam.geo;
  spec.options.seed = seed;
  spec.options.weights = wc.spec;
  return sssp::gen::generate<W>(spec);
}

/// Every family x weight case x seeds 1..num_seeds, from two sources per graph (vertex 0 and a
/// seed-dependent vertex). The
/// decrease-key adversary has Theta(n^2) edges, so its n is capped at dense_cap.
template <class W>
Stats sweep(std::size_t base_n, std::size_t step_n, std::uint64_t num_seeds,
            std::size_t dense_cap) {
  Stats st;
  std::uint64_t graphs = 0;
  for (const Family& fam : families()) {
    for (const WeightCase& wc : weight_cases<W>()) {
      for (std::uint64_t seed = 1; seed <= num_seeds; ++seed) {
        std::size_t n = base_n + step_n * static_cast<std::size_t>(seed);  // even
        if (fam.type == "adversarial_decrease_key") n = std::min(n, dense_cap);
        const EdgeList<W> el = generate<W>(fam, wc, n, seed);
        ++graphs;
        const std::size_t vn = el.num_vertices;
        const std::string label = fam.label + " " + wc.label + " seed=" + std::to_string(seed);
        for (const Vertex s : {Vertex{0}, static_cast<Vertex>((seed * 7919U) % vn)}) {
          check_instance(label, el, s, st);
        }
      }
    }
  }
  CHECK(graphs == families().size() * weight_cases<W>().size() * num_seeds);
  return st;
}

template <class W>
const char* type_name() {
  return WeightTraits<W>::is_exact ? "uint64" : "double";
}

}  // namespace

// ------------------------------------------------------------------ tests

TEMPLATE_TEST_CASE("P4-DIFF-ALL: Mode A equals three oracles on every family x 20 seeds",
                   "[differential][paper][sssp][p4-diff-all]", std::uint64_t, double) {
  using W = TestType;
  constexpr std::uint64_t kSeeds = 20;
  const Stats st = sweep<W>(28, 4, kSeeds, 60);  // n = 32 .. 108
  st.print(std::string("sweep ") + type_name<W>());
  CHECK(st.mismatches == 0);
  CHECK(st.runs == families().size() * weight_cases<W>().size() * kSeeds * 2 * 2);
  CHECK(st.reasons.size() == 2);   // both variants recorded their own reason
  CHECK(st.contracted_runs > 0);   // contraction or unreached vertices changed N
  CHECK(st.required_throws == 0);  // generator weights never overflow
}

TEMPLATE_TEST_CASE("P4-DIFF-ALL: Mode A on larger graphs of every family x 20 seeds",
                   "[differential][paper][sssp][p4-diff-all][slow]", std::uint64_t, double) {
  using W = TestType;
  constexpr std::uint64_t kSeeds = 20;
  const Stats st = sweep<W>(220, 20, kSeeds, 200);  // n = 240 .. 620
  st.print(std::string("large sweep ") + type_name<W>());
  CHECK(st.mismatches == 0);
  CHECK(st.runs == families().size() * weight_cases<W>().size() * kSeeds * 2 * 2);
}

TEMPLATE_TEST_CASE("P4-DIFF-ALL: boundary catalogues from every source",
                   "[differential][paper][sssp][p4-diff-all]", std::uint64_t, double) {
  using W = TestType;
  Stats st;
  for (const auto& c : sssp::test::required_graphs<W>()) {
    const EdgeList<W> el{c.n, c.edges};
    for (Vertex s = 0; s < c.n; s += (c.n > 40 ? 37 : 1)) check_instance(c.name, el, s, st);
  }
  for (const auto& c : sssp::gen::small_graph_suite<W>()) {
    for (Vertex s = 0; s < c.graph.num_vertices; ++s) check_instance(c.name, c.graph, s, st);
  }
  const W one{1};
  const W zero{0};
  const std::vector<std::pair<std::string, EdgeList<W>>> extra{
      {"n=1 no edges", {1, {}}},
      {"n=1 loops", {1, {{0, 0, zero}, {0, 0, W{5}}, {0, 0, one}}}},
      {"isolated vertices", {6, {{1, 2, one}, {4, 5, W{2}}}}},
      {"three components",
       {9,
        {{0, 1, one},
         {1, 2, W{2}},
         {3, 4, zero},
         {4, 5, one},
         {6, 7, W{3}},
         {7, 8, zero},
         {8, 6, one}}}},
      {"all zero, cycles",
       {5, {{0, 1, zero}, {1, 2, zero}, {2, 0, zero}, {3, 4, zero}, {2, 3, zero}, {4, 4, zero}}}},
      {"zero parallel to positive",
       {3, {{0, 1, W{7}}, {0, 1, zero}, {1, 2, W{4}}, {2, 1, zero}, {0, 2, W{9}}}}},
      {"multi edges both orientations",
       {4,
        {{0, 1, W{9}},
         {1, 0, W{3}},
         {0, 1, W{3}},
         {1, 2, W{5}},
         {2, 1, W{2}},
         {2, 3, one},
         {3, 2, one},
         {0, 3, W{8}}}}},
      {"self loops everywhere",
       {4,
        {{0, 0, one},
         {0, 1, W{2}},
         {1, 1, zero},
         {1, 2, W{2}},
         {2, 2, W{9}},
         {2, 3, zero},
         {3, 3, one}}}},
      {"zero class not containing source",
       {6,
        {{0, 1, W{4}},
         {1, 2, zero},
         {2, 3, zero},
         {3, 1, zero},
         {0, 4, W{2}},
         {4, 3, W{1}},
         {3, 5, W{6}}}}},
      {"zero edge to an unreachable part", {5, {{0, 1, one}, {2, 3, zero}, {3, 4, zero}}}},
  };
  for (const auto& [name, el] : extra) {
    for (Vertex s = 0; s < el.num_vertices; ++s) check_instance(name, el, s, st);
  }
  // L_N jumps (N + 2 a power of two) and their neighbours, on the contracted size N: a path of
  // N classes, each class a zero-weight star of 3 vertices.
  for (unsigned b = 2; b <= 9; ++b) {
    for (const std::uint64_t big_n : {(U64{1} << b) - 3, (U64{1} << b) - 2, (U64{1} << b) - 1}) {
      if (big_n == 0) continue;
      EdgeList<W> el{static_cast<std::size_t>(3 * big_n), {}};
      for (std::uint64_t i = 0; i < big_n; ++i) {
        const auto c = static_cast<Vertex>(3 * i);
        el.edges.push_back({c, c + 1, zero});
        el.edges.push_back({c + 2, c, zero});
        if (i + 1 < big_n) el.edges.push_back({c + 1, c + 3, W{static_cast<W>(i % 3 + 1)}});
      }
      check_instance("L_N jump N=" + std::to_string(big_n), el, 0, st);
      check_instance("L_N jump N=" + std::to_string(big_n), el,
                     static_cast<Vertex>(el.num_vertices - 1), st);
    }
  }
  st.print(std::string("boundary ") + type_name<W>());
  CHECK(st.mismatches == 0);
}

TEST_CASE("P4-DIFF-ALL: an empty graph or an invalid source throws std::out_of_range",
          "[differential][paper][sssp][p4-diff-all]") {
  const auto empty64 = CsrGraph<U64>::from_edges(0, {});
  const auto empty_d = CsrGraph<double>::from_edges(0, {});
  const auto g64 = CsrGraph<U64>::from_edges(3, {{0, 1, 1}});
  const auto gd = CsrGraph<double>::from_edges(3, {{0, 1, 1.0}});
  for (const PaperVariant v : sssp::paper::kPaperVariants) {
    CHECK_THROWS_AS(sssp::paper::paper_sssp(empty64, 0, PaperOptions{v}), std::out_of_range);
    CHECK_THROWS_AS(sssp::paper::paper_sssp(empty_d, 0, PaperOptions{v}), std::out_of_range);
    CHECK_THROWS_AS(sssp::paper::paper_sssp(g64, 3, PaperOptions{v}), std::out_of_range);
    CHECK_THROWS_AS(sssp::paper::paper_sssp(gd, 3, PaperOptions{v}), std::out_of_range);
    CHECK_THROWS_AS(sssp::paper::paper_sssp(g64, kNoVertex, PaperOptions{v}), std::out_of_range);
  }
}

TEMPLATE_TEST_CASE("P4-DIFF-ALL: every graph on 4 vertices with weights {0, 1, 2}",
                   "[differential][paper][sssp][p4-diff-all]", std::uint64_t, double) {
  using W = TestType;
  const std::vector<W> weights{W{0}, W{1}, W{2}};
  Stats st;
  std::uint64_t graphs = 0;
  sssp::test::for_each_small_graph<W>(4, std::span<const W>(weights), [&](const auto& edges) {
    ++graphs;
    const EdgeList<W> el{4, edges};
    for (Vertex s = 0; s < 4; ++s) check_instance("exhaustive n=4", el, s, st);
  });
  st.print(std::string("exhaustive n=4 ") + type_name<W>());
  CHECK(graphs == sssp::test::count_small_graphs(4, 3));
  CHECK(st.mismatches == 0);
}

TEMPLATE_TEST_CASE("P4-DIFF-ALL: every graph on 5 vertices with weights {0, 1}",
                   "[differential][paper][sssp][p4-diff-all][slow]", std::uint64_t, double) {
  using W = TestType;
  const std::vector<W> weights{W{0}, W{1}};
  Stats st;
  std::uint64_t graphs = 0;
  sssp::test::for_each_small_graph<W>(5, std::span<const W>(weights), [&](const auto& edges) {
    ++graphs;
    // Every labelled graph is enumerated, so source 0 covers every (graph, source) pair up to
    // isomorphism.
    check_instance("exhaustive n=5", EdgeList<W>{5, edges}, 0, st);
  });
  st.print(std::string("exhaustive n=5 ") + type_name<W>());
  CHECK(graphs == sssp::test::count_small_graphs(5, 2));
  CHECK(st.mismatches == 0);
}

TEST_CASE("P4-DIFF-ALL: uint64 weights near the infinity sentinel",
          "[differential][paper][sssp][p4-diff-all][overflow]") {
  constexpr U64 kMax = std::numeric_limits<U64>::max();  // the sentinel 2^64 - 1
  constexpr U64 k63 = U64{1} << 63U;
  Stats st;
  const std::vector<std::pair<std::string, EdgeList<U64>>> cases{
      {"single edge 2^64-2", {2, {{0, 1, kMax - 1}}}},
      {"single edge at the sentinel", {2, {{0, 1, kMax}}}},
      {"2^63 + 2^63", {3, {{0, 1, k63}, {1, 2, k63}}}},
      {"2^63 + (2^63 - 2)", {3, {{0, 1, k63}, {1, 2, k63 - 2}}}},
      {"2^63 + (2^63 - 1)", {3, {{0, 1, k63}, {1, 2, k63 - 1}}}},
      {"2^63 star", {4, {{0, 1, k63}, {0, 2, k63}, {0, 3, k63 - 1}}}},
      {"2^63 star from a leaf", {4, {{0, 1, k63}, {0, 2, k63}, {0, 3, 1}}}},
      {"heavy chord off the shortest path", {3, {{0, 1, 1}, {1, 2, 1}, {0, 2, kMax}}}},
      {"heavy parallel edge", {2, {{0, 1, kMax}, {0, 1, 3}, {1, 0, kMax - 1}}}},
      {"heavy edge back to the source", {3, {{0, 1, k63}, {1, 2, 1}, {2, 0, k63}}}},
      {"zero class behind a 2^63 edge", {4, {{0, 1, k63}, {1, 2, 0}, {2, 3, 0}, {3, 1, 0}}}},
      {"sentinel loop", {2, {{0, 0, kMax}, {0, 1, 5}, {1, 1, kMax}}}},
      {"2^62 chain of four",
       {5,
        {{0, 1, U64{1} << 62U},
         {1, 2, U64{1} << 62U},
         {2, 3, U64{1} << 62U},
         {3, 4, (U64{1} << 62U) - 1}}}},
      {"2^62 chain of four, last +0",
       {5,
        {{0, 1, U64{1} << 62U},
         {1, 2, U64{1} << 62U},
         {2, 3, U64{1} << 62U},
         {3, 4, U64{1} << 62U}}}},
  };
  for (const auto& [name, el] : cases) {
    for (Vertex s = 0; s < el.num_vertices; ++s) check_instance(name, el, s, st);
  }
  // Random: weights in [2^60, 2^62] on small trees and sparse graphs, so that some reachable
  // distances exceed the sentinel and some do not.
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    const std::size_t n = 4 + seed % 5;
    EdgeList<U64> tree{n,
                       sssp::test::random_tree_edges<U64>(n, U64{1} << 60U, U64{1} << 62U, seed)};
    EdgeList<U64> rnd{n, sssp::test::random_edges<U64>(n, 2 * n, U64{1} << 60U, kMax, seed)};
    for (Vertex s = 0; s < n; ++s) {
      check_instance("random huge tree seed=" + std::to_string(seed), tree, s, st);
      check_instance("random huge graph seed=" + std::to_string(seed), rnd, s, st);
    }
  }
  st.print("uint64 near sentinel");
  CHECK(st.mismatches == 0);
  CHECK(st.required_throws > 0);  // the cases above do exercise the overflow contract
}

TEST_CASE("P4-DIFF-ALL: double absorption, subnormal and huge weights",
          "[differential][paper][sssp][p4-diff-all][overflow]") {
  constexpr double kMax = std::numeric_limits<double>::max();
  constexpr double kDen = std::numeric_limits<double>::denorm_min();
  constexpr double k53 = 0x1p53;
  Stats st;
  const std::vector<std::pair<std::string, EdgeList<double>>> cases{
      {"2^53 + 1 absorbed", {3, {{0, 1, k53}, {1, 2, 1.0}}}},
      {"2^53 + 1 + 1 vs 2^53 + 2", {4, {{0, 1, k53}, {1, 2, 1.0}, {2, 3, 1.0}, {0, 3, k53 + 2}}}},
      {"2^53 + 2 + 1 rounds up", {3, {{0, 1, k53 + 2}, {1, 2, 1.0}}}},
      {"1 + 2^53 (order)", {3, {{0, 1, 1.0}, {1, 2, k53}, {0, 2, k53 + 2}}}},
      {"2^60 chain of ones",
       {6,
        {{0, 1, 0x1p60},
         {1, 2, 1.0},
         {2, 3, 1.0},
         {3, 4, 1.0},
         {4, 5, 1.0},
         {0, 5, 0x1p60 + 256}}}},
      {"subnormal path", {4, {{0, 1, kDen}, {1, 2, kDen}, {2, 3, 3 * kDen}, {0, 3, 4 * kDen}}}},
      {"subnormal vs one", {4, {{0, 1, 1.0}, {1, 2, kDen}, {2, 3, kDen}, {0, 3, 1.0}}}},
      {"tiny and huge", {4, {{0, 1, 1e-300}, {1, 2, 1e300}, {2, 3, 1e-300}, {0, 3, 1e300}}}},
      {"single DBL_MAX", {2, {{0, 1, kMax}}}},
      {"DBL_MAX + DBL_MAX", {3, {{0, 1, kMax}, {1, 2, kMax}}}},
      {"DBL_MAX/2 twice is finite", {3, {{0, 1, kMax / 2}, {1, 2, kMax / 2}}}},
      {"heavy chord off the shortest path", {3, {{0, 1, kMax}, {1, 2, 1.0}, {0, 2, 1.0}}}},
      {"heavy edge back to the source", {3, {{0, 1, kMax}, {1, 2, kMax}, {0, 2, kMax}}}},
      {"zero class behind DBL_MAX", {4, {{0, 1, kMax}, {1, 2, 0.0}, {2, 3, 0.0}}}},
  };
  for (const auto& [name, el] : cases) {
    for (Vertex s = 0; s < el.num_vertices; ++s) check_instance(name, el, s, st);
  }
  // Random: topology from the er / tree generators, weights redrawn here as (1 + U) * 2^e with
  // e uniform in [-1000, 1000] (the generators cap weights at 2^53; this needs ~2000-bit exact
  // integers in the reference), and wide-range weights over [1, 2^62] from the generator.
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    sssp::gen::GraphSpec spec;
    spec.type = seed % 2 == 0 ? "er" : "tree";
    spec.n = 30;
    spec.avg_degree = 3.0;
    spec.options.seed = seed;
    auto a = sssp::gen::generate<double>(spec);
    sssp::test::SplitMix64 rng(seed);
    for (auto& e : a.edges) {
      e.weight = std::ldexp(1.0 + rng.unit(), static_cast<int>(rng.between(0, 2000)) - 1000);
    }
    spec.options.weights.dist = sssp::gen::WeightDist::kWideRange;
    spec.options.weights.lo = 1.0;
    spec.options.weights.hi = 1.0;
    spec.options.weights.wide_log2 = 62;
    const auto b = sssp::gen::generate<double>(spec);
    for (const Vertex s : {Vertex{0}, Vertex{17}}) {
      check_instance("weights 2^[-1000,1000] seed=" + std::to_string(seed), a, s, st);
      check_instance("wide_range 2^62 seed=" + std::to_string(seed), b, s, st);
    }
  }
  st.print("double extremes");
  CHECK(st.mismatches == 0);
  CHECK(st.required_throws > 0);
  CHECK(st.inexact_double_vertices > 0);  // rounding did occur, so the tolerance was exercised
}

TEST_CASE("P4-DIFF-ALL: paper-faithful scan base case at N in the thousands",
          "[differential][paper][sssp][p4-diff-all][slow]") {
  Stats st;
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    sssp::gen::GraphSpec spec;
    spec.type = seed % 2 == 0 ? "road" : "geometric";
    spec.n = 2500;
    spec.avg_degree = 5.0;
    spec.options.seed = seed;
    spec.options.weights.dist = sssp::gen::WeightDist::kZeroMixed;
    spec.options.weights.zero_probability = 0.2;
    spec.options.weights.lo = 1;
    spec.options.weights.hi = 1000;
    check_instance("large " + spec.type + " seed=" + std::to_string(seed),
                   sssp::gen::generate<U64>(spec), static_cast<Vertex>(seed * 13), st);
    check_instance("large " + spec.type + " seed=" + std::to_string(seed),
                   sssp::gen::generate<double>(spec), static_cast<Vertex>(seed * 13), st);
  }
  st.print("large N");
  CHECK(st.mismatches == 0);
}
