// Lemma 5.9 (Sec 5.3, txt l.1152) of arXiv 2611.34911v1, exhaustively on small graphs
// [P4-PROPERTY-LEMMAS].
//
// Statement: in the local radial state representation, fix W > 0 and keep the block MST joins,
// the state-to-owner joins at theta and the pivot corridor joins of weight <= W; the undirected
// components restricted to the states are exactly the strongly connected components of the state
// arcs of length <= W. The strict version (< W) holds likewise (l.1172).
//
// Already covered elsewhere (not repeated here): test_local_joins_property.cpp checks the same
// equality, and I0 / O0 / I* / O* against brute force, on seeded random graph families.
// New here:
//   * every connected simple graph on N <= 4 vertices with weights from a multi-scale set (and,
//     tagged [slow], N = 5 with two weight sets), for every separator size b in 1 .. N + 1;
//   * every threshold W at which a partition can change, plus one W above all of them;
//   * the SCCs come from a bit-matrix transitive closure (Warshall), not Kosaraju, and the
//     join components from a separate union-find, both written here;
//   * besides the library's join records, the lemma is checked on reference joins rebuilt from
//     the definition (every within-window arc at its length, theta by a fixpoint bottleneck
//     iteration, the corridors at their weights). Components of the within-window arcs at a
//     threshold equal those of any block MST, so the reference needs no MST;
//   * seeded random graphs on 6..9 vertices through the shrinking driver (shrink_support.hpp):
//     a failure is reported with its seed and a 1-minimal edge list.
//
// The input precondition of the owner construction (connected, positive weights) is enforced by
// skipping disconnected graphs, which the enumeration also produces.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "property/paper/shrink_support.hpp"
#include "sssp/graph.hpp"
#include "support/enumerate.hpp"
#include "support/prng.hpp"

using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
using sssp::paper::RadialArcKind;
using sssp::paper::StateId;
namespace sh = sssp::test::shrink;

namespace {

template <class W>
constexpr W kInf = sssp::WeightTraits<W>::infinity();

struct Coverage {
  std::uint64_t instances = 0;
  std::uint64_t thresholds = 0;
  std::uint64_t multi_pivot = 0;      // instances with >= 2 pivots
  std::uint64_t macro_sccs = 0;       // threshold SCCs containing a macro arc
  std::uint64_t nontrivial_sccs = 0;  // threshold SCCs with >= 2 states
  std::uint64_t strict_differs = 0;   // thresholds where <= W and < W partitions differ
  std::uint64_t finite_theta_nonpivot = 0;
  std::uint64_t far_transfers = 0;  // far transfer arcs
  std::uint64_t near_transfers = 0;
};

// Partition labels of the states: SCCs of the arcs accepted by `keep`, by a transitive closure on
// bit rows (S <= 64 states).
template <class W, class Keep>
std::vector<int> scc_by_closure(const sssp::paper::RadialStateGraph<W>& rs, Keep keep) {
  const std::size_t s = rs.num_states();
  std::vector<std::uint64_t> reach(s, 0);
  for (std::size_t x = 0; x < s; ++x) reach[x] = std::uint64_t{1} << x;
  for (const auto& a : rs.graph.arcs) {
    if (keep(a.length)) reach[a.tail] |= std::uint64_t{1} << a.head;
  }
  for (std::size_t k = 0; k < s; ++k) {
    for (std::size_t i = 0; i < s; ++i) {
      if ((reach[i] >> k) & 1U) reach[i] |= reach[k];
    }
  }
  std::vector<int> label(s, -1);
  int next = 0;
  for (std::size_t x = 0; x < s; ++x) {
    if (label[x] >= 0) continue;
    for (std::size_t y = x; y < s; ++y) {
      if (((reach[x] >> y) & 1U) && ((reach[y] >> x) & 1U)) label[y] = next;
    }
    ++next;
  }
  return label;
}

struct Join {
  std::size_t a;
  std::size_t b;
};

// Union-find labels of vertices 0..n-1 over the joins whose weight is accepted by `keep`,
// restricted to the first `states` vertices.
template <class W, class Keep>
std::vector<std::size_t> join_labels(std::size_t n, std::size_t states,
                                     const std::vector<std::pair<Join, W>>& joins, Keep keep) {
  std::vector<std::size_t> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = i;
  const auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  for (const auto& [j, w] : joins) {
    if (keep(w)) p[find(j.a)] = find(j.b);
  }
  std::vector<std::size_t> out(states);
  for (std::size_t x = 0; x < states; ++x) out[x] = find(x);
  return out;
}

template <class A, class B>
bool same_partition(const std::vector<A>& a, const std::vector<B>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    for (std::size_t j = i + 1; j < a.size(); ++j) {
      if ((a[i] == a[j]) != (b[i] == b[j])) return false;
    }
  }
  return true;
}

// theta(z) = max(I*(z), O*(z)) from the definition (l.1115-1116): I* is the least, over macro
// arcs x -> y and admissible paths (every arc kind but macro) from y to z, of the largest arc
// length on the macro arc and the path; O* the same from z to a macro tail. Fixpoint iteration
// (Bellman-Ford style), independent of the library's block sweeps.
template <class W>
std::vector<W> theta_reference(const sssp::paper::RadialStateGraph<W>& rs) {
  const std::size_t s = rs.num_states();
  std::vector<W> in(s, kInf<W>);
  std::vector<W> out(s, kInf<W>);
  const auto& arcs = rs.graph.arcs;
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (rs.arc_kind[a] != RadialArcKind::kMacro) continue;
    in[arcs[a].head] = std::min(in[arcs[a].head], arcs[a].length);
    out[arcs[a].tail] = std::min(out[arcs[a].tail], arcs[a].length);
  }
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t a = 0; a < arcs.size(); ++a) {
      if (rs.arc_kind[a] == RadialArcKind::kMacro) continue;
      const auto& arc = arcs[a];
      if (in[arc.tail] != kInf<W>) {
        const W v = std::max(in[arc.tail], arc.length);
        if (v < in[arc.head]) {
          in[arc.head] = v;
          changed = true;
        }
      }
      if (out[arc.head] != kInf<W>) {
        const W v = std::max(out[arc.head], arc.length);
        if (v < out[arc.tail]) {
          out[arc.tail] = v;
          changed = true;
        }
      }
    }
  }
  std::vector<W> theta(s);
  for (std::size_t z = 0; z < s; ++z) theta[z] = std::max(in[z], out[z]);
  return theta;
}

// Checks Lemma 5.9 on one instance; std::nullopt if it holds or the input is outside the
// precondition (disconnected), otherwise a description of the first violation.
template <class W>
std::optional<std::string> lemma59_violation(const sh::Case<W>& c, Coverage* cov) {
  if (c.b == 0 || !sh::connected(c.n, c.edges)) return std::nullopt;
  const auto g = CsrGraph<W>::from_edges(c.n, c.edges);
  const auto od = sssp::paper::owner_decomposition(g, c.b);
  const auto rs = sssp::paper::local_radial_states(g, od);
  const auto j = sssp::paper::compute_joins(od, rs);
  const std::size_t s = rs.num_states();
  if (s > 64) return "more than 64 states: the closure oracle needs S <= 64";
  if (j.num_states != s) return "num_states differs from the radial state graph";
  const std::size_t np = od.pivots.size();
  const std::size_t nj = s + np;

  // The library's joins.
  std::vector<std::pair<Join, W>> lib;
  for (const auto& r : j.joins) {
    if (r.a >= nj || r.b >= nj) return "join endpoint out of range";
    lib.push_back({{r.a, r.b}, r.weight});
  }
  // Reference joins from the definition.
  const auto theta = theta_reference(rs);
  std::vector<std::pair<Join, W>> ref;
  const auto& arcs = rs.graph.arcs;
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    if (rs.arc_kind[a] == RadialArcKind::kWithinWindow) {
      ref.push_back({{arcs[a].tail, arcs[a].head}, arcs[a].length});
    }
  }
  for (std::size_t z = 0; z < s; ++z) {
    if (theta[z] != kInf<W>) ref.push_back({{z, s + od.owner[rs.projection[z]]}, theta[z]});
  }
  for (const auto& cor : od.corridors) ref.push_back({{s + cor.lower, s + cor.upper}, cor.weight});

  // Every W at which some partition can change, and one above all of them.
  std::set<W> values;
  for (const auto& a : arcs) {
    if (W{0} < a.length) values.insert(a.length);
  }
  for (const auto& [jn, w] : lib) values.insert(w);
  for (const auto& [jn, w] : ref) values.insert(w);
  values.erase(kInf<W>);
  W top{1};
  for (const W w : values) top = std::max(top, w);
  values.insert(top + top);

  if (cov != nullptr) {
    ++cov->instances;
    if (np >= 2) ++cov->multi_pivot;
    for (const auto k : rs.arc_kind) {
      cov->far_transfers += k == RadialArcKind::kFarTransfer ? 1 : 0;
      cov->near_transfers += k == RadialArcKind::kNearTransfer ? 1 : 0;
    }
    for (std::size_t z = 0; z < s; ++z) {
      if (theta[z] != kInf<W> && !od.is_pivot(rs.projection[z])) ++cov->finite_theta_nonpivot;
    }
  }
  for (const W w : values) {
    const auto le = [w](W x) { return !(w < x); };
    const auto lt = [w](W x) { return x < w; };
    const auto s_le = scc_by_closure(rs, le);
    const auto s_lt = scc_by_closure(rs, lt);
    const struct {
      const char* name;
      const std::vector<std::pair<Join, W>>* joins;
    } sources[] = {{"library joins", &lib}, {"reference joins", &ref}};
    for (const auto& src : sources) {
      const bool ok_le = same_partition(s_le, join_labels(nj, s, *src.joins, le));
      const bool ok_lt = same_partition(s_lt, join_labels(nj, s, *src.joins, lt));
      if (!ok_le || !ok_lt) {
        std::ostringstream os;
        os << src.name << ": components differ from the SCCs at W = " << w << " ("
           << (ok_le ? "" : "<= W ") << (ok_lt ? "" : "< W") << "), states " << s << ", pivots "
           << np;
        return os.str();
      }
    }
    if (cov != nullptr) {
      ++cov->thresholds;
      if (!same_partition(s_le, s_lt)) ++cov->strict_differs;
      std::map<int, std::size_t> size;
      for (const int l : s_le) ++size[l];
      for (const auto& [l, k] : size) cov->nontrivial_sccs += k > 1 ? 1 : 0;
      std::set<int> with_macro;
      for (std::size_t a = 0; a < arcs.size(); ++a) {
        if (rs.arc_kind[a] == RadialArcKind::kMacro && le(arcs[a].length) &&
            s_le[arcs[a].tail] == s_le[arcs[a].head]) {
          with_macro.insert(s_le[arcs[a].tail]);
        }
      }
      cov->macro_sccs += with_macro.size();
    }
  }
  return std::nullopt;
}

// Integer weights as W; for double the same values are quarter-integers (x / 4), so every
// comparison is exact and the weight pattern is the same.
template <class W>
std::vector<W> scaled(const std::vector<std::uint64_t>& xs) {
  std::vector<W> out;
  for (const auto x : xs) {
    if constexpr (std::is_floating_point_v<W>) {
      out.push_back(static_cast<W>(x) / W{4});
    } else {
      out.push_back(static_cast<W>(x));
    }
  }
  return out;
}

template <class W>
void run_exhaustive(std::size_t n, const std::vector<W>& weights, std::size_t max_b,
                    Coverage& cov) {
  std::uint64_t graphs = 0;
  sssp::test::for_each_small_graph<W>(n, weights, [&](const std::vector<Edge<W>>& edges) {
    ++graphs;
    if (!sh::connected(n, edges)) return;
    for (std::size_t b = 1; b <= max_b; ++b) {
      const sh::Case<W> c{n, edges, b};
      const auto v = lemma59_violation(c, &cov);
      if (v) {
        sh::Case<W> m = c;
        const auto r = sh::shrink_case(
            m, [](const sh::Case<W>& x) { return lemma59_violation<W>(x, nullptr); });
        FAIL("Lemma 5.9 fails on " << sh::describe(c) << ": " << *v
                                   << "\n  minimal: " << sh::describe(m) << ": " << r.message);
      }
    }
  });
  CHECK(graphs == sssp::test::count_small_graphs(n, weights.size()));
}

void require_coverage(const Coverage& cov) {
  WARN("instances " << cov.instances << ", thresholds " << cov.thresholds
                    << ", multi-pivot instances " << cov.multi_pivot << ", SCCs with a macro arc "
                    << cov.macro_sccs << ", nontrivial SCCs " << cov.nontrivial_sccs
                    << ", strict != non-strict " << cov.strict_differs
                    << ", finite theta at nonpivots " << cov.finite_theta_nonpivot
                    << ", far / near transfer arcs " << cov.far_transfers << " / "
                    << cov.near_transfers);
  CHECK(cov.multi_pivot > 0);
  CHECK(cov.macro_sccs > 0);
  CHECK(cov.nontrivial_sccs > 0);
  CHECK(cov.strict_differs > 0);
  CHECK(cov.finite_theta_nonpivot > 0);
  CHECK(cov.far_transfers > 0);
  CHECK(cov.near_transfers > 0);
}

}  // namespace

TEMPLATE_TEST_CASE("Lemma 5.9 holds on every connected graph with N <= 4 and every b",
                   "[property][lemma59]", std::uint64_t, double) {
  using W = TestType;
  // 1 and 2 give ties and two thresholds; 1500 lies beyond K^2 = 256 (b = 1) and K^2 = 1024
  // (b = 2, 3) times 1, so far transfers and several near components occur.
  Coverage cov;
  for (std::size_t n = 1; n <= 4; ++n) run_exhaustive<W>(n, scaled<W>({1, 2, 1500}), n + 1, cov);
  require_coverage(cov);
}

TEMPLATE_TEST_CASE("Lemma 5.9 holds on every connected graph with N = 4 and four weights [slow]",
                   "[property][lemma59][slow]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  run_exhaustive<W>(4, scaled<W>({1, 2, 3, 1500}), 5, cov);
  require_coverage(cov);
}

TEST_CASE("Lemma 5.9 holds on every connected graph with N = 5 and N = 6 [slow]",
          "[property][lemma59][slow]") {
  using W = std::uint64_t;
  Coverage cov;
  run_exhaustive<W>(5, {1, 1500}, 6, cov);  // every b; 3^10 graphs
  run_exhaustive<W>(6, {1}, 3, cov);        // all-equal weights (ties only); 2^15 graphs
  require_coverage(cov);
}

TEMPLATE_TEST_CASE("Lemma 5.9 on seeded random graphs with shrinking", "[property][lemma59]",
                   std::uint64_t, double) {
  using W = TestType;
  const auto weights = scaled<W>({1, 2, 3, 5, 40, 1500, 90000});
  const auto generate = [&](std::uint64_t seed) {
    sssp::test::SplitMix64 rng(seed);
    sh::Case<W> c;
    c.n = 6 + rng.below(4);
    c.b = 1 + rng.below(3);
    for (Vertex v = 1; v < c.n; ++v) {  // a random spanning tree keeps the graph connected
      c.edges.push_back({static_cast<Vertex>(rng.below(v)), v, weights[rng.below(weights.size())]});
    }
    const std::size_t extra = rng.below(2 * c.n);
    for (std::size_t i = 0; i < extra; ++i) {
      const auto u = static_cast<Vertex>(rng.below(c.n));
      const auto v = static_cast<Vertex>(rng.below(c.n));
      if (u != v) c.edges.push_back({u, v, weights[rng.below(weights.size())]});
    }
    return c;
  };
  Coverage cov;
  const auto out = sh::run_seeded<W>(
      1, 300, generate, [&](const sh::Case<W>& c) { return lemma59_violation(c, &cov); }, weights);
  INFO(out.report());
  CHECK_FALSE(out.failing_seed.has_value());
  CHECK(out.cases == 300);
  require_coverage(cov);
}
