// End-to-end overflow test of the Mode A SSSP paths (P4-OVERFLOW-E2E) on uint64 inputs whose
// auxiliary quantities leave 64 bits although the distances fit, and on the neighbouring inputs
// whose distances do not fit.
//
// Two paths are checked on every input and source:
//   * the public entry sssp::paper::paper_sssp (Alg 8, Sec 8.3; both variants of ADR 0008);
//   * Algorithm 4 (local_hierarchy, Sec 5.4) followed by Evaluate (Alg 3, Sec 4.3) with the
//     single initial offer 0 at the native entry of the source; the distance of v is the least
//     emission over the copies of v (Lemma 4.6, Lemma 5.10).
//
// Overflow rule (include/sssp/paper.hpp, decided 2026-11-27): std::overflow_error may be thrown
// only if a tentative path length d[u] + w of a reachable u reaches the sentinel infinity(); it
// must not be thrown because an auxiliary quantity (mass S(C), expiry, slot ends, root and
// forest distances, gaps) reaches it. Evaluate is stricter: it skips a relaxation value beyond
// the range and throws only if a state that received nothing else never emits (D-23). Hence,
// with d the exact distances:
//   * every d(v) < infinity()                  -> Evaluate returns d (no exception);
//   * no arc with d[u] + w >= infinity()       -> paper_sssp returns d (no exception);
//   * some d(v) >= infinity()                  -> both throw std::overflow_error;
//   * otherwise paper_sssp may throw std::overflow_error or return d (checked either way).
// The first point holds Evaluate to its own rule: a relaxation value beyond the range is a
// tentative length of a walk that is not shortest, which Evaluate skips. Each test case also
// states how many sources fall under each point (Tally, counted by hand), so that the exact
// checks are not vacuous: with the base-case Dijkstra of Alg 8 l.5, an arc back over a heavy
// edge (e.g. 1 -> 0 on one edge of 2^63: 2^63 + 2^63) already allows paper_sssp to throw.
//
// Exact reference: tests/python/reference/exact_sssp.py. On a path the distance is the sum of
// the edges between the endpoints, |pos(v) - pos(s)| with pos the prefix sums; this closed form
// is evaluated here in 128 bits and was compared with the reference output for every input and
// source below when the test was written (equal everywhere). Inputs that are not paths use the
// literal reference output (and the expected_dist of their regression records).
//
// Tree verification: paper_sssp returns its own tree; Evaluate returns emissions only, so the
// test builds a witness tree from the emitted distances (for each v an incident edge that is
// tight under them, from a strictly closer neighbour; all weights are positive) and passes the
// distances with that tree to verify_sssp, which checks the tree and the optimality certificate
// on every arc independently of both paths, and the exact reference.
//
// Inputs (all uint64):
//   p62    the path of four edges of 2^62, every source (2 is the reported one); sources 0 and
//          4 reach 2^64 (overflow), 1 - 3 fit (Sec 4.1 mass 5 * 2^62 and root distance 2^64 are
//          auxiliary).
//   q3     the minimised child-order counterexample 20261127-child-order-overflow, three edges
//          of q = ceil(2^64 / 5), every source; q5 the same with five edges (5 q = 2^64 + 4).
//   mass   regression 20261127-evaluate-mass-overflow (one edge of 2^63) and its neighbours:
//          two edges of 2^63 (2^64), 2^63 + (2^63 - 1) (exactly the sentinel), and
//          2^63 + (2^63 - 2) (the largest representable distance).
//   relax  regression 20261127-evaluate-relax-overflow, every source.
//   heavy  the path of k vertices whose first four edges are 2^62 and the rest 1. At k = 209953
//          (the first radial size of Alg 4, L = 18) with the paper's constants: [slow]. The fast
//          version runs k = 40 through the radial branch with test thresholds c = 1, e = 1 and
//          several separators (not Mode A constants; the construction is the same code).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/paper.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"

namespace {

using W = std::uint64_t;
using sssp::CsrGraph;
using sssp::Edge;
using sssp::Vertex;
__extension__ typedef unsigned __int128 U128;  // exact sums of 64-bit lengths (test only)

constexpr W kInf = sssp::WeightTraits<W>::infinity();
constexpr W kQ62 = W{1} << 62;
constexpr W k2p63 = W{1} << 63;
constexpr W kQ5 = 3689348814741910324ULL;  // ceil(2^64 / 5); 5 kQ5 = 2^64 + 4

std::string u128_str(U128 x) {
  if (x == 0) return "0";
  std::string s;
  while (x != 0) {
    s.insert(s.begin(), static_cast<char>('0' + static_cast<int>(x % 10)));
    x /= 10;
  }
  return s;
}

std::vector<Edge<W>> path_edges(const std::vector<W>& lengths) {
  std::vector<Edge<W>> e;
  e.reserve(lengths.size());
  for (Vertex v = 0; v < lengths.size(); ++v) e.push_back({v, v + 1, lengths[v]});
  return e;
}

std::vector<W> heavy_prefix_lengths(std::size_t k) {
  std::vector<W> len(k - 1, W{1});
  for (std::size_t i = 0; i < 4; ++i) len[i] = kQ62;
  return len;
}

/// Closed form on a path: d(s, v) = |pos(v) - pos(s)|, pos the 128-bit prefix sums.
std::vector<U128> path_distances(const std::vector<W>& lengths, Vertex s) {
  std::vector<U128> pos(lengths.size() + 1, 0);
  for (std::size_t i = 0; i < lengths.size(); ++i) pos[i + 1] = pos[i] + lengths[i];
  std::vector<U128> d(pos.size());
  for (std::size_t v = 0; v < pos.size(); ++v)
    d[v] = pos[v] >= pos[s] ? pos[v] - pos[s] : pos[s] - pos[v];
  return d;
}

std::vector<U128> literal(const std::vector<W>& d) {
  return {d.begin(), d.end()};
}

struct Classified {
  bool fits = true;            // every exact distance is below the sentinel
  bool tentative_hit = false;  // some arc u -> v has d[u] + w >= sentinel
};

Classified classify(const CsrGraph<W>& g, const std::vector<U128>& exact) {
  Classified c;
  for (Vertex u = 0; u < g.num_vertices(); ++u) {
    if (exact[u] >= kInf) c.fits = false;
    for (const auto& a : g.arcs(u)) {
      if (exact[u] + a.weight >= kInf) c.tentative_hit = true;
    }
  }
  return c;
}

std::vector<W> narrow(const std::vector<U128>& d) {
  std::vector<W> r(d.size());
  for (std::size_t i = 0; i < d.size(); ++i) r[i] = static_cast<W>(d[i]);
  return r;
}

/// Witness tree for distances `dist`: an incident edge (u, v) with dist[u] < dist[v] and
/// dist[v] - dist[u] == w; none if no such edge exists (verify_sssp then reports it).
sssp::SsspResult<W> with_witness_tree(const CsrGraph<W>& g, Vertex s, std::vector<W> dist) {
  sssp::SsspResult<W> r;
  const std::size_t n = g.num_vertices();
  r.parent.assign(n, sssp::kNoVertex);
  r.parent_edge.assign(n, sssp::kNoEdge);
  for (Vertex v = 0; v < n; ++v) {
    if (v == s || dist[v] == kInf) continue;
    for (const auto& a : g.arcs(v)) {
      if (dist[a.to] < dist[v] && dist[v] - dist[a.to] == a.weight) {
        r.parent[v] = a.to;
        r.parent_edge[v] = a.edge;
        break;
      }
    }
  }
  r.dist = std::move(dist);
  return r;
}

/// Alg 3 from the single offer 0 at the native entry of s; least emission per vertex.
std::vector<W> evaluate_from(const CsrGraph<W>& g, const sssp::paper::LocalHierarchyResult<W>& lh,
                             Vertex s) {
  sssp::paper::SeededOrderedHierarchy<W> seeded{
      lh.hierarchy, {{s, W{0}}}, {{lh.native_entry[s], W{0}}}};
  const auto valid = sssp::paper::validate_seeded_ordered_hierarchy(seeded, g);
  INFO(valid.summary());
  REQUIRE(valid.ok());
  const auto r = sssp::paper::evaluate(seeded);
  std::vector<W> best(g.num_vertices(), kInf);
  for (const auto& [x, v] : r.emissions) {
    const Vertex p = lh.hierarchy.projection[x];
    if (v < best[p]) best[p] = v;
  }
  return best;
}

void check_distances(const CsrGraph<W>& g, Vertex s, const std::vector<W>& got,
                     const std::vector<W>& expected) {
  std::size_t mismatches = 0;
  Vertex first_bad = sssp::kNoVertex;
  for (Vertex v = 0; v < g.num_vertices(); ++v) {
    if (got[v] != expected[v]) {
      if (mismatches++ == 0) first_bad = v;
    }
  }
  INFO("source " << s << ", first mismatching vertex " << first_bad);
  CHECK(mismatches == 0);
}

/// Which branch of the overflow rule each source was checked under (counted per source), so
/// that a test states how many sources are held to an exact result (see the file comment).
struct Tally {
  std::size_t eval_exact = 0;    // Evaluate must return the exact distances
  std::size_t eval_throw = 0;    // Evaluate must throw std::overflow_error
  std::size_t paper_strict = 0;  // paper_sssp must return the exact distances
  std::size_t paper_may = 0;     // paper_sssp may throw std::overflow_error, else exact
  std::size_t paper_throw = 0;   // paper_sssp must throw std::overflow_error

  friend bool operator==(const Tally&, const Tally&) = default;
};

std::ostream& operator<<(std::ostream& os, const Tally& t) {
  return os << "{eval_exact " << t.eval_exact << ", eval_throw " << t.eval_throw
            << ", paper_strict " << t.paper_strict << ", paper_may " << t.paper_may
            << ", paper_throw " << t.paper_throw << "}";
}

/// The Alg 4 -> Alg 3 path for source s against the exact distances.
void check_evaluate(const CsrGraph<W>& g, const sssp::paper::LocalHierarchyResult<W>& lh, Vertex s,
                    const std::vector<U128>& exact, Tally& tally) {
  const Classified c = classify(g, exact);
  INFO("Evaluate, source " << s << ", fits " << c.fits);
  if (!c.fits) {
    ++tally.eval_throw;
    CHECK_THROWS_AS(evaluate_from(g, lh, s), std::overflow_error);
    return;
  }
  ++tally.eval_exact;
  std::vector<W> got;
  REQUIRE_NOTHROW(got = evaluate_from(g, lh, s));
  const std::vector<W> expected = narrow(exact);
  check_distances(g, s, got, expected);
  const auto r = with_witness_tree(g, s, got);
  const auto report = sssp::verify_sssp(g, s, r, &expected);
  INFO(report.summary());
  CHECK(report.ok);
}

/// The public entry for source s against the exact distances, both variants.
void check_paper_sssp(const CsrGraph<W>& g, Vertex s, const std::vector<U128>& exact,
                      Tally& tally) {
  const Classified c = classify(g, exact);
  if (!c.fits) {
    ++tally.paper_throw;
  } else if (c.tentative_hit) {
    ++tally.paper_may;
  } else {
    ++tally.paper_strict;
  }
  for (const auto variant : sssp::paper::kPaperVariants) {
    INFO("paper_sssp " << sssp::paper::algorithm_name(variant) << ", source " << s << ", fits "
                       << c.fits << ", tentative reaches the sentinel " << c.tentative_hit);
    const sssp::paper::PaperOptions opt{variant};
    if (!c.fits) {
      CHECK_THROWS_AS(sssp::paper::paper_sssp(g, s, opt), std::overflow_error);
      continue;
    }
    sssp::SsspResult<W> r;
    if (c.tentative_hit) {
      // Allowed to throw (a tentative path length reaches the sentinel); if it returns, the
      // result must be exact.
      try {
        r = sssp::paper::paper_sssp(g, s, opt);
      } catch (const std::overflow_error&) {
        continue;
      }
    } else {
      REQUIRE_NOTHROW(r = sssp::paper::paper_sssp(g, s, opt));
    }
    const std::vector<W> expected = narrow(exact);
    check_distances(g, s, r.dist, expected);
    const auto report = sssp::verify_sssp(g, s, r, &expected);
    INFO(report.summary());
    CHECK(report.ok);
    CHECK(r.diag.fallbacks.size() == 1);  // ADR 0008: one base-case record per run
  }
}

/// Alg 4 with the paper's constants; checks the branch that the size selects.
sssp::paper::LocalHierarchyResult<W> mode_a_hierarchy(const CsrGraph<W>& g) {
  sssp::Diagnostics diag;
  sssp::paper::LocalHierarchyResult<W> lh;
  REQUIRE_NOTHROW(lh = sssp::paper::local_hierarchy(g, &diag));
  const bool radial =
      sssp::paper::local_construction_applies(g.num_vertices(), sssp::paper::kPaperLocalParams);
  if (radial) {
    CHECK(diag.fallbacks.empty());
  } else {
    REQUIRE(diag.fallbacks.size() == 1);
    CHECK(diag.fallbacks[0].reason ==
          sssp::paper::local_base_case_reason(sssp::paper::kPaperLocalParams));
  }
  return lh;
}

/// Both paths for every source in `sources` with exact(s) the reference distances.
template <class ExactFn>
Tally check_all(const CsrGraph<W>& g, const std::vector<Vertex>& sources, ExactFn exact) {
  const auto lh = mode_a_hierarchy(g);
  Tally tally;
  for (const Vertex s : sources) {
    const auto d = exact(s);
    check_evaluate(g, lh, s, d, tally);
    check_paper_sssp(g, s, d, tally);
  }
  return tally;
}

std::vector<Vertex> every_vertex(std::size_t n) {
  std::vector<Vertex> v(n);
  for (Vertex i = 0; i < n; ++i) v[i] = i;
  return v;
}

Tally check_path(const std::vector<W>& lengths) {
  const auto g = CsrGraph<W>::from_edges(lengths.size() + 1, path_edges(lengths));
  return check_all(g, every_vertex(g.num_vertices()),
                   [&](Vertex s) { return path_distances(lengths, s); });
}

}  // namespace

TEST_CASE("P4-OVERFLOW-E2E: the closed form of the path matches the exact reference",
          "[overflow][paper]") {
  // Spot values printed by exact_sssp.py (p62 from 0, q5 from 0, heavy k = 209953 from 2).
  const auto p62 = path_distances({kQ62, kQ62, kQ62, kQ62}, 0);
  CHECK(u128_str(p62[4]) == "18446744073709551616");
  const auto q5 = path_distances(std::vector<W>(5, kQ5), 0);
  CHECK(u128_str(q5[5]) == "18446744073709551620");
  const auto heavy = path_distances(heavy_prefix_lengths(209953), 2);
  CHECK(u128_str(heavy[209952]) == "9223372036854985756");
  CHECK(u128_str(heavy[0]) == "9223372036854775808");
}

// The expected tallies below are counted by hand from the exact distances (which sources have
// a distance >= 2^64 - 1, and which have an arc u -> v with d[u] + w >= 2^64 - 1); they state
// how many runs are held to the exact result, so that no branch is vacuous.

TEST_CASE("P4-OVERFLOW-E2E: path of four edges of 2^62, every source", "[overflow][paper]") {
  // Sources 1, 2, 3 fit (at most 3 * 2^62); 0 and 4 reach 2^64. From 1 (and 3) the arc
  // 4 -> 3 (0 -> 1) gives 3 * 2^62 + 2^62 = 2^64: only source 2 binds paper_sssp.
  CHECK(check_path({kQ62, kQ62, kQ62, kQ62}) == Tally{3, 2, 1, 2, 2});
}

TEST_CASE("P4-OVERFLOW-E2E: child-order counterexample, three edges of ceil(2^64 / 5)",
          "[overflow][paper]") {
  // Every d[u] + w is at most 3q + q = 4q < 2^64 - 1.
  CHECK(check_path(std::vector<W>(3, kQ5)) == Tally{4, 0, 4, 0, 0});
}

TEST_CASE("P4-OVERFLOW-E2E: five edges of ceil(2^64 / 5), the end-to-end distance is 2^64 + 4",
          "[overflow][paper]") {
  // Sources 0 and 5 reach 5q = 2^64 + 4; from 1 (4) the arc 5 -> 4 (0 -> 1) gives 5q.
  CHECK(check_path(std::vector<W>(5, kQ5)) == Tally{4, 2, 2, 2, 2});
}

TEST_CASE("P4-OVERFLOW-E2E: regression evaluate-mass-overflow and its boundary",
          "[overflow][paper][regression]") {
  // On each of these paths some arc back towards the source gives d[u] + w >= 2^64 - 1 for
  // every source, so paper_sssp may throw everywhere; Evaluate must not, where d fits.
  SECTION("one edge of 2^63 (the record)") {
    CHECK(check_path({k2p63}) == Tally{2, 0, 0, 2, 0});
  }
  SECTION("2^63 + 2^63 = 2^64 from an end") {
    CHECK(check_path({k2p63, k2p63}) == Tally{1, 2, 0, 1, 2});
  }
  SECTION("2^63 + (2^63 - 1) is exactly the sentinel from an end") {
    CHECK(check_path({k2p63, k2p63 - 1}) == Tally{1, 2, 0, 1, 2});
  }
  SECTION("2^63 + (2^63 - 2) = 2^64 - 2 is the largest finite distance") {
    CHECK(check_path({k2p63, k2p63 - 2}) == Tally{3, 0, 0, 3, 0});
  }
}

TEST_CASE("P4-OVERFLOW-E2E: regression evaluate-relax-overflow, every source",
          "[overflow][paper][regression]") {
  const auto g = CsrGraph<W>::from_edges(4, {{0, 1, 1}, {1, 2, kInf - 1}, {1, 3, 1}, {3, 2, 1}});
  // exact_sssp.py output per source (source 0 is the record's expected_dist). The heavy edge
  // gives d[u] + w >= 2^64 - 1 from every source (d(1) >= 1 or d(2) >= 2).
  const std::vector<std::vector<W>> ref{{0, 1, 3, 2}, {1, 0, 2, 1}, {3, 2, 0, 1}, {2, 1, 1, 0}};
  CHECK(check_all(g, every_vertex(4), [&](Vertex s) { return literal(ref[s]); }) ==
        Tally{4, 0, 0, 4, 0});
}

TEST_CASE("P4-OVERFLOW-E2E: heavy-prefix path k = 40 through the radial branch (test thresholds)",
          "[overflow][paper]") {
  // Alg 4 with c = 1, e = 1 (k = 40 > L_40 = 6) instead of the paper's 2 L^4, so that the
  // radial branch (Sec 5.1-5.4) runs at a small size; local_assembly_applies needs b < k.
  const std::size_t k = 40;
  const auto lengths = heavy_prefix_lengths(k);
  const auto g = CsrGraph<W>::from_edges(k, path_edges(lengths));
  const std::size_t b = GENERATE(2, 3, 5, 8, 13);
  sssp::paper::LocalHierarchyParams p;
  p.threshold_coeff = 1;
  p.threshold_exponent = 1;
  p.separator = b;
  REQUIRE(sssp::paper::local_construction_applies(k, p));
  sssp::Diagnostics diag;
  sssp::paper::LocalHierarchyResult<W> lh;
  INFO("separator " << b);
  REQUIRE_NOTHROW(lh = sssp::paper::local_hierarchy(g, &diag, p));
  CHECK(diag.fallbacks.empty());
#if SSSP_DIAGNOSTICS
  CHECK(diag.counters.at("local_hierarchy_radial") == 1);
#endif
  // Sources 1, 2, 3 fit (at most 3 * 2^62 + 35); every other source reaches 2^64 at vertex 0
  // or 39.
  Tally tally;
  for (Vertex s = 0; s < k; ++s) check_evaluate(g, lh, s, path_distances(lengths, s), tally);
  CHECK(tally == Tally{3, 37, 0, 0, 0});
}

TEST_CASE("P4-OVERFLOW-E2E: heavy-prefix path k = 40, public entry, every source",
          "[overflow][paper]") {
  // Base case of Alg 4 at k = 40. paper_sssp is bound at source 2 only (from 1 and 3 an arc
  // back over a heavy edge gives 3 * 2^62 + 2^62 = 2^64).
  CHECK(check_path(heavy_prefix_lengths(40)) == Tally{3, 37, 1, 2, 37});
}

TEST_CASE("P4-OVERFLOW-E2E: heavy-prefix path k = 209953, Mode A radial branch, Evaluate",
          "[slow][overflow][paper]") {
  const std::size_t k = 209953;
  const auto lengths = heavy_prefix_lengths(k);
  const auto g = CsrGraph<W>::from_edges(k, path_edges(lengths));
  REQUIRE(sssp::paper::local_construction_applies(k, sssp::paper::kPaperLocalParams));
  // Source 2 is the reported case (fits, largest distance 2^63 + 209948); source 3 fits as
  // well (paper_sssp may throw there: the arc 0 -> 1 gives 3 * 2^62 + 2^62); source 0 reaches
  // 2^64 at vertex 4 and beyond.
  CHECK(check_all(g, {2, 3, 0}, [&](Vertex s) { return path_distances(lengths, s); }) ==
        Tally{2, 1, 1, 1, 1});
}
