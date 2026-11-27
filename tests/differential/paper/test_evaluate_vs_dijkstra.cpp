// Differential tests of Evaluate (Alg 3 / Alg 2, Sec 4.3) against oracles that share no code
// with it (tests/support/paper/eval_oracle.hpp): delta of Def 4.2 from a physical Dijkstra, the
// protected set from a tight-arc BFS, dist_D from a state-graph Dijkstra, and the I1-I8 checker
// of tests/support/paper/eval_invariants.hpp attached through the src-side hook type.
//
// Claims checked (Lemma 4.6 and Def 4.2):
//   * every protected state emits exactly delta(pi(x)) (= dist_D(x) for protected states);
//   * every emission is >= delta(pi(x)); unprotected states are not compared with dist_D
//     (Q-29: they may emit less than dist_D through supplemental offers, or nothing);
//   * at every physical vertex served by a protected state the least emission is delta(v);
//   * the invariant checker reports no violation;
//   * supplemental offers on / off (D-14) give identical protected emissions.
// Weights: std::uint64_t (exact). The double runs use integer-valued lengths below 2^53, where
// every sum and every subtraction of Evaluate (origin rule) and of child_order (gaps, D-15) is
// exact, so the comparison is exact equality (Tolerance::exact()): a tolerance would only
// hide a real error here.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "sssp/verify.hpp"
#include "support/floyd_warshall.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/paper/seeded_instances.hpp"

using sssp::WeightTraits;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::EvalResult;
using sssp::paper::StateId;
using sssp::test::SplitMix64;
using namespace sssp::test::paper;

// Counters are compiled out with SSSP_DIAGNOSTICS=0 (release preset): assertions on them run
// only when they are counted.
constexpr bool kCounted = SSSP_DIAGNOSTICS != 0;
#define CHECK_COUNTER(...)    \
  do {                        \
    if constexpr (kCounted) { \
      CHECK(__VA_ARGS__);     \
    }                         \
  } while (false)

namespace {

using U = std::uint64_t;

template <class W>
struct Run {
  EvalResult<W> result;
  InvariantReport report;
};

template <class W>
Run<W> run_checked(const SeededInstance<W>& inst, const EvalReference<W>& ref, bool supplemental) {
  CheckerOptions copt;
  copt.supplemental_offers = supplemental;
  InvariantChecker<W> checker(inst.seeded, ref, copt);
  ForwardingObserver<EvalHooks<W>, W> adapter(checker);
  EvalOptions<W> opt;
  opt.supplemental_offers = supplemental;
  opt.hooks = &adapter;
  Run<W> r{sssp::paper::evaluate<W>(inst.seeded, opt), {}};
  r.report = checker.report();
  return r;
}

/// Lemma 4.6 against the oracle; returns the number of protected states.
template <class W>
std::size_t check_lemma46(const SeededInstance<W>& inst, const EvalReference<W>& ref,
                          const EvalResult<W>& r) {
  const auto& h = inst.seeded.hierarchy;
  const std::size_t n = h.num_states();
  REQUIRE(r.value.size() == n);
  std::size_t prot = 0;
  for (StateId x = 0; x < n; ++x) {
    INFO("state " << x);
    if (ref.is_protected[x]) {
      ++prot;
      REQUIRE(ref.state_distance[x] == ref.state_target[x]);
      REQUIRE(r.value[x] == ref.state_target[x]);
    }
    if (r.value[x] != WeightTraits<W>::infinity()) REQUIRE(r.value[x] >= ref.state_target[x]);
  }
  std::size_t emitted_protected = 0;
  for (const auto& [x, v] : r.emissions) {
    REQUIRE(r.value[x] == v);
    if (ref.is_protected[x]) ++emitted_protected;
  }
  REQUIRE(emitted_protected == ref.num_protected());
  const auto minima = physical_minima<W>(inst.num_physical, h.projection, r.emissions);
  const auto served = served_targets<W>(inst.num_physical, h.projection, ref);
  for (std::size_t v = 0; v < inst.num_physical; ++v) {
    INFO("physical vertex " << v);
    if (served[v]) REQUIRE(minima[v] == ref.target[v]);
  }
  return prot;
}

/// Counter totals over several runs (for branch coverage claims and the l.799-801 claim).
struct Totals {
  std::map<std::string, std::uint64_t> sum;
  void add(const sssp::Diagnostics& d) {
    for (const auto& [k, v] : d.counters) sum[k] += v;
  }
  [[nodiscard]] std::uint64_t operator[](const std::string& k) const {
    const auto it = sum.find(k);
    return it == sum.end() ? 0 : it->second;
  }
  /// Branch counters of Alg 2 / Alg 3, for the test log.
  [[nodiscard]] std::string branches() const {
    std::ostringstream os;
    for (const char* k :
         {"eval_opens", "eval_origin_second_rule", "eval_scans", "eval_clipped_slot_visits",
          "eval_resumed_slot_visits", "eval_failed_run_visits", "eval_registrations",
          "eval_unopened_not_due", "eval_leaf_failed_checks", "eval_supplemental_offers",
          "eval_supplemental_decreases", "eval_completed_nodes", "eval_expired_nodes",
          "eval_rejected_offers", "eval_root_leaf_emissions", "eval_roots_skipped",
          "eval_origin_rounding_fallbacks"}) {
      os << k << "=" << (*this)[k] << " ";
    }
    return os.str();
  }
};

template <class W>
void run_family(const std::vector<SeededInstance<W>>& insts, bool supplemental, Totals& totals) {
  for (const auto& inst : insts) {
    INFO(inst.name << " seed " << inst.seed << " supplemental " << supplemental);
    const auto val = sssp::paper::validate_seeded_ordered_hierarchy<W>(inst.seeded, inst.physical);
    INFO(val.summary());
    REQUIRE(val.ok());
    const auto ref = compute_eval_reference<W>(inst.seeded, inst.physical);
    const auto r = run_checked<W>(inst, ref, supplemental);
    INFO(r.report.summary());
    REQUIRE(r.report.ok());
    check_lemma46<W>(inst, ref, r.result);
    totals.add(r.result.diag);
  }
}

template <class W>
std::vector<SeededInstance<W>> base_cases(const std::string& suffix) {
  std::vector<SeededInstance<W>> out;
  for (const std::uint64_t seed : {1ULL, 2ULL, 3ULL, 4ULL, 5ULL}) {
    for (auto& inst : base_case_catalogue<W>(seed)) {
      const auto& nm = inst.name;
      if (nm.size() >= suffix.size() &&
          nm.compare(nm.size() - suffix.size(), suffix.size(), suffix) == 0) {
        out.push_back(std::move(inst));
      }
    }
  }
  return out;
}

template <class W>
std::vector<SeededInstance<W>> zero_arc_family(std::size_t count) {
  std::vector<SeededInstance<W>> out;
  for (std::size_t k = 0; k < count; ++k) {
    SplitMix64 rng(0x2e40ULL + k);
    ReferenceParams p;
    p.num_states = 2 + rng.below(25);
    p.num_physical = 1 + rng.below(p.num_states);
    p.num_arcs = rng.below(3 * p.num_states + 1);
    p.max_len = rng.below(2) == 0 ? 2 : 9;
    p.zero_arcs = 1 + rng.below(p.num_states + 1);
    p.num_sources = 1 + rng.below(3);
    p.extra_physical_edges = (k % 3 == 0) ? rng.below(p.num_physical + 1) : 0;
    p.extra_offers = rng.below(4);
    auto inst = zero_arc_instance<W>(p, 0x2e40ULL + k);
    inst.name = "zero_arcs/" + std::to_string(k);
    out.push_back(std::move(inst));
  }
  return out;
}

}  // namespace

TEST_CASE("evaluate vs dijkstra: BaseCaseSingleSource", "[evaluate][differential]") {
  Totals t;
  const auto insts = base_cases<U>("/single");
  REQUIRE(!insts.empty());
  run_family<U>(insts, true, t);
  for (const auto& inst : insts) {
    INFO(inst.name << " seed " << inst.seed);
    const auto ref = compute_eval_reference<U>(inst.seeded, inst.physical);
    const auto r = sssp::paper::evaluate<U>(inst.seeded);
    const auto& proj = inst.seeded.hierarchy.projection;
    const auto minima = physical_minima<U>(inst.num_physical, proj, r.emissions);
    // Every reachable vertex is served in the base case (its native state is protected).
    CHECK(minima == ref.target);
    const auto dj = sssp::baseline::dijkstra_binary_heap<U>(inst.physical, 0);
    CHECK(minima == dj.dist);
  }
}

TEST_CASE("evaluate vs dijkstra: MultiSource", "[evaluate][differential]") {
  Totals t;
  const auto bases = base_cases<U>("/multi");
  run_family<U>(bases, true, t);
  std::vector<SeededInstance<U>> refs;
  for (auto& inst : reference_catalogue<U>(7000, 300)) {
    if (inst.seeded.sources.size() > 1) refs.push_back(std::move(inst));
  }
  REQUIRE(refs.size() > 50);
  run_family<U>(refs, true, t);
  // The oracle's delta against min_i c_i + Floyd-Warshall on the physical graph.
  for (const std::vector<SeededInstance<U>>* fam : {&bases, &std::as_const(refs)}) {
    for (const auto& inst : *fam) {
      if (inst.num_physical > 32) continue;
      const auto fw = sssp::test::floyd_warshall<U>(inst.num_physical, inst.physical.edges());
      const auto delta = target_distances<U>(inst.physical, inst.seeded.sources);
      for (std::size_t v = 0; v < inst.num_physical; ++v) {
        U best = WeightTraits<U>::infinity();
        for (const auto& src : inst.seeded.sources) {
          const U d = fw[src.vertex][v];
          if (d != WeightTraits<U>::infinity() && src.offset + d < best) best = src.offset + d;
        }
        REQUIRE(delta[v] == best);
      }
    }
  }
}

TEST_CASE("evaluate vs dijkstra: RandomDirectedStateGraph", "[evaluate][differential]") {
  Totals t;
  run_family<U>(reference_catalogue<U>(1, 2000), true, t);
  // Branches reached by the catalogue (otherwise the test would not exercise them).
  CHECK_COUNTER(t["eval_origin_second_rule"] > 0);
  CHECK_COUNTER(t["eval_supplemental_decreases"] > 0);
  CHECK_COUNTER(t["eval_failed_run_visits"] > 0);
  CHECK_COUNTER(t["eval_root_leaf_emissions"] > 0);
  CHECK_COUNTER(t["eval_roots_skipped"] > 0);
  // Sec 4.3 l.799-801: no boundary cuts a local slot while events remain beyond it, so in
  // Sec 4 alone a clipped slot is never resumed in a later portion.
  CHECK_COUNTER(t["eval_resumed_slot_visits"] == 0);
  WARN("RandomDirectedStateGraph counter totals: " << t.branches());
}

TEST_CASE("evaluate vs dijkstra: RandomDirectedStateGraph with integer-valued doubles",
          "[evaluate][differential]") {
  Totals t;
  run_family<double>(reference_catalogue<double>(1, 600), true, t);
  // Same instances in both weight types give bit-identical emissions.
  const auto u = reference_catalogue<U>(1, 600);
  const auto d = reference_catalogue<double>(1, 600);
  for (std::size_t k = 0; k < u.size(); ++k) {
    const auto ru = sssp::paper::evaluate<U>(u[k].seeded);
    const auto rd = sssp::paper::evaluate<double>(d[k].seeded);
    REQUIRE(ru.value.size() == rd.value.size());
    for (std::size_t x = 0; x < ru.value.size(); ++x) {
      if (ru.value[x] == WeightTraits<U>::infinity()) {
        REQUIRE(rd.value[x] == WeightTraits<double>::infinity());
      } else {
        REQUIRE(rd.value[x] == static_cast<double>(ru.value[x]));
      }
    }
  }
}

TEST_CASE("evaluate vs dijkstra: SupplementalOffersOffEquivalent",
          "[evaluate][differential][d14]") {
  Totals on;
  Totals off;
  auto insts = reference_catalogue<U>(500, 1000);
  for (auto& b : base_cases<U>("/single")) insts.push_back(std::move(b));
  run_family<U>(insts, true, on);
  run_family<U>(insts, false, off);
  for (const auto& inst : insts) {
    INFO(inst.name << " seed " << inst.seed);
    const auto ref = compute_eval_reference<U>(inst.seeded, inst.physical);
    EvalOptions<U> no;
    no.supplemental_offers = false;
    const auto a = sssp::paper::evaluate<U>(inst.seeded);
    const auto b = sssp::paper::evaluate<U>(inst.seeded, no);
    for (StateId x = 0; x < inst.seeded.hierarchy.num_states(); ++x) {
      if (ref.is_protected[x]) REQUIRE(a.value[x] == b.value[x]);
    }
  }
  CHECK_COUNTER(off["eval_supplemental_offers"] == 0);
  CHECK_COUNTER(on["eval_supplemental_offers"] > 0);
  CHECK_COUNTER(off["eval_expired_nodes"] > 0);  // without the offers, unreached copies expire
  WARN("supplemental on: " << on.branches());
  WARN("supplemental off: " << off.branches());
}

TEST_CASE("evaluate vs dijkstra: ZeroLengthArcsAcyclic", "[evaluate][differential]") {
  Totals t;
  const auto insts = zero_arc_family<U>(600);
  std::size_t zero_arcs = 0;
  for (const auto& inst : insts) {
    for (const auto& a : inst.seeded.hierarchy.graph.arcs) zero_arcs += a.length == 0 ? 1 : 0;
  }
  REQUIRE(zero_arcs > 0);
  run_family<U>(insts, true, t);
  run_family<U>(insts, false, t);
}

TEST_CASE("evaluate vs dijkstra: I1-I8 hold on 10^4 random hierarchies with N <= 60",
          "[evaluate][differential][invariants]") {
  Totals t;
  std::size_t protected_total = 0;
  for (std::size_t k = 0; k < 10000; ++k) {
    SplitMix64 rng(0x4e7aULL + k);
    ReferenceParams p;
    p.num_states = 2 + rng.below(59);
    p.num_physical = 1 + rng.below(p.num_states);
    p.num_arcs = rng.below(4 * p.num_states + 1);
    p.max_len = 1 + rng.below(16);
    p.num_sources = 1 + rng.below(3);
    p.extra_physical_edges = rng.below(2) == 0 ? 0 : rng.below(p.num_physical + 1);
    p.extra_offers = rng.below(5);
    const auto inst = reference_instance<U>(p, 0x4e7aULL + k);
    INFO("instance " << k);
    const auto ref = compute_eval_reference<U>(inst.seeded, inst.physical);
    const bool supplemental = (k % 4) != 3;
    const auto r = run_checked<U>(inst, ref, supplemental);
    INFO(r.report.summary());
    REQUIRE(r.report.ok());
    protected_total += check_lemma46<U>(inst, ref, r.result);
    t.add(r.result.diag);
  }
  CHECK(protected_total > 0);
  CHECK_COUNTER(t["eval_resumed_slot_visits"] == 0);
  WARN("10^4 random hierarchies, counter totals: " << t.branches());
}

namespace {

/// The instance with every length multiplied by f in double (one rounding per length). Rounding
/// is monotone, so Def 4.1 still holds: each link and arc of the reference hierarchy is a
/// physical edge with the identical rounded value, and link <= height is preserved.
SeededInstance<double> scaled(const SeededInstance<U>& u, double f) {
  const auto sc = [f](U w) { return static_cast<double>(w) * f; };
  SeededInstance<double> d;
  d.name = u.name + "/x" + std::to_string(f);
  d.seed = u.seed;
  d.num_physical = u.num_physical;
  std::vector<sssp::Edge<double>> edges;
  for (const auto& e : u.physical.edges()) edges.push_back({e.u, e.v, sc(e.weight)});
  d.physical = sssp::CsrGraph<double>::from_edges(u.num_physical, std::move(edges));
  const auto& hu = u.seeded.hierarchy;
  auto& hd = d.seeded.hierarchy;
  hd.graph.num_states = hu.graph.num_states;
  for (const auto& a : hu.graph.arcs) hd.graph.arcs.push_back({a.tail, a.head, sc(a.length)});
  hd.projection = hu.projection;
  for (const auto& nd : hu.nodes) {
    sssp::paper::HierarchyNode<double> m;
    m.parent = nd.parent;
    m.children = nd.children;
    m.state = nd.state;
    m.height = sc(nd.height);
    for (const auto& lk : nd.links) m.links.push_back({lk.x, lk.y, sc(lk.length), lk.arc});
    hd.nodes.push_back(std::move(m));
  }
  hd.leaf_of_state = hu.leaf_of_state;
  hd.roots = hu.roots;
  hd.c0 = hu.c0;
  for (const auto& s : u.seeded.sources) d.seeded.sources.push_back({s.vertex, sc(s.offset)});
  for (const auto& o : u.seeded.offers) d.seeded.offers.push_back({o.state, sc(o.value)});
  return d;
}

}  // namespace

TEST_CASE("evaluate vs dijkstra: non-integer doubles stay within a rounding tolerance (D-15)",
          "[evaluate][differential][double]") {
  // Lengths k * 0.1 are not representable, so the sums of Evaluate, the masses S, the gaps g of
  // order_children (materialised from a signed pair, D-15) and the origin b - j nu (D-16) are
  // rounded. Reference: the EXACT integer instance (protected set and delta in uint64), scaled.
  // Tolerance (test-only, documented, deliberately not tight): every real Evaluate compares is
  // a sum, or for g a difference of two sums, of at most K = 4 (N + M + L + H) rounded terms
  // (N states, M arcs, L links, H nodes: a key is an offer plus a simple path, plus
  // supplemental increments S + g whose tour charges use each link and arc at most twice,
  // Lemma 4.4), so its error is at most gamma_K (Tolerance::for_path_hops(K)) times the sum of
  // the magnitudes of its terms, which is at most 4 * total + max delta, total = the sum of all
  // arc, link and height values and source offsets.
  std::size_t compared = 0;
  Totals t;
  for (const auto& u : reference_catalogue<U>(3000, 400)) {
    for (const double f : {0.1, 0.3}) {
      const auto d = scaled(u, f);
      INFO(d.name << " seed " << d.seed);
      const auto ref = compute_eval_reference<U>(u.seeded, u.physical);  // exact
      const auto r = sssp::paper::evaluate<double>(d.seeded);
      t.add(r.diag);
      const auto& h = u.seeded.hierarchy;
      double total = 0;
      std::size_t links = 0;
      for (const auto& a : d.seeded.hierarchy.graph.arcs) total += a.length;
      for (const auto& nd : d.seeded.hierarchy.nodes) {
        total += nd.height;
        for (const auto& lk : nd.links) total += lk.length;
        links += nd.links.size();
      }
      for (const auto& s : d.seeded.sources) total += s.offset;
      double max_delta = 0;
      for (StateId x = 0; x < h.num_states(); ++x) {
        if (ref.state_target[x] != WeightTraits<U>::infinity()) {
          max_delta = std::max(max_delta, static_cast<double>(ref.state_target[x]) * f);
        }
      }
      const auto k = 4 * (h.num_states() + h.graph.arcs.size() + links + h.nodes.size());
      const double slack = sssp::Tolerance::for_path_hops(k).slack(4.0 * total + max_delta);
      for (StateId x = 0; x < h.num_states(); ++x) {
        if (ref.state_target[x] == WeightTraits<U>::infinity()) continue;
        const double delta = static_cast<double>(ref.state_target[x]) * f;
        INFO("state " << x << " value " << r.value[x] << " delta " << delta << " slack " << slack);
        if (ref.is_protected[x]) {
          ++compared;
          REQUIRE(r.value[x] != WeightTraits<double>::infinity());
          REQUIRE(std::abs(r.value[x] - delta) <= slack);
        }
        if (r.value[x] != WeightTraits<double>::infinity()) REQUIRE(r.value[x] >= delta - slack);
      }
    }
  }
  CHECK(compared > 0);
  WARN("non-integer doubles: " << compared << " protected states compared; " << t.branches());
}
