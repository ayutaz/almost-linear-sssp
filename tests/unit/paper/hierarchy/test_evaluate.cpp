// Unit tests of Evaluate (Alg 3) and Scan (Alg 2), Sec 4.3 of arXiv 2611.34911 v1, on hand
// examples whose traces were derived by hand from the paper (expected values are NOT taken
// from an implementation run).
//
// Hand derivations (Sec 4.1 normalisation, Sec 4.2 child order, then Alg 3 / Alg 2):
//
// chain: states a = 0, b = 1, c = 2; arcs 0 a->b 1, 1 b->a 1, 2 b->c 4, 3 c->b 4; pi = id on
//   the physical path 0 -1- 1 -4- 2; source (0, 0), offer a = 0. Threshold-SCC hierarchy:
//   P = {a, b} (h 1), R = {P, c} (h 4), links (a, b, 1) and (b, c, 4).
//   Sec 4.1: S(P) = 1 + 1 = 2, S(R) = 4 + 4 + 2 = 10; 4 S(P) >= h(R) links P and R, table
//   1, 2, 4: nu(P) = 1, nu(R) = 4. No arc is shorter than the unit at its LCA: no light arcs.
//   Sec 4.2: F = a - b - c, preorder a, b, c. At P the DFS roots are a then b, both finish
//   at once, marks 0 and D_F(a, b) = 1; the rank order is the REVERSE finishing order [b, a],
//   g(b, a) = 1. At R: provisional reps a (for P) and c; finishing [P, c], ranks [c, P],
//   g(c, P) = D_F(a, c) = 5, rep(P) = a. Leaf layout c, b, a.
//   Alg 3: R opens at D0 = 0, origin 0, expiry 10. Slot [0, 4): the run {c, P} has minimum 0
//   at a, child P: left subrun {c}, supplemental offer at c = max{4, 0 + S(P) + 5} = 7. The run
//   {c} (7) is not due; P opens at 0, 0 + S(P) = 2 < 4 so origin 0 (first rule), deadline
//   min(10, 2) = 2. P slot [0, 1): run {b, a} splits at a, supplemental at b =
//   max{1, 0 + 0 + 1} = 1; {b} is not due; a emits 0, relaxes a->b (1). P slot [1, 2): {b}
//   splits, b emits 1, relaxes b->a (2, settled) and b->c (5). P completed: retired. R slot
//   [4, 8): {c} (5) splits, c emits 5, relaxes c->b (9). R completed.
//
// align: as chain with lengths 3 and 8 and the source (2, 0) at c (offer c = 0).
//   S(P) = 6, S(R) = 22, table 3, 6, 12: nu(P) = 3, nu(R) = 6; ranks [c, P] with
//   g(c, P) = 3 + 8 = 11 and [b, a] with g(b, a) = 3; layout c, b, a.
//   R slot [0, 6): run {c, P} minimum 0 at c (rank 0): right subrun {P}, supplemental at
//   rep(P) = a: max{6, 0 + 0 + 11} = 11; c emits 0, relaxes c->b (8). {P} (8) not due.
//   R slot [6, 12): {P} splits at b; P opens at D0 = 8, 8 + 6 = 14 >= 12, so the second origin
//   rule: largest 12 - 3 j <= 8 is 6 (j = 2); deadline min(22, 14) = 14. P slot [6, 9): run
//   {b, a} minimum 8 at b: right subrun {a}, supplemental at a = max{9, 8 + 0 + 3} = 11; b
//   emits 8, relaxes b->a (11) and b->c (16, settled). P slot [9, 12): {a} emits 11, relaxes
//   a->b (14, settled). P completed and retired; R completed.
//
// pair: two states without arcs under one node R (h 1, link (a, b, 1)), physical edge 0 -1- 1,
//   offer a = 0. S(R) = 2, nu(R) = 1; ranks [b, a], g = 1, layout b, a. b is not protected
//   (no state path reaches it). With supplemental offers: slot [0, 1) splits at a,
//   supplemental at b = max{1, 0 + 0 + 1} = 1, a emits 0; slot [1, 2): b emits 1. Without them
//   (D-14): b keeps key infinity, the slots [1, 2) and [2, 3) (2 <= expiry 2) pass, R retires
//   with b unemitted.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"

using sssp::CsrGraph;
using sssp::WeightTraits;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::evaluate;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
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
struct Fixture {
  SeededOrderedHierarchy<W> s;
  CsrGraph<W> phys;
};

template <class W>
Fixture<W> path_fixture(W ab, W bc, sssp::Vertex source_vertex, StateId offer_state) {
  StateGraph<W> d;
  d.num_states = 3;
  d.arcs = {{0, 1, ab}, {1, 0, ab}, {1, 2, bc}, {2, 1, bc}};
  Fixture<W> f;
  f.s.hierarchy = sssp::paper::build_threshold_scc_hierarchy<W>(d, {0, 1, 2});
  f.phys = sssp::paper::physical_graph_of_arcs<W>(f.s.hierarchy.graph, {0, 1, 2}, 3);
  f.s.sources = {{source_vertex, W{0}}};
  f.s.offers = {{offer_state, W{0}}};
  return f;
}

template <class W>
Fixture<W> chain() {
  return path_fixture<W>(W{1}, W{4}, 0, 0);
}
template <class W>
Fixture<W> align() {
  return path_fixture<W>(W{3}, W{8}, 2, 2);
}

/// Two arc-free states under one node of height 1 (see the file comment).
template <class W>
Fixture<W> pair_fixture() {
  StateGraph<W> d;
  d.num_states = 2;
  auto h = sssp::paper::make_leaf_forest<W>(d, {0, 1});
  sssp::paper::add_internal_node<W>(h, W{1}, {0, 1}, {{0, 1, W{1}}});
  sssp::paper::collect_roots(h);
  Fixture<W> f;
  f.s.hierarchy = std::move(h);
  f.phys = CsrGraph<W>::from_edges(2, {{0, 1, W{1}}});
  f.s.sources = {{0, W{0}}};
  f.s.offers = {{0, W{0}}};
  return f;
}

/// Runs evaluate with a RecordingObserver attached through the src-side hook type.
template <class W>
struct Recorded {
  sssp::paper::EvalResult<W> result;
  std::vector<EvalEvent<W>> events;
};

template <class W>
Recorded<W> record(const SeededOrderedHierarchy<W>& s, bool supplemental = true) {
  RecordingObserver<W> rec;
  ForwardingObserver<EvalHooks<W>, W> adapter(rec);
  EvalOptions<W> opt;
  opt.supplemental_offers = supplemental;
  opt.hooks = &adapter;
  Recorded<W> r{evaluate<W>(s, opt), {}};
  r.events = rec.events;
  return r;
}

/// Runs evaluate with the I1-I8 checker attached and returns its report.
template <class W>
InvariantReport checked(const Fixture<W>& f, bool supplemental = true) {
  const auto ref = compute_eval_reference<W>(f.s, f.phys);
  CheckerOptions copt;
  copt.supplemental_offers = supplemental;
  InvariantChecker<W> checker(f.s, ref, copt);
  ForwardingObserver<EvalHooks<W>, W> adapter(checker);
  EvalOptions<W> opt;
  opt.supplemental_offers = supplemental;
  opt.hooks = &adapter;
  (void)evaluate<W>(f.s, opt);
  return checker.report();
}

const char* kind_name(EvalEventKind k) {
  static constexpr const char* kNames[] = {"initial", "relax",      "supplemental",
                                           "open",    "scan_begin", "scan_end",
                                           "emit",    "retire",     "finish"};
  return kNames[static_cast<std::size_t>(k)];
}

/// Event equality ignoring the node id (ids of the normalised forest are not hand-derived).
template <class W>
bool same_event(const EvalEvent<W>& a, const EvalEvent<W>& b) {
  if (a.kind != b.kind) return false;
  switch (a.kind) {
    case EvalEventKind::kInitialOffer:
    case EvalEventKind::kSupplementalOffer:
    case EvalEventKind::kEmit:
      return a.state == b.state && a.value == b.value;
    case EvalEventKind::kRelaxOffer:
      return a.arc == b.arc && a.value == b.value;
    case EvalEventKind::kOpen:
      return a.states == b.states && a.value == b.value && a.origin == b.origin &&
             a.unit == b.unit && a.mass == b.mass;
    case EvalEventKind::kScanBegin:
      return a.states == b.states && a.lo == b.lo && a.hi == b.hi && a.slot_lo == b.slot_lo &&
             a.slot_hi == b.slot_hi && a.deadline == b.deadline;
    case EvalEventKind::kScanEnd:
    case EvalEventKind::kRetire:
      return a.states == b.states;
    case EvalEventKind::kFinish:
      return true;
  }
  return false;
}

template <class W>
void require_trace(const std::vector<EvalEvent<W>>& got, const std::vector<EvalEvent<W>>& want) {
  const std::size_t n = std::min(got.size(), want.size());
  for (std::size_t i = 0; i < n; ++i) {
    INFO("event " << i << ": got " << kind_name(got[i].kind) << ", want "
                  << kind_name(want[i].kind));
    REQUIRE(same_event(got[i], want[i]));
  }
  REQUIRE(got.size() == want.size());
}

// States in leaf-layout order (hand-derived: c, b, a).
const std::vector<StateId> kR{2, 1, 0};
const std::vector<StateId> kP{1, 0};

template <class W>
std::vector<EvalEvent<W>> chain_trace() {
  return {
      ev_initial<W>(0, 0),
      ev_open<W>(0, kR, 0, 0, 4, 10),
      ev_scan_begin<W>(0, kR, 0, 4, 0, 4, 10),
      ev_supplemental<W>(2, 7),       // left subrun {c}: max{4, 0 + S(P) + g(c, P)} = max{4, 2 + 5}
      ev_open<W>(0, kP, 0, 0, 1, 2),  // first origin rule: 0 + 2 < 4
      ev_scan_begin<W>(0, kP, 0, 1, 0, 1, 2),
      ev_supplemental<W>(1, 1),  // left subrun {b}: max{1, 0 + 0 + 1}
      ev_emit<W>(0, 0),
      ev_relax<W>(0, 1),
      ev_scan_end<W>(0, kP),
      ev_scan_begin<W>(0, kP, 1, 2, 1, 2, 2),
      ev_emit<W>(1, 1),
      ev_relax<W>(1, 2),
      ev_relax<W>(2, 5),
      ev_scan_end<W>(0, kP),
      ev_retire<W>(0, kP),
      ev_scan_end<W>(0, kR),
      ev_scan_begin<W>(0, kR, 4, 8, 4, 8, 10),
      ev_emit<W>(2, 5),
      ev_relax<W>(3, 9),
      ev_scan_end<W>(0, kR),
      ev_finish<W>(),
  };
}

template <class W>
std::vector<EvalEvent<W>> align_trace() {
  return {
      ev_initial<W>(2, 0),
      ev_open<W>(0, kR, 0, 0, 6, 22),
      ev_scan_begin<W>(0, kR, 0, 6, 0, 6, 22),
      ev_supplemental<W>(0, 11),  // right subrun {P}, at rep(P) = a: max{6, 0 + 0 + 11}
      ev_emit<W>(2, 0),
      ev_relax<W>(3, 8),
      ev_scan_end<W>(0, kR),
      ev_scan_begin<W>(0, kR, 6, 12, 6, 12, 22),
      ev_open<W>(0, kP, 8, 6, 3, 6),  // second origin rule: 8 + 6 >= 12, 12 - 2 * 3 = 6
      ev_scan_begin<W>(0, kP, 6, 9, 6, 9, 14),
      ev_supplemental<W>(0, 11),  // right subrun {a}: max{9, 8 + 0 + 3}
      ev_emit<W>(1, 8),
      ev_relax<W>(1, 11),
      ev_relax<W>(2, 16),
      ev_scan_end<W>(0, kP),
      ev_scan_begin<W>(0, kP, 9, 12, 9, 12, 14),
      ev_emit<W>(0, 11),
      ev_relax<W>(0, 14),
      ev_scan_end<W>(0, kP),
      ev_retire<W>(0, kP),
      ev_scan_end<W>(0, kR),
      ev_finish<W>(),
  };
}

const std::vector<StateId> kPair{1, 0};

template <class W>
std::vector<EvalEvent<W>> pair_trace_with_supplemental() {
  return {
      ev_initial<W>(0, 0),
      ev_open<W>(0, kPair, 0, 0, 1, 2),
      ev_scan_begin<W>(0, kPair, 0, 1, 0, 1, 2),
      ev_supplemental<W>(1, 1),
      ev_emit<W>(0, 0),
      ev_scan_end<W>(0, kPair),
      ev_scan_begin<W>(0, kPair, 1, 2, 1, 2, 2),
      ev_emit<W>(1, 1),
      ev_scan_end<W>(0, kPair),
      ev_finish<W>(),
  };
}

template <class W>
std::vector<EvalEvent<W>> pair_trace_without_supplemental() {
  return {
      ev_initial<W>(0, 0),
      ev_open<W>(0, kPair, 0, 0, 1, 2),
      ev_scan_begin<W>(0, kPair, 0, 1, 0, 1, 2),
      ev_emit<W>(0, 0),
      ev_scan_end<W>(0, kPair),
      ev_scan_begin<W>(0, kPair, 1, 2, 1, 2, 2),
      ev_scan_end<W>(0, kPair),
      ev_scan_begin<W>(0, kPair, 2, 3, 2, 3, 2),  // slot containing the expiry 2 (inclusive)
      ev_scan_end<W>(0, kPair),
      ev_retire<W>(0, kPair),  // root expired with b unemitted
      ev_finish<W>(),
  };
}

[[maybe_unused]] U counter(const sssp::Diagnostics& d, const std::string& name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

}  // namespace

TEST_CASE("evaluate: chain hand example, first origin rule and both split sides",
          "[evaluate][alg2][alg3]") {
  const auto f = chain<U>();
  const auto r = record<U>(f.s);
  require_trace<U>(r.events, chain_trace<U>());
  CHECK(r.result.value == std::vector<U>{0, 1, 5});
  CHECK(checked<U>(f).summary() == "OK");
}

TEST_CASE("evaluate: align hand example, second origin rule b - j nu (Sec 4.3 l.794)",
          "[evaluate][alg2][alg3]") {
  const auto f = align<U>();
  const auto r = record<U>(f.s);
  require_trace<U>(r.events, align_trace<U>());
  CHECK(r.result.value == std::vector<U>{11, 8, 0});
  CHECK(checked<U>(f).summary() == "OK");
}

TEST_CASE("evaluate: integer-valued doubles give the same hand traces", "[evaluate]") {
  const auto c = record<double>(chain<double>().s);
  require_trace<double>(c.events, chain_trace<double>());
  CHECK(c.result.value == std::vector<double>{0, 1, 5});
  const auto a = record<double>(align<double>().s);
  require_trace<double>(a.events, align_trace<double>());
  CHECK(checked<double>(align<double>()).summary() == "OK");
}

TEST_CASE("evaluate: hand-counted Lemma 4.7 counters on the chain example", "[evaluate][lemma47]") {
  const auto r = evaluate<U>(chain<U>().s);
  const auto& d = r.diag;
  CHECK_COUNTER(counter(d, "eval_calls") == 1);
  CHECK_COUNTER(counter(d, "eval_states") == 3);
  CHECK_COUNTER(counter(d, "eval_arcs") == 4);
  CHECK_COUNTER(counter(d, "eval_initial_offers") == 1);
  CHECK_COUNTER(counter(d, "eval_opens") == 2);
  CHECK_COUNTER(counter(d, "eval_root_slots") == 2);
  CHECK_COUNTER(counter(d, "eval_child_slot_visits") == 2);
  CHECK_COUNTER(counter(d, "eval_scans") == 4);
  CHECK_COUNTER(counter(d, "eval_run_visits") == 6);
  CHECK_COUNTER(counter(d, "eval_failed_run_visits") == 2);
  CHECK_COUNTER(counter(d, "eval_registrations") == 4);
  CHECK_COUNTER(counter(d, "eval_interval_splits") == 2);
  CHECK_COUNTER(counter(d, "eval_child_visits") == 4);
  CHECK_COUNTER(counter(d, "eval_leaf_failed_checks") == 0);
  CHECK_COUNTER(counter(d, "eval_emissions") == 3);
  CHECK_COUNTER(counter(d, "eval_relaxations") == 4);
  CHECK_COUNTER(counter(d, "eval_relax_decreases") == 1);
  CHECK_COUNTER(counter(d, "eval_rejected_offers") == 2);
  CHECK_COUNTER(counter(d, "eval_supplemental_offers") == 2);
  CHECK_COUNTER(counter(d, "eval_supplemental_decreases") == 2);
  CHECK_COUNTER(counter(d, "eval_completed_nodes") == 2);
  CHECK_COUNTER(counter(d, "eval_expired_nodes") == 0);
  CHECK_COUNTER(counter(d, "eval_origin_second_rule") == 0);
  CHECK_COUNTER(counter(d, "eval_clipped_slot_visits") == 0);
  CHECK_COUNTER(counter(d, "eval_resumed_slot_visits") == 0);
  CHECK(r.diag.fallbacks.empty());

  const auto a = evaluate<U>(align<U>().s);
  CHECK_COUNTER(counter(a.diag, "eval_origin_second_rule") == 1);
  CHECK_COUNTER(counter(a.diag, "eval_origin_steps") == 2);  // 12 -> 9 -> 6
}

TEST_CASE(
    "evaluate: an unprotected state emits a value >= delta or nothing; without supplemental "
    "offers the root expires (D-14, Q-29)",
    "[evaluate][alg2][d14]") {
  const auto f = pair_fixture<U>();
  REQUIRE(sssp::paper::validate_seeded_ordered_hierarchy<U>(f.s, f.phys).ok());

  const auto with = record<U>(f.s, true);
  require_trace<U>(with.events, pair_trace_with_supplemental<U>());
  CHECK(with.result.value == std::vector<U>{0, 1});
  CHECK(checked<U>(f, true).summary() == "OK");

  const auto without = record<U>(f.s, false);
  require_trace<U>(without.events, pair_trace_without_supplemental<U>());
  CHECK(without.result.value == std::vector<U>{0, WeightTraits<U>::infinity()});
  CHECK_COUNTER(counter(without.result.diag, "eval_expired_nodes") == 1);
  CHECK_COUNTER(counter(without.result.diag, "eval_supplemental_offers") == 0);
  CHECK(checked<U>(f, false).summary() == "OK");
}

TEST_CASE(
    "evaluate: one-state roots emit in root-quotient order, roots without an offer are "
    "skipped (Alg 3 l.5-7)",
    "[evaluate][alg3]") {
  StateGraph<U> d;
  d.num_states = 4;  // state 3 has no arcs and no offer
  d.arcs = {{0, 1, 2}, {1, 2, 3}, {0, 2, 7}};
  Fixture<U> f;
  f.s.hierarchy = sssp::paper::build_threshold_scc_hierarchy<U>(d, {0, 1, 2, 3});
  f.phys = sssp::paper::physical_graph_of_arcs<U>(f.s.hierarchy.graph, {0, 1, 2, 3}, 4);
  f.s.sources = {{0, 0}};
  f.s.offers = {{0, 0}};
  const auto r = record<U>(f.s);
  require_trace<U>(r.events, {
                                 ev_initial<U>(0, 0),
                                 ev_emit<U>(0, 0),
                                 ev_relax<U>(0, 2),
                                 ev_relax<U>(2, 7),
                                 ev_emit<U>(1, 2),
                                 ev_relax<U>(1, 5),
                                 ev_emit<U>(2, 5),
                                 ev_finish<U>(),
                             });
  CHECK(r.result.value == std::vector<U>{0, 2, 5, WeightTraits<U>::infinity()});
  CHECK_COUNTER(counter(r.result.diag, "eval_roots_skipped") == 1);
  CHECK(checked<U>(f).summary() == "OK");
}

TEST_CASE("evaluate: only the least initial offer per state is posted (Lemma 4.7 input scan)",
          "[evaluate][alg3]") {
  auto f = chain<U>();
  f.s.offers = {{2, 9}, {0, 3}, {0, 0}};
  const auto r = record<U>(f.s);
  std::vector<EvalEvent<U>> initial;
  for (const auto& e : r.events) {
    if (e.kind == EvalEventKind::kInitialOffer) initial.push_back(e);
  }
  require_trace<U>(initial, {ev_initial<U>(0, 0), ev_initial<U>(2, 9)});
  CHECK(r.result.value == std::vector<U>{0, 1, 5});
  CHECK_COUNTER(counter(r.result.diag, "eval_initial_offers") == 3);
  CHECK(checked<U>(f).summary() == "OK");
}

TEST_CASE("evaluate: emissions list and values agree; no hooks gives the same result",
          "[evaluate]") {
  const auto f = align<U>();
  const auto hooked = record<U>(f.s).result;
  const auto plain = evaluate<U>(f.s);
  CHECK(plain.value == hooked.value);
  REQUIRE(plain.emissions.size() == 3);
  CHECK(plain.emissions[0] == std::pair<StateId, U>{2, 0});
  CHECK(plain.emissions[1] == std::pair<StateId, U>{1, 8});
  CHECK(plain.emissions[2] == std::pair<StateId, U>{0, 11});
}

TEST_CASE("evaluate: empty hierarchy and hierarchy without offers", "[evaluate]") {
  SeededOrderedHierarchy<U> empty;
  const auto r0 = evaluate<U>(empty);
  CHECK(r0.value.empty());
  CHECK(r0.emissions.empty());

  auto f = chain<U>();
  f.s.offers.clear();
  const auto r1 = record<U>(f.s);
  CHECK(r1.result.value == std::vector<U>(3, WeightTraits<U>::infinity()));
  require_trace<U>(r1.events, {ev_finish<U>()});
}

TEST_CASE("evaluate: invalid initial offers are rejected", "[evaluate]") {
  auto f = chain<U>();
  f.s.offers = {{3, 0}};
  CHECK_THROWS_AS(evaluate<U>(f.s), std::out_of_range);
  f.s.offers = {{0, WeightTraits<U>::infinity()}};
  CHECK_THROWS_AS(evaluate<U>(f.s), std::invalid_argument);

  auto g = chain<double>();
  g.s.offers = {{0, std::numeric_limits<double>::quiet_NaN()}};
  CHECK_THROWS_AS(evaluate<double>(g.s), std::invalid_argument);
  g.s.offers = {{0, -1.0}};
  CHECK_THROWS_AS(evaluate<double>(g.s), std::invalid_argument);
  g.s.offers = {{0, std::numeric_limits<double>::infinity()}};
  CHECK_THROWS_AS(evaluate<double>(g.s), std::invalid_argument);
}
