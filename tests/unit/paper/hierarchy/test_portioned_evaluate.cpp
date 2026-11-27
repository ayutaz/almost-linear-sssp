// Unit tests of PortionedEvaluator: the cell-by-cell visit of Lemma 7.7 (Sec 7.4 of arXiv
// 2611.34911 v1, l.1835-1841 and l.1871-1909, Figure 2) on hand examples. Expected traces are
// derived by hand from the paper and from the hand derivations of test_evaluate.cpp (the same
// fixtures), not taken from an implementation run.
//
// align (see test_evaluate.cpp): states a = 0, b = 1, c = 2; arcs 0 a->b 3, 1 b->a 3, 2 b->c 8,
//   3 c->b 8; offer c = 0. Normalised: P = {a, b} with nu 3, S 6; R = {c, P} with nu 6, S 22;
//   ranks [c, P] and [b, a]; layout c, b, a. Bulk Evaluate: R slots [0, 6), [6, 12); P opens at
//   8 with origin 6 (second rule), slots [6, 9), [9, 12); emissions c 0, b 8, a 11.
//
// align cut at 10 (cells [0, 10), [10, +inf)):
//   cell 1: R opens at its first due key 0 (Lemma 7.7 "Open the root on its first due key"),
//     origin 0, cap c = min(+inf, 0 + 22) = 22. Root slot [0, 6) lies inside the cell: c emits 0
//     (supplemental 11 at a, relaxation c->b 8). Root slot [6, 12) meets the cell in [6, 10):
//     {P} splits at b (8), P opens at 8, origin 6, deadline 14; P slot [6, 9) is inside the
//     portion: b emits 8 (supplemental 11 at a, relaxations 11 and 16). P slot [9, 12) meets the
//     portion in [9, 10): a (11) is not due; the slot is clipped and stays P's current slot. P
//     is not retired (14 >= 10, "No temporary cut retires a node", l.1891). Root slot [6, 12)
//     is clipped too and stays current.
//   cell 2: root slot [6, 12) resumes with the portion [10, 12) (l.1840-1841, l.1883-1885); the
//     cursor restarts at R's first live piece, the registered child P; P's slot [9, 12)
//     resumes with the portion [10, 12): {a} splits at 11, a emits 11 (relaxation 14 to the
//     settled b). P and R complete.
//   One child slot and one root slot are clipped, and each is resumed once.
//
// align with deadline 9 (cells [0, 10), [10, +inf)): cap min(9, 22) = 9; P's deadline is
//   min(9, 14) = 9; in [9, 10) the key 11 of a is above the cap; P covers its deadline (9 < 10)
//   and retires with a unemitted; R then has no pieces left and completes after one cell.
//
// pair without supplemental offers (see test_evaluate.cpp): R = {a, b}, nu 1, S 2, offer a = 0,
//   b never gets a key. Cells of width 1: [0, 1) emits a; [1, 2) and [2, 3) scan R with nothing
//   due; after [2, 3) the cap 2 is covered (2 < 3, "including equality") and R retires.
//
// single: one state, no arcs, offer 5. A singleton root emits its due minimum once in its cell
//   (l.1835); with deadline 3 it is retired once a cell ends after 3.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/primitives/split_findmin.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"

using sssp::CsrGraph;
using sssp::WeightTraits;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::PortionedEvaluator;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using namespace sssp::test::paper;

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
Fixture<W> align() {
  return path_fixture<W>(W{3}, W{8}, 2, 2);
}

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

template <class W>
Fixture<W> single_fixture(W offer) {
  StateGraph<W> d;
  d.num_states = 1;
  auto h = sssp::paper::make_leaf_forest<W>(d, {0});
  sssp::paper::collect_roots(h);
  Fixture<W> f;
  f.s.hierarchy = std::move(h);
  f.phys = CsrGraph<W>::from_edges(1, {});
  f.s.sources = {{0, offer}};
  f.s.offers = {{0, offer}};
  return f;
}

/// Two states and no internal node: two roots (Lemma 7.7 requires one).
template <class W>
Fixture<W> two_roots() {
  StateGraph<W> d;
  d.num_states = 2;
  auto h = sssp::paper::make_leaf_forest<W>(d, {0, 1});
  sssp::paper::collect_roots(h);
  Fixture<W> f;
  f.s.hierarchy = std::move(h);
  f.phys = CsrGraph<W>::from_edges(2, {});
  return f;
}

template <class W>
NormalizedHierarchy<W> prepare(const SeededOrderedHierarchy<W>& s) {
  NormalizedHierarchy<W> nh = sssp::paper::normalize<W>(s.hierarchy);
  sssp::paper::order_children<W>(nh);
  return nh;
}

template <class W>
W inf() {
  return WeightTraits<W>::infinity();
}

using Emissions = std::vector<std::pair<StateId, U>>;

template <class W>
std::vector<std::pair<StateId, W>> to_vec(std::span<const std::pair<StateId, W>> s) {
  return {s.begin(), s.end()};
}

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
    INFO("event " << i);
    REQUIRE(same_event(got[i], want[i]));
  }
  REQUIRE(got.size() == want.size());
}

const std::vector<StateId> kR{2, 1, 0};
const std::vector<StateId> kP{1, 0};

/// Hand-derived trace of align cut at 10 (file comment).
template <class W>
std::vector<EvalEvent<W>> align_cut_trace() {
  return {
      ev_initial<W>(2, 0),
      // cell [0, 10)
      ev_open<W>(0, kR, 0, 0, 6, 22),
      ev_scan_begin<W>(0, kR, 0, 6, 0, 6, 22),
      ev_supplemental<W>(0, 11),
      ev_emit<W>(2, 0),
      ev_relax<W>(3, 8),
      ev_scan_end<W>(0, kR),
      ev_scan_begin<W>(0, kR, 6, 10, 6, 12, 22),  // root slot clipped at the cell end
      ev_open<W>(0, kP, 8, 6, 3, 6),
      ev_scan_begin<W>(0, kP, 6, 9, 6, 9, 14),
      ev_supplemental<W>(0, 11),
      ev_emit<W>(1, 8),
      ev_relax<W>(1, 11),
      ev_relax<W>(2, 16),
      ev_scan_end<W>(0, kP),
      ev_scan_begin<W>(0, kP, 9, 10, 9, 12, 14),  // child slot clipped at the cell end
      ev_scan_end<W>(0, kP),
      ev_scan_end<W>(0, kR),
      // cell [10, +inf)
      ev_scan_begin<W>(0, kR, 10, 12, 6, 12, 22),  // root slot resumed
      ev_scan_begin<W>(0, kP, 10, 12, 9, 12, 14),  // child slot resumed from P's first piece
      ev_emit<W>(0, 11),
      ev_relax<W>(0, 14),
      ev_scan_end<W>(0, kP),
      ev_retire<W>(0, kP),
      ev_scan_end<W>(0, kR),
      ev_finish<W>(),
  };
}

[[maybe_unused]] U counter(const sssp::Diagnostics& d, const std::string& name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

}  // namespace

TEMPLATE_TEST_CASE("portioned: a cell boundary clips a child slot, which resumes (Lemma 7.7)",
                   "[portioned][lemma7.7]", std::uint64_t, double) {
  using W = TestType;
  const auto f = align<W>();
  const auto nh = prepare<W>(f.s);
  RecordingObserver<W> rec;
  ForwardingObserver<EvalHooks<W>, W> adapter(rec);
  EvalOptions<W> opt;
  opt.hooks = &adapter;
  PortionedEvaluator<W> ev(nh, opt);
  REQUIRE(ev.num_states() == 3);
  ev.initial_offer(2, W{0});
  CHECK_FALSE(ev.opened());
  const auto c1 = to_vec<W>(ev.visit_cell(W{0}, W{10}));
  CHECK(c1 == std::vector<std::pair<StateId, W>>{{2, W{0}}, {1, W{8}}});
  CHECK(ev.opened());
  CHECK_FALSE(ev.finished());
  const auto c2 = to_vec<W>(ev.visit_cell(W{10}, inf<W>()));
  CHECK(c2 == std::vector<std::pair<StateId, W>>{{0, W{11}}});
  CHECK(ev.finished());
  const auto r = ev.take_result();
  CHECK(r.value == std::vector<W>{11, 8, 0});
  require_trace<W>(rec.events, align_cut_trace<W>());
  CHECK_COUNTER(counter(r.diag, "eval_cells") == 2);
  CHECK_COUNTER(counter(r.diag, "eval_resumed_slot_visits") == 1);
  CHECK_COUNTER(counter(r.diag, "eval_clipped_slot_visits") == 1);
  CHECK_COUNTER(counter(r.diag, "eval_root_resumed_slot_visits") == 1);
  CHECK_COUNTER(counter(r.diag, "eval_root_clipped_slot_visits") == 1);
  CHECK_COUNTER(counter(r.diag, "eval_root_slots") == 2);
  CHECK_COUNTER(counter(r.diag, "eval_scans") == 6);
}

TEST_CASE("portioned: the cut run satisfies I1-I8 under the Sec 4 checker", "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  const auto ref = compute_eval_reference<U>(f.s, f.phys);
  InvariantChecker<U> checker(f.s, ref);
  ForwardingObserver<EvalHooks<U>, U> adapter(checker);
  EvalOptions<U> opt;
  opt.hooks = &adapter;
  PortionedEvaluator<U> ev(nh, opt);
  ev.initial_offer(2, 0);
  for (U lo = 0; lo < 30; ++lo) (void)ev.visit_cell(lo, lo + 1);
  CHECK(ev.finished());
  (void)ev.take_result();
  CHECK(checker.report().summary() == "OK");
}

TEMPLATE_TEST_CASE("portioned: one cell [0, +inf) reproduces the bulk trace event by event",
                   "[portioned]", std::uint64_t, double) {
  using W = TestType;
  for (const auto& f : {align<W>(), path_fixture<W>(W{1}, W{4}, 0, 0), pair_fixture<W>()}) {
    const auto nh = prepare<W>(f.s);
    RecordingObserver<W> bulk_rec;
    ForwardingObserver<EvalHooks<W>, W> bulk_adapter(bulk_rec);
    EvalOptions<W> bopt;
    bopt.hooks = &bulk_adapter;
    const auto bulk = sssp::paper::evaluate_normalized<W>(nh, f.s.offers, bopt);

    RecordingObserver<W> rec;
    ForwardingObserver<EvalHooks<W>, W> adapter(rec);
    EvalOptions<W> opt;
    opt.hooks = &adapter;
    PortionedEvaluator<W> ev(nh, opt);
    for (const auto& o : f.s.offers) ev.initial_offer(o.state, o.value);
    (void)ev.visit_cell(W{0}, inf<W>());
    CHECK(ev.finished());
    const auto r = ev.take_result();
    CHECK(r.value == bulk.value);
    CHECK(r.emissions == bulk.emissions);
    require_trace<W>(rec.events, bulk_rec.events);
  }
}

TEST_CASE("portioned: an inclusive atom deadline caps the root (Sec 7.4 l.1837-1838)",
          "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  PortionedEvaluator<U> ev(nh);
  ev.initial_offer(2, 0);
  ev.set_deadline(9);
  const auto c1 = to_vec<U>(ev.visit_cell(0, 10));
  CHECK(c1 == Emissions{{2, 0}, {1, 8}});
  CHECK(ev.finished());
  CHECK(to_vec<U>(ev.visit_cell(10, 20)).empty());
  const auto r = ev.take_result();
  CHECK(r.value == std::vector<U>{inf<U>(), 8, 0});
  // A later, smaller deadline before opening keeps the smaller one.
  PortionedEvaluator<U> ev2(nh);
  ev2.initial_offer(2, 0);
  ev2.set_deadline(30);
  ev2.set_deadline(9);
  ev2.set_deadline(40);
  (void)ev2.visit_cell(0, 10);
  CHECK(ev2.finished());
  CHECK(ev2.value()[0] == inf<U>());
}

TEST_CASE("portioned: the cap is covered including equality (Sec 7.4 l.1841)", "[portioned]") {
  const auto f = pair_fixture<U>();
  const auto nh = prepare<U>(f.s);
  EvalOptions<U> opt;
  opt.supplemental_offers = false;
  PortionedEvaluator<U> ev(nh, opt);
  ev.initial_offer(0, 0);
  CHECK(to_vec<U>(ev.visit_cell(0, 1)) == Emissions{{0, 0}});
  CHECK(to_vec<U>(ev.visit_cell(1, 2)).empty());
  CHECK_FALSE(ev.finished());
  CHECK(to_vec<U>(ev.visit_cell(2, 3)).empty());  // expiry 2 < 3: covered
  CHECK(ev.finished());
  CHECK(ev.value() == std::vector<U>{0, inf<U>()});
  CHECK_FALSE(ev.external_offer(1, 7));  // retired: later proposals are ignored (Lemma 7.8)
}

TEST_CASE("portioned: an empty cell leaves an unopened root unopened (Lemma 7.7)", "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  PortionedEvaluator<U> ev(nh);
  ev.initial_offer(2, 5);
  CHECK(to_vec<U>(ev.visit_cell(0, 2)).empty());
  CHECK_FALSE(ev.opened());
  // Before the root opens, cells need not be consecutive (Alg 7 visits an atom only from its
  // registration on); c = 5 is due in [4, 7), c->b gives 13, not due in R's slot [5, 11).
  CHECK(to_vec<U>(ev.visit_cell(4, 7)) == Emissions{{2, 5}});
  CHECK(ev.opened());
  CHECK_FALSE(ev.finished());
}

TEST_CASE("portioned: a singleton root emits its due minimum once or retires at its deadline",
          "[portioned]") {
  const auto f = single_fixture<U>(5);
  const auto nh = prepare<U>(f.s);
  PortionedEvaluator<U> ev(nh);
  ev.initial_offer(0, 5);
  CHECK(to_vec<U>(ev.visit_cell(0, 4)).empty());
  CHECK(to_vec<U>(ev.visit_cell(4, 8)) == Emissions{{0, 5}});
  CHECK(ev.finished());

  PortionedEvaluator<U> capped(nh);
  capped.initial_offer(0, 5);
  capped.set_deadline(3);
  CHECK(to_vec<U>(capped.visit_cell(0, 4)).empty());
  CHECK(capped.finished());
  CHECK(capped.value()[0] == inf<U>());
}

TEST_CASE("portioned: external offers between cells (Lemma 7.7 arrival hypothesis)",
          "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  PortionedEvaluator<U> ev(nh);
  ev.initial_offer(2, 0);
  (void)ev.visit_cell(0, 10);  // c 0, b 8 emitted; a has key 11
  // An effective decrease that arrives before its cell and not below the current cell is
  // accepted: a = 10 lies in the next cell [10, 20).
  CHECK(ev.external_offer(0, 10));
  CHECK_FALSE(ev.external_offer(0, 12));  // no decrease
  CHECK_FALSE(ev.external_offer(1, 3));   // b is settled: ignored even below the cell
  // An effective decrease inside an already visited cell violates the hypothesis.
  CHECK_THROWS_AS(ev.external_offer(0, 9), std::invalid_argument);
  CHECK(to_vec<U>(ev.visit_cell(10, 20)) == Emissions{{0, 10}});
  CHECK(ev.finished());
  const auto r = ev.take_result();
  CHECK(r.value == std::vector<U>{10, 8, 0});
  CHECK_COUNTER(counter(r.diag, "eval_external_offers") == 4);
  CHECK_COUNTER(counter(r.diag, "eval_external_decreases") == 1);
}

TEST_CASE("portioned: preconditions are checked", "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  SECTION("one root (Lemma 7.7)") {
    const auto two = two_roots<U>();
    const auto nh2 = prepare<U>(two.s);
    CHECK_THROWS_AS(PortionedEvaluator<U>(nh2), std::invalid_argument);
  }
  SECTION("a normalised and ordered hierarchy") {
    auto raw = sssp::paper::normalize<U>(f.s.hierarchy);
    CHECK_THROWS_AS(PortionedEvaluator<U>(raw), std::invalid_argument);
  }
  SECTION("offers") {
    PortionedEvaluator<U> ev(nh);
    CHECK_THROWS_AS(ev.initial_offer(3, 0), std::out_of_range);
    CHECK_THROWS_AS(ev.initial_offer(0, inf<U>()), std::invalid_argument);
    CHECK_THROWS_AS(ev.external_offer(3, 0), std::out_of_range);
    ev.initial_offer(2, 0);
    (void)ev.visit_cell(0, 4);
    CHECK_THROWS_AS(ev.initial_offer(1, 9), std::logic_error);
  }
  SECTION("cells") {
    PortionedEvaluator<U> ev(nh);
    ev.initial_offer(2, 0);
    CHECK_THROWS_AS(ev.visit_cell(3, 3), std::invalid_argument);
    (void)ev.visit_cell(0, 4);
    CHECK_THROWS_AS(ev.visit_cell(3, 8), std::invalid_argument);  // overlaps
    CHECK_THROWS_AS(ev.visit_cell(5, 8), std::invalid_argument);  // gap after opening
    CHECK_NOTHROW(ev.visit_cell(4, 8));
  }
  SECTION("the first cell starts at or before every finite initial timestamp") {
    PortionedEvaluator<U> ev(nh);
    ev.initial_offer(2, 0);
    CHECK_THROWS_AS(ev.visit_cell(1, 4), std::logic_error);
  }
  SECTION("the deadline is fixed once the root opened") {
    PortionedEvaluator<U> ev(nh);
    ev.initial_offer(2, 0);
    (void)ev.visit_cell(0, 4);
    CHECK_THROWS_AS(ev.set_deadline(3), std::logic_error);
  }
  SECTION("no offers from inside a visit") {
    struct Reentrant final : EvalHooks<U> {
      PortionedEvaluator<U>* ev = nullptr;
      bool threw = false;
      void on_emit(StateId /*x*/, U v) override {
        try {
          (void)ev->external_offer(0, v + 100);
        } catch (const std::logic_error&) {
          threw = true;
        }
      }
    } hooks;
    EvalOptions<U> opt;
    opt.hooks = &hooks;
    PortionedEvaluator<U> ev(nh, opt);
    hooks.ev = &ev;
    ev.initial_offer(2, 0);
    (void)ev.visit_cell(0, 4);
    CHECK(hooks.threw);
  }
  SECTION("a shared structure must hold the whole layout") {
    sssp::paper::SplitFindmin<U> sf(4);
    CHECK_THROWS_AS(PortionedEvaluator<U>(nh, sf, 2), std::out_of_range);
  }
}

TEST_CASE("portioned: two atoms share one split-findmin layout (Alg 7 l.8)", "[portioned]") {
  const auto f = align<U>();
  const auto nh = prepare<U>(f.s);
  // Layout: atom A at [0, 3), atom B at [3, 6); one interval before any registration.
  sssp::paper::SplitFindmin<U> sf(6);
  PortionedEvaluator<U> a(nh, sf, 0);
  PortionedEvaluator<U> b(nh, sf, 3);
  for (StateId x = 0; x < 3; ++x) {
    CHECK(a.position(x) < 3);
    CHECK(b.position(x) == a.position(x) + 3);
  }
  a.initial_offer(2, 0);
  b.initial_offer(2, 4);
  // The family-level view: one query over the unsplit interval sees the least key of both.
  CHECK(sf.query(0).key == 0);
  CHECK(sf.query(0).index == a.position(2));
  std::vector<std::pair<StateId, U>> ea, eb;
  for (U lo = 0; lo < 40; lo += 5) {
    for (const auto& e : a.visit_cell(lo, lo + 5)) ea.push_back(e);
    for (const auto& e : b.visit_cell(lo, lo + 5)) eb.push_back(e);
  }
  CHECK(a.finished());
  CHECK(b.finished());
  CHECK(ea == Emissions{{2, 0}, {1, 8}, {0, 11}});
  CHECK(eb == Emissions{{2, 4}, {1, 12}, {0, 15}});
}
