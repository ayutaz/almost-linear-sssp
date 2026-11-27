// Self-tests of the Evaluate invariant checker (tests/support/paper/eval_invariants.hpp).
//
// Two hand-made hook traces follow Alg 2 / Alg 3 step by step on tiny threshold-SCC
// hierarchies (the normalisation values nu, S are worked out by hand from Sec 4.1):
//   chain:  a <-1-> b <-4-> c, pi = identity on the path 0 -1- 1 -4- 2, source (0, 0).
//           P = {a, b} (h 1), R = {P, c} (h 4). S(P) = 1 + 1 = 2, S(R) = 4 + 4 + 2 = 10;
//           4 S(P) >= h(R), so P and R share a unit table from h_min = 1: nu(P) = 1,
//           nu(R) = 4. delta = [0, 1, 5], every state protected.
//   align:  a <-3-> b <-8-> c, source (2, 0) at c. P = {a, b} (h 3), R = {P, c} (h 8).
//           S(P) = 6, S(R) = 22; table from 3: nu(P) = 3, nu(R) = 6. delta = [11, 8, 0].
//           P opens at D0 = 8 in R's slot [6, 12): D0 + S(P) = 14 >= 12, so the second
//           origin rule applies: the largest 12 - 3j <= 8 is 6 (Sec 4.3 l.794).
// Child order (Sec 4.2, derived in the test "the hand-trace fixtures have the Sec 4.2 child
// order" and checked against order_children()): no light arcs, so every child is a DFS root
// that finishes at once, in the F preorder a, b, c of its provisional representative; the rank
// order is the REVERSE finishing order: [c, P] at R with rep(P) = a, g(c, P) = D_F(a, c)
// (5 / 11), and [b, a] at P with g(b, a) = D_F(a, b) (1 / 3). The leaf layout is c, b, a.
// Supplemental offers (Sec 4.3 l.826, u = max{b, z + S(i) + g(i, j)} at x_j):
//   chain  R slot [0, 4): minimum 0 at a, i = P, left subrun {c}: c = max{4, 0 + 2 + 5} = 7;
//          P slot [0, 1): minimum 0 at a, i = a, left subrun {b}: b = max{1, 0 + 0 + 1} = 1.
//   align  R slot [0, 6): minimum 0 at c, i = c, right subrun {P}: rep(P) = a =
//          max{6, 0 + 0 + 11} = 11; P slot [6, 9): minimum 8 at b, i = b, right subrun {a}:
//          a = max{9, 8 + 0 + 3} = 11.
// Every broken variant changes one thing and must raise the named violation; the correct
// traces and plain state-graph Dijkstra traces must raise none. The checker itself does not see
// the child order, the representatives or the gaps (it checks I1-I8 and the protocol only), so
// the traces are tied to order_children() and to evaluate() by the tests that recompute the
// supplemental offers and replay Evaluate.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"

using sssp::CsrGraph;
using sssp::paper::NodeId;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateGraph;
using sssp::paper::StateId;
using namespace sssp::test::paper;

namespace {

using U = std::uint64_t;
using K = EvalViolationKind;
using Ev = EvalEvent<U>;

struct Fixture {
  SeededOrderedHierarchy<U> s;
  CsrGraph<U> phys;
  EvalReference<U> ref;
  std::uint32_t p = 0;  // node ids (messages only)
  std::uint32_t r = 0;
};

Fixture make_fixture(U ab, U bc, sssp::Vertex source_vertex, StateId offer_state) {
  StateGraph<U> d;
  d.num_states = 3;
  d.arcs = {{0, 1, ab}, {1, 0, ab}, {1, 2, bc}, {2, 1, bc}};
  Fixture f;
  f.s.hierarchy = sssp::paper::build_threshold_scc_hierarchy<U>(d, {0, 1, 2});
  f.phys = sssp::paper::physical_graph_of_arcs<U>(f.s.hierarchy.graph, {0, 1, 2}, 3);
  f.s.sources = {{source_vertex, 0}};
  f.s.offers = {{offer_state, 0}};
  f.ref = compute_eval_reference<U>(f.s, f.phys);
  const NodeId leaf_a = f.s.hierarchy.leaf_of_state[0];
  f.p = f.s.hierarchy.nodes[leaf_a].parent;
  f.r = f.s.hierarchy.nodes[f.p].parent;
  return f;
}

// Leaf-layout intervals (rank order: R = [c, P], P = [b, a]).
const std::vector<StateId> kP{1, 0};
const std::vector<StateId> kR{2, 1, 0};

// Correct trace of Evaluate on the chain fixture (see the file comment). Arc ids: 0 a->b,
// 1 b->a, 2 b->c, 3 c->b.
std::vector<Ev> chain_trace(const Fixture& f) {
  const auto p = f.p;
  const auto r = f.r;
  return {
      ev_initial<U>(0, 0),
      ev_open<U>(r, kR, 0, 0, 4, 10),           // root: origin = D0, c = D0 + S = 10
      ev_scan_begin<U>(r, kR, 0, 4, 0, 4, 10),  // slot [0, 4)
      ev_supplemental<U>(2, 7),                 // split at P, left subrun {c}: max{4, 0+2+5}
      ev_open<U>(p, kP, 0, 0, 1, 2),            // D0 + S(P) = 2 < 4: origin = D0
      ev_scan_begin<U>(p, kP, 0, 1, 0, 1, 2),   // c_P = min(10, 2)
      ev_supplemental<U>(1, 1),                 // split at a, left subrun {b}: max{1, 0+0+1}
      ev_emit<U>(0, 0),                         // a ({b} is passed first: key 1 not due)
      ev_relax<U>(0, 1),                        // a->b
      ev_scan_end<U>(p, kP),                    // b's key 1 is outside [0, 1)
      ev_scan_begin<U>(p, kP, 1, 2, 1, 2, 2),   // next slot of P
      ev_emit<U>(1, 1),                         // b
      ev_relax<U>(1, 2),                        // b->a (settled; still reported)
      ev_relax<U>(2, 5),                        // b->c
      ev_scan_end<U>(p, kP),
      ev_retire<U>(p, kP),    // P completed
      ev_scan_end<U>(r, kR),  // c's key 5 is outside [0, 4)
      ev_scan_begin<U>(r, kR, 4, 8, 4, 8, 10),
      ev_emit<U>(2, 5),   // c
      ev_relax<U>(3, 9),  // c->b
      ev_scan_end<U>(r, kR),
      ev_finish<U>(),
  };
}

// Correct trace on the alignment fixture (second origin rule).
std::vector<Ev> align_trace(const Fixture& f) {
  const auto p = f.p;
  const auto r = f.r;
  return {
      ev_initial<U>(2, 0),
      ev_open<U>(r, kR, 0, 0, 6, 22),
      ev_scan_begin<U>(r, kR, 0, 6, 0, 6, 22),
      ev_supplemental<U>(0, 11),  // split at c, right subrun {P}, at rep(P) = a: max{6, 0+0+11}
      ev_emit<U>(2, 0),           // c (then P's run is not due: key 8 at b)
      ev_relax<U>(3, 8),          // c->b
      ev_scan_end<U>(r, kR),
      ev_scan_begin<U>(r, kR, 6, 12, 6, 12, 22),
      ev_open<U>(p, kP, 8, 6, 3, 6),            // second rule: origin 6
      ev_scan_begin<U>(p, kP, 6, 9, 6, 9, 14),  // c_P = min(22, 8 + 6)
      ev_supplemental<U>(0, 11),                // split at b, right subrun {a}: max{9, 8+0+3}
      ev_emit<U>(1, 8),                         // b
      ev_relax<U>(1, 11),                       // b->a
      ev_relax<U>(2, 16),                       // b->c (settled)
      ev_scan_end<U>(p, kP),
      ev_scan_begin<U>(p, kP, 9, 12, 9, 12, 14),
      ev_emit<U>(0, 11),   // a
      ev_relax<U>(0, 14),  // a->b (settled)
      ev_scan_end<U>(p, kP),
      ev_retire<U>(p, kP),
      ev_scan_end<U>(r, kR),
      ev_finish<U>(),
  };
}

InvariantReport run(const Fixture& f, const std::vector<Ev>& trace, CheckerOptions opt = {}) {
  InvariantChecker<U> c(f.s, f.ref, opt);
  replay<U>(trace, c);
  return c.report();
}

// Index of the first event matching kind and state (or node for node events).
std::size_t find(const std::vector<Ev>& t, EvalEventKind k, StateId x, std::size_t from = 0) {
  for (std::size_t i = from; i < t.size(); ++i) {
    if (t[i].kind == k && (x == sssp::paper::kNoState || t[i].state == x)) return i;
  }
  FAIL("event not found");
  return 0;
}
constexpr StateId kAny = sssp::paper::kNoState;

Fixture chain() {
  return make_fixture(1, 4, 0, 0);
}
Fixture align() {
  return make_fixture(3, 8, 2, 2);
}

}  // namespace

TEST_CASE("eval invariants: the fixtures have the hand-computed oracle values",
          "[eval-invariants]") {
  const auto c = chain();
  CHECK(c.ref.state_target == std::vector<U>{0, 1, 5});
  CHECK(c.ref.is_protected == std::vector<bool>{true, true, true});
  const auto a = align();
  CHECK(a.ref.state_target == std::vector<U>{11, 8, 0});
  CHECK(a.ref.is_protected == std::vector<bool>{true, true, true});
}

namespace {
// The Sec 4.1 / 4.2 data the hand traces rely on, read from normalize() + order_children().
struct OrderedFixture {
  sssp::paper::NormalizedHierarchy<U> nh;
  NodeId r = 0;  // normalised node ids
  NodeId p = 0;
  NodeId leaf[3] = {0, 0, 0};
};

OrderedFixture ordered(const Fixture& f) {
  OrderedFixture o{sssp::paper::normalize(f.s.hierarchy)};
  sssp::paper::order_children(o.nh);
  o.r = o.nh.node_of_input[f.r];
  o.p = o.nh.node_of_input[f.p];
  for (StateId x = 0; x < 3; ++x) o.leaf[x] = o.nh.forest.leaf_of_state[x];
  return o;
}
}  // namespace

// Hand derivation from the paper (both fixtures have the same shape). Sec 4.1: F has the links
// a - b (P) and b - c (R), rooted at a, preorder a, b, c. No arc is shorter than the unit at its
// LCA (1 < 1 / 4 < 4 in chain, 3 < 3 / 8 < 6 in align), so neither node has light arcs.
// Sec 4.2 l.726-730: at R the provisional representatives are a (for P) and c; the DFS launches
// at P (preorder of a is smaller), P has no outgoing light arc and finishes at once, then c.
// Lemma 4.4 (l.735): the rank order is the REVERSE finishing order, [c, P]; rep(P) = a, the
// provisional representative of a DFS root; the charge between the two marks is the move along
// F between the successive DFS roots, D_F(a, c). At P the same gives [b, a] and D_F(a, b).
TEST_CASE("eval invariants: the hand-trace fixtures have the Sec 4.2 child order",
          "[eval-invariants]") {
  struct Expected {
    Fixture f;
    U nu_p, nu_r, mass_p, mass_r, gap_cp, gap_ba;
  };
  for (const auto& e :
       {Expected{chain(), 1, 4, 2, 10, 1 + 4, 1}, Expected{align(), 3, 6, 6, 22, 3 + 8, 3}}) {
    const auto o = ordered(e.f);
    const auto& nh = o.nh;
    REQUIRE(nh.children_ordered);
    CHECK(nh.forest.nodes[o.r].parent == sssp::paper::kNoNode);
    CHECK(nh.forest.nodes[o.p].parent == o.r);
    CHECK(nh.light(o.r).empty());
    CHECK(nh.light(o.p).empty());
    CHECK(nh.unit[o.p] == e.nu_p);
    CHECK(nh.unit[o.r] == e.nu_r);
    CHECK(nh.mass[o.p] == e.mass_p);
    CHECK(nh.mass[o.r] == e.mass_r);
    // R: ranks [c, P], rep(c) = c, rep(P) = a, g(c, P) = D_F(a, c).
    CHECK(nh.forest.nodes[o.r].children == std::vector<NodeId>{o.leaf[2], o.p});
    CHECK(nh.rep[o.leaf[2]] == StateId{2});
    CHECK(nh.rep[o.p] == StateId{0});
    CHECK(nh.gap_next[o.leaf[2]] == e.gap_cp);
    CHECK(nh.gap_next[o.p] == U{0});
    // P: ranks [b, a], g(b, a) = D_F(a, b).
    CHECK(nh.forest.nodes[o.p].children == std::vector<NodeId>{o.leaf[1], o.leaf[0]});
    CHECK(nh.rep[o.leaf[1]] == StateId{1});
    CHECK(nh.rep[o.leaf[0]] == StateId{0});
    CHECK(nh.gap_next[o.leaf[1]] == e.gap_ba);
    CHECK(nh.gap_next[o.leaf[0]] == U{0});
  }
}

namespace {
// One split of Alg 2 l.5 that creates subruns: the run is the whole child list of the node
// (the first split there), z is the run minimum at `min_state`, b the nominal slot end.
struct Split {
  bool at_root;  // R or P
  StateId min_state;
  U z;
  U b;
};

// The supplemental offers of Alg 2 l.6 / Sec 4.3 l.826 for these splits, with the child order,
// representatives, gaps and masses taken from order_children(), not from the hand trace:
// u = max{b, z + S(i) + g(i, j)} at x_j for the child j of each new subrun adjacent to i.
std::vector<Ev> supplemental_from_order(const OrderedFixture& o, const std::vector<Split>& splits) {
  const auto& nh = o.nh;
  std::vector<Ev> out;
  for (const auto& s : splits) {
    const NodeId x = s.at_root ? o.r : o.p;
    NodeId i = nh.forest.leaf_of_state[s.min_state];
    while (nh.forest.nodes[i].parent != x) i = nh.forest.nodes[i].parent;
    const auto& ch = nh.forest.nodes[x].children;
    std::size_t k = 0;
    while (ch[k] != i) ++k;
    const U base = s.z + nh.mass[i];
    if (k > 0) {  // left subrun, adjacent child ch[k - 1]
      const NodeId j = ch[k - 1];
      out.push_back(ev_supplemental<U>(nh.rep[j], std::max(s.b, base + nh.gap_next[j])));
    }
    if (k + 1 < ch.size()) {  // right subrun, adjacent child ch[k + 1]
      const NodeId j = ch[k + 1];
      out.push_back(ev_supplemental<U>(nh.rep[j], std::max(s.b, base + nh.gap_next[i])));
    }
  }
  return out;
}

// Every field but `node` (a free id, messages only).
bool same_event(const Ev& a, const Ev& b) {
  return a.kind == b.kind && a.state == b.state && a.arc == b.arc && a.states == b.states &&
         a.value == b.value && a.origin == b.origin && a.unit == b.unit && a.mass == b.mass &&
         a.lo == b.lo && a.hi == b.hi && a.slot_lo == b.slot_lo && a.slot_hi == b.slot_hi &&
         a.deadline == b.deadline;
}

// Hand-derived split points of the two traces (see the file comment).
const std::vector<Split> kChainSplits{{true, 0, 0, 4}, {false, 0, 0, 1}};
const std::vector<Split> kAlignSplits{{true, 2, 0, 6}, {false, 1, 8, 9}};
}  // namespace

TEST_CASE("eval invariants: hand-trace supplemental offers follow order_children",
          "[eval-invariants]") {
  struct Case {
    const char* name;
    Fixture f;
    std::vector<Ev> trace;
    std::vector<Split> splits;
  };
  const auto c = chain();
  const auto a = align();
  for (const auto& e : {Case{"chain", c, chain_trace(c), kChainSplits},
                        Case{"align", a, align_trace(a), kAlignSplits}}) {
    INFO(e.name);
    const auto want = supplemental_from_order(ordered(e.f), e.splits);
    std::vector<Ev> got;
    std::vector<U> slot_end;  // nominal end b of the scan each offer is posted in
    U b = 0;
    for (const auto& ev : e.trace) {
      if (ev.kind == EvalEventKind::kScanBegin) b = ev.slot_hi;
      if (ev.kind == EvalEventKind::kSupplementalOffer) {
        got.push_back(ev);
        slot_end.push_back(b);
      }
    }
    REQUIRE(got.size() == want.size());
    REQUIRE(got.size() == e.splits.size());  // one subrun per split in both fixtures
    for (std::size_t k = 0; k < got.size(); ++k) {
      INFO("supplemental offer " << k);
      CHECK(slot_end[k] == e.splits[k].b);
      CHECK(got[k].state == want[k].state);
      CHECK(got[k].value == want[k].value);
    }
  }
}

TEST_CASE("eval invariants: Evaluate reproduces the hand traces", "[eval-invariants]") {
  const auto c = chain();
  const auto a = align();
  for (const auto& [name, f, want] :
       {std::tuple{"chain", c, chain_trace(c)}, std::tuple{"align", a, align_trace(a)}}) {
    INFO(name);
    RecordingObserver<U> rec;
    ForwardingObserver<sssp::paper::EvalHooks<U>, U> adapter(rec);
    sssp::paper::EvalOptions<U> opt;
    opt.hooks = &adapter;
    (void)sssp::paper::evaluate<U>(f.s, opt);
    const auto& got = rec.events;
    for (std::size_t k = 0; k < std::min(got.size(), want.size()); ++k) {
      INFO("event " << k << ": got " << static_cast<int>(got[k].kind) << " state " << got[k].state
                    << " value " << got[k].value << ", want " << static_cast<int>(want[k].kind)
                    << " state " << want[k].state << " value " << want[k].value);
      REQUIRE(same_event(got[k], want[k]));
    }
    REQUIRE(got.size() == want.size());
  }
}

TEST_CASE("eval invariants: correct hand traces raise no violation", "[eval-invariants]") {
  const auto c = chain();
  const auto rc = run(c, chain_trace(c));
  INFO(rc.summary());
  CHECK(rc.ok());
  const auto a = align();
  const auto ra = run(a, align_trace(a));
  INFO(ra.summary());
  CHECK(ra.ok());
}

TEST_CASE("eval invariants: I1 detects offers below delta", "[eval-invariants]") {
  const auto f = chain();
  auto t = chain_trace(f);
  t[find(t, EvalEventKind::kSupplementalOffer, 2)].value = 4;  // delta(c) = 5
  const auto rep = run(f, t);
  CHECK(rep.has(K::kOfferBelowTarget));
  CHECK(rep.has(Invariant::kI1));
}

TEST_CASE("eval invariants: I2 detects double, wrong, low and missing emissions",
          "[eval-invariants]") {
  const auto f = chain();
  SECTION("double emission") {
    auto t = chain_trace(f);
    const auto i = find(t, EvalEventKind::kEmit, 1);
    t.insert(t.begin() + static_cast<std::ptrdiff_t>(i) + 1, ev_emit<U>(1, 1));
    const auto rep = run(f, t);
    CHECK(rep.has(K::kEmitTwice));
    CHECK(rep.has(Invariant::kI2));
  }
  SECTION("protected state emits a value above delta") {
    // Drop the relaxation b->c: c keeps the supplemental key 7 and emits it.
    auto t = chain_trace(f);
    t.erase(t.begin() +
            static_cast<std::ptrdiff_t>(
                find(t, EvalEventKind::kRelaxOffer, kAny, find(t, EvalEventKind::kEmit, 1)) + 1));
    t[find(t, EvalEventKind::kEmit, 2)].value = 7;
    const auto rep = run(f, t);
    INFO(rep.summary());
    CHECK(rep.has(K::kProtectedWrongValue));
    CHECK_FALSE(rep.has(K::kEmitNotKey));
  }
  SECTION("emission below delta") {
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kEmit, 2)].value = 3;
    const auto rep = run(f, t);
    CHECK(rep.has(K::kEmitBelowTarget));
  }
  SECTION("a protected state never emits") {
    auto t = chain_trace(f);
    const auto i = find(t, EvalEventKind::kEmit, 2);
    t.erase(t.begin() + static_cast<std::ptrdiff_t>(i),
            t.begin() + static_cast<std::ptrdiff_t>(i) + 2);  // emit c, relax c->b
    const auto rep = run(f, t);
    CHECK(rep.has(K::kProtectedNotEmitted));
    CHECK(rep.has(K::kKeyLeftInPortion));  // c's key 5 stays in [4, 8)
  }
}

TEST_CASE("eval invariants: I3 detects offers below the portion and unpadded supplements",
          "[eval-invariants]") {
  const auto f = chain();
  auto t = chain_trace(f);
  // During R's slot [4, 8): a supplemental offer 2 at b (delta 1, so I1 holds).
  const auto i = find(t, EvalEventKind::kScanBegin, kAny, find(t, EvalEventKind::kRetire, kAny));
  t.insert(t.begin() + static_cast<std::ptrdiff_t>(i) + 1, ev_supplemental<U>(1, 2));
  const auto rep = run(f, t);
  CHECK(rep.has(K::kOfferBelowPortion));
  CHECK(rep.has(K::kSupplementalNotPadded));
  CHECK_FALSE(rep.has(K::kOfferBelowTarget));
  CHECK(rep.has(Invariant::kI3));
}

TEST_CASE("eval invariants: I4 detects emissions outside the portion or node",
          "[eval-invariants]") {
  const auto f = chain();
  SECTION("c emitted during R's slot [0, 4)") {
    auto t = chain_trace(f);
    const auto i = find(t, EvalEventKind::kEmit, 2);
    auto em = t[i];
    t.erase(t.begin() + static_cast<std::ptrdiff_t>(i));
    const auto end_r = find(t, EvalEventKind::kRetire, kAny) + 1;  // still inside [0, 4)
    t.insert(t.begin() + static_cast<std::ptrdiff_t>(end_r), em);
    const auto rep = run(f, t);
    INFO(rep.summary());
    CHECK(rep.has(K::kEmitOutsidePortion));
  }
  SECTION("c emitted inside P's scan") {
    auto t = chain_trace(f);
    const auto i = find(t, EvalEventKind::kEmit, 2);
    t.erase(t.begin() + static_cast<std::ptrdiff_t>(i));
    const auto in_p =
        find(t, EvalEventKind::kRelaxOffer, kAny, find(t, EvalEventKind::kEmit, 1)) + 2;
    t.insert(t.begin() + static_cast<std::ptrdiff_t>(in_p), ev_emit<U>(2, 5));
    const auto rep = run(f, t);
    INFO(rep.summary());
    CHECK(rep.has(K::kEmitOutsideNode));
    CHECK(rep.has(K::kEmitOutsidePortion));  // 5 is not in [1, 2)
  }
  SECTION("a child portion outside the parent's") {
    auto t = chain_trace(f);
    auto& sb = t[find(t, EvalEventKind::kScanBegin, kAny, find(t, EvalEventKind::kOpen, kAny, 2))];
    sb.deadline = 11;  // parent's deadline is 10
    const auto rep = run(f, t);
    CHECK(rep.has(K::kPortionNotNested));
  }
}

TEST_CASE("eval invariants: I5 detects a due key left in the scanned portion",
          "[eval-invariants]") {
  const auto f = chain();
  auto t = chain_trace(f);
  // Skip b's emission in P's slot [1, 2) (and its relaxations); the scan still returns.
  const auto i = find(t, EvalEventKind::kEmit, 1);
  t.erase(t.begin() + static_cast<std::ptrdiff_t>(i),
          t.begin() + static_cast<std::ptrdiff_t>(i) + 3);
  const auto rep = run(f, t);
  INFO(rep.summary());
  CHECK(rep.has(K::kKeyLeftInPortion));
  CHECK(rep.has(Invariant::kI5));
}

TEST_CASE("eval invariants: I5 accepts keys beyond the deadline or of retired states",
          "[eval-invariants]") {
  const auto f = chain();
  auto t = chain_trace(f);
  // Retire P (with b unemitted) before P's second slot: b is then not "left" in I5's sense,
  // but I6 and I2 fire.
  const auto i = find(t, EvalEventKind::kScanBegin, kAny, find(t, EvalEventKind::kScanEnd, kAny));
  t.erase(t.begin() + static_cast<std::ptrdiff_t>(i),
          t.begin() + static_cast<std::ptrdiff_t>(find(t, EvalEventKind::kRetire, kAny)));
  const auto rep = run(f, t);
  INFO(rep.summary());
  CHECK_FALSE(rep.has(K::kKeyLeftInPortion));
  CHECK(rep.has(K::kRetiredProtectedUnemitted));
  CHECK(rep.has(K::kProtectedNotEmitted));
}

TEST_CASE("eval invariants: I6 detects early retirement, a too small mass and late deadlines",
          "[eval-invariants]") {
  const auto f = chain();
  SECTION("mass below the physical diameter") {
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kOpen, kAny)].mass = 1;  // delta(c) = 5 > 0 + 1
    const auto rep = run(f, t);
    CHECK(rep.has(K::kMassBelowDiameter));
    CHECK(rep.has(Invariant::kI6));
  }
  SECTION("deadline after the node's expiry") {
    auto t = chain_trace(f);
    const auto i = find(t, EvalEventKind::kOpen, kAny, 2);
    t[i + 1].deadline = 3;  // P's expiry is 0 + 2
    const auto rep = run(f, t);
    CHECK(rep.has(K::kDeadlineAfterExpiry));
  }
}

TEST_CASE("eval invariants: I7 detects origin and grid violations", "[eval-invariants]") {
  SECTION("first rule: D0 + S < b forces origin = D0") {
    const auto f = chain();
    auto t = chain_trace(f);
    auto& open_p = t[find(t, EvalEventKind::kOpen, kAny, 2)];
    open_p.origin = 1;  // also no longer <= D0
    CHECK(run(f, t).has(K::kOriginRule));
  }
  SECTION("second rule: origin must be the largest b - j nu <= D0") {
    const auto f = align();
    auto t = align_trace(f);
    const auto i = find(t, EvalEventKind::kOpen, kAny, 2);
    t[i].origin = 8;  // the first rule misapplied
    // Keep the slots consistent with the wrong origin so only the rule itself is wrong.
    t[i + 1].lo = 8;
    t[i + 1].hi = 11;
    t[i + 1].slot_lo = 8;
    t[i + 1].slot_hi = 11;
    const auto rep = run(f, t);
    INFO(rep.summary());
    CHECK(rep.has(K::kOriginRule));
  }
  SECTION("root origin differs from D0") {
    const auto f = chain();
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kOpen, kAny)].origin = 1;
    CHECK(run(f, t).has(K::kOriginRule));
  }
  SECTION("a slot not on the node's grid") {
    const auto f = align();
    auto t = align_trace(f);
    auto& sb =
        t[find(t, EvalEventKind::kScanBegin, kAny,
               find(t, EvalEventKind::kScanEnd, kAny, find(t, EvalEventKind::kOpen, kAny, 2)))];
    sb.slot_lo = 10;  // P's grid is 6 + 3k
    sb.slot_hi = 13;
    sb.lo = 10;
    const auto rep = run(f, t);
    CHECK(rep.has(K::kSlotNotOnGrid));
    CHECK(rep.has(Invariant::kI7));
  }
}

TEST_CASE("eval invariants: I8 detects opening above an unemitted protected delta",
          "[eval-invariants]") {
  // Open P at D0 = 1 although a (protected, delta 0) is pending below it.
  const auto f = chain();
  auto t = chain_trace(f);
  t[find(t, EvalEventKind::kOpen, kAny, 2)].value = 1;
  const auto rep = run(f, t);
  INFO(rep.summary());
  CHECK(rep.has(K::kOpenAboveProtected));
  CHECK(rep.has(K::kOpenNotMinimum));
  CHECK(rep.has(Invariant::kI8));
}

TEST_CASE("eval invariants: protocol checks", "[eval-invariants]") {
  const auto f = chain();
  SECTION("emission differs from the current key") {
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kEmit, 1)].value = 2;  // b's key is 1
    CHECK(run(f, t).has(K::kEmitNotKey));
  }
  SECTION("relaxation value differs from emitted(tail) + length") {
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kRelaxOffer, kAny, find(t, EvalEventKind::kEmit, 1)) + 1].value = 6;
    CHECK(run(f, t).has(K::kRelaxValue));
  }
  SECTION("unknown node") {
    auto t = chain_trace(f);
    t[find(t, EvalEventKind::kOpen, kAny, 2)].states = {0, 2};
    CHECK(run(f, t).has(K::kUnknownNode));
  }
  SECTION("unbalanced scans") {
    auto t = chain_trace(f);
    t.erase(t.end() - 2);  // last scan_end
    CHECK(run(f, t).has(K::kScanEndMismatch));
  }
  SECTION("supplemental offers while disabled") {
    CheckerOptions opt;
    opt.supplemental_offers = false;
    CHECK(run(f, chain_trace(f), opt).has(K::kSupplementalDisabled));
  }
  SECTION("emission outside any scan at a non-root leaf") {
    auto t = chain_trace(f);
    t.insert(t.begin() + 1, ev_emit<U>(0, 0));
    CHECK(run(f, t).has(K::kEmitOutsideScan));
  }
  SECTION("hooks after on_finish") {
    auto t = chain_trace(f);
    t.push_back(ev_emit<U>(0, 0));
    CHECK(run(f, t).has(K::kEventAfterFinish));
  }
}

TEST_CASE("eval invariants: a state-graph Dijkstra trace satisfies I1 and I2",
          "[eval-invariants]") {
  for (const auto& f : {chain(), align()}) {
    const auto rep = run(f, state_dijkstra_trace<U>(f.s, f.ref), CheckerOptions{true, true, true});
    INFO(rep.summary());
    CHECK(rep.ok());
  }
  // Without allow_unscoped_emits the same trace is a protocol violation, not an invariant one.
  const auto f = chain();
  const auto rep = run(f, state_dijkstra_trace<U>(f.s, f.ref));
  CHECK(rep.has(K::kEmitOutsideScan));
  CHECK_FALSE(rep.has(Invariant::kI1));
  CHECK_FALSE(rep.has(Invariant::kI2));
}

namespace {
// Stand-in for a src-side virtual hook type with the same member functions.
struct MockSrcHooks {
  MockSrcHooks() = default;
  MockSrcHooks(const MockSrcHooks&) = default;
  MockSrcHooks(MockSrcHooks&&) = default;
  MockSrcHooks& operator=(const MockSrcHooks&) = default;
  MockSrcHooks& operator=(MockSrcHooks&&) = default;
  virtual ~MockSrcHooks() = default;
  virtual void on_initial_offer(StateId, U) {}
  virtual void on_relax_offer(sssp::paper::ArcId, U) {}
  virtual void on_supplemental_offer(StateId, U) {}
  virtual void on_open(std::uint32_t, std::span<const StateId>, U, U, U, U) {}
  virtual void on_scan_begin(std::uint32_t, std::span<const StateId>, U, U, U, U, U) {}
  virtual void on_scan_end(std::uint32_t, std::span<const StateId>) {}
  virtual void on_emit(StateId, U) {}
  virtual void on_retire(std::uint32_t, std::span<const StateId>) {}
  virtual void on_finish() {}
};
}  // namespace

TEST_CASE("eval invariants: a src-side hook type forwards to the checker", "[eval-invariants]") {
  const auto f = chain();
  InvariantChecker<U> checker(f.s, f.ref);
  ForwardingObserver<MockSrcHooks, U> adapter(checker);
  MockSrcHooks& hooks = adapter;  // what an Evaluate in src/ would hold
  RecordingObserver<U> rec;
  replay<U>(chain_trace(f), rec);
  // Drive the src-side interface by hand from the recorded events.
  for (const auto& e : rec.events) {
    switch (e.kind) {
      case EvalEventKind::kInitialOffer:
        hooks.on_initial_offer(e.state, e.value);
        break;
      case EvalEventKind::kRelaxOffer:
        hooks.on_relax_offer(e.arc, e.value);
        break;
      case EvalEventKind::kSupplementalOffer:
        hooks.on_supplemental_offer(e.state, e.value);
        break;
      case EvalEventKind::kOpen:
        hooks.on_open(e.node, e.states, e.value, e.origin, e.unit, e.mass);
        break;
      case EvalEventKind::kScanBegin:
        hooks.on_scan_begin(e.node, e.states, e.lo, e.hi, e.slot_lo, e.slot_hi, e.deadline);
        break;
      case EvalEventKind::kScanEnd:
        hooks.on_scan_end(e.node, e.states);
        break;
      case EvalEventKind::kEmit:
        hooks.on_emit(e.state, e.value);
        break;
      case EvalEventKind::kRetire:
        hooks.on_retire(e.node, e.states);
        break;
      case EvalEventKind::kFinish:
        hooks.on_finish();
        break;
    }
  }
  INFO(checker.report().summary());
  CHECK(checker.report().ok());
  CHECK(checker.emissions().size() == 3);
  CHECK(rec.events.size() == chain_trace(f).size());
}
