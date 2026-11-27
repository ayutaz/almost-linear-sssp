#pragma once
// Observation hooks for Algorithm 2 (Scan) / Algorithm 3 (Evaluate) and a checker of the
// invariants I1-I8 of research notes 02 sec 2.8 (derived from Lemma 4.6 and Sec 4.3).
// Test-only: an Evaluate implementation in src/ must not include this file. It calls hooks of
// its own type (EvalOptions::checks, research notes 02 sec 8); tests bind that type to the
// checker below through ForwardingObserver or a ten-line adapter. See "For the Evaluate
// implementer" at the end of this comment.
//
// Invariants (oracle: delta and protected set from eval_oracle.hpp):
//   I1  every offer (initial, relaxation, supplemental) at x has value >= delta(pi(x));
//   I2  each state emits at most once, never below delta(pi(x)); every protected state emits,
//       and its value is exactly delta(pi(x)) (checked at on_finish for missing ones);
//   I3  while a portion [a', b') is being scanned, every new offer has value >= a';
//       a supplemental offer is padded to the NOMINAL right end b (A7 / Q-28);
//   I4  an emission inside Scan(X, I, c) is at a state below X with value in I and <= c;
//       a child's portion is inside its parent's portion, its deadline is <= the parent's;
//       a node opens at a value due in the parent's portion; the slots of a node never go
//       back in time;
//   I5  on return from Scan(X, I, c), no unemitted, unretired state below X has a key in
//       I intersected with (-inf, c];
//   I6  when X retires, every protected state below X has emitted; every state below X has
//       delta <= D0 + S(X) (S is a physical-diameter bound, Lemma 4.3); a scan's deadline is
//       <= the node's expiry D0 + S(X);
//   I7  origin rule (Sec 4.3, l.794): at a root, origin = D0; otherwise origin = D0 if
//       D0 + S < b, else the largest b - j nu <= D0, where [a, b) is the parent's NOMINAL slot
//       (A6 / Q-28); the origin lies in [a, b); D0 lies in the first slot; every nominal slot
//       is [origin + k nu, origin + (k+1) nu);
//   I8  when X opens at D0, every unemitted protected state below X has delta >= D0.
// Protocol checks (kind kProtocol) make the trace self-consistent: an emission equals the
// current key (the least offer reported so far), a relaxation offer equals emitted(tail) + l,
// a node is identified by exactly the leaf set of some node of the input hierarchy, scans nest,
// initial offers come first, and supplemental offers appear only if enabled.
//
// Not checked (needs data the hooks do not carry): the rank part of I3 (keys falling into
// children ranked below the cursor are >= b'), the cursor order of Alg 2 l.2 and the run list.
//
// For the Evaluate implementer
// ----------------------------
// Call the hooks below in this order and with these meanings; all arguments are plain types,
// so a hook type in src/ with the same member functions needs no test types.
//   on_initial_offer(x, v)       Alg 3 l.3, once per posted initial offer (before anything else).
//   on_relax_offer(arc, v)       every attempted decrease along arc `arc` of the INPUT state graph
//                                (s.hierarchy.graph.arcs index), v = emitted(tail) + l; report it
//                                even if the head is settled/retired or v does not decrease.
//   on_supplemental_offer(x, v)  Alg 2 l.6, v = max{b, z + S(i) + g(i, j)}, inside the scan.
//   on_open(node, states, D0, origin, nu, S)
//                                Alg 2 l.9 / Alg 3 l.9, internal nodes only. `states` are the
//                                states below the node (its leaf-layout interval); `node` is any
//                                id, used in messages only. For a root, origin = D0.
//   on_scan_begin(node, states, a', b', a, b, c)
//                                entry of Scan(X, I = [a', b'), nominal [a, b), c). For the root
//                                calls of Alg 3 l.10-11, [a', b') = [a, b).
//   on_scan_end(node, states)    return of that Scan (after Alg 2 l.12 retired its pieces).
//   on_emit(x, v)                Alg 2 l.11 / Alg 3 l.7, before relaxing the outgoing arcs.
//   on_retire(node, states)      any node or leaf being retired (Alg 2 l.10/l.12, expiry walk).
//                                One call on the expiring node covers every state below it.
//   on_finish()                  end of Evaluate.
// A normalised (contracted, Lemma 4.3) node has the same leaf set as the topmost input node it
// came from, so `states` identifies it in the input hierarchy without sharing NodeIds.
// EvalOptions{supplemental_offers} corresponds to CheckerOptions::supplemental_offers.

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "support/paper/eval_oracle.hpp"

namespace sssp::test::paper {

using sssp::paper::ArcId;
using sssp::paper::kNoArc;
using sssp::paper::kNoNode;
using sssp::paper::kNoState;
using sssp::paper::NodeId;

// ------------------------------------------------------------------ hook interface

/// Reference hook interface (virtual, all no-ops). See the header comment for the contract.
template <class W>
class EvalObserver {
 public:
  EvalObserver() = default;
  EvalObserver(const EvalObserver&) = default;
  EvalObserver(EvalObserver&&) noexcept = default;
  EvalObserver& operator=(const EvalObserver&) = default;
  EvalObserver& operator=(EvalObserver&&) noexcept = default;
  virtual ~EvalObserver() = default;

  virtual void on_initial_offer(StateId /*x*/, W /*value*/) {}
  virtual void on_relax_offer(ArcId /*arc*/, W /*value*/) {}
  virtual void on_supplemental_offer(StateId /*x*/, W /*value*/) {}
  virtual void on_open(std::uint32_t /*node*/, std::span<const StateId> /*states*/, W /*d0*/,
                       W /*origin*/, W /*unit*/, W /*mass*/) {}
  virtual void on_scan_begin(std::uint32_t /*node*/, std::span<const StateId> /*states*/, W /*lo*/,
                             W /*hi*/, W /*slot_lo*/, W /*slot_hi*/, W /*deadline*/) {}
  virtual void on_scan_end(std::uint32_t /*node*/, std::span<const StateId> /*states*/) {}
  virtual void on_emit(StateId /*x*/, W /*value*/) {}
  virtual void on_retire(std::uint32_t /*node*/, std::span<const StateId> /*states*/) {}
  virtual void on_finish() {}
};

/// Duck-typed form of the same contract, for an Evaluate that takes its observer as a template
/// parameter instead of a virtual base.
template <class O, class W>
concept EvalObserverFor =
    requires(O& o, StateId x, ArcId a, W v, std::uint32_t node, std::span<const StateId> states) {
      o.on_initial_offer(x, v);
      o.on_relax_offer(a, v);
      o.on_supplemental_offer(x, v);
      o.on_open(node, states, v, v, v, v);
      o.on_scan_begin(node, states, v, v, v, v, v);
      o.on_scan_end(node, states);
      o.on_emit(x, v);
      o.on_retire(node, states);
      o.on_finish();
    };

/// Adapter from a src-side virtual hook base `Hooks` whose member functions have the names and
/// argument lists of EvalObserver<W>: derives from Hooks and forwards every call to `target`.
/// If Hooks differs (a missing or renamed hook), `override` makes this fail to compile, which
/// is the intended early signal; write a dedicated adapter in that case.
template <class Hooks, class W>
class ForwardingObserver final : public Hooks {
 public:
  explicit ForwardingObserver(EvalObserver<W>& target) : target_(&target) {}

  void on_initial_offer(StateId x, W value) override { target_->on_initial_offer(x, value); }
  void on_relax_offer(ArcId arc, W value) override { target_->on_relax_offer(arc, value); }
  void on_supplemental_offer(StateId x, W value) override {
    target_->on_supplemental_offer(x, value);
  }
  void on_open(std::uint32_t node, std::span<const StateId> states, W d0, W origin, W unit,
               W mass) override {
    target_->on_open(node, states, d0, origin, unit, mass);
  }
  void on_scan_begin(std::uint32_t node, std::span<const StateId> states, W lo, W hi, W slot_lo,
                     W slot_hi, W deadline) override {
    target_->on_scan_begin(node, states, lo, hi, slot_lo, slot_hi, deadline);
  }
  void on_scan_end(std::uint32_t node, std::span<const StateId> states) override {
    target_->on_scan_end(node, states);
  }
  void on_emit(StateId x, W value) override { target_->on_emit(x, value); }
  void on_retire(std::uint32_t node, std::span<const StateId> states) override {
    target_->on_retire(node, states);
  }
  void on_finish() override { target_->on_finish(); }

 private:
  EvalObserver<W>* target_;
};

// ------------------------------------------------------------------ recorded traces

enum class EvalEventKind : std::uint8_t {
  kInitialOffer,
  kRelaxOffer,
  kSupplementalOffer,
  kOpen,
  kScanBegin,
  kScanEnd,
  kEmit,
  kRetire,
  kFinish,
};

/// One hook call. Field use per kind: offers/emit use state (or arc) and value; open uses node,
/// states, value = D0, origin, unit, mass; scan_begin uses node, states, lo, hi, slot_lo,
/// slot_hi, deadline; scan_end / retire use node, states.
template <class W>
struct EvalEvent {
  EvalEventKind kind{};
  StateId state = kNoState;
  ArcId arc = kNoArc;
  std::uint32_t node = 0;
  std::vector<StateId> states;
  W value{};
  W origin{};
  W unit{};
  W mass{};
  W lo{};
  W hi{};
  W slot_lo{};
  W slot_hi{};
  W deadline{};
};

template <class W>
EvalEvent<W> ev_initial(StateId x, W v) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kInitialOffer;
  e.state = x;
  e.value = v;
  return e;
}
template <class W>
EvalEvent<W> ev_relax(ArcId arc, W v) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kRelaxOffer;
  e.arc = arc;
  e.value = v;
  return e;
}
template <class W>
EvalEvent<W> ev_supplemental(StateId x, W v) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kSupplementalOffer;
  e.state = x;
  e.value = v;
  return e;
}
template <class W>
EvalEvent<W> ev_open(std::uint32_t node, std::vector<StateId> states, W d0, W origin, W unit,
                     W mass) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kOpen;
  e.node = node;
  e.states = std::move(states);
  e.value = d0;
  e.origin = origin;
  e.unit = unit;
  e.mass = mass;
  return e;
}
template <class W>
EvalEvent<W> ev_scan_begin(std::uint32_t node, std::vector<StateId> states, W lo, W hi, W slot_lo,
                           W slot_hi, W deadline) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kScanBegin;
  e.node = node;
  e.states = std::move(states);
  e.lo = lo;
  e.hi = hi;
  e.slot_lo = slot_lo;
  e.slot_hi = slot_hi;
  e.deadline = deadline;
  return e;
}
template <class W>
EvalEvent<W> ev_scan_end(std::uint32_t node, std::vector<StateId> states) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kScanEnd;
  e.node = node;
  e.states = std::move(states);
  return e;
}
template <class W>
EvalEvent<W> ev_emit(StateId x, W v) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kEmit;
  e.state = x;
  e.value = v;
  return e;
}
template <class W>
EvalEvent<W> ev_retire(std::uint32_t node, std::vector<StateId> states) {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kRetire;
  e.node = node;
  e.states = std::move(states);
  return e;
}
template <class W>
EvalEvent<W> ev_finish() {
  EvalEvent<W> e;
  e.kind = EvalEventKind::kFinish;
  return e;
}

/// Feeds a recorded trace to an observer.
template <class W>
void replay(std::span<const EvalEvent<W>> events, EvalObserver<W>& obs) {
  for (const auto& e : events) {
    switch (e.kind) {
      case EvalEventKind::kInitialOffer:
        obs.on_initial_offer(e.state, e.value);
        break;
      case EvalEventKind::kRelaxOffer:
        obs.on_relax_offer(e.arc, e.value);
        break;
      case EvalEventKind::kSupplementalOffer:
        obs.on_supplemental_offer(e.state, e.value);
        break;
      case EvalEventKind::kOpen:
        obs.on_open(e.node, e.states, e.value, e.origin, e.unit, e.mass);
        break;
      case EvalEventKind::kScanBegin:
        obs.on_scan_begin(e.node, e.states, e.lo, e.hi, e.slot_lo, e.slot_hi, e.deadline);
        break;
      case EvalEventKind::kScanEnd:
        obs.on_scan_end(e.node, e.states);
        break;
      case EvalEventKind::kEmit:
        obs.on_emit(e.state, e.value);
        break;
      case EvalEventKind::kRetire:
        obs.on_retire(e.node, e.states);
        break;
      case EvalEventKind::kFinish:
        obs.on_finish();
        break;
    }
  }
}

/// Records every hook call (for debugging a failing run, or to replay it into a checker).
template <class W>
class RecordingObserver final : public EvalObserver<W> {
 public:
  std::vector<EvalEvent<W>> events;

  void on_initial_offer(StateId x, W v) override { events.push_back(ev_initial<W>(x, v)); }
  void on_relax_offer(ArcId arc, W v) override { events.push_back(ev_relax<W>(arc, v)); }
  void on_supplemental_offer(StateId x, W v) override {
    events.push_back(ev_supplemental<W>(x, v));
  }
  void on_open(std::uint32_t node, std::span<const StateId> states, W d0, W origin, W unit,
               W mass) override {
    events.push_back(ev_open<W>(node, {states.begin(), states.end()}, d0, origin, unit, mass));
  }
  void on_scan_begin(std::uint32_t node, std::span<const StateId> states, W lo, W hi, W slot_lo,
                     W slot_hi, W deadline) override {
    events.push_back(
        ev_scan_begin<W>(node, {states.begin(), states.end()}, lo, hi, slot_lo, slot_hi, deadline));
  }
  void on_scan_end(std::uint32_t node, std::span<const StateId> states) override {
    events.push_back(ev_scan_end<W>(node, {states.begin(), states.end()}));
  }
  void on_emit(StateId x, W v) override { events.push_back(ev_emit<W>(x, v)); }
  void on_retire(std::uint32_t node, std::span<const StateId> states) override {
    events.push_back(ev_retire<W>(node, {states.begin(), states.end()}));
  }
  void on_finish() override { events.push_back(ev_finish<W>()); }
};

/// A plain state-graph Dijkstra written as a hook trace: the initial offers, then for every
/// state in dist_D order its emission followed by the relaxation of all its outgoing arcs, then
/// on_finish. It has no scans, so check it with CheckerOptions::allow_unscoped_emits. Under
/// Def 4.1 it satisfies I1 and I2 (protected states emit dist_D = delta).
template <class W>
std::vector<EvalEvent<W>> state_dijkstra_trace(const SeededOrderedHierarchy<W>& s,
                                               const EvalReference<W>& ref) {
  const auto& g = s.hierarchy.graph;
  std::vector<std::vector<ArcId>> out(g.num_states);
  for (ArcId i = 0; i < g.arcs.size(); ++i) out[g.arcs[i].tail].push_back(i);
  std::vector<EvalEvent<W>> ev;
  for (const auto& o : s.offers) ev.push_back(ev_initial<W>(o.state, o.value));
  for (const StateId x : ref.dijkstra_order) {
    const W d = ref.state_distance[x];
    ev.push_back(ev_emit<W>(x, d));
    for (const ArcId i : out[x]) {
      ev.push_back(ev_relax<W>(i, WeightTraits<W>::checked_add(d, g.arcs[i].length)));
    }
  }
  ev.push_back(ev_finish<W>());
  return ev;
}

// ------------------------------------------------------------------ violations

enum class Invariant : std::uint8_t { kI1, kI2, kI3, kI4, kI5, kI6, kI7, kI8, kProtocol };

enum class EvalViolationKind : std::uint8_t {
  // I1
  kOfferBelowTarget,  ///< an offer below delta(pi(x))
  // I2
  kEmitTwice,            ///< a state emitted more than once
  kEmitBelowTarget,      ///< an emission below delta(pi(x))
  kProtectedWrongValue,  ///< a protected state emitted a value != delta(pi(x))
  kProtectedNotEmitted,  ///< a protected state never emitted (at on_finish)
  // I3
  kOfferBelowPortion,      ///< an offer below the current portion's left end a'
  kSupplementalNotPadded,  ///< a supplemental offer below the nominal right end b
  // I4
  kEmitOutsidePortion,  ///< emission value not in [a', b') or above the deadline c
  kEmitOutsideNode,     ///< emission at a state not below the scanned node
  kPortionNotNested,    ///< child portion outside the parent portion / deadline above parent's
  kPortionOutsideSlot,  ///< portion [a', b') not inside the nominal slot [a, b), or empty
  kOpenNotDue,          ///< a node opened at a value not due in the parent's portion
  kSlotWentBack,        ///< a node's nominal slot starts before its previous one
  // I5
  kKeyLeftInPortion,  ///< key in I and <= c left below X at the end of Scan(X, I, c)
  // I6
  kRetiredProtectedUnemitted,  ///< node retired with an unemitted protected state below it
  kMassBelowDiameter,          ///< some state below X has delta > D0 + S(X)
  kDeadlineAfterExpiry,        ///< a scan deadline above the node's expiry D0 + S(X)
  // I7
  kOriginRule,               ///< origin differs from the rule of Sec 4.3 l.794
  kOriginOutsideParentSlot,  ///< origin outside the parent's nominal slot
  kSlotNotOnGrid,            ///< nominal slot not [origin + k nu, origin + (k+1) nu)
  // I8
  kOpenAboveProtected,  ///< opened at D0 > delta of an unemitted protected state below it
  // protocol
  kStateOutOfRange,       ///< state or arc id out of range
  kUnknownNode,           ///< `states` is not the leaf set of any input node
  kInitialOfferLate,      ///< initial offer after the evaluation started
  kInitialOfferUnknown,   ///< initial offer that is not one of the seeded offers
  kRelaxBeforeEmit,       ///< relaxation offer along an arc whose tail has not emitted
  kRelaxValue,            ///< relaxation offer != emitted(tail) + l
  kSupplementalDisabled,  ///< supplemental offer although they are disabled
  kSupplementalOutside,   ///< supplemental offer outside a scan or at a state not below X
  kEmitNotKey,            ///< emitted value != current key (least offer reported so far)
  kEmitAfterRetire,       ///< emission at a retired state
  kEmitOutsideScan,       ///< emission outside any scan at a state that is not a root leaf
  kOpenTwice,             ///< a node opened twice
  kOpenLeaf,              ///< on_open for a leaf
  kOpenNotMinimum,        ///< D0 != least key below the node
  kBadUnitOrMass,         ///< nu <= 0 or S < 0 / non-finite
  kScanUnopened,          ///< scan of a node that is not open (or already retired)
  kScanNotNested,         ///< scan / open of a node that is not strictly below the scanned node
  kScanEndMismatch,       ///< on_scan_end does not match the innermost scan / unclosed scans
  kEventAfterFinish,      ///< any hook after on_finish
};

[[nodiscard]] inline Invariant invariant_of(EvalViolationKind k) noexcept {
  using K = EvalViolationKind;
  switch (k) {
    case K::kOfferBelowTarget:
      return Invariant::kI1;
    case K::kEmitTwice:
    case K::kEmitBelowTarget:
    case K::kProtectedWrongValue:
    case K::kProtectedNotEmitted:
      return Invariant::kI2;
    case K::kOfferBelowPortion:
    case K::kSupplementalNotPadded:
      return Invariant::kI3;
    case K::kEmitOutsidePortion:
    case K::kEmitOutsideNode:
    case K::kPortionNotNested:
    case K::kPortionOutsideSlot:
    case K::kOpenNotDue:
    case K::kSlotWentBack:
      return Invariant::kI4;
    case K::kKeyLeftInPortion:
      return Invariant::kI5;
    case K::kRetiredProtectedUnemitted:
    case K::kMassBelowDiameter:
    case K::kDeadlineAfterExpiry:
      return Invariant::kI6;
    case K::kOriginRule:
    case K::kOriginOutsideParentSlot:
    case K::kSlotNotOnGrid:
      return Invariant::kI7;
    case K::kOpenAboveProtected:
      return Invariant::kI8;
    default:
      return Invariant::kProtocol;
  }
}

[[nodiscard]] inline const char* to_string(Invariant i) noexcept {
  static constexpr const char* kNames[] = {"I1", "I2", "I3", "I4",      "I5",
                                           "I6", "I7", "I8", "protocol"};
  return kNames[static_cast<std::size_t>(i)];
}

struct EvalViolation {
  EvalViolationKind kind;
  StateId state = kNoState;
  std::uint32_t node = 0;
  std::string message;

  [[nodiscard]] Invariant invariant() const noexcept { return invariant_of(kind); }
};

struct InvariantReport {
  std::vector<EvalViolation> violations;  ///< the first max_recorded violations
  std::size_t total = 0;

  [[nodiscard]] bool ok() const noexcept { return total == 0; }
  [[nodiscard]] bool has(EvalViolationKind k) const noexcept {
    return std::any_of(violations.begin(), violations.end(),
                       [k](const EvalViolation& v) { return v.kind == k; });
  }
  [[nodiscard]] bool has(Invariant i) const noexcept {
    return std::any_of(violations.begin(), violations.end(),
                       [i](const EvalViolation& v) { return v.invariant() == i; });
  }
  [[nodiscard]] std::size_t count(Invariant i) const noexcept {
    return static_cast<std::size_t>(
        std::count_if(violations.begin(), violations.end(),
                      [i](const EvalViolation& v) { return v.invariant() == i; }));
  }
  [[nodiscard]] std::string summary() const {
    if (ok()) return "OK";
    std::string s = std::to_string(total) + " violation(s)\n";
    for (const auto& v : violations) {
      s += std::string(to_string(v.invariant())) + ": " + v.message + "\n";
    }
    return s;
  }
};

struct CheckerOptions {
  bool supplemental_offers = true;    ///< EvalOptions::supplemental_offers
  bool check_relax_values = true;     ///< relaxation offer must equal emitted(tail) + l
  bool allow_unscoped_emits = false;  ///< accept emissions outside scans (flat traces)
  std::size_t max_recorded = 256;
};

// ------------------------------------------------------------------ checker

/// Checks I1-I8 (and the protocol) on the hook calls of one Evaluate run of `s`. Construct one
/// checker per run; the referenced hierarchy and reference must outlive it.
template <class W>
class InvariantChecker final : public EvalObserver<W> {
 public:
  InvariantChecker(const SeededOrderedHierarchy<W>& s, const EvalReference<W>& ref,
                   CheckerOptions opt = {})
      : s_(&s), ref_(&ref), opt_(opt) {
    const auto& h = s.hierarchy;
    const std::size_t n = h.num_states();
    key_.assign(n, inf());
    emitted_.assign(n, false);
    emitted_value_.assign(n, W{});
    retired_.assign(n, false);
    initial_values_.resize(n);
    for (const auto& o : s.offers) {
      if (o.state < n) initial_values_[o.state].push_back(o.value);
    }
    // Leaf sets of the input nodes (unique per node: internal nodes have >= 2 children).
    const std::size_t nn = h.nodes.size();
    tin_.assign(nn, 0);
    tout_.assign(nn, 0);
    node_state_.resize(nn);
    std::size_t clock = 0;
    std::vector<std::pair<NodeId, std::size_t>> stack;
    for (const NodeId r : h.roots) {
      stack.emplace_back(r, 0);
      tin_[r] = clock++;
      while (!stack.empty()) {
        auto& [v, i] = stack.back();
        if (i < h.nodes[v].children.size()) {
          const NodeId c = h.nodes[v].children[i++];
          tin_[c] = clock++;
          stack.emplace_back(c, 0);
        } else {
          tout_[v] = clock;
          stack.pop_back();
        }
      }
    }
    for (NodeId v = 0; v < nn; ++v) {
      std::vector<StateId> leaves;
      for (StateId x = 0; x < n; ++x) {
        if (under(v, x)) leaves.push_back(x);
      }
      by_leaves_.emplace(std::move(leaves), v);
    }
  }

  [[nodiscard]] const InvariantReport& report() const noexcept { return report_; }
  [[nodiscard]] bool emitted(StateId x) const { return emitted_.at(x); }
  [[nodiscard]] std::vector<std::pair<StateId, W>> emissions() const { return emissions_; }

  void on_initial_offer(StateId x, W value) override {
    if (!live() || !state_ok(x, "initial offer")) return;
    if (started_) {
      add(EvalViolationKind::kInitialOfferLate, x, 0, "initial offer after evaluation started");
    }
    const auto& vals = initial_values_[x];
    if (std::find(vals.begin(), vals.end(), value) == vals.end()) {
      add(EvalViolationKind::kInitialOfferUnknown, x, 0,
          "initial offer " + fmt(value) + " is not a seeded offer");
    }
    offer(x, value, "initial");
  }

  void on_relax_offer(ArcId arc, W value) override {
    if (!live()) return;
    const auto& g = s_->hierarchy.graph;
    if (arc >= g.arcs.size()) {
      add(EvalViolationKind::kStateOutOfRange, kNoState, 0, "relaxation along unknown arc");
      return;
    }
    started_ = true;
    const auto& a = g.arcs[arc];
    if (!emitted_[a.tail]) {
      add(EvalViolationKind::kRelaxBeforeEmit, a.tail, 0,
          "relaxation of arc " + std::to_string(arc) + " before its tail emitted");
    } else if (opt_.check_relax_values) {
      const W expect = sat_add(emitted_value_[a.tail], a.length);
      if (expect != value) {
        add(EvalViolationKind::kRelaxValue, a.head, 0,
            "relaxation of arc " + std::to_string(arc) + " offers " + fmt(value) + ", expected " +
                fmt(expect));
      }
    }
    offer(a.head, value, "relaxation");
  }

  void on_supplemental_offer(StateId x, W value) override {
    if (!live() || !state_ok(x, "supplemental offer")) return;
    started_ = true;
    if (!opt_.supplemental_offers) {
      add(EvalViolationKind::kSupplementalDisabled, x, 0, "supplemental offer while disabled");
    }
    if (scans_.empty()) {
      add(EvalViolationKind::kSupplementalOutside, x, 0, "supplemental offer outside a scan");
    } else {
      const Frame& f = scans_.back();
      if (!under(f.node, x)) {
        add(EvalViolationKind::kSupplementalOutside, x, f.id,
            "supplemental offer at a state not below the scanned node");
      }
      if (value < f.slot_hi) {
        add(EvalViolationKind::kSupplementalNotPadded, x, f.id,
            "supplemental offer " + fmt(value) + " below the nominal slot end " + fmt(f.slot_hi));
      }
    }
    offer(x, value, "supplemental");
  }

  void on_open(std::uint32_t id, std::span<const StateId> states, W d0, W origin, W unit,
               W mass) override {
    if (!live()) return;
    started_ = true;
    const NodeId v = resolve(id, states);
    if (v == kNoNode) return;
    const auto& h = s_->hierarchy;
    if (h.is_leaf(v)) {
      add(EvalViolationKind::kOpenLeaf, kNoState, id, "on_open for a leaf");
      return;
    }
    NodeState& ns = node_state_[v];
    if (ns.opened) {
      add(EvalViolationKind::kOpenTwice, kNoState, id, "node opened twice");
      return;
    }
    if (!(unit > W{0}) || !finite(unit) || !finite(mass) || mass < W{0}) {
      add(EvalViolationKind::kBadUnitOrMass, kNoState, id,
          "invalid unit " + fmt(unit) + " or mass " + fmt(mass));
      return;
    }
    ns.opened = true;
    ns.d0 = d0;
    ns.origin = origin;
    ns.unit = unit;
    ns.mass = mass;
    ns.expiry = sat_add(d0, mass);
    // Least key below the node and per-state checks.
    W least = inf();
    for (const StateId x : states) {
      if (!emitted_[x] && !retired_[x] && key_[x] < least) least = key_[x];
      if (emitted_[x]) {
        add(EvalViolationKind::kOpenNotMinimum, x, id, "a state below the node already emitted");
      }
      const W t = ref_->state_target[x];
      if (t != inf() && ns.expiry < t) {
        add(EvalViolationKind::kMassBelowDiameter, x, id,
            "delta " + fmt(t) + " > D0 + S = " + fmt(ns.expiry));
      }
      if (ref_->is_protected[x] && !emitted_[x] && t < d0) {
        add(EvalViolationKind::kOpenAboveProtected, x, id,
            "opened at D0 = " + fmt(d0) + " above protected delta " + fmt(t));
      }
    }
    if (least != d0) {
      add(EvalViolationKind::kOpenNotMinimum, kNoState, id,
          "D0 = " + fmt(d0) + " but the least key below the node is " + fmt(least));
    }
    if (scans_.empty()) {
      if (h.nodes[v].parent != kNoNode) {
        add(EvalViolationKind::kScanNotNested, kNoState, id, "non-root opened outside a scan");
      }
      if (origin != d0) {
        add(EvalViolationKind::kOriginRule, kNoState, id,
            "root origin " + fmt(origin) + " != D0 " + fmt(d0));
      }
    } else {
      const Frame& f = scans_.back();
      if (!strictly_below(v, f.node)) {
        add(EvalViolationKind::kScanNotNested, kNoState, id,
            "opened node is not below the scanned node");
      }
      if (!(f.lo <= d0 && d0 < f.hi && d0 <= f.deadline)) {
        add(EvalViolationKind::kOpenNotDue, kNoState, id,
            "D0 = " + fmt(d0) + " not due in the parent's portion");
      }
      check_origin(id, d0, origin, unit, mass, f.slot_lo, f.slot_hi);
    }
    if (!(origin <= d0 && d0 < sat_add(origin, unit))) {
      add(EvalViolationKind::kOriginRule, kNoState, id,
          "D0 = " + fmt(d0) + " not in the first slot [" + fmt(origin) + ", +nu)");
    }
  }

  void on_scan_begin(std::uint32_t id, std::span<const StateId> states, W lo, W hi, W slot_lo,
                     W slot_hi, W deadline) override {
    if (!live()) return;
    started_ = true;
    const NodeId v = resolve(id, states);
    if (v == kNoNode) return;
    NodeState& ns = node_state_[v];
    if (!ns.opened || ns.retired) {
      add(EvalViolationKind::kScanUnopened, kNoState, id, "scan of a node that is not open");
    }
    if (scans_.empty()) {
      if (s_->hierarchy.nodes[v].parent != kNoNode) {
        add(EvalViolationKind::kScanNotNested, kNoState, id, "non-root scanned outside a scan");
      }
    } else {
      const Frame& f = scans_.back();
      if (!strictly_below(v, f.node)) {
        add(EvalViolationKind::kScanNotNested, kNoState, id,
            "scanned node is not below the parent scan's node");
      }
      if (lo < f.lo || f.hi < hi || f.deadline < deadline) {
        add(EvalViolationKind::kPortionNotNested, kNoState, id,
            "portion [" + fmt(lo) + ", " + fmt(hi) + ") deadline " + fmt(deadline) +
                " not inside the parent's [" + fmt(f.lo) + ", " + fmt(f.hi) + ") deadline " +
                fmt(f.deadline));
      }
    }
    if (!(lo < hi) || lo < slot_lo || slot_hi < hi) {
      add(EvalViolationKind::kPortionOutsideSlot, kNoState, id,
          "portion [" + fmt(lo) + ", " + fmt(hi) + ") not a nonempty part of slot [" +
              fmt(slot_lo) + ", " + fmt(slot_hi) + ")");
    }
    if (ns.opened) {
      if (ns.expiry < deadline) {
        add(EvalViolationKind::kDeadlineAfterExpiry, kNoState, id,
            "deadline " + fmt(deadline) + " after expiry " + fmt(ns.expiry));
      }
      if (sat_add(slot_lo, ns.unit) != slot_hi || !on_grid(ns.origin, ns.unit, slot_lo)) {
        add(EvalViolationKind::kSlotNotOnGrid, kNoState, id,
            "slot [" + fmt(slot_lo) + ", " + fmt(slot_hi) + ") not on the grid origin " +
                fmt(ns.origin) + " + k * " + fmt(ns.unit));
      }
      if (ns.has_slot && slot_lo < ns.last_slot_lo) {
        add(EvalViolationKind::kSlotWentBack, kNoState, id,
            "slot " + fmt(slot_lo) + " before previous slot " + fmt(ns.last_slot_lo));
      }
      ns.has_slot = true;
      ns.last_slot_lo = slot_lo;
    }
    scans_.push_back(Frame{v, id, lo, hi, slot_lo, slot_hi, deadline});
  }

  void on_scan_end(std::uint32_t id, std::span<const StateId> states) override {
    if (!live()) return;
    const NodeId v = resolve(id, states);
    if (scans_.empty() || (v != kNoNode && scans_.back().node != v)) {
      add(EvalViolationKind::kScanEndMismatch, kNoState, id,
          "on_scan_end does not match the innermost scan");
      return;
    }
    const Frame f = scans_.back();
    scans_.pop_back();
    for (StateId x = 0; x < key_.size(); ++x) {
      if (!under(f.node, x) || emitted_[x] || retired_[x]) continue;
      const W k = key_[x];
      if (f.lo <= k && k < f.hi && k <= f.deadline) {
        add(EvalViolationKind::kKeyLeftInPortion, x, f.id,
            "key " + fmt(k) + " left in [" + fmt(f.lo) + ", " + fmt(f.hi) + ") with deadline " +
                fmt(f.deadline));
      }
    }
  }

  void on_emit(StateId x, W value) override {
    if (!live() || !state_ok(x, "emission")) return;
    started_ = true;
    if (emitted_[x]) {
      add(EvalViolationKind::kEmitTwice, x, 0, "state emitted twice");
      return;
    }
    if (retired_[x]) add(EvalViolationKind::kEmitAfterRetire, x, 0, "emission at a retired state");
    if (value != key_[x]) {
      add(EvalViolationKind::kEmitNotKey, x, 0,
          "emitted " + fmt(value) + " but the current key is " + fmt(key_[x]));
    }
    const W t = ref_->state_target[x];
    if (value < t) {
      add(EvalViolationKind::kEmitBelowTarget, x, 0,
          "emitted " + fmt(value) + " below delta " + fmt(t));
    }
    if (ref_->is_protected[x] && value != t) {
      add(EvalViolationKind::kProtectedWrongValue, x, 0,
          "protected state emitted " + fmt(value) + ", delta is " + fmt(t));
    }
    if (scans_.empty()) {
      const auto& h = s_->hierarchy;
      const NodeId leaf = h.leaf_of_state[x];
      if (!opt_.allow_unscoped_emits && h.nodes[leaf].parent != kNoNode) {
        add(EvalViolationKind::kEmitOutsideScan, x, 0,
            "emission outside any scan at a state that is not a root leaf");
      }
    } else {
      const Frame& f = scans_.back();
      if (!under(f.node, x)) {
        add(EvalViolationKind::kEmitOutsideNode, x, f.id,
            "emission at a state not below the scanned node");
      }
      if (!(f.lo <= value && value < f.hi && value <= f.deadline)) {
        add(EvalViolationKind::kEmitOutsidePortion, x, f.id,
            "emitted " + fmt(value) + " outside [" + fmt(f.lo) + ", " + fmt(f.hi) +
                ") with deadline " + fmt(f.deadline));
      }
    }
    emitted_[x] = true;
    emitted_value_[x] = value;
    emissions_.emplace_back(x, value);
  }

  void on_retire(std::uint32_t id, std::span<const StateId> states) override {
    if (!live()) return;
    started_ = true;
    const NodeId v = resolve(id, states);
    if (v == kNoNode) return;
    node_state_[v].retired = true;
    for (const StateId x : states) {
      if (emitted_[x]) continue;
      if (ref_->is_protected[x]) {
        add(EvalViolationKind::kRetiredProtectedUnemitted, x, id,
            "retired before protected state emitted (delta " + fmt(ref_->state_target[x]) + ")");
      }
      retired_[x] = true;
    }
  }

  void on_finish() override {
    if (!live()) return;
    finished_ = true;
    if (!scans_.empty()) {
      add(EvalViolationKind::kScanEndMismatch, kNoState, scans_.back().id,
          "evaluation finished inside a scan");
    }
    for (StateId x = 0; x < emitted_.size(); ++x) {
      if (ref_->is_protected[x] && !emitted_[x]) {
        add(EvalViolationKind::kProtectedNotEmitted, x, 0,
            "protected state never emitted (delta " + fmt(ref_->state_target[x]) + ")");
      }
    }
  }

 private:
  struct NodeState {
    bool opened = false;
    bool retired = false;
    bool has_slot = false;
    W d0{};
    W origin{};
    W unit{};
    W mass{};
    W expiry{};
    W last_slot_lo{};
  };
  struct Frame {
    NodeId node;
    std::uint32_t id;
    W lo, hi, slot_lo, slot_hi, deadline;
  };

  static constexpr W inf() noexcept { return WeightTraits<W>::infinity(); }
  static bool finite(W w) noexcept {
    if constexpr (std::is_floating_point_v<W>) {
      return std::isfinite(w);
    } else {
      return w != inf();
    }
  }
  /// a + b, saturating at infinity (the checker must not throw on a faulty trace).
  static W sat_add(W a, W b) noexcept {
    if (a == inf() || b == inf()) return inf();
    if constexpr (std::is_floating_point_v<W>) {
      return a + b;
    } else {
      return a >= inf() - b ? inf() : static_cast<W>(a + b);
    }
  }
  static std::string fmt(W w) {
    if (w == inf()) return "inf";
    if constexpr (std::is_floating_point_v<W>) {
      std::ostringstream os;
      os.precision(17);
      os << w;
      return os.str();
    } else {
      return std::to_string(w);
    }
  }

  [[nodiscard]] bool under(NodeId v, StateId x) const {
    const NodeId leaf = s_->hierarchy.leaf_of_state[x];
    return tin_[v] <= tin_[leaf] && tin_[leaf] < tout_[v];
  }
  [[nodiscard]] bool strictly_below(NodeId v, NodeId anc) const {
    return v != anc && tin_[anc] <= tin_[v] && tin_[v] < tout_[anc];
  }
  /// nu-grid membership: slot_lo = origin + k * unit for some k >= 0.
  [[nodiscard]] static bool on_grid(W origin, W unit, W slot_lo) {
    if (slot_lo < origin) return false;
    if constexpr (std::is_integral_v<W>) {
      return (slot_lo - origin) % unit == 0;
    } else {
      W t = origin;
      for (std::size_t step = 0; step < (std::size_t{1} << 24) && t < slot_lo; ++step) t += unit;
      return t == slot_lo;
    }
  }

  void check_origin(std::uint32_t id, W d0, W origin, W unit, W mass, W a, W b) {
    if (!(a <= origin && origin < b)) {
      add(EvalViolationKind::kOriginOutsideParentSlot, kNoState, id,
          "origin " + fmt(origin) + " outside the parent's slot [" + fmt(a) + ", " + fmt(b) + ")");
    }
    if (sat_add(d0, mass) < b) {
      if (origin != d0) {
        add(EvalViolationKind::kOriginRule, kNoState, id,
            "D0 + S < b, so origin must be D0 = " + fmt(d0) + ", got " + fmt(origin));
      }
      return;
    }
    // Largest b - j * unit <= d0 by successive subtraction (Sec 4.3: O(1 + S/nu) steps).
    if (!(d0 < b)) return;  // reported as kOpenNotDue
    W t = b;
    if constexpr (std::is_integral_v<W>) {
      const W j = (b - d0) / unit + ((b - d0) % unit == 0 ? 0 : 1);
      if (j > b / unit) {  // b - j * unit < 0: no such origin in an unsigned grid
        add(EvalViolationKind::kOriginOutsideParentSlot, kNoState, id,
            "the rule's origin b - j nu would be negative");
        return;
      }
      t = b - j * unit;
    } else {
      for (std::size_t step = 0; step < (std::size_t{1} << 24) && d0 < t; ++step) t -= unit;
    }
    if (origin != t) {
      add(EvalViolationKind::kOriginRule, kNoState, id,
          "D0 + S >= b, so origin must be the largest b - j nu <= D0 = " + fmt(t) + ", got " +
              fmt(origin));
    }
  }

  void offer(StateId x, W value, const char* what) {
    const W t = ref_->state_target[x];
    if (value < t) {
      add(EvalViolationKind::kOfferBelowTarget, x, 0,
          std::string(what) + " offer " + fmt(value) + " below delta " + fmt(t));
    }
    if (!scans_.empty() && value < scans_.back().lo) {
      add(EvalViolationKind::kOfferBelowPortion, x, scans_.back().id,
          std::string(what) + " offer " + fmt(value) + " below the portion start " +
              fmt(scans_.back().lo));
    }
    if (value < key_[x]) key_[x] = value;
  }

  NodeId resolve(std::uint32_t id, std::span<const StateId> states) {
    std::vector<StateId> sorted(states.begin(), states.end());
    std::sort(sorted.begin(), sorted.end());
    const auto it = by_leaves_.find(sorted);
    if (it == by_leaves_.end()) {
      add(EvalViolationKind::kUnknownNode, kNoState, id,
          "states (" + std::to_string(states.size()) + ") are not the leaf set of an input node");
      return kNoNode;
    }
    return it->second;
  }

  bool state_ok(StateId x, const char* what) {
    if (x < key_.size()) return true;
    add(EvalViolationKind::kStateOutOfRange, x, 0, std::string(what) + " at an unknown state");
    return false;
  }

  bool live() {
    if (!finished_) return true;
    add(EvalViolationKind::kEventAfterFinish, kNoState, 0, "hook called after on_finish");
    return false;
  }

  void add(EvalViolationKind k, StateId x, std::uint32_t node, std::string msg) {
    ++report_.total;
    if (report_.violations.size() >= opt_.max_recorded) return;
    if (x != kNoState) msg = "state " + std::to_string(x) + ": " + msg;
    msg = "node " + std::to_string(node) + ": " + msg;
    report_.violations.push_back(EvalViolation{k, x, node, std::move(msg)});
  }

  const SeededOrderedHierarchy<W>* s_;
  const EvalReference<W>* ref_;
  CheckerOptions opt_;
  InvariantReport report_;

  std::vector<W> key_;
  std::vector<bool> emitted_;
  std::vector<W> emitted_value_;
  std::vector<bool> retired_;
  std::vector<std::vector<W>> initial_values_;
  std::vector<std::pair<StateId, W>> emissions_;

  std::vector<std::size_t> tin_, tout_;
  std::map<std::vector<StateId>, NodeId> by_leaves_;
  std::vector<NodeState> node_state_;
  std::vector<Frame> scans_;
  bool started_ = false;
  bool finished_ = false;
};

static_assert(EvalObserverFor<InvariantChecker<std::uint64_t>, std::uint64_t>);
static_assert(EvalObserverFor<EvalObserver<double>, double>);

}  // namespace sssp::test::paper
