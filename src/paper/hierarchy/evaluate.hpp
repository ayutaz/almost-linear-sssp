#pragma once
// Evaluation of a seeded ordered hierarchy: Algorithm 2 (Scan) and Algorithm 3 (Evaluate) of
// Sec 4.3 (l.781-882), with the correctness claim of Lemma 4.6 and the visit accounting of
// Lemma 4.7.
//
// evaluate() runs, in this order:
//   Alg 3 l.2  normalize() (Sec 4.1, Lemma 4.3) and order_children() (Sec 4.2, Lemma 4.5); the
//              leaf layout is the DFS order of the leaves with children in rank order and the
//              roots in root_order, so every node, every child and every run of consecutive
//              children is an interval of the layout. One split-findmin structure
//              (Lemma 3.11, SplitFindmin) holds a key per layout position, all +infinity.
//   Alg 3 l.3  splits the layout at every root boundary and posts the least initial offer of
//              every state (an input scan keeps the minimum per state, Lemma 4.7).
//   Alg 3 l.4  processes the roots in root_order; a root without a finite key is skipped
//              (l.5), a one-state root emits its key (l.6-7), any other root opens at its least
//              key D0 with origin D0 and expiry c = D0 + S(X) (l.9) and runs Scan over the
//              successive width-nu(X) slots starting at D0 through the slot containing c, until
//              all its children finish (l.10-11).
//   Alg 2      Scan(X, I = [a', b'), nominal slot [a, b), inclusive deadline c) walks X's list
//              of pieces (registered children and maximal runs of unregistered children) from
//              its first piece in rank order (l.2). A run whose interval minimum z is due
//              (a' <= z < b', z <= c) is replaced by its nonempty left subrun, the child i of the
//              minimising leaf (level ancestor, Lemma 3.9) and its nonempty right subrun, with
//              at most two interval splits (l.3-5); each new subrun gets the supplemental offer
//              max{b, z + S(i) + g(i, j)} at the representative of its child j next to i
//              (l.6, Sec 4.3 l.824-826); the cursor resumes at the first new piece (l.7). A
//              registered internal child is opened at its least key if necessary, with the
//              two-case origin rule of Sec 4.3 l.794-797, and its live slots meeting I through
//              its deadline min(c, D0 + S(i)) are visited recursively (l.9-10). A registered leaf
//              emits its due key once and relaxes its outgoing arcs (l.11). Completed children
//              and children whose deadline is covered are retired and unlinked (l.12).
//
// Interpretations (open questions in docs/research/open-questions.md):
//   [Interpretation, Q-27 / research notes 02-A5] Slot advancement and resumption. Every
//     Scan(X, I, c) walks X's list from its first live piece ("restart each intersected slot at
//     its first live child rank", Lemma 7.7 l.1884-1885). A registered internal child keeps its
//     slot index across visits ("its clock persists across visits", Sec 4.3 l.861) and a visit
//     covers every slot J of the child with J meeting I and start(J) <= its deadline ("advances
//     through all its local slots intersecting the current visit", l.863-864; Alg 2 l.10). After
//     Scan(i, J n I, c_i) the slot index advances iff J ends inside I (end(J) <= b'); otherwise
//     J stays the child's current slot and the next portion resumes it, again from the child's
//     first live piece ("retain registrations, cursors' nominal slot identities and unapplied
//     offers", Lemma 7.7 l.1883-1884). The child is retired after the visit iff it completed or
//     its deadline is covered (c_i < b'). In Sec 4 alone a clipped slot is never resumed
//     ("No future boundary cuts a local slot while physically relevant events remain beyond that
//     boundary", l.799-801): the counter eval_resumed_slot_visits measures this.
//   [Interpretation, Q-28 / 02-A6, A7] The b of the origin rule and of the supplemental padding
//     is the right end of the parent's NOMINAL slot, not of the current portion (Lemma 7.7
//     l.1887, l.1890).
//   [Interpretation, Q-29 / 02-A9] Only protected states are promised delta (Lemma 4.6); other
//     states emit a value >= delta or nothing (they may retire unemitted).
//   [Interpretation, Q-16 / 02-A17] "Place settled leaves in singleton intervals" (l.859): a
//     registered leaf already is a singleton interval and is not touched after it emits; no key
//     is ever raised.
//
// Model of computation (Sec 2, l.292-311): lengths are copied, compared and added (counted into
// Diagnostics::ops). [Deviation candidate, proposed D-16] The second origin rule's b - j nu is
// materialised as ONE real by j successive subtractions (the paper: "as signed pairs of sums",
// l.795-796). Exact for integer W; for double each subtraction is rounded (integer-valued
// doubles below 2^53 are exact). With rounding an exact tie D0 = b - j nu can come out one ulp
// above D0, so the subtractions would leave the parent's slot, which exact arithmetic excludes
// (Sec 4.3 l.798-800): for double the origin is then D0 (the first rule's choice; the grid is not
// aligned to b, which the slot-resumption rule of Q-27 tolerates) and the counter
// eval_origin_rounding_fallbacks records it; for integer W the same situation means Def 4.1 is
// violated and std::logic_error is thrown. The gaps g(i, j) come materialised from order_children
// (proposed D-15, child_order.hpp); with a gap a few ulps short a supplemental offer can be
// a few ulps below delta for non-integer doubles.
//
// The split-findmin structure is the two-level SplitFindmin of Lemma 3.11 (not
// SplitFindminRef, so no D-3). LCA / level ancestors are ForestQueries (Lemma 3.9).

#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

/// Observation hooks of Evaluate, all no-ops by default. The names and argument lists are the
/// observer contract of the Sec 4 test support (EvalObserver): tests derive an adapter
/// (ForwardingObserver<EvalHooks<W>, W>) that feeds the I1-I8 invariant checker.
///
/// Call order and meaning (node ids are NodeIds of the normalised forest, messages only;
/// `states` is the node's interval of the leaf layout):
///   on_initial_offer(x, v)        Alg 3 l.3, the least initial offer of x, before anything else
///   on_relax_offer(arc, v)        every relaxation of an arc of the INPUT state graph (index into
///                                 graph.arcs), v = emitted(tail) + l, also when it is rejected
///   on_supplemental_offer(x, v)   Alg 2 l.6, v = max{b, z + S(i) + g(i, j)}
///   on_open(node, states, D0, origin, nu, S)   Alg 3 l.9 / Alg 2 l.9 (internal nodes)
///   on_scan_begin(node, states, a', b', a, b, c) / on_scan_end(node, states)   one Scan call
///   on_emit(x, v)                 Alg 3 l.7 / Alg 2 l.11, before the outgoing arcs are relaxed
///   on_retire(node, states)       a child disposed of after completion or expiry, or a root
///                                 that expired with pieces left; one call covers its subtree
///   on_finish()                   end of Evaluate
template <class W>
class EvalHooks {
 public:
  EvalHooks() = default;
  EvalHooks(const EvalHooks&) = default;
  EvalHooks(EvalHooks&&) noexcept = default;
  EvalHooks& operator=(const EvalHooks&) = default;
  EvalHooks& operator=(EvalHooks&&) noexcept = default;
  virtual ~EvalHooks() = default;

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

template <class W>
struct EvalOptions {
  /// Post the supplemental offers of Alg 2 l.6. false is the test-only switch D-14: Lemma 4.6
  /// does not use them (correctness is unchanged), only the visit bound of Lemma 4.7 is lost.
  bool supplemental_offers = true;
  /// Observation hooks; nullptr costs one pointer test per event.
  EvalHooks<W>* hooks = nullptr;
};

template <class W>
struct EvalResult {
  /// Emitted value per state; WeightTraits<W>::infinity() where the state never emitted.
  /// Protected states (Def 4.2) carry delta(pi(x)) (Lemma 4.6).
  std::vector<W> value;
  /// (state, value) in emission order.
  std::vector<std::pair<StateId, W>> emissions;
  /// Counters of normalize, order_children and Evaluate (names in evaluate.cpp), plus length
  /// operations in diag.ops. Evaluate has no fallback: diag.fallbacks stays empty.
  Diagnostics diag;
};

/// Alg 3 on a seeded ordered hierarchy that satisfies Def 4.1 / 4.2 (the caller's obligation;
/// validate_seeded_ordered_hierarchy checks it). Only s.offers are read; the sources define
/// delta for the caller's checks. Throws std::out_of_range for an offer at a non-existent state
/// and std::invalid_argument for an offer that is not finite and nonnegative, plus whatever
/// normalize() / order_children() throw on a malformed hierarchy.
template <class W>
EvalResult<W> evaluate(const SeededOrderedHierarchy<W>& s, const EvalOptions<W>& options = {});

/// The same on an already normalised and ordered hierarchy (nh.children_ordered must be true),
/// for callers that reuse one normalisation. Counters are added to `diag` if non-null; the
/// returned EvalResult::diag then stays empty.
template <class W>
EvalResult<W> evaluate_normalized(const NormalizedHierarchy<W>& nh,
                                  std::span<const StateOffer<W>> offers,
                                  const EvalOptions<W>& options = {}, Diagnostics* diag = nullptr);

}  // namespace sssp::paper
