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
// PortionedEvaluator (end of this file) runs the same Scan state machine cell by cell, as
// Lemma 7.7 and the atom visits of Sec 7.4 (l.1835-1841) require: only the root loop differs.
//
// Interpretations (open questions in docs/research/open-questions.md):
//   [Interpretation, Q-27 / research notes 02-A5] Slot advancement and resumption. Every
//     Scan(X, I, c) walks X's list from its first live piece ("restart each intersected slot at
//     its first live child rank", Lemma 7.7 l.1885). A registered internal child keeps its
//     slot index across visits ("its clock persists across visits", Sec 4.3 l.861) and a visit
//     covers every slot J of the child with J meeting I and start(J) <= its deadline ("advances
//     through all its local slots intersecting the current visit", l.863-864; Alg 2 l.10). After
//     Scan(i, J n I, c_i) the slot index advances iff J ends inside I (end(J) <= b'); otherwise
//     J stays the child's current slot and the next portion resumes it, again from the child's
//     first live piece ("retain registrations, cursors' nominal slot identities and unapplied
//     offers", Lemma 7.7 l.1884-1885). The child is retired after the visit iff it completed or
//     its deadline is covered (c_i < b'). In Sec 4 alone a clipped slot is never resumed with
//     exact arithmetic ("No future boundary cuts a local slot while physically relevant events
//     remain beyond that boundary", l.799-801): the counter eval_resumed_slot_visits measures
//     this; it is nonzero under the external cells of PortionedEvaluator and can be nonzero for
//     non-integer doubles, whose rounded slot grids need not align with the parent's (see the
//     origin fallback D-19 below).
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
// Diagnostics::ops). [Deviation D-19] The second origin rule's b - j nu is
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

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/primitives/split_findmin.hpp"
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
///   on_external_offer(x, v)       PortionedEvaluator only: an external offer of Lemma 7.7
///                                 (between cells), also when it is rejected. Not part of the
///                                 Sec 4 observer contract (tests/support/paper), which has no
///                                 external offers; ForwardingObserver leaves it a no-op.
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
  virtual void on_external_offer(StateId /*x*/, W /*value*/) {}
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

/// Cell-by-cell Evaluate of Lemma 7.7 (Sec 7.4, l.1835-1841 "Visiting an active atom",
/// l.1871-1909, Figure 2), the API OwnerDistances (Alg 7 l.13) calls once per marker cell of
/// an active atom. It runs the same Scan state machine (Alg 2) as evaluate(); only the root
/// loop differs:
///   * the hierarchy has exactly one root (Lemma 7.7; Alg 7 caps every atom to a single root,
///     Lemmas 5.13-5.14) and must be normalised and ordered (Sec 4.1-4.2, Alg 7 l.7);
///   * visit_cell([L, H)) visits one external cell. A singleton root emits its due minimum
///     once (l.1835). An internal root opens on its first due key D0 with origin D0 ("an empty
///     cell leaves an unopened root unopened"); its inclusive cap is c = min(deadline,
///     D0 + S(root)) (l.1837-1838). The cell then visits the successive nominal root slots
///     whose intersections with [L, H) contain times <= c, calling Scan(root, J n [L, H), c)
///     with the nominal slot J (l.1838-1839, l.1883-1884);
///   * [Interpretation, Q-27] a slot J (root or child) clipped at H stays the node's current
///     slot; the next cell resumes it, and every Scan walks the node's list from its first live
///     piece ("On resuming any clipped slot, start at its first live child rank with its
///     registrations, offers, and deadlines intact", l.1840-1841; l.1884-1885). Registrations,
///     cursors' nominal slot identities, unapplied offers, origins (with b the parent's NOMINAL
///     slot end, Q-28) and supplemental padding (nominal b, l.1890) are kept across cells. A
///     child is retired only when it completed or its deadline is covered by the current portion
///     (c_i < b'), so "No temporary cut retires a node" (l.1891);
///   * the root stops once all its children finish or the cap is covered including equality
///     (c < H, l.1841); then finished() is true and later cells are no-ops.
/// Emissions relax the hierarchy's own arcs inside the visit; the caller relaxes the outside
/// exits of each emitted state after visit_cell returns (l.1839-1840, "relax the allowed exits
/// on emission"; interatom arrivals belong to later cells, Lemma 7.8).
///
/// Offers: initial_offer() before the first cell (Alg 7 l.9; Lemma 7.7's Q), external_offer()
/// between cells (Lemma 7.7's E). The arrival hypothesis of Lemma 7.7 ("arrives before the visit
/// to its timestamp's cell, and is never inserted below the current cell's left endpoint") is
/// enforced: an external offer that would decrease a pending key below the end of the last
/// visited cell throws std::invalid_argument; one at a settled or retired state is ignored
/// (Lemma 7.8: "a settled or retired head ignores later proposals"). A finite pending key below
/// the first visited cell at an unopened root (and not above the deadline) throws
/// std::logic_error ("The first cell starts at or before every finite initial timestamp").
///
/// Cells: lo < hi with lo finite (hi may be +infinity for a last open cell); cells never
/// overlap or go back; once the root has opened they must be consecutive (lo == previous hi,
/// "Visit it in consecutive cells"). Before that, skipped cells count as empty visits (Alg 7
/// visits an atom only from its registration on). Lemma 7.7's work bound assumes cells of a grid
/// [j q, (j + 1) q) with q >= every unit; other cells are accepted (a slot is then scanned once
/// per cell it meets).
///
/// Layout: by default the evaluator owns a split-findmin structure (Lemma 3.11) over its own
/// leaf layout. The second constructor places the layout at [base, base + num_states()) of a
/// caller's structure, so that several atoms share "one leaf layout in family order" (Alg 7
/// l.8, Lemma 7.9) and the family can split its runs of atoms in the same structure (Alg 7
/// l.12). The evaluator then never touches positions outside its range and splits the
/// structure at base and base + num_states() on its first cell (no-ops if the family already
/// isolated the atom). Keys must change only through this evaluator's offers.
///
/// Counters (Diagnostics::counters; into `diag` if non-null, else into the result's diag):
/// those of evaluate() (eval_resumed_slot_visits, eval_clipped_slot_visits, ...) plus
/// eval_portioned_calls, eval_cells (C), eval_external_offers (E), eval_external_decreases,
/// eval_root_slots (first portion of a root slot), eval_root_resumed_slot_visits,
/// eval_root_clipped_slot_visits, eval_root_empty_cells (cell with an unopened root and no due
/// key). No fallback: diag.fallbacks stays empty.
///
/// Not reentrant: calling offers, set_deadline or visit_cell from a hook throws
/// std::logic_error. After an exception from visit_cell the evaluator must be discarded. `nh`
/// (and a shared structure) must outlive the evaluator.
template <class W>
class PortionedEvaluator {
 public:
  /// Throws std::invalid_argument unless nh is normalised and ordered with exactly one root
  /// (an empty hierarchy is accepted and finished at once).
  explicit PortionedEvaluator(const NormalizedHierarchy<W>& nh, const EvalOptions<W>& options = {},
                              Diagnostics* diag = nullptr);
  /// The same with the layout at [base, base + num_states()) of `shared` (std::out_of_range if
  /// it does not fit).
  PortionedEvaluator(const NormalizedHierarchy<W>& nh, SplitFindmin<W>& shared, std::uint32_t base,
                     const EvalOptions<W>& options = {}, Diagnostics* diag = nullptr);
  PortionedEvaluator(PortionedEvaluator&&) noexcept;
  PortionedEvaluator& operator=(PortionedEvaluator&&) noexcept;
  PortionedEvaluator(const PortionedEvaluator&) = delete;
  PortionedEvaluator& operator=(const PortionedEvaluator&) = delete;
  ~PortionedEvaluator();

  [[nodiscard]] std::size_t num_states() const noexcept;
  /// Layout position of state x in the split-findmin structure (base included).
  [[nodiscard]] std::uint32_t position(StateId x) const;

  /// Alg 7 l.9: an initial offer, before the first cell (std::logic_error afterwards).
  /// std::out_of_range for a non-existent state, std::invalid_argument for a value that is not
  /// finite and nonnegative.
  void initial_offer(StateId x, W value);
  /// Lemma 7.7: an external offer between cells. Returns true iff the pending key decreased.
  bool external_offer(StateId x, W value);
  /// The atom's inclusive deadline (Sec 7.4 l.1825 "z + L_Y"); the smaller of all calls is
  /// kept. Only while the root has not opened (std::logic_error afterwards). Default +infinity.
  void set_deadline(W deadline);

  /// Visits the cell [lo, hi). Returns the emissions of this cell in emission order; the span
  /// is valid until the next call of a non-const member.
  std::span<const std::pair<StateId, W>> visit_cell(W lo, W hi);

  /// The root has opened (internal root) or emitted (singleton root).
  [[nodiscard]] bool opened() const noexcept;
  /// The root completed or covered its cap; later cells do nothing.
  [[nodiscard]] bool finished() const noexcept;
  /// Emitted value per state so far (infinity where none).
  [[nodiscard]] const std::vector<W>& value() const noexcept;
  /// All emissions so far, in emission order.
  [[nodiscard]] const std::vector<std::pair<StateId, W>>& emissions() const noexcept;

  /// Calls on_finish and moves the result out (counters included when no `diag` was given).
  /// The evaluator is unusable afterwards (std::logic_error).
  EvalResult<W> take_result();

 private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace sssp::paper
