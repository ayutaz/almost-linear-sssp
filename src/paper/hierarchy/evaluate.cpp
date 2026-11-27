// Evaluate (Alg 3) and Scan (Alg 2), Sec 4.3 of arXiv 2611.34911 v1, and the cell-by-cell
// visit of Lemma 7.7 (Sec 7.4, PortionedEvaluator). See evaluate.hpp for the overview and the
// interpretations (Q-27, Q-28, Q-29, Q-16).
//
// Implementation notes
//   * One state machine (Evaluator) serves both entry points. evaluate() runs Alg 3's root
//     loop (evaluate_root); PortionedEvaluator runs the cell procedure of Sec 7.4 l.1835-1841
//     (visit_root_cell). Below the root both use the same Scan frames, child slot advance /
//     resumption (start_child_slot, after_child_slot) and retirement, so a clipped slot is
//     resumed by exactly the code that the bulk run would use for an unclipped one.
//   * Scan recursion is driven by an explicit stack of frames (one frame per active
//     Scan(X, I, c) call) so that the depth of the hierarchy does not bound the call stack.
//   * Layout positions are absolute positions of the split-findmin structure: [0, n) for an
//     own structure, [base, base + n) of a shared one (Alg 7 l.8); layout_ is indexed by
//     position - base.
//   * Every list piece of every opened node is an interval of the leaf layout, and the
//     split-findmin intervals are exactly the current pieces (plus unopened registered children,
//     which keep their own interval until visited, l.940-941). A run's minimum is therefore
//     query(first position of the run): the "sample leaf" of l.945-946.
//   * Counters (Diagnostics::counters, compiled out with SSSP_DIAGNOSTICS=0):
//       eval_calls, eval_states, eval_arcs, eval_initial_offers (Q), eval_offered_states
//       eval_roots_skipped (Alg 3 l.5), eval_root_leaf_emissions (l.6-7)
//       eval_opens, eval_origin_second_rule, eval_origin_steps (successive subtractions),
//       eval_origin_rounding_fallbacks (double only, D-19; see open())
//       eval_scans (all Scan calls), eval_root_slots (Alg 3 l.10), eval_child_slot_visits
//       (Alg 2 l.10), eval_clipped_slot_visits (slot not inside the portion),
//       eval_resumed_slot_visits (portion starting after the slot start)
//       eval_run_visits, eval_failed_run_visits (run passed: no due minimum), eval_registrations,
//       eval_interval_splits (successful split-findmin splits incl. root boundaries)
//       eval_child_visits (registered pieces visited), eval_leaf_failed_checks (registered leaf
//       not due), eval_unopened_not_due (registered internal child without a due key)
//       eval_emissions, eval_relaxations, eval_relax_decreases, eval_rejected_offers (offer at a
//       settled or retired state), eval_supplemental_offers, eval_supplemental_decreases
//       eval_completed_nodes, eval_expired_nodes, eval_retire_piece_visits, eval_retired_leaves
//       eval_sfm_* (SplitFindminCounters of the evaluator's own structure)
//     PortionedEvaluator only: eval_portioned_calls, eval_cells (Lemma 7.7's C),
//       eval_external_offers (E), eval_external_decreases, eval_root_resumed_slot_visits,
//       eval_root_clipped_slot_visits, eval_root_empty_cells; eval_root_slots counts the first
//       portion of each root slot and eval_initial_offers each initial_offer() call.

#include "paper/hierarchy/evaluate.hpp"

#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/model/length.hpp"
#include "paper/primitives/split_findmin.hpp"

namespace sssp::paper {

namespace {

constexpr std::uint32_t kNil = std::numeric_limits<std::uint32_t>::max();

enum class Status : std::uint8_t { kUnopened, kLive, kDone };

/// The LCA / level-ancestor structure of a normalised hierarchy (Lemma 3.9).
template <class W>
const ForestQueries& queries_of(const NormalizedHierarchy<W>& nh) {
  if (!nh.queries.has_value()) {
    throw std::invalid_argument("evaluate: hierarchy is not normalised (Sec 4.1)");
  }
  return *nh.queries;
}

/// One element of an opened node's list (Sec 4.3 l.811-818): a registered child, or a maximal
/// run of unregistered children given by a rank range [r0, r1) of the owner's children.
struct Piece {
  bool is_run = false;
  NodeId child = kNoNode;  ///< registered child (is_run == false)
  std::uint32_t r0 = 0;    ///< run: first rank
  std::uint32_t r1 = 0;    ///< run: one past the last rank
  std::uint32_t prev = kNil;
  std::uint32_t next = kNil;
};

/// One active Scan(X, I = [lo, hi), nominal [a, b), c) (Alg 2), plus the visit of a registered
/// internal child in progress (Alg 2 l.9-10).
template <class W>
struct Frame {
  NodeId x;
  W lo, hi, a, b, c;
  std::uint32_t p;         ///< cursor: current piece of x's list, kNil at the end
  NodeId child = kNoNode;  ///< child whose slots are being visited (piece p)
  W child_c{};             ///< its deadline min(c, D0 + S(child))
};

template <class W>
class Evaluator {
 public:
  /// `shared` == nullptr: an own split-findmin structure over the layout [0, n). Otherwise the
  /// layout occupies [base, base + n) of `shared` (PortionedEvaluator, Alg 7 l.8).
  Evaluator(const NormalizedHierarchy<W>& nh, const EvalOptions<W>& opt, Diagnostics* diag,
            EvalResult<W>& out, SplitFindmin<W>* shared = nullptr, std::uint32_t base = 0)
      : nh_(nh),
        f_(nh.forest),
        qs_(queries_of(nh)),
        opt_(opt),
        hooks_(opt.hooks),
        diag_(diag),
        ops_(diag != nullptr ? &diag->ops : nullptr),
        out_(out),
        n_(nh.forest.num_states()),
        base_(base) {
    if (shared == nullptr) {
      sf_ = &own_sf_.emplace(n_ == 0 ? 1 : n_);
    } else {
      if (std::size_t{base} + n_ > shared->size()) {
        throw std::out_of_range("PortionedEvaluator: layout does not fit the shared structure");
      }
      sf_ = shared;
    }
  }

  void run(std::span<const StateOffer<W>> offers) {
    bump("eval_calls");
    bump("eval_states", n_);
    bump("eval_arcs", f_.graph.arcs.size());
    bump("eval_initial_offers", offers.size());
    out_.value.assign(n_, inf());
    // Lemma 4.7: an input scan keeps the minimum initial offer at each state.
    std::vector<W> best(n_, inf());
    for (const auto& o : offers) {
      check_offer(o.state, o.value);
      if (less(o.value, best[o.state])) best[o.state] = o.value;
    }
    if (n_ == 0) {
      if (hooks_ != nullptr) hooks_->on_finish();
      return;
    }
    init_structures();

    // Alg 3 l.3: split the layout at every root boundary, post the least initial offers.
    for (const NodeId r : nh_.root_order) split(lo_[r]);
    for (StateId x = 0; x < n_; ++x) {
      if (best[x] == inf()) continue;
      bump("eval_offered_states");
      if (hooks_ != nullptr) hooks_->on_initial_offer(x, best[x]);
      sf_->decrease(pos_[x], best[x]);
    }

    // Alg 3 l.4: roots in topological order of the root quotient.
    for (const NodeId r : nh_.root_order) {
      if (f_.is_leaf(r)) {
        const StateId x = f_.nodes[r].state;
        const W key = sf_->key(pos_[x]);
        if (key == inf()) {  // l.5
          bump("eval_roots_skipped");
          continue;
        }
        bump("eval_root_leaf_emissions");
        emit(x, key);  // l.6-7
        status_[r] = Status::kDone;
        continue;
      }
      const SfmMin<W> m = sf_->query(lo_[r]);
      if (m.key == inf()) {  // l.5
        bump("eval_roots_skipped");
        continue;
      }
      evaluate_root(r, m.key);
    }
    if (hooks_ != nullptr) hooks_->on_finish();
    SSSP_COUNT(if (diag_ != nullptr) sf_->counters().add_to(diag_->counters, "eval_sfm_"));
  }

  // ------------------------------------------------------------ Lemma 7.7 (PortionedEvaluator)

  /// Sets up an evaluation driven cell by cell (Sec 7.4 l.1835-1841). One root (Lemma 7.7).
  void portioned_init() {
    bump("eval_portioned_calls");
    bump("eval_states", n_);
    bump("eval_arcs", f_.graph.arcs.size());
    out_.value.assign(n_, inf());
    if (n_ == 0) {
      finished_ = true;
      return;
    }
    if (nh_.root_order.size() != 1) {
      throw std::invalid_argument("PortionedEvaluator: Lemma 7.7 needs exactly one root");
    }
    root_ = nh_.root_order.front();
    init_structures();
  }

  [[nodiscard]] std::size_t num_states() const noexcept { return n_; }
  [[nodiscard]] std::uint32_t position(StateId x) const {
    if (x >= n_) throw std::out_of_range("PortionedEvaluator: non-existent state");
    return pos_[x];
  }
  [[nodiscard]] bool root_opened() const noexcept { return root_opened_; }
  [[nodiscard]] bool finished() const noexcept { return finished_; }

  /// Alg 7 l.9 (Lemma 7.7's Q initial offers), before the first cell.
  void initial_offer(StateId x, W v) {
    not_in_visit();
    if (visited_) {
      throw std::logic_error("PortionedEvaluator: initial offers precede the first cell");
    }
    check_offer(x, v);
    bump("eval_initial_offers");
    if (hooks_ != nullptr) hooks_->on_initial_offer(x, v);
    const std::uint32_t q = pos_[x];
    if (less(v, sf_->key(q))) sf_->decrease(q, v);
  }

  /// Lemma 7.7's E external offers, under its arrival hypothesis (evaluate.hpp).
  bool external_offer(StateId x, W v) {
    not_in_visit();
    check_offer(x, v);
    bump("eval_external_offers");
    if (hooks_ != nullptr) hooks_->on_external_offer(x, v);
    if (settled_[x] != 0 || retired_[x] != 0) {  // Lemma 7.8: "ignores later proposals"
      bump("eval_rejected_offers");
      return false;
    }
    const std::uint32_t q = pos_[x];
    if (!less(v, sf_->key(q))) return false;
    if (visited_ && less(v, last_hi_)) {
      throw std::invalid_argument(
          "PortionedEvaluator: external offer below a visited cell (Lemma 7.7 arrival "
          "hypothesis)");
    }
    sf_->decrease(q, v);
    bump("eval_external_decreases");
    return true;
  }

  /// Sec 7.4 l.1825, l.1837-1838: the atom's lasting inclusive deadline.
  void set_deadline(W c) {
    not_in_visit();
    if (!T::is_valid_weight(c) && c != inf()) {
      throw std::invalid_argument("PortionedEvaluator: deadline must be nonnegative");
    }
    if (root_opened_) {
      throw std::logic_error("PortionedEvaluator: the deadline is fixed once the root opened");
    }
    deadline_ = min(deadline_, c);
  }

  /// One external cell [lo, hi); returns the index of its first emission in out_.emissions.
  std::size_t visit_cell(W lo, W hi) {
    not_in_visit();
    if (!T::is_valid_weight(lo) || (!T::is_valid_weight(hi) && hi != inf()) || !less(lo, hi)) {
      throw std::invalid_argument("PortionedEvaluator: a cell is [lo, hi) with finite lo < hi");
    }
    if (visited_ && less(lo, last_hi_)) {
      throw std::invalid_argument("PortionedEvaluator: cells overlap or go back in time");
    }
    if (visited_ && root_opened_ && !finished_ && less(last_hi_, lo)) {
      throw std::invalid_argument(
          "PortionedEvaluator: an open root needs consecutive cells (Lemma 7.7)");
    }
    const std::size_t first = out_.emissions.size();
    bump("eval_cells");
    if (!finished_) {
      in_visit_ = true;
      const VisitGuard guard{in_visit_};
      if (!visited_) {
        // Isolate the atom's interval (Alg 7 l.12, "splitting the family's live runs"): no-ops
        // for an own structure or for an atom the family already split off.
        split(base_);
        split(static_cast<std::uint32_t>(base_ + n_));
      }
      visit_root_cell(lo, hi);
    }
    visited_ = true;
    last_hi_ = hi;
    return first;
  }

  void finish_portioned() {
    if (hooks_ != nullptr) hooks_->on_finish();
    SSSP_COUNT(if (diag_ != nullptr && own_sf_.has_value()) sf_->counters().add_to(diag_->counters,
                                                                                   "eval_sfm_"));
  }

 private:
  /// Clears the reentrancy flag when a visit ends, also by an exception.
  struct VisitGuard {
    bool& flag;
    explicit VisitGuard(bool& f) : flag(f) {}
    VisitGuard(const VisitGuard&) = delete;
    VisitGuard& operator=(const VisitGuard&) = delete;
    ~VisitGuard() { flag = false; }
  };

  using T = WeightTraits<W>;
  using L = model::Len<W>;

  static constexpr W inf() noexcept { return T::infinity(); }

  // ------------------------------------------------------------ lengths (Sec 2 model)

  [[nodiscard]] W add(W a, W b) const { return ops_.add(L::finite(a), L::finite(b)).value(); }
  [[nodiscard]] bool less(W a, W b) const {
    SSSP_COUNT(if (ops_.counters() != nullptr)++ ops_.counters()->comparisons);
    return a < b;
  }
  [[nodiscard]] W min(W a, W b) const { return less(b, a) ? b : a; }
  [[nodiscard]] W max(W a, W b) const { return less(a, b) ? b : a; }
  /// One subtraction of the origin rule ([Deviation D-19], evaluate.hpp).
  [[nodiscard]] W sub(W a, W b) const {
    SSSP_COUNT(if (ops_.counters() != nullptr)++ ops_.counters()->additions);
    return a - b;
  }
  /// Due in the window of frame F: a' <= z < b' and z <= c (Sec 4.3 l.786-789).
  [[nodiscard]] bool due(W z, const Frame<W>& F) const {
    return z != inf() && !less(z, F.lo) && less(z, F.hi) && !less(F.c, z);
  }

  void bump([[maybe_unused]] const char* name, [[maybe_unused]] std::uint64_t k = 1) {
    SSSP_COUNT(if (diag_ != nullptr) diag_->counters[name] += k);
  }

  void check_offer(StateId x, W v) const {
    if (x >= n_) throw std::out_of_range("evaluate: initial offer at a non-existent state");
    if (!T::is_valid_weight(v) || v == inf()) {
      throw std::invalid_argument("evaluate: initial offers must be finite and nonnegative");
    }
  }

  /// Layout, outgoing arcs and the per-state and per-node arrays (Alg 3 l.2-3).
  void init_structures() {
    build_layout();
    build_out_arcs();
    settled_.assign(n_, 0);
    retired_.assign(n_, 0);
    const std::size_t nn = nh_.num_nodes();
    status_.assign(nn, Status::kUnopened);
    slot_lo_.assign(nn, W{});
    slot_hi_.assign(nn, W{});
    expiry_.assign(nn, W{});
    head_.assign(nn, kNil);
  }

  void not_in_visit() const {
    if (in_visit_) throw std::logic_error("PortionedEvaluator: not reentrant (called from a hook)");
  }

  // ------------------------------------------------------------ layout (Alg 3 l.2)

  /// DFS order of the leaves with children in rank order and roots in root_order; every node
  /// occupies [lo_, hi_) (Sec 4.3 l.809, "hierarchy-leaf DFS layout respecting these
  /// topological child orders").
  void build_layout() {
    const std::size_t nn = nh_.num_nodes();
    lo_.assign(nn, 0);
    hi_.assign(nn, 0);
    pos_.assign(n_, 0);
    layout_.assign(n_, kNoState);
    std::uint32_t next = base_;
    std::vector<std::pair<NodeId, std::size_t>> stack;
    for (const NodeId r : nh_.root_order) {
      stack.emplace_back(r, 0);
      lo_[r] = next;
      while (!stack.empty()) {
        auto& [v, k] = stack.back();
        const auto& nd = f_.nodes[v];
        if (nd.children.empty()) {
          pos_[nd.state] = next;
          layout_[next - base_] = nd.state;
          ++next;
          hi_[v] = next;
          stack.pop_back();
          continue;
        }
        if (k < nd.children.size()) {
          const NodeId c = nd.children[k++];
          lo_[c] = next;
          stack.emplace_back(c, 0);
        } else {
          hi_[v] = next;
          stack.pop_back();
        }
      }
    }
    if (next != base_ + n_) {
      throw std::invalid_argument("evaluate: forest leaves are not the states");
    }
  }

  /// Outgoing arcs per state in arc order; state self-arcs are discarded (Sec 4.1 l.719).
  void build_out_arcs() {
    const auto& arcs = f_.graph.arcs;
    out_begin_.assign(n_ + 1, 0);
    for (const auto& a : arcs) {
      if (a.tail != a.head) ++out_begin_[a.tail + 1];
    }
    for (std::size_t s = 0; s < n_; ++s) out_begin_[s + 1] += out_begin_[s];
    out_arc_.assign(out_begin_[n_], 0);
    std::vector<std::size_t> at(out_begin_.begin(), out_begin_.end() - 1);
    for (ArcId i = 0; i < arcs.size(); ++i) {
      if (arcs[i].tail != arcs[i].head) out_arc_[at[arcs[i].tail]++] = i;
    }
  }

  [[nodiscard]] std::span<const StateId> states_of(NodeId v) const {
    return {layout_.data() + (lo_[v] - base_), layout_.data() + (hi_[v] - base_)};
  }

  void split(std::uint32_t pos) {
    if (sf_->split(pos)) bump("eval_interval_splits");
  }

  // ------------------------------------------------------------ lists

  std::uint32_t new_run(std::uint32_t r0, std::uint32_t r1) {
    Piece p;
    p.is_run = true;
    p.r0 = r0;
    p.r1 = r1;
    pieces_.push_back(p);
    return static_cast<std::uint32_t>(pieces_.size() - 1);
  }
  std::uint32_t new_child(NodeId c) {
    Piece p;
    p.child = c;
    pieces_.push_back(p);
    return static_cast<std::uint32_t>(pieces_.size() - 1);
  }
  /// Replaces piece p of x's list by the chain `seq` (in order), in its list position.
  void replace(NodeId x, std::uint32_t p, std::span<const std::uint32_t> seq) {
    const std::uint32_t before = pieces_[p].prev;
    const std::uint32_t after = pieces_[p].next;
    std::uint32_t last = before;
    for (const std::uint32_t q : seq) {
      pieces_[q].prev = last;
      if (last == kNil) {
        head_[x] = q;
      } else {
        pieces_[last].next = q;
      }
      last = q;
    }
    pieces_[last].next = after;
    if (after != kNil) pieces_[after].prev = last;
    pieces_[p].prev = pieces_[p].next = kNil;
  }
  void unlink(NodeId x, std::uint32_t p) {
    const std::uint32_t a = pieces_[p].prev;
    const std::uint32_t b = pieces_[p].next;
    if (a == kNil) {
      head_[x] = b;
    } else {
      pieces_[a].next = b;
    }
    if (b != kNil) pieces_[b].prev = a;
    pieces_[p].prev = pieces_[p].next = kNil;
  }
  /// All leaves settled or retired: nothing is left in the list (Sec 4.3 l.792).
  [[nodiscard]] bool completed(NodeId x) const { return head_[x] == kNil; }

  // ------------------------------------------------------------ offers and emissions

  /// Posts value v at state y unless y is settled or retired (Sec 4.3 l.859-860); returns true
  /// if the key decreased.
  bool post(StateId y, W v) {
    if (settled_[y] != 0 || retired_[y] != 0) {
      bump("eval_rejected_offers");
      return false;
    }
    const std::uint32_t q = pos_[y];
    if (!less(v, sf_->key(q))) return false;
    sf_->decrease(q, v);
    return true;
  }

  /// Emits state x at value v once and relaxes its outgoing arcs (Alg 2 l.11, Alg 3 l.7).
  void emit(StateId x, W v) {
    settled_[x] = 1;
    out_.value[x] = v;
    out_.emissions.emplace_back(x, v);
    bump("eval_emissions");
    if (hooks_ != nullptr) hooks_->on_emit(x, v);
    const auto& arcs = f_.graph.arcs;
    for (std::size_t k = out_begin_[x]; k < out_begin_[x + 1]; ++k) {
      const ArcId a = out_arc_[k];
      const W nv = add(v, arcs[a].length);
      bump("eval_relaxations");
      if (hooks_ != nullptr) hooks_->on_relax_offer(a, nv);
      if (post(arcs[a].head, nv)) bump("eval_relax_decreases");
    }
  }

  /// Alg 2 l.6: supplemental offer at state x (Sec 4.3 l.824-832).
  void supplemental(StateId x, W u) {
    bump("eval_supplemental_offers");
    if (hooks_ != nullptr) hooks_->on_supplemental_offer(x, u);
    if (post(x, u)) bump("eval_supplemental_decreases");
  }

  // ------------------------------------------------------------ opening, retirement

  /// Opens internal node x at its least key D0 (Sec 4.3 l.790-797). `root`: origin D0;
  /// otherwise the two-case rule with [a, b) the parent's NOMINAL slot (Q-28): D0 if
  /// D0 + S(x) < b, else the largest b - j nu(x) <= D0, by successive subtractions.
  void open(NodeId x, W d0, bool root, W a, W b) {
    const W nu = nh_.unit[x];
    const W s = nh_.mass[x];
    const W expiry = add(d0, s);
    W origin = d0;
    if (!root && !less(expiry, b)) {
      bump("eval_origin_second_rule");
      // Sec 4.3 l.798-800: the origin lies inside the parent's slot, so no subtraction may go
      // below a (t - nu >= a, i.e. t >= a + nu).
      const W floor = add(a, nu);
      W t = b;
      bool left_slot = false;
      while (less(d0, t)) {
        if (less(t, floor)) {
          left_slot = true;
          break;
        }
        t = sub(t, nu);
        bump("eval_origin_steps");
      }
      if (!left_slot) {
        origin = t;
      } else if constexpr (std::is_floating_point_v<W>) {
        // [Deviation D-19] Exact arithmetic cannot get here (l.798-800);
        // rounded subtractions can (an exact tie D0 = b - j nu computed as D0 + 1 ulp). Start at
        // D0 as in the first rule: the grid is then not aligned to b, which costs at most a
        // clipped slot and is safe under the slot-resumption rule (Q-27).
        bump("eval_origin_rounding_fallbacks");
        origin = d0;
      } else {
        throw std::logic_error("evaluate: origin rule left the parent's slot (Def 4.1 violated)");
      }
    }
    status_[x] = Status::kLive;
    slot_lo_[x] = origin;
    slot_hi_[x] = add(origin, nu);
    expiry_[x] = expiry;
    // "On first opening X, its whole child list is one run" (l.818-819).
    const auto d = static_cast<std::uint32_t>(f_.nodes[x].children.size());
    head_[x] = new_run(0, d);
    bump("eval_opens");
    if (hooks_ != nullptr) hooks_->on_open(x, states_of(x), d0, origin, nu, s);
  }

  /// Retires x after its deadline was covered (Alg 2 l.12, Sec 4.3 l.856-860): walks only the
  /// pieces still linked below x, recursing into registered children and marking the leaves
  /// of unregistered runs and unopened registered children.
  void retire(NodeId x) {
    bump("eval_expired_nodes");
    if (hooks_ != nullptr) hooks_->on_retire(x, states_of(x));
    std::vector<NodeId> stack{x};
    while (!stack.empty()) {
      const NodeId v = stack.back();
      stack.pop_back();
      if (!f_.is_leaf(v) && status_[v] == Status::kLive) {
        for (std::uint32_t p = head_[v]; p != kNil; p = pieces_[p].next) {
          bump("eval_retire_piece_visits");
          const Piece& pc = pieces_[p];
          if (pc.is_run) {
            const auto& ch = f_.nodes[v].children;
            mark_retired(lo_[ch[pc.r0]], hi_[ch[pc.r1 - 1]]);
          } else {
            stack.push_back(pc.child);
          }
        }
        head_[v] = kNil;
      } else if (status_[v] != Status::kDone) {
        // A registered leaf, or a registered child that never opened: its whole interval is
        // untouched (no state below it has emitted).
        mark_retired(lo_[v], hi_[v]);
      }
      status_[v] = Status::kDone;
    }
  }
  void mark_retired(std::uint32_t lo, std::uint32_t hi) {
    for (std::uint32_t q = lo; q < hi; ++q) {
      const StateId y = layout_[q - base_];
      if (settled_[y] == 0 && retired_[y] == 0) {
        retired_[y] = 1;
        bump("eval_retired_leaves");
      }
    }
  }

  // ------------------------------------------------------------ Alg 3 l.8-11

  void evaluate_root(NodeId r, W d0) {
    open(r, d0, /*root=*/true, W{}, W{});
    const W c = expiry_[r];  // l.9: c = D0 + S(X)
    // l.10: successive width-nu slots beginning at D0 through the slot containing c.
    for (W s = d0; !less(c, s);) {
      const W e = add(s, nh_.unit[r]);
      bump("eval_root_slots");
      slot_lo_[r] = s;
      slot_hi_[r] = e;
      push_scan(r, s, e, s, e, c);
      drive();
      if (completed(r)) break;  // l.11
      s = e;
    }
    if (completed(r)) {
      bump("eval_completed_nodes");
      status_[r] = Status::kDone;
    } else {
      retire(r);
    }
  }

  // ------------------------------------------------------------ Sec 7.4, Lemma 7.7

  /// Due in the cell [lo, hi) under the inclusive cap c.
  [[nodiscard]] bool due_in(W z, W lo, W hi, W c) const {
    return z != inf() && !less(z, lo) && less(z, hi) && !less(c, z);
  }

  /// Lemma 7.7: "The first cell starts at or before every finite initial timestamp", and
  /// external offers are "never inserted below the current cell's left endpoint".
  void check_not_below(W z, W lo) const {
    if (z != inf() && less(z, lo) && !less(deadline_, z)) {
      throw std::logic_error(
          "PortionedEvaluator: a pending key lies below the cell (Lemma 7.7 hypotheses)");
    }
  }

  /// "Visiting an active atom" (Sec 7.4 l.1835-1841) in the cell [lo, hi).
  void visit_root_cell(W lo, W hi) {
    const NodeId r = root_;
    // The cap is covered once the cell ends after it ("including equality"); after a last
    // cell [lo, +inf) nothing is left for later.
    const auto covered = [&](W c) { return less(c, hi) || hi == inf(); };
    if (f_.is_leaf(r)) {
      // l.1835: "At a singleton root, emit the due minimum once ... and retire the leaf".
      const StateId x = f_.nodes[r].state;
      const W key = sf_->key(pos_[x]);
      check_not_below(key, lo);
      if (due_in(key, lo, hi, deadline_)) {
        bump("eval_root_leaf_emissions");
        root_opened_ = true;
        status_[r] = Status::kDone;
        emit(x, key);
        finished_ = true;
      } else if (covered(deadline_)) {
        retire(r);
        finished_ = true;
      }
      return;
    }
    if (status_[r] == Status::kUnopened) {
      // l.1836-1837: "open on the first due key D0 and anchor its clock at D0"; Lemma 7.7: "an
      // empty cell leaves an unopened root unopened".
      const SfmMin<W> m = sf_->query(lo_[r]);
      check_not_below(m.key, lo);
      if (!due_in(m.key, lo, hi, deadline_)) {
        bump("eval_root_empty_cells");
        if (covered(deadline_)) {
          retire(r);
          finished_ = true;
        }
        return;
      }
      open(r, m.key, /*root=*/true, W{}, W{});
      root_opened_ = true;
      // l.1837-1838: "the smaller of the atom deadline and the root expiry D0 + S(root)".
      cap_ = min(deadline_, expiry_[r]);
    }
    // l.1838-1839: the successive nominal root slots whose intersections with the cell contain
    // times at or below c, each by Scan(root, J n [lo, hi), c) with the nominal slot J.
    const W nu = nh_.unit[r];
    for (;;) {
      const W jl = slot_lo_[r];
      const W jh = slot_hi_[r];
      if (!less(lo, jh)) throw std::logic_error("PortionedEvaluator: root slot behind the cell");
      const W plo = max(jl, lo);
      if (!less(jl, hi) || less(cap_, plo)) break;
      if (less(jl, lo)) {
        bump("eval_root_resumed_slot_visits");  // l.1840: resumed at the first live piece
      } else {
        bump("eval_root_slots");
      }
      if (less(hi, jh)) bump("eval_root_clipped_slot_visits");
      push_scan(r, plo, min(jh, hi), jl, jh, cap_);
      drive();
      // [Interpretation, Q-27] a slot clipped at hi stays current for the next cell.
      if (completed(r) || less(hi, jh)) break;
      slot_lo_[r] = jh;
      slot_hi_[r] = add(jh, nu);
    }
    // l.1841: "Stop once all children finish or the cap has been covered, including equality".
    if (completed(r)) {
      bump("eval_completed_nodes");
      status_[r] = Status::kDone;
      finished_ = true;
    } else if (covered(cap_)) {
      retire(r);
      finished_ = true;
    }
  }

  // ------------------------------------------------------------ Alg 2

  void push_scan(NodeId x, W lo, W hi, W a, W b, W c) {
    bump("eval_scans");
    if (hooks_ != nullptr) hooks_->on_scan_begin(x, states_of(x), lo, hi, a, b, c);
    frames_.push_back(Frame<W>{x, lo, hi, a, b, c, head_[x]});
  }

  /// Runs the innermost Scan calls until the stack is empty.
  void drive() {
    while (!frames_.empty()) {
      const std::size_t fi = frames_.size() - 1;
      if (frames_[fi].child != kNoNode) {
        after_child_slot(fi);
      } else if (frames_[fi].p == kNil) {
        // End of the list. Alg 2 l.12: completed children and children whose deadline is
        // covered were retired at their visit; runs and registered leaves carry only x's own
        // deadline, which x's parent (or Alg 3) handles.
        const NodeId x = frames_[fi].x;
        frames_.pop_back();
        if (hooks_ != nullptr) hooks_->on_scan_end(x, states_of(x));
      } else if (pieces_[frames_[fi].p].is_run) {
        visit_run(fi);
      } else {
        bump("eval_child_visits");
        const NodeId i = pieces_[frames_[fi].p].child;
        if (f_.is_leaf(i)) {
          visit_leaf(fi, i);
        } else {
          visit_internal(fi, i);
        }
      }
    }
  }

  /// Alg 2 l.3-7 on the run at the cursor.
  void visit_run(std::size_t fi) {
    bump("eval_run_visits");
    const Frame<W> F = frames_[fi];
    const std::uint32_t p = F.p;
    const std::uint32_t r0 = pieces_[p].r0;
    const std::uint32_t r1 = pieces_[p].r1;
    const auto& ch = f_.nodes[F.x].children;
    const SfmMin<W> m = sf_->query(lo_[ch[r0]]);  // l.4, Lemma 3.11
    if (!due(m.key, F)) {                         // passed for this portion (l.821)
      bump("eval_failed_run_visits");
      frames_[fi].p = pieces_[p].next;
      return;
    }
    // l.5: the child of X containing the minimising leaf, by a level-ancestor query (l.817).
    const NodeId leaf = f_.leaf_of_state[layout_[m.index - base_]];
    const NodeId i = qs_.level_ancestor(leaf, qs_.depth(F.x) + 1);
    const std::uint32_t r = nh_.rank[i];
    if (r < r0 || r >= r1 || ch[r] != i)
      throw std::logic_error("evaluate: run minimum outside run");
    split(lo_[i]);  // at most two interval splits isolate the states below i
    split(hi_[i]);
    bump("eval_registrations");
    std::uint32_t seq[3];
    std::size_t k = 0;
    const std::uint32_t left = r0 < r ? new_run(r0, r) : kNil;
    const std::uint32_t mid = new_child(i);
    const std::uint32_t right = r + 1 < r1 ? new_run(r + 1, r1) : kNil;
    if (left != kNil) seq[k++] = left;
    seq[k++] = mid;
    if (right != kNil) seq[k++] = right;
    replace(F.x, p, std::span<const std::uint32_t>(seq, k));
    // l.6: u = max{b, z + S(i) + g(i, j)} at x_j for the child j of each new subrun next to i,
    // padded with the NOMINAL slot end b (Q-28, Lemma 7.7 l.1890).
    if (opt_.supplemental_offers && (left != kNil || right != kNil)) {
      const W zs = add(m.key, nh_.mass[i]);
      if (left != kNil) {
        const NodeId j = ch[r - 1];
        supplemental(nh_.rep[j], max(F.b, add(zs, nh_.gap_next[j])));  // g(j, i) = g(i, j)
      }
      if (right != kNil) {
        const NodeId j = ch[r + 1];
        supplemental(nh_.rep[j], max(F.b, add(zs, nh_.gap_next[i])));
      }
    }
    // l.7: resume at the first piece just created.
    frames_[fi].p = left != kNil ? left : mid;
  }

  /// Alg 2 l.11 at a registered leaf: emit its due value once.
  void visit_leaf(std::size_t fi, NodeId i) {
    const std::uint32_t p = frames_[fi].p;
    const std::uint32_t next = pieces_[p].next;
    const StateId x = f_.nodes[i].state;
    const W key = sf_->key(pos_[x]);
    if (settled_[x] == 0 && retired_[x] == 0 && due(key, frames_[fi])) {
      unlink(frames_[fi].x, p);
      status_[i] = Status::kDone;
      emit(x, key);
    } else {
      bump("eval_leaf_failed_checks");
    }
    frames_[fi].p = next;
  }

  /// Alg 2 l.9-10 at a registered internal child: open it if necessary with its least due offer,
  /// then visit its live slots that meet the portion, through its inclusive deadline.
  void visit_internal(std::size_t fi, NodeId i) {
    const Frame<W> F = frames_[fi];
    if (status_[i] == Status::kUnopened) {
      const SfmMin<W> m = sf_->query(lo_[i]);
      if (!due(m.key, F)) {
        bump("eval_unopened_not_due");
        frames_[fi].p = pieces_[F.p].next;
        return;
      }
      open(i, m.key, /*root=*/false, F.a, F.b);
    }
    frames_[fi].child = i;
    frames_[fi].child_c = min(F.c, expiry_[i]);
    if (!start_child_slot(fi)) finish_child(fi);
  }

  /// [Interpretation, Q-27] Starts Scan(i, J n I, c_i) for the child's current slot J if J meets
  /// I and starts no later than c_i; returns false if there is no such slot.
  bool start_child_slot(std::size_t fi) {
    const Frame<W> F = frames_[fi];
    const NodeId i = F.child;
    const W jl = slot_lo_[i];
    const W jh = slot_hi_[i];
    // Slots of a live child are contiguous across its parent's portions (they advance only up to
    // the portion end), so the current slot always reaches into the portion.
    if (!less(F.lo, jh)) throw std::logic_error("evaluate: child slot behind the portion");
    if (!less(jl, F.hi) || less(F.child_c, jl)) return false;
    bump("eval_child_slot_visits");
    if (less(jl, F.lo)) bump("eval_resumed_slot_visits");
    if (less(F.hi, jh)) bump("eval_clipped_slot_visits");
    push_scan(i, max(jl, F.lo), min(jh, F.hi), jl, jh, F.child_c);
    return true;
  }

  /// Back in frame fi after Scan(i, J n I, c_i) returned.
  void after_child_slot(std::size_t fi) {
    const NodeId i = frames_[fi].child;
    if (!completed(i) && !less(frames_[fi].hi, slot_hi_[i])) {
      // J ended inside I: advance to the next slot (the clock persists, l.861).
      slot_lo_[i] = slot_hi_[i];
      slot_hi_[i] = add(slot_lo_[i], nh_.unit[i]);
      if (start_child_slot(fi)) return;
    }
    // Completed, or J is clipped at b' (it stays current for the next portion), or no further
    // slot meets I.
    finish_child(fi);
  }

  /// Alg 2 l.12 for the child at the cursor: retire it if it completed or covered its deadline
  /// (c_i < b'), and move the cursor on.
  void finish_child(std::size_t fi) {
    const Frame<W> F = frames_[fi];
    const NodeId i = F.child;
    const std::uint32_t next = pieces_[F.p].next;
    if (completed(i)) {
      bump("eval_completed_nodes");
      status_[i] = Status::kDone;
      if (hooks_ != nullptr) hooks_->on_retire(i, states_of(i));
      unlink(F.x, F.p);
    } else if (less(F.child_c, F.hi)) {
      retire(i);
      unlink(F.x, F.p);
    }
    frames_[fi].child = kNoNode;
    frames_[fi].p = next;
  }

  const NormalizedHierarchy<W>& nh_;
  const OrderedHierarchy<W>& f_;
  const ForestQueries& qs_;
  const EvalOptions<W>& opt_;
  EvalHooks<W>* hooks_;
  Diagnostics* diag_;
  model::LenOps<W> ops_;
  EvalResult<W>& out_;
  std::size_t n_;
  std::uint32_t base_;  ///< first layout position in *sf_
  std::optional<SplitFindmin<W>> own_sf_;
  SplitFindmin<W>* sf_ = nullptr;

  std::vector<std::uint32_t> lo_, hi_, pos_;
  std::vector<StateId> layout_;
  std::vector<std::size_t> out_begin_;
  std::vector<ArcId> out_arc_;
  std::vector<std::uint8_t> settled_, retired_;
  std::vector<Status> status_;
  std::vector<W> slot_lo_, slot_hi_, expiry_;
  std::vector<std::uint32_t> head_;
  std::vector<Piece> pieces_;
  std::vector<Frame<W>> frames_;

  // PortionedEvaluator state (Lemma 7.7).
  NodeId root_ = kNoNode;
  W deadline_ = inf();  ///< the atom's inclusive deadline (Sec 7.4 l.1825)
  W cap_ = inf();       ///< min(deadline, D0 + S(root)) once the root opened (l.1837-1838)
  W last_hi_{};         ///< end of the last visited cell
  bool visited_ = false;
  bool in_visit_ = false;
  bool root_opened_ = false;
  bool finished_ = false;
};

}  // namespace

template <class W>
EvalResult<W> evaluate_normalized(const NormalizedHierarchy<W>& nh,
                                  std::span<const StateOffer<W>> offers,
                                  const EvalOptions<W>& options, Diagnostics* diag) {
  if (!nh.children_ordered || !nh.queries.has_value()) {
    throw std::invalid_argument("evaluate: hierarchy is not normalised and ordered (Sec 4.1-4.2)");
  }
  EvalResult<W> out;
  Evaluator<W> ev(nh, options, diag, out);
  ev.run(offers);
  return out;
}

template <class W>
EvalResult<W> evaluate(const SeededOrderedHierarchy<W>& s, const EvalOptions<W>& options) {
  const auto& h = s.hierarchy;
  if (h.num_states() == 0) {
    if (!s.offers.empty()) {
      throw std::out_of_range("evaluate: initial offer at a non-existent state");
    }
    EvalResult<W> out;
    if (options.hooks != nullptr) options.hooks->on_finish();
    return out;
  }
  Diagnostics diag;
  // Alg 3 l.2: normalisation (Lemma 4.3) and child orders and gaps (Lemma 4.5).
  NormalizedHierarchy<W> nh = normalize<W>(h, &diag);
  order_children<W>(nh, &diag);
  EvalResult<W> out = evaluate_normalized<W>(nh, s.offers, options, &diag);
  out.diag = std::move(diag);
  return out;
}

template EvalResult<std::uint64_t> evaluate<std::uint64_t>(
    const SeededOrderedHierarchy<std::uint64_t>&, const EvalOptions<std::uint64_t>&);
template EvalResult<double> evaluate<double>(const SeededOrderedHierarchy<double>&,
                                             const EvalOptions<double>&);
template EvalResult<std::uint64_t> evaluate_normalized<std::uint64_t>(
    const NormalizedHierarchy<std::uint64_t>&, std::span<const StateOffer<std::uint64_t>>,
    const EvalOptions<std::uint64_t>&, Diagnostics*);
template EvalResult<double> evaluate_normalized<double>(const NormalizedHierarchy<double>&,
                                                        std::span<const StateOffer<double>>,
                                                        const EvalOptions<double>&, Diagnostics*);

// ------------------------------------------------------------------ PortionedEvaluator

template <class W>
class PortionedEvaluator<W>::Impl {
 public:
  Impl(const NormalizedHierarchy<W>& nh, SplitFindmin<W>* shared, std::uint32_t base,
       const EvalOptions<W>& options, Diagnostics* diag)
      : options_(options),
        ev_(nh, options_, diag != nullptr ? diag : &out_.diag, out_, shared, base) {
    ev_.portioned_init();
  }

  Evaluator<W>& ev() {
    if (taken_) throw std::logic_error("PortionedEvaluator: result already taken");
    return ev_;
  }
  [[nodiscard]] const Evaluator<W>& state() const noexcept { return ev_; }
  [[nodiscard]] const EvalResult<W>& out() const noexcept { return out_; }

  EvalResult<W> take() {
    ev().finish_portioned();
    taken_ = true;
    return std::move(out_);
  }

 private:
  EvalOptions<W> options_;
  EvalResult<W> out_;
  Evaluator<W> ev_;  // after options_ and out_, which it references
  bool taken_ = false;
};

namespace {

/// The hierarchy, once checked to be normalised and ordered (before Evaluator reads it).
template <class W>
const NormalizedHierarchy<W>* checked_for_portions(const NormalizedHierarchy<W>* nh) {
  if (!nh->children_ordered || !nh->queries.has_value()) {
    throw std::invalid_argument(
        "PortionedEvaluator: hierarchy is not normalised and ordered (Sec 4.1-4.2)");
  }
  return nh;
}

}  // namespace

template <class W>
PortionedEvaluator<W>::PortionedEvaluator(const NormalizedHierarchy<W>& nh,
                                          const EvalOptions<W>& options, Diagnostics* diag)
    : impl_(std::make_unique<Impl>(*checked_for_portions(&nh), nullptr, 0, options, diag)) {}

template <class W>
PortionedEvaluator<W>::PortionedEvaluator(const NormalizedHierarchy<W>& nh, SplitFindmin<W>& shared,
                                          std::uint32_t base, const EvalOptions<W>& options,
                                          Diagnostics* diag)
    : impl_(std::make_unique<Impl>(*checked_for_portions(&nh), &shared, base, options, diag)) {}

template <class W>
PortionedEvaluator<W>::PortionedEvaluator(PortionedEvaluator&&) noexcept = default;
template <class W>
PortionedEvaluator<W>& PortionedEvaluator<W>::operator=(PortionedEvaluator&&) noexcept = default;
template <class W>
PortionedEvaluator<W>::~PortionedEvaluator() = default;

template <class W>
std::size_t PortionedEvaluator<W>::num_states() const noexcept {
  return impl_->state().num_states();
}
template <class W>
std::uint32_t PortionedEvaluator<W>::position(StateId x) const {
  return impl_->state().position(x);
}
template <class W>
void PortionedEvaluator<W>::initial_offer(StateId x, W value) {
  impl_->ev().initial_offer(x, value);
}
template <class W>
bool PortionedEvaluator<W>::external_offer(StateId x, W value) {
  return impl_->ev().external_offer(x, value);
}
template <class W>
void PortionedEvaluator<W>::set_deadline(W deadline) {
  impl_->ev().set_deadline(deadline);
}
template <class W>
std::span<const std::pair<StateId, W>> PortionedEvaluator<W>::visit_cell(W lo, W hi) {
  const std::size_t first = impl_->ev().visit_cell(lo, hi);
  return std::span<const std::pair<StateId, W>>(impl_->out().emissions).subspan(first);
}
template <class W>
bool PortionedEvaluator<W>::opened() const noexcept {
  return impl_->state().root_opened();
}
template <class W>
bool PortionedEvaluator<W>::finished() const noexcept {
  return impl_->state().finished();
}
template <class W>
const std::vector<W>& PortionedEvaluator<W>::value() const noexcept {
  return impl_->out().value;
}
template <class W>
const std::vector<std::pair<StateId, W>>& PortionedEvaluator<W>::emissions() const noexcept {
  return impl_->out().emissions;
}
template <class W>
EvalResult<W> PortionedEvaluator<W>::take_result() {
  return impl_->take();
}

template class PortionedEvaluator<std::uint64_t>;
template class PortionedEvaluator<double>;

}  // namespace sssp::paper
