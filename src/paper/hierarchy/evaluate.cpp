// Evaluate (Alg 3) and Scan (Alg 2), Sec 4.3 of arXiv 2611.34911 v1. See evaluate.hpp for the
// overview and the interpretations (Q-27, Q-28, Q-29, Q-16).
//
// Implementation notes
//   * Scan recursion is driven by an explicit stack of frames (one frame per active
//     Scan(X, I, c) call) so that the depth of the hierarchy does not bound the call stack.
//   * Every list piece of every opened node is an interval of the leaf layout, and the
//     split-findmin intervals are exactly the current pieces (plus unopened registered children,
//     which keep their own interval until visited, l.940-941). A run's minimum is therefore
//     query(first position of the run): the "sample leaf" of l.945-946.
//   * Counters (Diagnostics::counters, compiled out with SSSP_DIAGNOSTICS=0):
//       eval_calls, eval_states, eval_arcs, eval_initial_offers (Q), eval_offered_states
//       eval_roots_skipped (Alg 3 l.5), eval_root_leaf_emissions (l.6-7)
//       eval_opens, eval_origin_second_rule, eval_origin_steps (successive subtractions),
//       eval_origin_rounding_fallbacks (double only, proposed D-16; see open())
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
//       eval_sfm_* (SplitFindminCounters of the evaluator's structure)

#include "paper/hierarchy/evaluate.hpp"

#include <cstddef>
#include <limits>
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
  Evaluator(const NormalizedHierarchy<W>& nh, const EvalOptions<W>& opt, Diagnostics* diag,
            EvalResult<W>& out)
      : nh_(nh),
        f_(nh.forest),
        qs_(*nh.queries),
        opt_(opt),
        hooks_(opt.hooks),
        diag_(diag),
        ops_(diag != nullptr ? &diag->ops : nullptr),
        out_(out),
        n_(nh.forest.num_states()),
        sf_(n_ == 0 ? 1 : n_) {}

  void run(std::span<const StateOffer<W>> offers) {
    bump("eval_calls");
    bump("eval_states", n_);
    bump("eval_arcs", f_.graph.arcs.size());
    bump("eval_initial_offers", offers.size());
    out_.value.assign(n_, inf());
    // Lemma 4.7: an input scan keeps the minimum initial offer at each state.
    std::vector<W> best(n_, inf());
    for (const auto& o : offers) {
      if (o.state >= n_) throw std::out_of_range("evaluate: initial offer at a non-existent state");
      if (!WeightTraits<W>::is_valid_weight(o.value) || o.value == inf()) {
        throw std::invalid_argument("evaluate: initial offers must be finite and nonnegative");
      }
      if (less(o.value, best[o.state])) best[o.state] = o.value;
    }
    if (n_ == 0) {
      if (hooks_ != nullptr) hooks_->on_finish();
      return;
    }
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

    // Alg 3 l.3: split the layout at every root boundary, post the least initial offers.
    for (const NodeId r : nh_.root_order) split(lo_[r]);
    for (StateId x = 0; x < n_; ++x) {
      if (best[x] == inf()) continue;
      bump("eval_offered_states");
      if (hooks_ != nullptr) hooks_->on_initial_offer(x, best[x]);
      sf_.decrease(pos_[x], best[x]);
    }

    // Alg 3 l.4: roots in topological order of the root quotient.
    for (const NodeId r : nh_.root_order) {
      if (f_.is_leaf(r)) {
        const StateId x = f_.nodes[r].state;
        const W key = sf_.key(pos_[x]);
        if (key == inf()) {  // l.5
          bump("eval_roots_skipped");
          continue;
        }
        bump("eval_root_leaf_emissions");
        emit(x, key);  // l.6-7
        status_[r] = Status::kDone;
        continue;
      }
      const SfmMin<W> m = sf_.query(lo_[r]);
      if (m.key == inf()) {  // l.5
        bump("eval_roots_skipped");
        continue;
      }
      evaluate_root(r, m.key);
    }
    if (hooks_ != nullptr) hooks_->on_finish();
    SSSP_COUNT(if (diag_ != nullptr) sf_.counters().add_to(diag_->counters, "eval_sfm_"));
  }

 private:
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
  /// One subtraction of the origin rule ([Deviation candidate, proposed D-16], evaluate.hpp).
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
    std::uint32_t next = 0;
    std::vector<std::pair<NodeId, std::size_t>> stack;
    for (const NodeId r : nh_.root_order) {
      stack.emplace_back(r, 0);
      lo_[r] = next;
      while (!stack.empty()) {
        auto& [v, k] = stack.back();
        const auto& nd = f_.nodes[v];
        if (nd.children.empty()) {
          pos_[nd.state] = next;
          layout_[next] = nd.state;
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
    if (next != n_) throw std::invalid_argument("evaluate: forest leaves are not the states");
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
    return {layout_.data() + lo_[v], layout_.data() + hi_[v]};
  }

  void split(std::uint32_t pos) {
    if (sf_.split(pos)) bump("eval_interval_splits");
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
    if (!less(v, sf_.key(q))) return false;
    sf_.decrease(q, v);
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
        // [Deviation candidate, proposed D-16] Exact arithmetic cannot get here (l.798-800);
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
      const StateId y = layout_[q];
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
    const SfmMin<W> m = sf_.query(lo_[ch[r0]]);  // l.4, Lemma 3.11
    if (!due(m.key, F)) {                        // passed for this portion (l.821)
      bump("eval_failed_run_visits");
      frames_[fi].p = pieces_[p].next;
      return;
    }
    // l.5: the child of X containing the minimising leaf, by a level-ancestor query (l.817).
    const NodeId leaf = f_.leaf_of_state[layout_[m.index]];
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
    const W key = sf_.key(pos_[x]);
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
      const SfmMin<W> m = sf_.query(lo_[i]);
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
  SplitFindmin<W> sf_;

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

}  // namespace sssp::paper
