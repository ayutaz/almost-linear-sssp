// Property tests of the Evaluate invariants I1-I8 (research notes 02 sec 2.8, derived from
// Sec 4.3, Alg 2-3 and Lemma 4.6 of arXiv 2611.34911 v1) with Catch2 generators and the
// deterministic shrinking driver of shrink_support.hpp [P4-PROPERTY-EVAL].
//
// Already covered elsewhere (not repeated here):
//   * tests/support/paper/eval_invariants.hpp (InvariantChecker) checks I1-I8 and the hook
//     protocol on one run; tests/unit/paper/hierarchy/test_eval_invariants.cpp shows that it
//     detects each violation kind on hand-made traces;
//   * tests/differential/paper/test_evaluate_vs_dijkstra.cpp runs that checker on 10^4 uint64
//     reference hierarchies (N <= 60), the base-case catalogue (11 fixed shapes x 5 seeds),
//     zero-length arcs, integer-valued doubles, and supplemental offers off (D-14);
//   * test_evaluate_property.cpp / test_portioned_evaluate_property.cpp count visits (Lemma 4.7,
//     Lemma 7.7), not I1-I8.
// Gaps those tests leave, closed here:
//   * no I1-I8 failure there is reported with a seed and a shrunk input. Here every case is a
//     pure function of a seed (state graph or physical graph as an edge list), and a failure
//     prints the seed, the original and a 1-minimal input (replay: sh::shrink_case on
//     generate(seed) with the same violation function);
//   * random base-case (Alg 4) hierarchies of random physical graphs, not only fixed shapes,
//     and non-integer doubles whose Sec 4 arithmetic is exact (lengths k/8: every sum, the
//     origin subtractions D-19 and the gaps D-15 are exact), so I1-I8 hold with exact equality;
//   * parts of I3 / I4 / I5 that InvariantChecker does not check ("Not checked" in its header).
//     ScanOrderObserver below checks them on the same hook calls:
//       I3 (rank part)  an arc relaxed inside Scan(X) from a child of X to a lower-ranked child
//                       of X offers a value >= the nominal slot end b ("A backward child arc of a
//                       slot forces any newly generated value to be at least the nominal right
//                       endpoint", proof of Lemma 4.6; backward arcs are not light, Lemma 4.5);
//       cursor order    within one Scan(X, I, c) the children of X acted on (opened, scanned,
//                       emitted, retired) never go back in rank (Alg 2 l.2 "in cursor order",
//                       l.7 "Resume at the first piece just created"; Sec 4.3 "Keep the cursor
//                       on the first new piece so all earlier ranks have their opportunity");
//       I4 (direct)     an emission inside Scan(X) is at a leaf CHILD of X (Alg 2 l.11; the
//                       checker only asks for a state below X), and the nominal slots of X at
//                       which states below X emit never go back (research notes 2.8 I4, which
//                       the checker covers only through the scan slots);
//       I5 (closure)    at the end of Evaluate every state that received a finite offer has
//                       emitted or was retired (Alg 3 l.11 runs each root until its children
//                       finish or its expiry is covered, and expiry retires the rest, Sec 4.3
//                       "Mark each newly retired leaf"; offers never go to earlier roots, H4).
//
// Inputs (positive lengths; Def 4.1 / 4.2 are checked with validate_seeded_ordered_hierarchy,
// and a case that fails them is a precondition miss, not a counterexample):
//   state-graph  Case.n states, Case.edges directed arcs, Case.b physical vertices (state x is
//                a copy of vertex x mod b); threshold-SCC reference hierarchy; physical graph =
//                the arc multigraph, optionally with one shortcut of half length per
//                odd-numbered arc k, from its tail's vertex to the head's vertex of arc k-1
//                (fewer protected states, supplemental offers matter, copies expire); an
//                ordinary source at state 0, a weighted source at state n/2, a tight offer
//                delta at state n-1 and a non-tight offer delta + 5 at state (n-1)/2.
//   base-case    Case.n vertices and Case.edges a physical multigraph; Alg 4 base-case
//                hierarchy; Case.b sources (1..3), source i at vertex 5i mod n, offset 2i.
// Each block of seeds asserts that the checks were exercised (protected states, emissions,
// backward relaxations, cursor advances, second-rule origins, supplemental offers iff enabled)
// and reports the counts with WARN. Inclusive deadlines below the portion end and expiries with
// a pending key are rare at the regular sizes, so three pinned state graphs that reach them (the
// shrinker's minima under two injected faults, see the test) run in the regular suite as well.

#include <algorithm>
#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "property/paper/shrink_support.hpp"
#include "sssp/graph.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/paper/seeded_instances.hpp"
#include "support/prng.hpp"

using sssp::Edge;
using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::ArcId;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::EvalResult;
using sssp::paper::NodeId;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::StateId;
using sssp::test::SplitMix64;
using namespace sssp::test::paper;
namespace sh = sssp::test::shrink;

namespace {

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

/// Length k >= 0 of the test grid: k itself for integers, k/8 for double (exact, dyadic).
template <class W>
W len(std::uint64_t k) {
  if constexpr (std::is_floating_point_v<W>) {
    return static_cast<W>(k) * W{0.125};
  } else {
    return static_cast<W>(k);
  }
}

/// Candidate weights for shrinking (smallest grid lengths).
template <class W>
std::vector<W> shrink_weights() {
  return {len<W>(1), len<W>(2), len<W>(4)};
}

// ------------------------------------------------------------------ extra hook checks

/// Checks the parts of I3 / I4 / I5 that InvariantChecker leaves out (file comment), on the
/// normalised and ordered hierarchy whose NodeIds the hooks report (evaluate_normalized).
template <class W>
class ScanOrderObserver final : public EvalObserver<W> {
 public:
  struct Coverage {
    std::size_t backward_relaxations = 0;  ///< relaxations to a lower-ranked child (I3 rank)
    std::size_t cursor_advances = 0;       ///< actions at a strictly higher rank in one Scan
    std::size_t retired_states = 0;        ///< states retired without emitting
    std::size_t retired_pending = 0;       ///< ... of which held a finite key
  };

  explicit ScanOrderObserver(const NormalizedHierarchy<W>& nh) : nh_(&nh) {
    const auto& f = nh.forest;
    const std::size_t n = f.num_states();
    key_.assign(n, kInf<W>);
    emitted_.assign(n, false);
    retired_.assign(n, false);
    const std::size_t nn = f.nodes.size();
    has_emit_slot_.assign(nn, false);
    emit_slot_.assign(nn, W{});
  }

  [[nodiscard]] std::size_t total() const noexcept { return total_; }
  [[nodiscard]] const Coverage& coverage() const noexcept { return cov_; }
  [[nodiscard]] std::string summary() const {
    std::string s = std::to_string(total_) + " scan-order violation(s)\n";
    for (const auto& m : messages_) s += m + "\n";
    return s;
  }

  void on_initial_offer(StateId x, W v) override { offer(x, v); }
  void on_supplemental_offer(StateId x, W v) override { offer(x, v); }

  void on_relax_offer(ArcId arc, W v) override {
    const auto& arcs = nh_->forest.graph.arcs;
    if (arc >= arcs.size()) return;  // protocol, reported by InvariantChecker
    const StateId x = arcs[arc].tail;
    const StateId y = arcs[arc].head;
    // I3 (rank part): at every scanned node X containing both ends in different children, an
    // arc into a lower-ranked child carries a value >= the nominal slot end b of X.
    for (const Frame& fr : frames_) {
      const NodeId cx = child_containing(fr.node, x);
      const NodeId cy = child_containing(fr.node, y);
      if (cx == kNoNode || cy == kNoNode || cx == cy) continue;
      if (nh_->rank[cy] < nh_->rank[cx]) {
        ++cov_.backward_relaxations;
        if (v < fr.slot_hi) {
          add("I3 (rank): arc " + std::to_string(arc) + " " + std::to_string(x) + "->" +
              std::to_string(y) + " into a lower-ranked child of node " + std::to_string(fr.node) +
              " offers " + fmt(v) + " below the nominal slot end " + fmt(fr.slot_hi));
        }
      }
    }
    offer(y, v);
  }

  void on_open(std::uint32_t node, std::span<const StateId> /*states*/, W /*d0*/, W /*origin*/,
               W /*unit*/, W /*mass*/) override {
    act(node, "open");
  }

  void on_scan_begin(std::uint32_t node, std::span<const StateId> /*states*/, W /*lo*/, W /*hi*/,
                     W slot_lo, W slot_hi, W /*deadline*/) override {
    act(node, "scan");
    frames_.push_back(Frame{node, slot_lo, slot_hi, false, 0});
  }

  void on_scan_end(std::uint32_t node, std::span<const StateId> /*states*/) override {
    if (!frames_.empty() && frames_.back().node == node) frames_.pop_back();
  }

  void on_emit(StateId x, W /*v*/) override {
    if (x >= emitted_.size()) return;
    const NodeId leaf = nh_->forest.leaf_of_state[x];
    if (!frames_.empty() && nh_->forest.nodes[leaf].parent != frames_.back().node) {
      add("I4 (direct): state " + std::to_string(x) + " emits inside Scan(node " +
          std::to_string(frames_.back().node) + ") but its leaf is not a child of that node");
    }
    act(leaf, "emit");
    // I4: per node, the nominal slots at which states below it emit never go back.
    for (const Frame& fr : frames_) {
      if (has_emit_slot_[fr.node] && fr.slot_lo < emit_slot_[fr.node]) {
        add("I4 (slots): state " + std::to_string(x) + " emits in slot " + fmt(fr.slot_lo) +
            " of node " + std::to_string(fr.node) + " after an emission in slot " +
            fmt(emit_slot_[fr.node]));
      }
      has_emit_slot_[fr.node] = true;
      emit_slot_[fr.node] = fr.slot_lo;
    }
    emitted_[x] = true;
  }

  void on_retire(std::uint32_t node, std::span<const StateId> states) override {
    act(node, "retire");
    for (const StateId x : states) {
      if (x >= emitted_.size() || emitted_[x] || retired_[x]) continue;
      retired_[x] = true;
      ++cov_.retired_states;
      if (key_[x] != kInf<W>) ++cov_.retired_pending;
    }
  }

  void on_finish() override {
    // I5 (closure): nothing finite is left pending at a live state.
    for (StateId x = 0; x < key_.size(); ++x) {
      if (key_[x] != kInf<W> && !emitted_[x] && !retired_[x]) {
        add("I5 (closure): state " + std::to_string(x) + " ends with pending key " + fmt(key_[x]) +
            ", neither emitted nor retired");
      }
    }
  }

 private:
  struct Frame {
    NodeId node;
    W slot_lo;
    W slot_hi;
    bool has_rank;
    std::uint32_t last_rank;
  };

  static std::string fmt(W w) {
    if (w == kInf<W>) return "inf";
    std::ostringstream os;
    os.precision(17);
    os << w;
    return os.str();
  }

  /// The child of X whose subtree contains state x, or kNoNode if x is not below X.
  [[nodiscard]] NodeId child_containing(NodeId x_node, StateId x) const {
    const auto& f = nh_->forest;
    if (x >= f.leaf_of_state.size()) return kNoNode;
    NodeId v = f.leaf_of_state[x];
    while (v != kNoNode && f.nodes[v].parent != x_node) v = f.nodes[v].parent;
    return v;
  }

  /// Cursor order: an action at a child of the innermost scanned node never goes back in rank.
  void act(NodeId node, const char* what) {
    if (frames_.empty() || node >= nh_->forest.nodes.size()) return;
    Frame& fr = frames_.back();
    if (nh_->forest.nodes[node].parent != fr.node) return;  // nesting: InvariantChecker
    const std::uint32_t r = nh_->rank[node];
    if (fr.has_rank && r < fr.last_rank) {
      add(std::string("cursor order: ") + what + " at rank " + std::to_string(r) + " (node " +
          std::to_string(node) + ") after rank " + std::to_string(fr.last_rank) +
          " in the same Scan of node " + std::to_string(fr.node));
    }
    if (fr.has_rank && r > fr.last_rank) ++cov_.cursor_advances;
    if (!fr.has_rank || r > fr.last_rank) fr.last_rank = r;
    fr.has_rank = true;
  }

  void offer(StateId x, W v) {
    if (x < key_.size() && !emitted_[x] && !retired_[x] && v < key_[x]) key_[x] = v;
  }

  void add(std::string m) {
    ++total_;
    if (messages_.size() < 16) messages_.push_back(std::move(m));
  }

  const NormalizedHierarchy<W>* nh_;
  std::vector<W> key_;
  std::vector<bool> emitted_;
  std::vector<bool> retired_;
  std::vector<bool> has_emit_slot_;
  std::vector<W> emit_slot_;
  std::vector<Frame> frames_;
  std::vector<std::string> messages_;
  std::size_t total_ = 0;
  Coverage cov_;
};

/// Forwards every hook call to two observers.
template <class W>
class Tee final : public EvalObserver<W> {
 public:
  Tee(EvalObserver<W>& a, EvalObserver<W>& b) : a_(&a), b_(&b) {}
  void on_initial_offer(StateId x, W v) override {
    a_->on_initial_offer(x, v);
    b_->on_initial_offer(x, v);
  }
  void on_relax_offer(ArcId arc, W v) override {
    a_->on_relax_offer(arc, v);
    b_->on_relax_offer(arc, v);
  }
  void on_supplemental_offer(StateId x, W v) override {
    a_->on_supplemental_offer(x, v);
    b_->on_supplemental_offer(x, v);
  }
  void on_open(std::uint32_t node, std::span<const StateId> states, W d0, W origin, W unit,
               W mass) override {
    a_->on_open(node, states, d0, origin, unit, mass);
    b_->on_open(node, states, d0, origin, unit, mass);
  }
  void on_scan_begin(std::uint32_t node, std::span<const StateId> states, W lo, W hi, W slot_lo,
                     W slot_hi, W deadline) override {
    a_->on_scan_begin(node, states, lo, hi, slot_lo, slot_hi, deadline);
    b_->on_scan_begin(node, states, lo, hi, slot_lo, slot_hi, deadline);
  }
  void on_scan_end(std::uint32_t node, std::span<const StateId> states) override {
    a_->on_scan_end(node, states);
    b_->on_scan_end(node, states);
  }
  void on_emit(StateId x, W v) override {
    a_->on_emit(x, v);
    b_->on_emit(x, v);
  }
  void on_retire(std::uint32_t node, std::span<const StateId> states) override {
    a_->on_retire(node, states);
    b_->on_retire(node, states);
  }
  void on_finish() override {
    a_->on_finish();
    b_->on_finish();
  }

 private:
  EvalObserver<W>* a_;
  EvalObserver<W>* b_;
};

/// Counts over the cases of one generator block, to show that the checks are not vacuous.
struct Coverage {
  std::size_t valid = 0;
  std::size_t rejected = 0;  ///< precondition (Def 4.1 / 4.2) not met
  std::size_t protected_states = 0;
  std::size_t emissions = 0;
  std::size_t supplemental = 0;
  std::size_t second_rule_origins = 0;  ///< I7: origin != D0 at a non-root
  std::size_t retired_states = 0;
  std::size_t retired_pending = 0;
  std::size_t backward_relaxations = 0;
  std::size_t cursor_advances = 0;
  std::size_t capped_scans = 0;  ///< Scan(X, I, c) with c below the portion end (deadline cap)

  [[nodiscard]] std::string str() const {
    std::ostringstream os;
    os << "valid " << valid << ", rejected " << rejected << ", protected " << protected_states
       << ", emissions " << emissions << ", supplemental offers " << supplemental
       << ", second-rule origins " << second_rule_origins << ", retired states " << retired_states
       << " (" << retired_pending << " with a finite key)"
       << ", backward relaxations " << backward_relaxations << ", cursor advances "
       << cursor_advances << ", deadline-capped scans " << capped_scans;
    return os.str();
  }
};

/// Counts the hook events the coverage report needs (supplemental offers, I7 second rule).
template <class W>
class EventCounter final : public EvalObserver<W> {
 public:
  std::size_t supplemental = 0;
  std::size_t second_rule = 0;
  std::size_t capped_scans = 0;
  void on_supplemental_offer(StateId /*x*/, W /*v*/) override { ++supplemental; }
  void on_scan_begin(std::uint32_t /*node*/, std::span<const StateId> /*states*/, W /*lo*/, W hi,
                     W /*slot_lo*/, W /*slot_hi*/, W deadline) override {
    if (deadline < hi) ++capped_scans;
  }
  void on_open(std::uint32_t /*node*/, std::span<const StateId> /*states*/, W d0, W origin,
               W /*unit*/, W /*mass*/) override {
    if (origin != d0) ++second_rule;
  }
};

struct Config {
  bool supplemental = true;
  /// Test of the reporting pipeline only: raise delta of every protected state in the oracle,
  /// so that every correct run violates I2 (shrinking self-test below).
  bool inflate_protected_targets = false;
};

/// The property: I1-I8 (InvariantChecker) and the extra checks (ScanOrderObserver) hold on one
/// evaluation, and evaluate() agrees with evaluate_normalized(). std::nullopt if the property
/// holds or the input misses the precondition (Def 4.1 / 4.2).
template <class W>
std::optional<std::string> violation(const SeededInstance<W>& inst, const Config& cfg,
                                     Coverage* cov) {
  const auto val = sssp::paper::validate_seeded_ordered_hierarchy<W>(inst.seeded, inst.physical);
  if (!val.ok()) {
    if (cov != nullptr) ++cov->rejected;
    return std::nullopt;
  }
  auto ref = compute_eval_reference<W>(inst.seeded, inst.physical);
  if (cfg.inflate_protected_targets) {
    for (StateId x = 0; x < ref.state_target.size(); ++x) {
      if (ref.is_protected[x]) ref.state_target[x] = ref.state_target[x] + len<W>(1);
    }
  }
  auto nh = sssp::paper::normalize<W>(inst.seeded.hierarchy);
  sssp::paper::order_children<W>(nh);

  CheckerOptions copt;
  copt.supplemental_offers = cfg.supplemental;
  InvariantChecker<W> checker(inst.seeded, ref, copt);
  ScanOrderObserver<W> order(nh);
  EventCounter<W> events;
  Tee<W> inner(order, events);
  Tee<W> tee(checker, inner);
  ForwardingObserver<EvalHooks<W>, W> adapter(tee);
  EvalOptions<W> opt;
  opt.supplemental_offers = cfg.supplemental;
  opt.hooks = &adapter;
  EvalResult<W> r;
  EvalResult<W> plain;
  try {
    r = sssp::paper::evaluate_normalized<W>(nh, inst.seeded.offers, opt);
    EvalOptions<W> popt;
    popt.supplemental_offers = cfg.supplemental;
    plain = sssp::paper::evaluate<W>(inst.seeded, popt);
  } catch (const std::exception& e) {
    return std::string("Evaluate threw on a valid seeded hierarchy: ") + e.what();
  }
  std::string msg;
  if (!checker.report().ok()) msg += checker.report().summary();
  if (order.total() != 0) msg += order.summary();
  if (r.emissions != checker.emissions()) msg += "result emissions differ from the hook trace\n";
  if (plain.value != r.value || plain.emissions != r.emissions) {
    msg += "evaluate() and evaluate_normalized() differ\n";
  }
  if (cov != nullptr) {
    ++cov->valid;
    cov->protected_states += ref.num_protected();
    cov->emissions += r.emissions.size();
    cov->supplemental += events.supplemental;
    cov->second_rule_origins += events.second_rule;
    cov->retired_states += order.coverage().retired_states;
    cov->retired_pending += order.coverage().retired_pending;
    cov->backward_relaxations += order.coverage().backward_relaxations;
    cov->cursor_advances += order.coverage().cursor_advances;
    cov->capped_scans += events.capped_scans;
  }
  if (msg.empty()) return std::nullopt;
  return msg;
}

// ------------------------------------------------------------------ inputs from a Case

/// State-graph family (file comment). Requires c.n >= 1.
template <class W>
SeededInstance<W> state_graph_instance(const sh::Case<W>& c, bool shortcuts) {
  const std::size_t n = std::max<std::size_t>(c.n, 1);
  sssp::paper::StateGraph<W> g;
  g.num_states = n;
  for (const auto& e : c.edges) g.arcs.push_back({e.u, e.v, e.weight});
  const std::size_t p = std::clamp<std::size_t>(c.b, 1, n);
  std::vector<Vertex> proj(n);
  for (std::size_t x = 0; x < n; ++x) proj[x] = static_cast<Vertex>(x % p);
  SeededInstance<W> inst;
  inst.name = shortcuts ? "state-graph/shortcuts" : "state-graph";
  inst.num_physical = p;
  const auto base = sssp::paper::physical_graph_of_arcs<W>(g, proj, p);
  std::vector<Edge<W>> edges(base.edges().begin(), base.edges().end());
  if (shortcuts) {
    // Shortcuts only shorten physical distances, so every arc and link stays a long enough
    // walk (Def 4.1 H1, H5); they make state arcs longer than delta differences, so fewer
    // states are protected and unprotected copies can expire.
    for (std::size_t k = 1; k < g.arcs.size(); k += 2) {
      const auto& a = g.arcs[k];
      W half = a.length / 2;
      if constexpr (!std::is_floating_point_v<W>) half = std::max<W>(half, 1);
      edges.push_back({proj[a.tail], proj[g.arcs[k - 1].head], half});
    }
  }
  inst.physical = sssp::CsrGraph<W>::from_edges(p, std::move(edges));
  inst.seeded.hierarchy = sssp::paper::build_threshold_scc_hierarchy<W>(g, proj);
  auto& s = inst.seeded;
  const auto source = [&](StateId x, W offset) {
    s.sources.push_back({proj[x], offset});
    s.offers.push_back({x, offset});  // feasible: delta(pi(x)) <= offset
  };
  source(0, W{0});
  if (n >= 3) source(static_cast<StateId>(n / 2), len<W>(3));
  if (n >= 2) {
    const auto target = target_distances<W>(inst.physical, s.sources);
    const W t = target[proj[n - 1]];
    if (t != kInf<W>) s.offers.push_back({static_cast<StateId>(n - 1), t});  // tight offer
    const auto y = static_cast<StateId>((n - 1) / 2);
    const W u = target[proj[y]];
    if (u != kInf<W>) s.offers.push_back({y, u + len<W>(5)});  // feasible, not tight
  }
  return inst;
}

/// Base-case family (file comment). Requires positive weights.
template <class W>
SeededInstance<W> base_case_from(const sh::Case<W>& c) {
  const std::size_t n = std::max<std::size_t>(c.n, 1);
  const std::size_t k = std::clamp<std::size_t>(c.b, 1, 3);
  std::vector<PhysicalSource<W>> sources;
  sources.reserve(k);
  for (std::size_t i = 0; i < k; ++i) {
    sources.push_back({static_cast<Vertex>((5 * i) % n), len<W>(2 * i)});
  }
  return base_case_instance<W>("base-case", n, c.edges, std::move(sources));
}

template <class W>
sh::Case<W> random_state_graph(std::uint64_t seed, std::size_t max_states) {
  SplitMix64 rng(seed);
  sh::Case<W> c;
  c.n = 2 + rng.below(max_states - 1);
  const std::size_t m = rng.below(3 * c.n + 1);
  const std::uint64_t max_k = rng.below(2) == 0 ? 2 : 12;  // tie-heavy or spread lengths
  for (std::size_t i = 0; i < m; ++i) {
    c.edges.push_back({static_cast<Vertex>(rng.below(c.n)), static_cast<Vertex>(rng.below(c.n)),
                       len<W>(rng.between(1, max_k))});
  }
  c.b = 1 + rng.below(c.n);
  return c;
}

template <class W>
sh::Case<W> random_physical(std::uint64_t seed, std::size_t max_vertices) {
  SplitMix64 rng(seed);
  sh::Case<W> c;
  c.n = 1 + rng.below(max_vertices);
  const std::size_t m = rng.below(3 * c.n + 1);
  const std::uint64_t max_k = rng.below(2) == 0 ? 2 : 20;
  for (std::size_t i = 0; i < m; ++i) {
    c.edges.push_back({static_cast<Vertex>(rng.below(c.n)), static_cast<Vertex>(rng.below(c.n)),
                       len<W>(rng.between(1, max_k))});
  }
  c.b = 1 + rng.below(3);
  return c;
}

enum class Family : std::uint8_t { kStateGraph, kStateGraphShortcuts, kBaseCase };

template <class W>
SeededInstance<W> instance_of(Family fam, const sh::Case<W>& c) {
  switch (fam) {
    case Family::kStateGraph:
      return state_graph_instance<W>(c, false);
    case Family::kStateGraphShortcuts:
      return state_graph_instance<W>(c, true);
    case Family::kBaseCase:
      break;
  }
  return base_case_from<W>(c);
}

template <class W>
sh::Case<W> generate(Family fam, std::uint64_t seed, std::size_t max_n) {
  return fam == Family::kBaseCase ? random_physical<W>(seed, max_n)
                                  : random_state_graph<W>(seed, max_n);
}

/// Runs one block of seeds, shrinks the first failure, and checks that the block was not
/// vacuous.
template <class W>
void run_block(Family fam, const Config& cfg, std::uint64_t first, std::uint64_t count,
               std::size_t max_n) {
  Coverage cov;
  const auto out = sh::run_seeded<W>(
      first, count, [&](std::uint64_t s) { return generate<W>(fam, s, max_n); },
      [&](const sh::Case<W>& c) { return violation<W>(instance_of<W>(fam, c), cfg, &cov); },
      shrink_weights<W>());
  INFO("family " << static_cast<int>(fam) << ", supplemental " << cfg.supplemental << ", seeds "
                 << first << " .. " << first + count - 1);
  INFO(out.report());
  REQUIRE_FALSE(out.failing_seed.has_value());
  INFO(cov.str());
  // The precondition rejects few generated cases, and the checks see the events they test.
  CHECK(cov.valid >= count * 9 / 10);
  CHECK(cov.protected_states > 0);
  CHECK(cov.emissions > 0);
  CHECK(cov.cursor_advances > 0);
  CHECK(cov.backward_relaxations > 0);
  // Retirement of unemitted states (expiry) is rare at these sizes (about one case in several
  // hundred, and never in the base case, where every reachable state is protected), so it is
  // reported, not required; I6 on retirement is exercised by the differential tests.
  WARN(cov.str());
  if (cfg.supplemental) {
    CHECK(cov.supplemental > 0);
  } else {
    CHECK(cov.supplemental == 0);
  }
  CHECK(cov.second_rule_origins > 0);
}

constexpr std::uint64_t kBlock = 250;

}  // namespace

TEMPLATE_TEST_CASE("evaluate invariants: I1-I8 and Scan order on seeded state-graph hierarchies",
                   "[evaluate][property][invariants]", std::uint64_t, double) {
  using W = TestType;
  const bool supplemental = GENERATE(true, false);
  const bool shortcuts = GENERATE(false, true);
  const std::uint64_t block = GENERATE(range<std::uint64_t>(0, 4));
  Config cfg;
  cfg.supplemental = supplemental;
  run_block<W>(shortcuts ? Family::kStateGraphShortcuts : Family::kStateGraph, cfg,
               0xe1a10000ULL + block * kBlock, kBlock, 40);
}

TEMPLATE_TEST_CASE("evaluate invariants: I1-I8 and Scan order on seeded base-case hierarchies",
                   "[evaluate][property][invariants]", std::uint64_t, double) {
  using W = TestType;
  const bool supplemental = GENERATE(true, false);
  const std::uint64_t block = GENERATE(range<std::uint64_t>(0, 4));
  Config cfg;
  cfg.supplemental = supplemental;
  run_block<W>(Family::kBaseCase, cfg, 0xe1a20000ULL + block * kBlock, kBlock, 20);
}

TEMPLATE_TEST_CASE("evaluate invariants: larger seeded hierarchies [slow]",
                   "[evaluate][property][invariants][slow]", std::uint64_t, double) {
  using W = TestType;
  const auto fam = GENERATE(Family::kStateGraph, Family::kStateGraphShortcuts, Family::kBaseCase);
  const bool supplemental = GENERATE(true, false);
  const std::uint64_t block = GENERATE(range<std::uint64_t>(0, 4));
  Config cfg;
  cfg.supplemental = supplemental;
  run_block<W>(fam, cfg, 0xe1a30000ULL + block * 2 * kBlock, 2 * kBlock, 60);
}

TEMPLATE_TEST_CASE("evaluate invariants: pinned inputs that reach the deadline and expiry branches",
                   "[evaluate][property][invariants]", std::uint64_t, double) {
  // State-graph cases (lengths in grid units, len<W>) that random inputs of the regular size
  // rarely produce: an emission bounded by an inherited inclusive deadline c < b' (Sec 4.3
  // "Permanent deadlines are inclusive") and a child that expires with a pending key (Alg 2
  // l.12). They are the minimal counterexamples the shrinker printed (supplemental offers off,
  // [slow] seeds 3785556746 and 3785557787, and 3785557915) when, during P4-PROPERTY-EVAL, the
  // deadline test of due() or the retirement of an expired child was removed from
  // evaluate.cpp; on the unmodified code they satisfy every check.
  using W = TestType;
  struct Pinned {
    std::size_t n;
    std::vector<std::array<std::uint64_t, 3>> arcs;
  };
  const std::vector<Pinned> pinned = {
      {58, {{29, 40, 1}, {54, 40, 1}, {33, 29, 1}, {40, 54, 12}, {54, 33, 4}}},
      {25,
       {{8, 17, 1},
        {15, 24, 3},
        {6, 20, 1},
        {20, 7, 1},
        {7, 8, 12},
        {17, 6, 4},
        {20, 15, 2},
        {24, 6, 12}}},
      {12,
       {{9, 4, 1},
        {4, 3, 1},
        {2, 3, 1},
        {6, 1, 1},
        {0, 2, 4},
        {1, 10, 1},
        {4, 11, 1},
        {3, 9, 4},
        {10, 3, 11},
        {11, 6, 1},
        {9, 0, 1}}},
  };
  Coverage cov;
  for (const auto& pc : pinned) {
    sh::Case<W> c;
    c.n = pc.n;
    c.b = 1;
    for (const auto& [u, v, k] : pc.arcs) {
      c.edges.push_back({static_cast<Vertex>(u), static_cast<Vertex>(v), len<W>(k)});
    }
    for (const Family fam : {Family::kStateGraph, Family::kStateGraphShortcuts}) {
      for (const bool supplemental : {true, false}) {
        Config cfg;
        cfg.supplemental = supplemental;
        INFO(sh::describe(c) << " family " << static_cast<int>(fam) << " supplemental "
                             << supplemental);
        const auto v = violation<W>(instance_of<W>(fam, c), cfg, &cov);
        INFO(v.value_or("ok"));
        CHECK_FALSE(v.has_value());
      }
    }
  }
  INFO(cov.str());
  CHECK(cov.valid == pinned.size() * 2 * 2);
  CHECK(cov.capped_scans > 0);
  CHECK(cov.retired_pending > 0);
}

TEST_CASE("evaluate invariants: a failure is reported with its seed and a reproducible minimum",
          "[evaluate][property][invariants][shrink]") {
  // Self-test of the reporting pipeline on the real generator, checker and shrinker: with an
  // oracle whose protected deltas are raised by one, every valid case violates I2 (a protected
  // state emits its true delta). The driver must stop at the first seed, shrink to the smallest
  // state graph that still has a protected state (one state, no arcs) and reproduce exactly
  // the same report from the same seed.
  using W = std::uint64_t;
  Config cfg;
  cfg.inflate_protected_targets = true;
  const auto gen = [](std::uint64_t s) { return random_state_graph<W>(s, 24); };
  const auto bad = [&](const sh::Case<W>& c) {
    return violation<W>(state_graph_instance<W>(c, false), cfg, nullptr);
  };
  const std::uint64_t first = 0xe1a40000ULL;
  const auto a = sh::run_seeded<W>(first, 50, gen, bad, shrink_weights<W>());
  INFO(a.report());
  REQUIRE(a.failing_seed.has_value());
  const std::uint64_t seed = a.failing_seed.value_or(0);  // engaged: REQUIRE above
  CHECK(seed == first);
  CHECK(a.shrink.steps > 0);
  CHECK(a.minimal.n == 1);
  CHECK(a.minimal.edges.empty());
  CHECK(a.report().find("seed " + std::to_string(first)) != std::string::npos);
  CHECK(a.shrink.message.find("I2") != std::string::npos);
  // Replay: the same seed regenerates the original case and the same minimum.
  CHECK(gen(seed) == a.original);
  sh::Case<W> replay = gen(seed);
  const auto again = sh::shrink_case(replay, bad, shrink_weights<W>());
  CHECK(replay == a.minimal);
  CHECK(again.message == a.shrink.message);
  CHECK(bad(a.minimal).has_value());
  const auto b = sh::run_seeded<W>(first, 50, gen, bad, shrink_weights<W>());
  CHECK(b.report() == a.report());
  // With the true oracle the same seeds pass.
  Config ok;
  const auto good = sh::run_seeded<W>(
      first, 50, gen,
      [&](const sh::Case<W>& c) {
        return violation<W>(state_graph_instance<W>(c, false), ok, nullptr);
      },
      shrink_weights<W>());
  INFO(good.report());
  CHECK_FALSE(good.failing_seed.has_value());
}
