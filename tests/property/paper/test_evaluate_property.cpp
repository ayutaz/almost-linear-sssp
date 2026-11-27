// Property tests of the visit accounting of Lemma 4.7 (Sec 4.3, l.922-948), observed through the
// EvalHooks of Evaluate and the returned EvalResult, and on the counters Evaluate records in
// Diagnostics::counters.
//
// Exact bounds stated in the proof of Lemma 4.7 / Sec 4.3:
//   * a child is registered at most once (it then stays a registered piece until unlinked), so
//     registrations <= non-root nodes; at most two interval splits per registration (l.817),
//     plus the root boundaries of Alg 3 l.3;
//   * at most one supplemental offer per new subrun, two per registration (l.824);
//   * one opening per internal node; one emission per state (l.784); every arc is relaxed at
//     most once, by its tail's emission ("O(M) original relaxations", l.935);
//   * retirement visits every still-linked piece once (l.932-934); pieces ever created are
//     <= one per registered child + two subruns per registration + one run per opening.
// Paper-shaped bounds (the paper gives O(.); the constants are measured and reported with WARN,
// the assertions use a fixed factor 2 that the arguments below justify):
//   * slots of a node C: its slots start at its origin > D0 - nu(C) and run through the slot
//     containing its expiry D0 + S(C), i.e. at most S(C)/nu(C) + 2 slots (l.807-808), and a
//     slot is cut into at most two portions (l.799-801, Lemma 7.7 item 2):
//       scans <= 2 * sum_C (S(C)/nu(C) + 2);
//   * failed visits of a subrun created next to child i with neighbour j: O(1 + (S(i) +
//     g(i, j))/nu(X)) (l.925-928); its offer u <= max{b, z + S(i) + g} is due by the slot
//     containing u, at most 2 + (S(i) + g)/nu(X) slots after its creation. Over the two sides of
//     every adjacent pair (i, j) of X:
//       failed run visits <= 2 * sum_X sum_(i,j) (4 + (S(i) + S(j) + 2 g(i, j))/nu(X));
//   * visits of a registered child i: O(1 + S(i)/nu(X)) (l.930-931): it lives at most
//     S(i) + nu(i) <= S(i) + nu(X)/2 after its registration slot, so
//       child visits <= 2 * sum_(i non-root) (2 + S(i)/nu(parent(i))).
// Lemma 4.3 and 4.5 make all three sums O(N); their ratio to N is reported as well.
//
// Counters are compiled out with SSSP_DIAGNOSTICS=0 (ci-release, release): assertions on them run
// only when they are counted. What the hooks and the result expose is checked in every build
// (Tally below): Scan calls (on_scan_begin / on_scan_end, nested and balanced), openings,
// emissions, relaxations, supplemental offers and retirements, each against the bounds above
// (supplemental offers against 2 * non-root nodes, the registration bound). Registrations,
// interval splits, retired leaves, retirement piece visits, failed run visits and registered-
// child visits have no hook, so their bounds (and failed/N in the growth test) are checked only
// with diagnostics on.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "sssp/graph.hpp"
#include "support/graphs.hpp"
#include "support/paper/seeded_instances.hpp"

using sssp::paper::ArcId;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::NodeId;
using sssp::paper::StateId;
using sssp::test::SplitMix64;
using namespace sssp::test::paper;

namespace {

using U = std::uint64_t;

[[maybe_unused]] std::uint64_t counter(const sssp::Diagnostics& d, const std::string& name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

/// The paper-shaped sums of the file comment for one hierarchy (test arithmetic in double;
/// only used to bound counts), plus the normalised forest the hooks' node ids refer to.
struct Shape {
  std::size_t states = 0;
  std::size_t arcs = 0;
  std::size_t nodes = 0;
  std::size_t internal = 0;
  std::size_t roots = 0;
  double slot_sum = 0;     ///< sum_C (S(C)/nu(C) + 2)
  double fail_sum = 0;     ///< sum_X sum_(i,j) (4 + (S(i) + S(j) + 2 g)/nu(X))
  double child_sum = 0;    ///< sum_(i non-root) (2 + S(i)/nu(parent))
  std::vector<bool> leaf;  ///< per normalised node
};

Shape shape_of(const sssp::paper::SeededOrderedHierarchy<U>& s) {
  auto nh = sssp::paper::normalize<U>(s.hierarchy);
  sssp::paper::order_children<U>(nh);
  Shape sh;
  sh.states = s.hierarchy.num_states();
  sh.arcs = s.hierarchy.graph.arcs.size();
  sh.nodes = nh.num_nodes();
  sh.roots = nh.forest.roots.size();
  sh.leaf.assign(sh.nodes, false);
  for (NodeId x = 0; x < nh.num_nodes(); ++x) {
    if (nh.is_leaf(x)) {
      sh.leaf[x] = true;
      continue;
    }
    ++sh.internal;
    const double nu = static_cast<double>(nh.unit[x]);
    sh.slot_sum += static_cast<double>(nh.mass[x]) / nu + 2.0;
    const auto& ch = nh.forest.nodes[x].children;
    for (std::size_t k = 0; k < ch.size(); ++k) {
      sh.child_sum += 2.0 + static_cast<double>(nh.mass[ch[k]]) / nu;
      if (k + 1 < ch.size()) {
        sh.fail_sum +=
            4.0 + (static_cast<double>(nh.mass[ch[k]]) + static_cast<double>(nh.mass[ch[k + 1]]) +
                   2.0 * static_cast<double>(nh.gap_next[ch[k]])) /
                      nu;
      }
    }
  }
  return sh;
}

/// Counts what Evaluate shows through EvalHooks; available in every build. Order violations
/// (a Scan of an unopened node, an unmatched on_scan_end, a relaxation that does not follow
/// its tail's emission) are recorded and checked by the caller.
class Tally final : public EvalHooks<U> {
 public:
  explicit Tally(const sssp::paper::SeededOrderedHierarchy<U>& s) : arcs_(s.hierarchy.graph.arcs) {}

  std::size_t scans = 0;
  std::size_t scan_ends = 0;
  std::size_t supplemental = 0;
  std::size_t supplemental_outside_scan = 0;
  std::size_t scans_of_unopened = 0;
  std::size_t scans_after_retire = 0;
  std::size_t unmatched_scan_ends = 0;
  std::size_t misplaced_relaxations = 0;  ///< not right after the tail's emission, or v != d + l
  std::vector<NodeId> opened;
  std::vector<NodeId> retired;
  std::vector<std::pair<StateId, U>> emitted;
  std::vector<ArcId> relaxed;
  std::vector<NodeId> open_scans;  ///< the Scan stack at the end (empty after on_finish)
  bool finished = false;

  void on_supplemental_offer(StateId /*x*/, U /*value*/) override {
    ++supplemental;
    if (open_scans.empty()) ++supplemental_outside_scan;  // Alg 2 l.6 runs inside a Scan
  }
  void on_open(std::uint32_t node, std::span<const StateId> /*states*/, U /*d0*/, U /*origin*/,
               U /*unit*/, U /*mass*/) override {
    opened.push_back(node);
    mark(is_opened_, node);
  }
  void on_scan_begin(std::uint32_t node, std::span<const StateId> /*states*/, U /*lo*/, U /*hi*/,
                     U /*slot_lo*/, U /*slot_hi*/, U /*deadline*/) override {
    ++scans;
    // Alg 3 l.9-10 / Alg 2 l.9-10: a node is scanned only after it opened, never after it was
    // retired (Alg 2 l.12 unlinks it).
    if (!marked(is_opened_, node)) ++scans_of_unopened;
    if (marked(is_retired_, node)) ++scans_after_retire;
    open_scans.push_back(node);
  }
  void on_scan_end(std::uint32_t node, std::span<const StateId> /*states*/) override {
    ++scan_ends;
    if (open_scans.empty() || open_scans.back() != node) {
      ++unmatched_scan_ends;
    } else {
      open_scans.pop_back();
    }
  }
  void on_emit(StateId x, U value) override { emitted.emplace_back(x, value); }
  void on_relax_offer(ArcId arc, U value) override {
    relaxed.push_back(arc);
    // Alg 2 l.11 / Alg 3 l.7: the arc is relaxed by its tail's emission, with v = delta + l.
    if (emitted.empty() || arc >= arcs_.size() || arcs_[arc].tail != emitted.back().first ||
        value != emitted.back().second + arcs_[arc].length) {
      ++misplaced_relaxations;
    }
  }
  void on_retire(std::uint32_t node, std::span<const StateId> /*states*/) override {
    retired.push_back(node);
    mark(is_retired_, node);
  }
  void on_finish() override { finished = true; }

 private:
  static void mark(std::vector<bool>& v, std::uint32_t i) {
    if (i >= v.size()) v.resize(std::size_t{i} + 1, false);
    v[i] = true;
  }
  static bool marked(const std::vector<bool>& v, std::uint32_t i) { return i < v.size() && v[i]; }

  const std::vector<sssp::paper::StateArc<U>>& arcs_;
  std::vector<bool> is_opened_;
  std::vector<bool> is_retired_;
};

template <class T>
bool all_distinct(std::vector<T> v) {
  std::sort(v.begin(), v.end());
  return std::adjacent_find(v.begin(), v.end()) == v.end();
}

struct Ratios {
  double scans = 0;
  double failed = 0;
  double child = 0;
  double scans_per_n = 0;
  double failed_per_n = 0;
  double slot_sum_per_n = 0;
  double fail_sum_per_n = 0;
  void max_with(const Ratios& o) {
    scans = std::max(scans, o.scans);
    failed = std::max(failed, o.failed);
    child = std::max(child, o.child);
    scans_per_n = std::max(scans_per_n, o.scans_per_n);
    failed_per_n = std::max(failed_per_n, o.failed_per_n);
    slot_sum_per_n = std::max(slot_sum_per_n, o.slot_sum_per_n);
    fail_sum_per_n = std::max(fail_sum_per_n, o.fail_sum_per_n);
  }
  [[nodiscard]] std::string str() const {
    std::ostringstream os;
    os << "max scans/slot_sum " << scans << ", scans/N " << scans_per_n << ", slot_sum/N "
       << slot_sum_per_n << ", fail_sum/N " << fail_sum_per_n;
#if SSSP_DIAGNOSTICS
    os << ", failed/fail_sum " << failed << ", child/child_sum " << child << ", failed/N "
       << failed_per_n;
#else
    os << " (failed and child visits not counted: SSSP_DIAGNOSTICS=0)";
#endif
    return os.str();
  }
};

/// Checks the exact and shaped bounds of the file comment on one instance; returns the ratios.
Ratios check_counts(const SeededInstance<U>& inst) {
  INFO(inst.name << " seed " << inst.seed);
  Tally t(inst.seeded);
  EvalOptions<U> opt;
  opt.hooks = &t;
  const auto r = sssp::paper::evaluate<U>(inst.seeded, opt);
  const Shape sh = shape_of(inst.seeded);
  const auto& arcs = inst.seeded.hierarchy.graph.arcs;

  // Hook- and result-based checks (every build).
  CHECK(t.finished);
  CHECK(t.scan_ends == t.scans);
  CHECK(t.unmatched_scan_ends == 0);
  CHECK(t.open_scans.empty());
  CHECK(t.scans_of_unopened == 0);
  CHECK(t.scans_after_retire == 0);
  // One opening per internal node (Alg 3 l.9, Alg 2 l.9).
  CHECK(t.opened.size() <= sh.internal);
  CHECK(all_distinct(t.opened));
  for (const NodeId x : t.opened) CHECK((x < sh.nodes && !sh.leaf[x]));
  // One emission per state (l.784), the same sequence as the result.
  CHECK(t.emitted.size() <= sh.states);
  CHECK(t.emitted == r.emissions);
  std::vector<StateId> emitted_states;
  emitted_states.reserve(t.emitted.size());
  for (const auto& e : t.emitted) emitted_states.push_back(e.first);
  CHECK(all_distinct(emitted_states));
  // Every arc relaxed at most once, by its tail's emission (l.935); every arc leaving an emitted
  // state to another state is relaxed (Alg 2 l.11).
  CHECK(t.relaxed.size() <= sh.arcs);
  CHECK(all_distinct(t.relaxed));
  CHECK(t.misplaced_relaxations == 0);
  {
    std::vector<bool> is_emitted(sh.states, false);
    for (const StateId x : emitted_states) {
      if (x < sh.states) is_emitted[x] = true;
    }
    std::vector<bool> is_relaxed(sh.arcs, false);
    for (const ArcId a : t.relaxed) {
      if (a < sh.arcs) is_relaxed[a] = true;
    }
    std::size_t missed = 0;
    for (ArcId a = 0; a < sh.arcs; ++a) {
      if (arcs[a].tail != arcs[a].head && is_emitted[arcs[a].tail] && !is_relaxed[a]) ++missed;
    }
    CHECK(missed == 0);
  }
  // At most two supplemental offers per registration (l.824), registrations <= non-root nodes.
  CHECK(t.supplemental <= 2 * (sh.nodes - sh.roots));
  CHECK(t.supplemental_outside_scan == 0);
  // Retirements (completed or expired internal nodes, l.856-860): each node at most once.
  CHECK(t.retired.size() <= sh.internal);
  CHECK(all_distinct(t.retired));
  for (const NodeId x : t.retired) CHECK((x < sh.nodes && !sh.leaf[x]));
  const auto scans = static_cast<double>(t.scans);
  CHECK(scans <= 2.0 * sh.slot_sum);

  double failed = 0;
  double child = 0;
#if SSSP_DIAGNOSTICS
  // Counter-based checks: registrations, interval splits, retired leaves, retirement piece
  // visits, failed run visits and child visits have no hook.
  const auto& d = r.diag;
  // The hooks and the counters count the same events.
  CHECK(counter(d, "eval_scans") == t.scans);
  CHECK(counter(d, "eval_opens") == t.opened.size());
  CHECK(counter(d, "eval_emissions") == t.emitted.size());
  CHECK(counter(d, "eval_relaxations") == t.relaxed.size());
  CHECK(counter(d, "eval_supplemental_offers") == t.supplemental);
  CHECK(counter(d, "eval_completed_nodes") + counter(d, "eval_expired_nodes") >= t.retired.size());
  const auto regs = counter(d, "eval_registrations");
  // Exact bounds.
  CHECK(regs <= sh.nodes - sh.roots);
  CHECK(counter(d, "eval_interval_splits") <= 2 * regs + (sh.roots - 1));
  CHECK(counter(d, "eval_supplemental_offers") <= 2 * regs);
  CHECK(counter(d, "eval_opens") <= sh.internal);
  CHECK(counter(d, "eval_emissions") <= sh.states);
  CHECK(counter(d, "eval_relaxations") <= sh.arcs);
  CHECK(counter(d, "eval_completed_nodes") + counter(d, "eval_expired_nodes") <= sh.internal);
  CHECK(counter(d, "eval_retired_leaves") <= sh.states);
  // Pieces ever created: one per registered child, two subruns per registration, one run per
  // opening; retirement visits each at most once (l.932-934).
  CHECK(counter(d, "eval_retire_piece_visits") <= 3 * (sh.nodes - sh.roots) + sh.internal);
  CHECK(counter(d, "eval_interval_splits") == counter(d, "eval_sfm_splits"));
  // Paper-shaped bounds without a hook.
  failed = static_cast<double>(counter(d, "eval_failed_run_visits"));
  child = static_cast<double>(counter(d, "eval_child_visits"));
  CHECK(failed <= 2.0 * sh.fail_sum);
  CHECK(child <= 2.0 * sh.child_sum);
#endif
  Ratios q;
  q.scans = sh.slot_sum > 0 ? scans / sh.slot_sum : 0;
  q.failed = sh.fail_sum > 0 ? failed / sh.fail_sum : 0;
  q.child = sh.child_sum > 0 ? child / sh.child_sum : 0;
  const auto n = static_cast<double>(std::max<std::size_t>(1, sh.states));
  q.scans_per_n = scans / n;
  q.failed_per_n = failed / n;
  q.slot_sum_per_n = sh.slot_sum / n;
  q.fail_sum_per_n = sh.fail_sum / n;
  return q;
}

}  // namespace

TEST_CASE("evaluate property: Lemma 4.7 counts on reference hierarchies",
          "[evaluate][property][lemma47]") {
  Ratios worst;
  for (const auto& inst : reference_catalogue<U>(9000, 1500)) worst.max_with(check_counts(inst));
  WARN("Lemma 4.7 measured constants, reference catalogue: " << worst.str());
}

TEST_CASE("evaluate property: Lemma 4.7 counts on base-case hierarchies",
          "[evaluate][property][lemma47]") {
  Ratios worst;
  for (const std::uint64_t seed : {11ULL, 12ULL, 13ULL}) {
    for (const auto& inst : base_case_catalogue<U>(seed)) worst.max_with(check_counts(inst));
  }
  WARN("Lemma 4.7 measured constants, base-case catalogue: " << worst.str());
}

TEST_CASE("evaluate property: per-state counts stay bounded as N grows",
          "[evaluate][property][lemma47]") {
  // Lemma 4.7: total work O(M + N tau(N) + Q), and Lemma 4.3 / 4.5 make the shaped sums O(N).
  // Divided by N, neither the counts nor the sums may grow with N: the worst ratio at the
  // largest N must stay within twice the worst ratio at the smallest N (or 1). Two families:
  // threshold-SCC hierarchies of random state graphs (few distinct lengths: wide nodes) and
  // Alg 4 base-case hierarchies of random physical graphs (many lengths: deep Kruskal trees).
  // Scans are counted by the hooks (every build), failed run visits by a counter (diagnostics
  // only), the shaped sums from the structure (every build).
  for (const bool base_case : {false, true}) {
    std::vector<Ratios> per_n;
    for (const std::size_t n : {50, 100, 200, 400, 800}) {
      Ratios worst;
      for (std::uint64_t seed = 0; seed < 6; ++seed) {
        const std::uint64_t s = n * 1000 + seed;
        SplitMix64 rng(s);
        SeededInstance<U> inst;
        if (base_case) {
          inst = base_case_instance<U>("growth_base/" + std::to_string(n), n,
                                       sssp::test::random_edges<U>(n, 3 * n, 1, 1000, s),
                                       {{0, 0}, {static_cast<sssp::Vertex>(rng.below(n)), 5}}, s);
        } else {
          ReferenceParams p;
          p.num_states = n;
          p.num_physical = n / 2 + rng.below(n / 2);
          p.num_arcs = 3 * n;
          p.max_len = seed % 2 == 0 ? 4 : 64;
          p.num_sources = 1 + rng.below(3);
          p.extra_physical_edges = seed % 3 == 0 ? n / 4 : 0;
          p.extra_offers = 3;
          inst = reference_instance<U>(p, s);
          inst.name = "growth/" + std::to_string(n);
        }
        worst.max_with(check_counts(inst));
      }
      WARN("Lemma 4.7 growth, " << (base_case ? "base case" : "reference") << ", N = " << n << ": "
                                << worst.str());
      per_n.push_back(worst);
    }
    const Ratios& lo = per_n.front();
    const Ratios& hi = per_n.back();
    CHECK(lo.scans_per_n > 0);  // the hooks saw Scan calls (not a vacuous ratio)
    CHECK(hi.scans_per_n <= 2.0 * std::max(1.0, lo.scans_per_n));
#if SSSP_DIAGNOSTICS
    CHECK(hi.failed_per_n <= 2.0 * std::max(1.0, lo.failed_per_n));
#endif
    CHECK(hi.slot_sum_per_n <= 2.0 * std::max(1.0, lo.slot_sum_per_n));
    CHECK(hi.fail_sum_per_n <= 2.0 * std::max(1.0, lo.fail_sum_per_n));
  }
}
