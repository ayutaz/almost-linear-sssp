#pragma once
// Test-only helpers for the PortionedEvaluator tests (Lemma 7.7, Sec 7.4 of arXiv 2611.34911
// v1): a pool of seeded ordered hierarchies with one root, cell partitions of the time line,
// a runner that visits the cells in order (posting external offers just before the cell of
// their timestamp), and the comparison with the bulk evaluate (Alg 3).
//
// Why the bulk run is the reference. The bulk evaluate is checked against Dijkstra and the
// I1-I8 checker in test_evaluate_vs_dijkstra.cpp. Cutting the time line at cell boundaries
// changes only the order in which Alg 2 meets the portions of a nominal slot: within one
// portion the cursor still walks the list in rank order, and a key that becomes due in a
// later portion was not due in the earlier one (arcs into lower ranks are >= the nominal slot
// end, supplemental offers are padded to it, Lemma 7.7 proof item 1). So the portioned run
// emits exactly the bulk (state, value) pairs, and its emission sequence is the bulk sequence
// stably sorted by the cell containing each value. The tests below assert exactly that.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/primitives/split_findmin.hpp"
#include "sssp/graph.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/seeded_instances.hpp"
#include "support/prng.hpp"

namespace sssp::test::paper::portioned {

using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::EvalResult;
using sssp::paper::NormalizedHierarchy;
using sssp::paper::PortionedEvaluator;
using sssp::paper::StateId;
using sssp::paper::StateOffer;

template <class W>
W inf() {
  return WeightTraits<W>::infinity();
}

template <class W>
NormalizedHierarchy<W> prepare(const sssp::paper::SeededOrderedHierarchy<W>& s) {
  NormalizedHierarchy<W> nh = sssp::paper::normalize<W>(s.hierarchy);
  sssp::paper::order_children<W>(nh);
  return nh;
}

template <class W>
struct Case {
  std::string name;
  SeededInstance<W> inst;
  NormalizedHierarchy<W> nh;
  bool integer_valued = true;  ///< every length and offer is an integer (exact in double)
};

/// A strongly connected random state graph (a random Hamiltonian cycle plus random arcs), so
/// its threshold-SCC hierarchy has one root. `tenths`: lengths k / 10 (double only), which
/// are not exact in binary.
template <class W>
SeededInstance<W> strongly_connected_instance(std::uint64_t seed, bool tenths = false) {
  SplitMix64 rng(seed);
  const std::size_t n = 2 + rng.below(28);
  const std::size_t phys = 1 + rng.below(n);
  const std::uint64_t max_len = rng.below(2) == 0 ? 3 : 25;
  auto len = [&]() -> W {
    const std::uint64_t k = rng.between(1, max_len);
    if constexpr (std::is_floating_point_v<W>) {
      if (tenths) return static_cast<W>(k) / W{10};
    }
    return static_cast<W>(k);
  };
  std::vector<StateId> perm(n);
  for (std::size_t i = 0; i < n; ++i) perm[i] = static_cast<StateId>(i);
  for (std::size_t i = n; i > 1; --i) std::swap(perm[i - 1], perm[rng.below(i)]);
  sssp::paper::StateGraph<W> g;
  g.num_states = n;
  for (std::size_t i = 0; i < n; ++i) g.arcs.push_back({perm[i], perm[(i + 1) % n], len()});
  const std::size_t extra = rng.below(3 * n + 1);
  for (std::size_t i = 0; i < extra; ++i) {
    g.arcs.push_back(
        {static_cast<StateId>(rng.below(n)), static_cast<StateId>(rng.below(n)), len()});
  }
  auto proj = detail::random_projection<W>(n, phys, rng);
  SeededInstance<W> inst;
  inst.name = std::string(tenths ? "tenths/" : "strong/") + std::to_string(seed);
  inst.seed = seed;
  inst.num_physical = phys;
  const CsrGraph<W> base = sssp::paper::physical_graph_of_arcs<W>(g, proj, phys);
  std::vector<Edge<W>> edges(base.edges().begin(), base.edges().end());
  const std::size_t shortcuts = seed % 3 == 0 ? rng.below(phys + 1) : 0;
  for (std::size_t i = 0; i < shortcuts; ++i) {
    edges.push_back(
        {static_cast<Vertex>(rng.below(phys)), static_cast<Vertex>(rng.below(phys)), len()});
  }
  inst.physical = CsrGraph<W>::from_edges(phys, std::move(edges));
  inst.seeded.hierarchy = sssp::paper::build_threshold_scc_hierarchy<W>(g, proj);
  const std::size_t sources = seed % 2 == 0 ? 1 : 1 + rng.below(3);
  detail::seed_offers(inst, sources, /*max_offset=*/12, /*extra_offers=*/rng.below(4), max_len,
                      rng);
  return inst;
}

/// Every seeded hierarchy of the test-support catalogues with exactly one root (Lemma 7.7),
/// plus strongly connected instances (and, for double, instances with lengths k / 10).
template <class W>
std::vector<Case<W>> one_root_pool() {
  std::vector<Case<W>> out;
  auto take = [&](SeededInstance<W> inst, bool integer_valued) {
    if (inst.seeded.hierarchy.num_states() == 0 || inst.seeded.offers.empty()) return;
    NormalizedHierarchy<W> nh = prepare<W>(inst.seeded);
    if (nh.root_order.size() != 1) return;
    Case<W> c;
    c.name = inst.name + "#" + std::to_string(inst.seed);
    c.inst = std::move(inst);
    c.nh = std::move(nh);
    c.integer_valued = integer_valued;
    out.push_back(std::move(c));
  };
  for (const std::uint64_t seed : {11ULL, 12ULL}) {
    for (auto& inst : base_case_catalogue<W>(seed)) take(std::move(inst), true);
  }
  for (auto& inst : reference_catalogue<W>(5000, 120)) take(std::move(inst), true);
  for (std::uint64_t k = 0; k < 40; ++k) {
    SplitMix64 rng(7000 + k);
    ReferenceParams p;
    p.num_states = 2 + rng.below(20);
    p.num_physical = 1 + rng.below(p.num_states);
    p.num_arcs = p.num_states + rng.below(3 * p.num_states + 1);
    p.max_len = rng.below(2) == 0 ? 2 : 12;
    p.num_sources = 1 + rng.below(3);
    p.extra_offers = rng.below(4);
    p.zero_arcs = 1 + rng.below(p.num_states);
    take(zero_arc_instance<W>(p, 7000 + k), true);
  }
  for (std::uint64_t k = 0; k < 120; ++k) take(strongly_connected_instance<W>(9000 + k), true);
  if constexpr (std::is_floating_point_v<W>) {
    for (std::uint64_t k = 0; k < 60; ++k) {
      take(strongly_connected_instance<W>(11000 + k, /*tenths=*/true), false);
    }
  }
  return out;
}

// ------------------------------------------------------------------ partitions

/// Cells [cuts[i], cuts[i + 1]) and the last cell [cuts.back(), +inf).
template <class W>
struct Partition {
  std::string name;
  std::vector<W> cuts;
};

template <class W>
std::size_t cell_of(W v, const std::vector<W>& cuts) {
  return static_cast<std::size_t>(std::upper_bound(cuts.begin(), cuts.end(), v) - cuts.begin()) - 1;
}

template <class W>
void normalise_cuts(std::vector<W>& cuts) {
  std::erase_if(cuts, [](W v) { return v == inf<W>(); });
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
}

template <class W>
W least_offer(const Case<W>& c) {
  W m = inf<W>();
  for (const auto& o : c.inst.seeded.offers) m = std::min(m, o.value);
  return m;
}

/// Every emission is at most D0 + S(root) <= max offer + S(root) (Alg 3 l.9).
template <class W>
W horizon(const Case<W>& c) {
  W m{};
  for (const auto& o : c.inst.seeded.offers) m = std::max(m, o.value);
  const auto r = c.nh.root_order.front();
  if (!c.nh.forest.is_leaf(r)) m = m + c.nh.mass[r];
  return m + W{1};
}

/// The largest unit of the hierarchy (Lemma 7.7: "units at most q").
template <class W>
W max_unit(const Case<W>& c) {
  W q{};
  for (sssp::paper::NodeId v = 0; v < c.nh.num_nodes(); ++v) {
    if (!c.nh.forest.is_leaf(v)) q = std::max(q, c.nh.unit[v]);
  }
  return q == W{} ? W{1} : q;
}

template <class W>
W random_point(SplitMix64& rng, W hi) {
  if constexpr (std::is_floating_point_v<W>) {
    return hi * static_cast<W>(rng.below(1U << 20)) / static_cast<W>(1U << 20);
  } else {
    return hi <= 1 ? W{0} : static_cast<W>(rng.below(hi));
  }
}

/// The grid [j q, (j + 1) q) of Lemma 7.7 from cell j0 on, built by successive additions.
template <class W>
std::vector<W> grid_cuts(W q, std::size_t j0, W until) {
  std::vector<W> cuts;
  W t{};
  for (std::size_t j = 0; j < j0; ++j) t = t + q;
  for (; t <= until; t = t + q) cuts.push_back(t);
  cuts.push_back(t);
  return cuts;
}

/// Cut points taken from a bulk trace: every scan and slot boundary, deadline, origin, D0
/// and emitted value (cells starting exactly at the events of the bulk run).
template <class W>
std::vector<W> boundary_cuts(const std::vector<EvalEvent<W>>& events) {
  std::vector<W> cuts{W{0}};
  for (const auto& e : events) {
    switch (e.kind) {
      case EvalEventKind::kScanBegin:
        cuts.insert(cuts.end(), {e.lo, e.hi, e.slot_lo, e.slot_hi, e.deadline});
        if constexpr (!std::is_floating_point_v<W>) cuts.push_back(e.deadline + 1);
        break;
      case EvalEventKind::kOpen:
        cuts.insert(cuts.end(), {e.value, e.origin});
        break;
      case EvalEventKind::kEmit:
        cuts.push_back(e.value);
        break;
      default:
        break;
    }
  }
  normalise_cuts(cuts);
  return cuts;
}

/// Partitions for one case: one cell; random cuts (few and many); cells at every boundary
/// of the bulk run; unit cells (integers); the Lemma 7.7 grids with q = 1, 2, 3 times the
/// largest unit, from 0 and from the cell of the least offer; random cuts after a first cell
/// starting exactly at the least offer.
template <class W>
std::vector<Partition<W>> partitions(const Case<W>& c, const std::vector<EvalEvent<W>>& bulk,
                                     std::uint64_t seed) {
  SplitMix64 rng(seed);
  const W b = horizon(c);
  std::vector<Partition<W>> out;
  out.push_back({"single", {W{0}}});
  for (const std::size_t k : {std::size_t{1}, std::size_t{4}, std::size_t{16}, std::size_t{64}}) {
    std::vector<W> cuts{W{0}};
    for (std::size_t i = 0; i < k; ++i) cuts.push_back(random_point(rng, b));
    normalise_cuts(cuts);
    out.push_back({"random" + std::to_string(k), cuts});
  }
  out.push_back({"boundaries", boundary_cuts(bulk)});
  if constexpr (!std::is_floating_point_v<W>) {
    if (b <= 4000) {
      std::vector<W> cuts;
      for (W t = 0; t <= b; ++t) cuts.push_back(t);
      out.push_back({"unit", cuts});
    }
  }
  const W q = max_unit(c);
  W qm{};
  for (int m = 1; m <= 3; ++m) {
    qm = qm + q;
    out.push_back({"grid" + std::to_string(m), grid_cuts(qm, 0, b)});
    std::size_t j0 = 0;
    for (W t = qm; t <= least_offer(c); t = t + qm) ++j0;
    out.push_back({"grid" + std::to_string(m) + "-late", grid_cuts(qm, j0, b)});
  }
  {
    const W lo = least_offer(c);
    std::vector<W> cuts{lo};
    for (int i = 0; i < 8; ++i) {
      const W p = random_point(rng, b);
      if (lo < p) cuts.push_back(p);
    }
    normalise_cuts(cuts);
    out.push_back({"from-least-offer", cuts});
  }
  return out;
}

// ------------------------------------------------------------------ running

template <class W>
struct Run {
  EvalResult<W> result;
  std::vector<std::vector<std::pair<StateId, W>>> per_cell;
  std::size_t cells = 0;
};

/// Visits the cells of `cuts` in order until the evaluator finishes. `external` offers are
/// posted just before the visit of the cell containing their value (Lemma 7.7: "arrives
/// before the visit to its timestamp's cell"); `initial` offers before the first visit.
template <class W>
Run<W> run_cells(PortionedEvaluator<W>& ev, std::span<const StateOffer<W>> initial,
                 std::span<const StateOffer<W>> external, const std::vector<W>& cuts) {
  for (const auto& o : initial) ev.initial_offer(o.state, o.value);
  std::vector<std::vector<StateOffer<W>>> due(cuts.size());
  for (const auto& o : external) due[cell_of(o.value, cuts)].push_back(o);
  Run<W> r;
  for (std::size_t k = 0; k < cuts.size() && !ev.finished(); ++k) {
    for (const auto& o : due[k]) (void)ev.external_offer(o.state, o.value);
    const W hi = k + 1 < cuts.size() ? cuts[k + 1] : inf<W>();
    const auto got = ev.visit_cell(cuts[k], hi);
    r.per_cell.emplace_back(got.begin(), got.end());
    ++r.cells;
  }
  REQUIRE(ev.finished());
  r.result = ev.take_result();
  return r;
}

template <class W>
Run<W> run_portioned(const NormalizedHierarchy<W>& nh, std::span<const StateOffer<W>> initial,
                     std::span<const StateOffer<W>> external, const std::vector<W>& cuts,
                     const EvalOptions<W>& opt = {}, W deadline = inf<W>()) {
  PortionedEvaluator<W> ev(nh, opt);
  if (deadline != inf<W>()) ev.set_deadline(deadline);
  return run_cells(ev, initial, external, cuts);
}

/// The bulk emission sequence stably sorted by the cell containing each value.
template <class W>
std::vector<std::pair<StateId, W>> cell_stable(std::vector<std::pair<StateId, W>> e,
                                               const std::vector<W>& cuts) {
  std::stable_sort(e.begin(), e.end(), [&](const auto& x, const auto& y) {
    return cell_of(x.second, cuts) < cell_of(y.second, cuts);
  });
  return e;
}

/// value and emissions equal the bulk run's (emissions in cell-stable order); every emission
/// of a cell lies in that cell.
template <class W>
void require_matches_bulk(const EvalResult<W>& bulk, const Run<W>& run,
                          const std::vector<W>& cuts) {
  REQUIRE(run.result.value == bulk.value);
  REQUIRE(run.result.emissions == cell_stable(bulk.emissions, cuts));
  for (std::size_t k = 0; k < run.per_cell.size(); ++k) {
    for (const auto& [x, v] : run.per_cell[k]) {
      INFO("state " << x << " value " << v << " cell " << k);
      REQUIRE(cell_of(v, cuts) == k);
    }
  }
}

/// Counts Scan calls per nominal slot (node, [a, b)) and resumed portions (a' > a).
template <class W>
class ScanTally final : public EvalHooks<W> {
 public:
  explicit ScanTally(sssp::paper::NodeId root) : root_(root) {}
  std::map<std::tuple<std::uint32_t, W, W>, std::size_t> scans;
  std::size_t total = 0;
  std::size_t resumed = 0;
  std::size_t resumed_child = 0;

  void on_scan_begin(std::uint32_t node, std::span<const StateId> /*states*/, W lo, W /*hi*/,
                     W slot_lo, W slot_hi, W /*deadline*/) override {
    ++scans[{node, slot_lo, slot_hi}];
    ++total;
    if (slot_lo < lo) {
      ++resumed;
      if (node != root_) ++resumed_child;
    }
  }

 private:
  sssp::paper::NodeId root_;
};

}  // namespace sssp::test::paper::portioned
