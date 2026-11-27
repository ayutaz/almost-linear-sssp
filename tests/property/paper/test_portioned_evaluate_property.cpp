// Property tests of the work claims of Lemma 7.7 (Sec 7.4 of arXiv 2611.34911 v1, proof item 2,
// l.1902-1907) for PortionedEvaluator, by counting (no timing):
//   * on the grid of cells [j q, (j + 1) q) with q at least every unit, "each local nominal
//     scan, including scans of persistent children and unsplit runs, occurs at most twice":
//     the Scan calls per nominal slot (node, [a, b)) are <= 2, and the nominal slots scanned are
//     exactly those of the bulk run (where each is scanned once: in Sec 4 alone no clipped slot
//     resumes, l.799-801);
//   * on arbitrary cells a nominal slot is scanned at most once more per cell boundary strictly
//     inside it than in the bulk run (bulk: once with exact arithmetic; with rounded doubles the
//     origin fallback D-19 may already split a slot);
//   * "a cursor resumes from the first live piece, so the charged unit is the same list-piece
//     visit as in Lemma 4.7": run visits and registered-child visits are each at most twice the
//     bulk run's on the same grids. (Failed run visits alone are not bounded that way: a run
//     that registers in one bulk visit can fail in the first portion and register in the
//     second, so 0 failed visits can become 1; the paper charges list-piece visits, l.1904.)
//   * clipped slots really are resumed on these inputs (root and child level).
// With rounded doubles (lengths k / 10) the bound "at most twice" is checked relative to the
// bulk run: at most one more scan per slot than there.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "differential/paper/portioned_cases.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "support/paper/eval_invariants.hpp"

using namespace sssp::test::paper;
using namespace sssp::test::paper::portioned;
using sssp::paper::EvalOptions;

namespace {

constexpr bool kCounted = SSSP_DIAGNOSTICS != 0;

[[maybe_unused]] std::uint64_t counter(const sssp::Diagnostics& d, const std::string& name) {
  const auto it = d.counters.find(name);
  return it == d.counters.end() ? 0 : it->second;
}

template <class W>
std::uint64_t piece_visits(const sssp::Diagnostics& d) {
  return counter(d, "eval_run_visits") + counter(d, "eval_child_visits");
}

template <class W>
std::size_t interior_cuts(const std::vector<W>& cuts, W a, W b) {
  std::size_t k = 0;
  for (const W t : cuts) k += (a < t && t < b) ? 1 : 0;
  return k;
}

}  // namespace

TEMPLATE_TEST_CASE("Lemma 7.7: each local nominal scan occurs at most twice on a grid q >= units",
                   "[portioned][property][lemma7.7]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  std::size_t grids = 0;
  std::size_t twice = 0;
  std::size_t resumed = 0;
  std::size_t resumed_child = 0;
  std::uint64_t counted_resumed = 0;
  for (const auto& c : pool) {
    INFO("case " << c.name);
    const auto root = c.nh.root_order.front();
    ScanTally<W> bulk_tally(root);
    EvalOptions<W> bopt;
    bopt.hooks = &bulk_tally;
    const auto bulk = sssp::paper::evaluate_normalized<W>(c.nh, c.inst.seeded.offers, bopt);
    if (c.integer_valued) {
      // Exact arithmetic: the bulk run never resumes a slot and scans each slot once.
      CHECK(bulk_tally.resumed == 0);
      for (const auto& [slot, n] : bulk_tally.scans) CHECK(n == 1);
    }
    sssp::Diagnostics bulk_diag;
    if constexpr (kCounted) {
      (void)sssp::paper::evaluate_normalized<W>(c.nh, c.inst.seeded.offers, {}, &bulk_diag);
    }
    const W q = max_unit(c);
    W qm{};
    for (int m = 1; m <= 4; ++m) {
      qm = qm + q;
      std::size_t j0 = 0;
      for (W t = qm; t <= least_offer(c); t = t + qm) ++j0;
      for (const std::size_t start : {std::size_t{0}, j0}) {
        const auto cuts = grid_cuts(qm, start, horizon(c));
        INFO("grid q = " << qm << " from cell " << start);
        ScanTally<W> tally(root);
        sssp::Diagnostics diag;
        EvalOptions<W> opt;
        opt.hooks = &tally;
        PortionedEvaluator<W> ev(c.nh, opt, &diag);
        const auto run = run_cells<W>(ev, c.inst.seeded.offers, {}, cuts);
        REQUIRE(run.result.value == bulk.value);
        ++grids;
        // The same nominal slots as the bulk run, each scanned at most twice.
        REQUIRE(tally.scans.size() == bulk_tally.scans.size());
        for (const auto& [slot, n] : tally.scans) {
          REQUIRE(bulk_tally.scans.contains(slot));
          REQUIRE(n <= (c.integer_valued ? 2 : bulk_tally.scans.at(slot) + 1));
          twice += n == 2 ? 1 : 0;
        }
        CHECK(tally.total <= 2 * bulk_tally.total);
        resumed += tally.resumed;
        resumed_child += tally.resumed_child;
        if constexpr (kCounted) {
          // The charged list-piece visits of Lemma 4.7 at most double (l.1903-1904).
          CHECK(counter(diag, "eval_run_visits") <= 2 * counter(bulk_diag, "eval_run_visits"));
          CHECK(counter(diag, "eval_child_visits") <= 2 * counter(bulk_diag, "eval_child_visits"));
          CHECK(piece_visits<W>(diag) <= 2 * piece_visits<W>(bulk_diag));
          CHECK(counter(diag, "eval_cells") == run.cells);
          CHECK(counter(diag, "eval_resumed_slot_visits") +
                    counter(diag, "eval_root_resumed_slot_visits") ==
                tally.resumed);
          counted_resumed += counter(diag, "eval_resumed_slot_visits");
        }
      }
    }
  }
  CAPTURE(pool.size(), grids, twice, resumed, resumed_child, counted_resumed);
  CHECK(twice > 0);
  CHECK(resumed > 0);
  CHECK(resumed_child > 0);
  if constexpr (kCounted) CHECK(counted_resumed == resumed_child);
}

TEMPLATE_TEST_CASE("Lemma 7.7: on arbitrary cells a slot is scanned once per cell it meets",
                   "[portioned][property][lemma7.7]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  std::size_t over_two = 0;
  for (std::size_t i = 0; i < pool.size(); ++i) {
    const auto& c = pool[i];
    INFO("case " << c.name);
    const auto root = c.nh.root_order.front();
    ScanTally<W> bulk_tally(root);
    EvalOptions<W> bopt;
    bopt.hooks = &bulk_tally;
    RecordingObserver<W> rec;
    ForwardingObserver<sssp::paper::EvalHooks<W>, W> adapter(rec);
    EvalOptions<W> ropt;
    ropt.hooks = &adapter;
    const auto bulk = sssp::paper::evaluate_normalized<W>(c.nh, c.inst.seeded.offers, ropt);
    (void)sssp::paper::evaluate_normalized<W>(c.nh, c.inst.seeded.offers, bopt);
    for (const auto& p : partitions(c, rec.events, 13 * i + 1)) {
      INFO("partition " << p.name);
      ScanTally<W> tally(root);
      EvalOptions<W> opt;
      opt.hooks = &tally;
      PortionedEvaluator<W> ev(c.nh, opt);
      const auto run = run_cells<W>(ev, c.inst.seeded.offers, {}, p.cuts);
      REQUIRE(run.result.value == bulk.value);
      REQUIRE(tally.scans.size() == bulk_tally.scans.size());
      for (const auto& [slot, n] : tally.scans) {
        REQUIRE(bulk_tally.scans.contains(slot));
        const auto& [node, a, b] = slot;
        REQUIRE(n <= bulk_tally.scans.at(slot) + interior_cuts(p.cuts, a, b));
        over_two += n > 2 ? 1 : 0;
      }
    }
  }
  CAPTURE(over_two);
  CHECK(over_two > 0);  // narrow cells do cut some slots more than once
}
