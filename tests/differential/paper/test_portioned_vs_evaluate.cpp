// Differential tests: PortionedEvaluator (Lemma 7.7, Sec 7.4) against the bulk evaluate (Alg 3,
// Sec 4.3) of arXiv 2611.34911 v1 on every one-root seeded hierarchy of the test catalogues.
// See portioned_cases.hpp for why the bulk run is the reference and what is compared.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "differential/paper/portioned_cases.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/primitives/split_findmin.hpp"
#include "support/paper/eval_invariants.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/prng.hpp"

using namespace sssp::test::paper;
using namespace sssp::test::paper::portioned;
using sssp::paper::EvalHooks;
using sssp::paper::EvalOptions;
using sssp::paper::StateOffer;
using sssp::test::SplitMix64;

namespace {

template <class W>
struct Bulk {
  sssp::paper::EvalResult<W> result;
  std::vector<EvalEvent<W>> events;
};

template <class W>
Bulk<W> bulk_run(const Case<W>& c, bool supplemental = true) {
  RecordingObserver<W> rec;
  ForwardingObserver<EvalHooks<W>, W> adapter(rec);
  EvalOptions<W> opt;
  opt.supplemental_offers = supplemental;
  opt.hooks = &adapter;
  Bulk<W> b{sssp::paper::evaluate_normalized<W>(c.nh, c.inst.seeded.offers, opt), {}};
  b.events = std::move(rec.events);
  return b;
}

}  // namespace

TEMPLATE_TEST_CASE("portioned vs evaluate: random and boundary partitions give the bulk emissions",
                   "[portioned][differential][lemma7.7]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  REQUIRE(pool.size() >= 150);
  std::size_t runs = 0;
  std::size_t checked = 0;
  std::size_t emissions = 0;
  for (std::size_t i = 0; i < pool.size(); ++i) {
    const auto& c = pool[i];
    INFO("case " << c.name);
    const auto bulk = bulk_run(c);
    const std::span<const StateOffer<W>> offers = c.inst.seeded.offers;
    for (const auto& p : partitions(c, bulk.events, 31 * i + 7)) {
      INFO("partition " << p.name << " with " << p.cuts.size() << " cuts");
      if (c.integer_valued) {
        // The same run under the Sec 4 invariant checker (I1-I8 and the hook protocol).
        const auto ref = compute_eval_reference<W>(c.inst.seeded, c.inst.physical);
        InvariantChecker<W> checker(c.inst.seeded, ref);
        ForwardingObserver<EvalHooks<W>, W> adapter(checker);
        EvalOptions<W> opt;
        opt.hooks = &adapter;
        const auto run = run_portioned<W>(c.nh, offers, {}, p.cuts, opt);
        require_matches_bulk(bulk.result, run, p.cuts);
        INFO(checker.report().summary());
        REQUIRE(checker.report().ok());
        ++checked;
      } else {
        const auto run = run_portioned<W>(c.nh, offers, {}, p.cuts);
        require_matches_bulk(bulk.result, run, p.cuts);
      }
      ++runs;
      emissions += bulk.result.emissions.size();
    }
  }
  CAPTURE(pool.size(), runs, checked, emissions);
  CHECK(runs >= 10 * pool.size());
}

TEMPLATE_TEST_CASE("portioned vs evaluate: supplemental offers off (D-14) as well",
                   "[portioned][differential]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  for (std::size_t i = 0; i < pool.size(); ++i) {
    const auto& c = pool[i];
    INFO("case " << c.name);
    const auto bulk = bulk_run(c, /*supplemental=*/false);
    EvalOptions<W> opt;
    opt.supplemental_offers = false;
    for (const auto& p : partitions(c, bulk.events, 17 * i + 3)) {
      INFO("partition " << p.name);
      const auto run = run_portioned<W>(c.nh, c.inst.seeded.offers, {}, p.cuts, opt);
      require_matches_bulk(bulk.result, run, p.cuts);
    }
  }
}

TEMPLATE_TEST_CASE("portioned vs evaluate: external offers that arrive before their cell",
                   "[portioned][differential][lemma7.7]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  std::size_t external = 0;
  for (std::size_t i = 0; i < pool.size(); ++i) {
    const auto& c = pool[i];
    INFO("case " << c.name);
    const auto bulk = bulk_run(c);
    SplitMix64 rng(1000 + i);
    for (const auto& p : partitions(c, bulk.events, 5 * i + 1)) {
      if (least_offer(c) < p.cuts.front()) continue;
      INFO("partition " << p.name);
      // Every offer but the least becomes external with probability 1/2 (the least stays
      // initial so that the first cell starts at or before it either way).
      std::vector<StateOffer<W>> ini, ext;
      const W lo = least_offer(c);
      bool kept_least = false;
      for (const auto& o : c.inst.seeded.offers) {
        const bool least = !kept_least && o.value == lo;
        kept_least = kept_least || least;
        (least || rng.below(2) == 0 ? ini : ext).push_back(o);
      }
      external += ext.size();
      const auto run = run_portioned<W>(c.nh, ini, ext, p.cuts);
      require_matches_bulk(bulk.result, run, p.cuts);
      if constexpr (SSSP_DIAGNOSTICS != 0) {
        const auto it = run.result.diag.counters.find("eval_external_offers");
        const std::size_t e = it == run.result.diag.counters.end() ? 0 : it->second;
        CHECK(e <= ext.size());  // offers after the evaluator finished are not posted
      }
    }
  }
  CAPTURE(external);
  CHECK(external > 0);
}

TEMPLATE_TEST_CASE("portioned vs evaluate: an inclusive deadline keeps exactly the emissions <= it",
                   "[portioned][differential]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  std::size_t cut_off = 0;
  for (std::size_t i = 0; i < pool.size(); ++i) {
    const auto& c = pool[i];
    INFO("case " << c.name);
    const auto bulk = bulk_run(c);
    SplitMix64 rng(2000 + i);
    for (const auto& p : partitions(c, bulk.events, 11 * i + 5)) {
      INFO("partition " << p.name);
      const W deadline = random_point(rng, horizon(c));
      CAPTURE(deadline);
      sssp::paper::EvalResult<W> want;
      want.value.assign(bulk.result.value.size(), inf<W>());
      for (const auto& [x, v] : bulk.result.emissions) {
        if (deadline < v) {
          ++cut_off;
          continue;
        }
        want.emissions.emplace_back(x, v);
        want.value[x] = v;
      }
      const auto run = run_portioned<W>(c.nh, c.inst.seeded.offers, {}, p.cuts, {}, deadline);
      require_matches_bulk(want, run, p.cuts);
    }
  }
  CAPTURE(cut_off);
  CHECK(cut_off > 0);
}

TEMPLATE_TEST_CASE("portioned vs evaluate: atoms side by side in one split-findmin layout",
                   "[portioned][differential]", std::uint64_t, double) {
  using W = TestType;
  const auto pool = one_root_pool<W>();
  // Groups of up to five atoms laid out consecutively (Alg 7 l.8), visited cell by cell in
  // layout order with the same partition, each compared with its own bulk run.
  for (std::size_t g = 0; g < pool.size(); g += 5) {
    const std::size_t end = std::min(pool.size(), g + 5);
    std::size_t total = 0;
    W b{};
    for (std::size_t i = g; i < end; ++i) {
      total += pool[i].nh.forest.num_states();
      b = std::max(b, horizon(pool[i]));
    }
    SplitMix64 rng(3000 + g);
    std::vector<W> cuts{W{0}};
    for (int k = 0; k < 12; ++k) cuts.push_back(random_point(rng, b));
    normalise_cuts(cuts);
    sssp::paper::SplitFindmin<W> sf(total);
    std::vector<PortionedEvaluator<W>> evs;
    std::uint32_t base = 0;
    for (std::size_t i = g; i < end; ++i) {
      evs.emplace_back(pool[i].nh, sf, base);
      base += static_cast<std::uint32_t>(pool[i].nh.forest.num_states());
    }
    for (std::size_t i = g; i < end; ++i) {
      for (const auto& o : pool[i].inst.seeded.offers) {
        evs[i - g].initial_offer(o.state, o.value);
      }
    }
    std::vector<std::vector<std::pair<sssp::paper::StateId, W>>> got(end - g);
    for (std::size_t k = 0; k < cuts.size(); ++k) {
      const W hi = k + 1 < cuts.size() ? cuts[k + 1] : inf<W>();
      for (std::size_t j = 0; j < evs.size(); ++j) {
        if (evs[j].finished()) continue;
        for (const auto& e : evs[j].visit_cell(cuts[k], hi)) got[j].push_back(e);
      }
    }
    for (std::size_t i = g; i < end; ++i) {
      INFO("case " << pool[i].name);
      REQUIRE(evs[i - g].finished());
      const auto bulk = bulk_run(pool[i]);
      CHECK(evs[i - g].value() == bulk.result.value);
      CHECK(got[i - g] == cell_stable(bulk.result.emissions, cuts));
    }
  }
}
