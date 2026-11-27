// Differential tests of the Sec 5.5 replacement through Evaluate (Alg 3, Sec 4.3):
//   * before: Evaluate on the fine threshold-SCC hierarchy the partition hierarchy was collapsed
//     from (same underlying state graph: D plus both orientations of every internal edge);
//   * after: Evaluate on substitute_atoms(partition) with the seeds redirected to the native
//     copies (replacement procedure, Lemma 5.14-5.15).
// Lemma 4.6 fixes the minimum emission at every served physical vertex to delta(v); Lemma 5.14
// says replacement keeps every tight witness and creates no underestimate, so the served
// vertices and their per-vertex minima agree before and after, and every emission stays at
// least delta. On undirected inputs with the identity projection every vertex is served, and
// the per-vertex minima equal a Dijkstra on the physical graph written in the test support.
// Both branches of Alg 4 are exercised: the base case (Mode A parameters) and the radial
// construction with Mode B parameters (threshold k > L_k^2, separator b; atoms of >= 37
// states, N >= 32 and b + 1 <= N as Sec 5.4 l.1206 requires).
// Weights: std::uint64_t and quarter-integer doubles (exact sums, exact comparison).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/local_hierarchy.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/prng.hpp"
#include "unit/paper/local/atom_test_support.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::Vertex;
using sssp::WeightTraits;
using sssp::paper::LocalHierarchyParams;
using sssp::paper::SeededOrderedHierarchy;
namespace at = sssp::test::atom;
namespace lt = sssp::test::local;
namespace tp = sssp::test::paper;

namespace {

struct Seen {
  std::size_t instances = 0;
  std::size_t full_dijkstra = 0;
  std::size_t radial = 0;
};

template <class W>
void compare(const at::PartitionInstance<W>& inst, const LocalHierarchyParams& lp,
             std::uint64_t seed, Seen& seen) {
  const auto sub = sssp::paper::substitute_atoms(inst.partition, nullptr, lp);
  const std::size_t np = inst.num_physical;
  ++seen.instances;
  for (const auto& a : inst.partition.atoms) {
    seen.radial += sssp::paper::local_construction_applies(a.states.size(), lp) ? 1 : 0;
  }
  sssp::test::SplitMix64 rng(seed * 7919 + 11);
  for (int round = 0; round < 3; ++round) {
    INFO("round " << round);
    auto [sources, offers] = at::make_seeds(inst, 1 + static_cast<std::size_t>(round % 2), rng);
    const SeededOrderedHierarchy<W> before{inst.fine, sources, offers};
    const SeededOrderedHierarchy<W> after{sub.hierarchy, sources,
                                          sssp::paper::redirect_offers<W>(sub, offers)};
    {
      const auto v1 = sssp::paper::validate_seeded_ordered_hierarchy(before, inst.physical);
      INFO(v1.summary());
      REQUIRE(v1.ok());
      const auto v2 = sssp::paper::validate_seeded_ordered_hierarchy(after, inst.physical);
      INFO(v2.summary());
      REQUIRE(v2.ok());
    }
    const auto ref_before = tp::compute_eval_reference(before, inst.physical);
    const auto ref_after = tp::compute_eval_reference(after, inst.physical);
    const auto served_before = tp::served_targets<W>(np, inst.fine.projection, ref_before);
    const auto served_after = tp::served_targets<W>(np, sub.hierarchy.projection, ref_after);
    const auto r_before = sssp::paper::evaluate(before);
    const auto r_after = sssp::paper::evaluate(after);
    const auto m_before = tp::physical_minima<W>(np, inst.fine.projection, r_before.emissions);
    const auto m_after = tp::physical_minima<W>(np, sub.hierarchy.projection, r_after.emissions);
    for (Vertex v = 0; v < np; ++v) {
      INFO("vertex " << v);
      CHECK(served_before[v] == served_after[v]);
      if (served_before[v]) {
        CHECK(m_before[v] == ref_before.target[v]);
        CHECK(m_after[v] == m_before[v]);
      }
    }
    for (const auto& [x, val] : r_after.emissions) {
      CHECK(!(val < ref_after.target[sub.hierarchy.projection[x]]));
    }
    if (inst.all_undirected && !inst.partition.atoms.empty() && np == inst.full_graph.num_states) {
      std::vector<std::pair<Vertex, W>> src;
      for (const auto& s : sources) src.emplace_back(s.vertex, s.offset);
      std::vector<sssp::Edge<W>> edges(inst.physical.edges().begin(), inst.physical.edges().end());
      const auto want = lt::dijkstra_naive<W>(np, edges, src);
      for (Vertex v = 0; v < np; ++v) {
        INFO("dijkstra vertex " << v);
        if (want[v] == WeightTraits<W>::infinity()) {
          CHECK(m_after[v] == WeightTraits<W>::infinity());
        } else {
          CHECK(served_after[v]);
          CHECK(m_after[v] == want[v]);
        }
      }
      ++seen.full_dijkstra;
    }
  }
}

}  // namespace

TEMPLATE_TEST_CASE("atom substitution keeps Evaluate's per-vertex minima (base-case atoms)", "",
                   std::uint64_t, double) {
  using W = TestType;
  Seen seen;
  for (const bool one_way : {false, true}) {
    for (const bool fold : {false, true}) {
      for (const bool connect : {true, false}) {
        for (std::uint64_t seed = 1; seed <= 10; ++seed) {
          at::PartitionParams p;
          p.num_states = 32;
          p.one_way = one_way ? 10 : 0;
          p.fold_projection = fold;
          p.connect = connect;
          p.singleton_percent = 20;
          p.internal_max = seed % 2 == 0 ? 3 : 30;
          p.outer_min = 1;
          p.outer_max = seed % 3 == 0 ? 60 : 20;
          INFO("one_way=" << one_way << " fold=" << fold << " connect=" << connect
                          << " seed=" << seed);
          compare(at::make_partition_instance<W>(p, seed), sssp::paper::kPaperLocalParams, seed,
                  seen);
        }
      }
    }
  }
  CHECK(seen.full_dijkstra > 0);
}

TEMPLATE_TEST_CASE("atom substitution keeps Evaluate's per-vertex minima (radial atoms)", "",
                   std::uint64_t, double) {
  using W = TestType;
  Seen seen;
  for (const bool one_way : {false, true}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t b : {std::size_t{2}, std::size_t{6}}) {
        for (std::uint64_t seed = 1; seed <= 2; ++seed) {
          at::PartitionParams p;
          p.num_states = 90;
          p.min_group = 37;
          p.max_group = 45;
          p.internal_min = 1;
          p.internal_max = wd == lt::Weights::kTies ? 3 : 1000;
          p.internal_exponential = wd == lt::Weights::kExponential;
          p.outer_min = wd == lt::Weights::kExponential ? (std::uint64_t{1} << 42) : 2000;
          p.outer_max = wd == lt::Weights::kExponential ? (std::uint64_t{1} << 43) : 5000;
          p.one_way = one_way ? 8 : 0;
          p.largest_first = true;
          p.atom_percent = 100;
          INFO("one_way=" << one_way << " " << lt::name(wd) << " b=" << b << " seed=" << seed);
          compare(at::make_partition_instance<W>(p, seed), LocalHierarchyParams{1, 2, b}, seed,
                  seen);
        }
      }
    }
  }
  CHECK(seen.radial > 0);
  CHECK(seen.full_dijkstra > 0);
}
