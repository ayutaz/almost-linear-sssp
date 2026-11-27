// Property tests of the Mode B assembly engine (task P5-EXP-LOCAL) on seeded graph families,
// with oracles written in the test support (Kosaraju SCCs, brute-force sums):
//   * with the paper's constants (Sec 5.4: c = 10, support N^16 sigma) Mode A and Mode B produce
//     the assemblies and counters recorded before the two implementations were unified
//     (P5-LOCAL-UNIFY, Q-59), pinned as hashes of every field; D-8 and D-11 likewise;
//   * D-8 (research note 02 §4.3, Lambda = 2N * 3hK^2, support end 2(3hK^2)^2 N^3 sigma):
//     Def 4.1 with c0 = 2 and the threshold-SCC characterisation of the proof of Lemma 5.10
//     (l.1271-1275) for every N >= 2 and b in {1, 2, 3, 4, ceil(log2 N), N, N + 1}; the scale
//     relations the D-8 proof sketch uses: an unaligned block has total block-MST length below
//     sigma_U / 2 (the replaced l.1208), every finite theta and every block-MST edge of an
//     aligned block lies in the support of its key (the replaced l.1216-1221);
//   * D-11 (exact-height assembly, Q-39): Def 4.1 with c0 = 1, every height a join weight,
//     every node an SCC of the arcs <= its height and its children SCCs of the arcs < it;
//   * grid against exact: every internal node of a grid hierarchy has the state set of an
//     internal node of the exact-height hierarchy, whose height H lies in [g, upper) (at most
//     upper for an inclusive bound): the grid only coarsens the exact heights by less than a
//     factor of two.

#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <type_traits>
#include <vector>

#include "experimental/paper_experimental/local_assembly.hpp"
#include "experimental/paper_experimental/params.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/assembly_core.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/threshold_joins.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper_experimental/local_mode_b_support.hpp"
#include "unit/paper/local/assembly_hash.hpp"
#include "unit/paper/local/clustered_graphs.hpp"
#include "unit/paper/local/local_test_support.hpp"

using sssp::CsrGraph;
using sssp::Diagnostics;
using sssp::WeightTraits;
using sssp::paper::BlockScaleKind;
using sssp::paper::JoinKind;
using sssp::paper::NodeId;
using sssp::paper::StateId;
namespace pe = sssp::paper_experimental;
namespace lt = sssp::test::local;
namespace mb = sssp::test::mode_b;

namespace {

template <class W>
struct Instance {
  sssp::paper::OwnerDecomposition<W> od;
  sssp::paper::RadialStateGraph<W> rs;
  sssp::paper::ThresholdJoins<W> joins;
};

template <class W>
Instance<W> build(const CsrGraph<W>& g, std::uint64_t b) {
  auto od = sssp::paper::owner_decomposition(g, static_cast<std::size_t>(b));
  auto rs = sssp::paper::local_radial_states(g, od);
  auto joins = sssp::paper::compute_joins(od, rs);
  return {std::move(od), std::move(rs), std::move(joins)};
}

struct Coverage {
  std::size_t aligned = 0, unaligned = 0, unattached = 0, separate_nodes = 0, global_nodes = 0;
};

/// D-8 scale relations (research note 02 §4.3) on a grid assembly with the Lambda formulas.
template <class W>
void check_lambda_scales(const Instance<W>& in, const sssp::paper::LocalAssembly<W>& a,
                         Coverage& cov) {
  const std::size_t n = in.od.num_vertices();
  const std::uint64_t b = in.od.separator;
  const unsigned k = in.rs.shell_ratio_log2;
  const auto& jr = in.joins.joins;
  const auto end_of = [&](std::uint32_t key) {
    return pe::lambda_support_end<W>(in.od.corridors[key].weight, n, b, k);
  };
  for (std::size_t u = 0; u < a.blocks.size(); ++u) {
    const auto& bs = a.blocks[u];
    W total{0};
    W heaviest{0};
    for (const auto& j : jr) {
      if (j.kind != JoinKind::kBlockMst || in.rs.block[j.a] != u) continue;
      total = total + j.weight;
      heaviest = std::max(heaviest, j.weight);
    }
    switch (bs.kind) {
      case BlockScaleKind::kAligned: {
        ++cov.aligned;
        const W sigma = in.od.corridors[bs.sigma_key].weight;
        CHECK(!(pe::lambda_scaled<W>(bs.s, n, b, k) < sigma));
        CHECK(!(end_of(bs.sigma_key) < heaviest));
        break;
      }
      case BlockScaleKind::kUnaligned: {
        ++cov.unaligned;
        const W sigma = in.od.corridors[bs.sigma_key].weight;
        CHECK(pe::lambda_scaled<W>(bs.s, n, b, k) < sigma);
        CHECK(total + total < sigma);  // total block-MST length < sigma_U / 2
        break;
      }
      case BlockScaleKind::kUnattached:
        ++cov.unattached;
        break;
      case BlockScaleKind::kSingleton:
        break;
    }
  }
  for (StateId z = 0; z < in.joins.num_states; ++z) {
    if (in.joins.theta[z] == WeightTraits<W>::infinity()) continue;
    CHECK(!(end_of(in.joins.b0_key[z]) < in.joins.theta[z]));
  }
  for (std::uint32_t c = 0; c < a.supports.size(); ++c) CHECK(a.supports[c].upper == end_of(c));
}

/// Every internal node of `grid` has the state set of an internal node of `exact` with height
/// H in [g, upper) (or [g, upper] for an inclusive bound).
template <class W>
void check_grid_coarsens_exact(const sssp::paper::LocalAssembly<W>& grid,
                               const sssp::paper::LocalAssembly<W>& exact, Coverage& cov) {
  const auto& gh = grid.result.hierarchy;
  const auto& eh = exact.result.hierarchy;
  std::map<std::vector<StateId>, W> exact_nodes;
  for (std::size_t v = 0; v < eh.nodes.size(); ++v) {
    if (eh.is_leaf(static_cast<NodeId>(v))) continue;
    exact_nodes.emplace(mb::states_below(eh, static_cast<NodeId>(v)), eh.nodes[v].height);
  }
  for (std::size_t v = 0; v < gh.nodes.size(); ++v) {
    if (gh.is_leaf(static_cast<NodeId>(v))) continue;
    if (grid.node_origin[v] == sssp::paper::NodeOrigin::kSeparateGrid) ++cov.separate_nodes;
    if (grid.node_origin[v] == sssp::paper::NodeOrigin::kGlobalGrid) ++cov.global_nodes;
    const auto it = exact_nodes.find(mb::states_below(gh, static_cast<NodeId>(v)));
    INFO("grid node " << v << " height " << gh.nodes[v].height);
    REQUIRE(it != exact_nodes.end());
    const W g = gh.nodes[v].height;
    const W upper = grid.node_upper[v];
    CHECK(!(it->second < g));
    if (grid.node_upper_inclusive[v] != 0) {
      CHECK(!(upper < it->second));
    } else {
      CHECK(it->second < upper);
    }
  }
}

template <class W>
void check_exact(const CsrGraph<W>& g, const Instance<W>& in,
                 const sssp::paper::LocalAssembly<W>& e) {
  const auto& h = e.result.hierarchy;
  REQUIRE(h.c0 == 1);
  const auto report = sssp::paper::validate_ordered_hierarchy(h, g);
  INFO(report.summary());
  REQUIRE(report.ok());
  std::vector<W> weights;
  for (const auto& j : in.joins.joins) weights.push_back(j.weight);
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.is_leaf(static_cast<NodeId>(v))) continue;
    CHECK(std::find(weights.begin(), weights.end(), h.nodes[v].height) != weights.end());
  }
  mb::require_threshold_sccs(h, e.node_upper, e.node_upper_inclusive);
}

// ------------------------------------------------------------------ pinned outputs
//
// P5-LOCAL-UNIFY (Q-59): Mode A and Mode B now call one shared assembly core
// (paper/local/assembly_core.hpp), so comparing the Mode B engine with Mode A field for field
// would compare the core with itself. The reference
// side of that comparison is replaced by the same values in fixed form: FNV-1a hashes of every
// field of the assemblies (unit/paper/local/assembly_hash.hpp) and of the counters, computed
// with the code before the refactor (two separate implementations, commit d2a75b1) on the
// instance lists below. This is not a change of expected values: the hashes are the outputs
// the old field-for-field test compared, recorded instead of recomputed.

/// kPaperRange: N in {32, 41}, where the Sec 5.4 formulas apply (l.1206). kSmall: N < 32,
/// accepted by D-8 and D-11 only.
enum class Range { kPaperRange, kSmall };

template <class W, class Fn>
void for_each_pinned_instance(Range range, Fn fn) {
  if (range == Range::kPaperRange) {
    for (const auto f : {lt::Family::kPath, lt::Family::kSpider, lt::Family::kRandomTree,
                         lt::Family::kSparse, lt::Family::kGrid, lt::Family::kComplete}) {
      for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
        for (const std::size_t n : {32U, 41U}) {
          const auto g = CsrGraph<W>::from_edges(n, lt::make_graph<W>(f, wd, n, 31 * n + 7));
          for (const std::uint64_t b : {1U, 2U, 3U, 4U, 6U, static_cast<unsigned>(n - 1)}) {
            INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " b=" << b);
            fn(build(g, b));
          }
        }
      }
    }
    // Clustered graphs: unaligned blocks under the paper's c = 10.
    if constexpr (std::is_same_v<W, std::uint64_t>) {
      for (const std::size_t n : {32U, 40U, 64U}) {
        for (const std::uint64_t b : {2U, 3U, 5U, 8U}) {
          INFO("clustered n=" << n << " b=" << b);
          fn(build(CsrGraph<W>::from_edges(n, lt::clustered(n, n + b)), b));
        }
      }
    }
    return;
  }
  for (const auto f : {lt::Family::kPath, lt::Family::kRandomTree, lt::Family::kSparse,
                       lt::Family::kGrid, lt::Family::kComplete}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (const std::size_t n : {2U, 3U, 5U, 8U, 13U, 21U}) {
        const auto g = CsrGraph<W>::from_edges(n, lt::make_graph<W>(f, wd, n, 1009 * n + 3));
        for (const std::uint64_t b : mb::separators(n)) {
          INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " b=" << b);
          fn(build(g, b));
        }
      }
    }
  }
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    for (const std::size_t n : {6U, 12U, 20U, 31U}) {
      for (const std::uint64_t b : {1U, 2U, 3U, 5U}) {
        INFO("clustered n=" << n << " b=" << b);
        fn(build(CsrGraph<W>::from_edges(n, lt::clustered(n, 3 * n + b)), b));
      }
    }
  }
}

/// Hash of all assemblies and of the accumulated counters that `assemble(instance, diag)`
/// produces on the instances of `ranges`.
struct Digest {
  std::uint64_t assemblies = 0;
  std::uint64_t counters = 0;
  std::size_t instances = 0;
};

template <class W, class Assemble>
Digest digest_of(std::initializer_list<Range> ranges, Assemble assemble) {
  lt::Fnv1a fa;
  Diagnostics diag;
  Digest d;
  for (const Range r : ranges) {
    for_each_pinned_instance<W>(r, [&](const Instance<W>& in) {
      lt::mix_assembly(fa, assemble(in, &diag));
      ++d.instances;
    });
  }
  lt::Fnv1a fc;
  lt::mix_counters(fc, diag.counters);
  CHECK(diag.fallbacks.empty());
  d.assemblies = fa.digest();
  d.counters = fc.digest();
  return d;
}

/// Hashes recorded before the refactor (see the comment above).
struct Pinned {
  Digest paper;   ///< Sec 5.4 with c = 10 and N^16 sigma, kPaperRange
  Digest lambda;  ///< D-8 with the grid, kPaperRange then kSmall
  Digest exact;   ///< D-11, kPaperRange then kSmall
};

template <class W>
Pinned pinned() {
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    return {{UINT64_C(0x15ded7254bec33ed), UINT64_C(0xd28d00acadb89ae6), 228},
            {UINT64_C(0x1319e72d0e13e570), UINT64_C(0x00fde01de5c552f3), 739},
            {UINT64_C(0x312d3a70d74a89c2), UINT64_C(0x936795cdfa2b8cc9), 739}};
  } else {
    return {{UINT64_C(0x8093d9f81fef05be), UINT64_C(0xd6f1a46753d31cd1), 216},
            {UINT64_C(0xa51f29ff08db7d57), UINT64_C(0x698a33677e0b1f3f), 711},
            {UINT64_C(0xac6fd3a89b4e41d2), UINT64_C(0xb93d2ad66c963ce6), 711}};
  }
}

void require_pinned(const Digest& got, const Digest& want) {
  CHECK(got.instances == want.instances);
  CHECK(got.assemblies == want.assemblies);
#if SSSP_DIAGNOSTICS
  CHECK(got.counters == want.counters);  // counters are compiled out of timed builds
#endif
}

}  // namespace

TEMPLATE_TEST_CASE(
    "Sec 5.4 assembly with the paper's constants equals the pinned pre-refactor "
    "output on families",
    "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  const auto want = pinned<W>().paper;
  // The shared core with the paper's rule (c = 10, N^16 sigma).
  require_pinned(digest_of<W>({Range::kPaperRange},
                              [](const Instance<W>& in, Diagnostics* diag) {
                                return sssp::paper::assembly::assemble(
                                    in.od, in.rs, in.joins,
                                    sssp::paper::assembly::Rule::paper(in.od.num_vertices()), diag);
                              }),
                 want);
  // Mode A.
  require_pinned(digest_of<W>({Range::kPaperRange},
                              [](const Instance<W>& in, Diagnostics* diag) {
                                return sssp::paper::assemble_local_hierarchy(in.od, in.rs, in.joins,
                                                                             diag);
                              }),
                 want);
  // Mode B with the paper's alignment and the grid.
  require_pinned(digest_of<W>({Range::kPaperRange},
                              [](const Instance<W>& in, Diagnostics* diag) {
                                return pe::assemble_local_hierarchy(in.od, in.rs, in.joins,
                                                                    pe::LocalAlignment::kPaper,
                                                                    pe::LocalAssembly::kGrid, diag);
                              }),
                 want);
}

TEMPLATE_TEST_CASE("D-8 and D-11 assemblies equal the pinned pre-refactor output on families",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  const auto variant = [](pe::LocalAlignment alignment, pe::LocalAssembly assembly) {
    return digest_of<W>(
        {Range::kPaperRange, Range::kSmall}, [=](const Instance<W>& in, Diagnostics* diag) {
          return pe::assemble_local_hierarchy(in.od, in.rs, in.joins, alignment, assembly, diag);
        });
  };
  require_pinned(variant(pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid),
                 pinned<W>().lambda);
  require_pinned(variant(pe::LocalAlignment::kLambda, pe::LocalAssembly::kExactHeight),
                 pinned<W>().exact);
  // D-11 does not use the alignment (Q-61).
  require_pinned(variant(pe::LocalAlignment::kPaper, pe::LocalAssembly::kExactHeight),
                 pinned<W>().exact);
}

TEMPLATE_TEST_CASE("D-8 and D-11: ordered hierarchies, threshold SCCs and scale relations",
                   "[paper_experimental][local]", std::uint64_t, double) {
  using W = TestType;
  Coverage cov;
  const auto run = [&](const CsrGraph<W>& g, std::uint64_t b) {
    const auto in = build(g, b);
    const auto grid = pe::assemble_local_hierarchy(
        in.od, in.rs, in.joins, pe::LocalAlignment::kLambda, pe::LocalAssembly::kGrid);
    const auto exact = pe::assemble_local_hierarchy(
        in.od, in.rs, in.joins, pe::LocalAlignment::kLambda, pe::LocalAssembly::kExactHeight);
    // D-8 with the grid: Def 4.1 (c0 = 2) and Lemma 5.10's invariant.
    REQUIRE(grid.result.hierarchy.c0 == 2);
    const auto report = sssp::paper::validate_ordered_hierarchy(grid.result.hierarchy, g);
    INFO(report.summary());
    REQUIRE(report.ok());
    mb::require_threshold_sccs(grid.result.hierarchy, grid.node_upper, grid.node_upper_inclusive);
    check_lambda_scales(in, grid, cov);
    // D-11.
    check_exact(g, in, exact);
    check_grid_coarsens_exact(grid, exact, cov);
    // The paper's Sec 5.4 formulas where they apply, against the same exact hierarchy.
    if (sssp::paper::local_assembly_applies(g.num_vertices(), static_cast<std::size_t>(b))) {
      const auto paper_grid = sssp::paper::assemble_local_hierarchy(in.od, in.rs, in.joins);
      check_grid_coarsens_exact(paper_grid, exact, cov);
    }
  };
  for (const auto f : {lt::Family::kPath, lt::Family::kRandomTree, lt::Family::kSparse,
                       lt::Family::kGrid, lt::Family::kComplete}) {
    for (const auto wd : {lt::Weights::kTies, lt::Weights::kUniform, lt::Weights::kExponential}) {
      for (std::size_t n = 2; n <= 40; n += (n < 12 ? 1 : 7)) {
        const auto g = CsrGraph<W>::from_edges(n, lt::make_graph<W>(f, wd, n, 1009 * n + 3));
        for (const std::uint64_t b : mb::separators(n)) {
          INFO(lt::name(f) << " " << lt::name(wd) << " n=" << n << " b=" << b);
          run(g, b);
        }
      }
    }
  }
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    for (const std::size_t n : {6U, 12U, 20U, 33U}) {
      for (const std::uint64_t b : {1U, 2U, 3U, 5U}) {
        INFO("clustered n=" << n << " b=" << b);
        run(CsrGraph<W>::from_edges(n, lt::clustered(n, 3 * n + b)), b);
      }
    }
    CHECK(cov.unaligned > 0);
  }
  CHECK(cov.aligned > 0);
  CHECK(cov.unattached > 0);
  CHECK(cov.global_nodes > 0);
}
