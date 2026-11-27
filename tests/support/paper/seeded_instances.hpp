#pragma once
// Seeded ordered hierarchies (Def 4.2) for testing Evaluate (Alg 3) independently of Sec 5.
// Test-only; every random choice comes from SplitMix64 with a fixed seed, so an instance is
// reproduced from (family, parameters, seed) on every platform.
//
// Families
//   reference_instance   threshold-SCC reference hierarchy (c0 = 1) of a random directed state
//                        graph; random projection onto fewer physical vertices (copies); the
//                        physical graph is the arc multigraph, optionally with extra random
//                        physical edges (shortcuts make state arcs strictly longer than the
//                        physical distance, so fewer states are protected and supplemental
//                        offers matter). Single source (s, 0) or several weighted sources.
//   zero_arc_instance    as above with extra length-0 arcs that never close a cycle through
//                        two or more states (Def 4.1 forbids those, Q-21); length-0 self-loops
//                        are included. Built on a copy where length 0 is replaced by 1 and
//                        positive lengths are even, which gives the same forest; the zero
//                        lengths are then restored (links keep their bound 1 >= 0, H5/H6 hold).
//   base_case_instance   Alg 4 base case (local_hierarchy_base_case) of a positive-weight
//                        physical graph: one state per vertex, both arc orientations; offers
//                        c_i at the native entry state of v_i, so every reachable vertex is
//                        served by a protected state.
// Initial offers (all feasible by construction, value >= delta(pi(x))):
//   * per source i, one state with image v_i gets offer c_i (if v_i has a state);
//   * `extra_offers` random states with finite delta get delta + slack, slack = 0 with
//     probability 1/2 (tight offers in the middle of the graph) and otherwise in [1, max_len].

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/hierarchy/reference_hierarchy.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/graph.hpp"
#include "support/graphs.hpp"
#include "support/paper/eval_oracle.hpp"
#include "support/prng.hpp"

namespace sssp::test::paper {

template <class W>
struct SeededInstance {
  std::string name;
  std::uint64_t seed = 0;
  std::size_t num_physical = 0;
  CsrGraph<W> physical;
  SeededOrderedHierarchy<W> seeded;
};

struct ReferenceParams {
  std::size_t num_states = 12;
  std::size_t num_physical = 8;  ///< <= num_states; states 0..num_physical-1 cover all vertices
  std::size_t num_arcs = 30;
  std::uint64_t min_len = 1;  ///< >= 1 (the reference hierarchy needs positive lengths)
  std::uint64_t max_len = 10;
  std::size_t num_sources = 1;    ///< 1 = ordinary source (s, 0)
  std::uint64_t max_offset = 10;  ///< source offsets in [0, max_offset] when num_sources > 1
  std::size_t extra_physical_edges = 0;
  std::size_t extra_offers = 2;
  std::size_t zero_arcs = 0;  ///< zero_arc_instance only
};

namespace detail {

template <class W>
W w_of(std::uint64_t x) {
  return static_cast<W>(x);
}

/// Sources and feasible offers for an already built hierarchy (see the header comment).
template <class W>
void seed_offers(SeededInstance<W>& inst, std::size_t num_sources, std::uint64_t max_offset,
                 std::size_t extra_offers, std::uint64_t max_len, SplitMix64& rng) {
  auto& s = inst.seeded;
  const auto& proj = s.hierarchy.projection;
  const std::size_t n = proj.size();
  if (n == 0) return;
  s.sources.clear();
  s.offers.clear();
  for (std::size_t i = 0; i < num_sources; ++i) {
    const auto x = static_cast<StateId>(rng.below(n));
    const W c = num_sources == 1 ? W{0} : w_of<W>(rng.between(0, max_offset));
    s.sources.push_back({proj[x], c});
    s.offers.push_back({x, c});  // feasible: delta(pi(x)) <= c + d(pi(x), pi(x)) = c
  }
  const auto target = target_distances<W>(inst.physical, s.sources);
  for (std::size_t k = 0; k < extra_offers; ++k) {
    const auto x = static_cast<StateId>(rng.below(n));
    const W t = target[proj[x]];
    if (t == WeightTraits<W>::infinity()) continue;
    const std::uint64_t slack =
        rng.below(2) == 0 ? 0 : rng.between(1, std::max<std::uint64_t>(1, max_len));
    s.offers.push_back({x, WeightTraits<W>::checked_add(t, w_of<W>(slack))});
  }
}

template <class W>
std::vector<Vertex> random_projection(std::size_t num_states, std::size_t num_physical,
                                      SplitMix64& rng) {
  std::vector<Vertex> proj(num_states);
  for (std::size_t x = 0; x < num_states; ++x) {
    proj[x] =
        x < num_physical ? static_cast<Vertex>(x) : static_cast<Vertex>(rng.below(num_physical));
  }
  return proj;
}

template <class W>
CsrGraph<W> physical_with_extras(const StateGraph<W>& g, const std::vector<Vertex>& proj,
                                 std::size_t num_physical, std::size_t extra, std::uint64_t min_len,
                                 std::uint64_t max_len, SplitMix64& rng) {
  const CsrGraph<W> base = sssp::paper::physical_graph_of_arcs<W>(g, proj, num_physical);
  std::vector<Edge<W>> edges(base.edges().begin(), base.edges().end());
  for (std::size_t i = 0; i < extra && num_physical > 0; ++i) {
    edges.push_back({static_cast<Vertex>(rng.below(num_physical)),
                     static_cast<Vertex>(rng.below(num_physical)),
                     w_of<W>(rng.between(std::max<std::uint64_t>(1, min_len), max_len))});
  }
  return CsrGraph<W>::from_edges(num_physical, std::move(edges));
}

}  // namespace detail

/// Threshold-SCC reference hierarchy of a random directed state graph (see header comment).
template <class W>
SeededInstance<W> reference_instance(const ReferenceParams& p, std::uint64_t seed) {
  SplitMix64 rng(seed);
  StateGraph<W> g;
  g.num_states = p.num_states;
  for (std::size_t i = 0; i < p.num_arcs && p.num_states > 0; ++i) {
    g.arcs.push_back({static_cast<StateId>(rng.below(p.num_states)),
                      static_cast<StateId>(rng.below(p.num_states)),
                      detail::w_of<W>(rng.between(p.min_len, p.max_len))});
  }
  auto proj = detail::random_projection<W>(p.num_states, p.num_physical, rng);
  SeededInstance<W> inst;
  inst.name = "reference";
  inst.seed = seed;
  inst.num_physical = p.num_physical;
  inst.physical = detail::physical_with_extras<W>(g, proj, p.num_physical, p.extra_physical_edges,
                                                  p.min_len, p.max_len, rng);
  inst.seeded.hierarchy = sssp::paper::build_threshold_scc_hierarchy<W>(g, proj);
  detail::seed_offers(inst, p.num_sources, p.max_offset, p.extra_offers, p.max_len, rng);
  return inst;
}

/// Reference hierarchy with p.zero_arcs extra length-0 arcs that are acyclic (forward in a
/// random permutation of the states) plus occasional length-0 self-loops.
template <class W>
SeededInstance<W> zero_arc_instance(const ReferenceParams& p, std::uint64_t seed) {
  SplitMix64 rng(seed);
  const std::size_t n = p.num_states;
  std::vector<std::size_t> rank(n);
  std::iota(rank.begin(), rank.end(), std::size_t{0});
  for (std::size_t i = n; i > 1; --i) std::swap(rank[i - 1], rank[rng.below(i)]);
  StateGraph<W> zeroed;  // the instance's real state graph
  zeroed.num_states = n;
  StateGraph<W> build;  // zero -> 1, positive lengths doubled (all positive, same forest)
  build.num_states = n;
  for (std::size_t i = 0; i < p.num_arcs && n > 0; ++i) {
    const auto u = static_cast<StateId>(rng.below(n));
    const auto v = static_cast<StateId>(rng.below(n));
    const std::uint64_t len = 2 * rng.between(p.min_len, p.max_len);
    zeroed.arcs.push_back({u, v, detail::w_of<W>(len)});
    build.arcs.push_back({u, v, detail::w_of<W>(len)});
  }
  for (std::size_t i = 0; i < p.zero_arcs && n > 1; ++i) {
    auto u = static_cast<StateId>(rng.below(n));
    auto v = static_cast<StateId>(rng.below(n));
    if (u != v && rank[u] > rank[v]) std::swap(u, v);  // forward in the permutation, or a loop
    zeroed.arcs.push_back({u, v, W{0}});
    build.arcs.push_back({u, v, W{1}});
  }
  auto proj = detail::random_projection<W>(n, p.num_physical, rng);
  SeededInstance<W> inst;
  inst.name = "zero_arcs";
  inst.seed = seed;
  inst.num_physical = p.num_physical;
  auto h = sssp::paper::build_threshold_scc_hierarchy<W>(build, proj);
  h.graph = zeroed;  // same arcs, same ids; only the zero lengths differ
  inst.seeded.hierarchy = std::move(h);
  inst.physical = detail::physical_with_extras<W>(
      zeroed, proj, p.num_physical, p.extra_physical_edges, 2 * p.min_len, 2 * p.max_len, rng);
  detail::seed_offers(inst, p.num_sources, p.max_offset, p.extra_offers, p.max_len, rng);
  return inst;
}

/// Alg 4 base case of the physical graph (n, edges) with the given sources; offers c_i at the
/// native entry state of v_i. Edge weights must be positive (the base case rejects zeros).
template <class W>
SeededInstance<W> base_case_instance(std::string name, std::size_t n, std::vector<Edge<W>> edges,
                                     std::vector<PhysicalSource<W>> sources,
                                     std::uint64_t seed = 0) {
  SeededInstance<W> inst;
  inst.name = std::move(name);
  inst.seed = seed;
  inst.num_physical = n;
  inst.physical = CsrGraph<W>::from_edges(n, std::move(edges));
  auto local = sssp::paper::local_hierarchy_base_case<W>(inst.physical);
  inst.seeded.hierarchy = std::move(local.hierarchy);
  for (const auto& src : sources) {
    inst.seeded.offers.push_back({local.native_entry.at(src.vertex), src.offset});
  }
  inst.seeded.sources = std::move(sources);
  return inst;
}

/// Base-case instances over the shapes of tests/support/graphs.hpp with positive weights:
/// single vertex, path, star, cycle, complete, tree, disconnected, grid and layered graphs
/// with all-equal weights (many ties), random multigraphs with weights in [1, 2] (ties) and
/// [1, 20]. Each shape appears with one ordinary source and with three weighted sources.
template <class W>
std::vector<SeededInstance<W>> base_case_catalogue(std::uint64_t seed) {
  struct Shape {
    std::string name;
    std::size_t n;
    std::vector<Edge<W>> edges;
  };
  const W one = W{1};
  std::vector<Shape> shapes;
  shapes.push_back({"single_vertex", 1, {}});
  shapes.push_back(
      {"path", 7, path_edges<W>(7, [](std::size_t i) { return static_cast<W>(i + 1); })});
  shapes.push_back({"star", 8, star_edges<W>(8)});
  shapes.push_back({"cycle_equal", 9, cycle_edges<W>(9, one)});
  shapes.push_back({"complete", 7, complete_edges<W>(7, 1, 20, seed + 1)});
  shapes.push_back({"tree", 15, random_tree_edges<W>(15, 1, 9, seed + 2)});
  shapes.push_back({"disconnected", 7, {{0, 1, W{2}}, {1, 2, W{2}}, {3, 4, one}, {5, 6, W{3}}}});
  shapes.push_back({"grid_equal", 16, grid_edges<W>(4, 4, W{3})});
  shapes.push_back({"layered_equal", 1 + 4 * 3, layered_edges<W>(4, 3, one)});
  shapes.push_back({"random_ties", 12, random_edges<W>(12, 30, 1, 2, seed + 3)});
  shapes.push_back({"random_wide", 12, random_edges<W>(12, 30, 1, 20, seed + 4)});
  // random_edges may produce self-loops; they are positive here, which the base case accepts.
  std::vector<SeededInstance<W>> out;
  SplitMix64 rng(seed);
  for (auto& sh : shapes) {
    out.push_back(base_case_instance<W>(sh.name + "/single", sh.n, sh.edges, {{0, W{0}}}, seed));
    std::vector<PhysicalSource<W>> multi;
    for (int i = 0; i < 3; ++i) {
      multi.push_back({static_cast<Vertex>(rng.below(sh.n)), static_cast<W>(rng.between(0, 6))});
    }
    out.push_back(base_case_instance<W>(sh.name + "/multi", sh.n, std::move(sh.edges),
                                        std::move(multi), seed));
  }
  return out;
}

/// `count` reference instances with seeds seed0, seed0 + 1, ...: alternating single / multiple
/// sources, with and without physical shortcuts, small and tie-heavy length ranges.
template <class W>
std::vector<SeededInstance<W>> reference_catalogue(std::uint64_t seed0, std::size_t count) {
  std::vector<SeededInstance<W>> out;
  for (std::size_t k = 0; k < count; ++k) {
    const std::uint64_t seed = seed0 + k;
    SplitMix64 rng(seed ^ 0x5eedULL);
    ReferenceParams p;
    p.num_states = 2 + rng.below(30);
    p.num_physical = 1 + rng.below(p.num_states);
    p.num_arcs = rng.below(4 * p.num_states + 1);
    p.max_len = rng.below(2) == 0 ? 2 : 12;  // tie-heavy or spread lengths
    p.num_sources = 1 + (k % 2) * rng.below(4);
    p.extra_physical_edges = (k % 3 == 0) ? rng.below(p.num_physical + 1) : 0;
    p.extra_offers = rng.below(4);
    auto inst = reference_instance<W>(p, seed);
    inst.name = "reference/" + std::to_string(seed);
    out.push_back(std::move(inst));
  }
  return out;
}

}  // namespace sssp::test::paper
