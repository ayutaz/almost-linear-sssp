#pragma once
// Test support for the Mode B local hierarchy (task P5-EXP-LOCAL): structural equality of
// hierarchies and assemblies, the threshold-SCC characterisation of the proof of Lemma 5.10
// (Sec 5.4, l.1271-1275) by a Kosaraju search written here, and Evaluate from a native entry
// against a Dijkstra written in the test support. Nothing here calls the code under test except
// through the arguments.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "experimental/paper_experimental/params.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/base_case.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"
#include "support/paper/eval_oracle.hpp"
#include "unit/paper/local/local_test_support.hpp"

namespace sssp::test::mode_b {

namespace pe = sssp::paper_experimental;
namespace sp = sssp::paper;

/// Mode B parameters: threshold k > coeff * L^exponent, constant separator b.
inline pe::PaperParams make_params(pe::LocalAlignment alignment, pe::LocalAssembly assembly,
                                   std::uint64_t b, pe::SizeThreshold threshold = {0, 0}) {
  pe::LocalParams l;
  l.threshold = threshold;
  l.separator = pe::SeparatorRule::constant(b);
  l.alignment = alignment;
  l.assembly = assembly;
  return {l, pe::OuterParams{}};
}

/// ceil(log2 n) for n >= 2 (1 for n <= 2), by integer doubling.
inline std::uint64_t ceil_log2(std::uint64_t n) {
  std::uint64_t k = 0;
  while ((std::uint64_t{1} << k) < n) ++k;
  return std::max<std::uint64_t>(k, 1);
}

/// The separator sizes of the acceptance criteria, b in {1, 2, 3, 4, ceil(log2 N)}, plus the
/// two degenerate extremes N and N + 1 (one pivot, no macro edges).
inline std::vector<std::uint64_t> separators(std::size_t n) {
  std::vector<std::uint64_t> b{1, 2, 3, 4, ceil_log2(n), n, n + 1};
  std::sort(b.begin(), b.end());
  b.erase(std::unique(b.begin(), b.end()), b.end());
  return b;
}

template <class W>
void require_same_hierarchy(const sp::OrderedHierarchy<W>& a, const sp::OrderedHierarchy<W>& b) {
  REQUIRE(a.c0 == b.c0);
  REQUIRE(a.graph.num_states == b.graph.num_states);
  REQUIRE(a.graph.arcs.size() == b.graph.arcs.size());
  for (std::size_t i = 0; i < a.graph.arcs.size(); ++i) {
    INFO("arc " << i);
    REQUIRE(a.graph.arcs[i].tail == b.graph.arcs[i].tail);
    REQUIRE(a.graph.arcs[i].head == b.graph.arcs[i].head);
    REQUIRE(a.graph.arcs[i].length == b.graph.arcs[i].length);
  }
  REQUIRE(a.projection == b.projection);
  REQUIRE(a.leaf_of_state == b.leaf_of_state);
  REQUIRE(a.roots == b.roots);
  REQUIRE(a.nodes.size() == b.nodes.size());
  for (std::size_t v = 0; v < a.nodes.size(); ++v) {
    INFO("node " << v);
    const auto& x = a.nodes[v];
    const auto& y = b.nodes[v];
    REQUIRE(x.parent == y.parent);
    REQUIRE(x.children == y.children);
    REQUIRE(x.state == y.state);
    REQUIRE(x.height == y.height);
    REQUIRE(x.links.size() == y.links.size());
    for (std::size_t i = 0; i < x.links.size(); ++i) {
      REQUIRE(x.links[i].x == y.links[i].x);
      REQUIRE(x.links[i].y == y.links[i].y);
      REQUIRE(x.links[i].length == y.links[i].length);
      REQUIRE(x.links[i].arc == y.links[i].arc);
    }
  }
}

template <class W>
void require_same_result(const sp::LocalHierarchyResult<W>& a,
                         const sp::LocalHierarchyResult<W>& b) {
  require_same_hierarchy(a.hierarchy, b.hierarchy);
  REQUIRE(a.native_entry == b.native_entry);
  REQUIRE(a.arc_edge == b.arc_edge);
}

template <class W>
void require_same_assembly(const sp::LocalAssembly<W>& a, const sp::LocalAssembly<W>& b) {
  require_same_result(a.result, b.result);
  REQUIRE(a.blocks.size() == b.blocks.size());
  for (std::size_t u = 0; u < a.blocks.size(); ++u) {
    INFO("block " << u);
    REQUIRE(a.blocks[u].kind == b.blocks[u].kind);
    REQUIRE(a.blocks[u].size == b.blocks[u].size);
    REQUIRE(a.blocks[u].s == b.blocks[u].s);
    REQUIRE(a.blocks[u].t == b.blocks[u].t);
    REQUIRE(a.blocks[u].z == b.blocks[u].z);
    REQUIRE(a.blocks[u].sigma_key == b.blocks[u].sigma_key);
  }
  REQUIRE(a.y_component == b.y_component);
  REQUIRE(a.node_upper == b.node_upper);
  REQUIRE(a.node_upper_inclusive == b.node_upper_inclusive);
  REQUIRE(a.node_origin == b.node_origin);
  REQUIRE(a.grid == b.grid);
  REQUIRE(a.grid_upper == b.grid_upper);
  REQUIRE(a.supports.size() == b.supports.size());
  for (std::size_t c = 0; c < a.supports.size(); ++c) {
    REQUIRE(a.supports[c].lower == b.supports[c].lower);
    REQUIRE(a.supports[c].upper == b.supports[c].upper);
    REQUIRE(a.supports[c].lower_index == b.supports[c].lower_index);
    REQUIRE(a.supports[c].upper_index == b.supports[c].upper_index);
  }
}

/// Counters of `d` without the parameter log ("params.*") and without the keys starting with
/// `skip_prefix` (the Mode B entry's own counters).
inline std::map<std::string, std::uint64_t> counters_without(const Diagnostics& d,
                                                             const std::string& skip_prefix) {
  std::map<std::string, std::uint64_t> out;
  for (const auto& [k, v] : d.counters) {
    if (k.rfind("params.", 0) == 0) continue;
    if (!skip_prefix.empty() && k.rfind(skip_prefix, 0) == 0) continue;
    out.emplace(k, v);
  }
  return out;
}

/// Kosaraju SCC ids of the state arcs accepted by `keep`.
template <class W, class Keep>
std::vector<std::size_t> scc(const sp::StateGraph<W>& d, Keep keep) {
  const std::size_t s = d.num_states;
  std::vector<std::vector<sp::StateId>> fwd(s), rev(s);
  for (const auto& a : d.arcs) {
    if (!keep(a.length)) continue;
    fwd[a.tail].push_back(a.head);
    rev[a.head].push_back(a.tail);
  }
  std::vector<std::size_t> order;
  std::vector<bool> seen(s, false);
  for (std::size_t r = 0; r < s; ++r) {
    if (seen[r]) continue;
    std::vector<std::pair<sp::StateId, std::size_t>> st{{static_cast<sp::StateId>(r), 0}};
    seen[r] = true;
    while (!st.empty()) {
      auto& [x, i] = st.back();
      if (i < fwd[x].size()) {
        const sp::StateId y = fwd[x][i++];
        if (!seen[y]) {
          seen[y] = true;
          st.emplace_back(y, 0);
        }
      } else {
        order.push_back(x);
        st.pop_back();
      }
    }
  }
  std::vector<std::size_t> comp(s, s);
  std::size_t nc = 0;
  for (std::size_t i = order.size(); i-- > 0;) {
    const std::size_t r = order[i];
    if (comp[r] != s) continue;
    std::vector<sp::StateId> st{static_cast<sp::StateId>(r)};
    comp[r] = nc;
    while (!st.empty()) {
      const sp::StateId x = st.back();
      st.pop_back();
      for (const sp::StateId y : rev[x]) {
        if (comp[y] == s) {
          comp[y] = nc;
          st.push_back(y);
        }
      }
    }
    ++nc;
  }
  return comp;
}

template <class W>
std::vector<sp::StateId> states_below(const sp::OrderedHierarchy<W>& h, sp::NodeId v) {
  std::vector<sp::StateId> out;
  std::vector<sp::NodeId> st{v};
  while (!st.empty()) {
    const sp::NodeId x = st.back();
    st.pop_back();
    if (h.is_leaf(x)) out.push_back(h.nodes[x].state);
    for (const sp::NodeId c : h.nodes[x].children) st.push_back(c);
  }
  std::sort(out.begin(), out.end());
  return out;
}

/// The SCC of `x` as a sorted state list.
inline std::vector<sp::StateId> members(const std::vector<std::size_t>& comp, sp::StateId x) {
  std::vector<sp::StateId> out;
  for (std::size_t y = 0; y < comp.size(); ++y) {
    if (comp[y] == comp[x]) out.push_back(static_cast<sp::StateId>(y));
  }
  return out;
}

/// Proof of Lemma 5.10 (l.1271-1275): below every internal node C lie exactly the states of one
/// SCC of the arcs shorter than its class upper bound (at most it, if inclusive), and below each
/// child the states of one SCC of the arcs shorter than h(C); the roots are the SCCs of D.
template <class W>
void require_threshold_sccs(const sp::OrderedHierarchy<W>& h, const std::vector<W>& upper,
                            const std::vector<std::uint8_t>& inclusive) {
  std::map<std::pair<W, int>, std::vector<std::size_t>> cache;
  const auto comp_of = [&](W bound, int incl) -> const std::vector<std::size_t>& {
    auto it = cache.find({bound, incl});
    if (it != cache.end()) return it->second;
    auto c = incl == 0 ? scc(h.graph, [&](W l) { return l < bound; })
                       : scc(h.graph, [&](W l) { return !(bound < l); });
    return cache.emplace(std::make_pair(bound, incl), std::move(c)).first->second;
  };
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    const auto node = static_cast<sp::NodeId>(v);
    if (h.is_leaf(node)) continue;
    INFO("node " << v << " height " << h.nodes[v].height << " upper " << upper[v]);
    const auto below = states_below(h, node);
    REQUIRE(members(comp_of(upper[v], inclusive[v]), below.front()) == below);
    const auto& low = comp_of(h.nodes[v].height, 0);
    for (const sp::NodeId c : h.nodes[v].children) {
      const auto cb = states_below(h, c);
      REQUIRE(members(low, cb.front()) == cb);
    }
  }
  const auto all = scc(h.graph, [](W) { return true; });
  for (const sp::NodeId r : h.roots) {
    const auto below = states_below(h, r);
    REQUIRE(members(all, below.front()) == below);
  }
}

/// Seeds the native entry of `src` with 0, checks Def 4.1 / 4.2, runs Evaluate and returns the
/// least emission per physical vertex; every vertex must equal the test-support Dijkstra and
/// every emission must be at least the distance of its image. Counters of Evaluate are added to
/// `eval_counters` if non-null.
template <class W>
std::vector<W> evaluate_and_check(const sp::LocalHierarchyResult<W>& lh, const CsrGraph<W>& g,
                                  const std::vector<Edge<W>>& edges, Vertex src,
                                  std::map<std::string, std::uint64_t>* eval_counters = nullptr) {
  const std::size_t n = g.num_vertices();
  const auto& h = lh.hierarchy;
  sp::SeededOrderedHierarchy<W> seeded{h, {{src, W{0}}}, {{lh.native_entry[src], W{0}}}};
  const auto valid = sp::validate_seeded_ordered_hierarchy(seeded, g);
  INFO(valid.summary());
  REQUIRE(valid.ok());
  const auto r = sp::evaluate(seeded);
  const auto want = local::dijkstra_naive<W>(n, edges, {{src, W{0}}});
  const auto got = test::paper::physical_minima<W>(n, h.projection, r.emissions);
  for (Vertex v = 0; v < n; ++v) {
    INFO("source " << src << " vertex " << v);
    REQUIRE(got[v] == want[v]);
  }
  for (const auto& [x, val] : r.emissions) REQUIRE(!(val < want[h.projection[x]]));
  if (eval_counters != nullptr) {
    for (const auto& [k, v] : r.diag.counters) (*eval_counters)[k] += v;
  }
  return got;
}

}  // namespace sssp::test::mode_b
