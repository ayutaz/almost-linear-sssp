// Mode B assembly of the local hierarchy (task P5-EXP-LOCAL): Sec 5.4 with the paper's or the
// D-8 alignment constants, or the D-11 exact-height assembly. See local_assembly.hpp.
//
// The grid assembly is the shared core (paper/local/assembly_core.hpp, task P5-LOCAL-UNIFY)
// with the constants chosen here. The exact-height assembly (D-11) is not a step of the paper,
// so it is implemented here; it shares the core's ClassMerger and link_and_order.

#include "experimental/paper_experimental/local_assembly.hpp"

#include <algorithm>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/local/assembly_core.hpp"

namespace sssp::paper_experimental {

namespace {

using paper::assembly::LengthScale;

constexpr const char* kContext = "paper_experimental::assemble_local_hierarchy";

/// D-8: Lambda = 2 * N * 3 * h * K^2, h = b + 1: integer multiplications by doubling and
/// adding, K^2 = 2^(2 k_log2) by doublings (Sec 2 l.297); each step saturates (Q-38).
LengthScale lambda_scale(std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  if (b == std::numeric_limits<std::uint64_t>::max()) return LengthScale::saturated();  // h
  return LengthScale{}.times(2).times(n).times(3).times(b + 1).doubled(2 * k_log2);
}

/// D-8: the support end factor 2 * (3 h K^2)^2 * N^3.
LengthScale lambda_support_scale(std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  if (b == std::numeric_limits<std::uint64_t>::max()) return LengthScale::saturated();
  LengthScale s;
  s.times(2);
  for (int i = 0; i < 2; ++i) s.times(3).times(b + 1).doubled(2 * k_log2);
  return s.times(n, 3);
}

/// D-11 (Q-39): no grid. All join records sorted by (weight, index), the records of one weight
/// merged simultaneously at height = upper = that weight (inclusive); then the post-processing
/// of Lemma 5.10 (l.1281-1293). Links are arcs of length <= the exact height, hence c0 = 1.
template <class W>
paper::LocalAssembly<W> assemble_exact_height(const paper::RadialStateGraph<W>& rs,
                                              const paper::ThresholdJoins<W>& joins,
                                              Diagnostics* diag) {
  using paper::NodeId;
  const std::size_t s = rs.num_states();
  const auto& jr = joins.joins;

  paper::LocalAssembly<W> out;
  out.result.hierarchy = paper::make_leaf_forest(rs.graph, rs.projection);
  auto& h = out.result.hierarchy;
  h.c0 = 1;
  out.result.native_entry = rs.native_inward;
  out.result.arc_edge = rs.arc_edge;
  out.node_upper.assign(s, W{0});
  out.node_upper_inclusive.assign(s, 0);
  out.node_origin.assign(s, paper::NodeOrigin::kExactLocal);

  // Join vertices: the states (node x is the leaf of state x), then one auxiliary vertex per
  // pivot without a node (Lemma 3.6, l.1281).
  std::vector<NodeId> node_of_vertex(joins.num_join_vertices(), paper::kNoNode);
  for (std::size_t x = 0; x < joins.num_states; ++x) node_of_vertex[x] = static_cast<NodeId>(x);
  paper::assembly::ClassMerger merger(
      std::move(node_of_vertex),
      [&jr](std::uint32_t j) { return std::pair<std::uint32_t, std::uint32_t>{jr[j].a, jr[j].b}; },
      [&out](W height, std::vector<NodeId> children, W upper, bool inclusive,
             paper::NodeOrigin origin) {
        const NodeId r =
            paper::add_internal_node(out.result.hierarchy, height, std::move(children), {});
        out.node_upper.push_back(upper);
        out.node_upper_inclusive.push_back(inclusive ? 1 : 0);
        out.node_origin.push_back(origin);
        return r;
      });

  std::vector<std::uint32_t> order(jr.size());
  for (std::size_t i = 0; i < jr.size(); ++i) order[i] = static_cast<std::uint32_t>(i);
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (jr[a].weight < jr[b].weight) return true;
    if (jr[b].weight < jr[a].weight) return false;
    return a < b;
  });
  std::size_t exact_classes = 0;
  for (std::size_t i = 0; i < order.size();) {
    std::size_t e = i + 1;
    while (e < order.size() && !(jr[order[i]].weight < jr[order[e]].weight)) ++e;
    const W w = jr[order[i]].weight;
    if (w == WeightTraits<W>::infinity()) {  // Sec 5.3: only finite joins
      throw std::logic_error(std::string(kContext) + ": an infinite join record");
    }
    merger.merge_class(std::span<const std::uint32_t>(order.data() + i, e - i), w, w, true,
                       paper::NodeOrigin::kExactLocal);
    ++exact_classes;
    i = e;
  }

  auto lo = paper::assembly::link_and_order<W>(h.nodes, h.graph.arcs, h.leaf_of_state,
                                               out.node_upper, out.node_upper_inclusive, kContext);
  h.roots = std::move(lo.roots);

  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["local_assemble_calls"]);
    SSSP_COUNT(diag->counters["local_assemble_blocks"] += rs.num_blocks);
    SSSP_COUNT(diag->counters["local_assemble_exact_classes"] += exact_classes);
    SSSP_COUNT(diag->counters["local_assemble_nodes"] += h.nodes.size() - s);
    SSSP_COUNT(diag->counters["local_assemble_links"] += lo.links);
    SSSP_COUNT(merger.dsu().add_counters(*diag));
  }
  static_cast<void>(exact_classes);
  return out;
}

}  // namespace

template <class W>
W lambda_scaled(W s, std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  return lambda_scale(n, b, k_log2).apply(s);
}

template <class W>
W lambda_support_end(W sigma, std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  return lambda_support_scale(n, b, k_log2).apply(sigma);
}

bool local_assembly_variant_applies(std::size_t n, std::size_t b, LocalAlignment alignment,
                                    LocalAssembly assembly) {
  if (alignment == LocalAlignment::kPaper && assembly == LocalAssembly::kGrid) {
    return paper::local_assembly_applies(n, b);  // Sec 5.4 l.1206
  }
  return n >= 2;  // D-8 / D-11: Lemmas 5.2-5.9 with N >= 2 (research note 02 §4.2)
}

template <class W>
paper::LocalAssembly<W> assemble_local_hierarchy(const paper::OwnerDecomposition<W>& od,
                                                 const paper::RadialStateGraph<W>& rs,
                                                 const paper::ThresholdJoins<W>& joins,
                                                 LocalAlignment alignment, LocalAssembly assembly,
                                                 Diagnostics* diag) {
  paper::assembly::check_inputs(od, rs, joins, kContext);
  const std::size_t n = od.num_vertices();
  if (!local_assembly_variant_applies(n, od.separator, alignment, assembly)) {
    throw std::invalid_argument(
        std::string(kContext) + ": N = " + std::to_string(n) +
        ", b = " + std::to_string(od.separator) +
        (alignment == LocalAlignment::kPaper && assembly == LocalAssembly::kGrid
             ? " violates N >= 32 and h = b + 1 <= N (Sec 5.4 l.1206)"
             : " violates N >= 2 (D-8 / D-11)"));
  }
  if (assembly == LocalAssembly::kExactHeight) {
    return assemble_exact_height(rs, joins, diag);  // D-11: the alignment is unused
  }
  // Step 3 test and step 4 support end: l.1206 / l.1216, or D-8's Lambda formulas.
  paper::assembly::Rule rule = paper::assembly::Rule::paper(n);
  if (alignment == LocalAlignment::kLambda) {
    rule.align = lambda_scale(n, od.separator, rs.shell_ratio_log2);
    rule.support_end = lambda_support_scale(n, od.separator, rs.shell_ratio_log2);
  }
  rule.context = kContext;
  return paper::assembly::assemble(od, rs, joins, rule, diag);
}

template std::uint64_t lambda_scaled<std::uint64_t>(std::uint64_t, std::uint64_t, std::uint64_t,
                                                    unsigned);
template double lambda_scaled<double>(double, std::uint64_t, std::uint64_t, unsigned);
template std::uint64_t lambda_support_end<std::uint64_t>(std::uint64_t, std::uint64_t,
                                                         std::uint64_t, unsigned);
template double lambda_support_end<double>(double, std::uint64_t, std::uint64_t, unsigned);
template paper::LocalAssembly<std::uint64_t> assemble_local_hierarchy<std::uint64_t>(
    const paper::OwnerDecomposition<std::uint64_t>&, const paper::RadialStateGraph<std::uint64_t>&,
    const paper::ThresholdJoins<std::uint64_t>&, LocalAlignment, LocalAssembly, Diagnostics*);
template paper::LocalAssembly<double> assemble_local_hierarchy<double>(
    const paper::OwnerDecomposition<double>&, const paper::RadialStateGraph<double>&,
    const paper::ThresholdJoins<double>&, LocalAlignment, LocalAssembly, Diagnostics*);

}  // namespace sssp::paper_experimental
