// Mode B assembly of the local hierarchy (task P5-EXP-LOCAL): Sec 5.4 with the paper's or the
// D-8 alignment constants, or the D-11 exact-height assembly. See local_assembly.hpp.
//
// The steps are those of the shared core (paper/local/assembly_core.hpp, task P5-LOCAL-UNIFY);
// this file only chooses the constants it passes.

#include "experimental/paper_experimental/local_assembly.hpp"

#include <limits>
#include <stdexcept>
#include <string>

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
  // Step 3 test and step 4 support end: l.1206 / l.1216, or D-8's Lambda formulas.
  paper::assembly::Rule rule = paper::assembly::Rule::paper(n);
  if (alignment == LocalAlignment::kLambda) {
    rule.align = lambda_scale(n, od.separator, rs.shell_ratio_log2);
    rule.support_end = lambda_support_scale(n, od.separator, rs.shell_ratio_log2);
  }
  rule.exact_height = assembly == LocalAssembly::kExactHeight;  // D-11: the alignment is unused
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
