#include "experimental/paper_experimental/sssp.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "paper/model/params.hpp"
#include "paper/model/random.hpp"
#include "paper/outer/base_case.hpp"
#include "paper/outer/contract_restore.hpp"

namespace sssp::paper_experimental {

paper::outer::OuterConstructionParams outer_construction_params(const PaperParams& params,
                                                                std::uint64_t n) {
  const LocalParams& local = params.local();
  const OuterParams& outer = params.outer();
  if (local.alignment != LocalAlignment::kPaper || local.assembly != LocalAssembly::kGrid) {
    throw std::invalid_argument(
        "paper_experimental: the atoms of the outer construction run the shared Alg 4; "
        "local.alignment = lambda (D-8) and local.assembly = exact-height (D-11) are unsupported");
  }
  const bool paper_local_separator =
      local.separator == SeparatorRule::power_of_l(paper::model::kLocalExponent);
  if (!paper_local_separator && local.separator.kind() != SeparatorRule::Kind::kConstant) {
    throw std::invalid_argument(
        "paper_experimental: local.separator = L_k^e with e != 4 is unsupported on the outer "
        "branch (one b for all atoms)");
  }
  if (outer.h_rule != HRule::kPaper) {
    throw std::invalid_argument(
        "paper_experimental: outer.h_rule = proof-minimum (D-9) is unsupported on the outer "
        "branch (H depends on the pivot count)");
  }

  paper::outer::OuterConstructionParams p;
  // Sec 7.1 l.1644 or D-7. b stays a power so that F comes from the exact h = b + 1 (Q-71).
  if (outer.separator != SeparatorRule::power_of_l(paper::model::kOuterExponent)) {
    p.separator = outer.separator.kind() == SeparatorRule::Kind::kConstant
                      ? paper::outer::OuterSeparator{outer.separator.parameter(), 1}
                      : paper::outer::OuterSeparator{
                            paper::model::L(n), static_cast<unsigned>(outer.separator.parameter())};
  }
  // Alg 5 l.3-6 or D-17, D-10; the paper's H (t is not needed).
  p.name_routes = params.name_routes_params(1);
  // Alg 4 l.2 and Sec 5.1 l.967 or D-7; independent of the atom size for the rules accepted
  // above (separator 0 keeps the shared component on b = L_k^4).
  p.local = params.local_size_params(n);
  return p;
}

template <class W>
SsspResult<W> paper_sssp(const CsrGraph<W>& g, Vertex s, const PaperParams& params,
                         std::uint64_t seed) {
  Diagnostics diag;
  params.record(diag);
  // Alg 8 l.2-3: contraction; N = number of retained vertices.
  const paper::outer::ZeroContraction<W> c = paper::outer::contract_zero_edges(g, s, &diag);
  const std::uint64_t big_n = c.graph.num_vertices();

  std::vector<W> contracted_dist;
  if (!params.outer_applies(big_n)) {
    // Alg 8 l.4-5 under the parameterised threshold: the recorded base case.
    SsspResult<W> base = paper::outer::outer_base_case(c.graph, c.source);
    paper::outer::absorb_diagnostics(diag, std::move(base.diag));
    contracted_dist = std::move(base.dist);
  } else {
    // Alg 8 l.8-15: the shared outer construction. Exceptions propagate (failure rule).
    const paper::outer::OuterConstructionParams p = outer_construction_params(params, big_n);
    paper::model::WordRng rng(seed);
    contracted_dist = paper::outer::outer_construction(c.graph, c.source, rng, p, &diag);
  }

  // Alg 8 l.6 / l.16-17: Lemma 8.5. It rejects labels that no tight original arc attains.
  SsspResult<W> r = paper::outer::restore_predecessor_tree(g, s, c, contracted_dist);
  paper::outer::absorb_diagnostics(diag, std::move(r.diag));
  r.diag = std::move(diag);
  return r;
}

template SsspResult<std::uint64_t> paper_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                              const PaperParams&, std::uint64_t);
template SsspResult<double> paper_sssp(const CsrGraph<double>&, Vertex, const PaperParams&,
                                       std::uint64_t);

}  // namespace sssp::paper_experimental
