#include "paper/sssp.hpp"

#include <cstdint>
#include <utility>
#include <vector>

#include "paper/model/random.hpp"
#include "paper/outer/base_case.hpp"
#include "paper/outer/contract_restore.hpp"
#include "paper/outer/shortest_paths.hpp"
#include "sssp/paper.hpp"

namespace sssp::paper {

template <class W>
std::vector<W> outer_construction(const outer::ZeroContraction<W>& c, std::uint64_t seed,
                                  Diagnostics& diag) {
  // Alg 8 l.8-15 with the paper's constants. The generator serves NameRoutes (l.9) only.
  model::WordRng rng(seed);
  return outer::outer_construction(c.graph, c.source, rng, outer::kPaperOuterParams, &diag);
}

template <class W>
SsspResult<W> paper_sssp(const CsrGraph<W>& g, Vertex s, PaperOptions options) {
  Diagnostics diag;
  // Alg 8 l.2: contract zero-edge components, discard loops, retain the source component, keep
  // original edge ids. Throws std::out_of_range for an invalid source.
  const outer::ZeroContraction<W> c = outer::contract_zero_edges(g, s, &diag);
  // Alg 8 l.3: N = number of retained vertices.
  const std::uint64_t big_n = c.graph.num_vertices();

  std::vector<W> contracted_dist;
  if (!takes_outer_construction(big_n)) {
    // Alg 8 l.4-5: base case on the contracted graph; it records the FallbackRecord with the
    // contracted size (N, m').
    SsspResult<W> base = options.variant == PaperVariant::kFaithful
                             ? outer::outer_base_case(c.graph, c.source)
                             : outer::outer_base_case_d1(c.graph, c.source);
    outer::absorb_diagnostics(diag, std::move(base.diag));
    contracted_dist = std::move(base.dist);
  } else {
    // Alg 8 l.8-15. Unreachable for N < 2^64 (Q-04).
    contracted_dist = outer_construction(c, options.seed, diag);
  }

  // Alg 8 l.6 (and l.16): Lemma 8.5, original distances and the first-discovery tree of finite
  // tight original arcs.
  SsspResult<W> r = outer::restore_predecessor_tree(g, s, c, contracted_dist);
  outer::absorb_diagnostics(diag, std::move(r.diag));
  r.diag = std::move(diag);
  return r;
}

template std::vector<std::uint64_t> outer_construction(const outer::ZeroContraction<std::uint64_t>&,
                                                       std::uint64_t, Diagnostics&);
template std::vector<double> outer_construction(const outer::ZeroContraction<double>&,
                                                std::uint64_t, Diagnostics&);
template SsspResult<std::uint64_t> paper_sssp(const CsrGraph<std::uint64_t>&, Vertex, PaperOptions);
template SsspResult<double> paper_sssp(const CsrGraph<double>&, Vertex, PaperOptions);

}  // namespace sssp::paper
