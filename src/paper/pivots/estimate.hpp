#pragma once
// Sec 6.2 "Estimating Distances in a Sample", Algorithm 6 Estimate(S, Q) and Lemmas 6.4-6.7
// (arXiv 2611.34911v1, txt l.1455-1600), Mode A.
//
// Input: a positive-weight, loopless sample graph S on t pivots (O(t) edges) and h demands with
// distinct endpoints. Output: a table of upper bounds and, per demand, a table index or a
// disconnection flag. Every table value has the form 4 Delta + a with Delta = 2^i a for a band
// anchor a (Alg 6 l.9); it is stored as the pair (a, i) so that it is never rounded or wrapped.
//
// Guarantees (Lemma 6.7): for every connected demand, d_S(x, y) <= Z deterministically (on every
// seed); E[Z | S] <= 256 u^2 d_S(x, y) with u = ceil(log2(t + 2)) (expectation over the
// estimator's randomness; measured in the tests, not assumed). Disconnected pairs are flagged.
//
// Constants (Sec 6.2). T = least power of two >= 4t, the top scale 2^l a >= t T^2 a and the
// window [a/T^2, T^2 a] enter Lemma 6.4 and Lemma 6.7(1), i.e. the correctness of Z >= d_S; they
// are fixed here and are NOT parameters, also in Mode B. H = least power of two >= (16t)^5 only
// enters E[Z] (Lemma 6.7(2); the proof needs H >= 4tT^2, research note Q-42), so Mode B may set
// it through EstimateParams::h_log2 (deviation D-9).
//
// Model of computation (Sec 2). Lengths are only compared, added and doubled. Power-of-two
// multiples (T a, T^2 w, 2^i a, H d, (H + j) Delta, 4 Delta + a) are formed by repeated doubling
// and addition, as the paper prescribes ("using shared doublings", "binary doubling and addition
// on the input lengths"). For integer weights these values are held exactly in multiword
// integers (research note Q-06 / Q-43: integer weights, exact power-of-two scalings); for double
// weights doubling is exact and additions round (IEEE-754).
//
// Randomness: only the caller's WordRng is used (uniform permutations and j in {0..H-1} by the
// rejection method of Sec 2), so a fixed seed reproduces the output on every platform.
// Las Vegas: the bound Z >= d_S does not depend on the random choices.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "paper/model/random.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::pivots {

/// A demand x-y with x != y (pivot indices).
struct Demand {
  Vertex x;
  Vertex y;
};

/// Index value of a demand whose endpoints lie in different sample components.
inline constexpr std::uint32_t kDisconnected = 0xffffffffu;

/// Table value Z = 4 Delta + anchor with Delta = 2^scale * anchor (Alg 6 l.9).
template <class W>
struct TableBound {
  W anchor;
  std::uint32_t scale;
};

/// Parameters of Estimate. Mode A uses the defaults (paper values). T, the top scale and the
/// window are deliberately not parameters (see the file comment).
struct EstimateParams {
  /// log2 H. nullopt: the paper's H = least power of two >= (16t)^5 (Sec 6.2, l.1491).
  /// Any other value is Mode B deviation D-9; Z >= d_S holds for every H >= 1.
  std::optional<unsigned> h_log2;
};

/// Integer constants of Sec 6.2 for t pivots: T = 2^tau (least power of two >= 4t), the top
/// scale index ell (least l with 2^l >= t T^2), and H = 2^eta.
struct EstimateConstants {
  unsigned tau = 0;
  unsigned ell = 0;
  unsigned eta = 0;
};

/// Throws std::invalid_argument for t = 0 and std::length_error for t >= 2^32.
EstimateConstants estimate_constants(std::size_t t, const EstimateParams& params = {});

template <class W>
struct EstimateResult {
  std::vector<TableBound<W>> table;
  /// Per demand: index into `table`, or kDisconnected.
  std::vector<std::uint32_t> index;
  /// Counters "estimate.*" (bands_processed, quotient_vertices, quotient_edges, vertex_scales,
  /// center_searches, settles, adjacency_scans, table_entries, disconnected_demands,
  /// max_windows_per_edge) and counted length operations.
  Diagnostics diag;
};

/// Optional record of the intermediate state, for invariant tests (Lemma 6.4, 6.6, 6.7, Q-41).
template <class W>
struct EstimateTrace {
  struct DemandTrace {
    W beta{};                            ///< weight of the merge (LCA in the merge forest)
    std::uint32_t band = kDisconnected;  ///< index into `bands`
    Vertex qx = kNoVertex;               ///< image of x in the band's quotient
    Vertex qy = kNoVertex;
  };
  struct ScaleTrace {
    std::uint32_t scale = 0;                 ///< i, Delta = 2^i a
    std::vector<std::uint32_t> order;        ///< center order (uniform permutation)
    std::vector<std::uint64_t> j;            ///< j in {0..H-1}, little-endian 64-bit words
    std::vector<Vertex> center;              ///< first covering center of each quotient vertex
    std::vector<std::uint32_t> table_index;  ///< table index of each vertex's part at this scale
  };
  struct BandTrace {
    W anchor{};
    std::vector<Vertex> image;            ///< per pivot: quotient vertex, or kNoVertex
    std::vector<Edge<W>> quotient_edges;  ///< retained edges on quotient vertices
    std::vector<std::uint32_t> retained;  ///< sample edge index of each quotient edge
    std::size_t num_quotient_vertices = 0;
    std::vector<std::uint32_t> component;  ///< quotient component of each quotient vertex
    std::vector<ScaleTrace> scales;        ///< i = ell (coarsest) down to 0
  };
  EstimateConstants constants;
  std::vector<DemandTrace> demand;
  std::vector<BandTrace> bands;  ///< bands with a nonempty demand group, increasing anchors
};

/// Algorithm 6. Throws std::out_of_range for an endpoint >= t, std::invalid_argument for a loop,
/// a nonpositive or non-finite weight (or the integer infinity sentinel) and for a demand with
/// equal endpoints; std::overflow_error if a sum of integer weights overflows.
template <class W>
EstimateResult<W> estimate(std::size_t t, std::span<const Edge<W>> sample,
                           std::span<const Demand> demands, model::WordRng& rng,
                           const EstimateParams& params = {}, EstimateTrace<W>* trace = nullptr);

// ---------------------------------------------------------- exact power-of-two scaled tests
// For integer W these are exact (multiword); for double W they use exact doublings and IEEE
// additions. They are what Estimate and NameRoutes use internally.

/// d <= Z = 4 * 2^scale * anchor + anchor.
template <class W>
bool bound_at_least(const TableBound<W>& z, W d);
/// lambda * 2^k >= Z (the median-pruning test P a* >= Z of Sec 6.1, with P = 2^k).
template <class W>
bool scaled_at_least(W lambda, unsigned k, const TableBound<W>& z);
/// lambda * 2^k >= d (route certificate test, Definition 6.1, with P = 2^k).
template <class W>
bool pow2_multiple_at_least(W lambda, unsigned k, W d);
/// The radius test of Sec 6.2: H d <= (H + j) Delta with H = 2^eta, Delta = 2^scale * anchor and
/// j given as little-endian words. Throws std::invalid_argument unless j < H.
template <class W>
bool within_radius(W d, W anchor, unsigned scale, unsigned eta, std::span<const std::uint64_t> j);

/// Output invariants of Lemma 6.7 against reference sample distances `dist` (t x t,
/// WeightTraits<W>::infinity() for unreachable pairs): sizes, every table anchor is a sample
/// weight, a demand is flagged iff its endpoints are disconnected, and d_S(x, y) <= Z for every
/// connected demand. Returns one message per violation (empty: all hold).
template <class W>
std::vector<std::string> check_estimate(std::size_t t, std::span<const Edge<W>> sample,
                                        std::span<const Demand> demands,
                                        const EstimateResult<W>& result,
                                        const std::vector<std::vector<W>>& dist);

}  // namespace sssp::paper::pivots
