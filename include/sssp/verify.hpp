#pragma once
// Independent verification of single-source shortest-path results.
//
// verify_sssp does not trust the algorithm that produced the result: it checks the
// shortest-path-tree structure and an optimality certificate directly against the graph.
// Without a reference, a result is correct iff
//   (1) dist[source] = 0 and the parent pointers form a tree rooted at source whose edges are
//       real graph edges with dist[v] = dist[parent(v)] + w(parent edge)  -> dist[v] is the
//       length of a real path, so dist[v] >= d*(v); and
//   (2) dist[v] <= dist[u] + w(u,v) on every arc                         -> dist[v] <= d*(v)
//       (induction along a shortest path from the source).
// Reachability (dist = infinity <=> unreachable) is checked by an unweighted BFS.
//
// Integer weights are checked exactly (Tolerance is ignored, overflow is reported).
// Floating-point weights are checked with a relative + absolute tolerance; see
// docs/design/numerics.md for its definition and justification.

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp {

/// Allowed deviation between two double distance values a and b:
///   |a - b| <= rel * max(|a|, |b|) + abs.
/// One-sided checks (a <= b) allow a <= b + rel * max(|a|, |b|) + abs.
struct Tolerance {
  double rel = 0.0;
  double abs = 0.0;

  /// No deviation: exact IEEE-754 equality / ordering.
  static constexpr Tolerance exact() noexcept { return {}; }

  /// gamma_k = k u / (1 - k u) with u = 2^-53: the relative error bound of a sum of k+1
  /// non-negative doubles evaluated in any order (Higham, Accuracy and Stability of Numerical
  /// Algorithms, 2nd ed., Sec. 4.2). Throws std::invalid_argument if k u >= 1/2.
  static Tolerance for_path_hops(std::uint64_t hops);

  /// for_path_hops(n - 1): every shortest path is simple, so it has at most n-1 edges.
  template <class W>
  static Tolerance for_graph(const CsrGraph<W>& g) {
    const std::size_t n = g.num_vertices();
    return for_path_hops(n == 0 ? 0 : n - 1);
  }

  [[nodiscard]] double slack(double magnitude) const noexcept { return rel * magnitude + abs; }
};

enum class FailureKind : std::uint8_t {
  kInvalidSource,          ///< source is not a vertex
  kSizeMismatch,           ///< dist/parent/parent_edge do not have n entries
  kInvalidDistance,        ///< NaN or negative distance
  kSourceDistance,         ///< dist[source] != 0
  kSourceParent,           ///< source has a parent
  kMissingParent,          ///< finite distance (v != source) without parent or parent edge
  kInvalidParent,          ///< parent vertex out of range
  kParentEdgeMissing,      ///< parent edge id out of range
  kParentEdgeMismatch,     ///< parent edge does not join parent[v] and v
  kTreeDistance,           ///< dist[v] != dist[parent] + w(parent edge) (within tolerance)
  kArithmeticOverflow,     ///< integer dist[parent] + w overflows
  kParentCycle,            ///< following parents never reaches the source
  kUnreachableHasParent,   ///< infinite distance but a parent is recorded
  kUnreachableButFinite,   ///< finite distance for a vertex BFS cannot reach
  kReachableButInfinite,   ///< infinite distance for a vertex BFS reaches
  kNotOptimal,             ///< some arc (u,v) has dist[v] > dist[u] + w (within tolerance)
  kReferenceSizeMismatch,  ///< reference has the wrong length
  kReferenceMismatch,      ///< dist[v] differs from reference[v] (within tolerance)
  kCount
};

[[nodiscard]] std::string_view to_string(FailureKind kind) noexcept;

struct Failure {
  FailureKind kind;
  Vertex vertex;  ///< offending vertex, or kNoVertex for whole-result failures
  std::string message;
};

struct VerificationReport {
  static constexpr std::size_t kMaxRecordedFailures = 256;

  bool ok = true;
  std::vector<Failure> failures;  ///< the first kMaxRecordedFailures failures, in check order
  std::size_t total_failures = 0;
  std::array<std::size_t, static_cast<std::size_t>(FailureKind::kCount)> counts{};

  void add(FailureKind kind, Vertex vertex, std::string message);
  [[nodiscard]] bool has(FailureKind kind) const noexcept { return count(kind) > 0; }
  [[nodiscard]] std::size_t count(FailureKind kind) const noexcept {
    return counts[static_cast<std::size_t>(kind)];
  }
  /// Human-readable multi-line summary ("OK" or counts per kind followed by messages).
  [[nodiscard]] std::string summary() const;
};

/// Verifies `result` as a single-source shortest-path solution on g from `source`.
/// If `reference` is non-null (it is a non-deduced parameter, so `nullptr` works), distances must
/// also agree with it. Instantiated for std::uint64_t (exact; `tol` ignored) and double.
template <class W>
VerificationReport verify_sssp(const CsrGraph<W>& g, Vertex source, const SsspResult<W>& result,
                               const std::vector<std::type_identity_t<W>>* reference = nullptr,
                               Tolerance tol = Tolerance::exact());

/// Writes {"n":..,"source":..,"dist":[..],"parent":[..],"parent_edge":[..]} on one line.
/// Doubles are printed with %.17g (round-trips exactly), infinity as null, kNoVertex/kNoEdge
/// as -1. Read by tests/python/reference/check_double.py for exact rational checking.
template <class W>
void write_result_json(std::ostream& os, Vertex source, const SsspResult<W>& result);

}  // namespace sssp
