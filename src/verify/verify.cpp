// Shortest-path-tree verification. See include/sssp/verify.hpp and docs/design/numerics.md.

#include "sssp/verify.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace sssp {
namespace {

constexpr double kUnitRoundoff = 0x1.0p-53;  // u for IEEE-754 binary64, round to nearest

template <class W>
std::string fmt(W x) {
  if (x == WeightTraits<W>::infinity()) {
    return "inf";
  }
  if constexpr (std::is_floating_point_v<W>) {
    std::array<char, 32> buf{};
    std::snprintf(buf.data(), buf.size(), "%.17g", static_cast<double>(x));
    return buf.data();
  } else {
    return std::to_string(x);
  }
}

std::string vtx(Vertex v) {
  return "vertex " + std::to_string(v);
}

/// a == b, exactly for integers, within `tol` for floating point.
template <class W>
bool same(W a, W b, const Tolerance& tol) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    if (a == b) {
      return true;  // covers inf == inf
    }
    if (std::isinf(a) || std::isinf(b)) {
      return false;
    }
    return std::fabs(a - b) <= tol.slack(std::max(std::fabs(a), std::fabs(b)));
  } else {
    return a == b;
  }
}

/// a <= b, exactly for integers, within `tol` for floating point.
template <class W>
bool at_most(W a, W b, const Tolerance& tol) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    if (a <= b) {
      return true;
    }
    if (std::isinf(a) || std::isinf(b)) {
      return false;
    }
    return a - b <= tol.slack(std::max(std::fabs(a), std::fabs(b)));
  } else {
    return a <= b;
  }
}

/// a + b for finite distance a and weight b; false on integer overflow (sum >= infinity()).
template <class W>
bool add(W a, W b, W& out) noexcept {
  if constexpr (std::is_integral_v<W>) {
    if (a >= WeightTraits<W>::infinity() - b) {
      return false;
    }
    out = static_cast<W>(a + b);
  } else {
    out = a + b;
  }
  return true;
}

std::vector<char> bfs_reachable(const auto& g, Vertex source) {
  std::vector<char> seen(g.num_vertices(), 0);
  std::deque<Vertex> queue{source};
  seen[source] = 1;
  while (!queue.empty()) {
    const Vertex u = queue.front();
    queue.pop_front();
    for (const auto& a : g.arcs(u)) {
      if (seen[a.to] == 0) {
        seen[a.to] = 1;
        queue.push_back(a.to);
      }
    }
  }
  return seen;
}

/// Reports one failure per parent cycle. Vertices whose chain ends at a missing/invalid parent
/// are not reported here (that failure is reported by the per-vertex checks).
void check_parent_cycles(const std::vector<Vertex>& parent, Vertex source,
                         VerificationReport& rep) {
  enum : char { kUnseen, kOnPath, kDone };
  const std::size_t n = parent.size();
  std::vector<char> state(n, kUnseen);
  std::vector<Vertex> path;
  for (Vertex start = 0; start < n; ++start) {
    if (state[start] != kUnseen) {
      continue;
    }
    path.clear();
    Vertex x = start;
    while (x < n && x != source && state[x] == kUnseen) {
      state[x] = kOnPath;
      path.push_back(x);
      x = parent[x];
    }
    if (x < n && x != source && state[x] == kOnPath) {
      rep.add(FailureKind::kParentCycle, x,
              vtx(x) + " lies on a parent cycle that never reaches the source");
    }
    for (const Vertex p : path) {
      state[p] = kDone;
    }
  }
}

}  // namespace

Tolerance Tolerance::for_path_hops(std::uint64_t hops) {
  const double ku = static_cast<double>(hops) * kUnitRoundoff;
  if (!(ku < 0.5)) {
    throw std::invalid_argument(
        "Tolerance::for_path_hops: k*u must be < 1/2 (k = " + std::to_string(hops) + ")");
  }
  return Tolerance{ku / (1.0 - ku), 0.0};
}

std::string_view to_string(FailureKind kind) noexcept {
  switch (kind) {
    case FailureKind::kInvalidSource:
      return "invalid_source";
    case FailureKind::kSizeMismatch:
      return "size_mismatch";
    case FailureKind::kInvalidDistance:
      return "invalid_distance";
    case FailureKind::kSourceDistance:
      return "source_distance";
    case FailureKind::kSourceParent:
      return "source_parent";
    case FailureKind::kMissingParent:
      return "missing_parent";
    case FailureKind::kInvalidParent:
      return "invalid_parent";
    case FailureKind::kParentEdgeMissing:
      return "parent_edge_missing";
    case FailureKind::kParentEdgeMismatch:
      return "parent_edge_mismatch";
    case FailureKind::kTreeDistance:
      return "tree_distance";
    case FailureKind::kArithmeticOverflow:
      return "arithmetic_overflow";
    case FailureKind::kParentCycle:
      return "parent_cycle";
    case FailureKind::kUnreachableHasParent:
      return "unreachable_has_parent";
    case FailureKind::kUnreachableButFinite:
      return "unreachable_but_finite";
    case FailureKind::kReachableButInfinite:
      return "reachable_but_infinite";
    case FailureKind::kNotOptimal:
      return "not_optimal";
    case FailureKind::kReferenceSizeMismatch:
      return "reference_size_mismatch";
    case FailureKind::kReferenceMismatch:
      return "reference_mismatch";
    case FailureKind::kCount:
      break;
  }
  return "unknown";
}

void VerificationReport::add(FailureKind kind, Vertex vertex, std::string message) {
  ok = false;
  ++total_failures;
  ++counts[static_cast<std::size_t>(kind)];
  if (failures.size() < kMaxRecordedFailures) {
    failures.push_back(Failure{kind, vertex, std::move(message)});
  }
}

std::string VerificationReport::summary() const {
  if (ok) {
    return "OK";
  }
  std::ostringstream os;
  os << total_failures << " failure(s):";
  for (std::size_t k = 0; k < counts.size(); ++k) {
    if (counts[k] != 0) {
      os << ' ' << to_string(static_cast<FailureKind>(k)) << '=' << counts[k];
    }
  }
  os << '\n';
  for (const auto& f : failures) {
    os << "  [" << to_string(f.kind) << "] " << f.message << '\n';
  }
  if (failures.size() < total_failures) {
    os << "  ... " << (total_failures - failures.size()) << " more not recorded\n";
  }
  return os.str();
}

template <class W>
VerificationReport verify_sssp(const CsrGraph<W>& g, Vertex source, const SsspResult<W>& result,
                               const std::vector<std::type_identity_t<W>>* reference,
                               Tolerance tol) {
  VerificationReport rep;
  const std::size_t n = g.num_vertices();
  const W inf = WeightTraits<W>::infinity();
  const auto& dist = result.dist;
  const auto& parent = result.parent;
  const auto& pedge = result.parent_edge;

  if (source >= n) {
    rep.add(
        FailureKind::kInvalidSource, source,
        "source " + std::to_string(source) + " is not a vertex (n = " + std::to_string(n) + ")");
    return rep;
  }
  if (dist.size() != n || parent.size() != n || pedge.size() != n) {
    rep.add(FailureKind::kSizeMismatch, kNoVertex,
            "expected " + std::to_string(n) + " entries, got dist=" + std::to_string(dist.size()) +
                " parent=" + std::to_string(parent.size()) +
                " parent_edge=" + std::to_string(pedge.size()));
    return rep;
  }

  // Values that later checks can rely on: finite, non-negative (or inf).
  std::vector<char> valid(n, 1);
  if constexpr (std::is_floating_point_v<W>) {
    for (Vertex v = 0; v < n; ++v) {
      if (std::isnan(dist[v]) || dist[v] < W{0}) {
        valid[v] = 0;
        rep.add(FailureKind::kInvalidDistance, v, vtx(v) + " has distance " + fmt(dist[v]));
      }
    }
  }

  if (valid[source] != 0 && dist[source] != WeightTraits<W>::zero()) {
    rep.add(FailureKind::kSourceDistance, source, "dist[source] = " + fmt(dist[source]));
  }
  if (parent[source] != kNoVertex || pedge[source] != kNoEdge) {
    rep.add(FailureKind::kSourceParent, source, "the source has a parent");
  }

  const auto reachable = bfs_reachable(g, source);

  // Per-vertex tree checks.
  for (Vertex v = 0; v < n; ++v) {
    if (v == source || valid[v] == 0) {
      continue;
    }
    if (dist[v] == inf) {
      if (parent[v] != kNoVertex || pedge[v] != kNoEdge) {
        rep.add(FailureKind::kUnreachableHasParent, v, vtx(v) + " has dist = inf but a parent");
      }
      if (reachable[v] != 0) {
        rep.add(FailureKind::kReachableButInfinite, v,
                vtx(v) + " is reachable from the source but has dist = inf");
      }
      continue;
    }
    if (reachable[v] == 0) {
      rep.add(FailureKind::kUnreachableButFinite, v,
              vtx(v) + " is unreachable but has dist = " + fmt(dist[v]));
    }
    const Vertex p = parent[v];
    const EdgeId e = pedge[v];
    if (p == kNoVertex || e == kNoEdge) {
      rep.add(FailureKind::kMissingParent, v,
              vtx(v) + " has finite dist = " + fmt(dist[v]) + " but no parent / parent edge");
      continue;
    }
    if (p >= n) {
      rep.add(FailureKind::kInvalidParent, v,
              vtx(v) + " has parent " + std::to_string(p) + " which is not a vertex");
      continue;
    }
    if (e >= g.num_edges()) {
      rep.add(FailureKind::kParentEdgeMissing, v,
              vtx(v) + " has parent edge " + std::to_string(e) +
                  " which does not exist (m = " + std::to_string(g.num_edges()) + ")");
      continue;
    }
    const auto& edge = g.edge(e);
    if (!((edge.u == p && edge.v == v) || (edge.u == v && edge.v == p))) {
      rep.add(FailureKind::kParentEdgeMismatch, v,
              vtx(v) + ": parent edge " + std::to_string(e) + " = {" + std::to_string(edge.u) +
                  "," + std::to_string(edge.v) + "} does not join parent " + std::to_string(p) +
                  " and " + std::to_string(v));
      continue;
    }
    if (valid[p] == 0 || dist[p] == inf) {
      rep.add(FailureKind::kTreeDistance, v,
              vtx(v) + ": parent " + std::to_string(p) + " has dist = " + fmt(dist[p]));
      continue;
    }
    W via{};
    if (!add(dist[p], edge.weight, via)) {
      rep.add(FailureKind::kArithmeticOverflow, v,
              vtx(v) + ": dist[parent] + w = " + fmt(dist[p]) + " + " + fmt(edge.weight) +
                  " overflows");
      continue;
    }
    if (!same(dist[v], via, tol)) {
      rep.add(FailureKind::kTreeDistance, v,
              vtx(v) + ": dist = " + fmt(dist[v]) + " but dist[parent " + std::to_string(p) +
                  "] + w = " + fmt(dist[p]) + " + " + fmt(edge.weight) + " = " + fmt(via));
    }
  }

  check_parent_cycles(parent, source, rep);

  // Optimality certificate: dist[v] <= dist[u] + w on every arc out of a finite vertex.
  for (Vertex u = 0; u < n; ++u) {
    if (valid[u] == 0 || dist[u] == inf) {
      continue;
    }
    for (const auto& a : g.arcs(u)) {
      const Vertex v = a.to;
      if (valid[v] == 0 || dist[v] == inf) {
        continue;  // already reported above
      }
      W via{};
      if (!add(dist[u], a.weight, via)) {
        continue;  // bound exceeds the range: vacuous
      }
      if (!at_most(dist[v], via, tol)) {
        rep.add(FailureKind::kNotOptimal, v,
                vtx(v) + ": dist = " + fmt(dist[v]) + " > dist[" + std::to_string(u) +
                    "] + w = " + fmt(dist[u]) + " + " + fmt(a.weight) + " = " + fmt(via) +
                    " (edge " + std::to_string(a.edge) + ")");
      }
    }
  }

  if (reference != nullptr) {
    if (reference->size() != n) {
      rep.add(FailureKind::kReferenceSizeMismatch, kNoVertex,
              "reference has " + std::to_string(reference->size()) + " entries, expected " +
                  std::to_string(n));
    } else {
      for (Vertex v = 0; v < n; ++v) {
        if (!same(dist[v], (*reference)[v], tol)) {
          rep.add(FailureKind::kReferenceMismatch, v,
                  vtx(v) + ": dist = " + fmt(dist[v]) + " but reference = " + fmt((*reference)[v]));
        }
      }
    }
  }
  return rep;
}

template <class W>
void write_result_json(std::ostream& os, Vertex source, const SsspResult<W>& result) {
  const auto id = [](auto x, auto none) {
    return x == none ? std::string("-1") : std::to_string(x);
  };
  os << "{\"n\":" << result.dist.size() << ",\"source\":" << source << ",\"dist\":[";
  for (std::size_t i = 0; i < result.dist.size(); ++i) {
    if (i != 0) {
      os << ',';
    }
    const W d = result.dist[i];
    os << (d == WeightTraits<W>::infinity() ? std::string("null") : fmt(d));
  }
  os << "],\"parent\":[";
  for (std::size_t i = 0; i < result.parent.size(); ++i) {
    os << (i != 0 ? "," : "") << id(result.parent[i], kNoVertex);
  }
  os << "],\"parent_edge\":[";
  for (std::size_t i = 0; i < result.parent_edge.size(); ++i) {
    os << (i != 0 ? "," : "") << id(result.parent_edge[i], kNoEdge);
  }
  os << "]}\n";
}

template VerificationReport verify_sssp(const CsrGraph<std::uint64_t>&, Vertex,
                                        const SsspResult<std::uint64_t>&,
                                        const std::vector<std::uint64_t>*, Tolerance);
template VerificationReport verify_sssp(const CsrGraph<double>&, Vertex, const SsspResult<double>&,
                                        const std::vector<double>*, Tolerance);
template void write_result_json(std::ostream&, Vertex, const SsspResult<std::uint64_t>&);
template void write_result_json(std::ostream&, Vertex, const SsspResult<double>&);

}  // namespace sssp
