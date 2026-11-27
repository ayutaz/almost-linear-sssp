#pragma once
// Common graph representation shared by every algorithm (baseline, paper, experimental),
// so that benchmarks compare algorithms on identical memory layouts.
//
// Graphs are UNDIRECTED with non-negative weights. Each undirected edge {u,v} is stored once
// in the edge list and as two arcs (u->v, v->u) in the CSR; a self-loop {u,u} is one arc.

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace sssp {

using Vertex = std::uint32_t;
using EdgeId = std::uint32_t;

inline constexpr Vertex kNoVertex = std::numeric_limits<Vertex>::max();
inline constexpr EdgeId kNoEdge = std::numeric_limits<EdgeId>::max();

/// Per-weight-type semantics. Integer weights are exact (overflow is detected);
/// floating-point weights are subject to IEEE-754 rounding (see docs/design/numerics.md).
template <class W, class Enable = void>
struct WeightTraits;

template <class W>
struct WeightTraits<W, std::enable_if_t<std::is_integral_v<W> && std::is_unsigned_v<W>>> {
  static constexpr bool is_exact = true;
  static constexpr W infinity() noexcept { return std::numeric_limits<W>::max(); }
  static constexpr W zero() noexcept { return W{0}; }
  static constexpr bool is_valid_weight(W) noexcept { return true; }
  /// Sum of two finite distances; throws if the result would reach infinity().
  static W checked_add(W a, W b) {
    if (a >= infinity() - b) throw std::overflow_error("integer distance overflow");
    return static_cast<W>(a + b);
  }
};

template <class W>
struct WeightTraits<W, std::enable_if_t<std::is_floating_point_v<W>>> {
  static constexpr bool is_exact = false;
  static constexpr W infinity() noexcept { return std::numeric_limits<W>::infinity(); }
  static constexpr W zero() noexcept { return W{0}; }
  static bool is_valid_weight(W w) noexcept { return std::isfinite(w) && w >= W{0}; }
  static W checked_add(W a, W b) noexcept { return a + b; }
};

template <class W>
struct Edge {
  Vertex u;
  Vertex v;
  W weight;
};

template <class W>
struct Arc {
  Vertex to;
  W weight;
  EdgeId edge;  ///< index into the original edge list
};

/// Validates an edge list; throws std::out_of_range for bad endpoints and
/// std::invalid_argument for negative / NaN / infinite weights.
template <class W>
void validate_edges(std::size_t n, std::span<const Edge<W>> edges) {
  if (n >= static_cast<std::size_t>(kNoVertex)) throw std::length_error("too many vertices");
  if (edges.size() >= static_cast<std::size_t>(kNoEdge)) throw std::length_error("too many edges");
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u >= n || e.v >= n) {
      throw std::out_of_range("edge " + std::to_string(i) + " has an endpoint out of range");
    }
    if constexpr (std::is_signed_v<W> || std::is_floating_point_v<W>) {
      if (!WeightTraits<W>::is_valid_weight(e.weight)) {
        throw std::invalid_argument("edge " + std::to_string(i) +
                                    " has a negative, NaN or infinite weight");
      }
    }
  }
}

template <class W>
class CsrGraph {
 public:
  using weight_type = W;

  CsrGraph() = default;

  static CsrGraph from_edges(std::size_t n, std::vector<Edge<W>> edges) {
    validate_edges<W>(n, edges);
    CsrGraph g;
    g.edges_ = std::move(edges);
    g.offsets_.assign(n + 1, 0);
    for (const auto& e : g.edges_) {
      ++g.offsets_[e.u + 1];
      if (e.u != e.v) ++g.offsets_[e.v + 1];
    }
    for (std::size_t v = 0; v < n; ++v) g.offsets_[v + 1] += g.offsets_[v];
    g.arcs_.resize(g.offsets_[n]);
    std::vector<std::size_t> cursor(g.offsets_.begin(), g.offsets_.end() - 1);
    for (std::size_t i = 0; i < g.edges_.size(); ++i) {
      const auto& e = g.edges_[i];
      const auto id = static_cast<EdgeId>(i);
      g.arcs_[cursor[e.u]++] = Arc<W>{e.v, e.weight, id};
      if (e.u != e.v) g.arcs_[cursor[e.v]++] = Arc<W>{e.u, e.weight, id};
    }
    return g;
  }

  [[nodiscard]] std::size_t num_vertices() const noexcept {
    return offsets_.empty() ? 0 : offsets_.size() - 1;
  }
  [[nodiscard]] std::size_t num_edges() const noexcept { return edges_.size(); }
  [[nodiscard]] std::size_t num_arcs() const noexcept { return arcs_.size(); }
  [[nodiscard]] std::size_t degree(Vertex v) const noexcept {
    return offsets_[v + 1] - offsets_[v];
  }
  [[nodiscard]] std::span<const Arc<W>> arcs(Vertex v) const noexcept {
    return {arcs_.data() + offsets_[v], arcs_.data() + offsets_[v + 1]};
  }
  [[nodiscard]] const Edge<W>& edge(EdgeId id) const noexcept { return edges_[id]; }
  [[nodiscard]] std::span<const Edge<W>> edges() const noexcept { return edges_; }

 private:
  std::vector<std::size_t> offsets_;
  std::vector<Arc<W>> arcs_;
  std::vector<Edge<W>> edges_;
};

}  // namespace sssp
