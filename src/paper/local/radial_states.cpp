// Radial state representation, Sec 5.2 (Def 5.4, Lemmas 5.5-5.7). See radial_states.hpp.

#include "paper/local/radial_states.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "paper/local/saturating.hpp"
#include "paper/model/ilog.hpp"
#include "paper/primitives/dsu.hpp"

namespace sssp::paper {

namespace {

/// Up to two copies of one phase of a vertex (native first).
struct Copies {
  std::array<StateId, 2> s{kNoState, kNoState};
  std::size_t size = 0;
  void push(StateId x) {
    if (x != kNoState) s[size++] = x;
  }
  [[nodiscard]] const StateId* begin() const { return s.data(); }
  [[nodiscard]] const StateId* end() const { return s.data() + size; }
};

template <class W>
class Builder {
 public:
  Builder(const CsrGraph<W>& g, const std::vector<W>& coordinate,
          const EdgeClassification& classification, unsigned k)
      : g_(g), edges_(g.edges()), n_(g.num_vertices()), k_(k) {
    out_.coordinate = coordinate;
    out_.classification = classification;
    out_.shell_ratio_log2 = k;
    positive_.resize(n_);
    for (Vertex v = 0; v < n_; ++v) positive_[v] = W{0} < coordinate[v] ? 1 : 0;
  }

  RadialStateGraph<W> run() {
    near_components_and_shells();
    create_states();
    create_arcs();
    create_blocks();
    return std::move(out_);
  }

 private:
  [[nodiscard]] const std::vector<W>& coord() const { return out_.coordinate; }

  /// Near radial edge (l.1030-1032): radial, both coordinates positive, ratio <= K^2.
  [[nodiscard]] bool near(EdgeId e) const {
    const auto& ed = edges_[e];
    if (out_.classification.edge_class[e] != EdgeClass::kRadial || ed.u == ed.v) return false;
    if (positive_[ed.u] == 0 || positive_[ed.v] == 0) return false;
    const bool u_low = !(coord()[ed.v] < coord()[ed.u]);
    const Vertex lo = u_low ? ed.u : ed.v;
    const Vertex hi = u_low ? ed.v : ed.u;
    return !(k2r_[lo] < coord()[hi]);
  }

  // l.1030-1038: near components, anchors a, and native shells a K^i <= r < a K^(i+1).
  void near_components_and_shells() {
    k2r_.resize(n_);
    for (Vertex v = 0; v < n_; ++v) k2r_[v] = sat::mul_pow2(coord()[v], 2 * k_);
    DisjointSets dsu(n_);
    std::vector<std::size_t> deg(n_ + 1, 0);
    for (EdgeId e = 0; e < edges_.size(); ++e) {
      if (!near(e)) continue;
      dsu.unite(edges_[e].u, edges_[e].v);
      ++deg[edges_[e].u + 1];
      ++deg[edges_[e].v + 1];
    }
    for (std::size_t v = 0; v < n_; ++v) deg[v + 1] += deg[v];
    std::vector<Vertex> nbr(deg[n_]);
    std::vector<std::size_t> pos(deg.begin(), deg.end() - 1);
    for (EdgeId e = 0; e < edges_.size(); ++e) {
      if (!near(e)) continue;
      nbr[pos[edges_[e].u]++] = edges_[e].v;
      nbr[pos[edges_[e].v]++] = edges_[e].u;
    }

    // Component ids in order of their smallest vertex; the anchor vertex is a least coordinate
    // (smallest id on ties).
    out_.near_component.assign(n_, kNoNearComponent);
    out_.shell.assign(n_, 0);
    std::vector<std::uint32_t> id_of_root(n_, kNoNearComponent);
    std::vector<Vertex> anchor_vertex;
    for (Vertex v = 0; v < n_; ++v) {
      if (positive_[v] == 0) continue;
      const auto r = dsu.find(v);
      if (id_of_root[r] == kNoNearComponent) {
        id_of_root[r] = static_cast<std::uint32_t>(anchor_vertex.size());
        anchor_vertex.push_back(v);
      }
      const auto c = id_of_root[r];
      out_.near_component[v] = c;
      if (coord()[v] < coord()[anchor_vertex[c]]) anchor_vertex[c] = v;
    }
    out_.component_anchor.resize(anchor_vertex.size());

    // Walk each component from its minimum; adjacent shells differ by at most two (l.1036),
    // so each step compares against O(1) cached powers a K^j.
    std::vector<std::uint8_t> seen(n_, 0);
    std::vector<Vertex> queue;
    std::vector<W> pw;
    for (std::uint32_t c = 0; c < anchor_vertex.size(); ++c) {
      const Vertex s = anchor_vertex[c];
      out_.component_anchor[c] = coord()[s];
      pw.assign(1, coord()[s]);
      const auto power = [&](std::size_t j) -> W {
        while (pw.size() <= j) pw.push_back(sat::mul_pow2(pw.back(), k_));
        return pw[j];
      };
      queue.assign(1, s);
      seen[s] = 1;
      for (std::size_t qi = 0; qi < queue.size(); ++qi) {
        const Vertex x = queue[qi];
        for (std::size_t t = deg[x]; t < deg[x + 1]; ++t) {
          const Vertex y = nbr[t];
          if (seen[y] != 0) continue;
          seen[y] = 1;
          const std::uint32_t ix = out_.shell[x];
          std::uint32_t j = ix >= 2 ? ix - 2 : 0;
          while (!(coord()[y] < power(j + 1))) ++j;
          if (coord()[y] < power(j) || j > ix + 2) {
            throw std::logic_error("build_radial_states: near shells differ by more than two");
          }
          out_.shell[y] = j;
          queue.push_back(y);
        }
      }
    }
  }

  // l.1039-1042: two-shell copies.
  void create_states() {
    out_.native_inward.assign(n_, kNoState);
    out_.lower_inward.assign(n_, kNoState);
    out_.native_outward.assign(n_, kNoState);
    out_.upper_outward.assign(n_, kNoState);
    const auto add = [&](Vertex v, Phase ph, std::uint32_t window) {
      const auto id = static_cast<StateId>(out_.projection.size());
      out_.projection.push_back(v);
      out_.phase.push_back(ph);
      out_.window.push_back(window);
      return id;
    };
    for (Vertex v = 0; v < n_; ++v) {
      if (positive_[v] == 0) {  // a pivot: one inward and one outward copy (l.1051)
        out_.native_inward[v] = add(v, Phase::kInward, 0);
        out_.native_outward[v] = add(v, Phase::kOutward, 0);
        continue;
      }
      const std::uint32_t i = out_.shell[v];
      out_.native_inward[v] = add(v, Phase::kInward, i);
      if (i >= 1) out_.lower_inward[v] = add(v, Phase::kInward, i - 1);
      out_.native_outward[v] = add(v, Phase::kOutward, i);
      out_.upper_outward[v] = add(v, Phase::kOutward, i + 1);
    }
    out_.graph.num_states = out_.projection.size();
  }

  [[nodiscard]] Copies inward(Vertex v) const {
    Copies c;
    c.push(out_.native_inward[v]);
    c.push(out_.lower_inward[v]);
    return c;
  }
  [[nodiscard]] Copies outward(Vertex v) const {
    Copies c;
    c.push(out_.native_outward[v]);
    c.push(out_.upper_outward[v]);
    return c;
  }
  /// The copy of positive vertex v in inward (outward) window j, or kNoState.
  [[nodiscard]] StateId inward_at(Vertex v, std::uint32_t j) const {
    if (j == out_.shell[v]) return out_.native_inward[v];
    if (j + 1 == out_.shell[v]) return out_.lower_inward[v];
    return kNoState;
  }
  [[nodiscard]] StateId outward_at(Vertex v, std::uint32_t j) const {
    if (j == out_.shell[v]) return out_.native_outward[v];
    if (j == out_.shell[v] + 1) return out_.upper_outward[v];
    return kNoState;
  }

  void arc(StateId t, StateId h, W len, RadialArcKind kind, EdgeId e) {
    out_.graph.arcs.push_back({t, h, len});
    out_.arc_kind.push_back(kind);
    out_.arc_edge.push_back(e);
  }

  // l.1042-1055: arcs.
  void create_arcs() {
    using K = RadialArcKind;
    for (EdgeId e = 0; e < edges_.size(); ++e) {
      const auto& ed = edges_[e];
      const EdgeClass cls = out_.classification.edge_class[e];
      if (cls == EdgeClass::kOmitted || ed.u == ed.v) continue;  // loops carry no arcs
      const W w = ed.weight;
      if (cls == EdgeClass::kMacro) {
        // Every copy of the tail resets to the native inward copy of the head (l.1053-1054).
        for (const auto [x, y] : {std::pair{ed.u, ed.v}, std::pair{ed.v, ed.u}}) {
          for (const StateId c : inward(x)) arc(c, out_.native_inward[y], w, K::kMacro, e);
          for (const StateId c : outward(x)) arc(c, out_.native_inward[y], w, K::kMacro, e);
        }
        continue;
      }
      const bool zu = positive_[ed.u] == 0;
      const bool zv = positive_[ed.v] == 0;
      if (zu && zv) {
        throw std::invalid_argument("build_radial_states: radial edge " + std::to_string(e) +
                                    " joins two zero-coordinate vertices");
      }
      if (zu || zv) {
        // Pivot radial edge (l.1051-1052; Q-32: from every inward copy of the positive end).
        const Vertex p = zu ? ed.u : ed.v;
        const Vertex v = zu ? ed.v : ed.u;
        for (const StateId c : inward(v)) arc(c, out_.native_inward[p], w, K::kPivot, e);
        arc(out_.native_outward[p], out_.native_outward[v], w, K::kPivot, e);
        continue;
      }
      if (near(e)) {
        for (const auto [x, y] : {std::pair{ed.u, ed.v}, std::pair{ed.v, ed.u}}) {
          // Inward copy at j: within window, or transfer down to a native shell i < j (l.1043).
          for (const StateId c : inward(x)) {
            const std::uint32_t j = out_.window[c];
            const StateId t = inward_at(y, j);
            if (t != kNoState) {
              arc(c, t, w, K::kWithinWindow, e);
            } else if (out_.shell[y] < j) {
              arc(c, out_.native_inward[y], w, K::kNearTransfer, e);
            }
          }
          // Outward copy at j: within window, or transfer up to a native shell i > j (l.1044).
          for (const StateId c : outward(x)) {
            const std::uint32_t j = out_.window[c];
            const StateId t = outward_at(y, j);
            if (t != kNoState) {
              arc(c, t, w, K::kWithinWindow, e);
            } else if (out_.shell[y] > j) {
              arc(c, out_.native_outward[y], w, K::kNearTransfer, e);
            }
          }
        }
        continue;
      }
      // Far positive radial edge (l.1045-1047).
      const bool u_high = coord()[ed.v] < coord()[ed.u];
      const Vertex hi = u_high ? ed.u : ed.v;
      const Vertex lo = u_high ? ed.v : ed.u;
      for (const StateId c : inward(hi)) arc(c, out_.native_inward[lo], w, K::kFarTransfer, e);
      for (const StateId c : outward(lo)) arc(c, out_.native_outward[hi], w, K::kFarTransfer, e);
    }
    // Phase switches (l.1049-1051): every inward copy to the native outward copy, length 0.
    for (Vertex v = 0; v < n_; ++v) {
      for (const StateId c : inward(v)) arc(c, out_.native_outward[v], W{0}, K::kSwitch, kNoEdge);
    }
  }

  // l.1041-1042: blocks are the components of the within-window arcs.
  void create_blocks() {
    const std::size_t s = out_.graph.num_states;
    DisjointSets dsu(s);
    for (std::size_t a = 0; a < out_.graph.arcs.size(); ++a) {
      if (out_.arc_kind[a] == RadialArcKind::kWithinWindow) {
        dsu.unite(out_.graph.arcs[a].tail, out_.graph.arcs[a].head);
      }
    }
    std::vector<std::uint32_t> id(s, std::numeric_limits<std::uint32_t>::max());
    out_.block.resize(s);
    for (StateId x = 0; x < s; ++x) {
      const auto r = dsu.find(x);
      if (id[r] == std::numeric_limits<std::uint32_t>::max()) {
        id[r] = static_cast<std::uint32_t>(out_.num_blocks++);
      }
      out_.block[x] = id[r];
    }
  }

  const CsrGraph<W>& g_;
  std::span<const Edge<W>> edges_;
  std::size_t n_;
  unsigned k_;
  std::vector<std::uint8_t> positive_;
  std::vector<W> k2r_;  ///< K^2 r_v (saturating), l.1032
  RadialStateGraph<W> out_;
};

}  // namespace

unsigned local_shell_ratio_log2(std::size_t hop_bound) {
  if (hop_bound == 0) throw std::invalid_argument("local_shell_ratio_log2: h must be >= 1");
  // 2^k >= 8h  <=>  k >= 3 + ceil(log2 h): integer word operations on h only.
  return 3 + model::ceil_log2(hop_bound);
}

template <class W>
RadialStateGraph<W> build_radial_states(const CsrGraph<W>& g, const std::vector<W>& coordinate,
                                        const EdgeClassification& classification,
                                        unsigned shell_ratio_log2, Diagnostics* diag) {
  const std::size_t n = g.num_vertices();
  const std::size_t m = g.num_edges();
  if (coordinate.size() != n) {
    throw std::invalid_argument("build_radial_states: one coordinate per vertex required");
  }
  if (classification.edge_class.size() != m || classification.macro_key.size() != m) {
    throw std::invalid_argument("build_radial_states: one class and key per edge required");
  }
  if (shell_ratio_log2 == 0) {
    throw std::invalid_argument("build_radial_states: the shell ratio K must be at least 2");
  }
  for (Vertex v = 0; v < n; ++v) {
    const W x = coordinate[v];
    bool ok = false;
    if constexpr (std::is_floating_point_v<W>) {
      ok = std::isfinite(x) && !(x < W{0});
    } else {
      ok = x != WeightTraits<W>::infinity();
    }
    if (!ok) {
      throw std::invalid_argument("build_radial_states: coordinate of vertex " + std::to_string(v) +
                                  " is not finite and nonnegative");
    }
  }
  auto rs = Builder<W>(g, coordinate, classification, shell_ratio_log2).run();
  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["radial_states_calls"]);
    SSSP_COUNT(diag->counters["radial_states_states"] += rs.graph.num_states);
    SSSP_COUNT(diag->counters["radial_states_arcs"] += rs.graph.arcs.size());
    SSSP_COUNT(diag->counters["radial_states_blocks"] += rs.num_blocks);
    SSSP_COUNT(diag->counters["radial_states_near_components"] += rs.component_anchor.size());
  }
  return rs;
}

template <class W>
EdgeClassification classify_local_edges(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                        unsigned macro_exponent) {
  const std::size_t n = g.num_vertices();
  const auto edges = g.edges();
  if (od.num_vertices() != n || od.reduced.size() != edges.size()) {
    throw std::invalid_argument("classify_local_edges: decomposition of a different graph");
  }
  // N^e sigma once per pivot-tree maximum (l.1026-1028), saturating (Q-38).
  std::vector<W> scaled(od.corridors.size());
  for (std::size_t c = 0; c < od.corridors.size(); ++c) {
    scaled[c] = sat::mul_int_pow(od.corridors[c].weight, n, macro_exponent);
  }
  EdgeClassification out{std::vector<EdgeClass>(edges.size(), EdgeClass::kOmitted),
                         std::vector<std::uint32_t>(edges.size(), kNoCorridor)};
  for (EdgeId e = 0; e < edges.size(); ++e) {
    const auto& ed = edges[e];
    if (od.reduced[e] == 0 || ed.u == ed.v) continue;
    out.edge_class[e] = EdgeClass::kRadial;
    const auto p = od.owner[ed.u];
    const auto q = od.owner[ed.v];
    if (p == q) continue;
    const std::uint32_t c = od.corridor_max(p, q);
    if (!(scaled[c] < ed.weight)) {  // w <= N^3 sigma(p, q): macro (l.1021)
      out.edge_class[e] = EdgeClass::kMacro;
      out.macro_key[e] = c;
    }
  }
  return out;
}

template <class W>
RadialStateGraph<W> local_radial_states(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                        Diagnostics* diag) {
  return build_radial_states(g, od.radius, classify_local_edges(g, od),
                             local_shell_ratio_log2(od.hop_bound), diag);
}

#define SSSP_INSTANTIATE_RADIAL(W)                                                               \
  template RadialStateGraph<W> build_radial_states<W>(const CsrGraph<W>&, const std::vector<W>&, \
                                                      const EdgeClassification&, unsigned,       \
                                                      Diagnostics*);                             \
  template EdgeClassification classify_local_edges<W>(const CsrGraph<W>&,                        \
                                                      const OwnerDecomposition<W>&, unsigned);   \
  template RadialStateGraph<W> local_radial_states<W>(const CsrGraph<W>&,                        \
                                                      const OwnerDecomposition<W>&, Diagnostics*);
SSSP_INSTANTIATE_RADIAL(std::uint64_t)
SSSP_INSTANTIATE_RADIAL(double)
#undef SSSP_INSTANTIATE_RADIAL

}  // namespace sssp::paper
