// Mode B assembly of the local hierarchy (task P5-EXP-LOCAL): Sec 5.4 with the paper's or the
// D-8 alignment constants, or the D-11 exact-height assembly. See local_assembly.hpp.
//
// Steps 1-10 follow src/paper/local/assemble.cpp (the step numbers of assemble.hpp) line for
// line; the only differences are AlignmentRule (step 3 test, step 4 support end) and the
// exact-height branch that replaces steps 1-8. Keep the two files in step: the property test
// "Mode B assembly engine with the paper's constants equals Mode A on families" compares them.

#include "experimental/paper_experimental/local_assembly.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paper/local/saturating.hpp"
#include "paper/primitives/dsu.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper_experimental {

namespace {

using paper::ArcId;
using paper::BlockScale;
using paper::BlockScaleKind;
using paper::DiArc;
using paper::Digraph;
using paper::DisjointSets;
using paper::ForestQueries;
using paper::GridSupport;
using paper::JoinKind;
using paper::kNoComponent;
using paper::kNoCorridor;
using paper::kNoNode;
using paper::kNoState;
using paper::NodeId;
using paper::NodeOrigin;
using paper::StateId;
namespace sat = paper::sat;

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

/// Stable counting sort of item indices 0..n-1 by key(i) in [0, k) (Lemma 3.1). Returns the
/// bucket offsets (size k + 1) and the sorted items.
template <class Key>
std::pair<std::vector<std::size_t>, std::vector<std::uint32_t>> bucket(std::size_t n, std::size_t k,
                                                                       Key key) {
  std::vector<std::size_t> begin(k + 1, 0);
  for (std::size_t i = 0; i < n; ++i) ++begin[key(i) + 1];
  for (std::size_t b = 0; b < k; ++b) begin[b + 1] += begin[b];
  std::vector<std::uint32_t> items(n);
  std::vector<std::size_t> pos(begin.begin(), begin.end() - 1);
  for (std::size_t i = 0; i < n; ++i) items[pos[key(i)]++] = static_cast<std::uint32_t>(i);
  return {std::move(begin), std::move(items)};
}

/// Class-wise simultaneous merges of join records (Sec 3.3, Lemma 3.6) building the hierarchy
/// on the states, as in Mode A: the disjoint sets range over all join vertices; rep_[root] is
/// the hierarchy node of the component's states, kNoNode for a component without states
/// (discarded, l.1281); a merged component with one nonempty part keeps that part's node.
template <class W>
class MergeBuilder {
 public:
  MergeBuilder(paper::LocalAssembly<W>& out, const paper::ThresholdJoins<W>& joins)
      : out_(out),
        joins_(joins),
        dsu_(joins.num_join_vertices()),
        rep_(joins.num_join_vertices(), kNoNode),
        stamp_of_(joins.num_join_vertices(), 0),
        group_of_(joins.num_join_vertices(), 0) {
    for (std::size_t x = 0; x < joins.num_states; ++x) rep_[x] = static_cast<NodeId>(x);
  }

  /// One class: the join records `js` merged simultaneously; every new component with >= 2
  /// nonempty parts becomes a node of the given height.
  void merge_class(std::span<const std::uint32_t> js, W height, W upper, bool inclusive,
                   NodeOrigin origin) {
    if (js.empty()) return;
    ++stamp_;
    pre_.clear();
    for (const std::uint32_t j : js) {
      for (const std::uint32_t v : {joins_.joins[j].a, joins_.joins[j].b}) {
        const auto r = dsu_.find(v);
        if (stamp_of_[r] != stamp_) {
          stamp_of_[r] = stamp_;
          pre_.push_back(r);
        }
      }
    }
    for (const std::uint32_t j : js) dsu_.unite(joins_.joins[j].a, joins_.joins[j].b);
    ++stamp_;
    group_root_.clear();
    gid_.assign(pre_.size(), 0);
    for (std::size_t i = 0; i < pre_.size(); ++i) {
      const auto r = dsu_.find(pre_[i]);
      if (stamp_of_[r] != stamp_) {
        stamp_of_[r] = stamp_;
        group_of_[r] = static_cast<std::uint32_t>(group_root_.size());
        group_root_.push_back(r);
      }
      gid_[i] = group_of_[r];
    }
    const std::size_t ng = group_root_.size();
    std::vector<std::size_t> begin(ng + 1, 0);
    for (std::size_t i = 0; i < pre_.size(); ++i) {
      if (rep_[pre_[i]] != kNoNode) ++begin[gid_[i] + 1];
    }
    for (std::size_t g = 0; g < ng; ++g) begin[g + 1] += begin[g];
    std::vector<NodeId> parts(begin[ng]);
    {
      std::vector<std::size_t> pos(begin.begin(), begin.end() - 1);
      for (std::size_t i = 0; i < pre_.size(); ++i) {
        if (rep_[pre_[i]] != kNoNode) parts[pos[gid_[i]]++] = rep_[pre_[i]];
      }
    }
    auto& h = out_.result.hierarchy;
    for (std::size_t g = 0; g < ng; ++g) {
      const std::size_t k = begin[g + 1] - begin[g];
      NodeId r = kNoNode;
      if (k == 1) {
        r = parts[begin[g]];
      } else if (k >= 2) {
        std::vector<NodeId> children(parts.begin() + static_cast<std::ptrdiff_t>(begin[g]),
                                     parts.begin() + static_cast<std::ptrdiff_t>(begin[g + 1]));
        r = paper::add_internal_node(h, height, std::move(children), {});
        out_.node_upper.push_back(upper);
        out_.node_upper_inclusive.push_back(inclusive ? 1 : 0);
        out_.node_origin.push_back(origin);
      }
      rep_[group_root_[g]] = r;
    }
  }

  [[nodiscard]] const DisjointSets& dsu() const noexcept { return dsu_; }

 private:
  paper::LocalAssembly<W>& out_;
  const paper::ThresholdJoins<W>& joins_;
  DisjointSets dsu_;
  std::vector<NodeId> rep_;
  std::vector<std::uint64_t> stamp_of_;
  std::vector<std::uint32_t> group_of_;
  std::uint64_t stamp_ = 0;
  std::vector<DisjointSets::Item> pre_;
  std::vector<DisjointSets::Item> group_root_;
  std::vector<std::uint32_t> gid_;
};

/// Index of the largest point <= w in pts[lo, hi), or hi if none (binary search).
template <class W>
std::size_t floor_index(const std::vector<W>& pts, std::size_t lo, std::size_t hi, W w) {
  const auto it = std::upper_bound(pts.begin() + static_cast<std::ptrdiff_t>(lo),
                                   pts.begin() + static_cast<std::ptrdiff_t>(hi), w);
  const auto i = static_cast<std::size_t>(it - pts.begin());
  return i == lo ? hi : i - 1;
}

[[noreturn]] void broken(const std::string& what) {
  throw std::logic_error("paper_experimental::assemble_local_hierarchy: " + what);
}

/// The two constants of Sec 5.4 that D-8 replaces.
template <class W>
struct AlignmentRule {
  LocalAlignment kind;
  std::uint64_t n;
  std::uint64_t b;
  unsigned k_log2;

  /// Step 3: s_U aligned for sigma_U (l.1206: s N^c >= sigma; D-8: s Lambda >= sigma).
  [[nodiscard]] bool aligned(W s, W sigma) const {
    const W scaled = kind == LocalAlignment::kPaper
                         ? sat::mul_int_pow(s, n, paper::kLocalAlignExponent)
                         : lambda_scaled(s, n, b, k_log2);
    return !(scaled < sigma);
  }
  /// Step 4: upper end of the support of key sigma (l.1216: N^16 sigma; D-8: 2(3hK^2)^2 N^3 sigma).
  [[nodiscard]] W support_end(W sigma) const {
    return kind == LocalAlignment::kPaper ? sat::mul_int_pow(sigma, n, paper::kLocalSupportExponent)
                                          : lambda_support_end(sigma, n, b, k_log2);
  }
};

struct GridStats {
  std::size_t y_nontrivial = 0;
  std::size_t global_classes = 0;
};

/// Steps 1-8 of the grid assembly (assemble.hpp), with `rule` for the alignment and support.
template <class W>
GridStats grid_merges(paper::LocalAssembly<W>& out, const paper::OwnerDecomposition<W>& od,
                      const paper::RadialStateGraph<W>& rs, const paper::ThresholdJoins<W>& joins,
                      const AlignmentRule<W>& rule, MergeBuilder<W>& mb) {
  const std::size_t s = rs.num_states();
  const std::size_t nb = rs.num_blocks;
  const std::size_t nkeys = od.corridors.size();
  const auto& jr = joins.joins;
  const auto& theta = joins.theta;
  GridStats stats;
  out.blocks.assign(nb, BlockScale<W>{});
  out.y_component.assign(s, kNoComponent);

  // Members and block-MST joins of each block (counting sorts, Lemma 3.1).
  const auto [mbegin, members] = bucket(s, nb, [&](std::size_t x) { return rs.block[x]; });
  std::vector<std::uint32_t> mst_ids;
  for (std::size_t i = 0; i < jr.size(); ++i) {
    if (jr[i].kind == JoinKind::kBlockMst) mst_ids.push_back(static_cast<std::uint32_t>(i));
  }
  const auto [ebegin, eorder] =
      bucket(mst_ids.size(), nb, [&](std::size_t i) { return rs.block[jr[mst_ids[i]].a]; });
  const auto block_msts = [&](std::size_t u) {
    std::vector<std::uint32_t> ids;
    ids.reserve(ebegin[u + 1] - ebegin[u]);
    for (std::size_t i = ebegin[u]; i < ebegin[u + 1]; ++i) ids.push_back(mst_ids[eorder[i]]);
    return ids;
  };

  // Step 1-2: s_U = a K^i for the lowest native shell i present (l.1184-1186); z_U and sigma_U
  // (l.1201-1203).
  std::vector<std::vector<W>> powers(rs.component_anchor.size());
  const auto power = [&](std::uint32_t comp, std::uint32_t shell) {
    auto& p = powers[comp];
    if (p.empty()) p.push_back(rs.component_anchor[comp]);
    while (p.size() <= shell) p.push_back(sat::mul_pow2(p.back(), rs.shell_ratio_log2));
    return p[shell];
  };
  std::vector<W> mst_max(nb, W{0});
  for (std::size_t u = 0; u < nb; ++u) {
    auto& bs = out.blocks[u];
    bs.size = mbegin[u + 1] - mbegin[u];
    if (bs.size == 1) continue;  // l.1183: singleton block, attachment join only
    const Vertex v0 = rs.projection[members[mbegin[u]]];
    const std::uint32_t comp = rs.near_component[v0];
    if (comp == paper::kNoNearComponent) broken("a nonsingleton block contains a pivot copy");
    std::uint32_t lowest = rs.shell[v0];
    StateId z = kNoState;
    for (std::size_t i = mbegin[u]; i < mbegin[u + 1]; ++i) {
      const StateId x = members[i];
      const Vertex v = rs.projection[x];
      if (rs.near_component[v] != comp) broken("a block spans two near components");
      lowest = std::min(lowest, rs.shell[v]);
      if (theta[x] != kInf<W> && (z == kNoState || theta[x] < theta[z])) z = x;
    }
    bs.s = power(comp, lowest);
    for (const std::uint32_t j : block_msts(u)) mst_max[u] = std::max(mst_max[u], jr[j].weight);
    bs.kind = BlockScaleKind::kUnattached;
    if (z != kNoState) {
      bs.z = z;
      bs.sigma_key = joins.b0_key[z];
      if (bs.sigma_key >= nkeys) broken("finite theta without a named corridor key");
      bs.kind = BlockScaleKind::kUnaligned;  // decided in step 3
    }
    // Lemma 5.8 (r <= theta) with s_U <= r: every attachment of U is at least s_U.
    for (std::size_t i = mbegin[u]; i < mbegin[u + 1]; ++i) {
      if (theta[members[i]] < bs.s) broken("theta below s_U in a nonsingleton block");
    }
  }

  // Step 3: alignment per named key by successive medians (l.1214-1218); the test itself is
  // rule.aligned (c = 10 of l.1206, or D-8's Lambda).
  {
    std::vector<std::uint32_t> attached;
    for (std::size_t u = 0; u < nb; ++u) {
      if (out.blocks[u].kind == BlockScaleKind::kUnaligned) {
        attached.push_back(static_cast<std::uint32_t>(u));
      }
    }
    const auto [kbegin, korder] = bucket(
        attached.size(), nkeys, [&](std::size_t i) { return out.blocks[attached[i]].sigma_key; });
    for (std::size_t c = 0; c < nkeys; ++c) {
      const W sigma = od.corridors[c].weight;
      std::vector<std::uint32_t> rem;
      for (std::size_t i = kbegin[c]; i < kbegin[c + 1]; ++i) rem.push_back(attached[korder[i]]);
      while (!rem.empty()) {
        std::vector<W> keys;
        keys.reserve(rem.size());
        for (const std::uint32_t u : rem) keys.push_back(out.blocks[u].s);
        const W m = paper::select_bfprt(keys, keys.size() / 2);
        const bool aligned = rule.aligned(m, sigma);
        std::vector<std::uint32_t> next;
        for (const std::uint32_t u : rem) {
          auto& bs = out.blocks[u];
          if (aligned && !(bs.s < m)) {
            bs.kind = BlockScaleKind::kAligned;
          } else if (!aligned && !(m < bs.s)) {
            bs.kind = BlockScaleKind::kUnaligned;
          } else {
            next.push_back(u);
          }
        }
        rem = std::move(next);
      }
    }
  }

  // Step 4: supports per named corridor key (Q-37) and the global doubling grid; the support
  // end is rule.support_end (N^16 sigma of l.1216, or D-8's 2(3hK^2)^2 N^3 sigma).
  out.supports.assign(nkeys, GridSupport<W>{});
  for (std::size_t c = 0; c < nkeys; ++c) {
    out.supports[c].lower = od.corridors[c].weight;
    out.supports[c].upper = rule.support_end(od.corridors[c].weight);  // Q-38
  }
  for (std::size_t u = 0; u < nb; ++u) {
    const auto& bs = out.blocks[u];
    if (bs.kind != BlockScaleKind::kAligned) continue;
    auto& sp = out.supports[bs.sigma_key];
    sp.lower = std::min(sp.lower, bs.s);
  }
  for (StateId z = 0; z < s; ++z) {  // l.1221: theta(z) lies in the support of its own key
    if (theta[z] == kInf<W>) continue;
    const std::uint32_t c = joins.b0_key[z];
    if (c >= nkeys) broken("finite theta without a named corridor key");
    if (theta[z] < out.supports[c].lower || out.supports[c].upper < theta[z]) {
      broken("theta outside the support of its named key");
    }
  }
  W wmax{0};  // [Supplement] no grid point above the largest finite join weight is needed
  for (const auto& j : jr) {
    if (j.weight != kInf<W>) wmax = std::max(wmax, j.weight);
  }
  if (nkeys > 0) {
    std::vector<std::uint32_t> by_lower(nkeys);
    for (std::uint32_t c = 0; c < nkeys; ++c) by_lower[c] = c;
    std::sort(by_lower.begin(), by_lower.end(), [&](std::uint32_t x, std::uint32_t y) {
      if (out.supports[x].lower < out.supports[y].lower) return true;
      if (out.supports[y].lower < out.supports[x].lower) return false;
      return x < y;
    });
    std::vector<std::pair<W, W>> comps;  // overlapping supports merged into [lo, hi]
    for (const std::uint32_t c : by_lower) {
      const auto& sp = out.supports[c];
      if (!comps.empty() && !(comps.back().second < sp.lower)) {
        comps.back().second = std::max(comps.back().second, sp.upper);
      } else {
        comps.emplace_back(sp.lower, sp.upper);
      }
    }
    for (std::size_t i = 0; i < comps.size(); ++i) {
      const W top = std::min(comps[i].second, wmax);
      W g = comps[i].first;
      for (;;) {
        out.grid.push_back(g);
        const W next = sat::add(g, g);
        if (top < next) {
          // Last class of the component; stop it at the next component start (l.1279).
          out.grid_upper.push_back(i + 1 < comps.size() ? std::min(next, comps[i + 1].first)
                                                        : next);
          break;
        }
        out.grid_upper.push_back(next);
        g = next;
      }
    }
    for (auto& sp : out.supports) {
      const auto lo = std::upper_bound(out.grid.begin(), out.grid.end(), sp.lower);
      sp.lower_index = static_cast<std::size_t>(lo - out.grid.begin()) - 1;  // grid[0] <= lower
      const auto hi = std::upper_bound(out.grid.begin(), out.grid.end(), sp.upper);
      sp.upper_index = static_cast<std::size_t>(hi - out.grid.begin());
    }
  }
  const auto classify = [&](W w, std::uint32_t c) {
    const auto& sp = out.supports[c];
    const std::size_t i = floor_index(out.grid, sp.lower_index, sp.upper_index, w);
    if (i == sp.upper_index || !(w < out.grid_upper[i])) {
      broken("a value outside the grid of its named key");
    }
    return i;
  };

  // Step 5: t_U (l.1232-1236).
  for (auto& bs : out.blocks) {
    if (bs.kind == BlockScaleKind::kSingleton) continue;
    bs.t = bs.kind == BlockScaleKind::kAligned ? out.grid[classify(bs.s, bs.sigma_key)] : bs.s;
  }

  std::vector<std::uint8_t> global(jr.size(), 1);

  // Step 6: Y components (l.1189-1200), each with its exact Kruskal merge tree.
  {
    DisjointSets ydsu(s);
    std::vector<std::uint32_t> light;
    for (std::size_t u = 0; u < nb; ++u) {
      const auto& bs = out.blocks[u];
      if (bs.kind == BlockScaleKind::kSingleton) continue;
      for (const std::uint32_t j : block_msts(u)) {
        if (jr[j].weight < bs.t) {
          ydsu.unite(jr[j].a, jr[j].b);
          light.push_back(j);
          global[j] = 0;
        }
      }
    }
    std::vector<std::uint32_t> y_of_root(s, kNoComponent);
    std::uint32_t ny = 0;
    std::vector<std::uint32_t> y_size;
    for (StateId x = 0; x < s; ++x) {
      if (out.blocks[rs.block[x]].kind == BlockScaleKind::kSingleton) continue;
      const auto r = ydsu.find(x);
      if (y_of_root[r] == kNoComponent) {
        y_of_root[r] = ny++;
        y_size.push_back(0);
      }
      out.y_component[x] = y_of_root[r];
      ++y_size[y_of_root[r]];
    }
    for (const std::uint32_t sz : y_size) stats.y_nontrivial += sz >= 2 ? 1 : 0;
    const auto [ybegin, yorder] =
        bucket(light.size(), ny, [&](std::size_t i) { return out.y_component[jr[light[i]].a]; });
    std::vector<std::uint32_t> edges;
    for (std::uint32_t y = 0; y < ny; ++y) {
      edges.clear();
      for (std::size_t i = ybegin[y]; i < ybegin[y + 1]; ++i) edges.push_back(light[yorder[i]]);
      std::sort(edges.begin(), edges.end(), [&](std::uint32_t a, std::uint32_t b) {
        if (jr[a].weight < jr[b].weight) return true;
        if (jr[b].weight < jr[a].weight) return false;
        return a < b;
      });
      for (std::size_t i = 0; i < edges.size();) {
        std::size_t e = i + 1;
        while (e < edges.size() && !(jr[edges[i]].weight < jr[edges[e]].weight)) ++e;
        const W w = jr[edges[i]].weight;
        mb.merge_class(std::span<const std::uint32_t>(edges.data() + i, e - i), w, w, true,
                       NodeOrigin::kExactLocal);
        i = e;
      }
    }
  }

  // Step 7: separate grids (l.1209-1213, l.1236) for unaligned and unattached blocks.
  std::vector<std::uint32_t> separate_attach;
  std::vector<W> separate_top(nb, W{0});
  for (std::size_t u = 0; u < nb; ++u) {
    const auto& bs = out.blocks[u];
    if (bs.kind != BlockScaleKind::kUnaligned && bs.kind != BlockScaleKind::kUnattached) continue;
    std::vector<W> pts{bs.s};
    while (!(mst_max[u] < sat::add(pts.back(), pts.back()))) {
      pts.push_back(sat::add(pts.back(), pts.back()));
    }
    std::vector<std::uint32_t> rest;
    for (const std::uint32_t j : block_msts(u)) {
      if (!(jr[j].weight < bs.t)) rest.push_back(j);
    }
    const auto [cbegin, corder] = bucket(rest.size(), pts.size(), [&](std::size_t i) {
      return floor_index(pts, 0, pts.size(), jr[rest[i]].weight);
    });
    std::vector<std::uint32_t> cls;
    for (std::size_t c = 0; c < pts.size(); ++c) {
      if (cbegin[c] == cbegin[c + 1]) continue;
      cls.clear();
      for (std::size_t i = cbegin[c]; i < cbegin[c + 1]; ++i) {
        cls.push_back(rest[corder[i]]);
        global[rest[corder[i]]] = 0;
      }
      const W upper = sat::add(pts[c], pts[c]);
      separate_top[u] = upper;
      mb.merge_class(cls, pts[c], upper, false, NodeOrigin::kSeparateGrid);
    }
    if (bs.kind == BlockScaleKind::kUnaligned)
      separate_attach.push_back(static_cast<std::uint32_t>(u));
  }

  // Step 8: global classes (l.1237-1238, l.1244-1250).
  std::vector<std::uint32_t> gl;
  std::vector<std::size_t> gl_class;
  for (std::size_t i = 0; i < jr.size(); ++i) {
    if (global[i] == 0) continue;
    std::uint32_t key = kNoCorridor;
    switch (jr[i].kind) {
      case JoinKind::kBlockMst:
        key = out.blocks[rs.block[jr[i].a]].sigma_key;
        if (out.blocks[rs.block[jr[i].a]].kind != BlockScaleKind::kAligned) {
          broken("a block-MST edge of a non-aligned block left for the global grid");
        }
        break;
      case JoinKind::kOwner:
        key = joins.b0_key[jr[i].a];
        break;
      case JoinKind::kCorridor:
        key = jr[i].ref;
        break;
    }
    gl.push_back(static_cast<std::uint32_t>(i));
    gl_class.push_back(classify(jr[i].weight, key));
  }
  for (const std::uint32_t u : separate_attach) {
    // Lemma 5.10 l.1276-1278: an unaligned block's local tree precedes its first attachment.
    const auto& bs = out.blocks[u];
    const W h0 = out.grid[classify(theta[bs.z], bs.sigma_key)];
    if (h0 < separate_top[u]) broken("a separate-grid class reaches the first attachment class");
  }
  {
    const auto [gbegin, gorder] =
        bucket(gl.size(), out.grid.size(), [&](std::size_t i) { return gl_class[i]; });
    std::vector<std::uint32_t> cls;
    for (std::size_t c = 0; c < out.grid.size(); ++c) {
      if (gbegin[c] == gbegin[c + 1]) continue;
      ++stats.global_classes;
      cls.clear();
      for (std::size_t i = gbegin[c]; i < gbegin[c + 1]; ++i) cls.push_back(gl[gorder[i]]);
      mb.merge_class(cls, out.grid[c], out.grid_upper[c], false, NodeOrigin::kGlobalGrid);
    }
  }
  return stats;
}

/// D-11 (Q-39): all join records by (weight, index), each weight merged simultaneously at
/// height = upper = that weight (inclusive). Returns the number of classes.
template <class W>
std::size_t exact_merges(const paper::ThresholdJoins<W>& joins, MergeBuilder<W>& mb) {
  const auto& jr = joins.joins;
  std::vector<std::uint32_t> order(jr.size());
  for (std::size_t i = 0; i < jr.size(); ++i) order[i] = static_cast<std::uint32_t>(i);
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (jr[a].weight < jr[b].weight) return true;
    if (jr[b].weight < jr[a].weight) return false;
    return a < b;
  });
  std::size_t classes = 0;
  for (std::size_t i = 0; i < order.size();) {
    std::size_t e = i + 1;
    while (e < order.size() && !(jr[order[i]].weight < jr[order[e]].weight)) ++e;
    const W w = jr[order[i]].weight;
    if (w == kInf<W>) broken("an infinite join record");  // Sec 5.3: only finite joins
    mb.merge_class(std::span<const std::uint32_t>(order.data() + i, e - i), w, w, true,
                   NodeOrigin::kExactLocal);
    ++classes;
    i = e;
  }
  return classes;
}

/// Steps 9-10 (l.1281-1293): links from the arcs at the lowest common ancestor and orders.
/// Returns the number of links.
template <class W>
std::size_t links_and_order(paper::LocalAssembly<W>& out, std::size_t s) {
  auto& h = out.result.hierarchy;
  const std::size_t nn = h.nodes.size();
  std::vector<Vertex> parent(nn);
  for (std::size_t v = 0; v < nn; ++v) {
    parent[v] = h.nodes[v].parent == kNoNode ? kNoVertex : h.nodes[v].parent;
  }
  const ForestQueries fq(parent);
  std::vector<std::uint32_t> rank_in_parent(nn, 0);
  for (std::size_t v = 0; v < nn; ++v) {
    const auto& ch = h.nodes[v].children;
    for (std::size_t i = 0; i < ch.size(); ++i)
      rank_in_parent[ch[i]] = static_cast<std::uint32_t>(i);
  }
  struct Assigned {
    NodeId node;
    std::uint32_t cx;
    std::uint32_t cy;
    ArcId arc;
  };
  std::vector<Assigned> assigned;
  std::vector<std::pair<Vertex, Vertex>> root_arcs;
  const auto& arcs = h.graph.arcs;
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    const StateId x = arcs[a].tail;
    const StateId y = arcs[a].head;
    if (x == y) continue;
    const std::optional<Vertex> l = fq.lca(x, y);
    if (!l) {
      root_arcs.emplace_back(fq.component(x), fq.component(y));
      continue;
    }
    const std::uint32_t d = fq.depth(*l) + 1;
    assigned.push_back({*l, rank_in_parent[fq.level_ancestor(x, d)],
                        rank_in_parent[fq.level_ancestor(y, d)], static_cast<ArcId>(a)});
  }
  const auto [abegin, aorder] =
      bucket(assigned.size(), nn, [&](std::size_t i) { return assigned[i].node; });
  std::size_t links = 0;
  for (std::size_t v = s; v < nn; ++v) {
    auto& node = h.nodes[v];
    const std::size_t k = node.children.size();
    std::vector<std::vector<std::pair<std::uint32_t, ArcId>>> adj(k);
    std::vector<DiArc> light;
    for (std::size_t i = abegin[v]; i < abegin[v + 1]; ++i) {
      const auto& as = assigned[aorder[i]];
      const W len = arcs[as.arc].length;
      const bool connecting =
          out.node_upper_inclusive[v] != 0 ? !(out.node_upper[v] < len) : len < out.node_upper[v];
      if (connecting) {
        adj[as.cx].emplace_back(as.cy, as.arc);
        adj[as.cy].emplace_back(as.cx, as.arc);
      }
      if (len < node.height) light.push_back({as.cx, as.cy});
    }
    // DFS over the child adjacency; one physical connector per discovery arc.
    std::vector<std::uint8_t> seen(k, 0);
    std::vector<std::uint32_t> stack{0};
    seen[0] = 1;
    while (!stack.empty()) {
      const std::uint32_t c = stack.back();
      stack.pop_back();
      for (const auto& [d, a] : adj[c]) {
        if (seen[d] != 0) continue;
        seen[d] = 1;
        stack.push_back(d);
        node.links.push_back({arcs[a].tail, arcs[a].head, arcs[a].length, a});
      }
    }
    if (node.links.size() + 1 != k) broken("children of a merge are not connected by its arcs");
    links += node.links.size();
    // Children in a topological order of the light arcs (shorter than the height).
    std::vector<Vertex> order;
    try {
      order = paper::topological_sort(Digraph::from_arcs(k, light));
    } catch (const std::invalid_argument&) {
      broken("light arcs between the children of a merge form a cycle");
    }
    std::vector<NodeId> sorted(k);
    for (std::size_t i = 0; i < k; ++i) sorted[i] = node.children[order[i]];
    node.children = std::move(sorted);
  }
  paper::collect_roots(h);
  {
    std::vector<std::uint32_t> root_rank(nn, 0);
    for (std::size_t i = 0; i < h.roots.size(); ++i) {
      root_rank[h.roots[i]] = static_cast<std::uint32_t>(i);
    }
    std::vector<DiArc> q;
    q.reserve(root_arcs.size());
    for (const auto& [rx, ry] : root_arcs) q.push_back({root_rank[rx], root_rank[ry]});
    std::vector<Vertex> order;
    try {
      order = paper::topological_sort(Digraph::from_arcs(h.roots.size(), q));
    } catch (const std::invalid_argument&) {
      broken("the roots are not the SCCs of the state graph");
    }
    std::vector<NodeId> sorted(h.roots.size());
    for (std::size_t i = 0; i < sorted.size(); ++i) sorted[i] = h.roots[order[i]];
    h.roots = std::move(sorted);
  }
  return links;
}

}  // namespace

template <class W>
W lambda_scaled(W s, std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  // Lambda = 2 * N * 3 * h * K^2, h = b + 1: integer multiplications by doubling and adding,
  // K^2 = 2^(2 k_log2) by doublings (Sec 2 l.297); each step saturates (Q-38).
  if (b == std::numeric_limits<std::uint64_t>::max()) return kInf<W>;  // h does not fit
  W x = sat::mul_int(s, 2);
  x = sat::mul_int(x, n);
  x = sat::mul_int(x, 3);
  x = sat::mul_int(x, b + 1);
  return sat::mul_pow2(x, 2 * k_log2);
}

template <class W>
W lambda_support_end(W sigma, std::uint64_t n, std::uint64_t b, unsigned k_log2) {
  // 2 * (3 h K^2)^2 * N^3 * sigma.
  if (b == std::numeric_limits<std::uint64_t>::max()) return kInf<W>;
  W x = sat::mul_int(sigma, 2);
  for (int i = 0; i < 2; ++i) {
    x = sat::mul_int(x, 3);
    x = sat::mul_int(x, b + 1);
    x = sat::mul_pow2(x, 2 * k_log2);
  }
  return sat::mul_int_pow(x, n, 3);
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
  const std::size_t n = od.num_vertices();
  const std::size_t s = rs.num_states();
  if (rs.num_vertices() != n || joins.num_states != s || joins.num_pivots != od.pivots.size() ||
      joins.theta.size() != s || joins.b0_key.size() != s || rs.block.size() != s) {
    throw std::invalid_argument(
        "paper_experimental::assemble_local_hierarchy: inputs of different constructions");
  }
  if (!local_assembly_variant_applies(n, od.separator, alignment, assembly)) {
    throw std::invalid_argument(
        "paper_experimental::assemble_local_hierarchy: N = " + std::to_string(n) +
        ", b = " + std::to_string(od.separator) +
        (alignment == LocalAlignment::kPaper && assembly == LocalAssembly::kGrid
             ? " violates N >= 32 and h = b + 1 <= N (Sec 5.4 l.1206)"
             : " violates N >= 2 (D-8 / D-11)"));
  }

  paper::LocalAssembly<W> out;
  out.result.hierarchy = paper::make_leaf_forest(rs.graph, rs.projection);
  auto& h = out.result.hierarchy;
  out.result.native_entry = rs.native_inward;
  out.result.arc_edge = rs.arc_edge;
  out.node_upper.assign(s, W{0});
  out.node_upper_inclusive.assign(s, 0);
  out.node_origin.assign(s, NodeOrigin::kExactLocal);

  MergeBuilder<W> mb(out, joins);
  GridStats stats;
  std::size_t exact_classes = 0;
  if (assembly == LocalAssembly::kGrid) {
    h.c0 = 2;  // Lemma 5.10
    const AlignmentRule<W> rule{alignment, n, od.separator, rs.shell_ratio_log2};
    stats = grid_merges(out, od, rs, joins, rule, mb);
  } else {
    h.c0 = 1;  // D-11: links are arcs of length <= the exact height
    exact_classes = exact_merges(joins, mb);
  }
  const std::size_t links = links_and_order(out, s);

  if (diag != nullptr) {
    const std::size_t nn = h.nodes.size();
    SSSP_COUNT(++diag->counters["local_assemble_calls"]);
    SSSP_COUNT(diag->counters["local_assemble_blocks"] += rs.num_blocks);
    if (assembly == LocalAssembly::kGrid) {
      std::size_t kinds[4] = {0, 0, 0, 0};
      for (const auto& bs : out.blocks) ++kinds[static_cast<std::size_t>(bs.kind)];
      SSSP_COUNT(diag->counters["local_assemble_aligned_blocks"] += kinds[1]);
      SSSP_COUNT(diag->counters["local_assemble_unaligned_blocks"] += kinds[2]);
      SSSP_COUNT(diag->counters["local_assemble_unattached_blocks"] += kinds[3]);
      SSSP_COUNT(diag->counters["local_assemble_y_components"] += stats.y_nontrivial);
      SSSP_COUNT(diag->counters["local_assemble_grid_points"] += out.grid.size());
      SSSP_COUNT(diag->counters["local_assemble_global_classes"] += stats.global_classes);
      static_cast<void>(kinds);
    } else {
      SSSP_COUNT(diag->counters["local_assemble_exact_classes"] += exact_classes);
    }
    SSSP_COUNT(diag->counters["local_assemble_nodes"] += nn - s);
    SSSP_COUNT(diag->counters["local_assemble_links"] += links);
    SSSP_COUNT(mb.dsu().add_counters(*diag));
    static_cast<void>(nn);
  }
  static_cast<void>(stats);
  static_cast<void>(exact_classes);
  static_cast<void>(links);
  return out;
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
