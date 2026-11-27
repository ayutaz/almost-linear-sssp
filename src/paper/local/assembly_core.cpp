// Shared core of the Sec 5.4 assembly and of the post-processing of Lemma 5.10. See
// assembly_core.hpp; the step numbers are those of assemble.hpp.

#include "paper/local/assembly_core.hpp"

#include <stdexcept>

#include "paper/primitives/forest_queries.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper::assembly {

namespace {

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

[[noreturn]] void broken(std::string_view context, const std::string& what) {
  throw std::logic_error(std::string(context) + ": " + what);
}

/// Exact merges (step 6, l.1189-1200): the records `ids` sorted by (weight, index), the records
/// of one weight merged simultaneously at height = upper = that weight (inclusive).
template <class W, class Merger>
void exact_merges(std::vector<std::uint32_t>& ids, const std::vector<JoinRecord<W>>& jr,
                  Merger& merger) {
  std::sort(ids.begin(), ids.end(), [&](std::uint32_t a, std::uint32_t b) {
    if (jr[a].weight < jr[b].weight) return true;
    if (jr[b].weight < jr[a].weight) return false;
    return a < b;
  });
  for (std::size_t i = 0; i < ids.size();) {
    std::size_t e = i + 1;
    while (e < ids.size() && !(jr[ids[i]].weight < jr[ids[e]].weight)) ++e;
    const W w = jr[ids[i]].weight;
    merger.merge_class(std::span<const std::uint32_t>(ids.data() + i, e - i), w, w, true,
                       NodeOrigin::kExactLocal);
    i = e;
  }
}

struct GridStats {
  std::size_t y_nontrivial = 0;
  std::size_t global_classes = 0;
};

/// Steps 1-8 of the grid assembly with the two scale tests of `rule`.
template <class W, class Merger>
GridStats grid_merges(LocalAssembly<W>& out, const OwnerDecomposition<W>& od,
                      const RadialStateGraph<W>& rs, const ThresholdJoins<W>& joins,
                      const Rule& rule, Merger& mb) {
  const std::string_view context = rule.context;
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

  // Step 1-2: s_U = a K^i for the lowest native shell i present (l.1184-1186), with a K^j cached
  // per near component and formed by doublings; z_U and sigma_U (l.1201-1203).
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
    if (comp == kNoNearComponent) broken(context, "a nonsingleton block contains a pivot copy");
    std::uint32_t lowest = rs.shell[v0];
    StateId z = kNoState;
    for (std::size_t i = mbegin[u]; i < mbegin[u + 1]; ++i) {
      const StateId x = members[i];
      const Vertex v = rs.projection[x];
      if (rs.near_component[v] != comp) broken(context, "a block spans two near components");
      lowest = std::min(lowest, rs.shell[v]);
      if (theta[x] != kInf<W> && (z == kNoState || theta[x] < theta[z])) z = x;
    }
    bs.s = power(comp, lowest);
    for (const std::uint32_t j : block_msts(u)) mst_max[u] = std::max(mst_max[u], jr[j].weight);
    bs.kind = BlockScaleKind::kUnattached;
    if (z != kNoState) {
      bs.z = z;
      bs.sigma_key = joins.b0_key[z];
      if (bs.sigma_key >= nkeys) broken(context, "finite theta without a named corridor key");
      bs.kind = BlockScaleKind::kUnaligned;  // decided in step 3
    }
    // Lemma 5.8 (r <= theta) with s_U <= r: every attachment of U is at least s_U. The Y
    // components and the separate grids below rely on it.
    for (std::size_t i = mbegin[u]; i < mbegin[u + 1]; ++i) {
      if (theta[members[i]] < bs.s) broken(context, "theta below s_U in a nonsingleton block");
    }
  }

  // Step 3: alignment per named key by successive medians (l.1214-1218); the scale of the test
  // is rule.align (N^c with c = 10 in the paper, l.1206).
  {
    std::vector<std::uint32_t> attached;
    for (std::size_t u = 0; u < nb; ++u) {
      if (out.blocks[u].kind == BlockScaleKind::kUnaligned) {
        attached.push_back(static_cast<std::uint32_t>(u));
      }
    }
    const auto [kbegin, korder] = bucket(
        attached.size(), nkeys, [&](std::size_t i) { return out.blocks[attached[i]].sigma_key; });
    std::vector<W> scales;
    for (std::size_t c = 0; c < nkeys; ++c) {
      scales.clear();
      for (std::size_t i = kbegin[c]; i < kbegin[c + 1]; ++i) {
        scales.push_back(out.blocks[attached[korder[i]]].s);
      }
      const auto aligned = aligned_by_medians<W>(scales, od.corridors[c].weight, rule.align);
      for (std::size_t i = kbegin[c]; i < kbegin[c + 1]; ++i) {
        if (aligned[i - kbegin[c]] != 0) {
          out.blocks[attached[korder[i]]].kind = BlockScaleKind::kAligned;
        }
      }
    }
  }

  // Step 4: supports per named corridor key (Q-37) and the global doubling grid; the support
  // end is rule.support_end (N^16 sigma in the paper, l.1216; saturating, Q-38).
  SupportGrid<W> sg;
  {
    std::vector<W> units(nkeys);
    for (std::size_t c = 0; c < nkeys; ++c) units[c] = od.corridors[c].weight;
    auto supports = key_supports<W>(units, rule.support_end);
    for (std::size_t u = 0; u < nb; ++u) {
      const auto& bs = out.blocks[u];
      if (bs.kind != BlockScaleKind::kAligned) continue;
      auto& sp = supports[bs.sigma_key];
      sp.lower = std::min(sp.lower, bs.s);
    }
    for (StateId z = 0; z < s; ++z) {  // l.1221: theta(z) lies in the support of its own key
      if (theta[z] == kInf<W>) continue;
      const std::uint32_t c = joins.b0_key[z];
      if (c >= nkeys) broken(context, "finite theta without a named corridor key");
      if (theta[z] < supports[c].lower || supports[c].upper < theta[z]) {
        broken(context, "theta outside the support of its named key");
      }
    }
    W wmax{0};  // [Supplement] no grid point above the largest finite join weight is needed
    for (const auto& j : jr) {
      if (j.weight != kInf<W>) wmax = std::max(wmax, j.weight);
    }
    sg = build_support_grid(std::move(supports), wmax);
  }
  // Classification of a value labelled with key c: binary search between the key's indices.
  const auto classify = [&](W w, std::uint32_t c) {
    const auto i = sg.classify(w, c);
    if (!i) broken(context, "a value outside the grid of its named key");
    return *i;
  };

  // Step 5: t_U (l.1232-1236).
  for (auto& bs : out.blocks) {
    if (bs.kind == BlockScaleKind::kSingleton) continue;
    bs.t = bs.kind == BlockScaleKind::kAligned ? sg.grid[classify(bs.s, bs.sigma_key)] : bs.s;
  }

  std::vector<std::uint8_t> global(jr.size(), 1);

  // Step 6: Y components (l.1189-1200): components of the block-MST edges below t_U, each with
  // its exact Kruskal merge tree (sorted edges, equal weights simultaneously).
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
      exact_merges(edges, jr, mb);
    }
  }

  // Step 7: separate grids (l.1209-1213, l.1236) for unaligned and unattached blocks.
  std::vector<std::uint32_t> separate_attach;  // blocks whose top class must precede theta(z_U)
  std::vector<W> separate_top(nb, W{0});
  for (std::size_t u = 0; u < nb; ++u) {
    const auto& bs = out.blocks[u];
    if (bs.kind != BlockScaleKind::kUnaligned && bs.kind != BlockScaleKind::kUnattached) continue;
    const std::vector<W> pts = doubling_points(bs.s, mst_max[u]);
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
      mb.merge_class(std::span<const std::uint32_t>(cls), pts[c], upper, false,
                     NodeOrigin::kSeparateGrid);
    }
    if (bs.kind == BlockScaleKind::kUnaligned)
      separate_attach.push_back(static_cast<std::uint32_t>(u));
  }

  // Step 8: global classes (l.1237-1238, l.1244-1250). Owner joins and corridor joins carry
  // their own key; the remaining block-MST edges of aligned blocks carry sigma_U (l.1235).
  std::vector<std::uint32_t> gl;
  std::vector<std::size_t> gl_class;
  for (std::size_t i = 0; i < jr.size(); ++i) {
    if (global[i] == 0) continue;
    std::uint32_t key = kNoCorridor;
    switch (jr[i].kind) {
      case JoinKind::kBlockMst:
        key = out.blocks[rs.block[jr[i].a]].sigma_key;
        if (out.blocks[rs.block[jr[i].a]].kind != BlockScaleKind::kAligned) {
          broken(context, "a block-MST edge of a non-aligned block left for the global grid");
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
    const W h0 = sg.grid[classify(theta[bs.z], bs.sigma_key)];
    if (h0 < separate_top[u]) {
      broken(context, "a separate-grid class reaches the first attachment class");
    }
  }
  {
    const auto [gbegin, gorder] =
        bucket(gl.size(), sg.grid.size(), [&](std::size_t i) { return gl_class[i]; });
    std::vector<std::uint32_t> cls;
    for (std::size_t c = 0; c < sg.grid.size(); ++c) {
      if (gbegin[c] == gbegin[c + 1]) continue;
      ++stats.global_classes;
      cls.clear();
      for (std::size_t i = gbegin[c]; i < gbegin[c + 1]; ++i) cls.push_back(gl[gorder[i]]);
      mb.merge_class(std::span<const std::uint32_t>(cls), sg.grid[c], sg.grid_upper[c], false,
                     NodeOrigin::kGlobalGrid);
    }
  }
  out.grid = std::move(sg.grid);
  out.grid_upper = std::move(sg.grid_upper);
  out.supports = std::move(sg.supports);
  return stats;
}

}  // namespace

template <class W>
std::vector<std::uint8_t> aligned_by_medians(std::span<const W> s, W sigma,
                                             const LengthScale& align) {
  std::vector<std::uint8_t> out(s.size(), 0);
  std::vector<std::uint32_t> rem(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) rem[i] = static_cast<std::uint32_t>(i);
  while (!rem.empty()) {
    std::vector<W> keys;
    keys.reserve(rem.size());
    for (const std::uint32_t i : rem) keys.push_back(s[i]);
    const W m = select_bfprt(keys, keys.size() / 2);
    const bool aligned = !(align.apply(m) < sigma);
    std::vector<std::uint32_t> next;
    for (const std::uint32_t i : rem) {
      if (aligned && !(s[i] < m)) {
        out[i] = 1;
      } else if (!aligned && !(m < s[i])) {
        out[i] = 0;
      } else {
        next.push_back(i);
      }
    }
    rem = std::move(next);
  }
  return out;
}

template <class W>
std::vector<GridSupport<W>> key_supports(std::span<const W> units, const LengthScale& end) {
  std::vector<GridSupport<W>> supports(units.size());
  for (std::size_t c = 0; c < units.size(); ++c) {
    supports[c].lower = units[c];
    supports[c].upper = end.apply(units[c]);
  }
  return supports;
}

template <class W>
SupportGrid<W> build_support_grid(std::vector<GridSupport<W>> supports, W cap) {
  SupportGrid<W> out;
  out.supports = std::move(supports);
  const std::size_t nkeys = out.supports.size();
  if (nkeys == 0) return out;
  std::vector<std::uint32_t> by_lower(nkeys);
  for (std::uint32_t c = 0; c < nkeys; ++c) by_lower[c] = c;
  std::sort(by_lower.begin(), by_lower.end(), [&](std::uint32_t x, std::uint32_t y) {
    if (out.supports[x].lower < out.supports[y].lower) return true;
    if (out.supports[y].lower < out.supports[x].lower) return false;
    return x < y;
  });
  // Merge overlapping supports into components [lo, hi].
  std::vector<std::pair<W, W>> comps;
  for (const std::uint32_t c : by_lower) {
    const auto& sp = out.supports[c];
    if (!comps.empty() && !(comps.back().second < sp.lower)) {
      comps.back().second = std::max(comps.back().second, sp.upper);
    } else {
      comps.emplace_back(sp.lower, sp.upper);
    }
  }
  for (std::size_t i = 0; i < comps.size(); ++i) {
    const W top = std::min(comps[i].second, cap);
    W g = comps[i].first;
    for (;;) {
      out.grid.push_back(g);
      const W next = sat::add(g, g);
      if (top < next || next == kInf<W>) {
        // Last class of the component; stop it at the next component start (l.1279).
        out.grid_upper.push_back(i + 1 < comps.size() ? std::min(next, comps[i + 1].first) : next);
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
  return out;
}

template <class W>
std::vector<W> doubling_points(W unit, W top) {
  std::vector<W> pts{unit};
  for (;;) {
    const W next = sat::add(pts.back(), pts.back());
    if (top < next || next == kInf<W>) break;
    pts.push_back(next);
  }
  return pts;
}

template <class W>
LinkOrder link_and_order(std::vector<HierarchyNode<W>>& nodes, std::span<const StateArc<W>> arcs,
                         std::span<const NodeId> leaf_of_state, std::span<const W> node_upper,
                         std::span<const std::uint8_t> node_upper_inclusive,
                         std::string_view context) {
  const std::size_t nn = nodes.size();
  std::vector<Vertex> parent(nn);
  for (std::size_t v = 0; v < nn; ++v) {
    parent[v] = nodes[v].parent == kNoNode ? kNoVertex : nodes[v].parent;
  }
  const ForestQueries fq(parent);
  std::vector<std::uint32_t> rank_in_parent(nn, 0);
  for (std::size_t v = 0; v < nn; ++v) {
    const auto& ch = nodes[v].children;
    for (std::size_t i = 0; i < ch.size(); ++i)
      rank_in_parent[ch[i]] = static_cast<std::uint32_t>(i);
  }
  // Step 9 (l.1282-1284): every arc goes to the lowest common ancestor of its endpoints' leaves.
  struct Assigned {
    NodeId node;
    std::uint32_t cx;
    std::uint32_t cy;
    ArcId arc;
  };
  std::vector<Assigned> assigned;
  std::vector<std::pair<Vertex, Vertex>> root_arcs;
  for (std::size_t a = 0; a < arcs.size(); ++a) {
    const NodeId x = leaf_of_state[arcs[a].tail];
    const NodeId y = leaf_of_state[arcs[a].head];
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
  LinkOrder out;
  for (std::size_t v = 0; v < nn; ++v) {
    auto& node = nodes[v];
    const std::size_t k = node.children.size();
    if (k == 0) continue;  // leaf
    std::vector<std::vector<std::pair<std::uint32_t, ArcId>>> adj(k);
    std::vector<DiArc> light;
    for (std::size_t i = abegin[v]; i < abegin[v + 1]; ++i) {
      const auto& as = assigned[aorder[i]];
      const W len = arcs[as.arc].length;
      const bool connecting =
          node_upper_inclusive[v] != 0 ? !(node_upper[v] < len) : len < node_upper[v];
      if (connecting) {
        adj[as.cx].emplace_back(as.cy, as.arc);
        adj[as.cy].emplace_back(as.cx, as.arc);
      }
      if (len < node.height) light.push_back({as.cx, as.cy});
    }
    // DFS over the child adjacency; one physical connector per discovery arc (l.1285-1292).
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
    if (node.links.size() + 1 != k) {
      broken(context, "children of a merge are not connected by its arcs");
    }
    out.links += node.links.size();
    // Step 10 (l.1292-1293): children in a topological order of the light arcs (shorter than
    // the height).
    std::vector<Vertex> order;
    try {
      order = topological_sort(Digraph::from_arcs(k, light));
    } catch (const std::invalid_argument&) {
      broken(context, "light arcs between the children of a merge form a cycle");
    }
    std::vector<NodeId> sorted(k);
    for (std::size_t i = 0; i < k; ++i) sorted[i] = node.children[order[i]];
    node.children = std::move(sorted);
  }
  // Roots (the SCCs of the whole state graph) in condensation order (l.1280).
  std::vector<NodeId> roots;
  for (std::size_t v = 0; v < nn; ++v) {
    if (nodes[v].parent == kNoNode) roots.push_back(static_cast<NodeId>(v));
  }
  std::vector<std::uint32_t> root_rank(nn, 0);
  for (std::size_t i = 0; i < roots.size(); ++i)
    root_rank[roots[i]] = static_cast<std::uint32_t>(i);
  std::vector<DiArc> q;
  q.reserve(root_arcs.size());
  for (const auto& [rx, ry] : root_arcs) q.push_back({root_rank[rx], root_rank[ry]});
  std::vector<Vertex> order;
  try {
    order = topological_sort(Digraph::from_arcs(roots.size(), q));
  } catch (const std::invalid_argument&) {
    broken(context, "the roots are not the SCCs of the state graph");
  }
  out.roots.resize(roots.size());
  for (std::size_t i = 0; i < roots.size(); ++i) out.roots[i] = roots[order[i]];
  return out;
}

template <class W>
void check_inputs(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                  const ThresholdJoins<W>& joins, std::string_view context) {
  const std::size_t s = rs.num_states();
  if (rs.num_vertices() != od.num_vertices() || joins.num_states != s ||
      joins.num_pivots != od.pivots.size() || joins.theta.size() != s || joins.b0_key.size() != s ||
      rs.block.size() != s) {
    throw std::invalid_argument(std::string(context) + ": inputs of different constructions");
  }
}

template <class W>
LocalAssembly<W> assemble(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                          const ThresholdJoins<W>& joins, const Rule& rule, Diagnostics* diag) {
  check_inputs(od, rs, joins, rule.context);
  const std::size_t s = rs.num_states();
  const auto& jr = joins.joins;

  LocalAssembly<W> out;
  out.result.hierarchy = make_leaf_forest(rs.graph, rs.projection);
  auto& h = out.result.hierarchy;
  out.result.native_entry = rs.native_inward;
  out.result.arc_edge = rs.arc_edge;
  out.node_upper.assign(s, W{0});
  out.node_upper_inclusive.assign(s, 0);
  out.node_origin.assign(s, NodeOrigin::kExactLocal);

  // Merges use the disjoint sets of Lemma 3.6 on the join vertices: the states (node x is the
  // leaf of state x), then one auxiliary vertex per pivot (no node).
  std::vector<NodeId> node_of_vertex(joins.num_join_vertices(), kNoNode);
  for (std::size_t x = 0; x < joins.num_states; ++x) node_of_vertex[x] = static_cast<NodeId>(x);
  ClassMerger mb(
      std::move(node_of_vertex),
      [&jr](std::uint32_t j) { return std::pair<std::uint32_t, std::uint32_t>{jr[j].a, jr[j].b}; },
      [&out](W height, std::vector<NodeId> children, W upper, bool inclusive, NodeOrigin origin) {
        const NodeId r = add_internal_node(out.result.hierarchy, height, std::move(children), {});
        out.node_upper.push_back(upper);
        out.node_upper_inclusive.push_back(inclusive ? 1 : 0);
        out.node_origin.push_back(origin);
        return r;
      });

  h.c0 = 2;  // Lemma 5.10
  const GridStats stats = grid_merges(out, od, rs, joins, rule, mb);

  // Steps 9-10: links and orders (l.1281-1293).
  auto lo = link_and_order<W>(h.nodes, h.graph.arcs, h.leaf_of_state, out.node_upper,
                              out.node_upper_inclusive, rule.context);
  h.roots = std::move(lo.roots);

  if (diag != nullptr) {
    const std::size_t nn = h.nodes.size();
    SSSP_COUNT(++diag->counters["local_assemble_calls"]);
    SSSP_COUNT(diag->counters["local_assemble_blocks"] += rs.num_blocks);
    std::size_t kinds[4] = {0, 0, 0, 0};
    for (const auto& bs : out.blocks) ++kinds[static_cast<std::size_t>(bs.kind)];
    SSSP_COUNT(diag->counters["local_assemble_aligned_blocks"] += kinds[1]);
    SSSP_COUNT(diag->counters["local_assemble_unaligned_blocks"] += kinds[2]);
    SSSP_COUNT(diag->counters["local_assemble_unattached_blocks"] += kinds[3]);
    SSSP_COUNT(diag->counters["local_assemble_y_components"] += stats.y_nontrivial);
    SSSP_COUNT(diag->counters["local_assemble_grid_points"] += out.grid.size());
    SSSP_COUNT(diag->counters["local_assemble_global_classes"] += stats.global_classes);
    static_cast<void>(kinds);
    SSSP_COUNT(diag->counters["local_assemble_nodes"] += nn - s);
    SSSP_COUNT(diag->counters["local_assemble_links"] += lo.links);
    SSSP_COUNT(mb.dsu().add_counters(*diag));
    static_cast<void>(nn);
  }
  static_cast<void>(stats);
  return out;
}

// Explicit instantiations for std::uint64_t.
template std::vector<std::uint8_t> aligned_by_medians<std::uint64_t>(std::span<const std::uint64_t>,
                                                                     std::uint64_t,
                                                                     const LengthScale&);
template std::vector<GridSupport<std::uint64_t>> key_supports<std::uint64_t>(
    std::span<const std::uint64_t>, const LengthScale&);
template SupportGrid<std::uint64_t> build_support_grid<std::uint64_t>(
    std::vector<GridSupport<std::uint64_t>>, std::uint64_t);
template std::vector<std::uint64_t> doubling_points<std::uint64_t>(std::uint64_t, std::uint64_t);
template LinkOrder link_and_order<std::uint64_t>(std::vector<HierarchyNode<std::uint64_t>>&,
                                                 std::span<const StateArc<std::uint64_t>>,
                                                 std::span<const NodeId>,
                                                 std::span<const std::uint64_t>,
                                                 std::span<const std::uint8_t>, std::string_view);
template void check_inputs<std::uint64_t>(const OwnerDecomposition<std::uint64_t>&,
                                          const RadialStateGraph<std::uint64_t>&,
                                          const ThresholdJoins<std::uint64_t>&, std::string_view);
template LocalAssembly<std::uint64_t> assemble<std::uint64_t>(
    const OwnerDecomposition<std::uint64_t>&, const RadialStateGraph<std::uint64_t>&,
    const ThresholdJoins<std::uint64_t>&, const Rule&, Diagnostics*);

// Explicit instantiations for double.
template std::vector<std::uint8_t> aligned_by_medians<double>(std::span<const double>, double,
                                                              const LengthScale&);
template std::vector<GridSupport<double>> key_supports<double>(std::span<const double>,
                                                               const LengthScale&);
template SupportGrid<double> build_support_grid<double>(std::vector<GridSupport<double>>, double);
template std::vector<double> doubling_points<double>(double, double);
template LinkOrder link_and_order<double>(std::vector<HierarchyNode<double>>&,
                                          std::span<const StateArc<double>>,
                                          std::span<const NodeId>, std::span<const double>,
                                          std::span<const std::uint8_t>, std::string_view);
template void check_inputs<double>(const OwnerDecomposition<double>&,
                                   const RadialStateGraph<double>&, const ThresholdJoins<double>&,
                                   std::string_view);
template LocalAssembly<double> assemble<double>(const OwnerDecomposition<double>&,
                                                const RadialStateGraph<double>&,
                                                const ThresholdJoins<double>&, const Rule&,
                                                Diagnostics*);

}  // namespace sssp::paper::assembly
