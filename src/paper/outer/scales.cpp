// Numerical scales and atoms of the outer construction, Sec 8.2 (Lemma 8.4). See scales.hpp.

#include "paper/outer/scales.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>

#include "paper/local/assemble.hpp"
#include "paper/local/assembly_core.hpp"
#include "paper/local/saturating.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/model/params.hpp"
#include "paper/pivots/name_routes.hpp"
#include "paper/primitives/dsu.hpp"

namespace sssp::paper::outer {

namespace {

template <class W>
constexpr W kInf = WeightTraits<W>::infinity();

[[noreturn]] void broken(const std::string& what) {
  throw std::logic_error("build_outer_scales: " + what);
}

// ---- multiword naturals for the exact h = b + 1 (S1, Q-71): little-endian 32-bit limbs
// without leading zeros. Integer word work on sizes (Sec 2 l.302-304), never on lengths.
using Nat = std::vector<std::uint32_t>;

Nat nat_from(std::uint64_t x) {
  Nat r;
  for (; x != 0; x >>= 32) r.push_back(static_cast<std::uint32_t>(x));
  return r;
}

Nat nat_mul(const Nat& a, const Nat& b) {
  if (a.empty() || b.empty()) return {};
  Nat r(a.size() + b.size(), 0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    std::uint64_t carry = 0;
    for (std::size_t j = 0; j < b.size(); ++j) {
      // a[i] b[j] + r[i + j] + carry <= (2^32 - 1)^2 + 2 (2^32 - 1) = 2^64 - 1: no wrap.
      const std::uint64_t cur = static_cast<std::uint64_t>(a[i]) * b[j] + r[i + j] + carry;
      r[i + j] = static_cast<std::uint32_t>(cur);
      carry = cur >> 32;
    }
    r[i + b.size()] = static_cast<std::uint32_t>(carry);
  }
  while (!r.empty() && r.back() == 0) r.pop_back();
  return r;
}

Nat nat_add(const Nat& a, const Nat& b) {
  Nat r;
  r.reserve(std::max(a.size(), b.size()) + 1);
  std::uint64_t carry = 0;
  for (std::size_t i = 0; i < std::max(a.size(), b.size()) || carry != 0; ++i) {
    std::uint64_t cur = carry;
    if (i < a.size()) cur += a[i];
    if (i < b.size()) cur += b[i];
    r.push_back(static_cast<std::uint32_t>(cur));
    carry = cur >> 32;
  }
  return r;
}

std::uint64_t nat_bit_length(const Nat& a) {
  if (a.empty()) return 0;
  std::uint64_t bits = 32 * static_cast<std::uint64_t>(a.size() - 1);
  for (std::uint32_t top = a.back(); top != 0; top >>= 1) ++bits;
  return bits;
}

/// Members of each class of a partition of 0..n-1 into k classes, by a stable counting sort
/// (Lemma 3.1): items[begin[c] .. begin[c + 1]) are the members of class c in increasing order.
struct Buckets {
  std::vector<std::size_t> begin;
  std::vector<std::uint32_t> items;

  [[nodiscard]] std::span<const std::uint32_t> of(std::size_t c) const {
    return {items.data() + begin[c], begin[c + 1] - begin[c]};
  }
};

template <class Key>
Buckets bucket(std::size_t n, std::size_t k, Key key) {
  auto [begin, items] = assembly::bucket(n, k, std::move(key));
  return {std::move(begin), std::move(items)};
}

/// [Q-67] The separate-grid table of one graph piece in input units: origin * 2^j for the
/// indices of the Sec 7.3 grid at or above the origin (the calendar grid is stored scaled by
/// 2^7 and starts 7 octaves below it), extended by kSeparateGridExtraOctaves doublings at the
/// top. A doubling that does not fit ends the table with the saturated sentinel.
template <class W>
std::vector<W> separate_table(const owner::GraphPiece<W>& piece) {
  if (piece.grid.size() <= owner::kGridOctavesBelow) broken("a calendar grid without its origin");
  const std::size_t count =
      piece.grid.size() - owner::kGridOctavesBelow + kSeparateGridExtraOctaves;
  std::vector<W> table;
  table.reserve(count);
  W x = piece.origin;
  for (std::size_t j = 0; j < count; ++j) {
    table.push_back(x);
    if (x == kInf<W>) break;
    x = sat::add(x, x);
  }
  return table;
}

}  // namespace

OuterScaleFactor outer_scale_factor_power_separator(unsigned log2_p, std::uint64_t base,
                                                    unsigned exponent, unsigned log2_k) {
  if (base == 0) {
    throw std::invalid_argument("outer_scale_factor_power_separator: b must be at least 1");
  }
  constexpr std::uint64_t kMaxBits = std::uint64_t{1} << 20;
  const Nat base_nat = nat_from(base);
  Nat b = nat_from(1);
  if (base > 1) {
    if (nat_bit_length(base_nat) * exponent > kMaxBits) {
      throw std::invalid_argument("outer_scale_factor_power_separator: separator too large");
    }
    for (unsigned i = 0; i < exponent; ++i) b = nat_mul(b, base_nat);
  }
  // l.2080: 2 P (2h+1) (8 K^2 h)^2 = 2^(7 + p + 4k) h^2 (2h+1); the least power of two strictly
  // above 2^c M is 2^(c + bitlen(M)).
  const Nat h = nat_add(b, nat_from(1));
  const Nat m = nat_mul(nat_mul(h, h), nat_add(nat_add(h, h), nat_from(1)));
  const std::uint64_t e = 7 + std::uint64_t{log2_p} + 4 * std::uint64_t{log2_k} + nat_bit_length(m);
  return OuterScaleFactor{
      e, e < 64 ? std::uint64_t{1} << e : std::numeric_limits<std::uint64_t>::max()};
}

OuterScaleFactor paper_outer_scale_factor(unsigned log2_p, std::uint64_t n) {
  // Sec 7.1 l.1644: b = L_N^64; Sec 8.1: K = 16.
  return outer_scale_factor_power_separator(log2_p, model::L(n), model::kOuterExponent,
                                            kOuterShellRatioLog2);
}

template <class W>
std::vector<W> outer_name_weights(const CsrGraph<W>& g,
                                  const owner::OuterEdgeClassification<W>& cls) {
  const std::size_t m = g.num_edges();
  if (cls.edge_class.size() != m || cls.name.size() != m) {
    throw std::invalid_argument("outer_name_weights: classification of a different graph");
  }
  // l.2077: the least physical weight of any macro with that name (Q-65: w_e, macros only).
  std::vector<W> a(cls.dictionary.size(), kInf<W>);
  const auto edges = g.edges();
  for (EdgeId e = 0; e < m; ++e) {
    if (cls.edge_class[e] != owner::OuterEdgeClass::kMacro) continue;
    if (cls.name[e] >= a.size()) {
      throw std::invalid_argument(
          "outer_name_weights: macro edge with a name outside the "
          "dictionary");
    }
    if (edges[e].weight < a[cls.name[e]]) a[cls.name[e]] = edges[e].weight;
  }
  return a;
}

template <class W>
std::vector<std::uint8_t> outer_alignment_by_median_pruning(const OuterScaleFactor& f, W a,
                                                            std::span<const W> s,
                                                            std::uint64_t* tests) {
  // l.2130-2131: F^2 s >= a is monotone in s, so the least aligned scale s* decides the group.
  const std::optional<W> star =
      pivots::median_prune<W>(s, [&](W x) { return outer_block_aligned(f, x, a); }, nullptr, tests);
  std::vector<std::uint8_t> flags(s.size(), 0);
  if (star) {
    for (std::size_t i = 0; i < s.size(); ++i) flags[i] = s[i] < *star ? 0 : 1;
  }
  return flags;
}

template <class W>
OuterScales<W> build_outer_scales(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                  const owner::OuterEdgeClassification<W>& cls,
                                  const owner::AuxGraph<W>& aux, const owner::Calendars<W>& cal,
                                  const OuterStates<W>& os, const OuterScaleFactor& f,
                                  Diagnostics* diag) {
  const auto& st = os.states;
  const auto& joins = os.joins;
  const auto& jr = joins.joins;
  const auto& theta = joins.theta;
  const std::size_t n = g.num_vertices();
  const std::size_t ns = st.num_states();
  const std::size_t nb = st.num_blocks;
  const std::size_t na = aux.num_vertices();
  if (od.num_vertices() != n || st.num_vertices() != n || aux.aux_vertex.size() != n ||
      od.radius.size() != n) {
    throw std::invalid_argument("build_outer_scales: inputs of different graphs");
  }
  if (cal.graph_piece.size() != na || cal.radius_index.size() != na) {
    throw std::invalid_argument("build_outer_scales: calendars of a different auxiliary graph");
  }
  if (joins.num_states != ns || theta.size() != ns || joins.b0_key.size() != ns ||
      st.block.size() != ns || st.projection.size() != ns || joins.num_pivots != od.pivots.size() ||
      os.num_ordinary_blocks > nb || st.arc_kind.size() != st.graph.arcs.size() ||
      st.arc_edge.size() != st.graph.arcs.size()) {
    throw std::invalid_argument("build_outer_scales: outer states of a different construction");
  }

  OuterScales<W> out;
  out.factor = f;

  // S0 (l.2077-2079): new names.
  out.name_weight = outer_name_weights(g, cls);
  const std::size_t nd = out.name_weight.size();

  // Keys of the grid (l.2081-2083): every route name of a finite join (singleton blocks
  // included), in increasing dictionary order, then every corridor of T_A (Q-65).
  out.name_key.assign(nd, kNoScaleKey);
  std::vector<std::uint8_t> used(nd, 0);
  for (StateId z = 0; z < ns; ++z) {
    if (theta[z] == kInf<W>) continue;
    const std::uint32_t d = joins.b0_key[z];
    // l.2078: a finite B0 stores the pointer of an attaining macro, whose name has a weight.
    if (d >= nd || out.name_weight[d] == kInf<W>) broken("finite theta without a macro name");
    used[d] = 1;
  }
  for (std::uint32_t d = 0; d < nd; ++d) {
    if (used[d] == 0) continue;
    out.name_key[d] = static_cast<std::uint32_t>(out.key_weight.size());
    out.key_weight.push_back(out.name_weight[d]);
  }
  out.num_name_keys = out.key_weight.size();
  for (const auto& c : od.corridors) out.key_weight.push_back(c.weight);
  const std::size_t nkeys = out.key_weight.size();

  // Members and block MST joins of every block (counting sorts, Lemma 3.1).
  const Buckets members = bucket(ns, nb, [&](std::size_t x) { return st.block[x]; });
  std::vector<std::uint32_t> mst_ids;
  for (std::uint32_t i = 0; i < jr.size(); ++i) {
    if (jr[i].kind == JoinKind::kBlockMst) mst_ids.push_back(i);
  }
  const Buckets block_mst =
      bucket(mst_ids.size(), nb, [&](std::size_t i) { return st.block[jr[mst_ids[i]].a]; });

  // S2 (Sec 8.1 l.1999, Sec 8.2 l.2083-2086): s_U, z_U and the name of every nonsingleton block.
  out.blocks.assign(nb, OuterBlockScale<W>{});
  std::vector<W> radii;
  for (std::size_t u = 0; u < nb; ++u) {
    auto& bs = out.blocks[u];
    const auto mem = members.of(u);
    bs.size = mem.size();
    if (bs.size <= 1) continue;  // l.2083-2085: only its own finite theta join
    if (u >= os.num_ordinary_blocks) broken("a special state in a nonsingleton block");
    radii.clear();
    for (const StateId x : mem) {
      const Vertex v = st.projection[x];
      if (od.is_pivot(v)) broken("a pivot copy in a nonsingleton block");
      radii.push_back(od.radius[v]);  // Q-35: the original radius, not J
    }
    bs.s = outer_block_scale(std::span<const W>(radii));
    for (const StateId x : mem) {
      const Vertex v = st.projection[x];
      if (od.radius[v] == bs.s && (bs.s_vertex == kNoVertex || v < bs.s_vertex)) bs.s_vertex = v;
      if (theta[x] != kInf<W> && (bs.z == kNoState || theta[x] < theta[bs.z])) bs.z = x;
      // Lemma 8.2 (r <= theta) with s_U <= r: the atoms and separate grids rely on it.
      if (theta[x] < bs.s) broken("theta below s_U in a nonsingleton block");
    }
    bs.kind = OuterBlockKind::kUnattached;
    if (bs.z != kNoState) {
      bs.name = joins.b0_key[bs.z];
      bs.key = out.name_key[bs.name];
      bs.kind = OuterBlockKind::kUnaligned;  // decided in S3
    }
  }

  // S3 (l.2086, l.2130-2131): alignment, grouped by name and decided by median pruning.
  {
    std::vector<std::uint32_t> attached;
    for (std::uint32_t u = 0; u < nb; ++u) {
      if (out.blocks[u].kind == OuterBlockKind::kUnaligned) attached.push_back(u);
    }
    const Buckets by_key = bucket(attached.size(), out.num_name_keys,
                                  [&](std::size_t i) { return out.blocks[attached[i]].key; });
    std::vector<W> scales;
    for (std::size_t k = 0; k < out.num_name_keys; ++k) {
      const auto group = by_key.of(k);
      if (group.empty()) continue;
      scales.clear();
      for (const std::uint32_t i : group) scales.push_back(out.blocks[attached[i]].s);
      const auto flags = outer_alignment_by_median_pruning<W>(
          f, out.key_weight[k], std::span<const W>(scales), &out.alignment_tests);
      for (std::size_t i = 0; i < group.size(); ++i) {
        if (flags[i] != 0) out.blocks[attached[group[i]]].kind = OuterBlockKind::kAligned;
      }
    }
  }

  // S4 (l.2081-2091): supports [a, F^3 a], lowered to s_U by the aligned blocks, and the grid.
  {
    std::vector<OuterAlignedScale<W>> aligned;
    for (const auto& bs : out.blocks) {
      if (bs.kind == OuterBlockKind::kAligned) aligned.push_back({bs.s, bs.key});
    }
    out.grid = build_outer_global_grid<W>(f, std::span<const W>(out.key_weight),
                                          std::span<const OuterAlignedScale<W>>(aligned));
  }
  const auto& points = out.grid.points;
  // Class of a value labelled with a key: binary search between the key's indices (l.2089).
  const auto global_class = [&](W w, std::uint32_t key, const char* what) {
    if (key >= nkeys) broken(std::string(what) + " without a grid key");
    const std::size_t i = outer_grid_class(out.grid, key, w);
    if (i == kNoGridClass) broken(std::string(what) + " outside the grid of its key");
    return i;
  };

  // S5 (l.2091-2093, l.2100-2104): the cutoffs t_U.
  out.piece_table.assign(cal.graph_pieces.size(), {});
  for (auto& bs : out.blocks) {
    if (bs.kind == OuterBlockKind::kSingleton) continue;
    if (bs.kind == OuterBlockKind::kAligned) {
      bs.t = outer_aligned_cutoff(out.grid, bs.key, bs.s);  // floor on its name's grid
      continue;
    }
    // Separate grid: the auxiliary grid pointer of the vertex attaining s_U (Q-67).
    const Vertex x = aux.aux_vertex[bs.s_vertex];
    if (x >= na || cal.graph_piece[x] >= cal.graph_pieces.size()) {
      broken("a block vertex without a graph piece");
    }
    bs.piece = cal.graph_piece[x];
    auto& table = out.piece_table[bs.piece];
    if (table.empty()) table = separate_table(cal.graph_pieces[bs.piece]);
    // grid[i] <= 2^7 r_v < 2 grid[i] with grid[7 + j] = 2^7 origin 2^j: the floor of r_v in
    // input units is table[i - 7].
    if (cal.radius_index[x] < owner::kGridOctavesBelow ||
        cal.radius_index[x] - owner::kGridOctavesBelow >= table.size()) {
      broken("a radius index outside the piece table");
    }
    bs.table_index = cal.radius_index[x] - owner::kGridOctavesBelow;
    const std::size_t width = std::min<std::size_t>(2, table.size() - bs.table_index);
    bs.t = outer_separate_cutoff(std::span<const W>(table).subspan(bs.table_index, width), bs.s);
    if (bs.t != table[bs.table_index]) broken("the grid pointer is not the floor of s_U");
  }

  // S6 (Sec 8.1 l.2000-2001; Lemma 8.4 item 3): the components of the block MST edges lighter
  // than t_U. A join below the cutoff is internal to an atom; the others are classified in S7.
  out.events.assign(jr.size(), OuterEvent<W>{});
  out.state_atom.assign(ns, kNoAtom);
  std::vector<std::uint32_t> light;  // MST joins below their cutoff
  {
    DisjointSets dsu(ns);
    for (const std::uint32_t j : mst_ids) {
      const auto& bs = out.blocks[st.block[jr[j].a]];
      if (jr[j].weight < bs.t) {
        dsu.unite(jr[j].a, jr[j].b);
        light.push_back(j);
      }
    }
    // A component has at least two states iff it contains a light edge.
    std::vector<std::uint8_t> has_light_edge(ns, 0);
    for (const std::uint32_t j : light) has_light_edge[dsu.find(jr[j].a)] = 1;
    std::vector<std::uint32_t> atom_of_root(ns, kNoAtom);
    for (StateId x = 0; x < ns; ++x) {
      if (out.blocks[st.block[x]].kind == OuterBlockKind::kSingleton) continue;
      const auto root = dsu.find(x);
      if (has_light_edge[root] == 0) continue;  // one state: a singleton part (Q-74)
      if (atom_of_root[root] == kNoAtom) {
        atom_of_root[root] = static_cast<std::uint32_t>(out.atoms.size());
        out.atoms.emplace_back();
        out.atoms.back().cutoff = out.blocks[st.block[x]].t;
        out.atom_block.push_back(st.block[x]);
      }
      out.state_atom[x] = atom_of_root[root];
      out.atoms[atom_of_root[root]].states.push_back(x);
    }
    if (diag != nullptr) SSSP_COUNT(dsu.add_counters(*diag));
  }
  const std::size_t natoms = out.atoms.size();
  {
    std::vector<std::uint32_t> local(ns, 0);
    for (const auto& atom : out.atoms) {
      for (std::size_t i = 0; i < atom.states.size(); ++i) {
        local[atom.states[i]] = static_cast<std::uint32_t>(i);
      }
    }
    // Internal graph H (Q-74): every within-window edge with both endpoints in the atom, once
    // per input edge (the arc from the smaller to the larger state). No other arc may lie
    // inside an atom: transfers and switches leave the window or the phase, pivot and special
    // arcs touch singleton blocks, and a macro joins two owner cells.
    const auto& arcs = st.graph.arcs;
    std::vector<std::uint32_t> inside;
    for (std::uint32_t a = 0; a < arcs.size(); ++a) {
      const std::uint32_t at = out.state_atom[arcs[a].tail];
      if (at == kNoAtom || at != out.state_atom[arcs[a].head]) continue;
      if (st.arc_kind[a] != RadialArcKind::kWithinWindow || arcs[a].tail == arcs[a].head) {
        broken("an arc inside an atom that is not a within-window edge");
      }
      if (arcs[a].tail < arcs[a].head) inside.push_back(a);
    }
    const Buckets atom_arcs = bucket(
        inside.size(), natoms, [&](std::size_t i) { return out.state_atom[arcs[inside[i]].tail]; });
    const Buckets atom_tree =
        bucket(light.size(), natoms, [&](std::size_t i) { return out.state_atom[jr[light[i]].a]; });
    // Slot of an input edge inside the atom being built (an atom projects injectively, so an
    // input edge occurs at most once in it); verified through the endpoints.
    std::vector<EdgeId> slot(g.num_edges(), kNoEdge);
    for (std::uint32_t y = 0; y < natoms; ++y) {
      auto& atom = out.atoms[y];
      for (const std::uint32_t i : atom_arcs.of(y)) {
        const std::uint32_t a = inside[i];
        if (st.arc_edge[a] >= slot.size()) broken("a within-window arc without an input edge");
        slot[st.arc_edge[a]] = static_cast<EdgeId>(atom.edges.size());
        atom.edges.push_back({local[arcs[a].tail], local[arcs[a].head], arcs[a].length});
      }
      for (const std::uint32_t i : atom_tree.of(y)) {
        const std::uint32_t j = light[i];
        const EdgeId e = jr[j].ref < arcs.size() ? st.arc_edge[jr[j].ref] : kNoEdge;
        const EdgeId id = e < slot.size() ? slot[e] : kNoEdge;
        const Vertex lo = local[std::min(jr[j].a, jr[j].b)];
        const Vertex hi = local[std::max(jr[j].a, jr[j].b)];
        if (id >= atom.edges.size() || atom.edges[id].u != lo || atom.edges[id].v != hi ||
            atom.edges[id].weight != jr[j].weight) {
          broken("a block MST edge without its within-window edge in the atom");
        }
        atom.spanning_tree.push_back(id);
        out.events[j] = {OuterEventGrid::kAtom, y, W{0}, atom.cutoff};
      }
      if (atom.spanning_tree.size() + 1 != atom.states.size()) {
        broken("the light MST edges do not span their atom");
      }
    }
  }

  // S7 (l.2091-2094): classes of the significant block MST edges and of the finite joins.
  std::vector<std::uint8_t> is_light(jr.size(), 0);
  for (const std::uint32_t j : light) is_light[j] = 1;
  std::vector<W> separate_top(nb, W{0});  // lower end of the top separate class of a block
  std::vector<std::uint8_t> has_separate(nb, 0);
  for (std::uint32_t j = 0; j < jr.size(); ++j) {
    if (is_light[j] != 0) continue;
    const W w = jr[j].weight;
    auto& ev = out.events[j];
    std::uint32_t key = kNoScaleKey;
    switch (jr[j].kind) {
      case JoinKind::kBlockMst: {
        const std::uint32_t u = st.block[jr[j].a];
        const auto& bs = out.blocks[u];
        if (bs.kind == OuterBlockKind::kAligned) {
          key = bs.key;  // global class in the index range of a_U
          break;
        }
        // Separate grid: the floor of w among the table points from t_U up (l.2100-2104).
        const auto& table = out.piece_table[bs.piece];
        const auto first = table.begin() + static_cast<std::ptrdiff_t>(bs.table_index);
        const auto it = std::upper_bound(first, table.end(), w);
        if (it == table.end()) broken("a block MST edge above the extended piece table");
        const auto pos = static_cast<std::size_t>(it - table.begin()) - 1;
        ev = {OuterEventGrid::kSeparate, pos - bs.table_index, table[pos], table[pos + 1]};
        if (separate_top[u] < ev.lower) separate_top[u] = ev.lower;
        has_separate[u] = 1;
        continue;
      }
      case JoinKind::kOwner:
        key = joins.b0_key[jr[j].a] < nd ? out.name_key[joins.b0_key[jr[j].a]] : kNoScaleKey;
        break;
      case JoinKind::kCorridor:
        key = static_cast<std::uint32_t>(out.num_name_keys + jr[j].ref);
        break;
    }
    const std::size_t i = global_class(w, key, "a join");
    ev = {OuterEventGrid::kGlobal, i, points[i], points[i + 1]};
  }
  // Unaligned blocks (l.2086-2087, l.2093-2094; Lemma 8.4 item 1): every MST edge is below
  // a_U / 2, and the block is complete before its first global class.
  std::vector<W> weights;
  for (std::uint32_t u = 0; u < nb; ++u) {
    const auto& bs = out.blocks[u];
    if (bs.kind != OuterBlockKind::kUnaligned) continue;
    weights.clear();
    for (const std::uint32_t i : block_mst.of(u)) weights.push_back(jr[mst_ids[i]].weight);
    if (first_edge_not_below_half(std::span<const W>(weights), out.key_weight[bs.key]) !=
        weights.size()) {
      broken("an unaligned block with an MST edge not below a_U / 2");
    }
    const W first_global = points[global_class(theta[bs.z], bs.key, "theta(z_U)")];
    if (has_separate[u] != 0 && first_global < sat::add(separate_top[u], separate_top[u])) {
      broken("a separate-grid class reaches the first attachment class");
    }
    if (first_global < bs.t) broken("the first attachment class lies below t_U");
  }

  if (diag != nullptr) {
    std::size_t kinds[4] = {0, 0, 0, 0};
    for (const auto& bs : out.blocks) ++kinds[static_cast<std::size_t>(bs.kind)];
    std::size_t grids[3] = {0, 0, 0};
    for (const auto& ev : out.events) ++grids[static_cast<std::size_t>(ev.grid)];
    std::size_t atom_states = 0;
    std::size_t atom_edges = 0;
    for (const auto& atom : out.atoms) {
      atom_states += atom.states.size();
      atom_edges += atom.edges.size();
    }
    SSSP_COUNT(++diag->counters["outer_scales.calls"]);
    SSSP_COUNT(diag->counters["outer_scales.scale_factor_log2"] = f.log2);
    SSSP_COUNT(diag->counters["outer_scales.name_keys"] += out.num_name_keys);
    SSSP_COUNT(diag->counters["outer_scales.blocks_singleton"] += kinds[0]);
    SSSP_COUNT(diag->counters["outer_scales.blocks_aligned"] += kinds[1]);
    SSSP_COUNT(diag->counters["outer_scales.blocks_unaligned"] += kinds[2]);
    SSSP_COUNT(diag->counters["outer_scales.blocks_unattached"] += kinds[3]);
    SSSP_COUNT(diag->counters["outer_scales.alignment_tests"] += out.alignment_tests);
    SSSP_COUNT(diag->counters["outer_scales.grid_points"] += points.size());
    SSSP_COUNT(diag->counters["outer_scales.grid_components"] += out.grid.component_begin.size());
    SSSP_COUNT(diag->counters["outer_scales.atoms"] += natoms);
    SSSP_COUNT(diag->counters["outer_scales.atom_states"] += atom_states);
    SSSP_COUNT(diag->counters["outer_scales.atom_edges"] += atom_edges);
    SSSP_COUNT(diag->counters["outer_scales.events_atom"] += grids[0]);
    SSSP_COUNT(diag->counters["outer_scales.events_separate"] += grids[1]);
    SSSP_COUNT(diag->counters["outer_scales.events_global"] += grids[2]);
    static_cast<void>(kinds);
    static_cast<void>(grids);
    static_cast<void>(atom_states);
    static_cast<void>(atom_edges);
  }
  return out;
}

template <class W>
std::vector<std::string> check_outer_scales(const CsrGraph<W>& g, const OwnerDecomposition<W>& od,
                                            const owner::OuterEdgeClassification<W>& cls,
                                            const owner::AuxGraph<W>& aux,
                                            const owner::Calendars<W>& cal,
                                            const OuterStates<W>& os, const OuterScales<W>& sc) {
  std::vector<std::string> bad;
  const auto fail = [&](const std::string& what) { bad.push_back(what); };
  const auto& st = os.states;
  const auto& joins = os.joins;
  const auto& jr = joins.joins;
  const auto& theta = joins.theta;
  const auto& f = sc.factor;
  const std::size_t ns = st.num_states();
  const std::size_t nb = st.num_blocks;
  const std::size_t nd = cls.dictionary.size();
  const auto edges = g.edges();
  if (sc.name_weight.size() != nd || sc.name_key.size() != nd || sc.blocks.size() != nb ||
      sc.state_atom.size() != ns || sc.events.size() != jr.size() ||
      sc.key_weight.size() != sc.num_name_keys + od.corridors.size() ||
      sc.grid.supports.size() != sc.key_weight.size() || sc.atom_block.size() != sc.atoms.size() ||
      sc.piece_table.size() != cal.graph_pieces.size()) {
    fail("sizes do not match the inputs");
    return bad;
  }
  const auto id = [](const char* what, std::size_t i) {
    return std::string(what) + " " + std::to_string(i) + ": ";
  };
  const auto& points = sc.grid.points;

  // S0: a <= w < 2 P (2h+1) a on every macro edge (Lemma 7.1), and a is attained.
  {
    std::vector<std::uint8_t> attained(nd, 0);
    for (EdgeId e = 0; e < edges.size(); ++e) {
      if (cls.edge_class[e] != owner::OuterEdgeClass::kMacro) continue;
      const std::uint32_t d = cls.name[e];
      if (d >= nd) {
        fail(id("edge", e) + "macro without a name");
        continue;
      }
      const W a = sc.name_weight[d];
      const W w = edges[e].weight;
      if (w < a) fail(id("edge", e) + "new name above the macro weight");
      if (w == a) attained[d] = 1;
      // 2 P a by doublings, then (2h + 1) by the binary expansion; saturating.
      std::uint64_t two_h_plus_one = std::numeric_limits<std::uint64_t>::max();
      if (od.hop_bound < (std::numeric_limits<std::uint64_t>::max() - 1) / 2) {
        two_h_plus_one = 2 * static_cast<std::uint64_t>(od.hop_bound) + 1;
      }
      const W bound = sat::mul_int(scale_pow2(a, std::uint64_t{cls.p_log2} + 1), two_h_plus_one);
      if (!(w < bound)) fail(id("edge", e) + "macro weight not below 2P(2h+1) times its name");
    }
    for (std::uint32_t d = 0; d < nd; ++d) {
      if ((sc.name_weight[d] != kInf<W>) != (attained[d] != 0)) {
        fail(id("name", d) + "weight not attained by a macro edge");
      }
    }
  }
  // Keys and supports: a <= B0 <= theta < F a, theta inside the support of its own key.
  for (StateId z = 0; z < ns; ++z) {
    if (theta[z] == kInf<W>) continue;
    const std::uint32_t d = joins.b0_key[z];
    if (d >= nd || sc.name_key[d] >= sc.num_name_keys) {
      fail(id("state", z) + "finite theta without a key");
      continue;
    }
    const W a = sc.name_weight[d];
    if (sc.key_weight[sc.name_key[d]] != a) fail(id("state", z) + "key weight is not the name");
    if (joins.b0[z] < a || theta[z] < joins.b0[z]) fail(id("state", z) + "a <= B0 <= theta fails");
    if (!(theta[z] < scale_pow2(a, f.log2))) fail(id("state", z) + "theta < F a fails");
    const auto& sp = sc.grid.supports[sc.name_key[d]];
    if (theta[z] < sp.lower || sp.upper < theta[z]) fail(id("state", z) + "theta outside support");
  }
  for (std::size_t c = 0; c < od.corridors.size(); ++c) {
    const auto& sp = sc.grid.supports[sc.num_name_keys + c];
    const W w = od.corridors[c].weight;
    if (sc.key_weight[sc.num_name_keys + c] != w) fail(id("corridor", c) + "key weight differs");
    if (w < sp.lower || sp.upper < w) fail(id("corridor", c) + "weight outside its support");
  }
  // Grid: strictly increasing; doubling inside a component.
  for (std::size_t c = 0; c < sc.grid.component_begin.size(); ++c) {
    const std::size_t lo = sc.grid.component_begin[c];
    const std::size_t hi =
        c + 1 < sc.grid.component_begin.size() ? sc.grid.component_begin[c + 1] : points.size();
    if (lo >= hi || hi > points.size()) {
      fail(id("grid component", c) + "empty or out of range");
      continue;
    }
    for (std::size_t i = lo; i + 1 < hi; ++i) {
      if (!(points[i] < points[i + 1]) || points[i + 1] != sat::add(points[i], points[i])) {
        fail(id("grid point", i) + "not doubled inside its component");
      }
    }
    // The next component starts above this one's last own point (at its sentinel, if shared).
    if (hi < points.size() && !(points[hi - 1] < points[hi])) {
      fail(id("grid component", c) + "next component does not start above it");
    }
  }

  // Blocks.
  std::vector<std::vector<StateId>> members(nb);
  for (StateId x = 0; x < ns; ++x) members[st.block[x]].push_back(x);
  std::vector<std::vector<std::uint32_t>> block_mst(nb);
  for (std::uint32_t j = 0; j < jr.size(); ++j) {
    if (jr[j].kind == JoinKind::kBlockMst) block_mst[st.block[jr[j].a]].push_back(j);
  }
  for (std::size_t u = 0; u < nb; ++u) {
    const auto& bs = sc.blocks[u];
    const auto& mem = members[u];
    if (bs.size != mem.size()) fail(id("block", u) + "size");
    if ((mem.size() <= 1) != (bs.kind == OuterBlockKind::kSingleton)) {
      fail(id("block", u) + "singleton kind");
      continue;
    }
    if (mem.size() <= 1) continue;
    if (u >= os.num_ordinary_blocks) fail(id("block", u) + "nonsingleton special block");
    W s = kInf<W>;
    Vertex sv = kNoVertex;
    StateId z = kNoState;
    for (const StateId x : mem) {
      const Vertex v = st.projection[x];
      if (od.is_pivot(v)) fail(id("block", u) + "contains a pivot copy");
      if (od.radius[v] < s || (od.radius[v] == s && v < sv)) {
        s = od.radius[v];
        sv = v;
      }
      if (theta[x] != kInf<W> && (z == kNoState || theta[x] < theta[z])) z = x;
    }
    if (bs.s != s || bs.s_vertex != sv) fail(id("block", u) + "s_U is not the least radius");
    if (!(W{0} < bs.s)) fail(id("block", u) + "s_U not positive");
    if (bs.z != z) fail(id("block", u) + "z_U does not minimise theta");
    for (const StateId x : mem) {
      if (theta[x] < bs.s) fail(id("block", u) + "theta below s_U");
    }
    if (!outer_cutoff_valid(bs.s, bs.t)) fail(id("block", u) + "s/2 < t <= s fails");
    if ((z == kNoState) != (bs.kind == OuterBlockKind::kUnattached)) {
      fail(id("block", u) + "attachment kind");
      continue;
    }
    if (z != kNoState) {
      if (bs.name != joins.b0_key[z] || bs.name >= nd || bs.key != sc.name_key[bs.name] ||
          bs.key >= sc.num_name_keys) {
        fail(id("block", u) + "name of z_U");
        continue;
      }
      const W a = sc.name_weight[bs.name];
      const bool aligned = !(scale_pow2(bs.s, 2 * f.log2) < a);
      if (aligned != (bs.kind == OuterBlockKind::kAligned)) fail(id("block", u) + "alignment");
      if (aligned) {
        if (!outer_aligned_scale_in_range(f, bs.s, a)) fail(id("block", u) + "s > F a");
        const auto& sp = sc.grid.supports[bs.key];
        if (bs.s < sp.lower) fail(id("block", u) + "support start above s_U");
        const auto it = std::lower_bound(points.begin(), points.end(), bs.t);
        const auto pos = static_cast<std::size_t>(it - points.begin());
        if (it == points.end() || *it != bs.t || pos < sp.lower_index || pos >= sp.upper_index) {
          fail(id("block", u) + "t_U is not a grid point of its name");
        }
        for (const StateId x : mem) {
          // Every finite theta of an aligned block lies in the support of the block's name.
          if (theta[x] != kInf<W> && (theta[x] < sp.lower || sp.upper < theta[x])) {
            fail(id("block", u) + "theta outside the support of the block name");
          }
        }
      } else {
        for (const std::uint32_t j : block_mst[u]) {
          if (!(sat::add(jr[j].weight, jr[j].weight) < a)) {
            fail(id("block", u) + "unaligned MST edge not below a/2");
          }
        }
        for (const StateId x : mem) {
          if (theta[x] != theta[z]) fail(id("block", u) + "unaligned block with two thetas");
        }
      }
    }
    if (bs.kind != OuterBlockKind::kAligned) {
      // Separate grid: t_U = origin 2^j of the piece of s_vertex.
      const Vertex x = aux.aux_vertex[bs.s_vertex];
      if (bs.piece != cal.graph_piece[x] || bs.piece >= cal.graph_pieces.size()) {
        fail(id("block", u) + "piece of the separate grid");
        continue;
      }
      if (aux.radius[x] != bs.s) fail(id("block", u) + "auxiliary radius differs from s_U");
      if (scale_pow2(cal.graph_pieces[bs.piece].origin, bs.table_index) != bs.t) {
        fail(id("block", u) + "t_U is not origin 2^j");
      }
      W top{0};
      bool any = false;
      for (const std::uint32_t j : block_mst[u]) {
        if (sc.events[j].grid != OuterEventGrid::kSeparate) continue;
        any = true;
        top = std::max(top, sc.events[j].lower);
        if (scale_pow2(bs.t, sc.events[j].index) != sc.events[j].lower) {
          fail(id("join", j) + "separate class is not t_U 2^j");
        }
      }
      if (any && z != kNoState && bs.key < sc.num_name_keys) {
        const std::size_t i = outer_grid_class(sc.grid, bs.key, theta[z]);
        if (i == kNoGridClass || points[i] < sat::add(top, top)) {
          fail(id("block", u) + "2 p_top <= g fails");
        }
      }
    }
  }

  // Atoms: the components of the block MST edges below t_U with at least two states.
  {
    DisjointSets dsu(ns);
    for (std::uint32_t j = 0; j < jr.size(); ++j) {
      if (jr[j].kind != JoinKind::kBlockMst) continue;
      const bool below = jr[j].weight < sc.blocks[st.block[jr[j].a]].t;
      if (below != (sc.events[j].grid == OuterEventGrid::kAtom)) {
        fail(id("join", j) + "atom class does not match w < t_U");
      }
      if (below) dsu.unite(jr[j].a, jr[j].b);
    }
    std::vector<std::uint32_t> count(ns, 0);
    for (StateId x = 0; x < ns; ++x) ++count[dsu.find(x)];
    for (StateId x = 0; x < ns; ++x) {
      const bool in_atom = count[dsu.find(x)] >= 2;
      if (in_atom != (sc.state_atom[x] != kNoAtom)) fail(id("state", x) + "atom membership");
    }
    for (std::uint32_t j = 0; j < jr.size(); ++j) {
      if (jr[j].kind != JoinKind::kBlockMst) continue;
      if ((dsu.find(jr[j].a) == dsu.find(jr[j].b)) !=
          (sc.state_atom[jr[j].a] != kNoAtom && sc.state_atom[jr[j].a] == sc.state_atom[jr[j].b])) {
        fail(id("join", j) + "atoms are not the light components");
      }
    }
  }
  std::vector<std::size_t> inside_arcs(sc.atoms.size(), 0);
  for (std::size_t a = 0; a < st.graph.arcs.size(); ++a) {
    const auto& arc = st.graph.arcs[a];
    const std::uint32_t y = sc.state_atom[arc.tail];
    if (y == kNoAtom || y != sc.state_atom[arc.head]) continue;
    if (st.arc_kind[a] != RadialArcKind::kWithinWindow) {
      fail(id("arc", a) + "inside an atom but not within a window");
    }
    ++inside_arcs[y];
  }
  for (std::size_t y = 0; y < sc.atoms.size(); ++y) {
    const auto& atom = sc.atoms[y];
    const std::size_t k = atom.states.size();
    if (k < 2) {
      fail(id("atom", y) + "fewer than two states");
      continue;
    }
    // Lemma 8.4: fewer than b vertices (one component of T - A).
    if (k >= od.separator) fail(id("atom", y) + "not fewer than b states");
    const auto& bs = sc.blocks[sc.atom_block[y]];
    if (atom.cutoff != bs.t) fail(id("atom", y) + "cutoff is not t_U");
    const Vertex v0 = st.projection[atom.states[0]];
    std::vector<Vertex> images;
    for (std::size_t i = 0; i < k; ++i) {
      const StateId x = atom.states[i];
      if (i > 0 && !(atom.states[i - 1] < x)) fail(id("atom", y) + "states not increasing");
      if (sc.state_atom[x] != y) fail(id("atom", y) + "state_atom differs");
      if (st.block[x] != sc.atom_block[y]) fail(id("atom", y) + "spans two blocks");
      const Vertex v = st.projection[x];
      images.push_back(v);
      // Lemma 8.4 item 3: one owner, one original radius, at least t_U.
      if (od.owner[v] != od.owner[v0]) fail(id("atom", y) + "two owners");
      if (od.radius[v] != od.radius[v0]) fail(id("atom", y) + "two radii");
      if (od.radius[v] < bs.s || bs.s < bs.t) fail(id("atom", y) + "r >= s_U >= t_U fails");
      if (theta[x] < bs.t) fail(id("atom", y) + "theta below the cutoff");
    }
    std::ranges::sort(images);
    if (std::ranges::adjacent_find(images) != images.end()) {
      fail(id("atom", y) + "projection not injective");
    }
    // H = the within-window edges inside the atom, once each (two arcs per edge).
    if (2 * atom.edges.size() != inside_arcs[y]) fail(id("atom", y) + "internal edge count");
    if (atom.spanning_tree.size() + 1 != k) fail(id("atom", y) + "spanning tree size");
    DisjointSets tree(k);
    for (const EdgeId e : atom.spanning_tree) {
      if (e >= atom.edges.size()) {
        fail(id("atom", y) + "tree edge out of range");
        continue;
      }
      const auto& ed = atom.edges[e];
      if (ed.u >= k || ed.v >= k || !tree.unite(ed.u, ed.v)) {
        fail(id("atom", y) + "spanning tree has a cycle");
      }
      if (!(ed.weight < atom.cutoff)) fail(id("atom", y) + "tree edge not below the cutoff");
    }
  }

  // Classes bracket their weights.
  for (std::uint32_t j = 0; j < jr.size(); ++j) {
    const auto& ev = sc.events[j];
    if (ev.grid == OuterEventGrid::kAtom) {
      if (jr[j].kind != JoinKind::kBlockMst) fail(id("join", j) + "atom class on a join");
      continue;
    }
    const W w = jr[j].weight;
    if (w < ev.lower || !(w < ev.upper)) fail(id("join", j) + "class does not bracket the weight");
    if (ev.upper != sat::add(ev.lower, ev.lower)) fail(id("join", j) + "class is not [g, 2g)");
    if (ev.grid == OuterEventGrid::kGlobal) {
      if (ev.index + 1 >= points.size() || points[ev.index] != ev.lower ||
          points[ev.index + 1] != ev.upper) {
        fail(id("join", j) + "global class is not a grid class");
      }
      if (jr[j].kind == JoinKind::kBlockMst &&
          sc.blocks[st.block[jr[j].a]].kind != OuterBlockKind::kAligned) {
        fail(id("join", j) + "global class for a block on a separate grid");
      }
    } else if (jr[j].kind != JoinKind::kBlockMst ||
               sc.blocks[st.block[jr[j].a]].kind == OuterBlockKind::kAligned) {
      fail(id("join", j) + "separate class outside a separate-grid block");
    }
  }
  return bad;
}

template <class W>
OuterMergePlan<W> build_outer_merge_plan(const OuterStates<W>& os, const OuterScales<W>& sc) {
  const auto& st = os.states;
  const auto& jr = os.joins.joins;
  const std::size_t ns = st.num_states();
  const std::size_t nb = st.num_blocks;
  if (sc.state_atom.size() != ns || sc.events.size() != jr.size() || sc.blocks.size() != nb) {
    throw std::invalid_argument("build_outer_merge_plan: scales of different outer states");
  }
  OuterMergePlan<W> plan;

  // Parts (Def 5.12): atoms and the states outside them, numbered by their least state.
  plan.part.assign(ns, 0);
  std::vector<std::uint32_t> part_of_atom(sc.atoms.size(), kNoAtom);
  for (StateId x = 0; x < ns; ++x) {
    const std::uint32_t y = sc.state_atom[x];
    if (y != kNoAtom && y >= sc.atoms.size()) {
      throw std::invalid_argument("build_outer_merge_plan: atom id out of range");
    }
    if (y != kNoAtom && part_of_atom[y] != kNoAtom) {
      plan.part[x] = part_of_atom[y];
      continue;
    }
    plan.part[x] = static_cast<std::uint32_t>(plan.part_state.size());
    if (y != kNoAtom) part_of_atom[y] = plan.part[x];
    plan.part_state.push_back(y == kNoAtom ? x : kNoState);
    plan.part_atom.push_back(y);
  }
  plan.num_parts = plan.part_state.size();

  // D: the outside arcs.
  plan.graph.num_states = ns;
  for (ArcId a = 0; a < st.graph.arcs.size(); ++a) {
    const auto& arc = st.graph.arcs[a];
    const std::uint32_t y = sc.state_atom[arc.tail];
    if (y != kNoAtom && y == sc.state_atom[arc.head]) continue;  // internal to an atom
    plan.graph.arcs.push_back(arc);
    plan.arc_of.push_back(a);
  }

  // Classes. Separate grids: joins bucketed by block, then by class index inside the block.
  std::vector<std::uint32_t> separate;
  std::vector<std::uint32_t> global;
  for (std::uint32_t j = 0; j < jr.size(); ++j) {
    if (sc.events[j].grid == OuterEventGrid::kSeparate) separate.push_back(j);
    if (sc.events[j].grid == OuterEventGrid::kGlobal) global.push_back(j);
  }
  const auto add_classes = [&](std::span<const std::uint32_t> js, std::size_t num_classes,
                               std::uint32_t block) {
    const Buckets by_class =
        bucket(js.size(), num_classes, [&](std::size_t i) { return sc.events[js[i]].index; });
    for (std::size_t c = 0; c < num_classes; ++c) {
      const auto items = by_class.of(c);
      if (items.empty()) continue;
      const auto& ev = sc.events[js[items.front()]];
      OuterMergeClass<W> cls{ev.grid, block, c, ev.lower, ev.upper, {}};
      cls.joins.reserve(items.size());
      for (const std::uint32_t i : items) cls.joins.push_back(js[i]);
      plan.classes.push_back(std::move(cls));
    }
  };
  {
    const Buckets by_block =
        bucket(separate.size(), nb, [&](std::size_t i) { return st.block[jr[separate[i]].a]; });
    std::vector<std::uint32_t> js;
    for (std::uint32_t u = 0; u < nb; ++u) {
      const auto items = by_block.of(u);
      if (items.empty()) continue;
      js.clear();
      std::size_t num_classes = 0;
      for (const std::uint32_t i : items) {
        js.push_back(separate[i]);
        num_classes = std::max(num_classes, sc.events[separate[i]].index + 1);
      }
      add_classes(js, num_classes, u);
    }
  }
  // Global classes in increasing grid index: all keys share the one sorted table (l.2089).
  add_classes(global, sc.grid.points.size(), kNoBlock);
  return plan;
}

template <class W>
OuterPartition<W> build_outer_partition(const OuterStates<W>& os, const OuterScales<W>& sc,
                                        Diagnostics* diag) {
  constexpr const char* kContext = "build_outer_partition";
  const auto& st = os.states;
  const auto& jr = os.joins.joins;
  const std::size_t ns = st.num_states();
  OuterMergePlan<W> plan = build_outer_merge_plan(os, sc);
  const std::size_t np = plan.num_parts;

  OuterPartition<W> out;
  auto& ph = out.hierarchy;
  ph.projection = st.projection;
  ph.atoms = sc.atoms;
  ph.c0 = 2;  // Q-68: height a, links below 2a (as Lemma 5.10)
  // Leaves: one per part (Def 5.12); node id = part id.
  ph.nodes.resize(np);
  for (std::size_t p = 0; p < np; ++p) ph.nodes[p].state = plan.part_state[p];
  ph.leaf_atom = plan.part_atom;
  out.node_upper.assign(np, W{0});
  out.node_grid.assign(np, OuterEventGrid::kAtom);
  out.leaf_of_state.assign(plan.part.begin(), plan.part.end());

  // Join vertices of the merges (Lemma 3.6): the parts, then one auxiliary vertex per pivot.
  // The states of an atom share one vertex, so the atom enters every merge as a single leaf.
  std::vector<NodeId> node_of_vertex(np + os.joins.num_pivots, kNoNode);
  for (std::size_t p = 0; p < np; ++p) node_of_vertex[p] = static_cast<NodeId>(p);
  const auto vertex = [&](std::uint32_t v) {
    return v < ns ? plan.part[v] : static_cast<std::uint32_t>(np + (v - ns));
  };
  const auto endpoints = [&](std::uint32_t j) {
    return std::pair<std::uint32_t, std::uint32_t>{vertex(jr[j].a), vertex(jr[j].b)};
  };
  const auto add_node = [&](W height, std::vector<NodeId> children, W upper, bool /*inclusive*/,
                            NodeOrigin origin) {
    const auto id = static_cast<NodeId>(ph.nodes.size());
    for (const NodeId c : children) ph.nodes[c].parent = id;
    HierarchyNode<W> node;
    node.height = height;
    node.children = std::move(children);
    ph.nodes.push_back(std::move(node));
    ph.leaf_atom.push_back(kNoAtom);
    out.node_upper.push_back(upper);
    out.node_grid.push_back(origin == NodeOrigin::kSeparateGrid ? OuterEventGrid::kSeparate
                                                                : OuterEventGrid::kGlobal);
    return id;
  };
  assembly::ClassMerger merger(std::move(node_of_vertex), endpoints, add_node);
  std::size_t separate_classes = 0;
  // l.2094: unions processed by their integer grid classes, in the order of the plan.
  for (const auto& cls : plan.classes) {
    const bool separate = cls.grid == OuterEventGrid::kSeparate;
    separate_classes += separate ? 1 : 0;
    merger.merge_class(std::span<const std::uint32_t>(cls.joins), cls.height, cls.upper, false,
                       separate ? NodeOrigin::kSeparateGrid : NodeOrigin::kGlobalGrid);
  }

  // Def 4.1 H2 and Def 5.12 on the forest: the order "separate grids first" is sound only if no
  // height decreases towards the root and no atom is cut above its parent (l.2093-2094).
  for (std::size_t v = 0; v < ph.nodes.size(); ++v) {
    const NodeId par = ph.nodes[v].parent;
    if (par == kNoNode) continue;
    const W below = v < np ? (ph.leaf_atom[v] == kNoAtom ? W{0} : ph.atoms[ph.leaf_atom[v]].cutoff)
                           : ph.nodes[v].height;
    if (ph.nodes[par].height < below) {
      throw std::logic_error(std::string(kContext) +
                             ": a height or an atom cutoff above the parent's height");
    }
  }

  // l.2096-2099: connecting trees by the child-quotient searches of Sec 5.4 and the child
  // order by the arcs lighter than each node.
  const std::vector<std::uint8_t> inclusive(ph.nodes.size(), 0);
  auto lo = assembly::link_and_order<W>(ph.nodes, std::span<const StateArc<W>>(plan.graph.arcs),
                                        std::span<const NodeId>(out.leaf_of_state),
                                        std::span<const W>(out.node_upper),
                                        std::span<const std::uint8_t>(inclusive), kContext);
  out.roots = std::move(lo.roots);
  ph.graph = std::move(plan.graph);
  out.arc_of = std::move(plan.arc_of);

  if (diag != nullptr) {
    SSSP_COUNT(++diag->counters["outer_scales.partition_calls"]);
    SSSP_COUNT(diag->counters["outer_scales.partition_parts"] += np);
    SSSP_COUNT(diag->counters["outer_scales.partition_nodes"] += ph.nodes.size() - np);
    SSSP_COUNT(diag->counters["outer_scales.partition_links"] += lo.links);
    SSSP_COUNT(diag->counters["outer_scales.partition_roots"] += out.roots.size());
    SSSP_COUNT(diag->counters["outer_scales.separate_classes"] += separate_classes);
    SSSP_COUNT(diag->counters["outer_scales.global_classes"] +=
               plan.classes.size() - separate_classes);
    SSSP_COUNT(merger.dsu().add_counters(*diag));
  }
  static_cast<void>(separate_classes);
  return out;
}

#define SSSP_INSTANTIATE_OUTER_SCALES(W)                                                          \
  template OuterPartition<W> build_outer_partition<W>(const OuterStates<W>&,                      \
                                                      const OuterScales<W>&, Diagnostics*);       \
  template OuterMergePlan<W> build_outer_merge_plan<W>(const OuterStates<W>&,                     \
                                                       const OuterScales<W>&);                    \
  template std::vector<W> outer_name_weights<W>(const CsrGraph<W>&,                               \
                                                const owner::OuterEdgeClassification<W>&);        \
  template std::vector<std::uint8_t> outer_alignment_by_median_pruning<W>(                        \
      const OuterScaleFactor&, W, std::span<const W>, std::uint64_t*);                            \
  template OuterScales<W> build_outer_scales<W>(                                                  \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, const owner::OuterEdgeClassification<W>&, \
      const owner::AuxGraph<W>&, const owner::Calendars<W>&, const OuterStates<W>&,               \
      const OuterScaleFactor&, Diagnostics*);                                                     \
  template std::vector<std::string> check_outer_scales<W>(                                        \
      const CsrGraph<W>&, const OwnerDecomposition<W>&, const owner::OuterEdgeClassification<W>&, \
      const owner::AuxGraph<W>&, const owner::Calendars<W>&, const OuterStates<W>&,               \
      const OuterScales<W>&);
SSSP_INSTANTIATE_OUTER_SCALES(std::uint64_t)
SSSP_INSTANTIATE_OUTER_SCALES(double)
#undef SSSP_INSTANTIATE_OUTER_SCALES

}  // namespace sssp::paper::outer
