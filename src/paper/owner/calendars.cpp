// Sec 7.3 (arXiv 2611.34911v1, txt l.1692-1779). See calendars.hpp for the statement, the
// scaling convention and the numbering.

#include "paper/owner/calendars.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/local/saturating.hpp"
#include "paper/model/ilog.hpp"
#include "paper/primitives/dsu.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper::owner {

namespace {

constexpr std::uint32_t kNil = std::numeric_limits<std::uint32_t>::max();

template <class W>
bool is_finite(W x) {
  return x != WeightTraits<W>::infinity();
}

/// 16 x by four doublings (saturating, so x <= 16 y stays exact, saturating.hpp).
template <class W>
W times16(W x) {
  return sat::mul_pow2(x, 4);
}

/// Sec 7.3 l.1693-1695: the closed intervals [r_u/4, 4R_u] and [r_v/4, 4R_v] meet.
template <class W>
bool overlaps(W ru, W Ru, W rv, W Rv) {
  return !(times16(Rv) < ru) && !(times16(Ru) < rv);
}

/// ceil(log2 h) for the window widths of Lemma 7.3 (an integer of the word model).
unsigned log_hop(std::size_t h) {
  return model::ceil_log2(std::max<std::uint64_t>(h, 1));
}

template <class W>
void require_aux(const AuxGraph<W>& aux) {
  const std::size_t na = aux.num_vertices();
  if (na == 0 || aux.graph.num_vertices() != na || aux.tree_parent.size() != na ||
      aux.tree_parent_edge.size() != na || aux.tree_order.size() != na || aux.radius.size() != na ||
      aux.owner_distance.size() != na || aux.tree_order.front() != kAuxRoot) {
    throw std::invalid_argument("calendars: inconsistent auxiliary graph");
  }
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex v = aux.tree_order[k];
    if (v >= na || v == kAuxRoot || aux.tree_parent[v] >= na ||
        aux.tree_parent_edge[v] >= aux.graph.num_edges()) {
      throw std::invalid_argument("calendars: T' is not a tree rooted at a");
    }
    if (!(W{0} < aux.radius[v])) {
      throw std::invalid_argument("calendars: zero radius off the root (input not positive)");
    }
  }
}

template <class W>
void require_matching(const AuxGraph<W>& aux, const Calendars<W>& cal) {
  const std::size_t na = aux.num_vertices();
  if (cal.overlap.size() != aux.graph.num_edges() || cal.tree_piece.size() != na ||
      cal.graph_piece.size() != na || cal.width.size() != na || cal.calendar.size() != na ||
      cal.radius_index.size() != na) {
    throw std::invalid_argument("find_atoms: calendars do not belong to the auxiliary graph");
  }
}

/// Adds `value` to the counter `name` (compiled out with SSSP_DIAGNOSTICS=0).
void add_counter(Diagnostics* diag, const char* name, std::uint64_t value) {
  static_cast<void>(name);
  static_cast<void>(value);
  if (diag == nullptr) return;
  SSSP_COUNT(diag->counters[name] += value);
}

/// One marker of the chain under construction: its value is mult * grid[grid] of the piece
/// (an integer position, never a real), linked in increasing order.
template <class W>
struct ChainNode {
  std::uint64_t mult;
  std::uint32_t grid;
  std::uint32_t next;
  W value;
};

/// Whether a marker created at grid index node_grid > j with multiplier mult sits at position t
/// of the width-grid[j] calendar, i.e. mult * 2^(node_grid - j) == t (integer test, no overflow).
bool at_position(std::uint64_t mult, std::uint32_t node_grid, std::uint32_t j, std::uint64_t t) {
  const std::uint32_t d = node_grid - j;
  if (d >= 64) return false;  // the position exceeds every 64-bit t (mult >= 1 off the head)
  const std::uint64_t low = (std::uint64_t{1} << d) - 1;
  return (t & low) == 0 && (t >> d) == mult;
}

}  // namespace

template <class W>
W calendar_scale(W x) {
  // Saturating doublings: an integer result that does not fit becomes infinity(); a double
  // overflows to +inf = infinity().
  const W s = sat::mul_pow2(x, kCalendarScaleLog2);
  if (!is_finite(s)) {
    throw std::overflow_error("calendars: 2^7 x does not fit the weight type");
  }
  return s;
}

template <class W>
Calendars<W> build_calendars(const AuxGraph<W>& aux, Diagnostics* diag) {
  require_aux(aux);
  const std::size_t na = aux.num_vertices();
  const auto edges = aux.graph.edges();
  Calendars<W> cal;

  // Scaled radii and owner distances (Lemma 7.3 l.1723-1724).
  std::vector<W> sr(na, W{0});
  std::vector<W> sR(na, W{0});
  for (Vertex v = 1; v < na; ++v) {
    sr[v] = calendar_scale(aux.radius[v]);
    sR[v] = calendar_scale(aux.owner_distance[v]);
  }

  // 1. Overlap of [r/4, 4R] on the edges of G' - a (l.1693-1695).
  cal.overlap.assign(edges.size(), 0);
  std::uint64_t overlap_edges = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == kAuxRoot || e.v == kAuxRoot || e.u == e.v) continue;
    if (overlaps(aux.radius[e.u], aux.owner_distance[e.u], aux.radius[e.v],
                 aux.owner_distance[e.v])) {
      cal.overlap[i] = 1;
      ++overlap_edges;
    }
  }

  // Tree pieces P: retained parent-child edges of T' - a; the root of P is the vertex whose
  // parent edge is not retained (the vertex of P closest to a, l.1695-1696).
  cal.tree_piece.assign(na, kNoPiece);
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex v = aux.tree_order[k];
    const Vertex p = aux.tree_parent[v];
    if (p != kAuxRoot && cal.overlap[aux.tree_parent_edge[v]] != 0) {
      cal.tree_piece[v] = cal.tree_piece[p];
    } else {
      cal.tree_piece[v] = static_cast<std::uint32_t>(cal.tree_piece_root.size());
      cal.tree_piece_root.push_back(v);
    }
  }

  // Graph pieces C: components of the overlapping edges (union-find, Lemma 3.6).
  DisjointSets dsu(na);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (cal.overlap[i] != 0) dsu.unite(edges[i].u, edges[i].v);
  }
  cal.graph_piece.assign(na, kNoPiece);
  std::vector<std::uint32_t> piece_of_root(na, kNil);
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex v = aux.tree_order[k];
    const auto root = dsu.find(v);
    if (piece_of_root[root] == kNil) {
      piece_of_root[root] = static_cast<std::uint32_t>(cal.graph_pieces.size());
      cal.graph_pieces.emplace_back();
    }
    cal.graph_piece[v] = piece_of_root[root];
  }
  for (const Vertex x : cal.tree_piece_root) {
    cal.tree_piece_graph_piece.push_back(cal.graph_piece[x]);
  }
  const std::size_t num_c = cal.graph_pieces.size();

  // Members of each C in tree order (one stable counting pass, Lemma 3.1).
  std::vector<std::uint32_t> c_offsets(num_c + 1, 0);
  for (Vertex v = 1; v < na; ++v) ++c_offsets[cal.graph_piece[v] + 1];
  for (std::size_t c = 0; c < num_c; ++c) c_offsets[c + 1] += c_offsets[c];
  std::vector<Vertex> c_members(na == 0 ? 0 : na - 1);
  {
    std::vector<std::uint32_t> cursor(c_offsets.begin(), c_offsets.end() - 1);
    for (std::size_t k = 1; k < na; ++k) {
      const Vertex v = aux.tree_order[k];
      c_members[cursor[cal.graph_piece[v]]++] = v;
    }
  }

  // 2. Grids (l.1696-1698): origin = least r, from 7 octaves below the origin (grid[0] = origin
  // in scaled units) to the first point >= 4 max R (scaled: 2^9 max R).
  for (std::size_t c = 0; c < num_c; ++c) {
    auto& gp = cal.graph_pieces[c];
    Vertex root = kNoVertex;
    W max_R{0};
    for (std::uint32_t k = c_offsets[c]; k < c_offsets[c + 1]; ++k) {
      const Vertex v = c_members[k];
      if (root == kNoVertex || aux.radius[v] < aux.radius[root]) root = v;
      max_R = std::max(max_R, aux.owner_distance[v]);
    }
    gp.grid_root = root;
    gp.origin = aux.radius[root];
    const W top = sat::mul_pow2(max_R, kCalendarScaleLog2 + kGridOctavesAbove);
    if (!is_finite(top)) throw std::overflow_error("calendars: 2^9 max R does not fit");
    W g = gp.origin;
    gp.grid.push_back(g);
    while (g < top) {
      g = sat::add(g, g);
      if (!is_finite(g)) throw std::overflow_error("calendars: a grid point does not fit");
      gp.grid.push_back(g);
    }
  }

  // Widths and upper markers (Lemma 7.3 l.1728-1735): BFS over the overlapping edges of each C
  // from its grid root, binary searches in windows around the parent's radius index.
  const unsigned c_h = log_hop(aux.hop_bound);
  const std::uint32_t radius_window = 4 + c_h;  // r_v <= 16 R_u <= 16 h r_u, and symmetrically
  const std::uint32_t upper_window = 1 + c_h;   // R_v <= h r_v
  cal.radius_index.assign(na, kNoPiece);
  cal.width_index.assign(na, kNoPiece);
  cal.upper_index.assign(na, kNoPiece);
  cal.width.assign(na, W{0});
  cal.upper.assign(na, W{0});
  std::uint64_t search_steps = 0;
  const auto place = [&](Vertex v, std::uint32_t lo, std::uint32_t hi) {
    const auto& grid = cal.graph_pieces[cal.graph_piece[v]].grid;
    // Greatest j in [lo, hi] with grid[j] <= 2^7 r_v.
    if (sr[v] < grid[lo]) throw std::logic_error("calendars: radius below its search window");
    while (lo < hi) {
      const std::uint32_t mid = lo + (hi - lo + 1) / 2;
      ++search_steps;
      if (grid[mid] <= sr[v]) {
        lo = mid;
      } else {
        hi = mid - 1;
      }
    }
    const std::uint32_t i = lo;
    if (i + 1 < grid.size() && grid[i + 1] <= sr[v]) {
      throw std::logic_error("calendars: radius above its search window (R > h r?)");
    }
    if (i < kWidthOffset) throw std::logic_error("calendars: radius below the grid");
    // Least k in [i + 1, i + 1 + c_h] with grid[k] > 2^7 R_v.
    std::uint32_t ulo = i + 1;
    std::uint32_t uhi = static_cast<std::uint32_t>(
        std::min<std::size_t>(grid.size() - 1, std::size_t{i} + upper_window));
    if (ulo > uhi || !(sR[v] < grid[uhi])) {
      throw std::logic_error("calendars: upper marker outside its window (R > h r)");
    }
    while (ulo < uhi) {
      const std::uint32_t mid = ulo + (uhi - ulo) / 2;
      ++search_steps;
      if (sR[v] < grid[mid]) {
        uhi = mid;
      } else {
        ulo = mid + 1;
      }
    }
    cal.radius_index[v] = i;
    cal.width_index[v] = i - kWidthOffset;
    cal.upper_index[v] = ulo;
    cal.width[v] = grid[i - kWidthOffset];
    cal.upper[v] = grid[ulo];
  };
  {
    std::vector<std::uint8_t> seen(na, 0);
    std::vector<Vertex> queue;
    for (std::size_t c = 0; c < num_c; ++c) {
      const Vertex root = cal.graph_pieces[c].grid_root;
      const auto top = static_cast<std::uint32_t>(cal.graph_pieces[c].grid.size() - 1);
      place(root, kGridOctavesBelow, kGridOctavesBelow);  // grid[7] = 2^7 origin
      seen[root] = 1;
      queue.assign(1, root);
      for (std::size_t k = 0; k < queue.size(); ++k) {
        const Vertex u = queue[k];
        const std::uint32_t iu = cal.radius_index[u];
        for (const auto& arc : aux.graph.arcs(u)) {
          if (cal.overlap[arc.edge] == 0 || seen[arc.to] != 0) continue;
          seen[arc.to] = 1;
          const std::uint32_t lo = iu > radius_window ? iu - radius_window : 0;
          const auto hi = static_cast<std::uint32_t>(
              std::min<std::size_t>(top, std::size_t{iu} + radius_window));
          place(arc.to, lo, hi);
          queue.push_back(arc.to);
        }
      }
    }
  }

  // 3. Width calendars and merged chains (l.1699-1700, Lemma 7.4).
  cal.calendar.assign(na, kNoCalendar);
  std::vector<ChainNode<W>> nodes;
  std::vector<std::vector<std::uint32_t>> node_markers;  // per calendar of C: chain node ids
  std::vector<std::uint32_t> max_upper;                  // per grid index of C
  std::vector<std::uint32_t> global_id;
  for (std::size_t c = 0; c < num_c; ++c) {
    auto& gp = cal.graph_pieces[c];
    max_upper.assign(gp.grid.size(), kNil);
    for (std::uint32_t k = c_offsets[c]; k < c_offsets[c + 1]; ++k) {
      const Vertex v = c_members[k];
      auto& m = max_upper[cal.width_index[v]];
      m = m == kNil ? cal.upper_index[v] : std::max(m, cal.upper_index[v]);
    }
    std::vector<std::uint32_t> calendar_of_index(gp.grid.size(), kNil);
    nodes.assign(1, ChainNode<W>{0, 0, kNil, W{0}});  // marker 0, shared by every calendar
    node_markers.clear();
    for (std::size_t jj = gp.grid.size(); jj-- > 0;) {  // coarse to fine
      if (max_upper[jj] == kNil) continue;
      const auto j = static_cast<std::uint32_t>(jj);
      const std::uint32_t span = max_upper[j] - j;
      if (span >= 63) throw std::length_error("calendars: too many cells in one calendar");
      const std::uint64_t cells = std::uint64_t{1} << span;  // max U+ / q
      const W q = gp.grid[j];
      calendar_of_index[j] = static_cast<std::uint32_t>(cal.calendars.size());
      gp.calendars.push_back(calendar_of_index[j]);
      cal.calendars.push_back(WidthCalendar<W>{static_cast<std::uint32_t>(c), j, q, {}});
      std::vector<std::uint32_t> mk;
      mk.reserve(static_cast<std::size_t>(cells) + 1);
      std::uint32_t prev = 0;
      mk.push_back(prev);
      std::uint32_t next = nodes[prev].next;
      // Lemma 7.4 l.1748-1750: walk the old chain only up to this calendar's end; every old
      // marker met there is a multiple of q and is reused, the others are inserted.
      for (std::uint64_t t = 1; t <= cells; ++t) {
        std::uint32_t id = 0;
        if (next != kNil && at_position(nodes[next].mult, nodes[next].grid, j, t)) {
          id = next;
          next = nodes[next].next;
        } else {
          if (nodes.size() >= kNil) throw std::length_error("calendars: too many markers");
          id = static_cast<std::uint32_t>(nodes.size());
          nodes.push_back(
              ChainNode<W>{t, j, next, WeightTraits<W>::checked_add(nodes[prev].value, q)});
          nodes[prev].next = id;
        }
        mk.push_back(id);
        prev = id;
      }
      node_markers.push_back(std::move(mk));
    }
    // Number the chain in increasing order.
    gp.first_marker = static_cast<std::uint32_t>(cal.marker_value.size());
    global_id.assign(nodes.size(), kNil);
    for (std::uint32_t id = 0; id != kNil; id = nodes[id].next) {
      if (cal.marker_value.size() >= kNil) throw std::length_error("calendars: too many markers");
      global_id[id] = static_cast<std::uint32_t>(cal.marker_value.size());
      cal.marker_value.push_back(nodes[id].value);
      cal.marker_graph_piece.push_back(static_cast<std::uint32_t>(c));
    }
    gp.end_marker = static_cast<std::uint32_t>(cal.marker_value.size());
    for (std::size_t k = 0; k < gp.calendars.size(); ++k) {
      auto& markers = cal.calendars[gp.calendars[k]].markers;
      markers.reserve(node_markers[k].size());
      for (const std::uint32_t id : node_markers[k]) markers.push_back(global_id[id]);
    }
    for (std::uint32_t k = c_offsets[c]; k < c_offsets[c + 1]; ++k) {
      const Vertex v = c_members[k];
      cal.calendar[v] = calendar_of_index[cal.width_index[v]];
    }
  }

  // 4. Interpiece marker arcs U+_u -> 16 q_v (l.1701-1703, Lemma 7.6 l.1766-1769).
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == kAuxRoot || e.v == kAuxRoot || e.u == e.v) continue;
    if (cal.graph_piece[e.u] == cal.graph_piece[e.v]) continue;
    Vertex tail = e.u;
    Vertex head = e.v;
    if (!(times16(aux.owner_distance[tail]) < aux.radius[head])) std::swap(tail, head);
    if (!(times16(aux.owner_distance[tail]) < aux.radius[head])) {
      throw std::logic_error("calendars: an edge between graph pieces overlaps");
    }
    const auto& from_cal = cal.calendars[cal.calendar[tail]];
    const auto& to_cal = cal.calendars[cal.calendar[head]];
    const std::uint64_t upper_pos = std::uint64_t{1}
                                    << (cal.upper_index[tail] - cal.width_index[tail]);
    cal.interpiece_arcs.push_back(MarkerArc{from_cal.markers[upper_pos], to_cal.markers[16],
                                            static_cast<EdgeId>(i), tail, head});
  }

  // The marker DAG and a topological order (l.1703-1704, Lemma 3.1).
  {
    std::vector<DiArc> arcs;
    arcs.reserve(cal.marker_value.size() + cal.interpiece_arcs.size());
    for (const auto& gp : cal.graph_pieces) {
      for (std::uint32_t m = gp.first_marker; m + 1 < gp.end_marker; ++m)
        arcs.push_back({m, m + 1});
    }
    for (const auto& a : cal.interpiece_arcs) arcs.push_back({a.from, a.to});
    const auto dag = Digraph::from_arcs(cal.marker_value.size(), arcs);
    const auto order = topological_sort(dag);
    cal.marker_order.assign(order.begin(), order.end());
  }

  add_counter(diag, "calendars.overlap_edges", overlap_edges);
  add_counter(diag, "calendars.tree_pieces", cal.tree_piece_root.size());
  add_counter(diag, "calendars.graph_pieces", num_c);
  add_counter(diag, "calendars.grid_points", cal.num_grid_points());
  add_counter(diag, "calendars.grid_search_steps", search_steps);
  add_counter(diag, "calendars.calendars", cal.calendars.size());
  add_counter(diag, "calendars.cells", cal.num_cells());
  add_counter(diag, "calendars.markers", cal.marker_value.size());
  add_counter(diag, "calendars.interpiece_arcs", cal.interpiece_arcs.size());
  return cal;
}

template <class W>
Atoms<W> find_atoms(const AuxGraph<W>& aux, const Calendars<W>& cal, Diagnostics* diag) {
  require_aux(aux);
  require_matching(aux, cal);
  const std::size_t na = aux.num_vertices();
  const auto edges = aux.graph.edges();

  // Defining edges w' < 2 max(q_u, q_v) on nonpivots (l.1705-1706, Lemma 7.5), compared as
  // 2^7 w' < 2 max(q) in the scaled units of the calendars (saturating: exact).
  DisjointSets dsu(na);
  std::vector<EdgeId> joins;
  std::uint64_t defining = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == kAuxRoot || e.v == kAuxRoot || e.u == e.v) continue;
    const W q = std::max(cal.width[e.u], cal.width[e.v]);
    if (!(sat::mul_pow2(e.weight, kCalendarScaleLog2) < sat::add(q, q))) continue;
    ++defining;
    // Lemma 7.5: such an edge cannot change the root bottleneck of the MST T'.
    if (aux.radius[e.u] != aux.radius[e.v]) {
      throw std::logic_error("find_atoms: a defining edge joins different radii (T' not an MST)");
    }
    if (dsu.unite(e.u, e.v)) joins.push_back(static_cast<EdgeId>(i));
  }

  Atoms<W> atoms;
  atoms.atom.assign(na, kNoAtom);
  std::vector<std::uint32_t> atom_of_root(na, kNil);
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex v = aux.tree_order[k];
    const auto root = dsu.find(v);
    if (atom_of_root[root] == kNil) {
      atom_of_root[root] = static_cast<std::uint32_t>(atoms.representative.size());
      atoms.representative.push_back(v);
    }
    atoms.atom[v] = atom_of_root[root];
  }
  const std::size_t num_atoms = atoms.representative.size();

  // Members in tree order and discovery-tree edges in edge order (stable counting passes).
  atoms.offsets.assign(num_atoms + 1, 0);
  for (Vertex v = 1; v < na; ++v) ++atoms.offsets[atoms.atom[v] + 1];
  for (std::size_t a = 0; a < num_atoms; ++a) atoms.offsets[a + 1] += atoms.offsets[a];
  atoms.vertices.assign(na - 1, kNoVertex);
  {
    std::vector<std::uint32_t> cursor(atoms.offsets.begin(), atoms.offsets.end() - 1);
    for (std::size_t k = 1; k < na; ++k) {
      const Vertex v = aux.tree_order[k];
      atoms.vertices[cursor[atoms.atom[v]]++] = v;
    }
  }
  atoms.tree_offsets.assign(num_atoms + 1, 0);
  for (const EdgeId i : joins) ++atoms.tree_offsets[atoms.atom[edges[i].u] + 1];
  for (std::size_t a = 0; a < num_atoms; ++a) atoms.tree_offsets[a + 1] += atoms.tree_offsets[a];
  atoms.tree_edges.assign(joins.size(), kNoEdge);
  {
    std::vector<std::uint32_t> cursor(atoms.tree_offsets.begin(), atoms.tree_offsets.end() - 1);
    for (const EdgeId i : joins) atoms.tree_edges[cursor[atoms.atom[edges[i].u]]++] = i;
  }
  atoms.length.assign(num_atoms, W{0});
  atoms.calendar.assign(num_atoms, kNoCalendar);
  std::uint64_t singletons = 0;
  std::uint64_t max_size = 0;
  for (std::size_t a = 0; a < num_atoms; ++a) {
    W len{0};
    for (std::uint32_t k = atoms.tree_offsets[a]; k < atoms.tree_offsets[a + 1]; ++k) {
      len = WeightTraits<W>::checked_add(len, edges[atoms.tree_edges[k]].weight);
    }
    atoms.length[a] = len;
    atoms.calendar[a] = cal.calendar[atoms.representative[a]];
    const std::uint64_t size = atoms.offsets[a + 1] - atoms.offsets[a];
    if (size == 1) ++singletons;
    max_size = std::max(max_size, size);
  }

  add_counter(diag, "atoms.atoms", num_atoms);
  add_counter(diag, "atoms.singletons", singletons);
  add_counter(diag, "atoms.defining_edges", defining);
  if (diag != nullptr) {
    SSSP_COUNT(diag->counters["atoms.max_size"] =
                   std::max(diag->counters["atoms.max_size"], max_size));
  }
  return atoms;
}

// ---------------------------------------------------------------------------------- checkers

namespace {

/// x == a + b exactly (integers), or within the recursive-summation bound (2M + 2) eps |x| on a
/// chain of M markers (double, see calendars.hpp).
template <class W>
bool sum_matches(W x, W a, W b, std::size_t chain) {
  if constexpr (std::is_integral_v<W>) {
    return a <= WeightTraits<W>::infinity() - b && x == static_cast<W>(a + b);
  } else {
    const W s = a + b;
    const W tol = static_cast<W>(2 * chain + 2) * std::numeric_limits<W>::epsilon() *
                  std::max(std::abs(x), std::abs(s));
    return std::abs(x - s) <= tol;
  }
}

/// x == y exactly (integers), or within (2M + 2) eps |x| (double).
template <class W>
bool value_matches(W x, W y, std::size_t chain) {
  return sum_matches(x, y, W{0}, chain);
}

}  // namespace

template <class W>
std::vector<std::string> check_calendars(const AuxGraph<W>& aux, const Calendars<W>& cal) {
  std::vector<std::string> out;
  const auto fail = [&out](const std::string& what) { out.push_back(what); };
  const std::size_t na = aux.num_vertices();
  const auto edges = aux.graph.edges();
  if (na == 0 || aux.radius.size() != na || aux.owner_distance.size() != na ||
      aux.tree_parent.size() != na || aux.tree_order.size() != na ||
      cal.overlap.size() != edges.size() || cal.tree_piece.size() != na ||
      cal.graph_piece.size() != na || cal.radius_index.size() != na ||
      cal.width_index.size() != na || cal.upper_index.size() != na || cal.width.size() != na ||
      cal.upper.size() != na || cal.calendar.size() != na ||
      cal.tree_piece_graph_piece.size() != cal.tree_piece_root.size() ||
      cal.marker_graph_piece.size() != cal.marker_value.size()) {
    fail("sizes do not match");
    return out;
  }
  const std::size_t num_p = cal.tree_piece_root.size();
  const std::size_t num_c = cal.graph_pieces.size();
  const std::size_t num_m = cal.marker_value.size();
  const auto r = [&](Vertex v) { return aux.radius[v]; };
  const auto R = [&](Vertex v) { return aux.owner_distance[v]; };

  // Overlap flags (l.1693-1695).
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    const bool off_root = e.u != kAuxRoot && e.v != kAuxRoot && e.u != e.v;
    const bool expect =
        off_root && !(sat::mul_pow2(R(e.v), 4) < r(e.u)) && !(sat::mul_pow2(R(e.u), 4) < r(e.v));
    if ((cal.overlap[i] != 0) != expect) fail("aux edge " + std::to_string(i) + ": overlap flag");
  }
  if (cal.tree_piece[kAuxRoot] != kNoPiece || cal.graph_piece[kAuxRoot] != kNoPiece ||
      cal.calendar[kAuxRoot] != kNoCalendar) {
    fail("root a carries a piece or calendar");
  }

  // Tree pieces: v is the root of its piece iff its parent is a or its parent edge is cut;
  // otherwise it shares its parent's piece. Lemma 7.3: R_x = r_x below a, 16 R_x <= 17 r_x at
  // other piece roots, R_v <= 17 R_u on retained edges.
  for (Vertex v = 1; v < na; ++v) {
    const std::string id = "vertex " + std::to_string(v) + ": ";
    const std::uint32_t P = cal.tree_piece[v];
    const std::uint32_t C = cal.graph_piece[v];
    if (P >= num_p || C >= num_c) {
      fail(id + "piece index out of range");
      return out;
    }
    if (cal.tree_piece_graph_piece[P] != C) fail(id + "tree piece not inside its graph piece");
    const Vertex p = aux.tree_parent[v];
    if (p >= na) {
      fail(id + "no parent in T'");
      return out;
    }
    const bool retained = p != kAuxRoot && cal.overlap[aux.tree_parent_edge[v]] != 0;
    if (retained) {
      if (cal.tree_piece[p] != P) fail(id + "retained tree edge between two tree pieces");
      if (cal.tree_piece_root[P] == v) fail(id + "piece root with a retained parent edge");
      if (sat::mul_int(R(p), 17) < R(v)) fail(id + "R_v > 17 R_u on a retained edge (Lemma 7.3)");
    } else {
      if (cal.tree_piece_root[P] != v) fail(id + "cut below its parent but not a piece root");
      if (p == kAuxRoot) {
        if (R(v) != r(v)) fail(id + "child of a with R != r (Lemma 7.3)");
      } else if (sat::mul_int(r(v), 17) < sat::mul_pow2(R(v), 4)) {
        fail(id + "piece root with 16 R > 17 r (Lemma 7.3)");
      }
    }
  }
  for (std::size_t P = 0; P < num_p; ++P) {
    const Vertex x = cal.tree_piece_root[P];
    if (x == kAuxRoot || x >= na || cal.tree_piece[x] != P) fail("tree piece root mismatch");
  }

  // Graph pieces: exactly the components of the overlapping edges.
  {
    DisjointSets dsu(na);
    for (std::size_t i = 0; i < edges.size(); ++i) {
      if (cal.overlap[i] != 0) dsu.unite(edges[i].u, edges[i].v);
    }
    std::vector<std::uint32_t> piece_of_root(na, kNil);
    std::vector<std::uint32_t> root_of_piece(num_c, kNil);
    for (Vertex v = 1; v < na; ++v) {
      const auto root = dsu.find(v);
      const std::uint32_t C = cal.graph_piece[v];
      if (piece_of_root[root] == kNil) piece_of_root[root] = C;
      if (root_of_piece[C] == kNil) root_of_piece[C] = root;
      if (piece_of_root[root] != C || root_of_piece[C] != root) {
        fail("vertex " + std::to_string(v) + ": graph piece is not an overlap component");
      }
    }
    for (std::size_t C = 0; C < num_c; ++C) {
      if (root_of_piece[C] == kNil) fail("graph piece " + std::to_string(C) + " is empty");
    }
  }

  // Grids (l.1696-1698).
  std::vector<W> min_r(num_c, WeightTraits<W>::infinity());
  std::vector<W> max_R(num_c, W{0});
  for (Vertex v = 1; v < na; ++v) {
    min_r[cal.graph_piece[v]] = std::min(min_r[cal.graph_piece[v]], r(v));
    max_R[cal.graph_piece[v]] = std::max(max_R[cal.graph_piece[v]], R(v));
  }
  for (std::size_t C = 0; C < num_c; ++C) {
    const auto& gp = cal.graph_pieces[C];
    const std::string id = "graph piece " + std::to_string(C) + ": ";
    if (gp.origin != min_r[C]) fail(id + "origin is not the least radius");
    if (gp.grid_root >= na || cal.graph_piece[gp.grid_root] != C || r(gp.grid_root) != gp.origin) {
      fail(id + "grid root is not a vertex of least radius");
    }
    if (gp.grid.size() <= kGridOctavesBelow || gp.grid.front() != gp.origin) {
      fail(id + "grid does not start 7 octaves below the origin");
      continue;
    }
    for (std::size_t j = 0; j + 1 < gp.grid.size(); ++j) {
      if (gp.grid[j + 1] != sat::add(gp.grid[j], gp.grid[j])) fail(id + "grid is not dyadic");
    }
    if (gp.grid[kGridOctavesBelow] != sat::mul_pow2(gp.origin, kCalendarScaleLog2)) {
      fail(id + "grid[7] is not the scaled origin");
    }
    if (gp.grid.back() < sat::mul_pow2(max_R[C], kCalendarScaleLog2 + kGridOctavesAbove)) {
      fail(id + "grid ends below 4 max R");
    }
  }

  // Widths and upper markers: 32 q <= r < 64 q and R < U+ <= 2R on grid values.
  for (Vertex v = 1; v < na; ++v) {
    const std::string id = "vertex " + std::to_string(v) + ": ";
    const auto& grid = cal.graph_pieces[cal.graph_piece[v]].grid;
    const std::uint32_t j = cal.width_index[v];
    const std::uint32_t k = cal.upper_index[v];
    if (j >= grid.size() || k >= grid.size() || cal.radius_index[v] != j + kWidthOffset ||
        cal.width[v] != grid[j] || cal.upper[v] != grid[k]) {
      fail(id + "width or upper marker is not the indexed grid value");
      continue;
    }
    const W sr = sat::mul_pow2(r(v), kCalendarScaleLog2);
    const W sR = sat::mul_pow2(R(v), kCalendarScaleLog2);
    const W q = cal.width[v];
    if (sr < sat::mul_pow2(q, 5) || !(sr < sat::mul_pow2(q, 6))) fail(id + "not 32 q <= r < 64 q");
    if (!(sR < cal.upper[v]) || sat::add(sR, sR) < cal.upper[v]) fail(id + "not R < U+ <= 2R");
  }

  // Calendars: one per occurring (C, q), coarse to fine within C, markers 0, q, ..., max U+.
  std::vector<std::vector<std::uint32_t>> expected_upper(num_c);
  for (std::size_t C = 0; C < num_c; ++C) {
    expected_upper[C].assign(cal.graph_pieces[C].grid.size(), kNil);
  }
  for (Vertex v = 1; v < na; ++v) {
    auto& row = expected_upper[cal.graph_piece[v]];
    if (cal.width_index[v] >= row.size()) continue;
    auto& m = row[cal.width_index[v]];
    m = m == kNil ? cal.upper_index[v] : std::max(m, cal.upper_index[v]);
  }
  std::vector<std::uint8_t> marker_used(num_m, 0);
  std::size_t calendars_expected = 0;
  for (std::size_t C = 0; C < num_c; ++C) {
    const auto& gp = cal.graph_pieces[C];
    const std::string id = "graph piece " + std::to_string(C) + ": ";
    if (gp.first_marker > gp.end_marker || gp.end_marker > num_m) {
      fail(id + "marker range out of bounds");
      continue;
    }
    const std::size_t chain = gp.end_marker - gp.first_marker;
    for (std::uint32_t m = gp.first_marker; m < gp.end_marker; ++m) {
      if (cal.marker_graph_piece[m] != C) fail(id + "marker of another graph piece in the chain");
      if (m + 1 < gp.end_marker && !(cal.marker_value[m] < cal.marker_value[m + 1])) {
        fail(id + "chain not strictly increasing at marker " + std::to_string(m));
      }
    }
    std::vector<std::uint32_t> listed;
    for (std::size_t jj = gp.grid.size(); jj-- > 0;) {
      if (expected_upper[C][jj] != kNil) listed.push_back(static_cast<std::uint32_t>(jj));
    }
    calendars_expected += listed.size();
    if (gp.calendars.size() != listed.size()) {
      fail(id + "calendars do not match the occurring widths");
      continue;
    }
    for (std::size_t k = 0; k < listed.size(); ++k) {
      const std::uint32_t cid = gp.calendars[k];
      if (cid >= cal.calendars.size()) {
        fail(id + "calendar id out of range");
        continue;
      }
      const auto& wc = cal.calendars[cid];
      const std::uint32_t j = listed[k];
      const std::string cname = "calendar " + std::to_string(cid) + ": ";
      if (wc.graph_piece != C || wc.grid_index != j || wc.width != gp.grid[j]) {
        fail(cname + "wrong graph piece or width (or not coarse to fine)");
        continue;
      }
      const std::uint32_t span = expected_upper[C][j] - j;
      if (span >= 63 || wc.markers.size() != (std::size_t{1} << span) + 1) {
        fail(cname + "cell count is not max U+ / q");
        continue;
      }
      for (std::size_t i = 0; i < wc.markers.size(); ++i) {
        const std::uint32_t m = wc.markers[i];
        if (m < gp.first_marker || m >= gp.end_marker || (i > 0 && m <= wc.markers[i - 1])) {
          fail(cname + "markers not increasing ids of the chain");
          break;
        }
        marker_used[m] = 1;
        if (i == 0 ? cal.marker_value[m] != W{0}
                   : !sum_matches(cal.marker_value[m], cal.marker_value[wc.markers[i - 1]],
                                  wc.width, chain)) {
          fail(cname + "marker " + std::to_string(i) + " is not i q");
        }
      }
      if (!value_matches(cal.marker_value[wc.markers.back()], gp.grid[expected_upper[C][j]],
                         chain)) {
        fail(cname + "calendar does not end at max U+");
      }
    }
  }
  if (calendars_expected != cal.calendars.size()) fail("calendars without a graph piece");
  for (std::size_t m = 0; m < num_m; ++m) {
    if (marker_used[m] == 0) fail("marker " + std::to_string(m) + " belongs to no calendar");
  }
  for (Vertex v = 1; v < na; ++v) {
    const std::uint32_t c = cal.calendar[v];
    if (c >= cal.calendars.size() || cal.calendars[c].graph_piece != cal.graph_piece[v] ||
        cal.calendars[c].grid_index != cal.width_index[v]) {
      fail("vertex " + std::to_string(v) + ": calendar is not (C_v, q_v)");
    }
  }
  if (!out.empty()) return out;

  // Interpiece arcs (l.1701-1703, Lemma 7.6): one per edge between graph pieces, in edge order.
  std::size_t next_arc = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto& e = edges[i];
    if (e.u == kAuxRoot || e.v == kAuxRoot || e.u == e.v) continue;
    if (cal.graph_piece[e.u] == cal.graph_piece[e.v]) continue;
    const std::string id = "aux edge " + std::to_string(i) + ": ";
    if (next_arc >= cal.interpiece_arcs.size() || cal.interpiece_arcs[next_arc].aux_edge != i) {
      fail(id + "edge between graph pieces without its marker arc");
      continue;
    }
    const auto& a = cal.interpiece_arcs[next_arc++];
    if (!((a.tail == e.u && a.head == e.v) || (a.tail == e.v && a.head == e.u))) {
      fail(id + "arc endpoints are not the edge's");
      continue;
    }
    if (!(sat::mul_pow2(R(a.tail), 4) < r(a.head))) fail(id + "not r_head > 16 R_tail");
    const auto& from_cal = cal.calendars[cal.calendar[a.tail]];
    const auto& to_cal = cal.calendars[cal.calendar[a.head]];
    const std::uint32_t span = cal.upper_index[a.tail] - cal.width_index[a.tail];
    if (span >= 63 || (std::size_t{1} << span) >= from_cal.markers.size() ||
        a.from != from_cal.markers[std::size_t{1} << span] || to_cal.markers.size() <= 16 ||
        a.to != to_cal.markers[16]) {
      fail(id + "arc is not U+_tail -> 16 q_head");
      continue;
    }
    const std::size_t chain = num_m;
    if (!value_matches(cal.marker_value[a.from], cal.upper[a.tail], chain) ||
        !value_matches(cal.marker_value[a.to], sat::mul_pow2(cal.width[a.head], 4), chain)) {
      fail(id + "arc marker values are not U+_tail and 16 q_head");
    }
    if (!(cal.marker_value[a.from] < cal.marker_value[a.to])) fail(id + "arc does not increase");
  }
  if (next_arc != cal.interpiece_arcs.size()) fail("marker arcs without an interpiece edge");

  // marker_order: a permutation along which every chain and interpiece arc points forward.
  if (cal.marker_order.size() != num_m) {
    fail("marker order is not a permutation");
    return out;
  }
  std::vector<std::size_t> pos(num_m, num_m);
  for (std::size_t k = 0; k < num_m; ++k) {
    const std::uint32_t m = cal.marker_order[k];
    if (m >= num_m || pos[m] != num_m) {
      fail("marker order is not a permutation");
      return out;
    }
    pos[m] = k;
  }
  for (const auto& gp : cal.graph_pieces) {
    for (std::uint32_t m = gp.first_marker; m + 1 < gp.end_marker; ++m) {
      if (pos[m + 1] < pos[m]) fail("marker order violates chain arc " + std::to_string(m));
    }
  }
  for (const auto& a : cal.interpiece_arcs) {
    if (a.from < num_m && a.to < num_m && pos[a.to] < pos[a.from]) {
      fail("marker order violates an interpiece arc");
    }
  }
  return out;
}

template <class W>
std::vector<std::string> check_atoms(const AuxGraph<W>& aux, const Calendars<W>& cal,
                                     const Atoms<W>& atoms) {
  std::vector<std::string> out;
  const auto fail = [&out](const std::string& what) { out.push_back(what); };
  const std::size_t na = aux.num_vertices();
  const auto edges = aux.graph.edges();
  const std::size_t num_a = atoms.num_atoms();
  if (na == 0 || cal.width.size() != na || cal.tree_piece.size() != na ||
      cal.graph_piece.size() != na || cal.calendar.size() != na || atoms.atom.size() != na ||
      atoms.offsets.size() != num_a + 1 || atoms.tree_offsets.size() != num_a + 1 ||
      atoms.length.size() != num_a || atoms.calendar.size() != num_a ||
      atoms.vertices.size() != na - 1 || atoms.offsets.back() != na - 1 ||
      atoms.tree_offsets.back() != atoms.tree_edges.size()) {
    fail("sizes do not match");
    return out;
  }
  if (atoms.atom[kAuxRoot] != kNoAtom) fail("the root a is in an atom");
  // b = h - 1 (h = b + 1, saturated).
  const std::size_t b = aux.hop_bound == 0 ? 0 : aux.hop_bound - 1;
  const auto defining = [&](std::size_t i) {
    const auto& e = edges[i];
    if (e.u == kAuxRoot || e.v == kAuxRoot || e.u == e.v) return false;
    const W q = std::max(cal.width[e.u], cal.width[e.v]);
    return sat::mul_pow2(e.weight, kCalendarScaleLog2) < sat::add(q, q);
  };

  // Membership lists against the per-vertex atom ids.
  std::vector<std::uint8_t> listed(na, 0);
  for (std::size_t a = 0; a < num_a; ++a) {
    if (atoms.offsets[a] >= atoms.offsets[a + 1]) {
      fail("atom " + std::to_string(a) + " is empty");
      return out;
    }
    for (const Vertex v : atoms.members(a)) {
      if (v == kAuxRoot || v >= na || listed[v] != 0 || atoms.atom[v] != a) {
        fail("atom " + std::to_string(a) + ": member list disagrees with the atom ids");
        return out;
      }
      listed[v] = 1;
    }
  }

  // Closed under defining edges.
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (defining(i) && atoms.atom[edges[i].u] != atoms.atom[edges[i].v]) {
      fail("aux edge " + std::to_string(i) + ": defining edge between two atoms");
    }
  }

  for (std::size_t a = 0; a < num_a; ++a) {
    const std::string id = "atom " + std::to_string(a) + ": ";
    const auto members = atoms.members(a);
    const Vertex rep = atoms.representative[a];
    const Vertex first = members.front();
    if (members.size() >= b) fail(id + "|Y| >= b (Lemma 7.5)");
    if (atoms.calendar[a] != cal.calendar[first]) fail(id + "calendar is not the members'");
    std::size_t tops = 0;
    for (const Vertex v : members) {
      if (cal.tree_piece[v] != cal.tree_piece[first] ||
          cal.graph_piece[v] != cal.graph_piece[first] ||
          cal.width_index[v] != cal.width_index[first] || aux.radius[v] != aux.radius[first]) {
        fail(id + "members differ in piece, width or radius (Lemma 7.5)");
      }
      const Vertex p = aux.tree_parent[v];
      if (p >= na || atoms.atom[p] != a) {
        ++tops;
        if (v != rep) fail(id + "representative is not the top of T'[Y]");
      }
    }
    if (tops != 1) fail(id + "T'[Y] is not connected");
    // Discovery tree: |Y| - 1 defining edges inside Y that connect Y.
    const std::size_t t_size = atoms.tree_offsets[a + 1] - atoms.tree_offsets[a];
    if (t_size + 1 != members.size()) {
      fail(id + "discovery tree does not have |Y| - 1 edges");
      continue;
    }
    DisjointSets dsu(na);
    W len{0};
    bool ok = true;
    for (std::uint32_t k = atoms.tree_offsets[a]; k < atoms.tree_offsets[a + 1]; ++k) {
      const EdgeId i = atoms.tree_edges[k];
      if (i >= edges.size() || !defining(i) || atoms.atom[edges[i].u] != a ||
          !dsu.unite(edges[i].u, edges[i].v)) {
        ok = false;
        break;
      }
      len = sat::add(len, edges[i].weight);
    }
    if (!ok) {
      fail(id + "discovery tree is not a forest of defining edges inside Y");
      continue;
    }
    if (len != atoms.length[a]) fail(id + "L_Y is not the discovery tree's weight");
    // L_Y <= 2 q |Y| (l.1708-1709), scaled.
    const W q = cal.width[first];
    if (sat::mul_int(sat::add(q, q), members.size()) <
        sat::mul_pow2(atoms.length[a], kCalendarScaleLog2)) {
      fail(id + "L_Y > 2 q |Y|");
    }
    // Sec 7.4 l.1793-1794: R_v <= R_rep + L_Y. For double, R_v (summed down T') and L_Y
    // (summed in edge order) round differently: allow (|Y| + 1) eps (R_rep + L_Y).
    if (rep < na) {
      W horizon = sat::add(aux.owner_distance[rep], atoms.length[a]);
      if constexpr (std::is_floating_point_v<W>) {
        horizon += static_cast<W>(members.size() + 1) * std::numeric_limits<W>::epsilon() * horizon;
      }
      for (const Vertex v : members) {
        if (horizon < aux.owner_distance[v]) fail(id + "R_v > R_rep + L_Y");
      }
    }
  }
  return out;
}

#define SSSP_CALENDARS_INSTANTIATE(W)                                                            \
  template W calendar_scale<W>(W);                                                               \
  template Calendars<W> build_calendars<W>(const AuxGraph<W>&, Diagnostics*);                    \
  template Atoms<W> find_atoms<W>(const AuxGraph<W>&, const Calendars<W>&, Diagnostics*);        \
  template std::vector<std::string> check_calendars<W>(const AuxGraph<W>&, const Calendars<W>&); \
  template std::vector<std::string> check_atoms<W>(const AuxGraph<W>&, const Calendars<W>&,      \
                                                   const Atoms<W>&);

SSSP_CALENDARS_INSTANTIATE(std::uint64_t)
SSSP_CALENDARS_INSTANTIATE(double)

#undef SSSP_CALENDARS_INSTANTIATE

}  // namespace sssp::paper::owner
