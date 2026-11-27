#pragma once
// Sec 7.3 "Width Classes, Calendars, and Atoms" (arXiv 2611.34911v1, txt l.1692-1779; Lemma 7.3
// l.1711, Lemma 7.4 l.1736, Lemma 7.5 l.1751, Lemma 7.6 l.1762), i.e. Alg 7 l.2-3, Mode A (Mode B
// reaches it through other separator sizes b; nothing here depends on b except the atom bound).
//
// Input: the auxiliary graph G' and its owner-spoke MST T' rooted at the common pivot a
// (AuxGraph, Sec 7.2), where r_v and R_v are the root bottleneck and the root distance in T'
// (l.1693) and h = b + 1 is the hop bound (R_v <= h r_v, Lemma 7.2).
//
// build_calendars (Alg 7 l.2, Lemma 7.6):
//   1. Overlap (l.1693-1695): an edge uv of G' - a overlaps iff the closed intervals [r_u/4, 4R_u]
//      and [r_v/4, 4R_v] meet, i.e. r_u <= 16 R_v and r_v <= 16 R_u (four doublings).
//      Tree pieces P: components of the overlapping parent-child edges of T' - a, rooted at the
//      vertex closest to a. Graph pieces C: components of all overlapping edges of G' - a
//      (union-find, Lemma 3.6). Every P lies in one C.
//   2. Grid of C (l.1696-1698): origin = least r in C; the dyadic grid runs from 7 octaves below
//      the origin (Q-46) up to the first point >= 4 max R (two octaves above the largest R).
//      Widths and upper markers (l.1698-1699, Lemma 7.3 l.1728-1735): from the parent in a BFS
//      tree of the overlapping edges of C rooted at a vertex of least radius, binary search the
//      grid index i_v of r_v (grid[i_v] <= r_v < grid[i_v + 1]) among the indices within
//      4 + ceil(log2 h) of the parent's index (r_v <= 16 R_u <= 16 h r_u and symmetrically), and
//      the index of U+_v among the 1 + ceil(log2 h) indices above i_v (R_v <= h r_v). Then
//      q_v = grid[i_v - 5] (32 q <= r < 64 q) and U+_v = the least grid value > R_v
//      (R < U+ <= 2R).
//   3. Width calendars (l.1699-1700, Lemma 7.4): one per occurring pair (C, q), with cells
//      [i q, (i + 1) q) from 0 up to the largest U+ of the pair, i.e. markers 0, q, ..., max U+.
//      The calendars of C are merged from coarse to fine into one increasing marker chain; a
//      finer calendar walks the chain only up to its own end and reuses every marker it meets
//      there (each such marker is a multiple of the finer width; l.1748-1750). Cell values are
//      formed by additions only.
//   4. Marker arcs (l.1701-1704, Lemma 7.6 l.1766-1775): for an edge of G' - a between two graph
//      pieces, label the ends with r_v > 16 R_u and add U+_u -> 16 q_v; together with the
//      consecutive-marker arcs of each chain they form a DAG whose arcs increase the value
//      (U+_u <= 2 R_u < r_v / 8 < 16 q_v). marker_order is a topological order (Lemma 3.1).
//
// find_atoms (Alg 7 l.3, l.1705-1709, Lemma 7.5): the components of the defining edges of
// G' - a, w'(uv) < 2 max(q_u, q_v); the edges that join two components in increasing edge order
// form the discovery tree, L_Y is its total auxiliary weight. Every nonroot vertex is in exactly
// one atom (singletons included).
//
// Model of computation (Sec 2): lengths are only compared, added, copied and doubled; grid
// indices, cell counts and marker positions are integers. Scaling (Lemma 7.3 l.1723-1724, "after
// scaling all quantities by a fixed power of two"): the grid starts 7 octaves below the origin,
// so every calendar value (grid points, q, U+, marker values) is stored as 2^7 times the length
// it stands for. With this scale grid[j] = origin * 2^j is formed by j doublings of the origin
// itself, grid[7] = 2^7 origin, and every marker is an integer multiple of a grid value, so the
// integer instantiation is exact. Compare an auxiliary length x with them through
// calendar_scale(x) = 2^7 x. For double the doublings are exact and the cell sums round (Q-06).
// Atom lengths L_Y are not scaled.
//
// Numbering (interpretation; the paper fixes none): tree pieces, graph pieces and atoms are
// numbered in the BFS order of T' (aux.tree_order) of their first vertex; a tree piece's root and
// an atom's representative are their first vertex in that order (the vertex closest to a). The
// calendars of a graph piece are numbered coarse to fine, graph piece by graph piece; the markers
// of a chain are consecutive ids in increasing value, graph piece by graph piece.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "paper/owner/aux_metric.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::owner {

inline constexpr std::uint32_t kNoPiece = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoCalendar = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoMarker = std::numeric_limits<std::uint32_t>::max();
inline constexpr std::uint32_t kNoAtom = std::numeric_limits<std::uint32_t>::max();

/// Sec 7.3 l.1698: the grid extends at least seven octaves below the origin (Q-46: taken as
/// exactly seven) and at least two above the largest R.
inline constexpr unsigned kGridOctavesBelow = 7;
inline constexpr unsigned kGridOctavesAbove = 2;
/// Lemma 7.3 l.1734: q_v is the grid value five indices below the radius index.
inline constexpr unsigned kWidthOffset = 5;
/// Calendar values are 2^kCalendarScaleLog2 times the lengths they stand for (see above).
inline constexpr unsigned kCalendarScaleLog2 = kGridOctavesBelow;

/// 2^7 x by seven doublings. Throws std::overflow_error if the result does not fit W (integers:
/// it would reach infinity(); double: it is not finite).
template <class W>
W calendar_scale(W x);

/// One graph piece C (Sec 7.3 l.1696-1698).
template <class W>
struct GraphPiece {
  W origin{};                            ///< least r in C (input units)
  Vertex grid_root = kNoVertex;          ///< first vertex (tree order) of radius `origin`
  std::vector<W> grid;                   ///< scaled: grid[j] = origin 2^j; grid[7] = 2^7 origin
  std::vector<std::uint32_t> calendars;  ///< calendar ids of C, coarse to fine
  std::uint32_t first_marker = 0;        ///< the chain of C is the marker ids [first, end)
  std::uint32_t end_marker = 0;
};

/// The width calendar of one pair (C, q) (l.1699-1700).
template <class W>
struct WidthCalendar {
  std::uint32_t graph_piece = kNoPiece;
  std::uint32_t grid_index = 0;  ///< q = grid[grid_index] of the graph piece
  W width{};                     ///< q (scaled)
  /// markers[i] = marker id of i q, i = 0 .. max U+ / q; cell i is [markers[i], markers[i+1]).
  std::vector<std::uint32_t> markers;

  [[nodiscard]] std::size_t num_cells() const noexcept {
    return markers.empty() ? 0 : markers.size() - 1;
  }
};

/// Interpiece marker arc U+_tail -> 16 q_head for the edge `aux_edge` of G' (l.1701-1703).
struct MarkerArc {
  std::uint32_t from = kNoMarker;
  std::uint32_t to = kNoMarker;
  EdgeId aux_edge = kNoEdge;
  Vertex tail = kNoVertex;  ///< the endpoint of smaller range: r_head > 16 R_tail
  Vertex head = kNoVertex;
};

/// build_calendars output. Aux-vertex-indexed vectors have aux.num_vertices() entries; the root
/// a has kNoPiece / kNoCalendar indices and zero widths.
template <class W>
struct Calendars {
  std::vector<std::uint8_t> overlap;        ///< per aux edge: 1 iff in G' - a and overlapping
  std::vector<std::uint32_t> tree_piece;    ///< per vertex: tree piece P
  std::vector<std::uint32_t> graph_piece;   ///< per vertex: graph piece C
  std::vector<std::uint32_t> radius_index;  ///< per vertex: i_v, grid[i_v] <= 2^7 r_v < 2 grid[i_v]
  std::vector<std::uint32_t> width_index;   ///< per vertex: i_v - 5
  std::vector<std::uint32_t> upper_index;   ///< per vertex: grid index of U+_v
  std::vector<W> width;                     ///< per vertex: q_v (scaled)
  std::vector<W> upper;                     ///< per vertex: U+_v (scaled)
  std::vector<std::uint32_t> calendar;      ///< per vertex: calendar of (C_v, q_v)

  std::vector<Vertex> tree_piece_root;                ///< per tree piece
  std::vector<std::uint32_t> tree_piece_graph_piece;  ///< per tree piece
  std::vector<GraphPiece<W>> graph_pieces;
  std::vector<WidthCalendar<W>> calendars;

  std::vector<W> marker_value;                    ///< per marker (scaled)
  std::vector<std::uint32_t> marker_graph_piece;  ///< per marker
  std::vector<MarkerArc> interpiece_arcs;         ///< increasing aux edge order
  std::vector<std::uint32_t> marker_order;        ///< a topological order of the marker DAG

  [[nodiscard]] std::size_t num_cells() const noexcept {
    std::size_t s = 0;
    for (const auto& c : calendars) s += c.num_cells();
    return s;
  }
  [[nodiscard]] std::size_t num_grid_points() const noexcept {
    std::size_t s = 0;
    for (const auto& c : graph_pieces) s += c.grid.size();
    return s;
  }
};

/// Alg 7 l.2 (Lemma 7.6). Throws std::invalid_argument if `aux` is inconsistent (sizes, T' not
/// rooted at a) or a nonroot vertex has r_v = 0 (the input is positive, l.1642),
/// std::overflow_error if 2^7 R or a grid point up to 2^10 max R does not fit W, and
/// std::logic_error if the binary searches leave their windows (the input violates R <= h r of
/// Lemma 7.2). Counters "calendars.*": overlap_edges, tree_pieces, graph_pieces, grid_points,
/// grid_search_steps (comparisons in the window searches), calendars, cells, markers,
/// interpiece_arcs.
template <class W>
Calendars<W> build_calendars(const AuxGraph<W>& aux, Diagnostics* diag = nullptr);

/// find_atoms output. Atom-indexed vectors have num_atoms() entries.
template <class W>
struct Atoms {
  std::vector<std::uint32_t> atom;          ///< per aux vertex (kNoAtom at the root)
  std::vector<std::uint32_t> offsets;       ///< members of atom i: vertices[offsets[i] ..]
  std::vector<Vertex> vertices;             ///< members, in T' BFS order within each atom
  std::vector<Vertex> representative;       ///< first member in T' BFS order
  std::vector<std::uint32_t> tree_offsets;  ///< discovery tree of atom i: tree_edges[...]
  std::vector<EdgeId> tree_edges;           ///< aux edges, increasing within each atom
  std::vector<W> length;                    ///< L_Y = total w' of the discovery tree
  std::vector<std::uint32_t> calendar;      ///< the common (C, q) calendar

  [[nodiscard]] std::size_t num_atoms() const noexcept { return representative.size(); }
  [[nodiscard]] std::span<const Vertex> members(std::size_t i) const {
    return {vertices.data() + offsets[i], vertices.data() + offsets[i + 1]};
  }
};

/// Alg 7 l.3 (Lemma 7.5). Throws std::invalid_argument if `cal` does not belong to `aux`, and
/// std::logic_error if a defining edge joins different radii (T' would not be an MST of G',
/// Lemma 7.2). Counters "atoms.*": atoms, singletons, defining_edges, max_size.
template <class W>
Atoms<W> find_atoms(const AuxGraph<W>& aux, const Calendars<W>& cal, Diagnostics* diag = nullptr);

/// Checks build_calendars' output against the definitions of Sec 7.3 and Lemmas 7.3, 7.4, 7.6:
/// overlap flags, tree pieces (rooted at the top, one per cut or root edge) and graph pieces
/// (components of the overlapping edges), the Lemma 7.3 root bounds R_x = r_x below a and
/// 16 R_x <= 17 r_x at other piece roots, R_v <= 17 R_u on retained tree edges, the grids
/// (origin, doubling, 7 octaves below, >= 4 max R on top), 32 q <= r < 64 q and R < U+ <= 2R,
/// one calendar per occurring (C, q) with markers 0, q, ..., max U+, increasing chains made of
/// calendar markers, one interpiece arc per edge between graph pieces with r_head > 16 R_tail
/// from U+_tail to 16 q_head, increasing arc values and a topological marker_order. Integer
/// values are compared exactly; for double, a marker is compared with the sum it stands for up
/// to the recursive-summation bound (2M + 2) eps |x| on a chain of M markers. Returns one
/// message per violation.
template <class W>
std::vector<std::string> check_calendars(const AuxGraph<W>& aux, const Calendars<W>& cal);

/// Checks find_atoms' output against Sec 7.3 l.1705-1709, Lemma 7.5 and Sec 7.4 l.1786-1794:
/// atoms are the components of the defining edges, each with < b = h - 1 vertices, one tree
/// piece, one graph piece, one width and one radius, a discovery tree of |Y| - 1 defining edges
/// spanning Y with L_Y its total weight and L_Y <= 2 q |Y|, T'[Y] connected with the
/// representative on top, and R_v <= R_rep + L_Y for every member (for double up to the
/// rounding bound (|Y| + 1) eps (R_rep + L_Y): R_v and L_Y are sums in different orders, Q-06).
/// Returns one message per violation.
template <class W>
std::vector<std::string> check_atoms(const AuxGraph<W>& aux, const Calendars<W>& cal,
                                     const Atoms<W>& atoms);

}  // namespace sssp::paper::owner
