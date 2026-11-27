#pragma once
// Text graph formats.
//
// 1. Project edge-list format ("*.txt"):
//      n m
//      u v w        (m lines, 0-based vertices, undirected edge {u,v})
//    Lines whose first non-blank character is '#' and blank lines are ignored anywhere.
//    Exactly m edge lines are required. Self-loops and parallel edges are kept as given.
//
// 2. 9th DIMACS Implementation Challenge shortest-path format ("*.gr"):
//      c <comment>
//      p sp n m     (exactly once, before the first arc; m = number of ARC lines)
//      a u v w      (m lines, 1-based vertices, directed arc u->v)
//    The graph is read as UNDIRECTED. DimacsArcMode chooses how arcs become edges, then
//    DuplicateEdgePolicy chooses what to do with remaining parallel edges (see below).
//
// Numbers are parsed and printed with std::from_chars / std::to_chars, so the result does not
// depend on the C locale. Integer weights (std::uint64_t) must be plain unsigned decimals
// (0 .. 2^64-1, no sign, no fraction, no exponent). Double weights accept any finite
// non-negative decimal or exponent form; -0 is normalised to +0; they are written in %.17g form
// (std::chars_format::general, precision 17), which round-trips every finite double bit-exactly.
//
// Every malformed input throws ParseError with the 1-based line number (0 = end of input or the
// input as a whole): negative / NaN / infinite / overflowing weights, out-of-range vertices
// (including DIMACS vertex 0), missing / duplicate / malformed header, too few or too many edge
// lines, non-numeric tokens, wrong token counts, n >= kNoVertex, m >= kNoEdge.

#include <cstddef>
#include <iosfwd>
#include <stdexcept>
#include <string>

#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"

namespace sssp::io {

class ParseError : public std::runtime_error {
 public:
  ParseError(std::size_t line, const std::string& reason);
  /// 1-based line of the offending input line; 0 = end of input / whole input.
  [[nodiscard]] std::size_t line() const noexcept { return line_; }

 private:
  std::size_t line_;
};

/// How DIMACS arcs become undirected edges.
enum class DimacsArcMode {
  /// Every arc line becomes one undirected edge. A symmetric pair u->v, v->u therefore yields two
  /// parallel edges; shortest-path distances are unaffected, but m doubles.
  kArcAsEdge,
  /// Arcs u->v and v->u with bit-identical weight are matched greedily in file order; each matched
  /// pair becomes ONE edge, emitted at the position and in the orientation of the first arc of the
  /// pair. Unmatched arcs (one-way, or asymmetric weights) become one edge each; a self-loop arc
  /// u->u always becomes one edge. This is the right mode for the DIMACS road networks, which
  /// list every road as two arcs, and it inverts write_dimacs exactly.
  kMergeReversePairs,
};

/// What to do with parallel edges (same unordered endpoint pair, including two self-loops on the
/// same vertex) that remain after DimacsArcMode has been applied.
enum class DuplicateEdgePolicy {
  kKeepAll,        ///< keep every edge (multigraph)
  kKeepMinWeight,  ///< keep one edge at its first position/orientation, with the minimum weight
  kReject,         ///< throw ParseError at the line of the first duplicate
};

struct DimacsReadOptions {
  DimacsArcMode arc_mode = DimacsArcMode::kArcAsEdge;
  DuplicateEdgePolicy duplicates = DuplicateEdgePolicy::kKeepAll;
};

// Instantiated for W = std::uint64_t and W = double.

template <class W>
EdgeList<W> read_edge_list(std::istream& in);
template <class W>
void write_edge_list(std::ostream& out, const EdgeList<W>& g);

template <class W>
EdgeList<W> read_dimacs(std::istream& in, const DimacsReadOptions& opts = {});
/// Writes every undirected edge as two arcs (u->v then v->u), a self-loop as one arc.
template <class W>
void write_dimacs(std::ostream& out, const EdgeList<W>& g);

/// File variants; throw std::runtime_error if the file cannot be opened or written.
template <class W>
EdgeList<W> read_edge_list_file(const std::string& path);
template <class W>
void write_edge_list_file(const std::string& path, const EdgeList<W>& g);
template <class W>
EdgeList<W> read_dimacs_file(const std::string& path, const DimacsReadOptions& opts = {});
template <class W>
void write_dimacs_file(const std::string& path, const EdgeList<W>& g);

}  // namespace sssp::io
