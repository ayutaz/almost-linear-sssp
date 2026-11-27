#pragma once
// Result and diagnostics types shared by all SSSP implementations.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "sssp/graph.hpp"

#ifndef SSSP_DIAGNOSTICS
#define SSSP_DIAGNOSTICS 1
#endif

#if SSSP_DIAGNOSTICS
/// Increment a diagnostic counter; compiled out when SSSP_DIAGNOSTICS=0 (timed runs).
#define SSSP_COUNT(expr) \
  do {                   \
    expr;                \
  } while (0)
#else
#define SSSP_COUNT(expr) \
  do {                   \
  } while (0)
#endif

namespace sssp {

/// Operation counts in the comparison-addition model plus memory-traffic proxies.
struct OpCounters {
  std::uint64_t vertex_scans = 0;  ///< vertices whose adjacency was scanned
  std::uint64_t edge_scans = 0;    ///< arcs inspected
  std::uint64_t relaxations = 0;   ///< successful distance improvements
  std::uint64_t comparisons = 0;   ///< weight/distance comparisons
  std::uint64_t additions = 0;     ///< weight/distance additions
  std::uint64_t queue_pushes = 0;
  std::uint64_t queue_pops = 0;
  std::uint64_t decrease_keys = 0;
};

/// A region of the input that an implementation handled with a direct search instead of the
/// paper's construction. Every fallback MUST be recorded; see CLAUDE.md "Fallback policy".
struct FallbackRecord {
  std::string reason;  ///< e.g. "outer: N <= L_N^64 (paper base case)"
  std::uint64_t vertices = 0;
  std::uint64_t edges = 0;
};

struct Diagnostics {
  OpCounters ops;
  /// Named execution-path counters, e.g. "outer_constructions", "local_hierarchies",
  /// "samplings", "slot_evaluations", "state_graph_vertices".
  std::map<std::string, std::uint64_t> counters;
  std::vector<FallbackRecord> fallbacks;
};

/// Single-source result: distances, a shortest-path tree (parent vertex + parent edge id), and
/// diagnostics. Unreachable vertices have dist = WeightTraits<W>::infinity() and kNoVertex parent.
template <class W>
struct SsspResult {
  std::vector<W> dist;
  std::vector<Vertex> parent;
  std::vector<EdgeId> parent_edge;
  Diagnostics diag;
};

}  // namespace sssp
