#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "sssp/graph_formats.hpp"
#include "sssp/graph_io.hpp"

using sssp::Edge;
using sssp::EdgeList;
using sssp::io::DimacsArcMode;
using sssp::io::DimacsReadOptions;
using sssp::io::DuplicateEdgePolicy;
using sssp::io::ParseError;

namespace {

template <class W>
EdgeList<W> read_el(const std::string& text) {
  std::istringstream in(text);
  return sssp::io::read_edge_list<W>(in);
}

template <class W>
EdgeList<W> read_gr(const std::string& text, DimacsReadOptions opts = {}) {
  std::istringstream in(text);
  return sssp::io::read_dimacs<W>(in, opts);
}

/// Returns the line number of the ParseError thrown by f, or -1 if none was thrown.
template <class F>
long parse_error_line(F&& f) {
  try {
    f();
  } catch (const ParseError& e) {
    return static_cast<long>(e.line());
  }
  return -1;
}

template <class W>
void require_edge(const Edge<W>& e, sssp::Vertex u, sssp::Vertex v, W w) {
  REQUIRE(e.u == u);
  REQUIRE(e.v == v);
  REQUIRE(e.weight == w);
}

template <class W>
bool same_bits(const EdgeList<W>& a, const EdgeList<W>& b) {
  if (a.num_vertices != b.num_vertices || a.edges.size() != b.edges.size()) return false;
  for (std::size_t i = 0; i < a.edges.size(); ++i) {
    const auto& x = a.edges[i];
    const auto& y = b.edges[i];
    if (x.u != y.u || x.v != y.v) return false;
    if constexpr (std::is_floating_point_v<W>) {
      if (std::bit_cast<std::uint64_t>(x.weight) != std::bit_cast<std::uint64_t>(y.weight)) {
        return false;
      }
    } else {
      if (x.weight != y.weight) return false;
    }
  }
  return true;
}

}  // namespace

// ------------------------------------------------------------------ edge-list format

TEST_CASE("edge-list format reads n, m and 0-based edges", "[graph][io][edgelist]") {
  const auto g = read_el<std::uint64_t>("3 2\n0 1 5\n1 2 7\n");
  REQUIRE(g.num_vertices == 3);
  REQUIRE(g.edges.size() == 2);
  require_edge<std::uint64_t>(g.edges[0], 0, 1, 5);
  require_edge<std::uint64_t>(g.edges[1], 1, 2, 7);
}

TEST_CASE("edge-list format accepts comments, blank lines and extra spaces",
          "[graph][io][edgelist]") {
  const auto g = read_el<double>("# header comment\n\n  3   2 \n# mid\n0 1 0.5\n\n2\t1  2.25\n");
  REQUIRE(g.num_vertices == 3);
  REQUIRE(g.edges.size() == 2);
  require_edge<double>(g.edges[0], 0, 1, 0.5);
  require_edge<double>(g.edges[1], 2, 1, 2.25);
}

TEST_CASE("edge-list format keeps self-loops, multi-edges and isolated vertices",
          "[graph][io][edgelist]") {
  const auto g = read_el<std::uint64_t>("5 3\n0 1 1\n0 1 2\n2 2 3\n");
  REQUIRE(g.num_vertices == 5);
  REQUIRE(g.edges.size() == 3);
  require_edge<std::uint64_t>(g.edges[2], 2, 2, 3);
}

TEST_CASE("edge-list format reads empty graphs", "[graph][io][edgelist]") {
  REQUIRE(read_el<double>("0 0\n").num_vertices == 0);
  const auto g = read_el<double>("4 0\n");
  REQUIRE(g.num_vertices == 4);
  REQUIRE(g.edges.empty());
}

TEST_CASE("edge-list writer output is the documented text", "[graph][io][edgelist]") {
  const EdgeList<std::uint64_t> g{3, {{0, 1, 5}, {2, 1, 7}}};
  std::ostringstream out;
  sssp::io::write_edge_list(out, g);
  REQUIRE(out.str() == "3 2\n0 1 5\n2 1 7\n");
}

TEST_CASE("edge-list round-trips uint64 weights including 2^64-1", "[graph][io][edgelist]") {
  const auto max = std::numeric_limits<std::uint64_t>::max();
  const EdgeList<std::uint64_t> g{4, {{0, 1, 0}, {1, 2, max}, {2, 3, 1ULL << 50}, {3, 3, 1}}};
  std::ostringstream out;
  sssp::io::write_edge_list(out, g);
  REQUIRE(same_bits(read_el<std::uint64_t>(out.str()), g));
}

TEST_CASE("edge-list round-trips double weights bit-exactly", "[graph][io][edgelist]") {
  const double two53 = 9007199254740992.0;
  const std::vector<double> ws{0.0,
                               0.1,
                               1.0 / 3.0,
                               1e-300,
                               5e-324,  // smallest subnormal
                               std::numeric_limits<double>::min(),
                               two53 + 2.0,
                               std::nextafter(1.0, 2.0),
                               std::numeric_limits<double>::max(),
                               123456789.123456789};
  EdgeList<double> g{2, {}};
  for (double w : ws) g.edges.push_back({0, 1, w});
  std::ostringstream out;
  sssp::io::write_edge_list(out, g);
  REQUIRE(same_bits(read_el<double>(out.str()), g));
}

TEST_CASE("double weights are written in %.17g form", "[graph][io][edgelist]") {
  const EdgeList<double> g{2, {{0, 1, 0.1}, {1, 0, 2.0}}};
  std::ostringstream out;
  sssp::io::write_edge_list(out, g);
  REQUIRE(out.str() == "2 2\n0 1 0.10000000000000001\n1 0 2\n");
}

TEST_CASE("negative zero is accepted and normalised to +0", "[graph][io][edgelist]") {
  const auto g = read_el<double>("2 1\n0 1 -0.0\n");
  REQUIRE(g.edges[0].weight == 0.0);
  REQUIRE_FALSE(std::signbit(g.edges[0].weight));
}

TEST_CASE("double type reads integers, exponents and decimals", "[graph][io][edgelist]") {
  const auto g = read_el<double>("2 3\n0 1 7\n0 1 1.5e3\n0 1 .25\n");
  REQUIRE(g.edges[0].weight == 7.0);
  REQUIRE(g.edges[1].weight == 1500.0);
  REQUIRE(g.edges[2].weight == 0.25);
}

TEST_CASE("edge-list format rejects invalid input with ParseError and line number",
          "[graph][io][edgelist][validation]") {
  using U = std::uint64_t;
  // negative weight
  REQUIRE(parse_error_line([] { read_el<U>("2 1\n0 1 -3\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<double>("2 1\n0 1 -0.5\n"); }) == 2);
  // NaN / infinity
  REQUIRE(parse_error_line([] { read_el<double>("2 1\n0 1 nan\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<double>("2 1\n0 1 inf\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<double>("2 1\n0 1 1e400\n"); }) == 2);
  // out-of-range vertex
  REQUIRE(parse_error_line([] { read_el<U>("2 2\n0 1 1\n# c\n1 2 1\n"); }) == 4);
  // too few edge lines (reported at EOF = 0)
  REQUIRE(parse_error_line([] { read_el<U>("3 2\n0 1 1\n"); }) == 0);
  // too many edge lines
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 1\n1 2 1\n"); }) == 3);
  // missing header / empty input
  REQUIRE(parse_error_line([] { read_el<U>(""); }) == 0);
  REQUIRE(parse_error_line([] { read_el<U>("# only a comment\n"); }) == 0);
  // malformed header / edge lines
  REQUIRE(parse_error_line([] { read_el<U>("3\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_el<U>("3 x\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_el<U>("3 1 9\n0 1 1\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 2 3\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 a 2\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 2x\n"); }) == 2);
  // integer type rejects fractional weights and overflow
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 1.5\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 18446744073709551616\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_el<U>("3 1\n0 1 +4\n"); }) == 2);
  // too many vertices for the 32-bit Vertex type
  REQUIRE(parse_error_line([] { read_el<U>("4294967295 0\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_el<U>("4294967296 0\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_el<U>("3 4294967295\n"); }) == 1);
}

TEST_CASE("ParseError message names the line and the problem", "[graph][io][validation]") {
  try {
    read_el<std::uint64_t>("2 1\n0 1 -3\n");
    FAIL("expected ParseError");
  } catch (const ParseError& e) {
    REQUIRE_THAT(e.what(), Catch::Matchers::ContainsSubstring("line 2"));
    REQUIRE_THAT(e.what(), Catch::Matchers::ContainsSubstring("negative"));
  }
  try {
    read_el<std::uint64_t>("2 1\n0 5 3\n");
    FAIL("expected ParseError");
  } catch (const ParseError& e) {
    REQUIRE_THAT(e.what(), Catch::Matchers::ContainsSubstring("out of range"));
  }
}

// ------------------------------------------------------------------ DIMACS format

namespace {
// Two undirected edges {1,2} w=3 and {2,3} w=4 written as symmetric arc pairs, plus comments.
const char* const kSymmetricGr =
    "c 9th DIMACS challenge style\n"
    "c\n"
    "p sp 3 4\n"
    "a 1 2 3\n"
    "a 2 1 3\n"
    "\n"
    "a 2 3 4\n"
    "a 3 2 4\n";
}  // namespace

TEST_CASE("DIMACS reader converts 1-based arcs to 0-based edges", "[graph][io][dimacs]") {
  const auto g = read_gr<std::uint64_t>(kSymmetricGr);
  REQUIRE(g.num_vertices == 3);
  REQUIRE(g.edges.size() == 4);  // kArcAsEdge (default): one edge per arc
  require_edge<std::uint64_t>(g.edges[0], 0, 1, 3);
  require_edge<std::uint64_t>(g.edges[1], 1, 0, 3);
  require_edge<std::uint64_t>(g.edges[2], 1, 2, 4);
  require_edge<std::uint64_t>(g.edges[3], 2, 1, 4);
}

TEST_CASE("DIMACS kMergeReversePairs merges symmetric arc pairs", "[graph][io][dimacs]") {
  const auto g = read_gr<std::uint64_t>(kSymmetricGr, {DimacsArcMode::kMergeReversePairs});
  REQUIRE(g.num_vertices == 3);
  REQUIRE(g.edges.size() == 2);
  require_edge<std::uint64_t>(g.edges[0], 0, 1, 3);
  require_edge<std::uint64_t>(g.edges[1], 1, 2, 4);
}

TEST_CASE("DIMACS kMergeReversePairs keeps asymmetric and one-way arcs", "[graph][io][dimacs]") {
  const std::string text =
      "p sp 4 6\n"
      "a 1 2 5\n"  // asymmetric pair: 5 vs 6 -> two edges
      "a 2 1 6\n"
      "a 3 4 1\n"  // one-way arc -> one edge
      "a 2 3 2\n"  // pair split by other arcs, reversed order -> one edge at first occurrence
      "a 4 4 9\n"  // self-loop arc -> one edge
      "a 3 2 2\n";
  const auto g = read_gr<std::uint64_t>(text, {DimacsArcMode::kMergeReversePairs});
  REQUIRE(g.edges.size() == 5);
  require_edge<std::uint64_t>(g.edges[0], 0, 1, 5);
  require_edge<std::uint64_t>(g.edges[1], 1, 0, 6);
  require_edge<std::uint64_t>(g.edges[2], 2, 3, 1);
  require_edge<std::uint64_t>(g.edges[3], 1, 2, 2);
  require_edge<std::uint64_t>(g.edges[4], 3, 3, 9);
}

TEST_CASE("DIMACS kMergeReversePairs pairs repeated arcs one-to-one", "[graph][io][dimacs]") {
  // three 1->2 arcs and one 2->1 arc (same weight): one pair + two unmatched arcs = 3 edges
  const std::string text = "p sp 2 4\na 1 2 1\na 1 2 1\na 2 1 1\na 1 2 1\n";
  const auto g = read_gr<double>(text, {DimacsArcMode::kMergeReversePairs});
  REQUIRE(g.edges.size() == 3);
}

TEST_CASE("DIMACS duplicate policies", "[graph][io][dimacs]") {
  // after merging: {1,2} w=5, {2,1} w=6 (parallel, asymmetric), {1,2} w=2, {3,3} twice
  const std::string text =
      "p sp 3 5\n"
      "a 1 2 5\n"
      "a 2 1 6\n"
      "a 2 1 2\n"
      "a 3 3 4\n"
      "a 3 3 1\n";
  SECTION("kKeepAll keeps every parallel edge") {
    const auto g = read_gr<std::uint64_t>(
        text, {DimacsArcMode::kMergeReversePairs, DuplicateEdgePolicy::kKeepAll});
    REQUIRE(g.edges.size() == 5);
  }
  SECTION("kKeepMinWeight keeps one edge per pair at its first position, with min weight") {
    const auto g = read_gr<std::uint64_t>(
        text, {DimacsArcMode::kMergeReversePairs, DuplicateEdgePolicy::kKeepMinWeight});
    REQUIRE(g.edges.size() == 2);
    require_edge<std::uint64_t>(g.edges[0], 0, 1, 2);
    require_edge<std::uint64_t>(g.edges[1], 2, 2, 1);
  }
  SECTION("kReject throws at the line of the first duplicate") {
    REQUIRE(parse_error_line([&] {
              read_gr<std::uint64_t>(
                  text, {DimacsArcMode::kMergeReversePairs, DuplicateEdgePolicy::kReject});
            }) == 3);
  }
  SECTION("kReject accepts a symmetric pair once merged") {
    REQUIRE(read_gr<std::uint64_t>(
                kSymmetricGr, {DimacsArcMode::kMergeReversePairs, DuplicateEdgePolicy::kReject})
                .edges.size() == 2);
  }
}

TEST_CASE("DIMACS writer emits both arcs per edge and one arc per self-loop",
          "[graph][io][dimacs]") {
  const EdgeList<std::uint64_t> g{3, {{0, 1, 3}, {2, 2, 4}}};
  std::ostringstream out;
  sssp::io::write_dimacs(out, g);
  const std::string s = out.str();
  REQUIRE_THAT(s, Catch::Matchers::ContainsSubstring("p sp 3 3\n"));
  REQUIRE_THAT(s, Catch::Matchers::ContainsSubstring("a 1 2 3\na 2 1 3\na 3 3 4\n"));
}

TEST_CASE("DIMACS write then read with kMergeReversePairs restores the edge list",
          "[graph][io][dimacs]") {
  const EdgeList<double> g{5,
                           {{0, 1, 0.1},
                            {1, 0, 0.1},  // parallel edge in the opposite direction
                            {2, 2, 1e-300},
                            {3, 1, 5e-324},
                            {1, 3, 7.0},
                            {0, 1, std::numeric_limits<double>::max()}}};
  std::ostringstream out;
  sssp::io::write_dimacs(out, g);
  std::istringstream in(out.str());
  const auto back = sssp::io::read_dimacs<double>(in, {DimacsArcMode::kMergeReversePairs});
  REQUIRE(same_bits(back, g));
}

TEST_CASE("DIMACS reader rejects invalid input with ParseError and line number",
          "[graph][io][dimacs][validation]") {
  using U = std::uint64_t;
  // vertex 0 is invalid (1-based)
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\na 0 1 1\n"); }) == 2);
  // vertex beyond n
  REQUIRE(parse_error_line([] { read_gr<U>("c x\np sp 2 1\na 1 3 1\n"); }) == 3);
  // negative / NaN weights
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\na 1 2 -1\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_gr<double>("p sp 2 1\na 1 2 nan\n"); }) == 2);
  // arc before the problem line
  REQUIRE(parse_error_line([] { read_gr<U>("a 1 2 1\np sp 2 1\n"); }) == 1);
  // duplicate problem line
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\np sp 2 1\na 1 2 1\n"); }) == 2);
  // missing problem line
  REQUIRE(parse_error_line([] { read_gr<U>("c only comments\n"); }) == 0);
  // wrong problem type / malformed problem line
  REQUIRE(parse_error_line([] { read_gr<U>("p max 2 1\na 1 2 1\n"); }) == 1);
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2\n"); }) == 1);
  // too few / too many arcs
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 2\na 1 2 1\n"); }) == 0);
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\na 1 2 1\na 2 1 1\n"); }) == 3);
  // unknown line type, malformed arc
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\nx 1 2 1\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\na 1 2\n"); }) == 2);
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 2 1\na 1 2 1.5\n"); }) == 2);
  // too many vertices
  REQUIRE(parse_error_line([] { read_gr<U>("p sp 4294967295 0\n"); }) == 1);
}

TEST_CASE("file variants round-trip and report unopenable files", "[graph][io][file]") {
  const std::string path = "sssp_test_graph_io_tmp.txt";
  const EdgeList<std::uint64_t> g{3, {{0, 1, 5}, {1, 2, 7}}};
  sssp::io::write_edge_list_file(path, g);
  REQUIRE(same_bits(sssp::io::read_edge_list_file<std::uint64_t>(path), g));
  sssp::io::write_dimacs_file(path, g);
  REQUIRE(same_bits(
      sssp::io::read_dimacs_file<std::uint64_t>(path, {DimacsArcMode::kMergeReversePairs}), g));
  std::remove(path.c_str());
  REQUIRE_THROWS_AS(sssp::io::read_edge_list_file<double>("/nonexistent/dir/graph.txt"),
                    std::runtime_error);
  REQUIRE_THROWS_AS(sssp::io::read_dimacs_file<double>("/nonexistent/dir/graph.gr"),
                    std::runtime_error);
  REQUIRE_THROWS_AS(sssp::io::write_edge_list_file("/nonexistent/dir/graph.txt", g),
                    std::runtime_error);
}
