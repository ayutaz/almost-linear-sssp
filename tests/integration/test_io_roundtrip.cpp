// Generators -> text formats -> parsers must reproduce every graph bit for bit.
#include <bit>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "sssp/generators.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/graph_io.hpp"

using sssp::EdgeList;

namespace {

template <class W>
bool bit_identical(const EdgeList<W>& a, const EdgeList<W>& b) {
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

/// One instance of every generator type with a weight distribution exercising the format.
template <class W>
std::vector<sssp::gen::NamedGraph<W>> corpus() {
  using sssp::gen::WeightDist;
  std::vector<sssp::gen::NamedGraph<W>> out;
  const WeightDist dists[] = {WeightDist::kUniformReal, WeightDist::kWideRange,
                              WeightDist::kExponential, WeightDist::kZeroMixed};
  std::size_t i = 0;
  for (const auto& info : sssp::gen::generator_catalog()) {
    sssp::gen::GraphSpec spec;
    spec.type = std::string(info.name);
    spec.n = 300;
    spec.options.seed = 100 + i;
    spec.options.weights.dist = dists[i % 4];
    spec.options.weights.hi = 1000.0;
    ++i;
    out.push_back({spec.type, sssp::gen::generate<W>(spec), 0, info.description.data()});
  }
  // configuration model with loops and multi-edges kept: parallel arcs of equal weight
  out.push_back(
      {"regular_keep",
       sssp::gen::random_regular<W>(
           100, 6, sssp::gen::MultiEdgePolicy::kKeep,
           {.weights = {.dist = WeightDist::kConstant, .lo = 3}, .seed = 1, .connected = false}),
       0, "multi-edges and self-loops"});
  for (auto& g : sssp::gen::small_graph_suite<W>()) out.push_back(std::move(g));
  return out;
}

std::string temp_path(const std::string& tag) {
  return (std::filesystem::temp_directory_path() / ("sssp_io_roundtrip_" + tag + ".txt")).string();
}

}  // namespace

TEMPLATE_TEST_CASE("edge-list format round-trips every generated graph bit for bit",
                   "[integration][io]", std::uint64_t, double) {
  for (const auto& g : corpus<TestType>()) {
    INFO(g.name);
    std::stringstream ss;
    sssp::io::write_edge_list(ss, g.graph);
    const auto back = sssp::io::read_edge_list<TestType>(ss);
    REQUIRE(bit_identical(back, g.graph));
  }
}

TEMPLATE_TEST_CASE("DIMACS format round-trips every generated graph with reverse-pair merging",
                   "[integration][io]", std::uint64_t, double) {
  sssp::io::DimacsReadOptions merge{.arc_mode = sssp::io::DimacsArcMode::kMergeReversePairs};
  for (const auto& g : corpus<TestType>()) {
    INFO(g.name);
    std::stringstream ss;
    sssp::io::write_dimacs(ss, g.graph);
    const auto back = sssp::io::read_dimacs<TestType>(ss, merge);
    REQUIRE(bit_identical(back, g.graph));
    // without merging every non-loop edge appears twice
    std::stringstream again;
    sssp::io::write_dimacs(again, g.graph);
    const auto doubled = sssp::io::read_dimacs<TestType>(again);
    std::size_t loops = 0;
    for (const auto& e : g.graph.edges) loops += (e.u == e.v);
    REQUIRE(doubled.edges.size() == 2 * g.graph.edges.size() - loops);
  }
}

TEMPLATE_TEST_CASE("file round-trip feeds an identical CSR graph", "[integration][io]",
                   std::uint64_t, double) {
  sssp::gen::GraphSpec spec;
  spec.type = "road";
  spec.rows = 20;
  spec.cols = 30;
  spec.options.seed = 4;
  const auto g = sssp::gen::generate<TestType>(spec);
  const std::string tag = std::is_integral_v<TestType> ? "u64" : "f64";
  const auto p1 = temp_path("el_" + tag);
  const auto p2 = temp_path("gr_" + tag);
  sssp::io::write_edge_list_file(p1, g);
  sssp::io::write_dimacs_file(p2, g);
  const auto a = sssp::to_csr(sssp::io::read_edge_list_file<TestType>(p1));
  const auto b = sssp::to_csr(sssp::io::read_dimacs_file<TestType>(
      p2, {.arc_mode = sssp::io::DimacsArcMode::kMergeReversePairs}));
  const auto c = sssp::to_csr(g);
  std::filesystem::remove(p1);
  std::filesystem::remove(p2);
  for (const auto* x : {&a, &b}) {
    REQUIRE(x->num_vertices() == c.num_vertices());
    REQUIRE(x->num_arcs() == c.num_arcs());
    for (sssp::Vertex v = 0; v < c.num_vertices(); ++v) {
      const auto xa = x->arcs(v);
      const auto ca = c.arcs(v);
      REQUIRE(xa.size() == ca.size());
      for (std::size_t i = 0; i < xa.size(); ++i) {
        REQUIRE(xa[i].to == ca[i].to);
        REQUIRE(xa[i].weight == ca[i].weight);
        REQUIRE(xa[i].edge == ca[i].edge);
      }
    }
  }
}
