// sssp_run: runs one baseline SSSP algorithm on a graph file and prints the result as JSON.
//
//   sssp_run --graph g.txt --algorithm dijkstra_binary --source 0 [--verify]
//
// The result is one line of sssp::write_result_json
//   {"n":..,"source":..,"dist":[..],"parent":[..],"parent_edge":[..]}
// (doubles as %.17g, unreachable as null, no parent as -1), readable by
// tests/python/reference/check_double.py. Exit status: 0 ok, 1 runtime error (unreadable or
// malformed graph, invalid source, integer overflow), 2 usage error, 3 --verify rejected the
// result.

#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

#include "sssp/baseline.hpp"
#include "sssp/graph.hpp"
#include "sssp/graph_formats.hpp"
#include "sssp/graph_io.hpp"
#include "sssp/verify.hpp"

namespace {

constexpr std::string_view kUsage =
    R"(usage: sssp_run --graph PATH --algorithm ALG --source S [options]
  --graph PATH                graph file
  --algorithm ALG             bellman_ford | dijkstra_binary | dijkstra_4ary |
                              dijkstra_fibonacci | dijkstra_radix (uint64 only)
  --source S                  source vertex (0-based)
  --format auto|edgelist|dimacs
                              input format (default auto: *.gr is DIMACS, anything else edge list)
  --weight-type uint64|double weight type (default uint64)
  --dimacs-arcs merge|edge    DIMACS: merge reverse arc pairs into one edge (default; inverts
                              sssp_gen --format dimacs) or read every arc as an edge
  --duplicates keep|min|reject
                              DIMACS: parallel edges after arc conversion (default keep)
  --out PATH                  write the result JSON to PATH instead of stdout
  --verify                    check the result with sssp::verify_sssp (shortest-path-tree
                              certificate, exact) and report "verify: OK" / failures on stderr
exit status: 0 ok, 1 runtime error, 2 usage error, 3 verification failed
)";

struct UsageError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

enum class Algorithm : std::uint8_t {
  kBellmanFord,
  kDijkstraBinary,
  kDijkstra4ary,
  kDijkstraFibonacci,
  kDijkstraRadix,
};

enum class Format : std::uint8_t { kAuto, kEdgeList, kDimacs };

struct Cli {
  std::string graph;
  bool has_algorithm = false;
  Algorithm algorithm = Algorithm::kBellmanFord;
  bool has_source = false;
  sssp::Vertex source = 0;
  Format format = Format::kAuto;
  bool use_double = false;
  sssp::io::DimacsReadOptions dimacs{sssp::io::DimacsArcMode::kMergeReversePairs,
                                     sssp::io::DuplicateEdgePolicy::kKeepAll};
  std::string out;
  bool verify = false;
};

Algorithm parse_algorithm(std::string_view a) {
  if (a == "bellman_ford") return Algorithm::kBellmanFord;
  if (a == "dijkstra_binary") return Algorithm::kDijkstraBinary;
  if (a == "dijkstra_4ary") return Algorithm::kDijkstra4ary;
  if (a == "dijkstra_fibonacci") return Algorithm::kDijkstraFibonacci;
  if (a == "dijkstra_radix") return Algorithm::kDijkstraRadix;
  throw UsageError("unknown algorithm: '" + std::string(a) + "'");
}

sssp::Vertex parse_source(const std::string& text) {
  std::uint64_t value = 0;
  const char* end = text.c_str() + text.size();
  const auto [ptr, ec] = std::from_chars(text.c_str(), end, value);
  if (text.empty() || ec != std::errc{} || ptr != end || value >= sssp::kNoVertex) {
    throw UsageError("invalid value for --source: '" + text + "'");
  }
  return static_cast<sssp::Vertex>(value);
}

Cli parse(int argc, char** argv) {
  Cli c;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    auto value = [&]() -> std::string_view {
      if (i + 1 >= argc) throw UsageError("missing value for " + std::string(a));
      return argv[++i];
    };
    if (a == "--help" || a == "-h") {
      std::cout << kUsage;
      std::exit(0);
    } else if (a == "--graph") {
      c.graph = std::string(value());
    } else if (a == "--algorithm") {
      c.algorithm = parse_algorithm(value());
      c.has_algorithm = true;
    } else if (a == "--source") {
      c.source = parse_source(std::string(value()));
      c.has_source = true;
    } else if (a == "--format") {
      const auto f = value();
      if (f == "auto") {
        c.format = Format::kAuto;
      } else if (f == "edgelist") {
        c.format = Format::kEdgeList;
      } else if (f == "dimacs") {
        c.format = Format::kDimacs;
      } else {
        throw UsageError("--format must be auto|edgelist|dimacs");
      }
    } else if (a == "--weight-type") {
      const auto t = value();
      if (t != "uint64" && t != "double") throw UsageError("--weight-type must be uint64|double");
      c.use_double = (t == "double");
    } else if (a == "--dimacs-arcs") {
      const auto m = value();
      if (m == "merge") {
        c.dimacs.arc_mode = sssp::io::DimacsArcMode::kMergeReversePairs;
      } else if (m == "edge") {
        c.dimacs.arc_mode = sssp::io::DimacsArcMode::kArcAsEdge;
      } else {
        throw UsageError("--dimacs-arcs must be merge|edge");
      }
    } else if (a == "--duplicates") {
      const auto d = value();
      if (d == "keep") {
        c.dimacs.duplicates = sssp::io::DuplicateEdgePolicy::kKeepAll;
      } else if (d == "min") {
        c.dimacs.duplicates = sssp::io::DuplicateEdgePolicy::kKeepMinWeight;
      } else if (d == "reject") {
        c.dimacs.duplicates = sssp::io::DuplicateEdgePolicy::kReject;
      } else {
        throw UsageError("--duplicates must be keep|min|reject");
      }
    } else if (a == "--out") {
      c.out = std::string(value());
    } else if (a == "--verify") {
      c.verify = true;
    } else {
      throw UsageError("unknown option: " + std::string(a));
    }
  }
  if (c.graph.empty()) throw UsageError("--graph is required");
  if (!c.has_algorithm) throw UsageError("--algorithm is required");
  if (!c.has_source) throw UsageError("--source is required");
  if (c.use_double && c.algorithm == Algorithm::kDijkstraRadix) {
    throw UsageError("dijkstra_radix requires --weight-type uint64");
  }
  return c;
}

bool is_dimacs(const Cli& c) {
  if (c.format != Format::kAuto) return c.format == Format::kDimacs;
  const std::string_view p = c.graph;
  return p.size() >= 3 && p.substr(p.size() - 3) == ".gr";
}

template <class W>
sssp::SsspResult<W> solve(const Cli& c, const sssp::CsrGraph<W>& g) {
  namespace b = sssp::baseline;
  switch (c.algorithm) {
    case Algorithm::kBellmanFord:
      return b::bellman_ford(g, c.source);
    case Algorithm::kDijkstraBinary:
      return b::dijkstra_binary_heap(g, c.source);
    case Algorithm::kDijkstra4ary:
      return b::dijkstra_4ary_heap(g, c.source);
    case Algorithm::kDijkstraFibonacci:
      return b::dijkstra_fibonacci_heap(g, c.source);
    case Algorithm::kDijkstraRadix:
      if constexpr (std::is_same_v<W, std::uint64_t>) {
        return b::dijkstra_radix_heap(g, c.source);
      } else {
        throw UsageError("dijkstra_radix requires --weight-type uint64");
      }
  }
  throw std::logic_error("unhandled algorithm");
}

template <class W>
int run(const Cli& c) {
  auto el = is_dimacs(c) ? sssp::io::read_dimacs_file<W>(c.graph, c.dimacs)
                         : sssp::io::read_edge_list_file<W>(c.graph);
  const auto g = sssp::to_csr(std::move(el));
  const auto result = solve<W>(c, g);

  if (c.out.empty()) {
    sssp::write_result_json(std::cout, c.source, result);
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("failed to write to stdout");
  } else {
    std::ofstream os(c.out);
    if (!os) throw std::runtime_error("cannot open " + c.out + " for writing");
    sssp::write_result_json(os, c.source, result);
    os.flush();
    if (!os) throw std::runtime_error("failed to write " + c.out);
  }

  if (c.verify) {
    const auto rep = sssp::verify_sssp(g, c.source, result);
    if (!rep.ok) {
      std::cerr << "verify: FAILED\n" << rep.summary();
      return 3;
    }
    std::cerr << "verify: OK\n";
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Cli c = parse(argc, argv);
    return c.use_double ? run<double>(c) : run<std::uint64_t>(c);
  } catch (const UsageError& e) {
    std::cerr << "sssp_run: " << e.what() << "\n" << kUsage;
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "sssp_run: " << e.what() << "\n";
    return 1;
  }
}
