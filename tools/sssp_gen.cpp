// sssp_gen: reproducible graph generator CLI.
//
//   sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt
//   sssp_gen --list
//
// Output is the project edge-list format (default) or DIMACS .gr; the same arguments always
// produce byte-identical output (see docs/design/graph_library.md).

#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "sssp/generators.hpp"
#include "sssp/graph_io.hpp"

namespace {

constexpr std::string_view kUsage = R"(usage: sssp_gen --type TYPE [options]
  --list                      list generator types and weight distributions
  --type TYPE                 generator (see --list)
  --n N                       vertices (default 1000)
  --avg-degree D              er, chung_lu (default 4)
  --degree D                  regular (default 3)
  --multi keep|drop|switch    regular: self-loop/multi-edge policy (default switch)
  --attach M                  ba (default 2)
  --rows R --cols C           grid, road (default floor(sqrt(n)) each)
  --exponent G                chung_lu power-law exponent > 2 (default 2.5)
  --jitter J --keep-prob P --shortcut-prob P --shortcut-radius R --int-scale S   road
  --weights DIST              uniform_int|uniform_real|log_uniform|exponential|constant|
                              zero_mixed|wide_range (default uniform_int)
  --lo X --hi X --mean X --zero-prob P --wide-log2 K   weight parameters
  --weight-type uint64|double (default uint64)
  --seed S                    64-bit seed (default 1)
  --connected                 join components with a random tree of edges
  --format edgelist|dimacs    output format (default edgelist)
  --out PATH                  output file (default stdout)
  --stats                     print n, m, components to stderr
)";

struct UsageError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

template <class T>
T parse_number(std::string_view flag, std::string_view text) {
  T value{};
  const auto* end = text.data() + text.size();
  const auto [ptr, ec] = std::from_chars(text.data(), end, value);
  if (ec != std::errc{} || ptr != end) {
    throw UsageError("invalid value for " + std::string(flag) + ": '" + std::string(text) + "'");
  }
  return value;
}

struct Cli {
  sssp::gen::GraphSpec spec;
  bool list = false;
  bool stats = false;
  bool use_double = false;
  bool dimacs = false;
  std::string out;
};

Cli parse(int argc, char** argv) {
  Cli c;
  auto& s = c.spec;
  auto& w = s.options.weights;
  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    auto value = [&]() -> std::string_view {
      if (i + 1 >= argc) throw UsageError("missing value for " + std::string(a));
      return argv[++i];
    };
    if (a == "--help" || a == "-h") {
      std::cout << kUsage;
      std::exit(0);
    } else if (a == "--list") {
      c.list = true;
    } else if (a == "--stats") {
      c.stats = true;
    } else if (a == "--connected") {
      s.options.connected = true;
    } else if (a == "--type") {
      s.type = std::string(value());
    } else if (a == "--n") {
      s.n = parse_number<std::size_t>(a, value());
    } else if (a == "--avg-degree") {
      s.avg_degree = parse_number<double>(a, value());
    } else if (a == "--degree") {
      s.degree = parse_number<std::size_t>(a, value());
    } else if (a == "--multi") {
      s.multi = sssp::gen::parse_multi_edge_policy(value());
    } else if (a == "--attach") {
      s.attach = parse_number<std::size_t>(a, value());
    } else if (a == "--rows") {
      s.rows = parse_number<std::size_t>(a, value());
    } else if (a == "--cols") {
      s.cols = parse_number<std::size_t>(a, value());
    } else if (a == "--exponent") {
      s.exponent = parse_number<double>(a, value());
    } else if (a == "--jitter") {
      s.road.jitter = parse_number<double>(a, value());
    } else if (a == "--keep-prob") {
      s.road.keep_probability = parse_number<double>(a, value());
    } else if (a == "--shortcut-prob") {
      s.road.shortcut_probability = parse_number<double>(a, value());
    } else if (a == "--shortcut-radius") {
      s.road.shortcut_radius = parse_number<int>(a, value());
    } else if (a == "--int-scale") {
      s.road.integer_scale = parse_number<double>(a, value());
    } else if (a == "--weights") {
      w.dist = sssp::gen::parse_weight_dist(value());
    } else if (a == "--lo") {
      w.lo = parse_number<double>(a, value());
    } else if (a == "--hi") {
      w.hi = parse_number<double>(a, value());
    } else if (a == "--mean") {
      w.mean = parse_number<double>(a, value());
    } else if (a == "--zero-prob") {
      w.zero_probability = parse_number<double>(a, value());
    } else if (a == "--wide-log2") {
      w.wide_log2 = parse_number<int>(a, value());
    } else if (a == "--weight-type") {
      const auto t = value();
      if (t != "uint64" && t != "double") throw UsageError("--weight-type must be uint64|double");
      c.use_double = (t == "double");
    } else if (a == "--seed") {
      s.options.seed = parse_number<std::uint64_t>(a, value());
    } else if (a == "--format") {
      const auto f = value();
      if (f != "edgelist" && f != "dimacs") throw UsageError("--format must be edgelist|dimacs");
      c.dimacs = (f == "dimacs");
    } else if (a == "--out") {
      c.out = std::string(value());
    } else {
      throw UsageError("unknown option: " + std::string(a));
    }
  }
  return c;
}

void print_list() {
  std::cout << "generator types:\n";
  for (const auto& g : sssp::gen::generator_catalog()) {
    std::cout << "  " << g.name << "\n      " << g.description
              << "\n      parameters: " << g.parameters << "\n";
  }
  std::cout << "weight distributions (--weights):\n";
  for (auto d : {sssp::gen::WeightDist::kUniformInt, sssp::gen::WeightDist::kUniformReal,
                 sssp::gen::WeightDist::kLogUniform, sssp::gen::WeightDist::kExponential,
                 sssp::gen::WeightDist::kConstant, sssp::gen::WeightDist::kZeroMixed,
                 sssp::gen::WeightDist::kWideRange}) {
    std::cout << "  " << sssp::gen::to_string(d) << "\n";
  }
}

template <class W>
void run(const Cli& c) {
  const auto g = sssp::gen::generate<W>(c.spec);
  if (c.stats) {
    std::cerr << "n=" << g.num_vertices << " m=" << g.edges.size()
              << " components=" << sssp::gen::count_components(g) << "\n";
  }
  auto write = [&](std::ostream& os) {
    if (c.dimacs) {
      sssp::io::write_dimacs(os, g);
    } else {
      sssp::io::write_edge_list(os, g);
    }
  };
  if (c.out.empty()) {
    write(std::cout);
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("failed to write to stdout");
  } else if (c.dimacs) {
    sssp::io::write_dimacs_file(c.out, g);
  } else {
    sssp::io::write_edge_list_file(c.out, g);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Cli c = parse(argc, argv);
    if (c.list) {
      print_list();
      return 0;
    }
    if (c.use_double) {
      run<double>(c);
    } else {
      run<std::uint64_t>(c);
    }
    return 0;
  } catch (const UsageError& e) {
    std::cerr << "sssp_gen: " << e.what() << "\n" << kUsage;
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "sssp_gen: " << e.what() << "\n";
    return 1;
  }
}
