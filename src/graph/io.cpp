#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "sssp/graph_io.hpp"

namespace sssp::io {

ParseError::ParseError(std::size_t line, const std::string& reason)
    : std::runtime_error(line == 0 ? "line 0 (end of input): " + reason
                                   : "line " + std::to_string(line) + ": " + reason),
      line_(line) {}

namespace {

constexpr std::size_t kMaxTokens = 5;

/// Whitespace-separated tokens of one line (at most kMaxTokens are kept; count is exact).
struct Tokens {
  std::string_view tok[kMaxTokens];
  std::size_t count = 0;
};

bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
}

Tokens tokenize(std::string_view s) {
  Tokens t;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && is_space(s[i])) ++i;
    if (i == s.size()) break;
    const std::size_t start = i;
    while (i < s.size() && !is_space(s[i])) ++i;
    if (t.count < kMaxTokens) t.tok[t.count] = s.substr(start, i - start);
    ++t.count;
  }
  return t;
}

std::string quoted(std::string_view s) {
  return "'" + std::string(s) + "'";
}

std::uint64_t parse_u64(std::string_view s, std::size_t line, const char* what) {
  if (!s.empty() && s.front() == '-') throw ParseError(line, std::string("negative ") + what);
  std::uint64_t value = 0;
  const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value, 10);
  if (ec == std::errc::result_out_of_range) {
    throw ParseError(line, std::string(what) + " " + quoted(s) + " overflows 64 bits");
  }
  if (ec != std::errc{} || ptr != s.data() + s.size()) {
    throw ParseError(line, std::string(what) + " " + quoted(s) + " is not an unsigned integer");
  }
  return value;
}

template <class W>
W parse_weight(std::string_view s, std::size_t line) {
  if constexpr (std::is_same_v<W, std::uint64_t>) {
    return parse_u64(s, line, "weight");
  } else {
    static_assert(std::is_same_v<W, double>);
    double value = 0.0;
    const auto [ptr, ec] =
        std::from_chars(s.data(), s.data() + s.size(), value, std::chars_format::general);
    if (ec == std::errc::result_out_of_range) {
      throw ParseError(line, "weight " + quoted(s) + " is out of the double range");
    }
    if (ec != std::errc{} || ptr != s.data() + s.size()) {
      throw ParseError(line, "weight " + quoted(s) + " is not a number");
    }
    if (std::isnan(value)) throw ParseError(line, "weight is NaN");
    if (std::isinf(value)) throw ParseError(line, "weight is infinite");
    if (value < 0.0) throw ParseError(line, "negative weight " + quoted(s));
    return value == 0.0 ? 0.0 : value;  // -0 -> +0
  }
}

std::size_t parse_vertex_count(std::string_view s, std::size_t line) {
  const auto n = parse_u64(s, line, "vertex count");
  if (n >= static_cast<std::uint64_t>(kNoVertex)) {
    throw ParseError(line, "vertex count " + quoted(s) + " exceeds the 32-bit Vertex limit");
  }
  return static_cast<std::size_t>(n);
}

/// Avoids trusting a huge header for the up-front allocation.
std::size_t reserve_hint(std::uint64_t m) {
  constexpr std::uint64_t kCap = std::uint64_t{1} << 24;
  return static_cast<std::size_t>(m < kCap ? m : kCap);
}

template <class W>
std::uint64_t weight_bits(W w) {
  if constexpr (std::is_floating_point_v<W>) {
    return std::bit_cast<std::uint64_t>(w);
  } else {
    return static_cast<std::uint64_t>(w);
  }
}

template <class W>
void append_weight(std::string& buf, W w) {
  char tmp[64];
  std::to_chars_result r{};
  if constexpr (std::is_floating_point_v<W>) {
    r = std::to_chars(tmp, tmp + sizeof tmp, w, std::chars_format::general, 17);
  } else {
    r = std::to_chars(tmp, tmp + sizeof tmp, w);
  }
  buf.append(tmp, r.ptr);
}

void append_u64(std::string& buf, std::uint64_t x) {
  char tmp[24];
  const auto r = std::to_chars(tmp, tmp + sizeof tmp, x);
  buf.append(tmp, r.ptr);
}

void flush_buffer(std::ostream& out, std::string& buf, bool force) {
  if (force || buf.size() >= (1U << 16)) {
    out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
    buf.clear();
  }
}

struct PairHash {
  std::size_t operator()(std::uint64_t k) const noexcept {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    return static_cast<std::size_t>(k);
  }
};

struct ArcKey {
  Vertex u;
  Vertex v;
  std::uint64_t w;
  bool operator==(const ArcKey&) const = default;
};

struct ArcKeyHash {
  std::size_t operator()(const ArcKey& k) const noexcept {
    const PairHash h;
    return h((std::uint64_t{k.u} << 32 | k.v) ^ h(k.w));
  }
};

template <class W>
void apply_duplicate_policy(EdgeList<W>& g, const std::vector<std::size_t>& lines,
                            DuplicateEdgePolicy policy) {
  if (policy == DuplicateEdgePolicy::kKeepAll) return;
  std::unordered_map<std::uint64_t, std::size_t, PairHash> first;  // {min,max} -> output index
  first.reserve(g.edges.size());
  std::vector<Edge<W>> kept;
  kept.reserve(g.edges.size());
  for (std::size_t i = 0; i < g.edges.size(); ++i) {
    const auto& e = g.edges[i];
    const Vertex a = e.u < e.v ? e.u : e.v;
    const Vertex b = e.u < e.v ? e.v : e.u;
    const std::uint64_t key = std::uint64_t{a} << 32 | b;
    const auto [it, inserted] = first.try_emplace(key, kept.size());
    if (inserted) {
      kept.push_back(e);
    } else if (policy == DuplicateEdgePolicy::kReject) {
      throw ParseError(lines[i], "duplicate edge {" + std::to_string(a + 1) + "," +
                                     std::to_string(b + 1) + "} (DuplicateEdgePolicy::kReject)");
    } else if (e.weight < kept[it->second].weight) {
      kept[it->second].weight = e.weight;
    }
  }
  g.edges = std::move(kept);
}

std::ifstream open_input(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in.is_open()) throw std::runtime_error("cannot open graph file for reading: " + path);
  return in;
}

std::ofstream open_output(const std::string& path) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out.is_open()) throw std::runtime_error("cannot open graph file for writing: " + path);
  return out;
}

void finish_output(std::ofstream& out, const std::string& path) {
  out.flush();
  if (!out) throw std::runtime_error("failed to write graph file: " + path);
}

}  // namespace

// ------------------------------------------------------------------ edge-list format

template <class W>
EdgeList<W> read_edge_list(std::istream& in) {
  EdgeList<W> g;
  bool have_header = false;
  std::uint64_t m = 0;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const Tokens t = tokenize(line);
    if (t.count == 0 || t.tok[0].front() == '#') continue;
    if (!have_header) {
      if (t.count != 2) throw ParseError(lineno, "header must be 'n m'");
      g.num_vertices = parse_vertex_count(t.tok[0], lineno);
      m = parse_u64(t.tok[1], lineno, "edge count");
      if (m >= static_cast<std::uint64_t>(kNoEdge)) {
        throw ParseError(lineno, "edge count exceeds the 32-bit EdgeId limit");
      }
      g.edges.reserve(reserve_hint(m));
      have_header = true;
      continue;
    }
    if (t.count != 3) throw ParseError(lineno, "edge line must be 'u v w'");
    if (g.edges.size() == m) {
      throw ParseError(lineno, "more edge lines than the " + std::to_string(m) + " declared");
    }
    const auto u = parse_u64(t.tok[0], lineno, "vertex");
    const auto v = parse_u64(t.tok[1], lineno, "vertex");
    if (u >= g.num_vertices || v >= g.num_vertices) {
      throw ParseError(lineno, "vertex out of range (n = " + std::to_string(g.num_vertices) +
                                   ", vertices are 0-based)");
    }
    const W w = parse_weight<W>(t.tok[2], lineno);
    g.edges.push_back(Edge<W>{static_cast<Vertex>(u), static_cast<Vertex>(v), w});
  }
  if (!have_header) throw ParseError(0, "missing 'n m' header");
  if (g.edges.size() != m) {
    throw ParseError(0, "expected " + std::to_string(m) + " edge lines, found " +
                            std::to_string(g.edges.size()));
  }
  return g;
}

template <class W>
void write_edge_list(std::ostream& out, const EdgeList<W>& g) {
  std::string buf;
  append_u64(buf, g.num_vertices);
  buf += ' ';
  append_u64(buf, g.edges.size());
  buf += '\n';
  for (const auto& e : g.edges) {
    append_u64(buf, e.u);
    buf += ' ';
    append_u64(buf, e.v);
    buf += ' ';
    append_weight(buf, e.weight);
    buf += '\n';
    flush_buffer(out, buf, false);
  }
  flush_buffer(out, buf, true);
}

// ------------------------------------------------------------------ DIMACS format

template <class W>
EdgeList<W> read_dimacs(std::istream& in, const DimacsReadOptions& opts) {
  EdgeList<W> g;
  std::vector<std::size_t> edge_line;  // source line of each edge (for kReject)
  std::unordered_map<ArcKey, std::uint32_t, ArcKeyHash> pending;  // unmatched arcs (merge mode)
  const bool merge = opts.arc_mode == DimacsArcMode::kMergeReversePairs;
  bool have_problem = false;
  std::uint64_t m = 0;
  std::uint64_t arcs = 0;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const Tokens t = tokenize(line);
    if (t.count == 0 || t.tok[0] == "c") continue;
    if (t.tok[0] == "p") {
      if (have_problem) throw ParseError(lineno, "duplicate problem line");
      if (t.count != 4 || t.tok[1] != "sp")
        throw ParseError(lineno, "problem line must be 'p sp n m'");
      g.num_vertices = parse_vertex_count(t.tok[2], lineno);
      m = parse_u64(t.tok[3], lineno, "arc count");
      if (m >= 2 * static_cast<std::uint64_t>(kNoEdge)) {
        throw ParseError(lineno, "arc count exceeds the 32-bit EdgeId limit");
      }
      g.edges.reserve(reserve_hint(merge ? m / 2 : m));
      edge_line.reserve(reserve_hint(merge ? m / 2 : m));
      have_problem = true;
      continue;
    }
    if (t.tok[0] != "a") throw ParseError(lineno, "unknown line type " + quoted(t.tok[0]));
    if (!have_problem) throw ParseError(lineno, "arc line before the 'p sp n m' problem line");
    if (t.count != 4) throw ParseError(lineno, "arc line must be 'a u v w'");
    if (arcs == m) {
      throw ParseError(lineno, "more arc lines than the " + std::to_string(m) + " declared");
    }
    ++arcs;
    const auto u1 = parse_u64(t.tok[1], lineno, "vertex");
    const auto v1 = parse_u64(t.tok[2], lineno, "vertex");
    if (u1 == 0 || v1 == 0) throw ParseError(lineno, "vertex 0 is invalid (DIMACS is 1-based)");
    if (u1 > g.num_vertices || v1 > g.num_vertices) {
      throw ParseError(lineno, "vertex out of range (n = " + std::to_string(g.num_vertices) + ")");
    }
    const W w = parse_weight<W>(t.tok[3], lineno);
    const auto u = static_cast<Vertex>(u1 - 1);
    const auto v = static_cast<Vertex>(v1 - 1);
    if (merge && u != v) {
      const auto rev = pending.find(ArcKey{v, u, weight_bits(w)});
      if (rev != pending.end() && rev->second > 0) {
        --rev->second;  // second arc of a symmetric pair: the edge already exists
        continue;
      }
      ++pending[ArcKey{u, v, weight_bits(w)}];
    }
    if (g.edges.size() + 1 >= static_cast<std::size_t>(kNoEdge)) {
      throw ParseError(lineno, "edge count exceeds the 32-bit EdgeId limit");
    }
    g.edges.push_back(Edge<W>{u, v, w});
    edge_line.push_back(lineno);
  }
  if (!have_problem) throw ParseError(0, "missing 'p sp n m' problem line");
  if (arcs != m) {
    throw ParseError(0,
                     "expected " + std::to_string(m) + " arc lines, found " + std::to_string(arcs));
  }
  apply_duplicate_policy(g, edge_line, opts.duplicates);
  return g;
}

template <class W>
void write_dimacs(std::ostream& out, const EdgeList<W>& g) {
  std::uint64_t arcs = 0;
  for (const auto& e : g.edges) arcs += (e.u == e.v) ? 1 : 2;
  std::string buf = "c undirected graph: each edge {u,v} is written as arcs u->v and v->u\np sp ";
  append_u64(buf, g.num_vertices);
  buf += ' ';
  append_u64(buf, arcs);
  buf += '\n';
  const auto arc = [&buf](Vertex a, Vertex b, W w) {
    buf += "a ";
    append_u64(buf, std::uint64_t{a} + 1);
    buf += ' ';
    append_u64(buf, std::uint64_t{b} + 1);
    buf += ' ';
    append_weight(buf, w);
    buf += '\n';
  };
  for (const auto& e : g.edges) {
    arc(e.u, e.v, e.weight);
    if (e.u != e.v) arc(e.v, e.u, e.weight);
    flush_buffer(out, buf, false);
  }
  flush_buffer(out, buf, true);
}

// ------------------------------------------------------------------ file variants

template <class W>
EdgeList<W> read_edge_list_file(const std::string& path) {
  auto in = open_input(path);
  return read_edge_list<W>(in);
}

template <class W>
void write_edge_list_file(const std::string& path, const EdgeList<W>& g) {
  auto out = open_output(path);
  write_edge_list(out, g);
  finish_output(out, path);
}

template <class W>
EdgeList<W> read_dimacs_file(const std::string& path, const DimacsReadOptions& opts) {
  auto in = open_input(path);
  return read_dimacs<W>(in, opts);
}

template <class W>
void write_dimacs_file(const std::string& path, const EdgeList<W>& g) {
  auto out = open_output(path);
  write_dimacs(out, g);
  finish_output(out, path);
}

#define SSSP_IO_INSTANTIATE(W)                                                            \
  template EdgeList<W> read_edge_list<W>(std::istream&);                                  \
  template void write_edge_list<W>(std::ostream&, const EdgeList<W>&);                    \
  template EdgeList<W> read_dimacs<W>(std::istream&, const DimacsReadOptions&);           \
  template void write_dimacs<W>(std::ostream&, const EdgeList<W>&);                       \
  template EdgeList<W> read_edge_list_file<W>(const std::string&);                        \
  template void write_edge_list_file<W>(const std::string&, const EdgeList<W>&);          \
  template EdgeList<W> read_dimacs_file<W>(const std::string&, const DimacsReadOptions&); \
  template void write_dimacs_file<W>(const std::string&, const EdgeList<W>&);

SSSP_IO_INSTANTIATE(std::uint64_t)
SSSP_IO_INSTANTIATE(double)

#undef SSSP_IO_INSTANTIATE

}  // namespace sssp::io
