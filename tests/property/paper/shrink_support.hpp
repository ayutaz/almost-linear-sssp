#pragma once
// Seeded property driver with deterministic shrinking, for the P4-PROPERTY lemma tests.
//
// Every random case is a pure function of one 64-bit seed: `generate(seed)` builds the input,
// `violation(input)` returns a message when the property fails and std::nullopt when it holds or
// when the input violates the property's precondition (a shrunk input that stops satisfying the
// precondition therefore never counts as a failure). On the first failing seed the input is
// shrunk greedily to a 1-minimal one: repeatedly delete one edge, then lower one weight to a
// smaller candidate, then drop isolated vertices, keeping each step only if the property still
// fails, until no single step does. Shrinking uses no randomness, so the reported seed always
// regenerates the same original input and the same minimal counterexample.
//
// The failure message names the seed, the original and the minimal input as an edge list, so a
// failure can be replayed with `shrink_case(generate(seed), violation)`.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::test::shrink {

/// A small undirected graph: vertices 0 .. n-1 and an edge list, plus the separator size b of
/// the owner construction (unused by properties that do not need it).
template <class W>
struct Case {
  std::size_t n = 0;
  std::vector<Edge<W>> edges;
  std::size_t b = 1;

  friend bool operator==(const Case& x, const Case& y) {
    if (x.n != y.n || x.b != y.b || x.edges.size() != y.edges.size()) return false;
    for (std::size_t i = 0; i < x.edges.size(); ++i) {
      const auto& e = x.edges[i];
      const auto& f = y.edges[i];
      if (e.u != f.u || e.v != f.v || !(e.weight == f.weight)) return false;
    }
    return true;
  }
};

template <class W>
std::string describe(const Case<W>& c) {
  std::ostringstream os;
  os.precision(17);
  os << "n=" << c.n << " b=" << c.b << " edges={";
  for (std::size_t i = 0; i < c.edges.size(); ++i) {
    const auto& e = c.edges[i];
    os << (i == 0 ? "" : ", ") << "{" << e.u << "," << e.v << "," << e.weight << "}";
  }
  os << "}";
  return os.str();
}

/// Removes the vertices without incident edges (keeping at least one) and renumbers the others
/// in increasing order.
template <class W>
Case<W> drop_isolated(const Case<W>& c) {
  std::vector<std::uint32_t> id(c.n, 0xffffffffu);
  for (const auto& e : c.edges) {
    id[e.u] = 0;
    id[e.v] = 0;
  }
  Case<W> out;
  out.b = c.b;
  for (std::size_t v = 0; v < c.n; ++v) {
    if (id[v] == 0) id[v] = static_cast<std::uint32_t>(out.n++);
  }
  if (out.n == 0) out.n = std::min<std::size_t>(c.n, 1);
  for (const auto& e : c.edges) out.edges.push_back({id[e.u], id[e.v], e.weight});
  return out;
}

struct ShrinkResult {
  std::string message;       ///< the property's message on the minimal input
  std::size_t steps = 0;     ///< accepted shrinking steps
  std::size_t attempts = 0;  ///< property evaluations during shrinking
};

/// Shrinks a failing case in place to a 1-minimal failing case (see the file comment).
/// `weights` lists the candidate replacement weights, tried from the smallest; a weight is only
/// replaced by a strictly smaller one, so the loop terminates. Requires violation(c) to fail.
template <class W, class Violation>
ShrinkResult shrink_case(Case<W>& c, Violation&& violation, const std::vector<W>& weights = {}) {
  ShrinkResult r;
  auto first = violation(std::as_const(c));
  r.message = first.value_or("");
  std::vector<W> cand = weights;
  std::sort(cand.begin(), cand.end());
  const auto accept = [&](Case<W>& next) {
    ++r.attempts;
    if (auto m = violation(std::as_const(next))) {
      c = std::move(next);
      r.message = *m;
      ++r.steps;
      return true;
    }
    return false;
  };
  for (bool progress = true; progress;) {
    progress = false;
    for (std::size_t i = 0; i < c.edges.size();) {
      Case<W> next = c;
      next.edges.erase(next.edges.begin() + static_cast<std::ptrdiff_t>(i));
      if (accept(next)) {
        progress = true;
      } else {
        ++i;
      }
    }
    for (std::size_t i = 0; i < c.edges.size(); ++i) {
      for (const W w : cand) {
        if (!(w < c.edges[i].weight)) break;
        Case<W> next = c;
        next.edges[i].weight = w;
        if (accept(next)) {
          progress = true;
          break;
        }
      }
    }
    {
      Case<W> next = drop_isolated(c);
      if (!(next == c) && accept(next)) progress = true;
    }
    if (c.b > 1) {
      for (std::size_t b = 1; b < c.b; ++b) {
        Case<W> next = c;
        next.b = b;
        if (accept(next)) {
          progress = true;
          break;
        }
      }
    }
  }
  return r;
}

/// Outcome of a seeded run: the first failing seed (if any), with its original and shrunk case.
template <class W>
struct SeededOutcome {
  std::uint64_t cases = 0;
  std::optional<std::uint64_t> failing_seed;
  Case<W> original;
  Case<W> minimal;
  ShrinkResult shrink;

  [[nodiscard]] std::string report() const {
    if (!failing_seed) return "no failure in " + std::to_string(cases) + " cases";
    std::ostringstream os;
    os << "seed " << *failing_seed << " fails: " << shrink.message
       << "\n  original: " << describe(original) << "\n  minimal (" << shrink.steps
       << " shrink steps): " << describe(minimal);
    return os.str();
  }
};

/// Runs `violation` on generate(seed) for seed = first .. first + count - 1, stops at the first
/// failure and shrinks it.
template <class W, class Generate, class Violation>
SeededOutcome<W> run_seeded(std::uint64_t first, std::uint64_t count, Generate&& generate,
                            Violation&& violation, const std::vector<W>& weights = {}) {
  SeededOutcome<W> out;
  for (std::uint64_t s = first; s < first + count; ++s) {
    Case<W> c = generate(s);
    ++out.cases;
    if (!violation(std::as_const(c))) continue;
    out.failing_seed = s;
    out.original = c;
    out.shrink = shrink_case(c, violation, weights);
    out.minimal = c;
    break;
  }
  return out;
}

/// True iff the edges connect all n vertices (n >= 1).
template <class W>
bool connected(std::size_t n, const std::vector<Edge<W>>& edges) {
  if (n == 0) return false;
  std::vector<std::size_t> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = i;
  const auto find = [&](std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  };
  std::size_t comps = n;
  for (const auto& e : edges) {
    const std::size_t a = find(e.u);
    const std::size_t b = find(e.v);
    if (a != b) {
      p[a] = b;
      --comps;
    }
  }
  return comps == 1;
}

}  // namespace sssp::test::shrink
