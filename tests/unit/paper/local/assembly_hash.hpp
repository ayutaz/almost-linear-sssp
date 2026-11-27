#pragma once
// Test support for task P5-LOCAL-UNIFY: a 64-bit FNV-1a hash over every field of a
// paper::LocalAssembly (hierarchy, maps, block scales, classes, grid, supports) and over a
// counter map. It pins the output of the Sec 5.4 assembly across a refactor: two assemblies
// with the same hash agree on every field unless FNV collides.
//
// Lengths enter by their value representation (std::uint64_t as is, double by its IEEE bit
// pattern), every other field as an unsigned 64-bit integer, and every vector with its size.
// The assembly performs only comparisons, additions and doublings on lengths and the test
// graphs carry integers or quarter-integers, so the hashes do not depend on the platform.

#include <bit>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/base_case.hpp"

namespace sssp::test::local {

class Fnv1a {
 public:
  void word(std::uint64_t x) noexcept {
    for (int i = 0; i < 8; ++i) {
      h_ ^= (x >> (8 * i)) & 0xffU;
      h_ *= 1099511628211ULL;
    }
  }

  template <class T>
  void value(T x) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
      word(std::bit_cast<std::uint64_t>(static_cast<double>(x)));
    } else if constexpr (std::is_enum_v<T>) {
      word(static_cast<std::uint64_t>(static_cast<std::underlying_type_t<T>>(x)));
    } else {
      word(static_cast<std::uint64_t>(x));
    }
  }

  template <class T>
  void values(const std::vector<T>& v) noexcept {
    word(v.size());
    for (const T& x : v) value(x);
  }

  void text(const std::string& s) noexcept {
    word(s.size());
    for (const char c : s) word(static_cast<unsigned char>(c));
  }

  [[nodiscard]] std::uint64_t digest() const noexcept { return h_; }

 private:
  std::uint64_t h_ = 14695981039346656037ULL;
};

template <class W>
void mix_hierarchy(Fnv1a& f, const sssp::paper::OrderedHierarchy<W>& h) {
  f.word(h.c0);
  f.word(h.graph.num_states);
  f.word(h.graph.arcs.size());
  for (const auto& a : h.graph.arcs) {
    f.value(a.tail);
    f.value(a.head);
    f.value(a.length);
  }
  f.values(h.projection);
  f.values(h.leaf_of_state);
  f.values(h.roots);
  f.word(h.nodes.size());
  for (const auto& node : h.nodes) {
    f.value(node.parent);
    f.values(node.children);
    f.value(node.state);
    f.value(node.height);
    f.word(node.links.size());
    for (const auto& l : node.links) {
      f.value(l.x);
      f.value(l.y);
      f.value(l.length);
      f.value(l.arc);
    }
  }
}

/// Every field of the assembly, in declaration order.
template <class W>
void mix_assembly(Fnv1a& f, const sssp::paper::LocalAssembly<W>& a) {
  mix_hierarchy(f, a.result.hierarchy);
  f.values(a.result.native_entry);
  f.values(a.result.arc_edge);
  f.word(a.blocks.size());
  for (const auto& b : a.blocks) {
    f.value(b.kind);
    f.value(b.size);
    f.value(b.s);
    f.value(b.t);
    f.value(b.z);
    f.value(b.sigma_key);
  }
  f.values(a.y_component);
  f.values(a.node_upper);
  f.values(a.node_upper_inclusive);
  f.values(a.node_origin);
  f.values(a.grid);
  f.values(a.grid_upper);
  f.word(a.supports.size());
  for (const auto& s : a.supports) {
    f.value(s.lower);
    f.value(s.upper);
    f.value(s.lower_index);
    f.value(s.upper_index);
  }
}

/// Counter names and values, in the map's (lexicographic) order.
inline void mix_counters(Fnv1a& f, const std::map<std::string, std::uint64_t>& counters) {
  f.word(counters.size());
  for (const auto& [name, value] : counters) {
    f.text(name);
    f.word(value);
  }
}

}  // namespace sssp::test::local
