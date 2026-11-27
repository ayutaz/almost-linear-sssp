// Unit tests for the internal union-find shared by src/graph (P2-UNIONFIND-DEDUP):
// sssp::detail::UnionFind in src/graph/union_find.hpp. Checked against an independent
// brute-force partition (every element carries a set label; a union relabels one whole set).

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "graph/union_find.hpp"
#include "sssp/generators.hpp"

namespace {

/// Brute-force partition oracle: O(n) per union, obviously correct.
class LabelPartition {
 public:
  explicit LabelPartition(std::size_t n) : label_(n) {
    for (std::size_t i = 0; i < n; ++i) label_[i] = i;
  }
  [[nodiscard]] bool same(std::size_t a, std::size_t b) const { return label_[a] == label_[b]; }
  void merge(std::size_t a, std::size_t b) {
    const std::size_t from = label_[b];
    const std::size_t to = label_[a];
    for (auto& l : label_) {
      if (l == from) l = to;
    }
  }
  [[nodiscard]] std::size_t set_size(std::size_t a) const {
    std::size_t s = 0;
    for (const auto l : label_) s += l == label_[a] ? 1 : 0;
    return s;
  }

 private:
  std::vector<std::size_t> label_;
};

}  // namespace

TEMPLATE_TEST_CASE("UnionFind starts with singletons", "[graph][union_find]", sssp::Vertex,
                   std::size_t) {
  sssp::detail::UnionFind<TestType> uf(5);
  for (TestType v = 0; v < 5; ++v) CHECK(uf.find(v) == v);
  sssp::detail::UnionFind<TestType> empty(0);
  (void)empty;
}

TEMPLATE_TEST_CASE("UnionFind::unite reports whether two different sets were merged",
                   "[graph][union_find]", sssp::Vertex, std::size_t) {
  sssp::detail::UnionFind<TestType> uf(6);
  CHECK(uf.unite(0, 1));
  CHECK(uf.find(0) == uf.find(1));
  CHECK_FALSE(uf.unite(1, 0));
  CHECK_FALSE(uf.unite(2, 2));
  CHECK(uf.unite(2, 3));
  CHECK(uf.find(0) != uf.find(2));
  CHECK(uf.unite(3, 1));
  CHECK(uf.find(0) == uf.find(2));
  CHECK_FALSE(uf.unite(0, 3));
  CHECK(uf.find(4) == 4);
  CHECK(uf.find(5) == 5);
}

TEMPLATE_TEST_CASE("UnionFind unites by size: the root of the larger set stays the root",
                   "[graph][union_find]", sssp::Vertex, std::size_t) {
  sssp::detail::UnionFind<TestType> uf(8);
  // {5, 6, 7} has root r; a singleton joined from either side hangs below r.
  REQUIRE(uf.unite(5, 6));
  REQUIRE(uf.unite(6, 7));
  const TestType r = uf.find(7);
  REQUIRE(uf.unite(0, 7));
  CHECK(uf.find(0) == r);
  REQUIRE(uf.unite(5, 1));
  CHECK(uf.find(1) == r);
  // Equal sizes: the first argument's root becomes the root.
  REQUIRE(uf.unite(2, 3));
  CHECK(uf.find(3) == 2);
  REQUIRE(uf.unite(4, 2));
  CHECK(uf.find(2) == 2);  // {2, 3} is larger than {4}
  CHECK(uf.find(4) == 2);
}

TEMPLATE_TEST_CASE("UnionFind matches a brute-force partition on random unions",
                   "[graph][union_find]", sssp::Vertex, std::size_t) {
  for (std::uint64_t seed = 1; seed <= 20; ++seed) {
    sssp::gen::Xoshiro256ss rng(seed);
    const std::size_t n = 1 + static_cast<std::size_t>(rng.uniform_below(80));
    sssp::detail::UnionFind<TestType> uf(n);
    LabelPartition oracle(n);
    for (std::size_t step = 0; step < 3 * n; ++step) {
      const auto a = static_cast<TestType>(rng.uniform_below(n));
      const auto b = static_cast<TestType>(rng.uniform_below(n));
      const bool different = !oracle.same(a, b);
      INFO("seed " << seed << " step " << step << " unite(" << a << ", " << b << ")");
      CHECK(uf.unite(a, b) == different);
      if (different) oracle.merge(a, b);
      // The partition (not the particular roots) must agree with the oracle.
      for (std::size_t x = 0; x < n; ++x) {
        const auto tx = static_cast<TestType>(x);
        CHECK((uf.find(tx) == uf.find(a)) == oracle.same(x, a));
        CHECK(oracle.same(uf.find(tx), x));  // the root is a member of x's set
      }
      // Union by size: the root of a's set is the root of every member.
      std::size_t members = 0;
      for (std::size_t x = 0; x < n; ++x) {
        members += uf.find(static_cast<TestType>(x)) == uf.find(a) ? 1 : 0;
      }
      CHECK(members == oracle.set_size(a));
    }
  }
}
