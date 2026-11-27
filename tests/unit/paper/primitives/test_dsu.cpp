// Sec 3.3, Lemma 3.6: disjoint sets with union by rank and full path compression.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <vector>

#include "paper/primitives/dsu.hpp"
#include "sssp/result.hpp"

using sssp::paper::DisjointSets;

namespace {

/// Naive model: every item stores its set label; a union relabels one whole set.
class NaiveSets {
 public:
  explicit NaiveSets(std::size_t k) : label_(k) {
    for (std::size_t i = 0; i < k; ++i) label_[i] = static_cast<std::uint32_t>(i);
  }
  bool same(std::uint32_t a, std::uint32_t b) const { return label_[a] == label_[b]; }
  bool unite(std::uint32_t a, std::uint32_t b) {
    const auto la = label_[a];
    const auto lb = label_[b];
    if (la == lb) return false;
    for (auto& l : label_) {
      if (l == lb) l = la;
    }
    return true;
  }

 private:
  std::vector<std::uint32_t> label_;
};

std::uint32_t floor_log2(std::uint64_t x) {
  std::uint32_t r = 0;
  while (x > 1) {
    x >>= 1;
    ++r;
  }
  return r;
}

/// Structural invariants from the proof of Lemma 3.6, checked without mutating `ds`.
void check_rank_invariants(const DisjointSets& ds) {
  const std::size_t k = ds.size();
  std::vector<std::uint64_t> members(k, 0);
  for (std::uint32_t x = 0; x < k; ++x) {
    std::uint32_t r = x;
    while (ds.parent(r) != r) {
      // Ranks increase strictly along parent chains.
      REQUIRE(ds.rank(ds.parent(r)) > ds.rank(r));
      r = ds.parent(r);
    }
    ++members[r];
  }
  std::size_t roots = 0;
  for (std::uint32_t x = 0; x < k; ++x) {
    // All ranks are at most floor(log2 k).
    REQUIRE(ds.rank(x) <= floor_log2(k));
    if (ds.parent(x) == x) {
      ++roots;
      // A root of rank r has at least 2^r original items.
      REQUIRE(members[x] >= (std::uint64_t{1} << ds.rank(x)));
    }
  }
  REQUIRE(roots == ds.num_sets());
}

}  // namespace

TEST_CASE("dsu starts with k singleton roots of rank zero", "[paper][dsu]") {
  DisjointSets ds(5);
  REQUIRE(ds.size() == 5);
  REQUIRE(ds.num_sets() == 5);
  for (std::uint32_t x = 0; x < 5; ++x) {
    REQUIRE(ds.parent(x) == x);
    REQUIRE(ds.rank(x) == 0);
    REQUIRE(ds.find(x) == x);
  }
}

TEST_CASE("dsu on a single item", "[paper][dsu]") {
  DisjointSets ds(1);
  REQUIRE(ds.find(0) == 0);
  REQUIRE_FALSE(ds.unite(0, 0));
  REQUIRE(ds.num_sets() == 1);
  REQUIRE(ds.rank(0) == 0);
}

TEST_CASE("dsu unite of a set with itself does nothing", "[paper][dsu]") {
  DisjointSets ds(4);
  REQUIRE(ds.unite(0, 1));
  REQUIRE_FALSE(ds.unite(1, 0));
  REQUIRE_FALSE(ds.unite(0, 1));
  REQUIRE_FALSE(ds.unite(2, 2));
  REQUIRE(ds.num_sets() == 3);
  REQUIRE(ds.find(0) == ds.find(1));
  REQUIRE(ds.find(2) != ds.find(0));
}

TEST_CASE("dsu links by rank and raises the survivor only on a tie", "[paper][dsu]") {
  DisjointSets ds(4);
  REQUIRE(ds.unite(0, 1));  // tie at rank 0
  const auto r01 = ds.find(0);
  REQUIRE(ds.rank(r01) == 1);
  REQUIRE(ds.unite(2, r01));  // rank 0 under rank 1: no promotion
  REQUIRE(ds.find(2) == r01);
  REQUIRE(ds.rank(r01) == 1);
  REQUIRE(ds.parent(2) == r01);
  REQUIRE(ds.unite(r01, 3));
  REQUIRE(ds.find(3) == r01);
  REQUIRE(ds.rank(r01) == 1);
  REQUIRE(ds.num_sets() == 1);
}

TEST_CASE("dsu find compresses every nonroot on the chain to the root", "[paper][dsu]") {
  // Build a binomial tree of rank 3 on 8 items; its deepest item has a chain of 3 edges.
  DisjointSets ds(8);
  for (std::uint32_t step = 1; step < 8; step <<= 1) {
    for (std::uint32_t i = 0; i + step < 8; i += 2 * step) REQUIRE(ds.unite(i, i + step));
  }
  check_rank_invariants(ds);
  // Find an item with chain length >= 2 before compression.
  std::uint32_t deep = 0;
  std::uint32_t best = 0;
  for (std::uint32_t x = 0; x < 8; ++x) {
    std::uint32_t len = 0;
    for (std::uint32_t r = x; ds.parent(r) != r; r = ds.parent(r)) ++len;
    if (len > best) {
      best = len;
      deep = x;
    }
  }
  REQUIRE(best >= 2);
  std::vector<std::uint32_t> chain;
  for (std::uint32_t r = deep; ds.parent(r) != r; r = ds.parent(r)) chain.push_back(r);
  const auto root = ds.find(deep);
  for (const auto x : chain) REQUIRE(ds.parent(x) == root);
  check_rank_invariants(ds);
}

TEST_CASE("dsu matches naive labels under random ops", "[paper][dsu]") {
  for (const std::uint32_t seed : {1u, 2u, 3u, 4u, 5u}) {
    for (const std::size_t k :
         {std::size_t{1}, std::size_t{2}, std::size_t{7}, std::size_t{64}, std::size_t{300}}) {
      std::mt19937_64 rng(seed * 1000003u + k);
      std::uniform_int_distribution<std::uint32_t> item(0, static_cast<std::uint32_t>(k - 1));
      DisjointSets ds(k);
      NaiveSets naive(k);
      std::size_t sets = k;
      for (int op = 0; op < 4 * static_cast<int>(k) + 20; ++op) {
        const auto a = item(rng);
        const auto b = item(rng);
        if (rng() % 2 == 0) {
          const bool merged = naive.unite(a, b);
          REQUIRE(ds.unite(a, b) == merged);
          if (merged) --sets;
        } else {
          REQUIRE((ds.find(a) == ds.find(b)) == naive.same(a, b));
        }
        REQUIRE(ds.num_sets() == sets);
      }
      check_rank_invariants(ds);
      for (std::uint32_t a = 0; a < k; ++a) {
        for (std::uint32_t b = 0; b < k && b < 40; ++b) {
          REQUIRE((ds.find(a) == ds.find(b)) == naive.same(a, b));
        }
      }
    }
  }
}

#if SSSP_DIAGNOSTICS
TEST_CASE("dsu counts finds and find traversals into diagnostics", "[paper][dsu][diagnostics]") {
  DisjointSets ds(4);
  REQUIRE(ds.find_calls() == 0);
  REQUIRE(ds.find_steps() == 0);
  (void)ds.find(2);  // root: no traversal
  REQUIRE(ds.find_calls() == 1);
  REQUIRE(ds.find_steps() == 0);
  REQUIRE(ds.unite(0, 1));  // two finds, both on roots
  REQUIRE(ds.find_calls() == 3);
  REQUIRE(ds.find_steps() == 0);
  const auto nonroot = ds.parent(0) == 0 ? 1u : 0u;
  (void)ds.find(nonroot);  // one parent-pointer traversal
  REQUIRE(ds.find_calls() == 4);
  REQUIRE(ds.find_steps() == 1);
  REQUIRE(ds.unions() == 1);

  sssp::Diagnostics diag;
  diag.counters["dsu.find_steps"] = 10;
  ds.add_counters(diag);
  REQUIRE(diag.counters.at("dsu.find_calls") == 4);
  REQUIRE(diag.counters.at("dsu.find_steps") == 11);
  REQUIRE(diag.counters.at("dsu.unions") == 1);
}
#endif
