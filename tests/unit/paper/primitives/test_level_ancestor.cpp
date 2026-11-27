// Sec 3.4, Lemma 3.8: level ancestors via macro-micro marking, jump pointers and ladders.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "paper/primitives/forest_queries.hpp"
#include "unit/paper/primitives/forest_gen.hpp"

using forest_gen::Parents;
using sssp::kNoVertex;
using sssp::Vertex;
using sssp::paper::LevelAncestor;
using sssp::paper::OrderedTree;

namespace {

std::uint32_t floor_log2(std::uint64_t x) {
  std::uint32_t r = 0;
  while (x > 1) {
    x >>= 1;
    ++r;
  }
  return r;
}

void check_all_queries(const Parents& parent, std::uint32_t block_size) {
  const OrderedTree tree(parent);
  const LevelAncestor la = block_size == 0 ? LevelAncestor(tree) : LevelAncestor(tree, block_size);
  const std::size_t k = parent.size();
  const std::uint32_t b = la.block_size();
  // Lemma 3.8: components left after deleting the marks have fewer than B vertices.
  REQUIRE(la.max_micro_tree_size() < b);
  // Proof of Lemma 3.8: each mark but the root pays for >= B residual vertices.
  REQUIRE(la.num_marks() <= k / b + 1);
  // Upward extensions at most double the total path array length.
  REQUIRE(la.ladder_entries() <= 2 * k);
  REQUIRE(la.jump_entries() <= la.num_marks() * (floor_log2(k) + 1));
  for (Vertex v = 0; v < k; ++v) {
    const auto dv = forest_gen::naive_depth(parent, v);
    REQUIRE(la.depth(v) == dv);
    for (std::uint32_t j = 0; j <= dv; ++j) {
      INFO("v " << v << " j " << j << " B " << b);
      REQUIRE(la.query(v, j) == forest_gen::naive_level_ancestor(parent, v, j));
    }
  }
}

}  // namespace

TEST_CASE("level ancestor on a hand-built tree", "[paper][la]") {
  // 0 -> {1, 2}, 1 -> {3, 4}, 4 -> {5}, 5 -> {7}, 2 -> {6}
  const Parents parent{kNoVertex, 0, 0, 1, 1, 4, 2, 5};
  const OrderedTree tree(parent);
  for (const std::uint32_t b : {0u, 2u, 3u, 4u, 5u}) {
    const LevelAncestor la = b == 0 ? LevelAncestor(tree) : LevelAncestor(tree, b);
    REQUIRE(la.query(7, 4) == 7);
    REQUIRE(la.query(7, 3) == 5);
    REQUIRE(la.query(7, 2) == 4);
    REQUIRE(la.query(7, 1) == 1);
    REQUIRE(la.query(7, 0) == 0);
    REQUIRE(la.query(6, 1) == 2);
    REQUIRE(la.query(3, 1) == 1);
    REQUIRE(la.query(0, 0) == 0);
  }
}

TEST_CASE("level ancestor of a single-vertex tree", "[paper][la]") {
  const OrderedTree tree(Parents{kNoVertex});
  const LevelAncestor la(tree);
  REQUIRE(la.query(0, 0) == 0);
  REQUIRE(la.depth(0) == 0);
  REQUIRE(la.num_marks() == 1);  // the root is always marked
}

TEST_CASE("level ancestor rejects depths below the query vertex", "[paper][la]") {
  const OrderedTree tree(forest_gen::path(6));
  const LevelAncestor la(tree);
  REQUIRE(la.query(5, 5) == 5);
  REQUIRE_THROWS_AS(la.query(5, 6), std::out_of_range);
  REQUIRE_THROWS_AS(la.query(0, 1), std::out_of_range);
}

TEST_CASE("level ancestor test block size must lie in [2, 8]", "[paper][la]") {
  const OrderedTree tree(forest_gen::path(5));
  REQUIRE_THROWS_AS(LevelAncestor(tree, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(LevelAncestor(tree, 9), std::invalid_argument);
}

TEST_CASE("level ancestor equals naive for all (v, j) on path, star, caterpillar and binary trees",
          "[paper][la]") {
  std::mt19937_64 rng(9);
  for (const std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{10},
                              std::size_t{40}, std::size_t{130}}) {
    for (const std::uint32_t b : {0u, 2u, 3u, 4u, 6u}) {
      check_all_queries(forest_gen::path(n), b);
      check_all_queries(forest_gen::star(n), b);
      check_all_queries(forest_gen::caterpillar(n, n / 2 + 1), b);
      check_all_queries(forest_gen::binary(n), b);
      check_all_queries(forest_gen::shuffle_labels(forest_gen::caterpillar(n, 5), rng), b);
    }
  }
}

TEST_CASE("level ancestor equals naive for all (v, j) on small random trees", "[paper][la]") {
  for (const std::uint64_t seed : {1u, 2u, 3u, 4u, 5u, 6u}) {
    std::mt19937_64 rng(seed + 100);
    for (std::size_t n = 1; n <= 80; n += 1 + n / 3) {
      for (const std::uint32_t b : {0u, 2u, 3u, 5u, 7u}) {
        check_all_queries(forest_gen::shuffle_labels(forest_gen::random_tree(n, rng), rng), b);
        check_all_queries(forest_gen::shuffle_labels(forest_gen::random_deep_tree(n, 2, rng), rng),
                          b);
      }
    }
  }
}

// ---------------------------------------------------------------- Lemma 3.9 (forests)

namespace {

template <class Q>
void check_forest_all_queries(const Q& q, const Parents& parent) {
  for (Vertex v = 0; v < parent.size(); ++v) {
    const auto dv = forest_gen::naive_depth(parent, v);
    REQUIRE(q.depth(v) == dv);
    REQUIRE(q.component(v) == forest_gen::naive_root(parent, v));
    for (std::uint32_t j = 0; j <= dv; ++j) {
      INFO("v " << v << " j " << j);
      REQUIRE(q.level_ancestor(v, j) == forest_gen::naive_level_ancestor(parent, v, j));
    }
    // Lemma 3.9: the requested depth must lie between zero and the original depth.
    REQUIRE_THROWS_AS(q.level_ancestor(v, dv + 1), std::out_of_range);
  }
}

}  // namespace

TEST_CASE("forest level ancestor uses original depths and stays inside the tree",
          "[paper][la][forest]") {
  // Trees {0 -> 1 -> 2}, {3 -> 4}; roots have depth 0 in the original forest.
  const Parents parent{kNoVertex, 0, 1, kNoVertex, 3};
  const sssp::paper::ForestQueries fq(parent);
  REQUIRE(fq.depth(0) == 0);
  REQUIRE(fq.depth(2) == 2);
  REQUIRE(fq.level_ancestor(2, 0) == 0);
  REQUIRE(fq.level_ancestor(2, 1) == 1);
  REQUIRE(fq.level_ancestor(4, 0) == 3);
  REQUIRE(fq.level_ancestor(3, 0) == 3);
  REQUIRE_THROWS_AS(fq.level_ancestor(3, 1), std::out_of_range);
  check_forest_all_queries(fq, parent);
  check_forest_all_queries(sssp::paper::BinaryLiftingRef(parent), parent);
}

TEST_CASE("forest level ancestor equals naive for all (v, j) in small random forests",
          "[paper][la][forest]") {
  for (const std::uint64_t seed : {41u, 42u, 43u, 44u, 45u}) {
    std::mt19937_64 rng(seed);
    for (std::size_t n = 1; n <= 90; n += 1 + n / 4) {
      for (const std::size_t trees : {std::size_t{1}, std::size_t{3}, n / 2 + 1}) {
        const auto parent =
            forest_gen::shuffle_labels(forest_gen::random_forest(n, trees, rng), rng);
        check_forest_all_queries(sssp::paper::ForestQueries(parent), parent);
        check_forest_all_queries(sssp::paper::ForestQueries(parent, 3), parent);
        check_forest_all_queries(sssp::paper::ForestQueries(parent, 5), parent);
        check_forest_all_queries(sssp::paper::BinaryLiftingRef(parent), parent);
      }
    }
  }
}
