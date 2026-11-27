// Sec 3.4, Lemma 3.7: lowest common ancestors from an Euler tour with +-1 block tables.
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
using sssp::paper::LcaTable;
using sssp::paper::OrderedTree;
using sssp::paper::paper_block_size;

namespace {

void check_all_pairs(const Parents& parent, std::uint32_t block_size) {
  const OrderedTree tree(parent);
  const LcaTable lca = block_size == 0 ? LcaTable(tree) : LcaTable(tree, block_size);
  REQUIRE(lca.tour_length() == 2 * parent.size() - 1);
  for (Vertex u = 0; u < parent.size(); ++u) {
    for (Vertex v = 0; v < parent.size(); ++v) {
      const auto expected = forest_gen::naive_lca(parent, u, v);
      REQUIRE(expected.has_value());
      INFO("u " << u << " v " << v << " B " << lca.block_size());
      REQUIRE(lca.lca(u, v) == *expected);
    }
  }
}

}  // namespace

TEST_CASE("lca block size follows Sec 3.4 B = max(2, floor(log2(k+1)/8))", "[paper][lca]") {
  REQUIRE(paper_block_size(1) == 2);
  REQUIRE(paper_block_size(1000) == 2);
  REQUIRE(paper_block_size((std::size_t{1} << 24) - 2) == 2);  // log2(k+1) just below 24
  REQUIRE(paper_block_size((std::size_t{1} << 24) - 1) == 3);  // log2(k+1) = 24
  REQUIRE(paper_block_size((std::size_t{1} << 32) - 2) == 3);
  REQUIRE(paper_block_size((std::size_t{1} << 32) - 1) == 4);
  REQUIRE(paper_block_size((std::size_t{1} << 40) - 1) == 5);
}

TEST_CASE("ordered tree rejects parent arrays that are not one rooted tree", "[paper][lca]") {
  REQUIRE_THROWS_AS(OrderedTree(Parents{}), std::invalid_argument);                      // no root
  REQUIRE_THROWS_AS(OrderedTree(Parents{kNoVertex, kNoVertex}), std::invalid_argument);  // two
  REQUIRE_THROWS_AS(OrderedTree(Parents{kNoVertex, 5}), std::out_of_range);
  REQUIRE_THROWS_AS(OrderedTree(Parents{kNoVertex, 2, 1}), std::invalid_argument);  // cycle
  REQUIRE_THROWS_AS(OrderedTree(Parents{kNoVertex, 1}), std::invalid_argument);     // self-loop
}

TEST_CASE("ordered tree records root, depths and children in increasing index order",
          "[paper][lca]") {
  //        2
  //      /   \.
  //     0     3
  //     |
  //     1
  const OrderedTree t(Parents{2, 0, kNoVertex, 2});
  REQUIRE(t.size() == 4);
  REQUIRE(t.root() == 2);
  REQUIRE(t.depth(2) == 0);
  REQUIRE(t.depth(0) == 1);
  REQUIRE(t.depth(3) == 1);
  REQUIRE(t.depth(1) == 2);
  REQUIRE(t.parent(1) == 0);
  REQUIRE(t.parent(2) == kNoVertex);
  const auto ch = t.children(2);
  REQUIRE(ch.size() == 2);
  REQUIRE(ch[0] == 0);
  REQUIRE(ch[1] == 3);
}

TEST_CASE("lca on a hand-built tree", "[paper][lca]") {
  // 0 -> {1, 2}, 1 -> {3, 4}, 4 -> {5}, 2 -> {6}
  const Parents parent{kNoVertex, 0, 0, 1, 1, 4, 2};
  const OrderedTree tree(parent);
  for (const std::uint32_t b : {0u, 2u, 3u, 4u, 5u}) {
    const LcaTable lca = b == 0 ? LcaTable(tree) : LcaTable(tree, b);
    REQUIRE(lca.lca(3, 5) == 1);
    REQUIRE(lca.lca(5, 3) == 1);
    REQUIRE(lca.lca(5, 6) == 0);
    REQUIRE(lca.lca(4, 5) == 4);
    REQUIRE(lca.lca(3, 3) == 3);
    REQUIRE(lca.lca(6, 2) == 2);
    REQUIRE(lca.lca(0, 5) == 0);
    REQUIRE(lca.lca(0, 0) == 0);
  }
}

TEST_CASE("lca of a single-vertex tree is the vertex", "[paper][lca]") {
  const OrderedTree tree(Parents{kNoVertex});
  const LcaTable lca(tree);
  REQUIRE(lca.tour_length() == 1);
  REQUIRE(lca.lca(0, 0) == 0);
}

TEST_CASE("lca equals naive on all pairs of path, star, caterpillar and binary trees",
          "[paper][lca]") {
  std::mt19937_64 rng(7);
  for (const std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{9},
                              std::size_t{33}, std::size_t{70}}) {
    for (const std::uint32_t b : {0u, 2u, 3u, 4u, 6u, 8u}) {
      check_all_pairs(forest_gen::path(n), b);
      check_all_pairs(forest_gen::star(n), b);
      check_all_pairs(forest_gen::caterpillar(n, n / 3 + 1), b);
      check_all_pairs(forest_gen::binary(n), b);
      check_all_pairs(forest_gen::shuffle_labels(forest_gen::caterpillar(n, 4), rng), b);
    }
  }
}

TEST_CASE("lca equals naive on all pairs of small random trees", "[paper][lca]") {
  for (const std::uint64_t seed : {1u, 2u, 3u, 4u, 5u, 6u}) {
    std::mt19937_64 rng(seed);
    for (std::size_t n = 1; n <= 48; n += 1 + n / 4) {
      for (const std::uint32_t b : {0u, 2u, 3u, 5u, 7u}) {
        check_all_pairs(forest_gen::shuffle_labels(forest_gen::random_tree(n, rng), rng), b);
        check_all_pairs(forest_gen::shuffle_labels(forest_gen::random_deep_tree(n, 3, rng), rng),
                        b);
      }
    }
  }
}

TEST_CASE("lca test block size must lie in [2, 8]", "[paper][lca]") {
  const OrderedTree tree(forest_gen::path(5));
  REQUIRE_THROWS_AS(LcaTable(tree, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(LcaTable(tree, 9), std::invalid_argument);
}
