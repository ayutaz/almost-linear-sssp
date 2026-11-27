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

// ---------------------------------------------------------------- Lemma 3.9 (forests)

namespace {

template <class Q>
void check_forest_all_pairs(const Q& q, const Parents& parent) {
  REQUIRE(q.size() == parent.size());
  for (Vertex u = 0; u < parent.size(); ++u) {
    REQUIRE(q.component(u) == forest_gen::naive_root(parent, u));
    REQUIRE(q.depth(u) == forest_gen::naive_depth(parent, u));
    for (Vertex v = 0; v < parent.size(); ++v) {
      INFO("u " << u << " v " << v);
      REQUIRE(q.lca(u, v) == forest_gen::naive_lca(parent, u, v));
    }
  }
}

}  // namespace

TEST_CASE("forest lca returns nullopt for vertices in different trees", "[paper][lca][forest]") {
  // Trees {0 -> 1, 0 -> 2}, {3 -> 4} and the isolated vertex 5.
  const Parents parent{kNoVertex, 0, 0, kNoVertex, 3, kNoVertex};
  const sssp::paper::ForestQueries fq(parent);
  const sssp::paper::BinaryLiftingRef ref(parent);
  for (const Vertex a : {0u, 1u, 2u}) {
    for (const Vertex b : {3u, 4u, 5u}) {
      REQUIRE_FALSE(fq.lca(a, b).has_value());
      REQUIRE_FALSE(fq.lca(b, a).has_value());
      REQUIRE_FALSE(ref.lca(a, b).has_value());
    }
  }
  REQUIRE(fq.lca(1, 2) == Vertex{0});
  REQUIRE(fq.lca(4, 3) == Vertex{3});
  REQUIRE(fq.lca(5, 5) == Vertex{5});
  REQUIRE(fq.component(4) == 3);
  REQUIRE(fq.component(5) == 5);
  REQUIRE(ref.lca(1, 2) == Vertex{0});
}

TEST_CASE("forest queries on the empty forest and on isolated vertices", "[paper][lca][forest]") {
  const sssp::paper::ForestQueries empty(Parents{});
  REQUIRE(empty.size() == 0);
  const Parents isolated(4, kNoVertex);
  check_forest_all_pairs(sssp::paper::ForestQueries(isolated), isolated);
  check_forest_all_pairs(sssp::paper::BinaryLiftingRef(isolated), isolated);
}

TEST_CASE("forest queries reject invalid parent arrays", "[paper][lca][forest]") {
  using sssp::paper::BinaryLiftingRef;
  using sssp::paper::ForestQueries;
  REQUIRE_THROWS_AS(ForestQueries(Parents{kNoVertex, 2}), std::out_of_range);  // 2 == n
  REQUIRE_THROWS_AS(ForestQueries(Parents{kNoVertex, 7}), std::out_of_range);
  REQUIRE_THROWS_AS(ForestQueries(Parents{kNoVertex, 2, 1}), std::invalid_argument);  // cycle
  REQUIRE_THROWS_AS(ForestQueries(Parents{0}), std::invalid_argument);                // self
  REQUIRE_THROWS_AS(BinaryLiftingRef(Parents{kNoVertex, 2}), std::out_of_range);
  REQUIRE_THROWS_AS(BinaryLiftingRef(Parents{kNoVertex, 2, 1}), std::invalid_argument);
  REQUIRE_THROWS_AS(BinaryLiftingRef(Parents{0}), std::invalid_argument);
}

TEST_CASE("forest lca equals naive on all pairs of small random forests", "[paper][lca][forest]") {
  for (const std::uint64_t seed : {31u, 32u, 33u, 34u, 35u}) {
    std::mt19937_64 rng(seed);
    for (std::size_t n = 1; n <= 45; n += 1 + n / 5) {
      for (const std::size_t trees : {std::size_t{1}, std::size_t{2}, n / 3 + 1, n}) {
        const auto parent =
            forest_gen::shuffle_labels(forest_gen::random_forest(n, trees, rng), rng);
        check_forest_all_pairs(sssp::paper::ForestQueries(parent), parent);
        check_forest_all_pairs(sssp::paper::ForestQueries(parent, 4), parent);
        check_forest_all_pairs(sssp::paper::BinaryLiftingRef(parent), parent);
      }
    }
  }
}
