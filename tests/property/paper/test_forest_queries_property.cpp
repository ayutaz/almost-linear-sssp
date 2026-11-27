// Sec 3.4, Lemma 3.9: differential property tests of the paper's forest queries (Lemmas 3.7,
// 3.8 on the augmented tree) against the binary-lifting reference on larger random forests.
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <vector>

#include "paper/primitives/forest_queries.hpp"
#include "unit/paper/primitives/forest_gen.hpp"

using forest_gen::Parents;
using sssp::Vertex;
using sssp::paper::BinaryLiftingRef;
using sssp::paper::ForestQueries;

namespace {

void compare_random_queries(const Parents& parent, std::mt19937_64& rng, std::size_t queries) {
  const ForestQueries fq(parent);
  const BinaryLiftingRef ref(parent);
  const std::size_t n = parent.size();
  REQUIRE(fq.size() == n);
  for (Vertex v = 0; v < n; ++v) {
    REQUIRE(fq.depth(v) == ref.depth(v));
    REQUIRE(fq.component(v) == ref.component(v));
  }
  std::uniform_int_distribution<Vertex> pick(0, static_cast<Vertex>(n - 1));
  for (std::size_t q = 0; q < queries; ++q) {
    const Vertex u = pick(rng);
    const Vertex v = pick(rng);
    REQUIRE(fq.lca(u, v) == ref.lca(u, v));
    const std::uint32_t j = std::uniform_int_distribution<std::uint32_t>(0, ref.depth(u))(rng);
    REQUIRE(fq.level_ancestor(u, j) == ref.level_ancestor(u, j));
    // Sec 5.4 usage: the child of the LCA towards u is LA(u, depth(lca) + 1).
    const auto w = ref.lca(u, v);
    if (w.has_value() && *w != u) {
      const auto child = fq.level_ancestor(u, ref.depth(*w) + 1);
      REQUIRE(ref.level_ancestor(child, ref.depth(*w)) == *w);
    }
  }
}

}  // namespace

TEST_CASE("forest queries agree with binary lifting on random forests",
          "[paper][forest][property]") {
  for (const std::uint64_t seed : {101u, 102u, 103u, 104u, 105u, 106u}) {
    std::mt19937_64 rng(seed);
    for (const std::size_t n : {std::size_t{500}, std::size_t{5000}, std::size_t{30000}}) {
      for (const std::size_t trees : {std::size_t{1}, std::size_t{7}, n / 10}) {
        INFO("seed " << seed << " n " << n << " trees " << trees);
        compare_random_queries(
            forest_gen::shuffle_labels(forest_gen::random_forest(n, trees, rng), rng), rng, 4000);
      }
    }
  }
}

TEST_CASE("forest queries agree with binary lifting on deep and degenerate shapes",
          "[paper][forest][property]") {
  for (const std::uint64_t seed : {201u, 202u, 203u}) {
    std::mt19937_64 rng(seed);
    const std::size_t n = 40000;
    compare_random_queries(forest_gen::shuffle_labels(forest_gen::path(n), rng), rng, 4000);
    compare_random_queries(forest_gen::shuffle_labels(forest_gen::star(n), rng), rng, 4000);
    compare_random_queries(forest_gen::shuffle_labels(forest_gen::caterpillar(n, 300), rng), rng,
                           4000);
    compare_random_queries(forest_gen::shuffle_labels(forest_gen::binary(n), rng), rng, 4000);
    compare_random_queries(forest_gen::shuffle_labels(forest_gen::random_deep_tree(n, 4, rng), rng),
                           rng, 4000);
  }
}
