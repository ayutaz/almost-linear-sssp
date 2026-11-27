#include <catch2/catch_test_macros.hpp>
#include <string_view>

#include "sssp/version.hpp"

TEST_CASE("library exposes version and target paper", "[version]") {
  REQUIRE(sssp::version_string() == std::string_view{"0.1.0"});
  REQUIRE(sssp::paper_arxiv_id() == std::string_view{"2611.34911v1"});
}
