// Sec 2 (txt l.302-306): integer logarithms are tabulated by the halving recurrence for
// floor(log2 i) starting at i = 1; ceil_log2 is the bit count used by uniform sampling (l.308).
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "paper/model/ilog.hpp"

using sssp::paper::model::ceil_log2;
using sssp::paper::model::floor_log2;
using sssp::paper::model::IlogTable;

namespace {
constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();
}  // namespace

TEST_CASE("ilog table equals bit_width-1 for 1..2^20", "[paper][model][ilog]") {
  constexpr std::uint64_t kRange = std::uint64_t{1} << 20;
  const IlogTable table(kRange);
  REQUIRE(table.range() == kRange);
  for (std::uint64_t i = 1; i <= kRange; ++i) {
    REQUIRE(table.floor_log2(i) == static_cast<unsigned>(std::bit_width(i) - 1));
  }
}

TEST_CASE("ilog table ceil_log2 matches the defining inequality", "[paper][model][ilog]") {
  const IlogTable table(5000);
  REQUIRE(table.ceil_log2(1) == 0);
  for (std::uint64_t i = 2; i <= 5000; ++i) {
    const unsigned c = table.ceil_log2(i);
    // 2^(c-1) < i <= 2^c
    REQUIRE((std::uint64_t{1} << (c - 1)) < i);
    REQUIRE(i <= (std::uint64_t{1} << c));
  }
}

TEST_CASE("ilog table rejects indices outside 1..K", "[paper][model][ilog]") {
  const IlogTable table(16);
  REQUIRE_THROWS_AS(table.floor_log2(0), std::out_of_range);
  REQUIRE_THROWS_AS(table.floor_log2(17), std::out_of_range);
  REQUIRE(table.floor_log2(16) == 4);
  const IlogTable empty(0);
  REQUIRE(empty.range() == 0);
  REQUIRE_THROWS_AS(empty.floor_log2(1), std::out_of_range);
}

TEST_CASE("word floor_log2 and ceil_log2 match definitions at power-of-two boundaries",
          "[paper][model][ilog]") {
  REQUIRE(floor_log2(1) == 0);
  REQUIRE(ceil_log2(1) == 0);
  REQUIRE(ceil_log2(2) == 1);
  REQUIRE(ceil_log2(3) == 2);
  for (unsigned k = 1; k < 64; ++k) {
    const std::uint64_t p = std::uint64_t{1} << k;
    REQUIRE(floor_log2(p) == k);
    REQUIRE(floor_log2(p - 1) == k - 1);
    REQUIRE(floor_log2(p + 1) == k);
    REQUIRE(ceil_log2(p) == k);
    REQUIRE(ceil_log2(p + 1) == k + 1);
    if (p - 1 >= 1) REQUIRE(ceil_log2(p - 1) == (k == 1 ? 0 : k));
  }
  REQUIRE(floor_log2(kMax) == 63);
  REQUIRE(ceil_log2(kMax) == 64);
}

TEST_CASE("word logarithms are undefined at zero", "[paper][model][ilog]") {
  REQUIRE_THROWS_AS(floor_log2(0), std::domain_error);
  REQUIRE_THROWS_AS(ceil_log2(0), std::domain_error);
}

TEST_CASE("word logarithms agree with the table on its range", "[paper][model][ilog]") {
  const IlogTable table(1u << 12);
  for (std::uint64_t i = 1; i <= table.range(); ++i) {
    REQUIRE(floor_log2(i) == table.floor_log2(i));
    REQUIRE(ceil_log2(i) == table.ceil_log2(i));
  }
}

TEST_CASE("word logarithms are usable in constant expressions", "[paper][model][ilog]") {
  STATIC_REQUIRE(floor_log2(1024) == 10);
  STATIC_REQUIRE(ceil_log2(1025) == 11);
}
