// Minimal JSON value used by the benchmark runner: --env-json is copied into every record
// verbatim (numbers keep their original token), and records are written as one JSON line.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

#include "sssp/bench_runner.hpp"

using sssp::bench::Json;
using sssp::bench::JsonError;

TEST_CASE("bench json: parse and dump round trip keeps order and number tokens", "[bench]") {
  const std::string text =
      R"({"b":1,"a":[true,false,null],"x":1.50,"e":-2E+3,"s":"t\"q\\\/\n","o":{}})";
  const Json j = Json::parse(text);
  REQUIRE(j.kind() == Json::Kind::kObject);
  REQUIRE(j.items().size() == 6);
  CHECK(j.items()[0].first == "b");
  CHECK(j.items()[1].first == "a");
  CHECK(j.find("x")->number_text() == "1.50");
  CHECK(j.find("e")->number_text() == "-2E+3");
  CHECK(j.find("s")->as_string() == "t\"q\\/\n");
  CHECK(j.find("a")->elements().size() == 3);
  CHECK(j.find("a")->elements()[0].as_bool());
  CHECK(j.find("a")->elements()[2].is_null());
  CHECK(j.find("missing") == nullptr);
  CHECK(j.dump() == R"({"b":1,"a":[true,false,null],"x":1.50,"e":-2E+3,"s":"t\"q\\/\n","o":{}})");
}

TEST_CASE("bench json: whitespace and unicode escapes", "[bench]") {
  const Json j = Json::parse(" { \"k\" :\t[ 1 , 2 ]\n, \"u\": \"\\u00e9\\ud83d\\ude00\" } ");
  CHECK(j.find("k")->elements().size() == 2);
  CHECK(j.find("u")->as_string() == "\xC3\xA9\xF0\x9F\x98\x80");
  // Control characters are escaped on output; other UTF-8 passes through unchanged.
  CHECK(Json::string("a\x01"
                     "b\xC3\xA9\t")
            .dump() == "\"a\\u0001b\xC3\xA9\\t\"");
}

TEST_CASE("bench json: duplicate keys keep the last value (Python json semantics)", "[bench]") {
  const Json j = Json::parse(R"({"k":1,"k":2})");
  REQUIRE(j.items().size() == 1);
  CHECK(j.find("k")->number_text() == "2");
}

TEST_CASE("bench json: malformed input is rejected", "[bench]") {
  for (const char* bad :
       {"",      "{",       "[1,]",        "{\"a\":1,}", "01",    "1.",       "-",
        "+1",    ".5",      "1e",          "nan",        "NaN",   "Infinity", "tru",
        "\"abc", "\"\\x\"", "\"\\ud800\"", "{1:2}",      "[1] 2", "\"a\nb\"", "{\"a\" 1}"}) {
    INFO(bad);
    CHECK_THROWS_AS(Json::parse(bad), JsonError);
  }
}

TEST_CASE("bench json: builders and numeric formatting", "[bench]") {
  CHECK(Json::integer(std::uint64_t{18446744073709551615ULL}).dump() == "18446744073709551615");
  CHECK(Json::integer(std::int64_t{-5}).dump() == "-5");
  CHECK(Json::real(4.0).dump() == "4");
  CHECK(Json::real(0.1).dump() == "0.1");
  CHECK(Json::real(2.5).dump() == "2.5");
  CHECK_THROWS_AS(Json::real(std::numeric_limits<double>::quiet_NaN()), JsonError);
  CHECK_THROWS_AS(Json::real(std::numeric_limits<double>::infinity()), JsonError);
  CHECK(Json::null().dump() == "null");
  CHECK(Json::boolean(true).dump() == "true");

  Json o = Json::object();
  o.set("z", Json::integer(std::uint64_t{1}));
  o.set("a", Json::array());
  o.set("z", Json::string("replaced"));  // replaces in place, keeps position
  o.find_mut("a")->push_back(Json::null());
  CHECK(o.dump() == R"({"z":"replaced","a":[null]})");
  // Parsed integer tokens are recognised as integers.
  CHECK(Json::parse("12").is_integer());
  CHECK_FALSE(Json::parse("12.0").is_integer());
  CHECK_FALSE(Json::parse("1e2").is_integer());
}
