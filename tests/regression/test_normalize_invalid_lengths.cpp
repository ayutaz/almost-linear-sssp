// Regression 20261127-normalize-invalid-lengths (record:
// cases/20261127-normalize-invalid-lengths.json), findings F-2, F-3 and A-4 of
// docs/reviews/P4-REVIEW-hierarchy.md.
//
// Def 4.1 (Sec 4, l.646-658) admits only nonnegative real arc lengths and links of length at most
// c0 h(C); infinity is a separate flag, not a length (Sec 2, l.292-293), so +inf (double) and the
// uint64 infinity sentinel 2^64 - 1 are not lengths either. normalize.hpp documents that invalid
// lengths raise std::invalid_argument. Expected values below come from that contract, not from a
// run of the implementation:
//   * normalize() throws std::invalid_argument for an arc (between two roots, a self-arc, or
//     below an LCA) or a link whose length is negative, NaN, +inf or the sentinel, and for a link
//     longer than c0 h(C) (H5);
//   * evaluate(), whose first step (Alg 3 l.2) is normalize(), throws the same exception and
//     returns promptly;
//   * valid lengths at the boundaries (-0.0, a link of exactly c0 h(C), masses that saturate)
//     are still accepted.
//
// Before the fix (review findings):
//   F-2  two singleton roots with a root arc of length -1.0: evaluate returned value 4 at the
//        head; NaN was ignored silently; +inf and the sentinel gave std::overflow_error (A-4).
//   F-3  one root of height 1 over two leaves, link (0, 1) of length NaN, +inf, 2^64 - 1 or 2^40
//        (> c0 h = 1), state 1 never offered: evaluate did not return within 10 s (Alg 3 l.10
//        walks width-1 slots up to the saturated mass). -3.0 gave std::logic_error.
//
// A hang must not stall ctest for its 600 s preset timeout: every call runs under a watchdog
// that ends the process with a failure after kDeadline.

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "sssp/graph.hpp"

namespace {

using sssp::paper::NodeId;
using sssp::paper::OrderedHierarchy;
using sssp::paper::SeededOrderedHierarchy;
using sssp::paper::StateGraph;
using U = std::uint64_t;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPosInf = std::numeric_limits<double>::infinity();
constexpr U kSentinel = sssp::WeightTraits<U>::infinity();
// Valid inputs of this size take microseconds, also under the sanitizers.
constexpr auto kDeadline = std::chrono::seconds(5);

/// Runs fn on this thread; if it has not returned after kDeadline, a watchdog prints the case
/// and ends the process with a failure (a hung call cannot be cancelled).
template <class Fn>
void within_deadline(const std::string& what, Fn&& fn) {
  std::mutex m;
  std::condition_variable cv;
  bool done = false;
  std::thread watchdog([&] {
    std::unique_lock lock(m);
    if (!cv.wait_for(lock, kDeadline, [&] { return done; })) {
      std::fprintf(stderr, "TIMEOUT: %s did not return within %lld s\n", what.c_str(),
                   static_cast<long long>(kDeadline.count()));
      std::fflush(stderr);
      std::_Exit(EXIT_FAILURE);
    }
  });
  const auto finish = [&] {
    {
      const std::lock_guard lock(m);
      done = true;
    }
    cv.notify_one();
    watchdog.join();
  };
  try {
    std::forward<Fn>(fn)();
  } catch (...) {
    finish();
    throw;
  }
  finish();
}

/// F-2 shape: two singleton roots, arc 0 -> 1 of length len (root quotient), offer 5 at state 0.
template <class W>
SeededOrderedHierarchy<W> two_roots(W len) {
  StateGraph<W> g;
  g.num_states = 2;
  g.arcs.push_back({0, 1, len});
  SeededOrderedHierarchy<W> s;
  s.hierarchy = sssp::paper::make_leaf_forest<W>(g, {0, 1});
  sssp::paper::collect_roots(s.hierarchy);
  s.sources.push_back({0, W{0}});
  s.offers.push_back({0, W{5}});
  return s;
}

/// F-3 shape: one root of height h over two leaves with link (0, 1) of length link; optional arc
/// 0 -> 1 of length arc; offer 0 at state 0 (without the arc, state 1 is never offered).
template <class W>
SeededOrderedHierarchy<W> pair_node(W h, W link, std::uint32_t c0 = 1, bool with_arc = false,
                                    W arc = W{0}) {
  StateGraph<W> g;
  g.num_states = 2;
  if (with_arc) g.arcs.push_back({0, 1, arc});
  SeededOrderedHierarchy<W> s;
  s.hierarchy = sssp::paper::make_leaf_forest<W>(g, {0, 1});
  s.hierarchy.c0 = c0;
  sssp::paper::add_internal_node<W>(s.hierarchy, h, {0, 1}, {{0, 1, link, sssp::paper::kNoArc}});
  sssp::paper::collect_roots(s.hierarchy);
  s.sources.push_back({0, W{0}});
  s.offers.push_back({0, W{0}});
  return s;
}

/// Self-arc 0 -> 0 of length len on two singleton roots; offer 0 at state 0.
template <class W>
SeededOrderedHierarchy<W> self_arc(W len) {
  StateGraph<W> g;
  g.num_states = 2;
  g.arcs.push_back({0, 0, len});
  SeededOrderedHierarchy<W> s;
  s.hierarchy = sssp::paper::make_leaf_forest<W>(g, {0, 1});
  sssp::paper::collect_roots(s.hierarchy);
  s.sources.push_back({0, W{0}});
  s.offers.push_back({0, W{0}});
  return s;
}

/// A non-root node with an invalid link: X = {1, 2} (height 1, link (1, 2) of length link) below
/// the root {0, X} (height 4, link (0, 1) of length 1). Offer 0 at state 0.
template <class W>
SeededOrderedHierarchy<W> nested(W link) {
  StateGraph<W> g;
  g.num_states = 3;
  SeededOrderedHierarchy<W> s;
  s.hierarchy = sssp::paper::make_leaf_forest<W>(g, {0, 1, 2});
  const NodeId x = sssp::paper::add_internal_node<W>(s.hierarchy, W{1}, {1, 2},
                                                     {{1, 2, link, sssp::paper::kNoArc}});
  sssp::paper::add_internal_node<W>(s.hierarchy, W{4}, {0, x}, {{0, 1, W{1}, sssp::paper::kNoArc}});
  sssp::paper::collect_roots(s.hierarchy);
  s.sources.push_back({0, W{0}});
  s.offers.push_back({0, W{0}});
  return s;
}

/// normalize() and evaluate() both reject the input with std::invalid_argument, in time.
template <class W>
void require_rejected(const std::string& what, const SeededOrderedHierarchy<W>& s) {
  INFO(what);
  within_deadline("normalize: " + what, [&] {
    CHECK_THROWS_AS(sssp::paper::normalize(s.hierarchy), std::invalid_argument);
  });
  within_deadline("evaluate: " + what,
                  [&] { CHECK_THROWS_AS(sssp::paper::evaluate(s), std::invalid_argument); });
}

/// Valid input: normalize() and evaluate() return, in time.
template <class W>
sssp::paper::EvalResult<W> require_accepted(const std::string& what,
                                            const SeededOrderedHierarchy<W>& s) {
  INFO(what);
  sssp::paper::EvalResult<W> r;
  within_deadline("normalize: " + what,
                  [&] { CHECK_NOTHROW(sssp::paper::normalize(s.hierarchy)); });
  within_deadline("evaluate: " + what, [&] { REQUIRE_NOTHROW(r = sssp::paper::evaluate(s)); });
  return r;
}

}  // namespace

TEST_CASE("20261127-normalize-invalid-lengths: F-2 root arcs with invalid lengths are rejected",
          "[regression][paper]") {
  require_rejected("double root arc -1.0", two_roots<double>(-1.0));
  require_rejected("double root arc NaN", two_roots<double>(kNaN));
  require_rejected("double root arc +inf", two_roots<double>(kPosInf));
  require_rejected("uint64 root arc 2^64 - 1 (sentinel)", two_roots<U>(kSentinel));
}

TEST_CASE("20261127-normalize-invalid-lengths: self-arcs with invalid lengths are rejected",
          "[regression][paper]") {
  // Self-arcs are discarded (Sec 4.1 l.719), but D still has nonnegative lengths only (Def 4.1).
  require_rejected("double self-arc -1.0", self_arc<double>(-1.0));
  require_rejected("double self-arc NaN", self_arc<double>(kNaN));
  require_rejected("double self-arc +inf", self_arc<double>(kPosInf));
  require_rejected("uint64 self-arc 2^64 - 1 (sentinel)", self_arc<U>(kSentinel));
}

TEST_CASE("20261127-normalize-invalid-lengths: arcs below an LCA with invalid lengths are rejected",
          "[regression][paper]") {
  require_rejected("double arc below LCA -1.0", pair_node<double>(4.0, 1.0, 1, true, -1.0));
  require_rejected("double arc below LCA NaN", pair_node<double>(4.0, 1.0, 1, true, kNaN));
  require_rejected("double arc below LCA +inf", pair_node<double>(4.0, 1.0, 1, true, kPosInf));
  require_rejected("uint64 arc below LCA 2^64 - 1", pair_node<U>(4, 1, 1, true, kSentinel));
}

TEST_CASE("20261127-normalize-invalid-lengths: F-3 root links with invalid lengths are rejected",
          "[regression][paper]") {
  require_rejected("double root link NaN", pair_node<double>(1.0, kNaN));
  require_rejected("double root link +inf", pair_node<double>(1.0, kPosInf));
  require_rejected("double root link -3.0", pair_node<double>(1.0, -3.0));
  require_rejected("double root link -0.5, arc 0 -> 1 of length 1",
                   pair_node<double>(1.0, -0.5, 1, true, 1.0));
  require_rejected("uint64 root link 2^64 - 1 (sentinel)", pair_node<U>(1, kSentinel));
}

TEST_CASE("20261127-normalize-invalid-lengths: F-3 links longer than c0 h(C) are rejected (H5)",
          "[regression][paper]") {
  require_rejected("uint64 h = 1, c0 = 1, link 2^40", pair_node<U>(1, U{1} << 40));
  require_rejected("uint64 h = 1, c0 = 1, link 2", pair_node<U>(1, 2));
  require_rejected("uint64 h = 1, c0 = 2, link 3", pair_node<U>(1, 3, 2));
  require_rejected("double h = 1, c0 = 1, link 2^40", pair_node<double>(1.0, 0x1p40));
  require_rejected("double h = 1, c0 = 3, link 3 + 2^-50",
                   pair_node<double>(1.0, 3.0 + 0x1p-50, 3));
}

TEST_CASE("20261127-normalize-invalid-lengths: invalid links below the root are rejected",
          "[regression][paper]") {
  require_rejected("double nested link NaN", nested<double>(kNaN));
  require_rejected("double nested link -1.0", nested<double>(-1.0));
  require_rejected("double nested link +inf", nested<double>(kPosInf));
  require_rejected("uint64 nested link 2^64 - 1", nested<U>(kSentinel));
  require_rejected("uint64 nested link 2 > c0 h = 1", nested<U>(2));
}

TEST_CASE("20261127-normalize-invalid-lengths: valid boundary lengths are still accepted",
          "[regression][paper]") {
  // -0.0 is the nonnegative real 0 (model::Len maps it to +0.0).
  const auto r1 = require_accepted("double root arc -0.0", two_roots<double>(-0.0));
  CHECK(r1.value[1] == 5.0);
  const auto r2 = require_accepted("double self-arc -0.0", self_arc<double>(-0.0));
  CHECK(r2.value[0] == 0.0);
  // A link of exactly c0 h(C) satisfies H5; the arc 0 -> 1 makes state 1 reachable.
  const auto r3 = require_accepted("uint64 h = 1, c0 = 1, link 1", pair_node<U>(1, 1, 1, true, 1));
  CHECK(r3.value[1] == 1);
  const auto r4 = require_accepted("uint64 h = 1, c0 = 2, link 2", pair_node<U>(1, 2, 2, true, 2));
  CHECK(r4.value[1] == 2);
  const auto r5 =
      require_accepted("double h = 1, c0 = 3, link 3", pair_node<double>(1.0, 3.0, 3, true, 3.0));
  CHECK(r5.value[1] == 3.0);
  // Root arcs of any finite length are valid (Def 4.1 has no upper bound on arcs).
  const auto r6 = require_accepted("uint64 root arc 2^64 - 7", two_roots<U>(kSentinel - 6));
  CHECK(r6.value[1] == kSentinel - 1);
  // c0 h(C) beyond the range of W bounds every finite link; the mass saturates (D-23).
  within_deadline("normalize: uint64 h = 2^63, c0 = 2, link 2^64 - 2", [] {
    CHECK_NOTHROW(sssp::paper::normalize(pair_node<U>(U{1} << 63, kSentinel - 1, 2).hierarchy));
  });
  within_deadline("normalize: double h = 2^1023, c0 = 2, link DBL_MAX", [] {
    CHECK_NOTHROW(sssp::paper::normalize(
        pair_node<double>(0x1p1023, std::numeric_limits<double>::max(), 2).hierarchy));
  });
}

TEST_CASE("20261127-normalize-invalid-lengths: the H5 check costs O(log c0), exact at the bound",
          "[regression][paper]") {
  // normalize.hpp promises O(M + N); c0 is an absolute constant of Def 4.1 (l.657), but
  // normalize() accepts every c0 <= 2^31 - 1. Forming c0 h by c0 - 1 additions would take 2^31
  // additions per node here; the binary expansion of c0 (Sec 2 l.297 doublings) takes 60.
  constexpr std::uint32_t kMaxC0 = std::numeric_limits<std::uint32_t>::max() / 2;
  const auto r =
      require_accepted("uint64 h = 1, c0 = 2^31 - 1, link 1", pair_node<U>(1, 1, kMaxC0, true, 1));
  CHECK(r.value[1] == 1);
  within_deadline("normalize: uint64 h = 1, c0 = 2^31 - 1, link 2^31 - 1", [] {
    CHECK_NOTHROW(sssp::paper::normalize(pair_node<U>(1, kMaxC0, kMaxC0).hierarchy));
  });
  within_deadline("normalize: double h = 1, c0 = 2^31 - 1, link 2^31 - 1", [] {
    CHECK_NOTHROW(sssp::paper::normalize(
        pair_node<double>(1.0, static_cast<double>(kMaxC0), kMaxC0).hierarchy));
  });
  require_rejected("uint64 h = 1, c0 = 2^31 - 1, link 2^31", pair_node<U>(1, U{1} << 31, kMaxC0));
  require_rejected("double h = 1, c0 = 2^31 - 1, link 2^31",
                   pair_node<double>(1.0, 0x1p31, kMaxC0));
}
