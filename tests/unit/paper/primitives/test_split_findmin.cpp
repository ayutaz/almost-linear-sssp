// Unit tests for Sec 3.5 "Interval Minima Under Splits and Decreases" (Lemmas 3.10, 3.11).
//
// Every implementation is checked against an independent naive model (linear scans over an
// explicit boundary vector). Agreement of each implementation with the model implies that all
// implementations agree with each other (ref == dyadic == blocked).
//
// Interface decisions tested here (docs/research A8, A9):
//  - ties are resolved by (infinite flag, value, original index), i.e. +inf keys tie on index;
//  - interval IDs: the initial interval has ID 0; the s-th successful split (s = 1, 2, ...)
//    gives the left part ID 2s-1 and the right part ID 2s (both sides get fresh IDs);
//  - split(pos) returns false for an existing boundary (pos = 0, pos = k, or a previous cut)
//    and throws std::out_of_range for pos > k.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "paper/primitives/split_findmin.hpp"

namespace {

using sssp::paper::IntervalBounds;
using sssp::paper::IntervalId;
using sssp::paper::SfmMin;

template <class Key>
constexpr Key kInf = sssp::WeightTraits<Key>::infinity();

/// Independent oracle: explicit interval starts, linear-scan minima, same ID allocation rule.
template <class Key>
class NaiveModel {
 public:
  explicit NaiveModel(std::vector<Key> keys)
      : keys_(std::move(keys)), id_of_start_(keys_.size(), 0), is_start_(keys_.size(), 0) {
    is_start_[0] = 1;
    bounds_.push_back({0, static_cast<std::uint32_t>(keys_.size())});
  }

  [[nodiscard]] std::uint32_t size() const { return static_cast<std::uint32_t>(keys_.size()); }
  [[nodiscard]] std::uint32_t start_of(std::uint32_t i) const {
    while (is_start_[i] == 0) --i;
    return i;
  }
  [[nodiscard]] std::uint32_t end_of(std::uint32_t i) const {
    ++i;
    while (i < size() && is_start_[i] == 0) ++i;
    return i;
  }
  bool split(std::uint32_t pos) {
    if (pos == 0 || pos >= size() || is_start_[pos] != 0) return false;
    const std::uint32_t lo = start_of(pos);
    const std::uint32_t hi = end_of(pos);
    is_start_[pos] = 1;
    id_of_start_[lo] = static_cast<IntervalId>(bounds_.size());
    bounds_.push_back({lo, pos});
    id_of_start_[pos] = static_cast<IntervalId>(bounds_.size());
    bounds_.push_back({pos, hi});
    return true;
  }
  void decrease(std::uint32_t i, Key key) {
    if (key < keys_[i]) keys_[i] = key;
  }
  [[nodiscard]] SfmMin<Key> query(std::uint32_t i) const {
    const std::uint32_t lo = start_of(i);
    const std::uint32_t hi = end_of(i);
    SfmMin<Key> best{keys_[lo], lo};
    for (std::uint32_t j = lo + 1; j < hi; ++j) {
      if (keys_[j] < best.key) best = {keys_[j], j};  // strict: first (smallest) index wins
    }
    return best;
  }
  [[nodiscard]] IntervalId interval_of(std::uint32_t i) const { return id_of_start_[start_of(i)]; }
  [[nodiscard]] IntervalBounds bounds(IntervalId id) const { return bounds_[id]; }
  [[nodiscard]] Key key(std::uint32_t i) const { return keys_[i]; }

 private:
  std::vector<Key> keys_;
  std::vector<IntervalId> id_of_start_;
  std::vector<std::uint8_t> is_start_;
  std::vector<IntervalBounds> bounds_;
};

template <class S, class Key>
void require_matches(const S& s, const NaiveModel<Key>& model) {
  REQUIRE(s.size() == model.size());
  for (std::uint32_t i = 0; i < model.size(); ++i) {
    const auto got = s.query(i);
    const auto want = model.query(i);
    INFO("position " << i);
    REQUIRE(got.index == want.index);
    REQUIRE(got.key == want.key);
    REQUIRE(s.interval_of(i) == model.interval_of(i));
    const auto gb = s.bounds(s.interval_of(i));
    const auto wb = model.bounds(model.interval_of(i));
    REQUIRE(gb.lo == wb.lo);
    REQUIRE(gb.hi == wb.hi);
    REQUIRE(s.key(i) == model.key(i));
  }
}

template <class Key>
Key random_key(std::mt19937_64& rng, std::uint64_t range) {
  return static_cast<Key>(std::uniform_int_distribution<std::uint64_t>(0, range)(rng));
}

/// A split position biased towards block edges of the blocked structures (multiples of the
/// outer and inner block sizes, and their neighbours), plus uniform positions.
std::uint32_t random_split_pos(std::mt19937_64& rng, std::uint32_t k) {
  const std::uint32_t b1 = sssp::paper::sfm_block_size(k);
  const std::uint32_t b2 = sssp::paper::sfm_block_size(b1);
  std::uniform_int_distribution<std::uint32_t> any(0, k);
  switch (std::uniform_int_distribution<int>(0, 3)(rng)) {
    case 0: {
      const std::uint32_t base = (any(rng) / b1) * b1;
      const std::uint32_t pos = base + std::uniform_int_distribution<std::uint32_t>(0, 2)(rng);
      return pos == 0 ? 0 : std::min(k, pos - 1);
    }
    case 1: {
      const std::uint32_t base = (any(rng) / b1) * b1;
      const std::uint32_t inner =
          (std::uniform_int_distribution<std::uint32_t>(0, b1)(rng) / b2) * b2;
      return std::min(k, base + inner);
    }
    default:
      return any(rng);
  }
}

/// Random interleaving of splits, decreases and queries; the structure is compared with the
/// naive model after every operation.
template <class S>
void run_random_sequence(std::uint32_t k, std::uint64_t seed, bool start_infinite) {
  using Key = typename S::key_type;
  std::mt19937_64 rng(seed);
  const std::uint64_t key_range = 6;  // small range: many ties
  std::vector<Key> init(k, kInf<Key>);
  if (!start_infinite) {
    for (auto& x : init) x = random_key<Key>(rng, key_range + 4);
  }
  S s = start_infinite ? S(static_cast<std::size_t>(k)) : S(std::span<const Key>(init));
  NaiveModel<Key> model(init);
  require_matches(s, model);
  const std::uint32_t ops = 6 * k + 40;
  for (std::uint32_t step = 0; step < ops; ++step) {
    INFO("k=" << k << " seed=" << seed << " step=" << step);
    const int kind = std::uniform_int_distribution<int>(0, 9)(rng);
    if (kind < 3) {
      const std::uint32_t pos = random_split_pos(rng, k);
      REQUIRE(s.split(pos) == model.split(pos));
    } else if (kind < 8) {
      const std::uint32_t i = std::uniform_int_distribution<std::uint32_t>(0, k - 1)(rng);
      const Key cur = model.key(i);
      Key next = random_key<Key>(rng, key_range);
      if (cur < next) next = cur;  // decrease only (Sec 3.5: a decrease lowers a key)
      s.decrease(i, next);
      model.decrease(i, next);
    } else {
      const std::uint32_t i = std::uniform_int_distribution<std::uint32_t>(0, k - 1)(rng);
      const auto got = s.query(i);
      const auto want = model.query(i);
      REQUIRE(got.index == want.index);
      REQUIRE(got.key == want.key);
    }
    require_matches(s, model);
  }
  // Finish by cutting everything: every interval becomes a singleton.
  for (std::uint32_t pos = 1; pos < k; ++pos) REQUIRE(s.split(pos) == model.split(pos));
  require_matches(s, model);
}

}  // namespace

#define SFM_DOUBLE_IMPLEMENTATIONS                                               \
  sssp::paper::SplitFindminRef<double>, sssp::paper::DyadicSplitFindmin<double>, \
      sssp::paper::OneLevelBlockedSplitFindmin<double>, sssp::paper::SplitFindmin<double>
#define SFM_IMPLEMENTATIONS                                                                    \
  sssp::paper::SplitFindminRef<std::uint64_t>, sssp::paper::DyadicSplitFindmin<std::uint64_t>, \
      sssp::paper::OneLevelBlockedSplitFindmin<std::uint64_t>,                                 \
      sssp::paper::SplitFindmin<std::uint64_t>, SFM_DOUBLE_IMPLEMENTATIONS

TEMPLATE_TEST_CASE("split-findmin: k=1 is one interval that cannot be split", "[split_findmin]",
                   SFM_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  TestType s(std::size_t{1});
  REQUIRE(s.size() == 1);
  REQUIRE(s.query(0).index == 0);
  REQUIRE(s.query(0).key == kInf<Key>);
  REQUIRE(s.interval_of(0) == 0);
  REQUIRE(s.bounds(0).lo == 0);
  REQUIRE(s.bounds(0).hi == 1);
  REQUIRE_FALSE(s.split(0));
  REQUIRE_FALSE(s.split(1));
  REQUIRE_THROWS_AS(s.split(2), std::out_of_range);
  s.decrease(0, Key{5});
  REQUIRE(s.query(0).key == Key{5});
  REQUIRE(s.query(0).index == 0);
  REQUIRE(s.key(0) == Key{5});
}

TEMPLATE_TEST_CASE("split-findmin: zero records are rejected", "[split_findmin]",
                   SFM_IMPLEMENTATIONS) {
  REQUIRE_THROWS_AS(TestType(std::size_t{0}), std::invalid_argument);
}

TEMPLATE_TEST_CASE("split-findmin: infinite keys return the smallest index of the interval",
                   "[split_findmin]", SFM_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  TestType s(std::size_t{10});
  REQUIRE(s.split(3));
  REQUIRE(s.split(7));
  REQUIRE(s.query(0).index == 0);
  REQUIRE(s.query(2).index == 0);
  REQUIRE(s.query(5).index == 3);
  REQUIRE(s.query(9).index == 7);
  REQUIRE(s.query(9).key == kInf<Key>);
  // One finite key beats every infinite key of its interval only.
  s.decrease(8, Key{100});
  REQUIRE(s.query(7).index == 8);
  REQUIRE(s.query(7).key == Key{100});
  REQUIRE(s.query(4).index == 3);
  REQUIRE(s.query(4).key == kInf<Key>);
}

TEMPLATE_TEST_CASE("split-findmin: equal keys are resolved by original index", "[split_findmin]",
                   SFM_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  const std::vector<Key> keys{9, 4, 7, 4, 4, 9};
  TestType s{std::span<const Key>(keys)};
  REQUIRE(s.query(5).index == 1);
  REQUIRE(s.query(5).key == Key{4});
  REQUIRE(s.split(2));  // [0,2) [2,6)
  REQUIRE(s.query(0).index == 1);
  REQUIRE(s.query(2).index == 3);
  s.decrease(5, Key{4});  // ties with 3 and 4; index 3 still wins
  REQUIRE(s.query(5).index == 3);
  s.decrease(2, Key{4});  // now index 2 wins the tie
  REQUIRE(s.query(5).index == 2);
  REQUIRE(s.split(5));  // [2,5) [5,6)
  REQUIRE(s.query(5).index == 5);
  REQUIRE(s.query(4).index == 2);
}

TEMPLATE_TEST_CASE("split-findmin: split at an existing boundary is rejected and changes nothing",
                   "[split_findmin]", SFM_IMPLEMENTATIONS) {
  TestType s(std::size_t{8});
  REQUIRE_FALSE(s.split(0));
  REQUIRE_FALSE(s.split(8));
  REQUIRE_THROWS_AS(s.split(9), std::out_of_range);
  REQUIRE(s.interval_of(5) == 0);
  REQUIRE(s.split(4));
  const auto left = s.interval_of(0);
  const auto right = s.interval_of(7);
  REQUIRE_FALSE(s.split(4));
  REQUIRE(s.interval_of(0) == left);
  REQUIRE(s.interval_of(7) == right);
  REQUIRE(s.bounds(left).hi == 4);
  REQUIRE(s.bounds(right).lo == 4);
}

TEMPLATE_TEST_CASE("split-findmin: both sides of a split receive fresh interval IDs",
                   "[split_findmin]", SFM_IMPLEMENTATIONS) {
  TestType s(std::size_t{12});
  REQUIRE(s.interval_of(11) == 0);
  REQUIRE(s.split(5));
  REQUIRE(s.interval_of(0) == 1);
  REQUIRE(s.interval_of(4) == 1);
  REQUIRE(s.interval_of(5) == 2);
  REQUIRE(s.split(9));  // splits [5,12)
  REQUIRE(s.interval_of(0) == 1);
  REQUIRE(s.interval_of(5) == 3);
  REQUIRE(s.interval_of(11) == 4);
  REQUIRE(s.bounds(3).lo == 5);
  REQUIRE(s.bounds(3).hi == 9);
  REQUIRE(s.bounds(4).lo == 9);
  REQUIRE(s.bounds(4).hi == 12);
  REQUIRE(s.split(1));  // splits [0,5)
  REQUIRE(s.interval_of(0) == 5);
  REQUIRE(s.interval_of(1) == 6);
  REQUIRE(s.bounds(6).hi == 5);
}

TEMPLATE_TEST_CASE(
    "split-findmin: exhaustive random operation sequences for k=1..70 match the "
    "naive model",
    "[split_findmin]", SFM_IMPLEMENTATIONS) {
  for (std::uint32_t k = 1; k <= 70; ++k) {
    for (std::uint64_t seed = 1; seed <= 3; ++seed) {
      run_random_sequence<TestType>(k, 1000 * k + seed, /*start_infinite=*/seed == 1);
    }
  }
}

TEMPLATE_TEST_CASE("split-findmin: longer random sequences on larger arrays match the naive model",
                   "[split_findmin]", SFM_IMPLEMENTATIONS) {
  for (const std::uint32_t k : {127U, 128U, 129U, 257U, 600U}) {
    run_random_sequence<TestType>(k, 77 + k, /*start_infinite=*/k % 2 == 0);
  }
}

TEMPLATE_TEST_CASE("split-findmin: double keys reject NaN", "[split_findmin]",
                   SFM_DOUBLE_IMPLEMENTATIONS) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> bad{1.0, nan, 2.0};
  REQUIRE_THROWS_AS(TestType{std::span<const double>(bad)}, std::invalid_argument);
  TestType s(std::size_t{3});
  REQUIRE_THROWS_AS(s.decrease(1, nan), std::invalid_argument);
  REQUIRE(s.query(1).key == kInf<double>);
}

TEST_CASE("split-findmin: block size follows Lemma 3.11", "[split_findmin]") {
  // B = max(2, ceil(log2(k+1))).
  REQUIRE(sssp::paper::sfm_block_size(1) == 2);
  REQUIRE(sssp::paper::sfm_block_size(3) == 2);
  REQUIRE(sssp::paper::sfm_block_size(4) == 3);
  REQUIRE(sssp::paper::sfm_block_size(7) == 3);
  REQUIRE(sssp::paper::sfm_block_size(8) == 4);
  REQUIRE(sssp::paper::sfm_block_size(70) == 7);
  REQUIRE(sssp::paper::sfm_block_size((std::uint64_t{1} << 32) - 1) == 32);
  REQUIRE(sssp::paper::sfm_block_size(std::uint64_t{1} << 32) == 33);
}

TEMPLATE_TEST_CASE("split-findmin: block-edge layouts of Lemma 3.11 give correct minima",
                   "[split_findmin]", SFM_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  // k = 64: B = 7, blocks [0,7) [7,14) ... [56,63) and the short final block [63,64).
  const std::uint32_t k = 64;
  std::vector<Key> keys(k);
  for (std::uint32_t i = 0; i < k; ++i) keys[i] = static_cast<Key>(1000 - 3 * i);
  TestType s{std::span<const Key>(keys)};
  NaiveModel<Key> model(keys);
  const auto cut = [&](std::uint32_t pos) {
    REQUIRE(s.split(pos) == model.split(pos));
    require_matches(s, model);
  };
  cut(7);   // between blocks
  cut(14);  // [7,14) is exactly one whole block
  cut(17);  // internal cut: [14,17) proper fragment
  cut(26);  // [17,26): fragment of block 2 + fragment of block 3, no whole block between
  cut(30);  // [26,30) spans two partial blocks
  cut(63);  // the short final block alone
  cut(45);  // internal cut inside a run of whole blocks [30,63)
  for (std::uint32_t i = 0; i < k; i += 5) {
    s.decrease(i, static_cast<Key>(i % 3));
    model.decrease(i, static_cast<Key>(i % 3));
    require_matches(s, model);
  }
  cut(1);
  cut(62);
  cut(15);
  require_matches(s, model);
}

TEST_CASE("split-findmin: two-level blocking uses B' = max(2, ceil(log2(B+1))) inside blocks",
          "[split_findmin]") {
  // Reading A7 of the research notes: the inner block size is derived from the outer block size
  // B, for every lower structure including a short final block.
  sssp::paper::SplitFindmin<std::uint64_t> s(std::size_t{100});
  REQUIRE(s.block_size() == 7);  // ceil(log2 101) = 7
  REQUIRE(s.num_blocks() == 15);
  for (std::uint32_t b = 0; b < s.num_blocks(); ++b) {
    REQUIRE(s.lower(b).block_size() == 3);  // max(2, ceil(log2 8)) = 3
  }
  REQUIRE(s.lower(14).size() == 2);  // final short block [98,100)
  sssp::paper::OneLevelBlockedSplitFindmin<std::uint64_t> one(std::size_t{100});
  REQUIRE(one.block_size() == 7);
  REQUIRE(one.lower(14).size() == 2);
}
