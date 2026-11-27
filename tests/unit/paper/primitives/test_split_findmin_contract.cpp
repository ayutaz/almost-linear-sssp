// Caller-facing contract of split-findmin (Sec 3.5, l.577-639) as used by Evaluate
// (Sec 4.3 Alg 2 l.841-842, l.859-862; proof of Lemma 4.7 l.939-947).
//
// Interval identifiers (open question Q-14, decided 2026-11-27): a split gives BOTH sides fresh
// identifiers and the old identifier becomes invalid; an identifier stays valid until the next
// split of its own interval. Numbering (header contract): the initial interval is 0, the s-th
// successful split gives 2s-1 to the left part and 2s to the right part, so identifiers are
// never reused. What Evaluate relies on (l.939-941, l.946-947: "stores the interval identifier
// at its chosen minimum", "a pointer to a sample leaf in its interval"):
//  * a stored (sample record, identifier) pair is current iff interval_of(sample) still returns
//    that identifier, and that happens iff no successful split has cut that interval since;
//  * a retired identifier is never returned again by interval_of for any record;
//  * decreases, queries and rejected splits never change any identifier;
//  * the minimum returned by query(i) lies in the interval of i.
// Everything is checked against bookkeeping written here (explicit boundary flags, linear
// scans), independent of the implementations and of the model in test_split_findmin.cpp.
//
// Key increases (open question Q-16, decision revised by the lead on 2026-11-27): Sec 3.5 only
// "lowers a key" (l.579), so decrease(i, key) with key above the current key is outside the
// interface. It throws std::invalid_argument in every build (with and without NDEBUG) and leaves
// the structure unchanged, so a caller that forgets to check settled/retired records before an
// update (Evaluate, l.862) fails loudly instead of keeping a stale minimum. Equal keys are not an
// increase and are accepted everywhere. The earlier decision (assert with assertions, silent
// no-op under NDEBUG) is superseded.
//
// Record indices (same revision, a [Supplement]: input checks are not paper operations):
// decrease, query, key and interval_of with i >= k throw std::out_of_range in every build and
// change nothing, like split(pos) with pos > k and bounds() of an identifier never issued.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
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

template <class Key>
constexpr Key kInf = sssp::WeightTraits<Key>::infinity();

// Own bookkeeping of the expected interval structure and identifiers.
struct IdBook {
  std::uint32_t k;
  std::vector<bool> is_start;          // is_start[p]: p is the first record of an interval
  std::vector<IntervalId> id;          // expected identifier of every record
  std::vector<IntervalBounds> issued;  // bounds of every identifier ever issued
  std::vector<bool> retired;           // identifier has been cut by a split
  std::uint32_t splits = 0;

  explicit IdBook(std::uint32_t n)
      : k(n), is_start(n, false), id(n, 0), issued{{0, n}}, retired{false} {
    is_start[0] = true;
  }
  [[nodiscard]] bool would_split(std::uint32_t pos) const {
    return pos > 0 && pos < k && !is_start[pos];
  }
  void split(std::uint32_t pos) {
    std::uint32_t lo = pos;
    while (!is_start[lo]) --lo;
    std::uint32_t hi = pos + 1;
    while (hi < k && !is_start[hi]) ++hi;
    retired[id[pos]] = true;
    is_start[pos] = true;
    ++splits;
    const IntervalId left = 2 * splits - 1;
    const IntervalId right = 2 * splits;
    for (std::uint32_t i = lo; i < pos; ++i) id[i] = left;
    for (std::uint32_t i = pos; i < hi; ++i) id[i] = right;
    issued.push_back({lo, pos});
    issued.push_back({pos, hi});
    retired.push_back(false);
    retired.push_back(false);
  }
};

template <class S>
void require_ids(const S& s, const IdBook& book) {
  for (std::uint32_t i = 0; i < book.k; ++i) {
    INFO("record " << i);
    const IntervalId got = s.interval_of(i);
    REQUIRE(got == book.id[i]);
    REQUIRE_FALSE(book.retired[got]);  // a retired identifier is never returned again
    REQUIRE(s.bounds(got) == book.issued[got]);
  }
}

struct Handle {
  std::uint32_t sample;
  IntervalId id;
};

template <class S>
void run_contract_sequence(std::uint32_t k, std::uint64_t seed) {
  using Key = typename S::key_type;
  std::mt19937_64 rng(seed);
  std::vector<Key> keys(k);
  for (auto& x : keys) {
    x = (rng() % 5 == 0) ? kInf<Key> : static_cast<Key>(rng() % 8);
  }
  S s{std::span<const Key>(keys)};
  IdBook book(k);
  std::vector<Handle> handles;
  require_ids(s, book);
  const std::uint32_t ops = 5 * k + 30;
  for (std::uint32_t step = 0; step < ops; ++step) {
    INFO("k=" << k << " seed=" << seed << " step=" << step);
    const auto kind = rng() % 10;
    if (kind < 3) {
      const auto pos = static_cast<std::uint32_t>(rng() % (k + 1));
      const bool expect = book.would_split(pos);
      const IntervalId before = pos < k ? s.interval_of(pos) : 0;
      REQUIRE(s.split(pos) == expect);
      if (expect) {
        REQUIRE(before == book.id[pos]);
        book.split(pos);
        REQUIRE(book.retired[before]);
      }
    } else if (kind < 6) {
      const auto i = static_cast<std::uint32_t>(rng() % k);
      Key next = static_cast<Key>(rng() % 8);
      if (keys[i] < next) next = keys[i];  // a decrease lowers a key (l.579)
      keys[i] = next;
      s.decrease(i, next);
    } else if (kind < 8) {
      // Evaluate-style handle: remember a sample record and the identifier of its interval.
      const auto i = static_cast<std::uint32_t>(rng() % k);
      handles.push_back({i, s.interval_of(i)});
    } else {
      // The minimum lies in the queried interval and is its least (key, index) record.
      const auto i = static_cast<std::uint32_t>(rng() % k);
      const auto got = s.query(i);
      REQUIRE(s.interval_of(got.index) == s.interval_of(i));
      const IntervalBounds b = book.issued[book.id[i]];
      std::uint32_t best = b.lo;
      for (std::uint32_t j = b.lo + 1; j < b.hi; ++j) {
        if (keys[j] < keys[best]) best = j;
      }
      REQUIRE(got.index == best);
      REQUIRE(got.key == keys[best]);
    }
    require_ids(s, book);
    // A handle is current iff its interval has not been split since it was taken.
    for (const auto& hd : handles) {
      INFO("handle sample " << hd.sample << " id " << hd.id);
      REQUIRE((s.interval_of(hd.sample) == hd.id) == !book.retired[hd.id]);
    }
  }
  for (std::uint32_t pos = 1; pos < k; ++pos) {
    REQUIRE(s.split(pos) == book.would_split(pos));
    if (book.would_split(pos)) book.split(pos);
  }
  require_ids(s, book);
  // Identifiers are never reused: 2s + 1 identifiers were issued, at most 2k - 1 in total.
  REQUIRE(book.splits == k - 1);
  REQUIRE(book.issued.size() == 2 * std::size_t{k} - 1);
  for (const auto& hd : handles) {
    REQUIRE((s.interval_of(hd.sample) == hd.id) == !book.retired[hd.id]);
  }
}

}  // namespace

#define SFM_CONTRACT_IMPLEMENTATIONS                                                             \
  sssp::paper::SplitFindminRef<std::uint64_t>, sssp::paper::DyadicSplitFindmin<std::uint64_t>,   \
      sssp::paper::OneLevelBlockedSplitFindmin<std::uint64_t>,                                   \
      sssp::paper::SplitFindmin<std::uint64_t>, sssp::paper::SplitFindminRef<double>,            \
      sssp::paper::DyadicSplitFindmin<double>, sssp::paper::OneLevelBlockedSplitFindmin<double>, \
      sssp::paper::SplitFindmin<double>

TEMPLATE_TEST_CASE("Q-14: a split retires the old interval id and issues two fresh ids",
                   "[split_findmin][q14]", SFM_CONTRACT_IMPLEMENTATIONS) {
  TestType s(std::size_t{10});
  REQUIRE(s.interval_of(3) == 0);
  const IntervalId old = s.interval_of(3);
  REQUIRE(s.split(4));
  // Neither side keeps the old identifier.
  for (std::uint32_t i = 0; i < 10; ++i) REQUIRE(s.interval_of(i) != old);
  REQUIRE(s.interval_of(0) == 1);
  REQUIRE(s.interval_of(9) == 2);
  // Splitting the right part leaves the left part's identifier valid.
  REQUIRE(s.split(7));
  for (std::uint32_t i = 0; i < 4; ++i) REQUIRE(s.interval_of(i) == 1);
  for (std::uint32_t i = 4; i < 10; ++i) REQUIRE(s.interval_of(i) != 2);
  REQUIRE(s.interval_of(4) == 3);
  REQUIRE(s.interval_of(9) == 4);
  // Decreases, queries and a rejected split change no identifier.
  s.decrease(5, typename TestType::key_type{1});
  (void)s.query(5);
  REQUIRE_FALSE(s.split(7));
  REQUIRE(s.interval_of(0) == 1);
  REQUIRE(s.interval_of(5) == 3);
  REQUIRE(s.interval_of(8) == 4);
}

TEMPLATE_TEST_CASE("Q-14: stored interval ids are current exactly until their interval is split",
                   "[split_findmin][q14][random]", SFM_CONTRACT_IMPLEMENTATIONS) {
  for (std::uint32_t k = 1; k <= 40; ++k) {
    for (std::uint64_t seed = 1; seed <= 3; ++seed) {
      run_contract_sequence<TestType>(k, 7919 * k + seed);
    }
  }
  for (const std::uint32_t k : {63U, 64U, 65U, 200U, 513U}) {
    run_contract_sequence<TestType>(k, 104729 + k);
  }
}

TEMPLATE_TEST_CASE("Q-16: a decrease to the current key is accepted and changes nothing",
                   "[split_findmin][q16]", SFM_CONTRACT_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  const std::vector<Key> keys = {Key{3}, Key{1}, kInf<Key>, Key{1}, Key{0}, Key{7}, Key{2}};
  TestType s{std::span<const Key>(keys)};
  REQUIRE(s.split(3));
  for (std::uint32_t i = 0; i < keys.size(); ++i) {
    const auto before = s.query(i);
    const IntervalId id = s.interval_of(i);
    s.decrease(i, keys[i]);
    REQUIRE(s.key(i) == keys[i]);
    REQUIRE(s.query(i) == before);
    REQUIRE(s.interval_of(i) == id);
  }
  // Ties stay resolved by original index after equal-key decreases.
  REQUIRE(s.query(0).index == 1);
  REQUIRE(s.query(6).index == 4);
}

namespace {

// Everything a caller can observe: keys, minima, identifiers and their bounds.
template <class S>
struct Observed {
  std::vector<typename S::key_type> keys;
  std::vector<sssp::paper::SfmMin<typename S::key_type>> mins;
  std::vector<IntervalId> ids;
  std::vector<IntervalBounds> bounds;
  friend bool operator==(const Observed&, const Observed&) = default;
};

template <class S>
Observed<S> observe(const S& s) {
  Observed<S> o;
  for (std::uint32_t j = 0; j < s.size(); ++j) {
    o.keys.push_back(s.key(j));
    o.mins.push_back(s.query(j));
    o.ids.push_back(s.interval_of(j));
    o.bounds.push_back(s.bounds(s.interval_of(j)));
  }
  return o;
}

}  // namespace

TEMPLATE_TEST_CASE("Q-16: raising a key throws in every build and changes nothing",
                   "[split_findmin][q16]", SFM_CONTRACT_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  const std::vector<Key> keys = {Key{4}, Key{2}, Key{9}, Key{2}, Key{5}, Key{1}, Key{6}, Key{3}};
  TestType s{std::span<const Key>(keys)};
  REQUIRE(s.split(5));
  s.decrease(2, Key{0});  // current minimum of [0, 5)
  struct Raise {
    std::uint32_t i;
    Key to;
  };
  // Raise an interval minimum (would change query), a non-minimum, a minimum to +inf, and the
  // last record of the array.
  const std::vector<Raise> raises = {
      {2, Key{3}}, {0, Key{8}}, {5, kInf<Key>}, {7, Key{4}}, {1, kInf<Key>}};
  const auto before = observe(s);
  for (const auto& r : raises) {
    INFO("record " << r.i << " raised from " << s.key(r.i) << " to " << r.to);
    REQUIRE(s.key(r.i) < r.to);
    REQUIRE_THROWS_AS(s.decrease(r.i, r.to), std::invalid_argument);
    REQUIRE(observe(s) == before);
  }
  // The structure stays usable: a genuine decrease and a split still work afterwards.
  s.decrease(7, Key{0});
  REQUIRE(s.query(5) == sssp::paper::SfmMin<Key>{Key{0}, 7});
  REQUIRE(s.query(0) == sssp::paper::SfmMin<Key>{Key{0}, 2});
  REQUIRE(s.split(2));
  REQUIRE(s.query(0) == sssp::paper::SfmMin<Key>{Key{2}, 1});
  REQUIRE(s.query(2) == sssp::paper::SfmMin<Key>{Key{0}, 2});
  if constexpr (std::floating_point<Key>) {
    // NaN keys are rejected the same way (std::invalid_argument, nothing changes).
    const auto mid = observe(s);
    REQUIRE_THROWS_AS(s.decrease(3, std::numeric_limits<Key>::quiet_NaN()), std::invalid_argument);
    REQUIRE(observe(s) == mid);
  }
}

TEMPLATE_TEST_CASE("Q-16: random raises after random splits throw and change nothing",
                   "[split_findmin][q16][random]", SFM_CONTRACT_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  for (std::uint32_t k = 1; k <= 40; ++k) {
    std::mt19937_64 rng(31 * k + 5);
    std::vector<Key> keys(k);
    for (auto& x : keys) x = (rng() % 4 == 0) ? kInf<Key> : static_cast<Key>(rng() % 6);
    TestType s{std::span<const Key>(keys)};
    for (std::uint32_t step = 0; step < k / 2; ++step) {
      (void)s.split(static_cast<std::uint32_t>(rng() % (k + 1)));
    }
    const auto before = observe(s);
    for (std::uint32_t i = 0; i < k; ++i) {
      if (!(keys[i] < kInf<Key>)) continue;
      INFO("k=" << k << " record " << i << " key " << keys[i]);
      REQUIRE_THROWS_AS(s.decrease(i, static_cast<Key>(keys[i] + Key{1})), std::invalid_argument);
      REQUIRE_THROWS_AS(s.decrease(i, kInf<Key>), std::invalid_argument);
      REQUIRE(observe(s) == before);
    }
  }
}

TEMPLATE_TEST_CASE("record indices beyond the array throw in every build and change nothing",
                   "[split_findmin][range]", SFM_CONTRACT_IMPLEMENTATIONS) {
  using Key = typename TestType::key_type;
  for (const std::uint32_t k : {1U, 2U, 7U, 10U, 13U, 64U, 65U}) {
    INFO("k=" << k);
    std::vector<Key> keys(k);
    for (std::uint32_t j = 0; j < k; ++j) keys[j] = static_cast<Key>((5 * j + 3) % 11);
    TestType s{std::span<const Key>(keys)};
    std::uint32_t splits = 0;
    for (const std::uint32_t pos : {k / 2, k / 3, k - 1}) {
      if (s.split(pos)) ++splits;
    }
    const auto before = observe(s);
    // k itself, positions inside a padded final block or dyadic padding, and far beyond.
    for (const std::uint32_t bad :
         {k, k + 1, k + 2, k + 7, 2 * k + 64, std::numeric_limits<std::uint32_t>::max()}) {
      INFO("index " << bad);
      REQUIRE_THROWS_AS(s.decrease(bad, Key{0}), std::out_of_range);
      REQUIRE_THROWS_AS(s.decrease(bad, kInf<Key>), std::out_of_range);
      REQUIRE_THROWS_AS((void)s.query(bad), std::out_of_range);
      REQUIRE_THROWS_AS((void)s.key(bad), std::out_of_range);
      REQUIRE_THROWS_AS((void)s.interval_of(bad), std::out_of_range);
      REQUIRE(observe(s) == before);
    }
    // Identifiers that were never issued and split positions beyond the array.
    REQUIRE_THROWS_AS((void)s.bounds(2 * splits + 1), std::out_of_range);
    REQUIRE_THROWS_AS(s.split(k + 1), std::out_of_range);
    REQUIRE(observe(s) == before);
    // The last valid index still works.
    s.decrease(k - 1, Key{0});
    REQUIRE(s.key(k - 1) == Key{0});
    REQUIRE(s.query(k - 1).key == Key{0});
  }
}
