#pragma once
// Interval minima under splits and decreases (split-findmin), paper Sec 3.5.
//
// A fixed array of k >= 1 records starts as one interval. split(pos) cuts the interval containing
// pos at the boundary pos into two nonempty intervals; decrease(i, key) lowers a key; query(i)
// returns the minimum record (key, original index) of the interval containing i.
//
// Shared interface decisions (docs/research notes, Sec 3.5 items A8-A10):
//  - Order: (infinite flag, value, original index), lexicographic. Because
//    WeightTraits<Key>::infinity() is the top element of Key, this equals the order on
//    (key, index); an interval whose keys are all infinite reports its smallest index (A9).
//  - Interval IDs: the initial interval has ID 0; the s-th successful split assigns ID 2s-1 to
//    the left part and 2s to the right part. Both sides get fresh IDs ("two new owner
//    descriptors", Sec 3.5 step 3); callers must not keep IDs across a split (A8). IDs are never
//    reused, so at most 2k-1 IDs exist; bounds() of a retired ID reports its last bounds.
//  - split(pos) returns false (and changes nothing) when pos is already a boundary, including
//    pos = 0 and pos = k; it throws std::out_of_range when pos > k.
//  - decrease(i, key) with key above the current key violates the precondition (Sec 3.5 only
//    lowers keys, A10): it asserts in debug builds and is ignored otherwise. NaN keys throw.
//
// Implementations:
//  - SplitFindminRef<Key>: boundary std::set + segment tree, O(log k) per operation (reference).
//  - DyadicSplitFindmin<Key>: Lemma 3.10, dyadic tree of maximal nodes.
//  - OneLevelBlockedSplitFindmin<Key>: Lemma 3.11, first blocking step (lower = dyadic).
//  - SplitFindmin<Key>: Lemma 3.11, two blocking steps (lower = first-step structure).

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper {

using IntervalId = std::uint32_t;

/// Half-open interval [lo, hi) of record positions.
struct IntervalBounds {
  std::uint32_t lo = 0;
  std::uint32_t hi = 0;
  friend bool operator==(const IntervalBounds&, const IntervalBounds&) = default;
};

/// Minimum record of an interval: its key and its original index.
template <class Key>
struct SfmMin {
  Key key;
  std::uint32_t index;
  friend bool operator==(const SfmMin&, const SfmMin&) = default;
};

/// Keys: unsigned integers or floating point, with WeightTraits<Key>::infinity() as +inf.
template <class Key>
concept SplitFindminKey = (std::unsigned_integral<Key> || std::floating_point<Key>) &&
                          requires { WeightTraits<Key>::infinity(); };

/// Strict order on records: (key, index) lexicographic, see the header comment.
template <class Key>
[[nodiscard]] constexpr bool sfm_less(const SfmMin<Key>& a, const SfmMin<Key>& b) noexcept {
  if (a.key < b.key) return true;
  if (b.key < a.key) return false;
  return a.index < b.index;
}

template <class Key>
[[nodiscard]] constexpr bool sfm_is_nan(Key key) noexcept {
  if constexpr (std::floating_point<Key>) {
    return key != key;
  } else {
    return false;
  }
}

/// ceil(log2(x)) for x >= 1, by doubling (integer word operations only, Sec 2).
[[nodiscard]] constexpr std::uint32_t sfm_ceil_log2(std::uint64_t x) noexcept {
  std::uint32_t c = 0;
  while (c < 64 && (std::uint64_t{1} << c) < x) ++c;
  return c;
}

/// Lemma 3.11 block size B = max(2, ceil(log2(k+1))).
[[nodiscard]] constexpr std::uint32_t sfm_block_size(std::uint64_t k) noexcept {
  const std::uint32_t c = sfm_ceil_log2(k + 1);
  return c < 2 ? 2 : c;
}

/// Diagnostic counters, compared with the bounds of Lemmas 3.10 and 3.11 by property tests.
/// Counted with SSSP_COUNT, so they stay zero when SSSP_DIAGNOSTICS=0.
struct SplitFindminCounters {
  std::uint64_t splits = 0;            ///< successful top-level splits
  std::uint64_t rejected_splits = 0;   ///< split calls at an existing boundary
  std::uint64_t init_scans = 0;        ///< records scanned at initialisation (dyadic level)
  std::uint64_t rescans = 0;           ///< record rescans for exposed dyadic nodes (Lemma 3.10)
  std::uint64_t exposed_nodes = 0;     ///< dyadic nodes exposed along cut paths
  std::uint64_t cut_path_steps = 0;    ///< dyadic cut-path descents
  std::uint64_t list_node_visits = 0;  ///< dyadic nodes relabelled in the two new lists
  std::uint64_t lower_splits = 0;      ///< blocked: direct split calls into lower structures
  std::uint64_t upper_splits = 0;      ///< blocked: direct split calls into the upper structure

  SplitFindminCounters& operator+=(const SplitFindminCounters& o) noexcept {
    splits += o.splits;
    rejected_splits += o.rejected_splits;
    init_scans += o.init_scans;
    rescans += o.rescans;
    exposed_nodes += o.exposed_nodes;
    cut_path_steps += o.cut_path_steps;
    list_node_visits += o.list_node_visits;
    lower_splits += o.lower_splits;
    upper_splits += o.upper_splits;
    return *this;
  }

  /// Adds the counters to a Diagnostics::counters-style map as "<prefix><name>".
  void add_to(std::map<std::string, std::uint64_t>& out, std::string_view prefix) const {
    const std::string p(prefix);
    out[p + "splits"] += splits;
    out[p + "rejected_splits"] += rejected_splits;
    out[p + "init_scans"] += init_scans;
    out[p + "rescans"] += rescans;
    out[p + "exposed_nodes"] += exposed_nodes;
    out[p + "cut_path_steps"] += cut_path_steps;
    out[p + "list_node_visits"] += list_node_visits;
    out[p + "lower_splits"] += lower_splits;
    out[p + "upper_splits"] += upper_splits;
  }
};

/// Common interface of every split-findmin implementation.
template <class S>
concept SplitFindminLike =
    SplitFindminKey<typename S::key_type> && std::constructible_from<S, std::size_t> &&
    std::constructible_from<S, std::span<const typename S::key_type>> &&
    requires(S s, const S cs, std::uint32_t i, IntervalId id, typename S::key_type key) {
      { cs.size() } -> std::convertible_to<std::size_t>;
      { s.split(i) } -> std::same_as<bool>;
      s.decrease(i, key);
      { cs.query(i) } -> std::same_as<SfmMin<typename S::key_type>>;
      { cs.interval_of(i) } -> std::same_as<IntervalId>;
      { cs.bounds(id) } -> std::same_as<IntervalBounds>;
      { cs.key(i) } -> std::same_as<typename S::key_type>;
      { cs.counters() } -> std::same_as<SplitFindminCounters>;
    };

namespace sfm_detail {

struct FromVector {};

[[nodiscard]] inline std::uint32_t checked_size(std::size_t k) {
  if (k == 0) throw std::invalid_argument("split-findmin needs k >= 1 records");
  if (k > (std::size_t{1} << 31)) throw std::length_error("split-findmin supports k <= 2^31");
  return static_cast<std::uint32_t>(k);
}

template <class Key>
void check_keys(std::span<const Key> keys) {
  for (const Key key : keys) {
    if (sfm_is_nan(key)) throw std::invalid_argument("split-findmin key is NaN");
  }
}

template <class Key>
void check_key(Key key) {
  if (sfm_is_nan(key)) throw std::invalid_argument("split-findmin key is NaN");
}

}  // namespace sfm_detail

// ------------------------------------------------------------------------------------------
// Reference implementation: interval starts in a std::set, range minima in a segment tree.
// Not from the paper; used as the differential oracle for the paper's structures.
// ------------------------------------------------------------------------------------------
template <SplitFindminKey Key>
class SplitFindminRef {
 public:
  using key_type = Key;

  explicit SplitFindminRef(std::size_t k)
      : SplitFindminRef(sfm_detail::FromVector{}, std::vector<Key>(sfm_detail::checked_size(k),
                                                                   WeightTraits<Key>::infinity())) {
  }

  explicit SplitFindminRef(std::span<const Key> keys)
      : SplitFindminRef(sfm_detail::FromVector{}, std::vector<Key>(keys.begin(), keys.end())) {}

  [[nodiscard]] std::size_t size() const noexcept { return k_; }

  bool split(std::uint32_t pos) {
    if (pos > k_) throw std::out_of_range("split position beyond the array");
    if (pos == 0 || pos == k_ || starts_.contains(pos)) {
      SSSP_COUNT(++counters_.rejected_splits);
      return false;
    }
    const auto it = std::prev(starts_.upper_bound(pos));
    const std::uint32_t lo = *it;
    const std::uint32_t hi = bounds_[id_of_start_[lo]].hi;
    starts_.insert(pos);
    id_of_start_[lo] = static_cast<IntervalId>(bounds_.size());
    bounds_.push_back({lo, pos});
    id_of_start_[pos] = static_cast<IntervalId>(bounds_.size());
    bounds_.push_back({pos, hi});
    SSSP_COUNT(++counters_.splits);
    return true;
  }

  void decrease(std::uint32_t i, Key key) {
    assert(i < k_);
    sfm_detail::check_key(key);
    assert(!(tree_[leaf_ + i].key < key) && "split-findmin decrease must not raise a key");
    if (tree_[leaf_ + i].key < key) return;
    std::size_t v = leaf_ + i;
    tree_[v].key = key;
    for (v >>= 1; v >= 1; v >>= 1) {
      const auto& a = tree_[2 * v];
      const auto& b = tree_[2 * v + 1];
      tree_[v] = sfm_less(b, a) ? b : a;
    }
  }

  [[nodiscard]] SfmMin<Key> query(std::uint32_t i) const {
    assert(i < k_);
    const IntervalBounds b = bounds_[interval_of(i)];
    // Iterative bottom-up range minimum over [b.lo, b.hi).
    SfmMin<Key> best{WeightTraits<Key>::infinity(), std::numeric_limits<std::uint32_t>::max()};
    for (std::size_t l = leaf_ + b.lo, r = leaf_ + b.hi; l < r; l >>= 1, r >>= 1) {
      if (l & 1U) {
        if (sfm_less(tree_[l], best)) best = tree_[l];
        ++l;
      }
      if (r & 1U) {
        --r;
        if (sfm_less(tree_[r], best)) best = tree_[r];
      }
    }
    return best;
  }

  [[nodiscard]] IntervalId interval_of(std::uint32_t i) const {
    assert(i < k_);
    return id_of_start_[*std::prev(starts_.upper_bound(i))];
  }

  [[nodiscard]] IntervalBounds bounds(IntervalId id) const { return bounds_.at(id); }

  [[nodiscard]] Key key(std::uint32_t i) const {
    assert(i < k_);
    return tree_[leaf_ + i].key;
  }

  [[nodiscard]] SplitFindminCounters counters() const noexcept { return counters_; }

 private:
  SplitFindminRef(sfm_detail::FromVector /*tag*/, std::vector<Key> keys)
      : k_(sfm_detail::checked_size(keys.size())), id_of_start_(k_, 0) {
    sfm_detail::check_keys<Key>(keys);
    leaf_ = 1;
    while (leaf_ < k_) leaf_ <<= 1;
    tree_.assign(2 * leaf_, SfmMin<Key>{WeightTraits<Key>::infinity(),
                                        std::numeric_limits<std::uint32_t>::max()});
    for (std::uint32_t i = 0; i < k_; ++i) tree_[leaf_ + i] = {keys[i], i};
    for (std::size_t v = leaf_ - 1; v >= 1; --v) {
      tree_[v] = sfm_less(tree_[2 * v + 1], tree_[2 * v]) ? tree_[2 * v + 1] : tree_[2 * v];
    }
    starts_.insert(0);
    bounds_.push_back({0, k_});
  }

  std::uint32_t k_;
  std::size_t leaf_ = 1;
  std::vector<SfmMin<Key>> tree_;  ///< 1-based heap layout; padding leaves never win
  std::set<std::uint32_t> starts_;
  std::vector<IntervalId> id_of_start_;
  std::vector<IntervalBounds> bounds_;
  SplitFindminCounters counters_;
};

extern template class SplitFindminRef<std::uint64_t>;
extern template class SplitFindminRef<double>;

// ------------------------------------------------------------------------------------------
// Sec 3.5, Lemma 3.10: dyadic representation.
//
// The array is padded to K = 2^D >= k positions and covered by a balanced binary tree of
// intervals (heap numbering: root 1, children 2v and 2v+1). Padding positions are ignored: a
// node's range is clipped to [lo, min(hi, k)) and nodes lying wholly in the padding are never
// used. Each live interval is represented by the ordered list of maximal tree nodes it
// contains; each listed node stores its minimum and its interval descriptor, each record
// points to its node, and the descriptor stores the minimum over its nodes.
// Initialisation and all splits cost O(k(1 + log k)); decrease, query and interval_of are O(1).
// ------------------------------------------------------------------------------------------
template <SplitFindminKey Key>
class DyadicSplitFindmin {
 public:
  using key_type = Key;

  explicit DyadicSplitFindmin(std::size_t k)
      : DyadicSplitFindmin(
            sfm_detail::FromVector{},
            std::vector<Key>(sfm_detail::checked_size(k), WeightTraits<Key>::infinity())) {}

  explicit DyadicSplitFindmin(std::span<const Key> keys)
      : DyadicSplitFindmin(sfm_detail::FromVector{}, std::vector<Key>(keys.begin(), keys.end())) {}

  [[nodiscard]] std::size_t size() const noexcept { return k_; }

  /// Sec 3.5 "At a split": at most one listed node straddles the cut. If one does, follow the
  /// cut path below it, expose the maximal nodes left and right of the cut and rescan their
  /// records; splice them into the list in place of the straddling node; cut the list and
  /// relabel both halves with fresh descriptors, recomputing their minima.
  bool split(std::uint32_t pos) {
    if (pos > k_) throw std::out_of_range("split position beyond the array");
    if (pos == 0 || pos == k_ || node_desc_[node_of_[pos]] != node_desc_[node_of_[pos - 1]]) {
      SSSP_COUNT(++counters_.rejected_splits);
      return false;
    }
    const std::uint32_t n = node_of_[pos];
    const Descriptor old = descs_[node_desc_[n]];
    std::uint32_t left_head = old.head;
    std::uint32_t left_tail = 0;
    std::uint32_t right_head = 0;
    if (node_lo_[n] == pos) {
      // The cut already lies between two listed nodes.
      left_tail = prev_[n];
      right_head = n;
      assert(left_tail != 0);
    } else {
      // n straddles the cut: walk the cut path (Lemma 3.10 proof: at most D steps).
      std::array<std::uint32_t, 32> left{};
      std::array<std::uint32_t, 32> right{};
      std::size_t nl = 0;
      std::size_t nr = 0;
      for (std::uint32_t x = n;;) {
        SSSP_COUNT(++counters_.cut_path_steps);
        const std::uint32_t a = 2 * x;
        const std::uint32_t b = 2 * x + 1;
        const std::uint32_t mid = node_lo_[b];  // == k_ when b lies wholly in the padding
        if (pos < mid) {
          if (mid < k_) right[nr++] = b;  // b lies wholly right of the cut
          x = a;
        } else {
          left[nl++] = a;  // a lies wholly left of the cut
          if (pos == mid) {
            right[nr++] = b;
            break;
          }
          x = b;
        }
      }
      // right[] was collected from the top of the path down, i.e. from right to left.
      std::reverse(right.begin(), right.begin() + static_cast<std::ptrdiff_t>(nr));
      for (std::size_t j = 0; j < nl; ++j) expose(left[j]);
      for (std::size_t j = 0; j < nr; ++j) expose(right[j]);
      // Splice left ++ right into the list in place of n.
      std::uint32_t last = prev_[n];
      const std::uint32_t after = next_[n];
      const auto link = [&](std::uint32_t v) {
        prev_[v] = last;
        if (last != 0) next_[last] = v;
        last = v;
      };
      for (std::size_t j = 0; j < nl; ++j) link(left[j]);
      for (std::size_t j = 0; j < nr; ++j) link(right[j]);
      next_[last] = after;
      if (after != 0) prev_[after] = last;
      if (old.head == n) left_head = left[0];
      left_tail = left[nl - 1];
      right_head = right[0];
    }
    // Split the ordered list at the cut and scan both halves.
    next_[left_tail] = 0;
    prev_[right_head] = 0;
    new_descriptor(old.lo, pos, left_head);
    new_descriptor(pos, old.hi, right_head);
    SSSP_COUNT(++counters_.splits);
    return true;
  }

  /// Sec 3.5: "On a decrease, compare with both stored minima."
  void decrease(std::uint32_t i, Key key) {
    assert(i < k_);
    sfm_detail::check_key(key);
    assert(!(keys_[i] < key) && "split-findmin decrease must not raise a key");
    if (keys_[i] < key) return;
    keys_[i] = key;
    const SfmMin<Key> rec{key, i};
    const std::uint32_t n = node_of_[i];
    if (sfm_less(rec, node_min_[n])) node_min_[n] = rec;
    Descriptor& d = descs_[node_desc_[n]];
    if (sfm_less(rec, d.min)) d.min = rec;
  }

  /// Sec 3.5: "A lookup from a record follows its node's descriptor pointer."
  [[nodiscard]] SfmMin<Key> query(std::uint32_t i) const {
    assert(i < k_);
    return descs_[node_desc_[node_of_[i]]].min;
  }

  [[nodiscard]] IntervalId interval_of(std::uint32_t i) const {
    assert(i < k_);
    return node_desc_[node_of_[i]];
  }

  [[nodiscard]] IntervalBounds bounds(IntervalId id) const {
    const Descriptor& d = descs_.at(id);
    return {d.lo, d.hi};
  }

  [[nodiscard]] Key key(std::uint32_t i) const {
    assert(i < k_);
    return keys_[i];
  }

  [[nodiscard]] SplitFindminCounters counters() const noexcept { return counters_; }

 private:
  struct Descriptor {
    std::uint32_t lo;
    std::uint32_t hi;
    std::uint32_t head;  ///< first listed node (the list is terminated by node 0)
    SfmMin<Key> min;
  };

  DyadicSplitFindmin(sfm_detail::FromVector /*tag*/, std::vector<Key> keys)
      : k_(sfm_detail::checked_size(keys.size())), keys_(std::move(keys)), node_of_(k_, 1) {
    sfm_detail::check_keys<Key>(keys_);
    const std::uint32_t depth = sfm_ceil_log2(k_);
    const std::size_t padded = std::size_t{1} << depth;
    node_lo_.assign(2 * padded, k_);
    node_hi_.assign(2 * padded, k_);
    for (std::uint32_t d = 0; d <= depth; ++d) {
      const std::size_t len = padded >> d;
      for (std::size_t j = 0; j < (std::size_t{1} << d); ++j) {
        const std::size_t v = (std::size_t{1} << d) + j;
        const std::size_t lo = j * len;
        if (lo >= k_) continue;  // wholly in the padding: ignored
        node_lo_[v] = static_cast<std::uint32_t>(lo);
        node_hi_[v] = static_cast<std::uint32_t>(std::min<std::size_t>(lo + len, k_));
      }
    }
    node_min_.assign(2 * padded, SfmMin<Key>{WeightTraits<Key>::infinity(), 0});
    node_desc_.assign(2 * padded, 0);
    prev_.assign(2 * padded, 0);
    next_.assign(2 * padded, 0);
    // Initially the single interval is represented by the root alone.
    SfmMin<Key> m{keys_[0], 0};
    for (std::uint32_t i = 1; i < k_; ++i) {
      if (sfm_less(SfmMin<Key>{keys_[i], i}, m)) m = {keys_[i], i};
    }
    SSSP_COUNT(counters_.init_scans += k_);
    node_min_[1] = m;
    descs_.reserve(2 * static_cast<std::size_t>(k_) - 1);
    descs_.push_back({0, k_, 1, m});
  }

  /// Point the records of a newly exposed node at it and compute its minimum.
  void expose(std::uint32_t v) {
    SSSP_COUNT(++counters_.exposed_nodes);
    const std::uint32_t lo = node_lo_[v];
    const std::uint32_t hi = node_hi_[v];
    SfmMin<Key> m{keys_[lo], lo};
    for (std::uint32_t i = lo; i < hi; ++i) {
      node_of_[i] = v;
      if (sfm_less(SfmMin<Key>{keys_[i], i}, m)) m = {keys_[i], i};
    }
    SSSP_COUNT(counters_.rescans += hi - lo);
    node_min_[v] = m;
  }

  /// Assign a fresh descriptor to the list starting at head and recompute its minimum.
  void new_descriptor(std::uint32_t lo, std::uint32_t hi, std::uint32_t head) {
    const auto id = static_cast<IntervalId>(descs_.size());
    SfmMin<Key> m = node_min_[head];
    for (std::uint32_t v = head; v != 0; v = next_[v]) {
      SSSP_COUNT(++counters_.list_node_visits);
      node_desc_[v] = id;
      if (sfm_less(node_min_[v], m)) m = node_min_[v];
    }
    descs_.push_back({lo, hi, head, m});
  }

  std::uint32_t k_;
  std::vector<Key> keys_;
  std::vector<std::uint32_t> node_of_;  ///< record -> its maximal listed node
  // Per tree node (heap numbering; index 0 is unused and terminates lists):
  std::vector<std::uint32_t> node_lo_;
  std::vector<std::uint32_t> node_hi_;  ///< clipped to k_
  std::vector<SfmMin<Key>> node_min_;
  std::vector<IntervalId> node_desc_;
  std::vector<std::uint32_t> prev_;
  std::vector<std::uint32_t> next_;
  std::vector<Descriptor> descs_;
  SplitFindminCounters counters_;
};

extern template class DyadicSplitFindmin<std::uint64_t>;
extern template class DyadicSplitFindmin<double>;

// ------------------------------------------------------------------------------------------
// Sec 3.5, Lemma 3.11: one blocking step.
//
// The array is divided into consecutive blocks of size B (default B = max(2, ceil(log2(k+1))),
// final block possibly short). A Lower structure per block holds the intersections of live
// intervals with the block; an Upper structure holds one proxy per block. While a block is
// whole (never cut internally) its proxy key is the block's least key; once cut internally
// the proxy is isolated in a singleton upper interval and ignored. A live interval (owner)
// consists of at most a proper left end fragment, one run of whole blocks (one upper interval)
// and a proper right end fragment; the owner descriptor stores its endpoints, its pieces and the
// minimum of their minima. Fragment and run descriptors point back to their owner.
//
// Upper is DyadicSplitFindmin in both blocking steps of the paper; Lower is either
// DyadicSplitFindmin (first step) or the first-step structure itself (second step).
// ------------------------------------------------------------------------------------------
template <SplitFindminLike Upper, SplitFindminLike Lower>
  requires std::same_as<typename Upper::key_type, typename Lower::key_type>
class BlockedSplitFindmin {
 public:
  using key_type = typename Upper::key_type;
  using upper_type = Upper;
  using lower_type = Lower;

  /// k records with key +inf, default block size B = max(2, ceil(log2(k+1))).
  explicit BlockedSplitFindmin(std::size_t k) : BlockedSplitFindmin(k, std::uint32_t{0}) {}
  /// As above with an explicit block size (0 selects the default).
  BlockedSplitFindmin(std::size_t k, std::uint32_t block_size)
      : BlockedSplitFindmin(
            sfm_detail::FromVector{},
            std::vector<key_type>(sfm_detail::checked_size(k), WeightTraits<key_type>::infinity()),
            block_size) {}

  explicit BlockedSplitFindmin(std::span<const key_type> keys)
      : BlockedSplitFindmin(keys, std::uint32_t{0}) {}
  BlockedSplitFindmin(std::span<const key_type> keys, std::uint32_t block_size)
      : BlockedSplitFindmin(sfm_detail::FromVector{},
                            std::vector<key_type>(keys.begin(), keys.end()), block_size) {}

  [[nodiscard]] std::size_t size() const noexcept { return k_; }
  [[nodiscard]] std::uint32_t block_size() const noexcept { return block_; }
  [[nodiscard]] std::uint32_t num_blocks() const noexcept { return nblocks_; }
  /// Lower structure of block b (for inspection and tests).
  [[nodiscard]] const Lower& lower(std::uint32_t b) const { return lower_.at(b); }

  /// Lemma 3.11 step 3.
  bool split(std::uint32_t pos) {
    if (pos > k_) throw std::out_of_range("split position beyond the array");
    if (pos == 0 || pos == k_) {
      SSSP_COUNT(++counters_.rejected_splits);
      return false;
    }
    const IntervalId o = owner_of(pos);
    if (owner_of(pos - 1) != o) {
      SSSP_COUNT(++counters_.rejected_splits);
      return false;
    }
    const std::uint32_t b = pos / block_;
    const std::uint32_t off = pos - b * block_;
    if (off == 0) {
      // "At a cut between blocks, split the upper whole run if necessary": only when block b
      // and block b-1 lie in the same upper run (otherwise b already starts an upper interval).
      if (whole_[b] != 0 && upper_.bounds(upper_.interval_of(b)).lo < b) upper_split(b);
    } else {
      // "At an internal cut, split the relevant lower interval."
      [[maybe_unused]] const bool cut = lower_[b].split(off);
      assert(cut);
      SSSP_COUNT(++counters_.lower_splits);
      if (whole_[b] != 0) {
        // "If the block was whole, isolate its proxy by splitting off its nonempty left and
        // right neighbors in the upper run."
        whole_[b] = 0;
        if (upper_.bounds(upper_.interval_of(b)).lo < b) upper_split(b);
        if (upper_.bounds(upper_.interval_of(b)).hi > b + 1) upper_split(b + 1);
      }
    }
    // "Form the two new owner descriptors from their endpoint indices."
    const Owner old = owners_[o];
    make_owner(old.lo, pos);
    make_owner(pos, old.hi);
    SSSP_COUNT(++counters_.splits);
    return true;
  }

  /// Lemma 3.11 step 2: update the lower structure; in a whole block also decrease the proxy if
  /// the block minimum decreased; refresh the owner's minimum from its pieces.
  void decrease(std::uint32_t i, key_type key) {
    assert(i < k_);
    sfm_detail::check_key(key);
    const std::uint32_t b = i / block_;
    const std::uint32_t off = i - b * block_;
    assert(!(lower_[b].key(off) < key) && "split-findmin decrease must not raise a key");
    if (lower_[b].key(off) < key) return;
    lower_[b].decrease(off, key);
    if (whole_[b] != 0 && key < upper_.key(b)) upper_.decrease(b, key);
    Owner& w = owners_[owner_of(i)];
    w.min = pieces_min(w);
  }

  /// Lemma 3.11 step 1: reach the owner through the lower descriptor (split block) or through
  /// the upper descriptor of the block's proxy (whole block); return the owner's minimum.
  [[nodiscard]] SfmMin<key_type> query(std::uint32_t i) const { return owners_[owner_of(i)].min; }

  [[nodiscard]] IntervalId interval_of(std::uint32_t i) const { return owner_of(i); }

  [[nodiscard]] IntervalBounds bounds(IntervalId id) const {
    const Owner& w = owners_.at(id);
    return {w.lo, w.hi};
  }

  [[nodiscard]] key_type key(std::uint32_t i) const {
    assert(i < k_);
    const std::uint32_t b = i / block_;
    return lower_[b].key(i - b * block_);
  }

  /// Own split counts plus the dyadic work (scans, rescans, exposures, list visits) of all
  /// nested structures. lower_splits / upper_splits count this level's direct calls only.
  [[nodiscard]] SplitFindminCounters counters() const {
    SplitFindminCounters nested = upper_.counters();
    for (const Lower& l : lower_) nested += l.counters();
    SplitFindminCounters c = counters_;
    c.init_scans += nested.init_scans;
    c.rescans += nested.rescans;
    c.exposed_nodes += nested.exposed_nodes;
    c.cut_path_steps += nested.cut_path_steps;
    c.list_node_visits += nested.list_node_visits;
    return c;
  }

 private:
  static constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

  /// Owner descriptor: endpoints, up to three pieces, minimum of their minima. A left fragment
  /// is represented by record lo, a right fragment by record hi-1, a run by its first proxy.
  struct Owner {
    std::uint32_t lo;
    std::uint32_t hi;
    std::uint32_t left_block;   ///< block of the proper left fragment, or kNone
    std::uint32_t run_first;    ///< first whole block of the run, or kNone
    std::uint32_t right_block;  ///< block of the proper right fragment, or kNone
    SfmMin<key_type> min;
  };

  BlockedSplitFindmin(sfm_detail::FromVector /*tag*/, const std::vector<key_type>& keys,
                      std::uint32_t block_size)
      : k_(sfm_detail::checked_size(keys.size())),
        block_(block_size != 0 ? block_size : sfm_block_size(k_)),
        nblocks_((k_ + block_ - 1) / block_),
        lower_(make_lowers(keys, block_)),
        upper_(std::span<const key_type>(proxy_keys(lower_))),
        whole_(nblocks_, 1),
        owner_up_(2 * static_cast<std::size_t>(nblocks_) - 1, 0),
        owner_low_(2 * static_cast<std::size_t>(block_) * nblocks_, 0) {
    owners_.reserve(2 * static_cast<std::size_t>(k_) - 1);
    make_owner(0, k_);
  }

  static std::vector<Lower> make_lowers(const std::vector<key_type>& keys, std::uint32_t block) {
    sfm_detail::check_keys<key_type>(keys);
    const auto k = static_cast<std::uint32_t>(keys.size());
    // Second blocking step (Lemma 3.11 proof): a blocked lower structure on a block of size B
    // uses block size B' = max(2, ceil(log2(B+1))). The paper does not say whether B' is
    // derived from B or from the block's own length; we derive it from B for every block,
    // including a short final block (research notes A7, an interpretation).
    const std::uint32_t inner = sfm_block_size(block);
    std::vector<Lower> lowers;
    lowers.reserve((k + block - 1) / block);
    for (std::uint32_t lo = 0; lo < k; lo += block) {
      const std::uint32_t len = std::min(block, k - lo);
      const std::span<const key_type> part(keys.data() + lo, len);
      if constexpr (std::constructible_from<Lower, std::span<const key_type>, std::uint32_t>) {
        lowers.emplace_back(part, inner);
      } else {
        lowers.emplace_back(part);
      }
    }
    return lowers;
  }

  /// Proxy keys: the least original key of each (initially whole) block.
  static std::vector<key_type> proxy_keys(const std::vector<Lower>& lowers) {
    std::vector<key_type> keys;
    keys.reserve(lowers.size());
    for (const Lower& l : lowers) keys.push_back(l.query(0).key);
    return keys;
  }

  void upper_split(std::uint32_t b) {
    [[maybe_unused]] const bool cut = upper_.split(b);
    assert(cut);
    SSSP_COUNT(++counters_.upper_splits);
  }

  [[nodiscard]] std::size_t low_slot(std::uint32_t b, IntervalId lower_id) const noexcept {
    // A lower structure on at most B records has at most 2B-1 interval IDs.
    return 2 * static_cast<std::size_t>(block_) * b + lower_id;
  }

  [[nodiscard]] IntervalId owner_of(std::uint32_t i) const {
    assert(i < k_);
    const std::uint32_t b = i / block_;
    if (whole_[b] != 0) return owner_up_[upper_.interval_of(b)];
    return owner_low_[low_slot(b, lower_[b].interval_of(i - b * block_))];
  }

  [[nodiscard]] SfmMin<key_type> fragment_min(std::uint32_t b, std::uint32_t i) const {
    const SfmMin<key_type> m = lower_[b].query(i - b * block_);
    return {m.key, b * block_ + m.index};
  }

  [[nodiscard]] SfmMin<key_type> run_min(std::uint32_t first) const {
    const SfmMin<key_type> u = upper_.query(first);  // least proxy (ties: smaller block)
    const SfmMin<key_type> m = lower_[u.index].query(0);
    assert(!(m.key < u.key) && !(u.key < m.key));
    return {u.key, u.index * block_ + m.index};
  }

  [[nodiscard]] SfmMin<key_type> pieces_min(const Owner& w) const {
    SfmMin<key_type> m{WeightTraits<key_type>::infinity(), kNone};
    const auto take = [&m](const SfmMin<key_type>& c) {
      if (sfm_less(c, m)) m = c;
    };
    if (w.left_block != kNone) take(fragment_min(w.left_block, w.lo));
    if (w.run_first != kNone) take(run_min(w.run_first));
    if (w.right_block != kNone) take(fragment_min(w.right_block, w.hi - 1));
    return m;
  }

  /// Create the owner of [lo, hi): identify its pieces from the endpoints, point their
  /// descriptors at the new owner and compute its minimum.
  void make_owner(std::uint32_t lo, std::uint32_t hi) {
    const auto id = static_cast<IntervalId>(owners_.size());
    Owner w{lo, hi, kNone, kNone, kNone, {}};
    const std::uint32_t bl = lo / block_;
    const std::uint32_t br = (hi - 1) / block_;
    if (whole_[bl] == 0) {
      w.left_block = bl;
      owner_low_[low_slot(bl, lower_[bl].interval_of(lo - bl * block_))] = id;
    } else {
      assert(lo == bl * block_);
    }
    if (br != bl && whole_[br] == 0) {
      w.right_block = br;
      owner_low_[low_slot(br, lower_[br].interval_of(hi - 1 - br * block_))] = id;
    }
    assert(whole_[br] == 0 || hi == std::min(k_, (br + 1) * block_));
    const std::uint32_t first = whole_[bl] != 0 ? bl : bl + 1;
    const bool has_run = bl == br ? whole_[bl] != 0 : (first < br || whole_[br] != 0);
    if (has_run) {
      w.run_first = first;
      owner_up_[upper_.interval_of(first)] = id;
    }
    w.min = pieces_min(w);
    owners_.push_back(w);
  }

  std::uint32_t k_;
  std::uint32_t block_;
  std::uint32_t nblocks_;
  std::vector<Lower> lower_;
  Upper upper_;
  std::vector<std::uint8_t> whole_;    ///< block never cut internally
  std::vector<IntervalId> owner_up_;   ///< upper interval ID -> owner
  std::vector<IntervalId> owner_low_;  ///< (block, lower interval ID) -> owner
  std::vector<Owner> owners_;          ///< owner descriptors = public interval IDs
  SplitFindminCounters counters_;
};

/// First blocking step of Lemma 3.11: O(k(1 + log log(k+4))) for initialisation and splits.
template <SplitFindminKey Key>
using OneLevelBlockedSplitFindmin =
    BlockedSplitFindmin<DyadicSplitFindmin<Key>, DyadicSplitFindmin<Key>>;

/// Lemma 3.11: dyadic upper structure, first-step blocked lower structures;
/// O(q + q' + k tau(k)) in total, worst-case O(1) query, decrease and descriptor lookup.
template <SplitFindminKey Key>
using SplitFindmin = BlockedSplitFindmin<DyadicSplitFindmin<Key>, OneLevelBlockedSplitFindmin<Key>>;

extern template class BlockedSplitFindmin<DyadicSplitFindmin<std::uint64_t>,
                                          DyadicSplitFindmin<std::uint64_t>>;
extern template class BlockedSplitFindmin<DyadicSplitFindmin<double>, DyadicSplitFindmin<double>>;
extern template class BlockedSplitFindmin<DyadicSplitFindmin<std::uint64_t>,
                                          OneLevelBlockedSplitFindmin<std::uint64_t>>;
extern template class BlockedSplitFindmin<DyadicSplitFindmin<double>,
                                          OneLevelBlockedSplitFindmin<double>>;

}  // namespace sssp::paper
