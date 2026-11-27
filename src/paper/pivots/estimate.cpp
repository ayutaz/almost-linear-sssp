// Sec 6.2, Algorithm 6 Estimate(S, Q), Lemmas 6.4-6.7 (arXiv 2611.34911v1, txt l.1455-1600).
// See estimate.hpp for the contract. Step comments cite the algorithm lines ("Alg 6 l.k") and the
// paragraphs "Bands and quotients" (l.1473-1488) and "Centers, partitions, and table indices"
// (l.1490-1508).
//
// Interpretations of points the paper leaves open (marked 【解釈】 at the code):
//   I-1 (Q-06, Q-43) Power-of-two multiples are exact: multiword integers for integer weights,
//       exact doublings for doubles. j in {0..H-1} is drawn as eta = log2 H random bits.
//   I-2 The paper scales every quotient edge weight by H; here the heap orders unscaled
//       distances (scaling by H is monotone and exact) and only the radius test forms H d by
//       eta doublings per settled vertex. Same comparisons, same O(log(t+2)) additions per
//       scaled number (Lemma 6.7(3)).
//   I-3 (Q-41) The part above the coarsest scale is the quotient component; Lemma 6.7(1) makes
//       the coarsest scale's refinement one part per component, which is checked.
//   I-4 Ties in the sample sort are broken by edge index (any order satisfies Sec 6.2).
//   I-5 Count arrays for grouping (parent part, center) are bounded by |V(Q_a)| of the band
//       instead of per component; the work per scale stays O(|V(Q_a)|).
//   I-6 The random draws of a band are made from the coarsest scale to the finest.
//   I-7 The binary merge forest (Kruskal with DisjointSets) is built here: the shared
//       simultaneous-merge builder (task P3-MERGE-EVENTS) produces non-binary nodes.
#include "paper/pivots/estimate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/model/ilog.hpp"
#include "paper/model/length.hpp"
#include "paper/primitives/dijkstra.hpp"
#include "paper/primitives/dsu.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "paper/primitives/search.hpp"

namespace sssp::paper::pivots {
namespace {

// ------------------------------------------------------------------ exact scaled lengths

/// Nonnegative integer of up to 512 bits. Only the operations of the comparison-addition model
/// are offered (copy, addition, doubling = self-addition, comparison); it holds power-of-two
/// multiples of integer lengths exactly (research note Q-06 / Q-43). 512 bits cover every value
/// formed for t < 2^32 (at most about 2^350); larger values throw std::overflow_error.
class WideUint {
 public:
  static constexpr std::size_t kWords = 8;

  WideUint() = default;
  explicit WideUint(std::uint64_t v) noexcept { w_[0] = v; }

  void add(const WideUint& o) {
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < kWords; ++i) {
      std::uint64_t s = 0;
      const bool c1 = __builtin_add_overflow(w_[i], o.w_[i], &s);
      const bool c2 = __builtin_add_overflow(s, carry, &s);
      w_[i] = s;
      carry = (c1 || c2) ? 1 : 0;
    }
    if (carry != 0) throw std::overflow_error("scaled length exceeds 512 bits");
  }
  void twice() { add(*this); }

  friend bool operator<(const WideUint& a, const WideUint& b) noexcept {
    for (std::size_t i = kWords; i-- > 0;) {
      if (a.w_[i] != b.w_[i]) return a.w_[i] < b.w_[i];
    }
    return false;
  }

 private:
  std::array<std::uint64_t, kWords> w_{};
};

template <class W>
struct Exact;

template <>
struct Exact<std::uint64_t> {
  using type = WideUint;
  static type of(std::uint64_t w) noexcept { return WideUint(w); }
  static void add(type& x, const type& y) { x.add(y); }
  static bool less(const type& a, const type& b) noexcept { return a < b; }
};

template <>
struct Exact<double> {
  using type = double;
  static type of(double w) noexcept { return w; }
  static void add(type& x, const type& y) {
    x = x + y;
    if (!std::isfinite(x)) throw std::overflow_error("scaled length overflows double");
  }
  static bool less(const type& a, const type& b) noexcept { return a < b; }
};

/// Counted arithmetic on scaled lengths: doublings and additions count as additions, tests as
/// comparisons (Diagnostics.ops).
template <class W>
class Arith {
 public:
  using X = typename Exact<W>::type;

  explicit Arith(OpCounters* ops = nullptr) noexcept : ops_(ops) {}

  [[nodiscard]] X of(W w) const noexcept { return Exact<W>::of(w); }
  /// x 2^k by k doublings (Sec 2: "multiplication of a real by 2^j takes j doublings").
  [[nodiscard]] X scaled(X x, unsigned k) const {
    for (unsigned i = 0; i < k; ++i) {
      const X y = x;
      Exact<W>::add(x, y);
    }
    count_add(k);
    return x;
  }
  [[nodiscard]] X scaled_len(W w, unsigned k) const { return scaled(of(w), k); }
  [[nodiscard]] X plus(X a, const X& b) const {
    Exact<W>::add(a, b);
    count_add(1);
    return a;
  }
  [[nodiscard]] bool less(const X& a, const X& b) const {
    count_cmp();
    return Exact<W>::less(a, b);
  }
  /// Z = 4 Delta + a = a 2^(i+2) + a (Alg 6 l.9).
  [[nodiscard]] X table_value(const TableBound<W>& z) const {
    return plus(scaled_len(z.anchor, z.scale + 2), of(z.anchor));
  }
  /// (H + j) Delta with Delta = 2^i a, H = 2^eta, by binary doubling and addition (Sec 6.2,
  /// l.1493-1494): sum Delta 2^b over the one bits b of H + j.
  [[nodiscard]] X radius_threshold(W a, unsigned i, unsigned eta,
                                   std::span<const std::uint64_t> j) const {
    X cur = scaled_len(a, i);
    X acc = of(W{0});
    for (unsigned b = 0; b < eta; ++b) {
      const std::size_t word = b >> 6U;
      if (word < j.size() && ((j[word] >> (b & 63U)) & 1U) != 0) acc = plus(acc, cur);
      cur = scaled(cur, 1);
    }
    return plus(acc, cur);  // the bit of H itself
  }

 private:
  void count_add([[maybe_unused]] std::uint64_t k) const {
    SSSP_COUNT(if (ops_ != nullptr) ops_->additions += k);
  }
  void count_cmp() const { SSSP_COUNT(if (ops_ != nullptr)++ ops_->comparisons); }

  OpCounters* ops_;
};

/// Number of 64-bit words of a j in {0..2^eta - 1}.
std::size_t j_words(unsigned eta) {
  return (static_cast<std::size_t>(eta) + 63U) >> 6U;
}

void check_j(unsigned eta, std::span<const std::uint64_t> j) {
  const std::size_t words = j_words(eta);
  for (std::size_t k = 0; k < j.size(); ++k) {
    std::uint64_t allowed = ~std::uint64_t{0};
    if (k >= words) {
      allowed = 0;
    } else if (k + 1 == words && (eta & 63U) != 0) {
      allowed = (std::uint64_t{1} << (eta & 63U)) - 1;
    }
    if ((j[k] & ~allowed) != 0) throw std::invalid_argument("within_radius: j must be < H");
  }
}

/// j uniform in {0..H-1}, H = 2^eta: eta unbiased bits (Sec 2: ceil(log2 H) bits; a power of two
/// never rejects).
std::vector<std::uint64_t> draw_j(unsigned eta, model::WordRng& rng) {
  std::vector<std::uint64_t> j(j_words(eta), 0);
  unsigned left = eta;
  for (auto& word : j) {
    const unsigned c = left < 64 ? left : 64;
    word = rng.next_bits(c);
    left -= c;
  }
  return j;
}

/// x * m for an integer x (WideUint) and a word m, by shift and add (integer size arithmetic).
WideUint mul_word(const WideUint& x, std::uint64_t m) {
  WideUint acc(0);
  WideUint cur = x;
  while (m != 0) {
    if ((m & 1U) != 0) acc.add(cur);
    m >>= 1U;
    if (m != 0) cur.twice();
  }
  return acc;
}

/// Sec 6.2 l.1491: least eta with 2^eta >= (16t)^5, in integer word arithmetic.
unsigned paper_h_log2(std::uint64_t t) {
  WideUint target(1);
  for (int k = 0; k < 5; ++k) target = mul_word(target, 16 * t);
  WideUint p(1);
  unsigned eta = 0;
  while (p < target) {
    p.twice();
    ++eta;
  }
  return eta;
}

constexpr unsigned kMaxHLog2 = 256;

// ------------------------------------------------------------------ input validation

template <class W>
void validate(std::size_t t, std::span<const Edge<W>> sample, std::span<const Demand> demands) {
  if (t >= (std::size_t{1} << 32)) throw std::length_error("estimate: t >= 2^32");
  if (sample.size() >= kDisconnected || demands.size() >= kDisconnected) {
    throw std::length_error("estimate: too many edges or demands");
  }
  for (std::size_t i = 0; i < sample.size(); ++i) {
    const auto& e = sample[i];
    if (e.u >= t || e.v >= t) {
      throw std::out_of_range("estimate: sample edge " + std::to_string(i) + " leaves 0..t-1");
    }
    if (e.u == e.v) throw std::invalid_argument("estimate: sample edge is a loop");
    if (!WeightTraits<W>::is_valid_weight(e.weight) || e.weight == WeightTraits<W>::infinity() ||
        !(W{0} < e.weight)) {
      throw std::invalid_argument("estimate: sample weights must be positive and finite");
    }
  }
  for (const auto& q : demands) {
    if (q.x >= t || q.y >= t) throw std::out_of_range("estimate: demand endpoint >= t");
    if (q.x == q.y) throw std::invalid_argument("estimate: demand endpoints must differ");
  }
}

// ------------------------------------------------------------------ sorting

/// Lemma 6.7(3): sort the sample by (weight, index) by repeatedly merging adjacent sorted runs
/// (bottom-up merge sort, O(t log(t+2)) comparisons). Returns edge indices in sorted order.
template <class W>
std::vector<std::uint32_t> sort_sample(std::span<const Edge<W>> s, const model::LenOps<W>& ops) {
  const std::size_t m = s.size();
  std::vector<std::uint32_t> a(m);
  std::iota(a.begin(), a.end(), 0U);
  std::vector<std::uint32_t> b(m);
  const auto less = [&](std::uint32_t x, std::uint32_t y) {
    return ops.less(model::EdgeKey<W>{s[x].weight, x}, model::EdgeKey<W>{s[y].weight, y});
  };
  for (std::size_t width = 1; width < m; width *= 2) {
    for (std::size_t lo = 0; lo < m; lo += 2 * width) {
      const std::size_t mid = std::min(lo + width, m);
      const std::size_t hi = std::min(lo + 2 * width, m);
      std::size_t i = lo;
      std::size_t k = mid;
      std::size_t o = lo;
      while (i < mid && k < hi) b[o++] = less(a[k], a[i]) ? a[k++] : a[i++];
      while (i < mid) b[o++] = a[i++];
      while (k < hi) b[o++] = a[k++];
    }
    a.swap(b);
  }
  return a;
}

// ------------------------------------------------------------------ first covering centers

template <class W>
struct CoverState {
  std::vector<W> best;  ///< best distance offered by earlier centers within R
  std::vector<std::uint32_t> best_stamp;
  std::vector<W> tent;  ///< tentative distances of the current search (timestamped)
  std::vector<std::uint32_t> tent_stamp;
  std::vector<std::uint32_t> done_stamp;
  std::uint32_t scale_id = 0;
  std::uint32_t search_id = 0;
  IndexedBinaryHeap<W> heap;

  void prepare(std::size_t n) {
    if (best.size() < n) {
      best.resize(n);
      best_stamp.resize(n, 0);
      tent.resize(n);
      tent_stamp.resize(n, 0);
      done_stamp.resize(n, 0);
    }
  }
};

struct CoverCounts {
  std::uint64_t searches = 0;
  std::uint64_t settles = 0;
  std::uint64_t scans = 0;
};

/// Sec 6.2 (l.1495-1500), Alg 6 l.7: one radius-limited Dijkstra search per center in
/// permutation order. A settled vertex whose distance does not strictly improve the stored value
/// is discarded (its adjacency is not scanned); the first improvement sets the permanent center
/// label. Lemma 6.6: the result is each vertex's first covering center in the order.
template <class W>
std::vector<Vertex> first_covering_centers(const CsrGraph<W>& q,
                                           std::span<const std::uint32_t> order,
                                           const typename Arith<W>::X& threshold, unsigned eta,
                                           const Arith<W>& ar, OpCounters* ops, CoverState<W>& st,
                                           CoverCounts& cnt) {
  const std::size_t n = q.num_vertices();
  std::vector<Vertex> center(n, kNoVertex);
  st.prepare(n);
  const std::uint32_t scale = ++st.scale_id;
  auto& heap = st.heap;
  for (const std::uint32_t c : order) {
    const std::uint32_t sid = ++st.search_id;
    ++cnt.searches;
    heap.reset(n);
    heap.set_counters(ops);
    st.tent[c] = W{0};
    st.tent_stamp[c] = sid;
    heap.push(c, W{0});
    while (!heap.empty()) {
      const auto [d, z] = heap.pop_min();
      st.done_stamp[z] = sid;
      ++cnt.settles;
      // Radius-limited: stop once the smallest value exceeds R, i.e. H d > (H + j) Delta.
      // 【解釈】I-2: only this test is H-scaled; the heap orders unscaled distances.
      if (ar.less(threshold, ar.scaled_len(d, eta))) break;
      // Discard the visit unless its distance is strictly smaller than the stored value.
      SSSP_COUNT(if (ops != nullptr)++ ops->comparisons);
      if (st.best_stamp[z] == scale && !(d < st.best[z])) continue;
      st.best[z] = d;
      st.best_stamp[z] = scale;
      if (center[z] == kNoVertex) center[z] = c;  // first improvement: permanent label
      ++cnt.scans;
      SSSP_COUNT(if (ops != nullptr)++ ops->vertex_scans);
      for (const auto& arc : q.arcs(z)) {
        const Vertex y = arc.to;
        SSSP_COUNT(if (ops != nullptr)++ ops->edge_scans);
        if (st.done_stamp[y] == sid) continue;
        SSSP_COUNT(if (ops != nullptr)++ ops->additions);
        const W nd = WeightTraits<W>::checked_add(d, arc.weight);
        if (st.tent_stamp[y] != sid) {
          st.tent_stamp[y] = sid;
          st.tent[y] = nd;
          heap.push(y, nd);
        } else {
          SSSP_COUNT(if (ops != nullptr)++ ops->comparisons);
          if (nd < st.tent[y]) {
            st.tent[y] = nd;
            heap.decrease(y, nd);
          }
        }
      }
    }
  }
  heap.set_counters(nullptr);
  return center;
}

}  // namespace

// ------------------------------------------------------------------ constants

EstimateConstants estimate_constants(std::size_t t, const EstimateParams& params) {
  if (t == 0) throw std::invalid_argument("estimate_constants: t must be positive");
  if (t >= (std::size_t{1} << 32)) throw std::length_error("estimate_constants: t >= 2^32");
  const auto tt = static_cast<std::uint64_t>(t);
  EstimateConstants c;
  c.tau = model::ceil_log2(4 * tt);          // T = least power of two >= 4t (l.1473)
  c.ell = 2 * c.tau + model::ceil_log2(tt);  // least l with 2^l >= t T^2 (l.1490)
  if (params.h_log2.has_value()) {           // Mode B, deviation D-9
    if (*params.h_log2 > kMaxHLog2) throw std::invalid_argument("estimate: h_log2 > 256");
    c.eta = *params.h_log2;
  } else {
    c.eta = paper_h_log2(tt);  // H = least power of two >= (16t)^5 (l.1491)
  }
  return c;
}

// ------------------------------------------------------------------ exact tests

template <class W>
bool bound_at_least(const TableBound<W>& z, W d) {
  const Arith<W> ar;
  return !ar.less(ar.table_value(z), ar.of(d));
}

template <class W>
bool scaled_at_least(W lambda, unsigned k, const TableBound<W>& z) {
  const Arith<W> ar;
  return !ar.less(ar.scaled_len(lambda, k), ar.table_value(z));
}

template <class W>
bool pow2_multiple_at_least(W lambda, unsigned k, W d) {
  const Arith<W> ar;
  return !ar.less(ar.scaled_len(lambda, k), ar.of(d));
}

template <class W>
bool within_radius(W d, W anchor, unsigned scale, unsigned eta, std::span<const std::uint64_t> j) {
  check_j(eta, j);
  const Arith<W> ar;
  return !ar.less(ar.radius_threshold(anchor, scale, eta, j), ar.scaled_len(d, eta));
}

// ------------------------------------------------------------------ Algorithm 6

template <class W>
EstimateResult<W> estimate(std::size_t t, std::span<const Edge<W>> sample,
                           std::span<const Demand> demands, model::WordRng& rng,
                           const EstimateParams& params, EstimateTrace<W>* trace) {
  validate<W>(t, sample, demands);
  using X = typename Arith<W>::X;
  EstimateResult<W> res;
  const std::size_t m = sample.size();
  const std::size_t h = demands.size();
  res.index.assign(h, kDisconnected);
  auto& cnt = res.diag.counters;
  cnt["estimate.calls"] += 1;
  if (trace != nullptr) {
    *trace = EstimateTrace<W>{};
    trace->demand.assign(h, {});
  }
  if (t < 2) {  // no loopless edge and no demand with distinct endpoints exists
    cnt["estimate.disconnected_demands"] += h;
    return res;
  }
  const EstimateConstants k = estimate_constants(t, params);
  if (trace != nullptr) trace->constants = k;
  OpCounters* ops = &res.diag.ops;
  const Arith<W> ar(ops);
  const model::LenOps<W> lops(ops);

  // Alg 6 l.2, "Bands and quotients": sort the sample by weight (【解釈】I-4: ties by index).
  const std::vector<std::uint32_t> order = sort_sample<W>(sample, lops);
  const auto weight_at = [&](std::size_t pos) { return sample[order[pos]].weight; };
  // T^2 w for every sample weight, computed once ("shared doublings", l.1486).
  std::vector<X> t2w(m);
  for (std::size_t pos = 0; pos < m; ++pos) t2w[pos] = ar.scaled_len(weight_at(pos), 2 * k.tau);

  // The least unassigned weight a starts a band containing the remaining weights below T a.
  std::vector<std::uint32_t> band_of_pos(m);
  std::vector<W> anchors;
  X anchor_t{};
  for (std::size_t pos = 0; pos < m; ++pos) {
    const W w = weight_at(pos);
    if (anchors.empty() || !ar.less(ar.of(w), anchor_t)) {
      anchors.push_back(w);
      anchor_t = ar.scaled_len(w, k.tau);
    }
    band_of_pos[pos] = static_cast<std::uint32_t>(anchors.size() - 1);
  }

  // 【解釈】I-7: binary merge forest built here. Kruskal sweep: each successful join creates a
  // binary parent labelled with its weight and band; the selected forest edges are kept
  // (l.1476-1479).
  DisjointSets dsu(t);
  std::vector<Vertex> top(t);
  std::iota(top.begin(), top.end(), Vertex{0});
  std::vector<Vertex> fparent(2 * t - 1, kNoVertex);
  std::vector<W> node_weight(2 * t - 1, W{0});
  std::vector<std::uint32_t> node_band(2 * t - 1, kDisconnected);
  std::vector<std::uint32_t> forest_pos;
  auto next = static_cast<Vertex>(t);
  for (std::size_t pos = 0; pos < m; ++pos) {
    const auto& e = sample[order[pos]];
    const auto ru = dsu.find(e.u);
    const auto rv = dsu.find(e.v);
    if (ru == rv) continue;
    fparent[top[ru]] = next;
    fparent[top[rv]] = next;
    node_weight[next] = e.weight;
    node_band[next] = band_of_pos[pos];
    dsu.unite(ru, rv);
    top[dsu.find(ru)] = next;
    ++next;
    forest_pos.push_back(static_cast<std::uint32_t>(pos));
  }
  fparent.resize(next);
  const ForestQueries merge_forest(fparent);

  // Alg 6 l.2-3: flag disconnected demands; group the others by the band of their merge (the
  // LCA in the merge forest, Lemma 3.9).
  std::vector<std::uint32_t> connected;
  std::vector<std::uint32_t> dband(h, kDisconnected);
  for (std::size_t q = 0; q < h; ++q) {
    const auto l = merge_forest.lca(demands[q].x, demands[q].y);
    if (!l.has_value()) continue;
    dband[q] = node_band[*l];
    if (trace != nullptr) trace->demand[q].beta = node_weight[*l];
    connected.push_back(static_cast<std::uint32_t>(q));
  }
  cnt["estimate.disconnected_demands"] += h - connected.size();
  const std::vector<std::uint32_t> by_band =
      counting_sort_stable(connected, anchors.size(), [&](std::uint32_t q) { return dband[q]; });

  // Contraction state: a component label at every original vertex, relabelling the smaller list
  // at each union (l.1484-1485).
  std::vector<Vertex> label(t);
  std::iota(label.begin(), label.end(), Vertex{0});
  std::vector<std::vector<Vertex>> members(t);
  for (Vertex v = 0; v < t; ++v) members[v].push_back(v);
  const auto unite_labels = [&](Vertex u, Vertex v) {
    Vertex a = label[u];
    Vertex b = label[v];
    if (a == b) throw std::logic_error("estimate: selected forest edge inside a component");
    if (members[a].size() < members[b].size()) std::swap(a, b);
    for (const Vertex x : members[b]) {
      label[x] = a;
      members[a].push_back(x);
    }
    std::vector<Vertex>().swap(members[b]);
  };

  std::vector<std::uint32_t> loc_stamp(t, 0);
  std::vector<Vertex> loc_id(t, kNoVertex);
  std::vector<std::uint32_t> windows(m, 0);
  std::uint32_t band_stamp = 0;
  std::size_t f = 0;   // next selected forest edge to contract
  std::size_t lo = 0;  // window cursors in the sorted sample
  std::size_t hi = 0;
  CoverState<W> cover;
  CoverCounts cc;
  std::uint64_t quotient_vertices = 0;
  std::uint64_t quotient_edges = 0;
  std::uint64_t vertex_scales = 0;
  std::uint64_t bands_processed = 0;

  // Alg 6 l.4: each band anchor with a nonempty demand group, in increasing order.
  for (std::size_t g = 0; g < by_band.size();) {
    const std::uint32_t band = dband[by_band[g]];
    std::size_t g_end = g;
    while (g_end < by_band.size() && dband[by_band[g_end]] == band) ++g_end;
    const W a = anchors[band];
    const X ax = ar.of(a);
    const X a_t2 = ar.scaled_len(a, 2 * k.tau);  // T^2 a

    // Alg 6 l.5: contract every edge of weight strictly below a / T^2 (test T^2 w < a): insert a
    // selected forest edge as soon as its weight falls below the cutoff.
    while (f < forest_pos.size() && ar.less(t2w[forest_pos[f]], ax)) {
      const auto& e = sample[order[forest_pos[f]]];
      unite_labels(e.u, e.v);
      ++f;
    }
    // Window [a / T^2, T^2 a] between two cursors: T^2 w >= a at the lower end, w <= T^2 a at
    // the upper end.
    while (lo < m && ar.less(t2w[lo], ax)) ++lo;
    hi = std::max(hi, lo);
    while (hi < m && !ar.less(a_t2, ar.of(weight_at(hi)))) ++hi;
    // Retain window edges with distinct contracted endpoints; temporary local identifiers for
    // the components encountered on retained edges.
    ++band_stamp;
    std::size_t nq = 0;
    const auto local = [&](Vertex lab) {
      if (loc_stamp[lab] != band_stamp) {
        loc_stamp[lab] = band_stamp;
        loc_id[lab] = static_cast<Vertex>(nq++);
      }
      return loc_id[lab];
    };
    const auto image = [&](Vertex v) {
      const Vertex lab = label[v];
      return loc_stamp[lab] == band_stamp ? loc_id[lab] : kNoVertex;
    };
    std::vector<Edge<W>> qedges;
    std::vector<std::uint32_t> retained;
    for (std::size_t pos = lo; pos < hi; ++pos) {
      const auto& e = sample[order[pos]];
      const Vertex lu = label[e.u];
      const Vertex lv = label[e.v];
      if (lu == lv) continue;
      const Vertex qu = local(lu);
      const Vertex qv = local(lv);
      qedges.push_back({qu, qv, e.weight});
      retained.push_back(order[pos]);
      ++windows[order[pos]];
    }
    const CsrGraph<W> qg = CsrGraph<W>::from_edges(nq, qedges);
    const Components comp = connected_components(qg);
    ++bands_processed;
    quotient_vertices += nq;
    quotient_edges += qedges.size();
    vertex_scales += static_cast<std::uint64_t>(nq) * (k.ell + 1);

    // Map the demand endpoints by their current component labels. Lemma 6.4: the images are
    // distinct, present and connected.
    std::vector<Vertex> qx(g_end - g);
    std::vector<Vertex> qy(g_end - g);
    for (std::size_t r = g; r < g_end; ++r) {
      const Demand& d = demands[by_band[r]];
      qx[r - g] = image(d.x);
      qy[r - g] = image(d.y);
      if (qx[r - g] == kNoVertex || qy[r - g] == kNoVertex || qx[r - g] == qy[r - g] ||
          comp.label[qx[r - g]] != comp.label[qy[r - g]]) {
        throw std::logic_error("estimate: demand images violate Lemma 6.4");
      }
    }

    typename EstimateTrace<W>::BandTrace* bt = nullptr;
    if (trace != nullptr) {
      trace->bands.emplace_back();
      bt = &trace->bands.back();
      bt->anchor = a;
      bt->image.resize(t);
      for (Vertex v = 0; v < t; ++v) bt->image[v] = image(v);
      bt->quotient_edges = qedges;
      bt->retained = retained;
      bt->num_quotient_vertices = nq;
      bt->component.assign(comp.label.begin(), comp.label.end());
      const auto bi = static_cast<std::uint32_t>(trace->bands.size() - 1);
      for (std::size_t r = g; r < g_end; ++r) {
        auto& dt = trace->demand[by_band[r]];
        dt.band = bi;
        dt.qx = qx[r - g];
        dt.qy = qy[r - g];
      }
    }

    // Alg 6 l.6-8: scales Delta = 2^i a, i = ell (coarsest) .. 0. At each scale an independent
    // uniform permutation and j; then refine every part by its center labels.
    // 【解釈】I-3 (Q-41): initial part of each vertex = its quotient component (above the coarsest
    // scale). 【解釈】I-6: draws from the coarsest scale to the finest.
    std::vector<std::uint32_t> prev_part(comp.label.begin(), comp.label.end());
    std::vector<Vertex> prev_node(comp.count, kNoVertex);
    std::vector<Vertex> tparent;            // partition forest of this band
    std::vector<std::uint32_t> node_table;  // table index of each part node
    std::vector<Vertex> verts(nq);
    std::iota(verts.begin(), verts.end(), Vertex{0});
    for (unsigned si = 0; si <= k.ell; ++si) {
      const unsigned i = k.ell - si;
      const std::vector<std::uint32_t> perm = model::random_permutation(nq, rng);
      const std::vector<std::uint64_t> j = draw_j(k.eta, rng);
      const X threshold = ar.radius_threshold(a, i, k.eta, j);
      const std::vector<Vertex> center =
          first_covering_centers<W>(qg, perm, threshold, k.eta, ar, ops, cover, cc);
      // Group equal (parent part, center) pairs by stable counting sorts (l.1503-1505).
      // 【解釈】I-5: count arrays bounded by |V(Q_a)|, not per component.
      const std::array<std::size_t, 2> ranges{prev_node.size(), nq};
      const std::vector<Vertex> sorted =
          radix_sort_stable(verts, ranges, [&](Vertex z, std::size_t c) -> std::size_t {
            return c == 0 ? prev_part[z] : center[z];
          });
      std::vector<std::uint32_t> part(nq);
      std::vector<Vertex> node;
      for (std::size_t r = 0; r < nq; ++r) {
        const Vertex z = sorted[r];
        const bool fresh = r == 0 || prev_part[z] != prev_part[sorted[r - 1]] ||
                           center[z] != center[sorted[r - 1]];
        if (fresh) {
          // A new part: its own index in the global table, storing 4 Delta + a (Alg 6 l.9).
          node.push_back(static_cast<Vertex>(tparent.size()));
          tparent.push_back(prev_node[prev_part[z]]);
          node_table.push_back(static_cast<std::uint32_t>(res.table.size()));
          res.table.push_back({a, i});
        }
        part[z] = static_cast<std::uint32_t>(node.size() - 1);
      }
      // Lemma 6.7(1) / Q-41: at the top scale (2^l >= t T^2) the first center covers its whole
      // quotient component, so each component is exactly one part.
      if (si == 0 && node.size() != comp.count) {
        throw std::logic_error("estimate: top-scale part differs from its quotient component");
      }
      if (bt != nullptr) {
        auto& sc = bt->scales.emplace_back();
        sc.scale = i;
        sc.order = perm;
        sc.j = j;
        sc.center = center;
        sc.table_index.resize(nq);
        for (Vertex z = 0; z < nq; ++z) sc.table_index[z] = node_table[node[part[z]]];
      }
      prev_part = std::move(part);
      prev_node = std::move(node);
    }
    // One distinct leaf per vertex below the finest parts; LCA preprocessing (Lemma 3.9).
    std::vector<Vertex> leaf(nq);
    for (Vertex z = 0; z < nq; ++z) {
      leaf[z] = static_cast<Vertex>(tparent.size());
      tparent.push_back(prev_node[prev_part[z]]);
    }
    const ForestQueries parts(tparent);
    // Alg 6 l.10: each demand of the band gets the index at the LCA of its leaves.
    for (std::size_t r = g; r < g_end; ++r) {
      const auto l = parts.lca(leaf[qx[r - g]], leaf[qy[r - g]]);
      if (!l.has_value() || *l >= node_table.size()) {
        throw std::logic_error("estimate: demand leaves have no common part");
      }
      res.index[by_band[r]] = node_table[*l];
    }
    g = g_end;
  }

  cnt["estimate.bands_processed"] += bands_processed;
  cnt["estimate.quotient_vertices"] += quotient_vertices;
  cnt["estimate.quotient_edges"] += quotient_edges;
  cnt["estimate.vertex_scales"] += vertex_scales;
  cnt["estimate.center_searches"] += cc.searches;
  cnt["estimate.settles"] += cc.settles;
  cnt["estimate.adjacency_scans"] += cc.scans;
  cnt["estimate.table_entries"] += res.table.size();
  std::uint64_t max_windows = 0;
  for (const auto c : windows) max_windows = std::max<std::uint64_t>(max_windows, c);
  auto& mw = cnt["estimate.max_windows_per_edge"];
  mw = std::max(mw, max_windows);
  return res;
}

// ------------------------------------------------------------------ invariant checker

template <class W>
std::vector<std::string> check_estimate(std::size_t t, std::span<const Edge<W>> sample,
                                        std::span<const Demand> demands,
                                        const EstimateResult<W>& result,
                                        const std::vector<std::vector<W>>& dist) {
  std::vector<std::string> v;
  if (result.index.size() != demands.size()) {
    v.push_back("index has " + std::to_string(result.index.size()) + " entries for " +
                std::to_string(demands.size()) + " demands");
    return v;
  }
  if (dist.size() != t) {
    v.push_back("reference distance matrix has the wrong size");
    return v;
  }
  std::vector<W> weights;
  weights.reserve(sample.size());
  for (const auto& e : sample) weights.push_back(e.weight);
  std::sort(weights.begin(), weights.end());
  for (std::size_t i = 0; i < result.table.size(); ++i) {
    if (!std::binary_search(weights.begin(), weights.end(), result.table[i].anchor)) {
      v.push_back("table entry " + std::to_string(i) + " has an anchor that is no sample weight");
    }
  }
  for (std::size_t q = 0; q < demands.size(); ++q) {
    const auto [x, y] = demands[q];
    if (x >= t || y >= t || dist[x].size() != t) {
      v.push_back("demand " + std::to_string(q) + " is out of range");
      continue;
    }
    const W d = dist[x][y];
    const bool disconnected = d == WeightTraits<W>::infinity();
    const std::uint32_t idx = result.index[q];
    const std::string who =
        "demand " + std::to_string(q) + " (" + std::to_string(x) + "," + std::to_string(y) + ")";
    if (idx == kDisconnected) {
      if (!disconnected) v.push_back(who + " is connected but flagged as disconnected");
      continue;
    }
    if (disconnected) {
      v.push_back(who + " is disconnected but has a table index");
      continue;
    }
    if (idx >= result.table.size()) {
      v.push_back(who + " has an index outside the table");
      continue;
    }
    if (!bound_at_least(result.table[idx], d)) v.push_back(who + ": Z < d_S (Lemma 6.7(1))");
  }
  return v;
}

// ------------------------------------------------------------------ instantiations

template bool bound_at_least<std::uint64_t>(const TableBound<std::uint64_t>&, std::uint64_t);
template bool scaled_at_least<std::uint64_t>(std::uint64_t, unsigned,
                                             const TableBound<std::uint64_t>&);
template bool pow2_multiple_at_least<std::uint64_t>(std::uint64_t, unsigned, std::uint64_t);
template bool within_radius<std::uint64_t>(std::uint64_t, std::uint64_t, unsigned, unsigned,
                                           std::span<const std::uint64_t>);
template EstimateResult<std::uint64_t> estimate<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>, std::span<const Demand>, model::WordRng&,
    const EstimateParams&, EstimateTrace<std::uint64_t>*);
template std::vector<std::string> check_estimate<std::uint64_t>(
    std::size_t, std::span<const Edge<std::uint64_t>>, std::span<const Demand>,
    const EstimateResult<std::uint64_t>&, const std::vector<std::vector<std::uint64_t>>&);
template bool bound_at_least<double>(const TableBound<double>&, double);
template bool scaled_at_least<double>(double, unsigned, const TableBound<double>&);
template bool pow2_multiple_at_least<double>(double, unsigned, double);
template bool within_radius<double>(double, double, unsigned, unsigned,
                                    std::span<const std::uint64_t>);
template EstimateResult<double> estimate<double>(std::size_t, std::span<const Edge<double>>,
                                                 std::span<const Demand>, model::WordRng&,
                                                 const EstimateParams&, EstimateTrace<double>*);
template std::vector<std::string> check_estimate<double>(std::size_t, std::span<const Edge<double>>,
                                                         std::span<const Demand>,
                                                         const EstimateResult<double>&,
                                                         const std::vector<std::vector<double>>&);

}  // namespace sssp::paper::pivots
