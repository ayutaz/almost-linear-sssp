// Definition 6.1, Sec 6.1, Algorithm 5 NameRoutes(G_P), Lemma 6.8 (arXiv 2611.34911v1,
// txt l.1367-1635). See name_routes.hpp for the contract. Step comments cite "Alg 5 l.k".
//
// Interpretations of points the paper leaves open (marked 【解釈】 at the code):
//   N-1 (Q-44) "At least half the edges of U qualify" is decided as 2 |qualified| >= |U|.
//   N-2 Among parallel edges of equal least weight the representative is the least index.
//   N-3 Loops are rejected: Definition 6.1 instances are loopless and Sec 6 lets the caller
//       discard loops before forming the instance.
//   N-4 The representatives come from the shared primitive edge_pairs (in increasing edge id);
//       U lists them in increasing (min, max) endpoint order (stable radix sort, Lemma 3.1).
//   N-5 Dictionary records are sorted by (weight, record index) with a merge sort; equal
//       weights are coalesced (Alg 5 l.16).
//   N-6 (Q-40) The dictionary size is recorded ("name_routes.dictionary_size"); Lemma 6.8 states
//       O(t log^2(N+2)), its proof gives O(t log^2(t+2)).
//   N-7 A certificate walk is a concatenation of instance edges, so check_route_names tests
//       dist_{G_P}(x, y) <= P lambda_e.
#include "paper/pivots/name_routes.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>

#include "paper/model/ilog.hpp"
#include "paper/model/length.hpp"
#include "paper/model/params.hpp"
#include "paper/primitives/edge_pairs.hpp"

namespace sssp::paper::pivots {
namespace {

constexpr std::uint32_t kNoRecord = 0xffffffffu;
constexpr unsigned kMaxPLog2 = 256;

template <class W>
void validate(const PivotRoutingInstance<W>& g, const NameRoutesParams& p) {
  if (g.t >= (std::size_t{1} << 32)) throw std::length_error("name_routes: t >= 2^32");
  if (g.edges.size() >= kNoRecord) throw std::length_error("name_routes: too many edges");
  for (std::size_t i = 0; i < g.edges.size(); ++i) {
    const auto& e = g.edges[i];
    if (e.u >= g.t || e.v >= g.t) {
      throw std::out_of_range("name_routes: edge " + std::to_string(i) + " leaves 0..t-1");
    }
    if (e.u == e.v) {
      // 【解釈】N-3: reject loops (the caller discards them before forming the instance).
      throw std::invalid_argument("name_routes: the instance must be loopless (Def 6.1)");
    }
    if (!WeightTraits<W>::is_valid_weight(e.weight) || e.weight == WeightTraits<W>::infinity() ||
        !(W{0} < e.weight)) {
      throw std::invalid_argument("name_routes: weights must be positive and finite");
    }
  }
  if (p.sample_factor == 0) throw std::invalid_argument("name_routes: sample_factor = 0");
  if (p.residual_factor < p.sample_factor) {
    throw std::invalid_argument("name_routes: residual_factor < sample_factor");
  }
  if (p.retry_cap.has_value() && *p.retry_cap == 0) {
    throw std::invalid_argument("name_routes: retry_cap must be positive");
  }
  if (p.p_log2.has_value() && *p.p_log2 > kMaxPLog2) {
    throw std::invalid_argument("name_routes: p_log2 > 256");
  }
}

/// Comparison sort of record indices by (weight, index): bottom-up merge of adjacent runs
/// (Lemma 6.8(4): the names are sorted by comparison).
template <class W>
std::vector<std::uint32_t> sort_records(const std::vector<W>& rec, const model::LenOps<W>& ops) {
  const std::size_t n = rec.size();
  std::vector<std::uint32_t> a(n);
  std::iota(a.begin(), a.end(), 0U);
  std::vector<std::uint32_t> b(n);
  const auto less = [&](std::uint32_t x, std::uint32_t y) {
    return ops.less(model::EdgeKey<W>{rec[x], x}, model::EdgeKey<W>{rec[y], y});
  };
  for (std::size_t width = 1; width < n; width *= 2) {
    for (std::size_t lo = 0; lo < n; lo += 2 * width) {
      const std::size_t mid = std::min(lo + width, n);
      const std::size_t hi = std::min(lo + 2 * width, n);
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

void add_diagnostics(Diagnostics& into, const Diagnostics& from) {
  for (const auto& [k, v] : from.counters) {
    if (k == "estimate.max_windows_per_edge") {
      into.counters[k] = std::max(into.counters[k], v);
    } else {
      into.counters[k] += v;
    }
  }
  into.ops.vertex_scans += from.ops.vertex_scans;
  into.ops.edge_scans += from.ops.edge_scans;
  into.ops.relaxations += from.ops.relaxations;
  into.ops.comparisons += from.ops.comparisons;
  into.ops.additions += from.ops.additions;
  into.ops.queue_pushes += from.ops.queue_pushes;
  into.ops.queue_pops += from.ops.queue_pops;
  into.ops.decrease_keys += from.ops.decrease_keys;
  into.fallbacks.insert(into.fallbacks.end(), from.fallbacks.begin(), from.fallbacks.end());
}

}  // namespace

NameRoutesConstants name_routes_constants(std::size_t t) {
  if (t >= (std::size_t{1} << 32)) throw std::length_error("name_routes_constants: t >= 2^32");
  NameRoutesConstants c;
  c.u = model::ceil_log2(static_cast<std::uint64_t>(t) + 2);  // Alg 5 l.2
  c.g = 4 * c.u + 3;
  // Alg 5 l.3: least power of two >= 8 * 256 * g * u^2 (integer word arithmetic; fits easily).
  const std::uint64_t target = std::uint64_t{8} * 256 * c.g * static_cast<std::uint64_t>(c.u) * c.u;
  c.p_log2 = model::ceil_log2(target);
  return c;
}

template <class W>
RouteNames<W> name_routes(const PivotRoutingInstance<W>& g, model::WordRng& rng,
                          const NameRoutesParams& params) {
  validate(g, params);
  const std::size_t t = g.t;
  const std::size_t m = g.edges.size();
  RouteNames<W> out;
  out.name.assign(m, 0);
  const NameRoutesConstants c = name_routes_constants(t);
  out.p_log2 = params.p_log2.value_or(c.p_log2);
  auto& cnt = out.diag.counters;
  cnt["name_routes.trials"] += 0;
  if (m == 0) return out;  // also every t < 2: the loopless instance has no edges (l.1422)
  const model::LenOps<W> lops(&out.diag.ops);
  const auto len = [](W w) { return model::Len<W>::finite(w); };

  // Alg 5 l.4: the least-weight representative of each distinct endpoint pair (Lemma 6.8(3)), by
  // the shared primitive (edge_pairs.hpp: two stable counting sorts of the integer pairs, Lemma
  // 3.1, one counted length comparison per further edge of a pair). 【解釈】N-2: ties by the least
  // index. Loops were rejected above, so every edge has a representative.
  const LightestEdgePairs pairs = lightest_edge_pairs<W>(t, g.edges, &out.diag.ops);
  const std::vector<EdgeId>& rep = pairs.representative;
  // 【解釈】N-4: U is a set; it is listed in increasing (min, max) endpoint order (stable radix
  // sort, Lemma 3.1), which fixes the samples a given seed draws (Alg 5 l.6).
  const std::array<std::size_t, 2> ranges{t, t};
  std::vector<std::uint32_t> residual =  // U, as representative edge ids
      radix_sort_stable(pairs.kept, ranges, [&](std::uint32_t o, std::size_t k) -> std::size_t {
        const auto& e = g.edges[o];
        return k == 0 ? std::min(e.u, e.v) : std::max(e.u, e.v);
      });
  cnt["name_routes.representatives"] += residual.size();

  std::vector<W> records;                        // dictionary records (unsorted)
  std::vector<std::uint32_t> ptr(m, kNoRecord);  // representative -> record
  const std::uint64_t tt = t;
  const std::uint64_t sample_size = model::capped_mul(params.sample_factor, tt, ~std::uint64_t{0});
  const std::uint64_t limit = model::capped_mul(params.residual_factor, tt, ~std::uint64_t{0});
  std::uint64_t trials = 0;
  std::uint64_t accepted = 0;
  std::uint64_t rejected = 0;
  std::uint64_t round_trials = 0;
  std::uint64_t max_round = 0;
  std::uint64_t consecutive = 0;
  std::uint64_t max_consecutive = 0;
  std::vector<std::uint32_t> qualified_group;  // per residual position: group id or kNoRecord
  // Sec 2 l.303-304, Q-51 decision (b): the sample draws of l.6 have ranges |U|, |U|-1, ...; one
  // table over 1..|U| (U only shrinks) supplies their bit counts in every trial, built in
  // O(|U|) once per call. Nothing is tabulated if the loop does not run.
  const model::IlogTable lg(residual.size() > limit ? residual.size() : 0);
  cnt["name_routes.ilog_table_entries"] += lg.build_steps();

  // Alg 5 l.5: while |U| > 32t.
  while (residual.size() > limit) {
    ++trials;
    ++round_trials;
    // l.6: a uniformly random 16t-element subset S of U (partial Fisher-Yates; U is a set, so
    // reordering it is harmless).
    partial_shuffle(std::span<std::uint32_t>(residual), static_cast<std::size_t>(sample_size), rng,
                    lg);
    std::vector<Edge<W>> sample;
    sample.reserve(static_cast<std::size_t>(sample_size));
    for (std::size_t k = 0; k < sample_size; ++k) sample.push_back(g.edges[residual[k]]);
    // l.7: Estimate on S for the endpoint pairs of all edges in U.
    std::vector<Demand> demands;
    demands.reserve(residual.size());
    for (const std::uint32_t o : residual) demands.push_back({g.edges[o].u, g.edges[o].v});
    const EstimateResult<W> est = estimate<W>(t, sample, demands, rng, params.estimate);
    add_diagnostics(out.diag, est.diag);

    // l.8: group all connected pairs by returned index; no qualified edges yet.
    std::vector<std::uint32_t> connected;
    for (std::size_t p = 0; p < residual.size(); ++p) {
      if (est.index[p] != kDisconnected) connected.push_back(static_cast<std::uint32_t>(p));
    }
    const std::vector<std::uint32_t> grouped = counting_sort_stable(
        connected, est.table.size(), [&](std::uint32_t p) { return est.index[p]; });
    qualified_group.assign(residual.size(), kNoRecord);
    std::vector<W> proposal;  // a* of each qualifying group
    std::size_t num_qualified = 0;
    // l.9-12: per nonempty group with bound Z, the smallest a* with P a* >= Z by median pruning;
    // qualify every edge of the group with a_e >= a* and propose lambda_e = a*.
    for (std::size_t i = 0; i < grouped.size();) {
      const std::uint32_t idx = est.index[grouped[i]];
      std::size_t k = i;
      std::vector<W> keys;
      while (k < grouped.size() && est.index[grouped[k]] == idx) {
        keys.push_back(g.edges[residual[grouped[k]]].weight);
        ++k;
      }
      const TableBound<W> z = est.table[idx];
      const std::optional<W> astar = median_prune<W>(
          keys, [&](W key) { return scaled_at_least(key, out.p_log2, z); }, &out.diag.ops);
      if (astar.has_value()) {
        const auto gid = static_cast<std::uint32_t>(proposal.size());
        proposal.push_back(*astar);
        for (std::size_t r = i; r < k; ++r) {
          if (!lops.less(len(g.edges[residual[grouped[r]]].weight), len(*astar))) {
            qualified_group[grouped[r]] = gid;
            ++num_qualified;
          }
        }
      }
      i = k;
    }
    // l.13: accept iff at least half the edges of U qualify. 【解釈】N-1 (Q-44):
    // 2 |qualified| >= |U|.
    if (2 * static_cast<std::uint64_t>(num_qualified) >= residual.size()) {
      // l.14: one record of weight a* per qualifying group; point each qualified edge to it and
      // remove those edges from U.
      std::vector<std::uint32_t> record_of(proposal.size());
      for (std::size_t q = 0; q < proposal.size(); ++q) {
        record_of[q] = static_cast<std::uint32_t>(records.size());
        records.push_back(proposal[q]);
      }
      std::size_t keep = 0;
      for (std::size_t p = 0; p < residual.size(); ++p) {
        if (qualified_group[p] != kNoRecord) {
          ptr[residual[p]] = record_of[qualified_group[p]];
        } else {
          residual[keep++] = residual[p];
        }
      }
      residual.resize(keep);
      ++accepted;
      max_round = std::max(max_round, round_trials);
      round_trials = 0;
      consecutive = 0;
    } else {
      // A rejected round discards its proposals; the next trial draws a fresh sample with fresh
      // estimator randomness (l.1421).
      ++rejected;
      ++consecutive;
      max_consecutive = std::max(max_consecutive, consecutive);
      if (params.retry_cap.has_value() && consecutive >= *params.retry_cap) {
        // Deviation D-10 (not in the paper): stop and let the residual certify itself below.
        cnt["name_routes.retry_cap_hits"] += 1;
        out.diag.fallbacks.push_back(FallbackRecord{
            "name_routes: retry cap of " + std::to_string(*params.retry_cap) +
                " consecutive rejected trials reached (deviation D-10, not in the paper); the "
                "residual edges certify their own weight",
            tt, residual.size()});
        break;
      }
    }
  }

  // l.15: a record with its own weight for each remaining representative.
  for (const std::uint32_t o : residual) {
    ptr[o] = static_cast<std::uint32_t>(records.size());
    records.push_back(g.edges[o].weight);
  }
  cnt["name_routes.terminal_edges"] += residual.size();
  cnt["name_routes.dictionary_records"] += records.size();

  // l.16: sort the records by weight, coalesce ties, and map the original edges through the
  // record pointers of their representatives. 【解釈】N-5: merge sort by (weight, record).
  const std::vector<std::uint32_t> order = sort_records(records, lops);
  std::vector<std::uint32_t> new_index(records.size());
  for (const std::uint32_t r : order) {
    if (out.dictionary.empty() || lops.less(len(out.dictionary.back()), len(records[r]))) {
      out.dictionary.push_back(records[r]);
    }
    new_index[r] = static_cast<std::uint32_t>(out.dictionary.size() - 1);
  }
  for (std::size_t o = 0; o < m; ++o) {
    const std::uint32_t r = ptr[rep[o]];
    if (r == kNoRecord) throw std::logic_error("name_routes: representative without a record");
    out.name[o] = new_index[r];
  }

  cnt["name_routes.trials"] += trials;
  cnt["name_routes.accepted_rounds"] += accepted;
  cnt["name_routes.rejected_trials"] += rejected;
  cnt["name_routes.max_trials_per_round"] = max_round;
  cnt["name_routes.max_consecutive_rejections"] = max_consecutive;
  cnt["name_routes.dictionary_size"] += out.dictionary.size();
  return out;
}

template <class W>
std::vector<std::string> check_route_names(const PivotRoutingInstance<W>& g,
                                           const RouteNames<W>& names,
                                           const std::vector<std::vector<W>>& dist) {
  std::vector<std::string> v;
  if (names.name.size() != g.edges.size()) {
    v.push_back("name has " + std::to_string(names.name.size()) + " entries for " +
                std::to_string(g.edges.size()) + " edges");
    return v;
  }
  if (dist.size() != g.t) {
    v.push_back("reference distance matrix has the wrong size");
    return v;
  }
  std::vector<W> weights;
  weights.reserve(g.edges.size());
  for (const auto& e : g.edges) weights.push_back(e.weight);
  std::sort(weights.begin(), weights.end());
  for (std::size_t i = 0; i < names.dictionary.size(); ++i) {
    if (i > 0 && !(names.dictionary[i - 1] < names.dictionary[i])) {
      v.push_back("dictionary is not strictly increasing at " + std::to_string(i));
    }
    if (!std::binary_search(weights.begin(), weights.end(), names.dictionary[i])) {
      v.push_back("dictionary entry " + std::to_string(i) + " is no input weight");
    }
  }
  for (std::size_t o = 0; o < g.edges.size(); ++o) {
    const auto& e = g.edges[o];
    const std::string who =
        "edge " + std::to_string(o) + " (" + std::to_string(e.u) + "," + std::to_string(e.v) + ")";
    const std::uint32_t idx = names.name[o];
    if (idx >= names.dictionary.size()) {
      v.push_back(who + " has a name outside the dictionary");
      continue;
    }
    const W lambda = names.dictionary[idx];
    if (e.weight < lambda) v.push_back(who + ": lambda_e > a_e");
    if (dist[e.u].size() != g.t) {
      v.push_back(who + ": reference row has the wrong size");
      continue;
    }
    const W d = dist[e.u][e.v];
    if (d == WeightTraits<W>::infinity() || !pow2_multiple_at_least(lambda, names.p_log2, d)) {
      v.push_back(who + ": dist_GP(x,y) > P lambda_e (Definition 6.1)");
    }
  }
  return v;
}

template RouteNames<std::uint64_t> name_routes<std::uint64_t>(
    const PivotRoutingInstance<std::uint64_t>&, model::WordRng&, const NameRoutesParams&);
template std::vector<std::string> check_route_names<std::uint64_t>(
    const PivotRoutingInstance<std::uint64_t>&, const RouteNames<std::uint64_t>&,
    const std::vector<std::vector<std::uint64_t>>&);
template RouteNames<double> name_routes<double>(const PivotRoutingInstance<double>&,
                                                model::WordRng&, const NameRoutesParams&);
template std::vector<std::string> check_route_names<double>(
    const PivotRoutingInstance<double>&, const RouteNames<double>&,
    const std::vector<std::vector<double>>&);

}  // namespace sssp::paper::pivots
