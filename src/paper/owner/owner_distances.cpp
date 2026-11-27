// Sec 7.4 and Algorithm 7 (arXiv 2611.34911v1, txt l.1780-1963). See owner_distances.hpp for
// the statement, the units and the interpretations.

#include "paper/owner/owner_distances.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "paper/hierarchy/child_order.hpp"
#include "paper/hierarchy/evaluate.hpp"
#include "paper/hierarchy/normalize.hpp"
#include "paper/local/atom_replacement.hpp"
#include "paper/local/saturating.hpp"
#include "paper/primitives/forest_queries.hpp"
#include "paper/primitives/split_findmin.hpp"

namespace sssp::paper::owner {

namespace {

constexpr std::uint32_t kNil = 0xffffffffu;

void count([[maybe_unused]] Diagnostics* diag, [[maybe_unused]] const char* name,
           [[maybe_unused]] std::uint64_t k = 1) {
  SSSP_COUNT(if (diag != nullptr) diag->counters[name] += k);
}

template <class W>
constexpr W inf() noexcept {
  return WeightTraits<W>::infinity();
}

/// a + b for finite scaled lengths; std::overflow_error if the sum does not fit W.
template <class W>
W add_fit(W a, W b) {
  const W s = sat::add(a, b);
  if (s == inf<W>()) throw std::overflow_error("owner_distances: a scaled length does not fit");
  return s;
}

/// x / 2^7 by seven halvings (see the units note in owner_distances.hpp). Every value passed
/// here is a sum of scaled lengths: for integers its low seven bits are zero.
template <class W>
W calendar_unscale(W x) {
  if (x == inf<W>()) return x;
  if constexpr (std::is_floating_point_v<W>) {
    for (unsigned i = 0; i < kCalendarScaleLog2; ++i) x /= W{2};
    return x;
  } else {
    if ((x & ((W{1} << kCalendarScaleLog2) - 1)) != 0) {
      throw std::logic_error("owner_distances: a scaled distance is not a multiple of 2^7");
    }
    return x >> kCalendarScaleLog2;
  }
}

/// Alg 7 l.5: a DFS preorder of T' (children in the BFS order of aux.tree_order).
template <class W>
std::vector<std::uint32_t> tree_preorder(const AuxGraph<W>& aux) {
  const std::size_t na = aux.num_vertices();
  std::vector<std::uint32_t> begin(na + 1, 0);
  for (std::size_t k = 1; k < na; ++k) ++begin[aux.tree_parent[aux.tree_order[k]] + 1];
  for (std::size_t v = 0; v < na; ++v) begin[v + 1] += begin[v];
  std::vector<Vertex> child(na == 0 ? 0 : na - 1);
  std::vector<std::uint32_t> at(begin.begin(), begin.end() - 1);
  for (std::size_t k = 1; k < na; ++k) {
    const Vertex v = aux.tree_order[k];
    child[at[aux.tree_parent[v]]++] = v;
  }
  std::vector<std::uint32_t> pre(na, kNil);
  std::vector<Vertex> stack{kAuxRoot};
  std::uint32_t next = 0;
  while (!stack.empty()) {
    const Vertex v = stack.back();
    stack.pop_back();
    pre[v] = next++;
    for (std::uint32_t i = begin[v + 1]; i > begin[v]; --i) stack.push_back(child[i - 1]);
  }
  if (next != na) throw std::invalid_argument("owner_distances: T' does not span G'");
  return pre;
}

/// Stable counting sort of `items` by key(item) < buckets (integer sorting, Alg 7 l.5).
template <class Key>
std::vector<std::uint32_t> counting_sort(const std::vector<std::uint32_t>& items,
                                         std::size_t buckets, Key key) {
  std::vector<std::uint32_t> start(buckets + 1, 0);
  for (const std::uint32_t x : items) ++start[key(x) + 1];
  for (std::size_t b = 0; b < buckets; ++b) start[b + 1] += start[b];
  std::vector<std::uint32_t> out(items.size());
  for (const std::uint32_t x : items) out[start[key(x)]++] = x;
  return out;
}

/// Lemma 7.7 check hook: Scan calls per local nominal slot (node, slot start).
template <class W>
class SlotScanCounter final : public EvalHooks<W> {
 public:
  void on_scan_begin(std::uint32_t node, std::span<const StateId> /*states*/, W /*lo*/, W /*hi*/,
                     W slot_lo, W /*slot_hi*/, W /*deadline*/) override {
    const std::uint32_t c = ++scans_[{node, slot_lo}];
    max_ = std::max(max_, c);
  }
  [[nodiscard]] std::uint32_t max_scans() const noexcept { return max_; }

 private:
  std::map<std::pair<std::uint32_t, W>, std::uint32_t> scans_;
  std::uint32_t max_ = 0;
};

/// An interatom exit u -> head of scaled length w'' (l.1818: listed at its tail u).
template <class W>
struct Exit {
  Vertex head;
  W length;
};

/// A contiguous run of unregistered atoms: family positions [first, last] (l.1801-1802).
struct Run {
  std::uint32_t first;
  std::uint32_t last;
};

/// The state of one Alg 7 run (l.6-14).
template <class W>
class OwnerDistances {
 public:
  OwnerDistances(const AuxGraph<W>& aux, Diagnostics* diag, const OwnerDistancesOptions& opt,
                 OwnerDistancesResult<W>& res)
      : aux_(aux),
        diag_(diag),
        opt_(opt),
        res_(res),
        cal_(res.calendars),
        atoms_(res.atoms),
        fam_(res.families),
        na_(aux.num_vertices()),
        num_atoms_(static_cast<std::uint32_t>(res.atoms.num_atoms())) {}

  void run() {
    scale_lengths();
    activation_lists();
    edge_scan();
    local_hierarchies();
    shared_layout();
    initial_offers();
    // l.10-13: markers in topological order.
    processed_.assign(cal_.marker_value.size(), 0);
    for (const std::uint32_t m : cal_.marker_order) {
      count(diag_, "owner_distances.markers");
      processed_[m] = 1;
      for (std::uint32_t r = act_begin_[m]; r < act_begin_[m + 1]; ++r) {
        const std::uint32_t f = act_family_[r];
        const auto& wc = cal_.calendars[fam_.calendar[f]];
        const W lo = cal_.marker_value[wc.markers[act_cell_[r]]];
        const W hi = cal_.marker_value[wc.markers[act_cell_[r] + 1]];
        register_due_atoms(f, lo, hi, m);  // l.12
        visit_active_atoms(f, lo, hi, m);  // l.13
      }
    }
    // l.14: the least emitted value at each vertex (zero at the root).
    for (Vertex v = 0; v < na_; ++v) {
      if (v != kAuxRoot && res_.scaled_distance[v] == inf<W>()) {
        throw std::logic_error(
            "owner_distances: a vertex was never emitted (Lemma 7.9 violated; for double see "
            "Q-55)");
      }
    }
    if (opt_.count_slot_scans) {
      for (const auto& h : hooks_)
        res_.max_slot_scans = std::max(res_.max_slot_scans, h->max_scans());
    }
  }

 private:
  // ------------------------------------------------------------------ setup

  /// Every length used below in calendar units (2^7 x, Lemma 7.3 l.1723-1724).
  void scale_lengths() {
    rs_.resize(na_);
    for (Vertex v = 0; v < na_; ++v) rs_[v] = calendar_scale(aux_.owner_distance[v]);
    ls_.resize(num_atoms_);
    for (std::uint32_t a = 0; a < num_atoms_; ++a) ls_[a] = calendar_scale(atoms_.length[a]);
    local_.assign(na_, kNil);
    for (std::uint32_t a = 0; a < num_atoms_; ++a) {
      const auto mem = atoms_.members(a);
      for (std::uint32_t i = 0; i < mem.size(); ++i) local_[mem[i]] = i;
    }
  }

  /// Alg 7 l.5 (l.1796-1800): one family record at the starting marker of every cell of the
  /// family's (C, q) calendar from cell 0 through the last cell with start <= horizon, or to the
  /// calendar's end.
  void activation_lists() {
    const std::size_t nm = cal_.marker_value.size();
    act_begin_.assign(nm + 1, 0);
    const auto cells_of = [&](std::uint32_t f) {
      const auto& wc = cal_.calendars[fam_.calendar[f]];
      std::uint32_t k = 0;
      while (k < wc.num_cells() && !(fam_.horizon[f] < cal_.marker_value[wc.markers[k]])) ++k;
      return k;
    };
    std::vector<std::uint32_t> cells(fam_.num_families());
    for (std::uint32_t f = 0; f < fam_.num_families(); ++f) {
      cells[f] = cells_of(f);
      const auto& wc = cal_.calendars[fam_.calendar[f]];
      for (std::uint32_t k = 0; k < cells[f]; ++k) ++act_begin_[wc.markers[k] + 1];
    }
    for (std::size_t m = 0; m < nm; ++m) act_begin_[m + 1] += act_begin_[m];
    act_family_.assign(act_begin_[nm], 0);
    act_cell_.assign(act_begin_[nm], 0);
    std::vector<std::uint32_t> at(act_begin_.begin(), act_begin_.end() - 1);
    for (std::uint32_t f = 0; f < fam_.num_families(); ++f) {
      const auto& wc = cal_.calendars[fam_.calendar[f]];
      for (std::uint32_t k = 0; k < cells[f]; ++k) {
        const std::uint32_t p = at[wc.markers[k]]++;
        act_family_[p] = f;
        act_cell_[p] = k;
      }
    }
    res_.activation_records = act_begin_[nm];
    count(diag_, "owner_distances.activation_records", act_begin_[nm]);
  }

  /// Alg 7 l.6 (l.1816-1818, l.1823-1824): one edge scan sorts every edge of G' into root
  /// edges, internal edges of one atom (chords included) and interatom exits at their tails;
  /// between graph pieces only the direction from the smaller range is kept.
  void edge_scan() {
    const auto edges = aux_.graph.edges();
    std::vector<Vertex> tail_only(edges.size(), kNoVertex);
    for (const auto& arc : cal_.interpiece_arcs) tail_only[arc.aux_edge] = arc.tail;
    internal_.assign(num_atoms_, {});
    internal_id_.assign(edges.size(), kNoEdge);
    exit_begin_.assign(na_ + 1, 0);
    std::vector<std::pair<Vertex, Exit<W>>> exits;
    for (EdgeId e = 0; e < edges.size(); ++e) {
      const auto& ed = edges[e];
      const W w = calendar_scale(ed.weight);
      if (ed.u == kAuxRoot || ed.v == kAuxRoot) {
        if (ed.u == ed.v) continue;  // a loop at a (build_aux_graph removes them)
        root_edges_.push_back({ed.u == kAuxRoot ? ed.v : ed.u, w});
        continue;
      }
      const std::uint32_t au = atoms_.atom[ed.u];
      const std::uint32_t av = atoms_.atom[ed.v];
      if (au == av) {
        internal_id_[e] = static_cast<EdgeId>(internal_[au].size());
        internal_[au].push_back({local_[ed.u], local_[ed.v], w});
        continue;
      }
      for (const auto& [x, y] : {std::pair{ed.u, ed.v}, std::pair{ed.v, ed.u}}) {
        if (tail_only[e] != kNoVertex && tail_only[e] != x) {
          count(diag_, "owner_distances.exits_wrong_direction");
          continue;
        }
        exits.push_back({x, Exit<W>{y, w}});
        ++exit_begin_[x + 1];
      }
    }
    for (std::size_t v = 0; v < na_; ++v) exit_begin_[v + 1] += exit_begin_[v];
    exit_.resize(exits.size());
    std::vector<std::uint32_t> at(exit_begin_.begin(), exit_begin_.end() - 1);
    for (const auto& [x, ex] : exits) exit_[at[x]++] = ex;
    count(diag_, "owner_distances.root_edges", root_edges_.size());
    count(diag_, "owner_distances.exits", exits.size());
  }

  /// Alg 7 l.6-7: LocalHierarchy on every G'[Y], capped at the width q (Lemma 5.13; the
  /// discovery tree has every edge < 2q, l.1829-1831), normalised (Lemma 4.3) and ordered
  /// (Lemma 4.5). Alg 4 records its own base-case fallback.
  void local_hierarchies() {
    nh_.reserve(num_atoms_);
    native_.resize(num_atoms_);
    std::size_t internal_edges = 0;
    for (std::uint32_t a = 0; a < num_atoms_; ++a) {
      const std::size_t k = atoms_.members(a).size();
      internal_edges += internal_[a].size();
      const auto h = CsrGraph<W>::from_edges(k, std::move(internal_[a]));
      std::vector<EdgeId> tree;
      for (std::uint32_t t = atoms_.tree_offsets[a]; t < atoms_.tree_offsets[a + 1]; ++t) {
        tree.push_back(internal_id_[atoms_.tree_edges[t]]);
      }
      const W q = cal_.width[atoms_.representative[a]];
      auto capped = cap_local_hierarchy(h, std::span<const EdgeId>(tree), q, diag_, opt_.local);
      res_.radial_atom[a] = local_construction_applies(k, opt_.local) ? 1 : 0;
      native_[a] = std::move(capped.local.native_entry);
      NormalizedHierarchy<W> n = normalize(capped.local.hierarchy, diag_);
      order_children(n, diag_);
      nh_.push_back(std::move(n));
    }
    internal_.clear();
    count(diag_, "owner_distances.atoms", num_atoms_);
    count(diag_, "owner_distances.internal_edges", internal_edges);
  }

  /// Alg 7 l.8: one leaf layout in family order, then representative preorder, then local
  /// child order (each evaluator's own layout); all keys infinite; a split at every family
  /// boundary.
  void shared_layout() {
    base_.assign(num_atoms_, 0);
    std::size_t total = 0;
    for (const std::uint32_t a : fam_.order) {
      base_[a] = static_cast<std::uint32_t>(total);
      total += nh_[a].forest.num_states();
    }
    count(diag_, "owner_distances.atom_states", total);
    atom_at_.assign(total, kNil);
    for (std::uint32_t a = 0; a < num_atoms_; ++a) {
      const std::size_t n = nh_[a].forest.num_states();
      for (std::size_t i = 0; i < n; ++i) atom_at_[base_[a] + i] = a;
    }
    sf_ = std::make_unique<SplitFindmin<W>>(total);
    for (std::uint32_t f = 0; f < fam_.num_families(); ++f) {
      sf_->split(base_[fam_.order[fam_.offsets[f]]]);
    }
    ev_.reserve(num_atoms_);
    for (std::uint32_t a = 0; a < num_atoms_; ++a) {
      EvalOptions<W> eo;
      if (opt_.count_slot_scans) {
        hooks_.push_back(std::make_unique<SlotScanCounter<W>>());
        eo.hooks = hooks_.back().get();
      }
      ev_.emplace_back(nh_[a], *sf_, base_[a], eo, diag_);
    }
    runs_.resize(fam_.num_families());
    active_.resize(fam_.num_families());
    for (std::uint32_t f = 0; f < fam_.num_families(); ++f) {
      runs_[f].push_back({fam_.offsets[f], fam_.offsets[f + 1] - 1});
    }
  }

  /// Alg 7 l.9: R_v at every native inward copy, and every direct offer from the root (Q-47).
  void initial_offers() {
    for (Vertex v = 0; v < na_; ++v) {
      if (v == kAuxRoot) continue;
      const std::uint32_t a = atoms_.atom[v];
      ev_[a].initial_offer(native_[a][local_[v]], rs_[v]);
      count(diag_, "owner_distances.initial_offers");
    }
    for (const auto& [v, w] : root_edges_) {
      const std::uint32_t a = atoms_.atom[v];
      ev_[a].initial_offer(native_[a][local_[v]], w);
      count(diag_, "owner_distances.root_offers");
    }
  }

  // ------------------------------------------------------------------ markers

  /// Alg 7 l.12 (l.1801-1815): register due atoms by splitting the family's live runs.
  void register_due_atoms(std::uint32_t f, W lo, W hi, std::uint32_t marker) {
    std::vector<Run> work = std::move(runs_[f]);
    runs_[f].clear();
    while (!work.empty()) {
      const Run run = work.back();
      work.pop_back();
      count(diag_, "owner_distances.run_checks");
      const SfmMin<W> m = sf_->query(base_[fam_.order[run.first]]);
      if (m.key != inf<W>() && m.key < lo) {
        throw std::logic_error(
            "owner_distances: a family run has a key below the current cell (Lemma 7.8 "
            "violated)");
      }
      if (!(m.key < hi)) {  // not due in [lo, hi) (an infinite key is never due)
        count(diag_, "owner_distances.run_checks_failed");
        runs_[f].push_back(run);
        continue;
      }
      // l.1802: the run splits at the atom i containing its first due minimum z.
      const W z = m.key;
      const std::uint32_t i = atom_at_[m.index];
      const std::uint32_t k = fam_.position[i];
      sf_->split(base_[i]);
      sf_->split(static_cast<std::uint32_t>(base_[i] + nh_[i].forest.num_states()));
      count(diag_, "owner_distances.registrations");
      res_.registration[i] = z;
      // l.1825: the registration time z with L_Y sets the lasting inclusive deadline z + L_Y.
      ev_[i].set_deadline(sat::add(z, ls_[i]));
      active_[f].push_back(i);
      // l.1803-1805: a padded proposal at the representative of each new subrun's boundary atom.
      if (k > run.first) {
        propose(fam_.order[k - 1], z, ls_[i], fam_.gap[k - 1], hi, marker);
        work.push_back({run.first, k - 1});
      }
      if (k < run.last) {
        propose(fam_.order[k + 1], z, ls_[i], fam_.gap[k], hi, marker);
        work.push_back({k + 1, run.last});
      }
    }
  }

  /// u = max{b0, z + L_i + g(i, j)} at the native inward copy of v = rep(j), discarded if
  /// u > R_v or if it does not lower the pending key (l.1803-1805).
  void propose(std::uint32_t j, W z, W li, W gap, W b0, std::uint32_t marker) {
    count(diag_, "owner_distances.supplemental_proposals");
    const Vertex v = atoms_.representative[j];
    const W u = std::max(b0, sat::add(sat::add(z, li), gap));
    if (rs_[v] < u) {
      count(diag_, "owner_distances.supplemental_capped");
      return;
    }
    if (ev_[j].external_offer(native_[j][local_[v]], u)) {
      count(diag_, "owner_distances.supplemental_decreases");
      log_arrival(ArrivalKind::kSupplemental, kNoVertex, v, u, marker);
    }
  }

  /// Alg 7 l.13 (l.1835-1841): visit every registered, unfinished atom of the family in the
  /// cell, then relax the allowed exits of its emissions.
  void visit_active_atoms(std::uint32_t f, W lo, W hi, std::uint32_t marker) {
    auto& act = active_[f];
    std::size_t keep = 0;
    for (const std::uint32_t a : act) {
      const auto emitted = ev_[a].visit_cell(lo, hi);
      ++res_.atom_cells[a];
      count(diag_, "owner_distances.atom_cell_visits");
      const auto mem = atoms_.members(a);
      for (const auto& [state, value] : emitted) {
        const Vertex u = mem[nh_[a].forest.projection[state]];
        count(diag_, "owner_distances.emissions");
        if (opt_.record_emissions) res_.emissions.emplace_back(u, value);
        res_.scaled_distance[u] = std::min(res_.scaled_distance[u], value);
        relax_exits(u, value, marker);
      }
      if (!ev_[a].finished()) act[keep++] = a;
    }
    act.resize(keep);
  }

  /// l.1821-1823 and Lemma 7.8: suppress the exits of u if the emitted value exceeds R_u;
  /// discard an offer to v that exceeds R_v.
  void relax_exits(Vertex u, W value, std::uint32_t marker) {
    if (rs_[u] < value) {
      count(diag_, "owner_distances.exits_suppressed");
      return;
    }
    for (std::uint32_t p = exit_begin_[u]; p < exit_begin_[u + 1]; ++p) {
      const Exit<W>& ex = exit_[p];
      count(diag_, "owner_distances.exit_relaxations");
      const W offer = sat::add(value, ex.length);
      if (rs_[ex.head] < offer) {
        count(diag_, "owner_distances.exits_capped");
        continue;
      }
      const std::uint32_t b = atoms_.atom[ex.head];
      if (ev_[b].external_offer(native_[b][local_[ex.head]], offer)) {
        count(diag_, "owner_distances.interatom_decreases");
        log_arrival(ArrivalKind::kInteratom, u, ex.head, offer, marker);
      }
    }
  }

  /// Lemma 7.8 arrival log: the starting marker of the head's cell containing `value` must not
  /// have been processed yet.
  void log_arrival(ArrivalKind kind, Vertex tail, Vertex head, W value, std::uint32_t marker) {
    if (!opt_.record_arrivals) return;
    ArrivalRecord<W> r;
    r.kind = kind;
    r.tail = tail;
    r.head = head;
    r.value = value;
    r.marker = marker;
    const auto& wc = cal_.calendars[cal_.calendar[head]];
    const auto begin = wc.markers.begin();
    const auto it = std::upper_bound(begin, wc.markers.end(), value, [&](W x, std::uint32_t mk) {
      return x < cal_.marker_value[mk];
    });
    if (it != begin && it != wc.markers.end()) {
      r.destination = *(it - 1);
      r.on_time = processed_[r.destination] == 0;
    }
    res_.arrivals.push_back(r);
  }

  const AuxGraph<W>& aux_;
  Diagnostics* diag_;
  const OwnerDistancesOptions& opt_;
  OwnerDistancesResult<W>& res_;
  const Calendars<W>& cal_;
  const Atoms<W>& atoms_;
  const Families<W>& fam_;
  std::size_t na_;
  std::uint32_t num_atoms_;

  std::vector<W> rs_;                 ///< 2^7 R_v
  std::vector<W> ls_;                 ///< 2^7 L_Y
  std::vector<std::uint32_t> local_;  ///< vertex -> index among its atom's members
  std::vector<std::uint32_t> act_begin_, act_family_, act_cell_;
  std::vector<std::vector<Edge<W>>> internal_;
  std::vector<EdgeId> internal_id_;
  std::vector<std::pair<Vertex, W>> root_edges_;
  std::vector<std::uint32_t> exit_begin_;
  std::vector<Exit<W>> exit_;
  std::vector<NormalizedHierarchy<W>> nh_;
  std::vector<std::vector<StateId>> native_;
  std::vector<std::uint32_t> base_;
  std::vector<std::uint32_t> atom_at_;
  std::unique_ptr<SplitFindmin<W>> sf_;
  std::vector<std::unique_ptr<SlotScanCounter<W>>> hooks_;
  std::vector<PortionedEvaluator<W>> ev_;
  std::vector<std::vector<Run>> runs_;
  std::vector<std::vector<std::uint32_t>> active_;
  std::vector<std::uint8_t> processed_;
};

}  // namespace

template <class W>
Families<W> build_families(const AuxGraph<W>& aux, const Calendars<W>& cal, const Atoms<W>& atoms,
                           Diagnostics* diag) {
  const std::size_t na = aux.num_vertices();
  if (cal.tree_piece.size() != na || cal.calendar.size() != na || atoms.atom.size() != na ||
      aux.tree_parent.size() != na || aux.tree_order.size() != na) {
    throw std::invalid_argument("build_families: calendars or atoms do not belong to G'");
  }
  Families<W> fam;
  fam.preorder = tree_preorder(aux);
  const std::uint32_t num_atoms = static_cast<std::uint32_t>(atoms.num_atoms());
  // l.1782-1783: representatives in T' preorder (one counting pass over the preorder), then
  // stable counting sorts by width calendar and by tree piece: each (P, q) class is contiguous
  // and keeps the preorder.
  std::vector<std::uint32_t> by_vertex(na, kNil);
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    by_vertex[fam.preorder[atoms.representative[a]]] = a;
  }
  std::vector<std::uint32_t> items;
  items.reserve(num_atoms);
  for (const std::uint32_t a : by_vertex) {
    if (a != kNil) items.push_back(a);
  }
  const auto rep = [&](std::uint32_t a) { return atoms.representative[a]; };
  items = counting_sort(items, cal.calendars.size(),
                        [&](std::uint32_t a) { return cal.calendar[rep(a)]; });
  items = counting_sort(items, cal.tree_piece_root.size(),
                        [&](std::uint32_t a) { return cal.tree_piece[rep(a)]; });
  fam.order = items;
  fam.family.assign(num_atoms, kNoFamily);
  fam.position.assign(num_atoms, 0);
  fam.gap.assign(num_atoms, W{0});
  fam.offsets.push_back(0);
  for (std::uint32_t k = 0; k < num_atoms; ++k) {
    const std::uint32_t a = fam.order[k];
    const Vertex r = rep(a);
    if (k == 0 || cal.tree_piece[r] != fam.tree_piece.back() ||
        cal.calendar[r] != fam.calendar.back()) {
      if (k != 0) fam.offsets.push_back(k);
      fam.tree_piece.push_back(cal.tree_piece[r]);
      fam.calendar.push_back(cal.calendar[r]);
      fam.horizon.push_back(W{0});
    }
    fam.family[a] = static_cast<std::uint32_t>(fam.tree_piece.size() - 1);
    fam.position[a] = k;
    // l.1786-1787: horizon max_Y (R_rep(Y) + L_Y), scaled.
    const W reach = add_fit(calendar_scale(aux.owner_distance[r]), calendar_scale(atoms.length[a]));
    fam.horizon.back() = std::max(fam.horizon.back(), reach);
  }
  if (num_atoms != 0) fam.offsets.push_back(num_atoms);
  // l.1785: consecutive gaps from tree distances, (R_x - R_lca) + (R_y - R_lca) (Lemma 3.9).
  if (num_atoms > 1) {
    const ForestQueries fq(std::span<const Vertex>(aux.tree_parent));
    for (std::uint32_t f = 0; f + 1 < fam.offsets.size(); ++f) {
      for (std::uint32_t k = fam.offsets[f]; k + 1 < fam.offsets[f + 1]; ++k) {
        const Vertex x = rep(fam.order[k]);
        const Vertex y = rep(fam.order[k + 1]);
        const auto lca = fq.lca(x, y);
        if (!lca.has_value()) throw std::invalid_argument("build_families: T' is not connected");
        const Vertex l = *lca;
        const W rx = calendar_scale(aux.owner_distance[x]);
        const W ry = calendar_scale(aux.owner_distance[y]);
        const W rl = calendar_scale(aux.owner_distance[l]);
        fam.gap[k] = add_fit(rx - rl, ry - rl);
      }
    }
  }
  count(diag, "owner_distances.families", fam.num_families());
  return fam;
}

template <class W>
OwnerDistancesResult<W> owner_distances(const AuxGraph<W>& aux, Diagnostics* diag,
                                        const OwnerDistancesOptions& options) {
  count(diag, "owner_distances.calls");
  OwnerDistancesResult<W> res;
  res.calendars = build_calendars(aux, diag);                          // l.2
  res.atoms = find_atoms(aux, res.calendars, diag);                    // l.3
  res.families = build_families(aux, res.calendars, res.atoms, diag);  // l.4-5
  const std::size_t na = aux.num_vertices();
  const std::size_t num_atoms = res.atoms.num_atoms();
  res.scaled_distance.assign(na, inf<W>());
  res.scaled_distance[kAuxRoot] = W{0};
  res.registration.assign(num_atoms, inf<W>());
  res.atom_cells.assign(num_atoms, 0);
  res.radial_atom.assign(num_atoms, 0);
  if (num_atoms != 0) {
    OwnerDistances<W> alg(aux, diag, options, res);
    alg.run();
  }
  res.distance.resize(na);
  for (Vertex v = 0; v < na; ++v) res.distance[v] = calendar_unscale(res.scaled_distance[v]);
  return res;
}

#define SSSP_OWNER_DISTANCES_INSTANTIATE(W)                                                        \
  template Families<W> build_families<W>(const AuxGraph<W>&, const Calendars<W>&, const Atoms<W>&, \
                                         Diagnostics*);                                            \
  template OwnerDistancesResult<W> owner_distances<W>(const AuxGraph<W>&, Diagnostics*,            \
                                                      const OwnerDistancesOptions&);

SSSP_OWNER_DISTANCES_INSTANTIATE(std::uint64_t)
SSSP_OWNER_DISTANCES_INSTANTIATE(double)

#undef SSSP_OWNER_DISTANCES_INSTANTIATE

}  // namespace sssp::paper::owner
