#pragma once
// Shared core of the Sec 5.4 assembly (arXiv 2611.34911v1, Sec 5.4 l.1182-1243, Alg 4 l.10-12,
// proof of Lemma 5.10 l.1265-1300), task P5-LOCAL-UNIFY (Q-59, ADR 0003, Q-04).
//
// One implementation of the support grid, the simultaneous merges, the links and the child
// order, with the constants of the two scale tests passed in:
//   * Mode A (paper/local/assemble.hpp) calls assemble() with Rule::paper(N): c = 10 of the
//     alignment test (l.1206) and the support end N^16 sigma (l.1216). These are the defaults;
//     nothing else in src/paper/ passes another rule.
//   * Mode B (experimental/paper_experimental/local_assembly.hpp) passes the Lambda formulas of
//     deviation D-8, or sets Rule::exact_height for deviation D-11.
//   * Sec 8.2 (outer scales) uses the same grid with another unit (the name weight a) and
//     another number of octaves (F^3 = 2^(3f)), research note 04 §2.5 and §4 M-6 / M-8. The
//     pieces it needs are callable on their own: key_supports + build_support_grid,
//     doubling_points, aligned_by_medians, ClassMerger, link_and_order. Sec 8.2 itself is not
//     implemented here.
//
// The step numbers in the comments are those of paper/local/assemble.hpp.
//
// Model of computation (Sec 2): lengths are copied, compared and added. Every scale factor is
// applied by integer multiplications (binary expansion) and doublings that saturate at
// WeightTraits<W>::infinity() (paper/local/saturating.hpp, Q-38, D-13).

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "paper/hierarchy/ordered_hierarchy.hpp"
#include "paper/local/assemble.hpp"
#include "paper/local/owner_decomposition.hpp"
#include "paper/local/radial_states.hpp"
#include "paper/local/saturating.hpp"
#include "paper/local/threshold_joins.hpp"
#include "paper/primitives/dsu.hpp"
#include "sssp/graph.hpp"
#include "sssp/result.hpp"

namespace sssp::paper::assembly {

// ------------------------------------------------------------------ scale factors

/// A positive integer factor applied to a length: a fixed sequence of integer multiplications
/// and doublings, performed in the order given (for double the order fixes the rounding of the
/// integer multiplications; doublings are exact). Examples: N^10 is int_power(N, 10); F^3 =
/// 2^(3f) is octaves(3 * f); D-8's Lambda = 2N * 3hK^2 is
/// LengthScale{}.times(2).times(N).times(3).times(h).doubled(2 * log2 K).
class LengthScale {
 public:
  /// The factor 1.
  LengthScale() = default;

  /// n^e: e integer multiplications by n (sat::mul_int_pow).
  [[nodiscard]] static LengthScale int_power(std::uint64_t n, unsigned e) {
    return LengthScale{}.times(n, e);
  }
  /// 2^k: k doublings (sat::mul_pow2). One octave is one doubling.
  [[nodiscard]] static LengthScale octaves(unsigned k) { return LengthScale{}.doubled(k); }
  /// A factor that does not fit any length type: apply() returns infinity.
  [[nodiscard]] static LengthScale saturated() {
    LengthScale s;
    s.saturated_ = true;
    return s;
  }

  /// Followed by e multiplications by the integer n.
  LengthScale& times(std::uint64_t n, unsigned e = 1) {
    if (e > 0) steps_.push_back({n, e, false});
    return *this;
  }
  /// Followed by k doublings.
  LengthScale& doubled(unsigned k) {
    if (k > 0) steps_.push_back({0, k, true});
    return *this;
  }

  /// x times the factor, saturating at WeightTraits<W>::infinity().
  template <class W>
  [[nodiscard]] W apply(W x) const {
    if (saturated_) return WeightTraits<W>::infinity();
    for (const Step& st : steps_) {
      x = st.pow2 ? sat::mul_pow2(x, st.count) : sat::mul_int_pow(x, st.factor, st.count);
    }
    return x;
  }

 private:
  struct Step {
    std::uint64_t factor;  ///< integer factor (unused for doublings)
    unsigned count;        ///< repetitions, or number of doublings
    bool pow2;
  };
  std::vector<Step> steps_;
  bool saturated_ = false;
};

/// The injected constants of the assembly.
struct Rule {
  /// Step 3 (l.1206): a block is aligned iff align(s_U) >= sigma_U. Paper: N^10.
  LengthScale align;
  /// Step 4 (l.1216): the support of a named key sigma ends at support_end(sigma). Paper: N^16.
  LengthScale support_end;
  /// NOT in the paper (deviation D-11, Q-39; Mode B only): skip steps 1-8 and merge all join
  /// records of one weight simultaneously at height = that weight (c0 = 1). `align` and
  /// `support_end` are then unused. Mode A never sets it.
  bool exact_height = false;
  /// Prefix of the exception messages.
  std::string context = "assemble_local_hierarchy";

  /// The paper's constants for N vertices: c = 10 (kLocalAlignExponent) and N^16
  /// (kLocalSupportExponent), with the grid.
  [[nodiscard]] static Rule paper(std::size_t n) {
    Rule r;
    r.align = LengthScale::int_power(n, kLocalAlignExponent);
    r.support_end = LengthScale::int_power(n, kLocalSupportExponent);
    return r;
  }
};

// ------------------------------------------------------------------ sorting and searching

/// Stable counting sort of item indices 0..n-1 by key(i) in [0, k) (Lemma 3.1). Returns the
/// bucket offsets (size k + 1) and the sorted items.
template <class Key>
std::pair<std::vector<std::size_t>, std::vector<std::uint32_t>> bucket(std::size_t n, std::size_t k,
                                                                       Key key) {
  std::vector<std::size_t> begin(k + 1, 0);
  for (std::size_t i = 0; i < n; ++i) ++begin[key(i) + 1];
  for (std::size_t b = 0; b < k; ++b) begin[b + 1] += begin[b];
  std::vector<std::uint32_t> items(n);
  std::vector<std::size_t> pos(begin.begin(), begin.end() - 1);
  for (std::size_t i = 0; i < n; ++i) items[pos[key(i)]++] = static_cast<std::uint32_t>(i);
  return {std::move(begin), std::move(items)};
}

/// Index of the largest point <= w in pts[lo, hi), or hi if none (binary search).
template <class W>
std::size_t floor_index(const std::vector<W>& pts, std::size_t lo, std::size_t hi, W w) {
  const auto it = std::upper_bound(pts.begin() + static_cast<std::ptrdiff_t>(lo),
                                   pts.begin() + static_cast<std::ptrdiff_t>(hi), w);
  const auto i = static_cast<std::size_t>(it - pts.begin());
  return i == lo ? hi : i - 1;
}

// ------------------------------------------------------------------ alignment (step 3)

/// Step 3 for the blocks of one named key (l.1214-1218): successive medians (BFPRT, Lemma 3.1)
/// with one test align(m) >= sigma per median m. If it holds, every s >= m is aligned,
/// otherwise every s <= m is unaligned. Returns 1 for the aligned entries of `s`, 0 otherwise;
/// the result equals the per-entry test because the scaling is monotone.
template <class W>
std::vector<std::uint8_t> aligned_by_medians(std::span<const W> s, W sigma,
                                             const LengthScale& align);

// ------------------------------------------------------------------ support grid (step 4)

/// The doubling grid over the merged supports and the supports with their index ranges.
/// Class i is [grid[i], grid_upper[i]) with grid_upper[i] <= 2 * grid[i].
template <class W>
struct SupportGrid {
  std::vector<W> grid;                   ///< class lower boundaries, increasing
  std::vector<W> grid_upper;             ///< matching class upper boundaries
  std::vector<GridSupport<W>> supports;  ///< per key, in the order given

  /// The class of a value labelled with `key`: binary search between the key's indices.
  /// nullopt if w lies in no class of that range.
  [[nodiscard]] std::optional<std::size_t> classify(W w, std::uint32_t key) const {
    const auto& sp = supports[key];
    const std::size_t i = floor_index(grid, sp.lower_index, sp.upper_index, w);
    if (i == sp.upper_index || !(w < grid_upper[i])) return std::nullopt;
    return i;
  }
};

/// The support [unit, end(unit)] of each key (l.1216), indices not yet set. `units` holds the
/// grid unit of each key: the named corridor weight sigma in Sec 5.4, the name weight a in
/// Sec 8.2. `end` is N^16 in Sec 5.4 and octaves(3f) in Sec 8.2. The caller may lower
/// `lower` afterwards (aligned blocks, l.1219) before build_support_grid.
template <class W>
std::vector<GridSupport<W>> key_supports(std::span<const W> units, const LengthScale& end);

/// Step 4 (l.1214-1224): sorts the supports by lower end (ties by key), merges overlapping or
/// touching ones into components and starts a doubling grid at the minimum of each component.
/// The last class of a component stops at the start of the next component (l.1279). Each key
/// gets lower_index (largest grid point <= lower) and upper_index (first grid point > upper,
/// grid.size() if none). [Supplement] No grid point above `cap` is generated (pass the largest
/// value that will be classified, or infinity for none); the grid also stops at the first
/// doubling that saturates. Requires 0 < lower <= upper for every support.
template <class W>
SupportGrid<W> build_support_grid(std::vector<GridSupport<W>> supports, W cap);

/// Points unit, 2 unit, 4 unit, ... of a separate grid (l.1209-1213): the last point p is the
/// first with top < 2p, so the classes [p, 2p) cover [unit, top]. Stops at a saturated doubling.
template <class W>
std::vector<W> doubling_points(W unit, W top);

// ------------------------------------------------------------------ merges (steps 6-8)

/// Class-wise simultaneous merges (Sec 3.3, Lemma 3.6). The disjoint sets range over join
/// vertices 0..k-1; a vertex either carries a hierarchy node (a state's leaf, an atom's leaf)
/// or kNoNode (auxiliary pivot vertex). After a class is merged, every new component with >= 2
/// node-carrying parts becomes a node via add_node; a component with one such part keeps that
/// node (unary compression) and one without any stays empty (discarded, l.1281).
///
/// Endpoints: j -> std::pair of the two join vertices of record j.
/// AddNode:   (W height, std::vector<NodeId> children, W upper, bool inclusive, NodeOrigin)
///            -> NodeId of the node it created.
/// Distinct join vertices must carry distinct nodes.
template <class Endpoints, class AddNode>
class ClassMerger {
 public:
  ClassMerger(std::vector<NodeId> node_of_vertex, Endpoints endpoints, AddNode add_node)
      : endpoints_(std::move(endpoints)),
        add_node_(std::move(add_node)),
        dsu_(node_of_vertex.size()),
        rep_(std::move(node_of_vertex)),
        stamp_of_(rep_.size(), 0),
        group_of_(rep_.size(), 0) {}

  /// One class: the records `js` merged simultaneously; new nodes get the given height, class
  /// upper bound (inclusive or not) and origin.
  template <class W>
  void merge_class(std::span<const std::uint32_t> js, W height, W upper, bool inclusive,
                   NodeOrigin origin) {
    if (js.empty()) return;
    ++stamp_;
    pre_.clear();
    for (const std::uint32_t j : js) {
      const auto [a, b] = endpoints_(j);
      for (const std::uint32_t v : {a, b}) {
        const auto r = dsu_.find(v);
        if (stamp_of_[r] != stamp_) {
          stamp_of_[r] = stamp_;
          pre_.push_back(r);
        }
      }
    }
    for (const std::uint32_t j : js) {
      const auto [a, b] = endpoints_(j);
      dsu_.unite(a, b);
    }
    ++stamp_;
    group_root_.clear();
    gid_.assign(pre_.size(), 0);
    for (std::size_t i = 0; i < pre_.size(); ++i) {
      const auto r = dsu_.find(pre_[i]);
      if (stamp_of_[r] != stamp_) {
        stamp_of_[r] = stamp_;
        group_of_[r] = static_cast<std::uint32_t>(group_root_.size());
        group_root_.push_back(r);
      }
      gid_[i] = group_of_[r];
    }
    const std::size_t ng = group_root_.size();
    std::vector<std::size_t> begin(ng + 1, 0);
    for (std::size_t i = 0; i < pre_.size(); ++i) {
      if (rep_[pre_[i]] != kNoNode) ++begin[gid_[i] + 1];
    }
    for (std::size_t g = 0; g < ng; ++g) begin[g + 1] += begin[g];
    std::vector<NodeId> parts(begin[ng]);
    {
      std::vector<std::size_t> pos(begin.begin(), begin.end() - 1);
      for (std::size_t i = 0; i < pre_.size(); ++i) {
        if (rep_[pre_[i]] != kNoNode) parts[pos[gid_[i]]++] = rep_[pre_[i]];
      }
    }
    for (std::size_t g = 0; g < ng; ++g) {
      const std::size_t k = begin[g + 1] - begin[g];
      NodeId r = kNoNode;
      if (k == 1) {
        r = parts[begin[g]];
      } else if (k >= 2) {
        std::vector<NodeId> children(parts.begin() + static_cast<std::ptrdiff_t>(begin[g]),
                                     parts.begin() + static_cast<std::ptrdiff_t>(begin[g + 1]));
        r = add_node_(height, std::move(children), upper, inclusive, origin);
      }
      rep_[group_root_[g]] = r;
    }
  }

  [[nodiscard]] const DisjointSets& dsu() const noexcept { return dsu_; }

 private:
  Endpoints endpoints_;
  AddNode add_node_;
  DisjointSets dsu_;
  std::vector<NodeId> rep_;  ///< per disjoint-set root: the node of its component, or kNoNode
  std::vector<std::uint64_t> stamp_of_;
  std::vector<std::uint32_t> group_of_;
  std::uint64_t stamp_ = 0;
  std::vector<DisjointSets::Item> pre_;
  std::vector<DisjointSets::Item> group_root_;
  std::vector<std::uint32_t> gid_;
};

// ------------------------------------------------------------------ links and order (9-10)

struct LinkOrder {
  std::vector<NodeId> roots;  ///< parentless nodes in condensation order (l.1280)
  std::size_t links = 0;      ///< links created
};

/// Steps 9-10 (l.1281-1293) on the forest `nodes` whose internal nodes have their children but
/// no links yet. Every arc of `arcs` whose endpoints lie in different leaves is assigned to the
/// pair of children it joins at the lowest common ancestor of the two leaves (Lemma 3.9). At
/// each internal node v the assigned arcs shorter than node_upper[v] (at most it, if
/// node_upper_inclusive[v] != 0) form the child adjacency and a DFS takes one arc per
/// discovered child as link; the children are then put in a topological order of the assigned
/// arcs shorter than the node's height. `leaf_of_state[x]` is the leaf below which state x lies
/// (the leaf of x in Sec 5.4, the leaf of its part in a partition hierarchy). node_upper and
/// node_upper_inclusive are indexed by node id. Throws std::logic_error (prefixed with
/// `context`) if the children of a node are not connected by its arcs, or if the light arcs
/// between children or the arcs between roots form a cycle.
template <class W>
LinkOrder link_and_order(std::vector<HierarchyNode<W>>& nodes, std::span<const StateArc<W>> arcs,
                         std::span<const NodeId> leaf_of_state, std::span<const W> node_upper,
                         std::span<const std::uint8_t> node_upper_inclusive,
                         std::string_view context);

// ------------------------------------------------------------------ the local assembly

/// Throws std::invalid_argument ("<context>: inputs of different constructions") unless the
/// three inputs belong to one construction.
template <class W>
void check_inputs(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                  const ThresholdJoins<W>& joins, std::string_view context);

/// Sec 5.4 and the post-processing of Lemma 5.10 (steps 1-10 of paper/local/assemble.hpp) with
/// the constants of `rule`. With Rule::paper(N) this is the paper's assembly. The caller checks
/// the preconditions of its constants (Mode A: local_assembly_applies). Throws
/// std::invalid_argument via check_inputs and std::logic_error if a scale relation the
/// construction relies on fails. Counters (if `diag` is non-null): those listed at
/// paper::assemble_local_hierarchy; with rule.exact_height only local_assemble_calls,
/// local_assemble_blocks, local_assemble_exact_classes, local_assemble_nodes,
/// local_assemble_links and the dsu.* counters.
template <class W>
LocalAssembly<W> assemble(const OwnerDecomposition<W>& od, const RadialStateGraph<W>& rs,
                          const ThresholdJoins<W>& joins, const Rule& rule,
                          Diagnostics* diag = nullptr);

}  // namespace sssp::paper::assembly
