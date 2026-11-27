#pragma once
// Ordered hierarchies (Sec 4, Definition 4.1, l.646-658) and seeded ordered hierarchies
// (Definition 4.2, l.659-672), plus an independent checker for every condition.
//
// Conditions of Definition 4.1, numbered as in the project's research notes:
//   H1  every state arc x->y of length l projects to a physical walk pi(x) ~> pi(y) of
//       length <= l;
//   H2  every internal node has >= 2 children and a positive height, nondecreasing towards
//       the root;
//   H3  at every internal node C, the arcs between distinct children of length < h(C)
//       ("light arcs") form an acyclic quotient on the children;
//   H4  the quotient of the roots on ALL inter-root arcs is acyclic;
//   H5  the children of every internal node C are connected by a tree of links; a link names
//       two endpoint states (in distinct children) and certifies a physical walk between
//       their images of length <= c0 * h(C), traversable in either direction;
//   H6  a link used solely as a certificate never becomes a relaxed arc.
//
// Design (read this before extending the types)
// ----------------------------------------------
// * The types below are the *input contract* of Sec 4 and are deliberately plain data:
//   producers (Alg 4 base case, Sec 5.4 assembly, Sec 5.5 atom replacement, test generators)
//   fill them, consumers (Sec 4.1 normalisation, Sec 4.2 child order, Alg 3 Evaluate) read
//   them. Nothing here depends on how a hierarchy was built.
// * Nodes are addressed by dense NodeId indices into OrderedHierarchy::nodes. Per-node data
//   added by later stages (Sec 4.1: unit nu(C), mass S(C), contraction map; Sec 4.2: rank in
//   parent, representative state rep(C), gap g(i,j)) belongs in PARALLEL ARRAYS indexed by
//   NodeId inside the later stage's own type (e.g. a NormalizedHierarchy holding a contracted
//   OrderedHierarchy plus `std::vector<W> unit, mass; std::vector<std::uint32_t> rank; ...`).
//   Do not add those fields to HierarchyNode: normalisation contracts nodes (Lemma 4.3), so its
//   node set differs from the input's, and keeping HierarchyNode minimal lets the same checker
//   validate the contracted forest as well.
// * H6 is enforced by type separation: evaluators relax only `graph.arcs`; links live in
//   HierarchyNode::links (a different type) and must never be inserted into `graph`. A link
//   may record the state arc it was taken from (`PhysicalLink::arc`, e.g. Lemma 5.10's DFS
//   discovery edges); kNoArc marks a certificate-only link (e.g. Lemma 5.13's zero-length
//   copy certificates). The checkable part of H6 is that a claimed backing arc exists, joins
//   the link's endpoints (either orientation) and is no longer than the link, so a consumer
//   may treat exactly those links as relaxable. [Interpretation, see validate_* below.]
// * Lengths are opaque reals (Sec 2, l.292-311): the code only copies, compares and adds W.
//   c0 is an absolute integer constant; c0 * h is formed by c0 - 1 additions.
// * Template definitions live in ordered_hierarchy.cpp and are explicitly instantiated for
//   std::uint64_t and double (the project's weight types).

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::paper {

using StateId = std::uint32_t;
using NodeId = std::uint32_t;
using ArcId = std::uint32_t;

inline constexpr StateId kNoState = std::numeric_limits<StateId>::max();
inline constexpr NodeId kNoNode = std::numeric_limits<NodeId>::max();
inline constexpr ArcId kNoArc = std::numeric_limits<ArcId>::max();
inline constexpr std::size_t kNoIndex = std::numeric_limits<std::size_t>::max();

/// Directed state arc tail -> head with a nonnegative length (Def 4.1, "finite directed graph D
/// with nonnegative arc lengths").
template <class W>
struct StateArc {
  StateId tail;
  StateId head;
  W length;
};

/// The state graph D. States are 0 .. num_states-1. Parallel arcs and self-loops are allowed
/// (normalisation drops self-loops, Sec 4.1 l.719).
template <class W>
struct StateGraph {
  std::size_t num_states = 0;
  std::vector<StateArc<W>> arcs;
};

/// A link of H5: endpoint states x, y and a certified upper bound `length` on a physical walk
/// between pi(x) and pi(y) (either direction). `arc` optionally names the state arc of D the
/// link was taken from; kNoArc means certificate-only (H6: never relaxed).
template <class W>
struct PhysicalLink {
  StateId x;
  StateId y;
  W length;
  ArcId arc = kNoArc;
};

/// One node of the rooted forest H. A leaf carries its state; an internal node carries its
/// children, its height h(C) and the links that connect its children (H5).
template <class W>
struct HierarchyNode {
  NodeId parent = kNoNode;
  std::vector<NodeId> children;        ///< empty for leaves
  StateId state = kNoState;            ///< leaves only
  W height{};                          ///< internal nodes only; must be > 0 (H2)
  std::vector<PhysicalLink<W>> links;  ///< internal nodes only; children - 1 links (H5)
};

/// Definition 4.1. `projection[x]` is pi(x), a vertex of the physical graph passed to the
/// checker. `leaf_of_state[x]` is the leaf node of state x. `roots` lists every parentless
/// node exactly once (use collect_roots after building). `c0` is the producer's declared
/// connection constant (1 for the Alg 4 base case and exact merges, <= 2 for Alg 4).
template <class W>
struct OrderedHierarchy {
  StateGraph<W> graph;
  std::vector<Vertex> projection;
  std::vector<HierarchyNode<W>> nodes;
  std::vector<NodeId> leaf_of_state;
  std::vector<NodeId> roots;
  std::uint32_t c0 = 1;

  [[nodiscard]] std::size_t num_states() const noexcept { return graph.num_states; }
  [[nodiscard]] bool is_leaf(NodeId v) const noexcept { return nodes[v].children.empty(); }
};

/// Weighted physical source (v_i, c_i) of Definition 4.2. An ordinary source s is (s, 0).
template <class W>
struct PhysicalSource {
  Vertex vertex;
  W offset;
};

/// Initial state offer (Def 4.2): a finite value; +infinity would mean "no offer" and is not
/// stored.
template <class W>
struct StateOffer {
  StateId state;
  W value;
};

/// Definition 4.2. Target distance delta(v) = min_i (c_i + d(v_i, v)) in the physical graph.
/// Every initial offer must be feasible (value >= delta(pi(x))). The requirement that each
/// requested target has a protected state is an input obligation the checker does not test.
template <class W>
struct SeededOrderedHierarchy {
  OrderedHierarchy<W> hierarchy;
  std::vector<PhysicalSource<W>> sources;
  std::vector<StateOffer<W>> offers;
};

// ------------------------------------------------------------------ builders

/// One leaf per state: node x is the leaf of state x, and every node is a root until
/// internal nodes are added. Throws std::invalid_argument if projection.size() !=
/// graph.num_states.
template <class W>
OrderedHierarchy<W> make_leaf_forest(StateGraph<W> graph, std::vector<Vertex> projection);

/// Appends an internal node with the given height, children and links, sets the children's
/// parent pointers and returns its id. Does not touch `roots` (call collect_roots at the end).
/// Throws std::invalid_argument if a child id is out of range or already has a parent.
template <class W>
NodeId add_internal_node(OrderedHierarchy<W>& h, W height, std::vector<NodeId> children,
                         std::vector<PhysicalLink<W>> links);

/// Sets h.roots to all parentless nodes in increasing id order.
template <class W>
void collect_roots(OrderedHierarchy<W>& h);

// ------------------------------------------------------------------ checker

/// Which clause of Def 4.1 / 4.2 a violation breaks. kStructure: the data does not describe a
/// forest whose leaves are exactly the states (the paper's conditions are then undefined).
enum class Condition : std::uint8_t { kStructure, kH1, kH2, kH3, kH4, kH5, kH6, kSeed };

enum class ViolationKind : std::uint8_t {
  // kStructure
  kStateCountMismatch,     ///< projection / leaf_of_state size differs from num_states
  kProjectionOutOfRange,   ///< pi(x) is not a vertex of the physical graph
  kArcEndpointOutOfRange,  ///< arc tail/head is not a state
  kNodeIndexOutOfRange,    ///< parent/child/root/leaf id is not a node
  kLeafMismatch,           ///< leaves are not exactly the states (one leaf per state)
  kParentChildMismatch,    ///< parent pointer and child list disagree
  kParentCycle,            ///< parent pointers do not reach a root
  kRootsMismatch,          ///< roots is not exactly the set of parentless nodes
  // kH1
  kArcLengthInvalid,  ///< negative, NaN or infinite arc length
  kArcWalkTooLong,    ///< d_phys(pi(x), pi(y)) > arc length
  // kH2
  kTooFewChildren,     ///< internal node with fewer than two children
  kNonPositiveHeight,  ///< height <= 0, NaN or infinite
  kHeightDecreases,    ///< internal child higher than its parent (reported at the child)
  // kH3
  kLightArcCycle,  ///< light-arc quotient on the children of `node` has a cycle
  // kH4
  kRootArcCycle,  ///< inter-root quotient has a cycle
  // kH5
  kLinkLengthInvalid,    ///< negative, NaN or infinite link length
  kLinkCount,            ///< number of links != number of children - 1
  kLinkEndpointOutside,  ///< a link endpoint is not a state below the node
  kLinkSameChild,        ///< both endpoints below the same child
  kLinkCycle,            ///< link joins two children already connected by earlier links
  kLinksNotSpanning,     ///< links leave the children disconnected
  kLinkTooLong,          ///< link length > c0 * h(C)
  kLinkWalkTooLong,      ///< d_phys(pi(x), pi(y)) > link length
  // kH6
  kLinkArcOutOfRange,  ///< claimed backing arc id is not an arc of D
  kLinkArcMismatch,    ///< backing arc does not join {x, y} or is longer than the link
  // kSeed (Def 4.2)
  kSourceOutOfRange,      ///< source vertex is not a physical vertex
  kSourceOffsetInvalid,   ///< source offset negative, NaN or infinite
  kOfferStateOutOfRange,  ///< offer names a non-existent state
  kOfferNotFinite,        ///< offer value is infinite/NaN (infinity means "no offer")
  kOfferInfeasible,       ///< offer value < delta(pi(x))
};

[[nodiscard]] Condition condition_of(ViolationKind kind) noexcept;
[[nodiscard]] const char* to_string(ViolationKind kind) noexcept;
[[nodiscard]] const char* to_string(Condition c) noexcept;

/// One violation with its location. Unused location fields hold kNo*.
struct Violation {
  ViolationKind kind;
  NodeId node = kNoNode;        ///< internal node (H2, H3, H5, H6) or offending node
  StateId state = kNoState;     ///< offending state (structure, seeds)
  ArcId arc = kNoArc;           ///< offending arc of D (H1, H3/H4 witness, structure)
  std::size_t link = kNoIndex;  ///< index into nodes[node].links (H5, H6) / offers (seeds)
  std::string message;

  [[nodiscard]] Condition condition() const noexcept { return condition_of(kind); }
};

struct ValidationReport {
  std::vector<Violation> violations;

  [[nodiscard]] bool ok() const noexcept { return violations.empty(); }
  [[nodiscard]] bool has(Condition c) const noexcept;
  [[nodiscard]] bool has(ViolationKind k) const noexcept;
  [[nodiscard]] std::size_t count(Condition c) const noexcept;
  /// One line per violation, for test diagnostics.
  [[nodiscard]] std::string summary() const;
};

/// Checks Definition 4.1 for `h` against the undirected physical graph `physical`, with the
/// connection constant h.c0. Each condition H1-H6 is checked independently and every
/// violation found is reported (not just the first). If the structure is malformed, H1-H6 are
/// skipped (they are undefined on a non-forest) and only kStructure violations are returned.
///
/// Physical walk bounds (H1, H5) are decided with an exact Dijkstra on `physical` (one search
/// per distinct source vertex; a direct edge or equal images short-circuit). This is a
/// checker, not part of the paper's algorithm. For double weights the Dijkstra sums are
/// rounded in path order, so a certificate computed by a different summation order can be
/// misjudged by one ulp; producers in this layer certify single edges, which is exact.
template <class W>
ValidationReport validate_ordered_hierarchy(const OrderedHierarchy<W>& h,
                                            const CsrGraph<W>& physical);

/// validate_ordered_hierarchy plus Definition 4.2: sources in range with valid offsets, every
/// initial offer names a state, is finite, and is feasible (value >= delta(pi(x))).
template <class W>
ValidationReport validate_seeded_ordered_hierarchy(const SeededOrderedHierarchy<W>& s,
                                                   const CsrGraph<W>& physical);

}  // namespace sssp::paper
