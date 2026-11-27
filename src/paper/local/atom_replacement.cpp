// Sec 5.5 (txt l.1302-1364): capping a local hierarchy (Lemma 5.13) and substituting atoms
// (replacement procedure l.1316-1319, Lemma 5.14-5.15). See atom_replacement.hpp.

#include "paper/local/atom_replacement.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "paper/local/saturating.hpp"
#include "paper/primitives/dsu.hpp"

namespace sssp::paper {

namespace {

constexpr std::uint32_t kNoPart = std::numeric_limits<std::uint32_t>::max();

template <class W>
bool positive_finite(W c) {
  using T = WeightTraits<W>;
  return T::is_valid_weight(c) && W{0} < c && c != T::infinity();
}

void count(Diagnostics* diag, const char* name, std::uint64_t by = 1) {
  if (diag != nullptr) SSSP_COUNT(diag->counters[name] += by);
  static_cast<void>(name);
  static_cast<void>(by);
}

/// Def 5.11: `tree` is a spanning tree of `g` (k - 1 edges, no cycle) and every edge is < 2c.
template <class W>
void check_spanning_tree(const CsrGraph<W>& g, std::span<const EdgeId> tree, W two_c) {
  const std::size_t k = g.num_vertices();
  if (k == 0) throw std::invalid_argument("cap_local_hierarchy: an atom has at least one state");
  if (tree.size() + 1 != k) {
    throw std::invalid_argument("cap_local_hierarchy: the spanning tree of H must have " +
                                std::to_string(k - 1) + " edges, got " +
                                std::to_string(tree.size()) + " (Def 5.11)");
  }
  DisjointSets dsu(k);
  for (const EdgeId id : tree) {
    if (id >= g.num_edges()) {
      throw std::invalid_argument("cap_local_hierarchy: spanning tree edge id out of range");
    }
    const auto& e = g.edge(id);
    if (!(e.weight < two_c)) {
      throw std::invalid_argument(
          "cap_local_hierarchy: a spanning tree edge is not shorter than 2c (Def 5.11)");
    }
    if (!dsu.unite(e.u, e.v)) {
      throw std::invalid_argument("cap_local_hierarchy: the given tree edges contain a cycle");
    }
  }
}

/// Lemma 5.13 on an Alg 4 result `local` of `internal`.
template <class W>
CappedLocalHierarchy<W> cap(LocalHierarchyResult<W> local, const CsrGraph<W>& internal,
                            std::span<const EdgeId> tree, W c, Diagnostics* diag) {
  const OrderedHierarchy<W>& in = local.hierarchy;
  const std::size_t num_nodes = in.nodes.size();
  const std::size_t num_states = in.num_states();
  // A node is "low" if it is a leaf or its height is strictly below c (l.1316).
  const auto low = [&](NodeId v) { return in.is_leaf(v) || in.nodes[v].height < c; };

  OrderedHierarchy<W> out = make_leaf_forest(in.graph, in.projection);
  out.c0 = in.c0;
  std::vector<NodeId> map(num_nodes, kNoNode);
  for (StateId x = 0; x < num_states; ++x) map[in.leaf_of_state[x]] = x;

  // l.1316-1317: the maximal low subtrees, including leaves not covered by any of them.
  std::vector<NodeId> retained;
  std::size_t removed = 0;
  for (NodeId v = 0; v < num_nodes; ++v) {
    if (!low(v)) {
      ++removed;
      continue;
    }
    const NodeId par = in.nodes[v].parent;
    if (par == kNoNode || !low(par)) retained.push_back(v);
  }

  // Copy every retained subtree (children before parents); remember the retained subtree of
  // every state for the connection step.
  std::vector<std::uint32_t> part(num_states, kNoPart);
  std::vector<std::pair<NodeId, bool>> stack;
  for (std::uint32_t i = 0; i < retained.size(); ++i) {
    stack.assign(1, {retained[i], false});
    while (!stack.empty()) {
      const auto [v, expanded] = stack.back();
      stack.pop_back();
      const auto& node = in.nodes[v];
      if (node.children.empty()) {
        part[node.state] = i;
        continue;
      }
      if (!expanded) {
        // H2 (heights nondecrease towards the root) makes every node below a low node low.
        if (!(node.height < c)) {
          throw std::logic_error(
              "cap_local_hierarchy: Alg 4 node of height >= c below a node "
              "of height < c (Def 4.1, H2)");
        }
        stack.push_back({v, true});
        for (const NodeId ch : node.children) stack.push_back({ch, false});
        continue;
      }
      std::vector<NodeId> children;
      children.reserve(node.children.size());
      for (const NodeId ch : node.children) children.push_back(map[ch]);
      map[v] = add_internal_node(out, node.height, std::move(children), node.links);
    }
  }

  CappedLocalHierarchy<W> res;
  res.removed_nodes = removed;
  if (retained.size() == 1) {
    // l.1335: a unique retained subtree already supplies the required root.
    res.root = map[retained[0]];
  } else {
    // l.1317: a common root of height c. Lemma 5.13 item 2 (l.1331-1334): zero-length
    // certificates join the copies of each old state, and each tree edge of H adds a
    // certificate between copies of its endpoints; a spanning tree of the quotient on the
    // retained subtrees is selected greedily. Each selected link is shorter than 2c.
    DisjointSets dsu(retained.size());
    std::vector<PhysicalLink<W>> links;
    std::vector<StateId> first(internal.num_vertices(), kNoState);
    for (StateId x = 0; x < num_states; ++x) {
      const Vertex v = in.projection[x];
      if (first[v] == kNoState) {
        first[v] = x;
      } else if (dsu.unite(part[first[v]], part[x])) {
        links.push_back({first[v], x, W{0}, kNoArc});
      }
    }
    for (const EdgeId id : tree) {
      const auto& e = internal.edge(id);
      const StateId x = local.native_entry[e.u];
      const StateId y = local.native_entry[e.v];
      if (dsu.unite(part[x], part[y])) links.push_back({x, y, e.weight, kNoArc});
    }
    if (dsu.num_sets() != 1) {
      throw std::logic_error(
          "cap_local_hierarchy: the certificates do not connect the retained "
          "subtrees (H is not connected)");
    }
    std::vector<NodeId> children;
    children.reserve(retained.size());
    for (const NodeId r : retained) children.push_back(map[r]);
    res.root = add_internal_node(out, c, std::move(children), std::move(links));
    out.c0 = std::max<std::uint32_t>(out.c0, 2);
    res.added_root = true;
    count(diag, "atom_cap_added_roots");
  }
  collect_roots(out);
  if (out.roots.size() != 1 || out.roots[0] != res.root) {
    throw std::logic_error("cap_local_hierarchy: the capped hierarchy does not have one root");
  }
  count(diag, "atom_cap_removed_nodes", removed);
  res.local.hierarchy = std::move(out);
  res.local.native_entry = std::move(local.native_entry);
  res.local.arc_edge = std::move(local.arc_edge);
  return res;
}

/// Structural checks of Def 5.12 (see the header). Returns the nodes in post-order (children
/// before parents), the atom of every state, the local index of every atom state and the leaf
/// of every atom.
struct PartitionLayout {
  std::vector<NodeId> post_order;
  std::vector<std::uint32_t> atom_of;
  std::vector<Vertex> local_of;
  std::vector<NodeId> atom_leaf;
};

template <class W>
PartitionLayout check_partition(const PartitionHierarchy<W>& ph) {
  const std::size_t n = ph.graph.num_states;
  const std::size_t num_nodes = ph.nodes.size();
  if (ph.projection.size() != n) {
    throw std::invalid_argument("substitute_atoms: projection size differs from the state count");
  }
  if (ph.leaf_atom.size() != num_nodes) {
    throw std::invalid_argument("substitute_atoms: leaf_atom size differs from the node count");
  }
  PartitionLayout lay{{},
                      std::vector<std::uint32_t>(n, kNoAtom),
                      std::vector<Vertex>(n, kNoVertex),
                      std::vector<NodeId>(ph.atoms.size(), kNoNode)};
  for (std::uint32_t a = 0; a < ph.atoms.size(); ++a) {
    const auto& states = ph.atoms[a].states;
    if (states.empty()) throw std::invalid_argument("substitute_atoms: an atom has no state");
    for (std::size_t i = 0; i < states.size(); ++i) {
      const StateId x = states[i];
      if (x >= n) throw std::out_of_range("substitute_atoms: atom state out of range");
      if (lay.atom_of[x] != kNoAtom) {
        throw std::invalid_argument("substitute_atoms: state " + std::to_string(x) +
                                    " lies in two atoms");
      }
      lay.atom_of[x] = a;
      lay.local_of[x] = static_cast<Vertex>(i);
    }
  }
  std::vector<std::uint8_t> singleton(n, 0);
  std::vector<std::uint32_t> refs(num_nodes, 0);
  for (NodeId v = 0; v < num_nodes; ++v) {
    const auto& node = ph.nodes[v];
    if (node.parent != kNoNode && node.parent >= num_nodes) {
      throw std::out_of_range("substitute_atoms: parent id out of range");
    }
    for (const NodeId ch : node.children) {
      if (ch >= num_nodes) throw std::out_of_range("substitute_atoms: child id out of range");
      if (ph.nodes[ch].parent != v) {
        throw std::invalid_argument("substitute_atoms: parent pointer and child list disagree");
      }
      ++refs[ch];
    }
    if (!node.children.empty()) {
      if (ph.leaf_atom[v] != kNoAtom) {
        throw std::invalid_argument("substitute_atoms: an atom leaf has children");
      }
      continue;
    }
    if (ph.leaf_atom[v] != kNoAtom) {
      const std::uint32_t a = ph.leaf_atom[v];
      if (a >= ph.atoms.size()) throw std::out_of_range("substitute_atoms: atom id out of range");
      if (node.state != kNoState || lay.atom_leaf[a] != kNoNode) {
        throw std::invalid_argument("substitute_atoms: atom " + std::to_string(a) +
                                    " needs exactly one leaf without a state");
      }
      lay.atom_leaf[a] = v;
      continue;
    }
    const StateId x = node.state;
    if (x >= n) throw std::out_of_range("substitute_atoms: leaf state out of range");
    if (lay.atom_of[x] != kNoAtom || singleton[x] != 0) {
      throw std::invalid_argument("substitute_atoms: state " + std::to_string(x) +
                                  " lies in two parts");
    }
    singleton[x] = 1;
  }
  for (std::uint32_t a = 0; a < ph.atoms.size(); ++a) {
    if (lay.atom_leaf[a] == kNoNode) {
      throw std::invalid_argument("substitute_atoms: atom " + std::to_string(a) + " has no leaf");
    }
  }
  for (StateId x = 0; x < n; ++x) {
    if (lay.atom_of[x] == kNoAtom && singleton[x] == 0) {
      throw std::invalid_argument("substitute_atoms: state " + std::to_string(x) +
                                  " lies in no part");
    }
  }
  // Post-order from the roots; every node must be reached exactly once.
  std::vector<std::pair<NodeId, bool>> stack;
  for (NodeId r = 0; r < num_nodes; ++r) {
    if (refs[r] != (ph.nodes[r].parent == kNoNode ? 0U : 1U)) {
      throw std::invalid_argument("substitute_atoms: parent pointer and child list disagree");
    }
    if (ph.nodes[r].parent != kNoNode) continue;
    stack.assign(1, {r, false});
    while (!stack.empty()) {
      const auto [v, expanded] = stack.back();
      stack.pop_back();
      if (expanded || ph.nodes[v].children.empty()) {
        lay.post_order.push_back(v);
        continue;
      }
      stack.push_back({v, true});
      for (const NodeId ch : ph.nodes[v].children) stack.push_back({ch, false});
    }
  }
  if (lay.post_order.size() != num_nodes) {
    throw std::invalid_argument("substitute_atoms: parent pointers contain a cycle");
  }
  // Def 5.11 / 5.12: an atom with a parent has cutoff at most the parent's height.
  for (std::uint32_t a = 0; a < ph.atoms.size(); ++a) {
    const NodeId par = ph.nodes[lay.atom_leaf[a]].parent;
    if (par != kNoNode && ph.nodes[par].height < ph.atoms[a].cutoff) {
      throw std::invalid_argument("substitute_atoms: atom " + std::to_string(a) +
                                  " has a cutoff above its parent's height (Def 5.12)");
    }
  }
  // [Interpretation, proposed Q-53] Internal edges belong to H, not to D.
  for (const auto& arc : ph.graph.arcs) {
    if (arc.tail >= n || arc.head >= n) {
      throw std::out_of_range("substitute_atoms: arc endpoint out of range");
    }
    if (lay.atom_of[arc.tail] != kNoAtom && lay.atom_of[arc.tail] == lay.atom_of[arc.head]) {
      throw std::invalid_argument(
          "substitute_atoms: an arc of D joins two states of one atom; internal edges belong to "
          "the atom's internal graph H");
    }
  }
  return lay;
}

}  // namespace

template <class W>
CappedLocalHierarchy<W> cap_local_hierarchy(const CsrGraph<W>& internal,
                                            std::span<const EdgeId> spanning_tree, W cutoff,
                                            Diagnostics* diag, const LocalHierarchyParams& p) {
  count(diag, "atom_cap_calls");
  if (!positive_finite(cutoff)) {
    throw std::invalid_argument(
        "cap_local_hierarchy: the cutoff must be positive and finite "
        "(Def 5.11)");
  }
  check_spanning_tree(internal, spanning_tree, sat::add(cutoff, cutoff));
  // Replacement procedure l.1316: run Algorithm 4 on H.
  return cap(local_hierarchy(internal, diag, p), internal, spanning_tree, cutoff, diag);
}

template <class W>
AtomSubstitution<W> substitute_atoms(const PartitionHierarchy<W>& ph, Diagnostics* diag,
                                     const LocalHierarchyParams& p) {
  count(diag, "atom_replace_calls");
  const PartitionLayout lay = check_partition(ph);
  const std::size_t n = ph.graph.num_states;
  const std::size_t num_atoms = ph.atoms.size();

  // Lemma 5.15: one local construction per atom, capped by Lemma 5.13.
  std::vector<CappedLocalHierarchy<W>> capped;
  capped.reserve(num_atoms);
  for (const auto& atom : ph.atoms) {
    const auto internal = CsrGraph<W>::from_edges(atom.states.size(), atom.edges);
    capped.push_back(cap_local_hierarchy(internal, std::span<const EdgeId>(atom.spanning_tree),
                                         atom.cutoff, diag, p));
  }
  count(diag, "atom_replace_atoms", num_atoms);

  // New states: singletons first (in state order), then each atom's local states as a block.
  AtomSubstitution<W> out;
  std::vector<StateId> singleton_new(n, kNoState);
  for (StateId x = 0; x < n; ++x) {
    if (lay.atom_of[x] != kNoAtom) continue;
    singleton_new[x] = static_cast<StateId>(out.old_state.size());
    out.old_state.push_back(x);
  }
  std::vector<StateId> offset(num_atoms);
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    offset[a] = static_cast<StateId>(out.old_state.size());
    for (const Vertex v : capped[a].local.hierarchy.projection) {
      out.old_state.push_back(ph.atoms[a].states[v]);
    }
  }
  const std::size_t total = out.old_state.size();
  if (total >= static_cast<std::size_t>(kNoState)) {
    throw std::length_error("substitute_atoms: too many states");
  }
  count(diag, "atom_replace_states", total);
  // l.1318: the entry of an atom state is its native inward copy.
  out.entry.resize(n);
  for (StateId x = 0; x < n; ++x) {
    const std::uint32_t a = lay.atom_of[x];
    out.entry[x] =
        a == kNoAtom ? singleton_new[x] : offset[a] + capped[a].local.native_entry[lay.local_of[x]];
  }
  // l.1319, Def 5.12: every copy keeps its old state's physical image.
  std::vector<Vertex> projection(total);
  for (StateId t = 0; t < total; ++t) projection[t] = ph.projection[out.old_state[t]];
  // Copies of every old state (bucketed, in increasing new id).
  std::vector<std::size_t> copy_start(n + 1, 0);
  for (const StateId x : out.old_state) ++copy_start[x + 1];
  for (std::size_t x = 0; x < n; ++x) copy_start[x + 1] += copy_start[x];
  std::vector<StateId> copies(total);
  {
    std::vector<std::size_t> cursor(copy_start.begin(), copy_start.end() - 1);
    for (StateId t = 0; t < total; ++t) copies[cursor[out.old_state[t]]++] = t;
  }

  // Arcs: the local arcs of every atom, then the old outside arcs.
  StateGraph<W> graph;
  graph.num_states = total;
  std::vector<ArcId> arc_offset(num_atoms);
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    arc_offset[a] = static_cast<ArcId>(graph.arcs.size());
    for (const auto& arc : capped[a].local.hierarchy.graph.arcs) {
      graph.arcs.push_back({offset[a] + arc.tail, offset[a] + arc.head, arc.length});
      out.old_arc.push_back(kNoArc);
    }
  }
  // l.1318-1319 (Lemma 5.14): an old arc x -> y enters y's native copy and is installed at every
  // copy of x. at_entry[i] is the copy whose tail is x's entry (used for backing arcs of links).
  std::vector<ArcId> at_entry(ph.graph.arcs.size(), kNoArc);
  std::uint64_t exit_copies = 0;
  for (ArcId i = 0; i < ph.graph.arcs.size(); ++i) {
    const auto& arc = ph.graph.arcs[i];
    const StateId head = out.entry[arc.head];
    for (std::size_t j = copy_start[arc.tail]; j < copy_start[arc.tail + 1]; ++j) {
      const StateId t = copies[j];
      if (t == out.entry[arc.tail]) at_entry[i] = static_cast<ArcId>(graph.arcs.size());
      graph.arcs.push_back({t, head, arc.length});
      out.old_arc.push_back(i);
      ++exit_copies;
    }
  }
  if (graph.arcs.size() >= static_cast<std::size_t>(kNoArc)) {
    throw std::length_error("substitute_atoms: too many arcs");
  }
  count(diag, "atom_replace_exit_copies", exit_copies);

  OrderedHierarchy<W> h = make_leaf_forest(std::move(graph), std::move(projection));
  h.c0 = ph.c0;
  // The capped hierarchies (their internal nodes were created children first).
  out.atom_root.resize(num_atoms);
  for (std::uint32_t a = 0; a < num_atoms; ++a) {
    const auto& ch = capped[a].local.hierarchy;
    std::vector<NodeId> map(ch.nodes.size(), kNoNode);
    for (NodeId u = 0; u < ch.nodes.size(); ++u) {
      const auto& node = ch.nodes[u];
      if (node.children.empty()) {
        map[u] = offset[a] + node.state;
        continue;
      }
      std::vector<NodeId> children;
      children.reserve(node.children.size());
      for (const NodeId c : node.children) children.push_back(map[c]);
      std::vector<PhysicalLink<W>> links;
      links.reserve(node.links.size());
      for (const auto& l : node.links) {
        links.push_back({offset[a] + l.x, offset[a] + l.y, l.length,
                         l.arc == kNoArc ? kNoArc : arc_offset[a] + l.arc});
      }
      map[u] = add_internal_node(h, node.height, std::move(children), std::move(links));
    }
    out.atom_root[a] = map[capped[a].root];
    h.c0 = std::max(h.c0, ch.c0);
  }
  // The partition forest, node for node (l.1335-1337: replacement duplicates endpoints only
  // inside a single old child; links name the entries of their endpoint states).
  out.partition_node.assign(ph.nodes.size(), kNoNode);
  for (const NodeId v : lay.post_order) {
    const auto& node = ph.nodes[v];
    if (ph.leaf_atom[v] != kNoAtom) {
      out.partition_node[v] = out.atom_root[ph.leaf_atom[v]];
      continue;
    }
    if (node.children.empty()) {
      out.partition_node[v] = singleton_new[node.state];  // leaf of state t is node t
      continue;
    }
    std::vector<NodeId> children;
    children.reserve(node.children.size());
    for (const NodeId c : node.children) children.push_back(out.partition_node[c]);
    std::vector<PhysicalLink<W>> links;
    links.reserve(node.links.size());
    for (const auto& l : node.links) {
      if (l.x >= n || l.y >= n)
        throw std::out_of_range("substitute_atoms: link state out of range");
      if (l.arc != kNoArc && l.arc >= ph.graph.arcs.size()) {
        throw std::out_of_range("substitute_atoms: link backing arc out of range");
      }
      links.push_back(
          {out.entry[l.x], out.entry[l.y], l.length, l.arc == kNoArc ? kNoArc : at_entry[l.arc]});
    }
    out.partition_node[v] =
        add_internal_node(h, node.height, std::move(children), std::move(links));
  }
  collect_roots(h);
  out.hierarchy = std::move(h);
  return out;
}

template <class W>
std::vector<StateOffer<W>> redirect_offers(const AtomSubstitution<W>& s,
                                           std::span<const StateOffer<W>> offers) {
  std::vector<StateOffer<W>> out;
  out.reserve(offers.size());
  for (const auto& o : offers) {
    if (o.state >= s.entry.size()) {
      throw std::out_of_range("redirect_offers: offer at a state that is not an underlying state");
    }
    out.push_back({s.entry[o.state], o.value});
  }
  return out;
}

template CappedLocalHierarchy<std::uint64_t> cap_local_hierarchy<std::uint64_t>(
    const CsrGraph<std::uint64_t>&, std::span<const EdgeId>, std::uint64_t, Diagnostics*,
    const LocalHierarchyParams&);
template CappedLocalHierarchy<double> cap_local_hierarchy<double>(const CsrGraph<double>&,
                                                                  std::span<const EdgeId>, double,
                                                                  Diagnostics*,
                                                                  const LocalHierarchyParams&);
template AtomSubstitution<std::uint64_t> substitute_atoms<std::uint64_t>(
    const PartitionHierarchy<std::uint64_t>&, Diagnostics*, const LocalHierarchyParams&);
template AtomSubstitution<double> substitute_atoms<double>(const PartitionHierarchy<double>&,
                                                           Diagnostics*,
                                                           const LocalHierarchyParams&);
template std::vector<StateOffer<std::uint64_t>> redirect_offers<std::uint64_t>(
    const AtomSubstitution<std::uint64_t>&, std::span<const StateOffer<std::uint64_t>>);
template std::vector<StateOffer<double>> redirect_offers<double>(
    const AtomSubstitution<double>&, std::span<const StateOffer<double>>);

}  // namespace sssp::paper
